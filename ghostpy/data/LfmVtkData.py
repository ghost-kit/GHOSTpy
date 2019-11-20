import numpy as np
import types as tp
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

import GpData as gpd
import ghostpy.algorithms.common as algc
import ghostpy.algorithms.convert as algx
import ghostpy.algorithms.geometry as algg

from scipy import spatial
from scipy.interpolate import LinearNDInterpolator as interp


class LfmVtkData(gpd.data):
    def __init__(self, filename=None, vector=None, name="vtkData"):
        """
        Initializes a GHOSTpy data object for a vtk data file
        :param filename: Name of the vtk file to read
        :param vector: Name of the vector to be read within the vtk data file
        :param name: common name for the reader. Default is "vtkData"  used to display the trace mode
        """
        assert isinstance(name, tp.StringType), "Must specify a string Name for the data Mode"
        self.name = name
        assert isinstance(filename, tp.StringType), "Filename must be a valid string"
        self.file = filename
        self.reader = self.vtk_xml_reader
        assert isinstance(vector, tp.StringType), "Vector Name must be a string"

        # Let us not save all of the intermediate work, as it may overwhelm our resources
        odata = self.reader.GetOutput()
        in_data = dsa.WrapDataObject(odata)
        points = in_data.Points
        bvalues = in_data.PointData.GetArray(vector)
        self.extents = self.reader.GetUpdateExtent()
        self.dims = np.array([0,0,0,0])
        self.dims[0] = self.extents[1]+1
        self.dims[1] = self.extents[3]+1
        self.dims[2] = self.extents[5]+1
        self.dims[3] = 3

        # print("Dims: {}".format(self.dims))

        self.gridX, self.gridY, self.gridZ = self.build_grid2(points, self.dims[:-1])
        # self.gridX = cv.cm_to_re(self.gridX)
        # self.gridY = cv.cm_to_re(self.gridY)
        # self.gridZ = cv.cm_to_re(self.gridZ)

        self.dataX, self.dataY, self.dataZ = self.build_grid2(bvalues, self.dims[:-1])

        blade_trees = []
        FLookup = []
        midY = int(self.dims[1]/2)  # Prevents calculating along the X axis
        midZ = int(self.dims[2]/2)  # Prevents calculating along the X axis
        for i in range(self.dims[0]):
            fanX = self.gridX[i, :, :]
            fanY = self.gridY[i, :, :]
            fanZ = self.gridZ[i, :, :]

            self.fanP = zip(fanX.ravel(), fanY.ravel(), fanZ.ravel())
            blade_trees.append(spatial.cKDTree(self.fanP, leafsize=1e6))

            posX = self.gridX[i, midY, midZ]
            posY = self.gridY[i, midY, midZ]
            posZ = self.gridZ[i, midY, midZ]

            FLookup.append(algx.__xyz_to_fan_angle__(xyz=[posX, posY, posZ]))

        self.FanLookup = np.array(FLookup)
        self.blades = np.array(blade_trees)

    def __str__(self):
        return "GHOSTpy {} Reader".format(self.name)

    @staticmethod
    def build_grid(data, dims):
        x = np.zeros(shape=dims)
        y = np.zeros_like(x)
        z = np.zeros_like(x)

        # print ("Shape of arrays: {}, {}, {}".format(np.shape(x), np.shape(y), np.shape(z)))

        count = 0
        for k in range(dims[2]):
            for j in range(dims[1]):
                for i in range(dims[0]):
                    # print ("index: {}".format(count))
                    # print(data[count])
                    x[i, j, k] = data[count][0]
                    y[i, j, k] = data[count][1]
                    z[i, j, k] = data[count][2]
                    count += 1

        assert False

    @staticmethod
    def build_grid2(data, dims):
        x, y, z = np.hsplit(data, 3)

        # The grids use Fortran ordering. Need to remember that.
        x = x.reshape(dims, order='F')
        y = y.reshape(dims, order='F')
        z = z.reshape(dims, order='F')

        return x, y, z

    @property
    def vtk_xml_reader(self):
        """
        Generates and returns the vtk reader object for the file associated with the data object
        :return: vtkObject for the reader
        """
        reader = vtk.vtkXMLStructuredGridReader()
        reader.SetFileName(self.file)
        reader.Update()
        assert isinstance(reader, vtk.vtkObject)
        return reader

    def get_name(self):
        return self.name

    def get_calc_boundary(self):
        return 3.0

    def get_trace_boundary(self):
        return 2.5

    def get_xyz(self, xyz):
        # print("Looking for: {} @ {} RE".format(xyz, algc.mag(xyz)))

        if algc.mag(xyz) < self.get_trace_boundary():
            return np.array([np.NaN, np.NaN, np.NaN])

        cell_coord = self.get_nhood(xyz)

        value = self.valuate_point(xyz, locA=cell_coord)
        if np.any(np.isnan(value)):
            return np.array([np.NaN, np.NaN, np.NaN])
        else:
            return self.valuate_point(xyz, locA=cell_coord)

    def get_rlp(self, rlp):
        """
        responds to list of rlp components and returns the values of data at the components
        :param rlp: array of (r,l,p) components to be probed
        :return: array of values coorosponding to the (r,l,p) components provided
        """
        xyz = algc.sphere_to_cart(r=rlp[0], lam=rlp[1], phi=rlp[2])
        return self.get_xyz(xyz=xyz)

    @staticmethod
    def find_nearest(array, value, axis=None):
        idx = (np.abs(value - array)).argmin(axis=axis)
        return idx

    def valuate_point(self, xyz, locA):

        # Need to validate the points.
        # This is to be done as follows:
        #   0) Does the point exist on the X axis?
        #   0A) Yes: interpolate and return coefficients
        #   0B) NO: Continue
        #   1) Are any points invalid?
        #   2) if YES: go to part 2.
        #   3) if NO: interpolate lines to the X,Y plane
        #   4) Interpolate X/Y plane to the X/Y lines
        #   5) Does the point exist on these lines?
        #   6) if YES: return the interpolation coefficients
        #   7) if NO: go to part 3
        #
        # Part 2:
        #   1) Does the point exist on any of the lines?
        #   2) if YES: interpolate and return the coefficients
        #   3) if NO:  Is the point outside of the grid?
        #   4) if YES:  return None (invalid point)
        #   5) if NO:  SEARCH for boundary (reverse faces individually)
        #   6) When Found: Return Coefficients
        #
        # Part 3:
        #   1) modify original boxing and retry
        #
        p0, p1, p2, p3, p4, p5, p6, p7 = locA

        x,y,z = xyz

        # check to see if we are on the x axis
        if np.isclose(y, 0, rtol=0, atol=1e-16) and np.isclose(z, 0, rtol=0, atol=1e-16):
            # print("We are on the X-axis, and should just interpolate along the line.")
            a = self.__get_raw_point__(p0)
            b = self.__get_raw_point__(p3)
            # print ("A: {}, B: {}".format(a,b))

            w = (x-a[0])/(b[0] - a[0])

            va = self.__get_raw_data__(p0)
            vb = self.__get_raw_data__(p3)

            # print("va: {}\nvb: {}".format(va,vb))
            # print("P1: {}, P2: {}".format(p0,p3))

            value =  va + (vb - va) * w

            # print("W: {}".format(w))
            # print("test: {}".format(a+(b-a)*w))

            # print("Value at {} = {} ({})".format(xyz,value, algc.mag(value)))

            return value

        else:
            # print("We are NOT on the x-axis, so we must be do more...")

            l0 = self.__get_raw_point__(p0)
            v0 = self.__get_raw_data__(p0)

            nanArray = np.array([np.NaN, np.NaN, np.NaN])

            if p1[1] >= 0:
                l1 = self.__get_raw_point__(p1)
                v1 = self.__get_raw_data__(p1)

            else:
                print ("Checking plus values:")
                l1 = nanArray
                v1 = nanArray

            if p2[1] >= 0 and p2[2] >= 0:
                l2 = self.__get_raw_point__(p2)
                v2 = self.__get_raw_data__(p2)
            else:
                l2 = nanArray
                v2 = nanArray

            if p3[2] >= 0:
                l3 = self.__get_raw_point__(p3)
                v3 = self.__get_raw_data__(p3)
            else:
                l3 = nanArray
                v3 = nanArray

            l4 = self.__get_raw_point__(p4)
            v4 = self.__get_raw_data__(p4)

            if p5[1] >= 0:
                l5 = self.__get_raw_point__(p5)
                v5 = self.__get_raw_data__(p5)
            else:
                l5 = nanArray
                v5 = nanArray

            if p6[1] >= 0 and p6[2] >= 0:
                l6 = self.__get_raw_point__(p6)
                v6 = self.__get_raw_data__(p6)
            else:
                l6 = nanArray
                v6 = nanArray

            if p7[2] >= 0:
                l7 = self.__get_raw_point__(p7)
                v7 = self.__get_raw_data__(p7)
            else:
                l7 = nanArray
                v7 = nanArray

            labels = np.array(["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"])
            dat = np.array([l0, l1, l2, l3, l4, l5, l6, l7])

            all_good_points = True
            # print("Checking to see if there are any bad points")
            if np.allclose(l1, nanArray, equal_nan=True):
                # print("Bad Point at p1")
                all_good_points = False
            if np.allclose(l2, nanArray, equal_nan=True):
                # print("Bad Point at p2")
                all_good_points = False
            if np.allclose(l3, nanArray, equal_nan=True):
                # print("Bad Point at p3")
                all_good_points = False
            if np.allclose(l5, nanArray, equal_nan=True):
                # print("Bad Point at p5")
                all_good_points = False
            if np.allclose(l6, nanArray, equal_nan=True):
                # print("Bad Point at p6")
                all_good_points = False
            if np.allclose(l7, nanArray, equal_nan=True):
                # print("Bad Point at p7")
                all_good_points = False


        if not all_good_points:

            if p2[2] >= self.dims[2] or p3[2] >= self.dims[2] or p6[2] >= self.dims[2] or p7[2] >= self.dims[2]:
                # print ("Dim 3 Out of Bounds: {}.  Maximum index is {}.".format(p2[2], self.dims[2]-1))
                return nanArray

            # print("Probably out of bounds are on the inner boundary")
            # assert False
            return nanArray

        else:
            # The ordering of the points in these faces is critical. if the order is changed (not just rotated) the
            #   consistency of the normals will also change.  NORMALS must be consistently all pointing out or all pointing
            #   in so we can quickly ascertain weather a point is inside or outside of the cell.

            faces = []
            faces.append(np.array([l1, l2, l3, l0]))  # mated with face 3
            faces.append(np.array([l7, l3, l2, l6]))  # mated with face 4
            faces.append(np.array([l7, l6, l5, l4]))  # mated with face 1
            faces.append(np.array([l4, l5, l1, l0]))  # mated with face 2
            faces.append(np.array([l3, l7, l4, l0]))  # mated with face 6
            faces.append(np.array([l5, l6, l2, l1]))  # mated with face 5
            faces = np.array(faces)

            # NOTE: We will use the newell function so we can 'rotate' the face array if we fail to get a normal on
            #   a non - n6 face.  If we cycle through all points and don't find one that works, we have a problem.
            #   it is important to note that the ordering cannot change, or the direction of the normal will change.
            #   we need the normal directions to be self consistent (either all pointing in or all pointing out)
            n1 = algg.newell_surf_normal(faces[0])
            n2 = algg.newell_surf_normal(faces[1])
            n3 = algg.newell_surf_normal(faces[2])
            n4 = algg.newell_surf_normal(faces[3])
            n5 = algg.newell_surf_normal(faces[4])
            n6 = algg.newell_surf_normal(faces[5])

            # The use of two opposing points for the normal calculations allows us to
            #   check the entire cell with only two surface points
            v0 = algg.vector_between(xyz, l0)
            v1 = algg.vector_between(xyz, l6)

            # the directions fo these results determine if the point is inside of the cell or not.
            #   any distance that has a sign opposite any of the others has a point outside of the cell
            #   this is easy to prove via contradiction.
            d1 = np.dot(v0, n1)
            d2 = np.dot(v1, n2)
            d3 = np.dot(v1, n3)
            d4 = np.dot(v0, n4)
            d5 = np.dot(v0, n5)
            d6 = np.dot(v1, n6)

            distances = np.array([d1, d2, d3, d4, d5, d6])
            dir_signs = np.sign(distances)

        if np.isnan(np.array([d1, d2, d3, d4, d5, d6])).sum() > 1:
            print ("Too Many NaN's.... Need to handle this here")

            print("Point 1: {}: {} = {}".format(p0, l0, v0))
            print("Point 2: {}: {} = {}".format(p1, l1, v1))
            print("Point 3: {}: {} = {}".format(p2, l2, v2))
            print("Point 4: {}: {} = {}".format(p3, l3, v3))
            print("Point 5: {}: {} = {}".format(p4, l4, v4))
            print("Point 6: {}: {} = {}".format(p5, l5, v5))
            print("Point 7: {}: {} = {}".format(p6, l6, v6))
            print("Point 8: {}: {} = {}".format(p7, l7, v7))

            # self.PlotCell(d1, d2, d3, d4, d5, d6, dat, faces, l0, l1, l2, l3, l4, l5, l6, l7, labels, n1, n2, n3, n4, n5, n6, xyz, show=True)

            assert False

        zeros = dir_signs == 0.0
        negs = np.argwhere(dir_signs == -1)
        pos = np.argwhere(dir_signs == 1)

        if np.any(zeros):
            print ("One of the faces contains the point needed. {}".format(np.argwhere(zeros)))
            print ("Zeros: {}".format(zeros))
            print ("Number of Zeros: {}".format(zeros.size))
            assert False

        elif (len(negs) > 0 and len(pos) == 0) or (len(pos) > 0 and len(negs) == 0):
            # print ("Within Bounds")

            fun = interp(points=np.array([l0, l1, l2, l3, l4, l5, l6, l7]), values=np.array([v0, v1, v2, v3, v4, v5, v6, v7]))

            # print ("Value for {}: {}".format(xyz, fun(xyz)))
            return fun(xyz).flatten()

        else:
            # print ("Out of Bounds. Must re-center cell.")

            # self.PlotCell(d1, d2, d3, d4, d5, d6, dat, faces, l0, l1, l2, l3, l4, l5, l6, l7, labels, n1, n2, n3, n4, n5, n6, xyz)

            # Face 0/2 shift on oX
            # Face 1/3 shift on oZ
            # Face 5/4 shift on oY
            face_pair_idx = {0: 2, 1: 3, 2: 0, 3: 1, 4: 5, 5: 4}
            face_comp_idx = {0: 0, 1: 2, 2: 0, 3: 2, 4: 1, 5: 1}
            face_point_idx = {0: [0, 1, 2, 3], 1: [2, 3, 6, 7], 2: [4, 5, 6, 7], 3: [0, 1, 4, 5], 4: [0, 3, 4, 7], 5: [1, 2, 5, 6]}

            fpoints = []
            fpoints.append([p1, p2, p3, p0])
            fpoints.append([p7, p3, p2, p6])
            fpoints.append([p7, p6, p5, p4])
            fpoints.append([p4, p5, p1, p0])
            fpoints.append([p3, p7, p4, p0])
            fpoints.append([p5, p6, p2, p1])

            if len(negs) > len(pos):
                # work from the pos list
                wlist = pos.flatten()

            else:
                # work from the negs list
                wlist = negs.flatten()

            # print ("work face index list: {}".format(wlist))
            # print ("Face Pairs Index list: {}".format(face_pair_idx))

            for idx in wlist:
                if fpoints[idx][0][face_comp_idx[idx]] < fpoints[face_pair_idx[idx]][0][face_comp_idx[idx]]:
                    # print ("Must Move index {} down".format(face_comp_idx[idx]))
                    adjuster = -2
                else:
                    # print ("Must Move index {} up.".format(face_comp_idx[idx]))
                    adjuster = 2

                # print (locA)

                for idx in face_point_idx[wlist[0]]:
                    # print ("Index: {}".format(idx))
                    # print ("component: {}".format(face_comp_idx[wlist[0]]))
                    locA[idx][face_comp_idx[wlist[0]]] += adjuster

                    # print (locA)

            return self.valuate_point(xyz, locA)

    def PlotCell(self, d1, d2, d3, d4, d5, d6, dat, faces, l1, l2, l3, l4, l5, l6, l7, l8, labels, n1, n2, n3, n4, n5, n6, xyz, show=False):
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d import proj3d
        from matplotlib.lines import Line2D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(10, 10), dpi=90)
        ax = fig.add_subplot(111, projection='3d')
        ax.plot([l1[0], l4[0]], [l1[1], l4[1]], [l1[2], l4[2]], 'k-')
        ax.plot([l4[0], l3[0]], [l4[1], l3[1]], [l4[2], l3[2]], 'k*-.')
        ax.plot([l3[0], l2[0]], [l3[1], l2[1]], [l3[2], l2[2]], 'k-')
        ax.plot([l2[0], l1[0]], [l2[1], l1[1]], [l2[2], l1[2]], 'k*-.')
        ax.plot([l5[0], l8[0]], [l5[1], l8[1]], [l5[2], l8[2]], 'k-')
        ax.plot([l8[0], l7[0]], [l8[1], l7[1]], [l8[2], l7[2]], 'k*-.')
        ax.plot([l7[0], l6[0]], [l7[1], l6[1]], [l7[2], l6[2]], 'k-')
        ax.plot([l6[0], l5[0]], [l6[1], l5[1]], [l6[2], l5[2]], 'k*-.')
        ax.plot([l1[0], l5[0]], [l1[1], l5[1]], [l1[2], l5[2]], 'k-')
        ax.plot([l2[0], l6[0]], [l2[1], l6[1]], [l2[2], l6[2]], 'k--')
        ax.plot([l4[0], l8[0]], [l4[1], l8[1]], [l4[2], l8[2]], 'k-')
        ax.plot([l3[0], l7[0]], [l3[1], l7[1]], [l3[2], l7[2]], 'k-')
        ax.scatter(l1[0], l1[1], l1[2], c='r', s=25)
        ax.scatter(l2[0], l2[1], l2[2], c='g', s=25)
        ax.scatter(l3[0], l3[1], l3[2], c='y')
        ax.scatter(l4[0], l4[1], l4[2], c='b')
        ax.scatter(l5[0], l5[1], l5[2], c='k')
        ax.scatter(l6[0], l6[1], l6[2], c='k')
        ax.scatter(l7[0], l7[1], l7[2], c='k')
        ax.scatter(l8[0], l8[1], l8[2], c='k')
        ax.scatter(xyz[0], xyz[1], xyz[2], c='y', s=25)
        ax.set_xlabel("X-Axis")
        ax.set_ylabel("Y-Axis")
        ax.set_zlabel("Z-Axis")
        # ==================================================
        # Verification Plotting only between double lines.
        ln1x = [l1[0], xyz[0]]
        ln1y = [l1[1], xyz[1]]
        ln1z = [l1[2], xyz[2]]
        ax.plot(ln1x, ln1y, ln1z, 'k-')
        ln7x = [l7[0], xyz[0]]
        ln7y = [l7[1], xyz[1]]
        ln7z = [l7[2], xyz[2]]
        ax.plot(ln7x, ln7y, ln7z, 'b-')
        # ===================================================
        # =================================================================================
        # Verification plotting ONLY between double lines
        #
        #
        # print ("d1: {}".format(d1))
        # print ("d2: {}".format(d2))
        # print ("d3: {}".format(d3))
        # print ("d4: {}".format(d4))
        # print ("d5: {}".format(d5))
        # print ("d6: {}".format(d6))
        dn1 = n1 * d1
        dn2 = n2 * d2
        dn3 = n3 * d3
        dn4 = n4 * d4
        dn5 = n5 * d5
        dn6 = n6 * d6
        ax.scatter((xyz + dn1)[0], (xyz + dn1)[1], (xyz + dn1)[2], c='r', marker='*')
        ax.scatter((xyz + dn2)[0], (xyz + dn2)[1], (xyz + dn2)[2], c='r', marker=">")
        ax.scatter((xyz + dn3)[0], (xyz + dn3)[1], (xyz + dn3)[2], c='b', marker="*")
        ax.scatter((xyz + dn4)[0], (xyz + dn4)[1], (xyz + dn4)[2], c='b', marker="<")
        ax.scatter((xyz + dn5)[0], (xyz + dn5)[1], (xyz + dn5)[2], c='r', marker="^")
        ax.scatter((xyz + dn6)[0], (xyz + dn6)[1], (xyz + dn6)[2], c='b', marker="v")
        # print ("V1: {}".format(v1))
        # print ("V2: {}".format(v2))
        fcolors = ['r', 'g', 'b', 'y', 'k', 'w']
        count = 0
        face_plots = []
        face_legends = []
        for face in faces:
            x, y, z = np.hsplit(face, 3)
            x = x.flatten()
            y = y.flatten()
            z = z.flatten()
            verts = [zip(x, y, z)]
            ax.add_collection3d(Poly3DCollection(verts, alpha=0.5, color=fcolors[count]))
            face_plots.append(Line2D([0], [0], linestyle="none", marker="s", alpha=0.5, markersize=10, markerfacecolor=fcolors[count]))
            face_legends.append("Face {}".format(count + 1))
            count += 1
        ax.legend(face_plots, face_legends, numpoints=1, loc='best')
        # print("Requested Point: {}".format(xyz))
        # print("Labels: {}".format(labels))
        # print("Locations: {}".format(dat))
        # print("Length of Labels: {}".format(len(labels)))
        for a in range(len(labels)):
            x, y, z = dat[a]
            lab = ax.text(x, y, z, " p{}".format(str(a + 1)), size=10, zorder=1, color='r')

        # plt.figtext(0.02, 0.01, "Face Distances from point {}:\n1: {}\n2: {}\n3: {}\n4: {}\n5: {}\n6: {}".format(xyz, d1, d2, d3, d4, d5, d6), size=10, zorder=1, color='k')

        # =================================================================================
        fig.savefig("points/point_{}.pdf".format(xyz))
        if show:
            plt.show()
        plt.close(fig)

    def get_nhood(self, xyz):
        # TODO: implement a quick tri-linear interpolator, even if it isn't correct.. I need to get the rest of this working!
        # TODO: Current Errors seen: Not getting proper grid with on-axis locations (-1 locations on the fan are not working correctly)
        # TODO: Some points are being miss-identified on the fan number.

        # TODO: Need to figure out the relationship between direction and which way to move the points to capture the correct cell.
        # TODO: Need to make sure the dot is always inside the box.

        # print("Getting Cell: {}".format(xyz))

        pc = self.get_nearest_point(xyz)
        # print ("Closest Point: {}".format(pc))

        if np.all(np.isnan(pc)):
            return np.array([pc, pc, pc, pc, pc, pc, pc, pc])

        pcXYZ = self.__get_raw_point__(pc)
        # print("Closest Point to {}: {}".format(xyz, pcXYZ))


        # Rotate to axis plane
        alpha = self.FanLookup[pc[0]]

        pcRot = algx.__rotate_x__(pcXYZ, algx.__xyz_to_fan_angle__(pcXYZ))
        xyzRot = algx.__rotate_x__(xyz, algx.__xyz_to_fan_angle__(xyz))
        xyzAlph = algx.__xyz_to_fan_angle__(xyz)

        # print("Requested Alpha: {}".format(xyzAlph))
        # print("Closest Alpha: {}".format(alpha))

        pcPol = algx.cart_to_polar([pcRot[0], pcRot[2]])
        xyzPol = algx.cart_to_polar([xyzRot[0], xyzRot[2]])

        # print("Closest Point Polar Coordinates: {}".format(pcPol))
        # print("Request Point Polar Coordinates: {}".format(xyzPol))

        if alpha <= xyzAlph:
            oX = np.array([1, 0, 0])
        else:
            oX = np.array([-1, 0, 0])

        if pcPol[0] <= xyzPol[0]:
            oZ = np.array([0, 0, 1])
        else:
            oZ = np.array([0, 0, -1])

        if pcPol[1] == 0:
            pcPol[1] = 2 * np.pi

        if pcPol[1] < xyzPol[1]:
            oY = np.array([0, -1, 0])
        else:
            oY = np.array([0, 1, 0])

        p1 = pc
        p2 = pc + oY
        p3 = pc + oY + oZ
        p4 = pc + oZ
        p5 = pc + oX
        p6 = pc + oX + oY
        p7 = pc + oX + oY + oZ
        p8 = pc + oX + oZ

        # This fix is for data-sets that repeat the first fan at the end.
        #  This will likely need to be removed once I move to pyLTR data reading.
        if p5[0] == -1:
            p5[0] = -2
        if p6[0] == -1:
            p6[0] = -2
        if p7[0] == -1:
            p7[0] = -2
        if p8[0] == -1:
            p8[0] = -2

        pts = np.array([p1, p2, p3, p4, p5, p6, p7, p8])

        return pts

    def __get_raw_data__(self, xyz):
        # print ("XYZ: {}".format(xyz))
        if np.any(np.isnan(xyz)):
            return xyz
        if xyz[0] >= self.dims[0]:
            xyz[0] -= self.dims[0]
            # print("adjusted x to wrap: {}".format(xyz))
        if xyz[1] >= self.dims[1] or xyz[2] >= self.dims[2]:
            # print(" y or z out of bounds")
            return np.array([np.NaN, np.NaN, np.NaN])
        return np.array([self.dataX[tuple(xyz)], self.dataY[tuple(xyz)], self.dataZ[tuple(xyz)]])

    def __get_raw_point__(self, xyz):
        # print ("XYZ: {}".format(xyz))
        if np.any(np.isnan(xyz)):
            return xyz

        if xyz[0] >= self.dims[0]:
            xyz[0] -= self.dims[0]
            # print("adjusted x to wrap: {}".format(xyz))
        if xyz[1] >= self.dims[1] or xyz[2] >= self.dims[2]:
            # print(" y or z out of bounds")
            return np.array([np.NaN, np.NaN, np.NaN])
        return np.array([self.gridX[tuple(xyz)], self.gridY[tuple(xyz)], self.gridZ[tuple(xyz)]])

    def get_nearest_point(self, xyz):

        # print ("Getting Nearest Point for {}".format(xyz))
        alpha_req = algx.__xyz_to_fan_angle__(xyz)

        alpha_idx = self.find_nearest(self.FanLookup, alpha_req)
        # print ("Nearest {} = Index: {}, Value: {}".format(alpha_req, alpha_idx, self.FanLookup[alpha_idx]))

        cpd, cpi = self.blades[alpha_idx].query(xyz)

        # print("CPI: {}".format(cpi))
        if cpi >= len(self.fanP):
            # print ("CPI out of bounds")
            return np.array([np.NaN, np.NaN, np.NaN])

        # print ("Closest Point Distance: {} @ index {}".format(cpd, cpi))
        cploc = self.fanP[cpi]

        # print ("Closest Point for {}: {}".format(xyz, cploc))

        # unravel the indices
        YZidx = np.unravel_index(cpi, self.dims[1:-1])
        # print("YZ: {}".format(YZidx))
        pret = np.array([alpha_idx, YZidx[0], YZidx[1]])

        return pret

# TODO NOTES:
# TODO: 1) grid the data set into blocks of interpolation areas
# TODO: 2) store the interpolators after they are created, but create first when first called
# TODO: 3) experiment with optimal size blocks
# TODO: 4) use a lookup table to quickly find the correct block.

