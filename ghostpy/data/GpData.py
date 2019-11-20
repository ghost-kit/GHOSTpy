class data:
    def get_name(self):
        """
        Return the name of the data mode
        :return:
        """
        pass

    def get_xyz(self, xyz):
        """
        should respond to a request for data at position X,Y,Z
        :param xyz: Cartesian X coordinate
                    Cartesian Y coordinate
                    Cartesian Z coordinate
        :return: Value of the data at (X,Y,Z)
        """
        pass

    def get_rlp(self, rlp):
        """
        should respond to a request for data from position radius, latitude, longitude (r, lambda, phi)
        :param rlp: spherical radius component
                    spherical lambda component (latitude)
                    spherical phi component (longitude)

        :return: value of data at (r,l,p)
        """
        pass

    def set_inner_boundary(self, re=1.0):
        """
        set the inner boundary for the data
            this will tell the system where to calculate the intersection.
        :param re: distance from center of the Earth for drift shell calculations
        :return: None
        """
        pass

    def get_calc_boundary(self):
        pass

    def get_trace_boundary(self):
        pass

    def get_reader(self):
        pass

    def get_integrator(self):
        pass

    def get_bisection_model(self):
        pass
