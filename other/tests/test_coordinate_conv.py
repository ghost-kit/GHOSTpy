import matplotlib.pyplot as plt
import numpy as np

from prototype import inv_common as ih

lats = np.linspace(start=0, stop=90, num=91)
lons = np.linspace(start=0, stop=360, num=1, endpoint=False)

data = {}
data2 = {}
data3 = {}
error = {}
error2 = {}
for lon in lons:
    data[lon] = {}
    data2[lon] = {}
    error[lon] = {}
    data3[lon] = {}
    error2[lon] = {}
    for lat in lats:
        theta = ih.deg_to_rad(lat)
        phi = ih.deg_to_rad(lon)
        loc = ih.get_location_from_theta_phi_r(theta=theta, phi=phi, r=1.0)
        data[lon][lat] = loc
        rtp = ih.get_r_theta_phi_from_location(loc)
        data2[lon][lat] = rtp
        error[lon][lat] = np.array([1.0, theta, phi]) - np.array(rtp)
        data3[lon][lat] = ih.get_lat(loc)
        error2[lon][lat] = theta - data3[lon][lat]

        print "Given: r:{} theta:{} phi:{}".format(1.0, theta, phi)
        print "Returned: r:{} theta:{} phi:{}".format(data2[lon][lat][0], data2[lon][lat][1], data2[lon][lat][2])
        print "Error: r:{} theta:{} phi:{}".format(error[lon][lat][0], error[lon][lat][1], error[lon][lat][2])
        print "===="

    X,Y = ih.dict_to_x_y(error[lon])
    X2,Y2 = ih.dict_to_x_y(error2[lon])

    r=[a[0] for a in Y]
    t=[a[1] for a in Y]
    p=[a[2] for a in Y]

    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.plot(X,t, label="get_r_t_p()")
    ax.set_xlabel("$\\theta$")
    ax.set_ylabel("$\\theta$ Error")
    ax.plot(X2,Y2, label="get_lat()")
    ax.set_title("$\\theta$ Conversion Error")

    plt.legend()
    plt.show()




