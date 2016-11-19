import numpy as np
import View


def SphereMeshGen(radii, dl):
# generates the xyz coodinate at the center of each cube
    N_Side = int(np.ceil(2 * radii / dl)) + 1  # number of segments at each side
    Side_Min = -radii
    Side_Max = radii

    Side = np.linspace(Side_Min, Side_Max, N_Side)
    dl_new = Side[1] - Side[0]
#    Side = Side[0:-1] + 0.5 * dl_new
    N_Cube = (Side.shape[0])**3
    xyz = np.zeros((N_Cube, 3))

    n = 0

    for x in Side:
        for y in Side:
            for z in Side:
                xyz[n] = [x, y, z]
                n += 1

    greater_radii = np.zeros(N_Cube, bool)
    for n in range(N_Cube):
        if np.linalg.norm(xyz[n]) > radii:
            greater_radii[n] = True

    xyz_sphere = xyz[np.logical_not(greater_radii)]
    return xyz_sphere, dl_new


if __name__ == "__main__":
    radii = input("Please input the radii")
    dl = input("Please input the mesh size")
    xyz, dl_new = SphereMeshGen(radii, dl)
    View.ViewPoints(xyz, 10)
    print "The new step size is %3.2f" % dl_new

