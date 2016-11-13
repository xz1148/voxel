import numpy as np

def SphereMeshGen(radii, dl)
# generates the xyz coodinate at the center of each cube
    N_Side = np.ceil(radii / dl)  # number of segments at each side
    N_Cube = N_Side ** 3

    Side_Min = -radii + 0.5 * dl
    Side_Max = radii + 0.5 * dl

    Side = np.linspace(Side_Min, Side_Max, N_Side)


    xyz = np.zeros(N_cube, 3)

    n = 0

    for x in N_Side:
        for y in N_Side:
            for z in N_Side:
                xyz[n] = [x, y, z]
                n += 1
    return xyz


if __name__ == "__main__":
    radii = 0.5
    dl = 0.1
    print SphereMeshGen(radii, dl)
