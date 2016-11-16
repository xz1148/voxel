import Mesh.MeshGen as MeshGen
import Mesh.View as View
import numpy as np
from Greens.Greens import *
order = 10
radii = 0.05
dl = 0.025
omega = 1e6 * 2 * np.pi
xyz_sphere, dl_new = MeshGen.SphereMeshGen(radii, dl)
size = 10
View.ViewPoints(xyz_sphere, size)
N_p = xyz_sphere.shape[0] #the number of points
sample, weight = np.polynomial.legendre.leggauss(order)


Gxx = np.asmatrix(np.zeros((N_p, N_p), complex))
Gyx = np.asmatrix(np.zeros((N_p, N_p), complex))
Gzx = np.asmatrix(np.zeros((N_p, N_p), complex))
Gxy = np.asmatrix(np.zeros((N_p, N_p), complex))
Gyy = np.asmatrix(np.zeros((N_p, N_p), complex))
Gzy = np.asmatrix(np.zeros((N_p, N_p), complex))
Gxz = np.asmatrix(np.zeros((N_p, N_p), complex))
Gyz = np.asmatrix(np.zeros((N_p, N_p), complex))
Gzz = np.asmatrix(np.zeros((N_p, N_p), complex))




# this is the matrix shows the field generated by x direction of current
direction = np.array([1,0,0])
for m in range(N_p):
    for n in range(N_p):
        if m == n:
            Gxx[m,n] = LocalE(omega, dl_new, sample, weight)
            Gyx[m,n] = 0.0 + 0.0j
            Gzx[m,n] = 0.0 + 0.0j
        else:
            Gxx[m,n] = RemoteE(omega, xyz_sphere[m, :], xyz_sphere[n, :],
                               dl_new, sample, weight, direction)[0]
            Gyx[m,n] = RemoteE(omega, xyz_sphere[m, :], xyz_sphere[n, :],
                               dl_new, sample, weight, direction)[1]
            Gzx[m,n] = RemoteE(omega, xyz_sphere[m, :], xyz_sphere[n, :],
                               dl_new, sample, weight, direction)[2]
Gx = np.concatenate((Gxx, Gyx, Gzx), axis=1)
# this is the matrix shows the field generated by x direction of current
direction = np.array([0,1,0])
for m in range(N_p):
    for n in range(N_p):
        if m == n:
            Gyy[m,n] = LocalE(omega, dl_new, sample, weight)
            Gxy[m,n] = 0.0 + 0.0j
            Gzy[m,n] = 0.0 + 0.0j
        else:
            Gxy[m,n] = RemoteE(omega, xyz_sphere[m, :], xyz_sphere[n, :],
                               dl_new, sample, weight, direction)[0]
            Gyy[m,n] = RemoteE(omega, xyz_sphere[m, :], xyz_sphere[n, :],
                               dl_new, sample, weight, direction)[1]
            Gzy[m,n] = RemoteE(omega, xyz_sphere[m, :], xyz_sphere[n, :],
                               dl_new, sample, weight, direction)[2]
Gy = np.concatenate((Gxy, Gyy, Gzy), axis=1)


# this is the matrix shows the field generated by x direction of current
direction = np.array([0,0,1])
for m in range(N_p):
    for n in range(N_p):
        if m == n:
            Gzz[m,n] = LocalE(omega, dl_new, sample, weight)
            Gyz[m,n] = 0.0 + 0.0j
            Gxz[m,n] = 0.0 + 0.0j
        else:
            Gxz[m,n] = RemoteE(omega, xyz_sphere[m, :], xyz_sphere[n, :],
                               dl_new, sample, weight, direction)[0]
            Gyz[m,n] = RemoteE(omega, xyz_sphere[m, :], xyz_sphere[n, :],
                               dl_new, sample, weight, direction)[1]
            Gzz[m,n] = RemoteE(omega, xyz_sphere[m, :], xyz_sphere[n, :],
                               dl_new, sample, weight, direction)[2]
Gz = np.concatenate((Gxz, Gyz, Gzz), axis=1)
G = np.vstack((Gx, Gy, Gz))
print G


