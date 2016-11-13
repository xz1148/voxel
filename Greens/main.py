import numpy as np
import EMConst
from scipy.integrate import nquad

eps0 = EMConst.eps0
mu0 = EMConst.mu0
PI = np.pi

def sample_points1D(x1, x2, sample):
    a = (x2 + x1) / 2.0
    b = (x2 - x1) / 2.0
    samples = sample*b + a
    return samples


def r_xyz(x, y, z):
    # this function calculates the distance form (x, y, z) to (0,0,0)
    r = np.sqrt(x**2+y**2+z**2)
    return r

def dGdr(k0, x, y, z):
    # dG / dr
    r0 = r_xyz(x, y, z)
    dGdr_r = (-1/(4*PI)) * (-1j*k0*np.exp(-1j*k0*r0))/r0 - np.exp(-1j*k0*r0)/r0**2
    return dGdr_r

def LocalA(dl, omega):
    # this function calculates the A vector contributed by a unit J( current )
    # input:
    # dl: the edge length of voxel
    k0 = omega * sqrt(mu0 * eps0)
    r_max = (3.0/(4.0*PI)) ** (1.0/3.0) * dl
    C0 = -mu0 / (4*PI*1j*k0)
    int_A_r_max = C0 * (r_max*np.exp(-1j*k0*r_max)+np.exp(-1j*k0*r_max)/(1j*k0))
    int_A_0 = C0 * (1/(1j*k0))
    A = int_A_r_max - int_A_0
    return A


def GradPhiIntegrand(dl, omega):
    # this funtion returns a function which used as integrand
    halfz = dl / 2
    k0 = omega * np.sqrt(mu0 * eps0)
    Charges = 1 / (1j * omega)  # equals to -rho / eps0, right hand side of Greens function
    def GradPhi(x, y):
        return -2*Charges*dGdr(k0,x,y,halfz) * (halfz/r_xyz(x,y,halfz))
    return GradPhi

def dblGaussianQuad(f, rangex, rangey, sample, weight):
    dl_x = rangex[1] - rangex[0]
    dl_y = rangey[1] - rangey[0]
    weight_x = weight / 2.0 * dl_x
    weight_y = weight / 2.0 * dl_y
    sample_x = sample_points1D(-rangex[0], rangex[0], sample)
    sample_y = sample_points1D(-rangey[0], rangey[0], sample)
    sample_xy = np.zeros((order, order), complex)
    weight_xy = np.zeros((order, order), float)
    for m in range(order):
        for n in range(order):
            sample_xy[m,n] = f(sample_x[m], sample_y[n])
            weight_xy[m,n] = weight_x[m] * weight_y[n]
    print sample_xy
    print weight_x
    print weight_xy
    Integral = np.sum(sample_xy * weight_xy)
    return Integral

if __name__ == "__main__":
    order = 5
    x, y, z = 0.0,0.0,0.0
    dl = 0.1
    omega = 1e6
    sample, weight = np.polynomial.legendre.leggauss(order)
#    integrand = GradPhiIntegrand(dl, omega)
    def integrand(x,y):
        return 1

    x_bound = [-dl/2, dl/2]

    y_bound = [-dl/2, dl/2]
    result1 = nquad(integrand, [x_bound, y_bound])
    result2 = dblGaussianQuad(integrand, x_bound, y_bound, sample, weight)
    print result1
    print result2
    #    def Integrand(x, y):
#        return LocalGradPhi(dl, omega, x, y)
#    print nquad(Integrand, [[-dl/2, dl/2],[-dl/2, dl/2]])
