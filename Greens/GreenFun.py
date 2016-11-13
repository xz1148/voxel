import numpy as np
from EMConst import eps0
from EMConst import mu0
from scipy.integrate import dblquad
import numpy.pi as PI

def r_xyz(x, y, z):
    # this function calculates the distance form (x, y, z) to (0,0,0)
    r = np.sqrt(x**2+y**2+z**2)
    return r

def dGdr(k0, x, y, z):
    # dG / dr
    r0 = r_xyz(x, y, z)
    dGdr_r = (-1/(4*PI)) * (-1j*k0*np.exp(-1j*k0*r0))/r0 - np.exp(-1j*k0*r0)/r0**2)
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


def LocalGradPhi(dl, omega):
    halfz = dl / 2
    k0 = omega * sqrt(mu0 * eps0)
    Charges = 1 / (1j * omega)  # equals to -rho / eps0, right hand side of Greens function
    GradPhi = dblquad(lambda x, y: -2*Charges*dGdr(k0,x,y,halfz) * (halfz/r_xyz(x,y,halfz)),
                      -dl/2, dl/2, -dl/2, dl/2)

if __name__ = "__main__"
    x, y, z = 0.0,0.0,0.0
    dl = 0.1
    omega = 1e6
    print LocalGradPhi(dl, omega)

