import numpy as np
import GreensLocal
import EMConst

eps0 = EMConst.eps0
mu0 = EMConst.mu0
PI = np.pi

def RemoteA(dl, omega, p0, p1):
    # this function calcualtes the remote A vector
    # input
    # dl: edge length
    # omega: angular frequency
    # p0: source coordinate
    # p1: field coordinate
    k0 = omega * np.sqrt(eps0*mu0)
    Volume = dl**3
    r = np.linalg.norm(p1 - p0)
    C0 = mu0 * np.exp(-1j*k0*r)/(4*PI*r) * Volume


def RemoteGradPhi(dl, omega, p0, p1, direction):
    # dir, the direction of the current
    # direction could either be
    # (1,0,0)(0,1,0)(0,0,1)
    k0 = omega * np.sqrt(eps0*mu0)
    if diection = [1,0,0]:
        # the charge is aligned in y-z plane

