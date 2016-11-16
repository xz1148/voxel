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
def r_xyz(p1, p0):
    # this function calculates the distance form (x, y, z) to (0,0,0)
    r = np.linalg.norm(p1 - p0)
    return r

def dGdr(k0, p1, p0):
    # dG / dr
    r0 = r_xyz(p1, p0)
    u = (p1 - p0) / r0
    dGdr_r = (-1/(4*PI)) * ((-1j*k0*np.exp(-1j*k0*r0))/r0 - np.exp(-1j*k0*r0)/r0**2)
    dGdr_xyz = dGdr_r*u
    return dGdr_xyz

def LocalA(dl, omega):
    # this function calculates the A vector contributed by a unit J( current )
    # input:
    # dl: the edge length of voxel
    k0 = omega * np.sqrt(mu0 * eps0)
    r_max = (3.0/(4.0*PI)) ** (1.0/3.0) * dl
    C0 = -mu0 / (4*PI*1j*k0)
    int_A_r_max = C0 * (r_max*np.exp(-1j*k0*r_max)+np.exp(-1j*k0*r_max)/(1j*k0))
    int_A_0 = C0 * (1/(1j*k0))
    A = int_A_r_max - int_A_0
    return A


def GradPhiIntegrand(omega, p1, p0):
    # this funtion returns a function which used as integrand
    # the integrand in xy plane
    k0 = omega * np.sqrt(mu0 * eps0)
    Charges = -1 / (1j * omega)  # equals to -rho / eps0, right hand side of Greens function
    def GradPhi_xy_source(x, y):
        p_source = np.array([p0[0]+x, p0[1]+y, p0[2]]) # the field point
        p_field = np.array([p1[0], p1[1], p1[2]]) # the field point
        dGdr_xyz = dGdr(k0, p_field, p_source)
        gradPhi = Charges*dGdr_xyz
        return gradPhi[0], gradPhi[1], gradPhi[2]
    return GradPhi_xy_source

def dblGaussianQuad(f, rangex, rangey, sample, weight):
    # this function calculates the double integration
    # f: function handle
    # rangex: the rangex of x
    # rangey: the rangey of y
    # sample and weight: the sample and weight given by gaussian quadrature
    dl_x = rangex[1] - rangex[0]
    dl_y = rangey[1] - rangey[0]
    weight_x = weight / 2.0 * dl_x
    weight_y = weight / 2.0 * dl_y
    sample_x = sample_points1D(-rangex[0], rangex[0], sample)
    sample_y = sample_points1D(-rangey[0], rangey[0], sample)
    sample_xy_x = np.zeros((order, order), complex)
    sample_xy_y = np.zeros((order, order), complex)
    sample_xy_z = np.zeros((order, order), complex)
    weight_xy = np.zeros((order, order), float)
    for m in range(order):
        for n in range(order):
            sample_xy_x[m,n], sample_xy_y[m,n], sample_xy_z[m,n] = f(sample_x[m], sample_y[n])
            weight_xy[m,n] = weight_x[m] * weight_y[n]
    Integral_x = np.sum(sample_xy_x * weight_xy)
    Integral_y = np.sum(sample_xy_y * weight_xy)
    Integral_z = np.sum(sample_xy_z * weight_xy)
    Integral = [Integral_x, Integral_y, Integral_z]
    return Integral

def GradPhi(omega, p1, p0, bound1, bound2, sample, weight, direction):
    # here direction means the normal vector of plane on which the source is
    # distributed, for instance, if the source is distributed on xy plane,
    # direction is [0, 0, 1]
    ux = np.array([1,0,0])
    uy = np.array([0,1,0])
    uz = np.array([0,0,1])
    if np.array_equal(direction, ux):
        M = np.matrix([[0.0,0.0,1.0],[0.0,1.0,0.0],[1.0,0.0,0.0]])
    elif np.array_equal(direction, uy):
        M = np.matrix([[0.0,0.0,1.0],[1.0,0.0,0.0],[0.0,1.0,0.0]])
    elif np.array_equal(direction, uz):
        M = np.matrix([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    p1_new = np.asarray((M*np.asmatrix(p1).T).T)[0]
    p0_new = np.asarray((M*np.asmatrix(p0).T).T)[0]
    integrand = GradPhiIntegrand(omega, p1_new, p0_new)
    gradPhi_new = dblGaussianQuad(integrand, bound1, bound2, sample, weight)
    gradPhi = np.asarray((np.linalg.inv(M)*np.asmatrix(gradPhi_new).T).T)[0]
    return gradPhi

def GradPhiRemote(omega, p1, p0, bound1, bound2):
    # because it calculates the remote field, the shape of charges can be
    # neglected
    area = (bound1[1] - bound1[0])*(bound2[1] - bound2[0])
    charges = -1/(1j * omega) * area
    k0 = omega * np.sqrt(eps0 * mu0)
    gradPhi = charges * dGdr(k0, p1, p0)
    return gradPhi

def ARemote(omega, dl, p1, p0):
    # this function calculates the remote A vector
    # the A vector can point to x and y and z directions
    k0 = omega*np.sqrt(mu0*eps0)
    Volume = dl**3
    r = r_xyz(p1, p0)
    G = -np.exp(-1j*k0*r) / (4*PI*r)
    A = -mu0*G*Volume
    return A
#def LocalE(dl, omega, sample, weight):
#    E_A = -1j*omega*LocalA(dl, omega)
#    rangex = [-dl/2, dl/2]
#    rangey = [-dl/2, dl/2]
#    E_GradPhi = -dblGaussianQuad(GradPhiIntegrand(dl,omega), rangex, rangey, sample, weight)
#    E = E_A + E_GradPhi
#    return E

if __name__ == "__main__":
    order = 10
    sample, weight = np.polynomial.legendre.leggauss(order)
    dl = 0.1
    omega = 1e6
    p1 = np.array([0.05, 0.0, 0.0])
    p0 = np.array([0.0, 0.0, 0.0])
    integrand = GradPhiIntegrand(omega, p1, p0)
    x_bound = [-dl/2, dl/2]
    y_bound = [-dl/2, dl/2]
    direction = np.array([1, 0, 0])
    k0 = omega * np.sqrt(eps0 * mu0)

    result2 = GradPhi(omega, p1, p0, x_bound, y_bound, sample, weight, direction)
    result3 = GradPhiRemote(omega, p1, p0, x_bound, y_bound)
    print result2

