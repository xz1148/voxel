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


def GradPhiIntegrand(dl, omega, p1, p0):
    # this funtion returns a function which used as integrand
    # the integrand in xy plane
    halfz = dl / 2
    k0 = omega * np.sqrt(mu0 * eps0)
    Charges = -1 / (1j * omega)  # equals to -rho / eps0, right hand side of Greens function
    def GradPhi_x(x, y):
        p_source = np.array([p0[0]+x, p0[1]+y, p0[2]]) # the field point
        p_field = np.array([p1[0], p1[1], p1[2]]) # the field point
        dGdr_x,_,_ = dGdr(k0, p1, p0)
        return Charges*dGdr_x
    def GradPhi_y(x, y):
        p_source = np.array([p0[0]+x, p0[1]+y, p0[2]]) # the field point
        p_field = np.array([p1[0], p1[1], p1[2]]) # the field point
        _,dGdr_y,_ = dGdr(k0, p1, p0)
        return Charges*dGdr_y
    def GradPhi_z(x, y):
        p_source = np.array([p0[0]+x, p0[1]+y, p0[2]]) # the field point
        p_field = np.array([p1[0], p1[1], p1[2]]) # the field point
        _,_,dGdr_z = dGdr(k0, p1, p0)
        return Charges*dGdr_z
    return GradPhi_x, GradPhi_y, GradPhi_z

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
    Integral = np.sum(sample_xy * weight_xy)
    return Integral

def LocalE(dl, omega, sample, weight):
    E_A = -1j*omega*LocalA(dl, omega)
    rangex = [-dl/2, dl/2]
    rangey = [-dl/2, dl/2]
    E_GradPhi = -dblGaussianQuad(GradPhiIntegrand(dl,omega), rangex, rangey, sample, weight)
    E = E_A + E_GradPhi
    return E

if __name__ == "__main__":
    order = 10
    sample, weight = np.polynomial.legendre.leggauss(order)
    dl = 0.1
    omega = 1e6
    integrand = GradPhiIntegrand(dl, omega)
    #def integrand(x,y):
    #    return 1
    x_bound = [-dl/2, dl/2]
    y_bound = [-dl/2, dl/2]
    result1 = nquad(integrand, [x_bound, y_bound])
    result2 = dblGaussianQuad(integrand, x_bound, y_bound, sample, weight)
    result3 = LocalE(dl, omega, sample, weight)
    print result3
    print result1
    print result2
    #    def Integrand(x, y):
#        return LocalGradPhi(dl, omega, x, y)
#    print nquad(Integrand, [[-dl/2, dl/2],[-dl/2, dl/2]])
