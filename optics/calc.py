from numpy import pi, cos, sin, arcsin, sqrt
import numpy as np

# Taken from Macleod - Thin Film Optical Filters (2001)
# Variable names
# n refractive index
# a incident angle
# B element of characteristic matrix, normalized E-field amplitude
# C element of characteristic matrix, normalized H-field amplitude
# Y = n * Y0 = C / B optical admittance
# Np = Y / cos(a) tilted optical admittance for p-polarization
# Ns = Y * cos(a) tilted optical admittance for s-polarization
eps0 = 8.85e-12
mu0 = 4e-7 * pi
Yfs = sqrt(eps0 / mu0) # fs for free-space

def _transfer_matrix(k0, n, d, a, pol):
    '''
    Outputs a 3d array of shape (N, 2, 2) where N is the number of wavevectors
    and each of the 2-by-2 sub-array associated with the particular wavevector.
    '''
    Y = n * Yfs
    if pol == "p":
        N = Y / cos(a)
    if pol == "s":
        N = Y * cos(a)
    # Invidual matrix element definition
    # cos(arg) and sin(arg) may overflow when metal thickness > 1 micron
    M00 = cos(k0 * n * d * cos(a))
    M01 = 1j * sin(k0 * n * d * cos(a)) / N
    M10 = 1j * sin(k0 * n * d * cos(a)) * N
    M11 = M00
    for arr in [n, M00, M01, M10, M11]:
        assert k0.shape == arr.shape,\
        "k0.shape " + str(k0.shape) + ", mismatch.shape " + str(arr.shape)
    M = np.empty((k0.shape[0], 2, 2), dtype=complex)
    M[:,0,0] = M00
    M[:,0,1] = M01
    M[:,1,0] = M10
    M[:,1,1] = M11
    return M

def _BandC_from_transfer_matrix(k0, n, d, n0, ns, a0, pol):
    '''
    k0 1d array of length #_of_wavelengths
    n 2d array of length (#_of_layers, #_of_wavelengths)
    d 1d array of length (#_of_layers)
    '''
    num_of_wavelengths = k0.shape[0]
    num_of_layers = d.shape[0]
    
    Y0 = n0 * Yfs
    Ys = ns * Yfs
    Y = n * Yfs
    a0 *= 2*pi/360

    # Initialize M as num_of_wavelengths-dimensional array of identity matrix
    M = np.zeros((num_of_wavelengths, 2, 2), complex)
    M[:,0,0], M[:,1,1] = 1, 1
    for i in range(num_of_layers):
        # Maybe possible to speed up here by taking out loop
        a = arcsin(n0 * sin(a0) / n[i])
        M = np.einsum('ijk,ikl->ijl', M, _transfer_matrix(k0, n[i], d[i], a, pol))
    B = M[:,0,0] + Ys * M[:,0,1]
    C = M[:,1,0] + Ys * M[:,1,1]
    return B, C

def TandR(k0, n, d, n0, ns, a0, pol):
    Y0 = n0 * Yfs
    Ys = ns * Yfs
    B, C = _BandC_from_transfer_matrix(k0, n, d, n0, ns, a0, pol)
    r = (Y0*B - C) / (Y0*B + C)
    R = (np.absolute(r))**2
    T = 4 * Y0 * Ys.real / (np.absolute(Y0*B + C))**2
    return T, R
    
def tanPSIandcosDELTA(k0, n, d, n0, ns, a0):
    Y0 = n0 * Yfs
    Bp, Cp = _BandC_from_transfer_matrix(k0, n, d, n0, ns, a0, pol="p")
    Bs, Cs = _BandC_from_transfer_matrix(k0, n, d, n0, ns, a0, pol="s")
    rp = (Y0*Bp - Cp) / (Y0*Bp + Cp)
    rs = (Y0*Bs - Cs) / (Y0*Bs + Cs)
    rho = rp / rs
    tanPSI = np.absolute(rho)
    DELTA = np.angle(rho)
    cosDELTA = cos(DELTA)
    return tanPSI, cosDELTA