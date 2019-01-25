from numpy import pi, cos, sin, sqrt
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
Y0 = sqrt(eps0 / mu0)

def interpolate(wl_raw, f_raw, wl):
    '''
    Precondition: wl is within the range of wl_raw
    f is a function of wl_raw.
    wl is the new domain that f is mapped onto via interpolation.
    All inputs are 1darrays.
    '''
    if f_raw.dtype == complex:
        return np.interp(wl, wl_raw, f_raw.real) + \
               1j * np.interp(wl, wl_raw, f_raw.imag)
    return np.interp(wl, wl_raw, f_raw)

def _transfer_matrix(k0, n, d, a, pol="p"):
    '''
    Outputs a 3d array of shape (N, 2, 2) where N is the number of wavevectors
    and each of the 2-by-2 sub-array associated with the particular wavevector.
    '''
    Y = n * Y0
    a *= 2*pi/360.
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

def TandR(k0, n, d, n0, ns, a):
    '''
    k0 1d array of length #_of_wavelengths
    n 2d array of length (#_of_layers, #_of_wavelengths)
    d 1d array of length (#_of_layers)
    '''
    num_of_wavelengths = k0.shape[0]
    num_of_layers = d.shape[0]
    
    Ys = ns * Y0
    Y = n * Y0

    # Initialize M as num_of_wavelengths-dimensional array of identity matrix
    M = np.zeros((num_of_wavelengths, 2, 2), complex)
    M[:,0,0], M[:,1,1] = 1, 1
    for i in range(num_of_layers):
        # Maybe possible to speed up here by taking out loop
        M = np.einsum('ijk,ikl->ijl',
                      M, _transfer_matrix(k0, n[i], d[i], a))
    B = M[:,0,0] + Ys * M[:,0,1]
    C = M[:,1,0] + Ys * M[:,1,1]

    r = (Y0*B - C) / (Y0*B + C)
    R = (np.absolute(r))**2
    T = 4 * Y0 * Ys.real / (np.absolute(Y0*B + C))**2

    return T, R
