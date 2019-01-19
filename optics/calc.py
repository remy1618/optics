from numpy import pi, cos, sin, sqrt
import numpy as np

eps0 = 8.85e-12
mu0 = 4e-7 * pi

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

def transmittance(k0, n, d, n0=1.0, ns=1.5):
    return _TR(k0, n, d, n0, ns)[0]

def reflectance(k0, n, d, n0=1.0, ns=1.5):
    return _TR(k0, n, d, n0, ns)[1]

def _transfer_matrix(k0, n, d):
    '''
    Outputs a 3d array of shape (N, 2, 2) where N is the number of wavevectors
    and each of the 2-by-2 sub-array associated with the particular wavevector.
    '''
    Y = n * sqrt(eps0 / mu0)
    # Invidual matrix element definition
    # cos(arg) and sin(arg) may overflow when metal thickness > 1 micron
    M00 = cos(k0 * n * d)
    M01 = 1j * sin(k0 * n * d) / Y
    M10 = 1j * sin(k0 * n * d) * Y
    M11 = M00
    for a in [n, M00, M01, M10, M11]:
        assert k0.shape == a.shape,\
        "k0.shape " + str(k0.shape) + ", mismatch.shape " + str(a.shape)
    M = np.empty((k0.shape[0], 2, 2), dtype=complex)
    M[:,0,0] = M00
    M[:,0,1] = M01
    M[:,1,0] = M10
    M[:,1,1] = M11
    return M

def TandR(k0, n, d, n0, ns):
    '''
    k0 1d array of length #_of_wavelengths
    n 2d array of length (#_of_layers, #_of_wavelengths)
    d 1d array of length (#_of_layers)
    '''
    num_of_wavelengths = k0.shape[0]
    num_of_layers = d.shape[0]
    
    Y0, Ys = n0 * sqrt(eps0 / mu0), ns * sqrt(eps0 / mu0)
    Y = n * sqrt(eps0 / mu0)

    M = np.zeros((num_of_wavelengths, 2, 2), complex)
    M[:,0,0], M[:,1,1] = 1, 1 # num_of_wavelengths-array of identity matrix
    t = np.empty(num_of_wavelengths, complex)
    r = np.empty(num_of_wavelengths, complex)
    # Maybe possible to speed up here by taking out loop
    for i in range(num_of_layers):
        M = np.einsum('ijk,ikl->ijl',
                      M, _transfer_matrix(k0, n[i], d[i]))
    # Denominator definition (to avoid multi-line equations)
    D = (Y0 * M[:,0,0] + Y0 * Ys * M[:,0,1] + M[:,1,0] + Ys * M[:,1,1])
    # Numerator for reflectance
    Nu = (Y0 * M[:,0,0] + Y0 * Ys * M[:,0,1] - M[:,1,0] - Ys * M[:,1,1])
    t = 2 * Y0 / D
    r = Nu / D
    T = (np.absolute(t))**2 * Ys / Y0
    R = (np.absolute(r))**2
    return T, R
