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
    f = np.empty(len(wl), f_raw.dtype)
    for i in range(len(wl)):
        j = 0
        while wl[i] > wl_raw[j]:
            j += 1
        # Now wl_raw[j - 1] < wl[i] <= wl_raw[j]
        # Linear interpolation y = y0 + (y1-y0)/(x1-x0) * (x-x0)
        f[i] = f_raw[j-1] + \
                   (f_raw[j] - f_raw[j-1]) / (wl_raw[j] - wl_raw[j-1]) * \
                   (wl[i] - wl_raw[j-1])
    return f

def transmittance(k0, n, d, n0=1.0, ns=1.5):
    return _TR(k0, n, d, n0, ns)[0]

def reflectance(k0, n, d, n0=1.0, ns=1.5):
    return _TR(k0, n, d, n0, ns)[1]

def _transfer_matrix(k0, n, d):
    assert isinstance(n, (complex, float, int)), "n is not a number but type " + str(type(n))
    assert isinstance(d, (complex, float, int)), "d is not a number but type " + str(type(d))
    Y = n * sqrt(eps0 / mu0)
    return np.array([
                     [cos(k0 * n * d), 1j * sin(k0 * n * d) / Y],
                     [1j * Y * sin(k0 * n * d), cos(k0 * n * d)]
                    ])

def _TR(k0, n, d, n0, ns):
    '''
    Calculates transmittance and reflectance at a single wavelength.
    n is an array of refractive index for each layer from top to bottom.
    d is an array of thickness for each layer from top to bottom.
    '''
    assert isinstance(n, (np.ndarray, list, tuple)), "n is not an array but type " + str(type(n))
    assert isinstance(d, (np.ndarray, list, tuple)), "d is not an array but type " + str(type(d))
    assert len(n) == len(d), "number of layers mismatch between n and d: " + "n = " + str(len(n)) + ", d = " + str(len(d))
    assert isinstance(n[0], (complex, float, int)), "values in n array not numbers but type " + str(type(n))
    assert isinstance(d[0], (complex, float, int)), "values in d array not numbers but type " + str(type(d))
    Y0, Ys = n0 * sqrt(eps0 / mu0), ns * sqrt(eps0 / mu0)
    Y = n * sqrt(eps0 / mu0)
    m = 1
    for i in range(len(n)):
        m = np.dot(m, _transfer_matrix(k0, n[i], d[i]))
    t = 2 * Y0 / (Y0 * m[0,0] + Y0 * Ys * m[0,1] + m[1,0] + Ys * m[1,1])
    T = (np.absolute(t))**2 * Ys / Y0
    r = (Y0 * m[0,0] + Y0 * Ys * m[0,1] - m[1,0] - Ys * m[1,1]) / \
        (Y0 * m[0,0] + Y0 * Ys * m[0,1] + m[1,0] + Ys * m[1,1])
    R = (np.absolute(r))**2
    return T, R
