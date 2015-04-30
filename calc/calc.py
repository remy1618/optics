import numpy as np

eps0 = 8.85e-12
mu0 = 4e-7 * pi

def interpolate(wl_raw, f, wl):
    '''
    f is a function of wl_raw.
    wl is the new domain that f is mapped onto via interpolation.
    All inputs are 1darrays.
    '''
    if wl[0] < wl_raw[0]:
        i = 1
        while wl[i] < wl_raw[0]:
            i += 1
        wl = wl[i:]
    if wl[-1] > wl_raw[-1]:
        i = -2
        while wl[i] > wl_raw[-1]:
            i -= 1
        wl = wl[:i + 1]
        
    f_new = np.empty(len(wl), f.dtype)
    for i in range(len(wl)):
        j = 0
        while wl[i] > wl_raw[j]:
            j += 1
        # Now wl_raw[j - 1] < wl[i] <= wl_raw[j]
        # Linear interpolation y = y0 + (y1-y0)/(x1-x0) * (x-x0)
        f_new[i] = f[j-1] + (f[j] - f[j-1]) / (wl_raw[j] - wl_raw[j-1]) * \
                   (wl[i] - wl_raw[j-1])
    return f_new


def reflectance(k0, n, d, n0=1.0, ns=1.5):
    return _RT(k0, n, d)[0]

def transmittance(k0, n, d, n0=1.0, ns=1.5):
    return _RT(k0, n, d)[1]

def _transfer_matrix(k0, n, d):
    assert isinstance(n, (complex, float, int)), "n is not a number"
    assert isinstance(d, (complex, float, int)), "d is not a number"
    Y = n * sqrt(eps0 / mu0)
    return np.array([
                     [cos(k0 * n * d), 1j * sin(k0 * n * d) / Y],
                     [1j * Y * sin(k0 * n * d), cos(k0 * n * d)]
                    ])

def _RT(k0, n, d, n0=1.0, ns=1.5):
    '''
    n an array of refractive index for each layer from top to bottom.
    d an array of thickness for each layer from top to bottom.
    '''
    assert isinstance(n, (np.ndarray, list, tuple)), "n is not an array"
    assert isinstance(d, (np.ndarray, list, tuple)), "d is not an array"
    assert len(n) == len(d), "number of layers mismatch between n and d"
    assert isinstance(n[0], (complex, float, int)), "values in n array not numbers"
    assert isinstance(d[0], (complex, float, int)), "values in d array not numbers"
    Y0, Ys = n0 * sqrt(eps0 / mu0), ns * sqrt(eps0 / mu0)
    Y = n * sqrt(eps0 / mu0)
    m = 1
    for i in range(len(n)):
        m = np.dot(m, _transfer_matrix(k0, n[i], d[i]))
    r = (Y0 * m[0,0] + Y0 * Ys * m[0,1] - m[1,0] - Ys * m[1,1]) / \
        (Y0 * m[0,0] + Y0 * Ys * m[0,1] + m[1,0] + Ys * m[1,1])
    R = (np.absolute(r))**2
    t = 2 * Y0 / (Y0 * m[0,0] + Y0 * Ys * m[0,1] + m[1,0] + Ys * m[1,1])
    T = (np.absolute(t))**2 * Ys / Y0
    return R, T
