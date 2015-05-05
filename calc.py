from numpy import pi, cos, sin, sqrt
import numpy as np

eps0 = 8.85e-12
mu0 = 4e-7 * pi

def interpolate(wl_raw, f_raw, wl):
    ''' Update documentation
    f is a function of wl_raw.
    wl is the new domain that f is mapped onto via interpolation.
    All inputs are 1darrays.
    '''

    # These 2 blocks might be unneeded now that structure is modified
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
    return wl, f

def TR_spectrum(struct, n0, ns):
    '''
    struct is a structure.MultiLayer
    '''
    if struct.unit == 'nm':
        wl = struct.wl / 1e9
    if struct.unit == 'micron':
        wl = struct.wl / 1e6
    k0 = 2 * pi / wl
    # room for np.ndarray optimization here
    temp_n = [layer.n for layer in struct._layers_list[::-1]]
    n = zip(*temp_n)
    d = [layer.d / 1e9 if struct.unit == 'nm' else \
         layer.d / 1e6 if struct.unit == 'micron' else None\
         for layer in struct._layers_list[::-1]]
    wl_len = len(struct.wl)
    T = np.empty(wl_len)
    R = np.empty(wl_len)
    for i in range(wl_len):
        T[i] = transmittance(k0[i], n[i], d, n0=n0, ns=ns)
        R[i] = reflectance(k0[i], n[i], d, n0=n0, ns=ns)
    return T, R

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
