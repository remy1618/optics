import numpy as np

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
