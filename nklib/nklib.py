import calc
import multilayer as ml

# Default wavelength range in nm
_wl_min = 305
_wl_max = 2600
_wl_step = 5


class Material(ml.Layer):
    '''
    Creates a customizable wl array and takes an arbitrary material's nk data
    and label to create ml.Layer.
    '''
    def __init__(self, thickness, mat_wl, mat_n, label,
                 wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                 unit="nm"):
        wl = np.arange(wl_min, wl_max, wl_step)
        n = calc.interpolate(mat_wl, mat_n, wl)
        ml.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class Ag(Material):
    '''
    Creates a Layer of silver. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, wl_min=wl_min, wl_max=wl_max,
                 wl_step=wl_step, unit="nm"):
        mat_wl = _Ag_wl
        mat_n = _Ag_nk
        label = 'Ag'
        Material.__init__(self, thickness, mat_wl, mat_n, label,
                          wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                          unit="nm")


class Al(Material):
    '''
    Creates a Layer of aluminum. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, wl_min=wl_min, wl_max=wl_max,
                 wl_step=wl_step, unit="nm"):
        mat_wl = _Al_wl
        mat_n = _Al_nk
        label = 'Al'
        Material.__init__(self, thickness, mat_wl, mat_n, label,
                          wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                          unit="nm")


class BK7(Material):
    '''
    Creates a Layer of BK7 glass. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, wl_min=wl_min, wl_max=wl_max,
                 wl_step=wl_step, unit="nm"):
        mat_wl = _BK7_wl
        mat_n = _BK7_nk
        label = 'BK7'
        Material.__init__(self, thickness, mat_wl, mat_n, label,
                          wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                          unit="nm")


class DLC3W(Material):
    '''
    Creates a Layer of DLC3W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, wl_min=wl_min, wl_max=wl_max,
                 wl_step=wl_step, unit="nm", extended=True):
        if extended:
            mat_wl = _DLC3W_ext_wl
            mat_n = _DLC3W_ext_nk
        else:
            mat_wl = _DLC3W_wl
            mat_n = _DLC3W_nk
        label = 'DLC3W'
        Material.__init__(self, thickness, mat_wl, mat_n, label,
                          wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                          unit="nm")


class DLC5W(Material):
    '''
    Creates a Layer of DLC5W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, wl_min=wl_min, wl_max=wl_max,
                 wl_step=wl_step, unit="nm", extended=True):
        if extended:
            mat_wl = _DLC5W_ext_wl
            mat_n = _DLC5W_ext_nk
        else:
            mat_wl = _DLC5W_wl
            mat_n = _DLC5W_nk
        label = 'DLC5W'
        Material.__init__(self, thickness, mat_wl, mat_n, label,
                          wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                          unit="nm")


class DLC10W(Material):
    '''
    Creates a Layer of DLC10W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, wl_min=wl_min, wl_max=wl_max,
                 wl_step=wl_step, unit="nm", extended=True):
        if extended:
            mat_wl = _DLC10W_ext_wl
            mat_n = _DLC10W_ext_nk
        else:
            mat_wl = _DLC10W_wl
            mat_n = _DLC10W_nk
        label = 'DLC10W'
        Material.__init__(self, thickness, mat_wl, mat_n, label,
                          wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                          unit="nm")


class DLC15W(Material):
    '''
    Creates a Layer of DLC15W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, wl_min=wl_min, wl_max=wl_max,
                 wl_step=wl_step, unit="nm", extended=True):
        if extended:
            mat_wl = _DLC15W_ext_wl
            mat_n = _DLC15W_ext_nk
        else:
            mat_wl = _DLC15W_wl
            mat_n = _DLC15W_nk
        label = 'DLC15W'
        Material.__init__(self, thickness, mat_wl, mat_n, label,
                          wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                          unit="nm")


class DLC20W(Material):
    '''
    Creates a Layer of DLC20W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, wl_min=wl_min, wl_max=wl_max,
                 wl_step=wl_step, unit="nm", extended=True):
        if extended:
            mat_wl = _DLC20W_ext_wl
            mat_n = _DLC20W_ext_nk
        else:
            mat_wl = _DLC20W_wl
            mat_n = _DLC20W_nk
        label = 'DLC20W'
        Material.__init__(self, thickness, mat_wl, mat_n, label,
                          wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                          unit="nm")


class DLC40W(Material):
    '''
    Creates a Layer of DLC40W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, wl_min=wl_min, wl_max=wl_max,
                 wl_step=wl_step, unit="nm", extended=True):
        if extended:
            mat_wl = _DLC40W_ext_wl
            mat_n = _DLC40W_ext_nk
        else:
            mat_wl = _DLC40W_wl
            mat_n = _DLC40W_nk
        label = 'DLC40W'
        Material.__init__(self, thickness, mat_wl, mat_n, label,
                          wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                          unit="nm")


class DLC60W(Material):
    '''
    Creates a Layer of DLC60W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, wl_min=wl_min, wl_max=wl_max,
                 wl_step=wl_step, unit="nm", extended=True):
        if extended:
            mat_wl = _DLC60W_ext_wl
            mat_n = _DLC60W_ext_nk
        else:
            mat_wl = _DLC60W_wl
            mat_n = _DLC60W_nk
        label = 'DLC60W'
        Material.__init__(self, thickness, mat_wl, mat_n, label,
                          wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                          unit="nm")


class ITO(Material):
    '''
    Creates a Layer of ITO. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, wl_min=wl_min, wl_max=wl_max,
                 wl_step=wl_step, unit="nm"):
        mat_wl = _ITO_wl
        mat_n = _ITO_nk
        label = 'ITO'
        Material.__init__(self, thickness, mat_wl, mat_n, label,
                          wl_min=wl_min, wl_max=wl_max, wl_step=wl_step,
                          unit="nm")


# Import data
_Ag_ndata = np.loadtxt("Ag_n.txt", skiprows=1)
_Ag_kdata = np.loadtxt("Ag_k.txt", skiprows=1)
_Ag_wl = Ag_ndata[:,0]
_Ag_nk = Ag_ndata[:,1] - 1j * Ag_kdata[:,1]

_Al_ndata = np.loadtxt("Al_n_Rakic1998.txt", skiprows=1)
_Al_kdata = np.loadtxt("Al_k_Rakic1998.txt", skiprows=1)
_Al_wl = Al_ndata[:,0]
_Al_nk = Al_ndata[:,1] - 1j * Al_kdata[:,1]

_BK7_ndata = np.loadtxt("BK7_n.txt", skiprows=1)
_BK7_kdata = np.loadtxt("BK7_k.txt", skiprows=1)
_BK7_wl = BK7_ndata[:,0]
_BK7_nk = BK7_ndata[:,0] - 1j * BK7_kdata[:,1]

_DLC_nkdata = [np.loadtxt("DLC{}W_nk.txt".format(n), skiprows=1) for \
               n in [3, 5, 10, 15, 20, 40, 60]]
_DLC_wl = [_DLC_nkdata[i][:,0] for i in range(7)]
_DLC_nk = [_DLC_nkdata[i][:,1] - 1j * _DLC_nkdata[i][:,2] for i in range(7)]
_DLC3W_wl, _DLC5W_wl, _DLC10W_wl, _DLC15W_wl, DLC20W_wl, DLC40W_wl, DLC60W_wl = _DLC_wl
_DLC3W_nk, _DLC5W_nk, _DLC10W_nk, _DLC15W_nk, DLC20W_nk, DLC40W_nk, DLC60W_nk = _DLC_nk

_DLC_ext_nkdata = [np.loadtxt("DLC{}W_extended_nk.txt".format(n), skiprows=1) \
                   for n in [3, 5, 10, 15, 20, 40, 60]]
_DLC_ext_wl = [_DLC_ext_nkdata[i][:,0] for i in range(7)]
_DLC_ext_nk = [_DLC_ext_nkdata[i][:,1] - 1j * _DLC_ext_nkdata[i][:,2] \
               for i in range(7)]
_DLC3W_ext_wl, _DLC5W_ext_wl, _DLC10W_ext_wl, _DLC15W_ext_wl, DLC20W_ext_wl, DLC40W_ext_wl, DLC60W_ext_wl = _DLC_ext_wl
_DLC3W_ext_nk, _DLC5W_ext_nk, _DLC10W_ext_nk, _DLC15W_ext_nk, DLC20W_ext_nk, DLC40W_ext_nk, DLC60W_ext_nk = _DLC_ext_nk

_ITO_ndata = np.loadtxt("ITO_n_Konig.txt", skiprows=1)
_ITO_kdata = np.loadtxt("ITO_k_Konig.txt", skiprows=1)
_ITO_wl = ITO_ndata[:,0]
_ITO_nk = ITO_ndata[:,0] - 1j * ITO_kdata[:,1]
