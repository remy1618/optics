import calc
import multilayer as ml

# Default wavelength range in nm
_wl_min = 305
_wl_max = 2600
_wl_step = 5

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
_BK7_nkdata = BK7_ndata[:,0] - 1j * BK7_kdata[:,1]

_ITO_ndata = np.loadtxt("ITO_n_Konig.txt", skiprows=1)
_ITO_kdata = np.loadtxt("ITO_k_Konig.txt", skiprows=1)
_ITO_wl = ITO_ndata[:,0]
_ITO_nkdata = ITO_ndata[:,0] - 1j * ITO_kdata[:,1]



##Ag_data = np.loadtxt("ag.txt")
##DLC3W_data = np.loadtxt("DLC3W.txt")
##DLC5W_data = np.loadtxt("DLC5W.txt")
##DLC10W_data = np.loadtxt("DLC10W.txt")
##DLC20W_data = np.loadtxt("DLC20W.txt")
##DLC40W_data = np.loadtxt("DLC40W.txt")
##DLC60W_data = np.loadtxt("DLC60W.txt")
##Al_data = np.loadtxt("Al.txt")
##BK7_data = np.loadtxt("BK7.txt")

#process all the stuff


class Material(ml.Layer):
    '''
    Creates a customizable wl array and takes an arbitrary material's nk data
    and label to feed into Layer.
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
