import os
import numpy as np
import calc
import structure as st


class Ag(st.Layer):
    '''
    Creates a Layer of silver. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _Ag_wl / 1e3 if unit == 'micron' else _Ag_wl
        n = _Ag_nk
        label = 'Ag'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class Al(st.Layer):
    '''
    Creates a Layer of aluminum. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _Al_wl / 1e3 if unit == 'micron' else _Al_wl
        n = _Al_nk
        label = 'Al'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class BK7(st.Layer):
    '''
    Creates a Layer of BK7 glass. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _BK7_wl / 1e3 if unit == 'micron' else _BK7_wl
        n = _BK7_nk
        label = 'BK7'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC3W(st.Layer):
    '''
    Creates a Layer of DLC3W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC3W_ext_wl if extended else _DLC3W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC3W_ext_nk if extended else _DLC3W_nk
        label = 'DLC3W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC5W(st.Layer):
    '''
    Creates a Layer of DLC5W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC5W_ext_wl if extended else _DLC5W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC5W_ext_nk if extended else _DLC5W_nk
        label = 'DLC5W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC10W(st.Layer):
    '''
    Creates a Layer of DLC10W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC10W_ext_wl if extended else _DLC10W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC10W_ext_nk if extended else _DLC10W_nk
        label = 'DLC10W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC15W(st.Layer):
    '''
    Creates a Layer of DLC15W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC15W_ext_wl if extended else _DLC15W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC15W_ext_nk if extended else _DLC15W_nk
        label = 'DLC15W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC20W(st.Layer):
    '''
    Creates a Layer of DLC20W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC20W_ext_wl if extended else _DLC20W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC20W_ext_nk if extended else _DLC20W_nk
        label = 'DLC20W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC40W(st.Layer):
    '''
    Creates a Layer of DLC40W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC40W_ext_wl if extended else _DLC40W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC40W_ext_nk if extended else _DLC40W_nk
        label = 'DLC40W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC60W(st.Layer):
    '''
    Creates a Layer of DLC60W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC40W_ext_wl if extended else _DLC40W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC40W_ext_nk if extended else _DLC40W_nk
        label = 'DLC60W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class ITO(st.Layer):
    '''
    Creates a Layer of ITO. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _ITO_wl / 1e3 if unit == 'micron' else _ITO_wl
        n = _ITO_nk
        label = 'ITO'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class ZnO(st.Layer):
    '''
    Creates a Layer of ZnO. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _ZnO_wl / 1e3 if unit == 'micron' else _ZnO_wl
        n = _ZnO_nk
        label = 'ZnO'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


# Import data
os.chdir(os.getcwd() + "\\opt_sim\\nklib_data") # possibly generalize

_Ag_ndata = np.loadtxt("Ag_n.txt", skiprows=1)
_Ag_kdata = np.loadtxt("Ag_k.txt", skiprows=1)
_Ag_wl = _Ag_ndata[:,0] * 1000
_Ag_nk = _Ag_ndata[:,1] - 1j * _Ag_kdata[:,1]

_Al_ndata = np.loadtxt("Al_n_Rakic1998.txt", skiprows=1)
_Al_kdata = np.loadtxt("Al_k_Rakic1998.txt", skiprows=1)
_Al_wl = _Al_ndata[:,0] * 1000
_Al_nk = _Al_ndata[:,1] - 1j * _Al_kdata[:,1]

_BK7_ndata = np.loadtxt("BK7_n.txt", skiprows=1)
_BK7_kdata = np.loadtxt("BK7_k.txt", skiprows=1)
_BK7_wl = _BK7_ndata[:,0] * 1000
# BK7 doesn't have the same number of data points for n and k
_wl_temp = _BK7_kdata[:,0] * 1000
_k_temp = _BK7_kdata[:,1]
_BK7_wl, _BK7_k = calc.interpolate(_wl_temp, _k_temp, _BK7_wl)
_BK7_nk = _BK7_ndata[:,1] - 1j * _BK7_k

_DLC_nkdata = [np.loadtxt("DLC{}W_nk.txt".format(n), skiprows=1) for \
               n in [3, 5, 10, 15, 20, 40, 60]]
_DLC_wl = [_DLC_nkdata[i][:,0] for i in range(7)]
_DLC_nk = [_DLC_nkdata[i][:,1] - 1j * _DLC_nkdata[i][:,2] for i in range(7)]
_DLC3W_wl, _DLC5W_wl, _DLC10W_wl, _DLC15W_wl, _DLC20W_wl, _DLC40W_wl, _DLC60W_wl = _DLC_wl
_DLC3W_nk, _DLC5W_nk, _DLC10W_nk, _DLC15W_nk, _DLC20W_nk, _DLC40W_nk, _DLC60W_nk = _DLC_nk

_DLC_ext_nkdata = [np.loadtxt("DLC{}W_extended_nk.txt".format(n), skiprows=1) \
                   for n in [3, 5, 10, 15, 20, 40, 60]]
_DLC_ext_wl = [_DLC_ext_nkdata[i][:,0] for i in range(7)]
_DLC_ext_nk = [_DLC_ext_nkdata[i][:,1] - 1j * _DLC_ext_nkdata[i][:,2] \
               for i in range(7)]
_DLC3W_ext_wl, _DLC5W_ext_wl, _DLC10W_ext_wl, _DLC15W_ext_wl, _DLC20W_ext_wl, _DLC40W_ext_wl, _DLC60W_ext_wl = _DLC_ext_wl
_DLC3W_ext_nk, _DLC5W_ext_nk, _DLC10W_ext_nk, _DLC15W_ext_nk, _DLC20W_ext_nk, _DLC40W_ext_nk, _DLC60W_ext_nk = _DLC_ext_nk
DLC_list = [DLC3W, DLC5W, DLC10W, DLC15W, DLC20W, DLC40W, DLC60W]

_ITO_ndata = np.loadtxt("ITO_n_Konig.txt", skiprows=1)
_ITO_kdata = np.loadtxt("ITO_k_Konig.txt", skiprows=1)
_ITO_wl = _ITO_ndata[:,0] * 1000
_ITO_nk = _ITO_ndata[:,1] - 1j * _ITO_kdata[:,1]

_ZnO_nkdata = np.loadtxt("ZnO_nk.txt", skiprows=1)
_ZnO_wl = _ZnO_nkdata[:,0] * 1000
_ZnO_nk = _ZnO_nkdata[:,1] - 1j * _ZnO_nkdata[:,2]

os.chdir("..")
