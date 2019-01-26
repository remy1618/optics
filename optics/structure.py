import numpy as np
import calc
import os.path
# Append path to support files
dir_path = os.path.dirname(os.path.realpath(__file__))
sup_path = dir_path + "\\support_files\\"

class Layer:

    def __init__(self, wl, n, d, label='untitled', unit='nm'):
        '''
        wl: wavelengths in nm (default) or micron
        n: wl-dependent refractive index of the layer
        d: thickness in nm (default) or micron of the layer
        label: a name for the layer
        unit: unit for distance, either in nm or micron
        '''
        if not isinstance(wl, (np.ndarray, list, tuple)):
            raise DataFormatException(
                "A range of wavelengths should be given as an array, list, "
                "or tuple.")
        if not isinstance(n, (np.ndarray, list, tuple)):
            raise DataFormatException(
                "Wavelength-dependent refractive index should be given as "
                "an array, list, or tuple.")
        if len(wl) != len(n):
            raise DataFormatException(
                "wl and n should have the same length.\nwl " + str(len(wl)) + \
                "\nn " + str(len(n)))
        if not unit in ['nm', 'micron']:
            raise DataFormatException(
                "unit should be given as 'nm' or 'micron'")      
        self.wl = wl
        self.n = n
        self.d = d
        self.label = label
        self.unit = unit

    def __repr__(self):
        '''
        Return a description of the layer.
        '''
        return '{} {} {}'.format(self.label, self.d, self.unit)
        

class MultiLayer:
    '''
    Initialize a multilayer structure. The structure is empty by default,
    unless the layers created from the Layer class are given in a list. The
    layers are added onto or removed from the multilayer structure using
    methods add_layer(), remove_layer(), or replace_layer().
    '''

    def __init__(self, layer_list=[], unit='nm', label=None,
                 min_wl=200, max_wl=2600, wl_step=5,
                 n0=1.0, ns=1.5):
        # Alternate constructor call for convenience
        if isinstance(layer_list, Layer):
            layer_list = [layer_list]
        if not unit in ['nm', 'micron']:
            raise DataFormatException(
                "unit should be given as 'nm' or 'micron'")
        for v in [min_wl, max_wl, wl_step]:
            if not isinstance(v, (int, float)):
                raise DataFormatException(
                    "wl bounds should be given as a number. {} given".format(
                    type(v)))
        
        self._layers_list = []
        self.unit = unit
        self.label = '' if not label else label
        self._user_label = False if not label else True
        self._min_wl = min_wl
        self._max_wl = max_wl
        self._wl_step = wl_step
        self.wl = np.arange(self._min_wl, self._max_wl, self._wl_step)
        self.wl_by_layer = ''
        self.n0 = n0
        self.ns = ns
        
        self.T = None
        self.R = None
        self.A = None
        self.T_color = None
        self.R_color = None
        self.tanPSI = None
        self.cosDELTA = None
        self._TR_calculated = False
        self._color_calculated = False
        self._tanPSIcosDELTA_calculated = False
        
        if layer_list:
            for layer in layer_list:
                self.add_layer(layer)
                self.wl_by_layer += layer.label + ' ' + str(layer.wl[0]) + \
                                    ' - ' + str(layer.wl[-1]) + ' ' + \
                                    layer.unit + '\n'
        
    def add_layer(self, layer, index=None):
        '''
        Add a layer to the top of multilayer if index is not given. Otherwise
        insert a layer indexed from 0.
        '''
        layer_wl = layer.wl * 1e3 if self.unit == 'nm' and \
                                     layer.unit == 'micron' else \
                   layer.wl / 1e3 if self.unit == 'micron' and \
                                     layer.unit == 'nm' else \
                   layer.wl
        self._min_wl = max(self._min_wl, layer_wl[0])
        self._max_wl = min(self._max_wl, layer_wl[-1])
        # If _max_wl is a a value such that _min_wl + N*_wl_step for some N,
        # add a small value to _max_wl so the last value will not be cut off
        # by np.arange
        if np.isclose((self._max_wl - self._min_wl) % self._wl_step, 0):
            self._max_wl += 0.01492873969 * self._wl_step
        self.wl = np.arange(self._min_wl, self._max_wl, self._wl_step)
                    
        index = len(self._layers_list) if not index else index
        self._layers_list.insert(index, layer)

        if not self._user_label:
            label_list = self.label.split('-') if self.label else []
            label_list.insert(index, layer.label + str(int(layer.d)))
            self.label = '-'.join(label_list)
        self._clear_calculations()
    
    def remove_layer(self, index=-1):
        '''
        Remove the top layer of the multilayer if index is not given. Otherwise
        the remove the layer indexed from 0.
        '''
        self._layers_list.pop(index)
        label_list = self.label.split('-')
        label_list.pop(index)
        if not self._user_label:
            self.label = '-'.join(label_list)
        self._clear_calculations()

    def replace_layer(self, index, layer):
        '''
        Replace the layer indexed from 0 with the newly given layer.
        '''
        self.remove_layer(index)
        self.add_layer(layer, index)

    def get_layer(self, index):
        '''
        Return the layer indexed from 0.
        '''
        return self._layers_list[index]

    def __repr__(self):
        '''
        Return a graphical representation of the structure.
        '''
        if not self._layers_list:
            return ''
        separator = '\n' + '-' * max([len(layer.__repr__())
                                      for layer in self._layers_list])
        structure_view = '\nTop'
        for layer in self._layers_list[::-1]:
            structure_view += separator + '\n' + layer.__repr__()
        structure_view += separator + '\n'
        return structure_view

    def _get_k0_n_d_from_self(self):
        num_layers = len(self._layers_list)
        num_wl = len(self.wl)
        # n a 2d array of layers (row) and wavelength (column)
        n = np.empty([num_layers, num_wl], complex)
        for i, L in enumerate(self._layers_list):
            L_wl = L.wl * 1e3 if self.unit == 'nm' and L.unit == 'micron' else\
                   L.wl / 1e3 if self.unit == 'micron' and L.unit == 'nm' else\
                   L.wl
            layer_n = np.interp(self.wl, L_wl, L.n)
            n[num_layers - 1 - i] = layer_n
        wl = self.wl / 1e9 if self.unit == 'nm' else self.wl / 1e6
        k0 = 2 * np.pi / wl
        # d has error: need to account for mismatchc of Layer.unit and ML.unit
        d = np.array([layer.d / 1e9 if layer.unit == 'nm' else layer.d / 1e6 \
                      for layer in self._layers_list[::-1]])
        return k0, n, d

    def calculate_TR(self, a0=0., pol="p"):
        if not self._layers_list:
            raise EmptyStructureException("Structure is empty.")
        k0, n, d = self._get_k0_n_d_from_self()
        self.T, self.R = calc.TandR(k0, n, d, n0=self.n0, ns=self.ns ,a0=a0, pol=pol)
        self.A = 1 - self.T - self.R
        self._TR_calculated = True

    def calculate_tanPSIcosDELTA(self, angle=53.):
        k0, n, d = self._get_k0_n_d_from_self()
        self.tanPSI, self.cosDELTA = calc.tanPSIandcosDELTA(k0, n, d, n0=self.n0, ns=self.ns, a0=angle)
        self._tanPSIcosDELTA_calculated = True

    def calculate_color(self):
        '''
        Generate a color tuple of (r, g, b) normalized to 255 for transmittance
        and reflectance.
        '''
        if not self._TR_calculated:
            self.calculate_TR()
            
        rgbdata = np.loadtxt(sup_path + "rgbdata.txt", skiprows=1)

        rgbwl = rgbdata[:,0]
        struct_wl = self.wl * 1e3 if self.unit == 'micron' else self.wl
        if struct_wl[0] > rgbwl[0] or struct_wl[-1] < rgbwl[-1]:
            raise Exception("The entire visible spectrum (390-830nm) must be available in the structure simulation.")
        
        rgb_sens = rgbdata[:,1:]
        rgb_sens /= np.sum(rgb_sens, 0) # Normalize sum to 1

        T_interp = np.interp(rgbwl, struct_wl, self.T)
        R_interp = np.interp(rgbwl, struct_wl, self.R)

        # Gives rgb in 0~1 float range
        self.T_color = np.sum(T_interp[:,np.newaxis] * rgb_sens, 0)
        self.R_color = np.sum(R_interp[:,np.newaxis] * rgb_sens, 0)

        # Gives rgb in 0~256 range
##        T_rgb = np.sum(T_interp[:,np.newaxis] * rgb_sens, 0) * 255
##        R_rgb = np.sum(R_interp[:,np.newaxis] * rgb_sens, 0) * 255
##        
##        self.T_color = tuple(T_rgb.astype(int))
##        self.R_color = tuple(R_rgb.astype(int))
                
    def _clear_calculations(self):
        self.T = None
        self.R = None
        self.A = None
        self.T_color = None
        self.R_color = None
        self.tanPSI = None
        self.cosDELTA = None
        self._TR_calculated = False
        self._color_calculated = False
        self._tanPSIcosDELTA_calculated = False

    def copy(self):
        return MultiLayer(layer_list=self._layers_list[:],
                          unit=self.unit, label=self.label,
                          min_wl=self._min_wl, max_wl=self._max_wl,
                          wl_step=self._wl_step, n0=self.n0,
                          ns=self.ns)



class DataFormatException(Exception):
    pass


class EmptyStructureException(Exception):
    pass
