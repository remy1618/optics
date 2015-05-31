import numpy as np
import calc

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
        self.wl_step = wl_step
        self.wl = None
        self.wl_by_layer = ''
        self.n0 = n0
        self.ns = ns
        self.T = None   # Do getter control with _TR_calculated
        self.R = None
        self.A = None
        self._TR_calculated = False
        
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
        if not self._layers_list:
            self.wl = layer_wl
        new_min_wl = max(self.wl[0], layer_wl[0])
        new_max_wl = min(self.wl[-1], layer_wl[-1])
        self.wl = np.arange(new_min_wl, new_max_wl, self.wl_step)
                    
        index = len(self._layers_list) if not index else index
        self._layers_list.insert(index, layer)

        if not self._user_label:
            label_list = self.label.split('-') if self.label else []
            label_list.insert(index, layer.label + str(layer.d))
            self.label = '-'.join(label_list)
        self._clear_TR()
    
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
        self._clear_TR()

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

    def calculate_TR(self):
        '''
        Calculate the reflectance and transmittance of the structure with the
        transfer matrix method from calc.py.
        '''
        if not self._layers_list:
            raise EmptyStructureException("Structure is empty.")

        # Room for np.ndarray optimization in the following code
        
        n_list = []    # Avoid side effect on Layer.n in the MultiLayer
        # interpolate seems costly. Unsure if obvious optimization exists.
        for L in self._layers_list:
            L_wl = L.wl * 1e3 if self.unit == 'nm' and L.unit == 'micron' else\
                   L.wl / 1e3 if self.unit == 'micron' and L.unit == 'nm' else\
                   L.wl
            layer_n = calc.interpolate(L_wl, L.n, self.wl)
            n_list.append(layer_n)
            
        wl = self.wl / 1e9 if self.unit == 'nm' else self.wl / 1e6
        k0 = 2 * np.pi / wl
        temp_n = [layer.n for layer in self._layers_list[::-1]]
        n = zip(*n_list[::-1])
        d = [layer.d / 1e9 if self.unit == 'nm' else layer.d / 1e6 \
             for layer in self._layers_list[::-1]]
        wl_len = len(self.wl)
        self.T = np.empty(wl_len)
        self.R = np.empty(wl_len)
        for i in range(wl_len):
            self.T[i] = calc.transmittance(k0[i], n[i], d,
                                           n0=self.n0, ns=self.ns)
            self.R[i] = calc.reflectance(k0[i], n[i], d,
                                         n0=self.n0, ns=self.ns)
        self.A = 1 - self.T - self.R
        self._TR_calculated = True

    def _clear_TR(self):
        self.T = None
        self.R = None
        self.A = None
        self._TR_calculated = False


class DataFormatException(Exception):
    pass


class EmptyStructureException(Exception):
    pass
