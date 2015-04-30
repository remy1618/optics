import numpy as np

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
                'A range of wavelengths should be given as an array, list, '
                'or tuple.')
        if not isinstance(n, (np.ndarray, list, tuple)):
            raise DataFormatException(
                'Wavelength-dependent refractive index should be given as '
                'an array, list, or tuple.')
        if len(wl) != len(n):
            print len(wl), len(n)
            raise DataFormatException(
                'wl and n should have the same length of data.')
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

    def __init__(self, layer_list=[], unit='nm', label=None):
        # Allow for the constructor call of MultiLayer(some_layer)
        # instead of MultiLayer([some_layer])
        if isinstance(layer_list, Layer):
            layer_list = [layer_list]
        if not unit in ['nm', 'micron']:
            raise DataFormatException(
                "unit should be given as 'nm' or 'micron'")
        self._layers_list = []
        self.wl_range = ''
        self._calculated = False
        self.unit = unit
        self.wl = None
        if not label:
            self.label = ''
            self._user_label = False
        else:
            self.label = label
            self._user_label = True
        if layer_list:
            if not all([layer_list[0].unit == layer_list[i].unit for i in range(len(layer_list))]):
                raise DataFormatException # or make conversion
            for layer in layer_list:
                self.add_layer(layer)
                self.wl_range += layer.label + ' ' + str(layer.wl[0]) + ' - ' \
                                 + str(layer.wl[-1]) + ' ' + layer.unit + '\n'
        self.T = None   #possibly do property: user try to index None
        self.R = None   #possibly do property: user try to index None
        #raise CalculateException in R/T_getter if not self._calculated
        

    def add_layer(self, layer_list, index=None):
        '''
        Add a layer to the top of multilayer if index is not given. Otherwise
        insert a layer indexed from 0.
        '''
        if isinstance(layer_list, Layer):
            layer_list = [layer_list]
        for layer in layer_list:
            if self.unit != layer.unit:
                raise DataFormatException(
                    'Unit mismatch between existing multilayer and new layers.')
        if not self._layers_list:
            self.wl = layer.wl
        if not all([np.allclose(self.wl, layer.wl) for layer in layer_list]):
            raise DataFormatException(
                'Wavelength range or interval of the new layer does not '
                'match that of the existing multilayer.')
        if isinstance(layer_list, Layer):
            layer_list = [layer_list]
        if index is None:
            for layer in layer_list:
                self._layers_list.append(layer)
                if not self._user_label:
                    if self.label:
                        self.label += '-'
                    self.label += layer.label + str(layer.d)
        else:
            for i, layer in enumerate(layer_list):
                self._layers_list.insert(index + i, layer)
                label_list = self.label.split('-')
                label_list.insert(index + i, layer.label + str(layer.d))
                if not self._user_label:
                    self.label = '-'.join(label_list)
        self._clear_calculation()
    
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
        self._clear_calculation()

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

    def calculate_TR(self, n0=1.0, ns=1.5):
        '''
        Calculate the reflectance and transmittance of the structure with the
        transfer matrix method.
        '''
        if not self._layers_list:
            raise EmptyStructureException('Structure is empty.')
        self.T, self.R = calc.TR_spectrum(self, n0, ns)
        self.A = 1 - self.T - self.R
        self._calculated = True

    def _clear_calculation(self):
        self.T = None
        self.R = None
        self.A = None
        self._calculated = False


class DataFormatException(Exception):
    pass


class EmptyStructureException(Exception):
    pass
