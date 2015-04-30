# Author: Remy Ko
# APD Group - Kherani

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
            raise DataFormatException(
                'wl and n should have the same lengths of data.')
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
                self.wl_range += layer.label + ' ' + layer.wl[0] + ' - ' + \
                                 layer.wl[-1] + ' ' + layer.unit + '\n'
        self.R = None   #possibly do property: user try to index None
        self.T = None   #possibly do property: user try to index None
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

    def view(self):
        '''
        Return a graphical representation of the structure.
        '''
        if not self._layers_list:
            return ''
        separator = '\n' + '-' * max([len(layer.view())
                                      for layer in self._layers_list])
        structure_view = '\nTop'
        for layer in self._layers_list[::-1]:
            structure_view += separator + '\n' + layer.view()
        structure_view += separator + '\n'
        print structure_view

    def calculate(self):
        #Make this automatic instead of getting the user to start it
        #(may depend on context)
        '''
        Calculate the reflectance and transmittance of the structure with the
        transfer matrix method.
        '''
        if not self._layers_list:
            raise EmptyStructureException('No layer exists.')
        if self.unit == 'nm':
            wl = self.wl / 1e9
        if self.unit == 'micron':
            wl = self.wl / 1e6
        k0 = 2 * pi / wl
        temp_n = [layer.n for layer in self._layers_list[::-1]]
        n = zip(*temp_n)
        d = [layer.d / 1e9 if self.unit == 'nm' else \
             layer.d / 1e6 if self.unit == 'micron' else None \
             for layer in self._layers_list[::-1]]
        wl_len = len(self.wl)
        self.R = np.empty(wl_len)
        self.T = np.empty(wl_len)
        for i in range(wl_len):
            self.R[i] = LayerStructure._R(k0[i], n[i], d)
            self.T[i] = LayerStructure._T(k0[i], n[i], d)
        self._calculated = True

    def _clear_calculation(self):
        self.R = None
        self.T = None
        self._calculated = False

    @staticmethod
    def _R(k0, n, d):
        return LayerStructure._RT(k0, n, d)[0]

    @staticmethod
    def _T(k0, n, d):
        return LayerStructure._RT(k0, n, d)[1]

    @staticmethod
    def _M(k0, n, d):
        assert isinstance(n, (complex, float, int)), "n is not a number"
        assert isinstance(d, (complex, float, int)), "d is not a number"
        Y = n * sqrt(eps0 / mu0)
        return np.array([
                         [cos(k0 * n * d), 1j * sin(k0 * n * d) / Y],
                         [1j * Y * sin(k0 * n * d), cos(k0 * n * d)]
                        ])

    @staticmethod
    def _RT(k0, n, d):
        '''For a single wavelength.
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
            m = np.dot(m, LayerStructure._M(k0, n[i], d[i]))
        r = (Y0 * m[0,0] + Y0 * Ys * m[0,1] - m[1,0] - Ys * m[1,1]) / \
            (Y0 * m[0,0] + Y0 * Ys * m[0,1] + m[1,0] + Ys * m[1,1])
        R = (np.absolute(r))**2
        t = 2 * Y0 / (Y0 * m[0,0] + Y0 * Ys * m[0,1] + m[1,0] + Ys * m[1,1])
        T = (np.absolute(t))**2 * Ys / Y0
        return R, T
