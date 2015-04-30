# Author: Remy Ko
# APD Group - Kherani

# Work in progress
# Last updated Nov 21, 2014

from numpy import pi, sqrt, sin, cos
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

n0 = 1.0
eps0 = 8.85e-12
mu0 = 4e-7 * pi
ns = 1.5

glass_1mm = np.loadtxt("glass_T_out.txt")
glass_6mm = np.loadtxt("clear_6_out.txt")
PET_3mil = np.loadtxt("PET_T_out.txt")


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

    def view(self):
        '''
        Return a description of the layer.
        '''
        return '{} {} {}'.format(self.label, self.d, self.unit)
        

class LayerStructure:
    '''
    Initialize a multilayer structure. The structure is empty by default,
    unless the layers created from the Layer class are given in a list. The
    layers are added onto or removed from the multilayer structure using
    methods add_layer(), remove_layer(), or replace_layer(). The reflectance
    and transmittance of the structure is calculated using the method
    calculate() which is implemented using the transfer matrix method.
    '''

    def __init__(self, layer_list=[], unit='nm', label=None):
        if isinstance(layer_list, Layer): #maybe keep this feature?
            layer_list = [layer_list]
        if not unit in ['nm', 'micron']:
            raise DataFormatException(
                "unit should be given as 'nm' or 'micron'")
        self._layers_list = []
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


def interpolate(wl_raw, f, wl):
    if wl[0] < wl_raw[0] or wl[-1] > wl_raw[-1]:
        raise BoundException("User-defined wavelength out of data range.")
    f_new = np.empty(len(wl), f.dtype)
    for i in range(len(wl)):    # does not work for case where wl range is out of bound of wl_raw range
        j = 0
        while wl[i] > wl_raw[j]:
            j += 1
        j_high = j
        j_low = j - 1
        # Linear interpolation y = y0 + (x-x0) * (y1-y0) / (x1-x0)
        value = f[j_low] + (wl[i] - wl_raw[j_low]) * (f[j_high] - f[j_low]) / (wl_raw[j_high] - wl_raw[j_low])
        f_new[i] = value
    return f_new


def plot(structure_list, min_wl=None, max_wl=None, curves='RT', title='Reflectance and Transmittance plots',
         substrate=None):
    '''
    Plot the reflectance and transmittance curve of given structures in the
    window. If given, min_wl and max_wl sets the wavelength range in the unit
    the LayerStructure instances are initialized in. Default is nm.
    Inputting 'R' or 'T' for the curves argument plots only the reflectance
    or transmittance, respectively. Default is both.
    '''
    if isinstance(structure_list, LayerStructure): #maybe keep this feature?
        structure_list = [structure_list]
    if not all([structure_list[0].unit == structure_list[i].unit \
                for i in range(len(structure_list))]):
        convert_to = input('Unit mismatch among structures. Enter 1 to '
                           'convert all structures to nm, 2 to convert all '
                           'structures to microns, or n to abort conversion.')
        #do conversion method in LayerStructure
    for structure in structure_list:
        if not structure._calculated:
            structure.calculate()

    unit = structure_list[0].unit
    plt.figure()
    for structure in structure_list:
        wl = _user_wl_range(structure, min_wl, max_wl)
        min_i, max_i = _user_wl_range_indices(structure, min_wl, max_wl)
        if 'R' in curves:
            R = structure.R[min_i:max_i]
            substrate_R = 0
            # gives wrong plot right now
            # if substrate ==  '1mm glass':
            #     substrate_R = 1 - interpolate(glass_1mm, wl) #ignores absorption
            # elif substrate == '6mm glass':
            #     substrate_R = 1 - interpolate(glass_6mm, wl)
            # elif substrate == '3mil PET':
            #     substrate_R = 1 - interpolate(PET_3mil, wl)
            R = R + substrate_R
            plt.plot(wl, R, label=structure.label + ' R')
        if 'T' in curves:
            T = structure.T[min_i:max_i]
            substrate_T = 1
            if substrate == '1mm glass':
                substrate_T = interpolate(glass_1mm[:,0], glass_1mm[:,1], wl)
            elif substrate == '6mm glass':
                substrate_T = interpolate(glass_6mm[:,0], glass_6mm[:,1], wl)
            elif substrate == '3mil PET':
                substrate_T = interpolate(PET_3mil[:,0], PET_3mil[:,1], wl)
            T = T * substrate_T
            plt.plot(wl, T, label=structure.label + ' T')
        if 'A' in curves:
            R = structure.R[min_i:max_i]
            T = structure.T[min_i:max_i]
            A = 1 - R - T
            plt.plot(wl, A, label=structure.label + ' A')
    plt.title(title)
    plt.xlabel('Wavelength ({})'.format(unit))
    if min_wl:
        plt.xlim(xmin=min_wl)
    if max_wl:
        plt.xlim(xmax=max_wl)
    plt.ylim(0, 1)
    plt.legend(loc=0)

def plot_n(layer_list, min_wl=None, max_wl=None, curves='nk'):
    '''
    Plot the complex refractive index of the given layers in the same window.
    If given, min_wl and max_wl sets the wavelength range in the unit
    the Layer instances are initialized in. Default is nm.
    '''
    if isinstance(layer_list, Layer): #maybe keep this feature?
        layer_list = [layer_list]
    if not all([layer_list[0].unit == layer_list[i].unit \
                for i in range(len(layer_list))]):
        convert_to = input('Unit mismatch among layers. Enter 1 to convert '
                           'all layers to nm, 2 to convert all layers to '
                           'microns, or n to abort conversion.')
        #do conversion method in LayerStructure

    unit = layer_list[0].unit
    plt.figure()
    for layer in layer_list:
        wl = _user_wl_range(layer, min_wl, max_wl)
        min_i, max_i = _user_wl_range_indices(layer, min_wl, max_wl)
        n = layer.n.real[min_i:max_i]
        k = -layer.n.imag[min_i:max_i]
        if 'n' in curves:
            plt.plot(wl, n, label=layer.label + ' n')
        if 'k' in curves:
            plt.plot(wl, k, label=layer.label + ' k')
    plt.title('Refractive index plots')
    plt.xlabel('Wavelength ({})'.format(unit))
    plt.legend(loc=0)
    
def show():
    plt.show()
    
def _user_wl_range(structure, min_wl, max_wl):
    return _custom_wl_range(structure, min_wl, max_wl)[0]

def _user_wl_range_indices(structure, min_wl, max_wl):
    return _custom_wl_range(structure, min_wl, max_wl)[1], \
           _custom_wl_range(structure, min_wl, max_wl)[2]

def _custom_wl_range(structure, min_wl, max_wl):
    '''
    Allow for user custom wavelength range.
    '''
    if not (isinstance(min_wl, (int, float)) and min_wl >= 0) and not \
       min_wl is None or \
       not (isinstance(max_wl, (int, float)) and max_wl >= 0) and not \
       max_wl is None or \
       max_wl < min_wl and not max_wl is None:  #maybe cleanup if possible
        raise BoundException('The given wavelength range is invalid.')

    min_index, index = 0, 0
    while index < len(structure.wl) and structure.wl[index] <= max_wl or \
          index < len(structure.wl) and not max_wl:
        if structure.wl[index] < min_wl:
            min_index = index + 1
        index += 1
    if not max_wl:
        index = None
    max_index = index
    return structure.wl[min_index:max_index], min_index, max_index


class BoundException(Exception):
    pass


class EmptyStructureException(Exception):
    pass


class DataFormatException(Exception):
    pass
