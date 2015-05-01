import matplotlib.pyplot as plt
import structure as st
import calc

def TR(structure_list, min_wl=None, max_wl=None, curves='TR'):
    '''
    Plot the reflectance and transmittance curve of given structures in the
    window. If given, min_wl and max_wl sets the wavelength range in the unit
    the LayerStructure instances are initialized in. Default is nm.
    Inputting 'R' or 'T' for the curves argument plots only the reflectance
    or transmittance, respectively. Default is both.
    '''
    if isinstance(structure_list, st.MultiLayer):
        structure_list = [structure_list]
    if not all([structure_list[0].unit == structure_list[i].unit \
                for i in range(len(structure_list))]):
        convert_to = input('Unit mismatch among structures. Enter 1 to '
                           'convert all structures to nm, 2 to convert all '
                           'structures to microns, or n to abort conversion.')
        #do conversion method in st.MultiLayer
    for structure in structure_list:
        if not structure._calculated:
            structure.calculate_TR()

    unit = structure_list[0].unit
    plt.figure()
    for structure in structure_list:
        wl = structure.wl
        R = structure.R
        T = structure.T
        if 'R' in curves:
            plt.plot(wl, R, label=structure.label + ' R')
        if 'T' in curves:
            plt.plot(wl, T, label=structure.label + ' T')
        if 'A' in curves:
            A = 1 - R - T
            plt.plot(wl, A, label=structure.label + ' A')
    plt.title('Reflectance and Transmittance plots')
    plt.xlabel('Wavelength ({})'.format(unit))
    if min_wl:
        plt.xlim(xmin=min_wl)
    if max_wl:
        plt.xlim(xmax=max_wl)
    plt.ylim(0, 1)
    plt.legend(loc=0)

def nk(layer_list, min_wl=None, max_wl=None, curves='nk'):
    '''
    Plot the complex refractive index of the given layers in the same window.
    If given, min_wl and max_wl sets the wavelength range in the unit
    the Layer instances are initialized in. Default is nm.
    '''
    if isinstance(layer_list, st.Layer):
        layer_list = [layer_list]
    if isinstance(layer_list, st.MultiLayer):
        layer_list = layer_list._layers_list
    if not all([layer_list[0].unit == layer_list[i].unit \
                for i in range(len(layer_list))]):
        convert_to = input('Unit mismatch among layers. Enter 1 to convert '
                           'all layers to nm, 2 to convert all layers to '
                           'microns, or n to abort conversion.')
        #do conversion method in st.MultiLayer

    unit = layer_list[0].unit
    plt.figure()
    for layer in layer_list:
        wl = layer.wl
        n = layer.n.real
        k = -layer.n.imag
        if 'n' in curves:
            plt.plot(wl, n, label=layer.label + ' n')
        if 'k' in curves:
            plt.plot(wl, k, label=layer.label + ' k')
    plt.title('Refractive index plots')
    plt.xlabel('Wavelength ({})'.format(unit))
    if min_wl:
        plt.xlim(xmin=min_wl)
    if max_wl:
        plt.xlim(xmax=max_wl)
    plt.legend(loc=0)
    
def show():
    plt.show()
