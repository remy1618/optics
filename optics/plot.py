import matplotlib.pyplot as plt
from skimage.io import imread
from skimage.color import rgb2hsv, hsv2rgb
import matplotlib.image as mpimg
from matplotlib.figure import figaspect
import numpy as np
from . import structure as st
import os.path
import warnings
# Append path to support files
dir_path = os.path.dirname(os.path.realpath(__file__))
sup_path = os.path.join(dir_path, "support_files")

def TR(structure_list, min_wl=None, max_wl=None, curves='TR', unit='nm',
       show_solar=False, legend=True, new_figure=True, grid=False):
    '''
    Plot the reflectance and transmittance curve of given structures in the
    window. If given, min_wl and max_wl sets the wavelength range in the unit
    the LayerStructure instances are initialized in. Default is nm.
    Inputting 'R' or 'T' for the curves argument plots only the reflectance
    or transmittance, respectively. Default is both.
    '''
    # Alternate call signature for convenience
    if isinstance(structure_list, st.MultiLayer):
        structure_list = [structure_list]
        
    if not unit in ['nm', 'micron']:
            raise DataFormatException(
                "unit should be given as 'nm' or 'micron'")
    for structure in structure_list:
        if not structure._TR_calculated:
            structure.calculate_TR()
    if not structure_list:
        return

    st_min_wl = None
    st_max_wl = None
    if new_figure:
        plt.figure()
    if show_solar:
        AM1p5_data = np.loadtxt(os.path.join(sup_path, "ASTMG173.txt"), skiprows=2)
        solar_wl = AM1p5_data[:,0]
        solar_intens = AM1p5_data[:,3]
        solar_intens /= max(solar_intens)
        plt.plot(solar_wl, solar_intens, label="AM1.5", color=(0.55,0.55,0.55))
    for structure in structure_list:
        wl = structure.wl * 1e3 if unit == 'nm' and structure.unit == 'micron' else\
             structure.wl / 1e3 if unit == 'micron' and structure.unit == 'nm' else\
             structure.wl
        T = structure.T
        R = structure.R
        if 'T' in curves:
            plt.plot(wl, T, label=structure.label + ' T')
        if 'R' in curves:
            plt.plot(wl, R, label=structure.label + ' R')
        if 'A' in curves:
            A = 1 - R - T
            plt.plot(wl, A, label=structure.label + ' A')
        st_min_wl = min(st_min_wl, wl[0]) if st_min_wl else wl[0]
        st_max_wl = max(st_max_wl, wl[-1]) if st_max_wl else wl[-1]
    
    plt.title('Reflectance and Transmittance plots')
    plt.xlabel('Wavelength ({})'.format(unit))
    if min_wl:
        plt.xlim(xmin=min_wl)
    else:
        plt.xlim(xmin=st_min_wl)
    if max_wl:
        plt.xlim(xmax=max_wl)
    else:
        plt.xlim(xmax=st_max_wl)
    plt.ylim(0, 1)
    if legend:
        plt.legend(loc=0)
    if grid:
        plt.grid(True)

def nk(layer_list, min_wl=None, max_wl=None, curves='nk', sep=False, unit='nm',
       grid=False):
    '''
    Plot the complex refractive index of the given layers in the same window.
    If given, min_wl and max_wl sets the wavelength range in the unit
    the Layer instances are initialized in. Default is nm.
    '''
    # Alternate call signature for convenience
    if isinstance(layer_list, st.Layer):
        layer_list = [layer_list]
    if isinstance(layer_list, st.MultiLayer):
        layer_list = layer_list._layers_list
    if not unit in ['nm', 'micron']:
            raise DataFormatException(
                "unit should be given as 'nm' or 'micron'")
    unit = layer_list[0].unit
            
    # Reduce redundant code merging everything below
    if curves == "nk" and sep:
        material_list = []
        plt.figure()
        plt.title("Real refractive index plot")
        plt.xlabel("Wavelength ({})".format(unit))
        if min_wl:
            plt.xlim(xmin=min_wl)
        if max_wl:
            plt.xlim(xmax=max_wl)
        for layer in layer_list:
            if layer.label in material_list:
                continue
            material_list.append(layer.label)
            wl = layer.wl * 1e3 if unit == 'nm' and layer.unit == 'micron' else\
                 layer.wl / 1e3 if unit == 'micron' and layer.unit == 'nm' else\
                 layer.wl
            n = layer.n.real
            plt.plot(wl, n, label=layer.label + ' n')
            if grid:
                plt.grid(True)
        plt.legend(loc=0)

        material_list = []
        plt.figure()
        plt.title("Imaginary refractive index plot")
        plt.xlabel("Wavelength ({})".format(unit))
        if min_wl:
            plt.xlim(xmin=min_wl)
        if max_wl:
            plt.xlim(xmax=max_wl)
        for layer in layer_list:
            if layer.label in material_list:
                continue
            material_list.append(layer.label)
            wl = layer.wl * 1e3 if unit == 'nm' and layer.unit == 'micron' else\
                 layer.wl / 1e3 if unit == 'micron' and layer.unit == 'nm' else\
                 layer.wl
            k = -layer.n.imag
            plt.plot(wl, k, label=layer.label + ' k')
            if grid:
                plt.grid(True)
        plt.legend(loc=0)
        return None

    material_list = []
    plt.figure()
    for layer in layer_list:
        if layer.label in material_list:
            continue
        material_list.append(layer.label)
        wl = layer.wl * 1e3 if unit == 'nm' and layer.unit == 'micron' else\
             layer.wl / 1e3 if unit == 'micron' and layer.unit == 'nm' else\
             layer.wl
        n = layer.n.real
        k = -layer.n.imag
        if 'n' in curves:
            plt.plot(wl, n, label=layer.label + ' n')
        if 'k' in curves:
            plt.plot(wl, k, label=layer.label + ' k')
        if grid:
                plt.grid(True)
    plt.title('Refractive index plots')
    plt.xlabel('Wavelength ({})'.format(unit))
    if min_wl:
        plt.xlim(xmin=min_wl)
    if max_wl:
        plt.xlim(xmax=max_wl)
    plt.legend(loc=0)   # this gives error at the moment
    return None

def view(structure, outdoor_lux=109870., indoor_lux=500.,
         outdoor_image=False, indoor_image=False, show_original=False):
    '''
    Default outdoor lux is for AM1.5 and indoor lux for office lighting.
    See https://en.wikipedia.org/wiki/Daylight and https://en.wikipedia.org/wiki/Lux
    '''
    if not structure._color_calculated:
        structure.calculate_color()

    if not outdoor_image:
        T_image_original = imread(os.path.join(sup_path, "test_outdoor.png")) / 256.
    else:
        T_image_original = imread(outdoor_image) / 256.
    if not indoor_image:
        R_image_original = imread(os.path.join(sup_path, "test_indoor.png")) / 256.
    else:
        R_image_original = imread(indoor_image) / 256.
    if not T_image_original.shape == R_image_original.shape:
        T_width, T_height = T_image_original.shape[0], T_image_original.shape[1]
        R_width, R_height = R_image_original.shape[0], R_image_original.shape[1]
        warnings.simplefilter("always")
        warnings.warn("Image dimensions do not match ({}x{} vs {}x{}). Images are \
cropped around their centres.".format(T_width, T_height, R_width, R_height),
Warning, stacklevel=2)
        width = min(T_width, R_width)
        height = min(T_height, R_height)
        T_image_original = T_image_original[(T_width-width)/2:(T_width+width)/2,
                                            (T_height-height)/2:(T_height+height)/2,:]
        R_image_original = R_image_original[(R_width-width)/2:(R_width+width)/2,
                                            (R_height-height)/2:(R_height+height)/2,:]

    # Assume outdoor photo taken at 109870 lux and indoor photo at 500 lux
    # Scaling for the perception of brightness is from Steven's Law
    T_brightness_ratio = np.cbrt(outdoor_lux / 109870.)
    R_brightness_ratio = np.cbrt(indoor_lux / 500.)
    T_image_hsv = rgb2hsv(T_image_original)
    R_image_hsv = rgb2hsv(R_image_original)
    T_image_hsv[:,:,2] *= T_brightness_ratio
    R_image_hsv[:,:,2] *= R_brightness_ratio
    T_image_hsv[T_image_hsv[:,:,2] > 1, 2] = 1.
    R_image_hsv[R_image_hsv[:,:,2] > 1, 2] = 1.
    T_image = hsv2rgb(T_image_hsv)
    R_image = hsv2rgb(R_image_hsv)

    T_filter = np.array(structure.T_color, float)
    R_filter = np.array(structure.R_color, float)
    T_image_coating = (T_image * T_filter)
    R_image_coating = (R_image * R_filter)
    overlay = T_image_coating + R_image_coating

    f, ax = plt.subplots(nrows=1, ncols=1)
    f.suptitle(structure.label)
    ax.set_title("Outdoor lux: {:d}    Indoor lux: {:d}".format(int(outdoor_lux), int(indoor_lux)))
    ax.imshow(overlay)
    ax.axis("off")

    if show_original:
        w, h = figaspect(2.)
        f2, ((ax1), (ax2)) = plt.subplots(nrows=2, ncols=1, sharex='col', figsize=(w,h))
        f2.subplots_adjust(left=0.03, right=0.97, hspace=0.0, wspace=0.0)
        ax1.axis("off")
        ax2.axis("off")
        ax1.set_aspect("equal")
        ax2.set_aspect("equal")
        ax1.imshow(T_image)
        ax2.imshow(R_image)    

def ellipsometry(structure_list, angle=51., min_eV=None, max_eV=None,
                 legend=True, grid=True):
    # Alternate call signature for convenience
    if isinstance(structure_list, st.MultiLayer):
        structure_list = [structure_list]
    for structure in structure_list:
        if not structure._tanPSIcosDELTA_calculated:
            structure.calculate_tanPSIcosDELTA(angle)
    if not structure_list:
        return

    st_min_eV = None
    st_max_eV = None
    plt.figure()
    plt.suptitle("Ellipsometry at {:.1f} degree incident angle".format(angle))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212, sharex=ax1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    for structure in structure_list:
        # Assumes only valid units are 'nm' and 'micron'
        E = 1240 / structure.wl[::-1] if structure.unit == 'nm' else\
            1.24 / structure.wl[::-1]
        tanPSI = structure.tanPSI[::-1]
        cosDELTA = structure.cosDELTA[::-1]
        st_min_eV = min(st_min_eV, E[0]) if st_min_eV else E[0]
        st_max_eV = min(st_max_eV, E[-1]) if st_max_eV else E[-1]

        ax1.set_title('tan(PSI)')
        ax1.plot(E, tanPSI, label=structure.label)
        ax2.set_title('cos(DELTA)')
        ax2.plot(E, cosDELTA, label=structure.label)

        ax2.set_xlabel('Energy (eV)')
        if min_eV:
            ax2.set_xlim(xmin=min_eV)
        else:
            ax2.set_xlim(xmin=st_min_eV)
        if max_eV:
            ax2.set_xlim(xmax=max_eV)
        else:
            ax2.set_xlim(xmax=st_max_eV)

        if legend:
            ax1.legend(loc=0)
            ax2.legend(loc=0)
        if grid:
            ax1.grid(True)
            ax2.grid(True)
        plt.tight_layout(rect=[0,0,1,0.95])

def show():
    plt.show()


class DataFormatException(Exception):
    pass
