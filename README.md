# opt_sim
Transfer matrix simulation of 1D optical multilayers

Dependecies: numpy, matplotlib

To use the library, the python file that imports it should reside in the directory that contains the package opt_sim.

Structure:

```
.
opt_sim/
	__init__.py
	calc.py
	nklib.py
	plot.py
	structure.py
	nklib_data\
		nk data text files
		DLC data extrapolation figures
```

Example code and brief documentation for usage of the library.

```
# opt_sim is a package that includes several modules. The call signatures in
# the code below shows the organization.
import opt_sim as opt

# Currently the materials available is
# Ag, Al, BK7, DLC3W, DLC5W, DLC10W, DLC15W, DLC20W, DLC40W, DLC60W, ITO

# Input thickness for each layer. Default unit is nm
Ag7 = opt.nklib.Ag(7)
DLC50 = opt.nklib.DLC5W(50)
# An example of using micron as unit
BK = opt.nklib.BK7(5000, unit='micron')

# Create a structure by making a list of Layers. Order is from bottom to top.
AgSCC = opt.structure.MultiLayer([DLC50, Ag7, DLC50])
plain_glass = opt.structure.MultiLayer([BK], label="plain glass")

# A structure with an explicit substrate can be made, noting to change
# the refractive index of the exiting mdeium ns to 1.0 from its default
# value of 1.5
AgSCC_glass = opt.structure.MultiLayer([BK, DLC50, Ag7, DLC50],
                                        ns=1.0, label="on glass")

# A graphical representation of the MultiLayer can be shown on the console
print AgSCC_glass

# TR plotting can take either a single MultiLayer or a list of MultiLayers
# as its argument. By default both TR are shown. With the curves keyword
# argument, a combination of T, R, and A can be plotted.
opt.plot.TR([AgSCC, AgSCC_glass])
opt.plot.TR(plain_glass, curves='TA')

# nk plotting is similar to TR plotting, switching out MultiLayers for Layers
# as its argument. Custom minimum and maximum wavelength value can be given
# for both plotting functions.
opt.plot.nk([DLC50, BK], curves='k', min_wl=500)

opt.plot.show()


# Example using least verbose call signatures by using import *
from opt_sim.structure import MultiLayer as ML
from opt_sim import *

Ag7 = Ag(7)
ITO30 = ITO(30)
ITO60 = ITO(60)
struct_AgITO = ML([ITO30, Ag7, ITO60, Ag7, ITO30])
print struct_AgITO

# If one material is more restricted in its data wavelength range, it can
# be explicitly examined with the wl_by_layer attribute of a MultiLayer
print struct_AgITO.wl_by_layer

TR(struct_AgITO)
show()


# Recommended call signatures
from opt_sim.structure import MultiLayer as ML
from opt_sim import nklib, plot

# Example making a DLC layer with a refractive index gradient. In opt_sim.nklib
# exists DLC_list that contains all 7 material types of DLC, ordering from
# DLC3W to DLC60W. The extended keyword argument refers to the extrapolation
# performed on the experimental DLC data to extend it from the visible range
# to the near infrared range. The option to not extrapolate the data is given.
# See \opt_sim\nklib_data for figures related to the extrapolation.
DLC_gradient = [DLC(20, extended=False) for DLC in nklib.DLC_list[::-1]]
struct_DLC = ML(DLC_gradient, label="DLC gradient 140nm")

# n and k can be plotted in separate windows if wished.
plot.nk(DLC_gradient, sep=True)
plot.TR(struct_DLC)
plot.show()

# Alternatively, the user can use matplotlib.pyplot (already imported as plt)
# to perform any plotting they wish. However, an explicit call to calculate
# T and R of the MultiLayer has to be called first.
struct = ML([Ag7, DLC50, Ag7])
struct.calculate_TR()
print struct.wl[::30]
print struct.T[::30]
print struct.R[::30]
# Code for plotting goes here...
```