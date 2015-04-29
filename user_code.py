from nklib import *
from multilayer import *

import opt_sim as opt

Ag7 = opt.nklib.Ag(7)
DLC50 = opt.nklib.DLC5W(50)

ml = opt.MultiLayer([DLC50, Ag7, DLC50])

opt.plot.TRspectrum(ml)
opt.plot.color(ml)
opt.plot.show()

print ml.wl_range
