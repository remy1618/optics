import opt_sim as opt

Ag7 = opt.nklib.Ag(7)
DLC50 = opt.nklib.DLC5W(50)

struct = opt.structure.MultiLayer([DLC50, Ag7, DLC50])
print struct

opt.plot.TR(struct)
opt.plot.nk(struct)
opt.plot.show()
