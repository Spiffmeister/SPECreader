using Pkg
Pkg.activate(".")
using SPECreader

SPEC = SPECreader.read_spec("testing\\G3V01L0Fi.002.sp.h5")