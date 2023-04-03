using Pkg
Pkg.activate(".")
using SPECreader

# SPEC = SPECreader.read_spec("testing\\G3V01L0Fi.002.sp.h5")
SPEC = SPECreader.read_spec("testing\\G1V02L0Fi.001.sp.h5")


SPECreader.GetMetric(SPEC,0.3,0.4,0.5,1)