# # An example of reading in a file
# 
#

using Pkg
Pkg.activate(".")
using SPECreader

# for this short example we'll read and plot

# First we'll read in the SPEC output file

speceq = SPECEquilibrium("testing/data/G3V01L0Fi.002.sp.h5")


# We can get the Fourier modes of the boundary of the equilibrium,

specbound = get_boundary(speceq)

# note that we also have the ability to read the boundary directly from the file if that's all we're interested in,

ReadBoundary("testing/data/G3V01L0Fi.002.sp.h5")

# One can also get the Fourier components of the axis,

specaxis = get_axis(speceq)

# To get RZ coordinates anywhere there is also a `get_RZ` function, with the last input being the volume number for the SPEC equilibrium.

boundary = [get_RZ(1.0,θ,ζ,speceq,1) for θ in 0:2π/100:2π, ζ in 0:2π/100:2π]

# The Poincar\'e data is not loaded into SPECEquilibrium object and must be loaded seperately,

poincare = ReadPoincare("testing/data/G3V01L0Fi.002.sp.h5")