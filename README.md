# SPECreader

This package is for reading and manipulating output from the [Stepped Pressure Equilibrium Code](https://github.com/PrincetonUniversity/SPEC) (SPEC).



## Example

Given the output file `G1V02L0Fi.001.sp.h5` it can be read by,

```julia
using Pkg
Pkg.activate(".")
using SPECreader

speceq = SPECEquilibrium("testing/data/G1V02L0Fi.001.sp.h5")
```







## Plotting


For plotting one could use either [Plots](https://docs.juliaplots.org/stable/) or [Makie](https://docs.makie.org/stable/). In the examples folder for plotting a mesh in 3D.




