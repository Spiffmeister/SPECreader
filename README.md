# SPECreader

<img alt="Static Badge" src="https://img.shields.io/badge/docs-dev-blue?label=docs&link=https%3A%2F%2Fspiffmeister.github.io%2FSPECreader%2Fdev%2F">

This package is for reading and manipulating output from the [Stepped Pressure Equilibrium Code](https://github.com/PrincetonUniversity/SPEC) (SPEC).

Currently the package is registered, so to add please use:
```julia
] add https://github.com/Spiffmeister/SPECreader.git
```


## Example

Given the output file `G1V02L0Fi.001.sp.h5` it can be read by,

```julia
using Pkg
Pkg.activate(".")
using SPECreader

speceq = SPECEquilibrium("testing/data/G1V02L0Fi.001.sp.h5")
```





## Plotting


For plotting one could use either [Plots](https://docs.juliaplots.org/stable/) or [Makie](https://docs.makie.org/stable/). The examples folder contains a script for plotting a mesh in 3D using Makie.




