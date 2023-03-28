# SPECreader

Julia version of [SPEC-field-reader](https://github.com/zhisong/SPEC-field-reader) or [pyoculus](https://github.com/zhisong/pyoculus) both by [Zhisong Qu](https://github.com/zhisong/) for reading an output from the [Stepped Pressure Equilibrium Code](https://github.com/PrincetonUniversity/SPEC).



## Example

Given the output file `G1V02L0Fi.001.sp.h5` it can be read by,

```julia
using SPECreader

SPEC = SPECreader.read_spec("G1V02L0Fi.001.sp.h5")
```




## Requirements

- `HDF5.jl`

