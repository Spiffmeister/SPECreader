module SPECreader

    using HDF5
    # using ZernikePolynomials
    # using FastChebInterp

    # Data Structures
    include("DataStructs.jl")
    # reading in HDF5
    include("IO.jl")

    export read_spec
end # module temp
