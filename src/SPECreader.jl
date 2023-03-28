"""
    SPECreader

Julia version of the SPEC-field-reader by Zhisong Qu https://github.com/zhisong/SPEC-field-reader
"""
module SPECreader

    using HDF5
    using LinearAlgebra
    # using ZernikePolynomials
    # using FastChebInterp

    # Data Structures
    include("DataStructs.jl")
    # reading in HDF5
    include("IO.jl")
    # Metrics
    include("Geometry.jl")
    # Fields


    export read_spec
end # module temp
