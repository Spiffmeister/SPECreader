"""
    SPECreader

Julia version of the SPEC-field-reader by Zhisong Qu https://github.com/zhisong/SPEC-field-reader
"""
module SPECreader

    using HDF5
    using LinearAlgebra
    using StaticArrays
    using NonlinearSolve
    using Optim

    # Data Structures
    include("DataStructs.jl")
    # reading in HDF5
    include("BasisFunctions.jl")
    # Metrics
    include("Geometry.jl")
    include("Fields.jl")
    # Fields
    # Basis functions for SPEC geometry


    export SPECEquilibrium
    export ReadBoundary, ReadPoincare
    
    export get_boundary, get_axis
    export get_RZ, find_sθζ

    export get_Bfield
    export field_line!

end