
using Revise
using Test

using Pkg
Pkg.activate(".")
using SPECreader


speceq = SPECEquilibrium("testing/data/G3V01L0Fi.002.sp.h5")




@testset "Test 1 volume RZ" begin
    RZ = (6.086373127567621, 0.0)
    @test all(get_RZ(0.0,0.0,0.0,speceq,1) .≈ RZ)

    RZ = (5.790604383896138, 0.047181582748737)
    @test all(get_RZ(0.2,0.3,0.4,speceq,1) .≈ RZ)

    RZ = (5.692641457375362, -0.116354800624682)
    @test all(get_RZ(-0.2,-0.3,-0.4,speceq,1) .≈ RZ)

    RZ = (5.182286170111236, 0.342000104406030)
    @test all(get_RZ(0.5,Float64(pi),pi/2,speceq,1) .≈ RZ)
end



@testset "Test 1 volume B" begin
    B = (0.0, -0.197332224079666, 0.141464914530581)
    @test all(get_Bfield(0.0,0.0,0.0,speceq) .≈ B)

    B = (0.007233436722876, -0.113344143560185, 0.160087751181934)
    @test all(get_Bfield(0.2,0.3,0.4,speceq) .≈ B)

    B = (-0.009184900352870, -0.061835475866279, 0.114093277649795)
    @test all(get_Bfield(-0.2,-0.3,-0.4,speceq) .≈ B)

    B = (0.065129059280715, -0.255911380467064, 0.284864657593101)
    @test all(get_Bfield(0.5,Float64(pi),pi/2,speceq) .≈ B)
end


