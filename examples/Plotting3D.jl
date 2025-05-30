using Pkg
Pkg.activate(".")
using SPECreader
using LinearAlgebra
using GLMakie, GeometryBasics




speceq = SPECEquilibrium("testing/data/G3V01L0Fi.002.sp.h5")



θ_rng = range(0.0,2π,101)
ζ_rng = range(0.0,2π,1001)

boundary_points_x = [get_RZ(1.0,θ,ζ,speceq,1)[1]*cos(ζ) for θ in θ_rng, ζ in ζ_rng]
boundary_points_z = [get_RZ(1.0,θ,ζ,speceq,1)[2] for θ in θ_rng, ζ in ζ_rng]
boundary_points_y = [norm(get_RZ(1.0,θ,ζ,speceq,1))*sin(ζ) for θ in θ_rng, ζ in ζ_rng]

# Convert to Point3fs
points = vec([Point3f(x,y,z) for (x,y,z) in zip(boundary_points_x,boundary_points_y,boundary_points_z)])

# Build faces for the plotting
_faces = decompose(QuadFace{GLIndex}, Tessellation(Rect(0, 0, 1, 1), size(boundary_points_z)))

# Compute the normals for the points
_normals = normalize.(points)

modB = [norm(get_Bfield(1.0,θ,ζ,speceq)) for θ in θ_rng, ζ in ζ_rng]

# Assign the colours to the faces
_color = FaceView(modB[:], _faces)

# Build the mesh from the stuff above
solidmesh = GeometryBasics.mesh(points, _faces, normal = _normals, color = _color)







f = Figure(size=(1200,600))
lscene = LScene(f[1,1], show_axis=false)
axf = lscene.scene
cam3d!(axf)

mesh!(axf,solidmesh, transparency=true)

