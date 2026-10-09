# KP505: one blade + hub (STEP solid, full scale, mm) -> 5-bladed propeller, model scale in metres,
# surface triangulation as ASCII STL for REEF3D (6DOF X 180).
#   python convert.py <step> <out.stl> [hmax_mm_model] [nelem_per_2pi]
import gmsh, sys, math, time
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

step, out = sys.argv[1], sys.argv[2]
hmax = float(sys.argv[3]) if len(sys.argv)>3 else 2.0       # max triangle size, model scale [mm]
ncurv = float(sys.argv[4]) if len(sys.argv)>4 else 40.0     # elements per 2 pi of curvature
Z, scale = 5, 1.0/31.6/1000.0                               # full scale mm -> model scale m (D 7.9 m -> 0.25 m)

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("kp505")
v = gmsh.model.occ.importShapes(step)
gmsh.model.occ.synchronize()
vol = [e for e in v if e[0]==3]
log("imported volumes", vol, "mass (full scale mm^3)", [gmsh.model.occ.getMass(*e) for e in vol])

# copies rotated about the shaft (x axis through the origin) by 72 deg
# the solid is the whole hub with one blade: the rotated copies minus the original leave their
# blades only (the hubs coincide), which are then fused with the original
blades = []
for k in range(1, Z):
    c = gmsh.model.occ.copy(vol)
    gmsh.model.occ.rotate(c, 0, 0, 0, 1, 0, 0, 2.0*math.pi*k/Z)
    b, _ = gmsh.model.occ.cut(c, vol, removeObject=True, removeTool=False)
    log("blade", k, b, "volume", sum(gmsh.model.occ.getMass(*e) for e in b))
    blades += b
fused, _ = gmsh.model.occ.fuse(vol, blades)
log("fused", fused)
gmsh.model.occ.dilate(fused, 0, 0, 0, scale, scale, scale)
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

vols = gmsh.model.getEntities(3)
mass = sum(gmsh.model.occ.getMass(*e) for e in vols)
bb = gmsh.model.getBoundingBox(-1, -1)
log("fused volumes", len(vols), "surfaces", len(gmsh.model.getEntities(2)), "volume [m^3]", mass)
log("bbox [m]", [round(x, 5) for x in bb])

# curvature-based surface mesh
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", ncurv)
gmsh.option.setNumber("Mesh.MeshSizeMax", hmax*1.0e-3)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.3*hmax*1.0e-3)
gmsh.option.setNumber("Mesh.Algorithm", 6)
log("meshing")
gmsh.model.mesh.generate(2)
log("meshed")

# ASCII STL
gmsh.option.setNumber("Mesh.Binary", 0)
gmsh.option.setNumber("Mesh.StlOneSolidPerSurface", 0)
gmsh.write(out)

ntri = 0
for (d, t) in gmsh.model.getEntities(2):
    types, tags, nodes = gmsh.model.mesh.getElements(d, t)
    for ty, tg in zip(types, tags):
        if ty == 2: ntri += len(tg)
log("triangles", ntri)
gmsh.finalize()
