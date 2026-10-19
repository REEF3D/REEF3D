# writes the regression cases of the horizontal diffusion of SFLOW on curvilinear grids, patches and blocks (P6a:
# A 212 1, A 260 3/4), of the Weber junction on three blocks (SFLOW, NHFLOW) and of the NHFLOW turbulence on the
# patches and the blocks (P6b: A 560): Regression/cases/sflow_2d_*, nhflow_3d_*
import os, sys, json, shutil, subprocess
T = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(T, "../../Regression/cases")

def fresh(d):
    """empty the case directory, keep its reference (regression.py bless)"""
    ref = os.path.join(d, "reference.json")
    keep = open(ref).read() if os.path.exists(ref) else None
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(d)
    if keep is not None:
        open(ref, "w").write(keep)

def derived(name, base, meta):
    d = os.path.join(C, name)
    os.makedirs(d, exist_ok=True)
    m = {"base": base}
    m.update(meta)
    json.dump(m, open(os.path.join(d, "case.json"), "w"), indent=1)

derived("sflow_2d_curv_channel_rot30_mixlen", "sflow_2d_curv_channel_rot30", {
    "ctrl_add": ["A 260 3", "W 2 0.01"],
    "tags": ["2d", "sflow", "curvilinear", "wetdry", "river", "turbulence"],
    "description": "As sflow_2d_curv_channel_rot30 (river reach with an emerged island, grid rotated by 30 deg) with the "
                   "horizontal diffusion (A 212 1) of W 2 and the mixing-length eddy viscosity (A 260 3): the stress "
                   "nu h (grad u + grad u^T) through the faces of the curvilinear cells; rotation-invariant to 1e-10 with "
                   "A 211 0 (sflow_curvilinear/t9).",
    "covers": ["sflow_ediff::curv_stress (orthogonal part, non-orthogonal part k.grad u, transposed part, no stress at "
               "walls and dry cells)", "sflow_ediff::algebraic_eddyv A 260 3 on the curvilinear grid (own_eddyv)"]})
derived("sflow_2d_mb_amr_parab", "sflow_2d_mb_amr", {
    "ctrl_add": ["A 260 4", "W 2 0.005"],
    "tags": ["2d", "sflow", "curvilinear", "multiblock", "amr", "turbulence"],
    "description": "As sflow_2d_mb_amr (four turned blocks around a vertex, two refinement levels) with the horizontal "
                   "diffusion of W 2 and the parabolic eddy viscosity (A 260 4) on every grid: the stress through the "
                   "faces of the blocks and of the patches from the ghost cells (both cells of a face the same flux).",
    "covers": ["sflow_amr: sflow_ediff of the patches and blocks with own_eddyv, the eddy viscosity of the first ghost "
               "layer", "sflow_f_logic: A 260 3/4 with REEFAMR (sflow_turb_void on level 0)"]})
derived("sflow_2d_mb_amr_mixlen_mpi2", "sflow_2d_mb_amr", {
    "np": 2, "ctrl_add": ["A 260 3"],
    "tags": ["2d", "sflow", "curvilinear", "multiblock", "amr", "turbulence", "mpi"],
    "description": "As sflow_2d_mb_amr with the mixing-length eddy viscosity (A 260 3) on 2 ranks: the velocity "
                   "gradients of the first ghost layer at the partition boundaries and in the blocks on other ranks.",
    "covers": ["sflow_ediff::algebraic_eddyv in the ghost cells of the partition and of the blocks"]})

# the Weber junction (tools/weber_junction): three blocks, W/dx 8, q* = 0.75, native outflow on block 0 (layout b),
# the parabolic eddy viscosity
for name, np_, extra in (("sflow_2d_mb_weber_junction", 1, []), ("sflow_2d_mb_weber_junction_mpi2", 2, [])):
    d = os.path.join(C, name)
    fresh(d)
    if np_ == 1:
        subprocess.run([sys.executable, os.path.join(T, "weber_case.py"), d, "--layout", "b", "--n", "8", "--np", "1",
                        "--tend", "1000", "--xctrl", "A 260 4"], check=True, stdout=subprocess.DEVNULL)
        os.remove(os.path.join(d, "blocks.json"))
        ctrl = open(os.path.join(d, "ctrl.txt")).read().splitlines()
        open(os.path.join(d, "ctrl.txt"), "w").write("\n".join(l for l in ctrl if not l.startswith("P 18")) + "\n")
        meta = {"np": 1, "steps": 200, "solver": "sflow", "files": ["blocks.dat"],
                "tags": ["2d", "sflow", "curvilinear", "multiblock", "river", "turbulence"],
                "description": "The 90 deg open-channel junction of Weber et al. (2001) on three blocks (W/dx 8): the "
                               "main channel downstream of the junction block 0 with the native outflow, the main "
                               "channel upstream and the junction block 1 with its inflow (G 54 1 1), the branch block 2 "
                               "with its inflow (G 54 2 1) on a part of block 1's +j face; a vertex of three blocks at "
                               "the downstream corner of the branch; the parabolic eddy viscosity (A 260 4).",
                "covers": ["the partial interface of block 1 and the vertex of three blocks with the inflow of two "
                           "blocks", "sflow_ediff::curv_stress across the block interfaces"]}
    else:
        meta = {"base": "sflow_2d_mb_weber_junction", "np": 2,
                "tags": ["2d", "sflow", "curvilinear", "multiblock", "river", "turbulence", "mpi"],
                "description": "As sflow_2d_mb_weber_junction on 2 ranks."}
    json.dump(meta, open(os.path.join(d, "case.json"), "w"), indent=1)
# P6b: NHFLOW turbulence on the blocks and the patches
derived("nhflow_3d_mb_quadrants_kw", "nhflow_3d_mb_quadrants", {
    "ctrl_add": ["A 560 2"],
    "tags": ["3d", "nhflow", "curvilinear", "multiblock", "turbulence", "quick"],
    "description": "As nhflow_3d_mb_quadrants (four turned blocks around a vertex, dam break) with k-omega (A 560 2) on "
                   "every block: the cells around the blocks take k, omega and the eddy viscosity of their neighbours, "
                   "the ring of the native grid those of the blocks; the implicit diffusion of the turbulence model (A 512 "
                   "0 -> 2) as one composite solve.  Four blocks = one block to the solver tolerance "
                   "(nhflow_multiblock/t2 with XCTRL='A 560 2').",
    "covers": ["nhflow_amr_turb.cpp: turb_objects, turb_fill0, turb_run, turb_save, turb_blocks_ini",
               "nhflow_amr::patch_objects / ini: A 512 0 with a turbulence model as A 512 2 (cdiff)"]})
derived("nhflow_3d_mb_amr_ke_mpi2", "nhflow_3d_mb_amr", {
    "np": 2, "ctrl_add": ["A 560 1"],
    "tags": ["3d", "nhflow", "curvilinear", "multiblock", "amr", "turbulence", "mpi"],
    "description": "As nhflow_3d_mb_amr (a refinement box over the vertex of four turned blocks) with k-epsilon (A 560 1) on "
                   "every block and patch, 2 ranks: the cells around the patches linear from the parent (minmod), the "
                   "covered cells the mean of their children.",
    "covers": ["nhflow_amr_turb.cpp: turb_fill (levels 1..), turb_restrict, the ring in the partition halos"]})
derived("nhflow_3d_amr_adaptive_kw", "nhflow_3d_amr_adaptive", {
    "ctrl_add": ["A 560 2"],
    "tags": ["3d", "nhflow", "amr", "turbulence", "quick"],
    "description": "As nhflow_3d_amr_adaptive (a regrid every step) with k-omega on every patch: the fresh patches take k, "
                   "omega and the eddy viscosity from an old patch of the same level or from their parent.",
    "covers": ["nhflow_amr_turb.cpp: turb_regrid, turb_free (patch_delete)"]})

# the Weber junction with NHFLOW and k-omega (W/dx 6, 4 layers): the vertex of block 1's +j face where the wall ends
# and the interface with the branch begins (the 2D wall lists of the own cells, reefamr::build_bc2D)
d = os.path.join(C, "nhflow_3d_mb_weber_junction_kw")
fresh(d)
subprocess.run([sys.executable, os.path.join(T, "weber_case.py"), d, "--solver", "nhflow", "--layout", "b", "--n", "6",
                "--nz", "4", "--np", "1", "--h0", "0.313", "--tend", "1000", "--xctrl", "A 560 2"], check=True,
               stdout=subprocess.DEVNULL)
os.remove(os.path.join(d, "blocks.json"))
ctrl = open(os.path.join(d, "ctrl.txt")).read().splitlines()
open(os.path.join(d, "ctrl.txt"), "w").write("\n".join(l for l in ctrl if not l.startswith("P 10") and not l.startswith("P 30")) + "\n")
json.dump({"np": 1, "steps": 150, "solver": "nhflow", "files": ["blocks.dat"],
           "tags": ["3d", "nhflow", "curvilinear", "multiblock", "river", "turbulence"],
           "description": "The 90 deg junction of Weber et al. (2001) on three blocks with NHFLOW and k-omega (W/dx 6, 4 "
                          "layers, bed friction A 519 1, free-surface omega A 567 3): the main channel downstream block 0 (native outflow), upstream and "
                          "junction block 1 (inflow G 54 1 1), the branch block 2 (inflow G 54 2 1) on a part of block 1's "
                          "+j face, both inflows with the log-law profile (B 61 2).  The "
                          "corner where block 1's wall ends and the interface begins: before P6b the wall ghost cells "
                          "there took the surface of the branch (a jet of 4 m/s).",
           "covers": ["reefamr::build_bc2D: the 2D wall lists of the own cells on multi-block grids",
                      "nhflow_amr_turb.cpp on blocks with open faces (G 54)",
                      "ioflow_amr_block::inflow_nhflow with B 61 2 (the log-law inflow on a face of a block)",
                      "nhflow_bcmom::wall_friction A 519 1 on the blocks",
                      "nhflow_komega_func::epsfsf A 567 3 on the blocks"]},
          open(os.path.join(d, "case.json"), "w"), indent=1)
print("written")
