# writes the regression cases of the horizontal diffusion of SFLOW on curvilinear grids, patches and blocks (P6a:
# A 212 1, A 260 3/4) and of the Weber junction on three blocks: Regression/cases/sflow_2d_*
import os, sys, json, shutil, subprocess
T = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(T, "../../Regression/cases")

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
    if os.path.exists(d):
        shutil.rmtree(d)
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
        os.makedirs(d)
        meta = {"base": "sflow_2d_mb_weber_junction", "np": 2,
                "tags": ["2d", "sflow", "curvilinear", "multiblock", "river", "turbulence", "mpi"],
                "description": "As sflow_2d_mb_weber_junction on 2 ranks."}
    json.dump(meta, open(os.path.join(d, "case.json"), "w"), indent=1)
print("written")
