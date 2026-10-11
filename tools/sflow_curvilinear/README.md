# SFLOW on curvilinear grids (G 50): tests

REEF3D::SFLOW runs on a boundary-fitted curvilinear grid in the horizontal plane (phase P1 of the
multi-block curvilinear concept). The grid is a logical box of `ni x nj` cells whose nodes have
arbitrary (x, y) positions.

## Running a case

1. Grid nodes, either
   - `curvgrid.dat` in the case folder (`G 50 2`): first line `ni nj`, then `x y zb` of node
     `(i,j)` for `j = 0..nj`, `i = 0..ni` (i fastest); `zb` is the bed level at the node, or
   - the `CURV` section of `DIVEMesh_Grid/grid-geometry.dat` (`G 50 1`).
   DIVEMesh's river corridor generator (`R 1 1`) writes both, and prints the DIVEMesh keys of the
   logical box to `DIVEMesh_Log/corridor.log`.
2. DIVEMesh grid = the logical box: `B 2 ni nj 1`, `B 10 0 ni*h 0 nj*h 0 h` with `h` about the mean
   spacing, `C 11 1` (inflow at i = 0), `C 14 2` (outflow at i = ni), `C 12 3`, `C 13 3` (banks).
3. `ctrl.txt`: `G 50 2` (or `1`). Cell bed levels are the mean of the four node bed levels.

Available on curvilinear grids so far: hydrostatic SFLOW (`A 220 0`), WENO or limited
reconstruction (`A 211`), RK2/RK3, Manning friction (`A 218`), wetting and drying, walls, discharge
inflow (`B 60`) along the inflow face normal, stage outflow (`F 50`, `F 62`), `F 72` boxes in
physical coordinates, several ranks, VTP output on the physical nodes. Refused with a message:
non-hydrostatic pressure and Boussinesq (`A 220 > 0`), diffusion and turbulence, mesh refinement
(`G 1`), wave generation, patch boundary conditions, sediment, 6DOF, SEASTATE, `I 11`.
Probes and other output that use x/y positions still refer to the logical box.

## Tests (`REEF3D_BIN`, `DIVEMESH_BIN` select the binaries; runs go to `./runs`)

| script | checks |
|---|---|
| `t1_dambreak_rotation.py` | dam break in a closed basin: Cartesian kernels, curvilinear kernels on the same nodes, and on the grid rotated by 90 deg and shifted agree per cell index to round-off (1.8e-15) |
| `t1b_channel_rotation.py` | river reach with discharge inflow and stage outflow: curvilinear = Cartesian on the same nodes (3.6e-15); rotated by 30 deg at UTM coordinates with an emerging island: equal to truncation level with WENO (the nonlinear weights act on the Cartesian components), to 8e-11 with a linear reconstruction (`A 211 0`) |
| `t2_lake_at_rest.py`, `t2b_variants.py` | lake at rest on a wavy, skewed grid: velocities 2.5e-15 fully wet; with an emerged island 1e-5 at the shoreline, the same as on the Cartesian grid (wet/dry treatment of SFLOW) |
| `t3_bend.py` | 90 deg bend, Manning friction: discharge through the outlet, superelevation 0.0099 m against U^2 B/(g R) = 0.0102 m; 2 ranks equal 1 rank to round-off |
| `make_regression_cases.py` | writes the `Regression/cases/sflow_2d_curv_*` cases |
| `read_curv.py` | reads the CURV section of a geometry file |
