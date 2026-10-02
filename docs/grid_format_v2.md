# DIVEMesh grid format v2 and the geometry core

DIVEMesh writes the grid in format v2. REEF3D only reads v2 and stops with a message
for older grid files (regenerate them with the current DIVEMesh).

## Files

* `DIVEMesh_Grid/grid-%06i.dat`, one per rank: cell flags, nodes, boundary and parallel
  surfaces, column flags, geodat bed level (`G 10`) and interpolated data (`D 10`).
* `DIVEMesh_Grid/grid-geometry.dat`, read by all ranks: the solid (`S`) and topo (`T`)
  entities as triangles (indexed mesh) with keyword, entity number, parameters,
  ray mode (`S 18`) and inversion (`S 9`, `T 9`).

Every file starts with a magic word, the format version and a byte-order marker,
followed by tagged sections (`tag[4] | int64 size | data`). Unknown sections are skipped.
The layout is documented in `src/gridfile_v2.h`.

Compared with v1:

* no `solid_dist`/`topo_dist` fields (16 bytes per cell) and no bed levels of the solids;
* cell and column flags are run-length coded;
* surface and parallel-surface lists are delta/varint coded;
* no ghost cell estimates and no unused header entries;
* the column data are written for every subdomain, also above the bottom layer (z-decomposition);
* the boundary surface count includes the plate surfaces (`S 201`).

Typical file sizes are 4 to 15 times smaller than v1.

## Geometry core

The ray-cast kernels of the 6DOF floating bodies were moved into a shared core:

| file | content |
|---|---|
| `geo_raycast_cart.cpp` | Cartesian inside/outside (ray parity) and grid-line distances (CFD, PTF, grid solids) |
| `geo_raycast_sigma.cpp` | vertical inside/outside and band distance on the sigma grid (NHFLOW) |
| `geo_raycast_column.cpp` | column crossings (FNPF) |
| `geo_primitive.cpp` | triangulated primitives (box, cylinders, wedges, hexahedron, sphere) |
| `geo_mesh.cpp` | triangle container, reader of `grid-geometry.dat` |
| `lexer_grid_solids.cpp` | solids and topography of the grid file for all modules |

Users of the core: the 6DOF objects (CFD, NHFLOW, FNPF), the NHFLOW solid forcing
(`A 581`-`A 590`) and the grid solids.

## Solids of the grid file in the modules

| module | `S` solids | `T` topography |
|---|---|---|
| CFD, PTF | signed distance field `solid` | signed distance field `topo` |
| NHFLOW | bed level (default) or immersed solids with `A 580 1` | bed level |
| SFLOW | `solidbed`, bed level | `topobed`, bed level |
| FNPF | bed level | bed level |

Columns fully blocked by solids stay inactive (column flags of DIVEMesh).

`A 580` (NHFLOW): 0 the solids raise the bed level (as before), 1 the solids are
immersed solids (direct forcing, as `A 581`-`A 590`) and the bed level only follows the
topography.

The geodat bed level (`G 10`) is united with the entities of its role (`G 9`). Before v2
the entities replaced it in the CFD distance field.

## Compatibility

The grid solids are built with the DIVEMesh ray conventions (shifted rays, faces on the
domain boundary as before), so the fields agree with v1. The exceptions are crossings that
fall exactly on a node, edge or vertex. Two v1 artefacts are not reproduced:

* in 2D, v1 measured the distance inside a solid to the domain wall at `ymax`;
* with an inverted STL (`S 9 2`), v1 clipped all distances to ±10 dx.
