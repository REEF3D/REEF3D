# MAST: app/

Standalone MAST executable: runs a MAST model without a fluid grid (engineering hydrodynamics, analytical waves, wind fields).

## Planned content

* `main` of the standalone executable
* used for fast design runs, load case batches and code checks without REEF3D's fluid solvers
* not part of the REEF3D build (the build files leave `src/mast/app/` out)

## May include from

everything except `reef3d`
