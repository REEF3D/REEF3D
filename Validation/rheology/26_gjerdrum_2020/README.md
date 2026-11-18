# 26 Gjerdrum quick-clay slide, 30.12.2020 (SFLOW, Herschel-Bulkley W 120 5)

Real terrain from the Kartverket laser projects on hoydedata.no (EPSG:25832, NN2000):

- pre-event: project 4241 "Gjerdrum Ullensaker Nannestad 5pkt 2020", flown 19.04.2020
- post-event: project 5001 "Gjerdrum 5pkt 2021", flown 08.01.2021 (covers the slide and the runout)

Both DTMs averaged to 5 m (bilinear from the ImageServer, then block mean onto the cells). The release is the
scar: inside the box 613150-613650 E, 6659700-6660500 N where post < pre - 1 m. The bed is the pre-event surface,
in the release area the post-event surface (the sliding surface with the material that stayed on it); the
initial level (F 73 1) is the pre-event surface in the release area. Mobilised volume on the 5 m grid
0.94 Mm³ (NGU 2021/023: 1.35 Mm³ mobilised, 0.47 Mm³ left over the sliding surface, 0.926 Mm³ deposited in
the outlet area). Domain 612300-615200 E, 6658700-6660700 N (local x - 612300, y - 6658700), 580 x 400 cells.

## Data and set-up

The DEMs are not part of the repository (5 m GeoTIFFs, 3 MB each). Download and set-up:

    python3 ../tools/fetch_hoydedata.py --list --bbox 612300 6658700 615200 6660700
    python3 ../tools/fetch_hoydedata.py 4241 gjerdrum_pre_2020_DTM_5m.tif  --bbox 611000 6658000 615500 6662500 --res 5
    python3 ../tools/fetch_hoydedata.py 5001 gjerdrum_post_2021_DTM_5m.tif --bbox 611000 6658000 615500 6662500 --res 5
    python3 ../tools/slide_terrain.py gjerdrum_pre_2020_DTM_5m.tif gjerdrum_post_2021_DTM_5m.tif . \
            --bbox 612300 6658700 615200 6660700 --dx 5 --release 613150 6659700 613650 6660500 --thr 1.0

This writes geo.dat (DiveMESH G 10 1), sflow-level.dat (SFLOW F 73 1) and reference.npz. The pre-event DTM has
126 cells without data in the south-east corner (filled from the post-event DTM, far from the flow).

## Runs

Instant release (the event was retrogressive over several minutes: the instant release is an upper bound of the
energy), Herschel-Bulkley with tau_y about the remoulded shear strength s_ur (0.5 kPa), K 15 Pa s^n, n 0.35,
rho 1900 kg/m³; series tau_y 0.25 / 0.5 / 1 / 2 kPa (W 125). 600 s.

    ../tools/run_case.sh <REEF3D> <DiveMESH> . <run>/tau500 "mpirun -np 2"
    python3 ../tools/rheo_eval.py terrainslide <run>/tau500 . results tau500

regions.json: boxes (EPSG:25832) for the volumes of the branches (Tistilbekken next to the scar, the valley to the east with the dam, the branch to
the south).

## Results (instant release, 600 s, 2 ranks, 31 - 49 min each)

F 73 placed 0.934 Mm³ in 3665 cells (0.938 Mm³ on the DEM grid; the DiveMESH interpolation of geo.dat).
Volume conserved to 4e-5. DEM deposit (post - pre > 0.5 m outside the release, patches of at least 1000 m²):
1.05 Mm³, 0.27 km², runout 1590 m from the release centre.

| tau_y (Pa) | runout (m) | IoU | hit rate | Tistilbekken (Mm³) | east valley (Mm³) | south (Mm³) | left in release (Mm³) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| DEM  | 1590 |  |  | 0.39 | 0.61 | 0.045 | (in the bed) |
| 250  | 2030 (east boundary) | 0.55 | 0.74 | 0.19 | 0.64 | 0.058 | 0.027 |
| 500  | 1960 | 0.56 | 0.77 | 0.26 | 0.57 | 0.041 | 0.048 |
| 1000 | 1500 | 0.50 | 0.70 | 0.41 | 0.40 | 0.022 | 0.088 |
| 2000 |  820 | 0.33 | 0.48 | 0.60 | 0.15 | 0.006 | 0.168 |

- The flow follows the ravines of the DEM: Tistilbekken, the junction, both channels around the ridge to the east
  and the start of the southern branch, at all tau_y; tau_y 0.5 - 1 kPa (s_ur about 0.5 kPa) bracket the runout.
- No single tau_y gives both branches: 1 kPa fills Tistilbekken (0.41 against 0.39 Mm³) and stops short in the east
  valley (0.40 against 0.61 Mm³), 0.5 kPa fills the east valley and empties Tistilbekken. With the instant release the
  whole mass leaves at once; in the retrogressive event the later blocks came down onto the deposit in the ravine.
  A staged release (later blocks after minutes) is the next test.
- False alarms (0.33 at 0.5 kPa): a thin sheet over the field west of the scar (the instant release spills over the
  rim of the scar) and the run beyond the end of the east deposit; with 0.25 kPa the material reaches the east boundary
  (x 615200) - extend the domain for weaker material.
- The southern branch: the DEM deposit runs 800 m further south along Tistilbekken (partly piped, the culverts are not
  in the bed); the model stops at the bend at 6659200 N.
- The remnant of 0.47 Mm³ over the sliding surface (NGU) is part of the bed (post-event surface in the release area);
  the model leaves 0.03 - 0.17 Mm³ more there.

![](results/26_gjerdrum_dem.png)
![](results/26_series_maps.png)
![](results/26_series_measures.png)
![](results/26_tau500.png)
