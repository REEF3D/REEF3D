<!-- Architect: Hans Bihs -->
# Reference data from public model test suites (second round, 2026-10-04)

Folder mapping used in this suite (file names unchanged):

| collected as | stored in |
|---|---|
| berkhoff_shoal/ | refdata/waves/berkhoff_shoal_funwave/ |
| conical_island/ | refdata/waves/conical_island/ |
| mase_kirby/ | refdata/waves/mase_kirby/ |
| synolakis_beach/ | refdata/waves/synolakis_runup_funwave/ |
| spheric_test2_kleefsman/ | refdata/dambreak/spheric_test2/ |
| porous_dambreak_lin1998/ | refdata/porous/lin1998/ |

The other data sets listed below (OSU shelf, composite beach, Monai valley, Ghia cavity) are not
used by a case yet and are not stored in the repository; they can be re-extracted from the listed
commits.


Collected 2026-10-04. Each data file has a `#` header that gives the reference, the source repo with commit and path, the column definitions with the script lines that define them, the sampling, and any caveats. Values are copied verbatim. The only changes are: header/comment lines dropped, delimiters normalised to single spaces, and `.mat` files converted to text by printing the stored doubles with Python `repr`.

Confidence levels: **A** = original instrument records (regular sampling). **B** = digitised, tabulated or derived by a third party. **C** = origin or meaning unclear.

## Source repositories

| Short | URL | Commit |
|---|---|---|
| FUNWAVE | https://github.com/fengyanshi/FUNWAVE-TVD | b4c322e7582035ee19df8e6409a3dfedaff1cb96 |
| NHWAVE | https://github.com/fengyanshi/NHWAVE | 286036b49872f7783de8cbcf574b018c314f3ed8 |
| AQUAgpusph | https://github.com/sanguinariojoe/aquagpusph | 9ca1e1b5e168b6e55ef0a171d481e769cfe7b4f8 |
| olaFlow | https://github.com/phicau/olaFlow | 9879096092a4fc13419a40bbe11c4bd1b87ddd1d |
| Basilisk mirror | https://github.com/comphy-lab/basilisk-C | 67ddf8a1936ea3a6662e0c6964bb18c1c1d6bbaa |
| cuIBM | https://github.com/barbagroup/cuIBM | 0b63f86c58e93e62f9dc720c08510cc88b10dd04 |
| sujaldave | https://github.com/sujaldave/2DLidDrivenCavityBenchmark | 53c82caee5af5ec9a8898d426f9449e6e582d869 |

## Files

| File | Content | Source (repo: path) | Conf. |
|---|---|---|---|
| berkhoff_shoal/section1..5.dat | Berkhoff, Booy & Radder (1982), elliptic shoal. Transverse sections at x = 1, 3, 5, 7, 9 m behind the shoal centre. Columns: y [m], amplitude a [mm]. H = 2a. H/H0 = a/23.2 | FUNWAVE: benchmarks/car_berkhoff_2d/postprocessing/sectionN.dat | B/C (regular 0.25/0.5 m spacing; values look like 23.2 mm x a 2-decimal ratio; tabulated or digitised, not documented) |
| berkhoff_shoal/section678.dat | Sections 6/7/8 along y = -2, 0, +2 m. Columns: -x [m], then a [mm] for sections 6, 7, 8 | FUNWAVE: same folder, section678.dat | B/C |
| berkhoff_shoal/FUNWAVE_DEPTH.F | FUNWAVE bathymetry generator (verbatim source code; includes the Qin Chen shoal formula) | FUNWAVE: car_berkhoff_2d/input/DEPTH.F | n/a |
| conical_island/briggs_case{A,B,C}_gauges_6_9_16_22.dat | Briggs et al. (1995) conical island, solitary wave with H/d = 0.045, 0.091 (0.096), 0.181. Columns: t [s], eta [m] at gauges 6, 9, 16, 22. dt = 0.04 s, t = 20 to 80 s | FUNWAVE: car_conical_island/work_case_X/X.txt (identical to the NTHMP Conical_GaugeABC.mat in NHWAVE) | A |
| conical_island/briggs_case{A,B,C}_runup.dat | Maximum run-up around the island. Columns: direction [deg], R [cm], R normalised (about R/32 cm). 24 points per case. The angle convention is NOT stated; the data suggest 270 deg = side facing the incident wave and 90 deg = lee side | NHWAVE: examples/inundation_benchmarks/Solitary_wave_on_a_conical_island/runup/conical_RU.mat | A/B |
| mase_kirby/mase_kirby_eta_hNNNmm.dat (12 files) | Mase & Kirby (1992), irregular waves on a 1:20 slope, offshore depth 47 cm. Measured eta [cm] at the gauge with still-water depth NNN mm. One column, dt = 0.05 s, 15000 samples | FUNWAVE: car_mase_kirby/postprocessing/r2dNNN.dat | A |
| mase_kirby/mase_kirby_skew_asym_data.dat | Skewness (row 1) and asymmetry (row 2) computed from the data, columns h = 47 ... 2.5 cm. The 13th column is never used and its meaning is unknown | FUNWAVE: car_mase_kirby/postprocessing/mkskew_man.out | B |
| osu_shelf_runup/osu_WG1..9.dat | Solitary wave (a = 0.39 m, h = 0.78 m) over a shelf with an island, OSU basin. Probably the 2009 workshop BM2 (Swigler 2009); the reference is NOT in the repo. Columns: t [s], eta [cm]. dt = 0.02 s | FUNWAVE: car_osu_runup/postprocess/WGk.txt | A (data). Reference: C |
| osu_shelf_runup/osu_ADV{A,B,C}_{U,V}.dat | ADV velocities at the same experiment. Columns: t [s], U or V [m/s]. Contains NaN. Not offset-corrected | FUNWAVE: car_osu_runup/postprocess/{U,V}_Velocity_AverageX.txt | A (data). Reference: C |
| composite_beach/composite_beach_case{A,B,C}_g4_g10.dat | Composite beach (Revere Beach), NTHMP BP5. Probable reference Briggs et al. (1996), not stated in the repo. Columns: t [s], eta [m] at gauges 4 to 10. dt = 0.05 s. Values are 0.001 ft multiples | FUNWAVE: sph_comp_beach/case_X/X.mat | A (data). Lab gauge positions: C |
| monai_valley/monai_gauges_5_7_9.dat | Monai valley 1:400 model (Matsuyama & Tanaka 2001), NTHMP BP7. Columns: t [s], eta [cm] at gauges 5, 7, 9. dt = 0.05 s | FUNWAVE: sph_monai_valley/postprocessing/A.mat | A (positions not established) |
| monai_valley/monai_incident_wave.dat | Monai incident boundary wave. Columns: t [s], eta [m] | FUNWAVE: sph_monai_valley/fft/MonaiValleyWave.txt | A/B |
| synolakis_beach/synolakis_H{0.0185,0.30}_tT.dat | Synolakis (1987), solitary wave on a 1:19.85 beach. Profiles x/d vs eta/d at t* = T | FUNWAVE: sph_sol_plane_mea/postprocessing/SW_01_*.mat | B |
| spheric_test2_kleefsman/spheric_test2_P1-P8_H1-H4.dat | Kleefsman et al. (2005), MARIN dam break with a box (SPHERIC Test 2). Columns: t [s], P1-P8 [Pa], then 4 water heights [m]. **Warning:** AQUAgpusph reads height columns 10-12 as H3, H2, H1, and arrival times confirm column 10 is the gauge nearest the reservoir | AQUAgpusph: examples/3D/spheric_testcase2_dambreak/src/templates/test_case_2_exp_data.dat | A (height naming: check paper) |
| porous_dambreak_lin1998/lin1998_CR35_tT.dat | Lin (1998) / Liu et al. (1999) dam break through a crushed-rock porous dam (D50 = 15.9 mm, n = 0.49). Free-surface profiles x [m], y [m] at t = 0, 0.35, 0.75, 1.15, 1.55, 1.95 s | olaFlow: tutorials/CR35_dambreak/labData/b1..6.dat | B (digitised) |
| beji_battjes_bar/bar_gaugeK.dat (K = 4..11) | Beji & Battjes (1993) / Luth et al. (1994) submerged bar, T = 2.02 s. Columns: t [s] (arbitrary origin), eta [cm]. Each file holds the digitised period twice | Basilisk mirror: basilisk-source/src/test/gauge-K | B (digitised) |
| ghia_cavity/ghia1982_centerlines_cuIBM.dat | Ghia et al. (1982), centre-line u and v for Re = 100, 1000, 3200, 5000, 10000 (numerical benchmark). Flag: Re = 3200 value at y = 0.4531 is -0.86636, an outlier | cuIBM: data/ghia_et_al_1982_lid_driven_cavity.dat | B |
| ghia_cavity/ghia_table_I_u.dat, ghia_table_II_v.dat, ghia_vortex_centers.dat | Ghia Tables I, II and V for Re = 100, 400, 1000, 5000 (independent transcription). One value differs from cuIBM: v at x = 0.5, Re = 1000 | sujaldave: literature/*.csv | B |

## Not extracted, with reasons

- **FUNWAVE car_mase_kirby/postprocessing/mketa.out** (2001 x 11, dt 0.01 s, m). comp_series.m plots it as if measured. However, it correlates with the measured r2d records only at r = 0.95 offshore, falling to 0.56 at h = 2.5 cm, and it does not match at any lag. It is most likely model output. Confidence C.
- **FUNWAVE car_mase_kirby/postprocessing/mkskew_kirby.out**: Kirby et al. model results.
- **FUNWAVE car_osu_runup/postprocess/tank.dat** (125 x 299): looks like a measured bathymetry survey, but its grid and coordinates are undocumented.
- **FUNWAVE car_sediment_lab_c2, car/sph_Nwave_statis, sph_sol_plane_ana, sph_nesting**: no measured data (model inputs, analytical solutions or statistics only).
- **NHWAVE current_benchmarks/Steady_flow_over_submerged_obstacle/vel_comp/SL_S*.DAT**: velocity time series for flow over a submerged conical mound (likely Lloyd & Stansby 1997). They are digitised, with non-monotonic time, and were not requested.
- **NHWAVE examples tingkirby, submerged_bar, berkhoff**: model inputs only, no measured data.
- **AQUAgpusph**: also contains SPHERIC Test 5 (examples/2D/spheric_testcase5_dambreak/src/doc/d38_*.csv) and Lobovsky et al. (2014) data. Not extracted (not requested).
- **Not found in any accessible public GitHub repo**: Ting & Kirby (1994/1995) wave height and setup; Whalin (1971) harmonic amplitudes; original (non-digitised) Beji-Battjes/Luth gauge records. OceanWave3D ships Whalin and Berkhoff input files only. gist.github.com and basilisk.fr are blocked by the proxy.
