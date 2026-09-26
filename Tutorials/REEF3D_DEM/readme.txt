REEF3D::DEM example cases (see docs/DEM.md)

1_CFD Settling Sphere Unresolved         3 mm glass bead, terminal velocity 0.351 m/s (Haider-Levenspiel)
2_CFD Settling Sphere Resolved           25 mm sphere, density ratio 1.25, d/dx = 5 (coarse, qualitative)
3_CFD Particle Pile Unresolved           250 spheres, boxes, cylinders and ellipsoids raining into a tank
4_NHFLOW Floating and Submerged Particles  floating boxes (surface-piercing large, unresolved small) and a submerged resolved sphere
5_NHFLOW Floes in Waves 2D               floating blocks in regular waves
6_CFD Box on Particles Distributed       large replicated box on distributed spheres across two ranks (dry)

Run DIVEMesh first, then REEF3D with 2 MPI ranks (M 10 2).
