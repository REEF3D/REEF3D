/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs
 *
This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include "control.h"
#include<math.h>

void control::ini_default()
{
    // Hydrodynamic Models
    A10=6;       // int hydrodynamic models: 2 SFLOW, 3 FNPF, 5 NHFLOW, 6 CFD, 7 Spectral
    
    // SFLOW
    A209=1;      // int interpolation sweeps for bed
    A210=3;		  // int time scheme for SFLOW velocities
    A211=4;		  // int convection scheme for SLOW velocities
    A212=-1;	  // int diffusion treatment for SFLOW velocities (-1: not set, becomes 1 with a turbulence model A 260 > 0, else 0)
    A213=1;      // int SFLOW WENO5 nonlinear weights (A 211 4): 0 WENO-JS, 1 WENO-Z (default), 2 TENO5
    A214=1;      // int convection for vertical velocity
    A225=1;      // int SFLOW Riemann solver: 1 HLL, 2 HLLC (contact-resolving tangential fluxes)
    A217=2;      // int slip or no-slip boundary conditions
    A218=0;      // int turn on roughness
    A219=1;      // int additional courant number constraint
    A220=2;		  // int non-hydrostatic pressure scheme for SFLOW
    A221=1;		  // int non-hydrostatic pressure in very shallow regions
    A222=1.0e-5; // double SFLOW TENO5 cutoff C_T (A 213 2)
    A223=0.5;    // double blending factor hydrostatic pressure gradient
    A224=1.159;  // double dispersion parameter alpha for A 220 3 (1: quadratic, 1.159: optimised)
    A240=1;      // int FSF algorithm SFLOW
    A242=0;          // int non-hydostatic pressure for shallow areas
    A243=2;      // int turn on wetting-drying
    A244=0.001; // double absolute wetting criterion value
    A246=1;      // int turn on breaking
    A247=0.6;    // double breaking parameter alpha
    A248=0;      // int turn on breaking persistence
    A249=0.3;    // double breaking persistence parameter beta
    A250=1.86;   // double viscosity breaking wave
    A226=0;      // int SFLOW free surface: 0 explicit (HLL), 1 semi-implicit theta method
    A227=0.55;   // double SFLOW semi-implicit: implicitness theta of the free surface
    A228=0.0;    // double SFLOW semi-implicit: maximum time step (s), 0 no cap
    A229=1;      // int SFLOW semi-implicit: wetting and drying iteration (Casulli 2009), 0 off
    A232=0.05;   // double SFLOW semi-implicit: depth below which the cell velocities of the momentum predictor are desingularised
    A231=0.0;    // double SFLOW semi-implicit: mass flux dissipation for transcritical flow, 0 off
    A251=0;    // double fsf-slope in x-dir
    A260=0;      // int turbulence model
    A261=0.267;  // double length scale factor
    A262=0.0667; // double parabolic turbulence model factor
    A263=10.0;   // double eddyv limiter factor set to high
    A264=3.6;    // double epsilon coefficient ce_gamma (Rastogi & Rodi 1978: nu_t = 0.077 u* h)


    // FNPF
    A309=1;      // int interpolation sweeps for bed
    A310=3;		  // int time scheme for FNPF velocities
    A311=5;		  // int convection scheme for FNPF velocities
    A312=2;      // int discretization for second-order gradient
    A313=3;      // int discretization for bed bc
    A314=2;      // int linear/non-linear
    A315=1;      // int eta gradient treatment: 0 legacy, 1 kinematic FSBC upwinded by dH/deta_x and symmetric at zero speed, 2 as 1 plus symmetric Ex,Ey for geometry (sigma metrics, breaking, wind, dynamic FSBC), 3 as 1 with the Godunov gradient (sonic points) in the kinematic FSBC
    A316=1;      // int wet-dry WENO5 near the shoreline: 0 legacy (zero gradient), 1 first-order fallback from wet neighbours
    A317=1;      // int bed metrics: 0 legacy (Bxx=Byy=0, left-biased Bx without wetting-drying), 1 bed curvature on, symmetric Bx
    A318=0;      // int WENO5 nonlinear weights for FNPF free-surface gradients: 0 WENO-JS, 1 WENO-Z, 2 TENO5
    A319=1.0e-5; // double TENO5 cutoff C_T for FNPF (A318 2)
    A320=1;		  // int order of Laplace equation
    A324=0;      // int FNPF Laplace: max Picard iterations for the explicit sigma cross-derivative terms (0 = lagged by one RK stage, legacy)
    A325=1.0e-3; // double FNPF Laplace: Picard stop when the change is below A325 times the change of the first solve (A324>0)
    A328=0;      // int FNPF Laplace: 1 sigma cross-derivatives implicit, 15-point (2D: 9-point) operator in REEFMG's BiCGStab, 7-point V-cycle preconditioner (N 10 1 only); 0 explicit
    A326=0.0;    // double static coastline (A343 1): coastline relaxation towards the running mean of eta, Fi over this time scale (s) instead of 0; 0 = relax to 0
    A327=1;      // int FNPF free-surface vertical velocity Fz: 0 legacy index-space one-sided stencil, 1 Fornberg weights on the sigma nodes ZN
    A329=1;      // int wave maker BC order
    
    
    A330=10;     // int dynamic wetting-drying (A343 2/3): minimum number of time steps a rewetted cell stays wet
    A331=1.0;    // double dynamic wetting-drying (A343 2/3): rewetting hysteresis margin as factor of A344
    A332=2;      // int dynamic wetting-drying (A343 2/3): width in cells of the front viscosity band (viscosity A346, needs A350 1), 0 off
    A334=1;      // int dynamic wetting-drying (A343 2/3): redistribute wet-dry clamp volume to wet neighbours: 1 on, 0 off
    A335=1.0e20; // double runup (A343 2): maximum bed elevation above still water level that can become wet
    A336=1;      // int dynamic wetting-drying (A343 2/3): eta gradients at the wet-dry front from wet cells only: 1 on, 0 off
    A337=1;      // int static coastline (A343 1): 0 legacy; 1 wet-only Fifsf clamp average, clamp volume redistribution (A334), wet-only eta gradients and Exx=0 at the coast (A336), viscosity A346 ramped inside the coastline band and in thin-film cells (WL < 10*A344); 2 as 1 but A346 only in thin-film cells (thin coastline)
    A338=1;      // int remove wet areas not connected to the main water body (coastline ini): 1 on, 0 off
    A339=1;      // int coastline initialization: 1 = fast sweeping, 0 = PDE reinitialization
    A341=0.0;    // double coastline damping distance factor for dxm
    A342=0.0;    // double coastline damping absolute distance
    A343=1;      // int wetting-drying: 0 off, 1 static coastline, 2 dynamic runup and rundown (no coastline damping), 3 dynamic rundown only
    A344=0.001; // double wetting criterion value
    A346=1.86;    // double viscosity damping within the coastline
    A347=1;     // int coastline relaxation for Fi and eta
    A348=1;     // int beach relaxation for Fi and eta
    A349=10.0;   // double initial coastline damping distance factor for dxm
    
    A350=0;      // int turn on breaking (which method)
    A351=3;      // int type of breaking detection (deep / shallow)
    A352=0;      // int additional filtering to viscosity based breaking
    A354=0.6;    // double breaking parameter alpha
    A355=1.25;   // double breaking parameter slope alpha
    A356=0.1;   // double breaking parameter slope beta
    A357=1;     // int breaking for Fi and eta
    
    A361=5;      // int breaking filter outer iter
    A362=2;      // int breaking filter inner iter
    A365=1.86;   // double viscosity breaking wave
    A366=1.0;   // double viscosity breaking wave alpha_eta
    A359=0.85;   // double breaking A350 4: onset threshold of B = u_s/C at the crest
    A360=8.0;    // double breaking A350 4: termination when the largest front-face slope is below this angle (deg)
    A364=1.0;    // double breaking A350 4: coefficient of the bore (hydraulic jump) dissipation g C H^3/(4 h_c)
    A367=1;      // int breaking A350 4: deep-water breaking strength b(Gamma) of Derakhti et al. (2018): 1 on, 0 off
    
    A370=0;      // int wind model
    A371_u=0.0;  // double wind velocity
    A371_dir=0.0;  // double wind direction
    A372=0;      // int wind forcing area
    A372_xs=A372_xe=A372_ys=A372_ye=0.0; // double wind forcing area
    A373=1;      // int wind forcing region
    A374=0;      // int wind decay
    A375=32.0;   // double Miles/Plant wind input coefficient beta, p_a = beta*rho_a*u*^2*deta/dn
    A376_s=0.5;  // double Jeffreys sheltering coefficient
    A376_sc=0.3; // double Jeffreys critical slope
    A376_c=0.0;  // double Jeffreys wave phase speed, <=0: use wave generation celerity
    A377=0.0;    // double reference depth for wind setup, <=0: mean still water depth
    A378=4;      // int low-pass filter passes on the wind forcing slope

    A380=0;      // int FNPF ice floes, 1: read ice_floes.dat
    A381=20.0;   // double ice lid stiffness k/(rho_w*g)
    A382=0.2;    // double ice lid damping ratio (of critical, per floe)
    A383=0.0055; // double ice-water drag coefficient
    A384_mu=0.3; // double ice contact friction coefficient
    A384_e=0.0;  // double ice contact restitution coefficient
    A385=1;      // int ice floes: domain boundaries act as walls
    A386=4;      // int ice footprint supersampling per cell direction
    A387=-1.0;   // double ice floe print interval in s, <0: follow fsf print settings
    A388=30;     // int ice contact solver iterations
    A389=2.0;    // double ice footprint edge taper half width in cells
    A390=0;      // int ice breaking, bits: 1 flexural, 2 contact splitting, 4 spalling (7 all)
    A391=0.5e6;  // double ice flexural strength [Pa]
    A392_K=1.0e5;// double ice fracture toughness K_IC [Pa m^0.5] for splitting
    A392_C=1.0;  // double splitting load coefficient, F = C*K_IC*h*sqrt(D)
    A393=-1.0;   // double minimum floe size after breaking [m], <0: 4 mean cells
    A394=1;      // int ice breaking check every n steps
    A395_dir=12; // int flexural cut directions
    A395_off=64; // int flexural beam bins (cut positions) per direction
    A396=0.05;   // double contact force averaging time for splitting [s]
    A397_E=5.0e9;// double ice Young's modulus [Pa] for the flexural stress, <=0: rigid-floe statics
    A397_nu=0.3; // double ice Poisson ratio
    A398=0.0;    // double ice crushing strength [Pa], contact force cap sigma_c*h*w, <=0: off
    A399_m=0.0;  // double Weibull modulus of the flexural strength, <=0: no scatter, no size effect
    A399_A=1.0;  // double Weibull reference area [m^2]
    A399_seed=1; // int Weibull random seed
    A400=1;      // int ice: max flexural cracks per floe and check
    A401=1;      // int ice: radial cracks per contact splitting event
    A402_L=-1.0; // double ice spall length [m], <0: floe thickness
    A402_C=3.0;  // double ice spall locality: chip chord <= C*chip depth
    A402_T=1.0;  // double ice rubble lid fade-out time [s]
    A403=1.0;    // double ice lid time step factor, dt <= A403/omega_lid; <0: old limit dt*1/sqrt(1+A381)

    
    // PTF
    
    
    // NHFLOW
    A509=1;      // int interpolation sweeps for bed
    A510=2;      // int NFHLOW time scheme
    A511=1;		// int NHFLOW HLL scheme
    A512=0;		// int NHFLOW diffusion
    A513=1;		// int NHFLOW wall diffusion
    A514=4;		// int NHFLOW reconstruction 
    A515=1;      // int Dirichlet wave BC type: 1 still-water flux, 2 actual-depth flux, 3 theoretical eta, 4 mass-conservative (still-water flux, U=UH/D, exact face flux)
    A516=1;      // int forcing WL/eta treatment
    A517=1;      // int preserve tangential velocities for forcing
    A518=2;      // int NHFLOW bed BC
    A519=0;      // int turn on bed roughness
    A520=2;		// int NFHLOW non-hydrostatic pressure scheme
    A521=0;		// int wetdry/forcing flux 
    A522=5.0;    // double p_alpha
    A524=2;      // int sigma gradients
    A526=1.1;    // double direct forcing factor
    A527=0;      // int NHFLOW WENO5 nonlinear weights (A514 4,5): 0 WENO-JS, 1 WENO-Z (default), 2 TENO5
    A528=1.0e-5; // double NHFLOW TENO5 cutoff C_T (A527 2)
    A529=0;      // int NHFLOW Poisson sigma cross terms: 0 explicit (P of the last solve), 1 implicit in the matrix, 15-point (2D: 9-point) operator in REEFMG's BiCGStab (N 10 1 only)
    A530=0;      // int NHFLOW Poisson sigma cross terms in porous cells (POR < 1): 0 off, 1 on
    A531=3.0;    // double Froude number limiter
    A532=1;      // int  Froude number limiter area
    A533=0;      // int  add veritcal velocity to CFL 
    A534=0;      // int NHFLOW hydrostatic semi-implicit free surface (theta method), 1 on
    A535=0.55;   // double NHFLOW semi-implicit: theta
    A536=0.0;    // double NHFLOW semi-implicit: maximum time step (s), 0: explicit gravity wave limit
    A537=1;      // int NHFLOW semi-implicit: wetting and drying iteration (Casulli 2009), 0 off
    A538=0.0;    // double NHFLOW semi-implicit: mass flux dissipation for transcritical flow, 0 off
    A539=0;      // int NHFLOW semi-implicit, curvilinear grid: non-orthogonal part of the face gradient of eta (explicit, deferred), 1 on
    A540=1;      // int NFHLOW wetdry scheme
    A544=0.001;  // double wetting & drying criterion
    A545=10.0;   // double deep criterion
    A546=0.05;   // double NHFLOW semi-implicit: depth below which the cell velocities of the momentum predictor are desingularised
    
    A550=0;      // int turn on breaking (which method)
    A551=3;      // int type of breaking detection (deep / shallow)
    A552=0;      // int additional filtering to viscosity based breaking

    A554=0.6;    // double breaking parameter alpha
    A555=1.25;   // double breaking parameter slope alpha
    A557=1.86;   // double viscosity breaking wave
    
    A560=0;      // int turbulence model
    A564=1;      // int eddyv limiter (NHFLOW: 1 realizability, 2 SST; >= 1 also the length-scale limit l <= kappa h)
    A565=0;      // int stabilization
    A566=0;      // int buoyancy
    A567=0;      // int fsf eps
    A568=1.0;      // double URANS C3 coefficient
    A561=0;        // int horizontal/vertical eddy viscosity split (HLES): 0 off, 1 nu_h = nu_t + nu_sgs (vertical-shear production), 2 nu_h = nu_sgs
    A562_cs=0.15;  // double horizontal Smagorinsky constant for A 561
    A562_nub=0.0;  // double background horizontal eddy viscosity for A 561
    A563=0.0;      // double QCR coefficient C_cr1 (0: off, Spalart 2000: 0.3)
    A559=0;        // int NHFLOW EARSM (Wallin & Johansson 2000) with k-omega (A 560 2/22): 0 off, 1 on
    
    A570=0;      // int wind model
    A571_u=0.0;  // double wind velocity
    A571_dir=0.0;  // double wind direction
    A572=0;      // int wind forcing area
    A572_xs=A572_xe=A572_ys=A572_ye=0.0; // double wind forcing area
    A573=1;      // int wind forcing region
    A574=0;      // int wind decay 
    A575=32.0;   // double Miles/Plant wind input coefficient beta
    A576_s=0.5;  // double Jeffreys sheltering coefficient
    A576_sc=0.3; // double Jeffreys critical slope
    A576_c=0.0;  // double Jeffreys wave phase speed, <=0: use wave generation celerity
    A578=4;      // int low-pass filter passes on the wind forcing slope
    A580=0;      // int DIVEMesh solids (S): 0 bed level, 1 immersed solids (direct forcing)
    
    A581=0;      // int solid box
    A583=0;      // int solid cylinder y
    A584=0;      // int solid cylinder z
    A585=0;      // int solid cylinder jacket member
    A586=0;      // int solid sphere
    A587=0;      // int solid wedge x
    A588=0;      // int solid wedge y
    A589=0;      // int solid wedge z
    
    A590=0;      // int STL
    A591=0;      // int scaling of stl geometry on/off (set by A 591)
    A591_x=A591_y=A591_z=1.0;  // double scaling of stl geometry
    A592=0;     // int translation on/off
    A592_dx=A592_dy=A592_dz=0.0;  // double translation of stl geometry
    A593=0;
    A593_x=A593_y=A593_z=A593_phi=A593_theta=A593_psi=0.0;
    A599=0;     // use dlm instead of df

    // Spectral
    A700=1;           // int mode: 1 nonstationary, 2 stationary
    A701=32;          // int number of frequencies (logarithmic spacing)
    A702_fmin=0.04;   // double lowest frequency [Hz]
    A702_fmax=1.0;    // double highest frequency [Hz]
    A703=36;          // int number of directions (full circle)
    A704=16;          // int tile size (cells) of the block-sparse action storage
    A705=0.05;        // double minimum water depth of an active (wet) cell [m]
    A706=60.0;        // double time step [s]
    A707=0;           // int iterations per time step (0: 1 nonstationary, 50 stationary)
    A708=0.001;       // double stationary convergence: max. relative change of Hs per iteration
    A709=100.0;       // double stationary convergence: percentage of the active cells that must meet A 708 (SWAN npnts, e.g. 99.5)
    A710=0;           // int initial spectrum: 0 zero, 1 parametric (B 85 spectrum, B 93 Hs Tp, B 88 gamma, B 130/131/134 spreading)
    A711=0;           // int boundary spectrum: 0 none, 1 parametric (as A 710 1), 2 SWAN 2D spectrum file seastate-boundary.spc, 3 time series of spectra at many locations (seastate-boundary.spc)
    A712_xm=1;        // int side x-: 0 open (no incoming waves), 1 boundary spectrum, 2 zero gradient
    A712_xp=0;        // int side x+ (as x-)
    A712_ym=0;        // int side y- (as x-)
    A712_yp=0;        // int side y+ (as x-)
    A713=1;           // int depth and current refraction
    A714=1;           // int frequency shift
    A716_w=0.0;       // double automatic fine direction sector around the parametric boundary spectrum: largest bin width there [deg] (0 off)
    A716_m=15.0;      // double ... margin on both sides of the boundary spreading [deg]
    A717=0;           // int computational direction sector (as SWAN CGRID SECTOR): 1 on, set by the key
    A717_th1=A717_th2=0.0; // double ... from this direction counter-clockwise to this one [deg], ccw from +x (propagation direction); the directions outside are not computed and hold no energy; on uniform directions only those of the sector are stored (Phase 12)
    A720=0;           // int prescribed current (stand-alone runs): 0 none, 1 U linear in x (A 721)
    A721_us=A721_ue=0.0; // double U at xs and at xe [m/s]
    A721_xs=A721_xe=0.0; // double xs, xe [m]
    A729=0;           // int drag of the wind input: 0 Wu (1982), 1 the met drag law W 207 (W 213, W 214)
    A730=0;           // int wind: 0 none, 1 uniform (A 731), 2 the met forcing W 200 (uniform, met-forcing.r3m, vortex, blend)
    A731_u10=0.0;     // double wind speed U10 [m/s]
    A731_dir=0.0;     // double wind direction [deg], direction the wind blows to, ccw from +x
    A732=0;           // int deep-water physics: 0 off, 1 Komen (wind input with wind, whitecapping, quadruplets A 733), 2 van der Westhuysen (Yan wind input, saturation-based whitecapping)
    A733=1;           // int quadruplets with A 732 1: 0 off, 1 DIA
    A734=1.5e-3;      // double linear wind growth (Cavaleri and Malanotte-Rizzoli), 0 off
    A735=0.1;         // double action density limiter with A 732 1 (SWAN, Ris 1997), 0 off
    A740=0;           // int depth-induced breaking: 0 off, 1 Battjes-Janssen
    A741_alpha=1.0;   // double breaking alpha
    A741_gamma=0.73;  // double breaking gamma
    A742=0;           // int bottom friction: 0 off, 1 JONSWAP
    A743=0.038;       // double JONSWAP friction coefficient [m^2/s^3]
    A744=0;           // int triads: 0 off, 1 LTA (Eldeberky)
    A745=0.05;        // double triads proportionality coefficient alpha_EB
    A736_cutfr=2.5;   // double triads: highest sum frequency / sigma_m01 (SWAN cutfr)
    A736_urcrit=0.63; // double triads: critical Ursell number of the biphase (SWAN urcrit; SWAN 40.x-41.31: 0.2)
    A736_urslim=0.1;  // double triads: lowest Ursell number with triads (SWAN urslim; SWAN 40.x-41.31: 0.01)
    A737=1;           // int with Battjes-Janssen breaking (A 740 1): 1 the total energy of a cell is limited to (gamma d)^2/4 as in SWAN, 0 off
    A746=10.0;        // double Roelvink breaking (A 740 2) exponent n
    A747=2.0;         // double surfbeat: maximum wave height to water depth ratio H/h (cap)
    A748=1;           // int surfbeat: roller 0 off, 1 on
    A749=0.1;         // double surfbeat: roller front slope beta
    A750=0;           // int coupling with SFLOW (A 10 2): 0 off, 1 SEASTATE <-> SFLOW (wave forcing, water level and currents back)
    A751=1;           // int wave forcing of SFLOW: 1 radiation stress, 2 vortex force (SFLOW velocity Eulerian, Stokes transport in the continuity)
    A752=0.0;         // double ramp-up time of the wave forcing [s], 0: none
    A753=2;           // int feedback SFLOW -> SEASTATE: 0 none, 1 water level, 2 water level and currents
    A760=0;           // int number of handover points (A 760 x y): 2D spectrum files for the wave generation (B 85 11)
    A761=90.0;        // double half-width of the handover direction sector around the mean direction [deg]
    A770=0;           // int surfbeat (wave-group action balance, single representative frequency): 0 off, 1 on
    A771=0.0;         // double surfbeat: representative period T_rep [s], 0: Tm-1,0 of the boundary spectrum
    A772=3600.0;      // double surfbeat: record length of the boundary time series [s] (repeats afterwards)
    A773=1;           // int surfbeat: seed of the random phases and directions
    A774=1;           // int surfbeat: long waves at the offshore boundary of the host: 0 absorbing only, 1 bound long wave in, absorbing out
    A780=0.0;         // double model start date-time YYYYMMDD.HHMMSS of the forcing files (A 711 3, A 730 2), 0: first time in the files
    A790=0;           // int bathymetry raster seastate-bathy.dat (stand-alone): 0 off (bed of the 2D grid, patches interpolated), 1 bed of all grids from the raster
    A791=0.0;         // double mesh refinement (G 1): refine where the water depth is below this value [m], 0 off
    A792=0;           // int mesh refinement (G 1): refine within this many cells of a dry or land cell (1-3), 0 off
    A793=0.0;         // double mesh refinement (G 1): refine where the depth changes by more than this fraction to a neighbour cell, 0 off
    A794=4;           // int mesh refinement (G 1): no refinement within this many level-0 cells of the sides with boundary spectra or zero gradient (A 712 1, 2)
    A795=0.0;         // double spectral sparsity (stationary): bins with less than this fraction of the energy of the cell are set to zero and not solved, 0 off
    A796=1;           // int geographic advection: 1 first-order upwind, 2 second-order upwind (SORDUP type)
    A715_th1=0.0;     // double fine direction sector from this direction ... [deg], ccw from +x
    A715_th2=0.0;     // double ... counter-clockwise to this direction [deg]
    A715_k=1;         // int ... the A 703 directions with centres in the sector divided into this many bins each (1 off)
    A718=0;           // int diffraction (phase-decoupled, Holthuijsen et al. 2003): 0 off, 1 from the total energy (as SWAN), 2 per frequency
    A719=0;           // int diffraction: smoothing steps of the energy (0: 0.4 (L/dx)^2 from the mean wavelength and the smallest cell)
    A723_kr=0.0;      // double reflecting coasts (A 723 Kr pown): wave-height reflection coefficient of the faces to land (0: off)
    A723_pown=0.0;    // double reflecting coasts: 0 specular, > 0 diffuse with cos^pown around the specular direction (as SWAN RDIFF)
    A724=3;           // int reflecting coasts: the coastline from the land faces within this many cells (0: the face itself, at most 3)
    A725=0;           // int number of obstacle structures (A 725 n type a b c): 1 d'Angremond (slope [deg], crest width), 2 per frequency from seastate-obstacle-n.dat, 3 porous (width, porosity, stone diameter)
    A726=0;           // int number of diffuse obstacle reflections (A 726 n pown)
    A754=0;           // int vegetation: 0 off, 1 Dalrymple / Suzuki et al. (2011) (as SWAN VEGETATION 1), 2 per frequency, Jacobsen et al. (2019) (SWAN 2)
    A755_h=0.0;       // double vegetation: height [m]
    A755_d=0.0;       // double vegetation: stem diameter [m]
    A755_n=0.0;       // double vegetation: stems per m^2
    A755_cd=1.0;      // double vegetation: drag coefficient
    A756=0;           // int vegetation: 1 stems per m^2 per cell from seastate-vegetation.dat (raster as seastate-bathy.dat)
    A722=0;           // int number of obstacles (A 722 xs ys xe ye Kt Kr zc): line obstacles, energy transmission Kt^2 (Kt < 0: Goda with crest level zc), reflection Kr^2
    A738=1;           // int stationary: source iterations per cell and sweep (max.), the source terms evaluated again from the spectrum just solved
    A739=1.0e-4;      // double ... until the estimated distance of the energy of the cell to its balance is below this fraction
    A798=1;           // int threads per MPI rank: the tiles of the sweeps, the patches of a refinement level, the loops over the cells
    A781=1;           // int MPI: 1 the sweeps pipelined across the ranks (the result of one rank; grids without refinement or with the composite sweep A 797 1), 0 each rank from the halo of the last exchange
    A799=0;           // int stationary convergence test: 0 change of Hs per iteration (A 708), 1 estimated distance of Hs to the solution
    A757_cds2=5.0e-5; // double A 732 2 (van der Westhuysen): coefficient Cds2 of the saturation-based whitecapping (SWAN cds2)
    A757_br=1.75e-3;  // double ... saturation threshold B_r (SWAN br)
    A762=0;           // int mesh refinement (G 1): refine where the water width (shortest wet run through the level-0 cell along x, y and the diagonals) is below this many cells of the level, 0 off
    A758_n=0;         // int FAS coarse-grid correction every this many iterations (nonstationary: of the step), 0 off; mesh refinement (G 1): level 0 the coarse grid, all levels corrected; uniform grid (Phase 12): a V-cycle over coarse grids of 2 x 2 cells
    A758_m=4;         // int ... coarse-grid iterations of each correction (uniform grid: m/2 before and m/2 after the next coarser grid, 4 m on the coarsest)
    A797=0;           // int mesh refinement (G 1): 0 sweeps level by level (V-cycle), 1 composite sweep (the finer cells within the sweep of their parent cell)
    A775=2;           // int surfbeat: geographic advection of the wave groups: 1 first-order upwind, 2 second order (van Leer)

    // NHFLOW Lagrangian particles
    L10=0;       // int particle tracking
    L11=1;       // int random seed
    L21=0;       // int point releases
    L22=0;       // int line releases
    L23=0;       // int box releases
    L31=0.0;     // double horizontal diffusivity Kh
    L32=0;       // int vertical diffusion: 0 off, 1 constant Kv, 2 Kv + nu_t/Sc
    L33=0.0;     // double vertical (background) diffusivity Kv
    L34=1.0;     // double turbulent Schmidt number
    L41=0;       // int windage wind given
    L41_u=0.0;   // double U10 for windage
    L41_dir=0.0; // double wind direction for windage
    L51=0;       // int beaching: 0 refloat when rewetted, 1 permanent
    L52=1;       // int bed contact: 0 reflect, 1 deposit sinking particles
    L53=0;       // int domain boundaries: 0 remove, 1 reflect
    L61=1.0;     // double print interval
    L62=1;       // int output: 1 vtp, 2 csv, 3 both


    // Boundary Conditions
	B10=0;			// int wall function velocities on/off
    B11=1;			// int wall function turbulence on/off
	B20=2;			// int slip or no-slip boundary condition for velocity
    B21=1;            // int preserve tangential velocities for forcing
    B22=2;            // int smooth vs rough wall law
    B23=1;            // int ghostcell extrapolation or reflective
	B29=0.5;		    // double gamma for gc image point
	B30=0;			// int type of pressure reference point
    B31=0.0;         // double pressure reference value
    B32=0;           // int pressure reference location
    B32_x=B32_y=B32_z=0.0; // double pressure reference location
    B33=1;           // int pressure gage virtual or inline
	B50=0.0001;		// double global wall roughness ks
	B51=-1.0;		// double global wall roughness ks
	B52=-1.0;		// double global wall roughness ks
	B53=-1.0;		// double global wall roughness ks
	B54=-1.0;		// double global wall roughness ks
	B55=-1.0;		// double global wall roughness ks
	B56=-1.0;		// double global wall roughness ks
    B57=0.0001;       // double global side wall roughness
	B60=0;            // int ioflow discharge
	B61=2;            // int plain or logarithmic inflow profile
    B62=0;            // int SFLOW/NHFLOW discharge inflow across the face: 0 as B 61, 1 by conveyance (depth average u ~ h^(2/3))
	B71=0;       // double distance for use relaxation method for fixed water level ini
	B75=1;		 // int type of outflow boundary conditions
    B76=1;      // int type of pressure inlet boundary condition
	B77=1;           // int outflow pressure controlled or free stream
	B81=0;		 // int focussed wave parameter
    B81_1=B81_2=0.0;   // double focused wave: focus location x, focus time
    B81_3=0.0;    // double unidirectional focused wave y is condisered 0
    B82=1;      // int type of focus point and time calculation
	B83=0.0025;  // double wave steepness parameter for focused breaking waves
    B84=1;       // Peak enhance method
	B85=1;		 // int PM or JONSWAP spectrum for irregular waves
	B86=10;		 // int number of regular waves for irregular wave generation
	B87=0;		 // int give ws and we for irregular wave generation
	B87_1=B87_2=0.0;  // double ws, we (B 87)
	B88=3.3;	 // double gamma for JONSWAP spectrum
	B89=0;            // int wave generation optimization
	B90=0;            // int iowave
	B91=0;          // int wave parameter wL
	B91_1=0.0;      // double wave amplitude
	B91_2=0.0;      // double wave length
	B92=0;            // int wave type
	B93=0;          // int wave parameter wT
	B93_1=0.0;      // double wave amplitude
	B93_2=0.0;      // double wave period
    B94=0;     // int set water depth for wave theory
    B94_wdt=0.0;    // double water depth for wave theory
    B95=0.0;        // double SFLOW and NHFLOW relaxation zones (B 98 2): 0 as before; >0 step-size consistent, reference step [s]: targets at the time of every stage output, relaxation factor r^(dt/B95) per stage; -1 the same with the reference step of the finest level (dt/2^G1 with G 7 1, else dt)
	B96_1=0.0;      // double dist1 for wave relax
	B96_2=0.0;      // double dist2 for wave relax
	B97=0;          // int NHFLOW numerical beach (B 99 1|2) relaxes to the inflow current (waves on a current)
	B98=0;          // int type of wave generation
	B99=0;			// int type of numerical beach
	B101=0;        // int ramp function wave geneartion
	B102=1.0;        // double factor ramp function wave generation
    B105=0;            // int wave generation origin changed
 	B105_1=0.0;        // double wave generation direction and line origin in Cartesian coordianate system
	B105_2=0.0;        // double wave generation line x origin
	B105_3=0.0;        // double wave generation line y origin
	B106=0;			// int read wave generation origin
	B107=0;			// int read numerical beach origin
    B108=0;        // int read wave generation  origin
    B500=B501=B502=B504=B520=B521=B524=0;   // int iowave redesign: wave sources (B 500-504), zones (B 520-524)
    B525=0;        // int iowave redesign: discharge of a clamped discharge edge (B 525 id Q t_ramp), repeatable
    B522=0;        // int multi-block grid (G 50 1): zones on a face of a block 1.. (B 522 id b: the B 521 edge of zone id is a face of block b), repeatable
    B513=0;          // int iowave redesign: vertical profile of the background current (B 513 id profile par), repeatable
    B505=0;          // int iowave: frame of a wave source (B 505 id frame; 0: B 92 generation frame, 1: global)
    B510=B511=B514=B523=0;   // int iowave redesign: tidal / current background (B 510-514), zone background (B 523)
    B515=0;          // int iowave: number of B 515 entries (counted by read_control)
    B530=-1;       // int iowave redesign: waves on the background, 1: on h_eff, 2: h_eff + Doppler, -1: auto (2 with a current, 1 with a level only)
    B530_N=10;     // int iowave redesign: k update interval [steps] of B 530
    B529=0;        // int iowave redesign: edge mass balance log every B 529 steps (0: off)
    B526=0.0;      // double iowave: low-pass time scale of the combined beach B 99 6 (0: peak period of the B 92 wave / 4)
    B527=0.0;      // double iowave: low-pass time scale of the combined generation B 98 6 (0: peak period of the B 92 wave / 4)
    B528=0.0;      // double iowave: reference period of the two-band absorption B 98 7 / B 99 7 (0: peak period of the B 92 wave)
    B534_d=B534_dx=B534_fl=B534_fh=0.0;   // double iowave: SIRW absorption at x- (B 98 8), 0: automatic from the reference period
    B535_d=B535_dx=B535_fl=B535_fh=0.0;   // double iowave: SIRW absorption at x+ (B 99 8), 0: automatic from the reference period
    B536=0;        // int iowave NHFLOW: wavemaker absorption target, 0 shallow water (legacy), 1 linear paddle theory
    B110=0;        // int read wave generation  origin
    B111_zs=0.0;	// double flap start
    B111_ze=0.0;	// double flap end
    B112_zs=0.0;	// double flap start
    B112_z2=0.0;	// double flap2 end/flap2 start
    B112_ze=1.0;	// double flag end
    B113=0;         // int 2nd-order wavemaker correction for piston/flap: 0 off, 1 sub+super, 2 sub, 3 super
    B113_J=10;      // int number of evanescent modes for the 2nd-order correction
    B114_fmin=0.0;  // double lower frequency of first-order components [Hz]
    B114_fmax=0.0;  // double upper frequency of first-order components [Hz], 0: automatic
    B114_f2min=0.0; // double lowest corrected 2nd-order frequency [Hz], 0: 0.1 x peak frequency
    B115=0;         // int activate vertical velocity component for flap wavemaker theory
    B116=1;         // int x or beta input for flap wavemaker theories
    B117=0.0;		  // double starting time shift for timeseries input
    B118=1.0;        // double wavemaker correction factor for NHFLOW
    B119=0;          // int FNPF paddle BC with 2nd-order Taylor terms (moving paddle) for piston/flap
    B120=-90.0;       // doubel delta t for wave generation
    B123=0.0;       // double flap AWA hinge location
    B125=0;         // int take 2D slice input for HDC
    B125_y=0.0;     // double 2D slice y-coor input for HDC
    B127=0;         // int turn of y-dir velociteis for HDC 
    B130=0;         // int directional spreading for irregular waves
    B131=0.0;       // double main direction for multidirectional irregular waves
    B132_s=-90.0;  // double start directional spreading
    B132_e= 90.0;  // double end directional spreading
    B133=1;         // int number of direction intervals for spreading function
    B134=-1.0;      // double shape parameter for spreading function
    B135=10.0;       // double peak value
    B136=1;         // int double summation method frequency vector
    B138=0;         // int seed number multidir waves
    B138_1=B138_2=0;  // int seeds (B 138)
    B139=0;         // int seed number wave spectrum
    B160=5;        // int number of vertical layers for 2D wave generation
    B170=1024;     // int number of Fourier modes for the generation of steady surface gravity waves
	B180=0;           // int gravity waves
    B181=0;         // int x-dir motion
	B181_1=0.0;     // double x-acceleration amplitude
	B181_2=0.0;      //double x-acceleration frequency
	B181_3=0.0;     // double wave phase change
    B182=0;         // int y-dir motion
	B182_1=0.0;     // double y-acceleration amplitude
	B182_2=0.0;      //double y-acceleration frequency
	B182_3=0.0;     // double wave phase change
    B183=0;         // int z-dir motion
	B183_1=0.0;     // double z-acceleration amplitude
	B183_2=0.0;      //double z-acceleration frequency
	B183_3=0.0;     // double wave phase change
	B191=0;			// int rotation around x-axis
	B191_1=0.0;		// double angle for rotation around x-axis
	B191_2=0.0;		// double frequency for rotation around x-axis
	B191_3=0.0;		// double y-coordinate for rotation around x-axis
	B191_4=0.0;		// double z-coordinate for rotation around x-axis
	B192=0;			// int rotation around y-axis
	B192_1=0.0;		// double angle for rotation around y-axis
	B192_2=0.0;		// double frequency forrotation around y-axis
	B192_3=0.0;		// double x-coordinate for rotation around y-axis
	B192_4=0.0;		// double z-coordinate for rotation around y-axis
    B194_s=-1.0e9; // double start rotation
	B194_e= 1.0e9; // double end rotation
    
    B200=0;			// int VRANS on/off -> assigned as 1 for VRANS Structure, 2 for Vegetation, 3 for Net interaction
    B201=0;          // int porosity
    B208=0;          // int print porous force
    B209=1.1;             // double eps for porous heaviside
    B201_n=1.0;         // double porosity n
    B201_d50=0.01;      // double porosity d50
    B201_alpha=0.0;     // double porosity alpha
    B201_beta=0.0;      // double porosity beta
    B202=0;             // int number of porous layers (NHFLOW), counted from the structure surface inwards
    B203=1;             // int VRANS surface: 1 faces resting on the bed do not define the porous interface / layer depth, 0 all faces
    B204=1;             // int NHFLOW VRANS drag: 0 explicit, 1 point-implicit with the stage update
    B210=0;      // int solid box
    B212=0;      // int solid cylinder y
    B213=0;      // int solid cylinder z
    B214=0;      // int solid cylinder jacket member
    B215=0;      // int solid sphere
    B217=0;      // int solid wedge x
    B218=0;      // int solid wedge y
    B219=0;      // int solid wedge z
    
    B230=0;      // int STL
    B231=0;      // int scaling of stl geometry on/off (set by B 231)
    B231_x=B231_y=B231_z=1.0;  // double scaling of stl geometry
    B232=0;     // int translation on/off
    B232_dx=B232_dy=B232_dz=0.0;  // double translation of stl geometry
    B233=0;
    B233_x=B233_y=B233_z=B233_phi=B233_theta=B233_psi=0.0;

    
	B240=0;			// int porous media
	B241=1;			// int porous media in x-direction
	B242=1;			// int porous media in y-direction
	B243=1;			// int porous media in z-direction
    B260=0.34;       // double C coefficient for VRANS
    B264=1.0e20;    // double KC number for VRANS
    B265=0;         // int NHFLOW VRANS inertia factor: 0 1/(1+C(1-n)/n^2), 1 n/(1+C(1-n)/n)
    B266=1;         // int CFD VRANS added mass: 0 C_m = B 260 in the porous zone, 1 C_m = B 260 (1-n)/n
    B268=1;         // int CFD VRANS drag: 0 explicit source, 1 point-implicit in the RK stages (momentum_rk)
    B270=0;         // int VRANS porous media box
    B274=0;         // int VRANS porous media vertical cylinder
    B281=0;         // int VRANS porous media wedge in x-direction
    B282=0;         // int VRANS porous media wedge in y-direction
    B291=0;         // int VRANS porous media plate in x-direction
    B295=1;         // int VRANS source treatment fo turbulence model
    B307=0;         // int VRANS interface smoothing interations
    B308=1;         // int porosity effects on fluid acceleration for vegetation
    B309=2.0;       // double Cm for vegetation
    B310=0;         // int VRANS vegetation box
    B311=0;         // int VRANS vegetation
    B321=0;         // int VRANS vegetation wedge in x-direction
    B322=0;         // int VRANS vegetation wedge in y-direction
    
    B411=0;        // int patchBC discharge
    B412=0;        // int patchBC pressure BC
    B413=0;        // int patchBC waterlevel
    B414=0;        // int patchBC perpendicular velocity
    B415=0;        // int patchBC velocity components
    B416=0;        // int patchBC horizontal inflow angle
    B417=0;        // int patchBC inflow normals
    B418=0;        // int patchBC outflow pressure condition
    B421=0;        // int patchBC hydrograph discharge
    B422=0;        // int patchBC hydrograph waterlevel
    B440=0;        // int patch BC inflow line
	B441=0;			// int rectangular inflow patch BC
    B442=0;			// int circular inflow patch BC

	// Concentration
	C1=10.0;		// double density concentration in water
	C2=0.0;			// double viscosity water + concentration
	C3=10.0;		// double density concentration in air
	C4=0.0;			// double viscosity air + concentration
	C10=0;          // int concentration transfer on/off
	C15=0;          // int concentration convection
	C20=0;          // int concentration diffusion
	C50_1=1.0; 		// double fill ration concentration area 1
	C50_2=0.0; 		// double fill ration concentration area 2
	C51=-1.0e7;		// double i-dir zero level set start
	C52=-1.0e7;		// double j-dir zero level set start
	C53=-1.0e7;		// double k-dir zero level set start
	C54=1.0e7;		// double i-dir zero level set end
	C55=1.0e7;		// double j-dir zero level set end
	C56=1.0e7;		// double k-dir zero level set end
	C57_1=0.0;      // a, plane
	C57_2=0.0;      // b
	C57_3=0.0;      // c
	C57_4=0.0;      // d
	C58_1=0.0;      // x0, sphere
	C58_2=0.0;      // y0
	C58_3=0.0;      // z0
	C58_4=0.0;      // r
	C75=0;			// int number of tiltboxes

    // Discretization
	D10=4;			// int convection scheme
	D11=2;			// int convection velocity scheme
    D12=0;          // int WENO5 nonlinear weights for momentum (D 10 4,5): 0 WENO-JS, 1 WENO-Z, 2 TENO5
    D13=1.0e-5;     // double TENO5 cutoff C_T (D 12 2)
	D20=2;			// int diffusion scheme
	D21=0;			// int print out implicit diffusion time and iterations
    D22=1;            // int diffusion wall boundary condition
    D23=2;            // int time accuracy of the implicit momentum diffusion (D 20 2): 1 first order, 2 second order
	D30=1;			// int pressure scheme

    // Free Surface
	F30=3;			    // int level set scheme
	F35=5;			    // int convection scheme for fsf
    F37=0;              // int WENO5 nonlinear weights for level set convection (F 35 4,5; F 305 4,5) and reini (F 40 3,23): 0 WENO-JS, 1 WENO-Z, 2 TENO5
    F38=1.0e-5;         // double TENO5 cutoff C_T (F 37 2)
	F40=3;			    // int reini scheme
	F42=-1.0;		// double maxlength
	F43=0.55;		// double factor for reini timestep
	F44=3;		        // int number reini time step
	F45=2.1;         // factor for calculation of epsi
	F46=0;            // int picard iteration for lsm or reini
	F47=10;            // int number of picard iterations
	F50=2;            // int bc phi, 1: inflow or 2: outflow
    F50_flag=0;       // int flag for lsm description
	F51=-1.0e20;		    // double i-dir zero level set start
	F52=-1.0e20;		    // double j-dir zero level set start
	F53=-1.0e20;		    // double k-dir zero level set start
	F54=1.0e20;		// double i-dir zero level set end
	F55=1.0e20;		// double j-dir zero level set end
	F56=1.0e20;		// double k-dir zero level set end
	F57_1=0.0;      // a, plane
	F57_2=0.0;      // b
	F57_3=0.0;      // c
	F57_4=0.0;      // d
	F58_1=0.0;      // x0, sphere
	F58_2=0.0;      // y0
	F58_3=0.0;      // z0
	F58_4=0.0;      // r
    F59_xm=0.0;      // xm, cylinder
	F59_ym=0.0;      // ym
	F59_zs=0.0;      // zs
	F59_ze=0.0;      // ze
    F59_r=0.0;      // r
    F60_flag=0;     // int F60 flag
	F60=-1.0e20;  // double ini z-dir
	F61=-1.0e20;  // double inflow  ini
	F62=-1.0e20;  // double outflow  ini
	F63=-1.0e20;  // double xstart phi interpolate with outflow h
	F64=0;			// int iterations for relaxing outflow waterlevel
    F65=0;            // int SFLOW river free-surface ini along i: 1 linear F 61 to F 62, 2 steady 1D flow of W 10 from F 62
    F68=0;            // int seaward free surface level 
	F70=0;             // int number of phi 1 ini boxes
	F71=0;             // int number of phi 2 ini boxes
	F72=0;             // int number of phi 1 ini regions
	F80=0;             // int time scheme VOF
    F85=0;             // int convection scheme VOF
    F88=10;             // int normal calculation method for VOF PLIC,2D MYC as default
    F90=4;             // int vof transport scheme for VOF PLIC, 2D COSMIC as default
    F91=0.1;          // rho threshold input for FCCRK3 PLIC
    F92=2;              // PLIC density method, fine vof-based as default
    F93=0.001;          // PLIC air phase threshold
    F94=0.999;          // PLIC water phase thresold
    F95=0.999;          // PLIC threshold for water in phasemarker functions
    F96=0;           // PLIC density softening from superfine to fine x direction start & end
    F96_xs=-1.0e20;
    F96_xe=-1.0e20;
    F97=0;           // PLIC density softening from superfine to fine y direction start & end
    F97_ys=-1.0e20;
    F97_ye=-1.0e20;
    F98=1; // Spray deletion algorithm
    
    F100=0; // Csaba's Similiarity runup case setup
    F101_nx=0.0; //vof ini for Csaba's similiarity runup case
    F101_ny=0.0;
    F101_nz=0.0;
    F101_s=0.0;
    F101_e=0.0;
    
    F112=0;            // int wedge x-dir
    F113=0;            // int wedge y-dir
    F114=0;            // int inv. wedge x-dir
    F115=0;            // int inv. wedge y-dir
    F150=0;         // int benchmark
    F151=0;         // int benchmark inverse sign of level set
    F300=0;             // int multiphase flow level set
    F305=5;             // int multiphase flow lsm convection
    F310=3;             // int multiphase flow reini
    F321=1.6;         // double epsi12
    F322=1.6;         // double epsi13
    F323=1.6;         // double epsi23
    F360=-1.0e20;  // double ini x-dir ls1
    F361=-1.0e20;  // double ini y-dir ls1
    F362=-1.0e20;  // double ini z-dir ls1
    F369=0;             // int number of phi 1 ini tiltboxes ls1
    F370=0;             // int number of phi 1 ini boxes ls1
    F371=0;             // int number of phi 2 ini boxes ls1
    F374=0;             // int number of pos ls1 ycyl
    F375=0;             // int number of neg ls1 ycyl
    F378=0;             // int number of pos ls1 sphere
    F379=0;             // int number of neg ls1 sphere
	F380=-1.0e20;  // double ini x-dir ls2
	F381=-1.0e20;  // double ini y-dir ls2
	F382=-1.0e20;  // double ini z-dir ls2
	F390=0;             // int number of phi 1 ini boxes ls2
	F391=0;             // int number of phi 2 ini boxes ls2
    F394=0;             // int number of pos ls2 ycyl
    F395=0;             // int number of neg ls2 ycyl
    F398=0;             // int number of pos ls2 sphere
    F399=0;             // int number of neg ls2 sphere

    // Grid
    // mesh refinement (SFLOW, FNPF, NHFLOW)
    G1=0;       // int SFLOW, FNPF and NHFLOW mesh refinement: number of refined levels
    G2=4;       // int SFLOW, FNPF, NHFLOW and CFD mesh refinement: regrid interval in time steps (0: static; FNPF and NHFLOW: with G 12, G 20, G 21, G 22 or G 23; CFD: with G 12, G 24, G 25 or G 26)
    G3=2;       // int SFLOW, FNPF, NHFLOW and CFD mesh refinement: buffer cells around flagged cells
    G4=8;       // int SFLOW, FNPF, NHFLOW and CFD mesh refinement: tile size in cells of the refined level
    G5=4;       // int SFLOW, FNPF, NHFLOW and CFD mesh refinement: regrids a refined tile is kept after its last flag (hysteresis)
    G6=0;       // int FNPF and NHFLOW mesh refinement: 1 doubles the sigma layers on every refined level
    G7=0;       // int SFLOW, NHFLOW and CFD mesh refinement: 0 one time step for all levels, 1 subcycling (every level two steps per step of the next coarser one)
    G10=0;      // int SFLOW, FNPF and NHFLOW mesh refinement: number of static refinement boxes
    G15=0;      // int CFD mesh refinement: number of static refinement boxes xs xe ys ye zs ze
    G11=0;      // int SFLOW, FNPF and NHFLOW mesh refinement: number of boxes without refinement
    G12=0;      // int SFLOW, FNPF, NHFLOW and CFD mesh refinement: refine around the moving body (SFLOW X 10 2/3, FNPF X 10 1, NHFLOW X 10 1/2, CFD X 10 1) and around the fixed solids (NHFLOW A 581-590)
    G12_r=0.5;  // double SFLOW, FNPF, NHFLOW and CFD mesh refinement: margin around the hull
    G13_L=0.0;  // double SFLOW, NHFLOW and CFD mesh refinement: length of the refined wake wedge from the bow (0: off)
    G13_a=19.47;// double SFLOW, NHFLOW and CFD mesh refinement: half angle of the wake wedge in degrees
    G20=0.0;    // double SFLOW, FNPF and NHFLOW mesh refinement: surface jump between neighbour cells that flags a cell (0: off)
    G21=0.0;    // double FNPF and NHFLOW mesh refinement: second difference of the surface along x or y that flags a cell (0: off)
    G22=0;      // int SFLOW, FNPF and NHFLOW mesh refinement: flag the shoreline (SFLOW: >0 on; NHFLOW and FNPF with G 30 1: the cells within this many cells (1-3) of the shoreline); 0 off
    G23=0;      // int NHFLOW mesh refinement with A 550 1, FNPF with A 350 1: 1 flags the breaking cells (the breaking zone on the finest grid)
    G24=0.0;    // double CFD mesh refinement: flag the cells with |phi| below this many interface thicknesses (psi of the flagging level); 0 off
    G25=0;      // int CFD mesh refinement: flag the cells within this many cells of a solid, the topography or a floating body; 0 off
    G26=0.0;    // double CFD mesh refinement: flag the cells with a vorticity above this value [1/s]; 0 off
    G30=0;      // int NHFLOW and FNPF mesh refinement: 1 patches may cover dry and shallow cells (wetting and drying in the patches; FNPF: the static coastline A 343 1 or the moving one A 343 2, 3), 0 they keep 4 level-0 cells away
    G31=1;      // int NHFLOW mesh refinement with A 512 2: 1 one implicit diffusion solve over the leaf cells of all grids (composite), 0 every grid solves its own
    G52=0;      // int curvilinear NHFLOW (G 50): 1 the i-walls are slip walls as the j-walls (normal velocity mirrored, the pressure correction Neumann), so that a wall acts the same in any index direction (set on multi-block grids); 0 the Cartesian x-walls
    G54=0;      // int multi-block grid (G 50 1): number of open boundaries on the blocks 1.. (G 54 b t v: block b, t 1: discharge v in through its -i face, t 2: water level v at its +i face, t 3: discharge time series discharge-<v>.dat through the -i face, t 4: water level time series waterlevel-<v>.dat at the +i face)
    G55=0;      // int multi-block grid (G 50 1): 0 the blocks 1.. of DIVEMesh as they are, m > 0 cut into pieces of at most (cells of the blocks 1..) / (m ranks) cells, so that they spread over the ranks (geo_mesh_split.cpp)
    G56=0;      // int multi-block grid (G 50 1): G 56 a dz places the blocks on the ranks and cuts them (G 55) by their load: a wet cell 1, a cell whose bed lies above the still water level (F 60, else F 62 with F 65) + dz the load a (Tagus, SI: a = 0.66 measured); the placement is improved by moves and swaps of blocks between the ranks; 0: by the number of cells, largest block first
    G56_a=1.0;  // double multi-block grid: G 56 a dz, a
    G56_z=0.0;  // double multi-block grid: G 56 a dz, dz (m)
    G53=4;      // int NHFLOW multi-block grid: coarse space across the blocks in the preconditioner of the composite solves, aggregates of this many columns in i and j (0: off)
    G51=2;      // int NHFLOW multi-block grid (G 50 1, DIVEMesh R 20 1): passes of the block-local V-cycles in the preconditioner of the pressure, the corrections of the neighbour blocks filled in between
    G50=0;      // int curvilinear horizontal grid (SFLOW): 0 Cartesian, 1 nodes of the CURV section of grid-geometry.dat, 2 nodes of curvgrid.dat
    G40=0;      // int mesh refinement on several ranks (NHFLOW, SFLOW, FNPF, CFD): 0 patches on the rank of the level-0 cells below, 1 placed for the load of the ranks, 2 test of 1
    G42=0.0;    // double CFD mesh refinement with floating bodies (X 10 1): a patch casts the bodies again when a vertex has moved more than G 42 times its cell size since its last cast (the stages in between keep that level set; their forcing velocity is that of the stage); 0 at every change of the pose
    G43=0;      // int mesh refinement (SFLOW, FNPF, NHFLOW, CFD): load statistics of the ranks every G 43 steps (busy and MPI time, patch and regrid time, predicted and measured load, the measured cost of a patch cell; REEF3D_<M>_AMR/REEF3D_<M>_AMR_load*.dat); 0 off
    G41=1.0;    // double mesh refinement on several ranks (NHFLOW, SFLOW, FNPF, CFD): the work of a patch cell against a level-0 cell in the placement of G 40 1, 2; 0 measured in the run (CFD: the timers of the ranks at every regrid, with a report of the measured load, and for the load map of a rebalance (G 46) the busy times of the ranks fitted to their level-0 and patch cells; NHFLOW, SFLOW, FNPF: the measured cost model, the busy time of the ranks fitted to their level-0 cells that are not solid, their patch cells and the parent cells they serve, at every regrid; a static refinement is placed once more at step 15)
    G45=0;      // int NHFLOW, SFLOW, FNPF and CFD mesh refinement: the load map every G 45 steps, the work per step of every level-0 column (CFD: cell, averaged over the window of G 47 steps) averaged over the run (its level-0 cells that are not solid plus the weight of a patch cell (G 41, or measured) times the patch cells over it), for the weighted partition of DIVEMesh (M 20 5, M 35 1); REEF3D_<M>_AMR/REEF3D_<M>_AMR_loadmap.dat; SEASTATE (stand-alone): the measured work of the cells per iteration of every level-0 column [microseconds], the patch cells over it included, REEF3D_SEASTATE_Log/REEF3D_SEASTATE_loadmap.dat; CFD with CPM (no refinement): a map of the cells, 1 per fluid cell plus the measured cost of the parcel updates in the cell, REEF3D_CFD_CPM_Particle/REEF3D_CFD_CPM_loadmap.dat; 0 off
    G46=0.0;    // double SFLOW, FNPF, NHFLOW, CFD, SEASTATE (stand-alone) on several ranks: rebalance when the measured load max/mean of the ranks (busy time, the MPI waits taken out; SEASTATE: the time of the solves of its cells) exceeds G 46 (SEASTATE: and the weighted partition lowers it by 10 % by its load map): checkpoint, load map, stop with exit code 3; the next start (tools/reef3d_rebalance.sh) decomposes with the weighted partition and the load map and continues (needs DIVEMesh M 50 1, 2); 0 off
    G47=100;    // int G 46: steps between the checks of the load (and at least before the first); SEASTATE stationary: iterations of the step
    G44=0.0;    // double NHFLOW, SFLOW and FNPF mesh refinement with G 40 1: weight of the patch phase in the placement (the cost of a placement is the maximum work of the ranks plus G 44 times their maximum patch work); 0 the total work only
    // grid
    G502=0;            // int sigma grid
    G505=0;             // int turn of topo and solid cells
    G540=3;         // int reini scheme for topo/solid/floating
    G541=2;          // int reini iterations

	// Heat
	H1=1.4e-7;      // thermal diffusivity water
	H2=2.216e-5;      // thermal diffusivity air
    H3=1;            // type of density calculation
    H4=0;           // int use beta coeff
    H4_beta1=0.0;   // double beta1
    H4_beta2=0.0;   // double beta2
    H9=1;           // int air-water assignment
	H10=0;          // int heat transfer on/off
    H15=5;          // int convection for heat transfer
	H50_1=20.0;     // double temperature 1
	H50_2=20.0;     // double temperature 2
	H51=-1.0e7;		    // double i-dir zero level set start
	H52=-1.0e7;		    // double j-dir zero level set start
	H53=-1.0e7;		    // double k-dir zero level set start
	H54=1.0e7;		// double i-dir zero level set end
	H55=1.0e7;		// double j-dir zero level set end
	H56=1.0e7;		// double k-dir zero level set end
	H57_1=0.0;      // a, plane
	H57_2=0.0;      // b
	H57_3=0.0;      // c
	H57_4=0.0;      // d
	H58_1=0.0;      // x0, sphere
	H58_2=0.0;      // y0
	H58_3=0.0;      // z0
	H58_4=0.0;      // r
    H61=0;          // int heat bc
    H62=0;          // int heat bc
    H63=0;          // int heat bc
    H64=0;          // int heat bc
    H65=0;          // int heat bc
    H66=0;          // int heat bc
    H61_T=0.0;      // double heat bc
    H62_T=0.0;      // double heat bc
    H63_T=0.0;      // double heat bc
    H64_T=0.0;      // double heat bc
    H65_T=0.0;      // double heat bc
    H66_T=0.0;      // double heat bc

    // Initialize
    I10=0;          // int initialize all
	I11=0;			// int initialize velocities with potential flow
	I12=0;          // int initialize pressure
	I13=0;          // int initialize turbulence
	I21=0;          // int set phase 2 velocities to zero after potential flow solver
	I30=0;			// int Fully intialize NWT
    I31=-1.0;       // double I 30 2: exposure mask smoothing length (< 0: automatic, 1 peak wave length)
    I32=32;         // int I 30 2: number of direction bins of the exposure masks
    I33=-1.0;       // double I 30 2: cell size of the coarse exposure mask grid (< 0: automatic)
    I35=15.0;       // double I 30 2: angular blur of the exposure masks across directions [deg] (0: off)
    I34=0.159155;   // double I 30 2: diffraction spreading coefficient, variance = I34*Lp*distance (1/(2 pi); 0: off)
	I40=0;			// int ini from state file
	I41=0;			// int ID of state file
    I44=0;          // int FNPF state with Fi
    I45_count=0;    // int I 45 <count> (count > 1): from the checkpoint REEF3D_Checkpoint/<count> instead of the latest (control::parse)
    I45=0;          // int SFLOW, FNPF, NHFLOW, CFD, SEASTATE (stand-alone): 1 restart from the latest checkpoint (REEF3D_Checkpoint/latest, P 48; it may hold a bare count) on any decomposition, <count> (> 1) from REEF3D_Checkpoint/<count> (SEASTATE stationary: the step continues at the iteration of the checkpoint)
    I50=0;			// double simtime ini
    I55=0.0;        // double reference pressure
	I56=0;          // int pressure above F56 set to zero
	I58_1=0.0;      // double vertical velocity for sphere initialization
	I58_2=0.0;      // double radius for sphere initialization
    I230=0;         // int read 2D flowfile
    I231=0.0;       // double starting x for flowfile
    I232=0.0;       // double starting y for flowfile
    I233=0.0;       // double starting z for flowfile
    I241=0.0;       // double delta t for flowfile
    I242=3;         // int FNPF coupling (B 92 62): time interpolation, 1 linear, 3 cubic

    // Numerics
	N10=1;			// int linear poisson solver: 1 REEFMG; hypre (10-39) only in builds with make HYPRE=1
	N11=11;         // int precondioner
    N12=0;         //  int REEFMG solver mode
    N13=0;         //  int REEFMG max number of coarse levels, 0: automatic
    N14=32;         //  int REEFMG V-cycle precision: 32 single, otherwise double
    N18=1;         //   int permuting axes for hypre
    N20=20;        //  int precondioner switch iter
    N22=0;         //  int precondioner switch SFLOW
    N15=0;         //  int REEFMG coarse-grid agglomeration: 1 on
    N23=0;         //  int precondioner switch FNPF
    N25=0;         //  int precondioner switch NHFLOW
    N26=1;         //  int precondioner switch CFD
	N40=3;			// int time scheme
	N41=1.0e+19; 	// double total time
	N43=1.0e-5;     // double stopping criteria convection-diffusion
	N44=1.0e-8;     // double stopping criteria pressure
	N45=1e8;		// max outer iter
	N46=250;		// int max number of solver iterations
	N47=0.3;		// doubel relaxation factor for time stepping
	N48=1;          // int adaptive timestepping
	N49=1.0;		// double max timestep or fixed timesteps
    N50=1;          // int adaptive timestepping method
	N61=500.0;      // double stopping criteria velocities

    // MPI
	M10=0;          // number of MPI processes
	M41=0;          // int M 40: 1 the weighted partition with the load map (DIVEMesh M 20 5, M 35 1) whatever DIVEMesh M 20 is; set by a rebalance restart (G 46)
	M40=0;          // int grid: 0 the grid files of the ranks from DIVEMesh, else its global grid file (DIVEMesh M 50 1, 2) decomposed by REEF3D; 1 always the global grid file, REEF3D decomposes for M 10 ranks with the decomposition keys of DIVEMesh (M 20, ...)

    // Printer
	P10=1;			 // int print file type
	P12=1;			 // int terminal print frequency
	P15=1;          // int print file numbering
    P16=1;          // int add timestamp (and run id) to paraview files
    P18=0;          // int LAGOON store: 1 instead of the vtu/vtp files, 2 as well as them
    P19=2;          // int print NHFLOW floating  based on vtp or vtu interval
	P20=-10;		// ith iteration file printed
    P21=0;          // int time averaged vtu print out
    P22=0.0;         // double start averging after transients
	P23=0;

			// int print test to vtu file
    P24=0;			// int print density to vtu file
    P25=0;			// int print solid to vtu file
	P26=0;			// int print cbed and conc to vtu file
	P27=0;			// int print topo to vtu file
	P28=0;			// int print fb to vtu file
	P29=0;			// int print walldist to vtu file
	P30=-1.0;       // double time between file printout in seconds
	P34=-1.0;       // double time between file printout in seconds for sediment
	P35=0;        	// int print for interval
    P37=1;			// int print restart state file
	P38=100;			// int print restartstate file each ith iteration
	P40=0;				// int print state file
	P41=1;			// int print state file each ith iteration
	P42=-1.0;			// double print state file each ith sec
    P43=0;             // int state print out selected area
    P44=1;             // int state print out components for FNPF
    P45=2;             // int print into single or continous state file
    P46=0;             // int print state iteration window
    P47=0;             // int print state time window
    P48_dt=0.0;        // double P 48 N dt: also a checkpoint every dt seconds of simulated time (at the first step at or past each multiple of dt); P 48 0 dt: by time only
    P49=-1;            // int the newest P 49 checkpoints are kept (REEF3D_Checkpoint/history), older ones removed; 0: all kept; not given: 2 (SEASTATE 1)
    P48=0;             // int SFLOW, FNPF, NHFLOW, CFD, SEASTATE (stand-alone): a checkpoint every P 48 steps (REEF3D_Checkpoint/<count>, for I 45 on any decomposition); SEASTATE stationary: every P 48 iterations of the step (REEF3D_Checkpoint/<count>_<iteration>); kept: P 49; 0 off
    P50=0;				// int wave theory wave gages
	P51=0;             // int print out wsf
	P52=0;            // int print out wsfline in x-dir
	P53=0;            // int print out wsfline for wave theory
	P54=10;			  // int ith iteration wsfline file  print out
	P55=-1.0;		  // double ith second wsfline files print out
	P56=0;            // int print out wsfline in y-dir
    P57=0;            // int add aditional info to WSF gage in FNPF
    P58=0;            // int print wave time series
    P59=0;              // int print breaking wave log FNPF
    P61=0;              // int print point probes
    P62=0;              // int print line probes
    P63=0;              // int print depth averaged point probe
    P64=0;              // int print pressure probes
    P65=0;              // int print velocity probes
    P66=0;              // int print velocity probes from wave theory
    P67=0;              // int print vertical velocity profile in x-direction
    P68=0;              // int print vertical velocity profile in y-direction
    P69=0;              // int print vertical velocity profile in  horiontal direction
    P71=0;           // int print viscosity to vtu
    P72=0;           // int print vof function
    P73=0;           // int print hx and hy for sflow vtp
    P74=0;           // int print omega nhflow
	P75=0;            // int print out vorticity vec
    P76=0;            // int print out bedload
    P77=0;            // int print out sediment parameters: 1
    P78=0;            // int print out sediment parameters: 2
	P79=0;            // int print out bed shear stress when running sediment transport
	P80=1;            // int force print out interval
    P81=0;            // int force print out
    P82=0;            // int add eddyv to viscous force
	P85=0;            // int ALE force print out for FNPF
    P88=0;            // int kinematics print out for FNPF
	P91=0.25;		  // double factor used in force calculation algorithm
    P92=0;           // int force from water or from water+air
	P101=0;			  // int print sloshing forces
    P110=0;           // int print significant wave height
    P111=0.0;         // double start averging after transients
    P112=0;           // int SFLOW maxima of the run (REEF3D_SFLOW_MAX)
    P113=0.0;         // double SFLOW maxima: start time (after the spin-up)
    P114=0.05;        // double SFLOW maxima: flooding depth
    P115=0;           // int SFLOW maxima at points
	P121=0;             // int bed level gages
	P122=0;             // int max bed level gages
	P123=0;             // int topoline in x-direction
	P124=0;             // int topoline in y-direction
	P125=0;             // int bed shear stress gages
	P126=0;             // int bed shear stress maxval
    P131=0;             // int max wetdry in vtp
    P133=0;             // int runup gage x-crossection
    P134=0;             // int runup gage y-crossection
    P140=0;             // int runup gage cylinder
    P141=0;             // int runup cylinder radius
    P144=0;           // int print NHFLOW depth-averaged velocity line probes
    P145=0;           // int print NHFLOW depth and time-averaged velocity line probes
    P146=0;           // int print NHFLOW time-averaged wsf line in x-direction
    P147=0;           // int print NHFLOW time-averaged wsf line in y-direction
    P148=0;           // int print NHFLOW time-averaged horizontal velocity profile
	P150=0;			  // int number of data points to read from grid file
	P151=1;			  // int type of data
	P152=4;			  // int type of boundary condition for data
    P160=0;            // int NHFLOW long-run statistics (P 160 t0 dt_out): harmonic analysis and time averages
    P160_ts=0.0;
    P160_dt=0.0;
    P161=1.0;          // double NHFLOW long-run statistics: Rayleigh criterion of the constituents
    P166=0;			  // int print discharge to terminal
	P167=0;			  // int discharge gages in x-direction
    P168=0;			  // int discharge gages in x-direction
	P180=0;			  // int print fsf
	P181=-10;		  // int ith iteration vtp (fsf, bed) printed
	P182=-1.0;       // double time between vtp (fsf, bed) file printout in seconds
    P183=-1.0;       // double time between vtp (fsf, bed) file printout in seconds for sediment time
    P184=0;       // int time between file printout in iterations
	P185=0;        	// int time between file printout in seconds
    P186=0;       // int time-averaged fsf vtp printouts
    P190=0;			  // int print topo
	P191=-10;		  // int ith iteration topo printed
	P192=-1.0;       // double time between topo file printout in seconds
    P194=0;       // int time between file printout in iterations
	P195=0;        	// int time between file printout in seconds
    P230=0;         // int print flowfile
    P240=0;         // int print potentialfile
    P250=0;             // int print CPM variables to vtu
    P310=0;             // int print breaking to fnpf vtu
    P311=0;             // int print breaking viscosity to fnpf vtu
	P351=0;             // int print out wsf lsm1
	P352=0;             // int print out wsf lsm2
    
    // DEM
    E10=0;              // int DEM on/off (input in dem.txt)
    E11=1;              // int fluid coupling: 0 dry, 1 unresolved, 2 resolved, 3 hybrid
    E12=5.0;            // double hybrid: resolved if d_eq/dx >= E12
    E13=1000;           // int max contact solver iterations
    E14=1.0e-5;         // double contact solver tolerance (relative velocity residual)
    E15=0.1;            // double print interval [s] for VTP and state file, 0 off
    E16=1;              // int min DEM substeps per fluid time step
    E17=0.5;            // double added mass coefficient (unresolved)
    E18=1.5;            // double resolved: Heaviside half width in cells
    E19=0.2;            // double penetration stabilisation factor
    E20=0.25;           // double max particle travel per DEM step, fraction of the smallest bounding radius
    E21=6;              // int contact manifold reduction: max points per pair and normal cluster, 0 off
    E22=1;              // int unresolved: fluid acceleration force (pressure gradient and added mass of the fluid acceleration), 0 off 1 on
    E23=2.0;            // double unresolved: kernel radius for the momentum exchange in cells (at least d_eq)
    E24=0.25;           // double distributed particles: max bounding radius as fraction of the smallest subdomain, larger ones are replicated; 0 replicates all
    E25=2.0;            // double unresolved: radius of the solid volume fraction kernel in equivalent diameters
    E26=0;              // int NHFLOW unresolved: 0 exchange at the centroid, 1 at the volume quadrature points (Peskin kernel)
    E27=0;              // int NHFLOW unresolved: fluid-side drag point-implicit (0 off, 1 on)
    E28=0;              // int NHFLOW: particle volume as porosity (0 off, 1 moving particles, 2 moving and fixed particles)
    E29=0;              // int NHFLOW unresolved: fixed particles exert drag on the fluid (0 off, 1 on)
    E30=1.0;            // double NHFLOW: relaxation time of the particle porosity [s] (E 28), 0 instantaneous

    // Particles
    Q10=0;              // int particle sediment model (S 10>0): 0 Eulerian bedload/exner, 1 CPM explicit Euler, 2 CPM RK2
    Q11=1;              // int particle movement model: 1 plain (empirical threshold), 2 MP-PIC
    Q12=2;              // int inter-particle stress model: 0 off, 1 Snider (2001), 2 packed bed (Johnson-Jackson contact pressure + mu(I) friction)
    Q13=1;              // int frictional shear stress for Q 12 2: 0 off, 1 on
    Q14=10;             // double continuum particle stress model parameter following Snider (2001) [Pressure]
    Q15=3.5;            // double continuum particle stress model parameter following Snider (2001)
    Q16=10e-7;          // double continuum particle stress model parameter following Snider (2001)
    Q17=0.6;            // double critical solid volume fraction
    Q18=0.5;            // double CFL coefficient for the particle sub-step based on the particle stress wave speed
    Q19=1;              // int grid-limited step: parcels enter a cell only if it has free volume (0 off, 1 on)
    Q26=0.5;            // double bed interface: solid fraction threshold relative to 1-S 24
    Q27=1;              // int number of smoothing passes for the solid volume fraction
    Q28=100;            // int maximum number of particle sub-steps per fluid time step
    Q32=-1.0;           // double packed bed: close packing, maximum solid fraction (<0: 1-S 24 + 0.035)
    Q33=0.05;           // double packed bed: contact pressure coefficient Fr [Pa] (Johnson & Jackson 1987), 0 off
    Q34=3.0;            // double packed bed: contact pressure exponent eta0
    Q35=5.0;            // double packed bed: contact pressure exponent eta1
    Q36=0.63;           // double packed bed: static friction coefficient mu_s
    Q37=1.13;           // double packed bed: dynamic friction coefficient mu_2
    Q38=0.6;            // double packed bed: inertial number scale I0
    Q24=0;              // int particles per cell
    Q25=1.25;           // double safety factor particle allocate
    Q29=1;              // int seed type
    Q30=0.58;           // double friction factor
    Q41=0;              // double default parcel factor
    Q45=-1.0;           // double solid fraction of the seeded region, default 1-S 24
    Q44=0;              // int all parcels fixed (e.g. flow through a fixed bed), 0 off, 1 on
    Q50=0;              // int fluid coupling of the parcels: 0 one-way, 1 two-way (drag reaction on the fluid, mixture continuity in the pressure equation)
    Q51=1;              // int drag law: 1 Andrews & O'Rourke (1996), 2 Gidaspow (Ergun below fluid fraction 0.8, Wen & Yu above)
    Q58=0;              // int sub-grid bedload layer (S 10 1): 0 off (resolved near-bed drag), 1 bedload layer, 2 bedload layer and pickup into the suspension, 3 bedload layer and Eulerian suspended load (hybrid), 4 as 3 without the onset criterion of the suspension (as the Eulerian model)
    Q59=6.5;            // double bedload layer: grain velocity coefficient, u_b = Q59 sqrt(R g d)(sqrt(theta) - 0.7 sqrt(theta_c))
    Q60=0.047;          // double bedload layer: critical Shields number theta_c of d50
    Q62=50.0;           // double bedload layer: mean hop length in grain diameters
    Q63=10.0;           // double bedload layer: relaxation time [s] of the bed level seen by the fluid, in bed time (Q 65 t)
    Q64=0.6;            // double bedload layer: thickness of the mobile surface layer of the bed in cells, the bed below is jammed
    Q65=1.0;            // double bedload layer: morphological factor, the pickup into the layer and the hops are Q 65 times faster
    Q66=0.0;            // double bedload layer: reference height of the bed shear stress in cells above the bed, 0: first cell
    Q67=0;              // int bedload layer: jammed bed only within the Mohr-Coulomb yield of the column (deep slope failures), 0 off, 1 on
    Q68=1;              // int bedload layer: bed shear stress, 0 log law, 1 max(log law, 0.3 rho k) with the turbulent kinetic energy, 2 sqrt(log law * 0.3 rho k) (two velocity scales)
    Q69=0;              // int bedload layer: seepage flow in the bed (Darcy), pore pressure and its force on the contact network, 0 off, 1 on
    Q70=0.0;            // double seepage (Q 69): excess pore pressure [Pa] at the bottom of the domain (upflow tests), 0: closed bottom
    Q72=2;              // int bedload layer: passes of a 1-2-1 filter over the columns of the bed level seen by the fluid, 0: off
    Q74=0;              // int bedload layer: 1 grains of the bed not released into the suspension see no resolved flow and no dispersion
    Q75=0.0;            // double bedload layer: sand slide of the bed surface beyond this slope (tan of the angle of repose), 0: off
    Q76=-1.0e20;        // double fixed (non-erodible) parcel bed upstream of x = Q76
    Q77=1.0e20;         // double fixed (non-erodible) parcel bed downstream of x = Q77
    Q79=1;              // int bedload layer: distance of the log law of the bed shear stress, 1 from the face under the first fluid cell (the wall of the k-omega wall law of the fluid at the bed), 0 from the bed level set (before patch 0039)
    Q78=0;              // int parcel merge and split: about Q78 parcels per cell; resting grains of the jammed bed (Q 58, S 10 1) merge in pairs in cells that hold more, and split again when they move or reach smaller cells; 0 off
    Q57=1;              // int Bagnold sheltering of the bed by the moving grains (S 10 1): 0 off, 1 on
    Q52=0;              // int turbulent dispersion of the parcels: 0 off, 1 random displacement with the eddy diffusivity
    Q56=0.0;            // double ride-over: blocked grains turn Q56 |u_blocked| upward (tan of the pivot angle), 0: off
    Q55=0.0;            // double dilatancy: theta_0(I) = theta_0 - Q55 I in sheared layers (Q 12 2), 0: off
    Q54=0.0;            // double lift coefficient C_L on the exposed grains (S 10 1 near-bed closure), 0: no lift
    Q53=1.0;            // double turbulent Schmidt number for Q 52
    Q61=0;              // point source x,y,z, iterations between particles
    Q73=0;              // int relax 
    Q110=0;             // int ini particle as box
    Q111=0;             // int ini particle as box for dummy particles
    Q180=0;             // int print vtp
    Q181=-10;           // int print vtp iter interval
    Q182=-1.0;          // double print vtp time interval

	// Sediment Transport
	S10=0;                  // int sediment transport module
	S11=0;                  // int bedload formula
	S12=0;                  // in Suspended Sediment, formula for boundary condition
	S13=10.0;               // double timestep for sediment transport
	S14=0.3;               // double relaxation timestep size for sediment transport
	S15=0;                  // int sediment time step: 0 S17*min(S13,CFL) 1 min(S17*dt,CFL) 2 S17*S13 3 S17*dt 4 S17*flow time since the last bed update
	S16=1;                  // int bed shear stress formulation
    S17=1.0;                  // double decoupling factor dtsed
    S18=0.823;                // double mu_d
	S19=1.0e+19; 			// double total time sediment
	S20=0.001;          // double sediment d50
	S21=3.0;          // double factor for d50 for calculation of ks in bedshear routine
    S22=2650.0;        // double sediment density
    S23=0.001;     // double sediment fall velocity
    S24=0.5;               // double porosity of sediment layer
    S25=1;     // int automatic sediment fall velocity
    S26_a=650.0;          // double alpha for VRANS sediment
    S26_b=2.2;            // double beta for VRANS sediment
    S27=1;              // int number of inner iterations
    S28=1;              // int use S21 in roughness BC
    S29=0;              // int sediment time step ramp up
    S30=0.047;          // double Shields parameter
    S31=2;              // int Exner formulation: 1 face flux qb*s_face, 2 upwinded cell flux vector qb*s (both conservative), 3 advective s*grad(qb) (not conservative)
    S32=1;              // int Exner discretisation: 1 FOU, 2 CDS, 4 WENO flux (conservative); 3 CDS-HJ (not conservative with S 31 1), 5 WENO-HJ (not conservative)
    S33=0;              // int non-equillibrium bedload 
    S34=1;              // int type of suspedned load D and E calculation
    S35=1.0;            // double Exner magnification factor
    S36=0;              // int alluvial/bedform/effective roughness
	S37=0;		        // int inertia bedshear stress
    S38=0;              // int shallow water bedshear calculation
    S39=100.0;            // double blend factor
    S40=0.1;            // double bedload adaptation length [m] (S33=1) or factor on water depth [-] (S33=4)
	S41=1;				// int type of sediment start criterion
	S42=1;				// int type of sediment interval criterion
    S43=1000;          // int number of water iteration, before sediment transport starts
    S44=10;            // int number of water timesteps between bed calculation
	S45=1.0;			// double flow simulation time, before sediment transport starts
	S46=1.0;			// double flow simulation time between bed calculation
	S47=1.0;			// double t/T, before sediment transport starts
	S48=0.1;			// int nt/T between bed calculation
    S49=200;            // int max number of non-equillibrium Gauss-Seidel sweeps (stops at convergence)
    S51=0;                  // int number of sediment fractions (count of S 51 lines), 0: single fraction d50=S20
    S52=0.0;                // double active layer thickness, <=0: 2*d90 of initial mixture
    S53=10.0;               // double initial substrate layer thickness
    S54=1;                  // int hiding/exposure function, 0: off, 1: Wu, Wang & Jia (2000)
    S55=0.6;                // double hiding/exposure exponent m
    S56=0;                  // int grain size for bed roughness ks=S21*d, 0: S20, 1: local d50, 2: local dm, 3: local d90
    S50=4;                  // int bc phi, 1: inflow fix or 2: outflow fix, 3: both fix
	S57=-1.0e20;        // double ini z-dir
    S58=0.0;            // double wave period for the wave friction bed shear stress S 16 9, <=0: period of the wave generation
    S59=0;              // int S 16 9 combination: 0 max(tau_log, tau_w), 1 max(tau_log, Soulsby DATA13), 2 Soulsby DATA13 (current from the running mean)
    S59_phi=0.0;        // double S 16 9 phase lead of the wave friction stress in deg
    S61=1;              // int reference height for cb
    S62=1;              // int type of suspended Exner solution
    S63=0.0;            // double S 16 10: reference height of the log law in cells above the bed, 0: first cell centre at least half a cell above the bed
    S71=-1.0e20;                 // int x start of erosion
    S72=1.0e20;          // int x end of erosion
    S73=0;       // double distance for use relaxation method for the sediment bed
    S74=0;           // int no erosion zone
    S76=0;           // int fixed-bed zones (S 76 xs xe ys ye): the bed held, the cells carry their bed load (supply of a reach)
    S75=0;       // double distance for use relaxation method for the sediment bed
    S77=0;          // int active sediment domain in x-direction
    S77_xs=-1.0e20; // double active sediment domain x_start
    S77_xe= 1.0e20; // double active sediment domain x_end
    S78=1;          // int inflow guard
    S80=0;                  // int type of slope reduction
    S81=35.0;              // double midphi for slope reduction
    S82=5.0;              // double delta phi for slope reduction
    S84=1;                  // int type of critical bed shear stress reduction limiters
    S85=1;                  // int bedload direction Koch and Flokstra
    S86=0;                  // int transverse bed slope effect on the bed-load direction: 1 f(theta) = S87 theta^S88 (d/h)^S89
    S87=9.0;                // double slope factor A_sh of f(theta) (Talmon et al. 1995: 9)
    S88=0.5;                // double Shields exponent B_sh of f(theta) (0.5)
    S89=0.3;                // double d/h exponent C_sh of f(theta) (0.3)
    S95=0;                  // int spiral-flow correction of the bed-load direction (depth-averaged SFLOW): 1 tan(delta) = alpha_I h/R_s
    S96=1.0;                // double spiral-flow intensity factor E_spir of alpha_I = E_spir 2/kappa^2 (1 - sqrt(g)/(kappa C))
    S97=1.0;                // double S 95 2: factor of the adaptation length of the spiral flow L_a = S97 (m-1) h/(2 kappa^2), m = kappa C/sqrt(g)
    S64=0.0;                // double SFLOW: suspended part f_s of the transport formula, relaxed over L_s = h|U|/(S65 ws) (0: off)
    S65=1.0;                // double S 64: ratio of the near-bed to the depth-averaged concentration in L_s
    S90=0;                  // int sandslide on/off
    S91=1;                  // int number of sandslide iterations
    S92=1.0;                // double sandslide correction factor
    S93=0.0;				// double delta phi for sandlide correciton
    S94=1;                // int sandslide intervall
	S100=3;					// int number of bed filter outer iterations (was hard-coded 3 in topovel2)
    S101=5;					// int number of bed filter inner iterations (was hard-coded 5 in topovel2)
    S102=0.0;              // double filter: fraction of the centre value kept per pass (0: 2dx wave removed)
    S103=0.0;              // double bed celerity CFL number for dtsed, 0: off
    S104=0.0;              // double limit of bed celerity x morphological factor / near-bed flow speed for dtsed, 0: off

    // Turbulence
	T10=0;			    // int turbulence model
	T12=5;              // int convection scheme
    T21=0;              // int type of LES filter
    T23=1.0;            // double URANS C3 coefficient
	T31=1.0/sqrt(3.0);  // double factor for limiter for eddy limiter
    T33=0;               // int kin source
    T34=1;               // int eddyv limiter
	T36=0;				// int explciti free surface dampong through dissipation
	T37=0.07;		    // double free-surface length scale y' for T36 / NHFLOW A567
    T38=1.6;            // double epsi fsf turbulence damping
    T39=0;              // blend fsf eddyv with sgs-eddyv
    T41=0;              // int RANS stabilization
    T42=0.05;           // double lambda2 factor (Larsen & Fuhrman stabilization)
    T44=1.0;            // double komega_func wall BC velocity factor for DF
    T45=0;              // int buoyancy term

    // Water Properties
	W1=998.2;		// double density water
	W2=1.004e-6;	// double viscosity water
	W3=1.205;		// double density air
	W4=1.41e-5;		// double viscosity air
	W5=0.0;			// double surface tension between phase 1 and phase 2
    W6=840.0;			// double density oil
	W7=3.0e-4;		// double viscosity oil
	W10=0.0;		// double discharge
    W11=0;         // int velocity inlet face 1
    W11_u=0.0;     // double u-velocity inlet face 1
    W11_v=0.0;     // double v-velocity inlet face 1
    W11_w=0.0;     // double w-velocity inlet face 1
    W12=0;         // int velocity inlet face 2
    W12_u=0.0;     // double u-velocity inlet face 2
    W12_v=0.0;     // double v-velocity inlet face 2
    W12_w=0.0;     // double w-velocity inlet face 2
    W13=0;         // int velocity inlet face 3
    W13_u=0.0;     // double u-velocity inlet face 3
    W13_v=0.0;     // double v-velocity inlet face 3
    W13_w=0.0;     // double w-velocity inlet face 3
    W14=0;         // int velocity inlet face 4
    W14_u=0.0;     // double u-velocity inlet face 4
    W14_v=0.0;     // double v-velocity inlet face 4
    W14_w=0.0;     // double w-velocity inlet face 4
    W15=0;         // int velocity inlet face 5
    W15_u=0.0;     // double u-velocity inlet face 5
    W15_v=0.0;     // double v-velocity inlet face 5
    W15_w=0.0;     // double w-velocity inlet face 5
    W16=0;         // int velocity inlet face 6
    W16_u=0.0;     // double u-velocity inlet face 6
    W16_v=0.0;     // double v-velocity inlet face 6
    W16_w=0.0;     // double w-velocity inlet face 6
	W20=0.0;		// double gi
	W21=0.0;		// double gj
	W22=0.0;		// double gk
    W29_x=0.0;		// double pressure gradient x-direction
    W29_y=0.0;		// double pressure gradient y-direction
    W29_z=0.0;		// double pressure gradient z-direction
	W30=0;			// int air compressibility on/off
	W31=20.0;		// double temperature for air compressibility in celsius
	W32=0;			// int isentropic air compressibility
	W33=1.4;		// double ratio of specific heats
	W34=-1.0;		// double interface band width
	W35=-1.0;		// double compensation layer width
	W36=101325.0;	// double reference pressure
    W41=0;         // int velocity source phase 1
    W50=0.0;        // double air inflow
    W50_air=0;      // int air inflow switch
    W90=0;           // int non-newtownian flow
    W95=0.001;       // double nu_0 maximum kinetic viscosity
    W96=1.0;         // double tau_0
    W97=0.00001;     // double K consitency parameter (viscosity)
    W98=1.0;         // int Herschel-Bulkley exponent (shear thinning/thickening)
    W101=0;          // int turn on Mohr-Coloumb
    W102_phi=30.0;   // double angle of repose
    W102_c=0.0;      // double coheasion
    W103=1.0;        // double MC transition factor
    W104=1.0;        // double shear rate dependent excess pore pressure factor
    W105_tau_00 = 0.0; // double grid-size dependency of the shear rate
    W105_C_total = 0.0; // double vol/vol concentration solids in total
    W105_C_kaolinite_chlorite = 0.0; // double vol/vol concentration of kaolinite_chlorite
    W105_C_illite = 0.0; // double vol/vol concentration of illite
    W105_C_montmorillonite = 0.0; // double vol/vol concentration of montmorillonite
    W106_b = 0.3;    // double consitency parameter k relation to tau_y
    W106_m_y = 0.2;  // double numerical parameter in Coulomb model [s]
    W107_mu_min = 0.0; // double minimal dynamic viscosity
    W107_delta = 0.0; // double internal friction angle
    W107_mu_0 = 0.0; // double maximal dynamic viscosity
    W108_a_2 = 0.3;  // double slurry phase volume contribution (only slurry and gravel)
    W110=1;          // int add rheology as source term or viscosity
    W111=1;          // int which pressure for MC
    W112=2.1;        // double threshold factor for pressure blening in W111 3
    W200=0;          // int atmospheric forcing (SFLOW, NHFLOW): 0 off, 1 uniform (W 201, W 202), 2 met file met-forcing.r3m, 3 vortex, 4 blend
    W201_u=0.0;      // double uniform wind speed U10 [m/s]
    W201_dir=0.0;    // double uniform wind direction [deg], counter-clockwise from the x axis, the direction the wind blows to
    W202=101325.0;   // double uniform air pressure [Pa]
    W203=0.0;        // double model start date-time YYYYMMDD.HHMMSS of the met forcing, 0: the first record of the file
    W204=0.0;        // double ramp time of the met forcing [s]
    W205=0;          // int geographic frame given (W 205 lon0 lat0 rot)
    W205_lon=W205_lat=W205_rot=0.0;  // double frame: transverse Mercator about lon0, lat0 [deg], x axis rot [deg] from east
    W206_x=W206_y=0.0;  // double model coordinates of (lon0, lat0) [m]
    W207=0;          // int drag law of the met forcing: 0 auto (1 for W 200 1/2, 5 for the vortex W 200 3), 1 Garratt, 2 Wu, 3 Smith and Banke, 4 Zijlema et al., 5 Garratt capped at W 213, 6 Zijlema et al. with the floor W 214
    W208=1;          // int air pressure gradient of the met forcing: 0 off, 1 on
    W209=101325.0;   // double reference air pressure of the inverse barometer [Pa]
    W210=1;          // int inverse barometer on the open edges and relaxation zones (B 520): 0 off, 1 on
    W211=0;          // int Coriolis force (SFLOW, NHFLOW): 0 off, 1 f-plane, 2 f per cell from the frame W 205
    W212=0.0;        // double f-plane Coriolis parameter [1/s], 0: 2 Omega sin(lat0) of W 205
    W213=2.5e-3;     // double cap of the drag coefficient of W 207 5 (saturation, Donelan et al. 2004)
    W214=1.5e-3;     // double floor of the drag coefficient of W 207 6 beyond the peak of Zijlema et al.
    W215=2;          // int vortex model of W 200 3: 1 Holland (1980), 2 GAHM, 3 Holland et al. (2010)
    W216=0.9;        // double vortex: reduction of the gradient wind to 10 m
    W217=0.8928;     // double vortex: 1-min to 10-min wind averaging factor
    W218=1;          // int vortex: inflow angle, 0 none, 1 ADCIRC (10 deg inside Rmax to 25 deg beyond 1.2 Rmax)
    W219=1;          // int vortex: translation, 0 none, 1 scaled by Vg/Vmax (ADCIRC), 2 full
    W220_n=0;        // int blend W 200 4: nest layer met-nest.r3m, 0 off, 1 on
    W220_v=1;        // int blend W 200 4: vortex layer met-vortex.trk, 0 off, 1 on
    W221_w=50000.0;  // double blend: taper width at the edges of the nest [m]
    W221_t=3600.0;   // double blend: taper after the first and before the last nest record [s]
    W222=1;          // int blend: basis of the vortex radii, 1 rmw of the track, 2 mean 34 kt radius, 3 metres
    W222_r1=3.0;     // double blend: pure vortex inside r1
    W222_r2=6.0;     // double blend: pure background beyond r2
    W223=1;          // int blend: vortex pressure, 0 as fitted, 1 deficit rescaled to the background on the ring r2
    W224=21600.0;    // double blend: taper of the vortex after the first and before the last track record [s]
    W225=1;          // int blend: taper, 1 cosine, 2 linear (ADCIRC NWS 29/30)
    W226=0;          // int blend: environment, 1 the background minus the vortex on the ring r2 (mean wind, pressure gradient) added to the vortex
    W227=0;          // int blend: nest in time, 0 linear at fixed points, 1 storm-relative (records shifted with the track centre)
    W_fb=0.0;        // double density of floating body

	// 6DOF
	X10=0;		// int turn 6DOF on
	X11_u=X11_v=X11_w=X11_p=X11_q=X11_r=1;		// int turn on degrees of freedom
    X12=1;      // int turn force calculation on
    X13=0;      // int NHFLOW floating body (X 10 1/2): water level in the pierced columns from the surrounding free surface (harmonic extension, 0 off)
    X14=1;      // int tangential velocity cfd
    X15=1;      // int tangential velocity nhflow (0 no slip, 1 tangential kept outside the body, 2 normal projection outside, 3 as 1 with the forcing inside the body only and the momentum fluxes at body faces from the fluid side, for bodies moving at speed, 4 normal projection outside ramped over the shell of 1, rotation invariant, 5 as 1 with the sharp switch at the hull surface instead of the shell ramp, 6 as 1 with the faces between body and fluid cells as reflective walls moving with the body in the HLL fluxes (mirror state of the fluid side about the body velocity; Cartesian grid, on the curvilinear grid the fluxes of 3), for wave reflection and wave loads, 7 as 6 with the forcing of 3 (inside the body only), for bodies at speed in waves)
    X16=0;          // int porous floating body (NHFLOW)
    X16_n=0.5;      // double porosity of the floating body
    X16_d50=0.01;   // double characteristic grain / element size of the floating body
    X16_alpha=1000.0; // double linear (Darcy) resistance coefficient
    X16_beta=1.1;   // double quadratic (Forchheimer) resistance coefficient
    X28=1;        // int FNPF: the added mass of the floating bodies (unit-mode psi solves, with several bodies the cross terms) every X28 time steps; in between the last one is used (1: every step)
    X30=1.0;      // double FNPF: tolerance of the psi solves of the body loads (pressure psi_0, added-mass unit modes) as a factor on N 44 (1: as the phi solves)
    X35=0;        // int FNPF: 1: the psi solves of the body loads start from the linear extrapolation in time of the last two solutions of their kind (psi_0 of the same stage, a unit mode) instead of the last one
    X36=0;        // int FNPF with mesh refinement (G 1, the zones G 12, not G 7 1): the unit-mode psi solves of the added mass of a body on its own patches only (the patches that reach into its zone, hull box grown by G 12), psi 0 at their edge and on level 0, instead of on all grids: 1 with the calibration (at the first added mass and after a new layout of the patches, every body also on all grids once; the differences of the added mass and of the cross terms are added to the window results), 2 without; the cross added mass from the bodies inside the patches only
    X37=0;        // int FNPF X 36 1: the calibration also every X 37 added-mass refreshes (0: at the first and after a new layout of the patches only)
    X29_d=X29_a=0.0;  // double FNPF with X28 > 1: the added mass again as soon as a body has moved by more than X29_d [m] or turned by more than X29_a [deg] since it was computed (0: off)
    X27=0;        // int FNPF: a free-surface column with a body node within X27 sigma layers below its free-surface node belongs to the footprint of the body (0: off)
    X17=0.5;        // double FNPF: relaxation of eta in the free-surface columns next to the footprint of a body, once per stage (0: off)
    X18=0.0;        // double FNPF: time constant [s] of the running mean of chi at body-fixed points, removed from the psi_0 free-surface data (0: off)
    X47=0.0;        // double FNPF: weight of the 5-point sawtooth filter (Longuet-Higgins and Cokelet) of eta and Fifsf in the free-surface columns within X47_n columns of the footprint of a body, once per stage (0: off; 1: the filtered value)
    X47_n=2;        // int FNPF: width of the filtered band at the footprint in columns (X 47)
    X19=1;        // int print out interval 6DOF log files
    X20=1;        // int number of floating body
    X21=1;        // int presribe homogeneous density floating body
    X21_d=900.0;  // double presribe homogeneous density floating body
    X22=0;        // int prescribe mass floating body
    X22_m=0;    // double prescribe mass floating body
    X23=0;        // int prescribe center of gravity
    X24=0;        // int prescribe moments of inertia
    X25_Cp=X25_Cq=X25_Cr=0.0;    // double damping rotation
    X26_Cu=X26_Cv=X26_Cw=0.0;    // double damping translational
    X31=4;        // int
    X32=1;        // int 
    X33=1;        // int 
    X34=0;        // int 
    X38=0;       // int type of viscous force calculation CFD
    X39=0;       // int type of viscous force calculation NHFLOW
    X40=3;		// int type of force calculation
	X41=0.6;    // double eps for continuous forcing heaviside
	X42=0.0;    // double distance for pressure force evaluation
	X43=1.0;    // double distance for shear stress evaluation
	X44=0.0;    // double viscosity in body
    X45=0;      // int type of lsm convection disc at fb
    X46=0;      // int 
    X48=0;      // int solid forcing level set
    X49=0;
    X50=1;      // int type of print out format for 6DOF structure
    X60=1;      // int type of print of force calculation
    X100=0;        // int delta x,y,z
    X100_x=X100_y=X100_z=0.0;
    X101=0;        // int ini Euler angles
    X101_phi=X101_theta=X101_psi=0.0;
    X102=0;        // int ini linear velocity
    X103=0;        // int ini angular velocity
    X103_p=X103_q=X103_r=0.0;
    X110=0;        // int rectangular box floating body
    X120=0;        // int sphere floating bod
    X120_rad=X120_xc=X120_yc=X120_zc=0.0;
    X131=0;        // int cylinder floating bod
    X131_rad=X131_h=X131_xc=X131_yc=X131_zc=0.0;
    X132=0;        // int cylinder floating bod
    X132_rad=X132_h=X132_xc=X132_yc=X132_zc=0.0;
    X133=0;        // int cylinder floating bod
    X133_rad=X133_h=X133_xc=X133_yc=X133_zc=0.0;
    X153=0;        // int symmetric wedge
    X163=0;        // int wedge
    X164=0;        // int hexahedron
    X170=0;        // int piston wavemaker
    X170_xs=X170_xe=X170_ys=X170_ye=X170_zs=X170_ze=0.0;
    X171=0;        // int flap wavemaker
    X171_xs=X171_xe=X171_ys=X171_ye=X171_zs=X171_z1=X171_ze=0.0;
    X172=0;        // int double-hinged flap wavemaker
    X172_xs=X172_xe=X172_ys=X172_ye=X172_zs=X172_z1=X172_z2=X172_ze=0.0;
    X180=0;        // int read .stl file for floating body geometry
    X181=1.0;   // double scale .stl geometry
    X181=0;     // int scale .stl geometry on/off
    X181_x=X181_y=X181_z=1.0;  // double scaling of stl geometry
    X182=0;     // int translation on/off
    X182_dx=X182_dy=X182_dz=0.0;  // double translation of stl geometry
    X183=0;
    X183_x=X183_y=X183_z=X183_phi=X183_theta=X183_psi=0.0;
    X184=1;     // int remeshing (X 185 4): close gaps and fill holes, 0 off, 1 harmonic patch, 2 biharmonic patch
    X185=4;     // int stl refinement
    X186=1.0;   // double refinement factor: target triangle size in fluid cells
    X187=30.0;  // double remeshing (X 185 4): feature angle in degrees
    X189=10;    // int remeshing (X 185 4): number of iterations
    X190=2;     // int remeshing (X 185 4): cell size, 1 grid metric (anisotropic), 2 geometric mean of dx,dy,dz, 3 smallest spacing
    X188=1;     // int ray cast algorithm
    X205=1;     // type of ramp up function
    X206=0;     // int ramp up velocity
    X206_ts=0.0;   // double ramp start
    X206_ts=0.0;   // double ramp start
    X207=0;     // int ramp up draft
    X207_ts=0.0;   // double ramp start
    X207_ts=0.0;   // double ramp start
	X210=0;		// int give fixed linear velocity
    X210_u=0.0; // double fixed u vel
    X210_v=0.0; // double fixed v vel
    X210_w=0.0; // double fixed w vel
	X211=0;		// int give fixed angular velocity
    X240=0;     // int read 6DOF motion file
    X241=0.0;     // double delta t for motion file
    X242_x=0.0;     // double delta x for motion file
    X242_y=0.0;     // double delta x for motion file
    X242_z=0.0;     // double delta x for motion file
    X243=0.0;     // double delta CoG for motion file
    X310=0;     // int mooring model
    X311=0;     // int number of simple taut mooring lines
    X312=0;     // int number of springs
    X313=0;     // int initial rotation of mooring end points with 6DOF body
    X314=0;     // int breaking mooring lines due to tension
    X315=0;     // int breaking mooring lines due to time
    X320=0;     // int net model
    X321=0;     // int number of nets
    X323_m=X323_d=X323_l=0.0;   // double dynamic net sinker properties
    X325_dt=0.001;   // double dynamic net time step
    X324=0;     // int number of nets
    X330=0;     // int impermeable membranes, geometry in membrane.dat
    X350=0;     // int ship module, parameters in ship.dat
	X325_relX=X325_relY=X325_relZ=0.01; // double dynamic net relaxation factors
	X400=0;         // sflow external pressure term
    X401_p0=0.0;    // sflow external pressure term p0
    X401_cl=2.0;    // sflow external pressure term cl
    X401_cb=16.0;   // sflow external pressure term cb
    X401_a=16.0;    // sflow external pressure term a
    X410=0;         // int etaval for draft
    X500=0;         // int PTO: 0 off, 1 explicit, 2 linearly implicit coupling
    X501=0;         // int PTO joint attachment point given
    X501_ax=X501_ay=0.0; // double PTO joint axis
    X501_az=1.0;
    X501_xa=X501_ya=X501_za=0.0; // double PTO joint attachment point, default CoG
    X502=0;         // int PTO linear damper
    X502_B=0.0;     // double damping coefficient [N s/m]
    X503=0;         // int PTO spring-damper
    X503_K=X503_B=0.0; // double stiffness [N/m], damping [N s/m]
    X504=0;         // int PTO end-stops
    X504_qmin=-1.0e20; // double lower stroke limit [m]
    X504_qmax=1.0e20;  // double upper stroke limit [m]
    X504_K=X504_C=0.0; // double end-stop stiffness [N/m], damping [N s/m]
    X505=0;         // int PTO generator limits
    X505_Fmax=X505_Pmax=0.0; // double force [N] and power [W] limit, 0: none
    X506=0;         // int PTO tuned passive damping (pto_hydro.dat)
    X506_T=0.0;     // double period [s], 0: wave input period (B 93)
    X506_C=0.0;     // double hydrostatic stiffness along the joint [N/m]
    X507=0;         // int PTO latching
    X507_t=0.0;     // double latching time [s]
    X507_K=X507_C=0.0; // double latch stiffness [N/m], damping [N s/m]
    X508=0;         // int PTO declutching
    X508_t=0.0;     // double declutching time [s]
    X510=0;         // int PTO drive-train losses
    X510_B=X510_Fc=0.0; // double loss damping [N s/m], loss Coulomb friction [N]
    X510_veps=0.0;  // double regularisation velocity of the loss friction [m/s], <= 0: 1e-3

	// Developer
	Y1=0;   // int turn on/off experimental screen force model
    Y2=0;   // int turn external moments on/off
    Y3=0;
    Y4=0;
    Y40=3;
	Y60=1;  // int require
    Y71=0;  // int turn on/off solid gcparax
    Y72=0;  // int turn on/off solid gcparax
    Y73=0;  // int turn on/off solid gcparax
    Y74=0;  // int turn on/off solid gcparax

	// FSI
	Z10=0;		// int turn FSI on
    Z11=0;      // int number of FSI beams (counted in read_control)
    Z12_ckx=Z12_cky=Z12_ckz=Z12_cdx=Z12_cdy=Z12_cdz=0.0;   // double fsi beam structural damping coefficients
    Z20=0;          // int flexible rod trees (soft corals, vegetation), 1: unresolved coupling, input rodtree.dat
    Z21=0.0;        // double rod-tree VTP print interval [s], 0: off
    Z30=0;          // int FEM solid structures (elastic, plastic, concrete damage, collapse), 1: coupling with CFD or NHFLOW, input fem.dat
    Z31=0.0;        // double FEM VTU print interval [s], 0: off
    Z40=0;          // int REEF3D::CORAL structural solver, 1: run the CORAL model of coral.txt with the fluid solver (CFD, NHFLOW, FNPF)
    
    // values of keys that also set their switch (read only with it): defined for the broadcast of
    // the control parameters to the other ranks (controlfill_send)
    A251_val=0.0;
    G501=0;      // (set by control::parse after the broadcast)
    B110_zs=B110_ze=0.0;
    F68_x=F68_y=F68_z=0.0;
    P43_xs=P43_xe=P43_ys=P43_ye=0.0;
    P46_is=P46_ie=0;
    P47_ts=P47_te=0;
    P101_xm=P101_ym=P101_zs=P101_ze=P101_r1=P101_r2=0.0;
    S29_ts=S29_te=S29_dts=S29_dte=0.0;
    X23_x=X23_y=X23_z=0.0;
    X24_Ix=X24_Iy=X24_Iz=0.0;
    X141=X142=X143=0;
    X153_xs=X153_xe=X153_ys=X153_ye=X153_zs=X153_ze=0.0;
    X206_te=X207_te=0.0;
    X211_p=X211_q=X211_r=0.0;
}
