/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs
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
--------------------------------------------------------------------*/

#include"lexer.h"

void lexer::ini_default()
{
    // Hydrodynamic Models
    A10=5;       // int turn on wave models
    A209=2;      // int interpolation sweeps for bed
	A210=3;		  // int time scheme for SFLOW velocities
	A211=4;		  // int convection scheme for SLOW velocities
	A212=0;		  // int diffusion treatment for SLOW velocities
    A214=1;      // int convection for vertical velocity
    A215=0;      // int conservative discretization
    A216=2;      // int convection velocity
    A217=2;      // int slip or no-slip boundary conditions
    A218=0;      // int turn on roughness
    A219=1;      // int additional courant number constraint
	A220=1;		  // int pressure scheme for SFLOW
	A221=1;		  // int include hydrostatic pressure for SFLOW
    A222=1;      // int turn on secondary pressure gradient
    A223=0.5;    // double blending factor pressure gradient
    A230=0;      // int turn on Boussinesq wave model
	A240=1;      // int FSF algorithm SFLOW
	A241=1;		  // int discretization of water level SFLOW
	A242=0;		  // int hydostatic pressure for shallow areas
    A243=1;      // int turn on wetting-drying
    A244_val=0.00005; // double absolute wetting criterion value
    A244=1;      // int absolute wetting criterion
    A245=0;      // int dx-based relative wetting citerion
    A245_val=0.001; // double dx-based relative wetting citerion value
    A246=1;      // int turn on breaking
    A247=0.6;    // double breaking parameter alpha
    A248=0;      // int turn on breaking persistence
    A249=0.3;    // double breaking persistence parameter beta
    A251=0;    // double fsf-slope in x-dir
    A260=0;      // int turbulence model
    A261=0.267;  // double length scale factor
    
    A300=1;      // int sigma or fixed grid
    A310=3;		  // int time scheme for FNPF velocities
    A311=5;		  // int convection scheme for FNPF velocities
    A312=2;      // int discretization for second-order gradient
    A313=3;      // int discretization for bed bc
    A320=1;		  // int order of Laplace equation
    A321=1;      // int boundary condition order for 4th-order Laplace equation
    A341=0.0;    // double coastline damping distance factor for dxm
    A342=0.0;    // double coastline damping absolute distance
    A343=1;      // int turn on wetting-drying
    A344=1;      // int absolute wetting criterion
    A344_val=0.00005; // double absolute wetting criterion value
    A345=0;      // int dx-based relative wetting citerion
    A345_val=0.001; // double dx-based relative wetting citerion value
    
    
    A350=0;      // int turn on breaking (which method)
    A351=0;      // int type of breaking detection (deep / shallow)
    
    A354=0.6;    // double breaking parameter alpha
    A355=1.25;   // double breaking parameter slope
    A356=0.3;    // double breaking persistence parameter beta
    
    A361=5;      // int breaking filter outer iter
    A362=2;      // int breaking filter inner iter
    A363=1;      // int breaking filter width
    A365=1.86;   // double viscosity breaking wave
    A368=0;      // int breaking waves in numerical beach
    A369=1.86;   // double viscosity relaxation breaking wave
    
    
    A410=1;      // int scheme eta       
    A440=1.6;    // double epsi for depth integration
    
    // Boundary Conditions
	B10=0;			// int wall laws velocities on/off
    B19=0;			// int direct forcing for solid walls
	B20=2;			// int slip or no-slip boundary condition for velocity
	B26=1;			// int boundary condition implementation level set method
	B28=0;			// int extend pressure along bed boundary
	B29=0.5;		// double gamma for gc image point
	B30=0;			// int inflow crossection via cbc
	B50=0.001;		// double global wall roughness ks
	B51=-1.0;		// double global wall roughness ks
	B52=-1.0;		// double global wall roughness ks
	B53=-1.0;		// double global wall roughness ks
	B54=-1.0;		// double global wall roughness ks
	B55=-1.0;		// double global wall roughness ks
	B56=-1.0;		// double global wall roughness ks
	B60=0;            // int ioflow discharge
	B61=2;            // int plain or logarithmic inflow profile
	B62=3;            // int plain or logarithmic rans ini
	B63=3;            // int zero, plain or logarithmic inflow profile for rans
	B64=1;            // int reduce eddy near inflow
	B65=1.0;         // double factor for log/plain velocity inflow profile
	B66_1=1.6;		// double factor 1 for inflow exting into air
	B66_2=3.6;		// double factor 2 for inflow exting into air
	B67=0;       // double distance for use relaxation method for fixed kin
	B68=0;       // double distance for use relaxation method for fixed omega
	B69=0;       // double distance for use relaxation method for fixed eddyv
	B70=0;       // double distance for use relaxation method for fixed water level
	B71=0;       // double distance for use relaxation method for fixed water level ini
    B74=3;		 // int crossection i-indice for the periodic boundary condition
	B75=1;		 // int type of outflow boundary conditions
	B76=0;             // int fixed pressure inlet/outlet for ioflow and iowave
	B77=1;           // int outflow pressure controlled or free stream
	B78=1;           // int fixed pressure outlet for ioflow and iowave
	B79=1.0;	 // double scale factor for pressure outflow
	B80=-1.0e9;	 // double height threshold for pressure at outlet
	B81=0;		 // int focussed wave parameter
    B81_3=0.0;    // double unidirectional focused wave y is condisered 0
    B82=1;      // int type of focus point and time calculation
	B83=0.0025;  // double wave steepness parameter for focused breaking waves
    B84=1;       // Peak enhance method
	B85=1;		 // int PM or JONSWAP spectrum for irregular waves
	B86=10;		 // int number of regular waves for irregular wave generation
	B87=0;		 // int give ws and we for irregular wave generation
	B88=3.3;	 // double gamma for JONSWAP spectrum
	B89=0;            // int wave generation optimization
	B90=0;            // int iowave
	B91=0;          // int wave parameter wL
	B91_1=0.0;      // double wave amplitude
	B91_2=0.0;      // double wave length
	B91_3=0.0;      // double wave phase change
	B92=0;            // int wave type
	B93=0;          // int wave parameter wT
	B93_1=0.0;      // double wave amplitude
	B93_2=0.0;      // double wave period
	B93_3=0.0;      // double wave phase change
	B96_1=0.0;      // double dist1 for wave relax
	B96_2=0.0;      // double dist2 for wave relax
	B96_3=0.0;      // double dist3 for wave relax
	B97=0.0;        // double wave direction in degree
	B98=0;          // int type of wave generation
	B99=0;			// int type of numerical beach
	B101=0;        // int ramp function wave geneartion
	B102=1.0;        // double factor ramp function wave generation
	B103=0.6;        // double factor1 for relaxation into air phase
	B104=3.0;        // double factor2 for relaxation into air phase
	B105_1=0.0;        // double wave generation direction and line origin in Cartesian coordianate system 
	B105_2=0.0;        // double wave generation line x origin
	B105_3=0.0;        // double wave generation line y origin 
	B106=0;			// int read wave generation origin
	B107=0;			// int read numerical beach origin
    B108=0;        // int read wave geenration  origin
	B109=0;			// int generate depth field
	B110=0;			// int read water depth
	B110_d=0.0;		// double water depth
    B111_zs=0.0;	// double flap start
    B111_ze=0.0;	// double flap end
    B112_zs=0.0;	// double flap start
    B112_z2=0.0;	// double flap2 end/flap2 start
	B112_ze=1.0;	// double flag end
    B115=0;         // int activate vertical velocity component for flap wavemaker theory
    B116=1;         // int x or beta input for flap wavemaker theories
    B117=0.0;		  // double starting time shift for timeseries input
	B118=6.0;        // double factor in relaxation function 1
	B119=3.5;        // double factor in relaxation function 2
    B120=-90.0;       // doubel delta t for wave generation
	B121=1;        // int air velocity on/off in wave relaxation zone
    B122=1.0;        // int air velocity on/off for active wave generation
    B123=0.0;       // double flap AWA hinge location
	B126=-1.0e20;  // doubel no paddle movement threshold for solitary waves
    B130=0;         // int directional spreading for irregular waves
    B131=0.0;       // double main direction for multidirectional irregular waves
    B132_s=-90.0;  // double start directional spreading
    B132_e= 90.0;  // double end directional spreading
    B133=1;         // int number of direction intervals for spreading function
    B134=-1.0;      // double shape parameter for spreading function
    B135=10.0;       // double peak value
    B136=1;         // int double summation method frequency vector 
	B140_1=0.0;    // doube x1 numerical beach
	B140_2=0.0;    // doube x2 numerical beach
	B140_3=0.0;    // doube beta numerical beach
    B160=5;        // int number of vertical layers for 2D wave generation
	B180=0;           // int gravity waves
	B181_1=0.0;     // double x-acceleration amplitude
	B181_2=0.0;      //double x-acceleration frequency
	B181_3=0.0;     // double wave phase change
	B182_1=0.0;     // double y-acceleration amplitude
	B182_2=0.0;      //double y-acceleration frequency
	B182_3=0.0;     // double wave phase change
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
	B210=0;			// int peridic boundary conditions
	B240=0;			// int porous media  
	B241=1;			// int porous media in x-direction
	B242=1;			// int porous media in y-direction
	B243=1;			// int porous media in z-direction
    B260=0.0;       // double C coefficient for VRANS
    B264=1.0e20;    // double KC number for VRANS
    B267=0.001;     // double d50 for VRANS
	B269=0;			 // int VRANS on/off
    B270=0;         // int VRANS porous media box
    B274=0;         // int VRANS porous media vertical cylinder
    B281=0;         // int VRANS porous media wedge in x-direction
    B291=0;         // int VRANS porous media plate in x-direction
    B308=1;         // int porosity effects on fluid acceleration for vegetation
    B309=2.0;       // double Cm for vegetation
    B310=0;         // int VRANS vegetation box
	
	// Concentration
	C1=10.0;		// double density concentration in water
	C2=0.0;			// double viscosity water + concentration	
	C3=10.0;		// double density concentration in air
	C4=0.0;			// double viscosity air + concentration
	C5=1.0;			// double Schmidt number
    C9=1;          // int only phase 1 concentration
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
	D20=1;			// int diffusion scheme
	D21=0;			// int print out implicit diffusion time and iterations
	D22=1;			// int turn on molecular diffusion in diffusion discretization for momentum equations
	D23=1;			// int turn on molecular diffusion in diffusion discretization for turbulence equations
    D24=1;          // int use RK timestep factor for implicit diffusion
	D29=1.0e-5;     // double stopping criteria implicit diffusion
	D30=1;			// int pressure scheme
	D31=1;			// int additional correction steps for the PISO solver
	D32=2;			// int density calculation poisson equation
    D38=0;          // int add hydrostatic pressure gradients to NSE 
	
    // Free Surface
	F10=2;			    // int free surface scheme
    F11=0;              // int mom-fsf coupled solver
	F19=1.0e-5;     	// double stopping criteria level set 
	F30=0;			    // int level set scheme
	F31=0;             // particle level set
	F32=64;			// number of particles per cell
	F33=0.5;		// factor for pls vec allocation
	F34=5000;		// printout iteration for pls
	F35=5;			    // int convection scheme for fsf
	F36=1;				// int RK3 scheme
	F39=0.5;			    // double reini constraint relaxation factor
	F40=0;			    // int reini scheme
	F41=1;          // int reini after ith iteration
	F42=-1.0;		// double maxlength
	F43=0.55;		// double factor for reini timestep
	F44=3;		        // int number reini time step
	F45=2.1;         // factor for calculation of epsi
	F46=0;            // int picard iteration for lsm or reini
	F47=10;            // int number of picard iterations
	F49=1;            // int no reinitialization for interface nodes
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
	F60=-1.0e20;  // double ini z-dir
	F61=-1.0e20;  // double inflow  ini
	F62=-1.0e20;  // double outflow  ini
	F63=-1.0e20;  // double xstart phi interpolate with outflow h
	F64=0;			// int fsf plane with angle on/off
	F64_xs=0.0;			// double xs
	F64_ys=0.0;			// double xs
	F64_zs=0.0;			// double xs
	F64_alpha=0.0;			// double alpha
	F70=0;             // int number of phi 1 ini boxes
	F71=0;             // int number of phi 2 ini boxes
	F72=0;             // int number of phi 1 ini regions
	F80=0;             // int time scheme VOF
	F84=1.0;             // double cgamma for vof compression
    F85=0;             // int convection scheme VOF
	F101=0;             // int air entrainment model
	F150=0;         // int benchmark
	F151=0;         // int benchmark inverse sign of level set
	F300=0;			 // int multiphase flow level set
	F305=5;			 // int multiphase flow lsm convection
	F310=3;			 // int multiphase flow reini
	F321=1.6;		 // double epsi12
	F322=1.6;		 // double epsi13
	F323=1.6;		 // double epsi23
	F350=0;			 // int multiphase flow fix level set inflow/outflow
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
    G2=0;          // int sigma grid
	G10=3;			// int xmargin inflow
	G11=3;			// int ymargin right
	G12=3;			// int zmargin bottom
	G20=3;			// int xmargin outflow
	G21=3;			// int ymargin left
	G22=3;			// int zmargin top
	G30=3;			// int extrapolated ghost cells
	G39=0;			// int read flag_solid
    G40=3;         // int reini scheme for topo
	G50=1;            // int look for geodat points in grid file
    G52=0;            // int number of smoothing iterations
	G60=0;            // int number of geotopo boxes
	G61=0;            // int number of geotopo wedges
	G81=0;            // int number of inflow boxes
	G95=0;            // int number of outflow boxes

	// Heat
	H1=1.4e-7;      // thermal diffusivity water
	H2=2.216e-5;      // thermal diffusivity air
	H10=0;          // int heat transfer on/off
	H50_1=20.0; // double temperature 1
	H50_2=20.0; // double temperature 2
	H51=0.0;		    // double i-dir zero level set start
	H52=0.0;		    // double j-dir zero level set start
	H53=0.0;		    // double k-dir zero level set start
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

    // Initialize
    I10=0;          // int initialize all
	I11=0;			// int initialize velocities with potential flow
	I12=0;          // int initialize pressure
	I13=0;          // int initialize turbulence
	I20=250;    // int number of potential flow iterations
	I21=0;          // int set phase 2 velocities to zero after potential flow solver
	I30=0;			// int Fully intialize NWT
	I40=0;			// int ini from state file
	I41=0;			// int ID of state file
    I55=0.0;        // double reference pressure
	I56=0;          // int pressure above F56 set to zero
	I58_1=0.0;      // double vertical velocity for sphere initialization
	I58_2=0.0;      // double radius for sphere initialization
    I230=0;         // int read 2D flowfile
    I231=0.0;       // double delta z for 2D/3D flowfile
    I232=0.0;       // double delta t for 2D/3D flowfile
    I240=0;         // int read 2D flowfile
    I241=0.0;       // double delta z for 2D/3D flowfile
    I242=0.0;       // double delta t for 2D/3D flowfile

    // Nunmerics
	N5=1;			// int true 2D
    N8=3;			// int linear solver velocities and scalars
	N9=1;         // int precondioner
	N10=14;			// int linear poisson solver
	N11=11;         // int precondioner
	N12=1;         // 
    N13=10;       // int inneriter 
    N14=5;          // int CGC gridnum
	N15=-1;			// int gcupdate during solver
	N16=1;			// int iteration limit for gcupdate gcupdate 
	N17=0.5;		// double relaxation factor
	N18=1.0e-10;	// double stopping criteria precon
	N40=1;			// int time scheme
	N41=1.0e+19; 	// double total time
	N42=3;			// int RK3 scheme
	N43=1.0e-6;     // double stopping criteria convection-diffusion
	N44=1.0e-8;     // double stopping criteria pressure
	N45=1e8;		// max outer iter
	N46=250;		// int max number of solver iterations
	N47=0.3;		// doubel relaxation factor for time stepping
	N48=1;          // int adaptive timestepping
	N49=1.0;		// double max timestep or fixed timesteps
	N50=1;			// int inner SIMPLE loop iterations
	N51=1.0e-5;     // double stopping criteria inner SIMPLE loop
	N52=1;          // int inner loop for turb and lsm
	N53=1.0;       // double pressure relaxation factor, implicit
	N54=1.0;       // double velocities relaxation factor, implicit
	N55=1.0;       // double turbulence relaxation factor, implicit
	N56=1.0;       // double levelset relaxation factor, implicit
	N57_1=0;       // int iteration start increasing N47
	N57_2=0;       // int iteration end increasing N47
	N57_3=0.3;     // double final N47
	N58=0;         // int enable timestep dependend relaxation factors
	N60=10;        // int maximum iteration of pjm correction
	N61=500.0;      // double stopping criteria velocities
	N71=2;          // order of gcepol extend
	N72=2;          // order of gcepol dirichlet
	N73=2;          // order of gcepol pressure

    // MPI
	M10=0;          // number of MPI processes

    // Printer
	P10=1;			 // int print file type
    P11=10;			 // int log print frequency
	P12=1;			 // int terminal print frequency
	P14=1;           // int print to folder
	P15=1;          // int print file numbering
	P17=0;			// int print paraview cutcell warning to console
	P18=2;			// int option for phi print out
	P19=1;			// int option for phi ccipol print out
	P20=-10;		// ith iteration file printed
	P23=0;			// int print test to vtu file
    P24=0;			// int print density to vtu file
    P25=0;			// int print solid to vtu file
	P26=0;			// int print cbed and conc to vtu file
	P27=0;			// int print topo to vtu file
	P28=0;			// int print fb to vtu file
	P29=0;			// int print walldist to vtu file
	P30=-1.0;       // double time between file printout in seconds
	P34=-1.0;       // double time between file printout in seconds for sediment
	P35=0;        	// int print for interval
	P40=0;				// int print state file
	P41=-10;				// int print state file each ith iteration
	P42=-1.0;				// double print state file each ith sec
	P50=0;				// int wave theory wave gages
	P51=0;             // int print out wsf
	P52=0;            // int print out wsf line in x-dir
	P53=0;            // int print out wsf line for wave theory
	P54=10;			  // int ith iteration wsf file  print out
	P55=-1.0;		  // double ith second wsf files print out
	P56=0;            // int print out wsf line in y-dir
	P59=0;			  // int print runup
	P61=0;			  // int print point probes
	P62=0;			  // int print line probes
	P66=0;			  // int print discharge to terminal
	P67=0;			  // int discharge gages in x-direction
    P71=0;           // int print viscosity to vtu
	P75=0;            // int print out vorticity vec
    P78=0;            // int print out velocity magnitude as scalar
	P79=0;            // int print out sediment shear stress when running sediment transport
	P81=0;            // int force print out
	P82_x=0.1;        // double reference length in x direction
	P82_y=0.1;        // double reference length in y direction
	P83=-1.6;		  // double epsi interface threshold for force calculation 
	P84=1.0;		  // double diameter for normalized vertical force calculation
	P85=0;            // int rectangular force print out
	P86_x=0.1;        // double reference length in x direction
	P86_y=0.1;        // double reference length in y direction
	P87=-1.6;		  // double epsi interface threshold for force calculation 
	P88=1.0;		  // double diameter for normalized vertical force calculation
	P89_cm=2.0;		  // double cm value 
	P89_cd=1.0;		  // double cd value 
	P90=0.0;		  // double obstacle truncation, distance between floor and obstacle
	P91=0.25;		  // double factor used in force calculation algorithm
    P92=0;           // int force from water or from water+air
    P93=0;             // int turn on Morison force calculation
	P101=0;			  // int print runup gages
	P121=0;             // int bed level gages
	P122=0;             // int max bed level gages
	P123=0;             // int topoline in x-direction
	P124=0;             // int topoline in y-direction
	P125=0;             // int bed shear stress gages
	P126=0;             // int bed shear stress maxval
	P150=0;			  // int number of data points to read from grid file
	P151=1;			  // int type of data
	P152=4;			  // int type of boundary condition for data
	P180=0;			  // int print fsf
	P181=10;		  // int ith iteration fsf printed
	P182=-1.0;       // double time between fsf file printout in seconds
    P191=0;         // print sloshing forces
    P210=0;			  // int print exportfile
	P211=1;		  // int ith iteration export printed
	P212=-1.0;       // double time between export file printout in seconds
    P230=0;         // int print flowfile
    P240=0;         // int print potentialfile
	P351=0;             // int print out wsf lsm1
	P352=0;             // int print out wsf lsm2
	
	// Sediment Transport
	S10=0;                  // int sediment transport module
	S11=0;                  // int bedload formula
	S12=0;                  // in Suspended Sediment, formula for boundary condition
	S13=10.0;               // double timestep for sediment transport
	S14=0.3;               // double relaxation timestep size for sediment transport
	S15=0;                  // int synchronize sediment time step with main solver
	S16=1;                  // int bed shear stress formulation
	S17=0;                  // int use inner topo iterations
	S18=100;                  // int maximum inner topo iterations
	S19=1.0e+19; 			// double total time sediment
	S20=0.001;          // double sediment d50
	S21=3.0;          // double factor for d50 for calculation of ks in bedshear routine
    S22=2650.0;        // double sediment density
    S23=0.000001;     // double sediment fall velocity
    S24=0.5;               // double porosity of sediment layer
    S25=30.0;             // double angle of respose
	S28=0.75;             // double factor for new fluid cell velocity interpolation
	S29=0.5;             // double factor for new solid cell velocity interpolation
    S30=0.047;          // double Shields parameter
    S31=1;		        // int exner scheme
	S37=2;		        // int number reini time step
	S38=0;            // int reinitialization for interface nodes
	S39=1;				// int order of time discretization for dz/dh
	S41=1;				// int type of sediment start criterion
	S42=1;				// int type of sediment interval criterion
    S43=1000;          // int number of water iteration, before sediment transport starts
    S44=10;            // int number of water timesteps between bed calculation
	S45=1.0;			// double flow simulation time, before sediment transport starts
	S46=1.0;			// double flow simulation time between bed calculation
	S47=1.0;			// double t/T, before sediment transport starts
	S48=0.1;			// int nt/T between bed calculation
    S50=4;                  // int bc phi, 1: inflow fix or 2: outflow fix, 3: both fix
	S57=-1.0e20;        // double ini z-dir
    S60=0;                  // int time stepping for suspended sediments
    S71=-1.0e20;                 // int x start of erosion
    S72=1.0e20;          // int x end of erosion
    S73=0;       // double distance for use relaxation method for the sediment bed
    S77=0;          // int active sediment domain in x-direction
    S77_xs=-1.0e20; // double active sediment domain x_start
    S77_xe= 1.0e20; // double active sediment domain x_end
    S80=0;                  // int type of slope reduction
    S81=35.0;              // double midphi for slope reduction
    S82=5.0;              // double delta phi for slope reduction
    S90=0;                  // int sandslide on/off
    S91=1;                  // int number of sandslide iterations
	S92=1.0;				// double safety factor sandslide correction
    S93=0.0;				// double delta phi for sandlide correciton
	S100=0;					// int number of bed filter outer iterations
    S101=0;					// int number of bed filter inner iterations
	S102=0;					// int number of bedload filter outer iterations
    S103=0;					// int number of bedload filter inner iterations
    S116=1.0;               // double factor for location of bed shear stress calculation
    S117=1.0;               // double factor for location of bed shear stress calculation

    // Turbulence
	T10=0;			    // int turbulence model
	T11=11;             // int time scheme for 2-eq turbulence models
	T12=5;              // int convection scheme
	T13=1.0e-5;     	// double stopping criteria turbulence
	T30=2;			    // int eddyv limiter
	T31=0.816;	        // double factor for limiter for eddy limiter in phase 1
	T32=0.816;	        // double factor for limiter for eddy limiter in phase 2
	T35=0.816; 			// double factor for limiter for eddy limiter near wall
	T36=0;				// int explciti free surface dampong through dissipation
	T37=0.07;		    // int damping coefficient for T36
    T38=1.6;            // double epsi fsf turbulence damping
    T39=1.0;            // double dirac function factor for fsf damping
	T40=1;              // int source term limiter
	T41=1;				// int source term definition for omega in kw model
	T42=10.0;			// double factor for limiter for production limiter
	T43=1;				// int Pk term with strain or vorticity    
	T51=0.0;            // double inflow turbulent intensity
	T52=0.00001;        // double inflow turbulent intensity
	T53=1.0;            // double factor for inflow turbulent dissipation

    // Water Properties
	W1=998.2;		// double density water
	W2=1.004e-6;	// double viscosity water
	W3=1.205;		// double density air
	W4=1.41e-5;		// double viscosity air
	W5=0.0;			// double surface tension between phase 1 and phase 2
	W6=840.0;			// double density oil
	W7=3.0e-4;		// double viscosity oil
	W10=0.0;		// double discharge
	W20=0.0;		// double gi
	W21=0.0;		// double gj
	W22=0.0;		// double gk
	W30=0;			// int air compressibility on/off
	W31=20.0;		// double temperature for air compressibility in celsius
    W90=0;           // int non-newtownian flow
    W95=0.001;       // double nu_0
	W96=1.0;         // double tau_0
    W97=0.00001;     // double K
    W98=1.0;         // int 
    W101=0;          // int turn on Mohr-Coloumb 
    W102_phi=30.0;   // double angle of repose
    W102_c=0.0;      // double c factor
    W103=1.0;        // double MC transition factor
    W104=1.0;        // double shear rate dependent excess pore pressure factor
    W110=1;          // int add rheology as source term or viscosity
    W111=1;           // int which pressure for MC
    W112=2.1;         // double threshold factor for pressure blening in W111 3
	
	// 6DOF
	X10=0;		// int turn 6DOF on 
	X11_u=X11_v=X11_w=X11_p=X11_q=X11_r=1;		// int turn on degrees of freedom
    X12=1;      // int turn force calculation on
	X13=0;      // int turn 6DOF algorithm wit quaternions on
	X18=0;		// int relaxation method solid velocities
	X19=1;		// int print out interval 6DOF log files
	X21=1;		// int presribe homogeneous density floating body
	X21_d=900.0;		// double presribe homogeneous density floating body
	X22=0;		// int prescribe mass floating body
	X22_m=0;	// double prescribe mass floating body
	X23=0;		// int prescribe center of gravity
	X24=0;		// int prescribe moments of inertia
	X25_Cp=X25_Cq=X25_Cr=0.0;	// double damping rotation
    X26_Ku=X26_Kv=X26_Kw=0.0;	// double spring coefficients
    X27=0;      // int spring origin given 
    X27_x=X27_y=X27_z=0.0;  // double spring origin
	X31=3;		// int boundary conditions for parallel velocity on floating body
	X32=1;		// int boundary conditions for orthogonal velocity on floating body
	X33=1;		// int boundary conditions for pressure on floating body
    X34=0;		// int boundary treatment for new solid velocity cells
	X38=0;		// int turn on direct forcing
    X40=1;		// int type of force calculation
	X41=1.75;		// double eps for lsm based force calculation
	X100=0;		// int delta x,y,z
	X100_x=X100_y=X100_z=0.0; 
	X101=0;		// int ini Euler angles
	X101_phi=X101_theta=X101_psi=0.0; 
	X102=0;		// int ini linear velocity
	X102_u=X102_v=X102_w=0.0; 
	X103=0;		// int ini angular velocity
	X103_p=X103_q=X103_r=0.0; 
	X110=0;		// int rectangular box floating body
	X120=0;		// int sphere floating bod
	X120_rad=X120_xc=X120_yc=X120_zc=0.0;
	X131=0;		// int cylinder floating bod
	X131_rad=X131_h=X131_xc=X131_yc=X131_zc=0.0;
	X132=0;		// int cylinder floating bod
	X132_rad=X132_h=X132_xc=X132_yc=X132_zc=0.0;
	X133=0;		// int cylinder floating bod
	X133_rad=X133_h=X133_xc=X133_yc=X133_zc=0.0;
    X153=0;		// int symmetric wedge
    X163=0;		// int wedge
    X164=0;		// int hexahedron
	X180=0;		// int read .stl file for floating body geometry
    X181=1.0;   // double scale .stl geometry
    X182=0;     // int translation on/off
    X182_x=X182_y=X182_z=0.0;  // double translation of stl geometry
    X183=0;
    X183_x=X183_y=X183_z=X183_phi=X183_theta=X183_psi=0.0;
	X210=0;		// int give fixed linear velocity
    X210_u=0.0; // double fixed u vel
    X210_v=0.0; // double fixed v vel
    X210_w=0.0; // double fixed w vel
	X211=0;		// int give fixed angular velocity
    X221=0;     // int read vec based motion file
    X311=0;     // int number of simple taut mooring lines
	X320=0;     // int number of nets
	
	// Developer 
	Y40=0;
    Y50=5;
	Y60=0;
    Y71=0;  // int turn on/off solid gcparax
    Y72=0;  // int turn on/off solid gcparax
    Y73=0;  // int turn on/off solid gcparax
    Y74=0;  // int turn on/off solid gcparax
    
	
	solveriter=0;
	mpirank=0;

	simtime=0.0;
	poissontime=0.0;
	pressval=0;
    alpha=0.0;
}
