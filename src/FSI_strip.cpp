/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"FSI_strip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

fsi_strip::fsi_strip(int num):nstrip(num),beam(num)
{}
    
fsi_strip::~fsi_strip(){}

void fsi_strip::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
	// Initialise parameter
	double gamma = p->X311_w[nstrip];			// specific weight [kg/m]
	double EA = p->X311_EA[nstrip];   	// stiffness [N]
    double L = p->X311_l[nstrip];      			// length of unstretched cable [m]
	double rho_c = p->X311_rho_c[nstrip];   		// density of material [kg/m3]
	double d_c = p->X311_d[nstrip];      		// diameter of the cable [m]
	double Ne = p->X311_H[nstrip];				// number of elements
 
    double A = gamma/rho_c;
    
    // Initialise beam
    iniBeam(Ne, EA/A, gamma/rho_c, rho_c, L, (EA/A)/(2.0*(1.0 + 0.5)), 1e-8, 1e-8, 1e-8);
    
    // Initialise material
    iniMaterial();

    // Initialise damping and compression effects
    iniDamping(10,0,0,0,0,0,true);

    // Meshing
    Eigen::VectorXd xIni = Eigen::VectorXd::Zero(Ne+1);   
    Eigen::VectorXd yIni = Eigen::VectorXd::Zero(Ne+1);   
    Eigen::VectorXd zIni = Eigen::VectorXd::Zero(Ne+1);   
	//mooring_Catenary *pcatenary;
	//pcatenary = new mooring_Catenary(line);
    //pcatenary->iniShape(p,a,pgc,xIni,yIni,zIni);
    meshBeam(xIni, yIni, zIni);

    // Initialise solver
    iniSolver();
	
	// Initialise communication 
//	ini_parallel(p, a, pgc);

    // Initialise print
	if(p->mpirank==0 && p->P14==1)
	{
		/*char str[1000];
		sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_mooring_force_%i.dat",line);
		eTout.open(str);
		eTout<<"time \t T"<<endl;*/	
	}

    // Initialise fields
   /* c_moor = Matrix3Xd::Zero(3,Ne+1); 
    c_moor_n = Matrix3Xd::Zero(3,Ne+1); 
    cdot_moor = Matrix3Xd::Zero(3,Ne+1); 
    cdot_moor_n = Matrix3Xd::Zero(3,Ne+1);
    cdotdot_moor = Matrix3Xd::Zero(3,Ne+1);

    getTransPos(c_moor); c_moor_n = c_moor;

    vector<double> three(3, 0);
	fluid_vel.resize(Ne+1, three);
	fluid_vel_n.resize(Ne+1, three);
	fluid_acc.resize(Ne+1, three);

    t_mooring = 0.0;
    t_mooring_n = 0.0;*/
}
	
void fsi_strip::start(lexer*,fdm*,ghostcell*){};
