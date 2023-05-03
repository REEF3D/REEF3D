/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2023 Tobias Martin

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

#include<sys/stat.h>
#include"mooring_dynamic.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void mooring_dynamic::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
	// Initialise parameter
	gamma = p->X311_w[line];			// specific weight [kg/m]
	double EA = p->X311_EA[line];   	// stiffness [N]
    L = p->X311_l[line];      			// length of unstretched cable [m]
	rho_c = p->X311_rho_c[line];   		// density of material [kg/m3]
	d_c = p->X311_d[line];      		// diameter of the cable [m]
	Ne = p->X311_H[line];				// number of elements
 
    A = gamma/rho_c;
    
    // Initialise beam
    iniBeam(Ne, EA/A, gamma/rho_c, rho_c, L, (EA/A)/(2.0*(1.0 + 0.5)), 1e-8, 1e-8, 1e-8);
    
    // Initialise material
    iniMaterial();

    // Initialise damping and compression effects
    iniDamping(10,0,0,0,0,0,false);

    // Meshing
    Eigen::VectorXd xIni = Eigen::VectorXd::Zero(Ne+1);   
    Eigen::VectorXd yIni = Eigen::VectorXd::Zero(Ne+1);   
    Eigen::VectorXd zIni = Eigen::VectorXd::Zero(Ne+1);   
	mooring_Catenary *pcatenary;
	pcatenary = new mooring_Catenary(line);
    pcatenary->iniShape(p,a,pgc,xIni,yIni,zIni);
    Eigen::Vector3d d0;  d0 << 1, 0, 0;
    meshBeam(xIni, yIni, zIni, d0);

    // Initialise solver
    iniSolver();
	
	// Initialise communication 
	ini_parallel(p, a, pgc);

    // Initialise print
	if(p->mpirank==0 && p->P14==1)
	{
		/*char str[1000];
		sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_mooring_force_%i.dat",line);
		eTout.open(str);
		eTout<<"time \t T"<<endl;*/	
	}

    // Initialise fields
    c_moor = Matrix3Xd::Zero(3,Ne+1); 
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
    t_mooring_n = 0.0;
    
    // Initialise breaking
    broken = false;
    breakTension = p->X314 > 0 ? p->X314_T[line]: 0.0;
    breakTime = p->X315 > 0 ? p->X315_t[line]: 0.0;
}

void mooring_dynamic::ini_parallel(lexer *p, fdm *a, ghostcell *pgc)
{
	p->Darray(xstart, p->mpi_size);
	p->Darray(xend, p->mpi_size);
	p->Darray(ystart, p->mpi_size);
	p->Darray(yend, p->mpi_size);
	p->Darray(zstart, p->mpi_size);
	p->Darray(zend, p->mpi_size);
	
	xstart[p->mpirank] = p->originx;
	ystart[p->mpirank] = p->originy;
	zstart[p->mpirank] = p->originz;
	xend[p->mpirank] = p->endx;
	yend[p->mpirank] = p->endy;
	zend[p->mpirank] = p->endz;
	
	for (int i = 0; i < p->mpi_size; i++)
	{
		MPI_Bcast(&xstart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
		MPI_Bcast(&xend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
		MPI_Bcast(&ystart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
		MPI_Bcast(&yend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
		MPI_Bcast(&zstart[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
		MPI_Bcast(&zend[i],1,MPI_DOUBLE,i,pgc->mpi_comm);
	}
}
