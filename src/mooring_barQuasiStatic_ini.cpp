/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2024 Tobias Martin

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
#include"mooring_barQuasiStatic.h"
#include"lexer.h"
#include"ghostcell.h"

void mooring_barQuasiStatic::initialize(lexer *p, ghostcell *pgc)
{		
	sigma = p->X311_H[line];

	double rho_f =1000.0;

	rho_c = p->X311_rho_c[line];
	w = p->X311_w[line]*9.81*(rho_c - rho_f)/rho_c;
	L = p->X311_l[line];

	p->Darray(x,sigma + 2); 
	p->Darray(y,sigma + 2);
	p->Darray(z,sigma + 2); 
	p->Darray(T,sigma + 2);
	p->Darray(Fb,sigma + 2);
	
	vector<double> three(3, 0);
	
	l0.resize(sigma+1, 0);
	l.resize(sigma+1, 0);
	A.resize(sigma+1, l);
	f.resize(sigma+1, three);
	B.resize(sigma+1, three);
		
	R.resize(sigma+1, three);
	v.resize(sigma+2, three);
	e.resize(sigma+1, three);
	c_coeff.resize(sigma+1, three);

	e_q.resize(sigma+1, three);
	e_d.resize(sigma+1, three);
	e_l.resize(sigma+1, three);

	if(p->mpirank==0 && p->P14==1)
	{
		char str[1000];
		sprintf(str,"./REEF3D_CFD_6DOF/REEF3D_6DOF_mooring_force_%i.dat",line);
		eTout.open(str);
		eTout<<"time \t T"<<endl;
	}		
	printtime = 0.0;
	
	
	// Initial element lengths
	l[0] = 0.0;
	l0[0] = 0.0;
	for (int j=0; j<sigma; j++)
	{
		l[j+1] = p->X311_l[line]/sigma;
		l0[j+1] = l[j+1];
	}
	
	
	// Filling system matrix A

	for (int j = 1; j < sigma + 1; j++)
	{		
		A[j-1][j-1] = -w;
		A[j-1][j]	= +w;
	}
		
	for (int j = 0; j < sigma; j++)
	{
		A[sigma][j] = 0.5*(l[j] + l[j+1]);
	}
	A[sigma][sigma] = 0.5*l[sigma];


	// Initial angles for direction vectors
	
	dx = p->X311_xe[line] - p->X311_xs[line];
	dy = p->X311_ye[line] - p->X311_ys[line];
	dz = p->X311_ze[line] - p->X311_zs[line];
	
	double magDist = sqrt(dx*dx+dy*dy+dz*dz);
		
	for (int j = 0; j < sigma + 1; j++)
	{		
		f[j][0] = dx/magDist;	
		f[j][1] = dy/magDist;	
		f[j][2] = dz/magDist;			
	}
	
	// Initialise communication 
	ini_parallel(p, pgc);

    // Initialise catenary
	pcatenary = new mooring_Catenary(line);

    // Initialise breaking
    broken = false;
    curr_time = 0.0;
    breakTension = p->X314 > 0 ? p->X314_T[line]: 0.0;
    breakTime = p->X315 > 0 ? p->X315_t[line]: 0.0;
}


void mooring_barQuasiStatic::ini_parallel(lexer *p, ghostcell *pgc)
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
