/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2020 Tobias Martin

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
#include"mooring_barDyn.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void mooring_barDyn::initialize(lexer *p, fdm *a, ghostcell *pgc)
{		
    //- Initialise line

    L = p->X311_l[line];                // Length 
    nl = p->X311_H[line];               // Number of elements    
    
    d_c = 0.0004; //p->X321_d[nNet];    // Diameter of twine
    rho_c = p->X311_rho_c[line];        // Density of material
    
    kappa = 0.01;                   // Elasticity constant
    C1_ = 1160.0;    
    C2_ = 37300.0;  
   

/*---------------------------------------*/

    // Initialise values
    al = L/nl;                          // Length of bars along height
    
    nf = nl;
    niK = nf - 1;                    // Number of inner knots
    nbK = 2;                       // Number of boundary knots
    nK = niK + nbK;                     // Total number of knots

    
    
    //- Initialise fields
    
    p->Darray(coupledField, nK, 4);		// fluid coupling matrix (velocity 1,2,3 + phi 4)
    p->Darray(coupledFieldn, nK, 4);
    
    p->Darray(l0, nf);			// initial bar length

    p->Iarray(Pi,nf);           // inner owner knots
    p->Iarray(Ni,nf);           // inner neighbour knots
    p->Darray(K,nK,3);          // knot coordinates

    p->Iarray(nfK, niK, 3);		          // Bars per Knot -> first entry is knot ID, then up to two bars
 
    x0_ = MatrixXd::Zero(nK,3); 
    x_ = MatrixXd::Zero(nK,3);
    xdot_ = MatrixXd::Zero(nK,3);
    xdotdot_ = MatrixXd::Zero(nK,3);

    top_xdot_ = MatrixXd::Zero(nbK,3);
    top_xdotdot_ = MatrixXd::Zero(nbK,3);

    mass_knot = VectorXd::Zero(nK);
    weight_knot = VectorXd::Zero(nK);
    sinker_knot = VectorXd::Zero(nK);
    added_mass = VectorXd::Zero(nK);
    forces_knot = MatrixXd::Zero(nK, 3);

    A_ = MatrixXd::Zero(nf,nf);   // Linear system matrix
    B_ = VectorXd::Zero(nf);      // Linear rhs 
    T_ = VectorXd::Zero(nf);      // Tension forces
    T_old = VectorXd::Zero(nf);
  

    // Get standard net coordinates in K ---> init with catenary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    int index = 0;
    for (int i = 0; i <= nl; i++)
    {
        for (int j = 0; j <= nl; j++)
        {
            K[index][0] = j;
            K[index][1] = i;
            K[index][2] = 0.0;

            index++;
        }
    }


    // Initialise inner owner and neighbour lists

    for (int index = 0; index < nf; index++)
    {
        Pi[index] = index;
        Ni[index] = index + 1;
    }


    // List of bars per knot nfK
    double tmp;
    bool bK;
    int indexK = 0;
    int indexbK = 0;
    
    for (int i = 0; i < niK; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            nfK[i][j] = -1;
        }
    }

    for (int i = 1; i < nK-1; i++)
    {
        nfK[indexK][0] = i;

        // Search adjointed bars
        index = 1;
        for (int j = 0; j < nf; j++)
        {
            if (Pi[j] == i || Ni[j] == i)
            {
                nfK[indexK][index] = j;
                index++;
            }
        }
        
        indexK++;
    }

   
    //- Initialise x and length of bars
    
    for (int j = 0; j < nK; j++)
    {
        double a = K[j][0];
        double b = K[j][1];
        double c = K[j][2];
       
        x_.row(j) << K[j][0], K[j][1], K[j][2];

        x_(j,0) += origin_x;
        x_(j,1) += origin_y;
        x_(j,2) += origin_z;		
    } 
    
    x0_ = x_;

    Vector3d fi;
    
    for (int i = 0; i < nf; i++)
    {
        fi(0) = K[Ni[i]][0] - K[Pi[i]][0];
        fi(1) = K[Ni[i]][1] - K[Pi[i]][1];
        fi(2) = K[Ni[i]][2] - K[Pi[i]][2];
        
        l0[i] = fi.norm();
    }


    // Mass lumping at each inner knot

    double l_solid = 0.0;

    index = 0;
    
    for (int i = 2; i < nK-1; i++)
    {
        // Calculate total length from bars adjoint to knot nfK[index]
        for (int j = 0; j < nf; j++)
        {
            if (Pi[j] == nfK[index][0] || Ni[j] == nfK[index][0])
            {
                l_solid += l0[j]/2.0;
            }
        }
        
        // Mass (in air)
        mass_knot(i) = rho_c*PI/4.0*d_c*d_c*l_solid; 
        
        // Weight (in water)
        weight_knot(i) = p->W1*PI/4.0*d_c*d_c*l_solid; 
        
        // Added mass assuming ca = 1.0
        added_mass(i) = p->W1*PI/4.0*d_c*d_c*1.0*l_solid; 
/* 
        // Add sinker to bottom row
        if (i >= nK - nd - 1)
        {
            mass_knot(i) += sinker_m;
            
            weight_knot(i) += p->W1*PI/4.0*sinker_d*sinker_d*sinker_l;

            added_mass(i) += p->W1*PI/4.0*sinker_d*sinker_d*sinker_l*1.0;
        }
*/     
        index++;
    }
    

    //- Initialise old variables   
    xnn_ = x_;    
    xn_ = x_;
    xdotnn_ = xdot_;    
    xdotn_ = xdot_;
    
  //  if (p->X325_dt == 0.0)
    {   
        dtnn_ = p->dt;
        dtn_ = p->dt;
        dt_ = p->dt;
    }
  /*  else
    {
        dtnn_ = p->X325_dt;
        dtn_ = p->X325_dt;
        dt_ = p->X325_dt;
    }*/
    t_net_n = 0.0;
    t_net = 0.0;


    //- Initialise printing

	if(p->mpirank==0 && p->P14==1)
	{
		char str[1000];
		sprintf(str,"./REEF3D_6DOF_Mooring/REEF3D_6DOF_mooring_force_%i.dat",line);
		eTout.open(str);
		eTout<<"time \t T"<<endl;
	}		
    print(p);
	printtime = 0.0;
	

	// Initialise communication 
	ini_parallel(p, a, pgc);
}


void mooring_barDyn::ini_parallel(lexer *p, fdm *a, ghostcell *pgc)
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
