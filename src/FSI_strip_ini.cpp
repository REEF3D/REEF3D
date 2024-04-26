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

#include"FSI_strip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void fsi_strip::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"FSI initialize"<<endl;
    
	// Initialise parameter
    double x_ini = p->Z11_x[nstrip]; // x-position of strip bottom
    double y_ini = p->Z11_y[nstrip]; // y-position of strip bottom
    double z_ini = p->Z11_z[nstrip]; // z-position of strip bottom
    double L = p->Z11_l[nstrip];     // Length of strip [m]
    W_el = p->Z11_w[nstrip];         // Width of strip [m]
    double T = p->Z11_t[nstrip];     // Thickness of strip [m]
	rho_s = p->Z11_rho[nstrip];      // Density of material [kg/m3]
	double E = p->Z11_e[nstrip];   	 // Young modulus [N/m^2]
	double Ix = p->Z11_ix[nstrip];   // X-moment of area [m^4]
	double Iy = p->Z11_iy[nstrip];   // Y-moment of area [m^4]
	double Iz = p->Z11_iz[nstrip];   // Z-moment of area [m^4]
	double Nu = p->Z11_nu[nstrip];   // Poisson ratio [-]
	Ne = p->Z11_n[nstrip];           // Number of elements

    thinStrip = false;
    if (p->Y2==1) thinStrip = true;

    gravity_vec << a->gi, a->gj, a->gk;
    rho_f = p->W1;
    A_el = W_el*T;

    // Initialise beam
    iniBeam(Ne, E, A_el, rho_s, L, E/(2.0*(1.0 + Nu)), Ix, Iy, Iz);

    // Initialise material
    iniMaterial();

    // Initialise damping and compression effects
    iniDamping(p->Z12_cdx,p->Z12_cdy,p->Z12_cdz,p->Z12_ckx,p->Z12_cky,p->Z12_ckz,true);

    // Meshing
    Eigen::Matrix3Xd ini_coord = Eigen::Matrix3Xd::Zero(3,Ne+1); 
    for (int n = 0; n < Ne+1; n++)
    {
        ini_coord.col(n) << x_ini, y_ini, z_ini + L/Ne*n;
    }
    
    Eigen::Vector3d d1,d2,d3;  
    d1 << 0, 0, 1;
    d2 << 0, 1, 0;
    d3 << 1, 0, 0;
    meshBeam(x_ini,y_ini,z_ini,d1,d2,d3);

    // Initialise solver
    iniSolver();
	
	// Initialise communication 
	ini_parallel(p, a, pgc);

    // Initialise cell size
    get_cellsize(p, a, pgc);

    // Initialise field
    getTransPos(x_el);
    getTransVel(xdot_el);
    getRotPos(q_el);
    getRotVel(qdot_el);
    
    t_strip = 0.0;
    t_strip_n = 0.0;

    F_el = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    P_el = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    P_el_n = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    M_el = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    I_el = Eigen::Matrix3Xd::Zero(3,Ne+2);   
    I_el_n = Eigen::Matrix3Xd::Zero(3,Ne+2);   

    // Initialise Lagrangian fields
    lagrangePoints.resize(Ne);  
    lagrangeVel.resize(Ne); 
    lagrangeVelCoup.resize(Ne); 
    lagrangeForceCoup.resize(Ne); 
    lagrangeArea.resize(Ne);    
    Xil.resize(Ne); 
    Xil_0.resize(Ne);
 
    l_el = L/Ne;
    int nl = ceil(l_el/dx_body);
    int nw = ceil(W_el/dx_body);
    double dl = l_el/nl;
    double dw = W_el/nw;

    for (int n = 0; n < Ne; n++)
    {
        lagrangePoints[n] = Eigen::Matrix3Xd::Zero(3,nl*nw);   
        lagrangeVel[n] = Eigen::MatrixXd::Zero(3,nl*nw);   
        lagrangeVelCoup[n] = Eigen::MatrixXd::Zero(3,nl*nw);   
        lagrangeForceCoup[n] = Eigen::MatrixXd::Zero(3,nl*nw);   
        lagrangeArea[n] = Eigen::VectorXd::Zero(nl*nw);   
        
        double l_0 = n*l_el;
        int ind = 0;
        for (int ii = 0; ii < nl; ii++)
        {
            for (int jj = 0; jj < nw; jj++)
            {
                lagrangePoints[n].col(ind) << x_ini, (y_ini - W_el/2.0 + 0.5*dw) + dw*jj, z_ini + l_0 + 0.5*dl + dl*ii;
                lagrangeArea[n](ind) = dw*dl;
                ind++;
            }
        }
    }

    // Initialise relative distance vectors
    for (int eI = 0; eI < Ne; eI++)
    {
        Xil[eI] = Eigen::Matrix3Xd::Zero(3,lagrangePoints[eI].cols());   
        Xil_0[eI] = Eigen::Matrix3Xd::Zero(3,lagrangePoints[eI].cols());   
        
        Eigen::Vector3d cg = (x_el.col(eI+1) + x_el.col(eI))/2.0;

        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            Xil_0[eI].col(pI) << lagrangePoints[eI].col(pI) - cg;
        }
    } 

    // Initialise print
    print_ini(p);
}

void fsi_strip::ini_parallel(lexer *p, fdm *a, ghostcell *pgc)
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

void fsi_strip::get_cellsize(lexer *p, fdm *a, ghostcell *pgc)
{
    Eigen::Vector3d coordI;
    coordI << p->Z11_t[nstrip]/2.0, p->Z11_w[nstrip]/2.0, p->Z11_l[nstrip]/2.0;
    
    if 
    (
        coordI(0) >= xstart[p->mpirank] && coordI(0) < xend[p->mpirank] &&
        coordI(1) >= ystart[p->mpirank] && coordI(1) < yend[p->mpirank] &&
        coordI(2) >= zstart[p->mpirank] && coordI(2) < zend[p->mpirank]
    )
    {
        int ii = p->posc_i(coordI(0));
        int jj = p->posc_j(coordI(1));
        int kk = p->posc_k(coordI(2));
        
        dx_body = p->DXN[ii + marge];
    }
    else
    {
         dx_body = 0.0;
    }

    dx_body = pgc->globalsum(dx_body);
}
	
