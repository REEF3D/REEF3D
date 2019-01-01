/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
#include"mooring_DGSEM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void mooring_DGSEM::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
	// Initialise parameter
	gamma = p->X311_w[line];			// specific weight [kg/m]
	EA = p->X311_EA[line];   			// stiffness times area
	L = p->X311_l[line];      			// length of unstretched cable [m]
	rho_c = p->X311_rho_c[line];   		// density of material [kg/m3]
	d_c = p->X311_d[line];      		// diameter of the cable [m]
	P = p->X311_P[line];				// polynomial order
	H = p->X311_H[line];				// number of elements
	relFac = p->X311_facT[line];		// Relaxation factor for mooring time step

	// Generate face coordinates
	facenodegen(p,0.0,L);

	// Get LGL quadrature points in standard element
	getqp(p);

	// Get Dr matrix
	getDr(p);
	
	// Get surface integrals
	getsurfInt(p);
	
	// Calculate coordinates of all nodes
	nodegen(p);
	
	// Get V and inverse of V
	getV(p);	
	
	// Initialise system
	iniCond(p,a,pgc);
	
	
	// Ini print
	if(p->mpirank==0 && p->P14==1)
	{
		char str[1000];
		sprintf(str,"./REEF3D_6DOF_Mooring/REEF3D_6DOF_mooring_force_%i.dat",line);
		eTout.open(str);
		eTout<<"time \t T"<<endl;	
	}
	printtime = 0.0;
}


void mooring_DGSEM::iniCond(lexer *p, fdm *a, ghostcell *pgc)
{
	t_mooring = 0.0;
	
	xs = p->X311_xs[line];
	ys = p->X311_ys[line];
	zs = p->X311_zs[line];
	
	r_xOn = p->X311_xe[line];
	r_yOn = p->X311_ye[line];
	r_zOn = p->X311_ze[line];
	v_xOn = 0.0;
	v_yOn = 0.0;
	v_zOn = 0.0;
	
	// Fields
	p->Darray(T, H, P+1);
	p->Darray(qmag, H, P+1);
	p->Darray(t_x, H, P+1);
	p->Darray(t_y, H, P+1);
	p->Darray(t_z, H, P+1);
	p->Darray(Fx, H, P+1);
	p->Darray(Fy, H, P+1);
	p->Darray(Fz, H, P+1);
	p->Darray(r_x, H, P+1);
	p->Darray(r_y, H, P+1);
	p->Darray(r_z, H, P+1);
	p->Darray(q_x, H, P+1);
	p->Darray(q_y, H, P+1);
	p->Darray(q_z, H, P+1);
	p->Darray(v_x, H, P+1);
	p->Darray(v_y, H, P+1);
	p->Darray(v_z, H, P+1);
	p->Darray(r_x1, H, P+1);
	p->Darray(r_y1, H, P+1);
	p->Darray(r_z1, H, P+1);
	p->Darray(q_x1, H, P+1);
	p->Darray(q_y1, H, P+1);
	p->Darray(q_z1, H, P+1);
	p->Darray(v_x1, H, P+1);
	p->Darray(v_y1, H, P+1);
	p->Darray(v_z1, H, P+1);
	p->Darray(r_x2, H, P+1);
	p->Darray(r_y2, H, P+1);
	p->Darray(r_z2, H, P+1);
	p->Darray(q_x2, H, P+1);
	p->Darray(q_y2, H, P+1);
	p->Darray(q_z2, H, P+1);
	p->Darray(v_x2, H, P+1);
	p->Darray(v_y2, H, P+1);
	p->Darray(v_z2, H, P+1);
	p->Darray(v_xn, H, P+1);
	p->Darray(v_yn, H, P+1);
	p->Darray(v_zn, H, P+1);
	p->Darray(vf, H, P+1, 3);
	p->Darray(vfn, H, P+1, 3);
	p->Darray(af, H, P+1, 3);
	
	// Fluxes
	p->Darray(fq_x, H, P+1); 
	p->Darray(fq_y, H, P+1); 
	p->Darray(fq_z, H, P+1);
	p->Darray(fv_x, H, P+1); 
	p->Darray(fv_y, H, P+1); 
	p->Darray(fv_z, H, P+1);
	p->Darray(dfr_x, H, 2);
	p->Darray(dfr_y, H, 2);
	p->Darray(dfr_z, H, 2);
	p->Darray(dfq_x, H, 2);
	p->Darray(dfq_y, H, 2);
	p->Darray(dfq_z, H, 2);
	p->Darray(dfv_x, H, 2);
	p->Darray(dfv_y, H, 2);
	p->Darray(dfv_z, H, 2);
	
	// RHS
	p->Darray(Lr_x, H, P+1);
	p->Darray(Lr_y, H, P+1);
	p->Darray(Lr_z, H, P+1);
	p->Darray(Lq_x, H, P+1);
	p->Darray(Lq_y, H, P+1);
	p->Darray(Lq_z, H, P+1);
	p->Darray(Lv_x, H, P+1);
	p->Darray(Lv_y, H, P+1);
	p->Darray(Lv_z, H, P+1);

	// Eigenvalues
	p->Darray(lambda_ct, H+1);
	p->Darray(lambda_cn, H+1);
	
	
	
	// Initial conditions

// Hanging catenary
/*	double xA, xB;
	if (line == 0)
	{
		xA = 3.85;
		xB = 5.85;
	}
	else
	{
		xA = 6.15;
		xB = 8.15;	
	}

	double zA = 0.3;
	double qc = gamma*(rho_c-1000)/rho_c*9.81;
	double va = qc*L/2;

	double H_ = 27.0;
	for (int i=1; i<100; i++)
	{
		double fH = -xB + xA + H_*L/EA + 2*H_/qc*asinh(qc*L/(2*H_));
		double dfH = L/EA + 2/qc*asinh(qc*L/(2*H_)) - L/(H_*sqrt(L*L*qc*qc/(4*H_*H_)+1));
		H_ -= fH/dfH;
	}
	
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < P+1; j++)
		{
			r_x[i][j] = xA + H_/EA*x[i][j] + H_/qc*(asinh((qc*x[i][j]-va)/H_) - asinh(-va/H_));
			q_x[i][j] = H_/EA + 1.0/sqrt(1.0 + ((qc*x[i][j]-va)/H_)*((qc*x[i][j]-va)/H_));
			r_z[i][j] = zA + 1/EA*(qc*x[i][j]*x[i][j]/2.0 - va*x[i][j]) + 1/qc*(sqrt(H_*H_+(qc*x[i][j]-va)*(qc*x[i][j]-va)) - sqrt(H_*H_+va*va));
			q_z[i][j] = (qc*x[i][j]-va)/EA + (qc*x[i][j]-va)/sqrt(H_*H_+(qc*x[i][j]-va)*(qc*x[i][j]-va));
		}
	}
*/

// Elastic catenary
	double FH, FV, d_xy;
	pcatenary = new mooring_Catenary(line);
	pcatenary->iniDyn(p, a, pgc, FH, FV);
	
	double dx = p->X311_xe[line] - p->X311_xs[line];			
	double dy = p->X311_ye[line] - p->X311_ys[line];	
	double alpha = atan(dy/dx);
	double w_ = p->X311_w[line]*9.81*(rho_c - 1000.0)/rho_c;
	double Lb = L - FV/w_;
	
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < P+1; j++)
		{
			if (x[i][j] <= Lb)
			{
				r_z[i][j] = zs;
				q_z[i][j] = 0.0;
				T[i][j] = fabs(FH);
					
				d_xy = x[i][j]*(1.0 + FH/EA);
					
				if (dx > 0)
				{
					r_x[i][j] = xs + d_xy*cos(alpha);
					q_x[i][j] = (1.0 + FH/EA)*cos(alpha);
				}
				else
				{
					r_x[i][j] = xs - d_xy*cos(alpha);
					q_x[i][j] = -(1.0 + FH/EA)*cos(alpha);
				}
					
				if (dy > 0)
				{
					r_y[i][j] = ys + d_xy*sin(alpha);
					q_y[i][j] = (1.0 + FH/EA)*sin(alpha);
				}
				else
				{
					r_y[i][j] = ys - d_xy*sin(alpha);
					q_y[i][j] = -(1.0 + FH/EA)*sin(alpha);
				}	
			}
			else
			{
				r_z[i][j] = zs + FH/w_*(sqrt(1+(w_*(x[i][j]-Lb)/FH)*(w_*(x[i][j]-Lb)/FH)) - 1) + w_*(x[i][j]-Lb)*(x[i][j]-Lb)/(2*EA);
				q_z[i][j] = w_*(x[i][j]-Lb)/EA + w_*(x[i][j]-Lb)/(FH*sqrt((w_*w_*(x[i][j]-Lb)*(x[i][j]-Lb))/(FH*FH) + 1));				
				T[i][j] = sqrt(FH*FH + (w_*(x[i][j] - Lb))*(w_*(x[i][j] - Lb)));
					
				d_xy = FH/w_*log(w_*(x[i][j]-Lb)/FH + sqrt(1+(w_*(x[i][j]-Lb)/FH)*(w_*(x[i][j]-Lb)/FH))) + FH*x[i][j]/EA;
			
				if (dx > 0)
				{
					r_x[i][j] = xs + Lb*cos(alpha) + d_xy*cos(alpha);
					q_x[i][j] = (FH / sqrt(FH*FH + w_*w_*(Lb-x[i][j])*(Lb-x[i][j])) + FH/EA)*cos(alpha);
				}
				else
				{
					r_x[i][j] = xs - Lb*cos(alpha) - d_xy*cos(alpha);
					q_x[i][j] = -(FH / sqrt(FH*FH + w_*w_*(Lb-x[i][j])*(Lb-x[i][j])) + FH/EA)*cos(alpha);
				}
				
				if (dy > 0)
				{
					r_y[i][j] = ys + Lb*sin(alpha) + d_xy*sin(alpha);
					q_y[i][j] = (FH / sqrt(FH*FH + w_*w_*(Lb-x[i][j])*(Lb-x[i][j])) + FH/EA)*sin(alpha);
				}
				else
				{
					r_y[i][j] = ys - Lb*sin(alpha) - d_xy*sin(alpha);
					q_y[i][j] = -(FH / sqrt(FH*FH + w_*w_*(Lb-x[i][j])*(Lb-x[i][j])) + FH/EA)*sin(alpha);
				}				
			}
		}	
	}

	delete pcatenary; 
	pcatenary = NULL;

// Benchmark vibration
/*
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			r_z[i][j] = 1.0*sin(PI/100*x[i][j]);  
			q_z[i][j] = PI/100*cos(PI/100*x[i][j]);
			r_x[i][j] = sqrt(x[i][j]*x[i][j]-1.0*sin(PI/100*x[i][j])*1.0*sin(PI/100*x[i][j]));
			q_x[i][j] = 1.1;
			
			r_y[i][j] = 0.0;
			q_y[i][j] = 0.0;
		}
	}
*/
// Benchmark linear shock
/*
	for (int i = 0; i < H; i++)
	{
		double x_S = 115/H*i;
		double x_E = 115/H*(i+1);
		
		for (int j = 0; j < (P+1); j++)
		{			
			r_y[i][j] = 0.0;
			q_y[i][j] = 0.0;
			r_z[i][j] = 0.0;  
			q_z[i][j] = 0.0;
			
			r_x[i][j] = x_S + (x_E - x_S)/P*(j);
			
			if (x[i][j] < 50)
			{
				q_x[i][j] = 1.1;
			}
			else
			{
				q_x[i][j] = 1.2;
			}
		}
	}
*/

// Benchmark nonlinear shock
/*
	for (int i = 0; i < H; i++)
	{
		double x_S = 1000/H*i;
		double x_E = 1000/H*(i+1);
		
		for (int j = 0; j < (P+1); j++)
		{			
			r_y[i][j] = 0.0;
			q_y[i][j] = 0.0;
			r_z[i][j] = 0.0;  
			q_z[i][j] = 0.0;
			
			r_x[i][j] = x_S + (x_E - x_S)/P*(j);
			q_x[i][j] = 1.1;
		}
	}
	q_x[H-1][P] = 1.2;
*/
	
	// Apply minMod slope limiter
	//limitFields(p,r_x,r_y,r_z,q_x,q_y,q_z,v_x,v_y,v_z);

	
	// Calculate CFL number
	if (P < 3)
	{
		cfl = (vx[1]-vx[0])/((2.0*P+1.0)*sqrt(EA/gamma));
	}
	else
	{
		cfl = (vx[1]-vx[0])/(P*P*sqrt(EA/gamma));
	}
	
	// Initialise communication 
	ini_parallel(p, a, pgc);
}


void mooring_DGSEM::ini_parallel(lexer *p, fdm *a, ghostcell *pgc)
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
