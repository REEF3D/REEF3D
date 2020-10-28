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

#include"mooring_DGSEM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void mooring_DGSEM::getL
(
	lexer *p, 
	double **& r_x_, double **& r_y_, double **& r_z_, 
	double **& q_x_, double **& q_y_, double **& q_z_, 
	double **& v_x_, double **& v_y_, double **& v_z_
)
{
	// Calculate tension forces and directions
	
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			qmag[i][j] = 
				sqrt
				(
					  q_x_[i][j]*q_x_[i][j] 
					+ q_y_[i][j]*q_y_[i][j] 
					+ q_z_[i][j]*q_z_[i][j]
				);
				
			T[i][j] = MAX(EA*(qmag[i][j] - 1.0), 0.0);
            
			t_x[i][j] = q_x_[i][j]/qmag[i][j];
			t_y[i][j] = q_y_[i][j]/qmag[i][j];
			t_z[i][j] = q_z_[i][j]/qmag[i][j];
		}
	}


	// Calculate maximum eigenvalues for LLF fluxes at interfaces
	
	lambda_ct[0] = sqrt(EA/gamma);
	lambda_ct[H] = sqrt(EA/gamma);
	
	lambda_cn[0] = sqrt(T[0][0]/(qmag[0][0]*gamma));
	lambda_cn[H] = sqrt(T[H-1][P]/(qmag[H-1][P]*gamma));
	
	for (int i = 1; i < H; i++)
	{
		lambda_ct[i] = sqrt(EA/gamma);
		
		lambda_cn[i] = 
			MAX
			(
				fabs(sqrt(T[i-1][P]/(qmag[i-1][P]*gamma))),
				fabs(sqrt(T[i][0]/(qmag[i][0]*gamma)))
			);
	}

	
	// Compute fluxes at each point
	
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			// fr_x[i][j] = 0.0; 
			// fr_y[i][j] = 0.0; 
			// fr_z[i][j] = 0.0; 
			fq_x[i][j] = -v_x_[i][j]/gamma; 
			fq_y[i][j] = -v_y_[i][j]/gamma; 
			fq_z[i][j] = -v_z_[i][j]/gamma;
			fv_x[i][j] = -T[i][j]*t_x[i][j]; 
			fv_y[i][j] = -T[i][j]*t_y[i][j]; 
			fv_z[i][j] = -T[i][j]*t_z[i][j];
		}
	}


	// Compute fluxes at cell boundaries using LLF fluxes
	
	for (int i = 0; i < H; i++)
	{
		if (i > 0)
		{
			// Left side
			dfr_x[i][0] = -lambda_ct[i]/2.0*(r_x_[i][0] - r_x_[i-1][P]);
			dfr_y[i][0] = -lambda_ct[i]/2.0*(r_y_[i][0] - r_y_[i-1][P]);
			dfr_z[i][0] = -lambda_ct[i]/2.0*(r_z_[i][0] - r_z_[i-1][P]);
			dfq_x[i][0] = +1.0/2.0*(fq_x[i-1][P] - fq_x[i][0]) - lambda_cn[i]/2.0*(q_x_[i][0] - q_x_[i-1][P]);
			dfq_y[i][0] = +1.0/2.0*(fq_y[i-1][P] - fq_y[i][0]) - lambda_cn[i]/2.0*(q_y_[i][0] - q_y_[i-1][P]);
			dfq_z[i][0] = +1.0/2.0*(fq_z[i-1][P] - fq_z[i][0]) - lambda_cn[i]/2.0*(q_z_[i][0] - q_z_[i-1][P]);				
			dfv_x[i][0] = +1.0/2.0*(fv_x[i-1][P] - fv_x[i][0]) - lambda_ct[i]/2.0*(v_x_[i][0] - v_x_[i-1][P]);
			dfv_y[i][0] = +1.0/2.0*(fv_y[i-1][P] - fv_y[i][0]) - lambda_ct[i]/2.0*(v_y_[i][0] - v_y_[i-1][P]);
			dfv_z[i][0] = +1.0/2.0*(fv_z[i-1][P] - fv_z[i][0]) - lambda_ct[i]/2.0*(v_z_[i][0] - v_z_[i-1][P]);	
		}
	
		if (i < (H-1))
		{
			// Right side
			dfr_x[i][1] = -lambda_ct[i+1]/2.0*(r_x_[i][P] - r_x_[i+1][0]);
			dfr_y[i][1] = -lambda_ct[i+1]/2.0*(r_y_[i][P] - r_y_[i+1][0]);
			dfr_z[i][1] = -lambda_ct[i+1]/2.0*(r_z_[i][P] - r_z_[i+1][0]);
			dfq_x[i][1] = +1.0/2.0*(fq_x[i][P] - fq_x[i+1][0]) - lambda_cn[i+1]/2.0*(q_x_[i][P] - q_x_[i+1][0]);
			dfq_y[i][1] = +1.0/2.0*(fq_y[i][P] - fq_y[i+1][0]) - lambda_cn[i+1]/2.0*(q_y_[i][P] - q_y_[i+1][0]);
			dfq_z[i][1] = +1.0/2.0*(fq_z[i][P] - fq_z[i+1][0]) - lambda_cn[i+1]/2.0*(q_z_[i][P] - q_z_[i+1][0]);
			dfv_x[i][1] = +1.0/2.0*(fv_x[i][P] - fv_x[i+1][0]) - lambda_ct[i+1]/2.0*(v_x_[i][P] - v_x_[i+1][0]);
			dfv_y[i][1] = +1.0/2.0*(fv_y[i][P] - fv_y[i+1][0]) - lambda_ct[i+1]/2.0*(v_y_[i][P] - v_y_[i+1][0]);
			dfv_z[i][1] = +1.0/2.0*(fv_z[i][P] - fv_z[i+1][0]) - lambda_ct[i+1]/2.0*(v_z_[i][P] - v_z_[i+1][0]);
		}
	}	
		
	
	// Boundary conditions for inflow and outflow fluxes
	getBoundaryFluxes(p,r_x_,r_y_,r_z_,q_x_,q_y_,q_z_,v_x_,v_y_,v_z_);


	// Calculate external forces
	getForces(r_z_,v_x_,v_y_,v_z_);

	
	// Compute fluxes
	
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			Lr_x[i][j] = rx*sIntdf(dfr_x[i],j) + v_x_[i][j]/gamma;
			Lr_y[i][j] = rx*sIntdf(dfr_y[i],j) + v_y_[i][j]/gamma;
			Lr_z[i][j] = rx*sIntdf(dfr_z[i],j) + v_z_[i][j]/gamma;
					
			Lq_x[i][j] = -rx*Drf(fq_x[i],j) + rx*sIntdf(dfq_x[i],j);
			Lq_y[i][j] = -rx*Drf(fq_y[i],j) + rx*sIntdf(dfq_y[i],j);
			Lq_z[i][j] = -rx*Drf(fq_z[i],j) + rx*sIntdf(dfq_z[i],j);
		
			Lv_x[i][j] = -rx*Drf(fv_x[i],j) + rx*sIntdf(dfv_x[i],j) + Fx[i][j];
			Lv_y[i][j] = -rx*Drf(fv_y[i],j) + rx*sIntdf(dfv_y[i],j) + Fy[i][j];
			Lv_z[i][j] = -rx*Drf(fv_z[i],j) + rx*sIntdf(dfv_z[i],j) + Fz[i][j];			
		}
	}
	
	// Correct source term with Dirichlet bc
	Lr_x[H-1][P] = Lr_x[H-1][P] - v_x_[H-1][P]/gamma + v_xO;
	Lr_y[H-1][P] = Lr_y[H-1][P] - v_y_[H-1][P]/gamma + v_yO;
	Lr_z[H-1][P] = Lr_z[H-1][P] - v_z_[H-1][P]/gamma + v_zO;	
}


double mooring_DGSEM::sIntdf(double *& df_i, int pj)
{
	double mult = 0.0;
	
	for (int i = 0; i < 2; i++)
	{
		mult += sInt[pj][i]*df_i[i];
	}
	
	return mult;
} 


double mooring_DGSEM::Drf(double *& f_i, int pj)
{
	double mult = 0.0;
	
	for (int i = 0; i < (P+1); i++)
	{
		mult += Dr[pj][i]*f_i[i];
	}
	
	return mult;
} 


void mooring_DGSEM::getForces
(
	double **& r_z_, double **& v_x_, double **& v_y_, double **& v_z_
)
{	
	double rho_f = 1000;
	double zg = 0.0;
	double cd_t = 0.5;
	double cd_n = 2.5;
	double cm_t = 0.0;
	double cm_n = 3.8;
	double Kg = 3.0e9;
	double nu = 0.01;
	double mu = 0.3;
	double xi = 1.0;
	
	double vx, vy, vz, vt, vnx, vny, vnz, vn_mag;
	double ax, ay, az, at, anx, any, anz;

	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < (P+1); j++)
		{
			// Gravity
			Fz[i][j] = -9.81*gamma*(rho_c - rho_f)/rho_c;
			
			// Drag
			vx = vf[i][j][0] - v_x_[i][j]/gamma;
			vy = vf[i][j][1] - v_y_[i][j]/gamma;
			vz = vf[i][j][2] - v_z_[i][j]/gamma;

			vt = vx*t_x[i][j] + vy*t_y[i][j] + vz*t_z[i][j];
	
			vnx = vx - vt*t_x[i][j];
			vny = vy - vt*t_y[i][j];
			vnz = vz - vt*t_z[i][j];
			vn_mag = sqrt(vnx*vnx + vny*vny + vnz*vnz);
    
			Fx[i][j] =  0.5*rho_f*d_c*sqrt(1.0 + (qmag[i][j] - 1.0))*(cd_t*fabs(vt)*vt*t_x[i][j] + cd_n*vn_mag*vnx);
			Fy[i][j] =  0.5*rho_f*d_c*sqrt(1.0 + (qmag[i][j] - 1.0))*(cd_t*fabs(vt)*vt*t_y[i][j] + cd_n*vn_mag*vny);
		//	Fz[i][j] += 0.5*rho_f*d_c*sqrt(1.0 + (qmag[i][j] - 1.0))*(cd_t*fabs(vt)*vt*t_z[i][j] + cd_n*vn_mag*vnz);

			// Added mass
			ax = af[i][j][0] - (v_x_[i][j] - v_xn[i][j])/(dt*gamma);
			ay = af[i][j][1] - (v_y_[i][j] - v_yn[i][j])/(dt*gamma);
			az = af[i][j][2] - (v_z_[i][j] - v_zn[i][j])/(dt*gamma);
			
			at = ax*t_x[i][j] + ay*t_y[i][j] + az*t_z[i][j];
			anx = ax - at*t_x[i][j];
			any = ay - at*t_y[i][j];
			anz = az - at*t_z[i][j];
			
			Fx[i][j] += rho_f*PI/4*d_c*d_c*(af[i][j][0] + cm_t*at*t_x[i][j] + cm_n*anx);
			Fy[i][j] += rho_f*PI/4*d_c*d_c*(af[i][j][1] + cm_t*at*t_y[i][j] + cm_n*any);
		//	Fz[i][j] += rho_f*PI/4*d_c*d_c*(af[i][j][2] + cm_t*at*t_z[i][j] + cm_n*anz);
			
			// Bottom
			if ((zg - r_z_[i][j]) >= 0.0)
			{
				Fx[i][j] += mu*(-9.81*gamma*(rho_c - rho_f)/rho_c)*tanh(PI*v_x_[i][j]/(gamma*nu))*v_x_[i][j]/sqrt(v_x_[i][j]*v_x_[i][j]);
				Fy[i][j] += mu*(-9.81*gamma*(rho_c - rho_f)/rho_c)*tanh(PI*v_y_[i][j]/(gamma*nu))*v_y_[i][j]/sqrt(v_y_[i][j]*v_y_[i][j]);
				Fz[i][j] += (Kg*d_c*(zg - r_z_[i][j]) - 2.0*xi*sqrt(Kg*gamma*d_c)*max(v_z_[i][j]/gamma,0.0));
			}
						
			Fx[i][j] = 0.0;
			Fy[i][j] = 0.0;
		}
	}
}














/*
	// Calculate maximum eigenvalues for LLF fluxes at interfaces	
	lambda_ct[0] = sqrt(EA/gamma);
	lambda_ct[H] = sqrt(EA/gamma);
	
	lambda_cn[0] = sqrt(T[0][0]/(qmag[0][0]*gamma));
	lambda_cn[H] = sqrt(T[H-1][P]/(qmag[H-1][P]*gamma));
	
	for (int i = 1; i < H; i++)
	{
		
		double ctM = sqrt(EA*10/gamma*exp(10*(qmag[i-1][P]-1)));
		double ctP = sqrt(EA*10/gamma*exp(10*(qmag[i][0]-1)));
		lambda_ct[i] = MAX(fabs(ctM),fabs(ctP));
		
		lambda_ct[i] = sqrt(EA/gamma);
	
		double cnM = sqrt(T[i-1][P]/(qmag[i-1][P]*gamma));
		double cnP = sqrt(T[i][0]/(qmag[i][0]*gamma));
		
		lambda_cn[i] = MAX(fabs(cnM),fabs(cnP));
	}

	lambda_ct[0] = sqrt(EA*10/gamma*exp(10*(qmag[0][0]-1)));
	lambda_ct[H] = sqrt(EA*10/gamma*exp(10*(qmag[H-1][P]-1)));
*/