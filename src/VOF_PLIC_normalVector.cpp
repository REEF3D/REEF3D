/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"


void VOF_PLIC::calcNormalFO(fdm* a, lexer* p)
{
    // 1st-order method 
    nx(i,j,k) = 
        (a->vof(i-1,j-1,k-1)+a->vof(i-1,j-1,k+1)+a->vof(i-1,j+1,k-1)
        +a->vof(i-1,j+1,k+1)+2.0*(a->vof(i-1,j-1,k)+a->vof(i-1,j+1,k)
        +a->vof(i-1,j,k-1)+a->vof(i-1,j,k+1))+4.0*a->vof(i-1,j,k)) 
        - (a->vof(i+1,j-1,k-1)+a->vof(i+1,j-1,k+1)+a->vof(i+1,j+1,k-1)
        +a->vof(i+1,j+1,k+1)+2.0*(a->vof(i+1,j-1,k)+a->vof(i+1,j+1,k)
        +a->vof(i+1,j,k-1)+a->vof(i+1,j,k+1))+4.0*a->vof(i+1,j,k));
                 
    ny(i,j,k) = 
        (a->vof(i-1,j-1,k-1)+a->vof(i-1,j-1,k+1)+a->vof(i+1,j-1,k-1)
        +a->vof(i+1,j-1,k+1)+2.0*(a->vof(i-1,j-1,k)+a->vof(i+1,j-1,k)
        +a->vof(i,j-1,k-1)+a->vof(i,j-1,k+1))+4.0*a->vof(i,j-1,k)) 
        - (a->vof(i-1,j+1,k-1)+a->vof(i-1,j+1,k+1)+a->vof(i+1,j+1,k-1)
        +a->vof(i+1,j+1,k+1)+2.0*(a->vof(i-1,j+1,k)+a->vof(i+1,j+1,k)
        +a->vof(i,j+1,k-1)+a->vof(i,j+1,k+1))+4.0*a->vof(i,j+1,k));

    nz(i,j,k) = 
        (a->vof(i-1,j-1,k-1)+a->vof(i-1,j+1,k-1)+a->vof(i+1,j-1,k-1)
        +a->vof(i+1,j+1,k-1)+2.0*(a->vof(i-1,j,k-1)+a->vof(i+1,j,k-1)
        +a->vof(i,j-1,k-1)+a->vof(i,j+1,k-1))+4.0*a->vof(i,j,k-1)) 
        - ( a->vof(i-1,j-1,k+1)+a->vof(i-1,j+1,k+1)+a->vof(i+1,j-1,k+1)
        +a->vof(i+1,j+1,k+1)+2.0*(a->vof(i-1,j,k+1)+a->vof(i+1,j,k+1)
        +a->vof(i,j-1,k+1)+a->vof(i,j+1,k+1))+4.0*a->vof(i,j,k+1));	
}



void VOF_PLIC::ininorVecLS(lexer* p)
{
	double *wn;
	double **dPN, **G, **invG;
	
	p->Darray(wn, 26);
	p->Darray(dPN, 26, 3);
	p->Darray(G, 3, 3);
	p->Darray(invG, 3, 3);
	
	p->Darray(nxCoeff, p->knox, p->knoy, p->knoz, 26);
	p->Darray(nyCoeff, p->knox, p->knoy, p->knoz, 26);
	p->Darray(nzCoeff, p->knox, p->knoy, p->knoz, 26);
	
	LOOP
	{
		//- Create distance stencil
		
		dPN[0][0] = p->DXP[IP]; 	dPN[0][1] = 0.0; 			dPN[0][2] = 0.0;
		dPN[1][0] = p->DXP[IP]; 	dPN[1][1] = p->DYP[JP]; 	dPN[1][2] = 0.0;
		dPN[2][0] = p->DXP[IP]; 	dPN[2][1] = -p->DYP[JM1]; 	dPN[2][2] = 0.0;
		dPN[3][0] = -p->DXP[IM1]; 	dPN[3][1] = 0.0; 			dPN[3][2] = 0.0;
		dPN[4][0] = -p->DXP[IM1]; 	dPN[4][1] = p->DYP[JP]; 	dPN[4][2] = 0.0;
		dPN[5][0] = -p->DXP[IM1]; 	dPN[5][1] = -p->DYP[JM1]; 	dPN[5][2] = 0.0;
		dPN[6][0] = 0.0; 			dPN[6][1] = p->DYP[JP]; 	dPN[6][2] = 0.0;
		dPN[7][0] = 0.0; 			dPN[7][1] = -p->DYP[JM1]; 	dPN[7][2] = 0.0;
		dPN[8][0] = 0.0; 			dPN[8][1] = 0.0; 			dPN[8][2] = p->DZP[KP];	
		dPN[9][0] = p->DXP[IP]; 	dPN[9][1] = 0.0; 			dPN[9][2] = p->DZP[KP];
		dPN[10][0] = p->DXP[IP]; 	dPN[10][1] = p->DYP[JP]; 	dPN[10][2] = p->DZP[KP];
		dPN[11][0] = p->DXP[IP]; 	dPN[11][1] = -p->DYP[JM1]; dPN[11][2] = p->DZP[KP];
		dPN[12][0] = -p->DXP[IM1]; dPN[12][1] = 0.0; 			dPN[12][2] = p->DZP[KP];
		dPN[13][0] = -p->DXP[IM1]; dPN[13][1] = p->DYP[JP]; 	dPN[13][2] = p->DZP[KP];
		dPN[14][0] = -p->DXP[IM1]; dPN[14][1] = -p->DYP[JM1]; dPN[14][2] = p->DZP[KP];
		dPN[15][0] = 0.0; 			dPN[15][1] = p->DYP[JP]; 	dPN[15][2] = p->DZP[KP];
		dPN[16][0] = 0.0; 			dPN[16][1] = -p->DYP[JM1]; dPN[16][2] = p->DZP[KP];	
		dPN[17][0] = 0.0; 			dPN[17][1] = 0.0; 			dPN[17][2] = -p->DZP[KM1];	
		dPN[18][0] = p->DXP[IP]; 	dPN[18][1] = 0.0; 			dPN[18][2] = -p->DZP[KM1];
		dPN[19][0] = p->DXP[IP]; 	dPN[19][1] = p->DYP[JP]; 	dPN[19][2] = -p->DZP[KM1];
		dPN[20][0] = p->DXP[IP]; 	dPN[20][1] = -p->DYP[JM1]; dPN[20][2] = -p->DZP[KM1];
		dPN[21][0] = -p->DXP[IM1]; dPN[21][1] = 0.0; 			dPN[21][2] = -p->DZP[KM1];
		dPN[22][0] = -p->DXP[IM1]; dPN[22][1] = p->DYP[JP]; 	dPN[22][2] = -p->DZP[KM1];
		dPN[23][0] = -p->DXP[IM1]; dPN[23][1] = -p->DYP[JM1]; dPN[23][2] = -p->DZP[KM1];
		dPN[24][0] = 0.0; 			dPN[24][1] = p->DYP[JP]; 	dPN[24][2] = -p->DZP[KM1];
		dPN[25][0] = 0.0; 			dPN[25][1] = -p->DYP[JM1]; dPN[25][2] = -p->DZP[KM1];
		
		
		//- Calculate weights and least-squares matrix
		
		for (int n = 0; n < 3; n++)
		{
			for (int m = 0; m < 3; m++)
			{
				G[n][m] = 0.0;
			}
		}
		
		for (int q = 0; q < 26; q++)
		{
			wn[q] = 1.0/sqrt(dPN[q][0]*dPN[q][0] + dPN[q][1]*dPN[q][1] + dPN[q][2]*dPN[q][2]);
			
			for (int n = 0; n < 3; n++)
			{
				for (int m = 0; m < 3; m++)
				{
					G[n][m] += wn[q]*wn[q]*dPN[q][n]*dPN[q][m];
				}
			}
		}

		
		//- Invert G
		
		double det = 
			 G[0][0]*(G[1][1]*G[2][2] - G[2][1]*G[1][2]) 
			- G[0][1]*(G[1][0]*G[2][2] - G[1][2]*G[2][0]) 
			+ G[0][2]*(G[1][0]*G[2][1] - G[1][1]*G[2][0]);

		double invdet = 1.0/det;

		invG[0][0] = (G[1][1]*G[2][2] - G[2][1]*G[1][2])*invdet;
		invG[0][1] = (G[0][2]*G[2][1] - G[0][1]*G[2][2])*invdet;
		invG[0][2] = (G[0][1]*G[1][2] - G[0][2]*G[1][1])*invdet;
		invG[1][0] = (G[1][2]*G[2][0] - G[1][0]*G[2][2])*invdet;
		invG[1][1] = (G[0][0]*G[2][2] - G[0][2]*G[2][0])*invdet;
		invG[1][2] = (G[1][0]*G[0][2] - G[0][0]*G[1][2])*invdet;
		invG[2][0] = (G[1][0]*G[2][0] - G[2][0]*G[1][1])*invdet;
		invG[2][1] = (G[2][0]*G[0][1] - G[0][0]*G[2][1])*invdet;
		invG[2][2] = (G[0][0]*G[1][1] - G[1][0]*G[0][1])*invdet;


		//- Calculate coefficients for normal vector calculation
		
		for (int q = 0; q < 26; q++)
		{
			nxCoeff[i][j][k][q] = wn[q]*wn[q]*(invG[0][0]*dPN[q][0] + invG[0][1]*dPN[q][1] + invG[0][2]*dPN[q][2]);
			nyCoeff[i][j][k][q] = wn[q]*wn[q]*(invG[1][0]*dPN[q][0] + invG[1][1]*dPN[q][1] + invG[1][2]*dPN[q][2]);
			nzCoeff[i][j][k][q] = wn[q]*wn[q]*(invG[2][0]*dPN[q][0] + invG[2][1]*dPN[q][1] + invG[2][2]*dPN[q][2]);
		}
	}


	//- Delete arrays
	
	p->del_Darray(wn, 26);
	p->del_Darray(dPN, 26, 3);
	p->del_Darray(G, 3, 3);
	p->del_Darray(invG, 3, 3);	
}

void VOF_PLIC::calcNormalLS(fdm* a, lexer* p)
{
	nx(i,j,k) = -1.0*(
		nxCoeff[i][j][k][0]*(a->vof(i+1,j,k) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][1]*(a->vof(i+1,j+1,k) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][2]*(a->vof(i+1,j-1,k) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][3]*(a->vof(i-1,j,k) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][4]*(a->vof(i-1,j+1,k) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][5]*(a->vof(i-1,j-1,k) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][6]*(a->vof(i,j+1,k) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][7]*(a->vof(i,j-1,k) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][8]*(a->vof(i,j,k+1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][9]*(a->vof(i+1,j,k+1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][10]*(a->vof(i+1,j+1,k+1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][11]*(a->vof(i+1,j-1,k+1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][12]*(a->vof(i-1,j,k+1) - a->vof(i,j,k))		
		+ nxCoeff[i][j][k][13]*(a->vof(i-1,j+1,k+1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][14]*(a->vof(i-1,j-1,k+1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][15]*(a->vof(i,j+1,k+1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][16]*(a->vof(i,j-1,k+1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][17]*(a->vof(i,j,k-1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][18]*(a->vof(i+1,j,k-1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][19]*(a->vof(i+1,j+1,k-1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][20]*(a->vof(i+1,j-1,k-1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][21]*(a->vof(i-1,j,k-1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][22]*(a->vof(i-1,j+1,k-1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][23]*(a->vof(i-1,j-1,k-1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][24]*(a->vof(i,j+1,k-1) - a->vof(i,j,k))
		+ nxCoeff[i][j][k][25]*(a->vof(i,j-1,k-1) - a->vof(i,j,k)));
	
	ny(i,j,k) = -1.0*(
		nyCoeff[i][j][k][0]*(a->vof(i+1,j,k) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][1]*(a->vof(i+1,j+1,k) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][2]*(a->vof(i+1,j-1,k) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][3]*(a->vof(i-1,j,k) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][4]*(a->vof(i-1,j+1,k) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][5]*(a->vof(i-1,j-1,k) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][6]*(a->vof(i,j+1,k) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][7]*(a->vof(i,j-1,k) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][8]*(a->vof(i,j,k+1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][9]*(a->vof(i+1,j,k+1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][10]*(a->vof(i+1,j+1,k+1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][11]*(a->vof(i+1,j-1,k+1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][12]*(a->vof(i-1,j,k+1) - a->vof(i,j,k))		
		+ nyCoeff[i][j][k][13]*(a->vof(i-1,j+1,k+1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][14]*(a->vof(i-1,j-1,k+1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][15]*(a->vof(i,j+1,k+1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][16]*(a->vof(i,j-1,k+1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][17]*(a->vof(i,j,k-1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][18]*(a->vof(i+1,j,k-1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][19]*(a->vof(i+1,j+1,k-1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][20]*(a->vof(i+1,j-1,k-1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][21]*(a->vof(i-1,j,k-1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][22]*(a->vof(i-1,j+1,k-1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][23]*(a->vof(i-1,j-1,k-1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][24]*(a->vof(i,j+1,k-1) - a->vof(i,j,k))
		+ nyCoeff[i][j][k][25]*(a->vof(i,j-1,k-1) - a->vof(i,j,k)));

	nz(i,j,k) = -1.0*(
		nzCoeff[i][j][k][0]*(a->vof(i+1,j,k) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][1]*(a->vof(i+1,j+1,k) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][2]*(a->vof(i+1,j-1,k) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][3]*(a->vof(i-1,j,k) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][4]*(a->vof(i-1,j+1,k) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][5]*(a->vof(i-1,j-1,k) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][6]*(a->vof(i,j+1,k) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][7]*(a->vof(i,j-1,k) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][8]*(a->vof(i,j,k+1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][9]*(a->vof(i+1,j,k+1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][10]*(a->vof(i+1,j+1,k+1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][11]*(a->vof(i+1,j-1,k+1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][12]*(a->vof(i-1,j,k+1) - a->vof(i,j,k))		
		+ nzCoeff[i][j][k][13]*(a->vof(i-1,j+1,k+1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][14]*(a->vof(i-1,j-1,k+1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][15]*(a->vof(i,j+1,k+1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][16]*(a->vof(i,j-1,k+1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][17]*(a->vof(i,j,k-1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][18]*(a->vof(i+1,j,k-1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][19]*(a->vof(i+1,j+1,k-1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][20]*(a->vof(i+1,j-1,k-1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][21]*(a->vof(i-1,j,k-1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][22]*(a->vof(i-1,j+1,k-1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][23]*(a->vof(i-1,j-1,k-1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][24]*(a->vof(i,j+1,k-1) - a->vof(i,j,k))
		+ nzCoeff[i][j][k][25]*(a->vof(i,j-1,k-1) - a->vof(i,j,k)));	
}


void VOF_PLIC::calcNormalWENO(fdm* a, lexer* p)
{
    //- WENO gradient scheme
    nx(i,j,k) = -normvec_x(a, a->vof);
    ny(i,j,k) = -normvec_y(a, a->vof);
    nz(i,j,k) = -normvec_z(a, a->vof);  
}


void VOF_PLIC::calcNormalPhi(fdm* a, lexer* p)
{
    nx(i,j,k) = 
        (a->phi(i-1,j-1,k-1)+a->phi(i-1,j-1,k+1)+a->phi(i-1,j+1,k-1)
        +a->phi(i-1,j+1,k+1)+2.0*(a->phi(i-1,j-1,k)+a->phi(i-1,j+1,k)
        +a->phi(i-1,j,k-1)+a->phi(i-1,j,k+1))+4.0*a->phi(i-1,j,k)) 
        - (a->phi(i+1,j-1,k-1)+a->phi(i+1,j-1,k+1)+a->phi(i+1,j+1,k-1)
        +a->phi(i+1,j+1,k+1)+2.0*(a->phi(i+1,j-1,k)+a->phi(i+1,j+1,k)
        +a->phi(i+1,j,k-1)+a->phi(i+1,j,k+1))+4.0*a->phi(i+1,j,k));
                 
    ny(i,j,k) = 
        (a->phi(i-1,j-1,k-1)+a->phi(i-1,j-1,k+1)+a->phi(i+1,j-1,k-1)
        +a->phi(i+1,j-1,k+1)+2.0*(a->phi(i-1,j-1,k)+a->phi(i+1,j-1,k)
        +a->phi(i,j-1,k-1)+a->phi(i,j-1,k+1))+4.0*a->phi(i,j-1,k)) 
        - (a->phi(i-1,j+1,k-1)+a->phi(i-1,j+1,k+1)+a->phi(i+1,j+1,k-1)
        +a->phi(i+1,j+1,k+1)+2.0*(a->phi(i-1,j+1,k)+a->phi(i+1,j+1,k)
        +a->phi(i,j+1,k-1)+a->phi(i,j+1,k+1))+4.0*a->phi(i,j+1,k));

    nz(i,j,k) = 
        (a->phi(i-1,j-1,k-1)+a->phi(i-1,j+1,k-1)+a->phi(i+1,j-1,k-1)
        +a->phi(i+1,j+1,k-1)+2.0*(a->phi(i-1,j,k-1)+a->phi(i+1,j,k-1)
        +a->phi(i,j-1,k-1)+a->phi(i,j+1,k-1))+4.0*a->phi(i,j,k-1)) 
        - ( a->phi(i-1,j-1,k+1)+a->phi(i-1,j+1,k+1)+a->phi(i+1,j-1,k+1)
        +a->phi(i+1,j+1,k+1)+2.0*(a->phi(i-1,j,k+1)+a->phi(i+1,j,k+1)
        +a->phi(i,j-1,k+1)+a->phi(i,j+1,k+1))+4.0*a->phi(i,j,k+1)); 
}


