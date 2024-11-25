/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Tobias Martin, Fabian Knoblauch
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

void VOF_PLIC::calcNormalFO(fdm* a, lexer* p, field4 voffield)
{   double nsum;
    // 1st-order method 
    nx(i,j,k) = 
        (voffield(i-1,j-1,k-1)+voffield(i-1,j-1,k+1)+voffield(i-1,j+1,k-1)
        +voffield(i-1,j+1,k+1)+2.0*(voffield(i-1,j-1,k)+voffield(i-1,j+1,k)
        +voffield(i-1,j,k-1)+voffield(i-1,j,k+1))+4.0*voffield(i-1,j,k)) 
        - (voffield(i+1,j-1,k-1)+voffield(i+1,j-1,k+1)+voffield(i+1,j+1,k-1)
        +voffield(i+1,j+1,k+1)+2.0*(voffield(i+1,j-1,k)+voffield(i+1,j+1,k)
        +voffield(i+1,j,k-1)+voffield(i+1,j,k+1))+4.0*voffield(i+1,j,k));
                 
    ny(i,j,k) = 
        (voffield(i-1,j-1,k-1)+voffield(i-1,j-1,k+1)+voffield(i+1,j-1,k-1)
        +voffield(i+1,j-1,k+1)+2.0*(voffield(i-1,j-1,k)+voffield(i+1,j-1,k)
        +voffield(i,j-1,k-1)+voffield(i,j-1,k+1))+4.0*voffield(i,j-1,k)) 
        - (voffield(i-1,j+1,k-1)+voffield(i-1,j+1,k+1)+voffield(i+1,j+1,k-1)
        +voffield(i+1,j+1,k+1)+2.0*(voffield(i-1,j+1,k)+voffield(i+1,j+1,k)
        +voffield(i,j+1,k-1)+voffield(i,j+1,k+1))+4.0*voffield(i,j+1,k));

    nz(i,j,k) = 
        (voffield(i-1,j-1,k-1)+voffield(i-1,j+1,k-1)+voffield(i+1,j-1,k-1)
        +voffield(i+1,j+1,k-1)+2.0*(voffield(i-1,j,k-1)+voffield(i+1,j,k-1)
        +voffield(i,j-1,k-1)+voffield(i,j+1,k-1))+4.0*voffield(i,j,k-1)) 
        - ( voffield(i-1,j-1,k+1)+voffield(i-1,j+1,k+1)+voffield(i+1,j-1,k+1)
        +voffield(i+1,j+1,k+1)+2.0*(voffield(i-1,j,k+1)+voffield(i+1,j,k+1)
        +voffield(i,j-1,k+1)+voffield(i,j+1,k+1))+4.0*voffield(i,j,k+1));	
        
        nsum=sqrt(nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k));
        nx(i,j,k)=nx(i,j,k)/nsum;
        ny(i,j,k)=ny(i,j,k)/nsum;
        nz(i,j,k)=nz(i,j,k)/nsum;
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

void VOF_PLIC::calcNormalLS(fdm* a, lexer* p, field4 voffield)
{   double nsum;
	nx(i,j,k) = -1.0*(
		nxCoeff[i][j][k][0]*(voffield(i+1,j,k) - voffield(i,j,k))
		+ nxCoeff[i][j][k][1]*(voffield(i+1,j+1,k) - voffield(i,j,k))
		+ nxCoeff[i][j][k][2]*(voffield(i+1,j-1,k) - voffield(i,j,k))
		+ nxCoeff[i][j][k][3]*(voffield(i-1,j,k) - voffield(i,j,k))
		+ nxCoeff[i][j][k][4]*(voffield(i-1,j+1,k) - voffield(i,j,k))
		+ nxCoeff[i][j][k][5]*(voffield(i-1,j-1,k) - voffield(i,j,k))
		+ nxCoeff[i][j][k][6]*(voffield(i,j+1,k) - voffield(i,j,k))
		+ nxCoeff[i][j][k][7]*(voffield(i,j-1,k) - voffield(i,j,k))
		+ nxCoeff[i][j][k][8]*(voffield(i,j,k+1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][9]*(voffield(i+1,j,k+1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][10]*(voffield(i+1,j+1,k+1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][11]*(voffield(i+1,j-1,k+1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][12]*(voffield(i-1,j,k+1) - voffield(i,j,k))		
		+ nxCoeff[i][j][k][13]*(voffield(i-1,j+1,k+1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][14]*(voffield(i-1,j-1,k+1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][15]*(voffield(i,j+1,k+1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][16]*(voffield(i,j-1,k+1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][17]*(voffield(i,j,k-1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][18]*(voffield(i+1,j,k-1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][19]*(voffield(i+1,j+1,k-1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][20]*(voffield(i+1,j-1,k-1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][21]*(voffield(i-1,j,k-1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][22]*(voffield(i-1,j+1,k-1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][23]*(voffield(i-1,j-1,k-1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][24]*(voffield(i,j+1,k-1) - voffield(i,j,k))
		+ nxCoeff[i][j][k][25]*(voffield(i,j-1,k-1) - voffield(i,j,k)));
	
	ny(i,j,k) = -1.0*(
		nyCoeff[i][j][k][0]*(voffield(i+1,j,k) - voffield(i,j,k))
		+ nyCoeff[i][j][k][1]*(voffield(i+1,j+1,k) - voffield(i,j,k))
		+ nyCoeff[i][j][k][2]*(voffield(i+1,j-1,k) - voffield(i,j,k))
		+ nyCoeff[i][j][k][3]*(voffield(i-1,j,k) - voffield(i,j,k))
		+ nyCoeff[i][j][k][4]*(voffield(i-1,j+1,k) - voffield(i,j,k))
		+ nyCoeff[i][j][k][5]*(voffield(i-1,j-1,k) - voffield(i,j,k))
		+ nyCoeff[i][j][k][6]*(voffield(i,j+1,k) - voffield(i,j,k))
		+ nyCoeff[i][j][k][7]*(voffield(i,j-1,k) - voffield(i,j,k))
		+ nyCoeff[i][j][k][8]*(voffield(i,j,k+1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][9]*(voffield(i+1,j,k+1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][10]*(voffield(i+1,j+1,k+1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][11]*(voffield(i+1,j-1,k+1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][12]*(voffield(i-1,j,k+1) - voffield(i,j,k))		
		+ nyCoeff[i][j][k][13]*(voffield(i-1,j+1,k+1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][14]*(voffield(i-1,j-1,k+1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][15]*(voffield(i,j+1,k+1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][16]*(voffield(i,j-1,k+1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][17]*(voffield(i,j,k-1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][18]*(voffield(i+1,j,k-1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][19]*(voffield(i+1,j+1,k-1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][20]*(voffield(i+1,j-1,k-1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][21]*(voffield(i-1,j,k-1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][22]*(voffield(i-1,j+1,k-1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][23]*(voffield(i-1,j-1,k-1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][24]*(voffield(i,j+1,k-1) - voffield(i,j,k))
		+ nyCoeff[i][j][k][25]*(voffield(i,j-1,k-1) - voffield(i,j,k)));

	nz(i,j,k) = -1.0*(
		nzCoeff[i][j][k][0]*(voffield(i+1,j,k) - voffield(i,j,k))
		+ nzCoeff[i][j][k][1]*(voffield(i+1,j+1,k) - voffield(i,j,k))
		+ nzCoeff[i][j][k][2]*(voffield(i+1,j-1,k) - voffield(i,j,k))
		+ nzCoeff[i][j][k][3]*(voffield(i-1,j,k) - voffield(i,j,k))
		+ nzCoeff[i][j][k][4]*(voffield(i-1,j+1,k) - voffield(i,j,k))
		+ nzCoeff[i][j][k][5]*(voffield(i-1,j-1,k) - voffield(i,j,k))
		+ nzCoeff[i][j][k][6]*(voffield(i,j+1,k) - voffield(i,j,k))
		+ nzCoeff[i][j][k][7]*(voffield(i,j-1,k) - voffield(i,j,k))
		+ nzCoeff[i][j][k][8]*(voffield(i,j,k+1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][9]*(voffield(i+1,j,k+1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][10]*(voffield(i+1,j+1,k+1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][11]*(voffield(i+1,j-1,k+1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][12]*(voffield(i-1,j,k+1) - voffield(i,j,k))		
		+ nzCoeff[i][j][k][13]*(voffield(i-1,j+1,k+1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][14]*(voffield(i-1,j-1,k+1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][15]*(voffield(i,j+1,k+1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][16]*(voffield(i,j-1,k+1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][17]*(voffield(i,j,k-1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][18]*(voffield(i+1,j,k-1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][19]*(voffield(i+1,j+1,k-1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][20]*(voffield(i+1,j-1,k-1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][21]*(voffield(i-1,j,k-1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][22]*(voffield(i-1,j+1,k-1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][23]*(voffield(i-1,j-1,k-1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][24]*(voffield(i,j+1,k-1) - voffield(i,j,k))
		+ nzCoeff[i][j][k][25]*(voffield(i,j-1,k-1) - voffield(i,j,k)));	
        
        nsum=sqrt(nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k));
        nx(i,j,k)=nx(i,j,k)/nsum;
        ny(i,j,k)=ny(i,j,k)/nsum;
        nz(i,j,k)=nz(i,j,k)/nsum;
}


void VOF_PLIC::calcNormalWENO(fdm* a, lexer* p, field4 voffield)
{
    //- WENO gradient scheme
    double nsum;
    nx(i,j,k) = -normvec_x(a, voffield);
    ny(i,j,k) = -normvec_y(a, voffield);
    nz(i,j,k) = -normvec_z(a, voffield);  
    nsum=sqrt(nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k));
    nx(i,j,k)=nx(i,j,k)/nsum;
    ny(i,j,k)=ny(i,j,k)/nsum;
    nz(i,j,k)=nz(i,j,k)/nsum;
    
}

void VOF_PLIC::calcNormalMassCentre(fdm* a, lexer* p, field4 voffield)
{
    double nsum,invvec_x,invvec_y, invvec_z; 
    double Vsum=0.0;
    double Vxsum=0.0;
    double Vysum=0.0;
    double Vzsum=0.0;
    
    // centre
    Vsum+=voffield(i,j,k)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
    // x*V = y*V = z*V = 0.0
    
    //north
    Vsum+=voffield(i+1,j,k)*p->DXN[IP1]*p->DYN[JP]*p->DZN[KP];
    Vxsum+=voffield(i+1,j,k)*p->DXN[IP1]*p->DYN[JP]*p->DZN[KP]   *p->DXP[IP];
    // y*V=0, z*V=0
    
    //south
    Vsum+=voffield(i-1,j,k)*p->DXN[IM1]*p->DYN[JP]*p->DZN[KP];
    Vxsum+=voffield(i-1,j,k)*p->DXN[IM1]*p->DYN[JP]*p->DZN[KP]   * -p->DXP[IM1];
     // y*V=0, z*V=0
     
     //top
     Vsum+=voffield(i,j,k+1)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP1];
     Vzsum+=voffield(i,j,k+1)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP1]  *p->DZP[KP];
     // x*V=0, y*V=0
     
     //bottom
     Vsum+=voffield(i,j,k-1)*p->DXN[IP]*p->DYN[JP]*p->DZN[KM1];
     Vzsum+=voffield(i,j,k-1)*p->DXN[IP]*p->DYN[JP]*p->DZN[KM1]  * -p->DZP[KM1];
     // x*V=0, y*V=0
     
     //north-top
     Vsum+=voffield(i+1,j,k+1)*p->DXN[IP1]*p->DYN[JP]*p->DZN[KP1];
     Vxsum+=voffield(i+1,j,k+1)*p->DXN[IP1]*p->DYN[JP]*p->DZN[KP1] *p->DXP[IP];
     Vzsum+=voffield(i+1,j,k+1)*p->DXN[IP1]*p->DYN[JP]*p->DZN[KP1] *p->DZP[KP];
     // y*V=0
     
     //north-bottom
     Vsum+=voffield(i+1,j,k-1)*p->DXN[IP1]*p->DYN[JP]*p->DZN[KM1];
     Vxsum+=voffield(i+1,j,k-1)*p->DXN[IP1]*p->DYN[JP]*p->DZN[KM1] *p->DXP[IP];
     Vzsum+=voffield(i+1,j,k-1)*p->DXN[IP1]*p->DYN[JP]*p->DZN[KM1] * -p->DZP[KM1];
     // y*V=0
     
     //south-top
     Vsum+=voffield(i-1,j,k+1)*p->DXN[IM1]*p->DYN[JP]*p->DZN[KP1];
     Vxsum+=voffield(i-1,j,k+1)*p->DXN[IM1]*p->DYN[JP]*p->DZN[KP1] * -p->DXP[IM1];
     Vzsum+=voffield(i-1,j,k+1)*p->DXN[IM1]*p->DYN[JP]*p->DZN[KP1] *p->DZP[KP];
     // y*V=0
     
     //south-bottom
     Vsum+=voffield(i-1,j,k-1)*p->DXN[IM1]*p->DYN[JP]*p->DZN[KM1];
     Vxsum+=voffield(i-1,j,k-1)*p->DXN[IM1]*p->DYN[JP]*p->DZN[KM1] * -p->DXP[IM1];
     Vzsum+=voffield(i-1,j,k-1)*p->DXN[IM1]*p->DYN[JP]*p->DZN[KM1] * -p->DZP[KM1];
     // y*V=0
     
     if(p->j_dir>0)
     {
         //west
         Vsum+=voffield(i,j+1,k)*p->DXN[IP]*p->DYN[JP1]*p->DZN[KP];
         Vysum+=voffield(i,j+1,k)*p->DXN[IP]*p->DYN[JP1]*p->DZN[KP] *p->DYP[JP];
         // x*V=0, z*V=0
         
         //north-west
         Vsum+=voffield(i+1,j+1,k)*p->DXN[IP1]*p->DYN[JP1]*p->DZN[KP];
         Vxsum+=voffield(i+1,j+1,k)*p->DXN[IP1]*p->DYN[JP1]*p->DZN[KP] *p->DXP[IP];
         Vysum+=voffield(i+1,j+1,k)*p->DXN[IP1]*p->DYN[JP1]*p->DZN[KP] *p->DYP[JP];
         // z*V=0
         
         //south-west
         Vsum+=voffield(i-1,j+1,k)*p->DXN[IM1]*p->DYN[JP1]*p->DZN[KP];
         Vxsum+=voffield(i-1,j+1,k)*p->DXN[IM1]*p->DYN[JP1]*p->DZN[KP] * -p->DXP[IM1];
         Vysum+=voffield(i-1,j+1,k)*p->DXN[IM1]*p->DYN[JP1]*p->DZN[KP] *p->DYP[JP];
         // z*V=0
         
         //west-top
         Vsum+=voffield(i,j+1,k+1)*p->DXN[IP]*p->DYN[JP1]*p->DZN[KP1];
         Vysum+=voffield(i,j+1,k+1)*p->DXN[IP]*p->DYN[JP1]*p->DZN[KP1] *p->DYP[JP];
         Vzsum+=voffield(i,j+1,k+1)*p->DXN[IP]*p->DYN[JP1]*p->DZN[KP1] *p->DZP[KP];
         // x*V=0
         
         //west-bottom
         Vsum+=voffield(i,j+1,k-1)*p->DXN[IP]*p->DYN[JP1]*p->DZN[KM1];
         Vysum+=voffield(i,j+1,k-1)*p->DXN[IP]*p->DYN[JP1]*p->DZN[KM1] *p->DYP[JP];
         Vzsum+=voffield(i,j+1,k-1)*p->DXN[IP]*p->DYN[JP1]*p->DZN[KM1] * -p->DZP[KM1];
         // x*V=0
         
         //north-west-top
         Vsum+=voffield(i+1,j+1,k+1)*p->DXN[IP1]*p->DYN[JP1]*p->DZN[KP1];
         Vxsum+=voffield(i+1,j+1,k+1)*p->DXN[IP1]*p->DYN[JP1]*p->DZN[KP1] *p->DXP[IP];
         Vysum+=voffield(i+1,j+1,k+1)*p->DXN[IP1]*p->DYN[JP1]*p->DZN[KP1] *p->DYP[JP];
         Vzsum+=voffield(i+1,j+1,k+1)*p->DXN[IP1]*p->DYN[JP1]*p->DZN[KP1] *p->DZP[KP];
         
         //north-west-bottom
         Vsum+=voffield(i+1,j+1,k-1)*p->DXN[IP1]*p->DYN[JP1]*p->DZN[KM1];
         Vxsum+=voffield(i+1,j+1,k-1)*p->DXN[IP1]*p->DYN[JP1]*p->DZN[KM1] *p->DXP[IP];
         Vysum+=voffield(i+1,j+1,k-1)*p->DXN[IP1]*p->DYN[JP1]*p->DZN[KM1] *p->DYP[JP];
         Vzsum+=voffield(i+1,j+1,k-1)*p->DXN[IP1]*p->DYN[JP1]*p->DZN[KM1] * -p->DZP[KM1];
         
         //south-west-top
         Vsum+=voffield(i-1,j+1,k+1)*p->DXN[IM1]*p->DYN[JP1]*p->DZN[KP1];
         Vxsum+=voffield(i-1,j+1,k+1)*p->DXN[IM1]*p->DYN[JP1]*p->DZN[KP1] * -p->DXP[IM1];
         Vysum+=voffield(i-1,j+1,k+1)*p->DXN[IM1]*p->DYN[JP1]*p->DZN[KP1] *p->DYP[JP];
         Vzsum+=voffield(i-1,j+1,k+1)*p->DXN[IM1]*p->DYN[JP1]*p->DZN[KP1] *p->DZP[KP];
         
         //south-west-bottom
         Vsum+=voffield(i-1,j+1,k-1)*p->DXN[IM1]*p->DYN[JP1]*p->DZN[KM1];
         Vxsum+=voffield(i-1,j+1,k-1)*p->DXN[IM1]*p->DYN[JP1]*p->DZN[KM1] * -p->DXP[IM1];
         Vysum+=voffield(i-1,j+1,k-1)*p->DXN[IM1]*p->DYN[JP1]*p->DZN[KM1] *p->DYP[JP];
         Vzsum+=voffield(i-1,j+1,k-1)*p->DXN[IM1]*p->DYN[JP1]*p->DZN[KM1] * -p->DZP[KM1];
         
         //east
         Vsum+=voffield(i,j-1,k)*p->DXN[IP]*p->DYN[JM1]*p->DZN[KP];
         Vysum+=voffield(i,j-1,k)*p->DXN[IP]*p->DYN[JM1]*p->DZN[KP] * -p->DYP[JM1];
         // x*V=0, z*V=0
         
         //north-east
         Vsum+=voffield(i+1,j-1,k)*p->DXN[IP1]*p->DYN[JM1]*p->DZN[KP];
         Vxsum+=voffield(i+1,j-1,k)*p->DXN[IP1]*p->DYN[JM1]*p->DZN[KP] *p->DXP[IP];
         Vysum+=voffield(i+1,j-1,k)*p->DXN[IP1]*p->DYN[JM1]*p->DZN[KP] * -p->DYP[JM1];
         // z*V=0
         
         //south-east
         Vsum+=voffield(i-1,j-1,k)*p->DXN[IM1]*p->DYN[JM1]*p->DZN[KP];
         Vxsum+=voffield(i-1,j-1,k)*p->DXN[IM1]*p->DYN[JM1]*p->DZN[KP] * -p->DXP[IM1];
         Vysum+=voffield(i-1,j-1,k)*p->DXN[IM1]*p->DYN[JM1]*p->DZN[KP] * -p->DYP[JM1];
         // z*V=0
         
         //east-top
         Vsum+=voffield(i,j-1,k+1)*p->DXN[IP]*p->DYN[JM1]*p->DZN[KP1];
         Vysum+=voffield(i,j-1,k+1)*p->DXN[IP]*p->DYN[JM1]*p->DZN[KP1] * -p->DYP[JM1];
         Vzsum+=voffield(i,j-1,k+1)*p->DXN[IP]*p->DYN[JM1]*p->DZN[KP1] *p->DZP[KP1];
         // x*V=0
         
         //east-bottom
         Vsum+=voffield(i,j-1,k-1)*p->DXN[IP]*p->DYN[JM1]*p->DZN[KM1];
         Vysum+=voffield(i,j-1,k-1)*p->DXN[IP]*p->DYN[JM1]*p->DZN[KM1] * -p->DYP[JM1];
         Vzsum+=voffield(i,j-1,k-1)*p->DXN[IP]*p->DYN[JM1]*p->DZN[KM1] * -p->DZP[KM1];
         // x*V=0
         
         //north-east-top
         Vsum+=voffield(i+1,j-1,k+1)*p->DXN[IP1]*p->DYN[JM1]*p->DZN[KP1];
         Vxsum+=voffield(i+1,j-1,k+1)*p->DXN[IP1]*p->DYN[JM1]*p->DZN[KP1] *p->DXP[IP];
         Vysum+=voffield(i+1,j-1,k+1)*p->DXN[IP1]*p->DYN[JM1]*p->DZN[KP1] * -p->DYP[JM1];
         Vzsum+=voffield(i+1,j-1,k+1)*p->DXN[IP1]*p->DYN[JM1]*p->DZN[KP1] *p->DZP[KP];
         
         //north-east-bottom
         Vsum+=voffield(i+1,j-1,k-1)*p->DXN[IP1]*p->DYN[JM1]*p->DZN[KM1];
         Vxsum+=voffield(i+1,j-1,k-1)*p->DXN[IP1]*p->DYN[JM1]*p->DZN[KM1] *p->DXP[IP];
         Vysum+=voffield(i+1,j-1,k-1)*p->DXN[IP1]*p->DYN[JM1]*p->DZN[KM1] * -p->DYP[JM1];
         Vzsum+=voffield(i+1,j-1,k-1)*p->DXN[IP1]*p->DYN[JM1]*p->DZN[KM1] * -p->DZP[KM1];
         
         //south-east-top
         Vsum+=voffield(i-1,j-1,k+1)*p->DXN[IM1]*p->DYN[JM1]*p->DZN[KP1];
         Vxsum+=voffield(i-1,j-1,k+1)*p->DXN[IM1]*p->DYN[JM1]*p->DZN[KP1] * -p->DXP[IM1];
         Vysum+=voffield(i-1,j-1,k+1)*p->DXN[IM1]*p->DYN[JM1]*p->DZN[KP1] * -p->DYP[JM1];
         Vzsum+=voffield(i-1,j-1,k+1)*p->DXN[IM1]*p->DYN[JM1]*p->DZN[KP1] *p->DZP[KP];
         
         //south-east-bottom
         Vsum+=voffield(i-1,j-1,k-1)*p->DXN[IM1]*p->DYN[JM1]*p->DZN[KM1];
         Vxsum+=voffield(i-1,j-1,k-1)*p->DXN[IM1]*p->DYN[JM1]*p->DZN[KM1] * -p->DXP[IM1];
         Vysum+=voffield(i-1,j-1,k-1)*p->DXN[IM1]*p->DYN[JM1]*p->DZN[KM1] * -p->DYP[JM1];
         Vzsum+=voffield(i-1,j-1,k-1)*p->DXN[IM1]*p->DYN[JM1]*p->DZN[KM1] * -p->DZP[KM1];
     }
     
     if(Vsum>1E-16)
        invvec_x=Vxsum/Vsum;
     else
         invvec_x=0.0;
         
     if(Vsum>1E-16)
        invvec_y=Vysum/Vsum;
     else
        invvec_y=0.0;
        
     if(Vsum>1E-16)
        invvec_z=Vzsum/Vsum;
     else
        invvec_z=0.0;
    
     nsum=sqrt(invvec_x*invvec_x+invvec_y*invvec_y+invvec_z*invvec_z);
     if(nsum>1E-16)
     {
        nx(i,j,k)=-invvec_x/nsum;
        ny(i,j,k)=-invvec_y/nsum;
        nz(i,j,k)=-invvec_z/nsum;
     }
     else
     {
         nx(i,j,k)=0.0;
         ny(i,j,k)=0.0;
         nz(i,j,k)=1.0;
         cout<<"nsum=0"<<endl;
     }
     
     //cout<<"nx: "<<nx(i,j,k)<<endl;
    // cout<<"ny: "<<ny(i,j,k)<<endl;
    // cout<<"nz: "<<nz(i,j,k)<<endl;
}

void VOF_PLIC::calcNormalPhi(fdm* a, lexer* p)
{
    double nsum;
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
        
    nsum=sqrt(nx(i,j,k)*nx(i,j,k)+ny(i,j,k)*ny(i,j,k)+nz(i,j,k)*nz(i,j,k));
    nx(i,j,k)=nx(i,j,k)/nsum;
    ny(i,j,k)=ny(i,j,k)/nsum;
    nz(i,j,k)=nz(i,j,k)/nsum;
}


void VOF_PLIC:: calcNormalWeymouth(fdm* a, lexer* p, field4 voffield)
{
    int dimswitch;
    int baseswitch;
    double n_max=0.0;
    double vof_max;
    double nsum;
    double n1m, n1c, n1p, n2m, n2c, n2p, n1_max, n2_max;
    double p1m, pcc, p1p, p2m, p2p;
    double n_x,n_y,n_z;
    
    //z-dir different schemes -> n_z=1.0)
    if(voffield(i,j,k+1)>voffield(i,j,k-1))
    {
        baseswitch=3;
        vof_max=voffield(i,j,k+1);
    }
    else if(voffield(i,j,k-1)>voffield(i,j,k-1))
    {
        baseswitch=-3;
        vof_max=voffield(i,j,k-1);
    }
    else
    {
        baseswitch=30;
        vof_max=voffield(i,j,k+1);
    }
    
    pcc=voffield(i,j,k-1)+voffield(i,j,k)+voffield(i,j,k+1);
    p1m=voffield(i-1,j,k-1)+voffield(i-1,j,k)+voffield(i-1,j,k+1);
    p1p=voffield(i+1,j,k-1)+voffield(i+1,j,k)+voffield(i+1,j,k+1);
    n1m=-(pcc-p1m)/(p->DXP[IM1]);
    n1c=-(p1p-p1m)/(p->DXP[IM1]+p->DXP[IP]);
    n1p=-(p1p-pcc)/(p->DXP[IP]);
    n1_max=max({fabs(n1m),fabs(n1c),fabs(n1p)});
    
    if(p->j_dir>0)
    {   
        p2m=voffield(i,j-1,k-1)+voffield(i,j-1,k)+voffield(i,j-1,k+1);
        p2p=voffield(i,j+1,k-1)+voffield(i,j+1,k)+voffield(i,j+1,k+1);
        n2m=-(pcc-p2m)/(p->DYP[JM1]);
        n2c=-(p2p-p2m)/(p->DYP[JM1]+p->DYP[JP]);
        n2p=-(p2p-pcc)/(p->DYP[JP]);
        n2_max=max({fabs(n2m),fabs(n2c),fabs(n2p)});
    }
    else
        n2_max=0.0;
    
    if(n1_max>=n2_max)
    {
        n_max=n1_max;
        dimswitch=1;
    }
    else
    {
        n_max=n2_max;
        dimswitch=2;
    }
    
    //x-dir different schemes -> n_x=1.0)
    if(voffield(i+1,j,k)>vof_max)
    {
        vof_max=voffield(i+1,j,k);
        baseswitch=1;
    }
    if(voffield(i-1,j,k)>vof_max)
    {
        vof_max=voffield(i-1,j,k);
        baseswitch=-1;
    }
    if(voffield(i+1,j,k)>voffield(i-1,j,k)-1E-06 && voffield(i+1,j,k)<voffield(i-1,j,k)+1E-06 && (baseswitch==1 || baseswitch==-1))
    {
        vof_max=voffield(i+1,j,k);
        baseswitch=10;
    }
    
    pcc=voffield(i-1,j,k)+voffield(i,j,k)+voffield(i+1,j,k);
    p1m=voffield(i-1,j,k-1)+voffield(i,j,k-1)+voffield(i+1,j,k-1);
    p1p=voffield(i-1,j,k+1)+voffield(i,j,k+1)+voffield(i+1,j,k+1);
    n1m=-(pcc-p1m)/(p->DZP[KM1]);
    n1c=-(p1p-p1m)/(p->DZP[KM1]+p->DZP[KP]);
    n1p=-(p1p-pcc)/(p->DZP[KP]);
    n1_max=max({fabs(n1m),fabs(n1c),fabs(n1p)});
    
    if(p->j_dir>0)
    {
        p2m=voffield(i-1,j-1,k)+voffield(i,j-1,k)+voffield(i+1,j-1,k);
        p2p=voffield(i-1,j+1,k)+voffield(i,j+1,k)+voffield(i+1,j+1,k);
        n2m=-(pcc-p2m)/(p->DYP[JM1]);
        n2c=-(p2p-p2m)/(p->DYP[JM1]+p->DYP[JP]);
        n2p=-(p2p-pcc)/(p->DYP[JP]);
        n2_max=max({fabs(n2m),fabs(n2c),fabs(n2p)});
    }
    else
        n2_max=0.0;
        
    if(n1_max>n_max)
    {
        n_max=n1_max;
        dimswitch=3;
    }
    
    if(n2_max>n_max)
    {
        n_max=n2_max;
        dimswitch=2;
    }
    
    if(p->j_dir>0)
    {   
        if(voffield(i,j+1,k)>vof_max)
        {
            vof_max=voffield(i,j+1,k);
            baseswitch=2;
        }
        if(voffield(i,j-1,k)>vof_max)
        {
            vof_max=voffield(i,j-1,k);
            baseswitch=-2;
        }
        if(voffield(i,j+1,k)>voffield(i,j-1,k)-1E-06 && voffield(i,j+1,k)<voffield(i,j-1,k)+1E-06 && (baseswitch==2 || baseswitch==-2))
        {
            vof_max=voffield(i,j+1,k);
            baseswitch=20;
        }
        
        pcc=voffield(i,j-1,k)+voffield(i,j,k)+voffield(i,j+1,k);
        p1m=voffield(i-1,j-1,k)+voffield(i-1,j,k)+voffield(i-1,j+1,k);
        p1p=voffield(i+1,j-1,k)+voffield(i+1,j,k)+voffield(i+1,j+1,k);
        p2m=voffield(i,j-1,k-1)+voffield(i,j,k-1)+voffield(i,j+1,k-1);
        p2p=voffield(i,j-1,k+1)+voffield(i,j,k+1)+voffield(i,j+1,k+1);
        n1m=-(pcc-p1m)/(p->DXP[IM1]);
        n1c=-(p1p-p1m)/(p->DXP[IM1]+p->DXP[IP]);
        n1p=-(p1p-pcc)/(p->DXP[IP]);
        n2m=-(pcc-p2m)/(p->DZP[KM1]);
        n2c=-(p2p-p2m)/(p->DZP[KM1]+p->DZP[KP]);
        n2p=-(p2p-pcc)/(p->DZP[KP]);
        n1_max=max({fabs(n1m),fabs(n1c),fabs(n1p)});
        n2_max=max({fabs(n2m),fabs(n2c),fabs(n2p)});
        
        if(n1_max>n_max)
        {
            n_max=n1_max;
            dimswitch=1;
        }
        if(n2_max>n_max)
        {
            n_max=n2_max;
            dimswitch=3;
        }
    }
    
    //if(baseswitch==10 || baseswitch==20 || baseswitch==30)
  //  {
        if(dimswitch==1)
        {
            if(voffield(i-1,j,k)>=voffield(i+1,j,k))
                baseswitch=-1;
            else
                baseswitch=1;
        }
        else if(dimswitch==2)
        {
            if(voffield(i,j-1,k)>=voffield(i,j+1,k))
                baseswitch=-2;
            else
                baseswitch=2;
        }
        else if(dimswitch==3)
        {
            if(voffield(i,j,k-1)>=voffield(i,j,k+1))
                baseswitch=-3;
            else
                baseswitch=3;
        }
        else
            cout<<"dimswitch fail"<<endl;
        
  //  }
    
    
    if(baseswitch==-1 || baseswitch==1)
    {
        n_x=1.0;
        pcc=voffield(i-1,j,k)+voffield(i,j,k)+voffield(i+1,j,k);
        p1m=voffield(i-1,j,k-1)+voffield(i,j,k-1)+voffield(i+1,j,k-1);
        p1p=voffield(i-1,j,k+1)+voffield(i,j,k+1)+voffield(i+1,j,k+1);
        /*if(voffield(i,j,k)>=0.5)
        {
            if(p1p>=p1m)
                n_z=-(p1p-pcc)/(p->DZP[KP]);
            else
                n_z=-(pcc-p1m)/(p->DZP[KM1]);
        }
        else
        {
            if(p1m<=p1p)
                n_z=-(pcc-p1m)/(p->DZP[KM1]);
            else
                n_z=-(p1p-pcc)/(p->DZP[KP]);
        }*/
        n_z=-(p1p-p1m)/(p->DZP[KP]+p->DZP[KM1]);
        
        if(p->j_dir>0)
        {   
            p2m=voffield(i-1,j-1,k)+voffield(i,j-1,k)+voffield(i+1,j-1,k);
            p2p=voffield(i-1,j+1,k)+voffield(i,j+1,k)+voffield(i+1,j+1,k);
           /* if(voffield(i,j,k)>=0.5)
            {
                if(p2p>=p2m)
                    n_y=-(p2p-pcc)/(p->DYP[JP]);
                else
                    n_y=-(pcc-p2m)/(p->DYP[JM1]);
            }
            else
            {
                if(p2m<=p2p)
                    n_y=-(pcc-p2m)/(p->DYP[JM1]);
                else
                    n_y=-(p2p-pcc)/(p->DYP[JP]);
            }*/
            n_y=-(p2p-p2m)/(p->DYP[JP]+p->DYP[IM1]);
        }
        else
            n_y=0.0;
        
        if(baseswitch==1)
        {
            n_z=-n_z;
            n_y=-n_y;
            n_x=-n_x;
        }
    }
    else if(baseswitch==3 || baseswitch==-3)
    {
        n_z=1.0;
        pcc=voffield(i,j,k-1)+voffield(i,j,k)+voffield(i,j,k+1);
        p1m=voffield(i-1,j,k-1)+voffield(i-1,j,k)+voffield(i-1,j,k+1);
        p1p=voffield(i+1,j,k-1)+voffield(i+1,j,k)+voffield(i+1,j,k+1);
       /* if(voffield(i,j,k)>=0.5)
        {
            if(p1p>=p1m)
                n_x=-(p1p-pcc)/(p->DXP[IP]);
            else
                n_x=-(pcc-p1m)/(p->DXP[IM1]);
        }
        else
        {
            if(p1m<=p1p)
                n_x=-(pcc-p1m)/(p->DXP[IM1]);
            else
                n_x=-(p1p-pcc)/(p->DXP[IP]);
        }
        */
        n_x=-(p1p-p1m)/(p->DXP[IP]+p->DXP[IM1]);
        
        if(p->j_dir>0)
        {
            p2m=voffield(i,j-1,k-1)+voffield(i,j-1,k)+voffield(i,j-1,k+1);
            p2p=voffield(i,j+1,k-1)+voffield(i,j+1,k)+voffield(i,j+1,k+1);
            /*if(voffield(i,j,k)>=0.5)
            {
                if(p2p>=p2m)
                    n_y=-(p2p-pcc)/(p->DYP[JP]);
                else
                    n_y=-(pcc-p2m)/(p->DYP[JM1]);
            }
            else
            {
                if(p2m<=p2p)
                    n_y=-(pcc-p2m)/(p->DYP[JM1]);
                else
                    n_y=-(p2p-pcc)/(p->DYP[JP]);
            }*/
            n_y=-(p2p-p2m)/(p->DYP[JP]+p->DYP[JM1]);
        }
        else
            n_y=0.0;
            
        if(baseswitch==3)
        {
            n_z=-n_z;
            n_y=-n_y;
            n_x=-n_x;
        }
        
    }
    else if(baseswitch==2 || baseswitch==-2)
    {
        n_y=1.0;
        pcc=voffield(i,j-1,k)+voffield(i,j,k)+voffield(i,j+1,k);
        p1m=voffield(i-1,j-1,k)+voffield(i-1,j,k)+voffield(i-1,j+1,k);
        p1p=voffield(i+1,j-1,k)+voffield(i+1,j,k)+voffield(i+1,j+1,k);
        p2m=voffield(i,j-1,k-1)+voffield(i,j,k-1)+voffield(i,j+1,k-1);
        p2p=voffield(i,j-1,k+1)+voffield(i,j,k+1)+voffield(i,j+1,k+1);
        /*if(voffield(i,j,k)>=0.5)
        {
            if(p1p>=p1m)
                n_x=-(p1p-pcc)/(p->DXP[IP]);
            else
                n_x=-(pcc-p1m)/(p->DXP[IM1]);
                
            if(p2p>=p2m)
                n_z=-(p2p-pcc)/(p->DZP[KP]);
            else
                n_z=-(pcc-p2m)/(p->DZP[KM1]);
        }
        
        n_z=-(p1p-p1m)/(p->DZP[KP]+p->DZP[KM1]);
        else
        {
            if(p1m<=p1p)
                n_x=-(pcc-p1m)/(p->DXP[IM1]);
            else
                n_x=-(p1p-pcc)/(p->DXP[IP]);
                
            if(p2m<=p2p)
                n_z=-(pcc-p2m)/(p->DZP[KM1]);
            else
                n_z=-(p2p-pcc)/(p->DZP[KP]);
        }*/
        
        n_x=-(p1p-p1m)/(p->DXP[IP]+p->DXP[IM1]);
        n_z=-(p2p-p2m)/(p->DZP[KP]+p->DZP[KM1]);
        
        
        if(baseswitch==2)
        {
            n_y=-n_y;
            n_x=-n_x;
            n_z=-n_z;
        }
    }
    
    nsum=sqrt(n_x*n_x+n_y*n_y+n_z*n_z);
    nx(i,j,k)=n_x/nsum;
    ny(i,j,k)=n_y/nsum;
    nz(i,j,k)=n_z/nsum;
    
    
    
}

void VOF_PLIC::calcNormalWang(fdm* a, lexer* p)
{
    double mx,my,mz;
    double mx_save,my_save,mz_save;
    double msum, nsum;
    double mincheck=1E06;
    
    for(int iloop=0; iloop<2;iloop++)
    {
        for(int jloop=0; jloop<2; jloop++)
        {
            for(int kloop=0; kloop<2; kloop++)
            {
                switch(iloop)
                {
                    case 0:
                        mx=(phistep(i,j,k)-phistep(i-1,j,k))/p->DXP[IM1];
                        break;
                    
                    case 1:
                        mx=(phistep(i+1,j,k)-phistep(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1]);
                        break;
                        
                    case 2:
                        mx=(phistep(i+1,j,k)-phistep(i,j,k))/p->DXP[IP];
                        break;
                }
                
                if(p->j_dir>0)
                {
                    switch(jloop)
                    {
                        case 0:
                            my=(phistep(i,j,k)-phistep(i,j-1,k))/p->DYP[JM1];
                            break;
                        
                        case 1:
                            my=(phistep(i,j+1,k)-phistep(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1]);
                            break;
                        
                        case 2:
                            my=(phistep(i,j+1,k)-phistep(i,j,k))/p->DYP[JP];
                            break;
                    }
                }
                else
                    my=0.0;
                
                switch(kloop)
                {
                    case 0:
                        mz=(phistep(i,j,k)-phistep(i,j,k-1))/p->DZP[KM1];
                        break;
                        
                    case 1:
                        mz=(phistep(i,j,k+1)-phistep(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]);
                        break;
                        
                    case 2:
                        mz=(phistep(i,j,k+1)-phistep(i,j,k))/p->DZP[KP];
                        break;
                }
                
                msum=sqrt(mx*mx+my*my+mz*mz);
                if(fabs(msum-1.0)<mincheck)
                {
                    mincheck=fabs(msum-1.0);
                    mx_save=mx;
                    my_save=my;
                    mz_save=mz;
                }
                
                
            }
        }
    }
    
    nsum=sqrt(mx_save*mx_save+my_save*my_save+mz_save*mz_save);
    nx(i,j,k)=-mx_save/nsum;
    ny(i,j,k)=-my_save/nsum;
    nz(i,j,k)=-mz_save/nsum;
}