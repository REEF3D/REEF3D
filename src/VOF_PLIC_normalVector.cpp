/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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


void VOF_PLIC::calcNormalFO(fdm* a, lexer* p, field& voffield)
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

void VOF_PLIC::calcNormalLS(fdm* a, lexer* p, field& voffield)
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


void VOF_PLIC::calcNormalWENO(fdm* a, lexer* p, field& voffield)
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

void VOF_PLIC::calcNormalMassCentre(fdm* a, lexer* p, field& voffield)
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


void VOF_PLIC:: calcNormalWeymouth(fdm* a, lexer* p, field& voffield)
{
    
}

void VOF_PLIC::calcNormalWang(fdm* a, lexer* p)
{

}

/*
void VOF_PLIC::calcNormalELVIRA2D(fdm* a, lexer* p, field& voffield)
{
    // Initialize arrays for 6 candidate normals
    double nx_candidates[6], nz_candidates[6];
    double min_error = 1e6;
    double nx_best = 0.0, nz_best = 0.0;
    double alpha_best = 0.0;
    // Column sums (3x3 stencil)
    double zsum_xm = voffield(i-1,j,k-1) + voffield(i-1,j,k) + voffield(i-1,j,k+1); // Left column
    double zsum_xc = voffield(i,j,k-1) + voffield(i,j,k) + voffield(i,j,k+1);       // Center column
    double zsum_xp = voffield(i+1,j,k-1) + voffield(i+1,j,k) + voffield(i+1,j,k+1); // Right column
    // Row sums
    double xsum_zm = voffield(i-1,j,k-1) + voffield(i,j,k-1) + voffield(i+1,j,k-1); // Bottom row
    double xsum_zc = voffield(i-1,j,k) + voffield(i,j,k) + voffield(i+1,j,k);       // Middle row
    double xsum_zp = voffield(i-1,j,k+1) + voffield(i,j,k+1) + voffield(i+1,j,k+1); // Top row
    // Compute candidate normals using central differences
    // x-direction candidates
    nx_candidates[0] = 1.0;
    nz_candidates[0] = (zsum_xp - zsum_xm)/(2.0*p->DXP[IP]); // Scale by dx for x-sweep
    nx_candidates[1] = -1.0;
    nz_candidates[1] = -(zsum_xp - zsum_xm)/(2.0*p->DXP[IP]);
    // z-direction candidates (vertical interface)
    nx_candidates[2] = (xsum_zp - xsum_zm)/(2.0*p->DZP[KP]); // Scale by dz for z-sweep
    nz_candidates[2] = 1.0;
    nx_candidates[3] = -(xsum_zp - xsum_zm)/(2.0*p->DZP[KP]);
    nz_candidates[3] = -1.0;
    // Mixed candidates (diagonal interfaces)
    nx_candidates[4] = 1.0;
    nz_candidates[4] = (voffield(i+1,j,k+1) - voffield(i-1,j,k+1)
                    - voffield(i+1,j,k-1) + voffield(i-1,j,k-1))
                    /(4.0*p->DXP[IP]); // Cross-derivative dx*dz
    nx_candidates[5] = -1.0;
    nz_candidates[5] = -(voffield(i+1,j,k+1) - voffield(i-1,j,k+1)
                        - voffield(i+1,j,k-1) + voffield(i-1,j,k-1))
                        /(4.0*p->DXP[IP]);
    // Test each candidate normal
    for(int m=0; m<6; ++m)
    {
        // Normalize candidate normal
        double mag = sqrt(nx_candidates[m]*nx_candidates[m] + nz_candidates[m]*nz_candidates[m]);
        nx_candidates[m] /= mag;
        nz_candidates[m] /= mag;
        // Calculate optimal alpha for this normal
        double alpha = calcAlphaFromInput(a, p,
                                        nx_candidates[m], 0.0, nz_candidates[m],
                                        p->DXP[IP], p->DYP[JP], p->DZP[KP],
                                        voffield(i,j,k));
        // Calculate reconstruction error with proper alpha
        double error = calcL2vofError2D(a, p, voffield,
                                      nx_candidates[m], 0.0, nz_candidates[m], alpha);
        // Store best normal and its alpha
        if(error < min_error)
        {
            min_error = error;
            nx_best = nx_candidates[m];
            nz_best = nz_candidates[m];
            alpha_best = alpha;
        }
    }
    // Set interface normal and alpha to best candidate
    nx(i,j,k) = nx_best;
    ny(i,j,k) = 0.0;
    nz(i,j,k) = nz_best;
    alpha(i,j,k) = alpha_best;
}
*/
void VOF_PLIC::calcNormalELVIRA2D(fdm* a, lexer* p, field& voffield)
{
    double n_y=0.0;
    double n_x,n_z,L2,nsum,r0,r0mod,recheck,voflocal;
    double nx_keep,ny_keep,nz_keep,r0_keep;
    double minL2=1E06;
    double zsum_xm,zsum_xc,zsum_xp;
    double xsum_zm,xsum_zc,xsum_zp;
    
    zsum_xm=voffield(i-1,j,k-1)*p->DZN[KM1]+voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k+1)*p->DZN[KP1];
    zsum_xc=voffield(i,j,k-1)*p->DZN[KM1]+voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k+1)*p->DZN[KP1];
    zsum_xp=voffield(i+1,j,k-1)*p->DZN[KM1]+voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k+1)*p->DZN[KP1];
    
    xsum_zm=voffield(i-1,j,k-1)*p->DXN[IM1]+voffield(i,j,k-1)*p->DXN[IP]+voffield(i+1,j,k-1)*p->DXN[IP1];
    xsum_zc=voffield(i-1,j,k)*p->DXN[IM1]+voffield(i,j,k)*p->DXN[IP]+voffield(i+1,j,k)*p->DXN[IP1];
    xsum_zp=voffield(i-1,j,k+1)*p->DXN[IM1]+voffield(i,j,k+1)*p->DXN[IP]+voffield(i+1,j,k+1)*p->DXN[IP1];
    
    //downwind in x
    n_z=1.0;
    n_x=-(zsum_xc-zsum_xm)/p->DXP[IM1];
    nsum=sqrt(n_x*n_x+n_z*n_z);
    n_z=n_z/nsum;
    n_x=n_x/nsum;
    r0=calcAlphaFromInput(a,p,n_x,n_y,n_z,p->DXN[IP],p->DYN[JP],p->DZN[KP],voffield(i,j,k));
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    //normvecdirtest
    n_z=-n_z;
    n_x=-n_x;
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    
    //central in x
    n_z=1.0;
    n_x=-(zsum_xp-zsum_xm)/(p->DXP[IM1]+p->DXP[IP]);
    nsum=sqrt(n_x*n_x+n_z*n_z);
    n_z=n_z/nsum;
    n_x=n_x/nsum;
    r0=calcAlphaFromInput(a,p,n_x,n_y,n_z,p->DXN[IP],p->DYN[JP],p->DZN[KP],voffield(i,j,k));
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    //normvecdirtest
    n_z=-n_z;
    n_x=-n_x;
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    
    //upwind in x
    n_z=1.0;
    n_x=-(zsum_xp-zsum_xc)/p->DXP[IP];
    nsum=sqrt(n_x*n_x+n_z*n_z);
    n_z=n_z/nsum;
    n_x=n_x/nsum;
    r0=calcAlphaFromInput(a,p,n_x,n_y,n_z,p->DXN[IP],p->DYN[JP],p->DZN[KP],voffield(i,j,k));
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    //normvecdirtest
    n_z=-n_z;
    n_x=-n_x;
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    
    //downwind in z
    n_x=1.0;
    n_z=-(xsum_zc-xsum_zm)/p->DZP[KM1];
    nsum=sqrt(n_x*n_x+n_z*n_z);
    n_z=n_z/nsum;
    n_x=n_x/nsum;
    r0=calcAlphaFromInput(a,p,n_x,n_y,n_z,p->DXN[IP],p->DYN[JP],p->DZN[KP],voffield(i,j,k));
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    //normvecdirtest
    n_z=-n_z;
    n_x=-n_x;
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    
    //central in z
    n_x=1.0;
    n_z=-(xsum_zp-xsum_zm)/(p->DZP[KM1]+p->DZP[KP]);
    nsum=sqrt(n_x*n_x+n_z*n_z);
    n_z=n_z/nsum;
    n_x=n_x/nsum;
    r0=calcAlphaFromInput(a,p,n_x,n_y,n_z,p->DXN[IP],p->DYN[JP],p->DZN[KP],voffield(i,j,k));
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    //normvecdirtest
    n_z=-n_z;
    n_x=-n_x;
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    
    //upwind in z
    n_x=1.0;
    n_z=-(xsum_zp-xsum_zc)/p->DZP[KP];
    nsum=sqrt(n_x*n_x+n_z*n_z);
    n_z=n_z/nsum;
    n_x=n_x/nsum;
    r0=calcAlphaFromInput(a,p,n_x,n_y,n_z,p->DXN[IP],p->DYN[JP],p->DZN[KP],voffield(i,j,k));
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    //normvecdirtest
    n_z=-n_z;
    n_x=-n_x;
    L2=calcL2vofError2D(a,p,voffield,n_x,n_y,n_z,r0);
    if(L2<minL2)
    {
        minL2=L2;
        nx_keep=n_x;
        ny_keep=n_y;
        nz_keep=n_z;
        r0_keep=r0;
    }
    
    nx(i,j,k)=nx_keep;
    ny(i,j,k)=ny_keep;
    nz(i,j,k)=nz_keep;
    alpha(i,j,k)=r0_keep;
}

double VOF_PLIC::calcL2vofError2D(fdm* a, lexer* p, field& voffield, double n_x, double n_y, double n_z, double r0)
{
    double ret=0.0;
    double r0mod,recheck,voflocal;
    
    //north
    r0mod=-(n_x*p->DXP[IP]-r0);
    recheck=0.5*(fabs(n_x)*p->DXN[IP1]+fabs(n_y)*p->DYN[JP]+fabs(n_z)*p->DZN[KP])-fabs(r0mod);
    if(recheck>0.0)
        voflocal=calculateVolume(n_x,n_y,n_z,p->DXN[IP1],p->DYN[JP],p->DZN[KP],r0mod);
    else if(r0mod>0.0)
        voflocal=1.0;
    else
        voflocal=0.0;
    ret+=(voffield(i+1,j,k)-voflocal)*(voffield(i+1,j,k)-voflocal);
    
    //south
    r0mod=-(n_x*(-p->DXP[IM1])-r0);
    recheck=0.5*(fabs(n_x)*p->DXN[IM1]+fabs(n_y)*p->DYN[JP]+fabs(n_z)*p->DZN[KP])-fabs(r0mod);
    if(recheck>0.0)
        voflocal=calculateVolume(n_x,n_y,n_z,p->DXN[IM1],p->DYN[JP],p->DZN[KP],r0mod);
    else if(r0mod>0.0)
        voflocal=1.0;
    else
        voflocal=0.0;
    ret+=(voffield(i-1,j,k)-voflocal)*(voffield(i-1,j,k)-voflocal);
    
    //top
    r0mod=-(n_z*p->DZP[IP]-r0);
    recheck=0.5*(fabs(n_x)*p->DXN[IP]+fabs(n_y)*p->DYN[JP]+fabs(n_z)*p->DZN[KP1])-fabs(r0mod);
    if(recheck>0.0)
        voflocal=calculateVolume(n_x,n_y,n_z,p->DXN[IP],p->DYN[JP],p->DZN[KP1],r0mod);
    else if(r0mod>0.0)
        voflocal=1.0;
    else
        voflocal=0.0;
    ret+=(voffield(i,j,k+1)-voflocal)*(voffield(i,j,k+1)-voflocal);
    
    //bottom
    r0mod=-(n_z*(-p->DZP[IM1])-r0);
    recheck=0.5*(fabs(n_x)*p->DXN[IP]+fabs(n_y)*p->DYN[JP]+fabs(n_z)*p->DZN[KM1])-fabs(r0mod);
    if(recheck>0.0)
        voflocal=calculateVolume(n_x,n_y,n_z,p->DXN[IP],p->DYN[JP],p->DZN[KM1],r0mod);
    else if(r0mod>0.0)
        voflocal=1.0;
    else
        voflocal=0.0;
    ret+=(voffield(i,j,k-1)-voflocal)*(voffield(i,j,k-1)-voflocal);
    
    //north-top
    r0mod=-(n_x*p->DXP[IP]+n_z*p->DZP[KP]-r0);
    recheck=0.5*(fabs(n_x)*p->DXN[IP1]+fabs(n_y)*p->DYN[JP]+fabs(n_z)*p->DZN[KP1])-fabs(r0mod);
    if(recheck>0.0)
        voflocal=calculateVolume(n_x,n_y,n_z,p->DXN[IP1],p->DYN[JP],p->DZN[KP1],r0mod);
    else if(r0mod>0.0)
        voflocal=1.0;
    else
        voflocal=0.0;
    ret+=(voffield(i+1,j,k+1)-voflocal)*(voffield(i+1,j,k+1)-voflocal);
    
    //south-top
    r0mod=-(n_x*(-p->DXP[IM1])+n_z*p->DZP[KP]-r0);
    recheck=0.5*(fabs(n_x)*p->DXN[IM1]+fabs(n_y)*p->DYN[JP]+fabs(n_z)*p->DZN[KP1])-fabs(r0mod);
    if(recheck>0.0)
        voflocal=calculateVolume(n_x,n_y,n_z,p->DXN[IM1],p->DYN[JP],p->DZN[KP1],r0mod);
    else if(r0mod>0.0)
        voflocal=1.0;
    else
        voflocal=0.0;
    ret+=(voffield(i-1,j,k+1)-voflocal)*(voffield(i-1,j,k+1)-voflocal);
    
    //north-bottom
    r0mod=-(n_x*p->DXP[IP]+n_z*(-p->DZP[KM1])-r0);
    recheck=0.5*(fabs(n_x)*p->DXN[IP1]+fabs(n_y)*p->DYN[JP]+fabs(n_z)*p->DZN[KM1])-fabs(r0mod);
    if(recheck>0.0)
        voflocal=calculateVolume(n_x,n_y,n_z,p->DXN[IP1],p->DYN[JP],p->DZN[KM1],r0mod);
    else if(r0mod>0.0)
        voflocal=1.0;
    else
        voflocal=0.0;
    ret+=(voffield(i+1,j,k-1)-voflocal)*(voffield(i+1,j,k-1)-voflocal);
    
    //south-bottom
    r0mod=-(n_x*(-p->DXP[IM1])+n_z*(-p->DZP[KM1])-r0);
    recheck=0.5*(fabs(n_x)*p->DXN[IM1]+fabs(n_y)*p->DYN[JP]+fabs(n_z)*p->DZN[KM1])-fabs(r0mod);
    if(recheck>0.0)
        voflocal=calculateVolume(n_x,n_y,n_z,p->DXN[IM1],p->DYN[JP],p->DZN[KM1],r0mod);
    else if(r0mod>0.0)
        voflocal=1.0;
    else
        voflocal=0.0;
    ret+=(voffield(i-1,j,k-1)-voflocal)*(voffield(i-1,j,k-1)-voflocal);
    ret=sqrt(ret);
    
    return ret;
}

void VOF_PLIC::calcNormalMYC3D(fdm* a,lexer* p, field& voffield)
{
    //First two Candidates by Centred Columns Scheme
    
    double zsum_xm, zsum_xp, zsum_yp, zsum_ym; 
    double xsum_zm, xsum_zp, xsum_yp, xsum_ym;
    double ysum_xm, ysum_xp, ysum_zm, ysum_zp;
    double vofsumup,vofsumdown,sign;
    double nx_Cz, ny_Cz, nz_Cz, nx_Cx, ny_Cx, nz_Cx, nx_Cy, ny_Cy, nz_Cy, nx_CC, nz_CC;
    double nsum;
    
    //Candidate CC1 height function is z dimension
    zsum_xm=voffield(i-1,j,k-1)*p->DZN[KM1]+voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k+1)*p->DZN[KP1];
    zsum_xp=voffield(i+1,j,k-1)*p->DZN[KM1]+voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k+1)*p->DZN[KP1];
    zsum_ym=voffield(i,j-1,k-1)*p->DZN[KM1]+voffield(i,j-1,k)*p->DZN[KP]+voffield(i,j-1,k+1)*p->DZN[KP1];
    zsum_yp=voffield(i,j+1,k-1)*p->DZN[KM1]+voffield(i,j+1,k)*p->DZN[KP]+voffield(i,j+1,k+1)*p->DZN[KP1];
    vofsumup=voffield(i-1,j,k+1)*p->DXN[IM1]*p->DYN[JP]+voffield(i,j,k+1)*p->DXN[IP]*p->DYN[JP]+voffield(i+1,j,k+1)*p->DXN[IP1]*p->DYN[JP]
            +voffield(i,j-1,k+1)*p->DXN[IP]*p->DYN[JM1]+voffield(i,j+1,k+1)*p->DXN[IP]*p->DYN[JP1];
    vofsumdown=voffield(i-1,j,k-1)*p->DXN[IM1]*p->DYN[JP]+voffield(i,j,k-1)*p->DXN[IP]*p->DYN[JP]+voffield(i+1,j,k-1)*p->DXN[IP1]*p->DYN[JP]
            +voffield(i,j-1,k-1)*p->DXN[IP]*p->DYN[JM1]+voffield(i,j+1,k-1)*p->DXN[IP]*p->DYN[JP1];
    
    if(vofsumup>vofsumdown)
        sign=-1.0;
    else
        sign=1.0;
    
    nz_Cz=1.0*sign;
    nx_Cz=-(zsum_xp-zsum_xm)/(p->DXP[IM1]+p->DXP[IP]);
    ny_Cz=-(zsum_yp-zsum_ym)/(p->DYP[JM1]+p->DYP[JP]);
    nsum=sqrt(nx_Cz*nx_Cz+ny_Cz*ny_Cz+nz_Cz*nz_Cz);
    if(nsum!=nsum)
        cout<<"nsumNAN Cz"<<endl;
    nz_Cz=nz_Cz/nsum;
    nx_Cz=nx_Cz/nsum;
    ny_Cz=nz_Cz/nsum;
    
    
    //Candidate CC2 height function is x dimension
    xsum_zm=voffield(i-1,j,k-1)*p->DXN[IM1]+voffield(i,j,k-1)*p->DXN[IP]+voffield(i+1,j,k-1)*p->DXN[IP1];
    xsum_zp=voffield(i-1,j,k+1)*p->DXN[IM1]+voffield(i,j,k+1)*p->DXN[IP]+voffield(i+1,j,k+1)*p->DXN[IP1];
    xsum_yp=voffield(i-1,j+1,k)*p->DXN[IM1]+voffield(i,j+1,k)*p->DXN[IP]+voffield(i+1,j+1,k)*p->DXN[IP1];
    xsum_ym=voffield(i-1,j-1,k)*p->DXN[IM1]+voffield(i,j-1,k)*p->DXN[IP]+voffield(i+1,j-1,k)*p->DXN[IP1];
    vofsumup=voffield(i+1,j,k-1)*p->DZN[KM1]*p->DYN[JP]+voffield(i+1,j,k)*p->DZN[KP]*p->DYN[JP]+voffield(i+1,j,k+1)*p->DZN[KP1]*p->DYN[JP]
            +voffield(i+1,j-1,k)*p->DZN[KP]*p->DYN[JM1]+voffield(i+1,j+1,k)*p->DZN[KP]*p->DYN[JP1];
    vofsumdown=voffield(i-1,j,k-1)*p->DZN[KM1]*p->DYN[JP]+voffield(i-1,j,k)*p->DZN[KP]*p->DYN[JP]+voffield(i-1,j,k+1)*p->DZN[KP1]*p->DYN[JP]
            +voffield(i-1,j-1,k)*p->DZN[KP]*p->DYN[JM1]+voffield(i-1,j+1,k)*p->DZN[KP]*p->DYN[JP1];
    
    if(vofsumup>vofsumdown)
        sign=-1.0;
    else
        sign=1.0;
        
    nx_Cx=1.0*sign;
    nz_Cx=-(xsum_zp-xsum_zm)/(p->DZP[KM1]+p->DZP[KP]);
    ny_Cx=-(xsum_yp-xsum_ym)/(p->DYP[JM1]+p->DYP[JP]);
    nsum=sqrt(nx_Cx*nx_Cx+ny_Cx*ny_Cx+nz_Cx*nz_Cx);
    if(nsum!=nsum)
        cout<<"nsumNAN Cx"<<endl;
    
    nx_Cx=nx_Cx/nsum;
    nz_Cx=nz_Cx/nsum;
    ny_Cx=ny_Cx/nsum;
    
    //Candidate CC3 height function is y dimension
    ysum_zm=voffield(i,j-1,k-1)*p->DYN[JM1]+voffield(i,j,k-1)*p->DYN[JP]+voffield(i,j+1,k-1)*p->DYN[JP1];
    ysum_zp=voffield(i,j-1,k+1)*p->DYN[JM1]+voffield(i,j,k+1)*p->DYN[JP]+voffield(i,j+1,k+1)*p->DYN[JP1];
    ysum_xm=voffield(i-1,j-1,k)*p->DYN[JM1]+voffield(i-1,j,k)*p->DYN[JP]+voffield(i-1,j+1,k)*p->DYN[JP1];
    ysum_xp=voffield(i+1,j-1,k)*p->DYN[JM1]+voffield(i+1,j,k)*p->DYN[JP]+voffield(i+1,j+1,k)*p->DYN[JP1];
    vofsumup=voffield(i,j+1,k)*p->DXN[IP]*p->DZN[KP]+voffield(i-1,j+1,k)*p->DXN[IM1]*p->DZN[KP]+voffield(i+1,j+1,k)*p->DXN[IP1]+p->DZN[KP]
            +voffield(i,j+1,k-1)*p->DXN[IP]*p->DZN[KM1]+voffield(i,j+1,k+1)*p->DXN[IP]*p->DZN[KP1];
    vofsumdown=voffield(i,j-1,k)*p->DXN[IP]*p->DZN[KP]+voffield(i-1,j-1,k)*p->DXN[IM1]*p->DZN[KP]+voffield(i+1,j-1,k)*p->DXN[IP1]+p->DZN[KP]
            +voffield(i,j-1,k-1)*p->DXN[IP]*p->DZN[KM1]+voffield(i,j-1,k+1)*p->DXN[IP]*p->DZN[KP1];
    
    if(vofsumup>vofsumdown)
        sign=-1.0;
    else
        sign=1.0;
    
    ny_Cy=1.0*sign;
    nx_Cy=-(ysum_xp-ysum_xm)/(p->DXP[IM1]+p->DXP[IP]);
    nz_Cy=-(ysum_zp-ysum_zm)/(p->DZP[KM1]+p->DZP[KP]);
    nsum=sqrt(nx_Cy*nx_Cy+ny_Cy*ny_Cy+nz_Cy*nz_Cy);
    if(nsum!=nsum)
        cout<<"nsumNAN Cy"<<endl;
        
    nx_Cy=nx_Cy/nsum;
    nz_Cy=nz_Cy/nsum;
    ny_Cy=ny_Cy/nsum;

    
// figure out which CC Candidate is used and inside compare tou Young
    if(((fabs(nz_Cz)>=fabs(nx_Cx) && fabs(nz_Cz)>=fabs(ny_Cy)) || (nx_Cx!=nx_Cx && ny_Cy!=ny_Cy)) && nz_Cz==nz_Cz)
    {
        nx(i,j,k)=nx_Cz;
        ny(i,j,k)=ny_Cz;
        nz(i,j,k)=nz_Cz;
            
        if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_Cz"<<endl;
        if(ny(i,j,k)!=ny(i,j,k))
                cout<<"NAN ny_Cz"<<endl;
        if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_Cz"<<endl;
                
    }
    else if(((fabs(nx_Cx)>=fabs(nz_Cz) && fabs(nx_Cx)>=fabs(ny_Cy)) || (nz_Cz!=nz_Cz && ny_Cy!=nz_Cy)) && nx_Cx==nx_Cx)
    {
        nx(i,j,k)=nx_Cx;
        ny(i,j,k)=ny_Cx;
        nz(i,j,k)=nz_Cx;
            
        if(nx(i,j,k)!=nx(i,j,k))
            cout<<"NAN nx_Cx"<<endl;
        if(ny(i,j,k)!=ny(i,j,k))
            cout<<"NAN ny_Cx"<<endl;
        if(nz(i,j,k)!=nz(i,j,k))
            cout<<"NAN nz_Cx"<<endl;
    }
    else if(ny_Cy==ny_Cy)
    {
        nx(i,j,k)=nx_Cy;
        ny(i,j,k)=ny_Cy;
        nz(i,j,k)=nz_Cy;
            
        if(nx(i,j,k)!=nx(i,j,k))
            cout<<"NAN nx_Cy"<<endl;
        if(ny(i,j,k)!=ny(i,j,k))
            cout<<"NAN ny_Cy"<<endl;
        if(nz(i,j,k)!=nz(i,j,k))
            cout<<"NAN nz_Cy"<<endl;
    }
    else
    {
        cout<<"All normal options NAN!";
        nx(i,j,k)=0.0;
        ny(i,j,k)=0.0;
        nz(i,j,k)=1.0;
    }
    return;
}

void VOF_PLIC::calcNormalMYC2D(fdm* a,lexer* p, field& voffield)
{
    //First two Candidates by Centred Columns Scheme
    
    double zsum_xm, zsum_xp, xsum_zm, xsum_zp;
    double vofsumup,vofsumdown,sign;
    double nx_Cz, nz_Cz, nx_Cx, nz_Cx, nx_CC, nz_CC;
    double ny_all=0.0;
    double nsum;
    
    //Candidate CC1 height function is z dimension
    zsum_xm=voffield(i-1,j,k-1)*p->DZN[KM1]+voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k+1)*p->DZN[KP1];
    zsum_xp=voffield(i+1,j,k-1)*p->DZN[KM1]+voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k+1)*p->DZN[KP1];
    vofsumup=voffield(i-1,j,k+1)*p->DXN[IM1]+voffield(i,j,k+1)*p->DXN[IP]+voffield(i+1,j,k+1)*p->DXN[IP1];
    vofsumdown=voffield(i-1,j,k-1)*p->DXN[IM1]+voffield(i,j,k-1)*p->DXN[IP]+voffield(i+1,j,k-1)*p->DXN[IP1];
    
    if(vofsumup>vofsumdown)
        sign=-1.0;
    else
        sign=1.0;
    
    nz_Cz=1.0*sign;
    nx_Cz=-(zsum_xp-zsum_xm)/(p->DXP[IM1]+p->DXP[IP]);
    nsum=sqrt(nx_Cz*nx_Cz+nz_Cz*nz_Cz);
    if(nsum!=nsum)
        cout<<"nsumNAN Cz"<<endl;
    nz_Cz=nz_Cz/nsum;
    nx_Cz=nx_Cz/nsum;
    
    
    //Candidate CC2 height function is x dimension
    xsum_zm=voffield(i-1,j,k-1)*p->DXN[IM1]+voffield(i,j,k-1)*p->DXN[IP]+voffield(i+1,j,k-1)*p->DXN[IP1];
    xsum_zp=voffield(i-1,j,k+1)*p->DXN[IM1]+voffield(i,j,k+1)*p->DXN[IP]+voffield(i+1,j,k+1)*p->DXN[IP1];
    vofsumup=voffield(i+1,j,k-1)*p->DZN[KM1]+voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k+1)*p->DZN[KP1];
    vofsumdown=voffield(i-1,j,k-1)*p->DZN[KM1]+voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k+1)*p->DZN[KP1] ;
    
    if(vofsumup>vofsumdown)
        sign=-1.0;
    else
        sign=1.0;
        
    nx_Cx=1.0*sign;
    nz_Cx=-(xsum_zp-xsum_zm)/(p->DZP[KM1]+p->DZP[KP]);
    nsum=sqrt(nx_Cx*nx_Cx+nz_Cx*nz_Cx);
    if(nsum!=nsum)
        cout<<"nsumNAN Cx"<<endl;
    
    nx_Cx=nx_Cx/nsum;
    nz_Cx=nz_Cx/nsum;
    
    //Third Candidate by Youngs as average of all 4 cell corners 1=i+1,k+1 2=i+1,k-1, 3=i-1,k-1, 4=i-1,k+1
    
    /*double nx1,nz1,nx2,nz2,nx3,nz3,nx4,nz4;
    double nx_CY, nz_CY;
    
    nx1=-((voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k+1)*p->DZN[KP1])/2.0
            -(voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k+1)*p->DZN[KP1])/2.0)
            /p->DXP[IP];
            
    nz1=-((voffield(i,j,k+1)*p->DXN[IP]+voffield(i+1,j,k+1))*p->DXN[IP1]/2.0
            -(voffield(i,j,k)*p->DXN[IP]+voffield(i+1,j,k)*p->DXN[IP1])/2.0)
            /p->DZP[IP];
            
    nx2=-((voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k-1)*p->DZN[KM1])/2.0
            -(voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k-1)*p->DZN[KM1])/2.0)
            /p->DXP[IP];
            
    nz2=-((voffield(i,j,k)*p->DXN[IP]+voffield(i+1,j,k)*p->DXN[IP1])/2.0
            -(voffield(i,j,k-1)*p->DXN[IP]+voffield(i+1,j,k-1)*p->DXN[IP1])/2.0)
            /p->DZP[KM1];
            
    nx3=-((voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k-1)*p->DZN[KM1])/2.0
            -(voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k-1)*p->DZN[KM1])/2.0)
            /p->DXP[IM1];
            
    nz3=-((voffield(i,j,k)*p->DXN[IP]+voffield(i-1,j,k)*p->DXN[IM1])/2.0
            -(voffield(i,j,k-1)*p->DXN[IP]+voffield(i-1,j,k-1)*p->DXN[IM1])/2.0)
            /p->DZP[KM1];
            
    nx4=-((voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k+1)*p->DZN[KP1])/2.0
            -(voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k+1)*p->DZN[KP1])/2.0)
            /p->DXP[IM1];
            
    nz4=-((voffield(i,j,k+1)*p->DXN[IP]+voffield(i-1,j,k+1)*p->DXN[IM1])/2.0
            -(voffield(i,j,k)*p->DXN[IP]+voffield(i-1,j,k)*p->DXN[IM1])/2.0)
            /p->DZP[KP];
            
    nsum=sqrt(nx1*nx1+nz1*nz1);
    if(nsum!=nsum)
        cout<<"nsumNAN n1"<<endl;
    nx1=nx1/nsum;
    nz1=nz1/nsum;
    
    nsum=sqrt(nx2*nx2+nz2*nz2);
    if(nsum!=nsum)
        cout<<"nsumNAN n2"<<endl;
    nx2=nx2/nsum;
    nz2=nz2/nsum;
    
    nsum=sqrt(nx3*nx3+nz3*nz3);
    if(nsum!=nsum)
        cout<<"nsumNAN n3"<<endl;
    nx3=nx3/nsum;
    nz3=nz3/nsum;
    
    nsum=sqrt(nx4*nx4+nz4*nz4);
    if(nsum!=nsum)
        cout<<"nsumNAN n4"<<endl;
    nx4=nx4/nsum;
    nz4=nz4/nsum;
    
    double divN=0.0;
    double vecsumX=0.0;
    double vecsumZ=0.0;
    
    if(nx1==nx1 && nz1==nz1)
    {
        vecsumX+=nx1;
        vecsumZ+=nz1;
        divN+=1.0;
    }
    
    if(nx2==nx2 && nz2==nz2)
    {
        vecsumX+=nx2;
        vecsumZ+=nz2;
        divN+=1.0;
    }
    
    if(nx3==nx3 && nz3==nz3)
    {
        vecsumX+=nx3;
        vecsumZ+=nz3;
        divN+=1.0;
    }
    
    if(nx4==nx4 && nz4==nz4)
    {
        vecsumX+=nx4;
        vecsumZ+=nz4;
        divN+=1.0;
    }
    
    nx_CY=vecsumX/divN;
    nz_CY=vecsumZ/divN;
    nsum=sqrt(nx_CY*nx_CY+nz_CY*nz_CY);
    if(nsum!=nsum)
        cout<<"nsumNAN CY"<<" divN="<<divN<<endl;
        
    nx_CY=nx_CY/nsum;
    nz_CY=nz_CY/nsum;
    
    //zcomp
    vofsumup=voffield(i-1,j,k+1)+voffield(i,j,k+1)+voffield(i+1,j,k+1);
    vofsumdown=voffield(i-1,j,k-1)+voffield(i,j,k-1)+voffield(i+1,j,k-1);
    
    if(vofsumup>vofsumdown)
        nz_CY=-fabs(nz_CY);
    else
        nz_CY=fabs(nz_CY);
    
    //xcomp
    vofsumup=voffield(i+1,j,k-1)+voffield(i+1,j,k)+voffield(i+1,j,k+1);
    vofsumdown=voffield(i-1,j,k-1)+voffield(i-1,j,k)+voffield(i-1,j,k+1);
    
    if(vofsumup>vofsumdown)
        nx_CY=-fabs(nx_CY);
    else
        nx_CY=fabs(nx_CY);*/
        
// figure out which CC Candidate is used and inside compare tou Young
    if((fabs(nz_Cz)>=fabs(nx_Cx) ||nx_Cx!=nx_Cx) && nz_Cz==nz_Cz)
    {
        nx_CC=nx_Cz;
        nz_CC=nz_Cz;
        
       // if(fabs(nz_CC)<fabs(nz_CY) || nz_CY!=nz_CY)
        {
            nx(i,j,k)=nx_CC;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CC;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_Cz"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_Cz"<<endl;
        }
       /* else
        {
            nx(i,j,k)=nx_CY;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CY;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_CY"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_CY"<<endl;
        }*/
    }
    else if(nx_Cx==nx_Cx)
    {
        nx_CC=nx_Cx;
        nz_CC=nz_Cx;
        
      //  if(fabs(nx_CC)<fabs(nx_CY) || nx_CY!=nx_CY)
        {
            nx(i,j,k)=nx_CC;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CC;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_Cx"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_Cx"<<endl;
        }
      /*  else
        {
            nx(i,j,k)=nx_CY;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CY;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_CY"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_CY"<<endl;
        }*/
    }
    else
    {
        cout<<"All normal options NAN!";
        nx(i,j,k)=0.0;
        ny(i,j,k)=0.0;
        nz(i,j,k)=1.0;
    }
    return;
}

void VOF_PLIC::calcNormalMYC2D_V2(fdm* a,lexer* p, field& voffield)
{
    //First two Candidates by Centred Columns Scheme
    
    double zsum_xm, zsum_xp, xsum_zm, xsum_zp;
    double vofsumup,vofsumdown,sign;
    double nx_Cz, nz_Cz, nx_Cx, nz_Cx, nx_CC, nz_CC;
    double ny_all=0.0;
    double nsum;
    
    //Candidate CC1 height function is z dimension
    zsum_xm=voffield(i-1,j,k-1)*p->DZN[KM1]+voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k+1)*p->DZN[KP1];
    zsum_xp=voffield(i+1,j,k-1)*p->DZN[KM1]+voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k+1)*p->DZN[KP1];
    vofsumup=voffield(i-1,j,k+1)*p->DXN[IM1]+voffield(i,j,k+1)*p->DXN[IP]+voffield(i+1,j,k+1)*p->DXN[IP1];
    vofsumdown=voffield(i-1,j,k-1)*p->DXN[IM1]+voffield(i,j,k-1)*p->DXN[IP]+voffield(i+1,j,k-1)*p->DXN[IP1];
    
    if(vofsumup>vofsumdown)
        sign=-1.0;
    else
        sign=1.0;
    
    nz_Cz=1.0*sign;
    nx_Cz=-(zsum_xp-zsum_xm)/(p->DXP[IM1]+p->DXP[IP]);
    nsum=sqrt(nx_Cz*nx_Cz+nz_Cz*nz_Cz);
    if(nsum!=nsum)
        cout<<"nsumNAN Cz"<<endl;
    nz_Cz=nz_Cz/nsum;
    nx_Cz=nx_Cz/nsum;
    
    
    //Candidate CC2 height function is x dimension
    xsum_zm=voffield(i-1,j,k-1)*p->DXN[IM1]+voffield(i,j,k-1)*p->DXN[IP]+voffield(i+1,j,k-1)*p->DXN[IP1];
    xsum_zp=voffield(i-1,j,k+1)*p->DXN[IM1]+voffield(i,j,k+1)*p->DXN[IP]+voffield(i+1,j,k+1)*p->DXN[IP1];
    vofsumup=voffield(i+1,j,k-1)*p->DZN[KM1]+voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k+1)*p->DZN[KP1];
    vofsumdown=voffield(i-1,j,k-1)*p->DZN[KM1]+voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k+1)*p->DZN[KP1] ;
    
    if(vofsumup>vofsumdown)
        sign=-1.0;
    else
        sign=1.0;
        
    nx_Cx=1.0*sign;
    nz_Cx=-(xsum_zp-xsum_zm)/(p->DZP[KM1]+p->DZP[KP]);
    nsum=sqrt(nx_Cx*nx_Cx+nz_Cx*nz_Cx);
    if(nsum!=nsum)
        cout<<"nsumNAN Cx"<<endl;
    
    nx_Cx=nx_Cx/nsum;
    nz_Cx=nz_Cx/nsum;
    
    //Third Candidate by Youngs as average of all 4 cell corners 1=i+1,k+1 2=i+1,k-1, 3=i-1,k-1, 4=i-1,k+1
    
    double nx1,nz1,nx2,nz2,nx3,nz3,nx4,nz4;
    double nx_CY, nz_CY;
    
    nx1=-((voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k+1)*p->DZN[KP1])/2.0
            -(voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k+1)*p->DZN[KP1])/2.0)
            /p->DXP[IP];
            
    nz1=-((voffield(i,j,k+1)*p->DXN[IP]+voffield(i+1,j,k+1))*p->DXN[IP1]/2.0
            -(voffield(i,j,k)*p->DXN[IP]+voffield(i+1,j,k)*p->DXN[IP1])/2.0)
            /p->DZP[IP];
            
    nx2=-((voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k-1)*p->DZN[KM1])/2.0
            -(voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k-1)*p->DZN[KM1])/2.0)
            /p->DXP[IP];
            
    nz2=-((voffield(i,j,k)*p->DXN[IP]+voffield(i+1,j,k)*p->DXN[IP1])/2.0
            -(voffield(i,j,k-1)*p->DXN[IP]+voffield(i+1,j,k-1)*p->DXN[IP1])/2.0)
            /p->DZP[KM1];
            
    nx3=-((voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k-1)*p->DZN[KM1])/2.0
            -(voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k-1)*p->DZN[KM1])/2.0)
            /p->DXP[IM1];
            
    nz3=-((voffield(i,j,k)*p->DXN[IP]+voffield(i-1,j,k)*p->DXN[IM1])/2.0
            -(voffield(i,j,k-1)*p->DXN[IP]+voffield(i-1,j,k-1)*p->DXN[IM1])/2.0)
            /p->DZP[KM1];
            
    nx4=-((voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k+1)*p->DZN[KP1])/2.0
            -(voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k+1)*p->DZN[KP1])/2.0)
            /p->DXP[IM1];
            
    nz4=-((voffield(i,j,k+1)*p->DXN[IP]+voffield(i-1,j,k+1)*p->DXN[IM1])/2.0
            -(voffield(i,j,k)*p->DXN[IP]+voffield(i-1,j,k)*p->DXN[IM1])/2.0)
            /p->DZP[KP];
            
    nsum=sqrt(nx1*nx1+nz1*nz1);
    if(nsum!=nsum)
        cout<<"nsumNAN n1"<<endl;
    nx1=nx1/nsum;
    nz1=nz1/nsum;
    
    nsum=sqrt(nx2*nx2+nz2*nz2);
    if(nsum!=nsum)
        cout<<"nsumNAN n2"<<endl;
    nx2=nx2/nsum;
    nz2=nz2/nsum;
    
    nsum=sqrt(nx3*nx3+nz3*nz3);
    if(nsum!=nsum)
        cout<<"nsumNAN n3"<<endl;
    nx3=nx3/nsum;
    nz3=nz3/nsum;
    
    nsum=sqrt(nx4*nx4+nz4*nz4);
    if(nsum!=nsum)
        cout<<"nsumNAN n4"<<endl;
    nx4=nx4/nsum;
    nz4=nz4/nsum;
    
    double divN=0.0;
    double vecsumX=0.0;
    double vecsumZ=0.0;
    
    if(nx1==nx1 && nz1==nz1)
    {
        vecsumX+=nx1;
        vecsumZ+=nz1;
        divN+=1.0;
    }
    
    if(nx2==nx2 && nz2==nz2)
    {
        vecsumX+=nx2;
        vecsumZ+=nz2;
        divN+=1.0;
    }
    
    if(nx3==nx3 && nz3==nz3)
    {
        vecsumX+=nx3;
        vecsumZ+=nz3;
        divN+=1.0;
    }
    
    if(nx4==nx4 && nz4==nz4)
    {
        vecsumX+=nx4;
        vecsumZ+=nz4;
        divN+=1.0;
    }
    
    nx_CY=vecsumX/divN;
    nz_CY=vecsumZ/divN;
    nsum=sqrt(nx_CY*nx_CY+nz_CY*nz_CY);
    if(nsum!=nsum)
        cout<<"nsumNAN CY"<<" divN="<<divN<<endl;
        
    nx_CY=nx_CY/nsum;
    nz_CY=nz_CY/nsum;
    
    //zcomp
  /*  vofsumup=voffield(i-1,j,k+1)+voffield(i,j,k+1)+voffield(i+1,j,k+1);
    vofsumdown=voffield(i-1,j,k-1)+voffield(i,j,k-1)+voffield(i+1,j,k-1);
    
    if(vofsumup>vofsumdown)
        nz_CY=-fabs(nz_CY);
    else
        nz_CY=fabs(nz_CY);
    
    //xcomp
    vofsumup=voffield(i+1,j,k-1)+voffield(i+1,j,k)+voffield(i+1,j,k+1);
    vofsumdown=voffield(i-1,j,k-1)+voffield(i-1,j,k)+voffield(i-1,j,k+1);
    
    if(vofsumup>vofsumdown)
        nx_CY=-fabs(nx_CY);
    else
        nx_CY=fabs(nx_CY);*/
        
// figure out which CC Candidate is used and inside compare tou Young
    if((fabs(nz_Cz)>=fabs(nx_Cx) ||nx_Cx!=nx_Cx) && nz_Cz==nz_Cz)
    {
        nx_CC=nx_Cz;
        nz_CC=nz_Cz;
        
        if((fabs(nz_CC)<fabs(nz_CY) ||nz_CY!=nz_CY) ||divN<3.5)
        {
            nx(i,j,k)=nx_CC;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CC;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_Cz"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_Cz"<<endl;
        }
        else
        {
            nx(i,j,k)=nx_CY;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CY;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_CY"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_CY"<<endl;
        }
    }
    else if(nx_Cx==nx_Cx)
    {
        nx_CC=nx_Cx;
        nz_CC=nz_Cx;
        
        if((fabs(nx_CC)<fabs(nx_CY) ||nx_CY!=nx_CY) ||divN<3.5)
        {
            nx(i,j,k)=nx_CC;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CC;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_Cx"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_Cx"<<endl;
        }
        else
        {
            nx(i,j,k)=nx_CY;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CY;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_CY"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_CY"<<endl;
        }
    }
    else
    {
        cout<<"All normal options NAN!";
        nx(i,j,k)=0.0;
        ny(i,j,k)=0.0;
        nz(i,j,k)=1.0;
    }
    return;
}

void VOF_PLIC::calcNormalMYC2D_V3(fdm* a,lexer* p, field& voffield)
{
    //First two Candidates by Centred Columns Scheme
    
    double zsum_xm, zsum_xp, xsum_zm, xsum_zp;
    double vofsumup,vofsumdown,sign;
    double nx_Cz, nz_Cz, nx_Cx, nz_Cx, nx_CC, nz_CC;
    double ny_all=0.0;
    double nsum;
    
    //Candidate CC1 height function is z dimension
    zsum_xm=voffield(i-1,j,k-1)*p->DZN[KM1]+voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k+1)*p->DZN[KP1];
    zsum_xp=voffield(i+1,j,k-1)*p->DZN[KM1]+voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k+1)*p->DZN[KP1];
    vofsumup=voffield(i-1,j,k+1)*p->DXN[IM1]+voffield(i,j,k+1)*p->DXN[IP]+voffield(i+1,j,k+1)*p->DXN[IP1];
    vofsumdown=voffield(i-1,j,k-1)*p->DXN[IM1]+voffield(i,j,k-1)*p->DXN[IP]+voffield(i+1,j,k-1)*p->DXN[IP1];
    
    if(vofsumup>vofsumdown)
        sign=-1.0;
    else
        sign=1.0;
    
    nz_Cz=1.0*sign;
    nx_Cz=-(zsum_xp-zsum_xm)/(p->DXP[IM1]+p->DXP[IP]);
    nsum=sqrt(nx_Cz*nx_Cz+nz_Cz*nz_Cz);
    if(nsum!=nsum)
        cout<<"nsumNAN Cz"<<endl;
    nz_Cz=nz_Cz/nsum;
    nx_Cz=nx_Cz/nsum;
    
    
    //Candidate CC2 height function is x dimension
    xsum_zm=voffield(i-1,j,k-1)*p->DXN[IM1]+voffield(i,j,k-1)*p->DXN[IP]+voffield(i+1,j,k-1)*p->DXN[IP1];
    xsum_zp=voffield(i-1,j,k+1)*p->DXN[IM1]+voffield(i,j,k+1)*p->DXN[IP]+voffield(i+1,j,k+1)*p->DXN[IP1];
    vofsumup=voffield(i+1,j,k-1)*p->DZN[KM1]+voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k+1)*p->DZN[KP1];
    vofsumdown=voffield(i-1,j,k-1)*p->DZN[KM1]+voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k+1)*p->DZN[KP1] ;
    
    if(vofsumup>vofsumdown)
        sign=-1.0;
    else
        sign=1.0;
        
    nx_Cx=1.0*sign;
    nz_Cx=-(xsum_zp-xsum_zm)/(p->DZP[KM1]+p->DZP[KP]);
    nsum=sqrt(nx_Cx*nx_Cx+nz_Cx*nz_Cx);
    if(nsum!=nsum)
        cout<<"nsumNAN Cx"<<endl;
    
    nx_Cx=nx_Cx/nsum;
    nz_Cx=nz_Cx/nsum;
    
    //Third Candidate by Youngs as average of all 4 cell corners 1=i+1,k+1 2=i+1,k-1, 3=i-1,k-1, 4=i-1,k+1
    
    double nx1,nz1,nx2,nz2,nx3,nz3,nx4,nz4;
    double nx_CY, nz_CY;
    
    nx1=-((voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k+1)*p->DZN[KP1])/2.0
            -(voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k+1)*p->DZN[KP1])/2.0)
            /p->DXP[IP];
            
    nz1=-((voffield(i,j,k+1)*p->DXN[IP]+voffield(i+1,j,k+1))*p->DXN[IP1]/2.0
            -(voffield(i,j,k)*p->DXN[IP]+voffield(i+1,j,k)*p->DXN[IP1])/2.0)
            /p->DZP[IP];
            
    nx2=-((voffield(i+1,j,k)*p->DZN[KP]+voffield(i+1,j,k-1)*p->DZN[KM1])/2.0
            -(voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k-1)*p->DZN[KM1])/2.0)
            /p->DXP[IP];
            
    nz2=-((voffield(i,j,k)*p->DXN[IP]+voffield(i+1,j,k)*p->DXN[IP1])/2.0
            -(voffield(i,j,k-1)*p->DXN[IP]+voffield(i+1,j,k-1)*p->DXN[IP1])/2.0)
            /p->DZP[KM1];
            
    nx3=-((voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k-1)*p->DZN[KM1])/2.0
            -(voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k-1)*p->DZN[KM1])/2.0)
            /p->DXP[IM1];
            
    nz3=-((voffield(i,j,k)*p->DXN[IP]+voffield(i-1,j,k)*p->DXN[IM1])/2.0
            -(voffield(i,j,k-1)*p->DXN[IP]+voffield(i-1,j,k-1)*p->DXN[IM1])/2.0)
            /p->DZP[KM1];
            
    nx4=-((voffield(i,j,k)*p->DZN[KP]+voffield(i,j,k+1)*p->DZN[KP1])/2.0
            -(voffield(i-1,j,k)*p->DZN[KP]+voffield(i-1,j,k+1)*p->DZN[KP1])/2.0)
            /p->DXP[IM1];
            
    nz4=-((voffield(i,j,k+1)*p->DXN[IP]+voffield(i-1,j,k+1)*p->DXN[IM1])/2.0
            -(voffield(i,j,k)*p->DXN[IP]+voffield(i-1,j,k)*p->DXN[IM1])/2.0)
            /p->DZP[KP];
            
    nsum=sqrt(nx1*nx1+nz1*nz1);
    if(nsum!=nsum)
        cout<<"nsumNAN n1"<<endl;
    nx1=nx1/nsum;
    nz1=nz1/nsum;
    
    nsum=sqrt(nx2*nx2+nz2*nz2);
    if(nsum!=nsum)
        cout<<"nsumNAN n2"<<endl;
    nx2=nx2/nsum;
    nz2=nz2/nsum;
    
    nsum=sqrt(nx3*nx3+nz3*nz3);
    if(nsum!=nsum)
        cout<<"nsumNAN n3"<<endl;
    nx3=nx3/nsum;
    nz3=nz3/nsum;
    
    nsum=sqrt(nx4*nx4+nz4*nz4);
    if(nsum!=nsum)
        cout<<"nsumNAN n4"<<endl;
    nx4=nx4/nsum;
    nz4=nz4/nsum;
    
    double divN=0.0;
    double vecsumX=0.0;
    double vecsumZ=0.0;
    
    if(nx1==nx1 && nz1==nz1)
    {
        vecsumX+=nx1;
        vecsumZ+=nz1;
        divN+=1.0;
    }
    
    if(nx2==nx2 && nz2==nz2)
    {
        vecsumX+=nx2;
        vecsumZ+=nz2;
        divN+=1.0;
    }
    
    if(nx3==nx3 && nz3==nz3)
    {
        vecsumX+=nx3;
        vecsumZ+=nz3;
        divN+=1.0;
    }
    
    if(nx4==nx4 && nz4==nz4)
    {
        vecsumX+=nx4;
        vecsumZ+=nz4;
        divN+=1.0;
    }
    
    nx_CY=vecsumX/divN;
    nz_CY=vecsumZ/divN;
    nsum=sqrt(nx_CY*nx_CY+nz_CY*nz_CY);
    if(nsum!=nsum)
        cout<<"nsumNAN CY"<<" divN="<<divN<<endl;
        
    nx_CY=nx_CY/nsum;
    nz_CY=nz_CY/nsum;
    
    //zcomp
    /*vofsumup=voffield(i-1,j,k+1)+voffield(i,j,k+1)+voffield(i+1,j,k+1);
    vofsumdown=voffield(i-1,j,k-1)+voffield(i,j,k-1)+voffield(i+1,j,k-1);
    
    if(vofsumup>vofsumdown)
        nz_CY=-fabs(nz_CY);
    else
        nz_CY=fabs(nz_CY);
    
    //xcomp
    vofsumup=voffield(i+1,j,k-1)+voffield(i+1,j,k)+voffield(i+1,j,k+1);
    vofsumdown=voffield(i-1,j,k-1)+voffield(i-1,j,k)+voffield(i-1,j,k+1);
    
    if(vofsumup>vofsumdown)
        nx_CY=-fabs(nx_CY);
    else
        nx_CY=fabs(nx_CY);*/
        
// figure out which CC Candidate is used and inside compare tou Young
    if((fabs(nz_Cz)>=fabs(nx_Cx) ||nx_Cx!=nx_Cx) && nz_Cz==nz_Cz)
    {
        nx_CC=nx_Cz;
        nz_CC=nz_Cz;
        
        if(fabs(nz_CC)<fabs(nz_CY) ||nz_CY!=nz_CY)
        {
            nx(i,j,k)=nx_CC;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CC;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_Cz"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_Cz"<<endl;
        }
        else
        {
            nx(i,j,k)=nx_CY;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CY;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_CY"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_CY"<<endl;
        }
    }
    else if(nx_Cx==nx_Cx)
    {
        nx_CC=nx_Cx;
        nz_CC=nz_Cx;
        
        if(fabs(nx_CC)<fabs(nx_CY) ||nx_CY!=nx_CY)
        {
            nx(i,j,k)=nx_CC;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CC;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_Cx"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_Cx"<<endl;
        }
        else
        {
            nx(i,j,k)=nx_CY;
            ny(i,j,k)=0.0;
            nz(i,j,k)=nz_CY;
            
            if(nx(i,j,k)!=nx(i,j,k))
                cout<<"NAN nx_CY"<<endl;
                
            if(nz(i,j,k)!=nz(i,j,k))
                cout<<"NAN nz_CY"<<endl;
        }
    }
    else
    {
        cout<<"All normal options NAN!";
        nx(i,j,k)=0.0;
        ny(i,j,k)=0.0;
        nz(i,j,k)=1.0;
    }
    return;
}