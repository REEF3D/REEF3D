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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_scalar_iweno.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_scalar_advec_CDS2.h"

nhflow_scalar_iweno::nhflow_scalar_iweno(lexer *p)
			:weno_nug_func(p), tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsi(1.0e-6)
{

    padvec = new nhflow_scalar_advec_CDS2(p);

}

nhflow_scalar_iweno::~nhflow_scalar_iweno()
{
}

void nhflow_scalar_iweno::start(lexer* p, fdm_nhf *d, double *F, int ipol, double *U, double *V, double *W)
{
    uf=vf=wf=0;

    DX=p->DXP;
    DY=p->DYP;
    DZ=p->DZP;
    

	count=0;

	LOOP
	{
        padvec->uadvec(ipol,U,iadvec,ivel2);
        padvec->vadvec(ipol,V,jadvec,jvel2);
        padvec->wadvec(ipol,W,kadvec,kvel2);
    

			if(iadvec>=0.0)
			{
            iqmin(p,d,F);
			is_min_x();
            weight_min_x();
			aij_south(p,d,F,d->L);
            }

			if(iadvec<0.0)
			{
            iqmax(p,d,F);
			is_max_x();
            weight_max_x();
			aij_north(p,d,F,d->L);
			}


			
			if(jadvec>=0.0 && p->j_dir==1)
			{
            jqmin(p,d,F);
			is_min_y();
            weight_min_y();
			aij_east(p,d,F,d->L);
			}

			if(jadvec<0.0 && p->j_dir==1)
			{
            jqmax(p,d,F);
			is_max_y();
            weight_max_y();
			aij_west(p,d,F,d->L);
			}



			if(kadvec>=0.0)
			{
            kqmin(p,d,F);
			is_min_z();
            weight_min_z();
			aij_bottom(p,d,F,d->L);
			}

			if(kadvec<0.0)
			{
            kqmax(p,d,F);
			is_max_z();
            weight_max_z();
			aij_top(p,d,F,d->L);
			}
                
		 ++count;
	}
}

void nhflow_scalar_iweno::aij_south(lexer* p, fdm_nhf *d, double *F, double *L)
{
	L[IJK]      -= F[Im3JK]*(-w3x*qfx[IP][uf][2][0]/DX[IM3])*iadvec
    
                 + F[Im2JK]*(w2x*qfx[IP][uf][1][1]/DX[IM2]
                           - w3x*(1.0 - qfx[IP][uf][2][0] - qfx[IP][uf][2][1])/DX[IM2]
                           + w3x*(qfx[IP][uf][2][0])/DX[IM3])*iadvec
                            
                 + F[Ip2JK]*(-w1x*qfx[IP][uf][0][1]/DX[IP1])*iadvec;  
                        
                           
	d->M.p[count] = (- w1x*(1.0 - qfx[IP][uf][0][0] + qfx[IP][uf][0][1])/DX[IP] 
                    - w2x*(qfx[IP][uf][1][0])/DX[IP]
                    + w1x*(qfx[IP][uf][0][0])/DX[IM1]
                    + w2x*(1.0 - qfx[IP][uf][1][0] + qfx[IP][uf][1][1])/DX[IM1]
                    + w3x*(qfx[IP][uf][2][1])/DX[IM1])*iadvec;
					 
	d->M.s[count] = (- w1x*qfx[IP][uf][0][0]/DX[IM1]
                    - w2x*(1.0 - qfx[IP][uf][1][0] + qfx[IP][uf][1][1])/DX[IM1]
                    - w3x*(qfx[IP][uf][2][1])/DX[IM1]
                    - w2x*(qfx[IP][uf][1][1])/DX[IM2]
                    + w3x*(1.0 - qfx[IP][uf][2][0] - qfx[IP][uf][2][1])/DX[IM2])*iadvec;
                           
	d->M.n[count] = (  w1x*qfx[IP][uf][0][1]/DX[IP1]
                    + w1x*(1.0 - qfx[IP][uf][0][0] + qfx[IP][uf][0][1])/DX[IP]
                    + w2x*(qfx[IP][uf][1][0])/DX[IP])*iadvec;
}

void nhflow_scalar_iweno::aij_north(lexer* p, fdm_nhf *d, double *F, double *L)
{
	L[IJK]      -= F[Im2JK]*(w3x*qfx[IP][uf][5][1]/DX[IM2])*iadvec
    
                 + F[Ip2JK]*(-w1x*qfx[IP][uf][3][1]/DX[IP2]
                            + w1x*(1.0 - qfx[IP][uf][3][0] - qfx[IP][uf][3][1])/DX[IP1]
                            - w2x*(qfx[IP][uf][4][1])/DX[IP1])*iadvec
                            
                 + F[Ip3JK]*(w1x*qfx[IP][uf][3][1]/DX[IP2])*iadvec;
                        
				 
	d->M.p[count] = (- w1x*(qfx[IP][uf][3][0])/DX[IP] 
                    - w2x*(1.0 - qfx[IP][uf][4][0] + qfx[IP][uf][4][1])/DX[IP] 
                    - w3x*(qfx[IP][uf][5][0])/DX[IP]
                    + w2x*(qfx[IP][uf][4][0])/DX[IM1]
                    + w3x*(1.0 - qfx[IP][uf][5][0] + qfx[IP][uf][5][1])/DX[IM1])*iadvec;
					 
	d->M.s[count] = (- w2x*qfx[IP][uf][4][0]/DX[IM1]
                    - w3x*(1.0 - qfx[IP][uf][5][0] + qfx[IP][uf][5][1])/DX[IM1]
                    - w3x*(qfx[IP][uf][5][1])/DX[IM2])*iadvec;
                           
	d->M.n[count] = (-w1x*(1.0 - qfx[IP][uf][3][0] - qfx[IP][uf][3][1])/DX[IP1]
                    + w2x*(qfx[IP][uf][4][1])/DX[IP1]
                    + w1x*(qfx[IP][uf][3][0])/DX[IP]
                    + w2x*(1.0 - qfx[IP][uf][4][0] + qfx[IP][uf][4][1])/DX[IP]
                    + w3x*qfx[IP][uf][5][0]/DX[IP])*iadvec;
}

void nhflow_scalar_iweno::aij_east(lexer* p, fdm_nhf *d, double *F, double *L)
{
	L[IJK]      += F[IJm3K]*(w3y*qfy[JP][vf][2][0]/DY[JM3])*jadvec
    
                 + F[IJm2K]*(-w2y*qfy[JP][vf][1][1]/DY[JM2]
                           + w3y*(1.0 - qfy[JP][vf][2][0] - qfy[JP][vf][2][1])/DY[JM2]
                           - w3y*(qfy[JP][vf][2][0])/DY[JM3])*jadvec
                            
                 + F[IJp2K]*(w1y*qfy[JP][vf][0][1]/DY[JP1])*jadvec;
                        
                           
	d->M.p[count] +=(-w1y*(1.0 - qfy[JP][vf][0][0] + qfy[JP][vf][0][1])/DY[JP] 
                    - w2y*(qfy[JP][vf][1][0])/DY[JP]
                    + w1y*(qfy[JP][vf][0][0])/DY[JM1]
                    + w2y*(1.0 - qfy[JP][vf][1][0] + qfy[JP][vf][1][1])/DY[JM1]
                    + w3y*(qfy[JP][vf][2][1])/DY[JM1])*jadvec;
					 
	d->M.e[count] = (-w1y*qfy[JP][vf][0][0]/DY[JM1]
                    - w2y*(1.0 - qfy[JP][vf][1][0] + qfy[JP][vf][1][1])/DY[JM1]
                    - w3y*(qfy[JP][vf][2][1])/DY[JM1]
                    + w2y*(-qfy[JP][vf][1][1])/DY[JM2]
                    + w3y*(1.0 - qfy[JP][vf][2][0] - qfy[JP][vf][2][1])/DY[JM2])*jadvec;
                           
	d->M.w[count] = ( w1y*qfy[JP][vf][0][1]/DY[JP1]
                    + w1y*(1.0 - qfy[JP][vf][0][0] + qfy[JP][vf][0][1])/DY[JP]
                    + w2y*(qfy[JP][vf][1][0])/DY[JP])*jadvec;
}

void nhflow_scalar_iweno::aij_west(lexer* p, fdm_nhf *d, double *F, double *L)
{
	L[IJK]      += F[IJm2K]*(-w3y*qfy[JP][vf][5][1]/DY[JM2])*jadvec
    
                 + F[IJp2K]*(w1y*qfy[JP][vf][3][1]/DY[JP2]
                           - w1y*(1.0 - qfy[JP][vf][3][0] - qfy[JP][vf][3][1])/DY[JP1]
                           + w2y*(qfy[JP][vf][4][1])/DY[JP1])*jadvec
                            
                 + F[IJp3K]*(-w1y*qfy[JP][vf][3][1]/DY[JP2])*jadvec;
                        
				 
	d->M.p[count] +=(-w1y*(qfy[JP][vf][3][0])/DY[JP] 
                    - w2y*(1.0 - qfy[JP][vf][4][0] + qfy[JP][vf][4][1])/DY[JP] 
                    - w3y*(qfy[JP][vf][5][0])/DY[JP]
                    + w2y*(qfy[JP][vf][4][0])/DY[JM1]
                    + w3y*(1.0 - qfy[JP][vf][5][0] + qfy[JP][vf][5][1])/DY[JM1])*jadvec;
					 
	d->M.e[count] = (-w2y*qfy[JP][vf][4][0]/DY[JM1]
                    - w3y*(1.0 - qfy[JP][vf][5][0] + qfy[JP][vf][5][1])/DY[JM1]
                    - w3y*(qfy[JP][vf][5][1])/DY[JM2])*jadvec;
                           
	d->M.w[count] = (-w1y*(1.0 - qfy[JP][vf][3][0] - qfy[JP][vf][3][1])/DY[JP1]
                    + w2y*(qfy[JP][vf][4][1])/DY[JP1]
                    + w1y*(qfy[JP][vf][3][0])/DY[JP]
                    + w2y*(1.0 - qfy[JP][vf][4][0] + qfy[JP][vf][4][1])/DY[JP]
                    + w3y*qfy[JP][vf][5][0]/DY[JP])*jadvec;
}

void nhflow_scalar_iweno::aij_bottom(lexer* p, fdm_nhf *d, double *F, double *L)
{
	L[IJK]     -= F[IJKm3]*(-w3z*qfz[KP][wf][2][0]/DZ[KM3])*kadvec
    
                + F[IJKm2]*( w2z*qfz[KP][wf][1][1]/DZ[KM2]
                             - w3z*(1.0 - qfz[KP][wf][2][0] - qfz[KP][wf][2][1])/DZ[KM2]
                             + w3z*(qfz[KP][wf][2][0])/DZ[KM3])*kadvec
                            
                + F[IJKp2]*(-w1z*qfz[KP][wf][0][1]/DZ[KP1])*kadvec;
                        
                           
	d->M.p[count]+= (-w1z*(1.0 - qfz[KP][wf][0][0] + qfz[KP][wf][0][1])/DZ[KP] 
                    - w2z*(qfz[KP][wf][1][0])/DZ[KP]
                    + w1z*(qfz[KP][wf][0][0])/DZ[KM1]
                    + w2z*(1.0 - qfz[KP][wf][1][0] + qfz[KP][wf][1][1])/DZ[KM1]
                    + w3z*(qfz[KP][wf][2][1])/DZ[KM1])*kadvec;
					 
	d->M.b[count] = (-w1z*qfz[KP][wf][0][0]/DZ[KM1]
                    - w2z*(1.0 - qfz[KP][wf][1][0] + qfz[KP][wf][1][1])/DZ[KM1]
                    - w3z*(qfz[KP][wf][2][1])/DZ[KM1]
                    - w2z*(qfz[KP][wf][1][1])/DZ[KM2]
                    + w3z*(1.0 - qfz[KP][wf][2][0] - qfz[KP][wf][2][1])/DZ[KM2])*kadvec;
                           
	d->M.t[count] = ( w1z*qfz[KP][wf][0][1]/DZ[KP1]
                    + w1z*(1.0 - qfz[KP][wf][0][0] + qfz[KP][wf][0][1])/DZ[KP]
                    + w2z*(qfz[KP][wf][1][0])/DZ[KP])*kadvec;
}

void nhflow_scalar_iweno::aij_top(lexer* p, fdm_nhf *d, double *F, double *L)
{
	L[IJK]     -= F[IJKm2]*(w3z*qfz[KP][wf][5][1]/DZ[KM2])*kadvec
    
                + F[IJKp2]*(- w1z*qfz[KP][wf][3][1]/DZ[KP2]
                             + w1z*(1.0 - qfz[KP][wf][3][0] - qfz[KP][wf][3][1])/DZ[KP1]
                             - w2z*(qfz[KP][wf][4][1])/DZ[KP1])*kadvec
                            
                + F[IJKp3]*(w1z*qfz[KP][wf][3][1]/DZ[KP2])*kadvec;
                        
				 
	d->M.p[count]+= (-w1z*(qfz[KP][wf][3][0])/DZ[KP] 
                    - w2z*(1.0 - qfz[KP][wf][4][0] + qfz[KP][wf][4][1])/DZ[KP] 
                    - w3z*(qfz[KP][wf][5][0])/DZ[KP]
                    + w2z*(qfz[KP][wf][4][0])/DZ[KM1]
                    + w3z*(1.0 - qfz[KP][wf][5][0] + qfz[KP][wf][5][1])/DZ[KM1])*kadvec;
					 
	d->M.b[count] = (-w2z*qfz[KP][wf][4][0]/DZ[KM1]
                    - w3z*(1.0 - qfz[KP][wf][5][0] + qfz[KP][wf][5][1])/DZ[KM1]
                    - w3z*(qfz[KP][wf][5][1])/DZ[KM2])*kadvec;
                           
	d->M.t[count] = (-w1z*(1.0 - qfz[KP][wf][3][0] - qfz[KP][wf][3][1])/DZ[KP1]
                    + w2z*(qfz[KP][wf][4][1])/DZ[KP1]
                    + w1z*(qfz[KP][wf][3][0])/DZ[KP]
                    + w2z*(1.0 - qfz[KP][wf][4][0] + qfz[KP][wf][4][1])/DZ[KP]
                    + w3z*qfz[KP][wf][5][0]/DZ[KP])*kadvec;
}

void nhflow_scalar_iweno::iqmin(lexer *p,fdm_nhf *d, double *F)
{	
	q1 = (F[Im2JK]-F[Im3JK])/DX[IM3];
	q2 = (F[Im1JK]-F[Im2JK])/DX[IM2];
	q3 = (F[IJK]-F[Im1JK])/DX[IM1];
	q4 = (F[Ip1JK]-F[IJK])/DX[IP];
	q5 = (F[Ip2JK]-F[Ip1JK])/DX[IP1];
}

void nhflow_scalar_iweno::jqmin(lexer *p,fdm_nhf *d, double *F)
{
	q1 = (F[IJm2K]-F[IJm3K])/DY[JM3];
	q2 = (F[IJm1K]-F[IJm2K])/DY[JM2];
	q3 = (F[IJK]-F[IJm1K])/DY[JM1];
	q4 = (F[IJp1K]-F[IJK])/DY[JP];
	q5 = (F[IJp2K]-F[IJp1K])/DY[JP1];
}

void nhflow_scalar_iweno::kqmin(lexer *p,fdm_nhf *d, double *F)
{
	q1 = (F[IJKm2]-F[IJKm3])/DZ[KM3];
	q2 = (F[IJKm1]-F[IJKm2])/DZ[KM2];
	q3 = (F[IJK]-F[IJKm1])/DZ[KM1];
	q4 = (F[IJKp1]-F[IJK])/DZ[KP];
	q5 = (F[IJKp2]-F[IJKp1])/DZ[KP1];
}

void nhflow_scalar_iweno::iqmax(lexer *p,fdm_nhf *d, double *F)
{
    q1 = (F[Im1JK]-F[Im2JK])/DX[IM2];
	q2 = (F[IJK]-F[Im1JK])/DX[IM1];
	q3 = (F[Ip1JK]-F[IJK])/DX[IP];
	q4 = (F[Ip2JK]-F[Ip1JK])/DX[IP1];
	q5 = (F[Ip3JK]-F[Ip2JK])/DX[IP2];
}

void nhflow_scalar_iweno::jqmax(lexer *p,fdm_nhf *d, double *F)
{
	q1 = (F[IJm1K]-F[IJm2K])/DY[JM2];
	q2 = (F[IJK]-F[IJm1K])/DY[JM1];
	q3 = (F[IJp1K]-F[IJK])/DY[JP];
	q4 = (F[IJp2K]-F[IJp1K])/DY[JP1];
	q5 = (F[IJp3K]-F[IJp2K])/DY[JP2];
}

void nhflow_scalar_iweno::kqmax(lexer *p,fdm_nhf *d, double *F)
{
	q1 = (F[IJKm1]-F[IJKp2])/DZ[KM2];
	q2 = (F[IJK]-F[IJKp1])/DZ[KM1];
	q3 = (F[IJKp1]-F[IJK])/DZ[KP];
	q4 = (F[IJKp2]-F[IJKp1])/DZ[KP1];
	q5 = (F[IJKp3]-F[IJKp2])/DZ[KP2];
}
