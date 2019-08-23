/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"iweno_hj_nug.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS2_vrans.h"


iweno_hj_nug::iweno_hj_nug(lexer *p)
			:weno_nug_func(p), tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsi(1.0e-6),deltin (1.0/p->DXM)
{
    if(p->B269==0)
    pflux = new flux_HJ_CDS2(p);
    
    if(p->B269>=1)
    pflux = new flux_HJ_CDS2_vrans(p);
}

iweno_hj_nug::~iweno_hj_nug()
{
}

void iweno_hj_nug::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    uf=vf=wf=0;
    
    if(ipol==1)
    wenoloop1(p,a,b,ipol,uvel,vvel,wvel);

    if(ipol==2)
    wenoloop2(p,a,b,ipol,uvel,vvel,wvel);

    if(ipol==3)
    wenoloop3(p,a,b,ipol,uvel,vvel,wvel);

    if(ipol==4)
    wenoloop4(p,a,b,ipol,uvel,vvel,wvel);

    if(ipol==5)
    wenoloop4(p,a,b,ipol,uvel,vvel,wvel);
}

void iweno_hj_nug::wenoloop1(lexer *p, fdm *a, field& f, int ipol, field& uvel, field& vvel, field& wvel)
{
    DX=p->DXN;
    DY=p->DYP;
    DZ=p->DZP;
    
    uf=1;
    
	count=0;
    
	ULOOP
	{
        pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
    

			if(iadvec>=0.0)
			{
            iqmin(p,a,f);
			is_min_x();
            weight_min_x();
			aij_south(p,a,f,a->F);
			}

			if(iadvec<0.0)
			{
            iqmax(p,a,f);
			is_max_x();
            weight_max_x();
			aij_north(p,a,f,a->F);
			}


			
			if(jadvec>=0.0)
			{
            jqmin(p,a,f);
			is_min_y();
            weight_min_y();
			aij_east(p,a,f,a->F);
			}

			if(jadvec<0.0)
			{
            jqmax(p,a,f);
			is_max_y();
            weight_max_y();
			aij_west(p,a,f,a->F);
			}



			if(kadvec>=0.0)
			{
            kqmin(p,a,f);
			is_min_z();
            weight_min_z();
			aij_bottom(p,a,f,a->F);
			}

			if(kadvec<0.0)
			{
            kqmax(p,a,f);
			is_max_z();
            weight_max_z();
			aij_top(p,a,f,a->F);
			}
		 ++count;
	}
}

void iweno_hj_nug::wenoloop2(lexer *p, fdm *a, field& f, int ipol, field& uvel, field& vvel, field& wvel)
{
    
    DX=p->DXP;
    DY=p->DYN;
    DZ=p->DZP;
    
    vf=1;
    
	count=0;

	VLOOP
	{
        pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
    

			if(iadvec>=0.0)
			{
            iqmin(p,a,f);
			is_min_x();
            weight_min_x();
			aij_south(p,a,f,a->G);
			}

			if(iadvec<0.0)
			{
            iqmax(p,a,f);
			is_max_x();
            weight_max_x();
			aij_north(p,a,f,a->G);
			}


			
			if(jadvec>=0.0)
			{
            jqmin(p,a,f);
			is_min_y();
            weight_min_y();
			aij_east(p,a,f,a->G);
			}

			if(jadvec<0.0)
			{
            jqmax(p,a,f);
			is_max_y();
            weight_max_y();
			aij_west(p,a,f,a->G);
			}



			if(kadvec>=0.0)
			{
            kqmin(p,a,f);
			is_min_z();
            weight_min_z();
			aij_bottom(p,a,f,a->G);
			}

			if(kadvec<0.0)
			{
            kqmax(p,a,f);
			is_max_z();
            weight_max_z();
			aij_top(p,a,f,a->G);
			}
		 ++count;
	}
}

void iweno_hj_nug::wenoloop3(lexer *p, fdm *a, field& f, int ipol, field& uvel, field& vvel, field& wvel)
{
    DX=p->DXP;
    DY=p->DYP;
    DZ=p->DZN;
    
    wf=1;
    
	count=0;

	WLOOP
	{
        pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
    

			if(iadvec>=0.0)
			{
            iqmin(p,a,f);
			is_min_x();
            weight_min_x();
			aij_south(p,a,f,a->H);
			}

			if(iadvec<0.0)
			{
            iqmax(p,a,f);
			is_max_x();
            weight_max_x();
			aij_north(p,a,f,a->H);
			}


			
			if(jadvec>=0.0)
			{
            jqmin(p,a,f);
			is_min_y();
            weight_min_y();
			aij_east(p,a,f,a->H);
			}

			if(jadvec<0.0)
			{
            jqmax(p,a,f);
			is_max_y();
            weight_max_y();
			aij_west(p,a,f,a->H);
			}



			if(kadvec>=0.0)
			{
            kqmin(p,a,f);
			is_min_z();
            weight_min_z();
			aij_bottom(p,a,f,a->H);
			}

			if(kadvec<0.0)
			{
            kqmax(p,a,f);
			is_max_z();
            weight_max_z();
			aij_top(p,a,f,a->H);
			}
		 ++count;
	}
}

void iweno_hj_nug::wenoloop4(lexer *p, fdm *a, field& f, int ipol, field& uvel, field& vvel, field& wvel)
{
    DX=p->DXP;
    DY=p->DYP;
    DZ=p->DZP;
    
	count=0;

	LOOP
	{
        pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
    

			if(iadvec>=0.0)
			{
            iqmin(p,a,f);
			is_min_x();
            weight_min_x();
			aij_south(p,a,f,a->L);
            
            
            
            
            
           is_south(f);
			alpha_calc();
			weight_calc();
            
            Ft    = (w1*third)*iadvec*deltin*f(i-3,j,k)
				 - (w2*sixth + w1*1.5)*iadvec*deltin*f(i-2,j,k)
				 + (w3*sixth)*iadvec*deltin*f(i+2,j,k);
				 
            Mp = (-w3*0.5 + w2*0.5 + w1*elvsix)*iadvec*deltin;
                             
            Ms = (-w3*third -w2 - w1*3.0)*iadvec*deltin;
            Mn = (w3 + w2*third)*iadvec*deltin;
            
            //if(p->mpirank==0)
            //cout<<"- M.p "<<a->M.p[count]<<" Mp "<<Mp<<" M.s "<<a->M.s[count]<<" Ms "<<Ms<<" M.n "<<a->M.n[count]<<" Mn "<<Mn<<" L: "<<a->L(i,j,k)<<" Ft: "<<Ft<<" ___  q1: "<<q1<<" q2: "<<q2<<" q3: "<<q3<<" q4: "<<q4<<" q5: "<<q5<<endl;
			}

			if(iadvec<0.0)
			{
            iqmax(p,a,f);
			is_max_x();
            weight_max_x();
			aij_north(p,a,f,a->L);
            
            
            
            
            is_north(f);
			alpha_calc();
			weight_calc();
            
            
            Ft   = -(w3*sixth)*iadvec*deltin*f(i-2,j,k)
				-	(-w1*1.5 - w2*sixth)*iadvec*deltin*f(i+2,j,k)
				-   (w1*third)*iadvec*deltin*f(i+3,j,k);
					
            Mp = (-w1*elvsix - w2*0.5 + w3*0.5)*iadvec*deltin;
                             
            Ms = (-w2*third - w3)*iadvec*deltin;
            Mn = (w1*3.0 + w2 + w3*third)*iadvec*deltin;
            
            if(p->mpirank==0)
            cout<<"+ M.p "<<a->M.p[count]<<" Mp "<<Mp<<" M.s "<<a->M.s[count]<<" Ms "<<Ms<<" M.n "<<a->M.n[count]<<" Mn "<<Mn<<" L: "<<a->L(i,j,k)<<" Ft: "<<Ft<<" ___  q1: "<<q1<<" q2: "<<q2<<" q3: "<<q3<<" q4: "<<q4<<" q5: "<<q5<<endl;
			}


			
			if(jadvec>=0.0)
			{
            jqmin(p,a,f);
			is_min_y();
            weight_min_y();
			aij_east(p,a,f,a->L);
			}

			if(jadvec<0.0)
			{
            jqmax(p,a,f);
			is_max_y();
            weight_max_y();
			aij_west(p,a,f,a->L);
			}



			if(kadvec>=0.0)
			{
            kqmin(p,a,f);
			is_min_z();
            weight_min_z();
			aij_bottom(p,a,f,a->L);
			}

			if(kadvec<0.0)
			{
            kqmax(p,a,f);
			is_max_z();
            weight_max_z();
			aij_top(p,a,f,a->L);
			}
            
        //cout<<"M.p "<<a->M.p[count]<<" M.s "<<a->M.s[count]<<" M.n "<<a->M.s[count]<<" M.b "<<a->M.b[count]<<" M.t "<<a->M.t[count]<<" L: "<<a->L(i,j,k)<<endl;
            
		 ++count;
	}
}

void iweno_hj_nug::aij_south(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k)    -= f(i-3,j,k)*(-w3x*qfx[IP][uf][2][0]/DX[IM3])*iadvec
    
                + f(i-2,j,k)*(w2x*qfx[IP][uf][1][1]/DX[IM2]
                           - w3x*(1.0 - qfx[IP][uf][2][0] - qfx[IP][uf][2][1])/DX[IM2]
                           + w3x*(qfx[IP][uf][2][0])/DX[IM3])*iadvec
                            
                + f(i+2,j,k)*(-w1x*qfx[IP][uf][0][1]/DX[IP1])*iadvec;  
                        
                           
	a->M.p[count] = (- w1x*(1.0 - qfx[IP][uf][0][0] + qfx[IP][uf][0][1])/DX[IP] 
                    - w2x*(qfx[IP][uf][1][0])/DX[IP]
                    + w1x*(qfx[IP][uf][0][0])/DX[IM1]
                    + w2x*(1.0 - qfx[IP][uf][1][0] + qfx[IP][uf][1][1])/DX[IM1]
                    + w3x*(qfx[IP][uf][2][1])/DX[IM1])*iadvec;
					 
	a->M.s[count] = (- w1x*qfx[IP][uf][0][0]/DX[IM1]
                    - w2x*(1.0 - qfx[IP][uf][1][0] + qfx[IP][uf][1][1])/DX[IM1]
                    - w3x*(qfx[IP][uf][2][1])/DX[IM1]
                    - w2x*(qfx[IP][uf][1][1])/DX[IM2]
                    + w3x*(1.0 - qfx[IP][uf][2][0] - qfx[IP][uf][2][1])/DX[IM2])*iadvec;
                           
	a->M.n[count] = (  w1x*qfx[IP][uf][0][1]/DX[IP1]
                    + w1x*(1.0 - qfx[IP][uf][0][0] + qfx[IP][uf][0][1])/DX[IP]
                    + w2x*(qfx[IP][uf][1][0])/DX[IP])*iadvec;
}

void iweno_hj_nug::aij_north(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k)    -= f(i-2,j,k)*(w3x*qfx[IP][uf][5][1]/DX[IM2])*iadvec
    
                + f(i+2,j,k)*(-w1x*qfx[IP][uf][3][1]/DX[IP2]
                            + w1x*(1.0 - qfx[IP][uf][3][0] - qfx[IP][uf][3][1])/DX[IP1]
                            - w2x*(qfx[IP][uf][4][1])/DX[IP1])*iadvec
                            
                + f(i+3,j,k)*(w1x*qfx[IP][uf][3][1]/DX[IP2])*iadvec;
                        
				 
	a->M.p[count] = (- w1x*(qfx[IP][uf][3][0])/DX[IP] 
                    - w2x*(1.0 - qfx[IP][uf][4][0] + qfx[IP][uf][4][1])/DX[IP] 
                    - w3x*(qfx[IP][uf][5][0])/DX[IP]
                    + w2x*(qfx[IP][uf][4][0])/DX[IM1]
                    + w3x*(1.0 - qfx[IP][uf][5][0] + qfx[IP][uf][5][1])/DX[IM1])*iadvec;
					 
	a->M.s[count] = (- w2x*qfx[IP][uf][4][0]/DX[IM1]
                    - w3x*(1.0 - qfx[IP][uf][5][0] + qfx[IP][uf][5][1])/DX[IM1]
                    - w3x*(qfx[IP][uf][5][1])/DX[IM2])*iadvec;
                           
	a->M.n[count] = (-w1x*(1.0 - qfx[IP][uf][3][0] - qfx[IP][uf][3][1])/DX[IP1]
                    + w2x*(qfx[IP][uf][4][1])/DX[IP1]
                    + w1x*(qfx[IP][uf][3][0])/DX[IP]
                    + w2x*(1.0 - qfx[IP][uf][4][0] + qfx[IP][uf][4][1])/DX[IP]
                    + w3x*qfx[IP][uf][5][0]/DX[IP])*iadvec;
}

void iweno_hj_nug::aij_east(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k)    += f(i,j-3,k)*(w3y*qfy[JP][vf][2][0]/DY[JM3])*jadvec
    
                + f(i,j-2,k)*(-w2y*qfy[JP][vf][1][1]/DY[JM2]
                           + w3y*(1.0 - qfy[JP][vf][2][0] - qfy[JP][vf][2][1])/DY[JM2]
                           - w3y*(qfy[JP][vf][2][0])/DY[JM3])*jadvec
                            
                + f(i,j+2,k)*(w1y*qfy[JP][vf][0][1]/DY[JP1])*jadvec;
                        
                           
	a->M.p[count] +=(-w1y*(1.0 - qfy[JP][vf][0][0] + qfy[JP][vf][0][1])/DY[JP] 
                    - w2y*(qfy[JP][vf][1][0])/DY[JP]
                    + w1y*(qfy[JP][vf][0][0])/DY[JM1]
                    + w2y*(1.0 - qfy[JP][vf][1][0] + qfy[JP][vf][1][1])/DY[JM1]
                    + w3y*(qfy[JP][vf][2][1])/DY[JM1])*jadvec;
					 
	a->M.e[count] = (-w1y*qfy[JP][vf][0][0]/DY[JM1]
                    - w2y*(1.0 - qfy[JP][vf][1][0] + qfy[JP][vf][1][1])/DY[JM1]
                    - w3y*(qfy[JP][vf][2][1])/DY[JM1]
                    + w2y*(-qfy[JP][vf][1][1])/DY[JM2]
                    + w3y*(1.0 - qfy[JP][vf][2][0] - qfy[JP][vf][2][1])/DY[JM2])*jadvec;
                           
	a->M.w[count] = ( w1y*qfy[JP][vf][0][1]/DY[JP1]
                    + w1y*(1.0 - qfy[JP][vf][0][0] + qfy[JP][vf][0][1])/DY[JP]
                    + w2y*(qfy[JP][vf][1][0])/DY[JP])*jadvec;
}

void iweno_hj_nug::aij_west(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k)    += f(i,j-2,k)*(-w3y*qfy[JP][vf][5][1]/DY[JM2])*jadvec
    
                + f(i,j+2,k)*(w1y*qfy[JP][vf][3][1]/DY[JP2]
                           - w1y*(1.0 - qfy[JP][vf][3][0] - qfy[JP][vf][3][1])/DY[JP1]
                           + w2y*(qfy[JP][vf][4][1])/DY[JP1])*jadvec
                            
                + f(i,j+3,k)*(-w1y*qfy[JP][vf][3][1]/DY[JP2])*jadvec;
                        
				 
	a->M.p[count] +=(-w1y*(qfy[JP][vf][3][0])/DY[JP] 
                    - w2y*(1.0 - qfy[JP][vf][4][0] + qfy[JP][vf][4][1])/DY[JP] 
                    - w3y*(qfy[JP][vf][5][0])/DY[JP]
                    + w2y*(qfy[JP][vf][4][0])/DY[JM1]
                    + w3y*(1.0 - qfy[JP][vf][5][0] + qfy[JP][vf][5][1])/DY[JM1])*jadvec;
					 
	a->M.e[count] = (-w2y*qfy[JP][vf][4][0]/DY[JM1]
                    - w3y*(1.0 - qfy[JP][vf][5][0] + qfy[JP][vf][5][1])/DY[JM1]
                    - w3y*(qfy[JP][vf][5][1])/DY[JM2])*jadvec;
                           
	a->M.w[count] = (-w1y*(1.0 - qfy[JP][vf][3][0] - qfy[JP][vf][3][1])/DY[JP1]
                    + w2y*(qfy[JP][vf][4][1])/DY[JP1]
                    + w1y*(qfy[JP][vf][3][0])/DY[JP]
                    + w2y*(1.0 - qfy[JP][vf][4][0] + qfy[JP][vf][4][1])/DY[JP]
                    + w3y*qfy[JP][vf][5][0]/DY[JP])*jadvec;
}

void iweno_hj_nug::aij_bottom(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k)    -= f(i,j,k-3)*(-w3z*qfz[KP][wf][2][0]/DZ[KM3])*kadvec
    
                + f(i,j,k-2)*( w2z*qfz[KP][wf][1][1]/DZ[KM2]
                             - w3z*(1.0 - qfz[KP][wf][2][0] - qfz[KP][wf][2][1])/DZ[KM2]
                             + w3z*(qfz[KP][wf][2][0])/DZ[KM3])*kadvec
                            
                + f(i,j,k+2)*(-w1z*qfz[KP][wf][0][1]/DZ[KP1])*kadvec;
                        
                           
	a->M.p[count]+= (-w1z*(1.0 - qfz[KP][wf][0][0] + qfz[KP][wf][0][1])/DZ[KP] 
                    - w2z*(qfz[KP][wf][1][0])/DZ[KP]
                    + w1z*(qfz[KP][wf][0][0])/DZ[KM1]
                    + w2z*(1.0 - qfz[KP][wf][1][0] + qfz[KP][wf][1][1])/DZ[KM1]
                    + w3z*(qfz[KP][wf][2][1])/DZ[KM1])*kadvec;
					 
	a->M.b[count] = (-w1z*qfz[KP][wf][0][0]/DZ[KM1]
                    - w2z*(1.0 - qfz[KP][wf][1][0] + qfz[KP][wf][1][1])/DZ[KM1]
                    - w3z*(qfz[KP][wf][2][1])/DZ[KM1]
                    - w2z*(qfz[KP][wf][1][1])/DZ[KM2]
                    + w3z*(1.0 - qfz[KP][wf][2][0] - qfz[KP][wf][2][1])/DZ[KM2])*kadvec;
                           
	a->M.t[count] = ( w1z*qfz[KP][wf][0][1]/DZ[KP1]
                    + w1z*(1.0 - qfz[KP][wf][0][0] + qfz[KP][wf][0][1])/DZ[KP]
                    + w2z*(qfz[KP][wf][1][0])/DZ[KP])*kadvec;
}

void iweno_hj_nug::aij_top(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k)    -= f(i,j,k-2)*(w3z*qfz[KP][wf][5][1]/DZ[KM2])*kadvec
    
                + f(i,j,k+2)*(- w1z*qfz[KP][wf][3][1]/DZ[KP2]
                             + w1z*(1.0 - qfz[KP][wf][3][0] - qfz[KP][wf][3][1])/DZ[KP1]
                             - w2z*(qfz[KP][wf][4][1])/DZ[KP1])*kadvec
                            
                + f(i,j,k+3)*(w1z*qfz[KP][wf][3][1]/DZ[KP2])*kadvec;
                        
				 
	a->M.p[count]+= (-w1z*(qfz[KP][wf][3][0])/DZ[KP] 
                    - w2z*(1.0 - qfz[KP][wf][4][0] + qfz[KP][wf][4][1])/DZ[KP] 
                    - w3z*(qfz[KP][wf][5][0])/DZ[KP]
                    + w2z*(qfz[KP][wf][4][0])/DZ[KM1]
                    + w3z*(1.0 - qfz[KP][wf][5][0] + qfz[KP][wf][5][1])/DZ[KM1])*kadvec;
					 
	a->M.b[count] = (-w2z*qfz[KP][wf][4][0]/DZ[KM1]
                    - w3z*(1.0 - qfz[KP][wf][5][0] + qfz[KP][wf][5][1])/DZ[KM1]
                    - w3z*(qfz[KP][wf][5][1])/DZ[KM2])*kadvec;
                           
	a->M.t[count] = (-w1z*(1.0 - qfz[KP][wf][3][0] - qfz[KP][wf][3][1])/DZ[KP1]
                    + w2z*(qfz[KP][wf][4][1])/DZ[KP1]
                    + w1z*(qfz[KP][wf][3][0])/DZ[KP]
                    + w2z*(1.0 - qfz[KP][wf][4][0] + qfz[KP][wf][4][1])/DZ[KP]
                    + w3z*qfz[KP][wf][5][0]/DZ[KP])*kadvec;
}


void iweno_hj_nug::iqmin(lexer *p,fdm *a, field& f)
{	
	q1 = (f(i-2,j,k)-f(i-3,j,k))/DX[IM3];
	q2 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q3 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q4 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
	q5 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

void iweno_hj_nug::jqmin(lexer *p,fdm *a, field& f)
{
	q1 = (f(i,j-2,k)-f(i,j-3,k))/DY[JM3];
	q2 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q3 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q4 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
	q5 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

void iweno_hj_nug::kqmin(lexer *p,fdm *a, field& f)
{
	q1 = (f(i,j,k-2)-f(i,j,k-3))/DZ[KM3];
	q2 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q3 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q4 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q5 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}

void iweno_hj_nug::iqmax(lexer *p,fdm *a, field& f)
{
    q1 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q2 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q3 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
	q4 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
	q5 = (f(i+3,j,k)-f(i+2,j,k))/DX[IP2];
}

void iweno_hj_nug::jqmax(lexer *p,fdm *a, field& f)
{
	q1 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q2 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q3 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
	q4 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
	q5 = (f(i,j+3,k)-f(i,j+2,k))/DY[JP2];
}

void iweno_hj_nug::kqmax(lexer *p,fdm *a, field& f)
{
	q1 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q2 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q3 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q4 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
	q5 = (f(i,j,k+3)-f(i,j,k+2))/DZ[KP2];
}





void iweno_hj_nug::is_south(field& b)
{

	is1 = tttw*pow( ( -b(i-3,j,k) + 3.0*b(i-2,j,k) -3.0*b(i-1,j,k) + b(i,j,k)),2.0)
		+ fourth*pow( (-b(i-3,j,k) + 5.0*b(i-2,j,k) - 7.0*b(i-1,j,k) + 3.0*b(i,j,k)),2.0);

	is2 = tttw*pow( ( -b(i-2,j,k) + 3.0*b(i-1,j,k) -3.0*b(i,j,k) + b(i+1,j,k)),2.0)
		+ fourth*pow( (-b(i-2,j,k) + b(i-1,j,k) + b(i,j,k) - b(i+1,j,k)),2.0);

	is3 = tttw*pow( (      -b(i-1,j,k) + 3.0*b(i,j,k) - 3.0*b(i+1,j,k) + b(i+2,j,k)),2.0)
		+ fourth*pow( (-3.0*b(i-1,j,k) + 7.0*b(i,j,k) - 5.0*b(i+1,j,k) + b(i+2,j,k)),2.0);
}

void iweno_hj_nug::is_north(field& b)
{

	is1 = tttw*pow( (b(i+3,j,k) - 3.0*b(i+2,j,k) + 3.0*b(i+1,j,k) - b(i,j,k)),2.0)
		+ fourth*pow( (b(i+3,j,k) - 5.0*b(i+2,j,k) + 7.0*b(i+1,j,k) - 3.0*b(i,j,k)),2.0);

	is2 = tttw*pow( (b(i+2,j,k) - 3.0*b(i+1,j,k) + 3.0*b(i,j,k) - b(i-1,j,k)),2.0)
		+ fourth*pow( (b(i+2,j,k) - b(i+1,j,k) - b(i,j,k) + b(i-1,j,k)),2.0);

	is3 = tttw*pow( (b(i+1,j,k) - 3.0*b(i,j,k) + 3.0*b(i-1,j,k) - b(i-2,j,k)),2.0)
		+ fourth*pow( (3.0*b(i+1,j,k) - 7.0*b(i,j,k) + 5.0*b(i-1,j,k) - b(i-2,j,k)),2.0);
}

void iweno_hj_nug::is_east(field& b)
{
	is1 = tttw*pow( (-b(i,j-3,k) + 3.0*b(i,j-2,k) - 3.0*b(i,j-1,k) + b(i,j,k)),2.0)
		+ fourth*pow( (-b(i,j-3,k) + 5.0*b(i,j-2,k) - 7.0*b(i,j-1,k) + 3.0*b(i,j,k)),2.0);

	is2 = tttw*pow( (-b(i,j-2,k) + 3.0*b(i,j-1,k) - 3.0*b(i,j,k) + b(i,j+1,k)),2.0)
		+ fourth*pow( (-b(i,j-2,k) + b(i,j-1,k) +b(i,j,k) - b(i,j+1,k)),2.0);

	is3 = tttw*pow( (-b(i,j-1,k) + 3.0*b(i,j,k) - 3.0*b(i,j+1,k) + b(i,j+2,k)),2.0)
		+ fourth*pow( (-3.0*b(i,j-1,k) + 7.0*b(i,j,k) - 5.0*b(i,j+1,k) + b(i,j+2,k)),2.0);
}

void iweno_hj_nug::is_west(field& b)
{
	is1 = tttw*pow( (b(i,j+3,k) - 3.0*b(i,j+2,k) + 3.0*b(i,j+1,k) - b(i,j,k)),2.0)
		+ fourth*pow( (b(i,j+3,k) - 5.0*b(i,j+2,k) + 7.0*b(i,j+1,k) - 3.0*b(i,j,k)),2.0);

	is2 = tttw*pow( (b(i,j+2,k) - 3.0*b(i,j+1,k) + 3.0*b(i,j,k) - b(i,j-1,k)),2.0)
		+ fourth*pow( (b(i,j+2,k) - b(i,j+1,k) - b(i,j,k) + b(i,j-1,k)),2.0);

	is3 = tttw*pow( (b(i,j+1,k) - 3.0*b(i,j,k) + 3.0*b(i,j-1,k) - b(i,j-2,k)),2.0)
		+ fourth*pow( (3.0*b(i,j+1,k) - 7.0*b(i,j,k) + 5.0*b(i,j-1,k) - b(i,j-2,k)),2.0);
}

void iweno_hj_nug::is_bottom(field& b)
{
	is1 = tttw*pow( (-b(i,j,k-3) + 3.0*b(i,j,k-2) -3.0*b(i,j,k-1) + b(i,j,k)),2.0)
		+ fourth*pow( (-b(i,j,k-3) + 5.0*b(i,j,k-2) - 7.0*b(i,j,k-1) + 3.0*b(i,j,k)),2.0);

	is2 = tttw*pow( (-b(i,j,k-2) + 3.0*b(i,j,k-1) -3.0*b(i,j,k) + b(i,j,k+1)),2.0)
		+ fourth*pow( (-b(i,j,k-2) + b(i,j,k-1) +b(i,j,k) - b(i,j,k+1)),2.0);

	is3 = tttw*pow( (-b(i,j,k-1) + 3.0*b(i,j,k) - 3.0*b(i,j,k+1) + b(i,j,k+2)),2.0)
		+ fourth*pow( (-3.0*b(i,j,k-1) + 7.0*b(i,j,k) - 5.0*b(i,j,k+1) + b(i,j,k+2)),2.0);
}

void iweno_hj_nug::is_top(field& b)
{
	is1 = tttw*pow( (b(i,j,k+3) - 3.0*b(i,j,k+2) + 3.0*b(i,j,k+1) - b(i,j,k)),2.0)
		+ fourth*pow( (b(i,j,k+3) - 5.0*b(i,j,k+2) + 7.0*b(i,j,k+1) - 3.0*b(i,j,k)),2.0);

	is2 = tttw*pow( (b(i,j,k+2) - 3.0*b(i,j,k+1) + 3.0*b(i,j,k) - b(i,j,k-1)),2.0)
		+ fourth*pow( (b(i,j,k+2) - b(i,j,k+1) - b(i,j,k) + b(i,j,k-1)),2.0);

	is3 = tttw*pow( (b(i,j,k+1) - 3.0*b(i,j,k) + 3.0*b(i,j,k-1) - b(i,j,k-2)),2.0)
		+ fourth*pow( (3.0*b(i,j,k+1) - 7.0*b(i,j,k) + 5.0*b(i,j,k-1) - b(i,j,k-2)),2.0);
}

void iweno_hj_nug::alpha_calc()
{
	alpha1=tenth/pow(epsi+is1,2.0);
	alpha2=sixten/pow(epsi+is2,2.0);
	alpha3=treten/pow(epsi+is3,2.0);
}

void iweno_hj_nug::weight_calc()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);
}
