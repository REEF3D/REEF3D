/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"nhflow_idiff.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"solver.h"

nhflow_idiff::nhflow_idiff(lexer* p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
}

nhflow_idiff::~nhflow_idiff()
{
}


void nhflow_idiff::diff_u(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *UHdiff, double *UH_in, double *UH, double *VH, double *WH, double alpha)
{/*
	starttime=pgc->timer();
	double visc_ddy_p,visc_ddy_m,visc_ddz_p,visc_ddz_m;
    
    LOOP
    UHdiff[IJK] = UH_in[IJK];
    
    pgc->start1(p,diff,gcval_u);

    count=0;

	LOOP // 
	{
	ev_ijk=a->eddyv(i,j,k);
	ev_ip_j_k=a->eddyv(i+1,j,k);
	ev_i_jm_k=a->eddyv(i,j-1,k);
	ev_i_jp_k=a->eddyv(i,j+1,k);
	ev_i_j_km=a->eddyv(i,j,k-1);
	ev_i_j_kp=a->eddyv(i,j,k+1);
	
	visc_ijk=a->visc(i,j,k);
	visc_ip_j_k=a->visc(i+1,j,k);
	visc_i_jm_k=a->visc(i,j-1,k);
	visc_i_jp_k=a->visc(i,j+1,k);
	visc_i_j_km=a->visc(i,j,k-1);
	visc_i_j_kp=a->visc(i,j,k+1);
	
	visc_ddy_p = 0.25*(visc_ijk+ev_ijk + visc_ip_j_k+ev_ip_j_k + visc_i_jp_k+ev_i_jp_k + a->visc(i+1,j+1,k)+a->eddyv(i+1,j+1,k));    
    
	visc_ddy_m = 0.25*(visc_i_jm_k+ev_i_jm_k  +a->visc(i+1,j-1,k)+a->eddyv(i+1,j-1,k) + visc_ijk+ev_ijk + visc_ip_j_k+ev_ip_j_k);
    
	visc_ddz_p = 0.25*(visc_ijk+ev_ijk + visc_ip_j_k+ev_ip_j_k + visc_i_j_kp+ev_i_j_kp + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1));
    
	visc_ddz_m = 0.25*(visc_i_j_km+ev_i_j_km + a->visc(i+1,j,k-1)+a->eddyv(i+1,j,k-1) + visc_ijk+ev_ijk + visc_ip_j_k+ev_ip_j_k);

	a->M.p[count] =  2.0*(visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP1]*p->DXP[IP])
				   + 2.0*(visc_ijk+ev_ijk)/(p->DXN[IP]*p->DXP[IP])
				   + visc_ddy_p/(p->DYP[JP]*p->DYN[JP])
				   + visc_ddy_m/(p->DYP[JM1]*p->DYN[JP])
				   + visc_ddz_p/(p->DZP[KP]*p->DZN[KP])
				   + visc_ddz_m/(p->DZP[KM1]*p->DZN[KP])
				   + CPOR1/(alpha*p->dt);
				  
	a->rhsvec.V[count] +=  ((v(i+1,j,k)-v(i,j,k))*visc_ddy_p - (v(i+1,j-1,k)-v(i,j-1,k))*visc_ddy_m)/(p->DXP[IP]*p->DYN[JP])
						 + ((w(i+1,j,k)-w(i,j,k))*visc_ddz_p - (w(i+1,j,k-1)-w(i,j,k-1))*visc_ddz_m)/(p->DXP[IP]*p->DZN[KP])

						 + (CPOR1*u_in(i,j,k))/(alpha*p->dt);
                         
	 
	 a->M.s[count] = -2.0*(visc_ijk+ev_ijk)/(p->DXN[IP]*p->DXP[IP]);
	 a->M.n[count] = -2.0*(visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP1]*p->DXP[IP]);
	 
	 a->M.e[count] = -visc_ddy_m/(p->DYP[JM1]*p->DYN[JP]);
	 a->M.w[count] = -visc_ddy_p/(p->DYP[JP]*p->DYN[JP]);
	 
	 a->M.b[count] = -visc_ddz_m/(p->DZP[KM1]*p->DZN[KP]);
	 a->M.t[count] = -visc_ddz_p/(p->DZP[KP]*p->DZN[KP]);
	 
     ++count;
	 }
     
    n=0;
    LOOP
	{
        if(p->wet[IJ]==1 && p->deep[IJ]==1 && d->breaking(i,j)==0)
        {
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            
            d->M.p[n]  =  (CPORNH*PORVALNH)/(p->W1*p->DXP[IP]*p->DXN[IP])
                        + (CPORNHm*PORVALNHm)/(p->W1*p->DXP[IM1]*p->DXN[IP])
                        
                        + (CPORNH*PORVALNH)/(p->W1*p->DYP[JP]*p->DYN[JP])*p->y_dir
                        + (CPORNHm*PORVALNHm)/(p->W1*p->DYP[JM1]*p->DYN[JP])*p->y_dir
                        
                        + (sigxyz2*CPORNH*PORVALNH)/(p->W1*p->DZP[KM1]*p->DZN[KP])
                        + (sigxyz2*CPORNHm*PORVALNHm)/(p->W1*p->DZP[KM1]*p->DZN[KM1]);


            d->M.n[n] = -(CPORNH*PORVALNH)/(p->W1*p->DXP[IP]*p->DXN[IP]);
            d->M.s[n] = -(CPORNHm*PORVALNHm)/(p->W1*p->DXP[IM1]*p->DXN[IP]);

            d->M.w[n] = -(CPORNH*PORVALNH)/(p->W1*p->DYP[JP]*p->DYN[JP])*p->y_dir;
            d->M.e[n] = -(CPORNHm*PORVALNHm)/(p->W1*p->DYP[JM1]*p->DYN[JP])*p->y_dir;

            d->M.t[n] = -(sigxyz2*CPORNH*PORVALNH)/(p->W1*p->DZP[KM1]*p->DZN[KP])     
                        - CPORNH*PORVALNH*p->sigxx[FIJK]/(p->W1*(p->DZN[KP]+p->DZN[KM1]));
                        
            d->M.b[n] = -(sigxyz2*CPORNHm*PORVALNHm)/(p->W1*p->DZP[KM1]*p->DZN[KM1]) 
                        + CPORNH*PORVALNH*p->sigxx[FIJK]/(p->W1*(p->DZN[KP]+p->DZN[KM1]));
            
            
            if(p->D33==0)
            d->rhsvec.V[n] +=  CPORNH*PORVALNH*2.0*p->sigx[FIJK]*(P[FIp1JKp1] - P[FIm1JKp1] - P[FIp1JKm1] + P[FIm1JKm1])
                            /(p->W1*(p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))
                        
                            + CPORNH*PORVALNH*2.0*p->sigy[FIJK]*(P[FIJp1Kp1] - P[FIJm1Kp1] - P[FIJp1Km1] + P[FIJm1Km1])
                            /(p->W1*(p->DYP[JP]+p->DYP[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;

            if(p->D33==1)
            {
            d->M.sb[n] = -CPORNH*PORVALNH*2.0*p->sigx[FIJK]/(p->W1*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]));
            d->M.st[n] =  CPORNH*PORVALNH*2.0*p->sigx[FIJK]/(p->W1*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]));
            d->M.nb[n] =  CPORNH*PORVALNH*2.0*p->sigx[FIJK]/(p->W1*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]));
            d->M.nt[n] = -CPORNH*PORVALNH*2.0*p->sigx[FIJK]/(p->W1*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]));
            
            d->M.eb[n] = -CPORNH*PORVALNH*2.0*p->sigy[FIJK]/(p->W1*(p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
            d->M.et[n] =  CPORNH*PORVALNH*2.0*p->sigy[FIJK]/(p->W1*(p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
            d->M.wb[n] =  CPORNH*PORVALNH*2.0*p->sigy[FIJK]/(p->W1*(p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
            d->M.wt[n] = -CPORNH*PORVALNH*2.0*p->sigy[FIJK]/(p->W1*(p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
            }
        }
        
        if(p->wet[IJ]==0 || p->deep[IJ]==0 || p->flag7[FIJK]<0 || d->breaking(i,j)==1)
        {
        d->M.p[n]  =  1.0;


        d->M.n[n] = 0.0;
        d->M.s[n] = 0.0;

        d->M.w[n] = 0.0;
        d->M.e[n] = 0.0;

        d->M.t[n] = 0.0;
        d->M.b[n] = 0.0;
        
        d->rhsvec.V[n] =  0.0;
        }
	
	++n;
	}
    
    
    
    
    
    
    
    
    n=0;
	ULOOP
	{
		if(p->flag1[Im1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.s[n]*u(i-1,j,k);
		a->M.s[n] = 0.0;
		}
		
		if(p->flag1[Ip1JK]<0)
		{
		a->rhsvec.V[n] -= a->M.n[n]*u(i+1,j,k);
		a->M.n[n] = 0.0;
		}
		
		if(p->flag1[IJm1K]<0)
		{
		a->rhsvec.V[n] -= a->M.e[n]*u(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag1[IJp1K]<0)
		{
		a->rhsvec.V[n] -= a->M.w[n]*u(i,j+1,k);
		a->M.w[n] = 0.0;
		}
		
		if(p->flag1[IJKm1]<0)
		{
		a->rhsvec.V[n] -= a->M.b[n]*u(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag1[IJKp1]<0)
		{
		a->rhsvec.V[n] -= a->M.t[n]*u(i,j,k+1);
		a->M.t[n] = 0.0;
		}

	++n;
	}
	
	psolv->start(p,a,pgc,diff,a->rhsvec,1);
    
	
    pgc->start1(p,diff,gcval_u);
    
    
	time=pgc->timer()-starttime;
	p->uiter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"udiffiter: "<<p->uiter<<"  udifftime: "<<setprecision(3)<<time<<endl;*/
}
