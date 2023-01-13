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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"idiff_IMEX_2D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"
#include"density_f.h"
#include"density_df.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
#include"density_rheo.h"

idiff_IMEX_2D::idiff_IMEX_2D(lexer* p){}

idiff_IMEX_2D::idiff_IMEX_2D(lexer* p, heat* pheat, concentration* pconc)
{
    if((p->F80==0||p->A10==5) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && (p->X10==0 || p->X13!=2))
	pd = new density_f(p);
    
    if((p->F80==0||p->A10==5) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && (p->X10==1 || p->X13!=2))  
	pd = new density_df(p);

	if(p->F80==0 && p->H10==0 && p->W30==1 && p->W90==0)
	pd = new density_comp(p);
	
	if(p->F80==0 && p->H10>0 && p->W90==0)
	pd = new density_heat(p,pheat);
	
	if(p->F80==0 && p->C10>0 && p->W90==0)
	pd = new density_conc(p,pconc);
    
    if(p->F80>0 && p->H10==0 && p->W30==0 && p->W90==0)
	pd = new density_vof(p);
    
    if(p->F30>0 && p->H10==0 && p->W30==0 && p->W90>0)
    pd = new density_rheo(p);

    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
}

idiff_IMEX_2D::~idiff_IMEX_2D()
{
}

void idiff_IMEX_2D::diff_u(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	double visc_ddz_p,visc_ddz_m;
    
    pgc->start1(p,u,gcval_u);

	count=0;
    if(p->i_dir==1)
    {
        ULOOP
        {
            a->rhsvec.V[count] = 0.0; 
            
            ev_ijk=a->eddyv(i,j,k);
            ev_ip_j_k=a->eddyv(i+1,j,k);
            ev_i_j_km=a->eddyv(i,j,k-1);
            ev_i_j_kp=a->eddyv(i,j,k+1);
            
            visc_ijk=a->visc(i,j,k);
            visc_ip_j_k=a->visc(i+1,j,k);
            visc_i_j_km=a->visc(i,j,k-1);
            visc_i_j_kp=a->visc(i,j,k+1);
            
            visc_ddz_p = (visc_ijk+ev_ijk + visc_ip_j_k+ev_ip_j_k + visc_i_j_kp+ev_i_j_kp + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1))*0.25;
            visc_ddz_m = (a->visc(i,j,k-1)+a->eddyv(i,j,k-1) + a->visc(i+1,j,k-1)+a->eddyv(i+1,j,k-1) + visc_ijk+ev_ijk + visc_ip_j_k+ev_ip_j_k)*0.25;
                
            a->M.p[count] =  2.0*(visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP1]*p->DXP[IP])   // changed
                           + 2.0*(visc_ijk+ev_ijk)/(p->DXN[IP]*p->DXP[IP])          // changed
                           + visc_ddz_p/(p->DZP[KP]*p->DZN[KP])
                           + visc_ddz_m/(p->DZP[KM1]*p->DZN[KP])
                           + CPOR1/(alpha*p->dt);
                          
             a->M.s[count] = -2.0*(visc_ijk+ev_ijk)/(p->DXN[IP]*p->DXP[IP]);          // changed
             a->M.n[count] = -2.0*(visc_ip_j_k+ev_ip_j_k)/(p->DXN[IP1]*p->DXP[IP]);   // changed
             
             a->M.b[count] = -visc_ddz_m/(p->DZP[KM1]*p->DZN[KP]);
             a->M.t[count] = -visc_ddz_p/(p->DZP[KP]*p->DZN[KP]);
            
             a->rhsvec.V[count] += 
                        ((w(i+1,j,k)-w(i,j,k))*visc_ddz_p - (w(i+1,j,k-1)-w(i,j,k-1))*visc_ddz_m)/(p->DXP[IP]*p->DZN[KP])  
                        + a->F(i,j,k) - PORVAL1*(a->press(i+1,j,k)-a->press(i,j,k))/(p->DXP[IP]*pd->roface(p,a,1,0,0));
             
             ++count;
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
        
        psolv->start(p,a,pgc,u,a->rhsvec,1);
    }
	
    pgc->start1(p,u,gcval_u);
    
	time=pgc->timer()-starttime;
	p->uiter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && p->count%p->P12==0)
	cout<<"udiffiter: "<<p->uiter<<"  udifftime: "<<setprecision(3)<<time<<endl;
}


void idiff_IMEX_2D::diff_v(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
}

void idiff_IMEX_2D::diff_w(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	starttime=pgc->timer();
	
	double visc_ddx_p,visc_ddx_m;
	
    pgc->start3(p,w,gcval_w);

	count=0;
    if(p->k_dir==1)
    {
        WLOOP
        {
            a->rhsvec.V[count] = 0.0; 
            
            ev_ijk=a->eddyv(i,j,k);
            ev_im_j_k=a->eddyv(i-1,j,k);
            ev_ip_j_k=a->eddyv(i+1,j,k);
            ev_i_j_kp=a->eddyv(i,j,k+1);
            
            visc_ijk=a->visc(i,j,k);
            visc_im_j_k=a->visc(i-1,j,k);
            visc_ip_j_k=a->visc(i+1,j,k);
            visc_i_j_kp=a->visc(i,j,k+1);
            
            visc_ddx_p = (visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp + visc_ip_j_k+ev_ip_j_k + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1))*0.25;
            visc_ddx_m = (visc_im_j_k+ev_im_j_k + a->visc(i-1,j,k+1)+a->eddyv(i-1,j,k+1) + visc_ijk+ev_ijk + visc_i_j_kp+ev_i_j_kp)*0.25;
            
            a->M.p[count] = 2.0*(visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP1]*p->DZP[KP])       // changed
                          + 2.0*(visc_ijk+ev_ijk)/(p->DZN[KP]*p->DZP[KP])              // changed
                          + visc_ddx_p/(p->DXP[IP]*p->DXN[IP])
                          + visc_ddx_m/(p->DXP[IM1]*p->DXN[IP])
                          + CPOR3/(alpha*p->dt);
                          
             a->M.s[count] = -visc_ddx_m/(p->DXP[IM1]*p->DXN[IP]);
             a->M.n[count] = -visc_ddx_p/(p->DXP[IP]*p->DXN[IP]);
             
             a->M.b[count] = -2.0*(visc_ijk+ev_ijk)/(p->DZN[KP]*p->DZP[KP]);           // changed
             a->M.t[count] = -2.0*(visc_i_j_kp+ev_i_j_kp)/(p->DZN[KP1]*p->DZP[KP]);    // changed
            
            a->rhsvec.V[count] += 
                  ((u(i,j,k+1)-u(i,j,k))*visc_ddx_p - (u(i-1,j,k+1)-u(i-1,j,k))*visc_ddx_m)/(p->DZP[KP]*p->DXN[IP])  
                + a->H(i,j,k) - PORVAL3*(a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1));

            ++count;
        }
        
        n=0;
        WLOOP
        {
            if(p->flag3[Im1JK]<0)
            {
            a->rhsvec.V[n] -= a->M.s[n]*w(i-1,j,k);
            a->M.s[n] = 0.0;
            }
            
            if(p->flag3[Ip1JK]<0)
            {
            a->rhsvec.V[n] -= a->M.n[n]*w(i+1,j,k);
            a->M.n[n] = 0.0;
            }
            
            if(p->flag3[IJm1K]<0)
            {
            a->rhsvec.V[n] -= a->M.e[n]*w(i,j-1,k);
            a->M.e[n] = 0.0;
            }
            
            if(p->flag3[IJp1K]<0)
            {
            a->rhsvec.V[n] -= a->M.w[n]*w(i,j+1,k);
            a->M.w[n] = 0.0;
            }
            
            if(p->flag3[IJKm1]<0)
            {
            a->rhsvec.V[n] -= a->M.b[n]*w(i,j,k-1);
            a->M.b[n] = 0.0;
            }
            
            if(p->flag3[IJKp1]<0)
            {
            a->rhsvec.V[n] -= a->M.t[n]*w(i,j,k+1);
            a->M.t[n] = 0.0;
            }

        ++n;
        }
        
        
        psolv->start(p,a,pgc,w,a->rhsvec,3);
    }
    
    
    pgc->start3(p,w,gcval_w);
	
	time=pgc->timer()-starttime;
	p->witer=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"wdiffiter: "<<p->witer<<"  wdifftime: "<<setprecision(3)<<time<<endl;
}


void idiff_IMEX_2D::diff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, double, double){}
void idiff_IMEX_2D::idiff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, double, double){}

