/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"outflow_pressure.h"
#include"poisson_pcorr.h"
#include"bc_noflux.h"
#include<mpi.h>
#include"lexer.h"
#include"fdm.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_df.h"
#include"density_sf.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
#include"density_rheo.h"
#include"density_pst.h"

poisson_pcorr::poisson_pcorr(lexer *p, heat *&pheat, concentration *&pconc) 
{
    if(p->F80==0 && p->F300==0 && p->W90==0)
    {
        if(p->W30==0 && p->C10==0 && p->H10==0)
        {
        if(p->X10==0 && p->Q10==0)
        pd = new density_f(p);
        
        if(p->X10==0 && p->Q10>=1)
        pd = new density_f(p);
        
        if(p->X10==1)  
        pd = new density_df(p);        
        }
        
        if(p->H10==0 && p->W30==1)
        pd = new density_comp(p);
        
        if(p->H10>0 && p->C10==0)
        pd = new density_heat(p,pheat);
        
        if(p->C10>0 && p->H10==0)
        pd = new density_conc(p,pconc);
    }
    
    if(p->F80>0 && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0)
	pd = new density_vof(p);
    
    if(p->F30>0 && p->H10==0 && p->W30==0  && p->F300==0 && p->W90>0)
    pd = new density_rheo(p);
    
    if(p->F300>=1)
    pd = new density_rheo(p);
}

poisson_pcorr::~poisson_pcorr()
{
}

void poisson_pcorr::start(lexer* p, fdm *a, field &press)
{
    a->M.reset();

    n=0;
    LOOP
	{
	a->M.p[n]  =  (CPOR1*PORVAL1)/(pd->roface(p,a,1,0,0)*p->DXP[IP]*p->DXN[IP])
                + (CPOR1m*PORVAL1m)/(pd->roface(p,a,-1,0,0)*p->DXP[IM1]*p->DXN[IP])
                
                + (CPOR2*PORVAL2)/(pd->roface(p,a,0,1,0)*p->DYP[JP]*p->DYN[JP])*p->y_dir
                + (CPOR2m*PORVAL2m)/(pd->roface(p,a,0,-1,0)*p->DYP[JM1]*p->DYN[JP])*p->y_dir
                
                + (CPOR3*PORVAL3)/(pd->roface(p,a,0,0,1)*p->DZP[KP]*p->DZN[KP])
                + (CPOR3m*PORVAL3m)/(pd->roface(p,a,0,0,-1)*p->DZP[KM1]*p->DZN[KP]);


   	a->M.n[n] = -(CPOR1*PORVAL1)/(pd->roface(p,a,1,0,0)*p->DXP[IP]*p->DXN[IP]);
	a->M.s[n] = -(CPOR1m*PORVAL1m)/(pd->roface(p,a,-1,0,0)*p->DXP[IM1]*p->DXN[IP]);

	a->M.w[n] = -(CPOR2*PORVAL2)/(pd->roface(p,a,0,1,0)*p->DYP[JP]*p->DYN[JP])*p->y_dir;
	a->M.e[n] = -(CPOR2m*PORVAL2m)/(pd->roface(p,a,0,-1,0)*p->DYP[JM1]*p->DYN[JP])*p->y_dir;

	a->M.t[n] = -(CPOR3*PORVAL3)/(pd->roface(p,a,0,0,1)*p->DZP[KP]*p->DZN[KP]);
	a->M.b[n] = -(CPOR3m*PORVAL3m)/(pd->roface(p,a,0,0,-1)*p->DZP[KM1]*p->DZN[KP]);
	
	++n;
	}
    

    // boundaries where the normal velocity is prescribed and not corrected by the projection
    // (walls, lid, bed, solids, inflow, wave generation, velocity patches): Neumann for the
    // pressure correction, i.e. the coefficient of the ghost cell is dropped from the diagonal
    bc_noflux_mask(p,noflux,BC_NOFLUX_WALLS|BC_NOFLUX_INFLOW);

    int ndirichlet=0;          // remaining boundary faces with a fixed pressure correction
    int pin_n=-1;              // row and side of one converted face, to pin the level if needed
    int pin_side=0;
    double pin_coef=0.0;

    n=0;
	LOOP
	{
        int nf = noflux[IJK];
        double *coef[6] = {&a->M.s[n],&a->M.w[n],&a->M.e[n],&a->M.n[n],&a->M.b[n],&a->M.t[n]};   // cs = 1..6
        bool ghost[6] = {p->flag4[Im1JK]<0 && (i+p->origin_i>0 || p->periodic1==0),
                         p->flag4[IJp1K]<0 && (j+p->origin_j<p->gknoy-1 || p->periodic2==0) && p->j_dir==1,
                         p->flag4[IJm1K]<0 && (j+p->origin_j>0 || p->periodic2==0) && p->j_dir==1,
                         p->flag4[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0),
                         p->flag4[IJKm1]<0 && (k+p->origin_k>0 || p->periodic3==0),
                         p->flag4[IJKp1]<0 && (k+p->origin_k<p->gknoz-1 || p->periodic3==0)};

        for(int cs=1; cs<=6; ++cs)
        {
            if(!ghost[cs-1])
            continue;

            if(nf & (1<<(cs-1)))
            {
                // prefer a lid/top face (cs 6) for pinning, otherwise the first one
                if(pin_n<0 || (cs==6 && pin_side!=6))
                {
                pin_n=n;
                pin_side=cs;
                pin_coef=*coef[cs-1];
                }

                a->M.p[n] += *coef[cs-1];
                *coef[cs-1] = 0.0;
            }
            else
            ++ndirichlet;
        }
    ++n;
    }

    // closed domain (no inflow/outflow/open boundary anywhere): the Neumann problem is singular,
    // so one face (on the lowest rank that has one) keeps the fixed value as a pressure level
    int nd_glob=0;
    MPI_Allreduce(&ndirichlet,&nd_glob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(nd_glob==0)
    {
        int myrank = (pin_n>=0) ? p->mpirank : 1<<30;
        int pinrank = 1<<30;
        MPI_Allreduce(&myrank,&pinrank,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

        if(p->mpirank==pinrank && pin_n>=0)
        {
        // undo the conversion: the face below sets M.x*press(ghost) into the rhs (press ghost = 0)
        double *coef[6] = {&a->M.s[pin_n],&a->M.w[pin_n],&a->M.e[pin_n],&a->M.n[pin_n],&a->M.b[pin_n],&a->M.t[pin_n]};
        a->M.p[pin_n] -= pin_coef;
        *coef[pin_side-1] = pin_coef;
        }
    }

    n=0;
	LOOP
	{
        // inflow
		if(p->flag4[Im1JK]<0 && (i+p->origin_i>0 || p->periodic1==0))
		{
		a->rhsvec.V[n] -= a->M.s[n]*press(i-1,j,k);
		a->M.s[n] = 0.0;
		}
        
        /*
        if(p->flag4[Im1JK]<0 &&  p->IO[Im1JK]==1)
		{
        pval=(p->fsfin - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
        //cout<<"FSFIN: "<<p->fsfin<<endl;
		a->rhsvec.V[n] -= a->M.s[n]*(-a->press(i,j,k)+pval);
		a->M.s[n] = 0.0;
		}*/
        
        // AWA inflow
        /*if(p->flag4[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0) && (p->IO[Ip1JK]==2 && p->B90==1 && p->B99>2))
        {
        pval=(p->fsfout - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
        
        a->rhsvec.V[n] += a->M.s[n]*(-a->press(i,j,k)+pval);
        
        a->M.s[n] = 0.0;
        }*/
        
        /*
        if(p->flag4[Im1JK]<0 && (i+p->origin_i>0 || p->periodic1==0) && p->IO[Im1JK]==1)
		{
             pval=a->press(i,j,k);
             
		//a->rhsvec.V[n] -= a->M.s[n]*(-a->press(i,j,k)+pval);
        a->rhsvec.V[n] -= a->M.s[n]*press(i-1,j,k);
		a->M.s[n] = 0.0;
		}*/
		
        // outflow (IO 2): the ghost pressure is outflow_pressure(), the value pressure_io has put
        // into the ghost cells before the projection, so the correction there is its difference
        if(p->flag4[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0) && p->IO[Ip1JK]==2)
		{
		a->rhsvec.V[n] -= a->M.n[n]*(outflow_pressure(p,a,i,j,k) - a->press(i+1,j,k));
		a->M.n[n] = 0.0;
		}

		if(p->flag4[Ip1JK]<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0))
		{
		a->rhsvec.V[n] -= a->M.n[n]*press(i+1,j,k);
		a->M.n[n] = 0.0;
		}
        
    // ----
		if(p->flag4[IJm1K]<0 && (j+p->origin_j>0 || p->periodic2==0) && p->j_dir==1)
		{
		a->rhsvec.V[n] -= a->M.e[n]*press(i,j-1,k);
		a->M.e[n] = 0.0;
		}
		
		if(p->flag4[IJp1K]<0 && (j+p->origin_j<p->gknoy-1 || p->periodic2==0) && p->j_dir==1)
		{
		a->rhsvec.V[n] -= a->M.w[n]*press(i,j+1,k);
		a->M.w[n] = 0.0;
		}
		
		if(p->flag4[IJKm1]<0 && (k+p->origin_k>0 || p->periodic3==0))
		{
		a->rhsvec.V[n] -= a->M.b[n]*press(i,j,k-1);
		a->M.b[n] = 0.0;
		}
		
		if(p->flag4[IJKp1]<0 && (k+p->origin_k<p->gknoz-1 || p->periodic3==0))
		{
		a->rhsvec.V[n] -= a->M.t[n]*press(i,j,k+1);
		a->M.t[n] = 0.0;
		}
	++n;
	}
    
  
}
