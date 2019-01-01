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

#include"fnpf_fg_laplace_cds2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"

fnpf_fg_laplace_cds2::fnpf_fg_laplace_cds2() 
{
}

fnpf_fg_laplace_cds2::~fnpf_fg_laplace_cds2()
{
}

void fnpf_fg_laplace_cds2::start(lexer* p, fdm *a, ghostcell *pgc, solver *psolv, field &f)
{
	n=0;
    LOOP
	{
	a->M.p[n]  =  1.0/(p->DXP[IP]*p->DXN[IP])*p->x_dir 
                + 1.0/(p->DXP[IM1]*p->DXN[IP])*p->x_dir 
                
                + 1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir 
                + 1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir 
                
                + 1.0/(p->DZP[KP]*p->DZN[KP])*p->z_dir
                + 1.0/(p->DZP[KM1]*p->DZN[KP])*p->z_dir;

   	a->M.n[n] = -1.0/(p->DXP[IP]*p->DXN[IP])*p->x_dir;
	a->M.s[n] = -1.0/(p->DXP[IM1]*p->DXN[IP])*p->x_dir;

	a->M.w[n] = -1.0/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
	a->M.e[n] = -1.0/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;

	a->M.t[n] = -1.0/(p->DZP[KP]*p->DZN[KP])*p->z_dir;
	a->M.b[n] = -1.0/(p->DZP[KM1]*p->DZN[KP])*p->z_dir;

	a->rhsvec.V[n] = 0.0;
	
	++n;
	}
    
    
    n=0;
	LOOP
	{
        if(p->flag4[IJK]>0)
        {
    
            if(p->flag4[Im1JK]<0)
            {
            a->rhsvec.V[n] -= a->M.s[n]*f(i-1,j,k);
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0)
            {
            a->rhsvec.V[n] -= a->M.n[n]*f(i+1,j,k);
            a->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0)
            {
            a->rhsvec.V[n] -= a->M.e[n]*f(i,j-1,k);
            a->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0)
            {
            a->rhsvec.V[n] -= a->M.w[n]*f(i,j+1,k);
            a->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0)
            {
            a->rhsvec.V[n] -= a->M.b[n]*f(i,j,k-1);
            a->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0)
            {
            a->rhsvec.V[n] -= a->M.t[n]*f(i,j,k+1);
            a->M.t[n] = 0.0;
            }
        }

	++n;
	}
    
    
    double starttime=pgc->timer();
    psolv->start(p,a,pgc,a->Fi,a->xvec,a->rhsvec,5,250,p->N44);
    double endtime=pgc->timer();
    pgc->start4(p,a->Fi,250);
    
    p->poissoniter=p->solveriter;
    p->poissontime=endtime-starttime;
	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	cout<<"Fi_iter: "<<p->solveriter<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;
    
    /*
    n=0;
	LOOP
	{
		if(fabs(a->phi(i,j,k))<p->dx*1.6)
        {
            
        a->M.p[n]=1.0e30;
        a->rhsvec.V[n] = -1.0e30*f(i,j,k);
        }

	++n;
	}
    */
    
    
    /*
    double ddx,ddy,ddz;
    
    n=0;
    LOOP
	{
    ddx = 0.5*pow(p->DXP[IP],2.0) + 0.5*pow(p->DXP[IM1],2.0);
    ddy = 0.5*pow(p->DYP[JP],2.0) + 0.5*pow(p->DYP[JM1],2.0);
    ddz = 0.5*pow(p->DZP[KP],2.0) + 0.5*pow(p->DZP[KM2],2.0);

	a->M.p[n] =  2.0/ddx*p->x_dir + 2.0/ddy*p->y_dir + 2.0/ddz*p->z_dir;

   	a->M.n[n] = -1.0/ddx*p->x_dir;
	a->M.s[n] = -1.0/ddx*p->x_dir;

	a->M.w[n] = -1.0/ddy*p->y_dir;
	a->M.e[n] = -1.0/ddy*p->y_dir;

	a->M.t[n] = -1.0/ddz*p->z_dir;
	a->M.b[n] = -1.0/ddz*p->z_dir;
    

	a->rhsvec.V[n] = 0.0;
	++n;
	}
    */
}

