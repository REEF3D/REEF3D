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

#include"ptf_laplace_cds2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"

ptf_laplace_cds2::ptf_laplace_cds2(lexer *p, fdm *a, ghostcell *pgc) : bc(p)
{
    // bc ini
    SLICELOOP4
	bc(i,j) = 0;
    
    pgc->gcsl_start4int(p,bc,50);
    
    if(p->B98>=3)
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];

    bc(i-1,j) = 1;
    }
    
    if(p->B99>=3)
    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];

    bc(i+1,j) = 2;
    }
}

ptf_laplace_cds2::~ptf_laplace_cds2()
{
}

void ptf_laplace_cds2::start(lexer* p, fdm *a, ghostcell *pgc, solver *psolv, field &f)
{
	n=0;
    FLUIDLOOP
	{
        
    if(p->flag4[IJK]>0)
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
    }
    
    
    if(p->flag4[IJK]<0)
    {
	a->M.p[n]  =  1.0;

   	a->M.n[n] = 0.0;
	a->M.s[n] = 0.0;

	a->M.w[n] = 0.0;
	a->M.e[n] = 0.0;

	a->M.t[n] = 0.0;
	a->M.b[n] = 0.0;

	a->rhsvec.V[n] = 0.0;
    }
	
	++n;
	}
    
    
    n=0;
	FLUIDLOOP
	{
        if(p->flag4[IJK]>0)
        {
            
            // south
            if(p->flag4[Im1JK]<AIR && bc(i-1,j)==0)
            {
            a->M.p[n] += a->M.s[n];
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[Im1JK]<AIR && bc(i-1,j)==1)
            {
            a->rhsvec.V[n] += a->M.s[n]*a->u(i-1,j,k)*p->DXP[IP];
            a->M.p[n] += a->M.s[n];
            a->M.s[n] = 0.0;
            }
            
            // north
            if(p->flag4[Ip1JK]<AIR && bc(i+1,j)==0)
            {
            a->M.p[n] += a->M.n[n];
            a->M.n[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<AIR && bc(i+1,j)==2)
            {
            a->rhsvec.V[n] -= a->M.n[n]*a->u(i+1,j,k)*p->DXP[IP1];
            a->M.p[n] += a->M.n[n];
            a->M.n[n] = 0.0;
            }
            
            // east
            if(p->flag4[IJm1K]<AIR)
            {
            a->M.p[n] += a->M.e[n];
            a->M.e[n] = 0.0;
            }
            
            // west
            if(p->flag4[IJp1K]<AIR)
            {
            a->M.p[n] += a->M.w[n];
            a->M.w[n] = 0.0;
            }
            
            // FSFBC
            if(p->flag4[IJKp1]==AIR)
            {
            // old
            //a->rhsvec.V[n] -= a->M.t[n]*f(i,j,k+1);
            //a->M.t[n] = 0.0;
            
            // new
            double lsv0,lsv1;
            
            lsv0 = fabs(a->phi(i,j,k));
            lsv1 = fabs(a->phi(i,j,k+1));

            lsv0 = fabs(lsv0)>1.0e-6?lsv0:1.0e20;
            
   
            a->rhsvec.V[n] -= a->M.t[n]*f(i,j,k+1)*(1.0 + lsv1/lsv0);
            a->M.p[n] -= a->M.t[n]*lsv1/lsv0;
            a->M.t[n] = 0.0;
            }
    
            // KBEDBC
            if(p->flag4[IJKm1]<AIR)
            {
            a->M.p[n] += a->M.b[n];
            a->M.b[n] = 0.0;
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
	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"Fi_iter: "<<p->solveriter<<"  Fi_time: "<<setprecision(3)<<p->poissontime<<endl;

}

