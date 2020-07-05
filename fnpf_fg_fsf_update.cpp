/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"fnpf_fg_fsf_update.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"onephase.h"
#include"slice.h"

fnpf_fg_fsf_update::fnpf_fg_fsf_update(lexer *p, fdm *a, ghostcell *pgc) 
{
    gcval_u = 10;
    gcval_v = 11;
    gcval_w = 12;
}

fnpf_fg_fsf_update::~fnpf_fg_fsf_update()
{
    
}

void fnpf_fg_fsf_update::fsfupdate(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, onephase *poneph, slice &eta)
{
    // update phi
    FLUIDLOOP
    a->phi(i,j,k) = eta(i,j) + p->phimean - p->pos_z();

    pgc->start4(p,a->phi,50);
    
    // update onephase
    poneph->update(p,a,pgc,pflow);
}

void fnpf_fg_fsf_update::etaloc(lexer *p, fdm *a, ghostcell *pgc)
{
    // find k location for eta, i.e. zero level set
    SLICELOOP4
    a->etaloc(i,j) = 0;
    
    LOOP
    a->etaloc(i,j) = MAX(a->etaloc(i,j),k);
}

void fnpf_fg_fsf_update::etaloc_sig(lexer *p, fdm *a, ghostcell *pgc)
{
}

void fnpf_fg_fsf_update::fsfbc_sig(lexer *p, fdm *a, ghostcell *pgc, slice &Fifsf, field &Fi)
{
}

void fnpf_fg_fsf_update::fsfbc(lexer *p, fdm *a, ghostcell *pgc, slice &Fifsf, field &Fi)
{
    AIRLOOP
    Fi(i,j,k)=0.0; 
    
    
    double lsv,val,dx;
    double nx,ny,nz,dnorm;
    double xp,yp,zp;
    
    double fac = 0.6;
    
    
    // insert DFSFBC into current Fi field
    AIRLOOP
    {
    lsv = a->phi(i,j,k);
    
        if(lsv>-3.1*(1.0/3.0)*(p->DXP[IP] + p->DYP[JP] + p->DZP[KP]))
        {
         nx = (a->phi(i+1,j,k)-a->phi(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
         ny = (a->phi(i,j+1,k)-a->phi(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]);
         nz = (a->phi(i,j,k+1)-a->phi(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);  

        dnorm = sqrt(nx*nx + ny*ny + nz*nz);
        
        nx/=dnorm;
        ny/=dnorm;
        nz/=dnorm;
        
        xp = p->pos_x() + nx*(fabs(lsv)+fac*p->DXP[IP]);
        yp = p->pos_y() + ny*(fabs(lsv)+fac*p->DYP[JP]);
        zp = p->pos_z() + nz*(fabs(lsv)+fac*p->DZP[KP]);
        
        dx = sqrt(pow(nx*p->DXP[IP],2.0) + pow(ny*p->DYP[JP],2.0) + pow(p->DZP[KP],dx));
        
        val = p->ccipol4(Fi, xp, yp, zp);
        
        Fi(i,j,k) =  (Fifsf(i,j)-val)*(fabs(lsv)/(fac*dx)) + Fifsf(i,j);
        
        //cout<<"Fi_epol: "<<Fi(i,j,k)<<" Fifsf: "<<Fifsf(i,j)<<" val: "<<val<<endl;
        }
    }
    
    /*
    // new: vertical extrapolation
    double lsv1,lsv2,lsv3;
    FILOOP4
    {
    lsv = fabs(a->phi(i,j,k))>1.0e-10?a->phi(i,j,k):1.0e10;
    val = Fi(i,j,k);
    
    lsv1 = a->phi(i,j,k+1);
    lsv2 = a->phi(i,j,k+2);
    lsv3 = a->phi(i,j,k+3);
        
    Fi(i,j,k+1) = (Fifsf(i,j)-val)*(fabs(lsv1)/(lsv)) + Fifsf(i,j);
    Fi(i,j,k+2) = (Fifsf(i,j)-val)*(fabs(lsv2)/(lsv)) + Fifsf(i,j); 
    Fi(i,j,k+3) = (Fifsf(i,j)-val)*(fabs(lsv3)/(lsv)) + Fifsf(i,j);
    }*/
    
    
    /*FILOOP4
    {
        Fi(i,j,k+1) = Fifsf(i,j);
        Fi(i,j,k+2) = Fifsf(i,j);  
        Fi(i,j,k+3) = Fifsf(i,j);
    }*/
}

void fnpf_fg_fsf_update::fsfepol(lexer *p, fdm *a, ghostcell *pgc, slice &eta, field &Fi)
{
}

void fnpf_fg_fsf_update::velcalc(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    double H,phival;
    double epsi = 1.6*p->DXM;
    
    ULOOP
    {
        phival = 0.5*(a->phi(i,j,k) + a->phi(i+1,j,k));
        
        if(phival>epsi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(fabs(phival)<=epsi)
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
	
    
    a->u(i,j,k) =(f(i+1,j,k)-f(i,j,k))/p->DXP[IP];
    }
    
    VLOOP
    {
        phival = 0.5*(a->phi(i,j,k) + a->phi(i,j+1,k));
        
        if(phival>epsi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(fabs(phival)<=epsi)
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
        
	a->v(i,j,k) =(f(i,j+1,k)-f(i,j,k))/p->DYP[JP];
    }
    
    WLOOP
    {
        phival = 0.5*(a->phi(i,j,k) + a->phi(i,j,k+1));
        
        if(phival>epsi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(fabs(phival)<=epsi)
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
        
	a->w(i,j,k) =(f(i,j,k+1)-f(i,j,k))/p->DZP[KP];
    }
    
    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

}


void fnpf_fg_fsf_update::velcalc_sig(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
}
