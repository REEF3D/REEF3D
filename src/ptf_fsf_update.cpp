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
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"ptf_fsf_update.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"onephase.h"
#include"slice.h"

ptf_fsf_update::ptf_fsf_update(lexer *p, fdm *a, ghostcell *pgc) 
{
    gcval_u = 10;
    gcval_v = 11;
    gcval_w = 12;
}

ptf_fsf_update::~ptf_fsf_update()
{

}

void ptf_fsf_update::fsfupdate(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, onephase *poneph, slice &eta)
{
    // update phi
    FLUIDLOOP
    a->phi(i,j,k) = eta(i,j) + p->phimean - p->pos_z();

    pgc->start4(p,a->phi,50);
    
    // update onephase
    poneph->update(p,a,pgc,pflow);
}

void ptf_fsf_update::etaloc(lexer *p, fdm *a, ghostcell *pgc)
{
    // find k location for eta, i.e. zero level set
    SLICELOOP4
    a->etaloc(i,j) = 0;
    
    LOOP
    a->etaloc(i,j) = MAX(a->etaloc(i,j),k);
}

void ptf_fsf_update::fsfbc(lexer *p, fdm *a, ghostcell *pgc, slice &Fifsf, field &Fi)
{
    AIRLOOP
    Fi(i,j,k)=0.0; 
    
    double lsv0,lsv1,lsv2,lsv3;
    double fival,lsval,dx,dist;
    
    if(p->A323==1)
    FILOOP4
    {
    Fi(i,j,k+1) =  Fifsf(i,j);
    Fi(i,j,k+2) =  Fifsf(i,j);
    Fi(i,j,k+3) =  Fifsf(i,j);
    }
    
    
    if(p->A323==2)
    FILOOP4
    {
    lsv0 = a->phi(i,j,k);
    lsv1 = a->phi(i,j,k+1);
    lsv2 = a->phi(i,j,k+2);
    lsv3 = a->phi(i,j,k+3);
    
    lsv0 = fabs(lsv0)>1.0e-6?lsv0:1.0e20;
    
    fival = Fi(i,j,k);
 

        Fi(i,j,k+1) =  ((Fifsf(i,j)-fival)/(fabs(lsv0)))*fabs(lsv1) + Fifsf(i,j);
        Fi(i,j,k+2) =  ((Fifsf(i,j)-fival)/(fabs(lsv0)))*fabs(lsv2) + Fifsf(i,j);
        Fi(i,j,k+3) =  ((Fifsf(i,j)-fival)/(fabs(lsv0)))*fabs(lsv3) + Fifsf(i,j);
        
        //if(p->mpirank==2)
        //cout<<"F_k: "<<fival<<" Fifsf: "<<Fifsf(i,j)<<" F_k+1: "<<Fi(i,j,k+1)<<"  | lsv0: "<<lsv0<<" lsv1: "<<lsv1<<endl;
    }
    
    
    double x0,x1,x2,y0,y1,y2;
    double x,y;
    
    
    if(p->A323==3)
    FILOOP4
    {
    
    x0 = -fabs(a->phi(i,j,k-1));
    x1 = -fabs(a->phi(i,j,k));
    x2 = 0.0;
    
    y0 = Fi(i,j,k-1);
    y1 = Fi(i,j,k);
    y2 = Fifsf(i,j);
    
        x = fabs(a->phi(i,j,k+1));
        
        Fi(i,j,k+1) =   ((x-x1)/(x0-x1)) * ((x-x2)/(x0-x2)) * y0
                      + ((x-x0)/(x1-x0)) * ((x-x2)/(x1-x2)) * y1
                      + ((x-x0)/(x2-x0)) * ((x-x1)/(x2-x1)) * y2;
        
        x = fabs(a->phi(i,j,k+2));
        Fi(i,j,k+2) =   ((x-x1)/(x0-x1)) * ((x-x2)/(x0-x2)) * y0
                     + ((x-x0)/(x1-x0)) * ((x-x2)/(x1-x2)) * y1
                     + ((x-x0)/(x2-x0)) * ((x-x1)/(x2-x1)) * y2;
                       
        x = fabs(a->phi(i,j,k+3));
        Fi(i,j,k+3) =   ((x-x1)/(x0-x1)) * ((x-x2)/(x0-x2)) * y0
                     + ((x-x0)/(x1-x0)) * ((x-x2)/(x1-x2)) * y1
                     + ((x-x0)/(x2-x0)) * ((x-x1)/(x2-x1)) * y2;
        
        //cout<<"F_k: "<<Fi(i,j,k)<<" Fifsf: "<<Fifsf(i,j)<<" F_k+1: "<<Fi(i,j,k+1)<<"  | x1: "<<x1<<" x: "<<x<<endl;
    }
    
    /*
    if(i+p->origin_i>0)
    FILOOP4
    {
    
    x0 = -fabs(a->phi(i,j,k));
    x1 = 0.0;
    
    y0 = Fi(i,j,k);
    y1 = Fifsf(i,j);
    


 
        x = fabs(a->phi(i,j,k+1));
        
        Fi(i,j,k+1) =   ((x-x1)/(x0-x1))   * y0
                      + ((x-x0)/(x1-x0))  * y1;
        
        x = fabs(a->phi(i,j,k+2));
        Fi(i,j,k+2) =   ((x-x1)/(x0-x1))  * y0
                     + ((x-x0)/(x1-x0))  * y1;
                       
        x = fabs(a->phi(i,j,k+3));
        Fi(i,j,k+3) =   ((x-x1)/(x0-x1)) * y0
                     + ((x-x0)/(x1-x0)) * y1;
        
        cout<<"x0: "<<x0<<" x1: "<<x1<<" y0: "<<y0<<" y1: "<<y1<<"  | Fi(i,j,k+1): "<<Fi(i,j,k+1)<<endl;
    }
    */
    /*
    
    //fill pos[]
	for(m=0;m<=orderdir-3;m++)
	pos[m]=-dx*double(orderdir-m-2);

	pos[orderdir-2]=0.0;
	pos[orderdir-1]=dist;

	for(m=0;m<margin;m++)
	x[m]=dx*double(m+2-ys);

    //fill y[]
	if(cs==6)
	for(m=0;m<=orderdir-2;m++)
	y[m]=f(i,j,k-orderdir+m+2);

	y[orderdir-1]=wallvalue;

		

	for(q=0; q<margin; ++q)
	{
	    y[orderdir+q]=0.0;

		for(m=0;m<orderdir;m++)
		{
			weight=1.0;
			for(n=0;n<orderdir;++n)
			{
			if(m!=n)
			weight*=(x[q]-pos[n])/(pos[m]-pos[n]+1.0e-20);
			}
		y[orderdir+q]+=weight*y[m];
		}
	}*/
}

void ptf_fsf_update::fsfbc0(lexer *p, fdm *a, ghostcell *pgc, slice &Fifsf, field &Fi)
{
    AIRLOOP
    Fi(i,j,k)=0.0; 

    FILOOP4
    {
    Fi(i,j,k+1) =  Fifsf(i,j);
    Fi(i,j,k+2) =  Fifsf(i,j);
    Fi(i,j,k+3) =  Fifsf(i,j);
    }
    
}

void ptf_fsf_update::fsfbc1(lexer *p, fdm *a, ghostcell *pgc, slice &Fifsf, field &Fi)
{
    AIRLOOP
    Fi(i,j,k)=0.0; 
    
    double lsv0,lsv1,lsv2,lsv3,lsv4;
    double fival0,fival1,lsval,dx,dist;

    FILOOP4
    {
    lsv0 = p->ZP[KM1];
    lsv1 = p->ZP[KP];
    lsv2 = p->ZP[KP1];
    lsv3 = p->ZP[KP2];
    lsv4 = p->ZP[KP3];
 
    
    fival0 = Fi(i,j,k-1);
    fival1 = Fi(i,j,k);
 
    Fi(i,j,k+1) =  ((fival1-fival0)/(lsv1-lsv0))*(lsv2-lsv1) + fival1;
    Fi(i,j,k+2) =  ((fival1-fival0)/(lsv1-lsv0))*(lsv3-lsv1) + fival1;
    Fi(i,j,k+3) =  ((fival1-fival0)/(lsv1-lsv0))*(lsv4-lsv1) + fival1;
    }
    
}

void ptf_fsf_update::fsfepol(lexer *p, fdm *a, ghostcell *pgc, slice &eta, field &Fi)
{
}

void ptf_fsf_update::velcalc(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    double H,phival;
    double epsi = 0.6*p->DXM;
    
    UFLUIDLOOP
    {
        epsi = 1.6*p->DZP[KP];
        
        phival = 0.5*(a->phi(i,j,k) + a->phi(i+1,j,k));
        
        if(phival>epsi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(fabs(phival)<=epsi)
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));

    a->u(i,j,k) = H*(f(i+1,j,k)-f(i,j,k))/p->DXP[IP];
    
    
    if(i+p->origin_i==0)
    a->u(i,j,k)=0.0;
    
    }
    
    VFLUIDLOOP
    {
        epsi = 1.6*p->DZP[KP];
        
        phival = 0.5*(a->phi(i,j,k) + a->phi(i,j+1,k));
        
        if(phival>epsi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(fabs(phival)<=epsi)
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
        
	a->v(i,j,k) = H*(f(i,j+1,k)-f(i,j,k))/p->DYP[JP];
    
    if(i+p->origin_i==0)
    a->v(i,j,k)=0.0;
    }
    
    WFLUIDLOOP
    {
        epsi = 1.6*p->DZP[KP];
        
        phival = 0.5*(a->phi(i,j,k) + a->phi(i,j,k+1));
        
        if(phival>epsi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(fabs(phival)<=epsi)
		H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
 
	a->w(i,j,k) = H*(f(i,j,k+1)-f(i,j,k))/p->DZP[KP];
    
    if(i+p->origin_i==0)
    a->w(i,j,k)=0.0;
    }
    
    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

}
