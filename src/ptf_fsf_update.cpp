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
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"ptf_fsf_update.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"onephase.h"
#include"slice.h"
#include"reinifluid_RK3.h"


ptf_fsf_update::ptf_fsf_update(lexer *p, fdm *a, ghostcell *pgc)
{
    gcval_u = 10;
    gcval_v = 11;
    gcval_w = 12;
    
    //preini = new reinifluid_RK3(p,1);
}

ptf_fsf_update::~ptf_fsf_update()
{

}

void ptf_fsf_update::fsfupdate(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, onephase *poneph, slice &eta)
{
    // update phi
    FLUIDLOOP
    a->phi(i,j,k) = eta(i,j) + p->phimean - p->pos_z();
    
    //preini->start(a,p, a->phi, pgc, pflow);

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

void ptf_fsf_update::fsfbc(lexer *p, fdm *a, ghostcell *pgc, slice &Fifsf, field &Fi, slice &eta)
{
    AIRLOOP
    Fi(i,j,k)=0.0;

    double lsv0,lsv1,lsv2,lsv3;
    double fival,lsval,dx,dist;
// ------

    
    if(p->A323==1)
    FILOOP4
    {
    Fi(i,j,k+1) =  Fifsf(i,j);
    Fi(i,j,k+2) =  Fifsf(i,j);
    Fi(i,j,k+3) =  Fifsf(i,j);
    }
    
    if(p->A323==2 || p->A323==4)
    FILOOP4
    {
    lsv0 = fabs(a->phi(i,j,k));
    lsv1 = fabs(a->phi(i,j,k+1));
    lsv2 = fabs(a->phi(i,j,k+2));
    lsv3 = fabs(a->phi(i,j,k+3));

    lsv0 = (fabs(lsv0)>1.0e-6?lsv0:1.0e20);
    //lsv0 += 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;

    fival = Fi(i,j,k);


        Fi(i,j,k+1) =  ((Fifsf(i,j)-fival)/(lsv0))*lsv1 + Fifsf(i,j);
        Fi(i,j,k+2) =  ((Fifsf(i,j)-fival)/(lsv0))*lsv2 + Fifsf(i,j);
        Fi(i,j,k+3) =  ((Fifsf(i,j)-fival)/(lsv0))*lsv3 + Fifsf(i,j);

        //if(p->mpirank==2)
        //cout<<"F_k: "<<fival<<" Fifsf: "<<Fifsf(i,j)<<" F_k+1: "<<Fi(i,j,k+1)<<"  | lsv0: "<<lsv0<<" lsv1: "<<lsv1<<endl;
    }
    
    double zpos_fs,zpos_1,zpos_2,zpos_3;
    
    if(p->A323==5 || p->A323==6)
    FILOOP4
    {
        zpos_fs=eta(i,j)-(p->ZP[KP]-p->F60);
        zpos_1=p->DZP[KP];
        zpos_2=p->ZP[KP2]-p->ZP[KP];
        zpos_3=p->ZP[KP3]-p->ZP[KP];
        
        zpos_fs = (fabs(zpos_fs)>1.0e-6?zpos_fs:1.0e20);
        
        Fi(i,j,k+1)=(Fifsf(i,j)-Fi(i,j,k))/zpos_fs*zpos_1+Fi(i,j,k);
        Fi(i,j,k+2)=(Fifsf(i,j)-Fi(i,j,k))/zpos_fs*zpos_2+Fi(i,j,k);
        Fi(i,j,k+3)=(Fifsf(i,j)-Fi(i,j,k))/zpos_fs*zpos_3+Fi(i,j,k);
    }
    
    double teta,Z_t,Z_b,Fi_p,Fi_b,Fi_fsf,a_epol,b_epol,c_epol;
    
    if(p->A323==7)
    {   
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k)));
        lsv1 = fabs(a->phi(i,j,k+1));
        lsv2 = fabs(a->phi(i,j,k+2));
        lsv3 = fabs(a->phi(i,j,k+3));
        Z_t=fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k));
        Z_b=fabs(a->phi(i,j,k-1))-fabs(a->phi(i,j,k));
        
        Fi_p=Fi(i,j,k);
        Fi_b=Fi(i,j,k-1);
        Fi_fsf=Fifsf(i,j);

    }
    
    
    double x0,x1,x2,y0,y1,y2;
    double x,y;
    double denom1,denom2,denom3,denom4,denom5,denom6;
    
    if(p->A323==3)
    FILOOP4
    {

    x0 = -fabs(a->phi(i,j,k-1));
    x1 = -fabs(a->phi(i,j,k));
    x2 = 0.0;

    y0 = Fi(i,j,k-1);
    y1 = Fi(i,j,k);
    y2 = Fifsf(i,j);
    
    denom1 = fabs(x0-x1)>1.0e-6?(x0-x1):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
    denom2 = fabs(x1-x0)>1.0e-6?(x1-x0):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
    denom3 = fabs(x2-x0)>1.0e-6?(x2-x0):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
    
    denom4 = fabs(x0-x2)>1.0e-6?(x0-x2):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
    denom5 = fabs(x1-x2)>1.0e-6?(x1-x2):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
    denom6 = fabs(x2-x1)>1.0e-6?(x2-x1):1.0e20 + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.000001;
    

        x = fabs(a->phi(i,j,k+1));
        Fi(i,j,k+1) =   ((x-x1)/denom1) * ((x-x2)/denom4) * y0
                      + ((x-x0)/denom2) * ((x-x2)/denom5) * y1
                      + ((x-x0)/denom3) * ((x-x1)/denom6) * y2;

        x = fabs(a->phi(i,j,k+2));
        Fi(i,j,k+2) =   ((x-x1)/denom1) * ((x-x2)/denom4) * y0
                      + ((x-x0)/denom2) * ((x-x2)/denom5) * y1
                      + ((x-x0)/denom3) * ((x-x1)/denom6) * y2;

        x = fabs(a->phi(i,j,k+3));
        Fi(i,j,k+3) =   ((x-x1)/denom1) * ((x-x2)/denom4) * y0
                      + ((x-x0)/denom2) * ((x-x2)/denom5) * y1
                      + ((x-x0)/denom3) * ((x-x1)/denom6) * y2;

        //cout<<"F_k: "<<Fi(i,j,k)<<" Fifsf: "<<Fifsf(i,j)<<" F_k+1: "<<Fi(i,j,k+1)<<"  | x1: "<<x1<<" x: "<<x<<endl;
        
    //if(i+p->origin_i==0)
    //Fi(i-1,j,k+1) = Fi(i,j,k+1);
    }
    
    if(p->A323>=11)
    FILOOP4
    {
            double dist_s_1,dist_s_2,dist_s_3,dist_s_4,dist_s_5;
            double dist_n_1,dist_n_2,dist_n_3,dist_n_4,dist_n_5;
            double dist_b_1,dist_b_2,dist_b_3,dist_b_4,dist_b_5;
            double dist_tot_1,dist_tot_2,dist_tot_3,dist_tot_4,dist_tot_5;
            
            dist_s_1=sqrt(p->DXP[IM1]*p->DXP[IM1]+a->phi(i-1,j,k+1)*a->phi(i-1,j,k+1));
            dist_s_2=sqrt(p->DXP[IM1]*p->DXP[IM1]+a->phi(i-1,j,k+2)*a->phi(i-1,j,k+2));
            dist_s_3=sqrt(p->DXP[IM1]*p->DXP[IM1]+a->phi(i-1,j,k+3)*a->phi(i-1,j,k+3));
            dist_s_4=sqrt(p->DXP[IM1]*p->DXP[IM1]+a->phi(i-1,j,k+4)*a->phi(i-1,j,k+4));
            dist_s_5=sqrt(p->DXP[IM1]*p->DXP[IM1]+a->phi(i-1,j,k+5)*a->phi(i-1,j,k+5));
            
            dist_n_1=sqrt(p->DXP[IP]*p->DXP[IP]+a->phi(i+1,j,k+1)*a->phi(i-1,j,k+1));
            dist_n_2=sqrt(p->DXP[IP]*p->DXP[IP]+a->phi(i+1,j,k+2)*a->phi(i-1,j,k+2));
            dist_n_3=sqrt(p->DXP[IP]*p->DXP[IP]+a->phi(i+1,j,k+3)*a->phi(i-1,j,k+3));
            dist_n_4=sqrt(p->DXP[IP]*p->DXP[IP]+a->phi(i+1,j,k+4)*a->phi(i-1,j,k+4));
            dist_n_5=sqrt(p->DXP[IP]*p->DXP[IP]+a->phi(i+1,j,k+5)*a->phi(i-1,j,k+5));
            
            dist_b_1=fabs(a->phi(i,j,k+1));
            dist_b_2=fabs(a->phi(i,j,k+2));
            dist_b_3=fabs(a->phi(i,j,k+3));
            dist_b_4=fabs(a->phi(i,j,k+4));
            dist_b_5=fabs(a->phi(i,j,k+5));
            
            if(a->phi(i-1,j,k) > a->phi(i,j,k) && a->phi(i+1,j,k) > a->phi(i,j,k))
            { 
                dist_tot_1=dist_s_1+dist_n_1+dist_b_1;
                dist_tot_2=dist_s_2+dist_n_2+dist_b_2;
                dist_tot_3=dist_s_3+dist_n_3+dist_b_3;
                dist_tot_4=dist_s_4+dist_n_4+dist_b_4;
                dist_tot_5=dist_s_5+dist_n_5+dist_b_5;    
            
           Fi(i,j,k+1)=Fifsf(i-1,j)*dist_s_1/dist_tot_1+Fifsf(i,j)*dist_b_1/dist_tot_1+Fifsf(i+1,j)*dist_n_1/dist_tot_1;
        Fi(i,j,k+2)=Fifsf(i-1,j)*dist_s_2/dist_tot_2+Fifsf(i,j)*dist_b_2/dist_tot_2+Fifsf(i+1,j)*dist_n_2/dist_tot_2;
            Fi(i,j,k+3)=Fifsf(i-1,j)*dist_s_3/dist_tot_3+Fifsf(i,j)*dist_b_3/dist_tot_3+Fifsf(i+1,j)*dist_n_3/dist_tot_3;
            Fi(i,j,k+4)=Fifsf(i-1,j)*dist_s_4/dist_tot_4+Fifsf(i,j)*dist_b_4/dist_tot_4+Fifsf(i+1,j)*dist_n_4/dist_tot_4;
            Fi(i,j,k+5)=Fifsf(i-1,j)*dist_s_5/dist_tot_5+Fifsf(i,j)*dist_b_5/dist_tot_5+Fifsf(i+1,j)*dist_n_5/dist_tot_5;
            
                
            }
            
            else if(a->phi(i-1,j,k) > a->phi(i,j,k) && a->phi(i+1,j,k) <= a->phi(i,j,k))
            {
                
                dist_tot_1=dist_s_1+dist_b_1;
                dist_tot_2=dist_s_2+dist_b_2;
                dist_tot_3=dist_s_3+dist_b_3;
                dist_tot_4=dist_s_4+dist_b_4;
                dist_tot_5=dist_s_5+dist_b_5;
                
                Fi(i,j,k+1)=Fifsf(i-1,j)*dist_s_1/dist_tot_1+Fifsf(i,j)*dist_b_1/dist_tot_1;
                Fi(i,j,k+2)=Fifsf(i-1,j)*dist_s_2/dist_tot_2+Fifsf(i,j)*dist_b_2/dist_tot_2;
                Fi(i,j,k+3)=Fifsf(i-1,j)*dist_s_3/dist_tot_3+Fifsf(i,j)*dist_b_3/dist_tot_3;
                Fi(i,j,k+4)=Fifsf(i-1,j)*dist_s_4/dist_tot_4+Fifsf(i,j)*dist_b_4/dist_tot_4;
                Fi(i,j,k+5)=Fifsf(i-1,j)*dist_s_5/dist_tot_5+Fifsf(i,j)*dist_b_5/dist_tot_5;
            }
            
            else if(a->phi(i-1,j,k) <= a->phi(i,j,k) && a->phi(i+1,j,k) > a->phi(i,j,k))
            {
                dist_tot_1=dist_n_1+dist_b_1;
                dist_tot_2=dist_n_2+dist_b_2;
                dist_tot_3=dist_n_3+dist_b_3;
                dist_tot_4=dist_n_4+dist_b_4;
                dist_tot_5=dist_n_5+dist_b_5;    
                
                Fi(i,j,k+1)=Fifsf(i,j)*dist_b_1/dist_tot_1+Fifsf(i+1,j)*dist_n_1/dist_tot_1;
                Fi(i,j,k+2)=Fifsf(i,j)*dist_b_2/dist_tot_2+Fifsf(i+1,j)*dist_n_2/dist_tot_2;
                Fi(i,j,k+3)=Fifsf(i,j)*dist_b_3/dist_tot_3+Fifsf(i+1,j)*dist_n_3/dist_tot_3;
                Fi(i,j,k+4)=Fifsf(i,j)*dist_b_4/dist_tot_4+Fifsf(i+1,j)*dist_n_4/dist_tot_4;
                Fi(i,j,k+5)=Fifsf(i,j)*dist_b_5/dist_tot_5+Fifsf(i+1,j)*dist_n_5/dist_tot_5;
            }
            
            else if(a->phi(i-1,j,k) <= a->phi(i,j,k) && a->phi(i+1,j,k) <= a->phi(i,j,k))
            {
                Fi(i,j,k+1) =  Fifsf(i,j);
                Fi(i,j,k+2) =  Fifsf(i,j);
                Fi(i,j,k+3) =  Fifsf(i,j);
            }
            
            else
                cout<<"Confusion of da highest orda!"<<endl;
    }
    
    

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


    //if(i+p->origin_i==0)
    //a->u(i,j,k)=0.0;

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

    //if(i+p->origin_i==0)
    //a->v(i,j,k)=0.0;
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

    //if(i+p->origin_i==0)
    //a->w(i,j,k)=0.0;
    }

    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

}
