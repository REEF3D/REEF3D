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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void ioflow_f::inflow_nhflow(lexer *p, fdm_nhf *d,ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    if(p->B60>0 || p->count==0)
    {
        if(p->B61==1)
        {
        inflow_plain_nhflow(p,d,pgc,U,V,W,UH,VH,WH);

            //if((p->count==0&&p->I11==1)||p->B60==3||p->B60==4)
            //outflow_log(p,a,pgc,u,v,w);
        }

        if(p->B61==2 || p->B61==4 || p->B61==5)
        {
        inflow_log_nhflow(p,d,pgc,U,V,W,UH,VH,WH);

            //if((p->count==0&&p->I11==1)||p->B60==3||p->B60==4)
            //outflow_log(p,a,pgc,u,v,w);
        }
    }
}

void ioflow_f::inflow_plain_nhflow(lexer *p, fdm_nhf *d,ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
        
        if(p->wet[IJ]==1 && p->DF[IJK]>0)
        {
        U[Im1JK]=p->Ui;
        U[Im2JK]=p->Ui;
        U[Im3JK]=p->Ui;
        }
        
        if(p->wet[IJ]==0 || p->DF[IJK]<0)
        {
        U[Im1JK]=0.0;
        U[Im2JK]=0.0;
        U[Im3JK]=0.0;
        }
		
        V[Im1JK]=0.0;
        V[Im2JK]=0.0;
        V[Im3JK]=0.0;
		
        W[Im1JK]=0.0;
        W[Im2JK]=0.0;
        W[Im3JK]=0.0;
        
        if(p->wet[IJ]==1 && p->DF[IJK]>0)
        {
        UH[Im1JK]=(d->eta(i,j)+d->depth(i,j))*p->Ui;
        UH[Im2JK]=(d->eta(i,j)+d->depth(i,j))*p->Ui;
        UH[Im3JK]=(d->eta(i,j)+d->depth(i,j))*p->Ui;
        }
        
        if(p->wet[IJ]==0 || p->DF[IJK]<0)
        {
        UH[Im1JK]=0.0;
        UH[Im2JK]=0.0;
        UH[Im3JK]=0.0;
        }
		
        VH[Im1JK]=0.0;
        VH[Im2JK]=0.0;
        VH[Im3JK]=0.0;
		
        WH[Im1JK]=0.0;
        WH[Im2JK]=0.0;
        WH[Im3JK]=0.0;
    }
}

void ioflow_f::inflow_log_nhflow(lexer *p, fdm_nhf *d,ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    double hmax=-1.0e20;
    double dmax=-1.0e20;
    double hmin=+1.0e20;
    
    double zcoor;

    double depth, ks, H, B, M, I;
    double tau, shearvel;
    const double visc = p->W2;
    double ratio;

    // water depth
    for(n=0;n<p->gcin_count;++n)
    {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];
        
        hmin=MIN(hmin,d->WL(i,j));
        hmax=MAX(hmax,d->WL(i,j));
    }
    hmax=pgc->globalmax(hmax);
    hmin=pgc->globalmin(hmin);
    dmax=pgc->globalmax(dmax);

  
    // bed shear stress and bed shear velocity
        if(p->S10==0)
        ks=p->B50;
        
        if(p->S10>0)
        ks=p->S20*p->S21;
        
        H=B=hmax;
        M=26.0/pow(ks,(1.0/6.0));
        I=pow(p->Ui/(M*pow(H,(2.0/3.0))),2.0);
        tau=(9.81*H*I*1000.0);
		
		if(p->mpirank==0 && p->count==0)
		cout<<"I   "<<I<<endl;
		
		shearvel = p->Ui/(2.5*log((11.0*H/ks)));

        for(n=0;n<p->gcin_count;n++)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];
        
        zcoor = p->ZSP[IJK]-d->bed(i,j);
        
        //cout<<"zcoor: "<<MAX(30.0*MIN(zcoor,dmax)/ks,1.0)<<endl;
        //cout<<"Uvel: "<<shearvel*2.5*log(MAX(30.0*MIN(zcoor,hmax)/ks,1.0))<<endl;
            
            if(p->wet[IJ]==1 && p->DF[IJK]>0)
            U[Im1JK]=U[Im2JK]=U[Im3JK] = shearvel*2.5*log(MAX(30.0*MIN(zcoor,hmax)/ks,1.0));
            
            if(p->wet[IJ]==0 || p->DF[IJK]<0)
            U[Im1JK]=U[Im2JK]=U[Im3JK] = 0.0;
            
        }


    // calculate discharge and correct velocities
    for(int q=0; q<5; ++q)
    {
    Qin_nhf(p,d,pgc);

	if(p->B60==1)
    ratio = p->W10/p->Qi;
	
	if(p->B60==2||p->B60==4)
	ratio = hydrograph_ipol(p,pgc,hydro_in,hydro_in_count)/p->Qi;

	if(fabs(p->Qi)<1.0e-20)
	ratio=1.0;
    
    //if(p->mpirank==0)
    //cout<<"W10: "<<p->W10<<" Qi: "<<p->Qi<<" ratio: "<<ratio<<endl;

        for(n=0;n<p->gcin_count;++n)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];
        

        U[Im1JK]*=ratio;
        U[Im2JK]*=ratio;
        U[Im3JK]*=ratio;
        }
    }
	
	for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
        V[Im1JK]=0.0;
        V[Im2JK]=0.0;
        V[Im3JK]=0.0;
		
        W[Im1JK]=0.0;
        W[Im2JK]=0.0;
        W[Im3JK]=0.0;
        
        if(p->wet[IJ]==1 && p->DF[IJK]>0)
        {
        UH[Im1JK] = U[Im1JK]*d->WL(i,j);
        UH[Im2JK] = U[Im2JK]*d->WL(i,j);
        UH[Im3JK] = U[Im3JK]*d->WL(i,j);
        }
        
        if(p->wet[IJ]==0 || p->DF[IJK]<0)
        {
        UH[Im1JK]=0.0;
        UH[Im2JK]=0.0;
        UH[Im3JK]=0.0;
        }
        
        VH[Im1JK]=0.0;
        VH[Im2JK]=0.0;
        VH[Im3JK]=0.0;
		
        WH[Im1JK]=0.0;
        WH[Im2JK]=0.0;
        WH[Im3JK]=0.0;
    }
}

void ioflow_f::rkinflow_nhflow(lexer *p, fdm_nhf *d,ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        U[Im3JK]=U[Im2JK]=U[Im1JK]=d->U[Im1JK];
        V[Im3JK]=V[Im2JK]=V[Im1JK]=d->V[Im1JK];
        W[Im3JK]=W[Im2JK]=W[Im1JK]=d->W[Im1JK];
        
        UH[Im3JK]=UH[Im2JK]=UH[Im1JK]=d->UH[Im1JK];
        VH[Im3JK]=VH[Im2JK]=VH[Im1JK]=d->VH[Im1JK];
        WH[Im3JK]=WH[Im2JK]=WH[Im1JK]=d->WH[Im1JK];
    }

}

void ioflow_f::rkinflow_nhflow(lexer *p, fdm_nhf *d,ghostcell *pgc, double *F, double *G)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
    //cout<<G[Im1JK]<<" "<<G[IJK]<<endl;

        F[Im3JK]=F[Im2JK]=F[Im1JK]=G[Im1JK];
    }
}

void ioflow_f::wavegen_precalc_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    
}

void ioflow_f::wavegen_precalc_ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    
}