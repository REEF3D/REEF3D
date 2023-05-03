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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"patchBC_interface.h"

void ioflow_f::inflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    int q;
    
    if(p->B60>0 || p->count==0)
    {
        if(p->B61==1)
        {
        inflow_plain(p,a,pgc,u,v,w);

            if((p->count==0&&p->I11==1)||p->B60==3||p->B60==4)
            outflow_log(p,a,pgc,u,v,w);
        }

        if(p->B61==2 || p->B61==4 || p->B61==5)
        {
        inflow_log(p,a,pgc,u,v,w);

            if((p->count==0&&p->I11==1)||p->B60==3||p->B60==4)
            outflow_log(p,a,pgc,u,v,w);
        }

        if(p->B61==3)
        {
        inflow_water(p,a,pgc,u,v,w);

            if(p->B60==3||p->B60==4)
            outflow_water(p,a,pgc,u,v,w);
        }
        
        if(p->B75==3)
        outflow_corresponding(p,a,pgc,u,v,w);\


        for(q=0;q<4;++q)
        for(n=0;n<p->gcin_count;++n)
        {
        i=p->gcin[n][0]+q;
        j=p->gcin[n][1];
        k=p->gcin[n][2];

        if(a->phi(i,j,k)<0.0)
        a->eddyv(i,j,k)=MIN(a->eddyv(i,j,k),1.0e-4);
        }
        pgc->start4(p,a->eddyv,24);
        
    }
    
    pBC->patchBC_ioflow(p,a,pgc,u,v,w);

}

void ioflow_f::inflow_plain(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
        if(a->topo(i,j,k)>0.0)
        {
        u(i-1,j,k)=p->Ui;
        u(i-2,j,k)=p->Ui;
        u(i-3,j,k)=p->Ui;
		
		v(i-1,j,k)=0.0;
        v(i-2,j,k)=0.0;
        v(i-3,j,k)=0.0;
		
		w(i-1,j,k)=0.0;
        w(i-2,j,k)=0.0;
        w(i-3,j,k)=0.0;
        
            // Air inflow
            if(p->W50_air==1 && a->phi(i,j,k)<-0.6*p->DXM)
            {
            u(i-1,j,k)+=p->W50;
            u(i-2,j,k)+=p->W50;
            u(i-3,j,k)+=p->W50;
            }
        }
        
        if(a->topo(i,j,k)<=0.0)
        u(i-1,j,k)=u(i-2,j,k)=u(i-3,j,k)=0.0;
    }
}

void ioflow_f::inflow_log(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    double hmax=-1.0e20;
    double dmax=-1.0e20;
    double hmin=+1.0e20;

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
        
        if(a->phi(i,j,k)>0.0)
        {
        hmin=MIN(hmin,p->ZN[KP]);
        hmax=MAX(hmax,p->ZN[KP1]);
        dmax=MAX(dmax,walldin[n]);
        }
    }
    hmax=pgc->globalmax(hmax);
    hmin=pgc->globalmin(hmin);
    dmax=pgc->globalmax(dmax);


    depth=hmax-hmin;
    
    // bed shear stress and bed shear velocity
        if(p->S10==0)
        ks=p->B50;
        
        if(p->S10>0)
        ks=p->S20*p->S21;
        
        H=B=depth;
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
        
            if(a->topo(i,j,k)>0.0)
            {
            u(i-1,j,k)=u(i-2,j,k)=u(i-3,j,k)= shearvel*2.5*log(MAX(30.0*MIN(walldin[n],dmax)/ks,1.0));
            
                // Air inflow
                if(p->W50_air==1 && a->phi(i,j,k)<-0.6*p->DXM)
                {
                u(i-1,j,k)+=p->W50;
                u(i-2,j,k)+=p->W50;
                u(i-3,j,k)+=p->W50;
                }
            }
            
            if(a->topo(i,j,k)<=0.0)
            u(i-1,j,k)=u(i-2,j,k)=u(i-3,j,k)=0.0;
        }


    // calculate discharge and correct velocities
    for(int q=0; q<5; ++q)
    {
    Qin(p,a,pgc);

	if(p->B60==1)
    ratio = p->W10/p->Qi;
	
	if(p->B60==2||p->B60==4)
	ratio = hydrograph_ipol(p,pgc,hydro_in,hydro_in_count)/p->Qi;

	if(fabs(p->Qi)<1.0e-20)
	ratio=1.0;


        for(n=0;n<p->gcin_count;++n)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];
        
        if(a->topo(i,j,k)>0.0)
        {
        u(i-1,j,k)*=ratio;
        u(i-2,j,k)*=ratio;
        u(i-3,j,k)*=ratio;
        }
        }
    }

    if((p->B61==4) && p->count>0)
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
	
		if(a->phi(i-1,j,k)<-epsi1*p->DXM && a->phi(i-1,j,k)>=-epsi2*p->DXM && a->topo(i,j,k)>0.0)
        {
        fac=1.0 - fabs(a->phi(i-1,j,k))/((epsi2-epsi1)*p->DXM);
        u(i-1,j,k)=u(i-1,j,k)*fac;
        u(i-2,j,k)=u(i-2,j,k)*fac;
        u(i-3,j,k)=u(i-3,j,k)*fac;
        }

        if(a->phi(i,j,k)<-epsi2*p->DXM && a->topo(i,j,k)>0.0)
        {
        u(i-1,j,k)=0.0;
        u(i-2,j,k)=0.0;
        u(i-3,j,k)=0.0;
        }

        if(a->phi(i,j,k)<-epsi2*p->DXM)
        pgc->dirichlet_ortho(p,u,p->DXM,10,1,1);
		
    }
	
	for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
		
		v(i-1,j,k)=0.0;
        v(i-2,j,k)=0.0;
        v(i-3,j,k)=0.0;
		
		w(i-1,j,k)=0.0;
        w(i-2,j,k)=0.0;
        w(i-3,j,k)=0.0;
    }
}

void ioflow_f::inflow_water(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        if(a->phi(i-1,j,k)>=0.0 && a->topo(i,j,k)>0.0)
        {
        u(i-1,j,k)=p->Ui;
        u(i-2,j,k)=p->Ui;
        u(i-3,j,k)=p->Ui;
        }

        if(a->phi(i-1,j,k)<-epsi1*p->DXM && a->phi(i-1,j,k)>=-epsi2*p->F45*p->DXM && a->topo(i,j,k)>0.0)
        {
        fac=1.0 - fabs(a->phi(i-1,j,k))/((epsi2-epsi1)*p->F45*p->DXM);
        u(i-1,j,k)=p->Ui*fac;
        u(i-2,j,k)=p->Ui*fac;
        u(i-3,j,k)=p->Ui*fac;
        }


        if(a->phi(i-1,j,k)<-epsi2*p->F45*p->DXM && a->topo(i,j,k)>0.0)
        {
        u(i-1,j,k)=0.0;
        u(i-2,j,k)=0.0;
        u(i-3,j,k)=0.0;
        }
        
        // Air inflow
        if(p->W50_air==1 && a->phi(i,j,k)<-0.6*p->DXM)
        {
        u(i-1,j,k)+=p->W50;
        u(i-2,j,k)+=p->W50;
        u(i-3,j,k)+=p->W50;
        }
        
        if(a->topo(i,j,k)<=0.0)
        u(i-1,j,k)=u(i-2,j,k)=u(i-3,j,k)=0.0;

        if(a->phi(i-1,j,k)<-epsi2*p->DXM)
        pgc->dirichlet_ortho(p,u,p->DXM,10,1,1);
    }
	
	for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
		
		v(i-1,j,k)=0.0;
        v(i-2,j,k)=0.0;
        v(i-3,j,k)=0.0;
		
		w(i-1,j,k)=0.0;
        w(i-2,j,k)=0.0;
        w(i-3,j,k)=0.0;
    }
}

void ioflow_f::rkinflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    //inflow(p,a,pgc,u,v,w);
    
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
    u(i-1,j,k) = u(i-2,j,k) = u(i-3,j,k) = a->u(i-1,j,k);
    v(i-1,j,k) = v(i-2,j,k) = v(i-3,j,k) = a->v(i-1,j,k);
    w(i-1,j,k) = w(i-2,j,k) = w(i-3,j,k) = a->w(i-1,j,k);
    
    //u(i-1,j,k) = u(i-2,j,k) = u(i-3,j,k) = u(i,j,k);
    //v(i-1,j,k) = v(i-2,j,k) = v(i-3,j,k) = v(i,j,k);
    //w(i-1,j,k) = w(i-2,j,k) = w(i-3,j,k) = w(i,j,k);
    }
    
    pBC->patchBC_rkioflow(p,a,pgc,u,v,w);
}

void ioflow_f::flowfile(lexer *p, fdm* a, ghostcell* pgc, turbulence *pturb)
{
}

void ioflow_f::inflow_fnpf(lexer *p, fdm_fnpf*, ghostcell *pgc, double *Fi, double *Uin,slice &Fifsf, slice &eta)
{

}

void ioflow_f::rkinflow_fnpf(lexer *p, fdm_fnpf*, ghostcell *pgc, slice &frk, slice &f)
{
    
}
