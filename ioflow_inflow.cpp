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
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"

void ioflow_f::inflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    int q;
    
    if(p->B60<5 || p->count==0)
    {
        if(p->B61==1)
        {
        inflow_plain(p,a,pgc,u,v,w);

            if(p->count==0||p->B60==3||p->B60==4)
            outflow_plain(p,a,pgc,u,v,w);
        }

        if(p->B61==2 || p->B61==4 || p->B61==5)
        {
        inflow_log(p,a,pgc,u,v,w);

            if(p->B60==3||p->B60==4)
            outflow_log(p,a,pgc,u,v,w);
        }

        if(p->B61==3)
        {
        inflow_water(p,a,pgc,u,v,w);

            if(p->B60==3||p->B60==4)
            outflow_water(p,a,pgc,u,v,w);
        }


        if(p->B64==1)
        {
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
    }

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
        }
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
        hmin=MIN(hmin,p->pos_z());
        hmax=MAX(hmax,a->phi(i,j,k));
        dmax=MAX(dmax,walldin[n]);
        }
    }
    hmax=pgc->globalmax(hmax);
    hmin=pgc->globalmin(hmin);
    dmax=pgc->globalmax(dmax);


    depth=hmax;//-hmin;
	
	
    // bed shear stress and bed shear velocity
        ks=p->B50;
        H=B=depth+0.5*p->dx;
        M=26.0/pow(ks,(1.0/6.0));
        I=pow(p->Ui/(M*pow(H,(2.0/3.0))),2.0);
        tau=(9.81*H*I*1000.0);
        //shearvel= sqrt(fabs(tau/1000.0));
		
		if(p->mpirank==0 && p->count==0)
		cout<<"I   "<<I<<endl;
		
		shearvel = p->Ui/(2.5*log((11.0*H/ks)));

        for(n=0;n<p->gcin_count;n++)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];
        
            if(a->topo(i,j,k)>0.0)
            u(i-1,j,k)=u(i-2,j,k)=u(i-3,j,k)= shearvel*2.5*log(MAX(30.0*MIN(walldin[n],dmax)/ks,1.0));
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
	
	//if(p->mpirank==0)
	//cout<<"RATIO: "<<ratio<<" Ui: "<<p->Ui<<" Qi: "<<p->Qi<<" HGQ: "<<hydrograph_ipol(p,a,pgc)<<endl;
	

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

    if((p->B61==4 || p->B62==5) && p->count>0)
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
	
		if(a->phi(i-1,j,k)<-epsi1*p->dx && a->phi(i-1,j,k)>=-epsi2*p->dx && a->topo(i,j,k)>0.0)
        {
        fac=1.0 - fabs(a->phi(i-1,j,k))/((epsi2-epsi1)*p->dx);
        u(i-1,j,k)=u(i-1,j,k)*fac;
        u(i-2,j,k)=u(i-2,j,k)*fac;
        u(i-3,j,k)=u(i-3,j,k)*fac;
        }

        if(a->phi(i,j,k)<-epsi2*p->dx && a->topo(i,j,k)>0.0)
        {
        u(i-1,j,k)=0.0;
        u(i-2,j,k)=0.0;
        u(i-3,j,k)=0.0;
        }

        if(a->phi(i,j,k)<-epsi2*p->dx)
        pgc->dirichlet_ortho(p,u,p->dx,10,1,1);
		
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

        if(a->phi(i-1,j,k)<-epsi1*p->dx && a->phi(i-1,j,k)>=-epsi2*p->F45*p->dx && a->topo(i,j,k)>0.0)
        {
        fac=1.0 - fabs(a->phi(i-1,j,k))/((epsi2-epsi1)*p->F45*p->dx);
        u(i-1,j,k)=p->Ui*fac;
        u(i-2,j,k)=p->Ui*fac;
        u(i-3,j,k)=p->Ui*fac;
        }


        if(a->phi(i-1,j,k)<-epsi2*p->F45*p->dx && a->topo(i,j,k)>0.0)
        {
        u(i-1,j,k)=0.0;
        u(i-2,j,k)=0.0;
        u(i-3,j,k)=0.0;
        }

        if(a->phi(i-1,j,k)<-epsi2*p->dx)
        pgc->dirichlet_ortho(p,u,p->dx,10,1,1);
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
}

void ioflow_f::flowfile(lexer *p, fdm* a, ghostcell* pgc, turbulence *pturb)
{
}

void ioflow_f::inflow_fnpf(lexer *p, fdm_fnpf*, ghostcell *pgc, double *Fi, double *Uin,slice &Fifsf, slice &eta)
{

}