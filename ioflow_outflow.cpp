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

void ioflow_f::outflow_plain(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    for(n=0;n<p->gcout_count;n++)
    {
    i=p->gcout[n][0]-1;
    j=p->gcout[n][1];
    k=p->gcout[n][2];
	
		u(i,j,k)=p->Uo;
        u(i+1,j,k)=p->Uo;
        u(i+2,j,k)=p->Uo;
        u(i+3,j,k)=p->Uo;
    }
}

void ioflow_f::outflow_log(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    double hmax=-1.0e20;
    double dmax=-1.0e20;
    double hmin=+1.0e20;

    double depth, ks, H, B, M, I;
    double tau, shearvel;
    const double visc = p->W2;
    double ratio;

    // water depth
    for(n=0;n<p->gcout_count;++n)
    {
        i=p->gcout[n][0]-1;
        j=p->gcout[n][1];
        k=p->gcout[n][2];
        if(a->phi(i,j,k)>0.0)
        {
        hmin=MIN(hmin,p->pos_z());
        hmax=MAX(hmax,a->phi(i,j,k));
        dmax=MAX(dmax,walldout[n]);
        }
    }
    hmax=pgc->globalmax(hmax);
    hmin=pgc->globalmin(hmin);
    dmax=pgc->globalmax(dmax);

    depth=hmax;

    // bed shear stress and bed shear velocity
        ks=p->B50;
        H=B=depth+0.5*p->dx;
        M=26.0/pow(ks,(1.0/6.0));
        I=pow(p->Uo/(M*pow(H,(2.0/3.0))),2.0);
        tau=(9.81*H*I*1000.0);
        shearvel= sqrt(fabs(tau/1000.0));


        for(n=0;n<p->gcout_count;n++)
        {
        i=p->gcout[n][0]-1;
        j=p->gcout[n][1];
        k=p->gcout[n][2];

            u(i+1,j,k)=u(i+2,j,k)=u(i+3,j,k)= (1.0-p->B65)*p->Uo + p->B65*(shearvel*2.5*log(MAX(30.0*MIN(walldout[n],dmax)/ks,1.0)));
        }


    // calculate discharge and correct velocities
    for(int q=0; q<5; ++q)
    {
    Qout(p,a,pgc);

    if(p->B60==1)
    ratio = p->W10/p->Qo;
	
	if(p->B60==2)
	ratio = hydrograph_ipol(p,pgc,hydro_in,hydro_in_count)/p->Qi;
	
	if(p->B60==3||p->B60==4)
	ratio = hydrograph_ipol(p,pgc,hydro_out,hydro_out_count)/p->Qo;

        for(n=0;n<p->gcout_count;++n)
        {
        i=p->gcout[n][0]-1;
        j=p->gcout[n][1];
        k=p->gcout[n][2];

        u(i+1,j,k)*=ratio;
        u(i+2,j,k)*=ratio;
        u(i+3,j,k)*=ratio;
        }
    }

    if(p->B61==4 && p->count>0)
    for(n=0;n<p->gcout_count;n++)
    {
    i=p->gcout[n][0]-1;
    j=p->gcout[n][1];
    k=p->gcout[n][2];

        if(a->phi(i,j,k)<-1.0*p->F45*p->dx)
        {
        u(i+1,j,k)=0.0;
        u(i+2,j,k)=0.0;
        u(i+3,j,k)=0.0;
        }
    }
}

void ioflow_f::outflow_water(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    for(n=0;n<p->gcout_count;n++)
    {
    i=p->gcout[n][0]-1;
    j=p->gcout[n][1];
    k=p->gcout[n][2];

        if(a->phi(i,j,k)>=-p->B66_1*p->dx)
        {
        u(i+1,j,k)=p->Uo;
        u(i+2,j,k)=p->Uo;
        u(i+3,j,k)=p->Uo;
        }

        if(a->phi(i-1,j,k)<-p->B66_1*p->dx && a->phi(i-1,j,k)>=-p->B66_2*p->dx)
        {
        fac=1.0 - fabs(a->phi(i-1,j,k))/((p->B66_2-p->B66_1)*p->dx);
        u(i+1,j,k)=p->Uo*fac;
        u(i+2,j,k)=p->Uo*fac;
        u(i+3,j,k)=p->Uo*fac;
        }


        if(a->phi(i-1,j,k)<-p->B66_2*p->dx)
        {
        u(i+1,j,k)=0.0;
        u(i+2,j,k)=0.0;
        u(i+3,j,k)=0.0;
        }
    }
}



