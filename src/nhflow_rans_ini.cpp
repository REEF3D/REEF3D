/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_rans_io.h"
#include"fdm_nhf.h"
#include"lexer.h"
#include"ghostcell.h"

void nhflow_rans_io::ini(lexer* p, fdm_nhf *d, ghostcell* pgc)
{
	/*gcval_kin=20;
	gcval_eps=30;
	gcval_edv=24;
	
	if(p->B60>=1)
	uref=p->Ui;
	
	if(p->B90>0)
	uref=0.01;
	
    if(fabs(uref)<1.0e-6)
    uref=0.01;

    plain_wallfunc(p,a,pgc);*/
}


void nhflow_rans_io::plain_wallfunc(lexer* p, fdm_nhf *d, ghostcell* pgc)
{
    /*double hmax=-1.0e20;
    double hmin=+1.0e20;

    // water depth
    LOOP
    if(a->phi(i,j,k)>0.0)
    {
        hmin=MIN(hmin,p->pos_z());
        hmax=MAX(hmax,p->pos_z());
    }

    hmax=pgc->globalmax(hmax);
    hmin=pgc->globalmin(hmin);

    depth=hmax-hmin;

	tau_calc(a,p,hmax);


	LOOP
	{
	a->eddyv(i,j,k)=sqrt(tau)*0.11*depth;

	kin(i,j,k)=0.5*kinbed;

	if(p->T10==1 || p->T10==11 || p->T10==21)
	eps(i,j,k)=(0.09*kin(i,j,k)*kin(i,j,k))/(a->eddyv(i,j,k)+1.0e-20);

	if(p->T10==2 || p->T10==12 || p->T10==22)
	eps(i,j,k)=(kin(i,j,k))/(a->eddyv(i,j,k));

	if(p->T10==3 || p->T10==13)
	eps(i,j,k)=(kin(i,j,k))/(a->eddyv(i,j,k));
	}

	GC4LOOP
	if(p->gcb4[n][4]==21 || p->gcb4[n][4]==5)
	{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];

        kin(i,j,k)=kinbed;

        if(p->T10==1 || p->T10==11)
        {
        eps(i,j,k)=(pow(0.09,0.75)*pow(kin(i,j,k),1.5))/(0.5*0.4*p->DXM);
        a->eddyv(i,j,k) = p->cmu*kin(i,j,k)*kin(i,j,k)/eps(i,j,k);
        }

        if(p->T10==2 || p->T10==12)
        {
        eps(i,j,k)=pow(kin(i,j,k),0.5)/(0.5*0.4*p->DXM*pow(0.09,0.25));
        a->eddyv(i,j,k) = kin(i,j,k)/eps(i,j,k);
        }

        if(p->T10==3 || p->T10==13)
        {
        eps(i,j,k)=pow(kin(i,j,k),0.5)/(0.5*0.4*p->DXM*pow(0.09,0.25));
        a->eddyv(i,j,k) = kin(i,j,k)/eps(i,j,k);
        }
	}

	pgc->start4(p,kin,20);
	pgc->start4(p,eps,30);
	pgc->start4(p,a->eddyv,24);*/

}

void nhflow_rans_io::inflow(lexer* p, fdm_nhf *d, ghostcell* pgc)
{
       /* GC4LOOP
        if(p->gcb4[n][4]==1)
        {
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];

		eps(i-1,j,k)=eps(i,j,k);
		eps(i-2,j,k)=eps(i,j,k);
		eps(i-2,j,k)=eps(i,j,k);

		kin(i-1,j,k)=kin(i,j,k);
		kin(i-2,j,k)=kin(i,j,k);
		kin(i-3,j,k)=kin(i,j,k);
        }
        
        GC4LOOP
        if(p->gcb4[n][4]==2)
        {
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];

		eps(i+1,j,k)=eps(i,j,k);
		eps(i+2,j,k)=eps(i,j,k);
		eps(i+2,j,k)=eps(i,j,k);

		kin(i+1,j,k)=kin(i,j,k);
		kin(i+2,j,k)=kin(i,j,k);
		kin(i+3,j,k)=kin(i,j,k);
        }*/
}

void nhflow_rans_io::tau_calc(fdm_nhf *d, lexer* p, double maxwdist)
{
	/*ks=p->B50;	
	
	H=B=depth+p->DXM;
	M=26.0/pow((ks/3.0),(1.0/6.0));
	I=pow(uref/(M*pow(H,(2.0/3.0))),2.0);
	tau=(9.81*H*I);
	kinbed = tau/sqrt(0.09);*/
}
