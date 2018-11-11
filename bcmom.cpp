/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"bcmom.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"

bcmom::bcmom(lexer* p):surftens(p),roughness(p),kappa(0.4)
{
	dx=p->dx;

	if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;


    bckin=0;
	if(p->T10>0 || p->T10<20)
	bckin=1;
	
	if(p->B27==1)
	wallfunc_type=1;
	
	if(p->B27==2)
	wallfunc_type=2;
	
	if(p->T10==0 || p->T10>=31)
	wallfunc_type=2;
}

bcmom::~bcmom()
{
}

void bcmom::bcmom_start(fdm* a, lexer* p,ghostcell *pgc, turbulence *pturb,field& b,int gcval)
{
	int q;

	if(gcval==10&&p->B10!=0&&(p->B5==1||p->B5==2))
	{
	    QGC1LOOP
		if((p->gcb1[q][4]==5 || p->gcb1[q][4]==21 || p->gcb1[q][4]==22 || p->gcb1[q][4]==41 || p->gcb1[q][4]==42 || p->gcb1[q][4]==43) && p->gcb1[q][3]!=1 && p->gcb1[q][3]!=4)
		wall_law_u(a,p,pturb,b,p->gcb1[q][0], p->gcb1[q][1], p->gcb1[q][2], p->gcb1[q][3], p->gcb1[q][4], p->gcd1[q]);
	}

	if(gcval==11&&p->B10!=0&&(p->B5==1||p->B5==2))
	{
		QGC2LOOP
		if((p->gcb2[q][4]==5 || p->gcb2[q][4]==21 || p->gcb2[q][4]==22 || p->gcb2[q][4]==41 || p->gcb2[q][4]==42 || p->gcb2[q][4]==43) && p->gcb2[q][3]!=2 && p->gcb2[q][3]!=3)
		wall_law_v(a,p,pturb,b,p->gcb2[q][0], p->gcb2[q][1], p->gcb2[q][2], p->gcb2[q][3], p->gcb2[q][4], p->gcd2[q]);
	}

	if(gcval==12&&p->B10!=0&&(p->B5==1||p->B5==2))
	{
		QGC3LOOP
		if((p->gcb3[q][4]==5 || p->gcb3[q][4]==21 || p->gcb3[q][4]==22 || p->gcb3[q][4]==41 || p->gcb3[q][4]==42 || p->gcb3[q][4]==43) && p->gcb3[q][3]!=5 && p->gcb3[q][3]!=6)
		wall_law_w(a,p,pturb,b,p->gcb3[q][0], p->gcb3[q][1], p->gcb3[q][2], p->gcb3[q][3], p->gcb3[q][4], p->gcd3[q]);

	}
	surface_tension(a,p,a->phi,gcval);
}

void bcmom::wall_law_u(fdm* a,lexer* p, turbulence *pturb,field& b,int ii,int jj,int kk,int cs,int bc,double dist)
{
	i=ii;
	j=jj;
	k=kk;

	if(p->B23==1)
    dist=MAX(0.1*p->dx,dist);
    
    if(p->B23==2)
    dist=0.5*p->dx;
	
	
	ks=ks_val(p,a,ii,jj,kk,cs,bc);

//rough
	if(p->B5==1)
	{
		if(30.0*dist<ks)
		dist=ks/30.0;

		uplus = (1.0/kappa)*log(30.0*(dist/ks));
	}

	if(p->B5==2)
	{
    ustar=sqrt(fabs((0.5*a->visc(i,j,k)+0.5*a->visc(i,j-1,k))*a->u(i,j,k)/dist));
    value=((9.0*dist*ustar)/(0.5*a->visc(i,j,k)+0.5*a->visc(i,j-1,k)));
	
	if(value>1.0)
	uplus=(1.0/kappa)*log(value);
	
	if(value<=1.0)
	uplus=1.0e20;
	}
	
	if(wallfunc_type==1)
	a->F(i,j,k) -=  (fabs(a->u(i,j,k))*pow(p->cmu,0.25)*pow(fabs(pturb->kinval(i,j,k)),0.5))/(uplus*dx);

	if(wallfunc_type==2)
	a->F(i,j,k) -= ((fabs(a->u(i,j,k))*a->u(i,j,k))/(uplus*uplus*dx));

	a->maxF=MAX(fabs(a->F(i,j,k)),a->maxF);

}

void bcmom::wall_law_v(fdm* a,lexer* p, turbulence *pturb,field& b,int ii,int jj,int kk,int cs,int bc,double dist)
{
	i=ii;
	j=jj;
	k=kk;

	if(p->B23==1)
    dist=MAX(0.1*p->dx,dist);
    
    if(p->B23==2)
    dist=0.5*p->dx;
	
	ks=ks_val(p,a,ii,jj,kk,cs,bc);

//rough

	if(p->B5==1)
	{
		if(30.0*dist<ks)
		dist=ks/30.0;

		uplus = (1.0/kappa)*log(30.0*(dist/ks));
	}

	if(p->B5==2)
	{
    ustar=sqrt(fabs((0.5*a->visc(i,j,k)+0.5*a->visc(i,j-1,k))*(0.5*a->ro(i,j,k)+0.5*a->ro(i,j-1,k))*a->v(i,j,k)/dist));
    value=((9.0*dist*ustar)/(0.5*a->visc(i,j,k)+0.5*a->visc(i,j-1,k)));
    
	if(value>1.0)
	uplus=(1.0/kappa)*log(value);
	
	if(value<=1.0)
	uplus=1.0e20;
	}

	if(wallfunc_type==1)
	a->G(i,j,k) -=  (fabs(a->v(i,j,k))*pow(p->cmu,0.25)*pow(fabs(pturb->kinval(i,j,k)),0.5))/(uplus*dx);

	if(wallfunc_type==2)
	a->G(i,j,k) -= ((fabs(a->v(i,j,k))*a->v(i,j,k))/(uplus*uplus*dx));
	
	a->maxG=MAX(fabs(a->G(i,j,k)),a->maxG);
}

void bcmom::wall_law_w(fdm* a,lexer* p, turbulence *pturb,field& b,int ii,int jj,int kk,int cs,int bc,double dist)
{
	i=ii;
	j=jj;
	k=kk;

	if(p->B23==1)
    dist=MAX(0.1*p->dx,dist);
    
    if(p->B23==2)
    dist=0.5*p->dx;
	
	ks=ks_val(p,a,ii,jj,kk,cs,bc);

//rough

	if(p->B5==1)
	{
		if(30.0*dist<ks)
		dist=ks/30.0;

		uplus = (1.0/kappa)*log(30.0*(dist/ks));
	}

	if(p->B5==2)
	{
    ustar=sqrt(fabs((0.5*a->visc(i,j,k)+0.5*a->visc(i,j-1,k))*(0.5*a->ro(i,j,k)+0.5*a->ro(i,j-1,k))*a->w(i,j,k)/dist));
    value=((9.0*dist*ustar)/(0.5*a->visc(i,j,k)+0.5*a->visc(i,j-1,k)));
    
	if(value>1.0)
	uplus=(1.0/kappa)*log(value);
	
	if(value<=1.0)
	uplus=1.0e20;
	}

    if(wallfunc_type==1)
	a->H(i,j,k) -=  (fabs(a->w(i,j,k))*pow(p->cmu,0.25)*pow(fabs(pturb->kinval(i,j,k)),0.5))/(uplus*dx);

	if(wallfunc_type==2)
	a->H(i,j,k) -= ((fabs(a->w(i,j,k))*a->w(i,j,k))/(uplus*uplus*dx));
	
    a->maxH=MAX(fabs(a->H(i,j,k)),a->maxH);
}




