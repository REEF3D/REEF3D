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

#include"komega.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans_f.h"
#include"vrans_v.h"

komega::komega(lexer* p, fdm* a, ghostcell *pgc) : rans_io(p,a),bc_komega(p)
{
    if(p->B269==0)
	pvrans = new vrans_v(p,a,pgc);
	
	if(p->B269==1)
	pvrans = new vrans_f(p,a,pgc);
}

komega::~komega()
{
}

void komega::clearfield(lexer *p, fdm*  a, field& b)
{
	LOOP
	b(i,j,k)=0.0;
}

void komega::isource(lexer *p, fdm* a)
{
	ULOOP
	a->F(i,j,k)=0.0;
}

void komega::jsource(lexer *p, fdm* a)
{
	VLOOP
	a->G(i,j,k)=0.0;
}

void komega::ksource(lexer *p, fdm* a)
{
	WLOOP
	a->H(i,j,k)=0.0;
}

void  komega::eddyvisc(lexer* p,fdm* a,field& kk,field& ee, ghostcell* pgc)
{
	double H;
	double epsi = 1.6*p->dx;
	double factor;
	
    if(p->T30==0)
	LOOP
	a->eddyv(i,j,k) = kk(i,j,k)/ee(i,j,k);


    if(p->T30==1)
	LOOP
    {
	a->eddyv(i,j,k) = MAX(MAX(kk(i,j,k)
					  /((ee(i,j,k))>(1.0e-20)?(ee(i,j,k)):(1.0e20)),0.0),
					  0.0001*a->visc(i,j,k));
	}

    if(p->T30==2)
	LOOP
    {
		if(a->phi(i,j,k)>epsi)
		H=1.0;

		if(a->phi(i,j,k)<-epsi)
		H=0.0;

		if(fabs(a->phi(i,j,k))<=epsi)
		H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));
		
		factor = H*p->T31 + (1.0-H)*p->T32;
		
	a->eddyv(i,j,k) = MAX(MIN(MAX(kk(i,j,k)
					  /((ee(i,j,k))>(1.0e-20)?(ee(i,j,k)):(1.0e20)),0.0), fabs(factor*kin(i,j,k))/strainterm(p,a)),
					  0.0001*a->visc(i,j,k));
	}
	
	if(p->T10==22)
	LOOP
	a->eddyv(i,j,k) = MIN(a->eddyv(i,j,k), p->dx*p->cmu*pow((kin(i,j,k)>(1.0e-20)?(kin(i,j,k)):(1.0e20)),0.5));

    pvrans->eddyv_func(p,a);

	pgc->start4(p,a->eddyv,29);
}

void  komega::kinsource(lexer* p,fdm* a,ghostcell *pgc,field& kk,field& ee)
{
    int count=0;
    
    
    LOOP
    {
    a->L(i,j,k)=0.0;
    a->rhsvec.V[count]=0.0;
    ++count;
    }    
    
    count=0;
    
    if(p->T40==0)
    LOOP
    {
	if(wallf(i,j,k)==0)
	a->rhsvec.V[count] = pk(p,a)
                  - p->cmu*kk(i,j,k)*ee(i,j,k);
				  
    a->maxK=MAX(fabs(a->rhsvec.V[count]),a->maxK);
    ++count;
    }
    

    if(p->T40==1)
	LOOP
	{
	if(wallf(i,j,k)==0)
	a->rhsvec.V[count] = pk(p,a)
                  - p->cmu*kk(i,j,k)*MAX(ee(i,j,k),0.0);

    a->maxK=MAX(fabs(a->rhsvec.V[count]),a->maxK);
    ++count;
	}


    if(p->T40==2 || p->T40==3)
	LOOP
	{
	if(wallf(i,j,k)==0)
	a->rhsvec.V[count] = MIN(pk(p,a), p->T42*p->cmu*MAX(kk(i,j,k),0.0)*MAX(ee(i,j,k),0.0))
                  - p->cmu*kin(i,j,k)*MAX(ee(i,j,k),0.0);

    a->maxK=MAX(fabs(a->rhsvec.V[count]),a->maxK);
    ++count;
	}
    
    pvrans->kw_source(p,a,kin);
}

void  komega::epssource(lexer* p,fdm* a, ghostcell *pgc, field& kk,field& ee)
{
	double epsi = 1.6*p->dx;
	double dirac;
    int count = 0;
    
    LOOP
    {
    a->L(i,j,k)=0.0;
    a->rhsvec.V[count]=0.0;
    ++count;
    }    
    
    count=0;
	
	if(p->T41==1)
    {
        if(p->T40==0)
        LOOP
        {
        a->rhsvec.V[count] =  kw_alpha * (ee(i,j,k)/kk(i,j,k))*pk(p,a)-  kw_beta*ee(i,j,k)*ee(i,j,k);

        a->maxE=MAX(fabs(a->rhsvec.V[count]),a->maxE);
        ++count;
        }


        if(p->T40==1 || p->T40==2)
        LOOP
        {
        a->rhsvec.V[count] =  kw_alpha * (MAX(ee(i,j,k),0.0)/(kk(i,j,k)>(1.0e-10)?(fabs(kk(i,j,k))):(1.0e20))) * pk(p,a) 
					-  kw_beta * ee(i,j,k)*MAX(ee(i,j,k),0.0);

        a->maxE=MAX(fabs(a->rhsvec.V[count]),a->maxE);
        ++count;
        }

        if(p->T40==3)
        LOOP
        {
        a->rhsvec.V[count] = kw_alpha * (MAX(ee(i,j,k),0.0)/(kk(i,j,k)>(1.0e-10)?(fabs(kk(i,j,k))):(1.0e20))) 
							   * MIN(pk(p,a), p->T42*p->cmu*MAX(kin(i,j,k),0.0))
							
					- kw_beta * ee(i,j,k)*MAX(ee(i,j,k),0.0);

        a->maxE=MAX(fabs(a->rhsvec.V[count]),a->maxE);
        ++count;
        }
    }
	
	if(p->T41==2)
    {
        if(p->T40==0)
        LOOP
        {
        a->rhsvec.V[count] =  kw_alpha * pk_w(p,a) -  kw_beta * ee(i,j,k)*ee(i,j,k);

        a->maxE=MAX(fabs(a->rhsvec.V[count]),a->maxE);
        ++count;
        }


        if(p->T40==1 || p->T40==2)
        LOOP
        {
        a->rhsvec.V[count] =  kw_alpha * pk_w(p,a) -  kw_beta * ee(i,j,k)*MAX(ee(i,j,k),0.0);

        a->maxE=MAX(fabs(a->rhsvec.V[count]),a->maxE);
        ++count;
        }

        if(p->T40==3)
        LOOP
        {
        a->rhsvec.V[count] = kw_alpha * MIN(pk_w(p,a), p->T42*p->cmu*MAX(kin(i,j,k),0.0)*MAX(ee(i,j,k),0.0)/MAX(a->eddyv(i,j,k),0.0001*a->visc(i,j,k)))
                      
					  - kw_beta * ee(i,j,k)*MAX(ee(i,j,k),0.0);

        a->maxE=MAX(fabs(a->rhsvec.V[count]),a->maxE);
        ++count;
        }
	}
    
    pvrans->omega_source(p,a,kin,eps);
}

void komega::epsfsf(lexer *p, fdm* a,ghostcell *pgc,field& kk,field& ee)
{
	double epsi = p->T38*p->dx;
	double dirac;
	
	if(p->T36>0)
	LOOP
	{
		if(fabs(a->phi(i,j,k))<epsi)
		dirac = (0.5/epsi)*(1.0 + cos((p->T39*PI*a->phi(i,j,k))/epsi));
		
		if(fabs(a->phi(i,j,k))>=epsi)
		dirac=0.0;
		
	if(dirac>0.0 && p->T36==1)
	ee(i,j,k) = dirac*2.5*pow(p->cmu,-0.25)*pow(fabs(kin(i,j,k)),0.5)*(1.0/p->T37); 
	
	if(dirac>0.0 && p->T36==2)
	ee(i,j,k) = dirac*2.5*pow(p->cmu,-0.25)*pow(fabs(kin(i,j,k)),0.5)*(1.0/p->T37 + 1.0/(a->walld(i,j,k))>1.0e-20?a->walld(i,j,k):1.0e20);
	}
}

void komega::rhs(lexer *p, fdm *a, ghostcell *pgc)
{
    n=0;
	if(p->D20<=1)
	LOOP
	{
	a->L(i,j,k) += a->rhsvec.V[n];
	a->rhsvec.V[n]=0.0;
	++n;
	}
}
	

