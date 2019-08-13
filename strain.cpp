/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"strain.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"

strain::strain(lexer *p, fdm *a)	: gradient(p),epsi(p->F45*p->dx),Pk(p)
{
	deltax=p->dx;
}

strain::~strain()
{
}

void strain::wallf_update(lexer *p, fdm *a, ghostcell *pgc, fieldint &wallf)
{
	int n;
	LOOP
	wallf(i,j,k)=0;
	
	GC4LOOP
	if((p->gcb4[n][4]==21 || p->gcb4[n][4]==22 || p->gcb4[n][4]==5 || p->gcb4[n][4]==41  || p->gcb4[n][4]==42 || p->gcb4[n][4]==43 || p->gcb4[n][4]==6 || p->gcb4[n][4]==7 || p->gcb4[n][4]==8))
	{
	i = p->gcb4[n][0];
	j = p->gcb4[n][1];
	k = p->gcb4[n][2];
	
	wallf(i,j,k)=1;
	}
}

void strain::Pk_update(lexer *p, fdm *a, ghostcell *pgc)
{
	int n;
	
	LOOP
    {
	s11 = pudx(p,a);
	s22 = pvdy(p,a);
	s33 = pwdz(p,a);
	s12 = (pudy(p,a) + pvdx(p,a));
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = (pvdz(p,a) + pwdy(p,a));

    Pk(i,j,k) = a->eddyv(i,j,k)*(2.0*s11*s11 + 2.0*s22*s22 + 2.0*s33*s33 + s12*s12 + s13*s13 + s23*s23);
    }	
}

double strain::sij(lexer *p, fdm *a, int ii, int jj)
{
	double s=0.0;

	if(ii==1 && jj==1)
	s = 2.0*pudx(p,a);

	if((ii==1 && jj==2) || (ii==2 && jj==1))
	s = pudy(p,a) + pvdx(p,a);

	if((ii==1 && jj==3) ||( ii==3 && jj==1))
	s = pudz(p,a) + pwdx(p,a);

	if(ii==2 && jj==2)
	s = pvdy(p,a);

	if((ii==2 && jj==3) || (ii==3 && jj==2))
	s = pvdz(p,a) + pwdy(p,a);

	if(ii==3 && jj==3)
	s = 2.0*pwdz(p,a);

	return 0.5*s;
}


double strain::qij(lexer *p, fdm *a, int ii, int jj)
{
	double q=0.0;

	if((ii==1 && jj==1) || (ii==2 && jj==2) || (ii==3 && jj==3))
	q = 0.0;

	if(ii==1 && jj==2)
	q = pudy(p,a) - pvdx(p,a);

	if(ii==2 && jj==1)
	q = -pudy(p,a) + pvdx(p,a);

	if(ii==1 && jj==3)
	q = pudz(p,a) - pwdz(p,a);

	if(ii==3 && jj==1)
	q = -pudx(p,a) + pwdz(p,a);

	if(ii==2 && jj==3)
	q = pvdz(p,a) - pwdy(p,a);

	if(ii==3 && jj==2)
	q = -pvdz(p,a) + pwdy(p,a);

	return q;
}

double strain::pk(lexer *p, fdm *a)
{ 
	return Pk(i,j,k);
}

double strain::pk_k(lexer *p, fdm *a)
{
	double pkterm=0.0;

    if(p->j_dir==1)
    {
	s11 = pudx(p,a);
	s22 = pvdy(p,a);
	s33 = pwdz(p,a);
	s12 = (pudy(p,a) + pvdx(p,a));
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = (pvdz(p,a) + pwdy(p,a));
    }
    
    if(p->j_dir==0)
    {
	s11 = pudx(p,a);
	s22 = 0.0;
	s33 = pwdz(p,a);
	s12 = 0.0;
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = 0.0;
    }

    pkterm = (2.0*s11*s11 + 2.0*s22*s22 + 2.0*s33*s33 + s12*s12 + s13*s13 + s23*s23);	

	return pkterm;
}

double strain::pk_w(lexer *p, fdm *a)
{
	double pkterm=0.0;

	if(p->j_dir==1)
    {
	s11 = pudx(p,a);
	s22 = pvdy(p,a);
	s33 = pwdz(p,a);
	s12 = (pudy(p,a) + pvdx(p,a));
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = (pvdz(p,a) + pwdy(p,a));
    }
    
    if(p->j_dir==0)
    {
	s11 = pudx(p,a);
	s22 = 0.0;
	s33 = pwdz(p,a);
	s12 = 0.0;
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = 0.0;
    }

    pkterm = (2.0*s11*s11 + 2.0*s22*s22 + 2.0*s33*s33 + s12*s12 + s13*s13 + s23*s23);

	return pkterm;
}

double strain::strainterm(lexer *p, fdm *a)
{
	double s=0.0;

	s11 = pudx(p,a);
	s22 = pvdy(p,a);
	s33 = pwdz(p,a);
	s12 = (pudy(p,a) + pvdx(p,a));
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = (pvdz(p,a) + pwdy(p,a));

    s = sqrt(s11*s11 + s22*s22 + s33*s33 + 0.5*s12*s12 + 0.5*s13*s13 + 0.5*s23*s23);

	return s;
}

double strain::strainplain(lexer *p, fdm *a)
{
	double s=0.0;

	s11 = pudx(p,a);
	s22 = pvdy(p,a);
	s33 = pwdz(p,a);
	s12 = (pudy(p,a) + pvdx(p,a));
	s13 = (pudz(p,a) + pwdx(p,a));
	s23 = (pvdz(p,a) + pwdy(p,a));

    s=fabs(s11)+fabs(s22)+fabs(s33)+0.5*fabs(s12)+0.5*fabs(s13)+0.5*fabs(s13);

	return s;
}
