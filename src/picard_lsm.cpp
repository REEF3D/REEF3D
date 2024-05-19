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

#include"picard_lsm.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

picard_lsm::picard_lsm(lexer *p) : gradient(p), epsi(p->F45*p->DXM)
{
}

picard_lsm::~picard_lsm()
{
}

void picard_lsm::volcalc(lexer *p, fdm *a, ghostcell *pgc, field& b)
{

    double H = 0.0;
    vol1=0.0;

    LOOP
	{
		if(b(i,j,k)>epsi)
		H=1.0;

		if(b(i,j,k)<-epsi)
		H=0.0;

		if(fabs(b(i,j,k))<=epsi)
		H=0.5*(1.0 + b(i,j,k)/epsi + (1.0/PI)*sin((PI*b(i,j,k))/epsi));

		vol1+=pow(p->DXM, 3.0)*H;
	}

	vol1 = pgc->globalsum(vol1);


	if(p->count==1)
	{
	inivol=vol1;
	netvol=vol1;
	}
}

void picard_lsm::volcalc2(lexer *p, fdm *a, ghostcell *pgc, field& b)
{
    double H = 0.0;
    vol2=0.0;

    LOOP
	{
		if(b(i,j,k)>epsi)
		H=1.0;

		if(b(i,j,k)<-epsi)
		H=0.0;

		if(fabs(b(i,j,k))<=epsi)
		H=0.5*(1.0 + b(i,j,k)/epsi + (1.0/PI)*sin((PI*b(i,j,k))/epsi));

		vol2+=pow(p->DXM, 3.0)*H;
	}

	vol2 = pgc->globalsum(vol2);

}


void picard_lsm::correct_ls(lexer *p, fdm *a, ghostcell *pgc, field& b)
{
    double r1,r2;
    double w;

    netvol=netvol + p->Qi*p->dt - p->Qo*p->dt;

    //cout<<p->mpirank<<"  NETVOL: "<<netvol<<"  INIVOL: "<<inivol<<endl;

    for(int n=0;n<p->F47;++n)
    {

    volcalc2(p,a,pgc,b);

	r1=pow((0.75/PI)*netvol,1.0/3.0);
	r2=pow((0.75/PI)*vol2,1.0/3.0);

    w=r1-r2;

    LOOP
    b(i,j,k)+=w;

    }



}
