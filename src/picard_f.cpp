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

#include"picard_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

picard_f::picard_f(lexer *p) : gradient(p), epsi(p->F45*p->DXM)
{
}

picard_f::~picard_f()
{
}

void picard_f::volcalc(lexer *p, fdm *a, ghostcell *pgc, field& b)
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
}

void picard_f::volcalc2(lexer *p, fdm *a, ghostcell *pgc, field& b)
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


void picard_f::correct_ls(lexer *p, fdm *a, ghostcell *pgc, field& b)
{
    double r1,r2;
    double w;
    double starttime, endtime;

    starttime=pgc->timer();

    for(int n=0;n<p->F47;++n)
    {

    volcalc2(p,a,pgc,b);

	r1=pow((0.75/PI)*vol1,1.0/3.0);
	r2=pow((0.75/PI)*vol2,1.0/3.0);

    w=r1-r2;

    //if(p->mpirank==0 &&(n==0||n==p->46-1))
    //cout<<"PICARD: "<<vol1<<" "<<vol2<<" "<<r1<<" "<<r2<<" "<<w<<endl;

    LOOP
    b(i,j,k)+=w;// *(vol1-vol2);

    }
    endtime=pgc->timer();

    //if(p->mpirank==0)
    //cout<<"PICARD TIME: "<<endtime-starttime<<endl;


}
