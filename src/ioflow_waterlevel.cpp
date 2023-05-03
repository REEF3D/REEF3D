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
#include"patchBC_interface.h"

void ioflow_f::fsfinflow(lexer *p, fdm *a, ghostcell *pgc)
{
    p->phimean=p->F56;

    if(p->F60>-1.0e20)
    {
    p->phimean=p->F60;
    p->phiout=p->F60;
    }

    if(p->F62>-1.0e20)
    p->phiout=p->F62;

    if(p->F61>-1.0e20)
    p->phiin=p->F62;

    count=0;
    zval=0.0;
    for(n=0;n<p->gcin_count;++n)
    if(p->gcin[n][3]>0)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        if(a->phi(i-1,j,k)>=0.0 && a->phi(i-1,j,k+1)<0.0)
        {
        zval+=-(a->phi(i-1,j,k)*p->DXM)/(a->phi(i-1,j,k+1)-a->phi(i-1,j,k)) + p->pos_z();
        ++count;
        }
    }

    count=pgc->globalisum(count);
    zval=pgc->globalsum(zval);

    if(count>0)
    {
    p->phimean=zval/double(count);

        if(p->F50==2 || p->F50==4)
        for(n=0;n<p->gcin_count;++n)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];

        a->phi(i-1,j,k)=p->phimean-p->pos_z();
        a->phi(i-2,j,k)=p->phimean-p->pos_z();
        a->phi(i-3,j,k)=p->phimean-p->pos_z();
        }
    p->phimean=pgc->globalmax(p->phimean);
    }

    // Outflow Water Level
    count=0;
    zval=0.0;
    for(n=0;n<p->gcout_count;++n)
    if(p->gcout[n][3]>0)
    {
    i=p->gcout[n][0];
    j=p->gcout[n][1];
    k=p->gcout[n][2];

        if(a->phi(i+1,j,k)>=0.0 && a->phi(i+1,j,k+1)<0.0)
        {
        zval+=-(a->phi(i+1,j,k)*p->DXM)/(a->phi(i+1,j,k+1)-a->phi(i+1,j,k)) + p->pos_z();
        ++count;
        }
    }

    count=pgc->globalisum(count);
    zval=pgc->globalsum(zval);

    if(count>0)
    {
    p->phiout=zval/double(count);
    p->phiout=pgc->globalmax(p->phiout);
    }
    
    
    pBC->patchBC_waterlevel(p,a,pgc,a->phi);
}

void ioflow_f::fsfrkout(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
        for(n=0;n<p->gcout_count;++n)
        {
        i=p->gcout[n][0];
        j=p->gcout[n][1];
        k=p->gcout[n][2];

        f(i+1,j,k)=a->phi(i+1,j,k);
        f(i+2,j,k)=a->phi(i+2,j,k);
        f(i+3,j,k)=a->phi(i+3,j,k);
        }
}

void ioflow_f::fsfrkin(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
        for(n=0;n<p->gcin_count;++n)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];

        f(i-1,j,k)=a->phi(i-1,j,k);
        f(i-2,j,k)=a->phi(i-2,j,k);
        f(i-3,j,k)=a->phi(i-3,j,k);
        }
}

void ioflow_f::fsfrkoutV(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
        for(int q=0;q<p->gcout_count;++q)
        {
        i=p->gcout[q][0];
        j=p->gcout[q][1];
        k=p->gcout[q][2];
        n=p->gcout[q][5];
		
		PCHECK
		{
        f.V[Ip1_J_K_4]=a->phi(i+1,j,k);
        f.V[Ip2_J_K_4]=a->phi(i+2,j,k);
        f.V[Ip3_J_K_4]=a->phi(i+3,j,k);
		}
        }
}

void ioflow_f::fsfrkinV(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
        for(int q=0;q<p->gcin_count;++q)
        {
        i=p->gcin[q][0];
        j=p->gcin[q][1];
        k=p->gcin[q][2];
        n=p->gcin[q][5];
		
		PCHECK
		{
        f.V[Im1_J_K_4]=a->phi(i+1,j,k);
        f.V[Im2_J_K_4]=a->phi(i+2,j,k);
        f.V[Im3_J_K_4]=a->phi(i+3,j,k);
		}
        }
}

void ioflow_f::fsfrkinVa(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
    for(int q=0;q<p->gcin4a_count;++q)
    {
        i=p->gcin4a[q][0];
        j=p->gcin4a[q][1];
        k=p->gcin4a[q][2];
        n=p->gcin4a[q][5];

        f.V[Im1_J_K_4a]=a->phi(i+1,j,k);
        f.V[Im2_J_K_4a]=a->phi(i+2,j,k);
        f.V[Im3_J_K_4a]=a->phi(i+3,j,k);

    }
}

void ioflow_f::fsfrkoutVa(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
    for(int q=0;q<p->gcout4a_count;++q)
    {
        i=p->gcout4a[q][0];
        j=p->gcout4a[q][1];
        k=p->gcout4a[q][2];
        n=p->gcout4a[q][5];
		
        f.V[Ip1_J_K_4a]=a->phi(i+1,j,k);
        f.V[Ip2_J_K_4a]=a->phi(i+2,j,k);
        f.V[Ip3_J_K_4a]=a->phi(i+3,j,k);
    }
}


double ioflow_f::wave_fsf(lexer *p, ghostcell *pgc, double x)
{
    double val=0.0;

    return val;
}
