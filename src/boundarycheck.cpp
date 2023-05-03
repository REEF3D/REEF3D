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

#include"boundarycheck.h"
#include"lexer.h"
#include"fdm.h"

boundarycheck::boundarycheck()
{
}

boundarycheck::~boundarycheck()
{
}

int boundarycheck::boundcheck(lexer *p, fdm *a, int ii, int jj, int kk, int margin)
{
    check=0;

    if(ii>=-margin && ii<p->knox+margin)
    if(jj>=-margin && jj<p->knoy+margin)
    if(kk>=-margin && kk<p->knoz+margin)
    check=1;

    if(check==1)
    if(p->flag4[(ii-p->imin-margin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0)
    if(p->flag4[(ii-p->imin+margin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin-margin)*p->kmax + kk-p->kmin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin+margin)*p->kmax + kk-p->kmin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin-margin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin+margin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0)
    check=0;

    return check;
}

int boundarycheck::boundcheck_ik(lexer *p, fdm *a, int ii, int jj, int kk, int margin)
{
    check=0;

    if(ii>=-margin && ii<p->knox+margin)
    if(kk>=-margin && kk<p->knoz+margin)
    check=1;

    if(check==1)
    if(p->flag4[(ii-p->imin-margin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0)
    if(p->flag4[(ii-p->imin+margin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin-margin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin+margin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0)
    check=0;

    return check;
}

int boundarycheck::positioncheck(lexer *p, fdm *a, double posx, double posy, double posz, int margin)
{
	int ii,jj,kk;	
    check=0;

    ii = p->posf_i(posx);
    jj = p->posf_j(posy);
    kk = p->posf_k(posz);
    
    if(ii>=-margin && ii<p->knox+margin)
    if(jj>=-margin && jj<p->knoy+margin)
    if(kk>=-margin && kk<p->knoz+margin)
    check=1;

    if(check==1)
    if(p->flag4[(ii-p->imin-margin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0)
    if(p->flag4[(ii-p->imin+margin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin-margin)*p->kmax + kk-p->kmin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin+margin)*p->kmax + kk-p->kmin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin-margin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin+margin]<0)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0)
    check=0;

    return check;
}

int boundarycheck::minboundcheck(lexer *p, int ii, int jj, int kk, int margin)
{
    check=1;

    if(ii>=0)
    if(jj>=0)
    if(kk>=0)
    check=1;

    return check;
}

int boundarycheck::maxboundcheck(lexer *p, int ii, int jj, int kk, int margin)
{
    check=1;

    if(ii<p->knox)
    if(jj<p->knoy)
    if(kk<p->knoz)
    check=1;

    return check;
}

int boundarycheck::ij_boundcheck(lexer *p, fdm *a, int ii, int jj, int margin)
{
    check=0;
    int c;

    if(ii>=0 && ii<p->knox)
    if(jj>=0 && jj<p->knoy)
    check=1;

    if(check==1)
    {
    check=0;
	
    for(c=0; c<p->knoz; ++c)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + c-p->kmin]>0 || p->flag5[IJK]<0)
    check=1;
    }

    return check;
}

int boundarycheck::ij_boundcheck_topo(lexer *p, fdm *a, int ii, int jj, int margin)
{
    check=0;
    int c;

    if(ii>=0 && ii<p->knox)
    if(jj>=0 && jj<p->knoy)
    check=1;

    if(check==1)
    {
    check=0;
	
    for(c=0; c<p->knoz; ++c)
    if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + c-p->kmin]>SOLID || p->flag5[IJK]<0)
    check=1;
    }

    return check;
}
