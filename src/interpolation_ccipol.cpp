/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"interpolation.h"
#include"field.h"
#include"slice.h"
#include"lexer.h"

double interpolation::ccipol1(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posc_i(xp)-1;
    j = p->posf_j(yp);
    k = p->posf_k(zp);

    // wa
    wa = (p->XN[IP2]-xp)/p->DXN[IP2];
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];
    
    if(p->j_dir==0)
    value = lint1_2D(f,i,j,k,wa,wb,wc);
    
    if(p->j_dir==1)
    value = lint1(f,i,j,k,wa,wb,wc);
    
    //cout<<" | i: "<<i<<" j: "<<j<<" k: "<<k<<" wa: "<<wa<<" wb: "<<wb<<" wc: "<<wc<<" xp: "<<xp<<" f: "<<f(i,j,k)<<" | ";

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol2(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posc_j(yp)-1;
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
    
    // wb
    wb = (p->YN[JP2]-yp)/p->DYN[JP2];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];

    if(p->j_dir==0)
    value = lint2_2D(f,i,j,k,wa,wb,wc);
    
    if(p->j_dir==1)
    value = lint2(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol3(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posc_k(zp)-1;
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZN[KP2]-zp)/p->DZN[KP2];
    

    if(p->j_dir==0)
    value = lint3_2D(f,i,j,k,wa,wb,wc);
    
    if(p->j_dir==1)
    value = lint3(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol4(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];

    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];

    // lint
    if(p->j_dir==0)
    value = lint4_2D(f,i,j,k,wa,wb,wc);
    
    if(p->j_dir==1)
    value = lint4(f,i,j,k,wa,wb,wc);
    
    //cout<<" | i: "<<i<<" j: "<<j<<" k: "<<k<<" wa: "<<wa<<" wb: "<<wb<<" wc: "<<wc<<" xp: "<<xp<<" yp: "<<yp<<" zp: "<<zp<<" | ";

    i=ii;
    j=jj;
    k=kk;
    
    return value;
}

double interpolation::ccipol4_c(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];

    // lint
    if(p->j_dir==0)
    value = lint4_2D(f,i,j,k,wa,wb,wc);
    
    if(p->j_dir==1)
    value = lint4c(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;
    
    return value;
}

double interpolation::ccipol4phi(fdm *a,field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];
    
    if(p->j_dir==0)
    value =  lint4phi_2D(a,f,i,j,k,wa,wb,wc);
    
    if(p->j_dir==1)
    value =  lint4phi(a,f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol4press(fdm *a,field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posf_k(zp);
    
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZN[KP];

    value =  lint4phi(a,f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol1_a(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posc_i(xp)-1;
    j = p->posf_j(yp);
    k = p->posf_k(zp);

    // wa
    wa = (p->XN[IP2]-xp)/p->DXN[IP2];
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];

    value = lint_a(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol2_a(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posc_j(yp)-1;
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
    
    // wb
    wb = (p->YN[JP2]-yp)/p->DYN[JP2];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];

    value = lint_a(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol3_a(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posc_k(zp)-1;
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZN[KP2]-zp)/p->DZN[KP2];

    value = lint_a(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol4a(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];

    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];

    
    if(p->j_dir==0)
    value = lint_a_2D(f,i,j,k,wa,wb,wc);
    
    if(p->j_dir==1)
    value = lint_a(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol4_b(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];

    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];


    value =  lint4b(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol4_kin(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];

    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];

    
    value = lint4kin(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol1c(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posc_i(xp)-1;
    j = p->posf_j(yp);
    k = p->posf_k(zp);

    // wa
    wa = (p->XN[IP2]-xp)/p->DXN[IP1];
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];
    
    if(p->j_dir==0)
    value = lint1_2D(f,i,j,k,wa,wb,wc);
    
    if(p->j_dir==1)
    value = lint1c(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol2c(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posc_j(yp)-1;
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
    
    // wb
    wb = (p->YN[JP2]-yp)/p->DYN[JP1];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];

    if(p->j_dir==0)
    value = lint2_2D(f,i,j,k,wa,wb,wc);
    
    if(p->j_dir==1)
    value = lint2c(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol3c(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posc_k(zp)-1;
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZN[KP2]-zp)/p->DZN[KP1];

    if(p->j_dir==0)
    value = lint3_2D(f,i,j,k,wa,wb,wc);
    
    if(p->j_dir==1)
    value = lint3c(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol4c(double *f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];

    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZP[KP];

    if(p->j_dir==0)
    value = lint4c(f,i,j,k,wa,wb,wc);
    
    if(p->j_dir==1)
    value = lint4c(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;
    

    return value;
}