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

double interpolation::ccipol4V(double *f, slice &WL, slice &bed, double xp, double yp, double zp)
{
    double wc1,wc2;
    
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posf_sig(i,j,zp);

    i = MAX(i,0);
    i = MIN(i,p->knox-1);
    
    j = MAX(j,0);
    j = MIN(j,p->knoy-1);
    
    k = MAX(k,0);
    k = MIN(k,p->knoz);
    
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];

    if(p->j_dir==0)
    j=0;
    
    //wc
    wc1 = ((p->ZSP[IJKp1]) - zp)/(p->ZSP[IJKp1]-p->ZSP[IJK]);
    
    wc2 = ((p->ZSP[Ip1JKp1]) - zp)/(p->ZSP[Ip1JKp1]-p->ZSP[Ip1JK]);
    
    if(p->j_dir==1)
    {
    wc3 = ((p->ZSP[IJp1Kp1]) - zp)/(p->ZSP[IJp1Kp1]-p->ZSP[IJp1K]);
    
    wc4 = ((p->ZSP[Ip1Jp1Kp1]) - zp)/(p->ZSP[Ip1Jp1Kp1]-p->ZSP[Ip1Jp1K]);
    }
        
    
    i = MAX(i,0);
    i = MIN(i,p->knox-1);
    
    j = MAX(j,0);
    j = MIN(j,p->knoy-1);
    
    k = MAX(k,0);
    k = MIN(k,p->knoz-1);
    
    
    wc = MAX(wc,0);
    wc = MIN(wc,1.0);

    if(p->j_dir==0)
    value = lint4V_2D(f,i,j,k,wa,wb,wc1,wc2);
    
    if(p->j_dir==1)
    value = lint4V(f,i,j,k,wa,wb,wc1,wc2,wc3,wc4);
    
    //if(value != value)
    //cout<<i<<" 4V "<<j<<" "<<k<<"   SIG: "<<value<<" "<<wc<<" "<<(p->ZSP[IJKp1]-zp)<<" | "<<(p->ZSN[FIJKp1]-p->ZSN[FIJK])<<" | "<<(p->ZSN[FIJK]-p->ZSN[FIJKm1])<<endl;

    i=ii;
    j=jj;
    k=kk;
    
    return value;
}


double interpolation::ccipol7P(double *f, slice &WL, slice &bed, double xp, double yp, double zp)
{
    double wc1,wc2;
    double ZSN_FIp1JK,ZSN_FIJK,ZSN_FIJKp1,ZSN_FIp1JKp1;
    double ZSN_FIp1Jp1K,ZSN_FIJp1K,ZSN_FIJp1Kp1,ZSN_FIp1Jp1Kp1;
    double WLval,bedval;
    double Z1,Z2;
    
    ii=i;
    jj=j;
    kk=k;
    
    // find ijk position for interpolation cell
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posc_sig(i,j,xp,yp,zp);
    
    // above FSF
    WLval = p->ccslipol4(WL,xp,yp); 
    bedval = p->ccslipol4(bed,xp,yp); 
    
    if((zp >= WLval + bedval) || k>p->knoz)
    {
    value=0.0;
    
    return value;
    }
    
    // in cell check
    /*if((xp<p->XP[IP] ) || (xp>p->XP[IP1]))
    cout<<"CELL xp: "<<zp<<" p->XP[IP]: "<<p->XP[IP]<<" p->XP[IP1]: "<<p->XP[IP1]<<" i: "<<i<<" j: "<<j<<" k: "<<k<<" knox: "<<p->knox<<endl;
    
    if((zp<p->ZSN[FIJK] && zp<p->ZSN[FIp1JK]) || (zp>p->ZSN[FIJKp1] && zp>p->ZSN[FIp1JKp1]))
    cout<<"CELL zp: "<<zp<<" ZN[KP]: "<<p->ZSN[FIJK]<<" ZSN[KP1]: "<<p->ZSN[FIJKp1]<<" i: "<<i<<" j: "<<j<<" k: "<<k<<" knox: "<<p->knox<<endl;
    */
    
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
        
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];

    if(p->j_dir==0)
    j=0;

    //wc
    ZSN_FIJK = p->ZN[KP]*WL(i,j) + bed(i,j);
    
    ZSN_FIp1JK = p->ZN[KP]*WL(i+1,j) + bed(i+1,j);
    
    ZSN_FIJKp1 = p->ZN[KP1]*WL(i,j) + bed(i,j);
    
    ZSN_FIp1JKp1 = p->ZN[KP1]*WL(i+1,j) + bed(i+1,j);
    
    if(p->count>2)
    {
    ZSN_FIJK = p->ZSN[FIJK];
    
    ZSN_FIp1JK = p->ZSN[FIp1JK];
    
    ZSN_FIJKp1 = p->ZSN[FIJKp1];
    
    ZSN_FIp1JKp1 = p->ZSN[FIp1JKp1];
    }
    
    
    //wc
    wc1 = (ZSN_FIJKp1 - zp)/(ZSN_FIJKp1-ZSN_FIJK);
    
    wc2 = (ZSN_FIp1JKp1 - zp)/(ZSN_FIp1JKp1-ZSN_FIp1JK);
    
    //wc1 = wc2 = ((p->ZSN[FIJKp1])-zp)/(p->DZN[KP]*WL(i,j));
    
    //cout<<"(ZSN_FIJKp1-ZSN_FIJK): "<<(ZSN_FIJKp1-ZSN_FIJK)<<" p->DZN[KP]*WL(i,j): "<<p->DZN[KP]*WL(i,j)<<endl;
    
    if(p->j_dir==1)
    {
    ZSN_FIJp1K = p->ZN[KP]*WL(i,j+1) + bed(i,j+1);
    
    ZSN_FIp1Jp1K = p->ZN[KP]*WL(i+1,j+1) + bed(i+1,j+1);
    
    ZSN_FIJp1Kp1 = p->ZN[KP1]*WL(i,j+1) + bed(i,j+1);
    
    ZSN_FIp1Jp1Kp1 = p->ZN[KP1]*WL(i+1,j+1) + bed(i+1,j+1);
    
    
    wc3 = (ZSN_FIJp1Kp1 - zp)/(ZSN_FIJp1Kp1-ZSN_FIJp1K);
    
    wc4 = (ZSN_FIp1Jp1Kp1 - zp)/(ZSN_FIp1Jp1Kp1-ZSN_FIp1Jp1K);
    
    
    //wc3 = wc4 = ((p->ZSN[FIJKp1])-zp)/(p->DZN[KP]*WL(i,j));
    }
    
    // ----------------------------------------------------------
    
    
    Z1 =  ((ZSN_FIp1JK - ZSN_FIJK)/p->DXP[IP])*(xp - p->XP[IP]) 
    
        + ((ZSN_FIJp1K - ZSN_FIJK)/p->DYP[JP])*(yp - p->YP[JP]) 
    
        + ZSN_FIJK;
    
    Z2 =  ((ZSN_FIp1JKp1- ZSN_FIJKp1)/p->DXP[IP])*(xp - p->XP[IP]) 
    
        + ((ZSN_FIJp1Kp1 - ZSN_FIJKp1)/p->DYP[JP])*(yp - p->YP[JP]) 
    
        + ZSN_FIJKp1;
        
    //if(wc1<0.0 || wc1>1.0 || wc2<0.0 || wc2>1.0)    
    //cout<<"wc: "<<(Z2-zp)/(Z2-Z1)<<" wc1: "<<wc1<<" wc2: "<<wc2<<" | wa: "<<wa<<endl;
        
    wc1 = wc2 = (Z2-zp)/(Z2-Z1);
    
    wc3 = wc4 = wc1;
    
    /*
    if(wc1<0.0 || wc1>1.0 || wc2<0.0 || wc2>1.0)
    { 
    cout<<" zp: "<<zp<<" Z1 "<<Z1<<" Z2 "<<Z2<<" | wa: "<<wa<<" wc2: "<<wc2<<" | p->DXP[IP]: "<<p->DXP[IP]<<" p->XP[IP1] - p->XP[IP]: "<<p->XP[IP1] - p->XP[IP]<<endl;
    cout<<"((ZSN_FIp1JKp1- ZSN_FIJKp1)/p->DXP[IP])*(p->XP[IP1] - p->XP[IP]) : "<<((ZSN_FIp1JKp1- ZSN_FIJKp1)/p->DXP[IP])*(p->XP[IP1] - p->XP[IP])+ ZSN_FIJKp1<<endl;  
    cout<<"CELL1 zp: "<<zp<<" ZN[KP]: "<<p->ZSN[FIJK]<<" ZSN[KP1]: "<<p->ZSN[FIJKp1]<<endl;
    cout<<"CELL2 zp_iP: "<<zp<<" ZN[KP]: "<<p->ZSN[FIp1JK]<<" ZSN[KP1]: "<<p->ZSN[FIp1JKp1]<<endl<<endl;
    }*/
    
    // ----------------------------------------------------------

    if(p->j_dir==0)
    value = lint7V_2D(f,i,j,k,wa,wb,wc1,wc2);
    
    if(p->j_dir==1)
    value = lint7V(f,i,j,k,wa,wb,wc1,wc2,wc3,wc4);
    
    
    //if(value != value)
    //if(wc1<0.0 || wc1>1.0 || wc2<0.0 || wc2>1.0) 
    //cout<<i<<" "<<j<<" 7P "<<zp<<" "<<k<<" | "<<p->ZSN[FIJK]<<" "<<p->ZSN[FIJKp1]<<" . "<<ZSN_FIJp1K<<" "<<ZSN_FIJp1Kp1<<"   SIG: "<<value<<" "<<wa<<" "<<wb<<" || "<<wc1<<" "<<wc2<<" | "<<(ZSN_FIJKp1 - zp)/(ZSN_FIJKp1-ZSN_FIJK)<<" "<<wc4<<endl;
    
    
    //if(wa<0.0 || wa>1.0) 
    //cout<<i<<" "<<j<<" 7P "<<zp<<" "<<k<<" | "<<p->ZSN[FIJK]<<" "<<p->ZSN[FIJKp1]<<" . "<<ZSN_FIJp1K<<" "<<ZSN_FIJp1Kp1<<"   wav: "<<wa<<" wb: "<<wb<<endl;
    
    i=ii;
    j=jj;
    k=kk;
    
    return value;
}


double interpolation::ccipol7V(double *f, slice &WL, slice &bed, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posc_sig(i,j,zp);

    i = MAX(i,0);
    i = MIN(i,p->knox-1);
    
    j = MAX(j,0);
    j = MIN(j,p->knoy-1);
    
    k = MAX(k,0);
    k = MIN(k,p->knoz);
    
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];

    if(p->j_dir==0)
    j=0;
    
    
    //wc
    wc1 = ((p->ZSN[FIJKp1]) - zp)/(p->ZSN[FIJKp1]-p->ZSN[FIJK]);
    
    wc2 = ((p->ZSN[FIp1JKp1]) - zp)/(p->ZSN[FIp1JKp1]-p->ZSN[FIp1JK]);

    
    if(p->j_dir==1)
    {
    wc3 = ((p->ZSN[FIJp1Kp1]) - zp)/(p->ZSN[FIJp1Kp1]-p->ZSN[FIJp1K]);
    
    wc4 = ((p->ZSN[FIp1Jp1Kp1]) - zp)/(p->ZSN[FIp1Jp1Kp1]-p->ZSN[FIp1Jp1K]);
    }
    
    i = MAX(i,0);
    i = MIN(i,p->knox-1);
    
    j = MAX(j,0);
    j = MIN(j,p->knoy-1);
    
    k = MAX(k,0);
    k = MIN(k,p->knoz);
    
    
    wc = MAX(wc,0);
    wc = MIN(wc,1.0);

    if(p->j_dir==0)
    value = lint7V_2D(f,i,j,k,wa,wb,wc1,wc2);
    
    if(p->j_dir==1)
    value = lint7V(f,i,j,k,wa,wb,wc1,wc2,wc3,wc4);
    
    //if(value != value)
    //cout<<i<<" 7V "<<j<<" "<<k<<"   SIG: "<<value<<" "<<wc<<" "<<(p->ZSP[IJKp1]-zp)<<" | "<<(p->ZSN[FIJKp1]-p->ZSN[FIJK])<<" | "<<(p->ZSN[FIJK]-p->ZSN[FIJKm1])<<endl;
    
    i=ii;
    j=jj;
    k=kk;
    
    return value;
}





double interpolation::ccipol7P_old(double *f, slice &WL, slice &bed, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posf_j(yp);
    k = p->posc_sig(i,j,zp);

    i = MAX(i,0);
    i = MIN(i,p->knox-1);
    
    j = MAX(j,0);
    j = MIN(j,p->knoy-1);
    
    k = MAX(k,0);
    k = MIN(k,p->knoz);
    
    // wa
    wa = (p->XP[IP1]-xp)/p->DXP[IP];
        
    // wb
    wb = (p->YP[JP1]-yp)/p->DYP[JP];

    if(p->j_dir==0)
    j=0;
    
    //wc
    wc = ((p->ZSN[FIJKp1])-zp)/(p->DZN[KP]*WL(i,j));
    
    i = MAX(i,0);
    i = MIN(i,p->knox-1);
    
    j = MAX(j,0);
    j = MIN(j,p->knoy-1);
    
    k = MAX(k,0);
    k = MIN(k,p->knoz);
    
    
    wc = MAX(wc,0);
    wc = MIN(wc,1.0);

    if(p->j_dir==0)
    value = lint7V_old_2D(f,i,j,k,wa,wb,wc);
    
    //if(p->j_dir==1)
    //value = lint7V(f,i,j,k,wa,wb,wc);
    
    if(zp > WL(i,j) + bed(i,j))
    value=0.0;
    
    //if(value != value)
    //cout<<i<<" 7P "<<j<<" "<<k<<"   SIG: "<<value<<" "<<wc<<" "<<(p->ZSP[IJKp1]-zp)<<" | "<<(p->ZSN[FIJKp1]-p->ZSN[FIJK])<<" | "<<(p->ZSN[FIJK]-p->ZSN[FIJKm1])<<endl;
    
    i=ii;
    j=jj;
    k=kk;
    
    return value;
}
