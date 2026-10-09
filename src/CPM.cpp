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

#include"CPM.h"
#include"lexer.h"
#include"ghostcell.h"

CPM::CPM(lexer *p, ghostcell *pgc) : P(p,pgc), bedch(p), tauGf(p), blTx(p), blTy(p), blGx(p), blGy(p), blH(p), blC(p), blCs(p), blNc(p), blSr(p), blMc(p), zbf(p), blZb(p), Tau(p), Ts(p),
                                               cellSum(p), Us(p), Vs(p), Ws(p), Pov(p), Fxy(p), Yr(p), Pse(p), Gsz(p), Pnos(p), Kc(p), KUx(p), KUy(p), KUz(p), Dsrc(p), dSx(p), dSy(p), dSz(p), Locc(p), Lout(p), Lin(p), Ltc(p), Lloc(p), LA(p), Lh0(p), Lh1(p), Lh2(p), Lh3(p),
                                               dPx(p),dPy(p),dPz(p),dTx(p),dTy(p),dTz(p),
                                               Kt(p),dKx(p),dKy(p),dKz(p),T0e(p),Tiso(p),rng(20261003ULL+p->mpirank),gauss(0.0,1.0)
{
    relax_ini(p);

    printcount=0;
    printtime=0.0;
    nsub=1;
    nrej_step=0;
    nclip_step=0;
    restored=0;
    logini=0;
    outvol=0.0;
    dtsub=0.0;
    cmax=0.0;
    
    // minimum cell size
    hmin=1.0e20;
    BASELOOP
    {
        hmin = MIN(hmin,p->DXN[IP]);
        hmin = MIN(hmin,p->DZN[KP]);
        
        if(p->j_dir==1)
        hmin = MIN(hmin,p->DYN[JP]);
    }
    hmin = pgc->globalmin(hmin);
    
    // open boundaries for the parcels (in- and outflow), from the grid boundary groups
    // sides: 1 -x, 2 +y, 3 -y, 4 +x, 5 -z, 6 +z
    for(int qn=0;qn<6;++qn)
    open_side[qn]=0;
    
    for(int qn=0;qn<p->gcb4_count;++qn)
    if(p->gcb4[qn][3]>=1 && p->gcb4[qn][3]<=6)
    if(p->gcb4[qn][4]==1 || p->gcb4[qn][4]==6 || p->gcb4[qn][4]==2 || p->gcb4[qn][4]==7 || p->gcb4[qn][4]==8)
    open_side[p->gcb4[qn][3]-1]=1;
    
    for(int qn=0;qn<6;++qn)
    open_side[qn] = pgc->globalimax(open_side[qn]);
    
    // vertical domain decomposition
    zsplit = pgc->globalimax((p->nb5>=0 || p->nb6>=0) ? 1 : 0);
    epsi=p->psi;
    
    // periodic sides are neither open nor walls
    perx = p->periodic1;
    pery = p->j_dir==1 ? p->periodic2 : 0;
    
    if(perx>0)
    open_side[0]=open_side[3]=0;
    
    if(pery>0)
    open_side[1]=open_side[2]=0;
    
    Ktmax=0.0;
    
    if(p->Q50==1 && p->Q11==2)
    coupled = this;

    // packed bed parameters
    theta_max = p->Q32>0.0 ? p->Q32 : (1.0-p->S24) + 0.035;
    theta_bed = p->Q26*(1.0-p->S24);
    theta_0 = 1.0-p->S24;
    
    Fr = p->Q33;
    eta0 = p->Q34;
    eta1 = p->Q35;
    mu_s = p->Q36;
    mu_2 = p->Q37;
    I0 = p->Q38;
    
    Vp = (1.0/6.0)*PI*pow(p->S20,3.0);

    if(p->mpirank==0 && p->Q11==2)
    {
        cout<<"CPM MP-PIC: time scheme "<<(p->Q10==2?"RK2":"EE1")<<"  stress model Q 12 "<<p->Q12<<"  friction Q 13 "<<p->Q13<<endl;
        
        if(p->Q12==2)
        cout<<"CPM packed bed: theta_bed "<<theta_bed<<" theta_0 "<<theta_0<<" theta_max "<<theta_max<<" Fr "<<Fr<<" eta0 "<<eta0<<" eta1 "<<eta1
            <<" mu_s "<<mu_s<<" mu_2 "<<mu_2<<" I0 "<<I0<<endl;
            
        if(theta_0>=theta_max)
        cout<<"CPM warning: 1-S 24 < theta_max is required, check Q 32 and S 24"<<endl;
    }
}
