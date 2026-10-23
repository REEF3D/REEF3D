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
#include<sstream>

CPM::CPM(lexer *p, ghostcell *pgc) : P(p,pgc), bedch(p), tauGf(p), blTx(p), blTy(p), blGx(p), blGy(p), blH(p), blC(p), blCs(p), blNc(p), blSr(p), blMc(p), zbf(p), blZb(p), zbt(p), Tau(p), Ts(p),
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

// ---------------------------------------------------------------------------------------------
// checkpoints (P 48, I 45)
// ---------------------------------------------------------------------------------------------
//  rank blob: the scalars of the model (several are counters of the rank), the random generator (Q 52
//  dispersion, Q 58 bedload layer) as its text state, one character per value, and the parcel slots.

#define CPM_INTS zbf_ini,zbl_ok,seep_ini,st_last,timestep_ini,nsub,zsplit,nrej_step,nclip_step,nit_step,restored, \
                 printcount,bag_count,blMc_ok,bl_nsleep,bl_nupd,bl_npick,bl_ndep,bl_nsus,bl_nset,bl_nslide,rs_nmerge,rs_nsplit
#define CPM_DOUBLES printtime,outvol,dtsub,cmax,Ktmax,shelter,expo_w,dbg_v0,dbg_ex,dbg_dv,dbg_sw

void CPM::ckpt_rank_state(lexer *p, std::vector<double> &v)
{
    v.push_back(1.0);   // version
    for(int x : {CPM_INTS})
    v.push_back(double(x));
    for(double x : {CPM_DOUBLES})
    v.push_back(x);

    std::ostringstream os;
    os<<rng<<' '<<gauss;
    const std::string r = os.str();
    v.push_back(double(r.size()));
    for(char c : r)
    v.push_back(double((unsigned char)c));

    P.ckpt_state(v);
}

bool CPM::ckpt_rank_restore(lexer *p, ghostcell *pgc, const std::vector<std::vector<double>> &pieces, bool same)
{
    // the scalars and the generator: this rank's piece (another decomposition: the piece of the same
    // number, else piece 0); the volume that left the domain is a sum over the ranks: rank 0 takes it all
    const size_t me = same ? 0 : (size_t(p->mpirank)<pieces.size() ? size_t(p->mpirank) : 0);
    auto head = [&](const std::vector<double> &v, size_t &k, bool set)
    {
        if(v.empty() || v[0]!=1.0)
        return false;
        k = 1;
        for(int *x : {&zbf_ini,&zbl_ok,&seep_ini,&st_last,&timestep_ini,&nsub,&zsplit,&nrej_step,&nclip_step,&nit_step,&restored,
                      &printcount,&bag_count,&blMc_ok,&bl_nsleep,&bl_nupd,&bl_npick,&bl_ndep,&bl_nsus,&bl_nset,&bl_nslide,&rs_nmerge,&rs_nsplit})
        {
            if(set) *x = int(v[k]);
            ++k;
        }
        for(double *x : {&printtime,&outvol,&dtsub,&cmax,&Ktmax,&shelter,&expo_w,&dbg_v0,&dbg_ex,&dbg_dv,&dbg_sw})
        {
            if(set) *x = v[k];
            ++k;
        }
        const size_t n = size_t(v[k++]);
        if(set)
        {
            std::string r(n,' ');
            for(size_t q=0; q<n; ++q)
            r[q] = char((unsigned char)(v[k+q]));
            std::istringstream is(r);
            is>>rng>>gauss;
        }
        k += n;
        return true;
    };

    if(pieces.empty() || me>=pieces.size())
    return false;

    size_t k=0;
    if(!head(pieces[me],k,true))
    return false;

    if(same)
    return P.ckpt_restore(p,pieces[0],k);

    // another decomposition: the parcels of every piece inside this subdomain
    if(p->mpirank==0)
    {
        double sum = 0.0;
        for(const std::vector<double> &v : pieces)
        {
            size_t j=0;
            if(head(v,j,false))
            sum += v[1+23+1];     // outvol: after the version and the 23 ints, the second double
        }
        outvol = sum;
    }
    else
    outvol = 0.0;

    const double gx = p->global_xmax, gy = p->global_ymax, gz = p->global_zmax;
    auto in = [](double v, double lo, double hi, double gmax){ return v>=lo && (v<hi || (hi>=gmax-1.0e-12 && v<=hi)); };
    auto mine = [&](double x, double y, double z)
    {
        return in(x,p->originx,p->endx,gx) && (p->j_dir==0 || in(y,p->originy,p->endy,gy)) && in(z,p->originz,p->endz,gz);
    };
    std::vector<std::vector<double>> parcels;
    for(const std::vector<double> &v : pieces)
    {
        size_t j=0;
        if(!head(v,j,false))
        return false;
        part::ckpt_parcels(v,j,mine,parcels);
    }
    P.ckpt_set(p,parcels);
    return true;
}
#undef CPM_INTS
#undef CPM_DOUBLES
