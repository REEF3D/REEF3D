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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"reefckpt.h"
#include"CPM.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<fstream>
#include<sys/stat.h>

/*--------------------------------------------------------------------
graded sediment (S 51 d fa fs)

The parcels carry their grain size (seed_diameter): each parcel holds grains of one fraction, all
parcels have the same volume, so the number of parcels of a fraction follows its volume share. S 20
is the reference grain size (parcel volume, mu(I) rheology) and should be the d50 of the mixture.

Composition of the bed surface: per bed column the topmost parcels of the bed (the ranking of
exposure_update) up to the volume of the surface layer

    V_a = theta_0 L_a A,    L_a = S 52, else 2 d90 of the seeded mixture, at least the coarsest grain

give the volume shares F_q of the fractions, and from them d50, dm and d90 of the surface (class
centres, log-linear interpolation, as the Eulerian mixture sediment_mixture::grain_stats). A column
holds a few parcels in L_a (L_a/h Q 24 on a 2D grid), so F_q of a single column is coarse; the log
gives the mean over the columns.

The composition of the surface acts on the bedload layer (Q 58, S 10 1) as in the Eulerian mixture:
  - pickup: the exposed parcels are those of the surface layer, a fraction takes its volume share F_q
    of the bed area (exposure_update), so its pickup rate is F_q V_m*(d_q) u_b(d_q)/L(d_q)
  - hiding and exposure (S 54 1, S 55): theta_c,q = Q 60 ((1-p_h)/p_h)^(-m),  p_h = sum_r F_r d_r/(d_q+d_r)
    (Wu, Wang & Jia 2000); S 54 0: theta_c,q = Q 60
  - grain roughness of the bed shear stress (S 56): ks = S 21 times S 20, or d50, dm, d90 of the surface
  - release into the resolved suspension (Q 58 2): per fraction with its settling velocity, its onset of
    suspension, its threshold and c_a,q = F_q c_a(d_q); the moving parcels of that fraction are released
  - hybrid suspension (Q 58 3 and 4): a concentration per fraction, see CPM_suspension.cpp
  - deposition: at random within the surface layer (well mixed), see deposit_depth
  - seeding: fa in the surface layer, fs below it

Output: REEF3D-CFD-CPM-Mixture.dat (surface grain sizes, surface share and volume per fraction),
the columns d50 dm d90 F_q of the bed files, the radius of the parcels.
--------------------------------------------------------------------*/

void CPM::mixture_ini(lexer *p)
{
    mx_nf = 0;
    
    double sum=0.0;
    for(int q=0;q<p->S51;++q)
    sum += MAX(p->S51_fa[q],0.0);
    
    if(p->S51<=0 || sum<=0.0)
    return;
    
    mx_nf = p->S51;
    mx_d.assign(p->S51_d, p->S51_d+mx_nf);
    
    mx_ord.resize(mx_nf);
    for(int q=0;q<mx_nf;++q)
    mx_ord[q] = q;
    
    std::sort(mx_ord.begin(), mx_ord.end(), [&](int q1, int q2){return mx_d[q1]<mx_d[q2];});
    
    std::vector<double> F(mx_nf);
    for(int q=0;q<mx_nf;++q)
    F[q] = MAX(p->S51_fa[q],0.0)/sum;
    
    double d50,dm,d90;
    mixture_sizes(F.data(),d50,dm,d90);
    
    mx_La = p->S52>0.0 ? p->S52 : 2.0*d90;
    mx_La = MAX(mx_La, mx_d[mx_ord[mx_nf-1]]);
    
    mx_F0 = F;
    mx_d50_0 = d50;
    mx_dm_0 = dm;
    mx_d90_0 = d90;
    
    // hybrid suspension (Q 58 3 and 4): concentration per fraction
    if(p->Q58>=3 && p->S10==1)
    for(int q=0;q<mx_nf;++q)
    {
        const double gmag = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
        const double R = (p->S22 - p->W1)/p->W1;
        const double d = mx_d[q];
        
        // settling velocity, Ferguson & Church (2004) as settling_velocity
        mx_ws.push_back(R*gmag*d*d/(18.0*p->W2 + sqrt(0.75*0.4*R*gmag*d*d*d)));
        
        if(q==0)
        for(int r=0;r<=mx_nf;++r)
        {
            mxFs.push_back(new slice4(p));
            
            for(int ii=-1;ii<p->knox+1;++ii)
            for(int jj=-1;jj<p->knoy+1;++jj)
            (*mxFs[r])(ii,jj) = 0.0;
        }
        
        mxC.push_back(new field4a(p));
        mxCbe.push_back(new slice4(p));
        mxSr.push_back(new slice4(p));
        mxMc.push_back(new slice4(p));
        
        const long long N = (long long)p->imax*p->jmax*p->kmax;
        
        for(long long m=0;m<N;++m)
        mxC[q]->V[m] = 0.0;
        
        for(int ii=-1;ii<p->knox+1;++ii)
        for(int jj=-1;jj<p->knoy+1;++jj)
        (*mxCbe[q])(ii,jj) = (*mxSr[q])(ii,jj) = (*mxMc[q])(ii,jj) = 0.0;
    }
    
    // accounts of the exchange with the suspension per fraction (Q 58 2 and higher)
    if(p->Q58>=2)
    for(int q=0;q<mx_nf;++q)
    {
        mxCs.push_back(new slice4(p));
        
        for(int ii=-1;ii<p->knox+1;++ii)
        for(int jj=-1;jj<p->knoy+1;++jj)
        (*mxCs[q])(ii,jj) = 0.0;
    }
    
    if(p->mpirank==0)
    {
        cout<<"CPM graded sediment: "<<mx_nf<<" fractions, d50 "<<d50<<" dm "<<dm<<" d90 "<<d90<<" m, surface layer "<<mx_La<<" m"<<endl;
        
        for(int q=0;q<mx_nf;++q)
        cout<<"  fraction "<<q+1<<"  d = "<<mx_d[q]<<" m   share = "<<F[q]<<endl;
        
        cout<<"  hiding/exposure: "<<(p->S54==1?"Wu, Wang & Jia (2000), m = ":"off")<<(p->S54==1?p->S55:0.0)
            <<"   roughness ks = S 21 x "<<(p->S56==1?"d50":(p->S56==2?"dm":(p->S56==3?"d90":"S 20")))<<(p->S56>0?" of the bed surface":"")<<endl;
        
        if(fabs(d50-p->S20) > 0.05*p->S20)
        cout<<"CPM warning: S 20 "<<p->S20<<" m is not the d50 of the S 51 mixture ("<<d50<<" m)"<<endl;
    }
}

// hiding and exposure factor of the critical Shields number of a grain d in the bed column (ic,jc):
// Wu, Wang & Jia (2000) with the composition of the surface layer of the column, as the Eulerian mixture
// (sediment_mixture::hiding); S 54 0: none
double CPM::mixture_hiding(lexer *p, int ic, int jc, double d)
{
    if(p->S54==0)
    return 1.0;
    
    const size_t c = size_t(ic+1)*size_t(p->knoy+2) + jc+1;
    const double *F = (c<mx_V.size() && mx_V[c]>0.0) ? &mx_F[c*mx_nf] : mx_F0.data();
    
    double ph=0.0;
    
    for(int r=0;r<mx_nf;++r)
    ph += F[r]*mx_d[r]/(d+mx_d[r]);
    
    ph = MIN(MAX(ph,1.0e-6),1.0-1.0e-6);
    
    return pow((1.0-ph)/ph, -p->S55);
}

// grain roughness of the bed column (ic,jc): S 21 times the grain size S 56 of its surface layer
// (0: S 20, 1: d50, 2: dm, 3: d90), as the Eulerian mixture
double CPM::mixture_ks(lexer *p, int ic, int jc)
{
    if(p->S56<=0)
    return p->S21*p->S20;
    
    ic = MAX(0, MIN(p->knox-1, ic));
    jc = MAX(0, MIN(p->knoy-1, jc));
    
    const size_t c = size_t(ic+1)*size_t(p->knoy+2) + jc+1;
    const bool ok = c<mx_V.size() && mx_V[c]>0.0;
    
    if(p->S56==1)
    return p->S21*(ok ? mx_d50[c] : mx_d50_0);
    
    if(p->S56==2)
    return p->S21*(ok ? mx_dm[c] : mx_dm_0);
    
    return p->S21*(ok ? mx_d90[c] : mx_d90_0);
}

// depth below the bed level at which a grain of the bedload layer is deposited: half a grain; graded sediment:
// at random within the surface layer L_a, which is taken as well mixed (Hirano 1971). With half a grain the
// fine grains came to rest above the coarse ones, and the ranking of the surface layer by height found a
// surface finer than the bed it was fed from (the 1 mm share of a 0.25/0.5/1 mm bed went from 0.30 to 0.21
// in 40 s of transport at theta = 0.2 in a periodic flume)
double CPM::deposit_depth(double d)
{
    if(mx_nf<=0)
    return 0.5*d;
    
    std::uniform_real_distribution<double> uni(0.0,1.0);
    
    return mx_La*uni(rng);
}

// fraction of a grain size: the nearest S 51 diameter
int CPM::mixture_fraction(double d) const
{
    int qm=0;
    
    for(int q=1;q<mx_nf;++q)
    if(fabs(d-mx_d[q]) < fabs(d-mx_d[qm]))
    qm = q;
    
    return qm;
}

// d50, dm and d90 of the volume shares F
void CPM::mixture_sizes(const double *F, double &d50, double &dm, double &d90) const
{
    dm = 0.0;
    for(int q=0;q<mx_nf;++q)
    dm += F[q]*mx_d[q];
    
    const double Pc[2] = {0.5,0.9};
    
    for(int m=0;m<2;++m)
    {
        // cumulative distribution at the class centres
        double cum=0.0, c0=0.0, val=mx_d[mx_ord[mx_nf-1]];
        bool found=false;
        
        for(int r=0;r<mx_nf && !found;++r)
        {
            const double c1 = cum + 0.5*F[mx_ord[r]];
            cum += F[mx_ord[r]];
            
            if(Pc[m]<=c1)
            {
                const double dr = mx_d[mx_ord[r]];
                
                if(r==0)
                val = dr;
                
                else
                {
                    const double dl = mx_d[mx_ord[r-1]];
                    val = exp(log(dl) + (Pc[m]-c0)/MAX(c1-c0,1.0e-20)*(log(dr)-log(dl)));
                }
                
                found = true;
            }
            
            c0 = c1;
        }
        
        if(m==0)
        d50 = val;
        
        else
        d90 = val;
    }
}

// composition of the surface layer of a bed column from its bed parcels b..e, sorted from the top
void CPM::mixture_column(lexer *p, size_t c, std::vector<std::pair<double,int>>::iterator b, std::vector<std::pair<double,int>>::iterator e, double area)
{
    const double vpar = P.ParcelFactor*Vp;
    double Va = theta_0*mx_La*area;
    double *F = &mx_F[c*mx_nf];
    double vol=0.0;
    
    for(auto q=b; q!=e && Va>0.0; ++q)
    {
        const int nn = q->second;
        const double vp = vpar*P.Wt[nn];
        const double f = MIN(1.0, Va/vp);
        
        F[mixture_fraction(P.D[nn])] += f*vp;
        vol += f*vp;
        Va -= f*vp;
    }
    
    mx_V[c] = vol;
    
    if(vol<=0.0)
    return;
    
    for(int q=0;q<mx_nf;++q)
    F[q] /= vol;
    
    mixture_sizes(F,mx_d50[c],mx_dm[c],mx_d90[c]);
}

// log of the mixture: grain sizes of the bed surface (mean over the bed columns, weighted with the volume of
// their surface layer), and per fraction its share of the surface, its volume, its volume in the bedload layer
// and its transport rate in the layer
void CPM::mixture_log(lexer *p, ghostcell *pgc)
{
    if(mx_nf<=0)
    return;
    
    const double vpar = P.ParcelFactor*Vp;
    const size_t nj = size_t(p->knoy+2);
    
    std::vector<double> Fs(mx_nf,0.0), V(mx_nf,0.0), Vb(mx_nf,0.0), Vs(mx_nf,0.0), Qb(mx_nf,0.0);
    
    // hybrid suspension: sediment volume of the fraction in its concentration
    if(!mxC.empty() && afdm!=nullptr)
    for(int q=0;q<mx_nf;++q)
    {
        for(i=0;i<p->knox;++i)
        for(j=0;j<p->knoy;++j)
        Vs[q] += susp_column(p,afdm,*mxC[q],i,j);
        
        Vs[q] = pgc->globalsum(Vs[q]);
    }

    double vs=0.0;
    
    if(mx_V.size()==size_t(p->knox+2)*nj)
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        const size_t c = size_t(i+1)*nj + j+1;
        
        vs += mx_V[c];
        
        for(int q=0;q<mx_nf;++q)
        Fs[q] += mx_V[c]*mx_F[c*mx_nf+q];
    }
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        const int q = mixture_fraction(P.D[n]);
        
        V[q] += vpar*P.Wt[n];
        
        if(P.Hop[n]>0.0)
        {
            Vb[q] += vpar*P.Wt[n];
            Qb[q] += vpar*P.Wt[n]*P.U[n];
        }
    }
    
    vs = pgc->globalsum(vs);
    
    for(int q=0;q<mx_nf;++q)
    {
        Fs[q] = pgc->globalsum(Fs[q])/MAX(vs,1.0e-300);
        V[q] = pgc->globalsum(V[q]);
        Vb[q] = pgc->globalsum(Vb[q]);
        Qb[q] = pgc->globalsum(Qb[q]);
    }
    
    if(p->mpirank!=0)
    return;
    
    double d50=0.0, dm=0.0, d90=0.0;
    
    if(vs>0.0)
    mixture_sizes(Fs.data(),d50,dm,d90);
    
    if(mixlogini==0)
    {
        mkdir("./REEF3D_CFD_CPM_Particle",0777);
        
        // a restart (I 45): continued from the checkpoint
        if(reefckpt_open(mixout,"./REEF3D_CFD_CPM_Particle/REEF3D-CFD-CPM-Mixture.dat",p->I45))
        {
            mixout<<"time\td50_surface[m]\tdm_surface[m]\td90_surface[m]";
            
            for(int q=0;q<mx_nf;++q)
            mixout<<"\tF_surface_"<<q+1;
            
            for(int q=0;q<mx_nf;++q)
            mixout<<"\tvolume_"<<q+1<<"[m3]";
            
            for(int q=0;q<mx_nf;++q)
            mixout<<"\tbedload_volume_"<<q+1<<"[m3]";
            
            for(int q=0;q<mx_nf;++q)
            mixout<<"\tqx_bedload_"<<q+1<<"[m2/s]";
            
            if(!mxC.empty())
            for(int q=0;q<mx_nf;++q)
            mixout<<"\tsuspended_volume_"<<q+1<<"[m3]";
            
            mixout<<endl;
        }
        
        mixlogini=1;
    }
    
    mixout<<p->simtime<<"\t"<<d50<<"\t"<<dm<<"\t"<<d90;
    
    for(int q=0;q<mx_nf;++q)
    mixout<<"\t"<<Fs[q];
    
    for(int q=0;q<mx_nf;++q)
    mixout<<"\t"<<V[q];
    
    for(int q=0;q<mx_nf;++q)
    mixout<<"\t"<<Vb[q];
    
    // transport rate of the fraction in the bedload layer per unit width, as qx_bedload of the CPM log
    const double Lx = p->global_xmax-p->global_xmin;
    const double Ly = p->j_dir==1 ? p->global_ymax-p->global_ymin : p->DYN[marge];
    
    for(int q=0;q<mx_nf;++q)
    mixout<<"\t"<<Qb[q]/(Lx*Ly);
    
    if(!mxC.empty())
    for(int q=0;q<mx_nf;++q)
    mixout<<"\t"<<Vs[q];
    
    mixout<<endl;
}

// hybrid suspension: the composition of the surface after the parcel step in slices of the grid, read back
// for the erosion of the next concentration step. A restart from a checkpoint finds it there: the composition
// from the parcels of the checkpoint is that after their last move, not the one the erosion used
void CPM::mixture_keep(lexer *p)
{
    if(mxFs.empty() || mx_V.empty())
    return;
    
    const size_t nj = size_t(p->knoy+2);
    
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        const size_t c = size_t(i+1)*nj + j+1;
        
        for(int q=0;q<mx_nf;++q)
        (*mxFs[q])(i,j) = mx_F[c*mx_nf+q];
        
        (*mxFs[mx_nf])(i,j) = mx_V[c];
    }
}

void CPM::mixture_restore(lexer *p)
{
    if(mxFs.empty())
    return;
    
    const size_t nj = size_t(p->knoy+2);
    const size_t nc = size_t(p->knox+2)*nj;
    
    mx_F.assign(nc*mx_nf, 0.0);
    mx_d50.assign(nc, 0.0);
    mx_dm.assign(nc, 0.0);
    mx_d90.assign(nc, 0.0);
    mx_V.assign(nc, 0.0);
    
    for(i=0;i<p->knox;++i)
    for(j=0;j<p->knoy;++j)
    {
        const size_t c = size_t(i+1)*nj + j+1;
        
        for(int q=0;q<mx_nf;++q)
        mx_F[c*mx_nf+q] = (*mxFs[q])(i,j);
        
        mx_V[c] = (*mxFs[mx_nf])(i,j);
        
        if(mx_V[c]>0.0)
        mixture_sizes(&mx_F[c*mx_nf],mx_d50[c],mx_dm[c],mx_d90[c]);
    }
}
