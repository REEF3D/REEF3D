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
    
    if(p->mpirank==0)
    {
        cout<<"CPM graded sediment: "<<mx_nf<<" fractions, d50 "<<d50<<" dm "<<dm<<" d90 "<<d90<<" m, surface layer "<<mx_La<<" m"<<endl;
        
        for(int q=0;q<mx_nf;++q)
        cout<<"  fraction "<<q+1<<"  d = "<<mx_d[q]<<" m   share = "<<F[q]<<endl;
        
        if(fabs(d50-p->S20) > 0.05*p->S20)
        cout<<"CPM warning: S 20 "<<p->S20<<" m is not the d50 of the S 51 mixture ("<<d50<<" m)"<<endl;
    }
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
// their surface layer), and per fraction its share of the surface, its volume and its volume in the bedload layer
void CPM::mixture_log(lexer *p, ghostcell *pgc)
{
    if(mx_nf<=0)
    return;
    
    const double vpar = P.ParcelFactor*Vp;
    const size_t nj = size_t(p->knoy+2);
    
    std::vector<double> Fs(mx_nf,0.0), V(mx_nf,0.0), Vb(mx_nf,0.0);
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
        Vb[q] += vpar*P.Wt[n];
    }
    
    vs = pgc->globalsum(vs);
    
    for(int q=0;q<mx_nf;++q)
    {
        Fs[q] = pgc->globalsum(Fs[q])/MAX(vs,1.0e-300);
        V[q] = pgc->globalsum(V[q]);
        Vb[q] = pgc->globalsum(Vb[q]);
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
    
    mixout<<endl;
}
