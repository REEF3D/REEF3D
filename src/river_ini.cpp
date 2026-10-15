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

#include"river_ini.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include<mpi.h>
#include<cmath>
#include<vector>
#include<algorithm>
#include<iostream>

river_ini::river_ini()
{
}


// the cells of the column below h connected to its lowest cell (neighbours in j): side pools behind
// higher ground stay dry
void river_ini::xs_run(const river_ini::xsec &c, double h, int &n0, int &n1)
{
    n0 = n1 = c.nmin;
    if(c.zb.empty() || h<=c.zb[c.nmin])
    {
        n1 = n0-1;
        return;
    }
    while(n0>0 && c.J[n0-1]==c.J[n0]-1 && h>c.zb[n0-1])
    --n0;
    while(n1<int(c.zb.size())-1 && c.J[n1+1]==c.J[n1]+1 && h>c.zb[n1+1])
    ++n1;
}

double river_ini::xs_area(const river_ini::xsec &c, double h)
{
    int n0,n1;
    xs_run(c,h,n0,n1);
    double A=0.0;
    for(int n=n0; n<=n1; ++n)
    A += (h-c.zb[n])*c.w[n];
    return A;
}

// unit-width conveyance q/sqrt(S) of a strip of depth d: Manning (law 0, SFLOW) kst d^(5/3); log law at the
// centre of the bottom layer (law 1, NHFLOW A 519 1, depth-uniform velocity) sqrt(g d) d u+(z0), z0 = dz0 d/2
double river_ini::strip(double d) const
{
    if(law==0)
    return kst*pow(d,5.0/3.0);

    const double ks = MAX(ksl,1.0e-4);
    const double z0 = MAX(0.5*dz0*d,ks/30.0);
    const double uplus = 2.5*MAX(1.0,log(30.0*z0/ks));
    return sqrt(g*d)*d*uplus;
}

// conveyance of the column (strips, law), and the conveyance-weighted length along of the column
double river_ini::xs_conv(const river_ini::xsec &c, double h, double kst, double &dsw) const
{
    int n0,n1;
    xs_run(c,h,n0,n1);
    double K=0.0, kd=0.0;
    dsw = 0.0;
    for(int n=n0; n<=n1; ++n)
    {
        const double k = (law==0) ? kst*pow(h-c.zb[n],5.0/3.0)*c.w[n] : strip(h-c.zb[n])*c.w[n];
        K += k;
        kd += k*c.ds[n];
    }
    if(K>0.0)
    dsw = kd/K;
    return K;
}

double river_ini::xs_zmin(const river_ini::xsec &c)
{
    double z=1.0e20;
    for(double v : c.zb)
    z = std::min(z,v);
    return z;
}

double river_ini::xs_dsmean(const river_ini::xsec &c)
{
    double s=0.0;
    for(double v : c.ds)
    s += v;
    return c.ds.empty() ? 0.0 : s/double(c.ds.size());
}

double river_ini::xs_energy(const river_ini::xsec &c, double h, double Q, double g)
{
    const double A = xs_area(c,h);
    return (A>1.0e-12) ? h + Q*Q/(2.0*g*A*A) : 1.0e20;
}

// critical level: minimum of the specific energy (golden section)
double river_ini::xs_critical(const river_ini::xsec &c, double Q, double g, double hmax)
{
    double a = xs_zmin(c)+1.0e-4, b = a+hmax;
    const double r = 0.5*(sqrt(5.0)-1.0);
    double x1 = b-r*(b-a), x2 = a+r*(b-a);
    double f1 = xs_energy(c,x1,Q,g), f2 = xs_energy(c,x2,Q,g);
    for(int it=0; it<100 && b-a>1.0e-5; ++it)
    {
        if(f1<f2)
        {
            b=x2; x2=x1; f2=f1;
            x1 = b-r*(b-a);
            f1 = xs_energy(c,x1,Q,g);
        }
        else
        {
            a=x1; x1=x2; f1=f2;
            x2 = a+r*(b-a);
            f2 = xs_energy(c,x2,Q,g);
        }
    }
    return 0.5*(a+b);
}

int river_ini::start(lexer *p, ghostcell *pgc, slice &bed, int verbose, int flaw)
{
    law = flaw;
    if(p->F65<1 || p->F65>2)
    return 0;

    if(p->geo_mb>0)
    {
        if(p->mpirank==0 && verbose==1)
        cout<<"!!! F 65: not available on multi-block grids, no river free-surface initialisation !!!"<<endl;
        return 0;
    }

    if(p->F62<-1.0e19 || (p->F65==1 && p->F61<-1.0e19) || (p->F65==2 && fabs(p->W10)<1.0e-12))
    {
        if(p->mpirank==0)
        cout<<"!!! F 65 "<<p->F65<<" needs "<<(p->F65==1?"F 61 and F 62":"F 62 and W 10")<<" !!!"<<endl;
        exit(1);
    }

    NI=p->gknox;
    NJ=p->gknoy;
    const long NC = long(NI)*long(NJ);
    g = fabs(p->W22)>1.0e-12 ? fabs(p->W22) : 9.81;
    Q = fabs(p->W10);
    const double ks = MAX(p->B50,1.0e-6);
    ksl = ks;

    // law 1: sigma height of the bottom layer
    if(law==1)
    {
    const int m0 = p->margin;
    const double zr = p->ZN[p->knoz+m0]-p->ZN[m0];
    dz0 = (zr>0.0) ? p->DZN[m0]/zr : 1.0/double(MAX(p->knoz,1));
    }
    kst = 20.0/pow(ks,1.0/6.0);
    const double hmax = 50.0;

    // global cross-sections: bed, width across and length along of every cell (0 width: no cell)
    std::vector<double> gz(NC,0.0), gw(NC,0.0), gs(NC,0.0);

    SLICELOOP4
    {
        const long q = long(i+p->origin_i)*NJ + long(j+p->origin_j);
        const double w = 0.5*(p->LEN1[Im1J] + p->LEN1[IJ]);
        gz[q] = bed(i,j);
        gw[q] = w;
        gs[q] = (w>1.0e-12) ? p->AREA2D[IJ]/w : 0.0;
    }

    MPI_Allreduce(MPI_IN_PLACE,gz.data(),int(NC),MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE,gw.data(),int(NC),MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE,gs.data(),int(NC),MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    xs.assign(NI,xsec());
    for(int I=0; I<NI; ++I)
    for(int J=0; J<NJ; ++J)
    {
        const long q = long(I)*NJ + J;
        if(gw[q]>1.0e-12)
        {
            xs[I].zb.push_back(gz[q]);
            xs[I].w.push_back(gw[q]);
            xs[I].ds.push_back(gs[q]);
            xs[I].J.push_back(J);
            if(xs[I].zb.back()<xs[I].zb[xs[I].nmin])
            xs[I].nmin = int(xs[I].zb.size())-1;
        }
    }

    // water level of the columns
    wl.assign(NI,p->F62);
    int ncrit=0;
    double hn_out=0.0;

    if(p->F65==1)
    {
        std::vector<double> s(NI,0.0);
        for(int I=1; I<NI; ++I)
        s[I] = s[I-1] + 0.5*(xs_dsmean(xs[I-1])+xs_dsmean(xs[I]));

        const double L = MAX(s[NI-1],1.0e-12);
        for(int I=0; I<NI; ++I)
        wl[I] = p->F61 + (p->F62-p->F61)*s[I]/L;
    }

    if(p->F65==2)
    {
        int last=NI-1;
        while(last>0 && xs[last].zb.empty())
        --last;

        const double hc_out = xs_critical(xs[last],Q,g,hmax);
        wl[last] = MAX(p->F62,hc_out);

        for(int I=last+1; I<NI; ++I)
        wl[I] = wl[last];

        for(int I=last-1; I>=0; --I)
        {
            if(xs[I].zb.empty())
            {
                wl[I] = wl[I+1];
                continue;
            }

            int Id=I+1;
            while(Id<last && xs[Id].zb.empty())
            ++Id;

            double dsd;
            const double Kd = xs_conv(xs[Id],wl[Id],kst,dsd);
            const double Sfd = (Kd>0.0) ? Q*Q/(Kd*Kd) : 0.0;
            const double Ed = xs_energy(xs[Id],wl[Id],Q,g);

            auto res = [&](double h)
            {
                double dsu;
                const double Ku = xs_conv(xs[I],h,kst,dsu);
                const double Sfu = (Ku>0.0) ? Q*Q/(Ku*Ku) : 1.0e20;
                const double dx = 0.5*(dsu+dsd)*double(Id-I);
                return xs_energy(xs[I],h,Q,g) - Ed - 0.5*(Sfu+Sfd)*dx;
            };

            const double hc = xs_critical(xs[I],Q,g,hmax);
            double lo=hc, hi=MAX(hc,wl[Id])+1.0;

            if(res(lo)>=0.0)
            {
                wl[I] = hc;
                ++ncrit;
                continue;
            }

            int it=0;
            while(res(hi)<0.0 && it<60)
            {
                hi += 2.0*(hi-lo);
                ++it;
            }

            for(it=0; it<80 && hi-lo>1.0e-5; ++it)
            {
                const double mid = 0.5*(lo+hi);
                if(res(mid)<0.0)
                lo=mid;
                else
                hi=mid;
            }
            wl[I] = 0.5*(lo+hi);
        }

        // normal-flow level at the outflow for the log: slope of the lowest bed over the last tenth
        int Iu = MAX(0,last-MAX(NI/10,2));
        double dl=0.0;
        for(int I=Iu; I<last; ++I)
        dl += 0.5*(xs_dsmean(xs[I])+xs_dsmean(xs[I+1]));
        const double S0 = (dl>0.0) ? MAX((xs_zmin(xs[Iu])-xs_zmin(xs[last]))/dl,1.0e-6) : 1.0e-3;

        double lo=xs_zmin(xs[last])+1.0e-4, hi=lo+hmax, dsn;
        for(int it=0; it<80; ++it)
        {
            const double mid=0.5*(lo+hi);
            if(xs_conv(xs[last],mid,kst,dsn)*sqrt(S0)<Q)
            lo=mid;
            else
            hi=mid;
        }
        hn_out = 0.5*(lo+hi);

        if(verbose==1 && p->mpirank==0)
        {
            cout<<"F 65 2: steady 1D river flow, Q "<<Q<<" m3/s, ks "<<ks<<" m"<<(law==1?" (log law, bottom layer)":"")<<": water level "<<wl[0]<<" m at the inflow, "
                <<wl[last]<<" m at the outflow ("<<ncrit<<" critical sections)"<<endl;
            cout<<"F 65 2: normal-flow level at the outflow "<<hn_out<<" m (bed slope "<<S0<<"), F 62 "<<p->F62<<" m, critical level "<<hc_out<<" m"<<endl;
        }
    }

    if(verbose==1 && p->F65==1 && p->mpirank==0)
    cout<<"F 65 1: water level linear along the river, "<<wl[0]<<" m at the inflow, "<<wl[NI-1]<<" m at the outflow"<<endl;

    // eta of the cells
    // the cells of a column connected to its lowest cell below the level
    inrun.assign(NC,0);
    for(int I=0; I<NI; ++I)
    {
        int n0,n1;
        xs_run(xs[I],wl[I],n0,n1);
        for(int n=n0; n<=n1; ++n)
        inrun[long(I)*NJ + xs[I].J[n]] = 1;
    }


    return 1;
}

double river_ini::eta(lexer *p, int ii, int jj, double bed) const
{
    const int I = ii+p->origin_i;
    const long q = long(I)*NJ + long(jj+p->origin_j);
    return (inrun[q]==1 ? wl[I] : MIN(wl[I],bed)) - p->wd;
}

void river_ini::velocity(lexer *p, int ii, int jj, double bed, double &u, double &v) const
{
    // Manning strip velocity along the i-direction of the cell
    const int I = ii+p->origin_i;
    const double d = wl[I] - bed;
    const int qc = (ii-p->imin)*p->jmax + (jj-p->jmin);
    const int qw = qc - p->jmax;

    u = 0.0;
    v = 0.0;

    if(d<=(law==0 ? p->A244 : p->A544) || xs[I].zb.empty() || inrun[long(I)*NJ + long(jj+p->origin_j)]==0)
    return;

    double dsw;
    const double K = xs_conv(xs[I],wl[I],kst,dsw);
    if(K<=0.0)
    return;

    const double um = (law==0) ? kst*pow(d,2.0/3.0)*Q/K : strip(d)/d*Q/K;

    double tx = p->SX1[qw] + p->SX1[qc];
    double ty = p->SY1[qw] + p->SY1[qc];
    const double tl = sqrt(tx*tx + ty*ty);
    if(tl<1.0e-12)
    return;
    tx/=tl;
    ty/=tl;

    u = um*tx;
    v = um*ty*p->y_dir;
}
