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

#include"seastate_f.h"
#include"seastate_amr.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<fstream>
#include<iomanip>
#include<sstream>
#include<string>
#include<vector>
#include<mpi.h>
#include<sys/stat.h>
#include<sys/types.h>

/*--------------------------------------------------------------------
Handover of 2D spectra to the phase-resolved wave generation (A 760)

For every point A 760 x y the spectrum of the cell that contains the
point is written to

  REEF3D_SEASTATE_Spectra/spectrum-file-2d_P<n>.dat

in the format of REEF3D's 2D spectrum input (B 85 11, file
spectrum-file-2d.dat, wave_lib_spectrum_file.cpp):

  fq_dir  theta_1 ... theta_M              [rad]
  omega_1 S(omega_1,theta_1) ... S(omega_1,theta_M)
  ...                                      [rad/s, m^2 s/rad^2]

theta is the direction of propagation, counter-clockwise from +x (the
convention of SEASTATE and of the wave generation, cos(beta) x +
sin(beta) y). The directions are the sector of half-width A 761 around
the grid direction closest to the mean direction, written unwrapped and
increasing (they may be < 0 or > 2 pi). The wave generation builds one
component per (omega, theta) with amplitude sqrt(2 S domega dbeta),
domega the forward difference of the omegas (the last one repeated),
dbeta the mean direction step. S is therefore written as
E(sig,theta) dsig/domega, so that every component carries exactly the
energy of its SEASTATE bin: the Hs of the generated waves equals the
Hs of the sector.

A companion file seastate-wavegen_P<n>.txt lists the location, the
cell, the water depth, Hs of the full spectrum and of the sector, Tp,
the mean direction, and the ctrl keys for the wave generation.

The files are rewritten after every step (the latest spectrum).
--------------------------------------------------------------------*/

void seastate_f::handover(lexer *p, ghostcell *pgc)
{
    const seastate_grid &g = *e->grid;
    const double pi = 3.14159265358979323846;

    if(p->mpirank==0)
    mkdir("./REEF3D_SEASTATE_Spectra",0777);

    pgc->globalmax(0.0);    // the folder exists before any rank writes

    for(int n=0; n<p->A760; ++n)
    {
    const double xp = p->A760_x[n], yp = p->A760_y[n];

    // the rank that owns the cell containing the point writes the file
    int owner = 0, ci=-1, cj=-1;

        ILOOP
        JLOOP
        if(p->XN[IP]<=xp && xp<p->XN[IP1] && p->YN[JP]<=yp && yp<p->YN[JP1])
        {
        owner = 1;
        ci = i;
        cj = j;
        }

    const int found = int(pgc->globalisum(owner));

        if(found==0)
        {
            if(p->mpirank==0)
            cout<<"SEASTATE handover point "<<n+1<<" ("<<xp<<", "<<yp<<") is outside the domain"<<endl;

            continue;
        }

        if(owner==0)
        continue;

    // mesh refinement: the finest grid that holds the point
    lexer *q = p;
    fdm_seastate *ee = e;

    if(pamr!=nullptr)
    pamr->locate(xp,yp,q,ee,ci,cj);

    i = ci;
    j = cj;

    char name[256];
    sprintf(name,"./REEF3D_SEASTATE_Spectra/spectrum-file-2d_P%i.dat",n+1);
    ofstream out(name);

    char iname[256];
    sprintf(iname,"./REEF3D_SEASTATE_Spectra/seastate-wavegen_P%i.txt",n+1);
    ofstream info(iname);

        if(ee->wet(i,j)==0 || ee->N->spec(i,j)==nullptr)
        {
        info<<"REEF3D::SEASTATE handover point "<<n+1<<" ("<<xp<<", "<<yp<<"): the cell is dry or land, no spectrum"<<endl;
        continue;
        }

    const float *N = ee->N->spec(i,j);

    seastate_param sp;
    sp.compute(g,N);

    // direction sector around the grid direction closest to the mean direction; with a fine direction
    // sector (A 715) the spectrum is written on uniform directions of the smallest bin width, linear
    // in direction between the bin centres
    const double dto = g.uniform ? g.dtheta : *std::min_element(g.dth.begin(),g.dth.end());
    const int nout = int(std::lround(2.0*pi/dto));
    const double t0 = g.uniform ? g.theta0 : g.theta[g.direction_bin(sp.dir*pi/180.0)];
    double dr = std::fmod(sp.dir*pi/180.0-g.theta0,2.0*pi);
    if(dr<0.0)
    dr += 2.0*pi;
    const int m0 = g.uniform ? int(std::lround(dr/g.dtheta)) % nout : 0;
    int nh = int(std::floor(p->A761*pi/180.0/dto + 1.0e-9));
    nh = std::max(nh,1);

    // energy density of frequency l at the written direction m
    auto Eout = [&](int l, int m)
    {
        if(g.uniform)
        {
        // full circle: nout = ndir; compact computational sector (A 717): nothing beyond the sector
        const int mm = ((m%nout)+nout)%nout;
        return mm<g.ndir ? g.sig[l]*double(N[g.bin(l,mm)]) : 0.0;
        }

        double t = std::fmod(t0 + double(m)*dto,2.0*pi);
        if(t<0.0)
        t += 2.0*pi;

        // bracketing bin centres (periodic)
        int b = int(std::upper_bound(g.theta.begin(),g.theta.end(),t) - g.theta.begin());
        const int m1 = (b==0) ? g.ndir-1 : b-1;
        const int m2 = (b==g.ndir) ? 0 : b;
        double d1 = g.theta[m1], d2 = g.theta[m2];
        if(d1>t) d1 -= 2.0*pi;
        if(d2<t) d2 += 2.0*pi;
        const double w = (d2>d1) ? (t-d1)/(d2-d1) : 0.0;

        return g.sig[l]*((1.0-w)*double(N[g.bin(l,m1)]) + w*double(N[g.bin(l,m2)]));
    };

    int ma, mb;
        if(2*nh+1>=nout)
        {
        ma = m0 - nout/2;
        mb = ma + nout - 1;
        }
        else
        {
        ma = m0 - nh;
        mb = m0 + nh;
        }

    // domega of the wave generation: forward difference, the last one repeated
    std::vector<double> dw(g.nsig);
    for(int l=0; l<g.nsig-1; ++l)
    dw[l] = g.sig[l+1]-g.sig[l];
    dw[g.nsig-1] = dw[g.nsig-2];

    out<<"fq_dir";
    out<<setprecision(10);
    for(int m=ma; m<=mb; ++m)
    out<<" "<<t0+double(m)*dto;
    out<<endl;

    double m0sec = 0.0;

        for(int l=0; l<g.nsig; ++l)
        {
        out<<g.sig[l];

            for(int m=ma; m<=mb; ++m)
            {
            const double E = Eout(l,m);

            out<<" "<<E*g.dsig[l]/dw[l];
            m0sec += E*g.dsig[l]*dto;
            }

        out<<endl;
        }

    out.close();

    const double hs_sec = 4.0*std::sqrt(std::max(m0sec,0.0));

    info<<"REEF3D::SEASTATE handover point "<<n+1<<endl;
    info<<"location "<<xp<<" "<<yp<<", cell centre "<<q->XP[IP]<<" "<<q->YP[JP]<<", cell size "<<q->DXN[IP]<<" m, simtime "<<p->simtime<<endl;
    info<<"water depth "<<ee->depth(i,j)<<" m"<<endl;
    info<<"Hs "<<sp.Hs<<" m (full spectrum), "<<hs_sec<<" m (written sector "<<(t0+double(ma)*dto)*180.0/pi<<" to "<<(t0+double(mb)*dto)*180.0/pi<<" deg)"<<endl;
    info<<"Tp "<<sp.Tp<<" s, Tm01 "<<sp.Tm01<<" s, mean direction "<<sp.dir<<" deg, spread "<<sp.spread<<" deg"<<endl;
    info<<endl;
    info<<"wave generation (FNPF/NHFLOW): copy spectrum-file-2d_P"<<n+1<<".dat to spectrum-file-2d.dat and use"<<endl;
    info<<"B 85 11"<<endl;
    info<<"B 92 31          (linear; 32 or 33 for second order, with B 89 0)"<<endl;
    info<<"B 93 "<<sp.Hs<<" "<<sp.Tp<<endl;
    info<<"B 94 "<<ee->depth(i,j)<<endl;
    info<<"B 130 1          (needed with B 85 11; the spreading comes from the file)"<<endl;
    info<<"directions are absolute (propagation, ccw from +x): rotate the SEASTATE frame into the phase-resolved frame if needed"<<endl;
    info.close();
    }
}

/*--------------------------------------------------------------------
Nesting output (A 785 1): 2D spectra at many points for the boundary of
a nested run

The points are read from seastate-spcout.dat (x y per line in model
coordinates, '#' starts a comment). Every A 786 seconds of model time
(0: every step; the first record at the start of the run) the spectra
at the points are appended to

  REEF3D_SEASTATE_SPC/seastate-spcout.spc

a nonstationary SWAN 2D spectral file: LOCATIONS (model coordinates),
AFREQ (the frequencies of the spectral grid), CDIR (its directions of
propagation, Cartesian), QUANT VaDens [m2/Hz/degr], then per record the
date-time line YYYYMMDD.HHMMSS (model time 0 + tref, the clock of the
forcing files: W 203 / A 780) and one FACTOR, ZERO or NODATA block per
point. The spectrum of a point is bilinear between the four cell
centres around it (level 0), over the wet cells with spectra (weights
renormalised); NODATA where none is wet or the point is outside.

A nested run with the same model frame (W 205, W 206) takes the file as
its seastate-boundary.spc (A 711 3, A 712 per side): the spectra of its
boundary cells come from the two nearest points, linear in time.
--------------------------------------------------------------------*/

namespace
{
    // YYYYMMDD.HHMMSS of t [s since 1970-01-01, UTC], proleptic Gregorian (H. Hinnant)
    std::string spc_stamp(double t)
    {
        const long s = long(std::floor(t + 0.5));
        long z = (s>=0 ? s : s-86399)/86400;
        const long sod = s - z*86400;
        z += 719468;
        const long era = (z>=0 ? z : z-146096)/146097;
        const long doe = z - era*146097;
        const long yoe = (doe - doe/1460 + doe/36524 - doe/146096)/365;
        const long doy = doe - (365*yoe + yoe/4 - yoe/100);
        const long mp = (5*doy + 2)/153;
        const long d = doy - (153*mp + 2)/5 + 1;
        const long m = mp + (mp<10 ? 3 : -9);
        const long y = yoe + era*400 + (m<=2 ? 1 : 0);

        char buf[64];
        snprintf(buf,sizeof(buf),"%04ld%02ld%02ld.%02ld%02ld%02ld",y,m,d,sod/3600,(sod/60)%60,sod%60);
        return std::string(buf);
    }
}

void seastate_f::spcout_ini(lexer *p, ghostcell *pgc)
{
    const seastate_grid &g = *e->grid;
    const double pi = 3.14159265358979323846;

    spc_x.clear();
    spc_y.clear();

    std::ifstream in("seastate-spcout.dat");
    std::string line;

    while(std::getline(in,line))
    {
    std::istringstream s(line.substr(0,line.find('#')));
    double x, y;
        if(s>>x>>y)
        {
        spc_x.push_back(x);
        spc_y.push_back(y);
        }
    }

    if(spc_x.empty())
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  A 785 1: no points in seastate-spcout.dat (x y per line, model coordinates)"<<endl<<endl;

    pgc->final(true);
    }

    // a restart appends to the file of the run it continues
    const bool append = restart && std::ifstream("./REEF3D_SEASTATE_SPC/seastate-spcout.spc").good();

    if(p->mpirank==0 && !append)
    {
    mkdir("./REEF3D_SEASTATE_SPC",0777);

    ofstream out("./REEF3D_SEASTATE_SPC/seastate-spcout.spc");
    out<<"SWAN   1                                Swan standard spectral file, version"<<endl;
    out<<"$   Data produced by REEF3D::SEASTATE (A 785 1), model coordinates"<<endl;
    out<<"TIME                                    time-dependent data"<<endl;
    out<<"     1                                  time coding option"<<endl;
    out<<"LOCATIONS                               locations in x-y-space"<<endl;
    out<<setw(6)<<spc_x.size()<<"                                  number of locations"<<endl;
    out<<setprecision(12);
        for(size_t n=0; n<spc_x.size(); ++n)
        out<<"   "<<spc_x[n]<<"   "<<spc_y[n]<<endl;
    out<<"AFREQ                                   absolute frequencies in Hz"<<endl;
    out<<setw(6)<<g.nsig<<"                                  number of frequencies"<<endl;
        for(int l=0; l<g.nsig; ++l)
        out<<"   "<<g.sig[l]/(2.0*pi)<<endl;
    out<<"CDIR                                    spectral Cartesian directions in degr"<<endl;
    out<<setw(6)<<g.ndir<<"                                  number of directions"<<endl;
        for(int m=0; m<g.ndir; ++m)
        {
        double d = std::fmod(g.theta[m]*180.0/pi,360.0);
        if(d<0.0)
        d += 360.0;
        out<<"   "<<d<<endl;
        }
    out<<"QUANT"<<endl;
    out<<"     1                                  number of quantities in table"<<endl;
    out<<"VaDens                                  variance densities in m2/Hz/degr"<<endl;
    out<<"m2/Hz/degr                              unit"<<endl;
    out<<"   -0.9900E+02                          exception value"<<endl;
    }

    if(p->mpirank==0)
    {
    cout<<"SEASTATE nesting output (A 785 1): "<<spc_x.size()<<" points of seastate-spcout.dat every "<<p->A786
        <<" s -> REEF3D_SEASTATE_SPC/seastate-spcout.spc"<<endl;
    }

    pgc->globalmax(0.0);    // the header exists before the first record

    // the first record at the start; after a restart the next one on the schedule of the run (times increase)
    spc_tnext = p->simtime;

    if(append && p->A786>0.0)
    spc_tnext = (std::floor(p->simtime/p->A786 + 1.0e-9) + 1.0)*p->A786;
}

void seastate_f::spcout(lexer *p, ghostcell *pgc, bool force)
{
    if(spc_x.empty())
    return;

    // due: the first record, then every A 786 s (0: every step)
    const double eps = 1.0e-6*std::max(p->dt,1.0);

    if(!force && p->simtime<spc_tnext-eps)
    return;

    if(p->A786>0.0)
    while(spc_tnext<=p->simtime+eps)
    spc_tnext += p->A786;

    const seastate_grid &g = *e->grid;
    const double pi = 3.14159265358979323846;
    const double fac = 2.0*pi*pi/180.0;     // E(sig,theta[rad]) -> E(f,theta[deg])
    const int nb = g.nsig*g.ndir;
    const int np = int(spc_x.size());

    // per point: nb densities and a flag (1 wet, 0.5 dry or land)
    std::vector<double> buf(size_t(np)*(nb+1),0.0), all;

    for(int n=0; n<np; ++n)
    {
    const double xp = spc_x[n], yp = spc_y[n];
    int ci=-1, cj=-1;

        ILOOP
        JLOOP
        if(p->XN[IP]<=xp && xp<p->XN[IP1] && p->YN[JP]<=yp && yp<p->YN[JP1])
        {
        ci = i;
        cj = j;
        }

        if(ci<0)
        continue;

    // the four cell centres around the point (ghost cells of the neighbour ranks included)
    const int i0 = (xp>=p->XP[ci+marge]) ? ci : ci-1;
    const int j0 = (p->j_dir==0 || yp>=p->YP[cj+marge]) ? cj : cj-1;
    const double wx = std::min(std::max((xp - p->XP[i0+marge])/(p->XP[i0+1+marge] - p->XP[i0+marge]),0.0),1.0);
    const double wy = (p->j_dir==0) ? 0.0 : std::min(std::max((yp - p->YP[j0+marge])/(p->YP[j0+1+marge] - p->YP[j0+marge]),0.0),1.0);

    double *b = &buf[size_t(n)*(nb+1)];
    double wsum = 0.0;

        for(int a=0; a<2; ++a)
        for(int c=0; c<2; ++c)
        {
        const double w = (a==0 ? 1.0-wx : wx)*(c==0 ? 1.0-wy : wy);

            if(w<=0.0)
            continue;

        const int ii = i0+a, jj = j0+c;

            if(e->wet(ii,jj)!=1 || e->N->spec(ii,jj)==nullptr)
            continue;

        const float *N = e->N->spec(ii,jj);

            for(int l=0; l<g.nsig; ++l)
            for(int m=0; m<g.ndir; ++m)
            b[l*g.ndir+m] += w*g.sig[l]*double(N[g.bin(l,m)])*fac;

        wsum += w;
        }

        if(wsum>0.0)
        {
        for(int k=0; k<nb; ++k)
        b[k] /= wsum;
        b[nb] = 1.0;
        }
        else
        {
        for(int k=0; k<nb; ++k)
        b[k] = 0.0;
        b[nb] = 0.5;
        }
    }

    if(p->mpirank==0)
    all.assign(buf.size(),0.0);

    MPI_Reduce(buf.data(),p->mpirank==0 ? all.data() : nullptr,int(buf.size()),MPI_DOUBLE,MPI_SUM,0,pgc->mpi_comm);

    if(p->mpirank!=0)
    return;

    ofstream out("./REEF3D_SEASTATE_SPC/seastate-spcout.spc",ios_base::app);
    out<<spc_stamp(tref + p->simtime)<<"                         date and time"<<endl;

    for(int n=0; n<np; ++n)
    {
    const double *b = &all[size_t(n)*(nb+1)];

        if(b[nb]<0.9)
        {
        out<<"NODATA"<<endl;
        continue;
        }

    double vmax = 0.0;
    for(int k=0; k<nb; ++k)
    vmax = std::max(vmax,b[k]);

        if(!(vmax>0.0))
        {
        out<<"ZERO"<<endl;
        continue;
        }

    const double factor = vmax/99999.0;
    out<<"FACTOR"<<endl;
    out<<"   "<<setprecision(8)<<std::scientific<<factor<<std::defaultfloat<<endl;

        for(int l=0; l<g.nsig; ++l)
        {
            for(int m=0; m<g.ndir; ++m)
            out<<" "<<setw(6)<<std::lround(b[l*g.ndir+m]/factor);
        out<<endl;
        }
    }
}
