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
Author: Widar W. Wang
--------------------------------------------------------------------*/

#include"wave_lib_spectrum.h"
#include"reef3d_mpmd.h"
#include"lexer.h"
#include"ghostcell.h"
#include"seastate_swan_spc.h"
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<ctime>
#include<filesystem>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<numeric>
#include<random>
#include<sstream>
#include<string>
#include<vector>
#include<sys/stat.h>
#include<sys/types.h>

/*--------------------------------------------------------------------
Wave components from a 2D spectrum: SWAN spectral file read directly
(B 85 12), and the 2D spectrum file of B 85 11 when B 143 or B 144 is
given.

SWAN file, read by rank 0 with seastate_swan_spc:
  SWAN.SPC or swan.spc in the case folder, otherwise the only *.spc
  B 140 n    location number, default 1
  B 141 x y  the location closest to x y, in the coordinates of the file
  B 142 n    record of a nonstationary file, default 1, -1 the last one
  Hs, Tp and the frequencies come from the file: B 86 and B 93 are not
  needed. The depth is not in the file: B 94 as usual.

Frame and limits, for B 85 11 with B 143 / B 144 as well:
  B 131      added to the directions [deg]: rotates the Cartesian
             direction of propagation of the spectrum into the frame of
             the model
  B 132 s e  sector kept, in the frame of the model, default -90 90: the
             waves that run into the domain from the generation at x = 0
  B 87 s e   frequency band kept [rad/s]; without B 87 the frequencies
             with wave lengths below 4 cells are dropped: the grid cannot
             hold them, in the generation zone they alias onto longer
             waves (the upper bins of a SWAN spectrum reach a few metres
             of wave length) and can blow up the free surface
  B 145 1    the rest scaled to the Hs of the full spectrum; default 0:
             what cannot enter the domain or stand on the grid is not
             generated and the log reports the Hs it leaves

Energy of a bin: a^2/2 = E(f,d) df dd with the trapezoidal widths of
the file, the integral SWAN and wavespectra use.

Components:
  B 143 N > 0  single summation (default, N = 512): N components of equal
               energy m0/N, each with its own frequency and direction.
               Frequencies: stratified in the cdf of E(w), E piecewise
               linear between the file's frequencies. Directions: the
               cdf of D(theta|w) inverted at a low-discrepancy sequence,
               so that every band of neighbouring frequencies holds the
               whole directional distribution; D blended between the two
               neighbouring frequencies, piecewise constant in the
               direction bins. sum a^2/2 = m0 exactly, and as no two
               components share a frequency the time-mean variance is
               the same at every point of the domain (double summation:
               it varies from point to point and does not average out).
  B 143 0      double summation: one component per bin of the file. The
               weakest bins that hold the fraction B 144 of the energy
               (default 0.001) are dropped and the rest scaled to keep m0.

  B 146 1      interpolate the spectrum to the bins of REEF3D first, then
               discretise it as for a theory spectrum; replaces B 143.
               E(w,theta) is bilinear between the bins of the file, so a
               file with few bins is enough (at the start only, on rank 0).
               B 136 4      B 86 components of equal energy, each with its
                            own frequency: B 86 must be a multiple of B 133.
                            Every B 133 neighbouring frequencies use all
                            B 133 equal energy directions of D(theta|w).
               B 136 1      B 86 x B 133 components, double summation
               B 136 2, 3   B 86 x B 133 components, single summation: the
                            frequency bin is split over the directions,
                            evenly (2) or at random (3)
               B 84 1 peak enhanced, 2 equal energy, 3 equal width bins of
               frequency, between B 87 or where the spectrum has energy.
               The B 133 directions (B 136 1 to 3) have equal width, in the
               sector B 132 where the spectrum has energy.

The components are built on rank 0 and broadcast (the random numbers
differ between the ranks with B 139 0).
--------------------------------------------------------------------*/

namespace
{
    const double spec2d_pi = 3.14159265358979323846;

    double spec2d_wavenumber(double w, double d)
    {
        const double g = 9.81;
        const double w2 = w*w;
        double k = fmax(w2/g, w/sqrt(g*d));

        for(int it=0; it<50; ++it)
        {
            const double kd = k*d;
            const double t = tanh(kd);
            const double ch = cosh(kd);
            const double dk = (g*k*t - w2)/(g*(t + kd/(ch*ch)));
            k -= dk;

            if(fabs(dk) < 1.0e-13*k)
            break;
        }

        return k;
    }

    // the SWAN file of the case folder
    std::string spec2d_find(std::string &error)
    {
        namespace fs = std::filesystem;

        for(const char *name : {"SWAN.SPC","swan.spc","SWAN.spc"})
        if(fs::exists(name))
        return name;

        std::vector<std::string> found;

        for(const auto &e : fs::directory_iterator("."))
        {
            std::string ext = e.path().extension().string();
            std::transform(ext.begin(),ext.end(),ext.begin(),[](unsigned char c){return std::tolower(c);});

            if(e.is_regular_file() && ext==".spc")
            found.push_back(e.path().filename().string());
        }

        if(found.size()==1)
        return found[0];

        error = found.empty() ? "no SWAN spectral file (SWAN.SPC, swan.spc or one *.spc) in the case folder"
                              : "several *.spc files in the case folder: name the one to read SWAN.SPC";
        return "";
    }

    // widths around the points, trapezoidal: half cells at the ends
    std::vector<double> spec2d_widths(const std::vector<double> &x)
    {
        const size_t n = x.size();
        std::vector<double> d(n,1.0);

        if(n>1)
        for(size_t i=0; i<n; ++i)
        d[i] = 0.5*(x[std::min(i+1,n-1)] - x[i==0 ? 0 : i-1]);

        return d;
    }

    // widths around sorted directions: periodic when the directions close the circle, a sector
    // (largest gap well above the others) has half bins at its ends
    std::vector<double> spec2d_widths_dir(const std::vector<double> &t)
    {
        const size_t n = t.size();
        std::vector<double> gap(n), d(n,2.0*spec2d_pi);

        if(n<2)
        return d;

        for(size_t i=0; i<n; ++i)
        {
            gap[i] = (i+1<n) ? t[i+1]-t[i] : t[0]+2.0*spec2d_pi-t[n-1];

            while(gap[i]<0.0)
            gap[i] += 2.0*spec2d_pi;
        }

        std::vector<double> g2(gap);
        std::nth_element(g2.begin(),g2.begin()+n/2,g2.end());
        const double med = g2[n/2];

        for(size_t i=0; i<n; ++i)
        {
            const double a = gap[(i+n-1)%n], b = gap[i];
            d[i] = 0.5*((a>1.5*med ? 0.0 : a) + (b>1.5*med ? 0.0 : b));
        }

        return d;
    }

    // E(w,theta) bilinear between the bins, on a fine even grid: frequencies between the first and
    // last bin with energy, one bin wider, directions over the circle. e is replaced by the energy of
    // the fine cells, with the sum of e kept
    void spec2d_refine(std::vector<double> &w, std::vector<double> &dwt, std::vector<double> &th, std::vector<double> &dth,
                       std::vector<double> &e, const int nwf, const int ntf)
    {
        const size_t nf = w.size(), nd = th.size();

        std::vector<size_t> o(nd);
        std::iota(o.begin(),o.end(),size_t(0));
        std::sort(o.begin(),o.end(),[&](size_t a, size_t b){return th[a]<th[b];});

        std::vector<double> ts(nd), gap(nd);
        double m0 = 0.0;

        for(size_t m=0; m<nd; ++m)
        ts[m] = th[o[m]];

        for(size_t m=0; m<nd; ++m)
        gap[m] = (m+1<nd) ? ts[m+1]-ts[m] : ts[0]+2.0*spec2d_pi-ts[nd-1];

        std::vector<double> g2(gap);
        std::nth_element(g2.begin(),g2.begin()+nd/2,g2.end());
        const double gmax = 1.5*g2[nd/2];

        // energy density, directions sorted
        std::vector<double> d(nf*nd,0.0);
        size_t l0 = nf, l1 = 0;

        for(size_t l=0; l<nf; ++l)
        for(size_t m=0; m<nd; ++m)
        {
            const double v = e[l*nd+o[m]];
            m0 += v;

            if(v>0.0)
            {
                d[l*nd+m] = v/(dwt[l]*dth[o[m]]);
                l0 = std::min(l0,l);
                l1 = std::max(l1,l);
            }
        }

        const double wa = w[l0>0 ? l0-1 : 0], wb = w[std::min(l1+1,nf-1)];
        const double hw = (wb-wa)/double(nwf-1), ht = 2.0*spec2d_pi/double(ntf);

        std::vector<double> wf(nwf), tf(ntf), ef(size_t(nwf)*ntf,0.0), row(nd);

        for(int i=0; i<nwf; ++i)
        wf[i] = wa + double(i)*hw;

        for(int j=0; j<ntf; ++j)
        tf[j] = -spec2d_pi + (double(j)+0.5)*ht;

        double sum = 0.0;

        for(int i=0; i<nwf; ++i)
        {
            size_t c = std::upper_bound(w.begin(),w.end(),wf[i]) - w.begin();
            c = std::min(std::max<size_t>(c,1),nf-1);

            const double x = fmin(fmax((wf[i]-w[c-1])/(w[c]-w[c-1]),0.0),1.0);

            for(size_t m=0; m<nd; ++m)
            row[m] = (1.0-x)*d[(c-1)*nd+m] + x*d[c*nd+m];

            for(int j=0; j<ntf; ++j)
            {
                // between two directions, the last gap closes the circle; nothing over a gap of a sector
                size_t m = std::upper_bound(ts.begin(),ts.end(),tf[j]) - ts.begin();
                m = (m==0) ? nd-1 : m-1;

                const size_t mn = (m+1)%nd;
                double dt = tf[j]-ts[m];

                if(dt<0.0)
                dt += 2.0*spec2d_pi;

                double v = 0.0;

                if(gap[m]<=gmax)
                v = row[m] + (row[mn]-row[m])*dt/gap[m];
                else if(dt<0.5*g2[nd/2])
                v = row[m];
                else if(gap[m]-dt<0.5*g2[nd/2])
                v = row[mn];

                ef[size_t(i)*ntf+j] = v;
                sum += v;
            }
        }

        if(sum>0.0)
        for(double &v : ef)
        v *= m0/sum;

        w.swap(wf);
        th.swap(tf);
        e.swap(ef);
        dwt = spec2d_widths(w);
        dth.assign(ntf,ht);
    }

    // the value where a rising sum first reaches the fraction q of its end
    double spec2d_quantile(const std::vector<double> &x, const std::vector<double> &cum, const double q)
    {
        const double v = q*cum.back();
        size_t c = std::lower_bound(cum.begin(),cum.end(),v) - cum.begin();
        c = std::min(std::max<size_t>(c,1),x.size()-1);

        const double h = cum[c]-cum[c-1];

        return h>0.0 ? x[c-1] + (x[c]-x[c-1])*(v-cum[c-1])/h : x[c];
    }
}

bool wave_lib_spectrum::spec2d_active(lexer *p)
{
    return p->B85==12 || (p->B85==11 && (p->B143>=0 || p->B144>=0.0 || p->B146==1));
}

// stop all ranks, the status comes from rank 0
void wave_lib_spectrum::spec2d_abort(lexer *p, int status, const std::string &msg)
{
    MPI_Bcast(&status,1,MPI_INT,0,reef3d_comm());

    if(status==0)
    return;

    if(p->mpirank==0)
    cout<<endl<<"!!! wave spectrum (B 85 "<<p->B85<<"): "<<msg<<" !!!"<<endl<<endl;

    MPI_Barrier(reef3d_comm());
    MPI_Finalize();
    exit(1);
}

void wave_lib_spectrum::spec2d_parameters(lexer *p)
{
    const double pi = spec2d_pi;
    const int rank = p->mpirank;

    // evaluated as a directional sea
    if(p->B130==0)
    p->B130 = 1;

    std::vector<double> w, dwt, th, dth, e;     // rad/s, rad, bin energy a^2/2 [m^2]
    std::string info;
    std::string error;
    int status = 0;

    // ---- the spectrum on rank 0 ------------------------------------------------
    if(p->B85==12 && rank==0)
    {
        seastate_swan_spc spc;
        seastate_swan_spc_select sel;
        sel.loc = p->B140;
        sel.time = p->B142;
        sel.nearest = p->B141>0;
        sel.x = p->B141_x;
        sel.y = p->B141_y;

        const std::string file = spec2d_find(error);

        if(file.empty() || !spc.read(file,error,sel))
        status = 1;
        else if(spc.nodata)
        {
        error = file+": location "+std::to_string(spc.loc)+" holds no spectrum in this record (ZERO or NODATA, a dry point?)";
        status = 1;
        }
        else
        {
            const size_t nf = spc.f.size(), nd = spc.dir.size();

            w.resize(nf);
            for(size_t l=0; l<nf; ++l)
            w[l] = 2.0*pi*spc.f[l];

            th.resize(nd);
            for(size_t m=0; m<nd; ++m)
            th[m] = spc.dir[m]*pi/180.0;

            dwt = spec2d_widths(w);
            dth = spec2d_widths_dir(th);

            e.resize(nf*nd);
            for(size_t l=0; l<nf; ++l)
            for(size_t m=0; m<nd; ++m)
            e[l*nd+m] = spc.E[l*nd+m]*(dwt[l]/(2.0*pi))*(dth[m]*180.0/pi);

            std::ostringstream s;
            s<<file<<": location "<<spc.loc<<" of "<<std::max<size_t>(spc.lx.size(),1)
             <<" ("<<spc.x<<", "<<spc.y<<"), record "<<spc.time;

            if(!spc.date.empty())
            s<<" ("<<spc.date<<")";

            s<<", "<<nf<<" frequencies "<<spc.f.front()<<" - "<<spc.f.back()<<" Hz, "<<nd<<" directions";
            info = s.str();
        }
    }

    if(p->B85==11)
    {
        // spectrum-file-2d.dat, all ranks; the energy of a bin as the wave generation of B 85 11 has it
        spectrum_file_2d_read(p);

        if(rank==0)
        {
            const int nf = ptnum_freq_2d, nd = ptnum_dir_2d;
            const double db = (nd>1) ? fabs(dir_2d[nd-1]-dir_2d[0])/double(nd-1) : 2.0*pi;

            w.assign(freq_2d,freq_2d+nf);
            th.assign(dir_2d,dir_2d+nd);
            dwt = spec2d_widths(w);
            dth.assign(nd,db);

            e.resize(size_t(nf)*nd);
            for(int l=0; l<nf; ++l)
            {
                const double dwf = (l<nf-1) ? freq_2d[l+1]-freq_2d[l] : freq_2d[l]-freq_2d[l-1];

                for(int m=0; m<nd; ++m)
                e[size_t(l)*nd+m] = spectrum_2d_1d[size_t(l)*nd+m]*dwf*db;
            }

            info = "spectrum-file-2d.dat: "+std::to_string(nf)+" frequencies, "+std::to_string(nd)+" directions";
        }
    }

    spec2d_abort(p,status,error);

    // ---- the components on rank 0 ----------------------------------------------
    std::vector<double> cw, cb, ce, cd;     // frequency, direction, energy, width
    double Hs_gen=0.0, Tp_file=0.0, ws_gen=0.0, we_gen=0.0;

    if(rank==0)
    {
        size_t nf = w.size(), nd = th.size();

        // frame of the model, directions in (-pi,pi]
        const double rot = p->B131*pi/180.0;
        const double sec_s = p->B132_s*pi/180.0, sec_e = p->B132_e*pi/180.0;

        for(size_t m=0; m<nd; ++m)
        th[m] = atan2(sin(th[m]+rot),cos(th[m]+rot));

        double m0_full=0.0, cx=0.0, sx=0.0;
        std::vector<double> rowd(nf,0.0);

        for(size_t l=0; l<nf; ++l)
        for(size_t m=0; m<nd; ++m)
        {
            const double v = e[l*nd+m];
            m0_full += v;
            rowd[l] += v/dwt[l];
            cx += v*cos(th[m]);
            sx += v*sin(th[m]);
        }

        if(!(m0_full>0.0))
        {
        status = 1;
        error = "the spectrum holds no energy";
        }
        else
        {
            const double dmean = atan2(sx,cx);
            const double r = sqrt(cx*cx+sx*sx)/m0_full;
            const double spread = sqrt(fmax(2.0*(1.0-r),0.0));
            const double wpk = w[std::max_element(rowd.begin(),rowd.end())-rowd.begin()];
            Tp_file = 2.0*pi/wpk;

            // the log format of B 85 11: the file can be used as spectrum-file-2d.dat
            {
                std::vector<int> o(nd);
                std::iota(o.begin(),o.end(),0);
                std::sort(o.begin(),o.end(),[&](int a, int b){return th[a]<th[b];});

                const double db = 2.0*pi/double(nd);

                mkdir("./REEF3D_Log-Wave",0777);
                ofstream out("./REEF3D_Log-Wave/REEF3D_wave-spectrum-2d.dat");
                out<<setprecision(10)<<"fq_dir";

                for(size_t m=0; m<nd; ++m)
                out<<" "<<th[o[m]];

                out<<endl;

                for(size_t l=0; l<nf; ++l)
                {
                    const double dwf = (l+1<nf) ? w[l+1]-w[l] : w[l]-w[l-1];
                    out<<w[l];

                    for(size_t m=0; m<nd; ++m)
                    out<<" "<<e[l*nd+o[m]]/(dwf*db);

                    out<<endl;
                }
            }

            // sector and band
            double m0_sec=0.0, m0_band=0.0;
            const double d0 = (p->B94==1) ? p->B94_wdt : p->phimean;
            const double kgrid = (p->DXM>0.0) ? 2.0*pi/(4.0*p->DXM) : 1.0e20;
            const double wgrid = sqrt(9.81*kgrid*tanh(fmin(kgrid*d0,20.0)));
            const double wa = p->B87==1 ? p->B87_1 : 0.0;
            const double wb = p->B87==1 ? p->B87_2 : wgrid;

            for(size_t l=0; l<nf; ++l)
            for(size_t m=0; m<nd; ++m)
            {
                double &v = e[l*nd+m];

                if(th[m]<sec_s || th[m]>sec_e)
                v = 0.0;

                m0_sec += v;

                if(w[l]<wa || w[l]>wb)
                v = 0.0;

                m0_band += v;
            }

            if(!(m0_band>0.0))
            {
            status = 1;
            error = "nothing of the spectrum is left in the sector B 132 and the band B 87: the mean direction is "
                  + std::to_string(dmean*180.0/pi)+" deg in the frame of the model, rotate with B 131";
            }
            else
            {
                double scale = 1.0;

                if(p->B145==1)
                scale = m0_full/m0_band;

                const bool bins = p->B146==1;
                const int nsingle = bins ? 0 : ((p->B143>=0) ? p->B143 : ((p->B85==12) ? 512 : 0));
                const int Nf = p->B86, Nd = p->B133;
                double repeat = 0.0, tbins = 0.0;

                if(bins)
                {
                    if(Nf<2 || Nd<2)
                    {
                    status = 1;
                    error = "B 146 1 uses the bins of B 86 and B 133: give B 86 and B 133, at least 2";
                    }
                    else if(p->B136==4 && Nf%Nd!=0)
                    {
                    status = 1;
                    error = "with B 136 4, B 86 ("+std::to_string(Nf)+") must be a multiple of B 133 ("+std::to_string(Nd)+"): the next ones are "
                          + std::to_string((Nf/Nd)*Nd)+" and "+std::to_string((Nf/Nd+1)*Nd);
                    }
                    else
                    {
                        tbins = MPI_Wtime();
                        spec2d_refine(w,dwt,th,dth,e,std::min(std::max(2048,8*Nf),8192),std::max(720,8*Nd));
                    }

                    nf = w.size();
                    nd = th.size();
                }

                if(status==0 && (nsingle>0 || bins))
                {
                    // energy density of the frequencies, linear in the cells between them
                    std::vector<double> row(nf,0.0), cum(nf,0.0);

                    for(size_t l=0; l<nf; ++l)
                    {
                        for(size_t m=0; m<nd; ++m)
                        row[l] += e[l*nd+m];

                        row[l] /= dwt[l];
                    }

                    for(size_t c=1; c<nf; ++c)
                    cum[c] = cum[c-1] + 0.5*(row[c-1]+row[c])*(w[c]-w[c-1]);

                    const double tot = cum[nf-1];

                    std::mt19937_64 gen(p->B139>0 ? (unsigned long long)p->B139 : (unsigned long long)time(0));
                    std::uniform_real_distribution<double> uni(0.0,1.0);

                    std::vector<double> mix(nd), mcum(nd+1);

                    // the frequency at u of the cdf of E(w), the direction at the share vs of D(theta|w)
                    auto component = [&](const double u, const double vs, double &cwn, double &cbn)
                    {
                        size_t c = std::upper_bound(cum.begin(),cum.end(),u) - cum.begin();
                        c = std::min(std::max<size_t>(c,1),nf-1);

                        // zero-energy cells have no width in the cdf
                        const double h = w[c]-w[c-1];
                        const double a = row[c-1], b = row[c];
                        const double du = u-cum[c-1];
                        double x;

                        if(fabs(b-a) > 1.0e-12*(a+b))
                        x = (-a + sqrt(fmax(a*a + 2.0*(b-a)*du/h,0.0)))/(b-a);
                        else
                        x = (a+b>0.0) ? du/(h*0.5*(a+b)) : 0.5;

                        x = fmin(fmax(x,0.0),1.0);
                        cwn = w[c-1] + x*h;

                        // D(theta|w), between the two rows with the weights of their share of E(w)
                        const double ra = (1.0-x)*a, rb = x*b;
                        const double wga = (ra+rb>0.0) ? ra/(ra+rb) : 0.5;
                        const double sa = a*dwt[c-1], sb = b*dwt[c];

                        mcum[0] = 0.0;
                        for(size_t m=0; m<nd; ++m)
                        {
                            const double da = sa>0.0 ? e[(c-1)*nd+m]/sa : 0.0;
                            const double db = sb>0.0 ? e[c*nd+m]/sb : 0.0;
                            mcum[m+1] = mcum[m] + wga*da + (1.0-wga)*db;
                        }

                        const double v = vs*mcum[nd];

                        size_t j = std::upper_bound(mcum.begin(),mcum.end(),v) - mcum.begin();
                        j = std::min(std::max<size_t>(j,1),nd) - 1;

                        const double pj = mcum[j+1]-mcum[j];
                        const double y = pj>0.0 ? (v-mcum[j])/pj : 0.5;

                        cbn = th[j] + (y-0.5)*dth[j];
                    };

                    if(!bins)
                    {
                        const double golden = 0.5*(sqrt(5.0)-1.0);
                        const double v0 = uni(gen);

                        cw.resize(nsingle);
                        cb.resize(nsingle);
                        ce.assign(nsingle,m0_band/double(nsingle));
                        cd.assign(nsingle,tot/double(nsingle));

                        for(int n=0; n<nsingle; ++n)
                        {
                            // stratified: one component in each interval of equal energy
                            const double u = tot*(double(n)+0.25+0.5*uni(gen))/double(nsingle);

                            component(u,fmod(v0 + double(n)*golden,1.0),cw[n],cb[n]);
                        }
                    }
                    else if(p->B136==4)
                    {
                        // equal energy in frequency and direction: the middle of each of the B 86 intervals of
                        // the cdf of E(w); the B 133 intervals of D(theta|w) in a new random order for every
                        // B 133 neighbouring frequencies, so that each band of them holds all the directions
                        cw.resize(Nf);
                        cb.resize(Nf);
                        ce.assign(Nf,m0_band/double(Nf));
                        cd.assign(Nf,tot/double(Nf));

                        std::vector<int> lvl(Nd);
                        std::iota(lvl.begin(),lvl.end(),0);

                        for(int n=0; n<Nf; ++n)
                        {
                            if(n%Nd==0)
                            std::shuffle(lvl.begin(),lvl.end(),gen);

                            component(tot*(double(n)+0.5)/double(Nf),(double(lvl[n%Nd])+0.5)/double(Nd),cw[n],cb[n]);
                        }
                    }
                    else
                    {
                        // the bins: edges of frequency and direction, the energy of the fine cells summed in them
                        std::vector<double> mar(nd,0.0), cmar(nd,0.0), cume(nf,0.0);

                        for(size_t l=0; l<nf; ++l)
                        for(size_t m=0; m<nd; ++m)
                        mar[m] += e[l*nd+m];

                        for(size_t m=0; m<nd; ++m)
                        cmar[m] = mar[m] + (m>0 ? cmar[m-1] : 0.0);

                        const double fa = p->B87==1 ? fmax(p->B87_1,w[0]) : spec2d_quantile(w,cum,1.0e-4);
                        const double fb = p->B87==1 ? fmin(p->B87_2,w[nf-1]) : spec2d_quantile(w,cum,1.0-1.0e-4);
                        const double ta = spec2d_quantile(th,cmar,1.0e-4), tb = spec2d_quantile(th,cmar,1.0-1.0e-4);
                        const double dbn = (tb-ta)/double(Nd);

                        std::vector<double> fe(Nf+1);

                        if(p->B84==2)
                        {
                            for(int n=0; n<=Nf; ++n)
                            fe[n] = spec2d_quantile(w,cum,double(n)/double(Nf));

                            fe[0] = w[0];
                            fe[Nf] = w[nf-1];
                        }
                        else
                        {
                            // B 84 1: the peak on an edge, as for a theory spectrum; B 84 3: equal widths
                            const double fp = w[std::max_element(row.begin(),row.end())-row.begin()];
                            int ns = int((fp-fa)/(fb-fa)*double(Nf)+0.5);

                            if(p->B84!=1 || fp<=fa || fp>=fb)
                            ns = 0;

                            ns = (ns>0) ? std::min(std::max(ns,1),Nf-1) : 0;

                            for(int n=0; n<=Nf; ++n)
                            fe[n] = (ns>0) ? (n<=ns ? fa + (fp-fa)*double(n)/double(ns) : fp + (fb-fp)*double(n-ns)/double(Nf-ns))
                                           : fa + (fb-fa)*double(n)/double(Nf);

                            if(ns>0)
                            repeat = 2.0*pi/fmin((fp-fa)/double(ns),(fb-fp)/double(Nf-ns));
                            else
                            repeat = 2.0*pi/((fb-fa)/double(Nf));
                        }

                        std::vector<double> eb(size_t(Nf)*Nd,0.0), ew(Nf,0.0), wbar(Nf,0.0);

                        for(size_t l=0; l<nf; ++l)
                        {
                            int n = int(std::upper_bound(fe.begin(),fe.end(),w[l]) - fe.begin()) - 1;
                            n = std::min(std::max(n,0),Nf-1);

                            for(size_t m=0; m<nd; ++m)
                            {
                                const double v = e[l*nd+m];

                                if(v>0.0)
                                {
                                    const int q = std::min(std::max(int(floor((th[m]-ta)/dbn)),0),Nd-1);
                                    eb[size_t(n)*Nd+q] += v;
                                    ew[n] += v;
                                    wbar[n] += v*w[l];
                                }
                            }
                        }

                        for(int n=0; n<Nf; ++n)
                        for(int q=0; q<Nd; ++q)
                        {
                            const double v = eb[size_t(n)*Nd+q];

                            if(v>0.0)
                            {
                                const double width = fe[n+1]-fe[n];
                                double wn = wbar[n]/ew[n];

                                if(p->B136==2)
                                wn = fe[n] + (double(q)+0.5)*width/double(Nd);

                                if(p->B136==3)
                                wn = fe[n] + (double(q)+uni(gen))*width/double(Nd);

                                cw.push_back(wn);
                                cb.push_back(ta + (double(q)+0.5)*dbn);
                                ce.push_back(v);
                                cd.push_back(width);
                            }
                        }
                    }
                }
                else if(status==0)
                {
                    // double summation: the bins, without the weakest ones
                    const double eps = (p->B144>=0.0) ? p->B144 : 0.001;

                    std::vector<size_t> o(nf*nd);
                    std::iota(o.begin(),o.end(),size_t(0));
                    std::sort(o.begin(),o.end(),[&](size_t a, size_t b){return e[a]>e[b];});

                    double sum=0.0;
                    size_t kept=0;

                    for(; kept<o.size() && e[o[kept]]>0.0; ++kept)
                    {
                        if(sum>=(1.0-eps)*m0_band)
                        break;

                        sum += e[o[kept]];
                    }

                    std::sort(o.begin(),o.begin()+kept);

                    for(size_t q=0; q<kept; ++q)
                    {
                        const size_t l = o[q]/nd, m = o[q]%nd;
                        cw.push_back(w[l]);
                        cb.push_back(th[m]);
                        ce.push_back(e[o[q]]*m0_band/sum);
                        cd.push_back(dwt[l]);
                    }
                }

                if(status==0)
                {
                double m0_gen=0.0;

                for(double &v : ce)
                {
                    v *= scale;
                    m0_gen += v;
                }

                Hs_gen = 4.0*sqrt(m0_gen);
                ws_gen = *std::min_element(cw.begin(),cw.end());
                we_gen = *std::max_element(cw.begin(),cw.end());

                // energy on wave lengths below 8 cells: it cannot be carried by the grid
                double m0_short=0.0;

                if(p->DXM>0.0)
                for(size_t n=0; n<cw.size(); ++n)
                if(2.0*pi/spec2d_wavenumber(cw[n],d0) < 8.0*p->DXM)
                m0_short += ce[n];

                const int prec = int(cout.precision());
                cout<<setprecision(5);
                cout<<"2D spectrum "<<info<<endl;
                cout<<"  Hs "<<4.0*sqrt(m0_full)<<" m, Tp "<<Tp_file<<" s, mean direction "<<dmean*180.0/pi
                    <<" deg in the model frame, spreading "<<spread*180.0/pi<<" deg"<<endl;
                cout<<"  Hs in the sector B 132 "<<4.0*sqrt(m0_sec)<<" m, also in the band "<<wa<<" - "<<(wb<1.0e10 ? wb : w[nf-1])<<" rad/s "
                    <<(p->B87==1 ? "(B 87) " : "(wave lengths from 4 cells) ")<<4.0*sqrt(m0_band)<<" m, generated "<<Hs_gen<<" m"<<endl;

                if(bins)
                {
                    double cg=0.0, sg=0.0;

                    for(size_t n=0; n<cw.size(); ++n)
                    {
                        cg += ce[n]*cos(cb[n]);
                        sg += ce[n]*sin(cb[n]);
                    }

                    const double rg = sqrt(cg*cg+sg*sg)/m0_gen;

                    cout<<"  interpolated to B 86 "<<Nf<<" and B 133 "<<Nd<<" (B 84 "<<p->B84<<", B 136 "<<p->B136<<") in "<<MPI_Wtime()-tbins<<" s: "<<cw.size()<<" components, mean direction "
                        <<atan2(sg,cg)*180.0/pi<<" deg, spreading "<<sqrt(fmax(2.0*(1.0-rg),0.0))*180.0/pi<<" deg"<<endl;

                    if(repeat>0.0 && repeat<p->N41)
                    cout<<"  !!! the frequency bins repeat the wave field after "<<repeat<<" s, shorter than N 41: more B 86, or B 84 2 !!!"<<endl;
                }
                else if(nsingle>0)
                cout<<"  "<<cw.size()<<" components, single summation"<<endl;
                else
                cout<<"  "<<cw.size()<<" components of "<<nf*nd<<" bins, double summation"<<endl;

                if(m0_sec<0.95*m0_full)
                cout<<"  !!! "<<100.0*(1.0-m0_sec/m0_full)<<" % of the energy runs out of the domain and is not generated (sector B 132); B 131 rotates, B 145 1 scales to the full Hs !!!"<<endl;

                if(m0_short>0.05*m0_gen)
                cout<<"  !!! "<<100.0*m0_short/m0_gen<<" % of the energy is on wave lengths below 8 cells; it is poorly resolved on this grid: refine, or limit the band with B 87 and scale with B 145 1 !!!"<<endl;

                cout<<setprecision(prec);
                }
            }
        }
    }

    spec2d_abort(p,status,error);

    // ---- all ranks -------------------------------------------------------------
    int ncomp = int(cw.size());
    MPI_Bcast(&ncomp,1,MPI_INT,0,reef3d_comm());

    cw.resize(ncomp);
    cb.resize(ncomp);
    ce.resize(ncomp);
    cd.resize(ncomp);

    MPI_Bcast(cw.data(),ncomp,MPI_DOUBLE,0,reef3d_comm());
    MPI_Bcast(cb.data(),ncomp,MPI_DOUBLE,0,reef3d_comm());
    MPI_Bcast(ce.data(),ncomp,MPI_DOUBLE,0,reef3d_comm());
    MPI_Bcast(cd.data(),ncomp,MPI_DOUBLE,0,reef3d_comm());

    double hdr[4] = {Hs_gen,Tp_file,ws_gen,we_gen};
    MPI_Bcast(hdr,4,MPI_DOUBLE,0,reef3d_comm());

    numcomp = ncomp;
    p->wN = ncomp;

    p->Darray(Si,numcomp);
    p->Darray(Sn,numcomp);
    p->Darray(Di,numcomp);
    p->Darray(Di_n,numcomp);
    p->Darray(wi,numcomp);
    p->Darray(dw,numcomp);
    p->Darray(Ai,numcomp);
    p->Darray(Li,numcomp);
    p->Darray(ki,numcomp);
    p->Darray(Ti,numcomp);
    p->Darray(ei,numcomp);
    p->Darray(beta,numcomp);
    p->Darray(beta_n,numcomp);
    p->Darray(cosbeta,numcomp);
    p->Darray(sinbeta,numcomp);

    // amplitudes_irregular: sqrt(2 S dw dbeta), S dw the energy of the component
    dbeta = 1.0;

    for(int n=0; n<numcomp; ++n)
    {
        wi[n] = cw[n];
        beta[n] = cb[n];
        cosbeta[n] = cos(cb[n]);
        sinbeta[n] = sin(cb[n]);
        dw[n] = cd[n];
        Si[n] = ce[n]/cd[n];
        Sn[n] = Si[n];
        Di[n] = 1.0;
        Di_n[n] = 1.0;
        beta_n[n] = cb[n];
        ki[n] = spec2d_wavenumber(cw[n],wD);
        Li[n] = 2.0*pi/ki[n];
    }

    // the sea state of the file replaces B 93
    p->wHs = hdr[0];
    p->wAs = 0.5*hdr[0];
    p->wA = p->wAs;
    p->wTp = hdr[1];
    p->wT = hdr[1];
    p->wwp = 2.0*pi/hdr[1];
    p->wLp = p->wTp*sqrt(9.81*wD);

    wp = p->wwp;
    ws = hdr[2];
    we = hdr[3];

    // peak wave speed
    p->wC = p->wwp/spec2d_wavenumber(p->wwp,wD);
}
