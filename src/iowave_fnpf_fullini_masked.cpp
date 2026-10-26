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

// I 30 2: full initialization of the FNPF wave field with exposure masks.
//
// As I 30 1, the whole domain is initialized with the wave theory of B 92, but every wave
// component is multiplied by an exposure mask for its propagation direction:
//
//  1. the wet/dry state of all cells is gathered on a coarse global grid (cell size I 33)
//  2. the component directions are grouped into I 32 direction bins; for each bin a sweep along the
//     propagation direction marks the cells that are reachable from the upwave boundary without
//     crossing land (geometric shadow)
//  3. the masks are smoothed over the wet cells with a normalized Gaussian-like filter of standard
//     deviation I 31 (diffraction into the shadow edges); land does not pull the masks down along
//     exposed coastlines
//  4. eta, Fifsf and Fi are summed over the components (B 89 decomposition: space part x time part),
//     each component weighted with the mask of its direction (linear between the two nearest bins)
//
// Wave theories without the component decomposition (B 92 other than 5, 31, 41, 51) use one mask
// for the main direction (B 105).
// The rest of the initialization (beach taper, coastline relaxation with A 349, velocities) is the
// same as for I 30 1.

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<vector>
#include<cmath>
#include<algorithm>
#include<fstream>
#include<iomanip>

namespace
{

struct mask_grid
{
    int nx=1, ny=1;
    double xs=0.0, ys=0.0, dx=1.0, dy=1.0;

    int id(int i, int j) const {return i*ny + j;}

    // bilinear interpolation of a cell centred field, clamped to the grid
    double interp(const std::vector<float> &f, double x, double y) const
    {
        double fi = (x-xs)/dx - 0.5;
        double fj = (y-ys)/dy - 0.5;

        fi = std::min(std::max(fi,0.0),double(nx-1));
        fj = std::min(std::max(fj,0.0),double(ny-1));

        const int i0 = std::min(int(fi),nx-2>0?nx-2:0);
        const int j0 = std::min(int(fj),ny-2>0?ny-2:0);
        const int i1 = std::min(i0+1,nx-1);
        const int j1 = std::min(j0+1,ny-1);
        const double wx = fi-i0;
        const double wy = fj-j0;

        return (1.0-wx)*(1.0-wy)*f[id(i0,j0)] + wx*(1.0-wy)*f[id(i1,j0)]
             + (1.0-wx)*wy*f[id(i0,j1)] + wx*wy*f[id(i1,j1)];
    }
};

// geometric exposure for the direction th (rad, global frame): sweep along the dominant axis of the
// direction; a cell is exposed when it is wet and the upwave point of the ray on the previous line
// (linear interpolation) is exposed; upwave points outside the domain are exposed (open boundary)
void exposure_sweep(const mask_grid &g, const std::vector<float> &W, std::vector<float> &L, double th)
{
    const double c = cos(th);
    const double s = sin(th);

    if(fabs(c)>=fabs(s))
    {
        const int dir = c>0.0 ? 1 : -1;
        const double shift = -(g.dx/g.dy)*s/fabs(c);   // upwave offset in j per line
        const int i0 = dir>0 ? 0 : g.nx-1;

        for(int j=0; j<g.ny; ++j)
        L[g.id(i0,j)] = W[g.id(i0,j)];

        for(int n=1; n<g.nx; ++n)
        {
            const int i = i0 + dir*n;
            const int ip = i - dir;

            for(int j=0; j<g.ny; ++j)
            {
                const double js = j + shift;
                double up = 1.0;

                if(js>=0.0 && js<=double(g.ny-1))
                {
                    const int ja = std::min(int(js),g.ny-1);
                    const int jb = std::min(ja+1,g.ny-1);
                    const double w = js-ja;
                    up = (1.0-w)*L[g.id(ip,ja)] + w*L[g.id(ip,jb)];
                }

                L[g.id(i,j)] = W[g.id(i,j)]*up;
            }
        }
    }

    else
    {
        const int dir = s>0.0 ? 1 : -1;
        const double shift = -(g.dy/g.dx)*c/fabs(s);   // upwave offset in i per line
        const int j0 = dir>0 ? 0 : g.ny-1;

        for(int i=0; i<g.nx; ++i)
        L[g.id(i,j0)] = W[g.id(i,j0)];

        for(int n=1; n<g.ny; ++n)
        {
            const int j = j0 + dir*n;
            const int jp = j - dir;

            for(int i=0; i<g.nx; ++i)
            {
                const double is = i + shift;
                double up = 1.0;

                if(is>=0.0 && is<=double(g.nx-1))
                {
                    const int ia = std::min(int(is),g.nx-1);
                    const int ib = std::min(ia+1,g.nx-1);
                    const double w = is-ia;
                    up = (1.0-w)*L[g.id(ia,jp)] + w*L[g.id(ib,jp)];
                }

                L[g.id(i,j)] = W[g.id(i,j)]*up;
            }
        }
    }
}

// one box filter pass of half width h along x (ax=1) or y (ax=0), via running sums
void box_pass(const mask_grid &g, std::vector<double> &f, int h, int ax, std::vector<double> &tmp)
{
    if(h<=0)
    return;

    const int n1 = ax ? g.ny : g.nx;     // lines
    const int n2 = ax ? g.nx : g.ny;     // points per line

    std::vector<double> cs(n2+1);

    for(int a=0; a<n1; ++a)
    {
        cs[0]=0.0;

        for(int b=0; b<n2; ++b)
        cs[b+1] = cs[b] + (ax ? f[g.id(b,a)] : f[g.id(a,b)]);

        for(int b=0; b<n2; ++b)
        {
            const int lo = std::max(b-h,0);
            const int hi = std::min(b+h,n2-1);
            const double v = (cs[hi+1]-cs[lo])/double(2*h+1);   // zero padding outside the domain

            if(ax)
            tmp[g.id(b,a)] = v;
            else
            tmp[g.id(a,b)] = v;
        }
    }

    f.swap(tmp);
}

// normalized smoothing over the wet cells: S(L W)/S(W), three box passes in each direction
// (standard deviation sqrt(h(h+1)) cells per direction); defined wherever wet cells are in reach
void smooth_masked(const mask_grid &g, const std::vector<float> &W, std::vector<float> &L, int hx, int hy)
{
    if(hx<=0 && hy<=0)
    return;

    const size_t N = size_t(g.nx)*g.ny;
    std::vector<double> A(N), B(N), tmp(N);

    for(size_t q=0; q<N; ++q)
    {
        A[q] = double(L[q])*W[q];
        B[q] = W[q];
    }

    for(int pass=0; pass<3; ++pass)
    {
        box_pass(g,A,hx,1,tmp);
        box_pass(g,A,hy,0,tmp);
        box_pass(g,B,hx,1,tmp);
        box_pass(g,B,hy,0,tmp);
    }

    for(size_t q=0; q<N; ++q)
    L[q] = B[q]>1.0e-6 ? float(std::min(std::max(A[q]/B[q],0.0),1.0)) : 0.0f;
}

}

void iowave::full_initialize_fnpf_masked(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    const double t0 = pgc->timer();

    if(p->mpirank==0)
    cout<<"full NWT initialize with exposure masks (I 30 2)"<<endl;

    p->wavetime = p->simtime;

    select_sources(nullptr);

    // ---------------------------------------------------------------- coarse grid
    double xmin=1.0e20, xmax=-1.0e20, ymin=1.0e20, ymax=-1.0e20;
    int nloc=0;

    SLICELOOP4
    {
        xmin = MIN(xmin,p->pos_x());
        xmax = MAX(xmax,p->pos_x());
        ymin = MIN(ymin,p->pos_y());
        ymax = MAX(ymax,p->pos_y());
        ++nloc;
    }

    xmin = pgc->globalmin(xmin);
    ymin = pgc->globalmin(ymin);
    xmax = pgc->globalmax(xmax);
    ymax = pgc->globalmax(ymax);

    // mean horizontal cell size from the cell count (DXM also contains the vertical spacing)
    const int ncell = pgc->globalisum(nloc);
    const double ext_x = MAX(xmax-xmin,1.0e-6);
    const double ext_y = MAX(ymax-ymin,1.0e-6);
    const double dxm = p->j_dir==0 ? ext_x/double(MAX(ncell-1,1)) : sqrt(ext_x*ext_y/double(MAX(ncell,1)));
    const double Lx = MAX(ext_x,dxm);
    const double Ly = MAX(ext_y,dxm);

    double dc = p->I33;

    if(dc<=0.0)
    dc = MAX(2.0*dxm, MAX(Lx,Ly)/1500.0);

    dc = MAX(dc,dxm);

    mask_grid g;
    g.nx = MAX(1,int(Lx/dc+0.5));
    g.ny = MAX(1,int(Ly/dc+0.5));
    g.dx = (Lx+dxm)/double(g.nx);
    g.dy = (Ly+dxm)/double(g.ny);
    g.xs = xmin - 0.5*dxm;
    g.ys = ymin - 0.5*dxm;

    if(p->j_dir==0)
    {
        g.ny = 1;
        g.dy = Ly+dxm;
    }

    const size_t N = size_t(g.nx)*g.ny;

    // ---------------------------------------------------------------- global wet mask
    std::vector<float> wsum(N,0.0f), csum(N,0.0f);

    SLICELOOP4
    {
        const int ci = std::min(std::max(int((p->pos_x()-g.xs)/g.dx),0),g.nx-1);
        const int cj = std::min(std::max(int((p->pos_y()-g.ys)/g.dy),0),g.ny-1);

        csum[g.id(ci,cj)] += 1.0f;

        if(p->wet[IJ]==1)
        wsum[g.id(ci,cj)] += 1.0f;
    }

    MPI_Allreduce(MPI_IN_PLACE,wsum.data(),int(N),MPI_FLOAT,MPI_SUM,pgc->mpi_comm);
    MPI_Allreduce(MPI_IN_PLACE,csum.data(),int(N),MPI_FLOAT,MPI_SUM,pgc->mpi_comm);

    std::vector<float> W(N,0.0f);

    for(size_t q=0; q<N; ++q)
    W[q] = (csum[q]>0.0f && wsum[q]>=0.5f*csum[q]) ? 1.0f : 0.0f;

    wsum.clear(); wsum.shrink_to_fit();
    csum.clear(); csum.shrink_to_fit();

    // ---------------------------------------------------------------- components and direction bins
    const int nc = decomp_ncomp(p);
    const bool decomp = nc>0;

    std::vector<double> th;     // global direction of each component (rad)

    if(decomp)
    {
        th.resize(nc);

        for(int n=0; n<nc; ++n)
        th[n] = gamma + decomp_dir(p,n);
    }

    double thmin=gamma, thmax=gamma;

    if(decomp)
    {
        thmin = *std::min_element(th.begin(),th.end());
        thmax = *std::max_element(th.begin(),th.end());
    }

    int Nb = MAX(1,p->I32);

    if(thmax-thmin<1.0e-6)
    Nb=1;

    std::vector<double> thb(Nb);

    for(int b=0; b<Nb; ++b)
    thb[b] = Nb==1 ? 0.5*(thmin+thmax) : thmin + b*(thmax-thmin)/double(Nb-1);

    // smoothing length: I 31, automatic 2 peak wave lengths (at least 2 coarse cells)
    double sig = p->I31;

    if(sig<0.0)
    {
        // peak wave length from the linear dispersion relation at the wave theory depth (B 94)
        const double T = p->wTp>0.0 ? p->wTp : p->wT;
        const double d = MAX(p->B94_wdt>0.0 ? p->B94_wdt : p->wd, 0.1);
        double Lp = 10.0*dxm;
        
        if(T>0.0)
        {
            const double w = 2.0*PI/T;
            double k = w*w/9.81;
            
            for(int it=0; it<50; ++it)
            k = w*w/(9.81*tanh(MAX(k*d,1.0e-10)));
            
            Lp = 2.0*PI/k;
        }
        
        sig = 2.0*Lp;
    }

    const int hx = sig>0.0 ? MAX(1,int(sig/g.dx+0.5)) : 0;
    const int hy = (sig>0.0 && g.ny>1) ? MAX(1,int(sig/g.dy+0.5)) : 0;

    // ---------------------------------------------------------------- masks, sampled at the local cells
    std::vector<float> L(N), Lmean(N,0.0f);
    std::vector<float> mloc(size_t(Nb)*MAX(nloc,1));
    std::vector<double> ebin(Nb,0.0);     // energy per bin (for the mean mask output)

    if(decomp)
    for(int n=0; n<nc; ++n)
    {
        const double a = wave_eta_space_cos(p,pgc,0.0,0.0,n);
        const double f = Nb==1 ? 0.0 : (th[n]-thmin)/(thmax-thmin)*(Nb-1);
        const int b0 = std::min(int(f),Nb-1);
        const int b1 = std::min(b0+1,Nb-1);
        const double w = f-b0;
        ebin[b0] += (1.0-w)*a*a;
        ebin[b1] += w*a*a;
    }

    double esum=0.0;
    for(int b=0; b<Nb; ++b)
    esum += ebin[b];

    double exposed_mean=0.0;

    for(int b=0; b<Nb; ++b)
    {
        exposure_sweep(g,W,L,thb[b]);
        smooth_masked(g,W,L,hx,hy);

        const double eb = esum>0.0 ? ebin[b]/esum : 1.0/double(Nb);

        for(size_t q=0; q<N; ++q)
        Lmean[q] += float(eb)*L[q];

        int q=0;
        SLICELOOP4
        {
            mloc[size_t(b)*nloc + q] = float(g.interp(L,p->pos_x(),p->pos_y()));
            ++q;
        }
    }

    // fraction of the wet area that is exposed (energy weighted)
    {
        double wa=0.0, ea=0.0;
        for(size_t q=0; q<N; ++q)
        if(W[q]>0.5f)
        {
            wa += 1.0;
            ea += Lmean[q];
        }
        exposed_mean = wa>0.0 ? ea/wa : 1.0;
    }

    // mean mask for checking: ESRI ASCII grid, rows from north (max y) to south
    if(p->mpirank==0)
    {
        std::ofstream out("REEF3D_FNPF_I30_mask.asc");
        out<<"ncols "<<g.nx<<"\nnrows "<<g.ny<<"\nxllcorner "<<setprecision(10)<<g.xs<<"\nyllcorner "<<g.ys
           <<"\ncellsize "<<g.dx<<"\nNODATA_value -1\n";
        out<<setprecision(4);
        for(int j=g.ny-1; j>=0; --j)
        {
            for(int i=0; i<g.nx; ++i)
            out<<(W[g.id(i,j)]>0.5f ? double(Lmean[g.id(i,j)]) : -1.0)<<(i<g.nx-1?" ":"");
            out<<"\n";
        }
    }

    L.clear(); L.shrink_to_fit();
    Lmean.clear(); Lmean.shrink_to_fit();
    W.clear(); W.shrink_to_fit();

    // component weights per local cell: linear between the two nearest direction bins
    std::vector<int> cb0, cb1;
    std::vector<double> cw;

    if(decomp)
    {
        cb0.resize(nc); cb1.resize(nc); cw.resize(nc);

        for(int n=0; n<nc; ++n)
        {
            const double f = Nb==1 ? 0.0 : (th[n]-thmin)/(thmax-thmin)*(Nb-1);
            cb0[n] = std::min(int(f),Nb-1);
            cb1[n] = std::min(cb0[n]+1,Nb-1);
            cw[n] = f-cb0[n];
        }
    }

    auto cmask = [&](int n, int q) -> double
    {
        if(!decomp)
        return mloc[q];

        return (1.0-cw[n])*mloc[size_t(cb0[n])*nloc + q] + cw[n]*mloc[size_t(cb1[n])*nloc + q];
    };

    // time parts
    std::vector<double> eTc, eTs, fTc, fTs;

    if(decomp)
    {
        eTc.resize(nc); eTs.resize(nc); fTc.resize(nc); fTs.resize(nc);

        for(int n=0; n<nc; ++n)
        {
            eTc[n] = wave_eta_time_cos(p,pgc,n);
            eTs[n] = wave_eta_time_sin(p,pgc,n);
            fTc[n] = wave_fi_time_cos(p,pgc,n);
            fTs[n] = wave_fi_time_sin(p,pgc,n);
        }
    }

    // local cell index of (i,j) in SLICELOOP4 order
    std::vector<int> qidx(size_t(p->imax)*p->jmax,-1);
    {
        int q=0;
        SLICELOOP4
        {
            qidx[IJ] = q;
            ++q;
        }
    }

    // ---------------------------------------------------------------- eta
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        xg = xgen(p);
        yg = ygen(p);
        const int q = qidx[IJ];

        double val=0.0;

        if(decomp)
        {
            for(int n=0; n<nc; ++n)
            {
                const double m = cmask(n,q);

                if(m>1.0e-4)
                val += m*(wave_eta_space_cos(p,pgc,xg,yg,n)*eTc[n] - wave_eta_space_sin(p,pgc,xg,yg,n)*eTs[n]);
            }
        }

        else
        val = cmask(0,q)*wave_eta(p,pgc,xg,yg);

        c->eta(i,j) = val;
    }

    for(int qn=0; qn<10;++qn)
    SLICELOOP4
    {
		db = distbeach(p);

		if(p->B99==1 || p->B99==2)
		{
            if(p->A10!=3 || p->A348==1 || p->A348==2)
            if(db<1.0e20)
            {
            if(p->wet[IJ]==1)
            c->eta(i,j) = relax4_nb(i,j)*c->eta(i,j);
            }
        }
    }

    // ---------------------------------------------------------------- Fifsf
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        xg = xgen(p);
        yg = ygen(p);
        z = c->eta(i,j);
        const int q = qidx[IJ];

        double val=0.0;

        if(decomp)
        {
            for(int n=0; n<nc; ++n)
            {
                const double m = cmask(n,q);

                if(m>1.0e-4)
                val += m*(wave_fi_space_cos(p,pgc,xg,yg,z,n)*fTc[n] - wave_fi_space_sin(p,pgc,xg,yg,z,n)*fTs[n]);
            }
        }

        else
        val = cmask(0,q)*wave_fi(p,pgc,xg,yg,z);

        c->Fifsf(i,j) = val;
    }

    for(int qn=0; qn<10;++qn)
    SLICELOOP4
    {
		db = distbeach(p);

		if(p->B99==1 || p->B99==2)
		{
            if(p->A10!=3 || p->A348==1 || p->A348==2)
            if(db<1.0e20)
            {
            if(p->wet[IJ]==1)
            c->Fifsf(i,j) = relax4_nb(i,j)*c->Fifsf(i,j);
            }
        }
    }

    // ---------------------------------------------------------------- Fi
    FLOOP
    if(p->wet[IJ]==1)
    {
        xg = xgen(p);
        yg = ygen(p);
        z = p->ZN[KP]*(c->eta(i,j) + p->wd - c->bed(i,j)) + c->bed(i,j)-p->phimean;
        const int q = qidx[IJ];

        double val=0.0;

        if(decomp)
        {
            for(int n=0; n<nc; ++n)
            {
                const double m = cmask(n,q);

                if(m>1.0e-4)
                val += m*(wave_fi_space_cos(p,pgc,xg,yg,z,n)*fTc[n] - wave_fi_space_sin(p,pgc,xg,yg,z,n)*fTs[n]);
            }
        }

        else
        val = cmask(0,q)*wave_fi(p,pgc,xg,yg,z);

        c->Fi[FIJK] = val;
    }

    for(int qn=0; qn<10;++qn)
    LOOP
    {
        db = distbeach(p);

        FKLOOP
        {
        if(p->A10!=3 || p->A348==1 || p->A348==3)
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            if(db<1.0e20)
            c->Fi[FIJK] = relax4_nb(i,j)*c->Fi[FIJK];
        }
        }
    }

    SLICELOOP4
    c->WL(i,j) = c->eta(i,j) + p->wd - c->bed(i,j);

    const double t1 = pgc->timer();

    if(p->mpirank==0)
    {
        cout<<" exposure masks: "<<g.nx<<" x "<<g.ny<<" cells of "<<g.dx<<" m, "<<Nb<<" direction bins ("
            <<thmin*180.0/PI<<" .. "<<thmax*180.0/PI<<" deg), smoothing "<<sig<<" m"<<endl;
        cout<<" components: "<<(decomp ? nc : 0)<<(decomp ? "" : " (no decomposition, one mask for B 105)")
            <<" . exposed wet area (energy weighted): "<<exposed_mean*100.0<<" %"
            <<" . mask file: REEF3D_FNPF_I30_mask.asc . time: "<<t1-t0<<" s"<<endl;
    }
}
