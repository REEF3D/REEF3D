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

//  REEFAMR: the load of the ranks.
//
//  load_step: once per step from the module; the load statistics (G 43, amr_load.h) with the work
//  predicted by the cost model of the placement.
//
//  Measured cost model (G 41 0, on several ranks; for the placement of G 40 1 and the load map):
//  the work of a rank per step is modelled as  a (x1 + wp x2 + ws x3)  with
//    x1  the level-0 cells that are not solid (and the level-0 blocks of the rank),
//    x2  the cells of its patches with their EXT cells and layers (as in place_patches),
//    x3  the parent blocks it serves to patches on other ranks (the block plans), in cells
//        of the parent level,
//  and busy, the wall-clock time per step without the MPI calls and the regrids (amr_clock), as the
//  measured work.  At every regrid the window since the last one gives one sample per rank (the
//  first window, with the initialisation, is left out); a, a wp and a ws are fitted by least
//  squares over the ranks and the last LWIN windows, with a ridge of 10 % of the diagonal towards
//  the last weights (directions the data do not resolve keep them).  wp in [0.1,16], ws in [0,16],
//  averaged with the last weights.  Without the MPI clock (-DREEF3D_NO_MPI_CLOCK) the model stays
//  fixed (wp = 1, ws = 0, x1 the solid-free cells).
//  A static refinement (no regrid) is placed once more at step 15 with the weights measured over
//  steps 5-15 (replace_due).
//
//  Load map (G 45 N): the work per step of every level-0 column, averaged over the steps of the run,
//    its level-0 cells that are not solid + wp (G 41, or measured) times the patch cells over it
//  (all levels, their EXT cells spread over the patch, the layers), by location (the rank that holds
//  a patch does not matter).  Every N steps rank 0 writes REEF3D_<M>_AMR/REEF3D_<M>_AMR_loadmap.dat
//  (GNX GNY, then a line of GNY values per i) for the weighted partition of DIVEMesh (M 20 5, M 35 1).
//  Not for multi-block grids.

#include"reefamr.h"
#include"lexer.h"
#include"ghostcell.h"
#include<mpi.h>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<sys/stat.h>

namespace
{
const size_t LWIN = 8;      // windows in the fit
}

void reefamr::load_counts(lexer *p, double &x1, double &x2, double &x3)
{
    const double kz = double(MAX(p->knoz,1));

    long act=0;
    for(int ii=0; ii<NX0; ++ii)
    for(int jj=0; jj<NY0; ++jj)
    if(flag0(ii+O0i,jj+O0j)>0)
    ++act;
    x1 = double(act)*kz;
    x2 = 0.0;
    x3 = 0.0;

    for(auto q : P)
    {
        if(q->lev==0)
        {
            x1 += double(q->nx)*double(q->ny)*double(MAX(q->pp->knoz,1));
            continue;
        }
        x2 += double(q->nx+2*EXT)*double(q->ny+2*EXT)*kz*double(vfac.empty() ? 1 : vfac[q->lev]);
    }

    for(int l=1; l<=maxlev && l<(int)bsrv.size(); ++l)
    x3 += double(bsrv[l].size())*kz*double(vfac.empty() ? 1 : vfac[l-1]);
}

// x1 of all ranks from the global solid flags (place 1) and the level-0 blocks
void reefamr::load_base(lexer *p)
{
    const int np = p->mpi_size;
    const double kz = double(MAX(p->knoz,1));
    lbase.assign(np,0.0);

    for(int r=0; r<np && r<(int)rbx0.size(); ++r)
    {
        long act=0;
        for(int I=rbx0[r]; I<=rbx1[r]; ++I)
        for(int J=rby0[r]; J<=rby1[r]; ++J)
        if(flag0(I,J)>0)
        ++act;
        lbase[r] = double(act)*kz;
    }

    for(int b=1; b<nblk; ++b)
    if(blk[b].rank>=0 && blk[b].rank<np)
    lbase[blk[b].rank] += double(blk[b].ni)*double(blk[b].nj)*kz;
}

void reefamr::load_window(lexer *p, bool use)
{
    if(!lmeas)
    return;

    const int np = p->mpi_size;

    if(use && lw_open && !lw_first && p->count>lw_c0)
    {
        const double now = MPI_Wtime();
        const double steps = double(p->count-lw_c0);
        const double busy = (now-lw_t0) - (amr_clock::mpi-lw_m0) - ((amr_clock::regrid-lw_r0)-(amr_clock::regrid_mpi-lw_rm0));

        double d[4];
        d[0] = busy/steps;
        load_counts(p,d[1],d[2],d[3]);

        vector<double> all((size_t)4*np);
        MPI_Allgather(d,4,MPI_DOUBLE,&all[0],4,MPI_DOUBLE,MPI_COMM_WORLD);

        // without the MPI clock the busy time is the wall-clock time: no measured model
        long calls = amr_clock::mpi_calls;
        long callsg = 0;
        MPI_Allreduce(&calls,&callsg,1,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
        if(callsg==0)
        {
            if(p->mpirank==0)
            cout<<par.name<<": G 41 0 needs the MPI clock (built with REEF3D_NO_MPI_CLOCK) -- the cost model of the placement stays fixed"<<endl;
            lmeas = false;
            lwp = 1.0;
            lws = 0.0;
            return;
        }

        vector<lsample> S(np);
        double sy=0.0, my=0.0;
        for(int r=0; r<np; ++r)
        {
            S[r] = {all[4*r],all[4*r+1],all[4*r+2],all[4*r+3]};
            sy += S[r].y;
            my = MAX(my,S[r].y);
        }
        lmeasL = (sy>0.0) ? my/(sy/np) : 1.0;

        lwin.push_back(S);
        if(lwin.size()>LWIN)
        lwin.pop_front();

        load_fit();
    }

    // a new window after the first one (the initialisation) leaves that one out
    if(lw_open)
    lw_first = false;

    lw_open = true;
    lw_t0 = MPI_Wtime();
    lw_m0 = amr_clock::mpi;
    lw_r0 = amr_clock::regrid;
    lw_rm0 = amr_clock::regrid_mpi;
    lw_c0 = p->count;
}

void reefamr::load_fit()
{
    double A[3][3] = {{0,0,0},{0,0,0},{0,0,0}}, b[3] = {0,0,0};
    double sy=0.0, sx=0.0;
    int n=0;

    for(const auto &S : lwin)
    for(const auto &s : S)
    {
        const double x[3] = {s.x1,s.x2,s.x3};
        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            A[i][j] += x[i]*x[j];
            b[i] += x[i]*s.y;
        }
        sy += s.y;
        sx += s.x1 + lwp*s.x2 + lws*s.x3;
        ++n;
    }

    if(n<2 || sy<=0.0 || sx<=0.0)
    return;

    // ridge towards the last weights
    const double a0 = sy/sx;
    const double t0[3] = {a0, a0*lwp, a0*lws};
    const double tr = A[0][0]+A[1][1]+A[2][2];
    for(int k=0; k<3; ++k)
    {
        const double lam = 0.1*A[k][k] + 1.0e-9*tr;
        A[k][k] += lam;
        b[k] += lam*t0[k];
    }

    // Gauss with partial pivoting
    double t[3];
    {
        double M[3][4];
        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            M[i][j] = A[i][j];
            M[i][3] = b[i];
        }
        for(int c=0; c<3; ++c)
        {
            int piv=c;
            for(int r=c+1; r<3; ++r)
            if(fabs(M[r][c])>fabs(M[piv][c]))
            piv=r;
            if(fabs(M[piv][c])<=1.0e-300)
            return;
            if(piv!=c)
            for(int j=0; j<4; ++j)
            std::swap(M[c][j],M[piv][j]);
            for(int r=c+1; r<3; ++r)
            {
                const double f = M[r][c]/M[c][c];
                for(int j=c; j<4; ++j)
                M[r][j] -= f*M[c][j];
            }
        }
        for(int i=2; i>=0; --i)
        {
            double v = M[i][3];
            for(int j=i+1; j<3; ++j)
            v -= M[i][j]*t[j];
            t[i] = v/M[i][i];
        }
    }

    if(!(t[0]>0.0))
    return;

    const double wp = MIN(MAX(t[1]/t[0],0.1),16.0);
    const double ws = MIN(MAX(t[2]/t[0],0.0),16.0);
    lwp = 0.5*(lwp+wp);
    lws = 0.5*(lws+ws);
    ++lfits;

    if(p0->mpirank==0)
    {
        const streamsize pr = cout.precision();
        cout<<par.name<<": measured cost model (G 41 0, "<<lwin.size()<<" windows): a patch cell "<<setprecision(3)<<lwp
            <<", a served parent cell "<<lws<<" level-0 cells; busy max/mean of the last window "<<lmeasL<<endl;
        cout.precision(pr);
    }
}

bool reefamr::replace_due(lexer *p)
{
    if(lmeas && par.place==1 && regrid_int==0 && maxlev>0 && p->count==15)
    {
        force_place = true;
        return true;
    }
    return false;
}

// patch cells per level-0 column of the rank box, current layout (the global patch table)
void reefamr::lmap_layout(lexer *p)
{
    const double kz = double(MAX(p->knoz,1));
    lmap_cur.assign((size_t)NX0*NY0,0.0);

    for(const reefamr_gpatch &G : GP)
    {
        const int l = G.lev;
        if(l<1)
        continue;
        const double nx = G.I1-G.I0+1, ny = G.J1-G.J0+1;
        const double f = (nx+2*EXT)*(ny+2*EXT)/(nx*ny)*kz*double(vfac.empty() ? 1 : vfac[l]);

        const int a0 = MAX(G.I0>>l,O0i), a1 = MIN(G.I1>>l,O0i+NX0-1);
        const int b0 = MAX(G.J0>>l,O0j), b1 = MIN(G.J1>>l,O0j+NY0-1);
        for(int I=a0; I<=a1; ++I)
        {
            const int ci = MIN(G.I1,((I+1)<<l)-1) - MAX(G.I0,I<<l) + 1;
            if(ci<=0)
            continue;
            for(int J=b0; J<=b1; ++J)
            {
                const int cj = MIN(G.J1,((J+1)<<l)-1) - MAX(G.J0,J<<l) + 1;
                if(cj>0)
                lmap_cur[(size_t)(I-O0i)*NY0+(J-O0j)] += f*double(ci)*double(cj);
            }
        }
    }
    lmap_dirty = false;
}

void reefamr::lmap_step(lexer *p)
{
    if(p->G45<=0 || nblk>0 || maxlev<1)
    return;

    if(lmap_dirty || lmap_cur.size()!=(size_t)NX0*NY0)
    lmap_layout(p);

    if(lmap_acc.size()!=lmap_cur.size())
    lmap_acc.assign(lmap_cur.size(),0.0);

    for(size_t n=0; n<lmap_cur.size(); ++n)
    lmap_acc[n] += lmap_cur[n];
    ++lmap_steps;

    if(p->count%p->G45==0)
    lmap_write(p);
}

void reefamr::lmap_write(lexer *p)
{
    const int np = p->mpi_size;
    const double kz = double(MAX(p->knoz,1));
    const double wp = lmeas ? lwp : ((p->G41>0.0) ? p->G41 : 1.0);
    const double st = double(MAX(lmap_steps,1L));

    // the work per step of the columns of the rank box
    vector<double> mine((size_t)NX0*NY0);
    for(int ii=0; ii<NX0; ++ii)
    for(int jj=0; jj<NY0; ++jj)
    {
        const size_t n = (size_t)ii*NY0+jj;
        mine[n] = ((flag0(ii+O0i,jj+O0j)>0) ? kz : 0.0) + wp*lmap_acc[n]/st;
    }

    vector<int> cnt(np), off(np,0);
    for(int r=0; r<np; ++r)
    cnt[r] = (rbx1[r]-rbx0[r]+1)*(rby1[r]-rby0[r]+1);
    for(int r=1; r<np; ++r)
    off[r] = off[r-1]+cnt[r-1];
    vector<double> all(p->mpirank==0 ? (size_t)(off[np-1]+cnt[np-1]) : 1);
    MPI_Gatherv(&mine[0],(int)mine.size(),MPI_DOUBLE,&all[0],&cnt[0],&off[0],MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(p->mpirank!=0)
    return;

    vector<double> map((size_t)GNX*GNY,0.0);
    for(int r=0; r<np; ++r)
    {
        const int ny = rby1[r]-rby0[r]+1;
        for(int I=rbx0[r]; I<=rbx1[r]; ++I)
        for(int J=rby0[r]; J<=rby1[r]; ++J)
        map[(size_t)I*GNY+J] = all[off[r]+(size_t)(I-rbx0[r])*ny+(J-rby0[r])];
    }

    string tag = par.name;
    for(auto &ch : tag)
    if(ch==' ')
    ch='_';
    const string dir = "./REEF3D_"+tag;
    mkdir(dir.c_str(),0777);
    ofstream out((dir+"/REEF3D_"+tag+"_loadmap.dat").c_str());
    out<<"# "<<par.name<<" load map (G 45): work per step of every level-0 column in level-0 cells, averaged over "<<lmap_steps
       <<" steps: the cells that are not solid + "<<wp<<" (a patch cell"<<(lmeas ? ", measured" : "")<<") times the patch cells over it"<<endl;
    out<<"# for the weighted partition of DIVEMesh: M 20 5, M 35 1 (copy as loadmap.dat or run DIVEMesh in this case folder)"<<endl;
    out<<GNX<<" "<<GNY<<endl;
    out<<setprecision(6);
    for(int I=0; I<GNX; ++I)
    {
        for(int J=0; J<GNY; ++J)
        out<<map[(size_t)I*GNY+J]<<(J+1<GNY ? " " : "");
        out<<"\n";
    }
}

void reefamr::load_step(lexer *p)
{
    // load map (G 45)
    lmap_step(p);

    // static refinement with the measured model: the window of the re-placement opens at step 5
    if(lmeas && par.place==1 && regrid_int==0 && maxlev>0 && p->count==5)
    load_window(p,false);

    if(p->G43<=0)
    return;

    if(!lstat.ready())
    lstat.configure(par.name);

    lstat.step(p,[&]()
    {
        amr_load_counts c;
        const double kz = double(MAX(p->knoz,1));
        double x1, x2, x3;
        load_counts(p,x1,x2,x3);

        c.n0 = double(NX0)*double(NY0)*kz;
        for(auto q : P)
        if(q->lev==0)
        c.n0 += double(q->nx)*double(q->ny)*double(MAX(q->pp->knoz,1));
        else
        ++c.npatch;
        c.n0a = x1;
        c.np = x2;

        // the work predicted by the cost model of the placement
        if(lmeas)
        c.wpred = x1 + lwp*x2 + lws*x3;
        else
        c.wpred = double(NX0)*double(NY0)*kz + ((p->G41>0.0) ? p->G41 : 1.0)*x2;
        return c;
    });
}
