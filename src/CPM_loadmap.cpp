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

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"amr_load.h"
#include<mpi.h>
#include<iomanip>
#include<sys/stat.h>

/*--------------------------------------------------------------------
Load map of CFD with CPM (G 45 N), for the weighted partition (DIVEMesh or REEF3D M 40: M 20 5, M 35 1;
REEF3D M 41 1), a map of the cells in the format of the load map of the CFD mesh refinement
(cfd_amr::write_loadmap):

  the work per step of every cell of the grid, in fluid cells:
      w(i,j,k) = f(i,j,k) + (u_s S(i,j,k) + u_a A(i,j,k))/(u_f N)
  f         1 for a fluid cell (flag4 > 0), else 0
  S, A      parcel updates in the cell over the N steps of the run: all parcels in every sub-step (S),
            and of those the ones that are awake (A), i.e. not sleeping grains of the jammed bed
  u_f       measured cost of a fluid cell per step: the busy time of the ranks outside CPM over their fluid cells
  u_s, u_a  measured cost of a parcel update and the extra cost of an awake one: least squares over the ranks,
            busy time of CPM = u_s S_r + u_a A_r (busy = wall - time in blocking MPI calls, mpi_clock.cpp);
            with one rank, or when the ranks do not separate the two, u_a = 2 u_s (CPM performance,
            validation 22: an awake update costs about three sleeping ones)

Every N steps rank 0 writes ./REEF3D_CFD_CPM_Particle/REEF3D_CFD_CPM_loadmap.dat: GNX GNY GNZ, then the
cells (i GNY + j) GNZ + k, ten per line. With G 46 (rebalance) the counts are kept without G 45 as well, and the
map is written when the run stops to rebalance (reefrun_ckpt, the hook loadmap of the driver), from the steps
since the start of the run; the restart decomposes with it (M 41 1). The weighted partition gives the active cells of DIVEMesh this
work (the parcels in cells that DIVEMesh counts as solid are not counted). Not for multi-block grids.
--------------------------------------------------------------------*/

void CPM::loadmap_substep(lexer *p)
{
    if(!loadmap_on(p))
    return;

    const size_t nc = size_t(p->knox)*p->knoy*p->knoz;

    if(lmS.size()!=nc)
    {
        lmS.assign(nc,0.0);
        lmA.assign(nc,0.0);
    }

    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        const int ic = MAX(0, MIN(p->knox-1, pci(p,n,P.X[n])));
        const int jc = p->j_dir==1 ? MAX(0, MIN(p->knoy-1, pcj(p,n,P.Y[n]))) : 0;
        const int kc = MAX(0, MIN(p->knoz-1, pck(p,n,P.Z[n])));
        const size_t c = (size_t(ic)*p->knoy + jc)*p->knoz + kc;

        lmS[c] += 1.0;

        if(n<int(lmAw.size()) && lmAw[n])
        lmA[c] += 1.0;
    }
}

void CPM::loadmap_step(lexer *p, fdm *a, ghostcell *pgc, double tcpm, double mcpm)
{
    if(!loadmap_on(p))
    return;

    const double now = MPI_Wtime();
    const double mpi = amr_clock::mpi;

    // the busy time of the step before this call (the first call opens the window)
    if(lm_open)
    {
        lmTstep += (now - lm_w0) - (mpi - lm_m0);
        lmTcpm += tcpm - mcpm;
        ++lm_steps;
    }

    lm_open = 1;
    lm_w0 = now;
    lm_m0 = mpi;

    if(p->G45>0 && lm_steps>0 && p->count%p->G45==0)
    loadmap_write(p,a,pgc);
}

void CPM::loadmap_write(lexer *p, fdm *a, ghostcell *pgc)
{
    const int GNX = p->gknox, GNY = p->j_dir==1 ? p->gknoy : 1, GNZ = p->gknoz;
    const size_t nc = size_t(p->knox)*p->knoy*p->knoz;

    if(lmS.size()!=nc)
    {
        lmS.assign(nc,0.0);
        lmA.assign(nc,0.0);
    }

    // the fluid cells
    std::vector<double> nf(nc,0.0);
    double nfr=0.0, Sr=0.0, Ar=0.0;

    LOOP
    if(p->flag4[IJK]>0)
    nf[(size_t(i)*p->knoy + j)*p->knoz + k] = 1.0;

    for(size_t c=0;c<nc;++c)
    {
        nfr += nf[c];
        Sr += lmS[c];
        Ar += lmA[c];
    }

    const double st = double(MAX(lm_steps,1L));
    const double tfl = MAX(0.0, lmTstep - lmTcpm);

    // the ranks: CPM busy time, parcel updates, fluid busy time and cells
    int np=1;
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    double mine[5] = {lmTcpm, Sr, Ar, tfl, nfr*st};
    std::vector<double> all(size_t(5)*np);
    MPI_Allgather(mine,5,MPI_DOUBLE,&all[0],5,MPI_DOUBLE,MPI_COMM_WORLD);

    double sC=0.0, sS=0.0, sA=0.0, sF=0.0, sN=0.0;
    double sSS=0.0, sSA=0.0, sAA=0.0, sCS=0.0, sCA=0.0;

    for(int r=0;r<np;++r)
    {
        const double C = all[5*r], S = all[5*r+1], A = all[5*r+2];
        sC += C;  sS += S;  sA += A;  sF += all[5*r+3];  sN += all[5*r+4];
        sSS += S*S;  sSA += S*A;  sAA += A*A;  sCS += C*S;  sCA += C*A;
    }

    // least squares C_r = u_s S_r + u_a A_r over the ranks
    double us=0.0, ua=0.0;
    int fit=0;
    const double det = sSS*sAA - sSA*sSA;

    if(np>=2 && det > 1.0e-6*sSS*sAA)
    {
        us = (sCS*sAA - sCA*sSA)/det;
        ua = (sCA*sSS - sCS*sSA)/det;
        fit = us>0.0 && ua>=0.0;
    }

    if(!fit)
    {
        // an awake update costs about three sleeping ones (validation 22): u_a = 2 u_s
        const double den = sS + 2.0*sA;
        us = den>0.0 ? sC/den : 0.0;
        ua = 2.0*us;
    }

    const double uf = sN>0.0 ? sF/sN : 0.0;
    const double fs = uf>0.0 ? us/uf/st : 0.0;
    const double fa = uf>0.0 ? ua/uf/st : 0.0;

    // the cells of this rank, gathered to rank 0 into the map of the grid (the ranks hold disjoint cells)
    std::vector<double> loc(nc);

    for(size_t c=0;c<nc;++c)
    loc[c] = nf[c] + fs*lmS[c] + fa*lmA[c];

    int box[6] = {p->origin_i, p->j_dir==1 ? p->origin_j : 0, p->origin_k, p->knox, p->knoy, p->knoz};
    std::vector<int> boxes(p->mpirank==0 ? size_t(6)*np : 0), cnt(np,0), dsp(np,0);
    MPI_Gather(box,6,MPI_INT,p->mpirank==0 ? &boxes[0] : nullptr,6,MPI_INT,0,MPI_COMM_WORLD);

    std::vector<double> buf;

    if(p->mpirank==0)
    {
        size_t tot=0;

        for(int r=0;r<np;++r)
        {
            cnt[r] = boxes[6*r+3]*boxes[6*r+4]*boxes[6*r+5];
            dsp[r] = int(tot);
            tot += size_t(cnt[r]);
        }

        buf.resize(MAX(tot,size_t(1)));
    }

    MPI_Gatherv(nc>0 ? &loc[0] : nullptr, int(nc), MPI_DOUBLE, p->mpirank==0 ? &buf[0] : nullptr, &cnt[0], &dsp[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(p->mpirank!=0)
    return;

    std::vector<double> map(size_t(GNX)*GNY*GNZ,0.0);

    for(int r=0;r<np;++r)
    {
        const int *bx = &boxes[6*r];
        size_t c = size_t(dsp[r]);

        for(int qi=0;qi<bx[3];++qi)
        for(int qj=0;qj<bx[4];++qj)
        for(int qk=0;qk<bx[5];++qk, ++c)
        {
            const int I = bx[0]+qi, J = bx[1]+qj, K = bx[2]+qk;

            if(I>=0 && I<GNX && J>=0 && J<GNY && K>=0 && K<GNZ)
            map[(size_t(I)*GNY + J)*GNZ + K] += buf[c];
        }
    }

    const double ps = uf>0.0 ? us/uf : 0.0, pa = uf>0.0 ? ua/uf : 0.0;

    const streamsize pr = cout.precision();
    cout<<"CPM load map (G 45, G 46): a parcel update costs "<<setprecision(3)<<ps<<" fluid cells, an awake one "<<ps+pa
        <<(fit ? " (measured)" : " (assumed: three sleeping ones)")<<", CPM "<<int(100.0*sC/MAX(sC+sF,1.0e-20)+0.5)<<" % of the busy time"<<endl;
    cout.precision(pr);

    mkdir("./REEF3D_CFD_CPM_Particle",0777);
    ofstream out("./REEF3D_CFD_CPM_Particle/REEF3D_CFD_CPM_loadmap.dat");
    out<<"# CFD CPM load map (G 45): work per step of every cell in fluid cells, averaged over "<<lm_steps
       <<" steps: 1 per fluid cell + "<<ps<<" per parcel update + "<<pa<<" per awake parcel update"<<(fit ? " (measured)" : " (assumed)")<<endl;
    out<<"# for the weighted partition: DIVEMesh or REEF3D M 40 with M 20 5, M 35 1 (REEF3D M 41 1)"<<endl;
    out<<GNX<<" "<<GNY<<" "<<GNZ<<endl;
    out<<setprecision(6);

    size_t m=0;

    for(double v : map)
    out<<v<<((++m)%10==0 ? "\n" : " ");

    out<<endl;
}
