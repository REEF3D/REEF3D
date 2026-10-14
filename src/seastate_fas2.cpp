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

#include"seastate_fas2.h"
#include"seastate_obstacle.h"
#include"seastate_implicit.h"
#include"seastate_exchange.h"
#include"seastate_store.h"
#include"seastate_grid.h"
#include"fdm_seastate.h"
#include"slice.h"
#include"slice4.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<mpi.h>

seastate_fas2::seastate_fas2(lexer *p, ghostcell *pgc, fdm_seastate *e, seastate_implicit *fs, seastate_exchange *fx, seastate_source *s,
                             const vector<float> *const Nside[4], slice *wU, slice *wD, int depth_, int maxdepth,
                             const vector<char> *fine_coast, const seastate_obstacle *fine_obs)
              : p0(p), pc(nullptr), ef(e), ec(nullptr), fsolv(fs), csolv(nullptr), fex(fx), cex(nullptr), src(s), wUf(wU), wDf(wD), depth(depth_)
{
    nbin = e->grid->nbin;
    ms = p->margin;
    kx = (p->knox+1)/2;
    ky = (p->knoy+1)/2;

    for(int k=0; k<4; ++k)
    Nsf[k] = (Nside!=nullptr) ? Nside[k] : nullptr;

    build_lexer(p,pgc);
    build_fields(p,fine_coast,fine_obs);

    csolv = new seastate_implicit(pc,ec);
    csolv->sources(src);
    csolv->sparsity(0.0);
    csolv->geographic_order(1);
    csolv->threads(p->A798);
    if(p->A700==2)
    csolv->source_iterations(p->A738,p->A739);

    bool sides=false;
    for(int k=0; k<4; ++k)
    if(Nsf[k]!=nullptr)
    sides = true;
    if(sides)
    {
    const vector<float> *sc[4] = {&Nsc[0],&Nsc[1],&Nsc[2],&Nsc[3]};
    csolv->boundary_sides(sc);
    }

    if(wUf!=nullptr && wDf!=nullptr)
    {
    wUc = new slice4(pc);
    wDc = new slice4(pc);
    csolv->wind_field(wUc,wDc);
    }

    cex = new seastate_exchange(pc,nbin,1);

    // Phase 15: obstacles and coasts (A 722, A 723, A 725) on the coarse grid, its own faces as on a refinement patch
    if(fine_obs!=nullptr && fine_obs->active())
    {
    cobs = new seastate_obstacle(*fine_obs);
    cobs->build(pc,ec,ci0,cj0,cgx,cgy);
    csolv->obstacles(cobs);
    }

    forcing();

    long nw = 0;
    for(int I=0; I<kx; ++I)
    for(int J=0; J<ky; ++J)
    if(ec->wet(I,J)==1)
    ++nw;
    ncoarse_wet = long(pgc->globalsum(double(nw)));

    // the next coarser grid, on all ranks or none
    const double more = (depth<maxdepth && (kx>1 || ky>1)) ? 1.0 : 0.0;
        if(pgc->globalmin(more)>0.5)
        {
        const vector<float> *sc[4] = {&Nsc[0],&Nsc[1],&Nsc[2],&Nsc[3]};
        child = new seastate_fas2(pc,pgc,ec,csolv,cex,src,sides ? sc : nullptr,wUc,wDc,depth+1,maxdepth,&coast,cobs);
        }
}

seastate_fas2::~seastate_fas2()
{
    delete child;
    delete cobs;
    delete csolv;
    delete cex;
    delete Nt;
    delete Tc;
    delete N0c;
    delete wUc;
    delete wDc;

    if(ec!=nullptr)
    {
    delete ec->N;
    delete ec->kw;
    delete ec->cg;
    delete ec;
    }

    if(pc!=nullptr)
    {
    const int m = marge;
    const int nxa = pc->knox+1+4*m;
    const int nya = pc->knoy+1+4*m;
    const int nsl = pc->imax*pc->jmax;

    pc->del_Darray(pc->XN,nxa); pc->del_Darray(pc->XP,nxa); pc->del_Darray(pc->DXN,nxa); pc->del_Darray(pc->DXP,nxa);
    pc->del_Darray(pc->YN,nya); pc->del_Darray(pc->YP,nya); pc->del_Darray(pc->DYN,nya); pc->del_Darray(pc->DYP,nya);
    pc->del_Iarray(pc->flagslice1,nsl);
    pc->del_Iarray(pc->flagslice2,nsl);
    pc->del_Iarray(pc->flagslice4,nsl);
    pc->del_Iarray(pc->wet,nsl);
    pc->del_Iarray(pc->wet_n,nsl);
    pc->del_Iarray(pc->deep,nsl);
    if(gcount[0]>0) pc->del_Iarray(pc->gcslpara1,gcount[0],6);
    if(gcount[1]>0) pc->del_Iarray(pc->gcslpara4,gcount[1],6);
    if(gcount[2]>0) pc->del_Iarray(pc->gcslpara3,gcount[2],6);
    if(gcount[3]>0) pc->del_Iarray(pc->gcslpara2,gcount[3],6);
    delete pc;
    }
}

// first fine index under coarse index I: interior cells 2I (2I+1 too where it exists), the halo cells beyond the
// ends of the rank under the fine halo cells
int seastate_fas2::child0(int I, int n, int nc) const
{
    if(I<nc)
    return 2*I;
    return n + 2*(I-nc);
}

int seastate_fas2::children(int I, int n, int nc, int c[2]) const
{
    const int a = child0(I,n,nc);
    int k=0;
    for(int d=0; d<2; ++d)
    {
    const int f = a+d;
    if(f<-ms || f>=n+ms)
    continue;
    if(I>=0 && I<nc && f>=n)     // an interior coarse cell has only interior children
    continue;
    if(I<0 && f>=0)
    continue;
    c[k++] = f;
    }
    return k;
}

long seastate_fas2::coarse_bytes() const
{
    long b = 0;
    for(const seastate_store *s : {ec->N,ec->kw,ec->cg,Nt,Tc,N0c})
    if(s!=nullptr)
    b += long(s->bytes());
    return b + (child!=nullptr ? child->coarse_bytes() : 0);
}

// --------------------------------------------------------------------- coarse grid
void seastate_fas2::build_lexer(lexer *p, ghostcell *pgc)
{
    pc = new lexer(*p,1);
    const int m = marge;

    pc->knox = kx;
    pc->knoy = ky;
    pc->knoz = p->knoz;
    pc->imin = -ms;
    pc->jmin = -ms;
    pc->imax = kx+2*ms;
    pc->jmax = ky+2*ms;
    pc->kmin = p->kmin;
    pc->kmax = p->kmax;
    pc->kmaxF = p->kmaxF;
    pc->i_dir = p->i_dir; pc->j_dir = p->j_dir; pc->k_dir = p->k_dir;
    pc->x_dir = p->x_dir; pc->y_dir = p->y_dir; pc->z_dir = p->z_dir;

    // nodes: every second fine node, the last node of the rank, beyond the ends the fine halo nodes
    const int nxa = kx+1+4*m, nya = ky+1+4*m;
    const int pnxa = p->knox+1+4*m, pnya = p->knoy+1+4*m;

    pc->Darray(pc->XN,nxa); pc->Darray(pc->XP,nxa); pc->Darray(pc->DXN,nxa); pc->Darray(pc->DXP,nxa);
    pc->Darray(pc->YN,nya); pc->Darray(pc->YP,nya); pc->Darray(pc->DYN,nya); pc->Darray(pc->DYP,nya);
    pc->ZN = p->ZN; pc->ZP = p->ZP; pc->DZN = p->DZN; pc->DZP = p->DZP;

    auto fnode = [&](const double *X, int n, int pna, int N, int nc)
    {
        const int F = (N<=nc) ? std::min(2*N,n) : n + 2*(N-nc);
        const int a = F+m;
        if(a<0)     return X[0] + a*(X[1]-X[0]);
        if(a>=pna)  return X[pna-1] + (a-pna+1)*(X[pna-1]-X[pna-2]);
        return X[a];
    };

    for(int a=0; a<nxa; ++a)
    pc->XN[a] = fnode(p->XN,p->knox,pnxa,a-m,kx);
    for(int a=0; a<nya; ++a)
    pc->YN[a] = fnode(p->YN,p->knoy,pnya,a-m,ky);

    for(int a=0; a<nxa-1; ++a) { pc->XP[a] = 0.5*(pc->XN[a]+pc->XN[a+1]); pc->DXN[a] = pc->XN[a+1]-pc->XN[a]; }
    for(int a=0; a<nya-1; ++a) { pc->YP[a] = 0.5*(pc->YN[a]+pc->YN[a+1]); pc->DYN[a] = pc->YN[a+1]-pc->YN[a]; }
    pc->XP[nxa-1] = pc->XP[nxa-2] + pc->DXN[nxa-2]; pc->DXN[nxa-1] = pc->DXN[nxa-2];
    pc->YP[nya-1] = pc->YP[nya-2] + pc->DYN[nya-2]; pc->DYN[nya-1] = pc->DYN[nya-2];
    for(int a=0; a<nxa-1; ++a) pc->DXP[a] = pc->XP[a+1]-pc->XP[a];
    for(int a=0; a<nya-1; ++a) pc->DYP[a] = pc->YP[a+1]-pc->YP[a];
    pc->DXP[nxa-1] = pc->DXP[nxa-2];
    pc->DYP[nya-1] = pc->DYP[nya-2];

    pc->DXM = 2.0*p->DXM;
    pc->DYM = 2.0*p->DYM;
    pc->DXD = p->DXD;
    pc->DYD = p->DYD;

    // flags: a coarse cell is fluid where a child is
    const int nsl = pc->imax*pc->jmax;
    pc->Iarray(pc->flagslice1,nsl);
    pc->Iarray(pc->flagslice2,nsl);
    pc->Iarray(pc->flagslice4,nsl);
    pc->Iarray(pc->wet,nsl);
    pc->Iarray(pc->wet_n,nsl);
    pc->Iarray(pc->deep,nsl);

    auto cij = [&](int I, int J) {return (I-pc->imin)*pc->jmax + (J-pc->jmin);};
    auto fij = [&](int ii, int jj) {return (ii-p->imin)*p->jmax + (jj-p->jmin);};

    for(int I=pc->imin; I<pc->imin+pc->imax; ++I)
    for(int J=pc->jmin; J<pc->jmin+pc->jmax; ++J)
    {
    int ci[2], cj[2];
    const int na = children(I,p->knox,kx,ci), nb = children(J,p->knoy,ky,cj);
    int f = -1;
        for(int a=0; a<na; ++a)
        for(int b=0; b<nb; ++b)
        {
        const int v = p->flagslice4[fij(ci[a],cj[b])];
        if(a==0 && b==0)
        f = v;
        if(v>0)
        f = std::max(f,v);
        }
    pc->flagslice4[cij(I,J)] = f;
    }

    for(int I=pc->imin; I<pc->imin+pc->imax; ++I)
    for(int J=pc->jmin; J<pc->jmin+pc->jmax; ++J)
    {
    const int q = cij(I,J);
    const int f4 = pc->flagslice4[q];
    pc->flagslice1[q] = f4;
    pc->flagslice2[q] = f4;
    if(f4>0 && I+1<pc->imin+pc->imax && pc->flagslice4[cij(I+1,J)]<0)
    pc->flagslice1[q] = pc->flagslice4[cij(I+1,J)];
    if(f4>0 && J+1<pc->jmin+pc->jmax && pc->flagslice4[cij(I,J+1)]<0)
    pc->flagslice2[q] = pc->flagslice4[cij(I,J+1)];
    pc->wet[q] = pc->wet_n[q] = pc->deep[q] = 0;
    }

    // ranks: the neighbours of the fine grid; the exchange lists are the boundary cells of the coarse grid
    pc->nb1 = p->nb1; pc->nb2 = p->nb2; pc->nb3 = p->nb3; pc->nb4 = p->nb4; pc->nb5 = p->nb5; pc->nb6 = p->nb6;
    pc->mpirank = p->mpirank;
    pc->mpi_size = p->mpi_size;
    pc->gcslpara1_count = pc->gcslpara2_count = pc->gcslpara3_count = pc->gcslpara4_count = 0;
    pc->gcslpara1 = pc->gcslpara2 = pc->gcslpara3 = pc->gcslpara4 = nullptr;

    // side 0 i- (nb1, gcslpara1), 1 i+ (nb4, gcslpara4), 2 j- (nb3, gcslpara3), 3 j+ (nb2, gcslpara2)
    const int nbs[4] = {p->nb1,p->nb4,p->nb3,p->nb2};
    int **lists[4] = {nullptr,nullptr,nullptr,nullptr};
    for(int s=0; s<4; ++s)
    if(nbs[s]>=0)
    {
    const int n = (s<2) ? ky : kx;
    gcount[s] = n;
    pc->Iarray(lists[s],n,6);
        for(int k=0; k<n; ++k)
        {
        lists[s][k][0] = (s==0) ? 0 : (s==1) ? kx-1 : k;
        lists[s][k][1] = (s==2) ? 0 : (s==3) ? ky-1 : k;
        }
    }
    pc->gcslpara1 = lists[0]; pc->gcslpara1_count = gcount[0];
    pc->gcslpara4 = lists[1]; pc->gcslpara4_count = gcount[1];
    pc->gcslpara3 = lists[2]; pc->gcslpara3_count = gcount[2];
    pc->gcslpara2 = lists[3]; pc->gcslpara2_count = gcount[3];

    // the neighbours must hold the same number of coarse cells along the shared side (rank boxes in a
    // tensor-product decomposition)
    {
    static const int opp[4] = {1,0,3,2};
    int bad = 0;
    MPI_Request req[8];
    int nreq = 0;
    int theirs[4] = {-1,-1,-1,-1};
        for(int s=0; s<4; ++s)
        if(nbs[s]>=0)
        {
        MPI_Irecv(&theirs[s],1,MPI_INT,nbs[s],950+opp[s],pgc->mpi_comm,&req[nreq++]);
        MPI_Isend(&gcount[s],1,MPI_INT,nbs[s],950+s,pgc->mpi_comm,&req[nreq++]);
        }
        if(nreq>0)
        MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);
        for(int s=0; s<4; ++s)
        if(nbs[s]>=0 && theirs[s]!=gcount[s])
        bad = 1;
    bad = int(pgc->globalmax(double(bad)));
        if(bad)
        {
            if(p->mpirank==0)
            cout<<endl<<"SEASTATE input error  --  A 758: the FAS coarse grid on a uniform grid needs neighbouring ranks with the same"
                <<" number of cells along their common side (a tensor-product decomposition)"<<endl<<endl;
        pgc->final(true);
        }
    }

    // global index of the coarse cells for the boundary tests of the solver: only the sides of the domain count
    ci0 = (p->origin_i==0) ? 0 : p->origin_i/2+1;
    cj0 = (p->origin_j==0) ? 0 : p->origin_j/2+1;
    cgx = (p->origin_i+p->knox>=p->gknox) ? ci0+kx : ci0+kx+1000000;
    cgy = (p->origin_j+p->knoy>=p->gknoy) ? cj0+ky : cj0+ky+1000000;
    pc->origin_i = ci0;
    pc->origin_j = cj0;
    pc->origin_k = p->origin_k;
    pc->gknox = cgx;
    pc->gknoy = cgy;
    pc->gknoz = p->gknoz;

    pc->vec2Dlength = nsl;
    pc->veclength = 0;
    pc->cellnum2D = kx*ky;
    pc->slicenum = nsl;
    pc->wd = p->wd;
    pc->dt = p->dt;
    pc->simtime = p->simtime;
    pc->count = p->count;
}

void seastate_fas2::build_fields(lexer *p, const vector<char> *fine_coast, const seastate_obstacle *fine_obs)
{
    ec = new fdm_seastate(pc);
    ec->grid = ef->grid;

    // environment: the mean of the wet children
    for(int I=pc->imin; I<pc->imin+pc->imax; ++I)
    for(int J=pc->jmin; J<pc->jmin+pc->jmax; ++J)
    {
    int ci[2], cj[2];
    const int na = children(I,p->knox,kx,ci), nb = children(J,p->knoy,ky,cj);
    double d=0.0, b=0.0, u=0.0, v=0.0;
    int n=0;
        for(int a=0; a<na; ++a)
        for(int c=0; c<nb; ++c)
        if(ef->wet(ci[a],cj[c])==1)
        {
        d += ef->depth(ci[a],cj[c]);
        b += ef->bed(ci[a],cj[c]);
        u += ef->U(ci[a],cj[c]);
        v += ef->V(ci[a],cj[c]);
        ++n;
        }
    const int q = (I-pc->imin)*pc->jmax + (J-pc->jmin);
    const bool w = n>0 && pc->flagslice4[q]>0;
    ec->depth(I,J) = w ? d/n : 0.0;
    ec->bed(I,J)   = w ? b/n : 0.0;
    ec->U(I,J)     = w ? u/n : 0.0;
    ec->V(I,J)     = w ? v/n : 0.0;
    ec->eta(I,J)   = 0.0;
    ec->wet(I,J)   = w ? 1 : 0;
    ec->wet0(I,J)  = ec->wet(I,J);
    }

    // tiles of the same area as on the fine grid (a coarse tile of A 704 cells would cover 4 times the area)
    int tile = p->A704;
    for(int d=0; d<depth; ++d)
    tile = std::max(tile/2,2);
    auto store = [&](int nb) {seastate_store *s = new seastate_store(pc->imin,pc->jmin,pc->imax,pc->jmax,nb,tile); s->build(ec->wet.V); s->fill(0.0f); return s;};

    ec->N  = store(nbin);
    ec->kw = store(ec->grid->nsig);
    ec->cg = store(ec->grid->nsig);
    Nt = store(nbin);
    Tc = store(nbin);
    if(p->A700==1)
    N0c = store(nbin);


    seastate_kinematics(pc,ec,0.0,ci0,cj0,cgx,cgy);

    // coarse cells over a dry cell of the original grid (no correction there): a dry child, or (the coarser grids)
    // a child that covers one
    coast.assign(size_t(kx)*ky,0);
    for(int I=0; I<kx; ++I)
    for(int J=0; J<ky; ++J)
    {
    int ci[2], cj[2];
    const int na = children(I,p->knox,kx,ci), nb = children(J,p->knoy,ky,cj);
        for(int a=0; a<na; ++a)
        for(int c=0; c<nb; ++c)
        if(ef->wet(ci[a],cj[c])!=1 || (fine_coast!=nullptr && (*fine_coast)[size_t(ci[a])*p->knoy+cj[c]]))
        coast[size_t(I)*ky+J] = 1;

        // Phase 15: a blocked face of the finer grid between the children (obstacle, structure)
        if(fine_obs!=nullptr && fine_obs->active())
        {
        if(na==2)
        for(int c=0; c<nb; ++c)
        if(fine_obs->east(ci[0],cj[c])!=nullptr)
        coast[size_t(I)*ky+J] = 1;
        if(nb==2)
        for(int a=0; a<na; ++a)
        if(fine_obs->north(ci[a],cj[0])!=nullptr)
        coast[size_t(I)*ky+J] = 1;
        }
    }
}

// boundary spectra per side cell (A 711 3): the mean of the fine side cells; wind field (A 730 2): the first wet child
void seastate_fas2::forcing()
{
    for(int s=0; s<4; ++s)
    {
    if(Nsf[s]==nullptr || Nsf[s]->empty())
    continue;
    const int nf = (s<2) ? p0->knoy : p0->knox;
    const int nc = (s<2) ? ky : kx;
    Nsc[s].assign(size_t(nc)*nbin,0.0f);
        for(int K=0; K<nc; ++K)
        {
        const int a = 2*K, b = std::min(2*K+1,nf-1);
        const float *A = Nsf[s]->data() + size_t(a)*nbin, *B = Nsf[s]->data() + size_t(b)*nbin;
        float *C = Nsc[s].data() + size_t(K)*nbin;
        for(int n=0; n<nbin; ++n)
        C[n] = 0.5f*(A[n]+B[n]);
        }
    }

    if(wUc!=nullptr)
    for(int I=pc->imin; I<pc->imin+pc->imax; ++I)
    for(int J=pc->jmin; J<pc->jmin+pc->jmax; ++J)
    {
    int ci[2], cj[2];
    const int na = children(I,p0->knox,kx,ci), nb = children(J,p0->knoy,ky,cj);
    bool done=false;
        for(int a=0; a<na && !done; ++a)
        for(int c=0; c<nb && !done; ++c)
        if(ef->wet(ci[a],cj[c])==1 || (a==na-1 && c==nb-1))
        {
        (*wUc)(I,J) = (*wUf)(ci[a],cj[c]);
        (*wDc)(I,J) = (*wDf)(ci[a],cj[c]);
        done = true;
        }
    }
}

// --------------------------------------------------------------------- correction
void seastate_fas2::correct(lexer *p, ghostcell *pgc, const seastate_store *N0, double rdt, const vector<float> &Nb, const int side[4],
                            bool refraction, bool fshift, int ncoarse)
{
    forcing();

    // 1: residual of the fine cells at the latest spectra (no solve), summed over the children of every coarse cell
    // into tau (residual_add)
    Tc->fill(0.0f);
    fsolv->residual_add(Tc,1);
        for(int q=0; q<4; ++q)
        fsolv->sweep(p,ef,q,N0,rdt,Nb,side,refraction,fshift);
    fsolv->residual_add(nullptr,0);

    // 2: restriction of N (N0): the means over the wet children; tau holds -R r (the mean of the residuals) for now
    Nt->fill(0.0f);
    if(N0c!=nullptr)
    N0c->fill(0.0f);

    for(int I=0; I<kx; ++I)
    for(int J=0; J<ky; ++J)
    {
    if(ec->wet(I,J)!=1)
    continue;
    float *nt = Nt->spec(I,J), *tc = Tc->spec(I,J);
    float *n0 = (N0c!=nullptr) ? N0c->spec(I,J) : nullptr;
    if(nt==nullptr || tc==nullptr)
    continue;

    int ci[2], cj[2];
    const int na = children(I,p->knox,kx,ci), nb = children(J,p->knoy,ky,cj);
    int n=0;
        for(int a=0; a<na; ++a)
        for(int c=0; c<nb; ++c)
        {
        if(ef->wet(ci[a],cj[c])!=1)
        continue;
        const float *f = ef->N->spec(ci[a],cj[c]);
        const float *f0 = (n0!=nullptr && N0!=nullptr) ? N0->spec(ci[a],cj[c]) : nullptr;
        if(f==nullptr)
        continue;
            for(int b=0; b<nbin; ++b)
            nt[b] += f[b];
            if(f0!=nullptr)
            for(int b=0; b<nbin; ++b)
            n0[b] += f0[b];
        ++n;
        }
    const float w = (n>0) ? 1.0f/float(n) : 0.0f;
        for(int b=0; b<nbin; ++b)
        {
        nt[b] *= w;
        tc[b] *= -w;
        }
        if(n0!=nullptr)
        for(int b=0; b<nbin; ++b)
        n0[b] *= w;
    }

    ec->N->copy_from(*Nt);
    cex->start(pc,pgc,*ec->N);

    // Phase 15: the coefficients of the coarse faces (Goda, d'Angremond, porous) from the restricted spectra
    if(cobs!=nullptr)
    cobs->update(pc,ec);

    // 3: the residual of the coarse problem at the restricted spectra added: tau = r_c(N_t) - R r
    csolv->residual_add(Tc,0);
        for(int q=0; q<4; ++q)
        csolv->sweep(pc,ec,q,N0c,rdt,Nb,side,refraction,fshift);
    csolv->residual_add(nullptr,0);

    // 4: the coarse problem: m iterations; with a coarser grid m/2, its correction (V-cycle), m/2; the coarsest
    // grid 4 m iterations
    auto coarse_iterations = [&](int n)
    {
        for(int it=0; it<n; ++it)
        for(int q=0; q<4; ++q)
        {
        csolv->sweep(pc,ec,q,N0c,rdt,Nb,side,refraction,fshift);
        cex->start(pc,pgc,*ec->N);
        }
    };

    csolv->fas_source(Tc,Nt,p0->A795>0.0);
        if(child!=nullptr)
        {
        const int pre = std::max(ncoarse/2,1);
        coarse_iterations(pre);
        child->correct(pc,pgc,N0c,rdt,Nb,side,refraction,fshift,ncoarse);
        coarse_iterations(std::max(ncoarse-pre,1));
        }
        else
        coarse_iterations(4*ncoarse);
    csolv->fas_source(nullptr,nullptr);

    // 5: the coarse change into the children (multiplicative within 1/4 and 4, additive where N_t vanishes), not
    // under coarse cells with a dry child
    vector<float> fac(nbin), add(nbin);
    double chg = 0.0;

    for(int I=0; I<kx; ++I)
    for(int J=0; J<ky; ++J)
    {
    if(ec->wet(I,J)!=1 || coast[size_t(I)*ky+J])
    continue;
    const float *nc = ec->N->spec(I,J), *nt = Nt->spec(I,J);
    if(nc==nullptr || nt==nullptr)
    continue;

    double et = 0.0, dc = 0.0;
        for(int b=0; b<nbin; ++b)
        {
        const double t = double(nt[b]), c = double(nc[b]);
        et += t;
        dc += std::fabs(c-t);
        fac[b] = 1.0f;
        add[b] = 0.0f;
            if(t>1.0e-30)
            fac[b] = float(std::min(std::max(c/t,0.25),4.0));
            else if(c>t)
            add[b] = float(c-t);
        }
    if(et>0.0)
    chg = std::max(chg,dc/et);

    int ci[2], cj[2];
    const int na = children(I,p->knox,kx,ci), nb = children(J,p->knoy,ky,cj);
        for(int a=0; a<na; ++a)
        for(int c=0; c<nb; ++c)
        {
        float *f = ef->N->spec(ci[a],cj[c]);
        if(f!=nullptr && ef->wet(ci[a],cj[c])==1)
        for(int b=0; b<nbin; ++b)
        f[b] = f[b]*fac[b] + add[b];
        }
    }
    change = pgc->globalmax(chg);

    fex->start(p,pgc,*ef->N);
}
