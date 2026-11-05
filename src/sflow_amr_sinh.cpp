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

#include"sflow_amr.h"
#include"sflow_momentum_SI.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"slice4.h"
#include<cmath>
#include<mpi.h>

//  The coupled eta-P system of the semi-implicit step with the quadratic non-hydrostatic pressure (A 226 1, A 220 2/3)
//  over the leaf cells of all grids of the mesh refinement (sflow_amr::si_step; the rows of one grid:
//  sflow_momentum_SI_nh_amr.cpp).
//
//  Unknowns: eta and P of the leaf cells. The operator: the leaf values into nE, nP of every grid, the covered cells the
//  mean of their children and the cells around the patches the values of the cells they stand for or the linear
//  interpolation of the coarser level (rsi_fill); h^3 lap(eta) and the face discharges on every grid; a coarse face
//  next to a finer patch takes the mean of the two fine face discharges (rsm), so that the mass rows are the
//  divergence of the discharges the face update takes (exact volume, as the hydrostatic eta system of
//  reefamr_si.cpp); then the rows of the leaf cells. FGMRES(30), relative tolerance 1e-6, right preconditioned: the
//  preconditioner of nh_solve on every grid (eta eliminated locally, one REEFMG2D V-cycle for psi, a second one for the
//  A 220 3 term), level 0 over its own cells with the mean residual of the children in the covered cells, every patch
//  over its own cells with 0 around it (additive). The wetting and drying iteration (A 229 1) runs over all grids.

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

typedef reefamr_comms_off comms_off;
}

void sflow_amr::si_nh_solve(lexer *p, ghostcell *pgc, rsi_result &res)
{
    const double t0 = MPI_Wtime();
    const vector<int> &G = rsi_grids();
    const int ng = (int)G.size();

    // the grids of the system, their leaf cells and the own cells (leaf and covered)
    struct gl { int g; lexer *q; fdm2D *b; sflow_momentum_SI *m; vector<int> leaf; int off; };
    vector<gl> L(ng);
    int nc = 0;

    for(int k=0; k<ng; ++k)
    {
        const int g = G[k];
        L[k].g = g;
        L[k].q = glex(g);
        L[k].b = (g<0) ? b0 : SP(g)->b;
        L[k].m = gsi(g);
        L[k].off = nc;

        lexer *q = L[k].q;
        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        if(q->flagslice4[lij(q,ii,jj)]>0 && rsi_leaf(g,ii,jj))
        {
        L[k].leaf.push_back(ii);
        L[k].leaf.push_back(jj);
        }

        nc += (int)L[k].leaf.size()/2;
    }

    // per grid on its own: f(grid entry), patches with the communication off
    auto each = [&](const std::function<void(gl&)> &f)
    {
        for(auto &c : L)
        {
            if(c.g<0)
            {
            f(c);
            continue;
            }

            comms_off guard(pgc);
            f(c);
        }
    };

    // the coefficients of every grid (the cells around the grids hold the filled WLq, phi, Lc0 of si_step)
    each([&](gl &c) { c.m->nh_assemble(c.q,c.b,pgc); c.m->nh_amr_cells(c.q,c.b); c.m->nh_amr_guess(c.q,c.b); });

    // own cells of grid g (the rows of its preconditioner): leaf and covered cells
    auto ownf = [&](int g) { return [this,g](int ii, int jj) { return rsi_own(g,ii,jj); }; };
    auto leaff = [&](int g) { return [this,g](int ii, int jj) { return rsi_leaf(g,ii,jj); }; };

    const bool m0nh_b = (psim0->Bdisp>0.0);
    const int N2 = 2*nc;
    std::vector<double> xv(N2), yv(N2), rv(N2), zv(N2);

    auto gather = [&](std::vector<double> &x, bool eta_p)
    {
        for(auto &c : L)
        for(size_t k=0; k<c.leaf.size(); k+=2)
        {
        const int ii=c.leaf[k], jj=c.leaf[k+1], n=c.off+(int)k/2;
        x[2*n]   = eta_p ? c.m->nE(ii,jj) : c.m->nZE(ii,jj);
        x[2*n+1] = eta_p ? c.m->nP(ii,jj) : c.m->nZP(ii,jj);
        }
    };

    auto scatter_fill = [&](const std::vector<double> &x)
    {
        for(auto &c : L)
        for(size_t k=0; k<c.leaf.size(); k+=2)
        {
        const int ii=c.leaf[k], jj=c.leaf[k+1], n=c.off+(int)k/2;
        c.m->nE(ii,jj) = x[2*n];
        c.m->nP(ii,jj) = x[2*n+1];
        }

        rsi_fill(9100,[&](int g) -> slice& { return gsi(g)->nE; });
        rsi_fill(9200,[&](int g) -> slice& { return gsi(g)->nP; });

        // A 220 3: h^3 lap(eta) next to the patch edges needs eta of the cells around the patch to second order
        if(m0nh_b)
        si_nh_quad([&](int g) -> slice& { return gsi(g)->nE; },false);
    };

    // the face discharges of the coarse faces next to the patches: the mean of the fine ones
    auto match = [&]()
    {
        rsm_begin(2);

        for(auto &c : L)
        {
            if(c.g<0 || P[c.g]->lev<1)
            continue;

            sflow_momentum_SI *m = c.m;
            for(int side=0; side<4; ++side)
            for(int r=0; r<rsm_count(c.g,side); ++r)
            {
                int dir,fi,fj;
                rsm_face(c.g,side,r,dir,fi,fj);
                double *s = rsm_slot(c.g,side,r);
                s[0] = (dir==1) ? m->nF1x(fi,fj) : m->nF1y(fi,fj);
                s[1] = (dir==1) ? m->nF2x(fi,fj) : m->nF2y(fi,fj);
            }
        }

        for(int l=maxlev; l>=1; --l)
        rsm_send(l,9300);

        for(auto &c : L)
        {
            sflow_momentum_SI *m = c.m;
            rsm_take(c.g,[&](int dir, int fi, int fj, const double *v, int sgn)
            {
                if(dir==1)
                {
                m->nF1x(fi,fj) = sgn*v[0];
                m->nF2x(fi,fj) = sgn*v[1];
                }
                else
                {
                m->nF1y(fi,fj) = sgn*v[0];
                m->nF2y(fi,fj) = sgn*v[1];
                }
            });
        }
    };

    // y = A x (affine: with the constant terms)
    auto apply = [&](const std::vector<double> &x, std::vector<double> &y, bool affine)
    {
        scatter_fill(x);
        each([&](gl &c) { c.m->nh_amr_lap(c.q,pgc); c.m->nh_amr_faces(c.q,affine); });
        match();

        for(auto &c : L)
        for(size_t k=0; k<c.leaf.size(); k+=2)
        {
        const int ii=c.leaf[k], jj=c.leaf[k+1], n=c.off+(int)k/2;
        c.m->nh_amr_row(c.q,c.b,ii,jj,affine,y[2*n],y[2*n+1]);
        }
    };

    auto dot = [&](const std::vector<double> &a, const std::vector<double> &c)
    {
        double s = 0.0;
        for(int k=0; k<N2; ++k)
        s += a[k]*c[k];
        return pgc->globalsum(s);
    };

    // z = M^-1 v: the residual of the leaf cells, the mean of the children in the covered cells
    auto precond = [&](const std::vector<double> &v, std::vector<double> &z)
    {
        for(auto &c : L)
        {
            lexer *q = c.q;
            for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
            for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
            {
            c.m->nRm(ii,jj) = 0.0;
            c.m->nRc(ii,jj) = 0.0;
            }

            for(size_t k=0; k<c.leaf.size(); k+=2)
            {
            const int ii=c.leaf[k], jj=c.leaf[k+1], n=c.off+(int)k/2;
            c.m->nRm(ii,jj) = v[2*n];
            c.m->nRc(ii,jj) = v[2*n+1];
            }
        }

        rsi_restrict(9400,[&](int g) -> slice& { return gsi(g)->nRm; });
        rsi_restrict(9500,[&](int g) -> slice& { return gsi(g)->nRc; });

        // level 0 first
        for(auto &c : L)
        if(c.g<0)
        c.m->nh_amr_prec(c.q,pgc,ownf(c.g));

        // the patches of level 1 (parent: the native grid on this rank): the level-0 correction interpolated into
        // the patch and its ring, the residual it leaves on the patch, the patch preconditioner on that (two-level,
        // multiplicative); finer levels: their preconditioner on the residual alone (additive)
        sflow_momentum_SI *m0 = psim0;
        lexer *q0 = p0;
        auto inq0 = [&](int a, int d) { return a>=q0->imin && a<q0->imin+q0->imax && d>=q0->jmin && d<q0->jmin+q0->jmax; };

        for(auto &c : L)
        {
            if(c.g<0)
            continue;

            comms_off guard(pgc);
            sflow_momentum_SI *m = c.m;
            lexer *q = c.q;
            reefamr_patch *P0 = P[c.g];
            const bool two = (P0->lev==1 && nblk==0);

            if(two)
            {
                // the interpolated level-0 correction on all cells of the patch
                for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
                for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
                {
                    double ze=0.0, zp=0.0;
                    const int I = ii-EXT+P0->I0, J = jj-EY+P0->J0;
                    const int ic = (I>>1) - O0i, jc = yparent(J) - O0j;

                    if(q->flagslice4[lij(q,ii,jj)]>0 && inq0(ic,jc) && q0->flagslice4[lij(q0,ic,jc)]>0)
                    {
                    const int ox = (I&1) ? 1 : -1, oy = (sy==1) ? ((J&1) ? 1 : -1) : 0;
                    ze = rsi_interp(q0,m0->nZE,ic,jc,ox,oy);
                    zp = rsi_interp(q0,m0->nZP,ic,jc,ox,oy);
                    }

                    m->nE(ii,jj) = ze;
                    m->nP(ii,jj) = zp;
                }

                m->nh_amr_lap(q,pgc);
                m->nh_amr_faces(q,false);

                for(size_t k=0; k<c.leaf.size(); k+=2)
                {
                    const int ii=c.leaf[k], jj=c.leaf[k+1];
                    double ym, yc;
                    m->nh_amr_row(q,c.b,ii,jj,false,ym,yc);
                    m->nRm(ii,jj) -= ym;
                    m->nRc(ii,jj) -= yc;
                }
            }

            m->nh_amr_prec(q,pgc,ownf(c.g));

            if(two)
            for(size_t k=0; k<c.leaf.size(); k+=2)
            {
                const int ii=c.leaf[k], jj=c.leaf[k+1];
                m->nZE(ii,jj) += m->nE(ii,jj);
                m->nZP(ii,jj) += m->nP(ii,jj);
            }
        }

        gather(z,false);
    };

    // the matrices of the preconditioners
    each([&](gl &c) { c.m->nh_amr_prec_setup(c.q,pgc,ownf(c.g)); });

    gather(xv,true);

    const int m = 30, maxit = 300;
    const double tol = 1.0e-6;
    std::vector<std::vector<double>> V(m+1,std::vector<double>(N2)), Z(m,std::vector<double>(N2));
    std::vector<double> H((m+1)*m,0.0), cs(m,0.0), sn(m,0.0), gg(m+1,0.0);

    res.iter = 0;
    res.wd_iter = 0;
    res.dry = 0;

    for(int wdit=0; wdit<(p->A229==1 ? 20 : 1); ++wdit)
    {
        ++res.wd_iter;

        // right-hand side: (rhsW - PW d, 0) - A(0)
        std::fill(zv.begin(),zv.end(),0.0);
        apply(zv,yv,true);

        for(auto &c : L)
        for(size_t k=0; k<c.leaf.size(); k+=2)
        {
        const int ii=c.leaf[k], jj=c.leaf[k+1], n=c.off+(int)k/2;
        rv[2*n]   = c.m->rhsW(ii,jj) - c.m->nPW(ii,jj)*c.b->depth(ii,jj) - yv[2*n];
        rv[2*n+1] = -yv[2*n+1];
        }

        const double bn = sqrt(dot(rv,rv));
        int it = 0;

        while(it<maxit)
        {
            apply(xv,yv,false);

            for(int k=0; k<N2; ++k)
            V[0][k] = rv[k] - yv[k];

            const double beta = sqrt(dot(V[0],V[0]));

            if(beta<=tol*MAX(bn,1.0e-30))
            break;

            for(int k=0; k<N2; ++k)
            V[0][k] /= beta;

            std::fill(gg.begin(),gg.end(),0.0);
            gg[0] = beta;

            int jj=0;
            for(; jj<m && it<maxit; ++jj, ++it)
            {
                precond(V[jj],Z[jj]);
                apply(Z[jj],V[jj+1],false);

                for(int l=0; l<=jj; ++l)
                {
                const double h = dot(V[jj+1],V[l]);
                H[l*m+jj] = h;
                for(int k=0; k<N2; ++k)
                V[jj+1][k] -= h*V[l][k];
                }

                const double h1 = sqrt(dot(V[jj+1],V[jj+1]));
                H[(jj+1)*m+jj] = h1;

                if(h1>0.0)
                for(int k=0; k<N2; ++k)
                V[jj+1][k] /= h1;

                for(int l=0; l<jj; ++l)
                {
                const double a0 = cs[l]*H[l*m+jj] + sn[l]*H[(l+1)*m+jj];
                H[(l+1)*m+jj] = -sn[l]*H[l*m+jj] + cs[l]*H[(l+1)*m+jj];
                H[l*m+jj] = a0;
                }

                const double a0 = H[jj*m+jj], a1 = H[(jj+1)*m+jj];
                const double rr = sqrt(a0*a0 + a1*a1);
                cs[jj] = rr>0.0 ? a0/rr : 1.0;
                sn[jj] = rr>0.0 ? a1/rr : 0.0;
                H[jj*m+jj] = rr;
                H[(jj+1)*m+jj] = 0.0;
                gg[jj+1] = -sn[jj]*gg[jj];
                gg[jj]   =  cs[jj]*gg[jj];

                if(fabs(gg[jj+1])<=tol*MAX(bn,1.0e-30))
                {
                ++jj; ++it;
                break;
                }
            }

            std::vector<double> yk(jj,0.0);
            for(int l=jj-1; l>=0; --l)
            {
            double s = gg[l];
            for(int q=l+1; q<jj; ++q)
            s -= H[l*m+q]*yk[q];
            yk[l] = s/H[l*m+l];
            }

            for(int l=0; l<jj; ++l)
            for(int k=0; k<N2; ++k)
            xv[k] += yk[l]*Z[l][k];

            if(jj>0 && fabs(gg[jj])<=tol*MAX(bn,1.0e-30))
            break;
        }

        res.iter += it;

        if(p->A229!=1)
        break;

        // the wet pattern of the leaf cells (si_eta2D::wd_update), over all grids
        scatter_fill(xv);
        int change = 0;
        for(auto &c : L)
        {
            if(c.g<0)
            change += c.m->nh_amr_wetdry(c.q,c.b,leaff(c.g));
            else
            {
            comms_off guard(pgc);
            change += c.m->nh_amr_wetdry(c.q,c.b,leaff(c.g));
            }
        }

        change = pgc->globalisum(change);

        if(change==0)
        break;

        // the wet pattern changes the preconditioner of the mass rows
        each([&](gl &c) { c.m->nh_amr_prec_setup(c.q,pgc,ownf(c.g)); });
    }

    for(auto &c : L)
    for(size_t k=0; k<c.leaf.size(); k+=2)
    if(c.m->nPW(c.leaf[k],c.leaf[k+1])<0.5)
    ++res.dry;

    res.dry = pgc->globalisum(res.dry);

    // the solution on all grids (covered cells and the cells around the grids filled), the face terms
    scatter_fill(xv);
    each([&](gl &c) { c.m->nh_amr_lap(c.q,pgc); });
    each([&](gl &c) { c.m->nh_amr_finish(c.q,c.b,pgc); });

    res.time = MPI_Wtime() - t0;

    for(auto &c : L)
    {
    c.m->nh_iter = res.iter;
    c.m->nh_time = res.time;
    }
}

// the cells around the patches that the coarser grid fills (kind 1, on this rank): quadratic interpolation of the
// coarse values instead of the linear one of rsi_fill / the flow fills. The second differences of the linear
// prolongation jump between the parent cells, so that h^3 lap(eta) of the A 220 3 term next to the patch edges was
// wrong by O(1) (0.5 eta'' with the refinement 2): pressure spikes at the patch edges, with patches moving with the
// waves (G 2 4) the Beji & Battjes bar blew up. Cells from another rank keep the linear value. wetonly: only where
// the parent and its neighbours are wet (the surface of the flow fills at a shoreline stays as it is)
void sflow_amr::si_nh_quad(const std::function<slice&(int)> &F, bool wetonly)
{
    for(int l=1; l<=maxlev && l<(int)lev.size(); ++l)
    for(int n : lev[l])
    {
        reefamr_patch *c = P[n];
        slice &X = F(n);

        for(const reefamr_fill &f : c->fill)
        {
            if(f.kind!=1)
            continue;

            lexer *q = glex(f.g);
            slice &S = F(f.g);
            const int a = f.si, d = f.sj;

            auto ok = [&](int u, int v) { return u>=q->imin && u<q->imin+q->imax && v>=q->jmin && v<q->jmin+q->jmax
                                                 && q->flagslice4[lij(q,u,v)]>0 && (!wetonly || q->wet[lij(q,u,v)]==1); };

            if(!ok(a,d))
            continue;

            const double c0 = S(a,d);
            auto val = [&](int u, int v) { return ok(u,v) ? S(u,v) : c0; };

            const double e = val(a+1,d), w = val(a-1,d);
            double v = c0 + 0.125*(e - w)*f.ox + (e - 2.0*c0 + w)/32.0;

            if(sy==1)
            {
                const double nn = val(a,d+1), ss = val(a,d-1);
                v += 0.125*(nn - ss)*f.oy + (nn - 2.0*c0 + ss)/32.0;

                if(ok(a+1,d+1) && ok(a+1,d-1) && ok(a-1,d+1) && ok(a-1,d-1))
                v += (S(a+1,d+1) - S(a+1,d-1) - S(a-1,d+1) + S(a-1,d-1))*f.ox*f.oy/64.0;
            }

            X(f.di,f.dj) = v;
        }
    }
}


// the face discharges of the coarse faces inside the patches (east and north face of a covered cell: the mean of the
// two fine faces on it), after the restriction of the cells: when a patch goes away, the uncovered cells start from
// the fine solution and their faces from the fine discharges (q^n of the next step), instead of the discharges of the
// coarse solution under the patch
void sflow_amr::si_nh_qrestrict()
{
    for(int l=maxlev; l>=1; --l)
    block_up(l,2,9700+l,
             [&](reefamr_patch *c, int id, int k, double *v)
             {
                 sflow_momentum_SI *m = gsi(id);
                 const int nby = c->ny>>sy;
                 const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);

                 // east face of the block: the fine x-faces (i0+1, j0..j0+sy); north face: the fine y-faces (i0..i0+1, j0+1);
                 // the faces on the box edge are matched (si_match)
                 v[0] = (i0+1>=EXT+c->nx-1) ? NAN : ((sy==0) ? m->qx(i0+1,j0) : 0.5*(m->qx(i0+1,j0) + m->qx(i0+1,j0+1)));
                 v[1] = (sy==0 || j0+1>=EY+c->ny-1) ? NAN : 0.5*(m->qy(i0,j0+1) + m->qy(i0+1,j0+1));
             },
             [&](const reefamr_block &B, int key, const double *v)
             {
                 sflow_momentum_SI *m = gsi(B.g);
                 if(v[0]==v[0])
                 m->qx(B.ic,B.jc) = v[0];
                 if(sy==1 && v[1]==v[1])
                 m->qy(B.ic,B.jc) = v[1];
             });
}
