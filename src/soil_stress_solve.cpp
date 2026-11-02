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

#include"soil_stress.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>
#include<algorithm>
#include<Eigen/Sparse>

using namespace std;

namespace
{
// shear strain energy of an edge with the cells A (lo, left), B (lo, right), C (hi, left), D (hi, right):
// gamma = du/dr (u of the rows lo: A-B, hi: C-D, across the rows) + dv/dc (v of the columns left: A-C,
// right: B-D, across the columns). A fluid or roller cell: free of shear. A row (column) of two walls:
// the wall is at the edge, u = 0 there, half the distance. de_r, de_c: the extent of the edge volume
// across the rows and the columns (zero for an edge without energy).
double shear(int tA, int tB, int tC, int tD,
             double ulo, double uhi, double hlo, double hhi, double hr,
             double vl, double vr, double hl, double hrr, double hc,
             double G, double &de_r, double &de_c)
{
    const int S = soil_stress::SOIL, W = soil_stress::WALL;
    de_r = de_c = 0.0;

    if(tA!=S && tB!=S && tC!=S && tD!=S)
    return 0.0;

    auto solidlike = [&](int t) {return t==S || t==W;};
    if(!solidlike(tA) || !solidlike(tB) || !solidlike(tC) || !solidlike(tD))
    return 0.0;

    double g1, g2;

    if(tA==W && tB==W)
    {
        de_r = 0.5*hhi;
        g1 = uhi/de_r;
    }
    else if(tC==W && tD==W)
    {
        de_r = 0.5*hlo;
        g1 = -ulo/de_r;
    }
    else
    {
        de_r = hr;
        g1 = (uhi-ulo)/de_r;
    }

    if(tA==W && tC==W)
    {
        de_c = 0.5*hrr;
        g2 = vr/de_c;
    }
    else if(tB==W && tD==W)
    {
        de_c = 0.5*hl;
        g2 = -vl/de_c;
    }
    else
    {
        de_c = hc;
        g2 = (vr-vl)/de_c;
    }

    return G*(g1+g2);
}
}

// the preconditioner: per rank (block Jacobi over the ranks) the Cholesky factors of the elastic block A
// (+ c0 rho M) and of the pressure block S + beta V + T/a0, beta V standing for B A^-1 B^T
// (beta = 1/(Ev C33), the constrained modulus); refactorised when the soil region changes or a0, c0 change
// by more than 20 %
struct soil_stress::factor
{
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > llt_u, llt_p;
    std::vector<int> lidx;          // local number of an unknown (block*N + q), -1 for none
    std::vector<int> lu, lp;        // the unknowns (block*N + q) of the local numbers
    std::vector<double> ctype_f;    // the soil region of the factors
    double a0f=0.0, c0f=0.0;
    bool ok_u=false, ok_p=false;
    Eigen::VectorXd ru, rp, zu, zp;
};

soil_stress::~soil_stress()
{
    delete fac;
}

// txz at x between the cells i, i+1 and z between k, k+1; dxe, dze: the extent of its volume
double soil_stress::edge_xz(int i, int j, int k, const vget &g, double &dxe, double &dze) const
{
    return shear(type(i,j,k),type(i+1,j,k),type(i,j,k+1),type(i+1,j,k+1),
                 U(g,0,i,j,k),U(g,0,i,j,k+1),DZN(k),DZN(k+1),DZP(k),
                 U(g,2,i,j,k),U(g,2,i+1,j,k),DXN(i),DXN(i+1),DXP(i),
                 Ev*Gv,dze,dxe);
}

// tyz at y between j, j+1 and z between k, k+1
double soil_stress::edge_yz(int i, int j, int k, const vget &g, double &dye, double &dze) const
{
    return shear(type(i,j,k),type(i,j+1,k),type(i,j,k+1),type(i,j+1,k+1),
                 U(g,1,i,j,k),U(g,1,i,j,k+1),DZN(k),DZN(k+1),DZP(k),
                 U(g,2,i,j,k),U(g,2,i,j+1,k),DYN(j),DYN(j+1),DYP(j),
                 Ev*Gv,dze,dye);
}

// txy at x between i, i+1 and y between j, j+1
double soil_stress::edge_xy(int i, int j, int k, const vget &g, double &dxe, double &dye) const
{
    return shear(type(i,j,k),type(i+1,j,k),type(i,j+1,k),type(i+1,j+1,k),
                 U(g,0,i,j,k),U(g,0,i,j+1,k),DYN(j),DYN(j+1),DYP(j),
                 U(g,1,i,j,k),U(g,1,i+1,j,k),DXN(i),DXN(i+1),DXP(i),
                 Ev*Gh,dye,dxe);
}

double soil_stress::strain_vol(int i, int j, int k, const vget &g) const
{
    double ev = (U(g,0,i,j,k)-U(g,0,i-1,j,k))/DXN(i) + (U(g,2,i,j,k)-U(g,2,i,j,k-1))/DZN(k);
    if(dim3)
    ev += (U(g,1,i,j,k)-U(g,1,i,j-1,k))/DYN(j);
    return ev;
}

void soil_stress::normal_stress(int i, int j, int k, const vget &g, double &sxx, double &syy, double &szz) const
{
    sxx = syy = szz = 0.0;
    if(type(i,j,k)!=SOIL)
    return;

    const double exx = (U(g,0,i,j,k)-U(g,0,i-1,j,k))/DXN(i);
    const double ezz = (U(g,2,i,j,k)-U(g,2,i,j,k-1))/DZN(k);
    const double eyy = dim3 ? (U(g,1,i,j,k)-U(g,1,i,j-1,k))/DYN(j) : 0.0;

    sxx = Ev*(C11*exx + C12*eyy + C13*ezz);
    syy = Ev*(C12*exx + C11*eyy + C13*ezz);
    szz = Ev*(C13*(exx+eyy) + C33*ezz);
}

// lumped mass of a face: half the volume of its soil cells
double soil_stress::face_mass(int c, int i, int j, int k) const
{
    int i2=i, j2=j, k2=k;
    if(c==0) ++i2;
    if(c==1) ++j2;
    if(c==2) ++k2;

    double m = 0.0;
    if(type(i,j,k)==SOIL)
    m += 0.5*DXN(i)*DYN(j)*DZN(k);
    if(type(i2,j2,k2)==SOIL)
    m += 0.5*DXN(i2)*DYN(j2)*DZN(k2);
    return rho*m;
}

// row of (A + c0 rho M) u - B^T p at the face c of the cell i,j,k
double soil_stress::mech_row(int c, int i, int j, int k, const vget &g) const
{
    double sL[3], sR[3], row = 0.0, t1, t2, d1, d2, e1, e2;

    if(c==0)
    {
        normal_stress(i,j,k,g,sL[0],sL[1],sL[2]);
        normal_stress(i+1,j,k,g,sR[0],sR[1],sR[2]);
        const double area = DYN(j)*DZN(k);
        row = (sL[0]-sR[0])*area;

        t1 = edge_xz(i,j,k-1,g,d1,e1);
        t2 = edge_xz(i,j,k,g,d2,e2);
        row += (t1*d1 - t2*d2)*DYN(j);

        if(dim3)
        {
            t1 = edge_xy(i,j-1,k,g,d1,e1);
            t2 = edge_xy(i,j,k,g,d2,e2);
            row += (t1*d1 - t2*d2)*DZN(k);
        }
        row -= (P(g,i,j,k)-P(g,i+1,j,k))*area;
    }

    if(c==1)
    {
        normal_stress(i,j,k,g,sL[0],sL[1],sL[2]);
        normal_stress(i,j+1,k,g,sR[0],sR[1],sR[2]);
        const double area = DXN(i)*DZN(k);
        row = (sL[1]-sR[1])*area;

        t1 = edge_xy(i-1,j,k,g,d1,e1);
        t2 = edge_xy(i,j,k,g,d2,e2);
        row += (t1*e1 - t2*e2)*DZN(k);

        t1 = edge_yz(i,j,k-1,g,d1,e1);
        t2 = edge_yz(i,j,k,g,d2,e2);
        row += (t1*d1 - t2*d2)*DXN(i);

        row -= (P(g,i,j,k)-P(g,i,j+1,k))*area;
    }

    if(c==2)
    {
        normal_stress(i,j,k,g,sL[0],sL[1],sL[2]);
        normal_stress(i,j,k+1,g,sR[0],sR[1],sR[2]);
        const double area = DXN(i)*DYN(j);
        row = (sL[2]-sR[2])*area;

        t1 = edge_xz(i-1,j,k,g,d1,e1);
        t2 = edge_xz(i,j,k,g,d2,e2);
        row += (t1*e1 - t2*e2)*DYN(j);

        if(dim3)
        {
            t1 = edge_yz(i,j-1,k,g,d1,e1);
            t2 = edge_yz(i,j,k,g,d2,e2);
            row += (t1*e1 - t2*e2)*DXN(i);
        }
        row -= (P(g,i,j,k)-P(g,i,j,k+1))*area;
    }

    if(model==2)
    row += c0*face_mass(c,i,j,k)*U(g,c,i,j,k);

    return row;
}

// transmissibility of the face of the cell i,j,k to the neighbour in direction dir (0 x-, 1 x+, 2 y-,
// 3 y+, 4 z-, 5 z+); nt: the type of the neighbour, 0 for an impermeable face
double soil_stress::flow_trans(int i, int j, int k, int &dir) const
{
    int ii=i, jj=j, kk=k;
    double area, dist, K;

    if(dir<2)
    {
        ii += dir==0 ? -1 : 1;
        area = DYN(j)*DZN(k);
        K = Kx;
        dist = dir==0 ? DXP(i-1) : DXP(i);
        if(type(ii,jj,kk)==FREE) dist = 0.5*DXN(i);
    }
    else if(dir<4)
    {
        jj += dir==2 ? -1 : 1;
        area = DXN(i)*DZN(k);
        K = Ky;
        dist = dir==2 ? DYP(j-1) : DYP(j);
        if(type(ii,jj,kk)==FREE) dist = 0.5*DYN(j);
    }
    else
    {
        kk += dir==4 ? -1 : 1;
        area = DXN(i)*DYN(j);
        K = Kz;
        dist = dir==4 ? DZP(k-1) : DZP(k);
        if(type(ii,jj,kk)==FREE) dist = 0.5*DZN(k);
    }

    const int t = type(ii,jj,kk);
    dir = t;
    if(t!=SOIL && t!=FREE)
    return 0.0;
    return K*area/dist;
}

// row of -B u - (S + T/a0) p at the soil cell i,j,k
double soil_stress::flow_row(int i, int j, int k, const vget &g) const
{
    const double V = DXN(i)*DYN(j)*DZN(k);
    const double pc = P(g,i,j,k);
    double row = -strain_vol(i,j,k,g)*V - stor[id(i,j,k)]*V*pc;

    static const int nbx[6] = {-1,1,0,0,0,0};
    static const int nby[6] = {0,0,-1,1,0,0};
    static const int nbz[6] = {0,0,0,0,-1,1};

    for(int dir=0;dir<6;++dir)
    {
        if(!dim3 && (dir==2 || dir==3))
        continue;
        int nt = dir;
        const double T = flow_trans(i,j,k,nt);
        if(T==0.0)
        continue;
        const double pn = nt==SOIL ? P(g,i+nbx[dir],j+nby[dir],k+nbz[dir]) : 0.0;
        row -= T/a0*(pc - pn);
    }
    return row;
}

void soil_stress::rhs(std::vector<double> &bb) const
{
    const vget gn{Xn.data(),-1,0};
    const vget gnn{Xnn.data(),-1,0};

    auto ijk = [&](int q, int &i, int &j, int &k)
    {
        i = q/(NY*NZ) - 1;
        j = (q/NZ)%NY - 1;
        k = q%NZ - 1;
    };

    // mechanics: the bed pressure on the free faces, inertia of the old time levels
    for(int c=0;c<3;++c)
    for(int q : own[c])
    {
        int i,j,k;
        ijk(q,i,j,k);
        int i2=i, j2=j, k2=k;
        double area;
        if(c==0) {++i2; area = DYN(j)*DZN(k);}
        else if(c==1) {++j2; area = DXN(i)*DZN(k);}
        else {++k2; area = DXN(i)*DYN(j);}

        double f = 0.0;
        if(type(i,j,k)==FREE)
        f += pb[id(i,j,k)]*area;
        if(type(i2,j2,k2)==FREE)
        f -= pb[id(i2,j2,k2)]*area;

        if(model==2)
        f -= face_mass(c,i,j,k)*(c1*Xn[c*N+q] + c2*Xnn[c*N+q]);

        bb[c*N+q] = f;
    }

    // storage: the old time levels, the bed pressure
    static const int nbx[6] = {-1,1,0,0,0,0};
    static const int nby[6] = {0,0,-1,1,0,0};
    static const int nbz[6] = {0,0,0,0,-1,1};

    for(int q : own[3])
    {
        int i,j,k;
        ijk(q,i,j,k);
        const double V = DXN(i)*DYN(j)*DZN(k);

        double g = V*(stor[q]*(a1*Xn[3*N+q] + a2*Xnn[3*N+q]) + a1*strain_vol(i,j,k,gn) + a2*strain_vol(i,j,k,gnn));

        for(int dir=0;dir<6;++dir)
        {
            if(!dim3 && (dir==2 || dir==3))
            continue;
            int nt = dir;
            const double T = flow_trans(i,j,k,nt);
            if(nt==FREE)
            g -= T*pb[id(i+nbx[dir],j+nby[dir],k+nbz[dir])];
        }
        bb[3*N+q] = g/a0;
    }
}

void soil_stress::apply(ghostcell *pgc, std::vector<double> &x, std::vector<double> &out)
{
    halo_all(pgc,x);
    const vget g{x.data(),-1,0};

    for(int c=0;c<4;++c)
    for(int q : own[c])
    {
        const int i = q/(NY*NZ) - 1;
        const int j = (q/NZ)%NY - 1;
        const int k = q%NZ - 1;
        out[c*N+q] = c<3 ? mech_row(c,i,j,k,g) : flow_row(i,j,k,g);
    }
}

void soil_stress::diagonal()
{
    for(int c=0;c<3;++c)
    for(int q : own[c])
    {
        const int i = q/(NY*NZ) - 1;
        const int j = (q/NZ)%NY - 1;
        const int k = q%NZ - 1;
        D[c*N+q] = mech_row(c,i,j,k,vget{nullptr,c,q});
    }
}

void soil_stress::factorize()
{
    if(fac==nullptr)
    fac = new factor;
    factor &F = *fac;

    const bool region = F.ctype_f!=ctype;
    const bool new_u = region || !F.ok_u || (model==2 && fabs(c0/F.c0f-1.0)>0.2);
    const bool new_p = region || !F.ok_p || fabs(a0/F.a0f-1.0)>0.2;

    if(!new_u && !new_p)
    return;

    auto ijk = [&](int q, int &i, int &j, int &k)
    {
        i = q/(NY*NZ) - 1;
        j = (q/NZ)%NY - 1;
        k = q%NZ - 1;
    };

    if(region || F.lidx.empty())
    {
        F.lidx.assign(size_t(4)*N,-1);
        F.lu.clear();
        F.lp.clear();
        for(int c=0;c<3;++c)
        for(int q : own[c])
        {
            F.lidx[c*N+q] = int(F.lu.size());
            F.lu.push_back(c*N+q);
        }
        for(int q : own[3])
        {
            F.lidx[3*N+q] = int(F.lp.size());
            F.lp.push_back(3*N+q);
        }
        F.ctype_f = ctype;
    }

    // a ghost of a serial periodic direction (periodic 1) is an unknown of the rank: its number
    auto local = [&](int c, int i, int j, int k) -> int
    {
        if(p->periodic1==1 && p->nb1==-2 && p->nb4==-2)
        {
            if(i<0) i += knox;
            if(i>=knox) i -= knox;
        }
        if(dim3 && p->periodic2==1 && p->nb3==-2 && p->nb2==-2)
        {
            if(j<0) j += knoy;
            if(j>=knoy) j -= knoy;
        }
        if(p->periodic3==1 && p->nb5==-2 && p->nb6==-2)
        {
            if(k<0) k += knoz;
            if(k>=knoz) k -= knoz;
        }
        if(i<0 || i>=knox || j<0 || j>=knoy || k<0 || k>=knoz)
        return -1;
        return F.lidx[c*N+id(i,j,k)];
    };

    const int nu = int(F.lu.size());
    const int np = int(F.lp.size());
    F.ru.resize(nu); F.zu.resize(nu);
    F.rp.resize(np); F.zp.resize(np);

    if(new_u)
    {
        std::vector<Eigen::Triplet<double> > trip;
        trip.reserve(size_t(nu)*(dim3?40:12));
        const int jw = dim3 ? 1 : 0;

        for(int r=0;r<nu;++r)
        {
            const int c = F.lu[r]/N, q = F.lu[r]%N;
            int i,j,k;
            ijk(q,i,j,k);
            for(int cc=0;cc<3;++cc)
            for(int i2=i-1;i2<=i+1;++i2)
            for(int j2=j-jw;j2<=j+jw;++j2)
            for(int k2=k-1;k2<=k+1;++k2)
            {
                const int q2 = id(i2,j2,k2);
                const int col = local(cc,i2,j2,k2);
                if(col<0 || col<r)
                continue;
                const double val = mech_row(c,i,j,k,vget{nullptr,cc,q2});
                if(val!=0.0)
                trip.emplace_back(col,r,val);   // lower triangle
            }
        }
        Eigen::SparseMatrix<double> A(nu,nu);
        A.setFromTriplets(trip.begin(),trip.end());
        F.llt_u.compute(A);
        F.ok_u = F.llt_u.info()==Eigen::Success;
        F.c0f = c0;
    }

    if(new_p)
    {
        static const int nbx[6] = {-1,1,0,0,0,0};
        static const int nby[6] = {0,0,-1,1,0,0};
        static const int nbz[6] = {0,0,0,0,-1,1};
        const double beta = 1.0/(Ev*C33);

        std::vector<Eigen::Triplet<double> > trip;
        trip.reserve(size_t(np)*4);
        for(int r=0;r<np;++r)
        {
            const int q = F.lp[r]%N;
            int i,j,k;
            ijk(q,i,j,k);
            const double V = DXN(i)*DYN(j)*DZN(k);
            double d = V*(stor[q] + beta);
            for(int dir=0;dir<6;++dir)
            {
                if(!dim3 && (dir==2 || dir==3))
                continue;
                int nt = dir;
                const double T = flow_trans(i,j,k,nt)/a0;
                if(T==0.0)
                continue;
                d += T;
                if(nt==SOIL)
                {
                    const int col = local(3,i+nbx[dir],j+nby[dir],k+nbz[dir]);
                    if(col>r)
                    trip.emplace_back(col,r,-T);
                }
            }
            trip.emplace_back(r,r,d);
        }
        Eigen::SparseMatrix<double> Pm(np,np);
        Pm.setFromTriplets(trip.begin(),trip.end());
        F.llt_p.compute(Pm);
        F.ok_p = F.llt_p.info()==Eigen::Success;
        F.a0f = a0;
    }
}

void soil_stress::precon(const std::vector<double> &x, std::vector<double> &out) const
{
    const factor &F = *fac;

    if(F.ok_u)
    {
        factor &G = *fac;
        for(size_t r=0;r<F.lu.size();++r)
        G.ru[r] = x[F.lu[r]];
        G.zu = F.llt_u.solve(G.ru);
        for(size_t r=0;r<F.lu.size();++r)
        out[F.lu[r]] = G.zu[r];
    }
    else
    for(int c=0;c<3;++c)
    for(int q : own[c])
    out[c*N+q] = x[c*N+q]/D[c*N+q];

    if(F.ok_p)
    {
        factor &G = *fac;
        for(size_t r=0;r<F.lp.size();++r)
        G.rp[r] = x[F.lp[r]];
        G.zp = F.llt_p.solve(G.rp);
        for(size_t r=0;r<F.lp.size();++r)
        out[F.lp[r]] = G.zp[r];
    }
    else
    for(int q : own[3])
    out[3*N+q] = x[3*N+q]/D[3*N+q];
}

double soil_stress::dot(ghostcell *pgc, const std::vector<double> &x, const std::vector<double> &z) const
{
    double s = 0.0;
    for(int c=0;c<4;++c)
    for(int q : own[c])
    s += x[c*N+q]*z[c*N+q];
    return pgc->globalsum(s);
}

// the coupled system with MINRES (Paige and Saunders), block-diagonal SPD preconditioner
void soil_stress::solve(lexer *p, ghostcell *pgc)
{
    iter = 0;
    resid = 0.0;
    if(soilcells==0)
    return;

    const size_t n4 = size_t(4)*N;
    for(auto *vec : {&b,&r1,&r2,&y,&v,&w,&w1,&w2})
    vec->assign(n4,0.0);

    // preconditioner: diag(A + c0 rho M) and diag(S + T/a0) + diag(B D_A^-1 B^T)
    diagonal();
    for(int c=0;c<3;++c)
    if(c!=1 || dim3)
    halo(pgc,D.data()+size_t(c)*N);

    for(int q : own[3])
    {
        const int i = q/(NY*NZ) - 1;
        const int j = (q/NZ)%NY - 1;
        const int k = q%NZ - 1;
        double d = -flow_row(i,j,k,vget{nullptr,3,q});

        auto add = [&](int c, int ii, int jj, int kk, double area)
        {
            const int qq = id(ii,jj,kk);
            if(int(kind[c][qq])==UNK)
            d += area*area/D[c*N+qq];
        };
        add(0,i,j,k,DYN(j)*DZN(k));
        add(0,i-1,j,k,DYN(j)*DZN(k));
        add(2,i,j,k,DXN(i)*DYN(j));
        add(2,i,j,k-1,DXN(i)*DYN(j));
        if(dim3)
        {
        add(1,i,j,k,DXN(i)*DZN(k));
        add(1,i,j-1,k,DXN(i)*DZN(k));
        }
        D[3*N+q] = d;
    }

    factorize();

    rhs(b);

    auto axpy = [&](std::vector<double> &out, double s, const std::vector<double> &x)
    {
        for(int c=0;c<4;++c)
        for(int q : own[c])
        out[c*N+q] += s*x[c*N+q];
    };

    // r1 = b - A x, y = M^-1 r1
    apply(pgc,X,y);
    for(int c=0;c<4;++c)
    for(int q : own[c])
    r1[c*N+q] = b[c*N+q] - y[c*N+q];

    precon(b,y);
    const double bnorm = sqrt(max(dot(pgc,b,y),0.0));

    precon(r1,y);
    double beta1 = sqrt(max(dot(pgc,r1,y),0.0));

    const double tol = p->S214_tol;
    const int maxit = p->S214_n;

    if(bnorm<1.0e-300)
    {
        for(int c=0;c<4;++c)
        for(int q : own[c])
        X[c*N+q] = 0.0;
        halo_all(pgc,X);
        return;
    }

    resid = beta1/bnorm;
    if(beta1<=tol*bnorm)
    {
        halo_all(pgc,X);
        return;
    }

    double oldb=0.0, beta=beta1, dbar=0.0, epsln=0.0, phibar=beta1, cs=-1.0, sn=0.0;
    r2 = r1;

    for(int itn=1; itn<=maxit; ++itn)
    {
        const double s = 1.0/beta;
        for(int c=0;c<4;++c)
        for(int q : own[c])
        v[c*N+q] = s*y[c*N+q];

        apply(pgc,v,y);

        if(itn>=2)
        axpy(y,-beta/oldb,r1);

        const double alfa = dot(pgc,v,y);
        axpy(y,-alfa/beta,r2);

        r1.swap(r2);
        for(int c=0;c<4;++c)
        for(int q : own[c])
        r2[c*N+q] = y[c*N+q];

        precon(r2,y);
        oldb = beta;
        beta = sqrt(max(dot(pgc,r2,y),0.0));

        const double oldeps = epsln;
        const double delta = cs*dbar + sn*alfa;
        const double gbar = sn*dbar - cs*alfa;
        epsln = sn*beta;
        dbar = -cs*beta;

        const double gamma = max(sqrt(gbar*gbar + beta*beta),1.0e-300);
        cs = gbar/gamma;
        sn = beta/gamma;
        const double phi = cs*phibar;
        phibar = sn*phibar;

        w1.swap(w2);
        w2.swap(w);
        for(int c=0;c<4;++c)
        for(int q : own[c])
        {
            const int n = c*N+q;
            w[n] = (v[n] - oldeps*w1[n] - delta*w2[n])/gamma;
            X[n] += phi*w[n];
        }

        iter = itn;
        resid = phibar/bnorm;

        if(phibar<=tol*bnorm || beta<1.0e-300)
        break;
    }

    iter_max = max(iter_max,iter);
    halo_all(pgc,X);
}
