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

#include"fnpf_ice.h"
#include"lexer.h"
#include"ghostcell.h"
#include"ice_contact.h"
#include<iomanip>
#include<algorithm>
#include<random>

// Ice breaking (A 390): a floe splits in two along a straight cut, at most once per floe and check.
//
// Flexural failure (A 390 1,3): along A 395 directions the net vertical load of the floe
//     q = phi*A*( p_lid - rho_i*g*h - rho_i*h*a_z(x,y) ),   a_z = a_G,z + (alpha x r)_z   (d'Alembert)
// is binned into a line load Q(s) (A 395 bins). The floe is a free-free beam on an elastic foundation
// along s (A 397, E > 0):
//     (B w'')'' + k w = Q,   B = D*b(s),  D = E h^3/(12(1-nu^2)),  k = rho_w*g*b(s),  b(s) = chord width
// w is the elastic deflection relative to the rigid-body motion: the rigid lid gives the loads of a rigid
// floe, where the floe bends the hydrostatic pressure relaxes by rho_w*g*w. Bending moment M = B w'',
// stress sigma = 6|M|/(h^2 b). Floes short compared with the flexural length (D/(rho_w g))^(1/4) recover
// rigid-floe statics, long floes break at a distance of the order of the flexural length.
// E <= 0: rigid-floe statics, M(s) = sum_{s_j > s} Q_j (s_j - s).
// Breaks when sigma >= sigma_f of the floe on the most stressed cut. sigma_f with Weibull scatter and
// size effect (A 399): sigma_f,i = sigma_f*(A_ref/A_i)^(1/m)*(-ln U)^(1/m)/Gamma(1+1/m), U uniform, new
// values for both pieces after a split.
// q sums to zero over a floe at rest and in rigid-body equilibrium, so there is no spurious moment.
//
// Contact splitting (A 390 2,3, 3D): a floe splits along the contact normal through the contact point
// when the pair normal force exceeds F_s = C*K_IC*h*sqrt(D) (A 392), D = 2*sqrt(A/pi); LEFM scaling of
// in-plane splitting (Bhat 1988, Lu et al. 2015), C to be calibrated. The force is the pair normal
// force of the non-smooth contact averaged over t_c (A 396), see contact_average: impact impulses are
// resolved within one or two steps, impulse/dt alone depends on dt and on the timing of the impact.
//
// Pieces with a caliper width below D_min (A 393) are not created. Mass, momentum and angular
// momentum are conserved: pieces keep the orientation and the rigid-body velocity field of the parent.

namespace
{
    double polygon_centroid(const vector<double> &x, const vector<double> &y, double &cx, double &cy)
    {
        const int n=int(x.size());
        double A=0.0;
        cx=cy=0.0;
        for(int q=0;q<n;++q)
        {
        const int q2=(q+1)%n;
        const double cr = x[q]*y[q2]-x[q2]*y[q];
        A  += 0.5*cr;
        cx += (x[q]+x[q2])*cr;
        cy += (y[q]+y[q2])*cr;
        }
        if(fabs(A)>0.0)
        {
        cx/=(6.0*A);
        cy/=(6.0*A);
        }
        return A;
    }

    double caliper(const vector<double> &x, const vector<double> &y)
    {
        const int n=int(x.size());
        double w=1.0e20;
        for(int q=0;q<n;++q)
        {
        const int q2=(q+1)%n;
        const double ex=x[q2]-x[q], ey=y[q2]-y[q];
        const double len=sqrt(ex*ex+ey*ey);
        if(len<=0.0)
        continue;
        double dm=0.0;
        for(int r=0;r<n;++r)
        dm = MAX(dm, fabs(-ey*(x[r]-x[q]) + ex*(y[r]-y[q]))/len);
        w = MIN(w,dm);
        }
        return w>1.0e19 ? 0.0 : w;
    }
}

void fnpf_ice::clip_halfplane(const vector<double> &x, const vector<double> &y, double nx, double ny, double s, double sgn,
                              vector<double> &ox, vector<double> &oy)
{
    // keep sgn*(n.x - s) <= 0 (Sutherland-Hodgman, one plane)
    ox.clear();
    oy.clear();
    const int n=int(x.size());

    for(int q=0;q<n;++q)
    {
        const int q2=(q+1)%n;
        const double d0 = sgn*(nx*x[q]  + ny*y[q]  - s);
        const double d1 = sgn*(nx*x[q2] + ny*y[q2] - s);

        if(d0<=0.0)
        {
        ox.push_back(x[q]);
        oy.push_back(y[q]);
        }

        if((d0<0.0 && d1>0.0) || (d0>0.0 && d1<0.0))
        {
        const double t = d0/(d0-d1);
        ox.push_back(x[q] + t*(x[q2]-x[q]));
        oy.push_back(y[q] + t*(y[q2]-y[q]));
        }
    }
}

double fnpf_ice::chord(const vector<double> &x, const vector<double> &y, double nx, double ny, double s)
{
    // length of the intersection of the line n.x = s with a convex polygon
    const int n=int(x.size());
    double tmin=1.0e20, tmax=-1.0e20;
    const double tx=-ny, ty=nx;

    for(int q=0;q<n;++q)
    {
        const int q2=(q+1)%n;
        const double d0 = nx*x[q]  + ny*y[q]  - s;
        const double d1 = nx*x[q2] + ny*y[q2] - s;

        if((d0<=0.0 && d1>0.0) || (d0>0.0 && d1<=0.0))
        {
        const double t = d0/(d0-d1);
        const double xi = x[q] + t*(x[q2]-x[q]);
        const double yi = y[q] + t*(y[q2]-y[q]);
        const double u = tx*xi + ty*yi;
        tmin = MIN(tmin,u);
        tmax = MAX(tmax,u);
        }
    }
    return tmax>tmin ? tmax-tmin : 0.0;
}

void fnpf_ice::strength(fnpf_ice_floe &fl)
{
    fl.sigf = sigf;
    
    if(wm>0.0 && fl.area>0.0)
    {
    uniform_real_distribution<double> U(1.0e-12,1.0);
    const double u = U(rng);
    fl.sigf = sigf*pow(wA/fl.area,1.0/wm)*pow(-log(u),1.0/wm)/tgamma(1.0+1.0/wm);
    }
}

int fnpf_ice::beam(int N, double ds, const double *B, const double *k, const double *Q, double *M)
{
    // free-free beam on an elastic foundation, finite volumes on N bins of width ds:
    //   V_{i+1/2} - V_{i-1/2} + k_i ds w_i = Q_i,   V_{i+1/2} = (M_{i+1}-M_i)/ds,
    //   M_i = B_i (w_{i-1} - 2 w_i + w_{i+1})/ds^2 (i = 1..N-2),  M_0 = M_{N-1} = 0 (free ends)
    // A = sum_r B_r/ds^3 s_r s_r^T + diag(k ds), symmetric positive definite, bandwidth 2
    if(N<3)
    return -1;
    
    vector<double> A(N*5,0.0);   // band storage: A[i*5 + (j-i+2)]
    auto a = [&](int i, int j) -> double& { return A[i*5 + (j-i+2)]; };
    
    for(int i=0;i<N;++i)
    a(i,i) += k[i]*ds;
    
    for(int r=1;r<N-1;++r)
    {
        const double c = B[r]/(ds*ds*ds);
        const int id[3] = {r-1,r,r+1};
        const double sv[3] = {1.0,-2.0,1.0};
        for(int p=0;p<3;++p)
        for(int q=0;q<3;++q)
        a(id[p],id[q]) += c*sv[p]*sv[q];
    }
    
    // banded Gaussian elimination, no pivoting (SPD)
    vector<double> w(Q,Q+N);
    for(int c=0;c<N;++c)
    {
        const double piv = a(c,c);
        if(fabs(piv)<1.0e-300)
        return -1;
        for(int i=c+1;i<=MIN(c+2,N-1);++i)
        {
            const double f = a(i,c)/piv;
            if(f==0.0)
            continue;
            for(int j=c;j<=MIN(c+2,N-1);++j)
            a(i,j) -= f*a(c,j);
            w[i] -= f*w[c];
        }
    }
    for(int c=N-1;c>=0;--c)
    {
        double sum=w[c];
        for(int j=c+1;j<=MIN(c+2,N-1);++j)
        sum -= a(c,j)*w[j];
        w[c] = sum/a(c,c);
    }
    
    M[0]=M[N-1]=0.0;
    for(int i=1;i<N-1;++i)
    M[i] = B[i]*(w[i-1]-2.0*w[i]+w[i+1])/(ds*ds);
    
    return 0;
}

void fnpf_ice::breaking_moments(lexer *p)
{
    // all ranks: line loads of the local footprint cells in bins along each direction, summed over ranks
    const int nbin = noff;
    const int nq = ndir*nbin;
    const size_t nf = floe.size();
    Mcut.assign(nq*nf + 2*ndir*nf,0.0);
    
    if(!(breakflag&1))
    return;
    
    // bin origin and width per floe and direction, from the planform extent (identical on all ranks)
    double *geo = &Mcut[nq*nf];
    
    for(size_t f=0; f<nf; ++f)
    {
        const fnpf_ice_floe &fl = floe[f];
        if(fl.type!=0)
        continue;
        
        for(int k=0; k<ndir; ++k)
        {
            const double th = PI*double(k)/double(ndir);
            const double nx = cos(th), ny = sin(th);
            double smin=1.0e20, smax=-1.0e20;
            
            for(size_t q=0; q<fl.wx.size(); ++q)
            {
            const double d = nx*(fl.wx[q]-fl.x[0]) + ny*(fl.wy[q]-fl.x[1]);
            smin = MIN(smin,d);
            smax = MAX(smax,d);
            }
            
            geo[2*(ndir*f+k)]   = smin;
            geo[2*(ndir*f+k)+1] = MAX(smax-smin,1.0e-12)/double(nbin);
        }
    }
    
    for(auto &e : cell)
    {
        const fnpf_ice_floe &fl = floe[e.f];
        
        if(fl.type!=0)
        continue;
        
        const double rx = e.xc - fl.x[0];
        const double ry = e.yc - fl.x[1];
        const double az = fl.acc[2] + fl.alp[0]*ry - fl.alp[1]*rx;
        const double qA = e.phi*e.area*(e.pl - fl.pw - fl.rho*fl.h*az);
        
        double *Q = &Mcut[nq*e.f];
        
        for(int k=0; k<ndir; ++k)
        {
            const double th = PI*double(k)/double(ndir);
            const double d = cos(th)*rx + sin(th)*ry;
            const double s0 = geo[2*(ndir*e.f+k)];
            const double ds = geo[2*(ndir*e.f+k)+1];
            const int ib = MAX(0, MIN(nbin-1, int(floor((d-s0)/ds))));
            Q[k*nbin+ib] += qA;
        }
    }
    
    if(p->mpi_size>1)
    MPI_Allreduce(MPI_IN_PLACE, Mcut.data(), int(nq*nf), MPI_DOUBLE, MPI_SUM, comm);
}

void fnpf_ice::breaking_decide(lexer *p)
{
    // rank 0: pick at most one cut per floe
    splits.clear();

    const size_t nf = floe.size();

    vector<split> best(nf);
    vector<double> ratio(nf,0.0);

    auto body_cut = [&](const fnpf_ice_floe &fl, double nx, double ny, double s, split &sp)
    {
        // world horizontal cut through the centre of mass frame -> body frame
        double bxn = fl.R[0][0]*nx + fl.R[1][0]*ny;
        double byn = fl.R[0][1]*nx + fl.R[1][1]*ny;
        const double len = sqrt(bxn*bxn + byn*byn);
        if(len<1.0e-12)
        return 0;
        sp.nbx = bxn/len;
        sp.nby = byn/len;
        sp.s = s/len;

        // both pieces at least D_min wide
        vector<double> ax,ay,bx2,by2;
        clip_halfplane(fl.bx,fl.by,sp.nbx,sp.nby,sp.s, 1.0,ax,ay);
        clip_halfplane(fl.bx,fl.by,sp.nbx,sp.nby,sp.s,-1.0,bx2,by2);
        if(ax.size()<3 || bx2.size()<3)
        return 0;
        if(is2D)
        {
        const double la = *max_element(ax.begin(),ax.end()) - *min_element(ax.begin(),ax.end());
        const double lb = *max_element(bx2.begin(),bx2.end()) - *min_element(bx2.begin(),bx2.end());
        return (la>=Dmin && lb>=Dmin) ? 1 : 0;
        }
        return (caliper(ax,ay)>=Dmin && caliper(bx2,by2)>=Dmin) ? 1 : 0;
    };

    // flexural: beam on elastic foundation (or rigid statics) along each direction
    if(breakflag&1)
    {
    const int nbin = noff;
    const int nq = ndir*nbin;
    const double *geo = &Mcut[nq*nf];
    vector<double> Bv(nbin),kv(nbin),Qv(nbin),Mv(nbin),bw(nbin);
    
    for(size_t f=0; f<nf; ++f)
    {
        const fnpf_ice_floe &fl = floe[f];
        if(fl.type!=0 || fl.Awet<=0.0)
        continue;
        
        // planform at the centre of mass level, relative to it
        const int nv = int(fl.bx.size());
        vector<double> px(nv),py(nv);
        for(int q=0;q<nv;++q)
        {
        px[q] = fl.R[0][0]*fl.bx[q] + fl.R[0][1]*fl.by[q];
        py[q] = fl.R[1][0]*fl.bx[q] + fl.R[1][1]*fl.by[q];
        }
        
        const double Dp = (Emod>0.0) ? Emod*fl.h*fl.h*fl.h/(12.0*(1.0-nu*nu)) : 0.0;
        
        for(int k=0; k<ndir; ++k)
        {
            const double th = PI*double(k)/double(ndir);
            const double nx = cos(th), ny = sin(th);
            const double s0 = geo[2*(ndir*f+k)];
            const double ds = geo[2*(ndir*f+k)+1];
            
            double bmax=0.0;
            for(int i=0;i<nbin;++i)
            {
            const double si = s0 + (i+0.5)*ds;
            bw[i] = is2D ? fl.width2D : chord(px,py,nx,ny,si);
            bmax = MAX(bmax,bw[i]);
            }
            if(bmax<=0.0)
            continue;
            
            for(int i=0;i<nbin;++i)
            {
            bw[i] = MAX(bw[i],1.0e-3*bmax);
            Qv[i] = Mcut[nq*f + k*nbin + i];
            Bv[i] = Dp*bw[i];
            kv[i] = rhow*g*bw[i];
            }
            
            if(Emod>0.0)
            {
            if(beam(nbin,ds,Bv.data(),kv.data(),Qv.data(),Mv.data())!=0)
            continue;
            }
            else
            {
            // rigid statics at the bin centres
            for(int i=0;i<nbin;++i)
            {
            double m=0.0;
            for(int j=i+1;j<nbin;++j)
            m += Qv[j]*(j-i)*ds;
            Mv[i]=m;
            }
            }
            
            for(int i=1;i<nbin-1;++i)
            {
                const double sig = 6.0*fabs(Mv[i])/(fl.h*fl.h*bw[i]);
                const double r = sig/MAX(fl.sigf,1.0e-20);
                
                if(r>=1.0 && r>ratio[f])
                {
                split sp;
                const double si = s0 + (i+0.5)*ds;
                if(body_cut(fl,nx,ny,si,sp))
                {
                sp.f = int(f);
                sp.mech = 1;
                sp.val = sig;
                sp.lim = fl.sigf;
                best[f] = sp;
                ratio[f] = r;
                }
                }
            }
        }
    }
    }
    
    // contact splitting, on the time-averaged pair normal force (contact_average)
    if((breakflag&2) && !is2D)
    for(const auto &kv : cforce)
    for(int side=0; side<2; ++side)
    {
        const int fid = side==0 ? kv.first.first : kv.first.second;
        if(fid<0 || fid>=int(nf))
        continue;
        
        const size_t f = size_t(fid);
        const fnpf_ice_floe &fl = floe[f];
        if(fl.type!=0)
        continue;
        
        const cavg &ca = kv.second;
        const double D = 2.0*sqrt(fl.area/PI);
        const double Fs = Csplit*KIC*fl.h*sqrt(D);
        const double r = ca.F/MAX(Fs,1.0e-20);
        
        if(r>=1.0 && r>ratio[f])
        {
        // cut along the contact normal through the contact point
        const double nx = -ca.ny, ny = ca.nx;
        const double s = nx*(ca.px-fl.x[0]) + ny*(ca.py-fl.x[1]);
        split sp;
        if(body_cut(fl,nx,ny,s,sp))
        {
        sp.f = int(f);
        sp.mech = 2;
        sp.val = ca.F;
        sp.lim = Fs;
        best[f] = sp;
        ratio[f] = r;
        }
        }
    }
    
    for(size_t f=0; f<nf; ++f)
    if(ratio[f]>=1.0)
    splits.push_back(best[f]);
}

void fnpf_ice::breaking_apply(lexer *p, ghostcell *pgc)
{
    int ns = (p->mpirank==0) ? int(splits.size()) : 0;

    if(p->mpi_size>1)
    {
    MPI_Bcast(&ns, 1, MPI_INT, 0, pgc->mpi_comm);

    vector<double> buf(7*ns);
    if(p->mpirank==0)
    for(int n=0;n<ns;++n)
    {
    const split &sp = splits[n];
    double *b=&buf[7*n];
    b[0]=sp.f; b[1]=sp.mech; b[2]=sp.nbx; b[3]=sp.nby; b[4]=sp.s; b[5]=sp.val; b[6]=sp.lim;
    }

    if(ns>0)
    MPI_Bcast(buf.data(), 7*ns, MPI_DOUBLE, 0, pgc->mpi_comm);

    if(p->mpirank>0)
    {
    splits.resize(ns);
    for(int n=0;n<ns;++n)
    {
    const double *b=&buf[7*n];
    splits[n].f=int(b[0]); splits[n].mech=int(b[1]); splits[n].nbx=b[2]; splits[n].nby=b[3];
    splits[n].s=b[4]; splits[n].val=b[5]; splits[n].lim=b[6];
    }
    }
    }

    const double t = p->simtime + p->dt;

    for(int n=0;n<ns;++n)
    {
        const split &sp = splits[n];
        const int nid = split_floe(p,size_t(sp.f),sp.nbx,sp.nby,sp.s,sp.mech);

        if(nid>=0)
        {
        ++nbreak;
        
        // the load that broke the floe is released
        const int pid = floe[sp.f].id;
        for(auto it=cforce.begin(); it!=cforce.end();)
        {
            if(it->first.first==pid || it->first.second==pid)
            it = cforce.erase(it);
            else
            ++it;
        }

        if(p->mpirank==0 && breakout.is_open())
        breakout<<setprecision(9)<<t<<" "<<floe[sp.f].id<<" "<<nid<<" "<<sp.mech<<" "<<sp.val<<" "<<sp.lim<<" "
                <<floe[sp.f].area<<" "<<floe.back().area<<endl;
        }
    }

    if(ns>0 && p->mpirank==0)
    cout<<"FNPF ice: "<<ns<<" floe(s) broken at t = "<<t<<", "<<nfloe<<" floes"<<endl;

    splits.clear();
}

int fnpf_ice::split_floe(lexer *p, size_t f, double nbx, double nby, double s, int mech)
{
    fnpf_ice_floe parent = floe[f];

    vector<double> ax,ay,bx,by;
    clip_halfplane(parent.bx,parent.by,nbx,nby,s, 1.0,ax,ay);
    clip_halfplane(parent.bx,parent.by,nbx,nby,s,-1.0,bx,by);

    if(ax.size()<3 || bx.size()<3)
    return -1;

    quat_to_matrix(parent.q,parent.R);

    auto make_piece = [&](const vector<double> &px, const vector<double> &py, fnpf_ice_floe &pc)
    {
        double cx,cy;
        polygon_centroid(px,py,cx,cy);

        pc = parent;
        pc.bx = px;
        pc.by = py;
        geometry(pc);   // re-centres on the piece centroid

        // world position and rigid-body velocity of the piece centre
        const double r[3] = {parent.R[0][0]*cx + parent.R[0][1]*cy,
                             parent.R[1][0]*cx + parent.R[1][1]*cy,
                             parent.R[2][0]*cx + parent.R[2][1]*cy};
        for(int a=0;a<3;++a)
        pc.x[a] = parent.x[a] + r[a];

        pc.v[0] = parent.v[0] + parent.w[1]*r[2] - parent.w[2]*r[1];
        pc.v[1] = parent.v[1] + parent.w[2]*r[0] - parent.w[0]*r[2];
        pc.v[2] = parent.v[2] + parent.w[0]*r[1] - parent.w[1]*r[0];

        pc.Fc[0]=pc.Fc[1]=0.0;
        pc.ncontact=0;
        kinematics(pc);
    };

    fnpf_ice_floe pa,pb;
    make_piece(ax,ay,pa);
    make_piece(bx,by,pb);
    
    // new strengths for both pieces (same random sequence on all ranks)
    strength(pa);
    strength(pb);

    pa.id = parent.id;
    pb.id = int(floe.size());

    floe[f] = pa;
    floe.push_back(pb);

    Yn.resize(floe.size());
    D1.resize(floe.size());
    D2.resize(floe.size());
    D3.resize(floe.size());
    ++nfloe;

    (void)mech;
    return pb.id;
}

void fnpf_ice::contact_average(lexer *p)
{
    // Rank 0, every step. The non-smooth contact gives impulses, an impact is resolved within one or two
    // steps and impulse/dt depends on dt and on where the impact falls within the step. For splitting
    // the pair normal force is averaged over t_c (A 396) with an exponential filter:
    //     F_avg <- F_avg + a*(F - F_avg),  a = min(1, dt/t_c)
    // An impact with impulse P gives F_avg ~ P/t_c, a sustained load gives F_avg = F.
    const double a = (tcon>0.0) ? MIN(1.0, p->dt/tcon) : 1.0;
    
    for(auto &kv : cforce)
    kv.second.F *= (1.0-a);
    
    for(const auto &rc : pcontact->records())
    {
        if(rc.a<0 || rc.a>=int(cmap.size()))
        continue;
        
        const int ida = floe[cmap[rc.a]].id;
        const int idb = (rc.b>=0 && rc.b<int(cmap.size())) ? floe[cmap[rc.b]].id : rc.b;
        
        cavg &ca = cforce[make_pair(ida,idb)];
        ca.F  += a*rc.Fn;
        ca.px = rc.px;
        ca.py = rc.py;
        ca.nx = rc.nx;
        ca.ny = rc.ny;
    }
    
    for(auto it=cforce.begin(); it!=cforce.end();)
    {
        if(it->second.F < 1.0e-6)
        it = cforce.erase(it);
        else
        ++it;
    }
}
