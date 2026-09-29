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

// Ice breaking (A 390, bits: 1 flexural, 2 contact splitting, 4 spalling): floes break along straight
// cuts, so all pieces stay convex. Per floe and check one mechanism acts, the most critical one; it may
// make several cuts, which are applied one after another to the piece each cut crosses.
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
// Several cracks per check: flexural up to A 400 parallel cracks along the governing direction, at the local
// stress maxima above the strength, at least D_min apart; splitting A 401 radial cracks through the contact
// point, fanned about the contact normal at -90 + 180 (k+1)/(A 401 + 1) deg: the floe breaks into wedges
// with their tips at the contact.
//
// Spalling (A 390 4, needs crushing A 398): at a crushing contact the crushed depth grows; when it reaches
// L_sp (A 402, <0: h) the tip in front of the cut parallel to the contact face, L_sp behind the face of the
// other body, breaks off. The chip is a new floe if it is at least D_min wide, otherwise it is cleared as
// rubble: it leaves the contact and stays in place as a type 3 chip whose lid pressure and surface
// depression fade out linearly over t_r (A 402), its planform area is logged. Only local chips
// spall: the chord of the cut must not exceed C_loc (A 402) times the chip depth, a straight cut across a
// wide floe face is not a spall; the floe keeps crushing then. Wedge tips after radial splitting, floe
// corners and small floes spall, which gives the crushing-spalling load cycles at structures.
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
        // loads of the first RK stage, i.e. of the accepted state t^n: the intermediate stages carry the
        // stiff lid modes with O(1) errors once omega_lid*dt ~ 1, the accepted state does not
        const double az = fl.acc1[2] + fl.alp1[0]*ry - fl.alp1[1]*rx;
        const double qA = e.phi*e.area*(e.pl1 - fl.pw - fl.rho*fl.h*az);
        
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

void fnpf_ice::world_polygon(fnpf_ice_floe &fl, vector<double> &px, vector<double> &py)
{
    // planform at the level of the centre of mass, world coordinates (as in the contact solve)
    quat_to_matrix(fl.q,fl.R);
    const int nv = int(fl.bx.size());
    px.resize(nv);
    py.resize(nv);
    for(int q=0; q<nv; ++q)
    {
    px[q] = fl.x[0] + fl.R[0][0]*fl.bx[q] + fl.R[0][1]*fl.by[q];
    py[q] = fl.x[1] + fl.R[1][0]*fl.bx[q] + fl.R[1][1]*fl.by[q];
    }
}

double fnpf_ice::piece_width(const vector<double> &x, const vector<double> &y) const
{
    if(x.size()<3)
    return 0.0;
    if(is2D)
    return *max_element(x.begin(),x.end()) - *min_element(x.begin(),x.end());
    return caliper(x,y);
}

int fnpf_ice::cut_valid(fnpf_ice_floe &fl, const split &sp, double &nbx, double &nby, double &sb, int &rubble)
{
    // world cut line n.x = s -> body frame of this floe, nb.xb = sb; checks the piece sizes.
    // Cracks (mech 1,2): both pieces at least D_min wide. Spall (mech 3): the remainder (n.x <= s)
    // at least D_min wide, the chip (n.x >= s) becomes a floe if it is at least D_min wide, else rubble.
    rubble = 0;
    quat_to_matrix(fl.q,fl.R);
    const double bxn = fl.R[0][0]*sp.nx + fl.R[1][0]*sp.ny;
    const double byn = fl.R[0][1]*sp.nx + fl.R[1][1]*sp.ny;
    const double len = sqrt(bxn*bxn + byn*byn);
    if(len<1.0e-12)
    return 0;
    nbx = bxn/len;
    nby = byn/len;
    sb  = (sp.s - sp.nx*fl.x[0] - sp.ny*fl.x[1])/len;
    
    vector<double> ax,ay,bx2,by2;
    clip_halfplane(fl.bx,fl.by,nbx,nby,sb, 1.0,ax,ay);
    clip_halfplane(fl.bx,fl.by,nbx,nby,sb,-1.0,bx2,by2);
    if(ax.size()<3 || bx2.size()<3)
    return 0;
    
    const double wa = piece_width(ax,ay);
    const double wb = piece_width(bx2,by2);
    
    if(sp.mech==3)
    {
    if(wa<Dmin)
    return 0;
    rubble = (wb<Dmin) ? 1 : 0;
    return 1;
    }
    
    return (wa>=Dmin && wb>=Dmin) ? 1 : 0;
}

void fnpf_ice::breaking_decide(lexer *p)
{
    // rank 0: the cuts of this check, per floe one mechanism (the most critical one):
    //   flexural   up to A 400 parallel cracks along the governing direction
    //   splitting  A 401 radial cracks through the contact point, fanned about the contact normal
    //   spalling   the crushed tip at a crushing contact, if no crack is due
    splits.clear();

    const size_t nf = floe.size();

    vector<vector<split>> cuts(nf);
    vector<double> ratio(nf,0.0);
    
    map<int,size_t> index;
    for(size_t f=0; f<nf; ++f)
    index[floe[f].id] = f;
    
    double nbx,nby,sb;
    int rb;
    
    auto line = [&](int f, int mech, double nx, double ny, double s, double val, double lim)
    {
        split sp;
        sp.f=f; sp.mech=mech; sp.nx=nx; sp.ny=ny; sp.s=s; sp.val=val; sp.lim=lim;
        return sp;
    };

    // flexural: beam on elastic foundation (or rigid statics) along each direction
    if(breakflag&1)
    {
    const int nbin = noff;
    const int nq = ndir*nbin;
    const double *geo = &Mcut[nq*nf];
    vector<double> Bv(nbin),kv(nbin),Qv(nbin),Mv(nbin),bw(nbin);
    vector<double> sig(ndir*nbin);
    
    for(size_t f=0; f<nf; ++f)
    {
        fnpf_ice_floe &fl = floe[f];
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
        int kbest=-1;
        fl.sigmax = 0.0;
        split best;
        
        for(int k=0; k<ndir; ++k)
        {
            const double th = PI*double(k)/double(ndir);
            const double nx = cos(th), ny = sin(th);
            const double s0 = geo[2*(ndir*f+k)];
            const double ds = geo[2*(ndir*f+k)+1];
            
            for(int i=0;i<nbin;++i)
            sig[k*nbin+i] = 0.0;
            
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
                sig[k*nbin+i] = 6.0*fabs(Mv[i])/(fl.h*fl.h*bw[i]);
                fl.sigmax = MAX(fl.sigmax, sig[k*nbin+i]);
                const double r = sig[k*nbin+i]/MAX(fl.sigf,1.0e-20);
                
                if(r>=1.0 && r>ratio[f])
                {
                const double si = s0 + (i+0.5)*ds;
                split sp = line(int(f),1,nx,ny,si + nx*fl.x[0] + ny*fl.x[1],sig[k*nbin+i],fl.sigf);
                if(cut_valid(fl,sp,nbx,nby,sb,rb))
                {
                best = sp;
                kbest = k;
                ratio[f] = r;
                }
                }
            }
        }
        
        if(kbest<0)
        continue;
        
        cuts[f].assign(1,best);
        
        // further cracks along the governing direction: local stress maxima above the strength,
        // strongest first, at least D_min apart
        if(ncrack>1)
        {
            const int k = kbest;
            const double th = PI*double(k)/double(ndir);
            const double nx = cos(th), ny = sin(th);
            const double s0 = geo[2*(ndir*f+k)];
            const double ds = geo[2*(ndir*f+k)+1];
            const double *sk = &sig[k*nbin];
            
            vector<int> peak;
            for(int i=1;i<nbin-1;++i)
            if(sk[i]>=fl.sigf && sk[i]>=sk[i-1] && sk[i]>=sk[i+1])
            peak.push_back(i);
            sort(peak.begin(),peak.end(),[&](int a, int b){return sk[a]>sk[b];});
            
            for(int i : peak)
            {
                if(int(cuts[f].size())>=ncrack)
                break;
                
                split sp = line(int(f),1,nx,ny,s0 + (i+0.5)*ds + nx*fl.x[0] + ny*fl.x[1],sk[i],fl.sigf);
                int ok=1;
                for(const auto &c : cuts[f])
                if(fabs(c.s-sp.s)<Dmin)
                ok=0;
                if(ok && cut_valid(fl,sp,nbx,nby,sb,rb))
                cuts[f].push_back(sp);
            }
        }
    }
    }
    
    // contact splitting, on the time-averaged pair normal force (contact_average):
    // radial cracks through the contact point, at the angles -90 + 180 (k+1)/(A 401 + 1) deg to the
    // contact normal, the one closest to the normal first (A 401 1: along the normal)
    if((breakflag&2) && !is2D)
    for(const auto &kv : cforce)
    for(int side=0; side<2; ++side)
    {
        const int fid = side==0 ? kv.first.first : kv.first.second;
        auto it = index.find(fid);
        if(fid<0 || it==index.end())
        continue;
        
        const size_t f = it->second;
        fnpf_ice_floe &fl = floe[f];
        if(fl.type!=0)
        continue;
        
        const cavg &ca = kv.second;
        const double D = 2.0*sqrt(fl.area/PI);
        const double Fs = Csplit*KIC*fl.h*sqrt(D);
        const double r = ca.F/MAX(Fs,1.0e-20);
        
        if(r>=1.0 && r>ratio[f])
        {
            vector<double> ang(nradial);
            for(int k=0;k<nradial;++k)
            ang[k] = (nradial==1) ? 0.0 : PI*(-0.5 + double(k+1)/double(nradial+1));
            sort(ang.begin(),ang.end(),[](double a, double b){return fabs(a)<fabs(b)-1.0e-12 || (fabs(fabs(a)-fabs(b))<=1.0e-12 && a<b);});
            
            vector<split> fan;
            for(double a : ang)
            {
                // crack direction: the contact normal turned by a, cut line normal perpendicular to it
                const double dx = cos(a)*ca.nx - sin(a)*ca.ny;
                const double dy = sin(a)*ca.nx + cos(a)*ca.ny;
                const double nx = -dy, ny = dx;
                split sp = line(int(f),2,nx,ny,nx*ca.px + ny*ca.py,ca.F,Fs);
                if(cut_valid(fl,sp,nbx,nby,sb,rb))
                fan.push_back(sp);
                else if(fan.empty())
                break;          // the first crack has to fit, as for a single crack
            }
            
            if(!fan.empty())
            {
            cuts[f] = fan;
            ratio[f] = r;
            }
        }
    }
    
    // spalling at crushing contacts: when the crushed depth reaches L_sp (A 402, <0: h), the tip of the
    // floe in front of the cut n.x = s_face - L_sp breaks off, n from the floe into the other body and
    // s_face the other body's face. Only local chips: the chord of the cut must not exceed C_loc times
    // the chip depth L_sp + crushed depth, a straight cut across a wide floe face is not a spall.
    if((breakflag&4) && p->A398>0.0)
    {
    vector<double> spr(nf,0.0);
    vector<split> spc(nf);
    vector<double> fx,fy,ox,oy;
    
    for(const auto &rc : pcontact->records())
    {
        if(rc.Fcap<=0.0 || rc.pen<=0.0 || rc.b<0)
        continue;
        if(rc.a>=int(cmap.size()) || rc.b>=int(cmap.size()))
        continue;
        
        for(int side=0; side<2; ++side)
        {
            const size_t f = size_t(cmap[side==0 ? rc.a : rc.b]);
            const size_t o = size_t(cmap[side==0 ? rc.b : rc.a]);
            fnpf_ice_floe &fl = floe[f];
            if(fl.type!=0 || ratio[f]>=1.0)
            continue;
            
            const double L = (Lspall>0.0) ? Lspall : fl.h;
            if(rc.pen<L || rc.pen/L<=spr[f])
            continue;
            
            const double nx = (side==0) ? rc.nx : -rc.nx;
            const double ny = (side==0) ? rc.ny : -rc.ny;
            
            world_polygon(fl,fx,fy);
            world_polygon(floe[o],ox,oy);
            
            double sface=1.0e20, sfront=-1.0e20;
            for(size_t q=0;q<ox.size();++q)
            sface = MIN(sface, nx*ox[q] + ny*oy[q]);
            for(size_t q=0;q<fx.size();++q)
            sfront = MAX(sfront, nx*fx[q] + ny*fy[q]);
            
            const double scut = sface - L;
            const double depth = sfront - scut;
            if(depth<=0.0)
            continue;
            
            if(chord(fx,fy,nx,ny,scut) > Cspall*depth)
            continue;
            
            split sp = line(int(f),3,nx,ny,scut,rc.pen,L);
            if(cut_valid(fl,sp,nbx,nby,sb,rb))
            {
            spc[f] = sp;
            spr[f] = rc.pen/L;
            }
        }
    }
    
    for(size_t f=0; f<nf; ++f)
    if(spr[f]>0.0 && cuts[f].empty())
    cuts[f].assign(1,spc[f]);
    }
    
    for(size_t f=0; f<nf; ++f)
    for(const auto &c : cuts[f])
    splits.push_back(c);
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
    b[0]=sp.f; b[1]=sp.mech; b[2]=sp.nx; b[3]=sp.ny; b[4]=sp.s; b[5]=sp.val; b[6]=sp.lim;
    }

    if(ns>0)
    MPI_Bcast(buf.data(), 7*ns, MPI_DOUBLE, 0, pgc->mpi_comm);

    if(p->mpirank>0)
    {
    splits.resize(ns);
    for(int n=0;n<ns;++n)
    {
    const double *b=&buf[7*n];
    splits[n].f=int(b[0]); splits[n].mech=int(b[1]); splits[n].nx=b[2]; splits[n].ny=b[3];
    splits[n].s=b[4]; splits[n].val=b[5]; splits[n].lim=b[6];
    }
    }
    }

    const double t = p->simtime + p->dt;
    
    // the cuts of one floe go one after another to the piece they cross: the floe itself or a piece
    // split off it in this check (same order and arithmetic on all ranks)
    map<int,vector<size_t>> family;
    int ncut=0, nspall=0, nrub=0;

    for(int n=0;n<ns;++n)
    {
        const split &sp = splits[n];
        vector<size_t> &fam = family[sp.f];
        if(fam.empty())
        fam.push_back(size_t(sp.f));
        
        for(size_t c : fam)
        {
            double nbx,nby,sb;
            int rubble;
            if(!cut_valid(floe[c],sp,nbx,nby,sb,rubble))
            continue;
            
            const int mode = (sp.mech==3) ? (rubble ? 2 : 1) : 0;
            double Aa=0.0, Ab=0.0;
            const int pid = floe[c].id;
            const int nid = split_floe(p,c,nbx,nby,sb,mode,Aa,Ab);
            
            if(nid==-1)
            continue;
            
            ++nbreak;
            if(sp.mech==3)
            ++nspall;
            else
            ++ncut;
            
            if(nid>=0)
            fam.push_back(floe.size()-1);
            else
            {
            ++nrub;
            Arubble += Ab;
            }
            
            // the load that broke the floe is released
            for(auto it=cforce.begin(); it!=cforce.end();)
            {
                if(it->first.first==pid || it->first.second==pid)
                it = cforce.erase(it);
                else
                ++it;
            }
            
            if(p->mpirank==0 && breakout.is_open())
            breakout<<setprecision(9)<<t<<" "<<pid<<" "<<nid<<" "<<sp.mech<<" "<<sp.val<<" "<<sp.lim<<" "
                    <<Aa<<" "<<Ab<<endl;
            break;
        }
    }

    if(ns>0 && p->mpirank==0 && ncut+nspall>0)
    {
    cout<<"FNPF ice: ";
    if(ncut>0)
    cout<<ncut<<" crack(s) ";
    if(nspall>0)
    cout<<nspall<<" spall(s) ("<<nrub<<" cleared as rubble, total rubble area "<<Arubble<<" m2) ";
    cout<<"at t = "<<t<<", "<<nfloe<<" floes"<<endl;
    }

    splits.clear();
}

int fnpf_ice::split_floe(lexer *p, size_t f, double nbx, double nby, double s, int mode, double &Aa, double &Ab)
{
    // mode 0: crack, two floes, new strengths for both
    // mode 1: spall, the chip (nb.x >= s) becomes a floe with a new strength, the remainder keeps its strength
    // mode 2: spall, the chip is cleared as rubble (removed), the remainder keeps its strength
    // returns the id of the new floe, -2 for a chip cleared as rubble (type 3, fading), -1 if no cut
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
    Aa = pa.area;
    Ab = pb.area;
    
    if(mode==2)
    {
    // rubble: fixed where it spalled, no contact, the lid fades out over t_r (prestep, footprint);
    // removing the chip at once releases its surface depression as a jet
    pa.id = parent.id;
    pb.id = int(floe.size());
    pb.type = 3;
    pb.t0 = p->simtime + p->dt;
    pb.fade = 1.0;
    for(int a=0;a<3;++a)
    pb.v[a] = pb.w[a] = 0.0;
    kinematics(pb);
    
    floe[f] = pa;
    floe.push_back(pb);
    
    Yn.resize(floe.size());
    D1.resize(floe.size());
    D2.resize(floe.size());
    D3.resize(floe.size());
    return -2;
    }
    
    // new strengths (same random sequence on all ranks)
    if(mode==0)
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
