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

// A 350 4: crest-tracking breaking model for FNPF
//
// Onset (Barthelemy et al. 2018, JFM 841; Derakhti et al. 2020, shallow water):
//   B = u_s / C > A 359 (0.85), u_s = horizontal surface velocity at the crest (= dFi_s/dx there,
//   the crest has deta/dx = 0), C = crest speed from tracking the crest position over time.
// Termination: largest front-face slope below tan(A 360 deg) (8 deg, Simon et al. 2018).
// Strength, per unit crest length and density (the larger of the two):
//   bore / hydraulic jump   eps = A 364 * g C H^3 / (4 h_c)   (discharge C h_t times the head loss H^3/(4 h_c h_t))
//   deep water (A 367 1)    eps = b C^5 / g, b = 0.034 (Gamma - 0.30)^2.5 (Derakhti, Banner & Kirby 2018),
//                           Gamma = T_b dB/dt at onset, T_b = L/C
//   H crest to trough ahead, h_c, h_t water depths at crest and trough.
// Filter: with A352 > 0 the filter of A350 1 (A361, A362) acts on eta and Fifsf in the breaking zones and,
// with A351 1/3, in the cells of the depth-induced criterion of A350 1 (A354) shoreward of the deepest onset.
// Coastline: the coastline / wet-dry front viscosity A346 of A350 1 is applied as well.
// Dissipation: the implicit damping of fnpf_fsfbc::damping with one viscosity per breaking crest on
// its front face (crest to trough ahead, sin^2 ramps at both ends). For a wave of permanent form
// (eta_t = -C eta_x, Fi_t = -C Fi_x) the damping terms dissipate
//   D = nu (fe + ff) C int S eta_x Fi_xx dx,   fe = A 366 (A 357 1, 3), ff = 1 (A 357 1, 2),
// and nu = eps / D-per-nu matches the dissipation of the zone to eps.
//
// The crests are tracked along x (2D, and in 3D along every grid row j: waves mainly in +-x).
// Extrema are gathered on all ranks, so the tracking is the same everywhere.
// Tracking and viscosity once per time step at t^n (first call of the step); the RK stages reuse vb.
//
// Surf zone settings: A352 0 and a coastline band A341 ~10 (A343 1). The A352 filter (A354 0.6) flattens
// the fronts before B reaches the threshold, and the breaking point is that of A350 1 again; without the
// filter and with the default band A341 2 the swash zone of Ting & Kirby goes unstable after 20-36 s.
// Ting & Kirby plunging (dx 0.05): onset at x = 13.2-14.6 m each wave (measured breaking point 14.3 m).

#include"fnpf_breaking.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<algorithm>
#include<sys/stat.h>

namespace
{
    // least-squares slope of (t, v) over the last n entries
    double bt_slope(const std::deque<std::pair<double,double>> &h, int n)
    {
        const int m = std::min((int)h.size(), n);

        if(m<2)
        return 0.0;

        const int s = (int)h.size()-m;
        double tm=0.0, vm=0.0;

        for(int q=s; q<(int)h.size(); ++q)
        {
            tm += h[q].first;
            vm += h[q].second;
        }
        tm/=double(m);
        vm/=double(m);

        double num=0.0, den=0.0;

        for(int q=s; q<(int)h.size(); ++q)
        {
            num += (h[q].first-tm)*(h[q].second-vm);
            den += (h[q].first-tm)*(h[q].first-tm);
        }

        return (den>0.0) ? num/den : 0.0;
    }

    // ramp of the viscosity zone, s in [0,1]: sin^2 over the first and last 10 %
    double bt_ramp(double s)
    {
        const double a = 0.1;

        if(s<0.0 || s>1.0)
        return 0.0;

        if(s<a)
        return pow(sin(0.5*PI*s/a),2.0);

        if(s>1.0-a)
        return pow(sin(0.5*PI*(1.0-s)/a),2.0);

        return 1.0;
    }
}

void fnpf_breaking::breaking_barthelemy(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    // A352 > 0: the breaking filter of A350 1 on the stage values in the breaking zones (every stage)
    if(p->count<=bt_count)
    {
        bt_filter(p,c,pgc,eta,eta_n,Fifsf,alpha);
        return;
    }

    bt_count = p->count;

    const double g = fabs(p->W22);
    const double t = p->simtime;
    const double dxm = p->DXM;

    // ------------------------------------------------------------
    // 1. local extrema of eta along x: sub-grid position and value (parabola through 3 cells),
    //    surface velocity dFi/dx at the extremum, water depth
    // ------------------------------------------------------------
    std::vector<double> loc;

    SLICELOOP4
    if(p->wet[IJ]==1 && p->wet[Im1J]==1 && p->wet[Ip1J]==1)
    {
        const double e0=c->eta(i-1,j), e1=c->eta(i,j), e2=c->eta(i+1,j);
        int type=0;

        if(e1>e0 && e1>=e2)
        type=1;

        if(e1<e0 && e1<=e2)
        type=-1;

        if(type==0)
        continue;

        const double x0=p->XP[IM1], x1=p->XP[IP], x2=p->XP[IP1];

        const double f01 = (e1-e0)/(x1-x0);
        const double f12 = (e2-e1)/(x2-x1);
        const double a = (f12-f01)/(x2-x0);

        double xc = x1;

        if(fabs(a)>1.0e-12)
        xc = 0.5*(x0+x1) - f01/(2.0*a);

        xc = MAX(xc, 0.5*(x0+x1));
        xc = MIN(xc, 0.5*(x1+x2));

        const double ec = e0 + f01*(xc-x0) + a*(xc-x0)*(xc-x1);

        const double F0=c->Fifsf(i-1,j), F1=c->Fifsf(i,j), F2=c->Fifsf(i+1,j);
        const double g01 = (F1-F0)/(x1-x0);
        const double g12 = (F2-F1)/(x2-x1);
        const double ag = (g12-g01)/(x2-x0);
        const double us = g01 + ag*(2.0*xc-x0-x1);

        loc.push_back(double(j+p->origin_j));
        loc.push_back(double(type));
        loc.push_back(xc);
        loc.push_back(ec);
        loc.push_back(us);
        loc.push_back(MAX(c->WL(i,j), 0.0));
    }

    // ------------------------------------------------------------
    // 2. gather on all ranks, sort by row and x
    // ------------------------------------------------------------
    const int nv=6;
    int nloc = (int)loc.size();
    int nproc=1;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    std::vector<int> cnt(nproc), dsp(nproc);
    MPI_Allgather(&nloc,1,MPI_INT,&cnt[0],1,MPI_INT,MPI_COMM_WORLD);

    int ntot=0;
    for(int q=0; q<nproc; ++q)
    {
        dsp[q]=ntot;
        ntot+=cnt[q];
    }

    std::vector<double> all(MAX(ntot,1));
    MPI_Allgatherv(loc.data(),nloc,MPI_DOUBLE,&all[0],&cnt[0],&dsp[0],MPI_DOUBLE,MPI_COMM_WORLD);

    struct ext {int row; int type; double x, eta, us, wl;};
    std::vector<ext> E;

    for(int q=0; q<ntot; q+=nv)
    E.push_back({int(all[q]+0.5), int(lround(all[q+1])), all[q+2], all[q+3], all[q+4], all[q+5]});

    std::sort(E.begin(),E.end(),[](const ext &A, const ext &B){ return (A.row!=B.row) ? A.row<B.row : A.x<B.x; });

    // ------------------------------------------------------------
    // 3. match the crests to the crests of the last call (same row, nearest to the predicted position)
    // ------------------------------------------------------------
    const double dtt = (bt_time>-1.0e19) ? t-bt_time : 0.0;

    std::vector<int> cidx;
    for(int m=0; m<(int)E.size(); ++m)
    if(E[m].type==1)
    cidx.push_back(m);

    struct pair_d {double d; int k; int m;};
    std::vector<pair_d> cand;

    if(dtt>0.0)
    for(int k=0; k<(int)bt_crests.size(); ++k)
    {
        const bt_crest &P = bt_crests[k];
        const double xp = P.x + (P.hasC ? P.C*dtt : 0.0);
        const double tol = P.hasC ? 2.0*dxm + 0.5*fabs(P.C)*dtt : 2.0*dxm + 1.5*sqrt(g*MAX(P.wl,1.0e-3))*dtt;

        for(int m : cidx)
        if(E[m].row==P.row && fabs(E[m].x-xp)<tol)
        cand.push_back({fabs(E[m].x-xp),k,m});
    }

    std::sort(cand.begin(),cand.end(),[](const pair_d &A, const pair_d &B){ return (A.d!=B.d) ? A.d<B.d : ((A.k!=B.k) ? A.k<B.k : A.m<B.m); });

    std::vector<int> kused(bt_crests.size(),0), mto(E.size(),-1);

    for(const pair_d &q : cand)
    if(kused[q.k]==0 && mto[q.m]<0)
    {
        kused[q.k]=1;
        mto[q.m]=q.k;
    }

    std::vector<bt_crest> N;

    for(int m : cidx)
    {
        bt_crest R;

        if(mto[m]>=0)
        R = bt_crests[mto[m]];
        else
        R.id = bt_nextid++;

        R.row = E[m].row;
        R.x   = E[m].x;
        R.eta = E[m].eta;
        R.us  = E[m].us;
        R.wl  = E[m].wl;

        // crest speed: least-squares slope of the sampled positions and the current one. A sample is
        // stored once the crest has moved about a quarter cell (with tiny time steps the step-to-step
        // displacement is far below the accuracy of the sub-grid crest position); up to 8 samples,
        // C from at least 3 samples + the current position (about one cell of travel)
        {
            const double cs = (R.hasC && fabs(R.C)>0.05) ? fabs(R.C) : MAX(sqrt(g*MAX(R.wl,1.0e-3)), 0.05);
            const double dts = 0.25*dxm/cs;
            const bool sample = (R.xh.empty() || t-R.xh.back().first>=dts);

            if(sample)
            {
                R.xh.push_back({t,R.x});
                while((int)R.xh.size()>8)
                R.xh.pop_front();
            }

            std::deque<std::pair<double,double>> hh = R.xh;

            if(!sample)
            hh.push_back({t,R.x});

            R.hasC = ((int)hh.size()>=4);
            R.C = R.hasC ? bt_slope(hh,9) : 0.0;

            if(R.hasC && fabs(R.C)<1.0e-6)
            R.hasC=false;

            R.B = R.hasC ? R.us/R.C : 0.0;

            if(R.hasC && sample)
            {
                R.Bh.push_back({t,R.B});
                while((int)R.Bh.size()>8)
                R.Bh.pop_front();
            }
        }

        N.push_back(R);
    }

    // ------------------------------------------------------------
    // 4. trough ahead / behind, onset
    // ------------------------------------------------------------
    std::vector<double> xa(N.size(),0.0), xb(N.size(),0.0), H(N.size(),0.0), ht(N.size(),0.0);
    std::vector<int> ok(N.size(),0);

    for(int k=0; k<(int)N.size(); ++k)
    {
        bt_crest &R = N[k];

        if(!R.hasC)
        {
            R.state=0;
            continue;
        }

        const int dir = (R.C>0.0) ? 1 : -1;
        int ta=-1, tb=-1;

        for(int m=0; m<(int)E.size(); ++m)
        if(E[m].row==R.row && E[m].type==-1)
        {
            const double d = (E[m].x-R.x)*dir;

            if(d>0.0 && (ta<0 || d<(E[ta].x-R.x)*dir))
            ta=m;

            if(d<0.0 && (tb<0 || -d<-(E[tb].x-R.x)*dir))
            tb=m;
        }

        if(ta<0)
        {
            R.state=0;
            continue;
        }

        // no breaking inside the wave generation and beach relaxation zones (B 96)
        if(p->B90>0 && (R.x < p->global_xmin + p->B96_1 || (p->B99>0 && R.x > p->global_xmax - p->B96_2)))
        {
            R.state=0;
            continue;
        }

        H[k]  = R.eta - E[ta].eta;
        ht[k] = MAX(E[ta].wl, 1.0e-3);
        xa[k] = MIN(R.x, E[ta].x);
        xb[k] = MAX(R.x, E[ta].x);
        ok[k] = (H[k]>2.0*c->wd_criterion) ? 1 : 0;

        if(ok[k]==0)
        {
            R.state=0;
            continue;
        }

        // onset
        if(R.state==0 && R.B>p->A359)
        {
            R.state = 1;
            R.ton = t;

            // still-water depth of the deepest onset so far: the surf zone for the filter (A352)
            bt_db = MAX(bt_db, R.wl-R.eta);

            // deep-water breaking strength from dB/dt at onset (D18)
            R.b = 0.0;

            if(p->A367==1)
            {
                const double L  = (tb>=0) ? fabs(E[ta].x-E[tb].x) : 2.0*fabs(E[ta].x-R.x);
                const double Tb = L/fabs(R.C);
                const double Gm = Tb*bt_slope(R.Bh,5);

                if(Gm>0.30)
                R.b = 0.034*pow(Gm-0.30,2.5);
            }
        }
    }

    // ------------------------------------------------------------
    // 5. front-face integrals of the breaking crests (all ranks), termination, viscosity
    // ------------------------------------------------------------
    const double fe = (p->A357==1 || p->A357==3) ? p->A366 : 0.0;
    const double ff = (p->A357==1 || p->A357==2) ? 1.0 : 0.0;

    std::vector<int> br;
    for(int k=0; k<(int)N.size(); ++k)
    if(N[k].state==1)
    br.push_back(k);

    const int nb = (int)br.size();
    std::vector<double> I(MAX(nb,1),0.0), smax(MAX(nb,1),0.0);

    for(int q=0; q<nb; ++q)
    {
        const int k = br[q];
        j = N[k].row - p->origin_j;

        if(j<0 || j>=p->knoy)
        continue;

        for(i=0; i<p->knox; ++i)
        if(p->flagslice4[IJ]>0 && p->XP[IP]>=xa[k] && p->XP[IP]<=xb[k])
        {
            const double S = bt_ramp((p->XP[IP]-xa[k])/MAX(xb[k]-xa[k],1.0e-12));
            const double ex = (c->eta(i+1,j)-c->eta(i-1,j))/(p->XP[IP1]-p->XP[IM1]);
            const double fxx = ((c->Fifsf(i+1,j)-c->Fifsf(i,j))/(p->XP[IP1]-p->XP[IP]) - (c->Fifsf(i,j)-c->Fifsf(i-1,j))/(p->XP[IP]-p->XP[IM1]))
                               /(0.5*(p->XP[IP1]-p->XP[IM1]));

            I[q] += S*ex*fxx*p->DXN[IP];
            smax[q] = MAX(smax[q], fabs(ex));
        }
    }

    if(nb>0)
    {
        MPI_Allreduce(MPI_IN_PLACE,&I[0],nb,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&smax[0],nb,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    }

    const double slope_end = tan(p->A360*PI/180.0);

    for(int q=0; q<nb; ++q)
    {
        bt_crest &R = N[br[q]];
        const int k = br[q];

        // termination (not in the step of onset)
        if(R.ton<t && smax[q]<slope_end)
        {
            R.state=0;
            R.nu=0.0;
            continue;
        }

        const double hc = MAX(R.wl, 1.0e-3);
        // moving bore: discharge through the jump C h_t, head loss H^3/(4 h_c h_t)
        const double eps_bore = p->A364*g*fabs(R.C)*pow(H[k],3.0)/(4.0*hc);
        const double eps_deep = R.b*pow(fabs(R.C),5.0)/g;
        R.eps = MAX(eps_bore, eps_deep);

        const double dpn = (fe+ff)*R.C*I[q];
        R.dpn = dpn;

        R.nu = (dpn>1.0e-12) ? R.eps/dpn : 0.0;
        // upper bound: diffusion over the front face within one passage of the crest, nu <= |C| L_zone
        R.nu = MIN(R.nu, fabs(R.C)*(xb[k]-xa[k]));
        R.nu = MAX(R.nu, 0.0);
    }

    // ------------------------------------------------------------
    // 6. viscosity field and flags
    // ------------------------------------------------------------
    SLICELOOP4
    {
        c->vb(i,j) = 0.0;
        c->breaking(i,j) = 0;
    }

    int nbreak=0;

    for(int q=0; q<nb; ++q)
    {
        const int k = br[q];

        if(N[k].state!=1)
        continue;

        ++nbreak;
        j = N[k].row - p->origin_j;

        if(j<0 || j>=p->knoy)
        continue;

        // viscosity on the front face, flags (filter A352) on the front face +- 2 cells, which
        // includes the crest cells
        for(i=0; i<p->knox; ++i)
        if(p->flagslice4[IJ]>0 && p->wet[IJ]==1 && p->XP[IP]>=xa[k]-2.0*dxm && p->XP[IP]<=xb[k]+2.0*dxm)
        {
            const double S = bt_ramp((p->XP[IP]-xa[k])/MAX(xb[k]-xa[k],1.0e-12));
            c->vb(i,j) = MAX(c->vb(i,j), N[k].nu*S);
            c->breaking(i,j) = 1;
        }
    }

    // coastline / wet-dry front viscosity A346 as with A350 1 (A343 1 with A337, A343 2/3 with A332)
    coast_viscosity(p,c,pgc);

    pgc->gcsl_start4(p, c->vb, 1);

    SLICELOOP4
    c->breaklog(i,j) = c->breaking(i,j);

    // crest log (P 59 1): crests with B > 0.5 or breaking
    if(p->P59==1 && p->mpirank==0)
    {
        if(!bt_log.is_open())
        {
            mkdir("./REEF3D_FNPF_Breaking_Log",0777);
            bt_log.open("./REEF3D_FNPF_Breaking_Log/REEF3D-FNPF-Breaking-Crests.dat");
            bt_log<<"# t row id x eta C us B state H h_c h_t b eps nu L dpn"<<endl;
        }

        for(int k=0; k<(int)N.size(); ++k)
        if(N[k].hasC && (N[k].B>0.5 || N[k].state==1))
        bt_log<<t<<" "<<N[k].row<<" "<<N[k].id<<" "<<N[k].x<<" "<<N[k].eta<<" "<<N[k].C<<" "<<N[k].us<<" "<<N[k].B<<" "<<N[k].state
              <<" "<<H[k]<<" "<<N[k].wl<<" "<<ht[k]<<" "<<N[k].b<<" "<<N[k].eps<<" "<<N[k].nu<<" "<<xb[k]-xa[k]<<" "<<N[k].dpn<<"\n";
    }

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking crests: "<<nbreak<<endl;

    bt_filter(p,c,pgc,eta,eta_n,Fifsf,alpha);

    bt_crests.swap(N);
    bt_time = t;
}

// A352 > 0: filter eta and Fifsf in the cells of the breaking zones (as A350 1 with A352 3).
// With A351 1/3 also the cells of the depth-induced criterion of A350 1 (eta_t > A354 sqrt(g h), +-2 cells),
// but only shoreward of the deepest breaking onset so far (still-water depth < bt_db): in the surf and swash
// zone the filter keeps the FNPF surface stable (Ting & Kirby stops at t = 20-36 s without it), offshore of
// the break point it would flatten the crests before B reaches the onset threshold.
void fnpf_breaking::bt_filter(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    if(p->A352<=0)
    return;

    if((p->A351==1 || p->A351==3) && p->count>1)
    {
        SLICELOOP4
        {
            bx(i,j)=0;

            if(p->wet[IJ]==1 && p->wd-c->bed(i,j)<bt_db && (eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A354*sqrt(9.81*c->WL(i,j)))
            bx(i,j)=2;
        }

        pgc->gcsl_start4int(p,bx,50);

        SLICELOOP4
        {
            int hit=0;

            for(int q=-2; q<=2; ++q)
            {
                if(bx(i+q,j)==2)
                hit=1;

                if(p->j_dir==1 && bx(i,j+q)==2)
                hit=1;
            }

            if(hit==1)
            c->breaking(i,j) = MAX(c->breaking(i,j), 2);
        }
    }

    SLICELOOP4
    if(c->breaking(i,j)>=1 && p->wet[IJ]==1)
    {
        filter(p,c,pgc,eta,p->A361,p->A362);
        filter(p,c,pgc,Fifsf,p->A361,p->A362);
    }
}
