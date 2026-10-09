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

#include"nhflow_fsf_body.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include<cmath>
#include<sys/stat.h>
#include"definitions.h"

// Water level in the columns pierced by a floating body (X 13 1, direct forcing X 10 1/2).
//
// The sigma grid spans bed to eta in every column, also where the hull pierces the free surface. There the
// upper part of the column is body, forced to the rigid velocity, and its "free surface" is not physical: the
// continuity update moves it with the hull. Its level enters the hydrostatic flux g(eta^2/2 + eta h), the
// Dirichlet condition P = 0 at the top of the column and the hydrostatic pressure on the hull in force_calc_stl
// (wd + eta(x,y) - z on the keel), so the restoring force and moment come from a level that heaves and tilts
// with the body, and P has to make up the difference with a lag.
//
// Between keel and eta the column holds body, not water: the water below the keel is fixed by the hull, so
// the level in a fully pierced column is free. After the water level update of every RK stage it is set to
// the harmonic extension of the surrounding free surface (Laplace over the footprint, Dirichlet data eta of
// the adjacent water columns, Neumann at walls), as fnpf_6DOF::footprint() does for FNPF. The flow, the
// pressure projection and the force calculation then all see the same level, that of the surrounding water.
//
//   footprint   Hb = smoothstep(-FB/(psi/2)) of the top cell of the column, psi the half width of the direct
//               forcing (A 526): Hb = 0 where the top cell centre is outside the hull (real water column,
//               never changed), 1 for a top cell half a forcing width inside. The ramp keeps the level continuous
//               when the hull moves across a column. The columns with 0 < Hb < 1 are the transition columns
//               at the waterline of the hull: their change is the upper bound of the change of real water.
//   solve       SOR over the list of footprint columns, warm started (columns entering the footprint start
//               from their current level), finite volume Laplacian (non-uniform grids), the SOR factor from
//               the footprint size. Without a footprint on more than one rank the sweeps run without
//               communication; otherwise eta_ext is exchanged every sweep and the residual checked every
//               ncheck sweeps.
//   update      eta = eta + Hb (eta_ext - eta), WL with it; detadt gets the change over alpha dt, so that the
//               grid velocity in omega_update is the one of the new level (the water below the keel does not
//               see a grid that moves with the hull). UH, VH, WH of the changed columns are scaled with the
//               new depth before velcalc (momentum()), so that the velocities are unchanged.
//
// Log (rank 0, every X 19 time steps): REEF3D_NHFLOW_6DOF/REEF3D_6DOF_pierced.dat
//   time, number of footprint columns, SOR sweeps since the last line, SOR factor, level change x area since
//   the last line for Hb = 1 (body volume, no water) and 0 < Hb < 1 (transition columns), the latter summed
//   over the run, total volume sum WL dA, water volume sum (1 - FHB) dV, wall time since the last line and
//   of the run.

nhflow_fsf_body::nhflow_fsf_body(lexer *p) : Hb(p), Hbold(p), etaext(p), ratio(p)
{
    mode = (p->X13==1 && (p->X10==1 || p->X10==2) && p->X16==0) ? 1 : 0;
    patch = false;
    told = false;

    gcval_eta = 50;

    if(p->F50==1)
	gcval_eta = 51;

    if(p->F50==2)
	gcval_eta = 52;

    if(p->F50==3)
	gcval_eta = 53;

    if(p->F50==4)
	gcval_eta = 54;

    nfoot = 0;
    itstep = 0;
    pending = false;
    omega = 1.0;
    dVcore = dVedge = dVedge_sum = 0.0;
    tstep = ttot = 0.0;

    SLICELOOP4
    {
    Hb(i,j) = 0.0;
    Hbold(i,j) = 0.0;
    etaext(i,j) = 0.0;
    ratio(i,j) = 1.0;
    }
}

nhflow_fsf_body::~nhflow_fsf_body()
{
    if(out.is_open())
    out.close();
}

void nhflow_fsf_body::level(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, double alpha)
{
    if(!told)
    {
        told = true;

        if(p->mpirank==0 && !patch)
        {
            if(mode==1)
            {
            cout<<"NHFLOW 6DOF X 13 1: water level in the pierced columns from the surrounding free surface"<<endl;

            if(p->G1>0)
            cout<<"NHFLOW 6DOF X 13 1: mesh refinement (G 1), only level 0 is treated"<<endl;
            }

            if(mode==0 && p->X13==1)
            cout<<"NHFLOW 6DOF X 13 1 ignored: needs X 10 1 or 2 and no porous body (X 16)"<<endl;
        }
    }

    pending = false;

    if(mode==0)
    return;

    starttime = pgc->timer();

    // footprint: smoothed Heaviside of the body level set in the top cell of the column
    fi.clear();
    fj.clear();

    k = p->knoz-1;

    SLICELOOP4
    {
        Hbold(i,j) = Hb(i,j);
        Hb(i,j) = 0.0;
        ratio(i,j) = 1.0;

        if(p->wet[IJ]==1)
        {
            double psi;

            if(p->j_dir==0)
            psi = p->A526*p->DXN[IP];

            if(p->j_dir==1)
            psi = p->A526*0.5*(p->DXN[IP] + p->DYN[JP]);

            // ramp over half the forcing width: the first column inside the hull is (nearly) fully treated,
            // also where the top cell centre is less than a forcing width above the keel
            const double dr = 0.5*psi;
            const double phi = -d->FB[IJK];

            if(phi >= dr)
            Hb(i,j) = 1.0;

            else if(phi > 0.0)
            {
            const double t = phi/dr;
            Hb(i,j) = t*t*(3.0 - 2.0*t);
            }
        }

        if(Hb(i,j)>0.0)
        {
        fi.push_back(i);
        fj.push_back(j);

        // column entering the footprint: start the extension from its current level
        if(Hbold(i,j)<=0.0)
        etaext(i,j) = d->eta(i,j);
        }
    }

    nfoot = pgc->globalisum(int(fi.size()));

    if(nfoot==0)
    {
    tstep += pgc->timer() - starttime;
    return;
    }

    pgc->gcsl_start4(p,Hb,1);
    pgc->gcsl_start4(p,etaext,1);

    solve(p,d,pgc);

    // new level, depth and grid velocity of the footprint columns
    const int nf = int(fi.size());

    for(int n=0; n<nf; ++n)
    {
        i = fi[n];
        j = fj[n];

        const double e0 = d->eta(i,j);
        const double e1 = e0 + Hb(i,j)*(etaext(i,j) - e0);
        const double dA = p->DXN[IP]*p->DYN[JP];

        if(Hb(i,j)>=1.0)
        dVcore += (e1 - e0)*dA;

        else
        dVedge += (e1 - e0)*dA;

        const double wl0 = WL(i,j);
        const double wl1 = e1 + d->depth(i,j);

        d->eta(i,j) = e1;
        WL(i,j) = wl1;

        ratio(i,j) = (wl0>p->A544 && wl1>p->A544) ? wl1/wl0 : 1.0;

        d->detadt(i,j) += (e1 - e0)/(alpha*p->dt);

        if(ratio(i,j)!=1.0)
        pending = true;
    }

    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4(p,WL,gcval_eta);
    pgc->gcsl_start4(p,d->detadt,1);

    tstep += pgc->timer() - starttime;
}

void nhflow_fsf_body::solve(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // harmonic extension of eta over the footprint (SOR, warm start)
    const int nf = int(fi.size());
    const int nrank = pgc->globalisum(nf>0 ? 1 : 0);
    const double tol = 1.0e-9*MAX(p->DXM,1.0e-12);
    const int itmax = 1000;
    const int ncheck = 4;

    // SOR factor of the model problem with N cells across the footprint
    const double N = (p->j_dir==1) ? sqrt(double(nfoot)) : double(nfoot);
    omega = 2.0/(1.0 + sin(PI/(N + 1.0)));
    omega = MIN(MAX(omega,1.0),1.9);

    auto add = [&](int ii, int jj, double w, double &s, double &ws)
    {
        // wall: Neumann (the face does not enter the Laplacian)
        if((ii<0 && p->nb1<0) || (ii>=p->knox && p->nb4<0))
        return;

        if(p->j_dir==1)
        if((jj<0 && p->nb3<0) || (jj>=p->knoy && p->nb2<0))
        return;

        if(p->wet[(ii-p->imin)*p->jmax + (jj-p->jmin)]==0)
        return;

        s  += w*(Hb(ii,jj)>0.0 ? etaext(ii,jj) : d->eta(ii,jj));
        ws += w;
    };

    int it;

    for(it=0; it<itmax; ++it)
    {
        double dmax = 0.0;

        for(int n=0; n<nf; ++n)
        {
            i = fi[n];
            j = fj[n];

            double s=0.0, ws=0.0;

            add(i+1,j,1.0/(p->DXN[IP]*p->DXP[IP]),s,ws);
            add(i-1,j,1.0/(p->DXN[IP]*p->DXP[IM1]),s,ws);

            if(p->j_dir==1)
            {
            add(i,j+1,1.0/(p->DYN[JP]*p->DYP[JP]),s,ws);
            add(i,j-1,1.0/(p->DYN[JP]*p->DYP[JM1]),s,ws);
            }

            if(ws>0.0)
            {
            const double de = omega*(s/ws - etaext(i,j));
            etaext(i,j) += de;
            dmax = MAX(dmax,fabs(de));
            }
        }

        if(nrank>1)
        pgc->gcsl_start4(p,etaext,1);

        if((it+1)%ncheck==0)
        {
            if(nrank>1)
            dmax = pgc->globalmax(dmax);

            if(dmax<tol)
            {
            ++it;
            break;
            }
        }
    }

    if(nrank<=1)
    pgc->gcsl_start4(p,etaext,1);

    itstep += pgc->globalimax(it);
}

void nhflow_fsf_body::momentum(lexer *p, double *UH, double *VH, double *WH)
{
    // velocities of the columns whose depth level() changed stay as they are
    if(mode==0 || !pending)
    return;

    const int nf = int(fi.size());

    for(int n=0; n<nf; ++n)
    {
        i = fi[n];
        j = fj[n];

        const double r = ratio(i,j);

        if(r!=1.0)
        KLOOP
        {
        UH[IJK] *= r;
        VH[IJK] *= r;
        WH[IJK] *= r;
        }
    }

    pending = false;
}

void nhflow_fsf_body::print(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
    // every X 19 steps; sweeps, volume changes and wall time are summed in between
    if(mode==0 || p->count%p->X19!=0)
    return;

    double Vt=0.0, Vw=0.0;

    SLICELOOP4
    if(p->wet[IJ]==1)
    Vt += WL(i,j)*p->DXN[IP]*p->DYN[JP];

    LOOP
    if(p->wet[IJ]==1)
    Vw += (1.0 - MIN(MAX(d->FHB[IJK],0.0),1.0))*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);

    Vt = pgc->globalsum(Vt);
    Vw = pgc->globalsum(Vw);

    const double dVc = pgc->globalsum(dVcore);
    const double dVe = pgc->globalsum(dVedge);
    const double ts  = pgc->globalmax(tstep);

    dVedge_sum += dVe;
    ttot += ts;

    if(p->mpirank==0)
    {
        if(!out.is_open())
        {
        mkdir("./REEF3D_NHFLOW_6DOF",0777);
        out.open("./REEF3D_NHFLOW_6DOF/REEF3D_6DOF_pierced.dat");
        out<<"time \t columns \t sweeps \t omega \t dV_core \t dV_trans \t dV_trans_sum \t V_total \t V_water \t t_step \t t_total"<<endl;
        }

        out<<p->simtime<<" \t "<<nfoot<<" \t "<<itstep<<" \t "<<omega<<" \t "<<dVc<<" \t "<<dVe<<" \t "<<dVedge_sum
           <<" \t "<<Vt<<" \t "<<Vw<<" \t "<<ts<<" \t "<<ttot<<endl;
    }

    dVcore = dVedge = 0.0;
    tstep = 0.0;
    itstep = 0;
}
