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
#include<cmath>

/*--------------------------------------------------------------------
Q 88 > 0: implicit Coulomb friction of the packed bed on the grid of CPM (step 2 of the deviatoric stress,
validation 36)

The parcels of the packed cells take their horizontal friction from a shear stress on the horizontal planes of the
grid, instead of the friction on their substrate:
    tau_xz = eta du/dz,  tau_yz = eta dv/dz,   eta = mu(I) P / sqrt(gamma^2 + gamma_r^2),
    gamma = |d(u,v)/dz|,  gamma_r = Q 88 [1/s]
with P = Tau, the packed-bed stress (the overburden, the normal stress on the horizontal planes), at least one grain
layer, and mu(I) of the friction. This is the Mohr-Coulomb yield on horizontal planes, the criterion of the column
yield ratio Yr of Q 67 and of the friction of the parcels inside the bed (contact normal = gravity), but implicit
and continuous through the column: below yield eta is large and the bed shears at about gamma_r at most, above
yield the stress is mu(I) P and the bed flows. The horizontal motion of the bed is then held by this friction, not
by the jam. The jam (Q 64, Q 67) still holds the vertical motion of the bed below the surface layer, and the
sleeping grains stay at rest.

Only the stress on horizontal planes: a column is not held sideways by its neighbours (sand has no tensile
strength, compression is taken by the isotropic stress Tau), and the vertical velocity stays with the isotropic
stress, the normal contact and the jam.

Per sub-step, after the parcel velocities (drag, forces, normal contact):
  1. the velocities of the parcels of the packed cells (not in the bedload layer) onto the grid (nearest cell), U*;
  2. cells: active (packed, theta >= theta_bed), rigid (solids, walls, mostly sleeping grains: at rest), open
     (theta < theta_bed: no stress through the face, the surface of the bed);
  3. eta from the shear rate of the parcel velocity at the start of the sub-step on the grid (lagged: with the
     shear rate of U* the acceleration of one sub-step would set the viscosity, and a slope below yield would creep
     at that rate; nothing is kept between the steps, so a restart is exact),
     then for u and v
         (theta rho_s (1 + dt D)/dt) (U - U*) = d/dz(eta dU/dz)
     with the point-implicit drag D of the parcels (U* = (U + dt F)/(1 + dt D)), so that the friction balances the
     full force F at rest. Tridiagonal per column, with a vertical domain decomposition conjugate gradients;
  4. the parcels of the packed cells take the grid velocity (PIC: interpolated from the active cells and from the
     rigid cells of their own column; the vertical velocity times 1 - Hjam), those of rigid cells stay at rest;
     the parcels of the open cells keep the friction on their substrate.
--------------------------------------------------------------------*/

void CPM::shear_alloc(lexer *p)
{
    if(sh_alloc)
    return;

    shU = new field4a(p); shV = new field4a(p); shW = new field4a(p);
    shF = new field4a(p); shJ = new field4a(p); shE = new field4a(p);
    shC = new field4a(p); shR = new field4a(p); shZ = new field4a(p);
    shP = new field4a(p); shA = new field4a(p); shX = new field4a(p);
    shUo = new field4a(p); shVo = new field4a(p); shI = new field4a(p);

    sh_alloc = true;
}

void CPM::shear_implicit(lexer *p, fdm *a, ghostcell *pgc, double dt)
{
    shear_alloc(p);

    field4a &Ug = *shU, &Vg = *shV, &Wg = *shW, &Wf = *shF, &Wj = *shJ, &Eta = *shE, &Cl = *shC, &Fi = *shI;
    field4a &Uo = *shUo, &Vo = *shVo;

    const double gmag = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    const double rog = (p->S22 - p->W1)*gmag;
    const double gr = p->Q88;

    // 1. the parcels onto the grid
    for(i=-1;i<p->knox+1;++i)
    for(j=-1;j<p->knoy+1;++j)
    for(k=-1;k<p->knoz+1;++k)
    Ug(i,j,k) = Vg(i,j,k) = Wg(i,j,k) = Wf(i,j,k) = Wj(i,j,k) = Fi(i,j,k) = Uo(i,j,k) = Vo(i,j,k) = 0.0;

    // (only the parcels of the packed cells: the grains of the open cells above the bed slide on their substrate;
    // nearest cell, so that the jammed grains of a neighbouring column do not hold a yielded one sideways)
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE && n<int(sh_hj.size()) && sh_hj[n]>=0.0)
    {
        i = pci(p,n,P.X[n]);
        j = pcj(p,n,P.Y[n]);
        k = pck(p,n,P.Z[n]);
        
        if(Ts(i,j,k)<theta_bed)
        continue;
        
        // the mean velocity of the grains of the cell, the sleeping ones at rest
        double w = P.ParcelFactor*Vp*P.Wt[n];
        
        if(sh_hj[n]>1.5)
        Wj(i,j,k) += w;
        
        else
        Wf(i,j,k) += w;
        
        Fi(i,j,k) += w*sh_fi[n];
        Ug(i,j,k) += w*P.URK1[n];
        Vg(i,j,k) += w*P.VRK1[n];
        Wg(i,j,k) += w*P.WRK1[n];
        
        // the velocity at the start of the sub-step, for the viscosity
        Uo(i,j,k) += w*P.U[n];
        Vo(i,j,k) += w*P.V[n];
    }

    // 2. classes: 1 active, 0 open, -1 rigid (solids, mostly sleeping grains)
    BASELOOP
    {
        double wt = Wf(i,j,k) + Wj(i,j,k);
        
        if(wt>1.0e-14)
        {
            Ug(i,j,k) /= wt;
            Vg(i,j,k) /= wt;
            Wg(i,j,k) /= wt;
            Fi(i,j,k) /= wt;
            Uo(i,j,k) /= wt;
            Vo(i,j,k) /= wt;
        }
        else
        {
            Ug(i,j,k) = Vg(i,j,k) = Wg(i,j,k) = Uo(i,j,k) = Vo(i,j,k) = 0.0;
            Fi(i,j,k) = 1.0;
        }

        double cl = 0.0;

        if(p->solidread>0 && a->solid(i,j,k)<0.0)
        cl = -1.0;

        else if(Ts(i,j,k)>=theta_bed)
        cl = Wf(i,j,k) >= Wj(i,j,k) ? 1.0 : -1.0;

        Cl(i,j,k) = cl;
    }

    pgc->start4a(p,Ug,1);
    pgc->start4a(p,Vg,1);
    pgc->start4a(p,Wg,1);
    pgc->start4a(p,Uo,1);
    pgc->start4a(p,Vo,1);
    pgc->start4a(p,Cl,1);

    auto cls = [&](int ii, int jj, int kk) -> double
    {
        if(wallcell(p,ii,jj,kk))
        return -1.0;

        return Cl(ii,jj,kk);
    };

    // 3a. the viscosity from the shear rate of the velocity at the start of the sub-step (lagged: with the shear
    //     rate of U* the acceleration of one sub-step alone would set the viscosity, and a slope below yield would
    //     creep at that rate). Vertical derivatives towards active and rigid neighbours (at rest), one-sided at
    //     open cells.
    
    auto val = [&](field &f, int kk) { return cls(i,j,kk)>0.5 ? f(i,j,kk) : 0.0; };
    
    auto ddz = [&](field &f) -> double
    {
        bool okm = cls(i,j,k-1)!=0.0, okp = cls(i,j,k+1)!=0.0;
        
        if(okm && okp)
        return (val(f,k+1)-val(f,k-1))/(p->DZP[KP]+p->DZP[KM1]);
        
        if(okp)
        return (val(f,k+1)-f(i,j,k))/p->DZP[KP];
        
        if(okm)
        return (f(i,j,k)-val(f,k-1))/p->DZP[KM1];
        
        return 0.0;
    };
    
    int nact = 0;
    
    BASELOOP
    {
        Eta(i,j,k) = 0.0;
        
        if(Cl(i,j,k)<0.5)
        continue;
        
        ++nact;
        
        double dudz = ddz(Uo);
        double dvdz = p->j_dir==1 ? ddz(Vo) : 0.0;
        double gamma = sqrt(dudz*dudz + dvdz*dvdz);
        
        double P = MAX(Tau(i,j,k), rog*MAX(Ts(i,j,k),theta_bed)*p->S20);
        double I = p->S20*gamma/sqrt(P/p->S22);
        double mu = mu_s + (mu_2-mu_s)/(I0/MAX(I,1.0e-10) + 1.0);
        
        Eta(i,j,k) = mu*P/sqrt(gamma*gamma + gr*gr);
    }
    
    pgc->start4a(p,Eta,1);
    
    nact = pgc->globalisum(nact);
    sh_nact = nact;
    
    // 3b. u and v: A U = b on the active cells (no active cell: no solve, the parcels of the rigid cells are still
    //     held in 4.)
    field4a &R = *shR, &Z = *shZ, &Pd = *shP, &Ap = *shA, &X = *shX;
    
    // coefficient of the neighbour kk per unit volume; a rigid neighbour goes to the diagonal only
    auto efac = [&](int kk, double hf, double &diag) -> double
    {
        double cn = cls(i,j,kk);
        
        if(cn==0.0)
        return 0.0;
        
        double e0 = Eta(i,j,k);
        double ef = cn>0.5 ? 2.0*e0*Eta(i,j,kk)/MAX(e0+Eta(i,j,kk),1.0e-300) : e0;
        double c = ef/(p->DZN[KP]*hf);
        
        diag += c;
        
        return cn>0.5 ? c : 0.0;
    };
    
    // y = A x on the active cells (x with valid halo)
    auto apply = [&](field4a &x, field4a &y, bool diagonly)
    {
        BASELOOP
        {
            y(i,j,k) = 0.0;
            
            if(Cl(i,j,k)<0.5)
            continue;
            
            double diag = Ts(i,j,k)*p->S22*Fi(i,j,k)/dt, off = 0.0, c;
            
            c = efac(k-1,p->DZP[KM1],diag); off += c*x(i,j,k-1);
            c = efac(k+1,p->DZP[KP],diag);  off += c*x(i,j,k+1);
            
            y(i,j,k) = diagonly ? diag : diag*x(i,j,k) - off;
        }
    };
    
    // the diagonal (Jacobi preconditioner)
    field4a &Dg = *shJ;
    apply(Ug,Dg,true);
    
    field4a *comp[2] = {&Ug,&Vg};
    int itsum = 0;
    
    for(int c=0; c<2 && nact>0; ++c)
    {
        if(c==1 && p->j_dir==0)
        continue;
        
        field4a &Us_ = *comp[c];
        
        // the columns within one subdomain: tridiagonal, direct
        if(zsplit==0)
        {
            const int nz = p->knoz;
            std::vector<double> lo(nz), di(nz), up(nz), rh(nz);
            
            ILOOP
            JLOOP
            {
                bool any = false;
                
                for(k=0; k<nz; ++k)
                {
                    if(Cl(i,j,k)>0.5)
                    {
                        any = true;
                        double m = Ts(i,j,k)*p->S22*Fi(i,j,k)/dt;
                        double d = m;
                        lo[k] = -efac(k-1,p->DZP[KM1],d);
                        up[k] = -efac(k+1,p->DZP[KP],d);
                        di[k] = d;
                        rh[k] = m*Us_(i,j,k);
                    }
                    else
                    {
                        lo[k] = up[k] = 0.0;
                        di[k] = 1.0;
                        rh[k] = Us_(i,j,k);
                    }
                }
                
                if(!any)
                continue;
                
                for(k=1; k<nz; ++k)
                {
                    double q = lo[k]/di[k-1];
                    di[k] -= q*up[k-1];
                    rh[k] -= q*rh[k-1];
                }
                
                rh[nz-1] /= di[nz-1];
                
                for(k=nz-2; k>=0; --k)
                rh[k] = (rh[k] - up[k]*rh[k+1])/di[k];
                
                for(k=0; k<nz; ++k)
                if(Cl(i,j,k)>0.5)
                Us_(i,j,k) = rh[k];
            }
            
            continue;
        }
        
        // with a vertical domain decomposition: conjugate gradients
        // b = (theta rho_s/dt) U*, start from x = U*
        double bnorm = 0.0;
        
        BASELOOP
        {
            X(i,j,k) = Cl(i,j,k)>0.5 ? Us_(i,j,k) : 0.0;
            
            if(Cl(i,j,k)>0.5)
            bnorm += pow(Ts(i,j,k)*p->S22*Fi(i,j,k)/dt*Us_(i,j,k),2.0);
        }
        
        pgc->start4a(p,X,1);
        bnorm = sqrt(pgc->globalsum(bnorm));
        
        apply(X,Ap,false);
        
        double rz = 0.0;
        
        BASELOOP
        {
            R(i,j,k) = Z(i,j,k) = Pd(i,j,k) = 0.0;
            
            if(Cl(i,j,k)<0.5)
            continue;
            
            R(i,j,k) = Ts(i,j,k)*p->S22*Fi(i,j,k)/dt*Us_(i,j,k) - Ap(i,j,k);
            Z(i,j,k) = R(i,j,k)/Dg(i,j,k);
            Pd(i,j,k) = Z(i,j,k);
            rz += R(i,j,k)*Z(i,j,k);
        }
        
        rz = pgc->globalsum(rz);
        
        int it;
        
        for(it=0; it<500; ++it)
        {
            double rn = 0.0;
            
            BASELOOP
            if(Cl(i,j,k)>0.5)
            rn += R(i,j,k)*R(i,j,k);
            
            rn = sqrt(pgc->globalsum(rn));
            
            if(rn <= 1.0e-8*MAX(bnorm,1.0e-300))
            break;
            
            pgc->start4a(p,Pd,1);
            apply(Pd,Ap,false);
            
            double pap = 0.0;
            
            BASELOOP
            if(Cl(i,j,k)>0.5)
            pap += Pd(i,j,k)*Ap(i,j,k);
            
            pap = pgc->globalsum(pap);
            
            if(fabs(pap)<1.0e-300)
            break;
            
            double alpha = rz/pap, rznew = 0.0;
            
            BASELOOP
            if(Cl(i,j,k)>0.5)
            {
                X(i,j,k) += alpha*Pd(i,j,k);
                R(i,j,k) -= alpha*Ap(i,j,k);
                Z(i,j,k) = R(i,j,k)/Dg(i,j,k);
                rznew += R(i,j,k)*Z(i,j,k);
            }
            
            rznew = pgc->globalsum(rznew);
            
            double beta = rznew/MAX(rz,1.0e-300);
            rz = rznew;
            
            BASELOOP
            if(Cl(i,j,k)>0.5)
            Pd(i,j,k) = Z(i,j,k) + beta*Pd(i,j,k);
        }
        
        itsum += it;
        
        BASELOOP
        if(Cl(i,j,k)>0.5)
        Us_(i,j,k) = X(i,j,k);
    }
    
    sh_iter = itsum;
    
    // the grid velocity: rigid cells at rest; open cells keep U* (not used: the interpolation in 4. skips them)
    BASELOOP
    if(Cl(i,j,k)<-0.5)
    Ug(i,j,k) = Vg(i,j,k) = Wg(i,j,k) = 0.0;
    
    pgc->start4a(p,Ug,1);
    pgc->start4a(p,Vg,1);
    pgc->start4a(p,Wg,1);
    
    // 4. the grains of the packed cells move with the grid (PIC): those of the active cells with the new velocity,
    //    those of the rigid cells stay; the grains of the open cells keep their friction (friction())
    const double gy = p->j_dir==1 ? 1.0 : 0.0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE && n<int(sh_hj.size()) && sh_hj[n]>=0.0 && sh_hj[n]<=1.0)
    {
        i = pci(p,n,P.X[n]);
        j = pcj(p,n,P.Y[n]);
        k = pck(p,n,P.Z[n]);
        
        if(Ts(i,j,k)<theta_bed)
        continue;
        
        double f = 1.0-sh_hj[n];
        double un=0.0, vn=0.0, wn=0.0;
        
        if(Cl(i,j,k)>0.5)
        {
            // interpolation from the packed cells only (active and rigid), not from the open cells above
            const int ii = pfi(p,n,P.X[n]), jj = p->j_dir==1 ? pfj(p,n,P.Y[n]) : 0, kk_ = pfk(p,n,P.Z[n]);
            const double wa = (p->XP[ii+1+marge]-P.X[n])/p->DXP[ii+marge];
            const double wb = p->j_dir==1 ? (p->YP[jj+1+marge]-P.Y[n])/p->DYP[jj+marge] : 1.0;
            const double wc = (p->ZP[kk_+1+marge]-P.Z[n])/p->DZP[kk_+marge];
            double ws = 0.0;
            
            for(int a1=0;a1<2;++a1)
            for(int b1=0;b1<(p->j_dir==1?2:1);++b1)
            for(int c1=0;c1<2;++c1)
            {
                const int i1=ii+a1, j1=jj+b1, k1=kk_+c1;
                
                const double c1s = cls(i1,j1,k1);
                
                // open cells not; rigid cells (at rest) only in the column of the parcel (the bed below)
                if(c1s==0.0 || (c1s<0.0 && (i1!=i || j1!=j)))
                continue;
                
                const double w = (a1==0?wa:1.0-wa)*(b1==0?wb:1.0-wb)*(c1==0?wc:1.0-wc);
                ws += w;
                
                if(c1s<0.0)
                continue;
                
                un += w*Ug(i1,j1,k1);
                vn += w*Vg(i1,j1,k1);
                wn += w*Wg(i1,j1,k1);
            }
            
            // the jam holds the vertical motion only (f = 1 - Hjam)
            if(ws>1.0e-12)
            {
                un *= 1.0/ws;
                vn *= gy/ws;
                wn *= f/ws;
            }
            else
            {
                un = Ug(i,j,k);
                vn = gy*Vg(i,j,k);
                wn = f*Wg(i,j,k);
            }
        }
        
        P.XRK1[n] += dt*(un-P.URK1[n]);
        P.YRK1[n] += dt*(vn-P.VRK1[n]);
        P.ZRK1[n] += dt*(wn-P.WRK1[n]);
        P.URK1[n] = un;
        P.VRK1[n] = vn;
        P.WRK1[n] = wn;
    }
}
