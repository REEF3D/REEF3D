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

#include"6DOF_obj_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"vrans_definitions.h"
#include"nhflow_bodyflux.h"
#include"nhflow_srchull.h"
#include<fstream>
#include<vector>
#include<cmath>

void sixdof_obj_nhflow::update_forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, int iter)
{
    // porous floating body: Darcy-Forchheimer resistance instead of rigid direct forcing
    if(p->X16==1)
    {
    update_forcing_nhflow_porous(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter);
    return;
    }
    
    forcing_kernel_nhflow(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter,false);
}

void sixdof_obj_nhflow::update_forcing_nhflow_wavemaker(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, int iter)
{
    // wavemaker (X 10 4): the same direct forcing (X 15 0/1/2) as the floating body, with the
    // paddle velocity of the stage as the target
    forcing_kernel_nhflow(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter,true);
}

void sixdof_obj_nhflow::target_velocity(lexer *p, bool wavemaker, double &uf, double &vf, double &wf)
{
    // target velocity of the direct forcing at the cell (i,j,k)
    if(wavemaker)
    {
        uf = wm_u(p,p->pos_x(),p->pos_z());
        vf = 0.0;
        wf = 0.0;
        return;
    }
    
    uf = u_fb(0) + u_fb(4)*(p->pos_z() - c_(2)) - u_fb(5)*(p->pos_y() - c_(1));
    vf = u_fb(1) + u_fb(5)*(p->pos_x() - c_(0)) - u_fb(3)*(p->pos_z() - c_(2));
    wf = u_fb(2) + u_fb(3)*(p->pos_y() - c_(1)) - u_fb(4)*(p->pos_x() - c_(0));
}

void sixdof_obj_nhflow::forcing_kernel_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, int iter, bool wavemaker)
{
    // X 13 2: no body cells in the water columns, no direct forcing
    // X 15 8: the hull is a sharp thin body (blocked links), no direct forcing
    if(d->hull!=nullptr || d->hullbody!=nullptr)
    return;
    
    // EXPERIMENT NHF_SRC (nhflow_srchull.h): the hull as mass sources, q = (1 + k1) Us d(chi)/dx, no forcing
    if(nhflow_srchull::active() && !wavemaker)
    {
        const int nn = p->imax*p->jmax*(p->kmax+2);
        
        if(int(nhflow_srchull::q.size())!=nn)
        nhflow_srchull::q.assign(nn,0.0);
        
        const double Us = nhflow_srchull::stream() - u_fb(0);
        const double fk = 1.0 + nhflow_srchull::k1();
        double Fx=0.0, Qp=0.0, Qn=0.0, Ffx=0.0, Ffz=0.0, Fqz=0.0;
        
        // NHF_SRC_FB = 1: feedback strength, 2: without the vertical term w d(chi)/dz; updated once per time step
        const bool fbk = nhflow_srchull::feedback()>0 && (iter==0 || nhflow_srchull::feedback()>=10);
        const double fbz = (nhflow_srchull::feedback()%10==2) ? 0.0 : 1.0;
        bool keep = nhflow_srchull::feedback()>0 && !fbk;
        static const double fbw = getenv("NHF_SRC_FBW")!=nullptr ? atof(getenv("NHF_SRC_FBW")) : 0.1;
        static const double wlf = getenv("NHF_SRC_WL")!=nullptr ? atof(getenv("NHF_SRC_WL")) : 0.0;
        static const bool nometric = getenv("NHF_SRC_NOMETRIC")!=nullptr && atoi(getenv("NHF_SRC_NOMETRIC"))>0;
        
        // NHF_SRC_RAMP = T [s]: the sources grow smoothly from 0 to their strength over the first T seconds. With the
        // transparent interior (cells inside forced to the stream) the hull then appears gradually out of the undisturbed
        // stream instead of at once: no impulsive start, much weaker sloshing of the tank (period 4 s in the Wigley tank)
        static const double tramp = getenv("NHF_SRC_RAMP")!=nullptr ? atof(getenv("NHF_SRC_RAMP")) : 0.0;
        const double ramp = (tramp>0.0 && p->simtime<tramp) ? 0.5*(1.0 - cos(PI*p->simtime/tramp)) : 1.0;
        static bool first = true;
        
        // NHF_SRC_STL = file: the sources from the hull geometry instead of the level set (Cases/WIGLEY/srcpts.py).
        // Lines "x y z a type": type 0 a hull facet, a = -n_x dA [m^2], source Us a in the cell that holds the point,
        // if the point is below the free surface of the column; type 1 a piece of the water line, a = n_x^3 dl/|n_h|
        // [m], source (Us^3/g) a in the top layer (the water line term, factor NHF_SRC_WL). Bilinear in x and y over
        // the four nearest columns; uniform Cartesian grid, fixed body. The smoothed level set of a thin stem gives
        // the wrong normal there, and the sources of the level set converge slowly (Wigley, Cases/WIGLEY).
        static const char *ptsfile = getenv("NHF_SRC_STL");
        static const bool totop = getenv("NHF_SRC_TOP")!=nullptr && atoi(getenv("NHF_SRC_TOP"))>0;
        
        if(ptsfile!=nullptr)
        {
            static std::vector<double> px,py,pz,pa;
            static std::vector<int> pt;
            static bool read = false;
            
            if(!read)
            {
                read = true;
                std::ifstream in(ptsfile);
                double a0,a1,a2,a3; int a4;
                
                while(in>>a0>>a1>>a2>>a3>>a4)
                {px.push_back(a0); py.push_back(a1); pz.push_back(a2); pa.push_back(a3); pt.push_back(a4);}
                
                if(p->mpirank==0)
                cout<<"NHF_SRC_STL: "<<px.size()<<" source points from "<<ptsfile<<endl;
            }
            
            LOOP
            nhflow_srchull::q[IJK] = 0.0;
            
            const double dxu = p->DXN[0+marge];
            const double dyu = p->DYN[0+marge];
            const double g = fabs(p->W22);
            
            for(size_t n=0;n<px.size();++n)
            {
                const double fx = (px[n] - p->XP[0+marge])/dxu;
                const double fy = (p->j_dir==0) ? 0.0 : (py[n] - p->YP[0+marge])/dyu;
                const int i0 = int(floor(fx));
                const int j0 = int(floor(fy));
                
                if(i0<-1 || i0>=p->knox || j0<-1 || j0>=p->knoy)
                continue;
                
                const double m = (pt[n]==0) ? fk*Us*pa[n] : wlf*fk*(Us*Us*Us/g)*pa[n];
                
                if(m==0.0)
                continue;
                
                for(int qi=0;qi<2;++qi)
                for(int qj=0;qj<((p->j_dir==0)?1:2);++qj)
                {
                    i = i0 + qi;
                    j = j0 + qj;
                    
                    if(i<0 || i>=p->knox || j<0 || j>=p->knoy)
                    continue;
                    
                    if(p->wet[IJ]!=1)
                    continue;
                    
                    const double w = (qi==0 ? 1.0-(fx-i0) : fx-i0)*((p->j_dir==0) ? 1.0 : (qj==0 ? 1.0-(fy-j0) : fy-j0));
                    
                    if(pt[n]==1)
                    k = p->knoz-1;
                    else
                    {
                        const double sg = (pz[n] - d->bed(i,j))/WL(i,j);
                        
                        // NHF_SRC_TOP=1: a facet above the free surface of the column (a wave trough at the hull) keeps
                        // its source, in the top layer; with the facets up to the still water line only (srcpts.py
                        // still) the sources then always sum to zero, as in the linear theory. Without it the wetted
                        // facets only: the bow wave wets freeboard facets, a trough dries facets, and the sources
                        // of -Us n_x no longer sum to zero (Wigley Fr 0.267: + 0.111 / - 0.103 m^3/s)
                        if(sg<p->ZN[0+marge])
                        continue;
                        
                        if(sg>=p->ZN[p->knoz+marge] && !totop)
                        continue;
                        
                        for(k=0;k<p->knoz-1;++k)
                        if(sg<p->ZN[KP1])
                        break;
                    }
                    
                    nhflow_srchull::q[IJK] += w*m/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j));
                }
            }
            
            keep = true;
        }
        
        // NHF_SRC_FORCE=1: transparent interior. The cells of the hull are forced to the undisturbed stream (NHF_SRC, 0, 0)
        // with the weight chi, instead of to the body velocity. The source sheet jumps the normal velocity by
        // sigma = -Us n_x, the water inside moves with the stream, u_n = Us n_x, so the water outside has u_n = 0:
        // no flow through the hull, by mass conservation, for any hull shape (no k1). Discretely: the blended field
        // u = u_e + chi (Us e_x - u_e) has div u = Us d(chi)/dx - u_e.grad(chi) = q exactly when u_e is tangent to
        // the hull. The velocity is continuous across the hull up to the perturbation of the stream (a no-slip body
        // has the jump Us there, and the reconstruction turns it into a numerical boundary layer and a wake).
        // NHF_SRC_FORCE=2: only the cells inside the source band (FB < -psi) are forced, fully. With the weight chi (1) the
        // direct forcing of every stage pins the whole band to the stream after a few steps, the sources then sit in
        // cells that cannot diverge, and the wake is back (2D ellipse: C_Dw 0.35, plain sources 0.00).
        static const int tfmode = getenv("NHF_SRC_FORCE")!=nullptr ? atoi(getenv("NHF_SRC_FORCE")) : 0;
        static const bool tforce = tfmode>0;
        
        // NHF_SRC_ANISO=1: the half width of the band from the grid spacing normal to the hull, psi = A 526 (|n_x| dx + |n_y| dy)
        // (weighted with the horizontal part of the normal; A 526 (dx + dy)/2 where the hull faces down), instead of
        // A 526 (dx + dy)/2 everywhere. On a grid stretched to the centre line (dy << dx) the sides of the hull are then as
        // sharp as dy allows. chi is one field, chi(x) = H(-FB(x)/psi(x)), for the sources and the forcing.
        // NHF_SRC_ANISO=2: also with the vertical part of the normal and the thickness of the layer,
        // psi = A 526 (|n_x| dx + |n_y| dy + |n_z| dz). With psi from dx, dy only the band below the hull (and above and below
        // a submerged body) is thinner than one layer wherever the layers are thicker than the horizontal cells (Wigley
        // L/112: psi 0.039 m, layers of 0.06-0.08 m at the keel): the forced cells (FB < -psi) then touch the cells of the
        // sources, which cannot diverge (as with NHF_SRC_FORCE=1). The 2D tests have square cells.
        static const int anisomode = getenv("NHF_SRC_ANISO")!=nullptr ? atoi(getenv("NHF_SRC_ANISO")) : 0;
        static const bool aniso = anisomode>0;
        
        auto IDX = [&](int ii, int jj, int kk)
        {
            return (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;
        };
        
        // NHF_SRC_PSI = psi [m]: the same band half width everywhere, for grids refined to the hull. The band of the cell
        // size (A 526 (dx + dy)/2, and the NHF_SRC_ANISO variants) changes from cell to cell on a stretched grid; the
        // surfaces chi = const are then not parallel to the hull, and the water outside is held tangent to those surfaces,
        // not to the hull (the "sources" of a changing band width, NHF_SRC_ANISO=2). Choose psi >= 2 cells normal to the
        // hull (2 dy at the sides, 2 dz at the bottom) and small against the half beam.
        static const double psifix = getenv("NHF_SRC_PSI")!=nullptr ? atof(getenv("NHF_SRC_PSI")) : 0.0;
        
        auto psiAt = [&](int ii, int jj, int kk)
        {
            if(psifix>0.0)
            return psifix;
            
            const double piso = (p->j_dir==0) ? p->A526*p->DXN[ii+marge] : p->A526*0.5*(p->DXN[ii+marge] + p->DYN[jj+marge]);
            
            if(!aniso)
            return piso;
            
            const int kp = MIN(kk+1,p->knoz-1);
            const int km = MAX(kk-1,0);
            const double gx = (d->FB[IDX(ii+1,jj,kk)] - d->FB[IDX(ii-1,jj,kk)])/(p->DXP[ii+marge] + p->DXP[ii-1+marge]);
            const double gy = (p->j_dir==0) ? 0.0 : (d->FB[IDX(ii,jj+1,kk)] - d->FB[IDX(ii,jj-1,kk)])/(p->DYP[jj+marge] + p->DYP[jj-1+marge]);
            const double gz = (kp>km) ? (d->FB[IDX(ii,jj,kp)] - d->FB[IDX(ii,jj,km)])/((p->ZP[kp+marge] - p->ZP[km+marge])*WL(ii,jj)) : 0.0;
            const double gh = sqrt(gx*gx + gy*gy);
            const double g2 = gh*gh + gz*gz;
            
            if(anisomode>=2)
            return (g2<1.0e-20) ? piso : p->A526*(fabs(gx)*p->DXN[ii+marge] + fabs(gy)*p->DYN[jj+marge] + fabs(gz)*p->DZN[kk+marge]*WL(ii,jj))/sqrt(g2);
            
            if(g2<1.0e-20 || gh<1.0e-10)
            return piso;
            
            const double ph = p->A526*(fabs(gx)*p->DXN[ii+marge] + fabs(gy)*p->DYN[jj+marge])/gh;
            const double w = gh*gh/g2;
            
            return w*ph + (1.0-w)*piso;
        };
        
        auto chiAt = [&](int ii, int jj, int kk)
        {
            const double ps = psiAt(ii,jj,kk);
            const double fbv = -d->FB[IDX(ii,jj,kk)];
            
            return fbv>ps ? 1.0 : (fbv<-ps ? 0.0 : 0.5*(1.0 + fbv/ps + (1.0/PI)*sin(PI*fbv/ps)));
        };
        
        // NHF_SRC_ANISO=3: as 2, but d(chi)/dx of the sources with the band width of the cell itself for its neighbours.
        // chi = H(-FB/psi(x)) with psi changing along the hull (normal, layer thickness, water depth) has a second part
        // in its derivative, -H' FB d(1/psi)/dx: sources of opposite sign on the two sides of the hull surface that are no
        // part of U n_x. They grow with A 526 (submerged Wigley double body, NHF_SRC_ANISO=2: C_wp 0.199 / 0.281 / 0.375e-3
        // for A 526 1.1 / 1.6 / 2.1, plain sources without any forcing 0.234, linear theory 0.30).
        auto chiP = [&](int ii, int jj, int kk, double ps)
        {
            const double fbv = -d->FB[IDX(ii,jj,kk)];
            
            return fbv>ps ? 1.0 : (fbv<-ps ? 0.0 : 0.5*(1.0 + fbv/ps + (1.0/PI)*sin(PI*fbv/ps)));
        };
        
        // NOT WORKING YET (19 Nov 2026, 2D tests Cases/forcing_tests mr_ellipse, mr_cyl): unstable with the full update;
        // relaxed (NHF_SRC_MREL 0.1-0.3) it runs, but the ellipse then leaks (speed over the crest 1.00 U, exact 1.23)
        // and the wave of the cylinder is unsteady and wrong. The divergence of the ghost field is not taken out of the
        // pressure projection (the sources of the other modes do that); kept as a starting point only.
        // NHF_SRC_FORCE=3: mirror hull, no sources. The cells inside the hull (FB < 0) are ghost cells of a free-slip
        // wall: the velocity at a probe point outside, X = X_cell + (d + d_p) n (d the distance of the cell to the hull,
        // d_p = max(d, NHF_SRC_PROBE cells), n the outward normal of the level set), is extended into the cell with its
        // tangential part unchanged and its normal part linear through zero at the hull:
        //      V = u_p - (1 + d/d_p) (u_p.n) n
        // The velocity is then continuous across the hull, with no flow through it by symmetry: no source band with
        // its jump U n_x of the normal velocity (which the HLL flux turns into a numerical boundary layer: the stern
        // wave system of the Wigley hull is then only 20-50 % effective and the hump of the wave resistance is 30-40 %
        // too low). Cells deeper than NHF_SRC_MDEPTH cells (default 2.6) get the stream. Fixed hull. With
        // NHF_SRC_WATER=1 (the ghost cells are water for the fluxes) and X 13 1 (level of the pierced columns).
        if(tfmode==3)
        {
            static const double mdepth = getenv("NHF_SRC_MDEPTH")!=nullptr ? atof(getenv("NHF_SRC_MDEPTH")) : 2.6;
            static const double mprobe = getenv("NHF_SRC_PROBE")!=nullptr ? atof(getenv("NHF_SRC_PROBE")) : 1.5;
            static const double mrel = getenv("NHF_SRC_MREL")!=nullptr ? atof(getenv("NHF_SRC_MREL")) : 0.3;
            
            pgc->start4V(p,U,10);
            pgc->start4V(p,V,11);
            pgc->start4V(p,W,12);
            
            auto zcc = [&](int ii, int jj, int kk)
            {
                return p->ZP[kk+marge]*WL(ii,jj) + d->bed(ii,jj);
            };
            
            // value of F at the height zq in the column (ii,jj): linear between the cell centres, constant beyond them
            auto colval = [&](double *F, int ii, int jj, double zq)
            {
                int k0 = 0;
                
                while(k0<p->knoz-2 && zcc(ii,jj,k0+1)<zq)
                ++k0;
                
                const double z0 = zcc(ii,jj,k0);
                const double z1 = zcc(ii,jj,k0+1);
                const double wz = (z1>z0) ? MIN(MAX((zq-z0)/(z1-z0),0.0),1.0) : 0.0;
                
                return (1.0-wz)*F[IDX(ii,jj,k0)] + wz*F[IDX(ii,jj,k0+1)];
            };
            
            double Ffx3=0.0;
            
            LOOP
            {
                nhflow_srchull::q[IJK] = 0.0;
                
                const double fb = d->FB[IJK];
                
                if(fb>=0.0 || p->wet[IJ]!=1)
                continue;
                
                d->FHB[IJK] = 1.0;
                
                const double dl = (p->j_dir==0) ? p->DXN[IP] : 0.5*(p->DXN[IP] + p->DYN[JP]);
                const double dist = -fb;
                
                double ut = fk*nhflow_srchull::stream();
                double vt = 0.0;
                double wt = 0.0;
                
                if(dist<mdepth*dl)
                {
                    const double dfs = (d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
                    double gx = (d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1]) + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*dfs;
                    double gy = (p->j_dir==0) ? 0.0 : (d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1]) + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*dfs;
                    double gz = dfs/WL(i,j);
                    const double gn = sqrt(gx*gx + gy*gy + gz*gz);
                    
                    if(gn>1.0e-10)
                    {
                        gx /= gn; gy /= gn; gz /= gn;
                        
                        const double dp = MAX(dist,mprobe*dl);
                        const double xq = p->XP[IP] + (dist+dp)*gx;
                        const double yq = p->YP[JP] + (dist+dp)*gy;
                        const double zq = zcc(i,j,k) + (dist+dp)*gz;
                        
                        int i0 = i;
                        while(i0<i+6 && p->XP[i0+1+marge]<=xq) ++i0;
                        while(i0>i-6 && p->XP[i0+marge]>xq) --i0;
                        
                        int j0 = j;
                        if(p->j_dir==1)
                        {
                        while(j0<j+6 && p->YP[j0+1+marge]<=yq) ++j0;
                        while(j0>j-6 && p->YP[j0+marge]>yq) --j0;
                        }
                        
                        const bool ok = i0>=-2 && i0+1<=p->knox+1 && (p->j_dir==0 || (j0>=-2 && j0+1<=p->knoy+1));
                        
                        if(ok)
                        {
                            const double wx = MIN(MAX((xq - p->XP[i0+marge])/(p->XP[i0+1+marge] - p->XP[i0+marge]),0.0),1.0);
                            const double wy = (p->j_dir==0) ? 0.0 : MIN(MAX((yq - p->YP[j0+marge])/(p->YP[j0+1+marge] - p->YP[j0+marge]),0.0),1.0);
                            
                            auto at = [&](double *F)
                            {
                                double v = (1.0-wx)*(1.0-wy)*colval(F,i0,j0,zq) + wx*(1.0-wy)*colval(F,i0+1,j0,zq);
                                
                                if(p->j_dir==1)
                                v += (1.0-wx)*wy*colval(F,i0,j0+1,zq) + wx*wy*colval(F,i0+1,j0+1,zq);
                                
                                return v;
                            };
                            
                            const double up = at(U);
                            const double vp = at(V);
                            const double wp = at(W);
                            const double un = up*gx + vp*gy + wp*gz;
                            const double fac = 1.0 + dist/dp;
                            
                            ut = up - fac*un*gx;
                            vt = vp - fac*un*gy;
                            wt = wp - fac*un*gz;
                        }
                    }
                }
                
                // the ghost values depend on the flow, the flow (through the pressure projection) on the ghost values:
                // the full update every stage is unstable (first test, 2D ellipse and cylinder: velocities of 100 m/s
                // after 0.3 s); relaxed by NHF_SRC_MREL (default 0.3)
                FX[IJK] += mrel*(ut - U[IJK])/(alpha[iter]*p->dt);
                FY[IJK] += mrel*(vt - V[IJK])/(alpha[iter]*p->dt);
                FZ[IJK] += mrel*(wt - W[IJK])/(alpha[iter]*p->dt);
                
                Ffx3 += p->W1*(ut - U[IJK])/(alpha[iter]*p->dt)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);
            }
            
            return;
        }
        
        // NHF_SRC_VMATCH=1: the stream inside the hull follows the outer flow. With a uniform stream inside, the water
        // beside the hull (faster amidships, slower at the ends, +-10 % under the waves of the hull itself) slips along
        // it, and the mixing in the band is a numerical boundary layer: a wake of 10 N behind the Wigley hull, five
        // times its wave resistance. No flow through the hull holds for ANY velocity field V inside, with the sources
        // q = div(chi V) (blended field u = u_e + chi (V - u_e), div u = div(chi V) when u_e is tangent). Here
        // V = (V_x(x, layer), 0, 0), V_x the mean x-velocity of the water in the shell 1.5 psi < FB < 3 psi of the same
        // station and sigma layer (2D: of the station), relaxed in time (NHF_SRC_VREL per step, default 0.02) and
        // smoothed over three stations; q = d(chi V_x)/dx at constant z, of which chi dV_x/dx are sources inside.
        static const int vmode = getenv("NHF_SRC_VMATCH")!=nullptr ? atoi(getenv("NHF_SRC_VMATCH")) : 0;
        static const bool vmatch = vmode==1;
        
        // NOT WORKING YET (19 Nov 2026, 2D tests Cases/forcing_tests v2_ellipse, v2_cyl): stable for the ellipse but with a
        // wake (u_min/U 0.86-0.92, uniform stream inside: 0.95-0.98); for the cylinder no proper flow around it (speed
        // over the crest 1.26 U, wave 1.3 mm instead of 11.7) and bursts of 30 U. The forcing inside is 1000 times larger
        // than with the uniform stream: the divergence of V as the pressure projection takes it (node based) is not the
        // q = div(chi V) of the central differences here, where V is not smooth (jump to the stream at NHF_SRC_MDEPTH,
        // medial axis). Kept as a starting point only.
        // NHF_SRC_VMATCH=2: the field inside is the tangential extension of the outer flow. For every cell of the hull and
        // of the band (FB < psi) the velocity u_p at the probe point on its normal, at the distance d_p = NHF_SRC_PROBE
        // cells (1.5) outside the hull, gives the target
        //      V = u_p - (1 + d/d_p) (u_p.n) n          (d = -FB: > 0 inside, < 0 in the outer half of the band)
        // tangential part as at the probe, normal part linear through zero AT the hull. V is relaxed in time (NHF_SRC_VREL
        // per step), the cells inside the band (FB < -psi) are forced to V, and the sources are q = div(chi V) at
        // constant z, so that the pressure projection takes the field inside as it is and the blended-field argument
        // (no flow through the hull for any V) holds. V is nearly tangent to the hull: q is small (chi div V), there
        // is no band of sources with a jump U n_x of the normal velocity. Cells deeper than NHF_SRC_MDEPTH cells (2.6)
        // and cells where the level set has no clear normal (|grad FB| < 0.7: thin ends) keep the stream (U,0,0): there
        // q = U d(chi)/dx as before. Fixed hull.
        static const bool vmir = vmode==2;
        static std::vector<double> VM1,VM2,VM3;
        static std::vector<double> THW;                // NHF_SRC_THIN: weight of the forcing in x and z, centre cells of the rows without a core
        static const double mdepth2 = getenv("NHF_SRC_MDEPTH")!=nullptr ? atof(getenv("NHF_SRC_MDEPTH")) : 2.6;
        static const double mprobe2 = getenv("NHF_SRC_PROBE")!=nullptr ? atof(getenv("NHF_SRC_PROBE")) : 1.5;
        static const double vrel = getenv("NHF_SRC_VREL")!=nullptr ? atof(getenv("NHF_SRC_VREL")) : 0.02;
        static std::vector<double> Vst;
        static long vstep = -1;
        static int vcall = 0;
        const int nxg = p->gknox;
        const int nzl = p->knoz;
        
        if(vmatch)
        {
            if(int(Vst.size())!=nxg*nzl)
            Vst.assign(nxg*nzl,fk*nhflow_srchull::stream());
            
            if(long(p->count)!=vstep)
            {
            vstep = long(p->count);
            vcall = 0;
            }
            
            ++vcall;
            
            // second call of the step: the reforcing of the first stage, velocities after the pressure projection
            if(vcall==2)
            {
                std::vector<double> acc(2*nxg*nzl,0.0);
                
                LOOP
                if(p->wet[IJ]==1)
                {
                    const double ps = psiAt(i,j,k);
                    const double fb = d->FB[IJK];
                    
                    if(fb>1.5*ps && fb<3.0*ps)
                    {
                    const int n = (i+p->origin_i)*nzl + k;
                    const double w = p->DYN[JP]*p->DZN[KP]*WL(i,j);
                    acc[2*n] += w*U[IJK];
                    acc[2*n+1] += w;
                    }
                }
                
                pgc->globalsum(acc.data(),2*nxg*nzl);
                
                std::vector<double> tg(nxg*nzl,fk*nhflow_srchull::stream());
                
                for(int ig=0; ig<nxg; ++ig)
                {
                    if(p->j_dir==0)
                    {
                        double sa=0.0, sw=0.0;
                        for(int kk=0; kk<nzl; ++kk)
                        {sa += acc[2*(ig*nzl+kk)]; sw += acc[2*(ig*nzl+kk)+1];}
                        
                        if(sw>0.0)
                        for(int kk=0; kk<nzl; ++kk)
                        tg[ig*nzl+kk] = sa/sw;
                    }
                    
                    else
                    for(int kk=0; kk<nzl; ++kk)
                    if(acc[2*(ig*nzl+kk)+1]>0.0)
                    tg[ig*nzl+kk] = acc[2*(ig*nzl+kk)]/acc[2*(ig*nzl+kk)+1];
                }
                
                const double vlo = 0.5*fk*fabs(nhflow_srchull::stream());
                const double vhi = 1.5*fk*fabs(nhflow_srchull::stream());
                
                for(int ig=0; ig<nxg; ++ig)
                for(int kk=0; kk<nzl; ++kk)
                {
                    const int im = MAX(ig-1,0);
                    const int ip = MIN(ig+1,nxg-1);
                    double t = 0.25*tg[im*nzl+kk] + 0.5*tg[ig*nzl+kk] + 0.25*tg[ip*nzl+kk];
                    t = (nhflow_srchull::stream()>=0.0) ? MIN(MAX(t,vlo),vhi) : MAX(MIN(t,-vlo),-vhi);
                    Vst[ig*nzl+kk] += vrel*(t - Vst[ig*nzl+kk]);
                }
            }
        }
        
        if(vmir)
        {
            const int nn3 = p->imax*p->jmax*(p->kmax+2);
            
            if(int(VM1.size())!=nn3)
            {
            VM1.assign(nn3,fk*nhflow_srchull::stream());
            VM2.assign(nn3,0.0);
            VM3.assign(nn3,0.0);
            }
            
            static long vstep2 = -1;
            static int vcall2 = 0;
            
            if(long(p->count)!=vstep2)
            {
            vstep2 = long(p->count);
            vcall2 = 0;
            }
            
            ++vcall2;
            
            if(vcall2==2)
            {
                pgc->start4V(p,U,10);
                pgc->start4V(p,V,11);
                pgc->start4V(p,W,12);
                
                auto zcc = [&](int ii, int jj, int kk)
                {
                    return p->ZP[kk+marge]*WL(ii,jj) + d->bed(ii,jj);
                };
                
                auto colval = [&](double *F, int ii, int jj, double zq)
                {
                    int k0 = 0;
                    
                    while(k0<p->knoz-2 && zcc(ii,jj,k0+1)<zq)
                    ++k0;
                    
                    const double z0 = zcc(ii,jj,k0);
                    const double z1 = zcc(ii,jj,k0+1);
                    const double wz = (z1>z0) ? MIN(MAX((zq-z0)/(z1-z0),0.0),1.0) : 0.0;
                    
                    return (1.0-wz)*F[IDX(ii,jj,k0)] + wz*F[IDX(ii,jj,k0+1)];
                };
                
                LOOP
                {
                    const double ps = psiAt(i,j,k);
                    const double fb = d->FB[IJK];
                    
                    double t1 = fk*nhflow_srchull::stream();
                    double t2 = 0.0;
                    double t3 = 0.0;
                    
                    if(fb<1.2*ps && p->wet[IJ]==1)
                    {
                        const double dl = (p->j_dir==0) ? p->DXN[IP] : 0.5*(p->DXN[IP] + p->DYN[JP]);
                        const double dist = -fb;
                        
                        if(dist<mdepth2*dl)
                        {
                            const int kp = MIN(k+1,p->knoz-1);
                            const int km = MAX(k-1,0);
                            const double dfs = (kp>km) ? (d->FB[IDX(i,j,kp)] - d->FB[IDX(i,j,km)])/(p->ZP[kp+marge] - p->ZP[km+marge]) : 0.0;
                            double gx = (d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1]) + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*dfs;
                            double gy = (p->j_dir==0) ? 0.0 : (d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1]) + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*dfs;
                            double gz = dfs/WL(i,j);
                            const double gn = sqrt(gx*gx + gy*gy + gz*gz);
                            
                            if(gn>0.7)
                            {
                                gx /= gn; gy /= gn; gz /= gn;
                                
                                const double dp = mprobe2*dl;
                                const double xq = p->XP[IP] + (dist+dp)*gx;
                                const double yq = p->YP[JP] + (dist+dp)*gy;
                                const double zq = zcc(i,j,k) + (dist+dp)*gz;
                                
                                int i0 = i;
                                while(i0<i+6 && p->XP[i0+1+marge]<=xq) ++i0;
                                while(i0>i-6 && p->XP[i0+marge]>xq) --i0;
                                
                                int j0 = j;
                                if(p->j_dir==1)
                                {
                                while(j0<j+6 && p->YP[j0+1+marge]<=yq) ++j0;
                                while(j0>j-6 && p->YP[j0+marge]>yq) --j0;
                                }
                                
                                const bool ok = i0>=-2 && i0+1<=p->knox+1 && (p->j_dir==0 || (j0>=-2 && j0+1<=p->knoy+1));
                                
                                if(ok)
                                {
                                    const double wx = MIN(MAX((xq - p->XP[i0+marge])/(p->XP[i0+1+marge] - p->XP[i0+marge]),0.0),1.0);
                                    const double wy = (p->j_dir==0) ? 0.0 : MIN(MAX((yq - p->YP[j0+marge])/(p->YP[j0+1+marge] - p->YP[j0+marge]),0.0),1.0);
                                    
                                    auto at = [&](double *F)
                                    {
                                        double v = (1.0-wx)*(1.0-wy)*colval(F,i0,j0,zq) + wx*(1.0-wy)*colval(F,i0+1,j0,zq);
                                        
                                        if(p->j_dir==1)
                                        v += (1.0-wx)*wy*colval(F,i0,j0+1,zq) + wx*wy*colval(F,i0+1,j0+1,zq);
                                        
                                        return v;
                                    };
                                    
                                    // the probe must lie in the water (thin ends, concave corners)
                                    if(at(d->FB)>0.5*dp)
                                    {
                                    const double up = at(U);
                                    const double vp = at(V);
                                    const double wp = at(W);
                                    const double un = up*gx + vp*gy + wp*gz;
                                    const double fac = 1.0 + dist/dp;
                                    
                                    t1 = up - fac*un*gx;
                                    t2 = vp - fac*un*gy;
                                    t3 = wp - fac*un*gz;
                                    }
                                }
                            }
                        }
                    }
                    
                    VM1[IJK] += vrel*(t1 - VM1[IJK]);
                    VM2[IJK] += vrel*(t2 - VM2[IJK]);
                    VM3[IJK] += vrel*(t3 - VM3[IJK]);
                }
                
                pgc->start4V(p,VM1.data(),1);
                pgc->start4V(p,VM2.data(),1);
                pgc->start4V(p,VM3.data(),1);
            }
        }
        
        // NHF_SRC_VMATCH=3: tangent stream. The field inside the hull is the stream spread over the local beam,
        //      V = Us (1, eta b_x + c_x, 0),   eta = (y - c)/b,
        // b(x,z) the half beam and c(x,z) the centre of the hull in the row of cells (x, layer) (in 2D: the half height and
        // the centre in the column, V = Us (1, 0, eta b_x + c_x)). V is tangent to the hull, V.n = 0 at eta = +-1: the
        // sources q = div(chi V) = chi Us b_x/b + V.grad(chi) are spread over the beam (in every row 2 Us b_x, the thin-ship
        // strength) and the band holds none: the normal velocity is continuous across the hull, 0 on both sides. With
        // the uniform stream (NHF_SRC_FORCE=2 alone) it jumps by Us n_x over the two cells of the band, and what the
        // flux scheme dissipates there is taken from the water that moves along the hull (loss of head at the hull,
        // 30-45 % of the wave energy missing; Cases/NOTES_gece_2026-11-19.md).
        // Rows that the hull crosses more than once (and rows next to such rows only) keep the uniform stream: there
        // q = Us d(chi)/dx as before. b not below NHF_SRC_TBMIN cells (0.5), |V_y| <= NHF_SRC_TVMAX Us (1.0). Fixed hull.
        // The pressure projection takes the difference u - chi V (nhflow_srchull::cv1..3, nhflow_pjm::rhs): it is free of
        // divergence in the discrete form of the projection itself, and 0 where the cells are forced.
        static const bool vtan = vmode==3;
        static const double tbmin = getenv("NHF_SRC_TBMIN")!=nullptr ? atof(getenv("NHF_SRC_TBMIN")) : 0.5;
        static const double tvmax = getenv("NHF_SRC_TVMAX")!=nullptr ? atof(getenv("NHF_SRC_TVMAX")) : 1.0;
        const bool vfield = vmir || vtan;
        
        if(vtan)
        {
            const int nn3 = p->imax*p->jmax*(p->kmax+2);
            
            if(int(VM1.size())!=nn3)
            {
            VM1.assign(nn3,fk*nhflow_srchull::stream());
            VM2.assign(nn3,0.0);
            VM3.assign(nn3,0.0);
            }
            
            const bool tz = (p->j_dir==0);              // across the hull: z in 2D, y in 3D
            const int nrow = tz ? nxg : nxg*nzl;
            std::vector<double> acc(4*nrow,0.0);        // per row: sum of the points into the hull, their number, out of it, number
            
            auto zcc = [&](int ii, int jj, int kk)
            {
                return p->ZP[kk+marge]*WL(ii,jj) + d->bed(ii,jj);
            };
            
            LOOP
            if(p->wet[IJ]==1)
            {
                if(tz && k==p->knoz-1)
                continue;
                
                const double f0 = d->FB[IJK];
                const double f1 = tz ? d->FB[IJKp1] : d->FB[IJp1K];
                
                if((f0>=0.0)==(f1>=0.0))
                continue;
                
                const double c0 = tz ? zcc(i,j,k) : p->YP[JP];
                const double c1 = tz ? zcc(i,j,k+1) : p->YP[JP1];
                const double cc = c0 + (c1-c0)*f0/(f0-f1);
                const int n = tz ? (i+p->origin_i) : (i+p->origin_i)*nzl + k;
                
                if(f0>=0.0)
                {acc[4*n] += cc; acc[4*n+1] += 1.0;}
                else
                {acc[4*n+2] += cc; acc[4*n+3] += 1.0;}
            }
            
            pgc->globalsum(acc.data(),4*nrow);
            
            // 0: no hull in the row, 1: one stretch of hull, 2: anything else
            auto rstate = [&](int ig, int kk)
            {
                if(ig<0 || ig>=nxg)
                return 0;
                
                const int n = tz ? ig : ig*nzl + kk;
                const int ne = int(acc[4*n+1]+0.5);
                const int nx = int(acc[4*n+3]+0.5);
                
                if(ne==0 && nx==0)
                return 0;
                
                if(ne==1 && nx==1 && acc[4*n+2]>acc[4*n])
                return 1;
                
                return 2;
            };
            
            auto rowb = [&](int ig, int kk) {const int n = tz ? ig : ig*nzl + kk; return 0.5*(acc[4*n+2] - acc[4*n]);};
            auto rowc = [&](int ig, int kk) {const int n = tz ? ig : ig*nzl + kk; return 0.5*(acc[4*n+2] + acc[4*n]);};
            
            const double us = fk*nhflow_srchull::stream();
            
            // NHF_SRC_LOCAL (1: u and w, 2: u only, 3: w only): rows without a forced core take the local flow as the stream.
            // Where the hull is thinner than the band (half beam b < psi: the ends, the keel; on the Wigley grid of
            // 19 Nov 31 % of the sources) no cell is forced, the water passes through and the sources are those of
            // the thin ship, Us b_x per side, whatever the flow does there. The hull needs v = u b_x + w b_z with the
            // local velocity: less where the flow is slow (bow, stern), and less where it rises or dives past the sloping
            // bottom (the hull is no longer transparent to the vertical motion of its own waves). With V = (u_r, ., w_r),
            // u_r, w_r the chi-weighted mean velocity of the row (relaxed over NHF_SRC_LTAU, 0.2 s), the sources of the
            // row are div(chi V): 2 (u_r b_x + w_r b_z) + 2 b div(u_r, w_r). Rows with a core keep Us (there the
            // forcing makes the flow tangent whatever V is); linear change over one cell of half beam.
            // Found on Wigley, Fn 0.289: C_wp 1.428e-3 with psi 0.043 m, 1.503e-3 with 0.060 m (experiment 0.98e-3).
            static const int lmode = getenv("NHF_SRC_LOCAL")!=nullptr ? atoi(getenv("NHF_SRC_LOCAL")) : 0;
            static const double ltau = getenv("NHF_SRC_LTAU")!=nullptr ? atof(getenv("NHF_SRC_LTAU")) : 0.2;
            static std::vector<double> rowU,rowW,OMv;
            const bool lrow = lmode>0 && !tz;
            
            // NHF_SRC_THIN=1 (20 Nov): rows without a forced core get one. Where the half beam is below psi (+ one cell) no cell
            // lies deeper than psi inside the hull, nothing is forced, and the water inside is not the stream V that the sources
            // div(chi V) are made for: through the hull surface then flows (u_in - V).n (bow: out of the hull, stern: into it,
            // where the water at the centre plane is slower than Us). The cells of such a row next to its centre plane
            // (inside the hull, |y - c| < 0.75 dy) are forced to V in x and z; v stays free, so that the cells can take the
            // divergence the projection asks for (chi < 1 there). Weight: 1 for b <= psi, 0 from b = psi + dy on.
            static const bool thin = getenv("NHF_SRC_THIN")!=nullptr && atoi(getenv("NHF_SRC_THIN"))>0;
            
            if(thin && !tz && int(THW.size())!=nn3)
            THW.assign(nn3,0.0);
            
            if(lrow)
            {
                std::vector<double> lac(3*nrow,0.0);
                
                LOOP
                if(p->wet[IJ]==1)
                {
                    const double ch = chiAt(i,j,k);
                    
                    if(ch>0.02)
                    {
                        const int n = (i+p->origin_i)*nzl + k;
                        lac[3*n] += ch;
                        lac[3*n+1] += ch*U[IJK];
                        lac[3*n+2] += ch*W[IJK];
                    }
                }
                
                pgc->globalsum(lac.data(),3*nrow);
                
                if(int(rowU.size())!=nrow)
                {
                rowU.assign(nrow,us);
                rowW.assign(nrow,0.0);
                }
                
                if(int(OMv.size())!=nn3)
                OMv.assign(nn3,0.0);
                
                const double al = MIN(1.0,p->dt/MAX(ltau,1.0e-6));
                
                for(int n=0; n<nrow; ++n)
                if(lac[3*n]>1.0e-12)
                {
                    if(lmode!=3)
                    rowU[n] += al*(lac[3*n+1]/lac[3*n] - rowU[n]);
                    
                    if(lmode!=2)
                    rowW[n] += al*(lac[3*n+2]/lac[3*n] - rowW[n]);
                }
            }
            
            LOOP
            {
                double t2 = 0.0;
                const int ig = i+p->origin_i;
                
                double v1 = us;
                double v3 = 0.0;
                
                if(lrow)
                {
                    double om = 0.0;
                    
                    if(p->wet[IJ]==1)
                    {
                        om = (rstate(ig,k)!=1) ? 1.0 : MIN(MAX((psiAt(i,j,k) + p->DYN[JP] - rowb(ig,k))/p->DYN[JP],0.0),1.0);
                        v1 = us + om*(rowU[ig*nzl+k] - us);
                        v3 = om*rowW[ig*nzl+k];
                    }
                    
                    OMv[IJK] = om;
                }
                
                if(thin && !tz)
                {
                    double tw = 0.0;
                    
                    if(p->wet[IJ]==1 && d->FB[IJK]<0.0 && rstate(ig,k)==1 && fabs(p->YP[JP]-rowc(ig,k))<0.75*p->DYN[JP])
                    tw = MIN(MAX((psiAt(i,j,k) + p->DYN[JP] - rowb(ig,k))/p->DYN[JP],0.0),1.0);
                    
                    THW[IJK] = tw;
                }
                
                if(d->FB[IJK]<1.5*psiAt(i,j,k) && p->wet[IJ]==1 && rstate(ig,k)==1)
                {
                    const int sm = rstate(ig-1,k);
                    const int sp = rstate(ig+1,k);
                    
                    if(sm!=2 || sp!=2)
                    {
                        const double b0 = rowb(ig,k);
                        const double c0 = rowc(ig,k);
                        
                        const double bm = (sm==2) ? b0 : (sm==1 ? rowb(ig-1,k) : 0.0);
                        const double cm = (sm==1) ? rowc(ig-1,k) : c0;
                        const double xm = (sm==2) ? p->XP[IP] : p->XP[IM1];
                        const double bp = (sp==2) ? b0 : (sp==1 ? rowb(ig+1,k) : 0.0);
                        const double cp = (sp==1) ? rowc(ig+1,k) : c0;
                        const double xp = (sp==2) ? p->XP[IP] : p->XP[IP1];
                        
                        const double bx = (bp-bm)/(xp-xm);
                        const double cx = (cp-cm)/(xp-xm);
                        const double dl = tz ? p->DZN[KP]*WL(i,j) : p->DYN[JP];
                        const double eta = MIN(MAX(((tz ? zcc(i,j,k) : p->YP[JP]) - c0)/MAX(b0,tbmin*dl),-1.3),1.3);
                        
                        t2 = MIN(MAX(v1*(eta*bx + cx),-tvmax*fabs(us)),tvmax*fabs(us));
                    }
                }
                
                VM1[IJK] = v1;
                VM2[IJK] = tz ? 0.0 : t2;
                VM3[IJK] = tz ? t2 : v3;
            }
            
            // NHF_SRC_VFLUX=1: the stream inside the hull follows the water level inside it. In the columns that the hull
            // pierces the level eta is that of the waves along the hull, and a horizontal stream (U,0,0) below a
            // sloping level is no steady flow: the continuity of the column makes up for U d(eta)/dx with a flow
            // through the hull bottom and the band (Wigley, Fn 0.316: level slope 0.05, 5 % of U; the sources of the
            // bow are 11-30 % above those of the still hull at Fn 0.25-0.41, the sinks of the stern are not, the
            // rest leaves the hull along its length). The hull then yields to its own waves. Here, with z_b the hull
            // bottom of the column, T0 its still draft and d = T0 + eta the depth of the stream,
            //      V_x = Us T0/d,     V_z = (V_x eta_x + V_y eta_y)(z - z_b)/d
            // the stream carries the flux Us T0 of the still hull whatever the level and is tangent to the level:
            // the sources of every column are those of the still hull, their sum 0. Columns of the band beside the
            // hull take T0 and z_b of the nearest pierced column (three cells around).
            static const bool vflux = getenv("NHF_SRC_VFLUX")!=nullptr && atoi(getenv("NHF_SRC_VFLUX"))>0;
            static std::vector<double> SF,SB;
            
            if(vflux && !tz)
            {
                if(int(SF.size())!=nn3)
                {
                SF.assign(nn3,0.0);
                SB.assign(nn3,0.0);
                }
                
                const int kt = p->knoz-1;
                
                SLICELOOP4
                {
                    double sc = 0.0;
                    double zb = 0.0;
                    
                    if(p->wet[IJ]==1 && d->FB[IDX(i,j,kt)]<0.0)
                    {
                        int kb = kt;
                        
                        while(kb>0 && d->FB[IDX(i,j,kb-1)]<0.0)
                        --kb;
                        
                        if(kb>0)
                        {
                            const double f0 = d->FB[IDX(i,j,kb-1)];
                            const double f1 = d->FB[IDX(i,j,kb)];
                            zb = zcc(i,j,kb-1) + (zcc(i,j,kb) - zcc(i,j,kb-1))*f0/(f0-f1);
                            
                            const double dd = d->bed(i,j) + WL(i,j) - zb;
                            const double t0 = dd - d->eta(i,j);
                            
                            if(dd>1.0e-6 && t0>0.0)
                            sc = MIN(MAX(t0/dd,0.5),2.0);
                        }
                    }
                    
                    for(int kk=0; kk<p->knoz; ++kk)
                    {
                    SF[IDX(i,j,kk)] = sc;
                    SB[IDX(i,j,kk)] = zb;
                    }
                }
                
                pgc->start4V(p,SF.data(),1);
                pgc->start4V(p,SB.data(),1);
                
                LOOP
                if(p->wet[IJ]==1 && d->FB[IJK]<1.5*psiAt(i,j,k))
                {
                    double sc = SF[IDX(i,j,kt)];
                    double zb = SB[IDX(i,j,kt)];
                    
                    if(sc<=0.0)
                    {
                        double best = 1.0e20;
                        
                        for(int di=-3; di<=3; ++di)
                        for(int dj=-3; dj<=3; ++dj)
                        {
                            const double sn = SF[IDX(i+di,j+dj,kt)];
                            const double r2 = double(di*di + dj*dj);
                            
                            if(sn>0.0 && r2<best)
                            {
                            best = r2;
                            sc = sn;
                            zb = SB[IDX(i+di,j+dj,kt)];
                            }
                        }
                    }
                    
                    if(sc>0.0)
                    {
                        const double dd = MAX(d->bed(i,j) + WL(i,j) - zb,1.0e-6);
                        const double ex = (d->eta(i+1,j) - d->eta(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
                        const double ey = (d->eta(i,j+1) - d->eta(i,j-1))/(p->DYP[JP] + p->DYP[JM1]);
                        const double zr = MIN(MAX((zcc(i,j,k) - zb)/dd,0.0),1.2);
                        
                        // rows of the local stream (NHF_SRC_LOCAL) keep their own velocity
                        const double om = lrow ? OMv[IJK] : 0.0;
                        const double s1 = (om>0.0) ? 1.0 + (1.0-om)*(sc-1.0) : sc;
                        
                        VM1[IJK] *= s1;
                        VM2[IJK] *= s1;
                        VM3[IJK] = (om>0.0) ? (1.0-om)*(VM1[IJK]*ex + VM2[IJK]*ey)*zr + om*VM3[IJK] : (VM1[IJK]*ex + VM2[IJK]*ey)*zr;
                    }
                }
            }
            
            pgc->start4V(p,VM1.data(),1);
            pgc->start4V(p,VM2.data(),1);
            pgc->start4V(p,VM3.data(),1);
            
            std::vector<double> &c1 = nhflow_srchull::cv1;
            std::vector<double> &c2 = nhflow_srchull::cv2;
            std::vector<double> &c3 = nhflow_srchull::cv3;
            
            if(int(c1.size())!=nn3)
            {
            c1.assign(nn3,0.0);
            c2.assign(nn3,0.0);
            c3.assign(nn3,0.0);
            }
            
            LOOP
            {
                const double ch = (p->wet[IJ]==1) ? ramp*chiAt(i,j,k) : 0.0;
                c1[IJK] = ch*VM1[IJK];
                c2[IJK] = ch*VM2[IJK];
                c3[IJK] = ch*VM3[IJK];
            }
            
            pgc->start4V(p,c1.data(),1);
            pgc->start4V(p,c2.data(),1);
            pgc->start4V(p,c3.data(),1);
        }
        
        auto Vin = [&](int ii, int kk)
        {
            if(!vmatch)
            return fk*nhflow_srchull::stream();
            
            const int ig = MIN(MAX(ii+p->origin_i,0),nxg-1);
            
            return Vst[ig*nzl + MIN(MAX(kk,0),nzl-1)];
        };
        
        LOOP
        {
            if(tforce)
            {
                const double psf = psiAt(i,j,k);
                const double fbv = -d->FB[IJK];
                double chi = fbv>psf ? 1.0 : (fbv<-psf ? 0.0 : 0.5*(1.0 + fbv/psf + (1.0/PI)*sin(PI*fbv/psf)));
                
                d->FHB[IJK] = MIN(d->FHB[IJK] + chi, 1.0);
                
                if(tfmode==2)
                chi = fbv>psf ? 1.0 : 0.0;
                
                // NHF_SRC_THIN: x and z in the centre cells of the rows without a core
                const double chit = (tfmode==2 && !THW.empty()) ? MAX(chi,THW[IJK]) : chi;
                
                // any uniform velocity V inside gives no flow through the hull if the sources are -V n_x: (1 + k1) Us for
                // both (NHF_SRC_K). k1 then only sets the slip between inside and outside, not the boundary condition
                // during the start-up ramp (NHF_SRC_RAMP) the sources are ramp div(chi V): the stream the interior is held at
                // must have the divergence ramp div V, i.e. U0 + ramp (V - U0) with the uniform stream U0 (no hull at ramp 0)
                const double us0 = fk*nhflow_srchull::stream();
                const double vtx = vfield ? (ramp<1.0 ? us0 + ramp*(VM1[IJK]-us0) : VM1[IJK]) : Vin(i,k);
                const double vty = vfield ? (ramp<1.0 ? ramp*VM2[IJK] : VM2[IJK]) : 0.0;
                const double vtz = vfield ? (ramp<1.0 ? ramp*VM3[IJK] : VM3[IJK]) : 0.0;

                FX[IJK] += chit*(vtx - U[IJK])/(alpha[iter]*p->dt);
                FY[IJK] += chi*(vty - V[IJK])/(alpha[iter]*p->dt);
                FZ[IJK] += chit*(vtz - W[IJK])/(alpha[iter]*p->dt);

                // force of the forcing on the water (momentum budget)
                Ffx += p->W1*chit*(vtx - U[IJK])/(alpha[iter]*p->dt)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);
                Ffz += p->W1*chit*(vtz - W[IJK])/(alpha[iter]*p->dt)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);
            }
            
            if(keep)
            {
                const double qv = nhflow_srchull::q[IJK];
                const double dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);
                Fx -= p->W1*qv*((1.0-nhflow_srchull::qv())*U[IJK] + nhflow_srchull::qv()*fk*nhflow_srchull::stream())*dV;
                Qp += MAX(qv,0.0)*dV;
                Qn += MIN(qv,0.0)*dV;
                continue;
            }
            
            const double psi = (psifix>0.0) ? psifix : ((p->j_dir==0) ? p->A526*p->DXN[IP] : p->A526*0.5*(p->DXN[IP] + p->DYN[JP]));
            
            auto Hs = [&](double fb)
            {
                if(-fb>psi)
                return 1.0;
                
                if(-fb<-psi)
                return 0.0;
                
                return 0.5*(1.0 + (-fb)/psi + (1.0/PI)*sin((PI*(-fb))/psi));
            };
            
            // sigma = -Us n_x on the hull (n outward) and grad(chi) = -n delta: q = +Us d(chi)/dx, sources at the bow
            // d(chi)/dx at constant z: the neighbours i-1, i+1 of a sigma layer are at other heights where the free
            // surface is inclined, d/dx|z = d/dx|sigma + sig_x d/dsigma (as the normal of the direct forcing). Without
            // the metric term the hull bottom gets sources Us d(chi)/dz times the slope of the layer: on Wigley the
            // sources and the sinks both grew by 17 % while the waves developed (0.1123 -> 0.132 m^3/s), on KCS by 20 %.
            // NHF_SRC_NOMETRIC=1: without it (all runs before 18 Nov 15:20)
            const double ps0 = (anisomode==3) ? psiAt(i,j,k) : 0.0;
            const double dcdx = (anisomode==3) ?
                                (chiP(i+1,j,k,ps0) - chiP(i-1,j,k,ps0))/(p->DXP[IP] + p->DXP[IM1])
                              + (nometric ? 0.0 : 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(chiP(i,j,MIN(k+1,p->knoz-1),ps0) - chiP(i,j,MAX(k-1,0),ps0))/(p->DZP[KP] + p->DZP[KM1]))
                              :
                                aniso ?
                                (chiAt(i+1,j,k) - chiAt(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1])
                              + (nometric ? 0.0 : 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(chiAt(i,j,MIN(k+1,p->knoz-1)) - chiAt(i,j,MAX(k-1,0)))/(p->DZP[KP] + p->DZP[KM1]))
                              :
                                (Hs(d->FB[Ip1JK]) - Hs(d->FB[Im1JK]))/(p->DXP[IP] + p->DXP[IM1])
                              + (nometric ? 0.0 : 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(Hs(d->FB[IJKp1]) - Hs(d->FB[IJKm1]))/(p->DZP[KP] + p->DZP[KM1]));
            double qv = ramp*fk*Us*dcdx;
            
            // NHF_SRC_VMATCH: q = d(chi V_x)/dx at constant z (fixed hull; a moving one: - u_body d(chi)/dx)
            if(vmatch)
            {
                const double cip = aniso ? chiAt(i+1,j,k) : Hs(d->FB[Ip1JK]);
                const double cim = aniso ? chiAt(i-1,j,k) : Hs(d->FB[Im1JK]);
                const double ckp = aniso ? chiAt(i,j,MIN(k+1,p->knoz-1)) : Hs(d->FB[IJKp1]);
                const double ckm = aniso ? chiAt(i,j,MAX(k-1,0)) : Hs(d->FB[IJKm1]);
                
                const double dcv = (cip*Vin(i+1,k) - cim*Vin(i-1,k))/(p->DXP[IP] + p->DXP[IM1])
                                 + (nometric ? 0.0 : 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(ckp*Vin(i,k+1) - ckm*Vin(i,k-1))/(p->DZP[KP] + p->DZP[KM1]));
                
                qv = ramp*(dcv - u_fb(0)*dcdx);
            }
            
            // NHF_SRC_VMATCH=2, 3: q = div(chi V) at constant z
            if(vfield)
            {
                const int kp = MIN(k+1,p->knoz-1);
                const int km = MAX(k-1,0);
                const double dsg = (kp>km) ? (p->ZP[kp+marge] - p->ZP[km+marge]) : 1.0;
                
                auto cv = [&](int ii, int jj, int kk, std::vector<double> &F)
                {
                    return chiAt(ii,jj,kk)*F[IDX(ii,jj,kk)];
                };
                
                const double d1s = (cv(i,j,kp,VM1) - cv(i,j,km,VM1))/dsg;
                const double d2s = (cv(i,j,kp,VM2) - cv(i,j,km,VM2))/dsg;
                const double d3s = (cv(i,j,kp,VM3) - cv(i,j,km,VM3))/dsg;
                
                double dv = (cv(i+1,j,k,VM1) - cv(i-1,j,k,VM1))/(p->DXP[IP] + p->DXP[IM1]) + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*d1s
                          + d3s/WL(i,j);
                
                if(p->j_dir==1)
                dv += (cv(i,j+1,k,VM2) - cv(i,j-1,k,VM2))/(p->DYP[JP] + p->DYP[JM1]) + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*d2s;
                
                qv = ramp*dv;
            }
            
            // NHF_SRC_WL (factor, 1: on): the water line term of the Neumann-Kelvin representation of a surface piercing
            // hull, 4 pi phi = int_S sigma G dS - (Us^2/g) oint_WL sigma n_x G dy (Brard 1972). The sources alone are a
            // hull with a free surface inside (the water flows through it); the line sources -(Us^2/g) sigma n_x t_y on
            // the water line (sinks at the bow, sources at the stern) make up for it in the linear theory (Wigley,
            // Fr 0.267: C_w 1.14e-3 -> 0.83e-3, experiment 0.85e-3; Cases/WIGLEY/nk_line.py). With sigma = -Us n_x,
            // t_y = n_x/|n_h| and grad(chi) = -n delta the flux per horizontal area of the column is
            //      -(Us^3/g) (d(chi)/dx)^3 / |grad(chi)|^2        in the top layer
            if(wlf!=0.0 && k==p->knoz-1)
            {
                const double cy = (p->j_dir==0) ? 0.0 : (Hs(d->FB[IJp1K]) - Hs(d->FB[IJm1K]))/(p->DYP[JP] + p->DYP[JM1]);
                const double cz = (Hs(d->FB[IJK]) - Hs(d->FB[IJKm1]))/(p->DZP[KM1]*WL(i,j));
                const double g2 = dcdx*dcdx + cy*cy + cz*cz;
                
                if(g2>1.0e-12)
                qv -= wlf*fk*(Us*Us*Us/fabs(p->W22))*dcdx*dcdx*dcdx/g2/(p->DZN[KP]*WL(i,j));
            }
            
            // NHF_SRC_FB: the strength from the flow itself. Across the sheet the normal velocity jumps by sigma,
            // u_n(d) = u_e,n - sigma chi(d); no flow through the hull, u_e,n = 0, needs sigma = -u_n/chi. Relaxed
            // (stable for every chi, error halved per update for a flat sheet), with q = sigma |grad chi| and
            // u_n |grad chi| = -u.grad(chi):
            //      q <- (1 - chi) q + (u - u_body).grad(chi)
            // First update from q = 0: Hogner, q = u.grad(chi). Fixed point: div u = u.grad(ln chi), the water
            // outside does not cross the hull whatever the flow inside; no k1.
            if(fbk)
            {
                const double chi = Hs(d->FB[IJK]);
                double uf,vf,wf;
                target_velocity(p,false,uf,vf,wf);
                
                const double dcdy = (p->j_dir==0) ? 0.0 : (Hs(d->FB[IJp1K]) - Hs(d->FB[IJm1K]))/(p->DYP[JP] + p->DYP[JM1]);
                
                double dcdz;
                if(k==0)
                dcdz = (Hs(d->FB[IJKp1]) - chi)/(p->DZP[KP]*WL(i,j));
                else
                if(k==p->knoz-1)
                dcdz = (chi - Hs(d->FB[IJKm1]))/(p->DZP[KM1]*WL(i,j));
                else
                dcdz = (Hs(d->FB[IJKp1]) - Hs(d->FB[IJKm1]))/((p->DZP[KP] + p->DZP[KM1])*WL(i,j));
                
                // under-relaxed (NHF_SRC_FBW, default 0.1): with the full update the cells at the outer edge of the band
                // (chi -> 0) are neutrally stable and the strength grows on the collocated grid (2D ellipse: in 0.2 s)
                const double qn = (1.0-chi)*nhflow_srchull::q[IJK] + (U[IJK]-uf)*dcdx + (V[IJK]-vf)*dcdy + fbz*(W[IJK]-wf)*dcdz;
                
                qv = first ? fk*Us*dcdx : (1.0-fbw)*nhflow_srchull::q[IJK] + fbw*qn;
            }
            
            if(p->wet[IJ]!=1)
            qv = 0.0;
            
            nhflow_srchull::q[IJK] = qv;
            
            const double dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);
            
            Fx -= p->W1*qv*((1.0-nhflow_srchull::qv())*U[IJK] + nhflow_srchull::qv()*fk*nhflow_srchull::stream())*dV;
            Qp += MAX(qv,0.0)*dV;
            Qn += MIN(qv,0.0)*dV;
        }
        
        first = false;
        
        // NHF_SRC_DUMP: the sources of the first stage as text (x y z q dV), one file per rank, for checks
        static bool dumped = getenv("NHF_SRC_DUMP")==nullptr;
        
        if(!dumped)
        {
            dumped = true;
            
            char name[200];
            snprintf(name,200,"./REEF3D_NHFLOW_6DOF/REEF3D_srchull_q_%i.dat",p->mpirank);
            std::ofstream dump(name);
            dump.precision(10);
            
            LOOP
            if(fabs(nhflow_srchull::q[IJK])>1.0e-12)
            dump<<p->pos_x()<<" "<<p->pos_y()<<" "<<p->pos_z()<<" "<<nhflow_srchull::q[IJK]<<" "<<p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j)<<" "<<d->FB[IJK]<<endl;
        }
        
        {
            Fx = pgc->globalsum(Fx);
            Qp = pgc->globalsum(Qp);
            Qn = pgc->globalsum(Qn);
            Ffx = pgc->globalsum(Ffx);
            Ffz = pgc->globalsum(Ffz);
            
            // momentum budget: the water gets the forcing, the momentum of the source water (+q u, -Fx here) and the
            // momentum of the level reset of the pierced columns (X 13 1). In a steady flow the opposite of their sum
            // is the momentum the flow has lost at the hull: the resistance R_x, whatever happens to the waves further
            // away. R_z: the vertical force on the hull without the buoyancy at rest is not in this budget.
            // The kernel runs twice per stage, for the forcing before and the reforcing after the pressure projection:
            // both forcings count, the source and the level terms once per stage. Every stage counts with its weight in
            // the step: the forcing F_s adds alpha_s dt F_s to the output of stage s, of which the step keeps
            // 1/2, 1/2 (SSP-RK2) or 1/6, 1/6, 2/3 (SSP-RK3) of dt F_s. (The first stage alone is not the mean: forcing
            // and level reset of the two stages trade places from step to step, Wigley with X 13 1.)
            // One line per time step, written when the next step begins.
            // (x13_adt is set only on the ranks with pierced columns: alpha dt of the stage here)
            const double Fx13 = pgc->globalsum(nhflow_srchull::x13_mx)/(alpha[iter]*p->dt);
            const double Q13 = pgc->globalsum(nhflow_srchull::x13_dv)/(alpha[iter]*p->dt);
            
            static long step_last = -1;
            static int stage_last = -1;
            static double t_acc=0.0, Fx_acc=0.0, Qp_acc=0.0, Qn_acc=0.0, Ff_acc=0.0, F13_acc=0.0, Q13_acc=0.0, Ffs_acc[3]={0.0,0.0,0.0};
            const bool newstep = (long(p->count)!=step_last);
            const bool newstage = newstep || iter!=stage_last;
            const bool write = newstep && step_last>=0;
            const double t_w=t_acc, Fx_w=Fx_acc, Qp_w=Qp_acc, Qn_w=Qn_acc, Ff_w=Ff_acc, F13_w=F13_acc, Q13_w=Q13_acc;
            const double Ffs_w[3]={Ffs_acc[0],Ffs_acc[1],Ffs_acc[2]};
            
            const double wrk2[3]={0.5,0.5,0.0};
            const double wrk3[3]={1.0/6.0,1.0/6.0,2.0/3.0};
            const double ws = (p->A510==3) ? wrk3[MIN(iter,2)] : wrk2[MIN(iter,2)];
            
            if(newstep)
            {
                step_last = long(p->count);
                t_acc = p->simtime; Qp_acc = Qp; Qn_acc = Qn; Fx_acc = 0.0; Ff_acc = 0.0; F13_acc = 0.0; Q13_acc = 0.0;
                Ffs_acc[0] = Ffs_acc[1] = Ffs_acc[2] = 0.0;
            }
            
            if(newstage)
            {
                stage_last = iter;
                Fx_acc += ws*Fx;
                F13_acc += ws*Fx13;
                Q13_acc += ws*Q13;
            }
            
            Ff_acc += ws*Ffx;
            Ffs_acc[MIN(iter,2)] += Ffx;
            
            const double Rx = -(Ff_w - Fx_w + F13_w);
            
            if(p->mpirank==0)
            {
                static std::ofstream out;
                
                if(!out.is_open())
                {
                    out.open("./REEF3D_NHFLOW_6DOF/REEF3D_srchull.dat");
                    out<<"time \t Fx_Lagally \t Q_sources \t Q_sinks \t F_forcing_x \t F_level_x \t R_x \t F_forcing_stage1 \t F_forcing_stage2 \t Q_level"<<endl;
                }
                
                if(write)
                out<<t_w<<" \t "<<Fx_w<<" \t "<<Qp_w<<" \t "<<Qn_w<<" \t "<<Ff_w<<" \t "<<F13_w<<" \t "<<Rx<<" \t "<<Ffs_w[0]<<" \t "<<Ffs_w[1]<<" \t "<<Q13_w<<endl;
            }
        }
        
        return;
    }
    
    // X 15 11: direct forcing with cell and face fractions and the body in the pressure projection
    // (the wavemaker, X 10 4, with the forcing of X 15 4)
    if(p->X15==11 && !wavemaker)
    {
    forcing_kernel_fraction(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter,wavemaker);
    return;
    }
    
    // Calculate forcing fields
    double H, uf, vf, wf;
    
    double du,dv,dw,udotn;
    double beta;   // 0.0 = free slip (tangential preserved), 1.0 = no slip
    
    if(p->X15==0)
    beta=1.0;
    
    if(p->X15>=1)
    beta=0.0;
    
    
    if(p->X15==0)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        target_velocity(p,wavemaker,uf,vf,wf);
         
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        FX[IJK] += H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += H*(wf - W[IJK])/(alpha[iter]*p->dt);
    }
    
    // X 15 3: as X 15 1, but only the cells inside the body are forced. A body moving through the
    // grid at speed (towed or self-propelled ship) otherwise pushes the water in the forcing band ahead
    // of the bow up to the body velocity; that water passes under the hull and is carried along with it.
    // X 15 5: as X 15 1 with the sharp switch at FB = 0 (full forcing for FB <= 0, normal-weighted outside),
    // the forcing before the shell ramp was introduced (9ce616dab)
    // X 15 6 / 7: the forcing of X 15 1 / 3; the body faces are reflective walls in the HLL fluxes (nhflow_HLL.cpp)
    if(p->X15==1 || p->X15==3 || p->X15==5 || p->X15==6 || p->X15==7)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        target_velocity(p,wavemaker,uf,vf,wf);
         
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        if((p->X15==3 || p->X15==7) && d->FB[IJK]>=0.0)
        continue;
        
    // Normal vectors calculation 
        if(p->geo_curv==1)
        fbnormal_curv(p,d,nx,ny);
        else
        {
		nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
            - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
        
        if(p->j_dir==0)
        ny = 0.0;
        
        if(p->j_dir==1)
		ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
            - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
        }
            
		nz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP]*WL(i,j) + p->DZP[KM1]*WL(i,j));

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;
        
        
        // full forcing inside the body, normal-weighted forcing (tangential velocity kept) outside,
        // ramped over a thin shell |FB| < dsh around the surface. A sharp switch at FB = 0 jumps
        // for cell centres on or next to the surface (a flat side through a row of cell centres):
        // the smallest rotation moves the row partly inside and partly outside, and the jump in
        // the forcing kicks a symmetric hull into sway and yaw.
        // X 15 5: the sharp switch at FB = 0 (2D cases and grid-aligned hulls, where the kick does not
        // arise; reproduces the floating-body motions of the code before the ramp)
        double sb;
        
        if(p->X15==5)
        sb = d->FB[IJK]<=0.0 ? 1.0 : 0.0;
        
        else
        {
        const double dsh = 0.25*Hpsi_nhflow(p,d);
        const double tb = MIN(MAX(0.5*(1.0 - d->FB[IJK]/dsh), 0.0), 1.0);
        sb = tb*tb*(3.0 - 2.0*tb);
        }
        
        FX[IJK] += (sb + (1.0-sb)*fabs(nx))*H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += (sb + (1.0-sb)*fabs(ny))*H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += (sb + (1.0-sb)*fabs(nz))*H*(wf - W[IJK])/(alpha[iter]*p->dt);
    
    }
    
    
    if(p->X15==2)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        target_velocity(p,wavemaker,uf,vf,wf);
        
        du = uf - U[IJK];
        dv = vf - V[IJK];
        dw = wf - W[IJK];
        
        if(d->FB[IJK] < 0.0)
        {
            // interior of the body: enforce the full rigid-body velocity
            FX[IJK] += H*du/(alpha[iter]*p->dt);
            FY[IJK] += H*dv/(alpha[iter]*p->dt);
            FZ[IJK] += H*dw/(alpha[iter]*p->dt);
        }
         
        else
        {
            // Normal vectors calculation 
            if(p->geo_curv==1)
            fbnormal_curv(p,d,nx,ny);
            else
            {
            nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
                - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
                
            if(p->j_dir==0)
            ny = 0.0;
            
            if(p->j_dir==1)
            ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
                - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
            }
                
            nz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP]*WL(i,j) + p->DZP[KM1]*WL(i,j));

            norm = sqrt(nx*nx + ny*ny + nz*nz);
                    
            if(norm > 1.0e-10)
            {
            nx /= norm;
            ny /= norm;
            nz /= norm;
            
            udotn = nx*du + ny*dv + nz*dw;
            
            FX[IJK] += H*(beta*du + (1.0-beta)*udotn*nx)/(alpha[iter]*p->dt);
            FY[IJK] += H*(beta*dv + (1.0-beta)*udotn*ny)/(alpha[iter]*p->dt);
            FZ[IJK] += H*(beta*dw + (1.0-beta)*udotn*nz)/(alpha[iter]*p->dt);
            }
            
            else
            {
            // degenerate gradient (flat level set): fall back to full forcing
            FX[IJK] += H*du/(alpha[iter]*p->dt);
            FY[IJK] += H*dv/(alpha[iter]*p->dt);
            FZ[IJK] += H*dw/(alpha[iter]*p->dt);
            }
        }
    
    }
    
    // X 15 4: normal projection outside the body (as X 15 2), ramped over the thin shell |FB| < dsh of X 15 1:
    //   F = H [sb du + (1 - sb)(du.n) n] / (alpha dt CPORNH)
    // sb = 1 inside the body (full rigid-body velocity), 0 outside (only the normal velocity is imposed,
    // the tangential velocity is kept). Unlike the component weighting |n_i| of X 15 1, the projection does
    // not depend on the orientation of the hull to the grid (rolled hull, bilges): the tangential velocity
    // outside is not damped. Unlike the sharp switch of X 15 2 at FB = 0, the forcing stays continuous when
    // the hull moves across cell centres. On grid-aligned faces X 15 1 and 4 are identical.
    // The forcing is divided by CPORNH as in nhflow_forcing::solid_forcing, so that the body velocity is
    // fully imposed also inside a porous structure (without porous media CPORNH = 1).
    if(p->X15==4 || (p->X15==11 && wavemaker))
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        if(H<=0.0)
        continue;
        
        target_velocity(p,wavemaker,uf,vf,wf);
        
        du = uf - U[IJK];
        dv = vf - V[IJK];
        dw = wf - W[IJK];
        
        const double dsh = 0.25*Hpsi_nhflow(p,d);
        const double tb = MIN(MAX(0.5*(1.0 - d->FB[IJK]/dsh), 0.0), 1.0);
        double sb = tb*tb*(3.0 - 2.0*tb);
        
        nx=ny=nz=0.0;
        udotn=0.0;
        
        if(sb<1.0)
        {
            // Normal vectors calculation 
            if(p->geo_curv==1)
            fbnormal_curv(p,d,nx,ny);
            else
            {
            nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
                - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
                
            if(p->j_dir==0)
            ny = 0.0;
            
            if(p->j_dir==1)
            ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
                - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
            }
                
            nz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP]*WL(i,j) + p->DZP[KM1]*WL(i,j));

            norm = sqrt(nx*nx + ny*ny + nz*nz);
            
            if(norm > 1.0e-10)
            {
            nx /= norm;
            ny /= norm;
            nz /= norm;
            
            udotn = nx*du + ny*dv + nz*dw;
            }
            
            // degenerate gradient (flat level set): fall back to full forcing
            else
            sb = 1.0;
        }
        
        const double fac = H/(alpha[iter]*p->dt*CPORNH);
        
        FX[IJK] += fac*(sb*du + (1.0-sb)*udotn*nx);
        FY[IJK] += fac*(sb*dv + (1.0-sb)*udotn*ny);
        FZ[IJK] += fac*(sb*dw + (1.0-sb)*udotn*nz);
    }
   
    pgc->start5V(p,d->FHB,50);
}
    
void sixdof_obj_nhflow::fbnormal_curv(lexer *p, fdm_nhf *d, double &gx, double &gy)
{
    // curvilinear grid (G 50): the horizontal part of -grad(FB) at the cell (i,j,k), the gradient along the sigma
    // level by the chain rule on the cell centres, then the sigma correction (sigx, sigy are physical)
    double fx,fy;
    p->geo_grad(i,j,d->FB[Ip1JK],d->FB[Im1JK],d->FB[IJp1K],d->FB[IJm1K],fx,fy);
    const double fs = (d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);

    gx = -fx - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*fs;
    gy = -fy - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*fs;
}

double sixdof_obj_nhflow::Hpsi_nhflow(lexer *p, fdm_nhf *d)
{
    // half width of the smoothed Heaviside of the direct forcing (as in Hsolidface_nhflow)
    if(p->geo_curv==1)
    return p->A526*p->geo_width(i,j);

    if(p->j_dir==0)
    return p->A526*p->DXN[IP];
    
    return p->A526*0.5*(p->DXN[IP] + p->DYN[JP]);
}

double sixdof_obj_nhflow::Hsolidface_nhflow(lexer *p, fdm_nhf *d, int aa, int bb, int cc)
{
    double psi, H, phival_fb;
    
    if(p->j_dir==0)
    psi = p->A526*(1.0/1.0)*(p->DXN[IP] + 0.0*p->DZN[KP]*p->WL[IJ]);
    
    if(p->j_dir==1)
    psi = p->A526*(1.0/2.0)*(p->DXN[IP] + p->DYN[JP] + 0.0*p->DZN[KP]*p->WL[IJ]);

    // curvilinear grid (G 50): the physical mean width of the cell
    if(p->geo_curv==1)
    psi = p->A526*p->geo_width(i,j);


    // Construct solid heaviside function
    phival_fb = d->FB[IJK];
    
    if(-phival_fb > psi)
    H = 1.0;

    if(-phival_fb < -psi)
    H = 0.0;

    if(fabs(phival_fb)<=psi)
    H = 0.5*(1.0 + (-phival_fb)/psi + (1.0/PI)*sin((PI*(-phival_fb))/psi));
    
    return H;
}

// X 15 11: direct forcing with the body in the pressure projection
//
//   fraction    kap of every cell: the part of the cell inside the body, from the signed distance FB of the
//               cell centre and its gradient, kap = 1/2 - FB/|(nx dx, ny dy, nz dz)| limited to 0..1 (exact
//               for the faces of the hull along the grid). It changes linearly while the surface moves through
//               the cell: no switching of cells, the volume of the body on the grid is that of the hull.
//   forcing     u = u* + tau (u_b - u*), tau = max(2 kap - 1, 0): only where the cell centre is inside the
//               hull, the rigid-body velocity in the cells of the body; continuous in the position of the hull.
//               The cells with their centre in the water are not forced: the body acts on them through the
//               fluxes of its faces and the pressure. A forcing of their normal velocity with the fraction
//               works against the projection, which has made the water flow around the hull: the pressure then
//               holds a part that only balances the forcing, and the pressure loads a mean force (a box free in
//               waves of 2 s drifted downwave with 4 times the drift force of its reflected and transmitted
//               waves, the mean force on the fixed box of 2 m in waves of 1.4 s was 63 % of that of its waves,
//               now 93 %). The reaction of the forcing as a load on the body did not help: the heave decay and
//               the mean force on the fixed box got worse.
//               Open: a box free in waves of 2 s (H 0.04 m) still gains about 0.002 m/s2 downwave, 3 times the
//               drift force of its waves (NHFLOW Floating validation 04).
//   projection  open fractions of the links of the Poisson equation: 1 - kap of the cell for the vertical
//               link through it (W of the cell), 1 - the fraction of the face from the mean FB of its two
//               cells for the horizontal links (d->MBX, MBY, MBZ; nhflow_membrane_beta.h, link mode). The
//               pressure correction then does not move water through the body, the pressure below the hull
//               follows from the normal velocity of the hull (its added mass is in the pressure), and the
//               Dirichlet condition at the top of a column covered by the body does not reach the water below.
//               The nodes inside the body hold no pressure (nhflow_membrane_row).
//   fluxes      the faces of the body carry the body through their closed part, in the divergence of the
//               projection, the continuity flux and the momentum fluxes (nhflow_bodyflux.h, nhflow_HLL.cpp)
//   level       the water level of the columns covered by the body is that of the surrounding free surface
//               (nhflow_fsf_body, as X 13 1)
//   The velocity is not forced again after the projection (the projected field keeps its continuity).
//   Loads: 6DOF_obj_nhflow_forces_stl.cpp, coupling: 6DOF_obj_nhflow_addedmass.cpp
//
//   Not yet: mesh refinement (G 1), curvilinear grids (G 50), the semi-implicit time step (A 510 11-13), several
//   bodies; checkpoints restart with the added mass and the acceleration of the first step after the restart
double sixdof_obj_nhflow::cellfraction_nhflow(lexer *p, fdm_nhf *d, slice &WL, double &ex, double &ey, double &ez)
{
    double gx,gy,gz;
    
    if(p->geo_curv==1)
    fbnormal_curv(p,d,gx,gy);
    else
    {
    gx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
        - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
    
    gy = 0.0;
    
    if(p->j_dir==1)
    gy = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
        - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
    }
    
    const double wl = MAX(WL(i,j),1.0e-6);
    
    gz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP]*wl + p->DZP[KM1]*wl);
    
    const double g = sqrt(gx*gx + gy*gy + gz*gz);
    
    const double dx = p->geo_curv==1 ? p->geo_width(i,j) : p->DXN[IP];
    const double dy = p->geo_curv==1 ? p->geo_width(i,j) : p->DYN[JP];
    const double dz = p->DZN[KP]*wl;
    
    // away from the surface (FB at the limit of its band) or a degenerate gradient: inside or outside
    if(g<0.1)
    {
    ex = ey = ez = 0.0;
    return d->FB[IJK]<0.0 ? 1.0 : 0.0;
    }
    
    ex = gx/g;
    ey = gy/g;
    ez = gz/g;
    
    // width of the cell in the direction of the normal. With the sum of the sides, |nx| dx + |ny| dy + |nz| dz (a
    // plane through a box cell), the cells at an edge of the hull, where the normal of the distance points along the
    // diagonal, stay open by a part although they lie inside: nodes with nearly closed links, whose pressure is
    // not determined by the water. With the length of (nx dx, ny dy, nz dz) the fraction is exact for the faces
    // along the grid and 1 in the cells at the edges
    const double w = sqrt(ex*ex*dx*dx + (p->j_dir==1 ? ey*ey*dy*dy : 0.0) + ez*ez*dz*dz);
    
    const double kap = MIN(MAX(0.5 - d->FB[IJK]/w, 0.0), 1.0);
    
    // no cells and links that are open or closed by less than 1 %: their pressure nodes are nearly without water
    return kap>0.99 ? 1.0 : (kap<0.01 ? 0.0 : kap);
}

void sixdof_obj_nhflow::fraction_ini(lexer *p, fdm_nhf *d)
{
    const int n = p->imax*p->jmax*(p->kmax+2);
    
    double **mb[] = {&d->MBETA,&d->MCHI,&d->MRCX,&d->MRCY,&d->MBX,&d->MBY,&d->MBZ};
    
    if(FNX==nullptr)
    {
    p->Darray(FNX,n);
    p->Darray(FNY,n);
    p->Darray(FNZ,n);
    }
    
    if(d->FBUX==nullptr)
    {
    p->Darray(d->FBUX,n);
    p->Darray(d->FBUY,n);
    p->Darray(d->FBWZ,n);
    p->Iarray(d->FBNEAR,p->imax*p->jmax);
    }
    
    // open fractions of the links (the mobilities of nhflow_membrane_beta.h, link mode) and the face corrections
    // of the last projection, 1 and 0 away from the body
    for(int q=0; q<7; ++q)
    if(*mb[q]==nullptr)
    {
        p->Darray(*mb[q],n);
        
        for(int m=0; m<n; ++m)
        (*mb[q])[m] = (q==1 || q==2 || q==3) ? 0.0 : 1.0;
    }
    
    d->MPROJ = 1;
}

void sixdof_obj_nhflow::forcing_kernel_fraction(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, int iter, bool wavemaker)
{
    // smallest open fraction of a link: the Poisson equation stays regular inside the body
    const double bmin = 1.0e-4;
    
    double uf,vf,wf,ex,ey,ez;
    
    fraction_ini(p,d);
    
    // acceleration of the time step (the loads of prescribed degrees of freedom, force_calc_stl)
    if(iter==0 && !fraction_reforce && p->simtime!=tfd_)
    {
        if(tfd_>=0.0 && p->simtime>tfd_)
        afd_ = (u_fb - ufd_)/(p->simtime - tfd_);
        
        ufd_ = u_fb;
        tfd_ = p->simtime;
    }
    
    LOOP
    {
        ex = ey = ez = 0.0;
        
        const double kap = (p->wet[IJ]==1) ? cellfraction_nhflow(p,d,WL,ex,ey,ez) : 0.0;
        
        FNX[IJK] = ex;
        FNY[IJK] = ey;
        FNZ[IJK] = ez;
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + kap, 1.0);
        d->MBETA[IJK] = MAX(1.0 - kap, bmin);
        d->MBZ[IJK] = d->MBETA[IJK];
        
        // forcing only where the cell centre is inside the hull
        if(kap<=0.5 || fraction_reforce)
        continue;
        
        target_velocity(p,wavemaker,uf,vf,wf);
        
        const double du = uf - U[IJK];
        const double dv = vf - V[IJK];
        const double dw = wf - W[IJK];
        
        const double tau = 2.0*kap - 1.0;
        const double fac = 1.0/(alpha[iter]*p->dt*CPORNH);

        // EXPERIMENT (ahmet_dev_ship, environment variable NHF_SLIP, not an input parameter yet): the cut cells
        // with their centre inside the hull (0.5 < kap < 1) hold water that slides along the hull. Forcing all
        // components towards the body velocity takes its tangential momentum (a no-slip wall for an inviscid
        // flow: a wake behind the body). 1: only the normal component is forced, 2: the cut cells are not forced
        static const int slip = getenv("NHF_SLIP")!=nullptr ? atoi(getenv("NHF_SLIP")) : 0;

        if(slip>=1 && kap<1.0)
        {
            if(slip==2)
            continue;

            const double dn = du*ex + dv*ey + dw*ez;

            FX[IJK] += fac*tau*dn*ex;
            FY[IJK] += fac*tau*dn*ey;
            FZ[IJK] += fac*tau*dn*ez;
            continue;
        }

        FX[IJK] += fac*tau*du;
        FY[IJK] += fac*tau*dv;
        FZ[IJK] += fac*tau*dw;
    }
    
    pgc->start5V(p,d->FHB,50);
    pgc->start4V(p,d->MBETA,1);
    pgc->start4V(p,d->MBZ,1);
    pgc->start4V(p,FNX,1);
    pgc->start4V(p,FNY,1);
    pgc->start4V(p,FNZ,1);
    
    // open fractions of the horizontal links
    LOOP
    {
        d->MBX[IJK] = 1.0;
        d->MBY[IJK] = 1.0;
        d->FBUX[IJK] = 0.0;
        d->FBUY[IJK] = 0.0;
        d->FBWZ[IJK] = 0.0;
        
        if(p->wet[IJ]==0)
        continue;
        
        const double dz = p->DZN[KP]*MAX(WL(i,j),1.0e-6);
        
        // normal velocity of the body at the centres of the x and y faces, vertical velocity at the cell centre
        {
        const double ry = p->pos_y() - c_(1);
        const double rx = p->pos_x() - c_(0);
        
        d->FBWZ[IJK] = u_fb(2) + u_fb(3)*ry - u_fb(4)*rx;
        
        if(p->wet[Ip1J]==1)
        d->FBUX[IJK] = u_fb(0) + u_fb(4)*(0.5*(p->ZSP[IJK] + p->ZSP[Ip1JK]) - c_(2)) - u_fb(5)*ry;
        
        if(p->j_dir==1 && p->wet[IJp1]==1)
        d->FBUY[IJK] = u_fb(1) + u_fb(5)*rx - u_fb(3)*(0.5*(p->ZSP[IJK] + p->ZSP[IJp1K]) - c_(2));
        }
        
        if(p->wet[Ip1J]==1 && p->flag4[Ip1JK]>0)
        d->MBX[IJK] = linkfraction_nhflow(d->FB[IJK],d->FB[Ip1JK],d->MBETA[IJK],d->MBETA[Ip1JK],
                                          (p->j_dir==1 ? 0.5*fabs(FNY[IJK]+FNY[Ip1JK])*p->DYN[JP] : 0.0) + 0.5*fabs(FNZ[IJK]+FNZ[Ip1JK])*dz, bmin);
        
        if(p->j_dir==1 && p->wet[IJp1]==1 && p->flag4[IJp1K]>0)
        d->MBY[IJK] = linkfraction_nhflow(d->FB[IJK],d->FB[IJp1K],d->MBETA[IJK],d->MBETA[IJp1K],
                                          0.5*fabs(FNX[IJK]+FNX[IJp1K])*p->DXN[IP] + 0.5*fabs(FNZ[IJK]+FNZ[IJp1K])*dz, bmin);
    }
    
    // columns near the body
    for(int q=0; q<p->imax*p->jmax; ++q)
    d->FBNEAR[q] = 0;
    
    LOOP
    if(d->MBETA[IJK]<1.0)
    d->FBNEAR[IJ] = 1;
    
    pgc->gcsl_start4Vint(p,d->FBNEAR,1);
    
    for(int qn=0; qn<nhflow_bodyflux::nnear; ++qn)
    {
        SLICELOOP4
        if(d->FBNEAR[IJ]==0 && (d->FBNEAR[Ip1J]==1 || d->FBNEAR[Im1J]==1 
           || (p->j_dir==1 && (d->FBNEAR[IJp1]==1 || d->FBNEAR[IJm1]==1))))
        d->FBNEAR[IJ] = 2;
        
        SLICELOOP4
        if(d->FBNEAR[IJ]==2)
        d->FBNEAR[IJ] = 1;
        
        pgc->gcsl_start4Vint(p,d->FBNEAR,1);
    }
    
    // no pressure at the nodes inside the body (between two cells of the body)
    LOOP
    if(k>0 && d->MBZ[IJK]<=bmin && d->MBZ[IJKm1]<=bmin)
    d->P[FIJK] = 0.0;
    
    pgc->start4V(p,d->MBX,1);
    pgc->start4V(p,d->MBY,1);
    pgc->start4V(p,d->FBUX,1);
    pgc->start4V(p,d->FBUY,1);
    pgc->start4V(p,d->FBWZ,1);
}

double sixdof_obj_nhflow::linkfraction_nhflow(double fa, double fb, double ba, double bb, double w, double bmin)
{
    // the face between two cells: FB at its centre from the two cell centres, w the change of FB over the face
    // (the components of the normal in the plane of the face times its sides). The open fraction of the face is the
    // part with FB > 0. A face parallel to the surface (w = 0) would switch when the surface passes it: the face
    // is open by at most twice the open fraction of its cells, so that it closes linearly while the surface
    // moves through the last half cell in front of it.
    if(ba>=1.0 && bb>=1.0)
    return 1.0;
    
    const double fm = 0.5*(fa + fb);
    const double af = w>1.0e-12 ? MIN(MAX(0.5 + fm/w, 0.0), 1.0) : (fm>0.0 ? 1.0 : 0.0);
    
    const double b = MIN(af, MIN(2.0*MIN(ba,bb), 1.0));
    
    return b>0.99 ? 1.0 : (b<0.01 ? bmin : b);
}
