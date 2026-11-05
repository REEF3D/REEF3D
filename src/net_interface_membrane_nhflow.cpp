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

#include"net_interface.h"
#include"net_membrane.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include"nhflow_membrane_beta.h"
#include"nhflow_thinbody.h"
#include"vrans_definitions.h"
#include<mpi.h>
#include<fstream>
#include<sstream>
#include<string>

// ctrl.txt: X 330 1 (A 520 1 or 2; moving membranes, structure rigid|flexible, switch A 520 2 to 1, see
// driver_logic_nhflow.cpp). membrane.dat (read by rank 0, broadcast):
//
//   # comment
//   membrane box       x0 x1 y0 y1 z_bottom z_top     starts a new membrane (in 2D y0, y1 are ignored)
//   membrane cylinder  xc yc R z_bottom z_top
//   membrane cylcone   xc yc R z_tip z_cone z_top  cylinder (z_cone to z_top) on a cone bottom (z_tip to z_cone), as the
//                                      closed flexible cage of Strand et al. (2013); the cone is the bag floor
//   name        text                   optional label
//   resistance  R_n [R_t]              hydraulic resistance [m/s], leakage u_n = (dp/rho)/R_n; default 1e4;
//                                      R_t default 0 (fixed membrane), R_n (moving membrane: the layer moves with it)
//   thickness   delta                  half width of the smeared layer [m]; default 1.5 max(dx,dy,dz)
//   mesh        h                      target triangle edge length [m]; default min(dx,dy) (fixed), max(dx,dy,dz) (moving),
//                                      1.5 delta (coupling iterated)
//   fill        dh                     initial inner water level above the outside level [m]; default 0
//   filling     lambda                 filling level V_water/V_bag (V_bag below the still water level): the inner
//                                      level starts lower by (1 - lambda) V_bag / A_waterplane, the flexible bag
//                                      then deflates to the inner volume; overrides fill
//   drain       T                      with filling: the bag starts full and the missing water is pumped out of its
//                                      interior over T [s] (quasi-static deflation instead of the sudden inner level drop)
//   print       dt                     vtp output interval [s] (REEF3D_NHFLOW_Membrane_VTP, with a .pvd
//                                      collection); default: NHFLOW print control P 30 / P 20; 0: off
//   floorpressure 0|1|3                static pressure below the floor: 0 uniform head difference,
//                                      1 local excess head (comparison only), 3 shape of the discrete
//                                      free-surface ramp, running average (tau), frozen at t = 2 tau;
//                                      default 3 (fixed membrane), 0 (moving membrane)
//   tau         t                      averaging time of floorpressure 3 [s]; default 2
//   projections n                      projection passes per stage (default 1: Rhie-Chow continuity flux)
//   mobility    layer|link             pressure coupling of the membrane. layer (default): the implicit porous factor
//                                      1/(1 + a K_n H) is the (isotropic) mobility of all cells of the smeared layer;
//                                      the layer fluid can then not be moved along the membrane by pressure and is held
//                                      with it (R_t = R_n for a moving membrane). link: only the links between two cell
//                                      centres (nodes, vertically) on opposite sides of the membrane get the porous-jump
//                                      mobility 1/(1 + a R_n/l); the layer fluid moves freely along the membrane
//                                      (R_t default 0), the normal resistance of the layer stays. Loads: momentum taken
//                                      out by the forcing + pressure difference across the blocked links.
//                                      Link mode: projections 1, coupling staggered (flexible), delta >= cell size
//               sharp                  only the links crossing the membrane are blocked (mobility 1/(1 + a R_n/l),
//                                      default R_n 1e5 m/s), no layer forcing: the cells on either side are free fluid.
//                                      Wall fluxes at the blocked faces, wall velocity at the blocked links in the
//                                      projection, cut cells below/above a floor, hydrostatic head below closed floors
//                                      from the outer free surface (no floorpressure); see nhflow_thinbody.h.
//                                      A 520 1, projections 1, structure fixed, rigid or flexible (flexible: staggered
//                                      coupling with an added mass, see addedmass, bending, smoothing). Loads: pressure
//                                      jump across the blocked links (no shear)
//   poisson     0|1                    membrane mobility in the pressure Poisson equation; default 1
//                                      (0 only to demonstrate the splitting leakage of the projection)
//
//   structure   fixed|rigid|flexible   fixed (default); rigid: moves with the floating body (X 10);
//                                      flexible: mass-spring membrane, top edge attached to the floating
//                                      body like the nets (X 320), or held in place without X 10
//                                      (coupling iterated, the default; with coupling staggered only a fixed or
//                                      prescribed collar motion, X 10 2 / X 11 2)
//   mass        m                      fabric mass per area [kg/m^2]; default 1
//   density     rho                    fabric density [kg/m^3] (buoyancy); default 1300
//   stiffness   Et                     membrane stiffness E t [N/m]; default 5e5
//   damping     zeta                   damping ratio of the edge dampers; default 0.1
//   compression f                      edge stiffness in compression as a fraction of E t (wrinkling); default 0.01
//   smoothing   n                      mobility sharp, structure flexible: smoothing passes of the membrane geometry that
//                                      the fluid sees (node displacement and velocity, neighbour average, removes folds
//                                      at the scale of the mesh); default 2, 0: off
//   bending     EI                     bending stiffness [N m] of the flexible membrane, a regularisation against folds
//                                      at the scale of the mesh (quadratic hinge model, flat rest state, Bergou et al.
//                                      2006); default 1e-3 with mobility sharp (the sharp walls need a membrane that is
//                                      smooth at the cell scale), 0 otherwise
//   sinker      w                      submerged weight along the floor edge [N/m]; default 0
//   attach      z                      nodes at or above z are attached; default the top edge z_top
//   addedmass   m_a                    mobility sharp, structure flexible: added mass per area [kg/m^2] of the
//                                      staggered coupling (the fabric is integrated with m + m_a against the load of
//                                      the last fluid step; stable for m_a above the added mass of the water, which
//                                      is ~ rho times the bag size); default 2 rho (L + H), L the radius (cylinder) or
//                                      half the smaller width (box), H the depth of the bag bottom below the still
//                                      water level. The extra inertia slows the dynamics, not the equilibrium
//   bodyaddedmass M                    added mass [kg] of the stabilised coupling to the floating body
//                                      (translation); default 2 rho V_bag (water of the bag below the still
//                                      water level) for a rigid membrane, 0 for a flexible one (its top
//                                      edge is coupled implicitly); 0: off. Rotations of a rigid bag: added
//                                      inertia 4 rho x inertia of that water about the centre of gravity,
//                                      scaled with M / (2 rho V_bag) when M is given
//   collar      D m EA EI [Cd [Ca]]    flexible membrane: its top edge is a floating pipe ring (no floating body,
//                                      X 10 0): diameter D [m], mass m [kg/m], axial and bending stiffness EA [N],
//                                      EI [N m^2], Morison drag / added-mass coefficients (default 1, 1).
//                                      Buoyancy from the local free surface, Froude-Krylov, added mass and drag
//                                      normal to the pipe axis; corotational bending (net_membrane_collar.cpp).
//                                      The collar centre line is the top edge of the bag (z_top). 3D only
//   mooring     xa ya za k T0          linear mooring spring [N/m], pretension T0 [N], from the anchor (xa,ya,za)
//                                      to the collar node nearest to it (horizontal distance); several lines
//   coupling    staggered|iterated [rtol [n]]
//                                      fluid-structure coupling of a flexible membrane. iterated (default): the
//                                      projection of every RK stage is repeated until the node velocities of fluid
//                                      and structure agree to rtol (default 1e-3), at most n iterations (default
//                                      50); IQN-ILS (net_membrane_coupling.cpp). staggered: the structure is
//                                      advanced once per time step with the implicit porous damper (stable, but an
//                                      extra inertia ~ rho R_n dt per area makes the dynamics depend on dt; not for
//                                      a freely floating collar)
//   couplingtol rtol [atol]            iterated: relative tolerance, absolute tolerance [m/s] (default 1e-5, rms)
//   couplingiter n                     iterated: maximum iterations per stage (default 50)
//   couplingreuse n                    iterated: converged stages whose secant information is reused (default 8)
//   couplingrelax w                    iterated: relaxation of the first iteration without history (default 0.5)
//   couplingrobin f                    iterated: scale of the Robin preconditioner (default 16); the converged
//                                      solution does not depend on it, the tolerances are divided by f
//   couplingcolumns n                  iterated: maximum number of IQN-ILS columns (default 100)
//   couplingfilter eps                 iterated: IQN-ILS QR filter (default 1e-2)
//   couplingqn ils|imvj                iterated: quasi-Newton update, IQN-ILS with reused columns (default) or
//                                      IQN-IMVJ with the inverse Jacobian carried over (n x n per RK stage)
//   couplinglog 0|1                    iterated: residual of every coupling iteration on screen (default 0)
//
// Parameter lines apply to the most recent 'membrane' line.

void net_interface::membrane_ini_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    string content;
    int len=0;

    if(p->mpirank==0)
    {
        ifstream f("membrane.dat");

        if(!f)
        len=-1;
        else
        {
            stringstream ss;
            ss<<f.rdbuf();
            content = ss.str();
            len = (int)content.size();
        }
    }

    MPI_Bcast(&len,1,MPI_INT,0,pgc->mpi_comm);

    if(len<0)
    {
        if(p->mpirank==0)
        cout<<"\n!!! X 330: membrane.dat not found in the case directory !!!\n"<<endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    content.resize(len);
    if(len>0)
    MPI_Bcast(&content[0],len,MPI_CHAR,0,pgc->mpi_comm);

    // parse
    vector<membrane_param> mp;
    istringstream is(content);
    string line;
    int lineno=0;
    bool error=false;

    while(getline(is,line))
    {
        ++lineno;

        const size_t hash = line.find('#');
        if(hash!=string::npos)
        line = line.substr(0,hash);

        istringstream ls(line);
        string key;

        if(!(ls>>key))
        continue;

        if(key=="membrane")
        {
            string shape;
            membrane_param m;
            ls>>shape;

            if(shape=="box")
            {
                m.shape=1;
                if(!(ls>>m.x0>>m.x1>>m.y0>>m.y1>>m.zb>>m.zt))
                error=true;
            }
            else if(shape=="cylinder")
            {
                m.shape=2;
                if(!(ls>>m.xc>>m.yc>>m.R>>m.zb>>m.zt))
                error=true;
            }
            else if(shape=="cylcone")
            {
                m.shape=3;
                if(!(ls>>m.xc>>m.yc>>m.R>>m.zb>>m.zc>>m.zt) || !(m.zb<m.zc && m.zc<m.zt))
                error=true;
            }
            else
            error=true;

            m.name = "membrane"+to_string(mp.size());
            mp.push_back(m);
        }
        else if(mp.empty())
        error=true;
        else if(key=="name")
        {
            if(!(ls>>mp.back().name))
            error=true;
        }
        else if(key=="resistance")
        {
            if(!(ls>>mp.back().Rn))
            error=true;
            ls>>mp.back().Rt;
            mp.back().Rngiven=1;
        }
        else if(key=="thickness")
        {
            if(!(ls>>mp.back().delta))
            error=true;
        }
        else if(key=="mesh")
        {
            if(!(ls>>mp.back().h))
            error=true;
        }
        else if(key=="fill")
        {
            if(!(ls>>mp.back().fill))
            error=true;
        }
        else if(key=="smoothing")
        {
            if(!(ls>>mp.back().smooth) || mp.back().smooth<0)
            error=true;
        }
        else if(key=="bending")
        {
            if(!(ls>>mp.back().EI) || mp.back().EI<0.0)
            error=true;
        }
        else if(key=="compression")
        {
            if(!(ls>>mp.back().compr) || mp.back().compr<0.0 || mp.back().compr>1.0)
            error=true;
        }
        else if(key=="drain")
        {
            if(!(ls>>mp.back().drain) || mp.back().drain<0.0)
            error=true;
        }
        else if(key=="filling")
        {
            if(!(ls>>mp.back().filling) || mp.back().filling<=0.0 || mp.back().filling>1.0)
            error=true;
        }
        else if(key=="floorpressure")
        {
            if(!(ls>>mp.back().floorp))
            error=true;
        }
        else if(key=="tau")
        {
            if(!(ls>>mp.back().tau) || mp.back().tau<=0.0)
            error=true;
        }
        else if(key=="projections")
        {
            if(!(ls>>mp.back().projections) || mp.back().projections<1)
            error=true;
        }
        else if(key=="mobility")
        {
            string mode;
            ls>>mode;
            
            if(mode=="layer")
            mp.back().link=0;
            else if(mode=="link")
            mp.back().link=1;
            else if(mode=="sharp")
            mp.back().link=2;
            else
            error=true;
        }
        else if(key=="poisson")
        {
            if(!(ls>>mp.back().poisson))
            error=true;
        }
        else if(key=="print")
        {
            if(!(ls>>mp.back().printdt))
            error=true;
        }
        else if(key=="structure")
        {
            string st;
            ls>>st;
            
            if(st=="fixed")
            mp.back().structure=0;
            else if(st=="rigid")
            mp.back().structure=1;
            else if(st=="flexible")
            mp.back().structure=2;
            else
            error=true;
        }
        else if(key=="mass")
        {
            if(!(ls>>mp.back().mA) || mp.back().mA<=0.0)
            error=true;
        }
        else if(key=="density")
        {
            if(!(ls>>mp.back().rhom) || mp.back().rhom<=0.0)
            error=true;
        }
        else if(key=="stiffness")
        {
            if(!(ls>>mp.back().EA) || mp.back().EA<=0.0)
            error=true;
        }
        else if(key=="damping")
        {
            if(!(ls>>mp.back().zeta) || mp.back().zeta<0.0)
            error=true;
        }
        else if(key=="sinker")
        {
            if(!(ls>>mp.back().sinker))
            error=true;
        }
        else if(key=="collar")
        {
            // collar D m EA EI [Cd [Ca]]
            mp.back().collar=1;
            
            if(!(ls>>mp.back().cD>>mp.back().cm>>mp.back().cEA>>mp.back().cEI)
               || mp.back().cD<=0.0 || mp.back().cm<0.0 || mp.back().cEA<=0.0 || mp.back().cEI<0.0)
            error=true;
            
            double v;
            if(ls>>v)
            {
                mp.back().cCd=v;
                
                if(ls>>v)
                mp.back().cCa=v;
            }
            
            if(mp.back().cCd<0.0 || mp.back().cCa<0.0)
            error=true;
        }
        else if(key=="mooring")
        {
            // mooring xa ya za k T0: linear spring from the anchor to the collar node nearest to it
            array<double,5> m;
            
            if(!(ls>>m[0]>>m[1]>>m[2]>>m[3]>>m[4]) || m[3]<0.0 || m[4]<0.0)
            error=true;
            else
            mp.back().moor.push_back(m);
        }
        else if(key=="addedmass")
        {
            if(!(ls>>mp.back().ma) || mp.back().ma<0.0)
            error=true;
        }
        else if(key=="bodyaddedmass")
        {
            if(!(ls>>mp.back().Mbody) || mp.back().Mbody<0.0)
            error=true;
        }
        else if(key=="coupling")
        {
            string st;
            ls>>st;
            
            if(st=="staggered")
            mp.back().coupling=0;
            else if(st=="iterated")
            {
                mp.back().coupling=1;
                
                double t;
                int n;
                
                if(ls>>t)
                {
                    mp.back().crtol=t;
                    
                    if(ls>>n)
                    mp.back().citer=n;
                }
            }
            else
            error=true;
            
            if(mp.back().crtol<=0.0 || mp.back().citer<1)
            error=true;
        }
        else if(key=="couplingtol")
        {
            if(!(ls>>mp.back().crtol) || mp.back().crtol<=0.0)
            error=true;
            
            ls>>mp.back().catol;
        }
        else if(key=="couplingiter")
        {
            if(!(ls>>mp.back().citer) || mp.back().citer<1)
            error=true;
        }
        else if(key=="couplingreuse")
        {
            if(!(ls>>mp.back().creuse) || mp.back().creuse<0)
            error=true;
        }
        else if(key=="couplingrobin")
        {
            if(!(ls>>mp.back().crobin) || mp.back().crobin<=0.0)
            error=true;
        }
        else if(key=="couplingcolumns")
        {
            if(!(ls>>mp.back().ccols) || mp.back().ccols<1)
            error=true;
        }
        else if(key=="couplingfilter")
        {
            if(!(ls>>mp.back().cfilt) || mp.back().cfilt<=0.0)
            error=true;
        }
        else if(key=="couplingqn")
        {
            string st;
            ls>>st;
            
            if(st=="ils")
            mp.back().cqn=0;
            else if(st=="imvj")
            mp.back().cqn=1;
            else
            error=true;
        }
        else if(key=="couplinglog")
        {
            if(!(ls>>mp.back().clog))
            error=true;
        }
        else if(key=="couplingrelax")
        {
            if(!(ls>>mp.back().crelax) || mp.back().crelax<=0.0 || mp.back().crelax>1.0)
            error=true;
        }
        else if(key=="attach")
        {
            if(!(ls>>mp.back().zattach))
            error=true;
        }
        else
        error=true;

        if(error)
        {
            if(p->mpirank==0)
            cout<<"\n!!! membrane.dat, line "<<lineno<<": cannot read '"<<line<<"' !!!\n"<<endl;
            MPI_Abort(pgc->mpi_comm,1);
        }
    }

    if(mp.empty())
    {
        if(p->mpirank==0)
        cout<<"\n!!! X 330: no membrane defined in membrane.dat !!!\n"<<endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    // mobility field for the pressure Poisson equation, 1 away from membranes
    if(d->MBETA==nullptr)
    p->Darray(d->MBETA,p->imax*p->jmax*(p->kmax+2));

    if(d->MCHI==nullptr)
    p->Darray(d->MCHI,p->imax*p->jmax*(p->kmax+2));

    if(d->MRCX==nullptr)
    p->Darray(d->MRCX,p->imax*p->jmax*(p->kmax+2));

    if(d->MRCY==nullptr)
    p->Darray(d->MRCY,p->imax*p->jmax*(p->kmax+2));

    for(int qn=0; qn<p->imax*p->jmax*(p->kmax+2); ++qn)
    {
    d->MBETA[qn]=1.0;
    d->MCHI[qn]=0.0;
    }

    // link mode: Rhie-Chow flux of the single projection (the converged wide divergence would average the cell
    // velocities across the blocked links)
    bool link=false, sharp=false;
    
    for(size_t m=0; m<mp.size(); ++m)
    if(mp[m].link>=1)
    {
        if(mp[m].link==2)
        sharp=true;
        
        link=true;
        
        if(mp[m].projections>1 && p->mpirank==0)
        cout<<"Membrane "<<m<<": mobility link uses projections 1"<<endl;
        
        mp[m].projections=1;
    }
    
    if(link)
    {
        p->Darray(d->MBX,p->imax*p->jmax*(p->kmax+2));
        p->Darray(d->MBY,p->imax*p->jmax*(p->kmax+2));
        p->Darray(d->MBZ,p->imax*p->jmax*(p->kmax+2));
        
        for(int qn=0; qn<p->imax*p->jmax*(p->kmax+2); ++qn)
        d->MBX[qn]=d->MBY[qn]=d->MBZ[qn]=1.0;
    }
    
    for(size_t m=0; m<mp.size(); ++m)
    d->MPROJ = MAX(d->MPROJ, mp[m].projections);

    // sharp mode: thin-body service (wall fluxes, projection right-hand side, cut cells, head below the floors)
    if(sharp && d->thinbody==nullptr)
    d->thinbody = new nhflow_thinbody(p,d,pgc);
    
    for(size_t m=0; m<mp.size(); ++m)
    {
        pmem.push_back(new net_membrane(m,mp[m]));
        
        if(mp[m].link==2)
        pmem.back()->tb_ = d->thinbody;
        
        pmem.back()->initialize_nhflow(p,d,pgc);
        pmem.back()->fill_nhflow(p,d,pgc);
    }
    
    // sharp mode: blocked links of the initial state, so that the first fluxes already see the walls
    if(d->thinbody!=nullptr)
    {
        membrane_links_nhflow(p,d,pgc,1.0);
        
        const int nl = d->thinbody->nlower(p,pgc);
        
        if(p->mpirank==0)
        cout<<"X 330 sharp: "<<nl<<" cells below closed floors"<<endl;
    }
}

void net_interface::membrane_links_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha)
{
    // 1. mobility field beta (also builds the membrane cell maps for this stage)
    for(int qn=0; qn<p->imax*p->jmax*(p->kmax+2); ++qn)
    d->MBETA[qn]=1.0;
    
    if(d->MBX!=nullptr)
    for(int qn=0; qn<p->imax*p->jmax*(p->kmax+2); ++qn)
    d->MBX[qn]=d->MBY[qn]=d->MBZ[qn]=1.0;
    
    if(d->thinbody!=nullptr)
    d->thinbody->begin(p);

    for(auto m : pmem)
    m->mobility_nhflow(p,d,pgc,alpha);

    pgc->start4V(p,d->MBETA,1);
    
    if(d->MBX!=nullptr)
    {
    pgc->start4V(p,d->MBX,1);
    pgc->start4V(p,d->MBY,1);
    pgc->start4V(p,d->MBZ,1);
    }
    
    if(d->thinbody!=nullptr)
    d->thinbody->finish(p,d,pgc);
}

void net_interface::membrane_forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha,
                                            double *UH, double *VH, double *WH, slice &WL)
{
    // 0. membrane positions for this stage (moving membranes)
    for(auto m : pmem)
    m->kinematics_nhflow(p,d,pgc);
    
    // 1. mobility field beta, blocked links (also builds the membrane cell maps for this stage)
    membrane_links_nhflow(p,d,pgc,alpha);
    
    // sharp mode: cells that crossed the body with the moving sigma grid
    if(d->thinbody!=nullptr)
    d->thinbody->side_change(p,d,UH,VH,WH,WL);

    // 2. static overpressure of the bag below its floor (prescribed pressure, see net_membrane)
    for(int qn=0; qn<p->imax*p->jmax*(p->kmax+2); ++qn)
    d->MCHI[qn]=0.0;

    for(auto m : pmem)
    m->static_pressure_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);

    pgc->start4V(p,d->MCHI,1);

    // 3. implicit porous-jump forcing
    //    Incremental pressure scheme (A 520 2): the predictor already contains the old non-hydrostatic
    //    pressure gradient, -a C/rho G P^n. For the projection to stay consistent, P^n has to act with
    //    the same mobility as the pressure correction (beta_h horizontal, beta vertical), not with the
    //    implicit forcing factor (I + a A)^-1 in the layer or 1 next to it. It is taken out before the
    //    forcing and put back with the correction mobility afterwards; without this the pressure
    //    increments accumulate in the layer and the scheme diverges.
    const bool incremental = (p->A520==2);

    if(incremental)
    membrane_pgrad(p,d,alpha,UH,VH,WH,WL,1);

    for(auto m : pmem)
    m->forcing_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);

    if(incremental)
    membrane_pgrad(p,d,alpha,UH,VH,WH,WL,-1);
    
    // sharp mode: vertical velocity of the cut cells (below / above a floor)
    if(d->thinbody!=nullptr)
    d->thinbody->cut_forcing(p,d,UH,VH,WH,WL);
}

void net_interface::membrane_pgrad(lexer *p, fdm_nhf *d, double alpha, double *UH, double *VH, double *WH, slice &WL, int mode)
{
    // Incremental scheme (A 520 2): the predictor contains the old pressure gradient -a C/rho G P^n
    // (nhflow_pjm_corr::upgrad, vpgrad, wpgrad: plain wide gradient). Next to the membrane it has to act
    // through the same operator as the pressure correction, i.e. with the face mobilities
    // (nhflow_membrane_gradx/grady) and beta on the vertical and sigma terms, not with the implicit
    // forcing factor (I + a A)^-1 in the layer or with mobility 1 next to it. Then P^n and PCORR act
    // like one full pressure P^{n+1} and the scheme is the projection of A 520 1.
    //   mode  1: take the plain gradient out before the forcing
    //   mode -1: put it back with the mobilities after the forcing
    double gx,gy,gz,sx,sy,a;
    const double *P = d->P;

    LOOP
    WETDRYDEEP
    {
        if(!nhflow_membrane_active(p,d,i,j,k))
        continue;

        a = alpha*p->dt*CPORNH;
        const double bv = nhflow_mbz(d,IJK);
        const double dPk = (P[FIJKp1]-P[FIJK])/p->DZN[KP];

        sx = 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*dPk;
        sy = 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*dPk;

        if(mode==1)
        {
            gx = (0.5*(P[FIp1JKp1]+P[FIp1JK])-0.5*(P[FIm1JKp1]+P[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1]) + sx;
            gy = p->j_dir==1 ? (0.5*(P[FIJp1Kp1]+P[FIJp1K])-0.5*(P[FIJm1Kp1]+P[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1]) + sy : 0.0;
            gz = dPk;
        }
        else
        {
            gx = -(nhflow_membrane_gradx(p,d,P,i,j,k) + bv*sx);
            gy = p->j_dir==1 ? -(nhflow_membrane_grady(p,d,P,i,j,k) + bv*sy) : 0.0;
            gz = -bv*dPk;
        }

        // wpgrad acts on WH without the water depth: dW = -a C/rho dP/(DZN WL)
        gx *= 1.0/p->W1;
        gy *= 1.0/p->W1;
        gz *= 1.0/(p->W1*MAX(WL(i,j),1.0e-20));

        d->U[IJK] += a*gx;
        UH[IJK]   += a*gx*WL(i,j);

        if(p->j_dir==1)
        {
        d->V[IJK] += a*gy;
        VH[IJK]   += a*gy*WL(i,j);
        }

        d->W[IJK] += a*gz;
        WH[IJK]   += a*gz*WL(i,j);
    }
}

void net_interface::membrane_reaction_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, slice &WL, bool finalize)
{
    for(auto m : pmem)
    m->reaction_nhflow(p,d,pgc,alpha,WL,finalize);
}

bool net_interface::membrane_iterated()
{
    for(auto m : pmem)
    if(m->iterated())
    return true;
    
    return false;
}

void net_interface::membrane_reforce_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    // strong coupling: forcing update for the node velocities of the next iteration
    for(auto m : pmem)
    m->reforce_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);
}

bool net_interface::membrane_couple_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, double alpha, slice &WL, int it)
{
    // strong coupling: loads of the projected velocity, structure, next node velocities; true when all converged
    bool conv=true;
    
    for(auto m : pmem)
    conv = m->couple_nhflow(p,d,pgc,iter,alpha,WL,it) && conv;
    
    return conv;
}

void net_interface::membrane_attach_nhflow(lexer *p, const Eigen::Vector3d &c, const Eigen::Matrix3d &R)
{
    for(auto m : pmem)
    m->attach_body(p,c,R);
}

void net_interface::membrane_body_nhflow(const Eigen::Vector3d &c, const Eigen::Matrix3d &R, const Eigen::Vector3d &v, const Eigen::Vector3d &w)
{
    for(auto m : pmem)
    m->set_body(c,R,v,w);
}

void net_interface::membraneForces_nhflow(lexer *p, const Eigen::Vector3d &c, const Eigen::Matrix3d &R, double &X, double &Y, double &Z, double &K, double &M, double &N)
{
    // load of all membranes on the floating body at its current position c and orientation R (flexible membranes:
    // linearised in the body motion since the last membrane step), moment about c; identical on all ranks
    X=Y=Z=K=M=N=0.0;
    
    for(auto m : pmem)
    {
        double x,y,z,k,mm,n;
        m->body_load(p,c,R,x,y,z,k,mm,n);
        
        X+=x; Y+=y; Z+=z;
        K+=k; M+=mm; N+=n;
    }
}

Eigen::Matrix3d net_interface::membrane_addedinertia_nhflow(lexer *p)
{
    Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
    
    for(auto m : pmem)
    I += m->body_addedinertia(p);
    
    return I;
}

double net_interface::membrane_addedmass_nhflow(lexer *p)
{
    double Ma=0.0;
    
    for(auto m : pmem)
    Ma += m->body_addedmass(p);
    
    return Ma;
}

Eigen::Matrix3d net_interface::membrane_stiffness_nhflow(lexer *p)
{
    // d F / d (body translation) of the flexible membranes attached to the body
    Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
    
    for(auto m : pmem)
    J += m->body_stiffness();
    
    return J;
}
