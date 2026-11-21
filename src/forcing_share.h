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

#ifndef FORCING_SHARE_H_
#define FORCING_SHARE_H_

#include"lexer.h"
#include"fdm.h"
#include"heaviside.h"
#include"increment.h"
#include<cmath>

// The share of a velocity face that the direct forcing takes to its target velocity (solids and bed
// of G 505 0, 6DOF bodies X 10), in one place for the forcing, the projection and the mesh refinement.
//
// The forcing sets u* = (1-w) u~ + w u_target at a face, w the share:
//   inside the surface (level set at the face < 0):   w = H
//   outside, tangential mode (B 21 1, X 14 > 0):      w = |n_d| H   (the normal component only)
// H the smoothed Heaviside of half width X 41 cells.
//
// X 46 and X 49: 1 the 6DOF bodies, 2 the bodies and the solids and the bed (the solids not tested yet).
// Not set in the control file: on (1) for bodies with the one-step scheme N 40 5 and the loads X 60 1
// (not with X 198 or the subcycled mesh refinement G 7 1),
// together with the loads from probe points off the hull (X 42 1.5, X 34 2) and the stabilised coupling
// X 197 1 (driver_logic_cfd); X 49 0 keeps the open projection.
//
// X 46 1: the share continuous across the surface.  The switch from |n_d| H to H at the level set 0
//         is a jump of up to 0.5 in the share of a face tangential to the surface; when a body moves,
//         each such face takes it once, a kick in the forcing and in the pressure of that step.
//         Here w = H (|n_d| + (1-|n_d|) s), s a cosine ramp from 0 at the surface to 1 at the depth
//         of the half width inside: outside unchanged, inside the tangential part comes in gradually.
//
// X 49 1: the projection does not correct the share w of a face (poisson_pcorr, pjm_corr): the
//         coefficient of the face and its velocity correction take the factor 1-w.  The forced
//         velocity is the blend above, so its divergence on the right-hand side is the divergence of
//         (1-w) u~ + w u_target already, and after the correction of the share 1-w the blend is free
//         of divergence.  A face with w = 1 keeps the target velocity exactly; the pressure cannot
//         move fluid through a body or a solid (the Darcy medium of permeability dt/rho of
//         solid_closed.h), in the air as little as in the water.  The share varies smoothly in time
//         as a body moves: no face is switched.

// the ramp of X 46 1 (phi the level set at the face, psi the half width)
inline double forcing_ramp(double phi, double psi)
{
    if(phi>=0.0)
    return 0.0;

    if(phi<=-psi)
    return 1.0;

    return 0.5*(1.0 - cos(PI*(-phi)/psi));
}

// the share of a face from its Heaviside H, |n_d| of the surface normal and the level set phi at
// the face (tangential mode; X 46 0: the caller keeps its own switch at phi = 0)
inline double forcing_weight_smooth(double H, double absn, double phi, double psi)
{
    return H*(absn + (1.0-absn)*forcing_ramp(phi,psi));
}

// |n_d| of the level set ls: central differences at cell (i0,j0,k0), or (face) the mean of the two
// cells of the face d (X 14 2, sixdof_obj_cfd::face_normal)
inline double forcing_absn(lexer *p, field &ls, int i0, int j0, int k0, int d, bool face)
{
    const int m = increment::marge;
    const int i1 = i0 + (d==0), j1 = j0 + (d==1), k1 = k0 + (d==2);

    double nx = (ls(i0+1,j0,k0) - ls(i0-1,j0,k0))/(2.0*p->DXN[i0+m]);
    double ny = (ls(i0,j0+1,k0) - ls(i0,j0-1,k0))/(2.0*p->DYN[j0+m]);
    double nz = (ls(i0,j0,k0+1) - ls(i0,j0,k0-1))/(2.0*p->DZN[k0+m]);

    if(face)
    {
    nx = 0.5*(nx + (ls(i1+1,j1,k1) - ls(i1-1,j1,k1))/(2.0*p->DXN[i1+m]));
    ny = 0.5*(ny + (ls(i1,j1+1,k1) - ls(i1,j1-1,k1))/(2.0*p->DYN[j1+m]));
    nz = 0.5*(nz + (ls(i1,j1,k1+1) - ls(i1,j1,k1-1))/(2.0*p->DZN[k1+m]));
    }

    const double norm = sqrt(nx*nx + ny*ny + nz*nz);

    if(norm<=1.0e-20)
    return 0.0;

    return fabs(d==0 ? nx : (d==1 ? ny : nz))/norm;
}

// the half width of the Heaviside at cell (i0,j0,k0) (ghostcell::Hsolidface, sixdof_obj_cfd::Hsolidface)
inline double forcing_psi(lexer *p, int i0, int j0, int k0)
{
    const int m = increment::marge;

    if(p->j_dir==0)
    return p->X41*0.5*(p->DXN[i0+m] + p->DZN[k0+m]);

    return p->X41*(1.0/3.0)*(p->DXN[i0+m] + p->DYN[j0+m] + p->DZN[k0+m]);
}

// the share of the face d (0 u, 1 v, 2 w) between cell (i0,j0,k0) and its neighbour in +d, as the
// forcing of the stage took it (ghostcell::solid_forcing, sixdof_obj_cfd::update_forcing).  With
// bodies (X 10) fbh1-3 hold the Heaviside of the bodies alone (sixdof_cfd::setup), the solids are
// taken from their level sets here; without bodies fbh1-3 are the solids
inline double forcing_share_face(lexer *p, fdm *a, int i0, int j0, int k0, int d, bool solids=true)
{
    const int i1 = i0 + (d==0), j1 = j0 + (d==1), k1 = k0 + (d==2);
    const double psi = forcing_psi(p,i0,j0,k0);

    double w = 0.0;

    // solids and bed (D 1 2: cut cells, walls of the projection without forcing, solid_aperture.h)
    if(solids && p->G505==0 && p->D1!=2 && (p->solidread>0 || p->topoforcing>0))
    {
        const double phs = 0.5*(a->solid(i0,j0,k0) + a->solid(i1,j1,k1));
        const double pht = 0.5*(a->topo(i0,j0,k0) + a->topo(i1,j1,k1));
        const bool both = p->topoforcing>0 && p->solidread>0;
        const double ph = both ? MIN(phs,pht) : (p->topoforcing>0 ? pht : phs);

        double Hs = heaviside(-ph,psi);

        if(p->B21==1 && Hs>0.0)
        {
            const bool topo_normal = p->topoforcing>0 && (p->solidread==0 || phs>=pht);
            field &ls = topo_normal ? (field&)a->topo : (field&)a->solid;

            if(p->X46==2)
            Hs = forcing_weight_smooth(Hs,forcing_absn(p,ls,i0,j0,k0,d,false),ph,psi);

            else if(ph>=0.0)
            Hs *= forcing_absn(p,ls,i0,j0,k0,d,false);
        }

        w += Hs;
    }

    // 6DOF bodies: the share the bodies took in their forcing of the stage (fbw1-3, summed over the
    // bodies).  Each body takes its share from its own level set; where two bodies meet (the leaves
    // of a mitre gate) the sum of their Heavisides is 1 but the sum of their shares is not, and the
    // level set of all bodies does not give the shares back.  Recomputed from it, a face in the
    // contact counted as closed while it was forced in part only: neither forced to the body
    // velocity nor corrected by the pressure (800 m/s within five steps at the gates of a lock)
    if(p->X10>0)
    {
        field &fbw = (d==0) ? (field&)a->fbw1 : (d==1 ? (field&)a->fbw2 : (field&)a->fbw3);
        w += fbw(i0,j0,k0);
    }

    return MAX(0.0,MIN(1.0,w));
}

// X 49: the open share 1-w of a face in the projection; below the limit the face is closed (a
// coefficient of that size only costs the condition of the matrix)
inline double forcing_open_face(lexer *p, fdm *a, int i0, int j0, int k0, int d)
{
    const double o = 1.0 - forcing_share_face(p,a,i0,j0,k0,d,p->X49==2);

    return o<1.0e-3 ? 0.0 : o;
}

// the face between cell (i0,j0,k0) and cell (i0+di,j0+dj,k0+dk)
inline double forcing_open(lexer *p, fdm *a, int i0, int j0, int k0, int di, int dj, int dk)
{
    if(di==1)  return forcing_open_face(p,a,i0,j0,k0,0);
    if(di==-1) return forcing_open_face(p,a,i0-1,j0,k0,0);
    if(dj==1)  return forcing_open_face(p,a,i0,j0,k0,1);
    if(dj==-1) return forcing_open_face(p,a,i0,j0-1,k0,1);
    if(dk==1)  return forcing_open_face(p,a,i0,j0,k0,2);
    return forcing_open_face(p,a,i0,j0,k0-1,2);
}

#endif
