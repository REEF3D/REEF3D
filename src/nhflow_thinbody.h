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


#ifndef NHFLOW_THINBODY_H_
#define NHFLOW_THINBODY_H_

#include"nhflow_convection.h"
#include"increment.h"
#include"slice4.h"
#include"lexer.h"

class lexer;
class fdm_nhf;
class ghostcell;
class slice;

using namespace std;

// Sharp thin bodies in NHFLOW: a service for bodies thinner than a cell (membranes X 330 'mobility sharp'; plates,
// porous sheets later). The body is not smeared over a layer of cells: only the links between two neighbouring cell
// centres (vertically: two nodes) that lie on opposite sides of the body are blocked, the cells on either side are
// ordinary fluid cells. The clients (net_membrane) mark in every stage
//   sideC       the side of the cell centres near the body (+1 / -1, 0 farther than one link length),
//   bx, by, bz  the blocked links: x-link (i,j,k)-(i+1,j,k), y-link (i,j,k)-(i,j+1,k), vertical link node k - node k+1
//               (through cell k, the "cut cell"); ux, uy, uz the wall velocity component along the link,
// together with the link mobilities d->MBX, MBY, MBZ of the projection (nhflow_membrane_beta.h). The service acts on
// NHFLOW through three hooks, without changes to the sigma grid:
//
// 1. Fluxes (flux_hook, the flux hook of the convection scheme, called before the divergence)
//    blocked faces    wall flux, separately for the two sides: hydrostatic flux g(eta^2/2 + eta h) of each side with
//                     its own head, advective flux with the wall velocity; continuity flux D u_wall (the face holds
//                     the left flux, the right cell gets the difference in its right-hand side)
//    near the body    the reconstruction (WENO5 by default; eta per column) would reach across the body. Faces within
//                     two cells of a column with a blocked link get a local flux: states from the cells on their own
//                     side (van Leer, slope 0 at a blocked link), central hydrostatic flux of the cell heads, HLL
//                     dissipation with the gravity wave speed for the continuity and the normal momentum, with the
//                     flow speed for the tangential momentum (as HLLC)
//    vertical faces   no momentum flux through a vertical face between cell centres on opposite sides; within two
//                     faces of it first-order upwind
//    below floors     closed bags: a column has one free surface (inside the bag), the water below the floor is
//                     outside water. Cells with an odd number of blocked vertical faces above them ("lower cells")
//                     use the hydrostatic head eta_L instead of eta: the outer free surface carried across the
//                     footprint (harmonic extension, Dirichlet eta at the columns around it). Faces next to lower cells:
//                     local flux with eta_L and flow-speed dissipation (no free surface below a floor), continuity the
//                     Rhie-Chow flux of the projection plus the rigid-lid correction (lid_correction: gradient of a
//                     lid potential, so that the lower cells of a footprint column have no net outflow; the inner
//                     free surface only sees the inner fluxes)
//    The non-hydrostatic pressure P then carries only the dynamic part; the static jump rho g (eta_L - eta) across
//    the floor is carried by the body.
// 2. Projection (projection_rhs in nhflow_pjm::rhs, nhflow_membrane_row, nhflow_membrane_gradx/grady): the face
//    velocity of a blocked link in the divergence is the wall velocity, not the average of the two cell velocities;
//    with the (near) zero link mobility the projected face flux through the body is the wall flux. Next to a floor
//    the vertical control volume of a pressure node spans two cell layers: the half face of a layer whose cells both
//    lie on the other side of the body than the node is closed (divergence and matrix), a horizontal node link
//    between nodes on opposite sides gets the blocked mobility, and the cell pressure of the velocity correction is
//    the node on the cell's own side (pcell). Without this the flow below the bag drives the water inside it.
//    Cells that crossed the body with the moving sigma grid get the velocity of their new side (side_change).
// 3. Cut cells (cut_forcing, before the projection): the vertical velocity of a cell with a blocked vertical node link
//    is the wall velocity (its W is the link velocity of the vertical Poisson control volumes). Pockets of the staircase
//    (cells with three or more blocked faces, e.g. below the rim of a cone floor) move with the body: in such a
//    corner cell the collocated pressure correction acts through one or two faces only and the cell velocity can grow
//    unchecked.
//
// Diffusion and turbulence: the implicit diffusion steps (momentum, k, epsilon/omega) get zero gradient across the
// blocked faces (matrix_walls), the faces at the body are walls with a wall function for the turbulence model and,
// with a turbulence model, wall friction in the momentum equations (nhflow_wall.h, nhflow_bcmom; ks = B57).
//
// Loads: total pressure jump across a blocked link (hydrostatic with eta / eta_L, plus P) times the face area
// (link_force). Shear on the body and the momentum of the cut-cell forcing are not included.
// Not covered: VRANS porosity at blocked faces, the incremental pressure scheme A 520 2, flexible membranes.

class nhflow_thinbody : public nhflow_flux_hook, public increment
{
public:
    nhflow_thinbody(lexer*, fdm_nhf*, ghostcell*);
    virtual ~nhflow_thinbody();
    
    // clients, every stage: begin() before marking, finish() after
    void begin(lexer*);
    void finish(lexer*, fdm_nhf*, ghostcell*);
    
    double *sideC;              // side of the cell centre, written by the clients
    double *sideN;              // side of the nodes (FIJK), written by the clients
    double *bx,*by,*bz;         // 1: blocked link
    double *ux,*uy,*uz;         // wall velocity along the blocked link
    
    // 1. fluxes
    void flux_hook(lexer*, fdm_nhf*, int, int, double*, double*) override;
    
    // 2. projection right-hand side (nhflow_pjm::rhs)
    void projection_rhs(lexer*, fdm_nhf*, double*, double*, double*, double);
    
    // 3. cut cells
    void cut_forcing(lexer*, fdm_nhf*, double*, double*, double*, slice&);
    
    // cells that changed side since the last stage (the sigma grid moves with the free surface, a cell centre next
    // to a floor can cross it): velocity of the neighbours on the new side, so that no momentum is carried through
    void side_change(lexer*, fdm_nhf*, double*, double*, double*, slice&);
    
    // walls for diffusion and turbulence (nhflow_wall.h): faces of cell (i,j,k) at the body, 0 x-, 1 x+, 2 y-,
    // 3 y+, 4 below, 5 above (blocked links; vertically the faces between cell centres on opposite sides)
    void wall_faces(lexer*, int, int, int, int*) const;
    
    // zero gradient across the blocked faces in d->M / d->rhsvec of an implicit step of F (momentum diffusion,
    // k, epsilon/omega): no diffusive or eddy-viscous exchange through the body
    void matrix_walls(lexer*, fdm_nhf*, const double*);
    
    // pressure of a cell for the horizontal velocity correction: the mean of its two nodes, next to a floor only the
    // node on the cell's own side (a cut layer would otherwise mix the pressure of both sides of the body)
    inline double pcell(lexer *p, const double *P, int ii, int jj, int kk) const
    {
        const int q = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;
        const int f = (ii-p->imin)*p->jmax*p->kmaxF + (jj-p->jmin)*p->kmaxF + kk-p->kmin;
        return pw[q]*P[f] + (1.0-pw[q])*P[f+1];
    }
    
    // nodes (i,j,k) and its +x / +y neighbour on opposite sides of the body (horizontal node link of the Poisson
    // equation crossing it)
    inline bool node_cut_x(lexer *p, int ii, int jj, int kk) const
    {
        const int f = (ii-p->imin)*p->jmax*p->kmaxF + (jj-p->jmin)*p->kmaxF + kk-p->kmin;
        return sideN[f]*sideN[f+p->jmax*p->kmaxF]<0.0;
    }
    
    inline bool node_cut_y(lexer *p, int ii, int jj, int kk) const
    {
        const int f = (ii-p->imin)*p->jmax*p->kmaxF + (jj-p->jmin)*p->kmaxF + kk-p->kmin;
        return sideN[f]*sideN[f+p->kmaxF]<0.0;
    }
    
    // node control volumes straddling a floor: the half face of layer l of the node row belongs to the other side when
    // both cells of that layer lie on the other side of the body than the node (side sN)
    inline double nside(lexer *p, int ii, int jj, int kk) const
    {
        return sideN[(ii-p->imin)*p->jmax*p->kmaxF + (jj-p->jmin)*p->kmaxF + kk-p->kmin];
    }
    
    inline bool other_side(int qa, int qb, double sN) const
    {
        return sN!=0.0 && sideC[qa]*sN<0.0 && sideC[qb]*sN<0.0;
    }
    
    // loads: force on the body [N] of the blocked link dir (0 x, 1 y, 2 vertical) of cell (i,j,k), in +dir
    double link_force(lexer*, fdm_nhf*, int, int, int, int);
    
    // hydrostatic head of a cell: eta, below a closed floor eta_L
    double head(lexer*, fdm_nhf*, int, int, int);
    
    int nlower(lexer*, ghostcell*);
    
private:
    void update_etaL(lexer*, fdm_nhf*, int);
    void lid_correction(lexer*, fdm_nhf*, double*, double*);
    
    ghostcell *pgc;
    int ncell;
    
    double *low;                // 1: lower cell (below a closed floor)
    double *pw;                 // weight of node k in the cell pressure (pcell), 0.5 away from the body
    double *side0;              // sideC of the previous stage
    double *fz;                 // 1: vertical momentum face between cells k and k+1 blocked (cell centres)
    slice4 etaL;                // hydrostatic head below closed floors
    slice4 fp;                  // 1: column with lower cells (footprint of a closed floor)
    slice4 cbx,cby;             // 1: column with a blocked link to its +x / +y neighbour
    slice4 phi,rL,Hx,Hy;        // rigid lid below closed floors: potential, lower outflow, open lower layer depth
    bool first, first_lid=true;
    int nlow, nlowg;
    
    int i,j,k;
};

#endif
