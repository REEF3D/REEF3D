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

#ifndef SEASTATE_OBSTACLE_H_
#define SEASTATE_OBSTACLE_H_

#include<vector>

class lexer;
class fdm_seastate;
class seastate_grid;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - line obstacles (A 722 xs ys xe ye Kt Kr zc)

A cell face is blocked by an obstacle when the obstacle segment crosses
the line between the two cell centres. Through a blocked face

  transmitted   Kt^2 of the energy flux leaving the upwind cell
                (the inflow of the downwind cell is multiplied by Kt^2)
  reflected     Kr^2 of the flux leaving the cell through the face
                re-enters the cell in the direction mirrored at the
                obstacle line, theta' = 2 alpha - theta (specular)
  dissipated    the rest, 1 - Kt^2 - Kr^2

Kt and Kr are wave-height coefficients (as SWAN OBSTACLE TRANS/REFL).
Kt < 0: transmission over a dam or low-crested breakwater after Goda
(as SWAN DAM GODA, alpha 2.6, beta 0.15) from the freeboard
Rc = zc - water level and the larger Hs of the two cells,

  Kt = 0.5 (1 - sin(pi/(2 alpha) (Rc/Hs + beta)))  for -beta-alpha < Rc/Hs < alpha-beta
  Kt = 1 below, 0 above.

Faces are stored for every cell of the rank including the ghost cells:
east(i,j) is the face between (i,j) and (i+1,j), north(i,j) the face
between (i,j) and (i,j+1).
--------------------------------------------------------------------*/

class seastate_obstacle
{
public:
    struct face
    {
        float kt2 = 1.0f;       // energy transmission
        float kr2 = 0.0f;       // energy reflection
        float alpha = 0.0f;     // direction of the obstacle line [rad]
        int obs = -1;           // obstacle (A 722 index)
    };

    explicit seastate_obstacle(lexer *p);

    bool active() const {return nobs>0;}

    // the faces crossed by the obstacles
    void build(lexer *p);

    // transmission after Goda (Kt < 0) from the present Hs of the cells on both sides
    void update(lexer *p, fdm_seastate *e);

    const face *east(int i, int j) const;
    const face *north(int i, int j) const;

    int faces_blocked() const {return int(faces.size());}

private:
    int nobs;
    int imin, jmin, ni, nj;
    std::vector<double> xs, ys, xe, ye, kt, kr, zc;
    std::vector<int> fe, fn;                        // index into faces, -1 none
    std::vector<face> faces;
    std::vector<int> fi, fj, fdir;                  // cell and side (0 east, 1 north) of each face
};

#endif
