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
class ghostcell;
class fdm_seastate;
class seastate_grid;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - line obstacles, structures and reflecting coasts
(A 722, A 723 - A 726)

Line obstacles (A 722 xs ys xe ye Kt Kr zc): a cell face is blocked by
an obstacle when the obstacle segment crosses the line between the two
cell centres. Through a blocked face

  transmitted   Kt^2 of the energy flux leaving the upwind cell
                (the inflow of the downwind cell is multiplied by Kt^2)
  reflected     Kr^2 of the flux leaving the cell through the face
                re-enters the cell in the direction mirrored at the
                obstacle line, theta' = 2 alpha - theta (specular), or
                (A 726 n pown) spread over the directions around it with
                cos^pown, as SWAN RDIFF
  dissipated    the rest, 1 - Kt^2 - Kr^2

Kt and Kr are wave-height coefficients (as SWAN OBSTACLE TRANS/REFL).
Kt < 0: transmission over a dam or low-crested breakwater after Goda
(as SWAN DAM GODA, alpha 2.6, beta 0.15) from the freeboard
Rc = zc - water level and the larger Hs of the two cells,

  Kt = 0.5 (1 - sin(pi/(2 alpha) (Rc/Hs + beta)))  for -beta-alpha < Rc/Hs < alpha-beta
  Kt = 1 below, 0 above.

Structures (A 725 n type a b c, obstacle n of A 722 counted from 1):
  type 1  rubble-mound breakwater after d'Angremond, Van der Meer and
          De Jong (1996), as SWAN DAM DANGREMOND: a the seaward slope
          [deg], b the crest width B [m], crest level zc of A 722;
          Kt = -0.4 Rc/Hs + 0.64 (B/Hs)^-0.31 (1 - exp(-0.5 xi)),
          0.075 <= Kt <= 0.9, for B/Hs < 8; Kt = -0.35 Rc/Hs +
          0.51 (B/Hs)^-0.65 (1 - exp(-0.41 xi)), 0.05 <= Kt <=
          0.93 - 0.006 B/Hs, for B/Hs > 12; linear in between; xi the
          breaker parameter tan(slope)/sqrt(Hs/L0p), L0p = g Tp^2/(2 pi);
          Hs and the discrete peak period Tp of the cell with the larger
          Hs; Kr of A 722
  type 2  Kt and Kr per frequency from seastate-obstacle-n.dat, lines
          f [Hz] Kt Kr ('$' comments), linear in f, constant beyond
          the ends (as SWAN TRANS1D, with a reflection per frequency)
  type 3  porous breakwater through the water column per frequency
          (Madsen 1974; Sollitt and Cross 1972): a the width B [m],
          b the porosity n, c the stone diameter D50 [m]. Inside the
          structure the waves obey the linear wave equations with the
          inertia s = 1 + 0.34 (1-n)/n and the linearised Forchheimer
          resistance g (a + b |q|) q (van Gent 1995: a = 1000 (1-n)^2
          nu/(n^3 g D50^2), b = 1.1 (1-n)/(n^3 g D50), q the discharge
          velocity) as f sig, f = n g (a + sqrt(8/pi) b q_rms)/sig
          (equivalent linearisation for Gaussian velocities), so that
          sig^2 (s + i f) = g k_s tanh(k_s h) (the progressive mode,
          continued from the real root without resistance) and the flux impedance
          relative to the water outside is gamma = n k/k_s. A slab of
          width B between the two plane-wave regions gives per frequency
            T = 4 gamma e^(i phi)/D,  R = (1 - gamma^2)(1 - e^(2 i phi))/D,
            D = (1 + gamma)^2 - (1 - gamma)^2 e^(2 i phi),  phi = k_s B,
          Kt^2 = |T|^2, Kr^2 = |R|^2, and q_rms is the rms of the depth-
          averaged discharge velocity in the structure for the waves of the
          cell with the larger Hs that travel towards the structure (fixed-
          point iteration); normal incidence, the same coefficients for all
          directions. A single-mode (plane-wave) model, meant for coarse
          rubble: with a strong resistance (fine stones, f of a few) the
          progressive and the first evanescent mode can exchange their roles
          and the coefficients jump between frequencies (validation 37)

Reflecting coasts (A 723 Kr pown): every face between an active cell
and a land or dry cell inside the domain reflects Kr^2 of the flux that
leaves the active cell through it, specular (pown 0) or diffuse, at the
coastline through the cell: the principal direction (total least squares)
of the midpoints of the faces between active and land cells within A 724
cells (default 3, at most 3: the ghost cells), averaged (as 2 alpha) over
the coastal cells within A 724 cells, so that a straight coast on the
staircase of cells reflects at its own direction (within about 1 deg);
A 724 0: the face itself. Faces with an obstacle keep the obstacle. The
coast faces are rebuilt when the active cells change (wetting and drying).

Faces are stored for every cell of the grid including the ghost cells:
east(i,j) is the face between (i,j) and (i+1,j), north(i,j) the face
between (i,j) and (i,j+1). The same object works on level 0 and on the
REEFAMR patches (build with the lexer and field of the grid).
--------------------------------------------------------------------*/

class seastate_obstacle
{
public:
    struct face
    {
        float kt2 = 1.0f;       // energy transmission
        float kr2 = 0.0f;       // energy reflection
        float alpha = 0.0f;     // direction of the obstacle line or the coastline [rad]
        int obs = -1;           // obstacle (A 722 index), -1: coast (A 723)
        int fq = -1;            // per-frequency coefficients (index into the pools), -1: kt2, kr2
        int wd = -1;            // diffuse reflection (index into the weights), -1: specular
    };

    seastate_obstacle(lexer *p, const seastate_grid *g);

    bool active() const {return nobs>0 || coast;}

    // the faces of the grid of lexer q (level 0 or a REEFAMR patch) and field e: obstacles, and (A 723)
    // coasts between active cells and land cells inside the domain; global index = local + (oi, oj),
    // the domain gnx x gny cells of the level; pgc: the coastline directions of the ghost cells from the
    // neighbouring ranks (level 0), nullptr: from the cells of the grid (patches)
    void build(lexer *q, fdm_seastate *e, int oi, int oj, int gnx, int gny, ghostcell *pgc=nullptr);

    // Goda, d'Angremond and porous coefficients from the present spectra; coasts rebuilt after a change
    // of the active cells
    void update(lexer *q, fdm_seastate *e);

    const face *east(int i, int j) const;
    const face *north(int i, int j) const;

    // energy transmission per frequency (nullptr: kt2 for all) and reflection of frequency l
    const float *kt2f(const face *f) const {return f->fq<0 ? nullptr : &fkt[size_t(f->fq)*nsig];}
    double kr2(const face *f, int l) const {return f->fq<0 ? double(f->kr2) : double(fkr[size_t(f->fq)*nsig+l]);}

    // diffuse reflection: weights of the directions -nd..nd around the specular one, nullptr: specular
    const float *diffuse(const face *f, int &nd) const
    {
        if(f->wd<0) {nd = 0; return nullptr;}
        nd = int(wdif[f->wd].size()/2);
        return wdif[f->wd].data();
    }

    int faces_blocked() const {return int(faces.size());}
    int faces_coast() const {return ncoast;}

private:
    const seastate_grid *g;
    int nobs, nsig;
    bool coast;
    double kr2c;
    int radius;
    int imin=0, jmin=0, ni=0, nj=0, oi=0, oj=0, gnx=0, gny=0;
    std::vector<double> xs, ys, xe, ye, kt, kr, zc;
    std::vector<int> type;                          // per obstacle: 0 Kt (Goda for Kt < 0), 1 d'Angremond, 2 table, 3 porous
    std::vector<double> pa, pb, pc;                 // structure parameters
    std::vector<std::vector<float>> tkt, tkr;       // type 2: Kt^2, Kr^2 per frequency
    std::vector<int> wobs;                          // per obstacle: diffuse weights, -1 specular
    int wcoast = -1;
    std::vector<std::vector<float>> wdif;           // diffuse weights, 2 nd + 1 each
    std::vector<int> fe, fn;                        // index into faces, -1 none
    std::vector<face> faces;
    std::vector<int> fi, fj, fdir;                  // cell and side (0 east, 1 north) of each face
    std::vector<float> fkt, fkr;                    // per-frequency coefficients
    std::vector<char> wmask;                        // active cells at the last build (coasts), ghost layers exchanged
    std::vector<char> wraw;                         // ... as the field had them (wetting and drying)
    int ncoast = 0;
    ghostcell *pgc = nullptr;
    int weights(double pown);
};

#endif
