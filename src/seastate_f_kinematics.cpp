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

#include"seastate_f.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_dispersion.h"
#include"lexer.h"
#include"ghostcell.h"

/*--------------------------------------------------------------------
Kinematics of the spectral grid: wave number k and group velocity cg
per active cell and frequency, depth and current gradients (central
differences, one-sided at the edge of the domain), rate of change of
the depth. Called after every change of
the environment (initially, and later by a host model every step).

The propagation velocities themselves are computed in the solver from
these quantities (plan Section 6: no per-bin velocity arrays):
  cx = cg cos(theta) + U,  cy = cg sin(theta) + V
  c_theta = sig/sinh(2kd) (sin(theta) dd/dx - cos(theta) dd/dy)
          + sin(theta)cos(theta)(dU/dx - dV/dy) + sin^2(theta) dV/dx - cos^2(theta) dU/dy
  c_sigma = k sig/sinh(2kd) (dd/dt + U dd/dx + V dd/dy)
          - cg k (cos^2(theta) dU/dx + sin(theta)cos(theta)(dU/dy + dV/dx) + sin^2(theta) dV/dy)
--------------------------------------------------------------------*/

// neighbour status: 1 active, 0 outside the domain (edge), -1 land or dry
static int status(lexer *p, int i, int j, sliceint &wet)
{
    if(wet(i,j)==1)
    return 1;

    const int gi = i + p->origin_i;
    const int gj = j + p->origin_j;

    if(gi<0 || gi>=p->gknox || gj<0 || gj>=p->gknoy)
    return 0;

    return -1;
}

void seastate_f::kinematics(lexer *p, ghostcell *pgc, double dt_depth)
{
    auto status = [&](lexer *pp, int ii, int jj) {return ::status(pp,ii,jj,e->wet);};

    const seastate_grid &g = *e->grid;

    IMALOOP
    JMALOOP
    {
        // k and cg
        float *kk = e->kw->spec(i,j);
        float *cc = e->cg->spec(i,j);

        if(kk!=nullptr)
        {
            for(int l=0; l<g.nsig; ++l)
            {
                if(e->wet(i,j)==1)
                {
                const double k = seastate_wavenumber(g.sig[l],e->depth(i,j));
                kk[l] = float(k);
                cc[l] = float(seastate_cg(g.sig[l],k,e->depth(i,j)));
                }
                else
                kk[l] = cc[l] = 0.0f;
            }
        }

    e->ddx(i,j)=e->ddy(i,j)=e->dUdx(i,j)=e->dUdy(i,j)=e->dVdx(i,j)=e->dVdy(i,j)=e->dddt(i,j)=0.0;
    e->refr(i,j)=0;

        // gradients: central differences; one-sided next to the edge of the domain; refraction and
        // frequency shift off next to land or dry cells (as SWAN)
        if(e->wet(i,j)==1 && i>-p->margin && i<p->knox+p->margin-1 && j>-p->margin && j<p->knoy+p->margin-1)
        {
        const int sw = status(p,i-1,j), se = status(p,i+1,j), ss = status(p,i,j-1), sn = status(p,i,j+1);

            if(sw>=0 && se>=0 && ss>=0 && sn>=0)
            {
            const int iw = sw==1 ? i-1 : i, ie = se==1 ? i+1 : i;
            const int js = ss==1 ? j-1 : j, jn = sn==1 ? j+1 : j;

                if(ie>iw)
                {
                const double rdx = 1.0/(p->XP[ie+marge]-p->XP[iw+marge]);
                e->ddx(i,j)  = (e->depth(ie,j)-e->depth(iw,j))*rdx;
                e->dUdx(i,j) = (e->U(ie,j)-e->U(iw,j))*rdx;
                e->dVdx(i,j) = (e->V(ie,j)-e->V(iw,j))*rdx;
                }

                if(jn>js)
                {
                const double rdy = 1.0/(p->YP[jn+marge]-p->YP[js+marge]);
                e->ddy(i,j)  = (e->depth(i,jn)-e->depth(i,js))*rdy;
                e->dUdy(i,j) = (e->U(i,jn)-e->U(i,js))*rdy;
                e->dVdy(i,j) = (e->V(i,jn)-e->V(i,js))*rdy;
                }

                if(dt_depth>0.0)
                e->dddt(i,j) = (e->depth(i,j)-e->depth_n(i,j))/dt_depth;

            e->refr(i,j) = 1;
            }
        }

    e->depth_n(i,j) = e->depth(i,j);
    }
}
