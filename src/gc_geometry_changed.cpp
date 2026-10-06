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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

// The flag rebuild (solid_forcing_flag_update, gcdf_update, gcb_velflagio) depends on the geometry
// only through these comparisons: solid, topo against 0 and against psi = -X41*dx (G5), fb against 0.
// geometry_changed() stores them per cell and returns true when any of them changed on any rank
// since the last call (always true on the first call), so the rebuild can be skipped otherwise.
bool ghostcell::geometry_changed(lexer *p, fdm *a)
{
    const int size = p->imax*p->jmax*p->kmax;
    int changed = 0;

    if((int)geo_sig.size()!=size)
    {
    geo_sig.assign(size,0xffff);
    changed=1;
    }

    double psi;

    BASELOOP
    {
    if(p->j_dir==0)
    psi = -p->X41*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]);

    if(p->j_dir==1)
    psi = -p->X41*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    const double s = a->solid(i,j,k);
    const double t = a->topo(i,j,k);
    const double f = a->fb(i,j,k);

    const uint16_t sig = (s<0.0) | (s>0.0)<<1 | (s<psi)<<2 | (s>=psi)<<3
                       | (t<0.0)<<4 | (t>0.0)<<5 | (t<psi)<<6 | (t>=psi)<<7
                       | (f<0.0)<<8 | (f>0.0)<<9;

    if(geo_sig[IJK]!=sig)
    {
    geo_sig[IJK]=sig;
    changed=1;
    }
    }

    return globalimax(changed)>0;
}
