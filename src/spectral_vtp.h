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

#ifndef SPECTRAL_VTP_H_
#define SPECTRAL_VTP_H_

#include"increment.h"
#include"vtp3D.h"

class lexer;
class fdm_spectral;
class ghostcell;
class slice;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::Spectral - VTP output of the integrated wave parameters on the
2D grid, folder REEF3D_SPECTRAL_VTP

  print interval: P 181 (iterations) or P 182 (seconds), as SFLOW
  point data: Hs, Tm01, Tm-10, Tp, dir, spread, depth, wet
  nodal values average the active (wet) neighbour cells only, so that
  land cells do not reduce the values along the coast
--------------------------------------------------------------------*/

class spectral_vtp : public increment, private vtp3D
{
public:
    spectral_vtp(lexer*,fdm_spectral*,ghostcell*);
    virtual ~spectral_vtp() = default;

    void start(lexer*,fdm_spectral*,ghostcell*);
    void print2D(lexer*,fdm_spectral*,ghostcell*);

private:
    void pvtp(lexer*,int);
    float node_wet(lexer*,fdm_spectral*,slice&);
    void write_scalar(lexer*,fdm_spectral*,slice&,ofstream&);

    char name[200];
    int n,iin,offset[200];
    float ffn;

    static const int nscalar = 8;
    static const char *scalar_name[nscalar];
};

#endif
