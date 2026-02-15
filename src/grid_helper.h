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

#ifndef GRID_HELPER_H_
#define GRID_HELPER_H_

#include"increment.h"

class lexer;

class grid_helper : public increment
{
public:
    grid_helper(lexer*);
    virtual ~grid_helper();

    // gcb
    void fillgcb1(lexer*);
    void fillgcb2(lexer*);
    void fillgcb3(lexer*);
    void fillgcb4a(lexer*);

    void fillgcb4_wall(lexer*);

    // dgc
    void make_dgc(lexer*);
    void fill_dgc1(lexer*);
    void fill_dgc2(lexer*);
    void fill_dgc3(lexer*);
    void fill_dgc4(lexer*);

private:
    int imin,imax,jmax,jmin,kmin,kmax;

    int *hgc;
    int **fgc;
};

#endif
