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

#ifndef PICARD_LSM_H_
#define PICARD_LSM_H_

#include"gradient.h"
#include"picard.h"

using namespace std;

class picard_lsm final : public gradient, public picard
{
public:

    picard_lsm(lexer *p);
    virtual ~picard_lsm();

    void volcalc(lexer*, fdm*, ghostcell*, field&) override final;
    void volcalc2(lexer*, fdm*, ghostcell*, field&);
    void correct_ls(lexer*, fdm*, ghostcell*, field&) override final;

private:
    double vol1,vol2,dvol2;
    // the reference and the net volume, the reference taken (1, first volcalc): one array for the
    // checkpoints (ckpt_registry)
    double vs[3] = {0.0,0.0,0.0};
    double &inivol = vs[0], &netvol = vs[1], &vol_ini = vs[2];
    unsigned long long ckid;
};

#endif

