/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"fnpf_fg.h"
#include"fnpf_fg_fsfbc.h"
#include"fnpf_fg_fsf_update.h"
#include"fnpf_fg_bed_update.h"
#include"slice4.h"

class fnpf_fg_laplace;
class field;

using namespace std;

#ifndef FNPF_FG_RK4_H_
#define FNPF_FG_RK4_H_

class fnpf_fg_RK4 : public fnpf_fg, public fnpf_fg_fsfbc, public fnpf_fg_fsf_update, public fnpf_fg_bed_update
{
public:
	fnpf_fg_RK4(lexer*, fdm*, ghostcell*);
	virtual ~fnpf_fg_RK4();
    
    virtual void start(lexer*, fdm*, ghostcell*, solver*, convection*, ioflow*, reini*,onephase*);
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*, reini*, onephase*);
    virtual void inidisc(lexer*, fdm*, ghostcell*);
    
    
private:

    int gcval,gcval_u,gcval_v,gcval_w;
    int gcval_eta,gcval_fifsf;
    int hypre_type;
    double starttime,endtime;

    slice4 erk1,erk2,erk3,erk;
    slice4 frk1,frk2,frk3,frk;

    fnpf_fg_laplace *plap;


};

#endif
