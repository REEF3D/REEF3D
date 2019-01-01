/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"fnpf_fg.h"
#include"fnpf_fg_fsfbc.h"
#include"slice4.h"

class fnpf_fg_laplace;
class fnpf_fg_fsf_update;
class fnpf_fg_bed_update;
class field;

using namespace std;

#ifndef FNPF_FG_RK3_H_
#define FNPF_FG_RK3_H_

class fnpf_fg_RK3 : public fnpf_fg, public fnpf_fg_fsfbc
{
public:
	fnpf_fg_RK3(lexer*, fdm*, ghostcell*);
	virtual ~fnpf_fg_RK3();
    
    virtual void start(lexer*, fdm*, ghostcell*, solver*, convection*, ioflow*, reini*,onephase*);
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*, reini*, onephase*);
    virtual void inidisc(lexer*, fdm*, ghostcell*);
    
private:

    int gcval_eta,gcval_fifsf,gcval;
    int hypre_type;
    double starttime,endtime;

    slice4 erk1,erk2;
    slice4 frk1,frk2;

    fnpf_fg_laplace *plap;
    fnpf_fg_fsf_update *pfsfupdate;
    fnpf_fg_bed_update *pbedupdate;

};

#endif
