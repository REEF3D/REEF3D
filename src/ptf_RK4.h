/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"ptf.h"
#include"ptf_fsfbc.h"
#include"slice4.h"

class ptf_laplace;
class field;
class ptf_fsf_update;
class ptf_bed_update;

using namespace std;

#ifndef PTF_RK4_H_
#define PTF_RK4_H_

class ptf_RK4 : public ptf, public ptf_fsfbc
{
public:
	ptf_RK4(lexer*, fdm*, ghostcell*);
	virtual ~ptf_RK4();
    
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

    ptf_laplace *plap;
    ptf_fsf_update *pfsfupdate;
    ptf_bed_update *pbedupdate;


};

#endif
