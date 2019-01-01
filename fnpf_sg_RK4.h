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

#include"fnpf_sg.h"
#include"fnpf_sg_ini.h"
#include"fnpf_sigma.h"
#include"slice4.h"

class fnpf_sg_laplace;
class fnpf_sg_fsfbc;
class field;

using namespace std;

#ifndef FNPF_SG_RK4_H_
#define FNPF_SG_RK4_H_

class fnpf_sg_RK4 : public fnpf_sg_ini, public fnpf_sigma
{
public:
	fnpf_sg_RK4(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_sg_RK4();
    
    virtual void start(lexer*, fdm_fnpf*, ghostcell*, solver*, convection*, ioflow*, reini*,onephase*);
    virtual void inidisc(lexer*, fdm_fnpf*, ghostcell*);
    
private:

    int gcval,gcval_u,gcval_v,gcval_w;
    int gcval_eta,gcval_fifsf;
    int hypre_type;
    double starttime,endtime;

    slice4 erk1,erk2,erk3,erk,en;
    slice4 frk1,frk2,frk3,frk;

    fnpf_sg_laplace *plap;
    fnpf_sg_fsfbc *pf;

};

#endif
