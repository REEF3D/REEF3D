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

#ifndef FNPF_RK2_H_
#define FNPF_RK2_H_

#include"fnpf.h"
#include"fnpf_ini.h"
#include"fnpf_sigma.h"
#include"slice4.h"

class fnpf_laplace;
class fnpf_fsf;
class field;

using namespace std;

class fnpf_RK2 : public fnpf_ini, public fnpf_sigma
{
public:
	fnpf_RK2(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_RK2();
    
    virtual void start(lexer*, fdm_fnpf*, ghostcell*, solver*, convection*, ioflow*, reini*,onephase*);
    virtual void inidisc(lexer*, fdm_fnpf*, ghostcell*, ioflow*, solver*);
    virtual void ini_wetdry(lexer*, fdm_fnpf*, ghostcell*);
    
private:

    int gcval,gcval_u,gcval_v,gcval_w;
    int gcval_eta,gcval_fifsf;
    int hypre_type;
    double starttime,endtime;

    slice4 erk1;
    slice4 frk1;

    fnpf_laplace *plap;
    fnpf_fsf *pf;
    
    int gcval_sl;
    double t0;

};

#endif
