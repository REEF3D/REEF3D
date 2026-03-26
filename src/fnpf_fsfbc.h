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

#ifndef FNPF_FSFBC_H_
#define FNPF_FSFBC_H_

#include"fnpf_breaking.h"
#include"sliceint4.h"

class fnpf_laplace;
class field;
class fnpf_convection;
class fnpf_ddx;
class fnpf_etadisc;
class solver2D;
class wind;

using namespace std;

class fnpf_fsfbc final : public fnpf_breaking
{
public:
	fnpf_fsfbc(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_fsfbc();
    
    void fsfdisc(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&) override final;
    void fsfdisc_ini(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&) override final;
    void kfsfbc(lexer*,fdm_fnpf*,ghostcell*) override final;
    void dfsfbc(lexer*,fdm_fnpf*,ghostcell*,slice&) override final;
    void fsfwvel(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&) override final;
    void wetdry(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&) override final;
    void coastline_eta(lexer*,fdm_fnpf*,ghostcell*,slice&) override final;
    void coastline_fi(lexer*,fdm_fnpf*,ghostcell*,slice&) override final;
    void damping(lexer*,fdm_fnpf*,ghostcell*,slice&,int,double) override final;
    
    
    void filter(lexer*, fdm_fnpf*,ghostcell*, slice&);

    fnpf_convection *pconvec;
    fnpf_convection *pconeta;
    fnpf_etadisc *pdf;
    fnpf_convection *pdx;
    fnpf_ddx *pddx;
    solver2D *psolv;
    wind *pwind;
    
    
private:
    double ivel,jvel,kvel;
    
    double visc;
    const double eps;


    int *temp;
    int gcval_eta,gcval_fifsf;
    
};

#endif
