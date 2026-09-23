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

#include "fnpf_breaking.h"
#include "slice4.h"
#include "fnpf_voiddisc.h"
#include "fnpf_cds2.h"
#include "fnpf_cds4.h"
#include "fnpf_cds6.h"
#include "fnpf_weno3.h"
#include "fnpf_weno5.h"
#include "fnpf_ddx_cds2.h"
#include "fnpf_ddx_cds4.h"
#include <optional>
#include <variant>

class fnpf_laplace;
class field;
class solver2D;
class wind;

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
    void coastline_eta(lexer*,fdm_fnpf*,ghostcell*,slice&) override final {};
    void coastline_fi(lexer*,fdm_fnpf*,ghostcell*,slice&) override final {};
    void coastline_fi_ini(lexer*,fdm_fnpf*,ghostcell*,slice&) override final {};
    void coastline_vel(lexer*,fdm_fnpf*,ghostcell*,double*) override final {};
    void damping(lexer*,fdm_fnpf*,ghostcell*,slice&,int,double) override final;

private:
    slice4 ef,df;

    std::variant<fnpf_voiddisc, fnpf_cds2, fnpf_cds4, fnpf_weno3, fnpf_weno5, fnpf_cds6> pconvec;
    std::variant<fnpf_ddx_cds2, fnpf_ddx_cds4> pddx;
    solver2D *psolv;
    wind *pwind;

    double ivel,jvel,kvel;

    double visc;

    int gcval_eta,gcval_fifsf;
};

#endif
