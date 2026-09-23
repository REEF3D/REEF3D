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

#ifndef FNPF_FSFBC_WD_H_
#define FNPF_FSFBC_WD_H_

#include "fnpf_breaking.h"
#include "sliceint4.h"
#include "slice4.h"
#include "fnpf_voiddisc.h"
#include "fnpf_cds2_wd.h"
#include "fnpf_cds4_wd.h"
#include "fnpf_cds6_wd.h"
#include "fnpf_weno3.h"
#include "fnpf_weno5.h"
#include "fnpf_weno5_wd.h"
#include "fnpf_cds4.h"
#include "fnpf_hires.h"
#include "fnpf_ddx_cds2.h"
#include "fnpf_ddx_cds4.h"
#include <optional>
#include <variant>

class fnpf_laplace;
class field;
class fnpf_coastline;
class solver2D;
class wind;

using namespace std;

class fnpf_fsfbc_wd final : public fnpf_breaking
{
public:
    fnpf_fsfbc_wd(lexer*, fdm_fnpf*, ghostcell*);
    virtual ~fnpf_fsfbc_wd();

    void fsfdisc(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&) override final;
    void fsfdisc_ini(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&) override final;
    void kfsfbc(lexer*,fdm_fnpf*,ghostcell*) override final;
    void dfsfbc(lexer*,fdm_fnpf*,ghostcell*,slice&) override final;
    void fsfwvel(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&) override final;
    void wetdry(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&) override final;
    void coastline_eta(lexer*,fdm_fnpf*,ghostcell*,slice&) override final;
    void coastline_fi(lexer*,fdm_fnpf*,ghostcell*,slice&) override final;
    void coastline_fi_ini(lexer*,fdm_fnpf*,ghostcell*,slice&) override final;
    void coastline_vel(lexer*,fdm_fnpf*,ghostcell*,double*) override final;
    void damping(lexer*,fdm_fnpf*,ghostcell*,slice&,int,double) override final;

    void coastline_Fz(lexer*,fdm_fnpf*,ghostcell*,slice&);

private:
    double rb3(lexer*,double);
    double rb4(lexer*,double);
    double rb5(lexer*,double);

    sliceint4 wetcoast;
    slice4 ef,df;

    std::variant<fnpf_voiddisc, fnpf_cds2_wd, fnpf_cds4_wd, fnpf_weno3, fnpf_weno5_wd, fnpf_cds6_wd> pconvec;
    std::optional<fnpf_weno5> pconeta; // eta discretisation next to fnpf_weno5_wd, otherwise pconvec is used
    std::variant<fnpf_hires, fnpf_cds4> pdx;
    std::variant<fnpf_ddx_cds2, fnpf_ddx_cds4> pddx;
    fnpf_coastline *pcoast;
    solver2D *psolv;
    wind *pwind;

    double ivel,jvel,kvel;

    double dist3,dist4,dist5,expinverse,db;

    double visc;

    int *temp;
    int gcval_eta,gcval_fifsf;

    int count_n;
    int coastline_count;

    static constexpr double eps = 1.0e-6;
};

#endif
