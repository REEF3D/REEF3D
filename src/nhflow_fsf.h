/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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


class convection;
class pressure;
class solver;
class fdm_nhf;
class lexer;
class field;
class ghostcell;
class fluid_update;
class heat;
class concentration;
class ioflow;
class slice;
class momentum;
class diffusion;
class poisson;
class turbulence;

using namespace std;

#ifndef NHFLOW_FSF_H_
#define NHFLOW_FSF_H_

class nhflow_fsf
{
public:    
    virtual void start(lexer*, fdm_nhf*, ghostcell*, ioflow*)=0;
    virtual void ini(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*)=0;
    
    virtual void rk2_step1(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double)=0;
    virtual void rk2_step2(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double)=0;
    
    virtual void rk3_step1(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double)=0;
    virtual void rk3_step2(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double)=0;
    virtual void rk3_step3(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double)=0;
    
    virtual void flux_update(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double)=0;
    
    virtual void kinematic_fsf(lexer*, fdm_nhf*, double*, double*, double*,slice&)=0;
    virtual void kinematic_bed(lexer*, fdm_nhf*, double*, double*, double*)=0;
    
    virtual void wetdry(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*, slice&)=0;
    virtual void wetdry_fluxes(lexer*, fdm_nhf*, ghostcell*,slice&,double*,double*,double*,double*,double*,double*)=0;
    
    virtual void ucorr(lexer*, fdm_nhf*, double*, slice&, double)=0;
    virtual void vcorr(lexer*, fdm_nhf*, double*, slice&, double)=0;


};

#endif
