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

#include"nhflow_fsf.h"
#include"increment.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"

class patchBC_interface;

using namespace std;

#ifndef NHFLOW_FSF_F_H_
#define NHFLOW_FSF_F_H_

class nhflow_fsf_f : public nhflow_fsf, public increment
{
public:
    nhflow_fsf_f(lexer*, fdm_nhf*, ghostcell*,ioflow*,patchBC_interface*);
	virtual ~nhflow_fsf_f();
    
    virtual void start(lexer*, fdm_nhf*, ghostcell*, ioflow*);
    virtual void ini(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*);
    
    virtual void rk2_step1(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double);
    virtual void rk2_step2(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double);
    
    virtual void rk3_step1(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double);
    virtual void rk3_step2(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double);
    virtual void rk3_step3(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double);
    
    virtual void flux_update(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double);
    
    virtual void kinematic_fsf(lexer*, fdm_nhf*, double*, double*, double*,slice&);
    virtual void kinematic_bed(lexer*, fdm_nhf*, double*, double*, double*);
    
    virtual void wetdry(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*, slice&);
    virtual void wetdry_fluxes(lexer*, fdm_nhf*, ghostcell*,slice&,double*,double*,double*,double*,double*,double*);
    
    virtual void breaking(lexer*, fdm_nhf*, ghostcell*,slice&, slice&, double);
    
    virtual void ucorr(lexer*, fdm_nhf*, double*, slice&, double);
    virtual void vcorr(lexer*, fdm_nhf*, double*, slice&, double);
    
    void update(lexer*,fdm_nhf*,ghostcell*,slice&);
    
private: 
    void filter(lexer*, fdm_nhf*, ghostcell*, slice&);
    
    double limiter(double, double);
    
    patchBC_interface *pBC;
    
    slice1 P;
    slice2 Q;
    slice4 K;
    int *temp;

    int gcval_phi,gcval_eta;
	double starttime;
    double phival,H;
	double d;
    double val, denom;
    double dfdx_min, dfdx_plus, dfdy_min, dfdy_plus;
    double detadx,detady;
    
    const double eps;

};

#endif
