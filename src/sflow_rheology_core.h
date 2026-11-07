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

#ifndef SFLOW_RHEOLOGY_CORE_H_
#define SFLOW_RHEOLOGY_CORE_H_

#include"sflow_rheology.h"
#include"increment.h"
#include"slice4.h"
#include"rheology_core.h"

// SFLOW adapter of the rheology core (W 120 > 0): single-layer depth-averaged flow of a
// rheological material (debris flow, mudflow, granular avalanche, landslide runout).
// - along-slope gravity W 20, W 21 as body force h g_x, h g_y (inclined reference plane: the bed is
//   the plane, W 22 its normal gravity); with W 20 = W 21 = 0 the bed slope drives the flow as usual
// - basal friction with the depth-integrated closures of the core, as a backward-Euler step on the
//   stage momenta after the explicit update (yield criterion: the layer stops when the driving force
//   of the stage is below the static threshold)
// - basal pore pressure ratio lambda (W 126): constant, or advected with the flow and relaxing to
//   W 127 lambda_h on the diffusion time scale 4 h^2/(pi^2 D)

class sflow_rheology_core final : public sflow_rheology, public increment
{
public:
    sflow_rheology_core(lexer*);
	virtual ~sflow_rheology_core() = default;
    
	void u_source(lexer*, fdm2D*, slice&, slice&, slice&) override final;
    void v_source(lexer*, fdm2D*, slice&, slice&, slice&) override final;
    
    void ini(lexer*, fdm2D*, ghostcell*) override final;
    void implicit(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double, bool) override final;
    
    slice4 lambda;
    
private:
    double gn(lexer*, fdm2D*);
    void pore_update(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&);
    
    rheo::basal_params q;
    slice4 lnew;
};

#endif
