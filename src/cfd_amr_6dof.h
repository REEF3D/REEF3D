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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef CFD_AMR_6DOF_H_
#define CFD_AMR_6DOF_H_

#include"6DOF.h"
#include"increment.h"
#include"fieldint4.h"
#include"field4a.h"

class sixdof_cfd;
class reinidisc;

using namespace std;

//  The floating bodies (6DOF_cfd, X 10 1) on a patch of the CFD mesh refinement (cfd_amr): the
//  sixdof object of the patch.  The bodies move on level 0 (sixdof_cfd: loads, equations of
//  motion, before the patches take their stage); a patch casts the hull on its own grid (level set
//  fb with the cells around the patch, ray-cast and reinitialisation workspace of the patch size)
//  and adds the direct forcing of the rigid-body velocity, as sixdof_cfd::start_cfd does it on
//  level 0.  With subcycling (G 7 1) cfd_amr sets the body to the time of the stage first.

class cfd_amr_6dof final : public sixdof, public increment
{
public:
    cfd_amr_6dof(lexer*, sixdof_cfd*);
    virtual ~cfd_amr_6dof();

    // the hull on the grid without forcing: fb, the flags of the solid boundary (patch set-up)
    void body(lexer*, fdm*, ghostcell*);

    void start_cfd(lexer*,fdm*,ghostcell*,int,field&,field&,field&,field&,field&,field&,bool) override final;
    void start_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool) override final {}
    void start_sflow(lexer*,fdm2D*,ghostcell*,int,slice&,slice&,slice&,slice&,slice&,slice&,slice&,bool) override final {}

    void reforce_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool) override final {}

    void ini(lexer*,ghostcell*) override final {}
    void initialize(lexer*, fdm*, ghostcell*) override final {}
    void initialize(lexer*, fdm2D*, ghostcell*) override final {}
    void initialize(lexer*, fdm_nhf*, ghostcell*) override final {}

    void isource(lexer*,fdm*,ghostcell*) override final {}
    void jsource(lexer*,fdm*,ghostcell*) override final {}
    void ksource(lexer*,fdm*,ghostcell*) override final {}

    void isource(lexer*,fdm_nhf*,ghostcell*,slice&) override final {}
    void jsource(lexer*,fdm_nhf*,ghostcell*,slice&) override final {}
    void ksource(lexer*,fdm_nhf*,ghostcell*,slice&) override final {}

    void isource2D(lexer*,fdm2D*,ghostcell*) override final {}
    void jsource2D(lexer*,fdm2D*,ghostcell*) override final {}

private:
    void cast(lexer*, fdm*, ghostcell*, int);
    void reset_fbh(lexer*, fdm*, ghostcell*);

    sixdof_cfd *b;
    fieldint4 cutl, cutr, fbio;
    field4a frk1, L, dt;
    reinidisc *prdisc;
};

#endif
