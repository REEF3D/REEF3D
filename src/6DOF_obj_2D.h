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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#ifndef SIXDOF_OBJ_2D_H_
#define SIXDOF_OBJ_2D_H_

#include"6DOF_obj.h"

//  6DOF body on a horizontal (2D) grid: level set of the waterline section and the pressure
//  patch of the hull (ship waves, X 400) for SFLOW and the NHFLOW ship-wave mode (X 10 3).

class sixdof_obj_2D : public sixdof_obj
{
public:
    
    sixdof_obj_2D(lexer*, ghostcell*, int);
	virtual ~sixdof_obj_2D();
    
    void initialize_shipwave(lexer*,ghostcell*,slice&,slice&);
    void update_position_2D(lexer*, ghostcell*,slice&);
    double Hsolidface_2D(lexer*, int,int);
    void updateForcing_box(lexer*, ghostcell*, slice&);
    void updateForcing_stl(lexer*, ghostcell*, slice&, slice&);
    void updateForcing_oned(lexer*, ghostcell*, slice&);
    void update_forcing_sflow(lexer*, ghostcell*, slice&, slice&, slice&, slice&, slice&, slice&, int);
    void solve_eqmotion_sflow(lexer*,ghostcell*,int,bool);
    void solve_eqmotion_oneway_sflow(lexer*,ghostcell*,int,bool);
    slice& amr_fs() {return fs;}

protected:
    
    void geometry_parameters_2D(lexer*, ghostcell*);
    void update_trimesh_2D(lexer*, ghostcell*);
    void ray_cast_2D(lexer*, ghostcell*);
    void ray_cast_2D_io_x(lexer*, ghostcell*,int,int);
    void ray_cast_2D_io_ycorr(lexer*, ghostcell*,int,int);
    void ray_cast_2D_x(lexer*, ghostcell*,int,int);
    void ray_cast_2D_y(lexer*, ghostcell*,int,int);
    void ray_cast_2D_z(lexer*, ghostcell*,int,int);
    void reini_2D(lexer*,ghostcell*,slice&);
    void disc_2D(lexer*,ghostcell*,slice&);
    void time_preproc_2D(lexer*);
    
    slice4 press,lrk1,lrk2,K,dts,fs,Ls,Bs,Rxmin,Rxmax,Rymin,Rymax,draft;
    sliceint5 cl,cr,fsio;
};

#endif
