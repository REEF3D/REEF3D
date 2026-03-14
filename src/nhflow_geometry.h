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

#ifndef NHFLOW_GEOMETRY_H_
#define NHFLOW_GEOMETRY_H_

#include"increment.h"
#include"slice4.h"
#include"vtp3D.h"

class lexer;
class fdm_nhf;
class ghostcell;
class slice;
class nhflow_reinidisc_fsf;

using namespace std;

class nhflow_geometry : public increment, private vtp3D
{
public:
	nhflow_geometry(lexer*, fdm_nhf*, ghostcell*);
	virtual ~nhflow_geometry();
    
    void ray_cast(lexer*, fdm_nhf*, ghostcell*);
    void reini_RK2(lexer*, fdm_nhf*, ghostcell*, double*);
    
    void objects_create(lexer*, ghostcell*);
    
    void objects_create_forcing(lexer*, ghostcell*);
    void objects_create_vrans(lexer*, ghostcell*);
    
private:
    void objects_allocate_forcing(lexer*, ghostcell*);
    void objects_allocate_vrans(lexer*, ghostcell*);
    
    void ray_cast_io(lexer*, fdm_nhf*, ghostcell*,int,int);
    void ray_cast_x(lexer*, fdm_nhf*, ghostcell*,int,int);
    void ray_cast_y(lexer*, fdm_nhf*, ghostcell*,int,int);
    void ray_cast_z(lexer*, fdm_nhf*, ghostcell*,int,int);
    void ray_cast_direct(lexer*, fdm_nhf*, ghostcell*,int,int);
    
    void box(lexer*, ghostcell*, int);
    void cylinder_y(lexer*, ghostcell*, int);
    void cylinder_z(lexer*, ghostcell*, int);
    void jacketmember(lexer*, ghostcell*, int);
    void sphere(lexer*, ghostcell*, int);
    void wedge_x(lexer*, ghostcell*, int);
    void wedge_y(lexer*, ghostcell*, int);
    void wedge_z(lexer*, ghostcell*, int);
    
    void read_stl(lexer*, ghostcell*);
    
    void rotation_tri(lexer*,double,double,double,double&,double&,double&, const double&, const double&, const double&);

    void create_triangle(double&,double&,double&,double&,double&,double&,double&,double&,double&,const double&,const double&,const double&);
    
    void rotation(double&,double&,double&,double,double,double);
    
	void rotate_triangle(lexer*,int,int);
    void rotation_ellipsoid(lexer*,int,double&,double&,double&,double,double,double);
    
    void angle_calc(double,double,double,double&,double&,double&);
    
    void print_vtp(lexer*);
    
    int *IO,*CR,*CL;
    double *FRK1,*dt,*L;
    
    double **tri_x,**tri_y,**tri_z,**tri_x0,**tri_y0,**tri_z0;
    
    int *tstart,*tend;
    int tricount;
    int entity_count, entity_sum;
    double xs,xe,ys,ye,zs,ze;
    double zmin,zmax;
    double starttime;
    
    int reiniter;
    int forcing_flag,solid_flag,floating_flag;
    int dlm_flag;
    
    const double epsi;
    
    nhflow_reinidisc_fsf *prdisc;

    
    double H,Ht, uf, vf, wf, ef;
    double efc;
	double nx, ny, nz,norm ;
	double phival_sf;
    double dirac;
    
    double phi,theta,psi;
    double xrot,yrot,zrot;
    
    double DSM;
    
    
    // -------------------------
    int box_num;
    double *box_xs,*box_xe,*box_ys,*box_ye,*box_zs,*box_ze;
    int cyly_num;
    double *cyly_xc,*cyly_zc,*cyly_ys,*cyly_ye,*cyly_r;
    int cylz_num;
    double *cylz_xc,*cylz_yc,*cylz_zs,*cylz_ze,*cylz_r;
    int jacket_num;
    double *jacket_xm1, *jacket_ym1, *jacket_zm1, *jacket_r1, *jacket_xm2, *jacket_ym2, *jacket_zm2, *jacket_r2;
    int sphere_num;
    double *sphere_xm,*sphere_ym,*sphere_zm,*sphere_r;
    int wedgex_num;
    double *wedgex_xs,*wedgex_xe,*wedgex_ys,*wedgex_ye,*wedgex_zs,*wedgex_ze;
    int wedgey_num;
    double *wedgey_xs,*wedgey_xe,*wedgey_ys,*wedgey_ye,*wedgey_zs,*wedgey_ze;
    int wedgez_num;
    double *wedgez_xs,*wedgez_xe,*wedgez_ys,*wedgez_ye,*wedgez_zs,*wedgez_ze;
    
    int stl_num;
    double stl_scale_x,stl_scale_y,stl_scale_z;
    double stl_trans_x,stl_trans_y,stl_trans_z;
    double stl_orig_x,stl_orig_y,stl_orig_z,stl_phi,stl_orig_theta,stl_orig_psi;
    int stl_invert,stl_dlm;
    // -------------------------
    
    
    

};

#endif
