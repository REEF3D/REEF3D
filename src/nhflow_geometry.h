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
#include"geo_raycast.h"
#include<vector>

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
    
    void ray_cast(lexer*, fdm_nhf*, ghostcell*, double*);
    void reini_RK2(lexer*, fdm_nhf*, ghostcell*, double*);
    
    void objects_create_forcing(lexer*, ghostcell*);
    void objects_create_vrans(lexer*, ghostcell*);
    
    void geometry_ini(lexer*, fdm_nhf*, ghostcell*);
    
    // VRANS: exposed surface faces (PORSTRUC distance, B 202 layer depth) and depth below them
    void layer_faces(lexer*, fdm_nhf*, ghostcell*);
    void layer_dist(lexer*, fdm_nhf*, double, const double*, double*);
    
    int forcing_flag,solid_flag,floating_flag;
    int dlm_flag;
    
    // mesh refinement (NHFLOW AMR, fixed solids A 581-590): the geometry of a patch casts its
    // solids within the patch box (amr_ray_box: the x and y extent of the patch arrays, as
    // nhflow_amr_6dof for the hulls; the lexer of a patch has the extent of the level-0 rank) and
    // prints no surface file and no messages
    bool amr_grid = false;
    void amr_ray_box(double x0, double x1, double y0, double y1) { amr_grid=true; amr_bx[0]=x0; amr_bx[1]=x1; amr_bx[2]=y0; amr_bx[3]=y1; }
    // the x and y extent of the surface triangles of the solids (entity n < solid_entities())
    int solid_entities() const { return solid_flag==1 ? entity_sum : 0; }
    void entity_box(int n, double *b) const;
    
    
private:
    void objects_allocate_forcing(lexer*, ghostcell*);
    void objects_allocate_vrans(lexer*, ghostcell*);
    
    // ray cast: geometry core kernels
    geo_raycast georay;
    
    // solids of the grid file as immersed solids (A 580 1)
    void grid_solids(lexer*, ghostcell*);
    int grid_solid_num, grid_solid_tri;
    vector<int> ent_raymode, ent_invert;
    
    // porous layers
    bool point_inside(lexer*, int, double, double, double);
    vector<int> tri_layer;
    double layer_dsmax;
    
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
    
    void print_vtp(lexer*,int);
    
    double NB;
    double amr_bx[4] = {0.0,0.0,0.0,0.0};
    
    int *IO,*CR,*CL;
    double *FRK1,*dt,*L;
    
    double DSM;
    
    double **tri_x,**tri_y,**tri_z;
    
    int *tstart,*tend;
    int tricount;
    int entity_count, entity_sum;
    double xs,xe,ys,ye,zs,ze;
    double zmin,zmax;
    double starttime;
    
    int reiniter;
    
    const double epsi;
    
    nhflow_reinidisc_fsf *prdisc;

    
    double phi,theta,psi;
    double xrot,yrot,zrot;
    
    
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
