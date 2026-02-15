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

#include"lexer.h"
#include"gridfile_v2.h"
#include"geo_mesh.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

// DIVEMesh grid, format v2 (see gridfile_v2.h)
//
// per rank: DIVEMesh_Grid/grid-%06i.dat   cell flags, nodes, boundary and parallel surfaces,
//                                         column flags, geodat bed level, data
// all ranks: DIVEMesh_Grid/grid-geometry.dat  solid (S) and topo (T) entities as triangles
//
// The solid and topo fields, the bed levels and the ghost cell estimates are built by
// REEF3D from the triangles (lexer::grid_solids, geometry core), not read.

void lexer::read_grid()
{
    int i,j,k,n,q;
    char name[100];
    
    gcwall_count=0;
    gcin_count=0;
    gcout_count=0;
    gcfsf_count=0;
    gcbed_count=0;
    
    gcpara1_count=gcpara2_count=gcpara3_count=gcpara4_count=gcpara5_count=gcpara6_count=0;
    gcparaco1_count=gcparaco2_count=gcparaco3_count=gcparaco4_count=gcparaco5_count=gcparaco6_count=0;
    
    surf_tot=0;
    
    const int padding = 6;
    snprintf(name,sizeof(name),"DIVEMesh_Grid/grid-%0*i.dat",padding,mpirank+1);
    
    gridv2::reader gf;
    
    if(!gf.open(name,gridv2::magic_grid))
    {
        cout<<endl;
        cout<<"!!! "<<gf.error<<" !!!"<<endl;
        cout<<"!!! the grid has to be generated with DIVEMesh, grid format v2; please check the manual!"<<endl<<endl<<endl<<endl;
        exit(1);
    }
    
    gridv2::section sc;
    
    // --------------------------------------------------------------------------------------------
    // HEAD
    if(!gf.find("HEAD",sc))
    gf.fail(this,"section HEAD missing");
    
    vector<int> iv;
    vector<double> dv;
    sc.get_list(iv,dv);
    gf.check(sc);
    
    if(!sc.good || iv.size()<62 || dv.size()<19)
    gf.fail(this,"section HEAD too short");
    
    q=0;
    const int DM_M10 = iv[q++];
    
    if(mpirank==0)
    if(DM_M10!=M10 || M10!=mpi_size || DM_M10!=mpi_size)
    {
        cout<<endl;
        cout<<"!!! Inconsistent M 10 parameter, needs to be the same in REEF3D and DIVEMesh !"<<endl;
        cout<<"mpi_size: "<<mpi_size<<" REEFD M10: "<<M10<<" DIVEMesh M10: "<<DM_M10<<endl;
        cout<<"!!! please check the manual!"<<endl<<endl<<endl<<endl;
        exit(1);
    }
    
    knox = iv[q++];
    knoy = iv[q++];
    knoz = iv[q++];
    
    gknox = iv[q++];
    gknoy = iv[q++];
    gknoz = iv[q++];
    
    origin_i = iv[q++];
    origin_j = iv[q++];
    origin_k = iv[q++];
    
    gcwall_count = iv[q++];
    
    gcpara1_count = iv[q++];
    gcpara2_count = iv[q++];
    gcpara3_count = iv[q++];
    gcpara4_count = iv[q++];
    gcpara5_count = iv[q++];
    gcpara6_count = iv[q++];
    
    gcparaco1_count = iv[q++];
    gcparaco2_count = iv[q++];
    gcparaco3_count = iv[q++];
    gcparaco4_count = iv[q++];
    gcparaco5_count = iv[q++];
    gcparaco6_count = iv[q++];
    
    gcslpara1_count = iv[q++];
    gcslpara2_count = iv[q++];
    gcslpara3_count = iv[q++];
    gcslpara4_count = iv[q++];
    
    gcslparaco1_count = iv[q++];
    gcslparaco2_count = iv[q++];
    gcslparaco3_count = iv[q++];
    gcslparaco4_count = iv[q++];
    
    nb1 = iv[q++];
    nb2 = iv[q++];
    nb3 = iv[q++];
    nb4 = iv[q++];
    nb5 = iv[q++];
    nb6 = iv[q++];
    
    mx = iv[q++];
    my = iv[q++];
    mz = iv[q++];
    
    bcside1 = iv[q++];
    bcside2 = iv[q++];
    bcside3 = iv[q++];
    bcside4 = iv[q++];
    bcside5 = iv[q++];
    bcside6 = iv[q++];
    
    periodic1 = iv[q++];
    periodic2 = iv[q++];
    periodic3 = iv[q++];
    
    periodicX1 = iv[q++];
    periodicX2 = iv[q++];
    periodicX3 = iv[q++];
    periodicX4 = iv[q++];
    periodicX5 = iv[q++];
    periodicX6 = iv[q++];
    
    i_dir = iv[q++];
    j_dir = iv[q++];
    k_dir = iv[q++];
    
    P150 = iv[q++];
    cms_flag = iv[q++];
    
    const int DM_marge = iv[q++];
    const int DM_rank = iv[q++];
    
    if(DM_marge!=marge)
    gf.fail(this,"node margin of the grid file differs from REEF3D");
    
    if(DM_rank!=mpirank+1)
    gf.fail(this,"grid file belongs to another rank");
    
    q=0;
    dx  = dv[q++];
    DXM = dv[q++];
    DYM = dv[q++];
    DZM = dv[q++];
    
    originx = dv[q++];
    originy = dv[q++];
    originz = dv[q++];
    endx = dv[q++];
    endy = dv[q++];
    endz = dv[q++];
    
    global_xmin = dv[q++];
    global_ymin = dv[q++];
    global_zmin = dv[q++];
    global_xmax = dv[q++];
    global_ymax = dv[q++];
    global_zmax = dv[q++];
    
    global_orig_x = dv[q++];
    global_orig_y = dv[q++];
    alpha_grid = dv[q++];
    
    // --------------------------------------------------------------------------------------------
    // geometry: all ranks
    gridgeo = new geo_mesh();
    gridgeo->read(this,"DIVEMesh_Grid/grid-geometry.dat");
    
    solidread = gridgeo->solidread;
    toporead = gridgeo->toporead;
    
    // estimated by lexer::grid_solids
    solid_gcb_est = topo_gcb_est = 0;
    solid_gcbextra_est = topo_gcbextra_est = tot_gcbextra_est = 0;
    
    gcb1_count=gcb2_count=gcb3_count=gcb4_count=gcb4a_count=gcb_fix=gcb_solid=gcb_topo=gcb_fb=gcwall_count;
    
    gcpara_sum=gcpara1_count+gcpara2_count+gcpara3_count+gcpara4_count+gcpara5_count+gcpara6_count;
    gcparaco_sum=gcparaco1_count+gcparaco2_count+gcparaco3_count+gcparaco4_count+gcparaco5_count+gcparaco6_count;
    
    grid::assign_margin();
    
    Iarray(flag4,imax*jmax*kmax);
    Darray(flag_solid,imax*jmax*kmax);
    Darray(flag_topo,imax*jmax*kmax);
    Darray(solidbed,imax*jmax);
    Darray(topobed,imax*jmax);
    Darray(geobed,imax*jmax);
    Darray(bed,imax*jmax);
    Iarray(wet,imax*jmax);
    Iarray(wet_n,imax*jmax);
    Iarray(deep,imax*jmax);
    Darray(depth,imax*jmax);
    Darray(WL,imax*jmax);
    Darray(data,imax*jmax);
    Iarray(flagslice1,imax*jmax);
    Iarray(flagslice2,imax*jmax);
    Iarray(flagslice4,imax*jmax);
    
    for(n=0;n<imax*jmax*kmax;++n)
    flag4[n]=-1;
    
    for(n=0;n<imax*jmax*kmax;++n)
    flag_solid[n]=0.0;
    
    for(n=0;n<imax*jmax;++n)
    {
    flagslice1[n]=-10;
    flagslice2[n]=-10;
    flagslice4[n]=-10;
    }
    
    if(gcb4_count>0)
    {
    Iarray(gcb1, gcb1_count,6);
    Iarray(gcb2, gcb2_count,6);
    Iarray(gcb3, gcb3_count,6);
    Iarray(gcb4, gcb4_count,6);
    Iarray(gcb4a, gcb4a_count,6);
    
    Darray(gcd1, gcb1_count);
    Darray(gcd2, gcb2_count);
    Darray(gcd3, gcb3_count);
    Darray(gcd4, gcb4_count);
    Darray(gcd4a, gcb4a_count);
    }
    
    Iarray(gcpara1, gcpara1_count,16);
    Iarray(gcpara2, gcpara2_count,16);
    Iarray(gcpara3, gcpara3_count,16);
    Iarray(gcpara4, gcpara4_count,16);
    Iarray(gcpara5, gcpara5_count,16);
    Iarray(gcpara6, gcpara6_count,16);
    
    Iarray(gcparaco1, gcparaco1_count,3);
    Iarray(gcparaco2, gcparaco2_count,3);
    Iarray(gcparaco3, gcparaco3_count,3);
    Iarray(gcparaco4, gcparaco4_count,3);
    Iarray(gcparaco5, gcparaco5_count,3);
    Iarray(gcparaco6, gcparaco6_count,3);
    
    gcbsl1_count=gcbsl2_count=gcbsl3_count=gcbsl4_count=gcbsl4a_count=1;
    
    Iarray(gcbsl1, gcbsl1_count,6);
    Iarray(gcbsl2, gcbsl2_count,6);
    Iarray(gcbsl3, gcbsl3_count,6);
    Iarray(gcbsl4, gcbsl4_count,6);
    Iarray(gcbsl4a, gcbsl4a_count,6);
    
    Iarray(gcslin, gcin_count,6);
    Iarray(gcslout, gcout_count,6);
    
    Iarray(gcslpara1, gcslpara1_count,2);
    Iarray(gcslpara2, gcslpara2_count,2);
    Iarray(gcslpara3, gcslpara3_count,2);
    Iarray(gcslpara4, gcslpara4_count,2);
    
    Iarray(gcslparaco1, gcslparaco1_count,4);
    Iarray(gcslparaco2, gcslparaco2_count,4);
    Iarray(gcslparaco3, gcslparaco3_count,4);
    Iarray(gcslparaco4, gcslparaco4_count,4);
    
    Darray(XN,knox+1+2*marge);
    Darray(YN,knoy+1+2*marge);
    Darray(ZN,knoz+1+2*marge);
    
    // --------------------------------------------------------------------------------------------
    // FLAG
    {
    vector<int> fl;
    
    if(!gf.find("FLAG",sc) || !sc.get_rle(fl,size_t(knox)*size_t(knoy)*size_t(knoz)))
    gf.fail(this,"section FLAG missing or inconsistent");
    
    size_t m=0;
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    for(k=0; k<knoz; ++k)
    flag4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin] = fl[m++];
    }
    
    // --------------------------------------------------------------------------------------------
    // NODE
    if(!gf.find("NODE",sc))
    gf.fail(this,"section NODE missing");
    
    for(i=-marge;i<knox+1+marge;++i)
    XN[IP]=sc.get_double();
    
    for(j=-marge;j<knoy+1+marge;++j)
    YN[JP]=sc.get_double();
    
    for(k=-marge;k<knoz+1+marge;++k)
    ZN[KP]=sc.get_double();
    
    if(!sc.good || !sc.done())
    gf.fail(this,"section NODE inconsistent");
    
    // --------------------------------------------------------------------------------------------
    // SURF: boundary surfaces
    if(!gf.find("SURF",sc))
    gf.fail(this,"section SURF missing");
    
    gcin_count=0;
    gcout_count=0;
    
    {
    vector<int> tb;
    
    if(!sc.get_table(tb,gcb4_count,5))
    gf.fail(this,"section SURF inconsistent");
    
    for(n=0; n<gcb4_count; ++n)
    {
        for(q=0; q<5; ++q)
        gcb4[n][q] = tb[5*n+q];     // i j k side group
        
        if(gcb4[n][4]==1 || gcb4[n][4]==6)
        ++gcin_count;
        
        if(gcb4[n][4]==2 || gcb4[n][4]==7 || gcb4[n][4]==8)
        ++gcout_count;
    }
    }
    
    if(!sc.good || !sc.done())
    gf.fail(this,"section SURF inconsistent");
    
    Iarray(gcin, gcin_count,6);
    Iarray(gcout, gcout_count,6);
    
    // --------------------------------------------------------------------------------------------
    // PARA: parallel surfaces
    {
    int **gcpara[6] = {gcpara1,gcpara2,gcpara3,gcpara4,gcpara5,gcpara6};
    const int num[6] = {gcpara1_count,gcpara2_count,gcpara3_count,gcpara4_count,gcpara5_count,gcpara6_count};
    
    if(!gf.find("PARA",sc))
    gf.fail(this,"section PARA missing");
    
    vector<int> tb;
    
    for(int d=0; d<6; ++d)
    {
        if(!sc.get_table(tb,num[d],3))
        gf.fail(this,"section PARA inconsistent");
        
        for(n=0; n<num[d]; ++n)
        {
        gcpara[d][n][0]=tb[3*n];
        gcpara[d][n][1]=tb[3*n+1];
        gcpara[d][n][2]=tb[3*n+2];
        gcpara[d][n][3]=1;
        }
    }
    
    if(!sc.good || !sc.done())
    gf.fail(this,"section PARA inconsistent");
    }
    
    // PACO: parallel corners
    {
    int **gcparaco[6] = {gcparaco1,gcparaco2,gcparaco3,gcparaco4,gcparaco5,gcparaco6};
    const int num[6] = {gcparaco1_count,gcparaco2_count,gcparaco3_count,gcparaco4_count,gcparaco5_count,gcparaco6_count};
    
    if(!gf.find("PACO",sc))
    gf.fail(this,"section PACO missing");
    
    vector<int> tb;
    
    for(int d=0; d<6; ++d)
    {
        if(!sc.get_table(tb,num[d],3))
        gf.fail(this,"section PACO inconsistent");
        
        for(n=0; n<num[d]; ++n)
        {
        gcparaco[d][n][0]=tb[3*n];
        gcparaco[d][n][1]=tb[3*n+1];
        gcparaco[d][n][2]=tb[3*n+2];
        }
    }
    
    if(!sc.good || !sc.done())
    gf.fail(this,"section PACO inconsistent");
    }
    
    // --------------------------------------------------------------------------------------------
    // SLFL: column flags
    {
    vector<int> fl;
    
    if(!gf.find("SLFL",sc) || !sc.get_rle(fl,size_t(knox)*size_t(knoy)))
    gf.fail(this,"section SLFL missing or inconsistent");
    
    size_t m=0;
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    flagslice4[(i-imin)*jmax + (j-jmin)] = fl[m++];
    }
    
    // SLPA: parallel slice surfaces
    {
    int **gcslpara[4] = {gcslpara1,gcslpara2,gcslpara3,gcslpara4};
    const int num[4] = {gcslpara1_count,gcslpara2_count,gcslpara3_count,gcslpara4_count};
    
    if(!gf.find("SLPA",sc))
    gf.fail(this,"section SLPA missing");
    
    vector<int> tb;
    
    for(int d=0; d<4; ++d)
    {
        if(!sc.get_table(tb,num[d],2))
        gf.fail(this,"section SLPA inconsistent");
        
        for(n=0; n<num[d]; ++n)
        {
        gcslpara[d][n][0]=tb[2*n];
        gcslpara[d][n][1]=tb[2*n+1];
        }
    }
    
    if(!sc.good || !sc.done())
    gf.fail(this,"section SLPA inconsistent");
    }
    
    // SLPC: parallel slice corners
    {
    int **gcslparaco[4] = {gcslparaco1,gcslparaco2,gcslparaco3,gcslparaco4};
    const int num[4] = {gcslparaco1_count,gcslparaco2_count,gcslparaco3_count,gcslparaco4_count};
    
    if(!gf.find("SLPC",sc))
    gf.fail(this,"section SLPC missing");
    
    vector<int> tb;
    
    for(int d=0; d<4; ++d)
    {
        if(!sc.get_table(tb,num[d],2))
        gf.fail(this,"section SLPC inconsistent");
        
        for(n=0; n<num[d]; ++n)
        {
        gcslparaco[d][n][0]=tb[2*n];
        gcslparaco[d][n][1]=tb[2*n+1];
        }
    }
    
    if(!sc.good || !sc.done())
    gf.fail(this,"section SLPC inconsistent");
    }
    
    // --------------------------------------------------------------------------------------------
    // GEOB: geodat bed level
    for(n=0;n<imax*jmax;++n)
    geobed[n]=global_zmin;
    
    if(gridgeo->geodat>0)
    {
        if(!gf.find("GEOB",sc))
        gf.fail(this,"section GEOB missing");
        
        for(i=0; i<knox; ++i)
        for(j=0; j<knoy; ++j)
        geobed[(i-imin)*jmax + (j-jmin)]=sc.get_double();
        
        if(!sc.good || !sc.done())
        gf.fail(this,"section GEOB inconsistent");
    }
    
    // DATA
    if(P150>0)
    {
        if(!gf.find("DATA",sc))
        gf.fail(this,"section DATA missing");
        
        for(i=0; i<knox; ++i)
        for(j=0; j<knoy; ++j)
        data[(i-imin)*jmax + (j-jmin)]=sc.get_double();
        
        if(!sc.good || !sc.done())
        gf.fail(this,"section DATA inconsistent");
    }
}
