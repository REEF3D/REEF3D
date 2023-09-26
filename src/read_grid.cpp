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

#include"lexer.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void lexer::read_grid()
{
	int i,n;
	int isurf,jsurf,ksurf,surfside,surfgroup,side1,side2,paraconum;
	int para_active;
    char name[100];
	int DM_M10;
    int iin;
    double ddn;

	gcwall_count=0;
	gcin_count=0;
	gcout_count=0;
    gcin6_count=0;
	gcout6_count=0;
	gcfsf_count=0;
	gcbed_count=0;
	gcpara1_count=0;
	gcpara2_count=0;
	gcpara3_count=0;
	gcpara4_count=0;
	gcpara5_count=0;
	gcpara6_count=0;
	gcparaco1_count=0;
	gcparaco2_count=0;
	gcparaco3_count=0;
	gcparaco4_count=0;
	gcparaco5_count=0;
	gcparaco6_count=0;
	surf_tot=0;

	if(mpirank<9)
	sprintf(name,"grid-00000%i.dat",mpirank+1);

	if(mpirank<99&&mpirank>8)
	sprintf(name,"grid-0000%i.dat",mpirank+1);

	if(mpirank<999&&mpirank>98)
	sprintf(name,"grid-000%i.dat",mpirank+1);

	if(mpirank<9999&&mpirank>998)
	sprintf(name,"grid-00%i.dat",mpirank+1);

	if(mpirank<99999&&mpirank>9998)
	sprintf(name,"grid-0%i.dat",mpirank+1);
	
	if(mpirank>99998)
	sprintf(name,"grid-%i.dat",mpirank+1);
    
    
    

// open file------------
	ifstream grid(name, ios_base::binary);
	
//read grid file-------------
    
    grid.read((char*)&iin, sizeof (int));
    DM_M10=iin;
	
	if(mpirank==0)
	if(DM_M10!=M10 || M10!=mpi_size || DM_M10!=mpi_size)
    {
    cout<<endl;
    cout<<"!!! Inconsistent M 10 parameter, needs to be the same in REEF3D and DIVEMesh !"<<endl;
    cout<<"mpi_size: "<<mpi_size<<" REEFD M10: "<<M10<<" DIVEMesh M10: "<<DM_M10<<endl;
    cout<<"!!! please check the manual!"<<endl<<endl<<endl<<endl;
    
    exit(0);
    }
    
    grid.read((char*)&iin, sizeof (int));
    knox=iin;
    grid.read((char*)&iin, sizeof (int));
    knoy=iin;
    grid.read((char*)&iin, sizeof (int));
    knoz=iin;
    
    grid.read((char*)&ddn, sizeof (double));
    dx=ddn;
    
    grid.read((char*)&ddn, sizeof (double));
    DRM=ddn;
    grid.read((char*)&ddn, sizeof (double));
    DSM=ddn;
    grid.read((char*)&ddn, sizeof (double));
    DTM=ddn;
    
    
    grid.read((char*)&ddn, sizeof (double));
    originx=ddn;
    grid.read((char*)&ddn, sizeof (double));
    originy=ddn;
    grid.read((char*)&ddn, sizeof (double));
    originz=ddn;
	
    grid.read((char*)&ddn, sizeof (double));
    endx=ddn;
    grid.read((char*)&ddn, sizeof (double));
    endy=ddn;
    grid.read((char*)&ddn, sizeof (double));
    endz=ddn;

    grid.read((char*)&ddn, sizeof (double));
    global_xmin=ddn;
    grid.read((char*)&ddn, sizeof (double));
    global_ymin=ddn;
    grid.read((char*)&ddn, sizeof (double));
    global_zmin=ddn;
    grid.read((char*)&ddn, sizeof (double));
    global_xmax=ddn;
    grid.read((char*)&ddn, sizeof (double));
    global_ymax=ddn;
    grid.read((char*)&ddn, sizeof (double));
    global_zmax=ddn;

    grid.read((char*)&iin, sizeof (int));
    gknox=iin;
    grid.read((char*)&iin, sizeof (int));
    gknoy=iin;
    grid.read((char*)&iin, sizeof (int));
    gknoz=iin;
    
    grid.read((char*)&iin, sizeof (int));
    origin_i=iin;
    grid.read((char*)&iin, sizeof (int));
    origin_j=iin;
    grid.read((char*)&iin, sizeof (int));
    origin_k=iin;
    
    grid.read((char*)&iin, sizeof (int));
    gcwall_count=iin;
    
    grid.read((char*)&iin, sizeof (int));
    gcpara1_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcpara2_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcpara3_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcpara4_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcpara5_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcpara6_count=iin;
    
    grid.read((char*)&iin, sizeof (int));
    gcparaco1_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcparaco2_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcparaco3_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcparaco4_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcparaco5_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcparaco6_count=iin;
    
    grid.read((char*)&iin, sizeof (int));
    gcslpara1_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslpara2_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslpara3_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslpara4_count=iin;

    grid.read((char*)&iin, sizeof (int));
    gcslparaco1_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslparaco2_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslparaco3_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslparaco4_count=iin;
    
    grid.read((char*)&iin, sizeof (int));
    nb1=iin;
    grid.read((char*)&iin, sizeof (int));
    nb2=iin;
    grid.read((char*)&iin, sizeof (int));
    nb3=iin;
    grid.read((char*)&iin, sizeof (int));
    nb4=iin;
    grid.read((char*)&iin, sizeof (int));
    nb5=iin;
    grid.read((char*)&iin, sizeof (int));
    nb6=iin;
    

	grid.read((char*)&iin, sizeof (int));
    mx=iin;
    grid.read((char*)&iin, sizeof (int));
    my=iin;
    grid.read((char*)&iin, sizeof (int));
    mz=iin;
    
    grid.read((char*)&iin, sizeof (int));
    mi=iin;
    grid.read((char*)&iin, sizeof (int));
    mj=iin;
    grid.read((char*)&iin, sizeof (int));
    mk=iin;
    
    grid.read((char*)&iin, sizeof (int));
    bcside1=iin;
    grid.read((char*)&iin, sizeof (int));
    bcside2=iin;
    grid.read((char*)&iin, sizeof (int));
    bcside3=iin;
    grid.read((char*)&iin, sizeof (int));
    bcside4=iin;
    grid.read((char*)&iin, sizeof (int));
    bcside5=iin;
    grid.read((char*)&iin, sizeof (int));
    bcside6=iin;
    
    grid.read((char*)&iin, sizeof (int));
    periodic1=iin;
    grid.read((char*)&iin, sizeof (int));
    periodic2=iin;
    grid.read((char*)&iin, sizeof (int));
    periodic3=iin;
    
    grid.read((char*)&iin, sizeof (int));
    periodicX1=iin;
    grid.read((char*)&iin, sizeof (int));
    periodicX2=iin;
    grid.read((char*)&iin, sizeof (int));
    periodicX3=iin;
    grid.read((char*)&iin, sizeof (int));
    periodicX4=iin;
    grid.read((char*)&iin, sizeof (int));
    periodicX5=iin;
    grid.read((char*)&iin, sizeof (int));
    periodicX6=iin;
    
    grid.read((char*)&iin, sizeof (int));
    i_dir=iin;
    grid.read((char*)&iin, sizeof (int));
    j_dir=iin;
    grid.read((char*)&iin, sizeof (int));
    k_dir=iin;
    
    grid.read((char*)&iin, sizeof (int));
    P150=iin;
    
    grid.read((char*)&iin, sizeof (int));
    solidread=iin; // solid
    grid.read((char*)&iin, sizeof (int));
    toporead=iin; // topo
    grid.read((char*)&iin, sizeof (int));
    solid_gcb_est=iin;
    grid.read((char*)&iin, sizeof (int));
    topo_gcb_est=iin;
    
    grid.read((char*)&iin, sizeof (int));
    solid_gcbextra_est=iin;
    grid.read((char*)&iin, sizeof (int));
    topo_gcbextra_est=iin;
    grid.read((char*)&iin, sizeof (int));
    tot_gcbextra_est=iin;
    
    grid.read((char*)&iin, sizeof (int));
    grid.read((char*)&iin, sizeof (int));
    grid.read((char*)&iin, sizeof (int));
    grid.read((char*)&iin, sizeof (int));
    grid.read((char*)&iin, sizeof (int));

    
// ---------------------------------------------------------------------------------------------------------------------	
// ---------------------------------------------------------------------------------------------------------------------	
   
    topo_gcb_est*=4;
		
	gcb1_count=gcb2_count=gcb3_count=gcb4_count=gcb4a_count=gcb_fix=gcb_solid=gcb_topo=gcb_fb=gcwall_count;
	
    gcpara_sum=gcpara1_count+gcpara2_count+gcpara3_count+gcpara4_count+gcpara5_count+gcpara6_count;
    gcparaco_sum=gcparaco1_count+gcparaco2_count+gcparaco3_count+gcparaco4_count+gcparaco5_count+gcparaco6_count;

	maxpara=maxparacount();
    assign_margin();
	
	Iarray(flag4,imax*jmax*kmax);
    
    //if(solidread==1)
	Darray(flag_solid,imax*jmax*kmax);
    
    //if(toporead==1)
    Darray(flag_topo,imax*jmax*kmax);
    
	Iarray(mgflag,imax*jmax*kmax);
	Darray(solidbed,imax*jmax);
    Darray(topobed,imax*jmax);
    Darray(bed,imax*jmax);
    Iarray(wet,imax*jmax);
    Iarray(wet_n,imax*jmax);
    Iarray(deep,imax*jmax);
    Darray(depth,imax*jmax);
	Darray(data,imax*jmax);
    Iarray(flagslice1,imax*jmax);
    Iarray(flagslice2,imax*jmax);
    Iarray(flagslice4,imax*jmax);
	Iarray(tpflagslice,imax*jmax);

	for(i=0;i<imax*jmax*kmax;++i)
	flag4[i]=-1;
	
	for(i=0;i<imax*jmax*kmax;++i)
	flag_solid[i]=0.0;
	
	for(i=0;i<imax*jmax;++i)
	{
	flagslice1[i]=-10;
	flagslice2[i]=-10;
	flagslice4[i]=-10;
	}
	
	if(gcb4_count>0)
	{
    Iarray(gcb1, gcb1_count,6);
    Iarray(gcb2, gcb2_count,6);
    Iarray(gcb3, gcb3_count,6);
    Iarray(gcb4, gcb4_count,6);
    Iarray(gcb4a, gcb4a_count,6);
    Iarray(gcb6, gcb4_count);
	
	Iarray(gcside4, gcb4_count);
	gcside4_size=gcb4_count;

    Darray(gcd1, gcb1_count);
    Darray(gcd2, gcb2_count);
    Darray(gcd3, gcb3_count);
    Darray(gcd4, gcb4_count);
    Darray(gcd4a, gcb4a_count);
	}
    
    if(periodic1==1||periodic2==1||periodic3==1)
    {
    gc4periodic_maxcount = 0;
    gc4periodic_maxcount = MAX(knox*knoy,knox*knoz);
    gc4periodic_maxcount = MAX(gc4periodic_maxcount,knoy*knoz);
    
    Iarray(gc4periodic_count,6);
    Iarray(gc4aperiodic_count,6);
    
    Iarray(gc4periodic,6,gc4periodic_maxcount);
    Iarray(gc4aperiodic,6,gc4periodic_maxcount);
    }
	
    Iarray(gcpara1, gcpara1_count,16);
    Iarray(gcpara2, gcpara2_count,16);
    Iarray(gcpara3, gcpara3_count,16);
    Iarray(gcpara4, gcpara4_count,16);
    Iarray(gcpara5, gcpara5_count,16);
    Iarray(gcpara6, gcpara6_count,16);

    Iarray(gcparaco1, gcparaco1_count,6);
    Iarray(gcparaco2, gcparaco2_count,6);
    Iarray(gcparaco3, gcparaco3_count,6);
    Iarray(gcparaco4, gcparaco4_count,6);
    Iarray(gcparaco5, gcparaco5_count,6);
    Iarray(gcparaco6, gcparaco6_count,6);
    
    
    // Slice allocation
    gcbsl1_count=gcbsl2_count=gcbsl3_count=gcbsl4_count=gcbsl4a_count=1;
    
    Iarray(gcbsl1, gcbsl1_count,6);
    Iarray(gcbsl2, gcbsl2_count,6);
    Iarray(gcbsl3, gcbsl3_count,6);
    Iarray(gcbsl4, gcbsl4_count,6);
    Iarray(gcbsl4a, gcbsl4a_count,6);
    
    Iarray(gcslin, gcin_count,6);
    Iarray(gcslout, gcout_count,6);
    
    Darray(gcdsl1, gcbsl1_count);
    Darray(gcdsl2, gcbsl2_count);
    Darray(gcdsl3, gcbsl3_count);
    Darray(gcdsl4, gcbsl4_count);
    Darray(gcdsl4a, gcbsl4a_count);
    
    Iarray(gcslpara1, gcslpara1_count,14);
    Iarray(gcslpara2, gcslpara2_count,14);
    Iarray(gcslpara3, gcslpara3_count,14);
    Iarray(gcslpara4, gcslpara4_count,14);
    
    Iarray(gcslparaco1, gcslparaco1_count,4);
    Iarray(gcslparaco2, gcslparaco2_count,4);
    Iarray(gcslparaco3, gcslparaco3_count,4);
    Iarray(gcslparaco4, gcslparaco4_count,4);
    
    Darray(XN,knox+1+4*marge);
    Darray(YN,knoy+1+4*marge);
    Darray(ZN,knoz+1+4*marge);
    
    Darray(RN,knox+1+4*marge);
    Darray(SN,knoy+1+4*marge);
    Darray(TN,knoz+1+4*marge);
    
	
// ---------------------------------------------------------------------------------------------------------------------		
// ---------------------------------------------------------------------------------------------------------------------	
	
//  Flag
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    for(k=0; k<knoz; ++k)
    {
    grid.read((char*)&iin, sizeof (int));
    flag4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=iin;
    }
	
	for(i=0;i<imax*jmax*kmax;++i)
    mgflag[i]=flag4[i];
    
// Nodes XYZ
    for(i=-marge;i<knox+1+marge;++i)
    {
    grid.read((char*)&ddn, sizeof (double));
    XN[IP]=ddn;
    }

    for(j=-marge;j<knoy+1+marge;++j)
    {
    grid.read((char*)&ddn, sizeof (double));
    YN[JP]=ddn;
    }
    
    for(k=-marge;k<knoz+1+marge;++k)
    {
    grid.read((char*)&ddn, sizeof (double));
    ZN[KP]=ddn;
    }
    
    // Nodes RST
    for(i=-marge;i<knox+1+marge;++i)
    {
    grid.read((char*)&ddn, sizeof (double));
    RN[IP]=ddn;
    }

    for(j=-marge;j<knoy+1+marge;++j)
    {
    grid.read((char*)&ddn, sizeof (double));
    SN[JP]=ddn;
    }
    
    for(k=-marge;k<knoz+1+marge;++k)
    {
    grid.read((char*)&ddn, sizeof (double));
    TN[KP]=ddn;
    }
    
//  Solid	
	if(solidread==1)
	for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    for(k=0; k<knoz; ++k)
    {
    grid.read((char*)&ddn, sizeof (double));
    flag_solid[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=ddn;
    }
    
//  Topo
	if(toporead==1)
	for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    for(k=0; k<knoz; ++k)
    {
    grid.read((char*)&ddn, sizeof (double));
    flag_topo[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=ddn;
    }
    
// Porous Structure
    
	
//  GC Surfaces
    gcin_count=0;
    gcout_count=0;
	for(i=0; i<gcb4_count; ++i)
	{     
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        surfside=iin;
        
        grid.read((char*)&iin, sizeof (int));
        surfgroup=iin;
    
			gcb4[i][0]=isurf;
			gcb4[i][1]=jsurf;
			gcb4[i][2]=ksurf;
			gcb4[i][3]=surfside;
			gcb4[i][4]=surfgroup;

			if(surfgroup==1 || surfgroup==6)
			++gcin_count;

			if(surfgroup==2 || surfgroup==7 || surfgroup==8)
			++gcout_count;
	}
    
    gcin4a_count=gcin_count;
	gcout4a_count=gcout_count;
    
    gcin6_count=gcin_count;
	gcout6_count=gcout_count;
    
	Iarray(gcin, gcin_count,6);
    Iarray(gcout, gcout_count,6);
    
    Iarray(gcin4a, gcin_count,6);
    Iarray(gcout4a, gcout_count,6);
    
    Iarray(gcin6, gcin6_count,6);
    Iarray(gcout6, gcout6_count,6);
    
    int count1=0;
    int count2=0;
    for(i=0; i<gcb4_count; ++i)
    {
        if(gcb4[i][4]==1 || gcb4[i][4]==6)
        {
        gcin6[count1][0] = gcb4[i][0];
        gcin6[count1][1] = gcb4[i][1];
        gcin6[count1][2] = gcb4[i][2];
        ++count1;
        }
        
        if(gcb4[i][4]==2 || gcb4[i][4]==7 || gcb4[i][4]==8)
        {
        gcout6[count2][0] = gcb4[i][0];
        gcout6[count2][1] = gcb4[i][1];
        gcout6[count2][2] = gcb4[i][2];
        ++count2;
        }
    }
    
//  Para Surfaces
	for(i=0; i<gcpara1_count; ++i)
	{ 
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

			gcpara1[i][0]=isurf;
			gcpara1[i][1]=jsurf;
			gcpara1[i][2]=ksurf;
			gcpara1[i][3]=1;
	}

	for(i=0; i<gcpara2_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

			gcpara2[i][0]=isurf;
			gcpara2[i][1]=jsurf;
			gcpara2[i][2]=ksurf;
			gcpara2[i][3]=1;
	}

	for(i=0; i<gcpara3_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

			gcpara3[i][0]=isurf;
			gcpara3[i][1]=jsurf;
			gcpara3[i][2]=ksurf;
			gcpara3[i][3]=1;
	}

	for(i=0; i<gcpara4_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

			gcpara4[i][0]=isurf;
			gcpara4[i][1]=jsurf;
			gcpara4[i][2]=ksurf;
			gcpara4[i][3]=1;
	}

	for(i=0; i<gcpara5_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

			gcpara5[i][0]=isurf;
			gcpara5[i][1]=jsurf;
			gcpara5[i][2]=ksurf;
			gcpara5[i][3]=1;
	}

	for(i=0; i<gcpara6_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

			gcpara6[i][0]=isurf;
			gcpara6[i][1]=jsurf;
			gcpara6[i][2]=ksurf;
			gcpara6[i][3]=1;
	}
	
//  Para Corners
	for(i=0; i<gcparaco1_count; ++i)
	{
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side1=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side2=iin;
        
        grid.read((char*)&iin, sizeof (int));
        paraconum=iin;

			gcparaco1[i][0]=isurf;
			gcparaco1[i][1]=jsurf;
			gcparaco1[i][2]=ksurf;
			gcparaco1[i][3]=side1;
			gcparaco1[i][4]=side2;
			gcparaco1[i][5]=paraconum;
	}

	for(i=0; i<gcparaco2_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side1=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side2=iin;
        
        grid.read((char*)&iin, sizeof (int));
        paraconum=iin;

			gcparaco2[i][0]=isurf;
			gcparaco2[i][1]=jsurf;
			gcparaco2[i][2]=ksurf;
			gcparaco2[i][3]=side1;
			gcparaco2[i][4]=side2;
			gcparaco2[i][5]=paraconum;
	}

	for(i=0; i<gcparaco3_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side1=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side2=iin;
        
        grid.read((char*)&iin, sizeof (int));
        paraconum=iin;

			gcparaco3[i][0]=isurf;
			gcparaco3[i][1]=jsurf;
			gcparaco3[i][2]=ksurf;
			gcparaco3[i][3]=side1;
			gcparaco3[i][4]=side2;
			gcparaco3[i][5]=paraconum;
	}

	for(i=0; i<gcparaco4_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side1=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side2=iin;
        
        grid.read((char*)&iin, sizeof (int));
        paraconum=iin;

			gcparaco4[i][0]=isurf;
			gcparaco4[i][1]=jsurf;
			gcparaco4[i][2]=ksurf;
			gcparaco4[i][3]=side1;
			gcparaco4[i][4]=side2;
			gcparaco4[i][5]=paraconum;
	}

	for(i=0; i<gcparaco5_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side1=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side2=iin;
        
        grid.read((char*)&iin, sizeof (int));
        paraconum=iin;

			gcparaco5[i][0]=isurf;
			gcparaco5[i][1]=jsurf;
			gcparaco5[i][2]=ksurf;
			gcparaco5[i][3]=side1;
			gcparaco5[i][4]=side2;
			gcparaco5[i][5]=paraconum;
	}

	for(i=0; i<gcparaco6_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side1=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side2=iin;
        
        grid.read((char*)&iin, sizeof (int));
        paraconum=iin;

			gcparaco6[i][0]=isurf;
			gcparaco6[i][1]=jsurf;
			gcparaco6[i][2]=ksurf;
			gcparaco6[i][3]=side1;
			gcparaco6[i][4]=side2;
			gcparaco6[i][5]=paraconum;
	}

	for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
    grid.read((char*)&iin, sizeof (int));
    flagslice4[(i-imin)*jmax + (j-jmin)]=iin;
    }
    
    //  Paraslice Surfaces
	for(i=0; i<gcslpara1_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

			gcslpara1[i][0]=isurf;
			gcslpara1[i][1]=jsurf;
	}
    
    for(i=0; i<gcslpara2_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

			gcslpara2[i][0]=isurf;
			gcslpara2[i][1]=jsurf;
	}
    
    for(i=0; i<gcslpara3_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

			gcslpara3[i][0]=isurf;
			gcslpara3[i][1]=jsurf;
	}
    
    for(i=0; i<gcslpara4_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

			gcslpara4[i][0]=isurf;
			gcslpara4[i][1]=jsurf;
	}
    
    //  Paraslice Surfaces
	for(i=0; i<gcslparaco1_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side1=iin;

			gcslparaco1[i][0]=isurf;
			gcslparaco1[i][1]=jsurf;
           gcslparaco1[i][3]=side1;
	}
    
    for(i=0; i<gcslparaco2_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side1=iin;

			gcslparaco2[i][0]=isurf;
			gcslparaco2[i][1]=jsurf;
            gcslparaco2[i][3]=side1;
	}
    
    for(i=0; i<gcslparaco3_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side1=iin;

			gcslparaco3[i][0]=isurf;
			gcslparaco3[i][1]=jsurf;
            gcslparaco3[i][3]=side1;
	}
    
    for(i=0; i<gcslparaco4_count; ++i)
	{
		grid.read((char*)&iin, sizeof (int));
        isurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;
        
        grid.read((char*)&iin, sizeof (int));
        side1=iin;
            
            gcslparaco4[i][0]=isurf;
            gcslparaco4[i][1]=jsurf;
            gcslparaco4[i][3]=side1;
	}
    
    
    // bed
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
    grid.read((char*)&ddn, sizeof (double));
    bed[(i-imin)*jmax + (j-jmin)]=ddn;
    }
    
	if(solidread>0)
	for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
    grid.read((char*)&ddn, sizeof (double));
    solidbed[(i-imin)*jmax + (j-jmin)]=ddn;
    }
    
    if(toporead>0)
	for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
    grid.read((char*)&ddn, sizeof (double));
    topobed[(i-imin)*jmax + (j-jmin)]=ddn;
    }
	
	if(P150>0)
	for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
    grid.read((char*)&ddn, sizeof (double));
    data[(i-imin)*jmax + (j-jmin)]=ddn;
    }
    
    gcin4a_count=gcin_count;
	gcout4a_count=gcout_count;
    

	grid.close();
}

void lexer::clear(char& b, int& j)
{
	b='a';
	j=0;
}
