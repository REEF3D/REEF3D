/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"fnpf_vtu3D.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fnpf_print_wsf.h"
#include"fnpf_print_wsfline.h"
#include"fnpf_print_wsfline_y.h"
#include"potentialfile_out.h"
#include<sys/stat.h>
#include<sys/types.h>

fnpf_vtu3D::fnpf_vtu3D(lexer* p, fdm_fnpf *c, ghostcell *pgc) 
{
    if(p->I40==0)
    {
	p->printtime=0.0;
	p->sedprinttime=0.0;
	p->fsfprinttime=0.0;
	p->probeprinttime=0.0;
	p->stateprinttime=0.0;
    p->exportprinttime=0.0;
    }
	
	p->Darray(printtime_wT,p->P35);
	
	for(int qn=0; qn<p->P35; ++qn)
	printtime_wT[qn]=p->P35_ts[qn]; 
	

	printcount=0;

	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_FNPF_VTU",0777);
    
    
    pwsf=new fnpf_print_wsf(p,c);
    
    pwsfline=new fnpf_print_wsfline(p,c,pgc);
    
    pwsfline_y=new fnpf_print_wsfline_y(p,c,pgc);
    
    if(p->P230>0)
    ppotentialfile = new potentialfile_out(p,c,pgc);
}

fnpf_vtu3D::~fnpf_vtu3D()
{
}

void fnpf_vtu3D::start(lexer* p, fdm_fnpf* c,ghostcell* pgc, ioflow *pflow)
{
    // Gages
	if(p->P51>0)
	pwsf->height_gauge(p,c,pgc,c->eta);
  
		// Print out based on iteration
        if(p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P10==1 && p->P20>0)
		{
        print_vtu(p,c,pgc);
		}
		
		// Print out based on time
        if((p->simtime>p->printtime && p->P30>0.0 && p->P34<0.0 && p->P10==1) || (p->count==0 &&  p->P30>0.0))
        {
        print_vtu(p,c,pgc);
		
        p->printtime+=p->P30;
        }
				
		// Print out based on time interval
		if(p->P10==1 && p->P35>0)
		for(int qn=0; qn<p->P35; ++qn)
		if(p->simtime>printtime_wT[qn] && p->simtime>=p->P35_ts[qn] && p->simtime<=(p->P35_te[qn]+0.5*p->P35_dt[qn]))
		{
		print_vtu(p,c,pgc);	
			
		printtime_wT[qn]+=p->P35_dt[qn];
		}
        
    // Gages
    if((p->P52>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P52>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline->start(p,c,pgc,pflow,c->eta);
    
    if((p->P56>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P56>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline_y->start(p,c,pgc,pflow,c->eta);
}

void fnpf_vtu3D::print_vtu(lexer* p, fdm_fnpf *c, ghostcell* pgc)
{
     //
    pgc->start4(p,c->u,110);
    pgc->start4(p,c->v,111);
	pgc->start4(p,c->w,112);
    pgc->gcsl_start4(p,c->WL,50);
    pgc->gcsl_start4(p,c->bed,50);
    pgc->start4(p,c->test,50);
	
	pgc->dgcpol(p,c->u,p->dgc1,p->dgc1_count,11);
	pgc->dgcpol(p,c->v,p->dgc2,p->dgc2_count,12);
	pgc->dgcpol(p,c->w,p->dgc3,p->dgc3_count,13);
    pgc->dgcpol(p,c->test,p->dgc4,p->dgc4_count,14);
    pgc->dgcpol(p,c->Fi4,p->dgc4,p->dgc4_count,14);
    pgc->dgcslpol(p,c->WL,p->dgcsl4,p->dgcsl4_count,14);
    pgc->dgcslpol(p,c->bed,p->dgcsl4,p->dgcsl4_count,14);
	
	c->u.ggcpol(p);
	c->v.ggcpol(p);
	c->w.ggcpol(p);
    c->Fi4.ggcpol(p);
    c->WL.ggcpol(p);
    c->test.ggcpol(p);
    
    i=-1;
    j=-1;
    if(i+p->origin_i==-1 && j+p->origin_j==-1 )
    c->WL(i,j) = c->WL(i+1,j+1);
    
    
    //----------
    
    if(p->mpirank==0)
    pvtu(p,pgc);

    name_iter(p,pgc);
	
	// Open File
	ofstream result;
	result.open(name, ios::binary);

    n=0;

	offset[n]=0;
	++n;
	
	// velocity
	offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)*3+4;
	++n;
	
	// scalars
	
		// pressure
	offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;
    
    // Fi
	offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;

    // test
    if(p->P23==1)
	{
	offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;
	}
    // elevation
	offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;
    
    if(p->P25==1)
	{
		// solid
	offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)+4;
	++n;
	}
	
	// Points
    offset[n]=offset[n-1]+4*(p->pointnum+p->ccptnum)*3+4;
    ++n;
	
	// Cells
    offset[n]=offset[n-1] + 4*p->tpcellnum*8 + 4*p->ccedgenum + 4;
    ++n;
    offset[n]=offset[n-1] + 4*(p->tpcellnum+p->ccellnum)+4;
    ++n;
	offset[n]=offset[n-1] + 4*(p->tpcellnum+p->ccellnum)+4;
    ++n;
	//---------------------------------------------

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<UnstructuredGrid>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<p->pointnum+p->ccptnum<<"\" NumberOfCells=\""<<p->tpcellnum+p->ccellnum<<"\">"<<endl;

    n=0;
    result<<"<PointData >"<<endl;
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	
	
    result<<"<DataArray type=\"Float32\" Name=\"pressure\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    
    if(p->A10==3)
	{
    result<<"<DataArray type=\"Float32\" Name=\"Fi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}
    
    if(p->P23==1)
	{
    result<<"<DataArray type=\"Float32\" Name=\"test\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}
	
    result<<"<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    
    if(p->P25==1)
	{
	result<<"<DataArray type=\"Float32\" Name=\"solid\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}
	result<<"</PointData>"<<endl;


    result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;

    result<<"<Cells>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"</Cells>"<<endl;

    result<<"</Piece>"<<endl;
    result<<"</UnstructuredGrid>"<<endl;

//----------------------------------------------------------------------------
    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";
	
	
//  Velocities
    iin=3*4*(p->pointnum+p->ccptnum);
	result.write((char*)&iin, sizeof (int));
    TPLOOP
	{
	ffn=float(p->ipol4(c->u));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ipol4(c->v));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ipol4(c->w));
	result.write((char*)&ffn, sizeof (float));
	}
    

	for(n=0;n<p->ccptnum;++n)
	{
	ffn=float(p->ccipol4(c->u,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2]));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ccipol4(c->v,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2]));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ccipol4(c->w,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2]));
	result.write((char*)&ffn, sizeof (float));

	}
	
//  Pressure
	iin=4*(p->pointnum+p->ccptnum);
	result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4press(c->press));
	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;++n)
	{
	ffn=float(0.0);//float(p->ccipol4press(a,c->press,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2]));
	result.write((char*)&ffn, sizeof (float));
	}
    
//  Fi
    iin=4*(p->pointnum+p->ccptnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
    ffn=float(c->Fi[FIJK]);//float(p->ipol4press(c->Fi4));
    
    if(k==-1 && j==-1)
	ffn=float(c->Fi[FIJp1Kp1]);//float(p->ipol4press(c->Fi4));
	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;++n)
	{
	ffn=float(0.0);//float(p->ccipol4(c->Fi4,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2]));
	result.write((char*)&ffn, sizeof (float));
	}
    
    if(p->P23==1)
	{
//  test
    iin=4*(p->pointnum+p->ccptnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4_a(c->test));
	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;++n)
	{
	ffn=float(p->ccipol4_a(c->test,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2]));
	result.write((char*)&ffn, sizeof (float));
	}
	}
	
//  elevation
	iin=4*(p->pointnum+p->ccptnum)*3;
	result.write((char*)&iin, sizeof (int));
    TPLOOP
	{
	ffn=float(p->ZN[KP1]*c->WL(i,j) + c->bed(i,j));
	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;++n)
	{
	ffn=float(p->ccpoint[n][2]);
	result.write((char*)&ffn, sizeof (float));
	}
    
    	
	if(p->P25==1)
	{
//  solid
    iin=4*(p->pointnum+p->ccptnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(1.0);//float(p->ipol4_a(c->solid));
	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;++n)
	{
	ffn=float(1.0);//float(p->ccipol4_a(c->solid,p->ccpoint[n][0],p->ccpoint[n][1],p->ccpoint[n][2]));
	result.write((char*)&ffn, sizeof (float));
	}
	}

//  XYZ
	double theta_y = p->B192_1*(PI/180.0);
	double omega_y = 2.0*PI*p->B192_2;

    if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    phase = omega_y*p->simtime;
    
	iin=4*(p->pointnum+p->ccptnum)*3;
	result.write((char*)&iin, sizeof (int));
    TPLOOP
	{
    zcoor = p->ZN[KP1]*c->WL(i,j) + c->bed(i,j); 
    
    if(c->wet(i,j)==0)
    zcoor=c->bed(i,j);
    
    if(i+p->origin_i==-1 && j+p->origin_j==-1 && c->wet(0,0)==1)
    zcoor = p->ZN[KP1]*c->WL(i,j) + c->bed(i,j); 
    
    ffn=float( (p->XN[IP1]-p->B192_3)*cos(theta_y*sin(phase)) - (zcoor-p->B192_4)*sin(theta_y*sin(phase)) + p->B192_3);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->YN[JP1]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float((p->XN[IP1]-p->B192_3)*sin(theta_y*sin(phase)) + (zcoor-p->B192_4)*cos(theta_y*sin(phase)) + p->B192_4);
	result.write((char*)&ffn, sizeof (float));
	}

	for(n=0;n<p->ccptnum;++n)
	{
    ffn=float(p->ccpoint[n][0]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ccpoint[n][1]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ccpoint[n][2]);
	result.write((char*)&ffn, sizeof (float));
	}

//  Connectivity
    iin=4*(p->tpcellnum)*8 + 4*p->ccedgenum;
    result.write((char*)&iin, sizeof (int));
    BASELOOP
    if(p->flag5[IJK]!=-20 && p->flag5[IJK]!=-30)
	{
	iin=int(c->nodeval(i-1,j-1,k-1)-1);
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i,j-1,k-1))-1;
	result.write((char*)&iin, sizeof (int));

    iin= int(c->nodeval(i,j,k-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i-1,j,k-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i-1,j-1,k))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i,j-1,k))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i,j,k))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i-1,j,k))-1;
	result.write((char*)&iin, sizeof (int));
	}


	for(n=0;n<p->ccellnum;++n)
	{
    iin=abs(c->pvccnode[n][0]);
	result.write((char*)&iin, sizeof (int));

	iin=abs(c->pvccnode[n][1]);
	result.write((char*)&iin, sizeof (int));

    iin=abs(c->pvccnode[n][2]);
	result.write((char*)&iin, sizeof (int));

	iin=abs(c->pvccnode[n][3]);
	result.write((char*)&iin, sizeof (int));

    if(c->ccedge[n]>4)
    {
	iin=abs(c->pvccnode[n][4]);
	result.write((char*)&iin, sizeof (int));

        if(c->ccedge[n]>5)
        {
        iin=abs(c->pvccnode[n][5]);
        result.write((char*)&iin, sizeof (int));

            if(c->ccedge[n]>6)
            {
            iin=abs(c->pvccnode[n][6]);
            result.write((char*)&iin, sizeof (int));

                if(c->ccedge[n]>7)
                {
                iin=abs(c->pvccnode[n][7]);
                result.write((char*)&iin, sizeof (int));
                }
            }
        }
    }

	}

//  Offset of Connectivity
    iin=4*(p->tpcellnum+p->ccellnum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->tpcellnum;++n)
	{
	iin=(n+1)*8;
	result.write((char*)&iin, sizeof (int));
	}

	for(n=0;n<p->ccellnum;++n)
	{
	iin+=c->ccedge[n];
	result.write((char*)&iin, sizeof (int));
	}

//  Cell types
    iin=4*(p->tpcellnum+p->ccellnum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->tpcellnum;++n)
	{
	iin=12;
	result.write((char*)&iin, sizeof (int));
	}

	for(n=0;n<p->ccellnum;++n)
	{
    if(c->ccedge[n]==4)
	iin=10;
	
	if(c->ccedge[n]==5)
	iin=14;

    if(c->ccedge[n]==6)
	iin=13;
	
	if(c->ccedge[n]==8)
	iin=12;

	result.write((char*)&iin, sizeof (int));
	}

	result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();
	
	++printcount;
	
	
	
	pgc->start1(p,c->u,114);
    pgc->start2(p,c->v,115);
	pgc->start3(p,c->w,116);

	pgc->dgcpol(p,c->u,p->dgc1,p->dgc1_count,11);
	pgc->dgcpol(p,c->v,p->dgc2,p->dgc2_count,12);
	pgc->dgcpol(p,c->w,p->dgc3,p->dgc3_count,13);
	
	c->u.ggcpol(p);
	c->v.ggcpol(p);
	c->w.ggcpol(p);
	
}
