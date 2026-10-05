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

#include"seastate_vtp.h"
#include"fdm_seastate.h"
#include"lexer.h"
#include"ghostcell.h"
#include"runlog.h"
#include<sys/stat.h>
#include<sys/types.h>

const char *seastate_vtp::scalar_name[seastate_vtp::nscalar] = {"Hs","Tm01","Tm-10","Tp","dir","spread","depth","wet"};

seastate_vtp::seastate_vtp(lexer *p, fdm_seastate *e, ghostcell *pgc)
{
    if(p->I40==0)
    p->printtime=0.0;

    p->printcount=0;

    if(p->mpirank==0)
    mkdir("./REEF3D_SEASTATE_VTP",0777);
}

void seastate_vtp::start(lexer *p, fdm_seastate *e, ghostcell *pgc)
{
    // print out based on iteration
    if((p->count%p->P181==0 && p->P182<0.0 && p->P10==1 && p->P181>0) || (p->count==0 && p->P182<0.0))
    print2D(p,e,pgc);

    // print out based on time
    if((p->simtime>p->printtime && p->P182>0.0 && p->P10==1) || (p->count==0 && p->P182>0.0))
    {
    print2D(p,e,pgc);

    p->printtime+=p->P182;
    }
}

float seastate_vtp::node_wet(lexer *p, fdm_seastate *e, slice &f)
{
    // node (i,j) is the corner of the cells (i,j), (i+1,j), (i,j+1), (i+1,j+1)
    double sum=0.0;
    int count=0;

    for(int qi=0; qi<2; ++qi)
    for(int qj=0; qj<2; ++qj)
    if(e->wet(i+qi,j+qj)==1)
    {
    sum += f(i+qi,j+qj);
    ++count;
    }

    return count>0 ? float(sum/double(count)) : 0.0f;
}

void seastate_vtp::write_scalar(lexer *p, fdm_seastate *e, slice &f, ofstream &result)
{
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));

    TPSLICELOOP
    {
    ffn=node_wet(p,e,f);
    result.write((char*)&ffn, sizeof(float));
    }
}

void seastate_vtp::print2D(lexer *p, fdm_seastate *e, ghostcell *pgc)
{
    int num = 0;
    if(p->P15==1)
    num = p->printcount;
    else if(p->P15==2)
    num = p->count;

    // fill the ghost cells for the nodal values
    pgc->gcsl_start4(p,e->Hs,50);
    pgc->gcsl_start4(p,e->Tm01,50);
    pgc->gcsl_start4(p,e->Tm10,50);
    pgc->gcsl_start4(p,e->Tp,50);
    pgc->gcsl_start4(p,e->dir,50);
    pgc->gcsl_start4(p,e->spread,50);

    if(p->mpirank==0)
    pvtp(p,num);

    // offsets
    n=0;
    offset[n]=0;
    ++n;

    // points
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D*3+sizeof(int);
    ++n;

    // scalars
    for(int q=0; q<nscalar; ++q)
    {
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
    ++n;
    }

    // cells
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum*3+sizeof(int);
    ++n;
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum+sizeof(int);
    ++n;

    sprintf(name,"./REEF3D_SEASTATE_VTP/REEF3D-SEASTATE-%08i-%06i.vtp",num,p->mpirank+1);
    ofstream result;
    result.open(name, ios::binary);

    vtp3D::beginning(p,result,p->pointnum2D,0,0,0,p->polygon_sum);

    n=0;
    vtp3D::points(result,offset,n);

    result<<"<PointData>\n";
    for(int q=0; q<nscalar; ++q)
    {
    result<<"<DataArray type=\"Float32\" Name=\""<<scalar_name[q]<<"\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
    result<<"</PointData>\n";

    vtp3D::polys(result,offset,n);

    vtp3D::ending(result);

    //----------------------------------------------------------------------------

    // XYZ: still water level
    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=p->XN[IP1];
    result.write((char*)&ffn, sizeof(float));

    ffn=p->YN[JP1];
    result.write((char*)&ffn, sizeof(float));

    ffn=p->wd;
    result.write((char*)&ffn, sizeof(float));
    }

    write_scalar(p,e,e->Hs,result);
    write_scalar(p,e,e->Tm01,result);
    write_scalar(p,e,e->Tm10,result);
    write_scalar(p,e,e->Tp,result);
    write_scalar(p,e,e->dir,result);
    write_scalar(p,e,e->spread,result);

    // depth: all fluid cells (sl_ipol4)
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=float(p->sl_ipol4(e->depth));
    result.write((char*)&ffn, sizeof(float));
    }

    // wet: fraction of active neighbour cells
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
    ffn=0.25f*float(e->wet(i,j)+e->wet(i+1,j)+e->wet(i,j+1)+e->wet(i+1,j+1));
    result.write((char*)&ffn, sizeof(float));
    }

    // connectivity
    iin=sizeof(int)*p->polygon_sum*3;
    result.write((char*)&iin, sizeof(int));
    SLICEBASELOOP
    {
    // triangle 1
    iin=int(e->nodeval(i-1,j-1))-1;
    result.write((char*)&iin, sizeof(int));

    iin=int(e->nodeval(i,j-1))-1;
    result.write((char*)&iin, sizeof(int));

    iin=int(e->nodeval(i,j))-1;
    result.write((char*)&iin, sizeof(int));

    // triangle 2
    iin=int(e->nodeval(i-1,j-1))-1;
    result.write((char*)&iin, sizeof(int));

    iin=int(e->nodeval(i,j))-1;
    result.write((char*)&iin, sizeof(int));

    iin=int(e->nodeval(i-1,j))-1;
    result.write((char*)&iin, sizeof(int));
    }

    // offset of connectivity
    iin=sizeof(int)*p->polygon_sum;
    result.write((char*)&iin, sizeof(int));
    for(n=0;n<p->polygon_sum;++n)
    {
    iin=(n+1)*3;
    result.write((char*)&iin, sizeof(int));
    }

    vtp3D::footer(result);

    result.close();

    ++p->printcount;
}

void seastate_vtp::pvtp(lexer *p, int num)
{
    sprintf(name,"./REEF3D_SEASTATE_VTP/REEF3D-SEASTATE-%08i.pvtp",num);

    ofstream result;
    result.open(name);

    vtp3D::beginningParallel(p,result);

    vtp3D::pointsParallel(result);

    result<<"<PPointData>\n";
    for(int q=0; q<nscalar; ++q)
    result<<"<PDataArray type=\"Float32\" Name=\""<<scalar_name[q]<<"\"/>\n";
    result<<"</PPointData>\n";

    char pname[200];
    for(n=0; n<p->M10; ++n)
    {
    sprintf(pname,"REEF3D-SEASTATE-%08i-%06i.vtp",num,n+1);
    result<<"<Piece Source=\""<<pname<<"\"/>\n";
    }

    vtp3D::endingParallel(result);

    result.close();

    if(p->plog)
    p->plog->written(p,num,"seastate","wave_parameters",name,p->M10);
}
