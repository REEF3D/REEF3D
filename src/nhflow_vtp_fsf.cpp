/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_vtp_fsf.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"sediment.h"
#include<sys/stat.h>
#include<sys/types.h>

nhflow_vtp_fsf::nhflow_vtp_fsf(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->I40==0)
        p->printtime=0.0;

    printcount=0;

    // Create Folder
    if(p->mpirank==0)
        mkdir("./REEF3D_NHFLOW_VTP_FSF",0777);

    // 3D
    gcval_eta = 55;
    gcval_fifsf = 60;
    // 2D
    if(p->j_dir==0)
    {
        gcval_eta = 155;
        gcval_fifsf = 160;
    }

    if(p->P131==1)
    {
        p->Iarray(wetmax,p->imax*p->jmax);

        SLICELOOP4
            wetmax[IJ] = 0;

        pgc->gcsl_start4Vint(p,wetmax,50);
    }
}

void nhflow_vtp_fsf::start(lexer *p, fdm_nhf *d, ghostcell* pgc, sediment *psed)
{
    print2D(p,d,pgc,psed);
}

void nhflow_vtp_fsf::print2D(lexer *p, fdm_nhf *d, ghostcell* pgc, sediment *psed)
{
    //pgc->gcsl_start4(p,d->eta,gcval_eta);
    //pgc->gcsl_start4(p,d->Fifsf,gcval_fifsf);
    //pgc->gcsl_start4(p,d->test2D,1);

    SLICELOOP4
    {
        if(d->breaking(i,j)>=1)
            d->breaking_print(i,j)=double(d->breaking(i,j));
        else if(d->breaking(i,j)==0)
            d->breaking_print(i,j)=0.0;
    }

    //pgd->gcsl_start4(p,d->breaking_print,50);

    int num = 0;
    if(p->P15==1)
        num = printcount;
    else if(p->P15==2)
        num = p->count;

    if(p->mpirank==0)
        pvtp(p,num);

    // offsets
    n=0;
    offset[n]=0;
    ++n;

    // Points
    offset[n]=offset[n-1]+sizeof(float)*p->pointnum2D*3+sizeof(int);
    ++n;

    // velocity
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)*3+sizeof(int);
    ++n;

    // Fifsf
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
    ++n;

    // WL
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
    ++n;

    // breaking
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
    ++n;

    // coastline
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
    ++n;

    // wetdry
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
    ++n;

    // test
    if(p->P23==1)
    {
        offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
        ++n;
    }

    // fb
    if(p->P28==1)
    {
        offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
        ++n;
    }

    // Hs
    if(p->P110==1)
    {
        offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
        ++n;
    }

    // wetdry_max
    if(p->P131==1)
    {
        offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+sizeof(int);
        ++n;
    }

    // Cells
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum*3+sizeof(int);
    ++n;
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum+sizeof(int);
    ++n;

    // Open File
    sprintf(name,"./REEF3D_NHFLOW_VTP_FSF/REEF3D-NHFLOW-FSF-%08i-%06i.vtp",num,p->mpirank+1);
    ofstream result;
    result.open(name, ios::binary);

    vtp3D::beginning(p,result,p->pointnum2D,0,0,0,p->polygon_sum);

    n=0;
    vtp3D::points(result,offset,n);

    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"eta\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"WL\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"breaking\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"coastline\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"wetdry\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    if(p->P23==1)
    {
        result<<"<DataArray type=\"Float32\" Name=\"test\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
    }
    if(p->P28==1)
    {
        result<<"<DataArray type=\"Float32\" Name=\"fb\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
    }
    if(p->P110==1)
    {
        result<<"<DataArray type=\"Float32\" Name=\"Hs\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
    }
    if(p->P131==1)
    {
        result<<"<DataArray type=\"Float32\" Name=\"wetdry_max\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
    }
    result<<"</PointData>\n";

    vtp3D::polys(result,offset,n);

    vtp3D::ending(result);

    //----------------------------------------------------------------------------

    float ffn;
    int iin;
    //  XYZ
    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=p->XN[IP1];
        result.write((char*)&ffn, sizeof(float));

        ffn=p->YN[JP1];
        result.write((char*)&ffn, sizeof(float));

        ffn=p->sl_ipol4eta(p->wet,d->eta,d->bed)+p->wd;
        result.write((char*)&ffn, sizeof(float));
    }

    //  Velocities
    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        k = p->knoz-1;

        if(p->j_dir==0)
        {
            jj=j;
            j=0;
            ffn=float(d->U[IJK]);
            j=jj;
        }
        else if(p->j_dir==1)
            ffn=float(0.5*(d->U[IJK]+d->U[IJp1K]));

        result.write((char*)&ffn, sizeof(float));

        if(p->j_dir==0)
        {
            jj=j;
            j=0;
            ffn=float(d->V[IJK]);
            j=jj;
        }
        else if(p->j_dir==1)
            ffn=float(0.5*(d->V[IJK]+d->V[IJp1K]));

        result.write((char*)&ffn, sizeof(float));

        if(p->j_dir==0)
        {
            jj=j;
            j=0;
            ffn=float(d->W[IJK]);
            j=jj;
        }
        else if(p->j_dir==1)
            ffn=float(0.5*(d->W[IJK]+d->W[IJp1K]));

        result.write((char*)&ffn, sizeof(float));
    }

    //  Eta
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4eta_wd(p->wet,d->eta,d->bed));
        result.write((char*)&ffn, sizeof(float));
    }

    //  WL
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(d->WL));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Breaking
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(d->breaking_print));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Coastline
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(d->coastline));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Wetdry
    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    TPSLICELOOP
    {
        ffn = 0.25*float((p->wet[IJ]+p->wet[Ip1J]+p->wet[IJp1]+p->wet[Ip1Jp1]));
        result.write((char*)&ffn, sizeof(float));
    }

    //  test
    if(p->P23==1)
    {
        iin=sizeof(float)*p->pointnum2D;
        result.write((char*)&iin, sizeof(int));
        TPSLICELOOP
        {
            ffn=float(p->sl_ipol4(d->test2D));
            result.write((char*)&ffn, sizeof(float));
        }
    }

    //  fb
    if(p->P28==1)
    {
        iin=sizeof(float)*p->pointnum2D;
        result.write((char*)&iin, sizeof(int));
        TPSLICELOOP
        {
            ffn=float(p->sl_ipol4(d->fs));
            result.write((char*)&ffn, sizeof(float));
        }
    }

    //  Hs
    if(p->P110==1)
    {
        iin=sizeof(float)*p->pointnum2D;
        result.write((char*)&iin, sizeof(int));
        TPSLICELOOP
        {
            ffn=float(p->sl_ipol4(d->Hs));
            result.write((char*)&ffn, sizeof(float));
        }
    }

    //  wetdry_max
    if(p->P131==1)
    {
        iin=sizeof(float)*p->pointnum2D;
        result.write((char*)&iin, sizeof(int));
        TPSLICELOOP
        {
            ffn = 0.25*float((wetmax[IJ]+wetmax[Ip1J]+wetmax[IJp1]+wetmax[Ip1Jp1]));
            result.write((char*)&ffn, sizeof(float));
        }
    }

    //  kin and eps
    //pnhfturb->print_2D(p,d,pgc,result,1);

    //  Connectivity
    iin=sizeof(int)*p->polygon_sum*3;
    result.write((char*)&iin, sizeof(int));
    SLICEBASELOOP
    {
        // Triangle 1
        iin=int(d->nodeval2D(i-1,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i,j))-1;
        result.write((char*)&iin, sizeof(int));

        // Triangle 2
        iin=int(d->nodeval2D(i-1,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i,j))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i-1,j))-1;
        result.write((char*)&iin, sizeof(int));
    }


    //  Offset of Connectivity
    iin=sizeof(int)*p->polygon_sum;
    result.write((char*)&iin, sizeof(int));
    for(n=0;n<p->polygon_sum;++n)
    {
        iin=(n+1)*3;
        result.write((char*)&iin, sizeof(int));
    }

    vtp3D::footer(result);

    result.close();

    ++printcount;
}

void nhflow_vtp_fsf::preproc(lexer *p, fdm_nhf *d, ghostcell* pgc)
{
    if(p->P131==1)
    {
        SLICELOOP4
            wetmax[IJ] = MAX(wetmax[IJ],p->wet[IJ]);

        pgc->gcsl_start4Vint(p,wetmax,50);
    }

}
