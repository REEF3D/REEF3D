/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include<fstream>
#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"

void sixdof_obj::print_normals_vtp(lexer *p, ghostcell *pgc)
{
    double x0,x1,x2,y0,y1,y2,z0,z1,z2;
    double xc,yc,zc;
    double at,bt,ct,st;
    double nx,ny,nz,norm;
    double n0,n1,n2;

    const double factor = 3.6;

    if(p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0) || (p->simtime>printtime && p->P30>0.0) || p->count==0))
    {
        printtimenormal+=p->P30;

        int num = 0;
        if(p->P15==1)
            num = printnormal_count;
        else if(p->P15==2)
            num = p->count;

        char path[300];
        if(p->A10==5)
            sprintf(path,"./REEF3D_NHFLOW_6DOF_Normals_VTP/REEF3D-6DOF-Normals-%i-%06i.vtp",n6DOF,num);
        else if(p->A10==6)
            sprintf(path,"./REEF3D_CFD_6DOF_Normals_VTP/REEF3D-6DOF-Normals-%i-%06i.vtp",n6DOF,num);

        ofstream result;
        result.open(path, ios::binary);

        // ---------------------------------------------------
        n=0;
        offset[n]=0;
        ++n;

        offset[n]=offset[n-1]+sizeof(float)*tricount*3*2 + sizeof(int);
        ++n;
        offset[n]=offset[n-1]+sizeof(int)*tricount*2 + sizeof(int);
        ++n;
        //---------------------------------------------

        vtp3D::beginning(p, result, tricount*2, 0, 0, 0, tricount);

        n=0;
        vtp3D::points(result, offset, n);

        vtp3D::polys(result, offset, n);

        vtp3D::ending(result);

        //----------------------------------------------------------------------------

        //  XYZ
        iin=sizeof(float)*tricount*3*2;
        result.write((char*)&iin, sizeof(int));
        for(n=0;n<tricount;++n)
        {
            x0 = tri_x[n][0];
            y0 = tri_y[n][0];
            z0 = tri_z[n][0];

            x1 = tri_x[n][1];
            y1 = tri_y[n][1];
            z1 = tri_z[n][1];

            x2 = tri_x[n][2];
            y2 = tri_y[n][2];
            z2 = tri_z[n][2];

            nx = (y1 - y0)*(z2 - z0) - (y2 - y0)*(z1 - z0);
            ny = (x2 - x0)*(z1 - z0) - (x1 - x0)*(z2 - z0);
            nz = (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0);

            norm = sqrt(nx*nx + ny*ny + nz*nz);

            nx /= norm > 1.0e-20 ? norm : 1.0e20;
            ny /= norm > 1.0e-20 ? norm : 1.0e20;
            nz /= norm > 1.0e-20 ? norm : 1.0e20;

            // Center of triangle
            xc = (x0 + x1 + x2)/3.0;
            yc = (y0 + y1 + y2)/3.0;
            zc = (z0 + z1 + z2)/3.0;

            ffn=xc;
            result.write((char*)&ffn, sizeof(float));

            ffn=yc;
            result.write((char*)&ffn, sizeof(float));

            ffn=zc;
            result.write((char*)&ffn, sizeof(float));

            if((p->j_dir==1 || fabs(ny)<0.001))
            {
                xc += nx*p->DXM*factor;
                yc += ny*p->DXM*factor;
                zc += nz*p->DXM*factor;
            }

            ffn=xc;
            result.write((char*)&ffn, sizeof(float));

            ffn=yc;
            result.write((char*)&ffn, sizeof(float));

            ffn=zc;
            result.write((char*)&ffn, sizeof(float));
        }

        //  Connectivity LINE
        int count=0;
        iin=sizeof(int)*tricount*2;
        result.write((char*)&iin, sizeof(int));
        for(n=0;n<tricount;++n)
        for(q=0;q<2;++q)
        {
            iin=count;
            result.write((char*)&iin, sizeof(int));
            ++count;
        }

        //  Offset of Connectivity
        iin=sizeof(int)*tricount;
        result.write((char*)&iin, sizeof(int));
        iin=0;
        for(n=0;n<tricount;++n)
        {
            iin+= 2;
            result.write((char*)&iin, sizeof(int));
        }

        vtp3D::footer(result);

        result.close();

        ++printnormal_count;
    }
}
