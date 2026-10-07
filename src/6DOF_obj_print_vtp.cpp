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

#include<fstream>
#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"
#include"runlog.h"
#include"lagoon_store.h"
#include<cctype>
#include<vector>

namespace
{
// the LAGOON body writer of this run (P 18), rank 0: one for all bodies
lagoon_bodies *lagoon_body_writer(lexer *p)
{
    static lagoon_bodies *writer = nullptr;
    if(writer==nullptr)
    {
        const std::string solver = p->A10==6 ? "CFD" : p->A10==3 ? "FNPF" : p->A10==2 ? "SFLOW" : "NHFLOW";
        std::string lower = solver;
        for(char &ch : lower)
            ch = char(std::tolower((unsigned char)ch));
        std::string run;
        if(p->plog)
            run = "{\"type\": \"run\", \"run\": " + lagoon_store::json_string(p->plog->id()) + "}";
        writer = new lagoon_bodies("./REEF3D_" + solver + ".lagoon", solver, lower + "_body",
                                   "REEF3D_" + solver + "_6DOF_VTP", run);
    }
    return writer;
}
}

void sixdof_obj::print_vtp(lexer *p, ghostcell *pgc)
{
    // print normals
    // print_normals_vtp(p,pgc);

    bool printflag=false;
    
    // SFFLOW
    if(p->A10==2)
    {
        if(((p->count%p->P181==0) && p->P182<0.0) 
            || (p->simtime>printtime && p->P182>0.0) 
            || (p->count==0 && p->P185==0))
            printflag=true;

        if(p->P185>0)
        for(int qn=0; qn<p->P185; ++qn)
        if(p->simtime>printtime_wT[qn] 
            && p->simtime>=p->P185_ts[qn] 
            && p->simtime<=(p->P185_te[qn]+0.5*p->P185_dt[qn]))
        {
            printflag=true;

            printtime_wT[qn]+=p->P185_dt[qn];
        }
    }
    
    // NHFLOW
    if((p->A10==5||p->A10==3))
    {
        if(p->P19==1)
        if(((p->count%p->P20==0) && p->P30<0.0) 
            || (p->simtime>printtime && p->P30>0.0) 
            || (p->count==0 && p->P35==0))
            printflag=true;
            
        if(p->P19==2)
        if(((p->count%p->P181==0) && p->P182<0.0) 
            || (p->simtime>printtime && p->P182>0.0) 
            || (p->count==0 && p->P185==0))
            printflag=true;
        

        if(p->P185>0)
        for(int qn=0; qn<p->P185; ++qn)
        if(p->simtime>printtime_wT[qn] 
            && p->simtime>=p->P185_ts[qn] 
            && p->simtime<=(p->P185_te[qn]+0.5*p->P185_dt[qn]))
        {
            printflag=true;

            printtime_wT[qn]+=p->P185_dt[qn];
        }
    }

    // CFD
    if(p->A10==6)
    {
        if(((p->count%p->P20==0) && p->P30<0.0) 
            || (p->simtime>printtime && p->P30>0.0) 
            || (p->count==0 && p->P35==0))
            printflag=true;

        if(p->P35>0)
        for(int qn=0; qn<p->P35; ++qn)
        if(p->simtime>printtime_wT[qn] 
            && p->simtime>=p->P35_ts[qn] 
            && p->simtime<=(p->P35_te[qn]+0.5*p->P35_dt[qn]))
        {
            printflag=true;

            printtime_wT[qn]+=p->P35_dt[qn];
        }
    }

    if(p->mpirank==0 && printflag)
    {
        if(p->A10==6)
        printtime+=p->P30;
        
        if(p->A10==2)
        printtime+=p->P182;
        
        if((p->A10==5||p->A10==3) && p->P19==1)
        printtime+=p->P30;
        
        if((p->A10==5||p->A10==3) && p->P19==2)
        printtime+=p->P182;

        int num=0;
        if(p->P15==1)
            num = p->printcount_sixdof;
        if(p->P15==2)
            num = p->count;
        if(num<0)
            num=0;

        char path[300];
        if(p->A10==2)
            sprintf(path,"./REEF3D_SFLOW_6DOF_VTP/REEF3D-6DOF-%i-%06i.vtp",n6DOF,num);
        else if((p->A10==5||p->A10==3))
            sprintf(path,(p->A10==3?"./REEF3D_FNPF_6DOF_VTP/REEF3D-6DOF-%i-%06i.vtp":"./REEF3D_NHFLOW_6DOF_VTP/REEF3D-6DOF-%i-%06i.vtp"),n6DOF,num);
        else if(p->A10==6)
            sprintf(path,"./REEF3D_CFD_6DOF_VTP/REEF3D-6DOF-%i-%06i.vtp",n6DOF,num);

        // P 18: the body in the LAGOON store, its mesh once and its motion (x = R x0 + c);
        // P 18 1 leaves out the VTP file then (P 18 2 writes it as well)
        bool stored = false;
        if(p->P18>0)
        {
            std::vector<double> x0(9*size_t(tricount)), x(9*size_t(tricount));
            for(n=0;n<tricount;++n)
            for(q=0;q<3;++q)
            {
                const size_t m = 9*size_t(n) + 3*q;
                x0[m] = tri_x0[n][q]; x0[m+1] = tri_y0[n][q]; x0[m+2] = tri_z0[n][q];
                x[m] = tri_x[n][q];   x[m+1] = tri_y[n][q];   x[m+2] = tri_z[n][q];
            }
            const double R[9] = {R_(0,0), R_(0,1), R_(0,2), R_(1,0), R_(1,1), R_(1,2), R_(2,0), R_(2,1), R_(2,2)};
            const double c[3] = {c_(0), c_(1), c_(2)};
            stored = lagoon_body_writer(p)->output(n6DOF, 3*tricount, x0.data(), x.data(), R, c, p->simtime, num);
        }
        if(stored && p->P18==1)
        {
            ++p->printcount_sixdof;
            return;
        }

        ofstream result;
        result.open(path, ios::binary);

        // ---------------------------------------------------
        n=0;
        offset[n]=0;
        ++n;

        offset[n]=offset[n-1]+sizeof(float)*tricount*3*3 + sizeof(int);
        ++n;
        offset[n]=offset[n-1]+sizeof(int)*tricount*3 + sizeof(int);
        ++n;
        offset[n]=offset[n-1]+sizeof(int)*tricount + sizeof(int);
        ++n;
        //---------------------------------------------

        vtp3D::beginning(p, result, tricount*3, 0, 0, 0, tricount);

        n=0;
        vtp3D::points(result, offset, n);

        vtp3D::polys(result, offset, n);

        vtp3D::ending(result);

        //----------------------------------------------------------------------------

        //  XYZ
        iin=sizeof(float)*tricount*3*3;
        result.write((char*)&iin, sizeof(int));
        for(n=0;n<tricount;++n)
        for(q=0;q<3;++q)
        {
            ffn=tri_x[n][q];
            result.write((char*)&ffn, sizeof(float));

            ffn=tri_y[n][q];
            result.write((char*)&ffn, sizeof(float));

            ffn=tri_z[n][q];
            result.write((char*)&ffn, sizeof(float));
        }

        //  Connectivity POLYGON
        int count=0;
        iin=sizeof(int)*tricount*3;
        result.write((char*)&iin, sizeof(int));
        for(n=0;n<tricount;++n)
        for(q=0;q<3;++q)
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
            iin+= 3;
            result.write((char*)&iin, sizeof(int));
        }

        vtp3D::footer(result);

        result.close();
        char sname[32];
        snprintf(sname,sizeof(sname),"body%i",n6DOF);
        if(p->plog)
        p->plog->written(p,num,sname,"floating_body",path,0);

        ++p->printcount_sixdof;
    }
}
