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

    if(p->P186>0 && p->mpirank==0)
        mkdir("./REEF3D_NHFLOW_VTP_TIME_AVG_FSF",0777);

    if(p->P186>0)
        initialize_avg_arrays(p);

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

        pgc->gcsl_start4Vint(p,wetmax,50);
    }
}

void nhflow_vtp_fsf::initialize_avg_arrays(lexer *p)
{
    if(avg_count == p->P186 && avg_npt == p->pointnum2D && !uavg.empty())
        return;

    avg_count = p->P186;
    avg_npt = p->pointnum2D;
    avg_started.assign(avg_count, false);
    avg_printed.assign(avg_count, false);
    avg_time.assign(avg_count, 0.0);
    avg_tstart.assign(avg_count, 0.0);
    avg_tend.assign(avg_count, 0.0);

    const size_t nacc = static_cast<size_t>(avg_count) * static_cast<size_t>(avg_npt);

    uavg.assign(nacc, 0.0);
    vavg.assign(nacc, 0.0);
    wavg.assign(nacc, 0.0);
    zavg.assign(nacc, 0.0);
    wlavg.assign(nacc, 0.0);
    etaavg.assign(nacc, 0.0);

    for(n=0; n<avg_count; ++n)
    {
        avg_tstart[n] = p->P186_tstart[n];
        avg_tend[n] = p->P186_tend[n];
    }
}

void nhflow_vtp_fsf::start_avg(lexer *p, fdm_nhf *d, ghostcell* pgc, sediment *psed)
{
    initialize_avg_arrays(p);

    if(avg_count<=0 || avg_npt<=0)
        return;

    const double eps = 1.0e-12;

    // Loop over all averaging intervals
    for(int q=0; q<avg_count; ++q)
    {
        if(avg_printed[q])
            continue;

        const double tbegin = avg_tstart[q];
        const double tend = avg_tend[q];
        const double t0 = p->simtime - p->dt_old;
        const double t1 = p->simtime;

        if(tend <= tbegin)
        {
            avg_printed[q] = true;
            continue;
        }

        if(t1 < tbegin - eps)
            continue;

        const double overlap_start = MAX(t0, tbegin);
        const double overlap_end = MIN(t1, tend);
        const double effective_dt = overlap_end - overlap_start;

        if(effective_dt > 0.0)
        {
            if(!avg_started[q])
                avg_started[q] = true;

            const size_t base = static_cast<size_t>(q) * static_cast<size_t>(avg_npt);
            int idx = 0;
            k = p->knoz - 1;
            TPSLICELOOP
            {
                double val_u = 0.0;
                double val_v = 0.0;
                double val_w = 0.0;
                double val_z = 0.0;
                double val_wl = 0.0;
                double val_eta = 0.0;

                if(p->j_dir)
                {
                    val_u = 0.5*(d->U[IJK]+d->U[IJp1K]);
                    val_v = 0.5*(d->V[IJK]+d->V[IJp1K]);
                    val_w = 0.5*(d->W[IJK]+d->W[IJp1K]);
                }
                else
                {
                    jj=j;
                    j=0;
                    val_u = d->U[IJK];
                    val_v = d->V[IJK];
                    val_w = d->W[IJK];
                    j=jj;
                }

                val_z = p->sl_ipol4eta(p->wet,d->eta,d->bed) + p->wd;
                val_wl = p->sl_ipol4(d->WL);
                val_eta = p->sl_ipol4eta_wd(p->wet,d->eta,d->bed);

                const size_t loc = base + static_cast<size_t>(idx);
                uavg[loc] += effective_dt * val_u;
                vavg[loc] += effective_dt * val_v;
                wavg[loc] += effective_dt * val_w;
                zavg[loc] += effective_dt * val_z;
                wlavg[loc] += effective_dt * val_wl;
                etaavg[loc] += effective_dt * val_eta;
                ++idx;
            }

            avg_time[q] += effective_dt;
        }

        if(t1 >= tend - eps && avg_time[q] > 0.0)
        {
            print2D_avg(p,d,pgc,psed,q);
            avg_printed[q] = true;
        }
    }
}

void nhflow_vtp_fsf::start(lexer *p, fdm_nhf *d, ghostcell* pgc, sediment *psed)
{
    if(p->P180==1)
        print2D(p,d,pgc,psed);
}

void nhflow_vtp_fsf::print2D(lexer *p, fdm_nhf *d, ghostcell* pgc, sediment *psed)
{
    SLICELOOP4
    {
        if(d->breaking(i,j)>=1)
            d->breaking_print(i,j)=double(d->breaking(i,j));
        else if(d->breaking(i,j)==0)
            d->breaking_print(i,j)=0.0;
    }

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
        ffn=float(p->Xout(p->XN[IP1],p->YN[JP1]));
        result.write((char*)&ffn, sizeof(float));

        ffn=float(p->Yout(p->XN[IP1],p->YN[JP1]));
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

void nhflow_vtp_fsf::print2D_avg(lexer *p, fdm_nhf *d, ghostcell* pgc, sediment *psed, int q)
{
    if(q<0 || q>=avg_count || avg_npt<=0 || avg_time[q] <= 0.0)
        return;

    int num = q + 1;

    if(p->mpirank==0)
        pvtp_avg(p,num);

    int offset_avg[20];
    n=0;
    offset_avg[n]=0;
    ++n;

    offset_avg[n]=offset_avg[n-1]+sizeof(float)*p->pointnum2D*3+sizeof(int);
    ++n;

    offset_avg[n]=offset_avg[n-1]+sizeof(float)*p->pointnum2D*3+sizeof(int);
    ++n;

    offset_avg[n]=offset_avg[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
    ++n;

    offset_avg[n]=offset_avg[n-1]+sizeof(float)*p->pointnum2D+sizeof(int);
    ++n;

    offset_avg[n]=offset_avg[n-1] + sizeof(int)*p->polygon_sum*3+sizeof(int);
    ++n;
    offset_avg[n]=offset_avg[n-1] + sizeof(int)*p->polygon_sum+sizeof(int);
    ++n;

    snprintf(name,sizeof(name),"./REEF3D_NHFLOW_VTP_TIME_AVG_FSF/REEF3D-NHFLOW-FSF-TIMEAVG-%08i-%06i.vtp",num,p->mpirank+1);
    ofstream result;
    result.open(name, ios::binary);

    vtp3D::beginning(p,result,p->pointnum2D,0,0,0,p->polygon_sum);

    n=0;
    vtp3D::points(result,offset_avg,n);

    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset_avg[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"eta\" format=\"appended\" offset=\""<<offset_avg[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"WL\" format=\"appended\" offset=\""<<offset_avg[n]<<"\"/>\n";
    ++n;
    result<<"</PointData>\n";

    vtp3D::polys(result,offset_avg,n);
    vtp3D::ending(result);

    float ffn;
    int iin;
    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    int idx=0;
    const size_t base = static_cast<size_t>(q) * static_cast<size_t>(avg_npt);
    TPSLICELOOP
    {
        ffn=float(p->Xout(p->XN[IP1],p->YN[JP1]));
        result.write((char*)&ffn, sizeof(float));

        ffn=float(p->Yout(p->XN[IP1],p->YN[JP1]));
        result.write((char*)&ffn, sizeof(float));

        ffn=float(wlavg[base + idx]/avg_time[q]);
        result.write((char*)&ffn, sizeof(float));
        ++idx;
    }

    iin=sizeof(float)*p->pointnum2D*3;
    result.write((char*)&iin, sizeof(int));
    idx=0;
    TPSLICELOOP
    {
        ffn=float(uavg[base + idx]/avg_time[q]);
        result.write((char*)&ffn, sizeof(float));
        ffn=float(vavg[base + idx]/avg_time[q]);
        result.write((char*)&ffn, sizeof(float));
        ffn=float(wavg[base + idx]/avg_time[q]);
        result.write((char*)&ffn, sizeof(float));
        ++idx;
    }

    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    idx=0;
    TPSLICELOOP
    {
        ffn=float(etaavg[base + idx]/avg_time[q]);
        result.write((char*)&ffn, sizeof(float));
        ++idx;
    }

    iin=sizeof(float)*p->pointnum2D;
    result.write((char*)&iin, sizeof(int));
    idx=0;
    TPSLICELOOP
    {
        ffn=float(wlavg[base + idx]/avg_time[q]);
        result.write((char*)&ffn, sizeof(float));
        ++idx;
    }

    iin=sizeof(int)*p->polygon_sum*3;
    result.write((char*)&iin, sizeof(int));
    SLICEBASELOOP
    {
        iin=int(d->nodeval2D(i-1,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i,j))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i-1,j-1))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i,j))-1;
        result.write((char*)&iin, sizeof(int));

        iin=int(d->nodeval2D(i-1,j))-1;
        result.write((char*)&iin, sizeof(int));
    }

    iin=sizeof(int)*p->polygon_sum;
    result.write((char*)&iin, sizeof(int));
    for(n=0;n<p->polygon_sum;++n)
    {
        iin=(n+1)*3;
        result.write((char*)&iin, sizeof(int));
    }

    vtp3D::footer(result);
    result.close();
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
