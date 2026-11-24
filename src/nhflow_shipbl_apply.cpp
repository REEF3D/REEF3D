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
Author: Ahmet Soydan
--------------------------------------------------------------------*/

#include"nhflow_shipbl.h"
#include"6DOF_obj_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"increment.h"
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<iostream>

// The hull triangles read from floating.stl are replaced by the displacement body of the hull (nhflow_shipbl): the
// solver then runs with the body as it does with a displacement body written by Cases/KCS/shipbl.py.
// Mass and inertia come from the control file (X 22 - X 24), as with the body file.
//
//   NHF_BL=1          the body from the double-body flow alone, from the start (single run, no first pass)
//   NHF_BL=2          the run starts with the hull as it is; the mean velocity of the run between NHF_BL_T0 and NHF_BL_T1
//                     at the distance NHF_BL_DN from the wall gives the edge velocity (shipbl.py edge), and at NHF_BL_T1
//                     the hull is replaced by its displacement body: the two runs of shipbl.py "all" in one
//   NHF_BL=3          as 2, the run starts with the body of NHF_BL=1
//   NHF_BL_L          length between perpendiculars [m]                                   (needed)
//   NHF_BL_XFP        x of the forward perpendicular in the tank, hull at even keel [m]   (needed)
//   NHF_BL_U          speed [m/s]                         default NHF_SRC
//   NHF_BL_T, _B      draft, beam [m]                     default from the triangles
//   NHF_BL_ZW         still water level [m]               default F 60
//   NHF_BL_SINK       sinkage of the triangles [m, + down], NHF_BL_TRIM trim [deg, + bow up] about midship at the
//                     still water line: the boundary layer is computed for the hull at even keel (as shipbl.py does)
//                     and the body is given the attitude of the hull
//   NHF_BL_XTR        trip at s = x/L                     default 0.05
//   NHF_BL_FSW=0      without the surface wake            default 1
//   NHF_BL_T0, _T1    averaging interval of the run [s]   (needed with NHF_BL=2, 3)
//   NHF_BL_NUP        number of updates of the body       default 1; the next ones after every further interval of the same
//                     length, each from the mean flow of its interval (the flow around the body of the update before:
//                     the stern flow settles as slowly as the waves, 3 L/U; early updates take most of the change)
//   NHF_BL_DN         distance of the sample points from the wall [m]     default NHF_SRC_PSI + 0.01 L
//   NHF_BL_DUMP=1     writes REEF3D_NHFLOW_6DOF/NHF_BL_dispbody.stl (NHF_BL=2, 3: NHF_BL_dispbody2.stl and the mean
//                     velocity at the sample points, NHF_BL_edge.bin, at NHF_BL_T1)
//
// Rank 0 computes (seconds; the dense matrix of the double body takes 8 N^2 bytes, N about 8000), the triangles
// are broadcast. Not written to the checkpoints: a restart begins again with the hull of floating.stl.

namespace
{
    struct blstate
    {
        int mode = 0;
        nhflow_shipbl *bl = nullptr;        // rank 0
        double U=0.0, xfp=0.0, yc=0.0, zw=0.0, xm=0.0, sink=0.0, ct=1.0, st=0.0;
        double t0=0.0, t1=0.0, w=0.0, dn=0.0, tsum=0.0;
        int ts=0, entity=0, M=0;
        bool done = false, ready = false;
        double tlast = 0.0;
        int nup=1, iup=0;
        std::vector<double> ph;             // sample points, hull frame at even keel (rank 0)
        std::vector<double> pt;             // sample points in the tank, M x (x, y - yc, z)
        std::vector<double> acc;            // M x (u v w weight) port, then starboard: time integrals of this rank
        std::vector<double> last;           // the same of the last time step (NHF_BL_DUMP)
    };
    blstate bs;

    double envd(const char *name, double def)
    {
        return getenv(name)!=nullptr ? atof(getenv(name)) : def;
    }

    // hull frame at even keel -> tank with the attitude of the hull (gen_kcs.py, shipbl.py to_tank)
    inline void to_tank(double hx, double hz, double &x, double &z)
    {
        const double xr = bs.xfp + hx - bs.xm;
        x = bs.xm + xr*bs.ct + hz*bs.st;
        z = bs.zw - xr*bs.st + hz*bs.ct - bs.sink;
    }

    void write_stl(const char *name, const std::vector<double> &buf)
    {
        FILE *fp = fopen(name,"w");
        if(fp==nullptr)
        return;
        fprintf(fp,"solid body\n");
        for(size_t n=0; n<buf.size()/9; ++n)
        {
            fprintf(fp," facet normal 0 0 0\n  outer loop\n");
            for(int q=0; q<3; ++q)
            fprintf(fp,"   vertex %.17g %.17g %.17g\n",buf[9*n+3*q],buf[9*n+3*q+1],buf[9*n+3*q+2]);      // exact: a run with this file equals the run that wrote it
            fprintf(fp,"  endloop\n endfacet\n");
        }
        fprintf(fp,"endsolid body\n");
        fclose(fp);
    }

    // body of the class (tank frame, even keel) -> triangles with the attitude of the hull
    void attitude(const std::vector<nhflow_shipbl::tri3> &body, std::vector<double> &buf)
    {
        buf.resize(9*body.size());
        for(size_t n=0; n<body.size(); ++n)
        for(int q=0; q<3; ++q)
        {
            to_tank(body[n][3*q]-bs.xfp,body[n][3*q+2]-bs.zw,buf[9*n+3*q],buf[9*n+3*q+2]);
            buf[9*n+3*q+1] = body[n][3*q+1];
        }
    }
}

void nhflow_shipbl_apply(lexer *p, ghostcell *pgc, double **&tri_x, double **&tri_y, double **&tri_z,
                         double **&tri_x0, double **&tri_y0, double **&tri_z0, int ts, int &tricount)
{
    const int mode = int(envd("NHF_BL",0.0));

    if(mode==0 || bs.mode!=0)
    return;

    const double L = envd("NHF_BL_L",0.0);
    const double xfp = envd("NHF_BL_XFP",-1.0e20);
    const double U = envd("NHF_BL_U",envd("NHF_SRC",0.0));

    if(L<=0.0 || xfp<-1.0e19 || U<=0.0)
    {
        if(p->mpirank==0)
        std::cout<<"NHF_BL: NHF_BL_L, NHF_BL_XFP and the speed (NHF_BL_U or NHF_SRC) are needed"<<std::endl;
        pgc->final(true);
    }

    bs.mode = mode;
    bs.U = U;
    bs.xfp = xfp;
    bs.zw = envd("NHF_BL_ZW",p->F60);
    bs.sink = envd("NHF_BL_SINK",0.0);
    const double trim = envd("NHF_BL_TRIM",0.0)*3.14159265358979323846/180.0;
    bs.xm = xfp + 0.5*L;
    bs.ct = cos(trim);
    bs.st = sin(trim);
    bs.ts = ts;
    bs.t0 = envd("NHF_BL_T0",-1.0);
    bs.t1 = envd("NHF_BL_T1",-1.0);
    bs.nup = std::max(int(envd("NHF_BL_NUP",1.0)),1);
    bs.w = bs.t1-bs.t0;
    bs.dn = envd("NHF_BL_DN",envd("NHF_SRC_PSI",-1.0e20) + 0.01*L);

    if(mode>=2 && (bs.t1<=bs.t0 || bs.t0<0.0 || bs.dn<=0.0))
    {
        if(p->mpirank==0)
        std::cout<<"NHF_BL="<<mode<<": NHF_BL_T0 < NHF_BL_T1 and NHF_BL_DN (or NHF_SRC_PSI) are needed"<<std::endl;
        pgc->final(true);
    }

    const double zw = bs.zw, sink = bs.sink, xm = bs.xm, ct = bs.ct, st = bs.st;

    int nnew=0;
    std::vector<double> buf;

    if(p->mpirank==0)
    {
        // hull to even keel:  x' = xm + xr cos + zr sin,  z' = zw - xr sin + zr cos - sink
        std::vector<nhflow_shipbl::tri3> hull(tricount-ts), body;
        double ylo=1.0e30, yhi=-1.0e30;
        for(int n=ts; n<tricount; ++n)
        for(int q=0; q<3; ++q)
        {
            const double X = tri_x[n][q]-xm, Z = tri_z[n][q]-zw+sink;
            hull[n-ts][3*q]   = xm + X*ct - Z*st;
            hull[n-ts][3*q+1] = tri_y[n][q];
            hull[n-ts][3*q+2] = zw + X*st + Z*ct;
            ylo = std::min(ylo,tri_y[n][q]);
            yhi = std::max(yhi,tri_y[n][q]);
        }

        bs.bl = new nhflow_shipbl;
        nhflow_shipbl &bl = *bs.bl;
        bl.L = L;
        bl.T = envd("NHF_BL_T",0.0);
        bl.B = envd("NHF_BL_B",0.0);
        bl.U = U;
        bl.nu = p->W2;
        bl.x_fp = xfp;
        bl.y_c = envd("NHF_BL_YC",0.5*(ylo+yhi));
        bl.z_w = zw;
        bl.xtr = envd("NHF_BL_XTR",0.05);
        bl.fsw = int(envd("NHF_BL_FSW",1.0))!=0;
        bs.yc = bl.y_c;

        bl.set_hull(hull);
        bl.build_grid();
        bl.double_body();

        if(mode>=2)
        bl.offset_points(bs.dn,bs.ph);

        if(mode!=2)
        {
            bl.solve();
            bl.body(body);
            attitude(body,buf);
            nnew = int(body.size());
        }

        std::cout<<"NHF_BL="<<mode<<": "<<(mode!=2 ? "the hull is replaced by its displacement body (double-body flow)" : "the run starts with the hull as it is")
                 <<", "<<tricount-ts<<" hull triangles, U "<<U<<" m/s, nu "<<p->W2<<", sinkage "<<sink<<" m, trim "<<envd("NHF_BL_TRIM",0.0)<<" deg";
        if(mode>=2)
        std::cout<<"; displacement body from the mean flow of t = "<<bs.t0<<" .. "<<bs.t1<<" s at "<<bs.dn<<" m from the wall, "<<bs.nup<<" update(s)";
        std::cout<<std::endl<<bl.log<<std::endl;
        bl.log.clear();

        if(mode!=2 && int(envd("NHF_BL_DUMP",0.0))!=0)
        write_stl("./REEF3D_NHFLOW_6DOF/NHF_BL_dispbody.stl",buf);

        if(mode==1)
        {
            delete bs.bl;
            bs.bl = nullptr;
        }
    }

    // sample points of the run
    if(mode>=2)
    {
        bs.M = int(bs.ph.size()/3);
        pgc->bcast_int(&bs.M,1);
        pgc->bcast_double(&bs.yc,1);
        bs.pt.resize(3*bs.M);
        if(p->mpirank==0)
        for(int m=0; m<bs.M; ++m)
        {
            to_tank(bs.ph[3*m],bs.ph[3*m+2],bs.pt[3*m],bs.pt[3*m+2]);
            bs.pt[3*m+1] = bs.ph[3*m+1];
        }
        pgc->bcast_double(bs.pt.data(),3*bs.M);
        bs.acc.assign(8*bs.M,0.0);
        bs.last.assign(8*bs.M,0.0);
    }

    if(mode==2)
    return;

    pgc->bcast_int(&nnew,1);
    if(p->mpirank!=0)
    buf.resize(9*nnew);
    pgc->bcast_double(buf.data(),9*nnew);

    p->Dresize(tri_x,tricount,ts+nnew,3,3);
    p->Dresize(tri_y,tricount,ts+nnew,3,3);
    p->Dresize(tri_z,tricount,ts+nnew,3,3);
    p->Dresize(tri_x0,tricount,ts+nnew,3,3);
    p->Dresize(tri_y0,tricount,ts+nnew,3,3);
    p->Dresize(tri_z0,tricount,ts+nnew,3,3);

    for(int n=0; n<nnew; ++n)
    for(int q=0; q<3; ++q)
    {
        tri_x[ts+n][q] = buf[9*n+3*q];
        tri_y[ts+n][q] = buf[9*n+3*q+1];
        tri_z[ts+n][q] = buf[9*n+3*q+2];
    }
    tricount = ts+nnew;
}

// NHF_BL=2, 3: called at the end of every time step (driver_loop_nhflow.cpp), p->simtime the time of the fields, p->dt the step
// done. The velocity of the cells is taken to the sample points as Cases/WIGLEY/edge.py does with the mean fields (linear in
// x and y between the cell centres, linear in z between the cell centres of the column, the top and bottom values beyond
// them), port and starboard. NHF_BL_DUMP=2: the cell centres instead of the velocity (test: the samples are the points).
void nhflow_shipbl_sample(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(bs.mode<2 || bs.done || bs.ready || p->simtime<bs.t0)
    return;

    const int M = bs.M;
    const int kn = p->knoz;
    const int marge = increment::marge;
    static const bool coords = int(envd("NHF_BL_DUMP",0.0))==2;

    auto idx = [&](int ii, int jj, int kk)
    {
        return (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;
    };

    auto sample = [&](double x, double y, double z, double *out)
    {
        if(!(x>=p->originx && x<p->endx))
        return false;
        if(p->j_dir==1 && !(y>=p->originy && y<p->endy))
        return false;

        int i0 = p->posc_i(x);
        i0 = std::max(std::min(i0,p->knox-1),0);
        if(x<p->XP[i0+marge])
        --i0;
        const double wx = (x-p->XP[i0+marge])/(p->XP[i0+1+marge]-p->XP[i0+marge]);

        int j0=0;
        double wy=0.0;
        if(p->j_dir==1)
        {
            j0 = p->posc_j(y);
            j0 = std::max(std::min(j0,p->knoy-1),0);
            if(y<p->YP[j0+marge])
            --j0;
            wy = (y-p->YP[j0+marge])/(p->YP[j0+1+marge]-p->YP[j0+marge]);
        }

        out[0]=out[1]=out[2]=0.0;
        for(int a=0; a<2; ++a)
        for(int b=0; b<(p->j_dir==1 ? 2 : 1); ++b)
        {
            const int ii=i0+a, jj=j0+b;
            const double wgt = (a==0 ? 1.0-wx : wx)*(b==0 ? 1.0-wy : wy);
            const double zlo = p->ZSP[idx(ii,jj,0)], zhi = p->ZSP[idx(ii,jj,kn-1)];
            if(!(zhi>zlo))
            return false;
            const double zz = std::min(std::max(z,zlo),zhi);
            int k0=0;
            while(k0<kn-2 && p->ZSP[idx(ii,jj,k0+1)]<zz)
            ++k0;
            const int n0 = idx(ii,jj,k0), n1 = idx(ii,jj,k0+1);
            const double wz = (zz-p->ZSP[n0])/(p->ZSP[n1]-p->ZSP[n0]);
            if(coords)
            {
                out[0] += wgt*p->XP[ii+marge];
                out[1] += wgt*p->YP[jj+marge];
                out[2] += wgt*((1.0-wz)*p->ZSP[n0] + wz*p->ZSP[n1]);
                continue;
            }
            out[0] += wgt*((1.0-wz)*d->U[n0] + wz*d->U[n1]);
            out[1] += wgt*((1.0-wz)*d->V[n0] + wz*d->V[n1]);
            out[2] += wgt*((1.0-wz)*d->W[n0] + wz*d->W[n1]);
        }
        return true;
    };

    for(int m=0; m<M; ++m)
    for(int side=0; side<2; ++side)
    {
        double v[3];
        const int n = 4*(side*M + m);
        bs.last[n]=bs.last[n+1]=bs.last[n+2]=bs.last[n+3]=0.0;
        if(!sample(bs.pt[3*m],bs.yc + (side==0 ? 1.0 : -1.0)*bs.pt[3*m+1],bs.pt[3*m+2],v))
        continue;
        for(int q=0; q<3; ++q)
        {
            bs.acc[n+q] += p->dt*v[q];
            bs.last[n+q] = v[q];
        }
        bs.acc[n+3] += p->dt;
        bs.last[n+3] = 1.0;
    }
    bs.tsum += p->dt;
    bs.tlast = p->simtime;

    if(p->simtime>=bs.t1)
    bs.ready = true;
}

// NHF_BL=2, 3: called in the last stage of a time step, before the position of the body is updated: after the last
// sample of an interval the hull is replaced by its displacement body.
void sixdof_obj_nhflow::shipbl_step(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(bs.mode<2 || bs.done || !bs.ready || n6DOF!=0)
    return;

    const int M = bs.M;
    const double tnew = bs.tlast;
    const int dumping = int(envd("NHF_BL_DUMP",0.0));
    bs.ready = false;

    // ------------------------------------------------ displacement body from the mean flow, the hull is replaced
    ++bs.iup;
    const bool lastup = bs.iup>=bs.nup;
    const double tsum = bs.tsum;
    bs.done = lastup;
    bs.t0 = bs.t1;
    bs.t1 += bs.w;
    pgc->globalsum(bs.acc.data(),8*M);
    if(dumping!=0)
    pgc->globalsum(bs.last.data(),8*M);

    int nnew=0;
    std::vector<double> buf;

    if(p->mpirank==0)
    {
        nhflow_shipbl &bl = *bs.bl;

        // port and starboard: u, w even, v odd. Points that the grid does not hold: the double body alone
        std::vector<double> vs(3*M), off(3);
        int nmiss=0;
        for(int m=0; m<M; ++m)
        {
            const double *a = &bs.acc[4*m], *b = &bs.acc[4*(M+m)];
            if(a[3]<0.5*tsum || b[3]<0.5*tsum)
            {
                ++nmiss;
                vs[3*m] = vs[3*m+1] = vs[3*m+2] = 1.0e30;
                continue;
            }
            vs[3*m]   = 0.5*(a[0]/a[3] + b[0]/b[3]);
            vs[3*m+1] = 0.5*(a[1]/a[3] - b[1]/b[3]);
            vs[3*m+2] = 0.5*(a[2]/a[3] + b[2]/b[3]);
        }

        if(dumping!=0 && lastup)
        {
            FILE *fp = fopen("./REEF3D_NHFLOW_6DOF/NHF_BL_edge.bin","wb");
            if(fp!=nullptr)
            {
                fwrite(bs.pt.data(),sizeof(double),bs.pt.size(),fp);        // M x (x, y - yc, z)
                fwrite(vs.data(),sizeof(double),vs.size(),fp);              // M x (u, v, w), mean, symmetric part
                fwrite(bs.last.data(),sizeof(double),bs.last.size(),fp);    // 2 M x (u, v, w, 1) of the last time step
                fclose(fp);
            }
        }

        bl.set_edge(bs.ph,vs,1.0e29);

        std::vector<nhflow_shipbl::tri3> body;
        bl.solve();
        bl.body(body);
        attitude(body,buf);
        nnew = int(body.size());

        std::cout<<"NHF_BL="<<bs.mode<<": t = "<<tnew<<" s, update "<<bs.iup<<" of "<<bs.nup<<", the hull is replaced by its displacement body from the mean flow of the last "
                 <<tsum<<" s ("<<M<<" sample points on each side, "<<nmiss<<" outside the grid)"<<std::endl<<bl.log<<std::endl;

        bl.log.clear();

        if(dumping!=0 && lastup)
        write_stl("./REEF3D_NHFLOW_6DOF/NHF_BL_dispbody2.stl",buf);

        if(lastup)
        {
            delete bs.bl;
            bs.bl = nullptr;
        }
    }

    bs.acc.assign(8*M,0.0);
    bs.tsum = 0.0;

    pgc->bcast_int(&nnew,1);
    if(p->mpirank!=0)
    buf.resize(9*nnew);
    pgc->bcast_double(buf.data(),9*nnew);

    const int ts = bs.ts, told = tricount;

    p->Dresize(tri_x,told,ts+nnew,3,3);
    p->Dresize(tri_y,told,ts+nnew,3,3);
    p->Dresize(tri_z,told,ts+nnew,3,3);
    p->Dresize(tri_x0,told,ts+nnew,3,3);
    p->Dresize(tri_y0,told,ts+nnew,3,3);
    p->Dresize(tri_z0,told,ts+nnew,3,3);

    for(int n=0; n<nnew; ++n)
    for(int q=0; q<3; ++q)
    {
        tri_x[ts+n][q] = buf[9*n+3*q];
        tri_y[ts+n][q] = buf[9*n+3*q+1];
        tri_z[ts+n][q] = buf[9*n+3*q+2];
    }
    tricount = ts+nnew;
    for(int e=0; e<entity_count; ++e)
    if(tstart[e]==ts)
    tend[e] = tricount;

    // orientation and remeshing as after reading (sixdof_obj::objects_create), then the body frame (the hull does not
    // move: the pose of the start). Mass and inertia stay
    geom.orient(p,pgc);

    if(p->X185>0 && p->X60!=2 && entity_count>0 && p->X170==0 && p->X171==0 && p->X172==0)
    geom.refine(p,pgc);

    if(p->mpirank==0)
    std::cout<<"NHF_BL: refined surface triangles: "<<tricount<<std::endl;

    for(int n=0; n<tricount; ++n)
    for(int q=0; q<3; ++q)
    {
        const Eigen::Vector3d r = R_.transpose()*Eigen::Vector3d(tri_x[n][q]-c_(0),tri_y[n][q]-c_(1),tri_z[n][q]-c_(2));
        tri_x0[n][q] = r(0);
        tri_y0[n][q] = r(1);
        tri_z0[n][q] = r(2);
    }
}
