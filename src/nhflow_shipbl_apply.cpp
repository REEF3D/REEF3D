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
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<vector>

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
//   NHF_BL_LOADS     the loads of the body (6DOF_obj_nhflow_forces_stl.cpp) are taken on the hull itself, not on its
//                     displacement body (0: on the body, as before; 2: on the kept hull triangles also before the first
//                     update of NHF_BL=2, a test). The body is 10 % larger than the hull (KCS, with a wake body behind the
//                     stern) and that volume counted as buoyancy. The hull as it was read is kept (start pose, outward
//                     winding) and follows the pose of the body; the hydrostatic part acts on it, the pressure of a hull
//                     triangle is taken at the displacement body, at the distance to it along the outward normal (the
//                     pressure across the boundary layer is the one at its edge), where the normal does not meet the
//                     body within NHF_BL_LMAX (default 0.02 L) at the hull.     default 1
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
        // pose of the hull (every rank): at the start of the run and at the end of the last time step, with the velocity
        // of the body there (u_fb). moved: the pose differs from the start or the hull has a velocity
        bool pose0 = false, moved = false;
        Eigen::Matrix3d R0 = Eigen::Matrix3d::Identity(), R = Eigen::Matrix3d::Identity();
        Eigen::Vector3d c0 = Eigen::Vector3d::Zero(), c = Eigen::Vector3d::Zero();
        Eigen::Matrix<double,6,1> ufb = Eigen::Matrix<double,6,1>::Zero();
        // loads on the hull itself (NHF_BL_LOADS, every rank): the hull as read (start pose, outward winding), for every
        // triangle the distance to the displacement body along the normal; swapped: the solver runs with the body
        int loads = 1;
        bool swapped = false, lstill = true;
        double L = 0.0;
        std::vector<double> hull, hoff;
        Eigen::Matrix3d lR = Eigen::Matrix3d::Identity();
        Eigen::Vector3d lc = Eigen::Vector3d::Zero();
    };
    blstate bs;

    double envd(const char *name, double def);

    // triangles (9 numbers each) in bins along x, rays against them (Moeller-Trumbore, both faces)
    struct raybins
    {
        const double *T = nullptr;
        int N = 0, nb = 1;
        double x0 = 0.0, h = 1.0, diag = 1.0;
        std::vector<int> start, id, stamp;
        int mark = 0;

        inline int bin(double x) const
        {
            const double r = (x-x0)/h;
            return r<=0.0 ? 0 : (r>=double(nb) ? nb-1 : int(r));
        }

        void build(const double *tri, int num, double hbin)
        {
            T = tri;
            N = num;
            double lo[3] = {1.0e30,1.0e30,1.0e30}, hi[3] = {-1.0e30,-1.0e30,-1.0e30};
            for(int n=0; n<3*N; ++n)
            for(int q=0; q<3; ++q)
            {
                lo[q] = std::min(lo[q],T[3*n+q]);
                hi[q] = std::max(hi[q],T[3*n+q]);
            }
            if(N==0)
            lo[0] = lo[1] = lo[2] = hi[0] = hi[1] = hi[2] = 0.0;
            diag = sqrt((hi[0]-lo[0])*(hi[0]-lo[0]) + (hi[1]-lo[1])*(hi[1]-lo[1]) + (hi[2]-lo[2])*(hi[2]-lo[2])) + 1.0e-12;
            x0 = lo[0];
            nb = std::max(1,std::min(8192,int((hi[0]-lo[0])/std::max(hbin,1.0e-12))+1));
            h = (hi[0]-lo[0])/double(nb) + 1.0e-12;
            start.assign(nb+1,0);
            for(int pass=0; pass<2; ++pass)
            {
                std::vector<int> pos(start.begin(),start.end()-1);
                for(int n=0; n<N; ++n)
                {
                    const double *a = &T[9*n];
                    const int b0 = bin(std::min(a[0],std::min(a[3],a[6]))), b1 = bin(std::max(a[0],std::max(a[3],a[6])));
                    for(int b=b0; b<=b1; ++b)
                    if(pass==0)
                    ++start[b+1];
                    else
                    id[pos[b]++] = n;
                }
                if(pass==0)
                {
                    for(int b=0; b<nb; ++b)
                    start[b+1] += start[b];
                    id.resize(start[nb]);
                }
            }
            stamp.assign(N,0);
            mark = 0;
        }

        // the ray o + t d (d a unit vector), tlo <= t <= thi: number of triangles met (not the triangle skip), the smallest t
        int cast(const double *o, const double *d, double tlo, double thi, int skip, double &tmin)
        {
            ++mark;
            int count = 0;
            tmin = 1.0e30;
            const double xa = o[0]+tlo*d[0], xb = o[0]+thi*d[0];
            const int b0 = bin(std::min(xa,xb)), b1 = bin(std::max(xa,xb));
            for(int b=b0; b<=b1; ++b)
            for(int k=start[b]; k<start[b+1]; ++k)
            {
                const int n = id[k];
                if(n==skip || stamp[n]==mark)
                continue;
                stamp[n] = mark;
                const double *a = &T[9*n];
                const double e1[3] = {a[3]-a[0],a[4]-a[1],a[5]-a[2]}, e2[3] = {a[6]-a[0],a[7]-a[1],a[8]-a[2]};
                const double pv[3] = {d[1]*e2[2]-d[2]*e2[1], d[2]*e2[0]-d[0]*e2[2], d[0]*e2[1]-d[1]*e2[0]};
                const double det = e1[0]*pv[0] + e1[1]*pv[1] + e1[2]*pv[2];
                if(fabs(det)<1.0e-300)
                continue;
                const double tv[3] = {o[0]-a[0],o[1]-a[1],o[2]-a[2]};
                const double u = (tv[0]*pv[0] + tv[1]*pv[1] + tv[2]*pv[2])/det;
                if(u<0.0 || u>1.0)
                continue;
                const double qv[3] = {tv[1]*e1[2]-tv[2]*e1[1], tv[2]*e1[0]-tv[0]*e1[2], tv[0]*e1[1]-tv[1]*e1[0]};
                const double v = (d[0]*qv[0] + d[1]*qv[1] + d[2]*qv[2])/det;
                if(v<0.0 || u+v>1.0)
                continue;
                const double t = (e2[0]*qv[0] + e2[1]*qv[1] + e2[2]*qv[2])/det;
                if(t<tlo || t>thi)
                continue;
                ++count;
                tmin = std::min(tmin,t);
            }
            return count;
        }
    };

    // unit normal of the winding and area of a triangle (9 numbers); false: no area
    inline bool tri_normal(const double *a, double *nv, double &area)
    {
        const double e1[3] = {a[3]-a[0],a[4]-a[1],a[5]-a[2]}, e2[3] = {a[6]-a[0],a[7]-a[1],a[8]-a[2]};
        nv[0] = e1[1]*e2[2]-e1[2]*e2[1];
        nv[1] = e1[2]*e2[0]-e1[0]*e2[2];
        nv[2] = e1[0]*e2[1]-e1[1]*e2[0];
        const double norm = sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);
        area = 0.5*norm;
        if(norm<1.0e-20)
        return false;
        nv[0] /= norm;
        nv[1] /= norm;
        nv[2] /= norm;
        return true;
    }

    // outward winding of the kept hull: a ray from a triangle along its normal leaves a closed surface through an even
    // number of other triangles (three rays, slightly apart, the majority counts). Returns the number of triangles turned
    int orient_hull(std::vector<double> &H)
    {
        const int N = int(H.size()/9);
        raybins rb;
        rb.build(H.data(),N,0.0);
        rb.build(H.data(),N,rb.diag/2000.0);
        const double tilt[3][2] = {{0.0,0.0},{0.02,0.011},{-0.013,0.017}};
        std::vector<char> turn(N,0);
        int nturn = 0;
        for(int n=0; n<N; ++n)
        {
            const double *a = &H[9*n];
            double nv[3], area;
            if(!tri_normal(a,nv,area))
            continue;
            const double o[3] = {0.31*a[0]+0.33*a[3]+0.36*a[6], 0.31*a[1]+0.33*a[4]+0.36*a[7], 0.31*a[2]+0.33*a[5]+0.36*a[8]};
            double t1[3] = {a[3]-a[0],a[4]-a[1],a[5]-a[2]};
            const double l1 = sqrt(t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2]);
            for(int q=0; q<3; ++q)
            t1[q] /= l1;
            const double t2[3] = {nv[1]*t1[2]-nv[2]*t1[1], nv[2]*t1[0]-nv[0]*t1[2], nv[0]*t1[1]-nv[1]*t1[0]};
            int out = 0;
            for(int r=0; r<3; ++r)
            {
                double d[3], tmin;
                for(int q=0; q<3; ++q)
                d[q] = nv[q] + tilt[r][0]*t1[q] + tilt[r][1]*t2[q];
                const double dl = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
                for(int q=0; q<3; ++q)
                d[q] /= dl;
                if(rb.cast(o,d,1.0e-9*rb.diag,rb.diag,n,tmin)%2==0)
                ++out;
            }
            if(out<2)
            {
                turn[n] = 1;
                ++nturn;
            }
        }
        for(int n=0; n<N; ++n)
        if(turn[n])
        for(int q=0; q<3; ++q)
        std::swap(H[9*n+3+q],H[9*n+6+q]);
        return nturn;
    }

    // distance from every hull triangle to the body (triangles buf, start pose) along the outward normal
    void hull_offsets(const std::vector<double> &buf, bool log)
    {
        const int N = int(bs.hull.size()/9);
        const double dmax = envd("NHF_BL_LMAX",0.02*bs.L), back = 0.001*bs.L;
        raybins rb;
        rb.build(buf.data(),int(buf.size()/9),0.25*dmax);
        bs.hoff.assign(N,0.0);
        double asum=0.0, vsum=0.0, omax=0.0;
        int nmiss=0;
        for(int n=0; n<N; ++n)
        {
            const double *a = &bs.hull[9*n];
            double nv[3], area, tmin;
            if(!tri_normal(a,nv,area))
            continue;
            const double o[3] = {(a[0]+a[3]+a[6])/3.0, (a[1]+a[4]+a[7])/3.0, (a[2]+a[5]+a[8])/3.0};
            if(rb.cast(o,nv,-back,dmax,-1,tmin)>0)
            bs.hoff[n] = std::max(tmin,0.0);
            else
            ++nmiss;
            asum += area;
            vsum += area*bs.hoff[n];
            omax = std::max(omax,bs.hoff[n]);
        }
        bs.swapped = true;
        if(log)
        std::cout<<"NHF_BL: the loads are taken on the hull itself ("<<N<<" triangles, "<<asum<<" m^2), the pressure at the displacement body: distance mean "
                 <<1.0e3*vsum/std::max(asum,1.0e-30)<<" mm, largest "<<1.0e3*omax<<" mm, volume between them "<<vsum<<" m^3, "<<nmiss<<" triangles whose normal does not meet the body"<<std::endl;
    }

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
    bs.L = L;
    bs.loads = int(envd("NHF_BL_LOADS",1.0));
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

    // the hull itself is kept for the loads (every rank)
    if(bs.loads!=0)
    {
        bs.hull.resize(9*(tricount-ts));
        for(int n=ts; n<tricount; ++n)
        for(int q=0; q<3; ++q)
        {
            bs.hull[9*(n-ts)+3*q]   = tri_x[n][q];
            bs.hull[9*(n-ts)+3*q+1] = tri_y[n][q];
            bs.hull[9*(n-ts)+3*q+2] = tri_z[n][q];
        }
        const int nturn = orient_hull(bs.hull);
        bs.hoff.assign(tricount-ts,0.0);
        double vol = 0.0;
        for(int n=0; n<tricount-ts; ++n)
        {
            const double *a = &bs.hull[9*n];
            vol += ((a[3]-a[0])*(a[7]-a[1]) - (a[4]-a[1])*(a[6]-a[0]))*(a[2]+a[5]+a[8])/6.0;
        }
        if(p->mpirank==0)
        std::cout<<"NHF_BL: hull kept for the loads, "<<tricount-ts<<" triangles, "<<nturn<<" turned to the outward winding, volume "<<vol<<" m^3"<<std::endl;
    }

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

    if(bs.loads!=0)
    hull_offsets(buf,p->mpirank==0);

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

    // A hull that moves (free sinkage and trim, motion in waves): the sample points are fixed to the hull, x = R b + c with
    // b = R0^T (x0 - c0) of the point x0 at the start pose, and the velocity is the one relative to the hull at the point,
    // v - (U_b + omega x r), turned back to the frame of the start pose, in which the boundary layer is computed
    const Eigen::Matrix3d Rrel = bs.R*bs.R0.transpose();
    const Eigen::Matrix3d Rback = Rrel.transpose();

    for(int m=0; m<M; ++m)
    for(int side=0; side<2; ++side)
    {
        double v[3];
        const int n = 4*(side*M + m);
        bs.last[n]=bs.last[n+1]=bs.last[n+2]=bs.last[n+3]=0.0;
        double xs = bs.pt[3*m], ys = bs.yc + (side==0 ? 1.0 : -1.0)*bs.pt[3*m+1], zs = bs.pt[3*m+2];
        Eigen::Vector3d r = Eigen::Vector3d::Zero();
        if(bs.moved)
        {
            r = Rrel*(Eigen::Vector3d(xs,ys,zs) - bs.c0);
            xs = bs.c(0) + r(0);
            ys = bs.c(1) + r(1);
            zs = bs.c(2) + r(2);
        }
        if(!sample(xs,ys,zs,v))
        continue;
        if(bs.moved && !coords)
        {
            const Eigen::Vector3d vr = Rback*Eigen::Vector3d(v[0] - (bs.ufb(0) + bs.ufb(4)*r(2) - bs.ufb(5)*r(1)),
                                                             v[1] - (bs.ufb(1) + bs.ufb(5)*r(0) - bs.ufb(3)*r(2)),
                                                             v[2] - (bs.ufb(2) + bs.ufb(3)*r(1) - bs.ufb(4)*r(0)));
            v[0] = vr(0);
            v[1] = vr(1);
            v[2] = vr(2);
        }
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
    if(bs.mode<1 || n6DOF!=0)
    return;

    // d = nullptr: only the pose and the velocity of the hull are taken (after the initialisation: the start pose; after
    // the position and velocity update of the last stage: the state at the end of the time step, used by the samples)
    if(d==nullptr)
    {
        if(!bs.pose0)
        {
            bs.R0 = R_;
            bs.c0 = c_;
            bs.pose0 = true;
        }
        bs.R = R_;
        bs.c = c_;
        bs.ufb = u_fb;
        const bool moved = (bs.R-bs.R0).cwiseAbs().maxCoeff()>1.0e-12 || (bs.c-bs.c0).cwiseAbs().maxCoeff()>1.0e-12
                           || bs.ufb.cwiseAbs().maxCoeff()>1.0e-12;
        if(moved && !bs.moved && p->mpirank==0 && bs.mode>=2)
        std::cout<<"NHF_BL: the hull moves, the sample points follow it and the velocity is taken relative to it (t = "<<p->simtime<<" s)"<<std::endl;
        if(moved)
        bs.moved = true;
        return;
    }

    if(bs.mode<2 || bs.done || !bs.ready)
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

    if(bs.loads!=0)
    hull_offsets(buf,p->mpirank==0);

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

    // orientation and remeshing as after reading (sixdof_obj::objects_create), then the body frame: the new body is
    // built at the start pose, b = R0^T (x0 - c0); the position update that follows takes it to the pose of the hull
    // now. Mass and inertia stay
    geom.orient(p,pgc);

    if(p->X185>0 && p->X60!=2 && entity_count>0 && p->X170==0 && p->X171==0 && p->X172==0)
    geom.refine(p,pgc);

    if(p->mpirank==0)
    std::cout<<"NHF_BL: refined surface triangles: "<<tricount<<std::endl;

    for(int n=0; n<tricount; ++n)
    for(int q=0; q<3; ++q)
    {
        const Eigen::Vector3d r = bs.moved ? Eigen::Vector3d(bs.R0.transpose()*Eigen::Vector3d(tri_x[n][q]-bs.c0(0),tri_y[n][q]-bs.c0(1),tri_z[n][q]-bs.c0(2)))
                                           : Eigen::Vector3d(R_.transpose()*Eigen::Vector3d(tri_x[n][q]-c_(0),tri_y[n][q]-c_(1),tri_z[n][q]-c_(2)));
        tri_x0[n][q] = r(0);
        tri_y0[n][q] = r(1);
        tri_z0[n][q] = r(2);
    }
}

// NHF_BL_LOADS: the triangles the loads of the body are taken on (6DOF_obj_nhflow_forces_stl.cpp). Returns the number of
// kept hull triangles (0: the triangles of the body as they are) and the first triangle of the hull in the body, ts; the
// pose R, c of the body now places them
int nhflow_shipbl_loads(const Eigen::Matrix3d &R, const Eigen::Vector3d &c, int &ts)
{
    if(bs.mode==0 || bs.hull.empty() || !(bs.loads==2 || (bs.loads==1 && bs.swapped)))
    return 0;

    ts = bs.ts;
    bs.lstill = !bs.pose0 || (R==bs.R0 && c==bs.c0);
    bs.lR = R*bs.R0.transpose();
    bs.lc = c;
    return int(bs.hull.size()/9);
}

// hull triangle n at the pose given to nhflow_shipbl_loads, and the distance to the displacement body along its normal
void nhflow_shipbl_loadtri(int n, double *vx, double *vy, double *vz, double &off)
{
    const double *h = &bs.hull[9*n];
    off = bs.hoff[n];
    for(int q=0; q<3; ++q)
    {
        if(bs.lstill)
        {
            vx[q] = h[3*q];
            vy[q] = h[3*q+1];
            vz[q] = h[3*q+2];
            continue;
        }
        const Eigen::Vector3d r = bs.lR*(Eigen::Vector3d(h[3*q],h[3*q+1],h[3*q+2]) - bs.c0) + bs.lc;
        vx[q] = r(0);
        vy[q] = r(1);
        vz[q] = r(2);
    }
}
