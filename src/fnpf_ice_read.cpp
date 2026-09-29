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

#include"fnpf_ice.h"
#include"lexer.h"
#include"ghostcell.h"
#include<sstream>
#include<string>
#include<random>
#include<algorithm>

// ice_floes.dat, one entry per line, '#' starts a comment, angles in degrees:
//
//   density  rho                               ice density for the entries that follow (default 917)
//   rect     xc yc h Lx Ly yaw                 rectangular floe
//   poly     xc yc h D n yaw                   regular n-gon, circumscribed diameter D
//   convex   h n x1 y1 ... xn yn               convex floe from world vertices (convex hull is taken)
//   field    xs xe ys ye h conc Dmin Dmax nmin nmax seed
//                                              random field of convex floes, n in [nmin,nmax] vertices,
//                                              diameters uniform in [Dmin,Dmax], placed without overlap
//                                              until the area concentration conc is reached
//   obstacle n x1 y1 ... xn yn                 fixed convex polygon, contact only, ice load is logged
//   disp     dz roll pitch                     initial offset of the previous floe from equilibrium
//   vel      u v                               initial drift velocity of the previous floe
//
// 2D flume (j_dir 0): only rect is meaningful, yc and yaw are ignored, Ly is the floe width
// the pressure acts on (per-unit-width results for Ly = 1).

void fnpf_ice::read(lexer *p, ghostcell *pgc)
{
    vector<double> buf;
    int size=0;

    if(p->mpirank==0)
    {
    read_file(p);
    serialize_geometry(buf);
    size = int(buf.size());
    }

    if(p->mpi_size>1)
    {
    MPI_Bcast(&size, 1, MPI_INT, 0, pgc->mpi_comm);

    if(p->mpirank>0)
    buf.resize(size);

    MPI_Bcast(buf.data(), size, MPI_DOUBLE, 0, pgc->mpi_comm);

    if(p->mpirank>0)
    deserialize_geometry(p,buf);
    }

    nfloe=nobst=0;
    for(size_t n=0; n<floe.size(); ++n)
    {
        fnpf_ice_floe &fl = floe[n];
        fl.id = int(n);

        if(fl.type==0)
        {
        ++nfloe;
        fl.pw   = fl.rho*g*fl.h;
        fl.klid = alpha*rhow*g;
        // damping: fraction zeta of critical for the floe mass per area on the lid spring
        fl.clid = zeta*2.0*sqrt(fl.klid*fl.rho*fl.h);
        }

        if(fl.type==1)
        ++nobst;
    }
}

void fnpf_ice::read_file(lexer *p)
{
    ifstream in("ice_floes.dat");

    if(!in.is_open())
    {
    cout<<"!!! FNPF ice (A 380): could not open ice_floes.dat !!!"<<endl;
    return;
    }

    double rho = 917.0;
    const double deg = PI/180.0;
    string line;
    int lnum=0;

    while(getline(in,line))
    {
        ++lnum;

        size_t hash = line.find('#');
        if(hash!=string::npos)
        line = line.substr(0,hash);

        istringstream ss(line);
        string key;

        if(!(ss>>key))
        continue;

        if(key=="density")
        {
        ss>>rho;
        }

        else if(key=="rect")
        {
        double xc,yc,h,Lx,Ly,yaw=0.0;
        ss>>xc>>yc>>h>>Lx>>Ly>>yaw;

        fnpf_ice_floe fl;
        fl.type=0;
        fl.h=h;
        fl.rho=rho;
        fl.bx = {-0.5*Lx, 0.5*Lx, 0.5*Lx, -0.5*Lx};
        fl.by = {-0.5*Ly,-0.5*Ly, 0.5*Ly,  0.5*Ly};
        geometry(fl);
        place(p,fl,xc,yc,yaw*deg,0.0,0.0,0.0);
        floe.push_back(fl);
        }

        else if(key=="poly")
        {
        double xc,yc,h,D,yaw=0.0;
        int nv;
        ss>>xc>>yc>>h>>D>>nv>>yaw;
        nv = MAX(nv,3);

        fnpf_ice_floe fl;
        fl.type=0;
        fl.h=h;
        fl.rho=rho;
        for(int q=0; q<nv; ++q)
        {
        fl.bx.push_back(0.5*D*cos(2.0*PI*q/double(nv)));
        fl.by.push_back(0.5*D*sin(2.0*PI*q/double(nv)));
        }
        geometry(fl);
        place(p,fl,xc,yc,yaw*deg,0.0,0.0,0.0);
        floe.push_back(fl);
        }

        else if(key=="convex" || key=="obstacle")
        {
        double h=0.0;
        int nv;

        if(key=="convex")
        ss>>h;

        ss>>nv;
        vector<double> xs(nv),ys(nv);
        for(int q=0; q<nv; ++q)
        ss>>xs[q]>>ys[q];

        add_polygon(p,xs,ys,h,key=="convex"?0:1);
        floe.back().rho = rho;
        if(floe.back().type==0)
        {
        geometry(floe.back());
        place(p,floe.back(),floe.back().x[0],floe.back().x[1],0.0,0.0,0.0,0.0);
        }
        }

        else if(key=="field")
        {
        double xs,xe,ys,ye,h,conc,Dmin,Dmax;
        int nmin,nmax,seed;
        ss>>xs>>xe>>ys>>ye>>h>>conc>>Dmin>>Dmax>>nmin>>nmax>>seed;

        size_t n0 = floe.size();
        floe_field(p,xs,xe,ys,ye,h,conc,Dmin,Dmax,nmin,nmax,seed);

        for(size_t n=n0; n<floe.size(); ++n)
        {
        floe[n].rho=rho;
        geometry(floe[n]);
        double roll,pitch,yaw;
        quat_to_euler(floe[n].q,roll,pitch,yaw);
        place(p,floe[n],floe[n].x[0],floe[n].x[1],yaw,0.0,0.0,0.0);
        }
        }

        else if(key=="disp")
        {
        double dz=0.0,roll=0.0,pitch=0.0;
        ss>>dz>>roll>>pitch;

        if(!floe.empty() && floe.back().type==0)
        {
        double r0,p0,yaw;
        quat_to_euler(floe.back().q,r0,p0,yaw);
        place(p,floe.back(),floe.back().x[0],floe.back().x[1],yaw,dz,roll*deg,pitch*deg);
        }
        }

        else if(key=="vel")
        {
        double u=0.0,v=0.0;
        ss>>u>>v;
        
        if(!floe.empty() && floe.back().type==0)
        {
        floe.back().v[0]=u;
        floe.back().v[1]=is2D ? 0.0 : v;
        }
        }
        
        else
        cout<<"FNPF ice: unknown entry '"<<key<<"' in ice_floes.dat, line "<<lnum<<endl;
    }

    in.close();
}

void fnpf_ice::add_polygon(lexer *p, vector<double> &xs, vector<double> &ys, double h, int type)
{
    convex_hull(xs,ys);

    fnpf_ice_floe fl;
    fl.type = type;
    fl.h = h;
    fl.bx = xs;
    fl.by = ys;

    // geometry() shifts the vertices to the centroid, the shift is the world position
    double A=0.0,cx=0.0,cy=0.0;
    const int nv=int(xs.size());
    for(int q=0; q<nv; ++q)
    {
    const int q2=(q+1)%nv;
    const double cr = xs[q]*ys[q2]-xs[q2]*ys[q];
    A  += 0.5*cr;
    cx += (xs[q]+xs[q2])*cr;
    cy += (ys[q]+ys[q2])*cr;
    }
    cx /= (6.0*A);
    cy /= (6.0*A);

    geometry(fl);

    fl.x[0]=cx;
    fl.x[1]=cy;
    fl.x[2]=wd;

    floe.push_back(fl);
}

void fnpf_ice::floe_field(lexer *p, double xs, double xe, double ys, double ye, double h, double conc,
                          double Dmin, double Dmax, int nmin, int nmax, int seed)
{
    mt19937 gen(seed);
    uniform_real_distribution<double> U(0.0,1.0);

    nmin = MAX(nmin,3);
    nmax = MAX(nmax,nmin);
    Dmax = MAX(Dmax,Dmin);

    const double atot = is2D ? (xe-xs) : (xe-xs)*(ye-ys);
    const double target = conc*atot;
    double acov=0.0;

    // candidates are only checked against floes of this field
    const size_t n0 = floe.size();
    int attempts=0;
    const int maxattempts=200000;

    while(acov<target && attempts<maxattempts)
    {
        ++attempts;

        const double D = Dmin + (Dmax-Dmin)*U(gen);

        fnpf_ice_floe fl;
        fl.type=0;
        fl.h=h;

        if(is2D)
        {
        fl.bx = {-0.5*D, 0.5*D, 0.5*D, -0.5*D};
        fl.by = {-0.5, -0.5, 0.5, 0.5};
        }
        else
        {
        const int nv = nmin + int((nmax-nmin+1)*U(gen)*0.999999);
        vector<double> ang(nv);
        for(auto &a : ang)
        a = 2.0*PI*U(gen);
        sort(ang.begin(),ang.end());

        vector<double> vx(nv),vy(nv);
        for(int q=0; q<nv; ++q)
        {
        const double r = 0.5*D*(0.8 + 0.2*U(gen));
        vx[q] = r*cos(ang[q]);
        vy[q] = r*sin(ang[q]);
        }
        convex_hull(vx,vy);

        if(vx.size()<3)
        continue;

        fl.bx=vx;
        fl.by=vy;
        }

        geometry(fl);
        
        // no slivers: at least 40% of the circumscribed disc and a width of half the diameter
        if(!is2D && (fl.area < 0.4*PI*fl.rbound*fl.rbound || fl.wmin < 0.5*D))
        continue;
        
        const double r = is2D ? 0.5*D : fl.rbound;
        const double gap = 0.02*D;

        if(xe-xs<=2.0*r || (!is2D && ye-ys<=2.0*r))
        continue;

        const double xc = xs + r + (xe-xs-2.0*r)*U(gen);
        const double yc = is2D ? 0.0 : ys + r + (ye-ys-2.0*r)*U(gen);

        // random yaw, stored in q for place()
        const double yaw = is2D ? 0.0 : 2.0*PI*U(gen);
        
        // world planform of the candidate
        const int nvc = int(fl.bx.size());
        vector<double> cx(nvc),cy(nvc);
        for(int q=0; q<nvc; ++q)
        {
        cx[q] = xc + cos(yaw)*fl.bx[q] - sin(yaw)*fl.by[q];
        cy[q] = yc + sin(yaw)*fl.bx[q] + cos(yaw)*fl.by[q];
        }
        
        // separating axis with a gap between convex polygons
        auto separated = [gap](const vector<double> &ax, const vector<double> &ay, const vector<double> &bx, const vector<double> &by)
        {
            for(int pass=0; pass<2; ++pass)
            {
                const vector<double> &px = pass ? bx : ax;
                const vector<double> &py = pass ? by : ay;
                const vector<double> &qx = pass ? ax : bx;
                const vector<double> &qy = pass ? ay : by;
                const int np = int(px.size());
                
                for(int q=0; q<np; ++q)
                {
                const int q2=(q+1)%np;
                const double ex = px[q2]-px[q], ey = py[q2]-py[q];
                const double len = sqrt(ex*ex+ey*ey);
                if(len<=0.0)
                continue;
                const double nx = ey/len, ny = -ex/len;
                double smin=1.0e20;
                for(size_t k=0; k<qx.size(); ++k)
                smin = MIN(smin, nx*(qx[k]-px[q]) + ny*(qy[k]-py[q]));
                if(smin>gap)
                return true;
                }
            }
            return false;
        };
        
        int ok=1;
        for(size_t n=n0; n<floe.size() && ok; ++n)
        {
            const double dx = floe[n].x[0]-xc;
            const double dy = is2D ? 0.0 : floe[n].x[1]-yc;
            const double rn = is2D ? 0.5*(*max_element(floe[n].bx.begin(),floe[n].bx.end()) - *min_element(floe[n].bx.begin(),floe[n].bx.end())) : floe[n].rbound;
            
            if(dx*dx+dy*dy >= (r+rn+gap)*(r+rn+gap))
            continue;
            
            if(is2D)
            {
            ok=0;
            continue;
            }
            
            // polygons of the placed floe
            const fnpf_ice_floe &fn = floe[n];
            const double yn = 2.0*atan2(fn.q[3],fn.q[0]);
            const int nvn = int(fn.bx.size());
            vector<double> px(nvn),py(nvn);
            for(int q=0; q<nvn; ++q)
            {
            px[q] = fn.x[0] + cos(yn)*fn.bx[q] - sin(yn)*fn.by[q];
            py[q] = fn.x[1] + sin(yn)*fn.bx[q] + cos(yn)*fn.by[q];
            }
            
            if(!separated(cx,cy,px,py))
            ok=0;
        }
        
        if(!ok)
        continue;
        
        fl.x[0]=xc;
        fl.x[1]=yc;
        fl.q[0]=cos(0.5*yaw);
        fl.q[1]=0.0;
        fl.q[2]=0.0;
        fl.q[3]=sin(0.5*yaw);

        acov += is2D ? D : fl.area;
        floe.push_back(fl);
    }

    if(p->mpirank==0)
    cout<<"FNPF ice: floe field with "<<floe.size()-n0<<" floes, concentration "<<acov/MAX(atot,1.0e-20)<<" (target "<<conc<<")"<<endl;
}

void fnpf_ice::geometry(fnpf_ice_floe &fl)
{
    const int nv = int(fl.bx.size());

    // orientation: counter-clockwise
    double A=0.0;
    for(int q=0; q<nv; ++q)
    {
    const int q2=(q+1)%nv;
    A += 0.5*(fl.bx[q]*fl.by[q2]-fl.bx[q2]*fl.by[q]);
    }

    if(A<0.0)
    {
    reverse(fl.bx.begin(),fl.bx.end());
    reverse(fl.by.begin(),fl.by.end());
    A=-A;
    }

    // centroid
    double cx=0.0,cy=0.0;
    for(int q=0; q<nv; ++q)
    {
    const int q2=(q+1)%nv;
    const double cr = fl.bx[q]*fl.by[q2]-fl.bx[q2]*fl.by[q];
    cx += (fl.bx[q]+fl.bx[q2])*cr;
    cy += (fl.by[q]+fl.by[q2])*cr;
    }
    cx /= (6.0*A);
    cy /= (6.0*A);

    for(int q=0; q<nv; ++q)
    {
    fl.bx[q]-=cx;
    fl.by[q]-=cy;
    }

    // second moments of area about the centroid
    double Jyy=0.0,Jxx=0.0,Jxy=0.0;   // int y^2, int x^2, int xy
    for(int q=0; q<nv; ++q)
    {
    const int q2=(q+1)%nv;
    const double x0=fl.bx[q], y0=fl.by[q], x1=fl.bx[q2], y1=fl.by[q2];
    const double cr = x0*y1-x1*y0;
    Jyy += (y0*y0 + y0*y1 + y1*y1)*cr/12.0;
    Jxx += (x0*x0 + x0*x1 + x1*x1)*cr/12.0;
    Jxy += (x0*y1 + 2.0*x0*y0 + 2.0*x1*y1 + x1*y0)*cr/24.0;
    }

    fl.area = A;
    fl.mass = fl.rho*fl.h*A;

    const double rh = fl.rho*fl.h;
    const double h2 = fl.mass*fl.h*fl.h/12.0;

    fl.Ib[0][0] = rh*Jyy + h2;
    fl.Ib[1][1] = rh*Jxx + h2;
    fl.Ib[2][2] = rh*(Jxx+Jyy);
    fl.Ib[0][1] = fl.Ib[1][0] = -rh*Jxy;
    fl.Ib[0][2] = fl.Ib[2][0] = 0.0;
    fl.Ib[1][2] = fl.Ib[2][1] = 0.0;

    fl.rbound=0.0;
    double xmin=1.0e20,xmax=-1.0e20;
    for(int q=0; q<nv; ++q)
    {
    fl.rbound = MAX(fl.rbound, sqrt(fl.bx[q]*fl.bx[q] + fl.by[q]*fl.by[q]));
    xmin = MIN(xmin,fl.bx[q]);
    xmax = MAX(xmax,fl.bx[q]);
    }

    fl.width2D = A/MAX(xmax-xmin,1.0e-20);
    
    // minimum caliper width: for a convex polygon the smallest over the edges of the largest vertex distance
    fl.wmin = 1.0e20;
    for(int q=0; q<nv; ++q)
    {
    const int q2=(q+1)%nv;
    const double ex = fl.bx[q2]-fl.bx[q], ey = fl.by[q2]-fl.by[q];
    const double len = sqrt(ex*ex+ey*ey);
    if(len<=0.0)
    continue;
    double dmax=0.0;
    for(int r=0; r<nv; ++r)
    dmax = MAX(dmax, fabs(-ey*(fl.bx[r]-fl.bx[q]) + ex*(fl.by[r]-fl.by[q]))/len);
    fl.wmin = MIN(fl.wmin,dmax);
    }
    if(fl.wmin>1.0e19)
    fl.wmin=0.0;
}

void fnpf_ice::place(lexer *p, fnpf_ice_floe &fl, double xc, double yc, double yaw, double dz, double roll, double pitch)
{
    if(is2D)
    {
    yc = 0.5*(p->global_ymin + p->global_ymax);
    yaw = 0.0;
    roll = 0.0;
    }

    fl.x[0] = xc;
    fl.x[1] = yc;
    fl.x[2] = (fl.type==0) ? wd - fl.rho*fl.h/rhow + 0.5*fl.h + dz : wd;

    // ZYX: q = qz(yaw)*qy(pitch)*qx(roll)
    const double cr=cos(0.5*roll),  sr=sin(0.5*roll);
    const double cp=cos(0.5*pitch), sp=sin(0.5*pitch);
    const double cyw=cos(0.5*yaw),  syw=sin(0.5*yaw);

    fl.q[0] = cyw*cp*cr + syw*sp*sr;
    fl.q[1] = cyw*cp*sr - syw*sp*cr;
    fl.q[2] = cyw*sp*cr + syw*cp*sr;
    fl.q[3] = syw*cp*cr - cyw*sp*sr;

    for(int q=0;q<3;++q)
    fl.v[q]=fl.w[q]=0.0;
}

void fnpf_ice::serialize_geometry(vector<double> &buf)
{
    buf.clear();
    buf.push_back(double(floe.size()));

    for(auto &fl : floe)
    {
    buf.push_back(double(fl.type));
    buf.push_back(fl.h);
    buf.push_back(fl.rho);
    buf.push_back(double(fl.bx.size()));
    for(auto v : fl.bx) buf.push_back(v);
    for(auto v : fl.by) buf.push_back(v);
    for(int q=0;q<3;++q) buf.push_back(fl.x[q]);
    for(int q=0;q<4;++q) buf.push_back(fl.q[q]);
    for(int q=0;q<3;++q) buf.push_back(fl.v[q]);
    }
}

void fnpf_ice::deserialize_geometry(lexer *p, vector<double> &buf)
{
    floe.clear();
    size_t pos=0;
    const int nf = int(buf[pos++]);

    for(int n=0; n<nf; ++n)
    {
    fnpf_ice_floe fl;
    fl.type = int(buf[pos++]);
    fl.h = buf[pos++];
    fl.rho = buf[pos++];
    const int nv = int(buf[pos++]);
    fl.bx.assign(buf.begin()+pos, buf.begin()+pos+nv); pos+=nv;
    fl.by.assign(buf.begin()+pos, buf.begin()+pos+nv); pos+=nv;
    for(int q=0;q<3;++q) fl.x[q]=buf[pos++];
    for(int q=0;q<4;++q) fl.q[q]=buf[pos++];
    for(int q=0;q<3;++q) fl.v[q]=buf[pos++];

    // vertices are already centred and counter-clockwise, geometry() leaves them unchanged
    geometry(fl);
    floe.push_back(fl);
    }
}

void fnpf_ice::convex_hull(vector<double> &xs, vector<double> &ys)
{
    // Andrew's monotone chain, counter-clockwise, collinear points removed
    const int n = int(xs.size());
    vector<pair<double,double>> pt(n);
    for(int q=0; q<n; ++q)
    pt[q] = make_pair(xs[q],ys[q]);

    sort(pt.begin(),pt.end());
    pt.erase(unique(pt.begin(),pt.end()),pt.end());

    const int m = int(pt.size());
    if(m<3)
    return;

    auto cross = [](const pair<double,double> &o, const pair<double,double> &a, const pair<double,double> &b)
    {
        return (a.first-o.first)*(b.second-o.second) - (a.second-o.second)*(b.first-o.first);
    };

    vector<pair<double,double>> H(2*m);
    int k=0;
    for(int q=0; q<m; ++q)
    {
    while(k>=2 && cross(H[k-2],H[k-1],pt[q])<=0.0) --k;
    H[k++]=pt[q];
    }
    for(int q=m-2, t=k+1; q>=0; --q)
    {
    while(k>=t && cross(H[k-2],H[k-1],pt[q])<=0.0) --k;
    H[k++]=pt[q];
    }
    H.resize(k-1);

    xs.resize(H.size());
    ys.resize(H.size());
    for(size_t q=0; q<H.size(); ++q)
    {
    xs[q]=H[q].first;
    ys[q]=H[q].second;
    }
}
