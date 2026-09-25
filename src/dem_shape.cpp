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

#include"dem_core.h"
#include<cmath>
#include<fstream>
#include<sstream>
#include<iostream>
#include<map>
#include<array>
#include<algorithm>
#include<cstring>

namespace
{
    const double dem_pi = 3.14159265358979323846;

    // Ericson, Real-Time Collision Detection, 5.1.5
    dem_vec closest_point_triangle(const dem_vec &p, const dem_vec &a, const dem_vec &b, const dem_vec &c)
    {
        dem_vec ab = b-a, ac = c-a, ap = p-a;
        double d1 = ab.dot(ap), d2 = ac.dot(ap);
        if(d1<=0.0 && d2<=0.0) return a;

        dem_vec bp = p-b;
        double d3 = ab.dot(bp), d4 = ac.dot(bp);
        if(d3>=0.0 && d4<=d3) return b;

        double vc = d1*d4 - d3*d2;
        if(vc<=0.0 && d1>=0.0 && d3<=0.0)
        return a + (d1/(d1-d3))*ab;

        dem_vec cp = p-c;
        double d5 = ab.dot(cp), d6 = ac.dot(cp);
        if(d6>=0.0 && d5<=d6) return c;

        double vb = d5*d2 - d1*d6;
        if(vb<=0.0 && d2>=0.0 && d6<=0.0)
        return a + (d2/(d2-d6))*ac;

        double va = d3*d6 - d5*d4;
        if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0)
        return b + ((d4-d3)/((d4-d3)+(d5-d6)))*(c-b);

        double denom = 1.0/(va+vb+vc);
        return a + ab*(vb*denom) + ac*(vc*denom);
    }
}

// ---------------------------------------------------------------------------------------------
// primitives
// ---------------------------------------------------------------------------------------------

void dem_shape::build_sphere(double r)
{
    type = DEM_SPHERE;
    dim = dem_vec(r,r,r);
    volume = 4.0/3.0*dem_pi*r*r*r;
    area = 4.0*dem_pi*r*r;
    rbound = r;
    deq = 2.0*r;
    sphericity = 1.0;
    double I = 0.4*volume*r*r;
    inertia_unit = dem_vec(I,I,I);

    icosphere(2);
    for(auto &v : vert)
    v *= r;

    nodes.clear();   // sphere contacts are analytic
}

void dem_shape::build_box(double lx, double ly, double lz)
{
    type = DEM_BOX;
    dim = dem_vec(0.5*lx,0.5*ly,0.5*lz);
    volume = lx*ly*lz;
    area = 2.0*(lx*ly + ly*lz + lx*lz);
    rbound = dim.norm();
    deq = pow(6.0*volume/dem_pi,1.0/3.0);
    sphericity = pow(dem_pi,1.0/3.0)*pow(6.0*volume,2.0/3.0)/area;
    inertia_unit = dem_vec(volume*(ly*ly+lz*lz)/12.0, volume*(lx*lx+lz*lz)/12.0, volume*(lx*lx+ly*ly)/12.0);

    vert.clear();
    tri.clear();
    for(int n=0; n<8; ++n)
    vert.push_back(dem_vec((n&1?1:-1)*dim(0),(n&2?1:-1)*dim(1),(n&4?1:-1)*dim(2)));

    int t[36] = {0,2,1, 1,2,3, 4,5,6, 5,7,6, 0,1,4, 1,5,4, 2,6,3, 3,6,7, 0,4,2, 2,4,6, 1,3,5, 3,7,5};
    tri.assign(t,t+36);
}

void dem_shape::build_cylinder(double r, double length)
{
    type = DEM_CYLINDER;
    double h = 0.5*length;
    dim = dem_vec(r,r,h);
    volume = dem_pi*r*r*length;
    area = 2.0*dem_pi*r*length + 2.0*dem_pi*r*r;
    rbound = sqrt(r*r+h*h);
    deq = pow(6.0*volume/dem_pi,1.0/3.0);
    sphericity = pow(dem_pi,1.0/3.0)*pow(6.0*volume,2.0/3.0)/area;
    double Ixx = volume*(3.0*r*r + length*length)/12.0;
    inertia_unit = dem_vec(Ixx,Ixx,0.5*volume*r*r);

    int ns=32;
    vert.clear();
    tri.clear();
    for(int n=0; n<ns; ++n)
    {
        double phi = 2.0*dem_pi*double(n)/double(ns);
        vert.push_back(dem_vec(r*cos(phi),r*sin(phi),-h));
        vert.push_back(dem_vec(r*cos(phi),r*sin(phi), h));
    }
    int cb = vert.size(); vert.push_back(dem_vec(0,0,-h));
    int ct = vert.size(); vert.push_back(dem_vec(0,0, h));

    for(int n=0; n<ns; ++n)
    {
        int b0 = 2*n, t0 = 2*n+1, b1 = 2*((n+1)%ns), t1 = 2*((n+1)%ns)+1;
        tri.insert(tri.end(),{b0,b1,t0, t0,b1,t1, cb,b1,b0, ct,t0,t1});
    }
}

void dem_shape::build_ellipsoid(double a, double b, double c, int res)
{
    vert.clear();
    tri.clear();

    int nt=24, np=48;
    for(int it=0; it<=nt; ++it)
    for(int ip=0; ip<np; ++ip)
    {
        double th = dem_pi*double(it)/double(nt);
        double ph = 2.0*dem_pi*double(ip)/double(np);
        vert.push_back(dem_vec(a*sin(th)*cos(ph), b*sin(th)*sin(ph), c*cos(th)));
    }

    for(int it=0; it<nt; ++it)
    for(int ip=0; ip<np; ++ip)
    {
        int v0 = it*np+ip, v1 = it*np+(ip+1)%np, v2 = (it+1)*np+ip, v3 = (it+1)*np+(ip+1)%np;
        if(it>0)    tri.insert(tri.end(),{v0,v2,v1});
        if(it<nt-1) tri.insert(tri.end(),{v1,v2,v3});
    }

    finish_mesh_shape(false);
    type = DEM_ELLIPSOID;
    dim = dem_vec(a,b,c);
    build_sdf_grid(res);
}

// ---------------------------------------------------------------------------------------------
// STL
// ---------------------------------------------------------------------------------------------

bool dem_shape::build_mesh(const string &file, double scale, int res)
{
    ifstream in(file, ios::binary);
    if(!in.is_open())
    {
        cout<<"DEM: cannot open STL file "<<file<<endl;
        return false;
    }

    in.seekg(0,ios::end);
    size_t fsize = in.tellg();
    in.seekg(0,ios::beg);

    vector<dem_vec> raw;

    char header[80];
    in.read(header,80);
    uint32_t ntri=0;
    in.read(reinterpret_cast<char*>(&ntri),4);

    if(in && fsize == 84 + size_t(ntri)*50)
    {
        // binary
        for(uint32_t n=0; n<ntri; ++n)
        {
            float buf[12];
            uint16_t attr;
            in.read(reinterpret_cast<char*>(buf),48);
            in.read(reinterpret_cast<char*>(&attr),2);
            for(int q=0; q<3; ++q)
            raw.push_back(dem_vec(buf[3+3*q],buf[4+3*q],buf[5+3*q])*scale);
        }
    }
    else
    {
        // ascii
        in.clear();
        in.seekg(0,ios::beg);
        string word;
        while(in>>word)
        if(word=="vertex")
        {
            double x,y,z;
            in>>x>>y>>z;
            raw.push_back(dem_vec(x,y,z)*scale);
        }
    }

    if(raw.size()<12 || raw.size()%3!=0)
    {
        cout<<"DEM: STL file "<<file<<" contains no valid triangles"<<endl;
        return false;
    }

    // weld vertices
    dem_vec lo = raw[0], hi = raw[0];
    for(auto &v : raw) {lo = lo.cwiseMin(v); hi = hi.cwiseMax(v);}
    double tolw = 1.0e-7*(hi-lo).norm();
    if(tolw<=0.0)
    {
        cout<<"DEM: STL file "<<file<<" is degenerate"<<endl;
        return false;
    }

    map<array<long long,3>,int> index;
    vert.clear();
    tri.clear();
    for(auto &v : raw)
    {
        array<long long,3> key = {llround(v(0)/tolw),llround(v(1)/tolw),llround(v(2)/tolw)};
        auto it = index.find(key);
        if(it==index.end())
        {
            index[key] = vert.size();
            tri.push_back(vert.size());
            vert.push_back(v);
        }
        else
        tri.push_back(it->second);
    }

    finish_mesh_shape();
    if(!(volume>0.0))
    {
        cout<<"DEM: STL file "<<file<<" does not enclose a volume"<<endl;
        return false;
    }
    type = DEM_MESH;
    build_sdf_grid(res);

    return true;
}

void dem_shape::mesh_mass_properties(dem_vec &com, dem_mat &J, double &vol)
{
    // D. Eberly, Polyhedral Mass Properties (Revisited)
    const double mult[10] = {1.0/6.0,1.0/24.0,1.0/24.0,1.0/24.0,1.0/60.0,1.0/60.0,1.0/60.0,1.0/120.0,1.0/120.0,1.0/120.0};
    double intg[10] = {0,0,0,0,0,0,0,0,0,0};

    auto sub = [](double w0, double w1, double w2, double &f1, double &f2, double &f3, double &g0, double &g1, double &g2)
    {
        double temp0 = w0+w1;
        f1 = temp0+w2;
        double temp1 = w0*w0;
        double temp2 = temp1 + w1*temp0;
        f2 = temp2 + w2*f1;
        f3 = w0*temp1 + w1*temp2 + w2*f2;
        g0 = f2 + w0*(f1+w0);
        g1 = f2 + w1*(f1+w1);
        g2 = f2 + w2*(f1+w2);
    };

    for(size_t t=0; t<tri.size(); t+=3)
    {
        const dem_vec &p0 = vert[tri[t]], &p1 = vert[tri[t+1]], &p2 = vert[tri[t+2]];
        dem_vec d = (p1-p0).cross(p2-p0);

        double f1x,f2x,f3x,g0x,g1x,g2x;
        double f1y,f2y,f3y,g0y,g1y,g2y;
        double f1z,f2z,f3z,g0z,g1z,g2z;
        sub(p0(0),p1(0),p2(0),f1x,f2x,f3x,g0x,g1x,g2x);
        sub(p0(1),p1(1),p2(1),f1y,f2y,f3y,g0y,g1y,g2y);
        sub(p0(2),p1(2),p2(2),f1z,f2z,f3z,g0z,g1z,g2z);

        intg[0] += d(0)*f1x;
        intg[1] += d(0)*f2x;
        intg[2] += d(1)*f2y;
        intg[3] += d(2)*f2z;
        intg[4] += d(0)*f3x;
        intg[5] += d(1)*f3y;
        intg[6] += d(2)*f3z;
        intg[7] += d(0)*(p0(1)*g0x + p1(1)*g1x + p2(1)*g2x);
        intg[8] += d(1)*(p0(2)*g0y + p1(2)*g1y + p2(2)*g2y);
        intg[9] += d(2)*(p0(0)*g0z + p1(0)*g1z + p2(0)*g2z);
    }

    for(int n=0; n<10; ++n)
    intg[n] *= mult[n];

    vol = intg[0];
    com = dem_vec(intg[1],intg[2],intg[3])/vol;

    J(0,0) = intg[5] + intg[6] - vol*(com(1)*com(1) + com(2)*com(2));
    J(1,1) = intg[4] + intg[6] - vol*(com(2)*com(2) + com(0)*com(0));
    J(2,2) = intg[4] + intg[5] - vol*(com(0)*com(0) + com(1)*com(1));
    J(0,1) = J(1,0) = -(intg[7] - vol*com(0)*com(1));
    J(1,2) = J(2,1) = -(intg[8] - vol*com(1)*com(2));
    J(0,2) = J(2,0) = -(intg[9] - vol*com(2)*com(0));
}

void dem_shape::finish_mesh_shape(bool rotate)
{
    dem_vec com;
    dem_mat J;
    double vol;

    mesh_mass_properties(com,J,vol);

    // inward oriented mesh: flip
    if(vol<0.0)
    {
        for(size_t t=0; t<tri.size(); t+=3)
        std::swap(tri[t+1],tri[t+2]);

        mesh_mass_properties(com,J,vol);
    }

    // principal axes (ellipsoids are built in their principal frame, keep the user's axis order)
    if(rotate)
    {
        Eigen::SelfAdjointEigenSolver<dem_mat> es(J);
        dem_mat A = es.eigenvectors();
        if(A.determinant()<0.0)
        A.col(2) *= -1.0;

        for(auto &v : vert)
        v = A.transpose()*(v - com);

        inertia_unit = es.eigenvalues();
    }
    else
    {
        for(auto &v : vert)
        v -= com;

        inertia_unit = J.diagonal();
    }

    volume = vol;

    area = 0.0;
    for(size_t t=0; t<tri.size(); t+=3)
    area += 0.5*((vert[tri[t+1]]-vert[tri[t]]).cross(vert[tri[t+2]]-vert[tri[t]])).norm();

    rbound = 0.0;
    for(auto &v : vert)
    rbound = std::max(rbound,v.norm());

    deq = pow(6.0*volume/dem_pi,1.0/3.0);
    sphericity = std::min(1.0,pow(dem_pi,1.0/3.0)*pow(6.0*volume,2.0/3.0)/area);
}

double dem_shape::mesh_distance(const dem_vec &y, double &winding) const
{
    double dmin = 1.0e20;
    double omega = 0.0;

    for(size_t t=0; t<tri.size(); t+=3)
    {
        const dem_vec &p0 = vert[tri[t]], &p1 = vert[tri[t+1]], &p2 = vert[tri[t+2]];

        double d = (closest_point_triangle(y,p0,p1,p2)-y).squaredNorm();
        dmin = std::min(dmin,d);

        // Van Oosterom & Strackee solid angle
        dem_vec a = p0-y, b = p1-y, c = p2-y;
        double la=a.norm(), lb=b.norm(), lc=c.norm();
        double num = a.dot(b.cross(c));
        double den = la*lb*lc + a.dot(b)*lc + a.dot(c)*lb + b.dot(c)*la;
        omega += 2.0*atan2(num,den);
    }

    winding = omega/(4.0*dem_pi);
    return sqrt(dmin);
}

void dem_shape::build_sdf_grid(int res)
{
    dem_vec lo = vert[0], hi = vert[0];
    for(auto &v : vert) {lo = lo.cwiseMin(v); hi = hi.cwiseMax(v);}

    res = std::max(res,8);
    gdx = (hi-lo).maxCoeff()/double(res);
    g0 = lo - dem_vec::Constant(3.0*gdx);

    for(int q=0; q<3; ++q)
    gn[q] = int(ceil((hi(q)-lo(q))/gdx)) + 7;

    sdfval.assign(size_t(gn[0])*gn[1]*gn[2],0.0f);

    if(tri.size()/3 > 20000)
    cout<<"DEM: warning, "<<tri.size()/3<<" triangles, SDF generation may take long"<<endl;

    for(int k=0; k<gn[2]; ++k)
    for(int j=0; j<gn[1]; ++j)
    for(int i=0; i<gn[0]; ++i)
    {
        dem_vec y = g0 + gdx*dem_vec(i,j,k);
        double w;
        double d = mesh_distance(y,w);
        sdfval[(size_t(k)*gn[1] + j)*gn[0] + i] = float(w>0.5 ? -d : d);
    }
}

double dem_shape::sdf_grid(const dem_vec &y) const
{
    dem_vec s = (y-g0)/gdx;
    dem_vec sc;
    for(int q=0; q<3; ++q)
    sc(q) = std::min(std::max(s(q),0.0),double(gn[q]-1)-1.0e-9);

    int i = int(sc(0)), j = int(sc(1)), k = int(sc(2));
    double fx = sc(0)-i, fy = sc(1)-j, fz = sc(2)-k;

    auto V = [&](int ii, int jj, int kk) {return double(sdfval[(size_t(kk)*gn[1] + jj)*gn[0] + ii]);};

    double c00 = V(i,j,k)*(1-fx)     + V(i+1,j,k)*fx;
    double c10 = V(i,j+1,k)*(1-fx)   + V(i+1,j+1,k)*fx;
    double c01 = V(i,j,k+1)*(1-fx)   + V(i+1,j,k+1)*fx;
    double c11 = V(i,j+1,k+1)*(1-fx) + V(i+1,j+1,k+1)*fx;
    double val = (c00*(1-fy) + c10*fy)*(1-fz) + (c01*(1-fy) + c11*fy)*fz;

    return val + gdx*(s-sc).norm();
}

dem_vec dem_shape::sdf_grid_grad(const dem_vec &y) const
{
    dem_vec s = (y-g0)/gdx;
    dem_vec sc;
    for(int q=0; q<3; ++q)
    sc(q) = std::min(std::max(s(q),0.0),double(gn[q]-1)-1.0e-9);

    if((s-sc).norm()>1.0e-9)
    return (s-sc).normalized();

    int i = int(sc(0)), j = int(sc(1)), k = int(sc(2));
    double fx = sc(0)-i, fy = sc(1)-j, fz = sc(2)-k;

    auto V = [&](int ii, int jj, int kk) {return double(sdfval[(size_t(kk)*gn[1] + jj)*gn[0] + ii]);};

    double gx = ((V(i+1,j,k)-V(i,j,k))*(1-fy) + (V(i+1,j+1,k)-V(i,j+1,k))*fy)*(1-fz)
              + ((V(i+1,j,k+1)-V(i,j,k+1))*(1-fy) + (V(i+1,j+1,k+1)-V(i,j+1,k+1))*fy)*fz;
    double gy = ((V(i,j+1,k)-V(i,j,k))*(1-fx) + (V(i+1,j+1,k)-V(i+1,j,k))*fx)*(1-fz)
              + ((V(i,j+1,k+1)-V(i,j,k+1))*(1-fx) + (V(i+1,j+1,k+1)-V(i+1,j,k+1))*fx)*fz;
    double gz = ((V(i,j,k+1)-V(i,j,k))*(1-fx) + (V(i+1,j,k+1)-V(i+1,j,k))*fx)*(1-fy)
              + ((V(i,j+1,k+1)-V(i,j+1,k))*(1-fx) + (V(i+1,j+1,k+1)-V(i+1,j+1,k))*fx)*fy;

    dem_vec g(gx,gy,gz);
    double nn = g.norm();
    return nn>1.0e-12 ? dem_vec(g/nn) : dem_vec(y.normalized());
}

// ---------------------------------------------------------------------------------------------
// signed distance
// ---------------------------------------------------------------------------------------------

double dem_shape::sdf(const dem_vec &y) const
{
    if(type==DEM_SPHERE)
    return y.norm() - dim(0);

    if(type==DEM_BOX)
    {
        dem_vec qv = y.cwiseAbs() - dim;
        return qv.cwiseMax(0.0).norm() + std::min(qv.maxCoeff(),0.0);
    }

    if(type==DEM_CYLINDER)
    {
        double dr = sqrt(y(0)*y(0)+y(1)*y(1)) - dim(0);
        double dz = fabs(y(2)) - dim(2);
        return std::min(std::max(dr,dz),0.0) + sqrt(std::max(dr,0.0)*std::max(dr,0.0) + std::max(dz,0.0)*std::max(dz,0.0));
    }

    return sdf_grid(y);
}

dem_vec dem_shape::sdf_grad(const dem_vec &y) const
{
    if(type==DEM_SPHERE)
    {
        double nn = y.norm();
        return nn>1.0e-14 ? dem_vec(y/nn) : dem_vec(dem_vec::UnitZ());
    }

    if(type==DEM_BOX || type==DEM_CYLINDER)
    {
        double h = 1.0e-6*rbound;
        dem_vec g;
        for(int q=0; q<3; ++q)
        {
            dem_vec e = dem_vec::Zero();
            e(q) = h;
            g(q) = sdf(y+e) - sdf(y-e);
        }
        double nn = g.norm();
        return nn>1.0e-14 ? dem_vec(g/nn) : dem_vec(dem_vec::UnitZ());
    }

    return sdf_grid_grad(y);
}

// ---------------------------------------------------------------------------------------------
// surface nodes and volume quadrature
// ---------------------------------------------------------------------------------------------

void dem_shape::make_nodes(double spacing)
{
    nodes.clear();

    if(type==DEM_SPHERE)
    return;

    spacing = std::max(spacing,1.0e-3*rbound);

    if(type==DEM_BOX)
    {
        // regular samples on every face, including edges and corners
        int nx = std::max(1,int(ceil(2.0*dim(0)/spacing)));
        int ny = std::max(1,int(ceil(2.0*dim(1)/spacing)));
        int nz = std::max(1,int(ceil(2.0*dim(2)/spacing)));

        for(int i=0; i<=nx; ++i)
        for(int j=0; j<=ny; ++j)
        for(int k=0; k<=nz; ++k)
        if(i==0 || i==nx || j==0 || j==ny || k==0 || k==nz)
        nodes.push_back(dem_vec(-dim(0)+2.0*dim(0)*i/nx, -dim(1)+2.0*dim(1)*j/ny, -dim(2)+2.0*dim(2)*k/nz));

        return;
    }

    if(type==DEM_CYLINDER)
    {
        double r = dim(0), h = dim(2);
        int nphi = std::max(8,int(ceil(2.0*3.14159265358979*r/spacing)));
        int nz = std::max(1,int(ceil(2.0*h/spacing)));
        int nr = std::max(1,int(ceil(r/spacing)));

        for(int k=0; k<=nz; ++k)
        for(int m=0; m<nphi; ++m)
        {
            double phi = 2.0*3.14159265358979*m/nphi;
            nodes.push_back(dem_vec(r*cos(phi),r*sin(phi),-h+2.0*h*k/nz));
        }

        for(int s=-1; s<=1; s+=2)
        {
            nodes.push_back(dem_vec(0,0,s*h));
            for(int ir=1; ir<nr; ++ir)
            {
                double rr = r*double(ir)/nr;
                int np = std::max(6,int(ceil(2.0*3.14159265358979*rr/spacing)));
                for(int m=0; m<np; ++m)
                {
                    double phi = 2.0*3.14159265358979*m/np;
                    nodes.push_back(dem_vec(rr*cos(phi),rr*sin(phi),s*h));
                }
            }
        }
        return;
    }

    // triangulated shapes: vertices, edge samples, face samples,
    // thinned to one node per voxel of half the node spacing (vertices first, so corners are kept)
    map<array<long long,3>,int> voxel;
    double hv = 0.5*spacing;
    auto addnode = [&](const dem_vec &y)
    {
        array<long long,3> key = {llround(y(0)/hv),llround(y(1)/hv),llround(y(2)/hv)};
        if(voxel.find(key)==voxel.end())
        {
            voxel[key]=1;
            nodes.push_back(y);
        }
    };

    for(auto &v : vert)
    addnode(v);

    map<pair<int,int>,int> edges;
    for(size_t t=0; t<tri.size(); t+=3)
    for(int e=0; e<3; ++e)
    {
        int a = tri[t+e], b = tri[t+(e+1)%3];
        edges[make_pair(std::min(a,b),std::max(a,b))] = 1;
    }

    for(auto &ed : edges)
    {
        const dem_vec &a = vert[ed.first.first], &b = vert[ed.first.second];
        int ns = int(floor((b-a).norm()/spacing));
        for(int s=1; s<ns; ++s)
        addnode(a + (b-a)*double(s)/double(ns));
    }

    for(size_t t=0; t<tri.size(); t+=3)
    {
        const dem_vec &a = vert[tri[t]], &b = vert[tri[t+1]], &c = vert[tri[t+2]];
        double lmax = std::max((b-a).norm(),std::max((c-b).norm(),(a-c).norm()));
        int ns = int(floor(lmax/spacing));
        for(int s1=1; s1<ns; ++s1)
        for(int s2=1; s1+s2<ns; ++s2)
        addnode(a + (b-a)*double(s1)/ns + (c-a)*double(s2)/ns);
    }
}

void dem_shape::make_quadrature(int n)
{
    qp.clear();
    qw.clear();

    n = std::max(n,2);
    dem_vec ext;
    if(type==DEM_SPHERE || type==DEM_BOX || type==DEM_CYLINDER)
    ext = dim;
    else
    {
        ext = dem_vec::Zero();
        for(auto &v : vert)
        ext = ext.cwiseMax(v.cwiseAbs());
    }

    double h = 2.0*ext.maxCoeff()/double(n);
    int nq[3];
    for(int q=0; q<3; ++q)
    nq[q] = std::max(1,int(ceil(2.0*ext(q)/h)));

    for(int i=0; i<nq[0]; ++i)
    for(int j=0; j<nq[1]; ++j)
    for(int k=0; k<nq[2]; ++k)
    {
        dem_vec y(-ext(0) + (i+0.5)*2.0*ext(0)/nq[0], -ext(1) + (j+0.5)*2.0*ext(1)/nq[1], -ext(2) + (k+0.5)*2.0*ext(2)/nq[2]);
        if(sdf(y)<0.0)
        qp.push_back(y);
    }

    if(qp.empty())
    qp.push_back(dem_vec::Zero());

    // centre the quadrature on the centroid (no spurious buoyancy torque for asymmetric shapes)
    dem_vec mean = dem_vec::Zero();
    for(auto &y : qp)
    mean += y;
    mean /= double(qp.size());
    for(auto &y : qp)
    y -= mean;

    qw.assign(qp.size(),volume/double(qp.size()));
}

void dem_shape::icosphere(int levels)
{
    double t = (1.0+sqrt(5.0))/2.0;
    vert = {dem_vec(-1,t,0),dem_vec(1,t,0),dem_vec(-1,-t,0),dem_vec(1,-t,0),
            dem_vec(0,-1,t),dem_vec(0,1,t),dem_vec(0,-1,-t),dem_vec(0,1,-t),
            dem_vec(t,0,-1),dem_vec(t,0,1),dem_vec(-t,0,-1),dem_vec(-t,0,1)};
    for(auto &v : vert) v.normalize();

    tri = {0,11,5, 0,5,1, 0,1,7, 0,7,10, 0,10,11, 1,5,9, 5,11,4, 11,10,2, 10,7,6, 7,1,8,
           3,9,4, 3,4,2, 3,2,6, 3,6,8, 3,8,9, 4,9,5, 2,4,11, 6,2,10, 8,6,7, 9,8,1};

    for(int l=0; l<levels; ++l)
    {
        map<pair<int,int>,int> mid;
        auto midpoint = [&](int a, int b)
        {
            auto key = make_pair(std::min(a,b),std::max(a,b));
            auto it = mid.find(key);
            if(it!=mid.end()) return it->second;
            int id = vert.size();
            vert.push_back(((vert[a]+vert[b])*0.5).normalized());
            mid[key]=id;
            return id;
        };

        vector<int> nt;
        for(size_t q=0; q<tri.size(); q+=3)
        {
            int a=tri[q], b=tri[q+1], c=tri[q+2];
            int ab=midpoint(a,b), bc=midpoint(b,c), ca=midpoint(c,a);
            nt.insert(nt.end(),{a,ab,ca, b,bc,ab, c,ca,bc, ab,bc,ca});
        }
        tri = nt;
    }
}
