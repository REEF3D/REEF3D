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

#include"bc_zone.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<cmath>

// ---------------------------------------------------------------------
// bc_zone
// ---------------------------------------------------------------------

bc_zone::bc_zone(int iid, bc_method m, double x_s, double y_s, double x_e, double y_e, double dd, double f)
                : id(iid), method(m), xs(x_s), ys(y_s), xe(x_e), ye(y_e), d(dd), fac(f)
{
    // quadrilateral as iowave::distbeach_ini / distgen_ini
    double dx,dy,l1,l2,n1x,n1y,n2x,n2y;

    dx = xe-xs;
    dy = ye-ys;

    n1x = -dy;
    n1y =  dx;

    n2x =  dy;
    n2y = -dx;

    l1 = sqrt(pow(n1x,2.0) + pow(n1y,2.0));
    l2 = sqrt(pow(n2x,2.0) + pow(n2y,2.0));

    l1 = fabs(l1)>1.0e-20?l1:1.0e20;
    l2 = fabs(l2)>1.0e-20?l2:1.0e20;

    P1[0] = xs + (n1x/l1)*d;
    P1[1] = ys + (n1y/l1)*d;

    P2[0] = xs + (n2x/l2)*d;
    P2[1] = ys + (n2y/l2)*d;

    P3[0] = xe + (n1x/l1)*d;
    P3[1] = ye + (n1y/l1)*d;

    P4[0] = xe + (n2x/l2)*d;
    P4[1] = ye + (n2y/l2)*d;
}

bool bc_zone::inside(double x0, double y0) const
{
    return intriangle(P1[0],P1[1],P3[0],P3[1],P2[0],P2[1],x0,y0)==1
        || intriangle(P3[0],P3[1],P4[0],P4[1],P2[0],P2[1],x0,y0)==1;
}

double bc_zone::line_dist(double x0, double y0) const
{
    double denom = sqrt(pow(ye-ys,2.0) + pow(xe-xs,2.0));
    denom = denom>1.0e-20?denom:1.0e20;

    return fabs((ye-ys)*x0 - (xe-xs)*y0 + xe*ys - ye*xs)/denom;
}

void bc_zone::box(double &bxs, double &bxe, double &bys, double &bye) const
{
    bxs = fmin(fmin(P1[0],P2[0]),fmin(P3[0],P4[0]));
    bxe = fmax(fmax(P1[0],P2[0]),fmax(P3[0],P4[0]));
    bys = fmin(fmin(P1[1],P2[1]),fmin(P3[1],P4[1]));
    bye = fmax(fmax(P1[1],P2[1]),fmax(P3[1],P4[1]));
}

// exact copy of iowave::intriangle
int bc_zone::intriangle(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double x0,double y0)
{
	double Px,Py,Pz;
	double Qx,Qy,Qz;
	double PQx,PQy,PQz;
	double Mx,My;
	double u,v,w;

    Px = x0;
    Py = y0;
    Pz = -10.0;

    Qx = x0;
    Qy = y0;
    Qz = +10.0;

    PQx = Qx-Px;
    PQy = Qy-Py;
    PQz = Qz-Pz;

    Mx = PQy*Pz - PQz*Py;
    My = PQz*Px - PQx*Pz;

    u = PQz*(Cx*By - Cy*Bx)
      + Mx*(Cx-Bx) + My*(Cy-By);

    v = PQz*(Ax*Cy - Ay*Cx)
      + Mx*(Ax-Cx) + My*(Ay-Cy);

    w = PQz*(Bx*Ay - By*Ax)
      + Mx*(Bx-Ax) + My*(By-Ay);

    int check=1;
    if(u==0.0 && v==0.0 && w==0.0)
    check = 0;

    if((u>=0.0 && v>=0.0 && w>=0.0) || (u<0.0 && v<0.0 && w<0.0) && check==1)
    return 1;

    else
    return 0;
}

// ---------------------------------------------------------------------
// bc_zone_set
// ---------------------------------------------------------------------

// Translator: iowave has completed the B 107 / B 108 lists from B 96 when
// none were given (constructor), so the lists are always the zones.
bc_zone_set bc_zone_set::from_legacy(lexer *p, ghostcell *pgc)
{
    bc_zone_set zs;

    const double fac = (p->B99==1) ? 2.0 : 1.0;

    for(int n=0; n<p->B108; ++n)
    zs.relax.emplace_back(n+1,bc_method::relax,p->B108_xs[n],p->B108_ys[n],p->B108_xe[n],p->B108_ye[n],p->B108_d[n],1.0);

    for(int n=0; n<p->B107; ++n)
    zs.beach.emplace_back(n+1,bc_method::beach,p->B107_xs[n],p->B107_ys[n],p->B107_xe[n],p->B107_ye[n],p->B107_d[n],fac);

    zs.read_input(p,pgc);

    return zs;
}

void bc_zone_set::read_input(lexer *p, ghostcell *pgc)
{
    // per zone: id method priority edge nsrc src... | s0 s1 width
    std::vector<int> iv;
    std::vector<double> dv;
    int n=0, ni=0;

    if(p->mpirank==0)
    {
        std::vector<int> id,method,prio,edge;
        std::vector<double> s0,s1,width;
        std::vector<std::vector<int>> src;

        auto find = [&](int k)->int
        {
            for(size_t q=0; q<id.size(); ++q)
            if(id[q]==k)
            return (int)q;
            return -1;
        };

        auto fail = [&](const std::string &msg)
        {
            std::cout<<std::endl<<"!!! bc_zone: "<<msg<<" !!!"<<std::endl<<std::endl;
            std::exit(1);
        };

        std::ifstream f("ctrl.txt");
        std::string line;

        for(int pass=0; pass<2; ++pass)
        {
            f.clear();
            f.seekg(0);

            while(std::getline(f,line))
            {
                std::istringstream ls(line);
                std::string c;
                int key,k;

                if(!(ls>>c) || c!="B" || !(ls>>key))
                continue;

                if(pass==0 && key==520)
                {
                    int m,pr;
                    if(!(ls>>k>>m>>pr))
                    fail("B 520 needs: id method priority");
                    if(find(k)>=0)
                    fail("zone id "+std::to_string(k)+" defined twice in B 520");
                    if(m!=1 && m!=2)
                    fail("B 520 method is 1 (relaxation) or 2 (beach)");

                    id.push_back(k); method.push_back(m); prio.push_back(pr); edge.push_back(0);
                    s0.push_back(0.0); s1.push_back(0.0); width.push_back(0.0);
                    src.emplace_back();
                }

                if(pass==1 && (key==521 || key==524))
                {
                    if(!(ls>>k))
                    fail("B "+std::to_string(key)+" needs a zone id");
                    int q = find(k);
                    if(q<0)
                    fail("B "+std::to_string(key)+" refers to zone "+std::to_string(k)+", which has no B 520");

                    if(key==521 && !(ls>>edge[q]>>s0[q]>>s1[q]>>width[q]))
                    fail("B 521 needs: id edge s0 s1 width");

                    if(key==524)
                    {
                        int s;
                        if(!(ls>>s))
                        fail("B 524 needs: id source");
                        src[q].push_back(s);
                    }
                }
            }
        }

        n = (int)id.size();

        for(int q=0; q<n; ++q)
        {
            if(edge[q]<1 || edge[q]>4)
            fail("zone "+std::to_string(id[q])+" needs a B 521 edge (1: x-, 2: x+, 3: y-, 4: y+)");
            if(width[q]<=0.0)
            fail("zone "+std::to_string(id[q])+": the B 521 width must be positive");

            iv.push_back(id[q]); iv.push_back(method[q]); iv.push_back(prio[q]); iv.push_back(edge[q]);
            iv.push_back((int)src[q].size());
            for(int s : src[q])
            iv.push_back(s);

            dv.push_back(s0[q]); dv.push_back(s1[q]); dv.push_back(width[q]);
        }

        ni = (int)iv.size();
    }

    pgc->bcast_int(&n,1);

    if(n==0)
    return;

    pgc->bcast_int(&ni,1);
    iv.resize(ni);
    dv.resize(3*n);
    pgc->bcast_int(iv.data(),ni);
    pgc->bcast_double(dv.data(),3*n);

    const double fac = (p->B99==1) ? 2.0 : 1.0;
    const double ext = 10.0*p->DXM;
    int pos=0;

    for(int q=0; q<n; ++q)
    {
        const int k=iv[pos], m=iv[pos+1], pr=iv[pos+2], e=iv[pos+3], ns=iv[pos+4];
        pos+=5;

        const double a=dv[3*q], b=dv[3*q+1], w=dv[3*q+2];
        const bool whole = b<=a;
        double xs,ys,xe,ye;

        if(e==1 || e==2)
        {
            xs = xe = (e==1) ? p->xcoormin : p->xcoormax;
            ys = whole ? p->ycoormin-ext : p->ycoormin+a;
            ye = whole ? p->ycoormax+ext : p->ycoormin+b;
        }
        else
        {
            ys = ye = (e==3) ? p->ycoormin : p->ycoormax;
            xs = whole ? p->xcoormin-ext : p->xcoormin+a;
            xe = whole ? p->xcoormax+ext : p->xcoormin+b;
        }

        bc_zone z(k, m==1 ? bc_method::relax : bc_method::beach, xs,ys,xe,ye,w, m==1 ? 1.0 : fac);
        z.priority = pr;
        z.user = true;

        for(int s=0; s<ns; ++s)
        z.sources.push_back(iv[pos+s]);
        pos+=ns;

        if(m==1)
        relax.push_back(z);
        else
        beach.push_back(z);
    }
}

const bc_zone* bc_zone_set::relax_zone_at(double x0, double y0) const
{
    const bc_zone *best=nullptr;

    for(const bc_zone &z : relax)
    if(z.inside(x0,y0) && (best==nullptr || z.priority>best->priority))
    best=&z;

    return best;
}

bool bc_zone_set::has_sources() const
{
    for(const bc_zone &z : relax)
    if(!z.sources.empty())
    return true;

    for(const bc_zone &z : beach)
    if(!z.sources.empty())
    return true;

    return false;
}

bool bc_zone_set::user_relax() const
{
    for(const bc_zone &z : relax)
    if(z.user)
    return true;

    return false;
}

bool bc_zone_set::user_beach() const
{
    for(const bc_zone &z : beach)
    if(z.user)
    return true;

    return false;
}

// former iowave::rb1_ext
double bc_zone_set::relax_weight(double x0, double y0) const
{
    double xdist, r=0.0, x=-1.0e20;
    const double dist=1.0e20;
    int test_all=0;

    for(const bc_zone &z : relax)
    {
        if(z.inside(x0,y0))
        {
        test_all=1;

        xdist = MIN(z.line_dist(x0,y0),dist);

        x = MAX(1.0-xdist/z.d,x);
        x = MAX(x,0.0);
        }

        if(test_all==1)
        r = 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);
    }

    if(test_all==0)
    r=1.0;

    return r;
}

// former iowave::rb1_flag
int bc_zone_set::relax_flag(double x0, double y0) const
{
    int flag=0;

    for(const bc_zone &z : relax)
    if(z.inside(x0,y0))
    flag=1;

    return flag;
}

// former iowave::rb3_ext
double bc_zone_set::beach_weight(double x0, double y0) const
{
    double x, r=0.0, dist2;
    const double dist=1.0e20;
    int test_all=0, count=0;

    for(const bc_zone &z : beach)
    if(z.inside(x0,y0))
    {
        test_all=1;

        x = MIN(z.line_dist(x0,y0),dist);

        dist2 = z.d;
        x=(dist2-fabs(x))/(dist2*z.fac);
        x=MAX(x,0.0);

        r += 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);
        ++count;
    }

    if(test_all==0)
    r=1.0;

    if(test_all==1)
    r/=double(count);

	return r;
}

// former iowave::distgen_calc
double bc_zone_set::relax_dist(double x0, double y0) const
{
    double dist=1.0e20;

    for(const bc_zone &z : relax)
    if(z.inside(x0,y0))
    dist = MIN(z.line_dist(x0,y0),dist);

    return dist;
}

// former iowave::distbeach_calc
double bc_zone_set::beach_dist(double x0, double y0) const
{
    double dist=1.0e20;

    for(const bc_zone &z : beach)
    if(z.inside(x0,y0))
    dist = MIN(z.line_dist(x0,y0),dist);

    return dist;
}

void bc_zone_set::norefine_boxes(lexer *p, std::vector<double> &fbox) const
{
    // old input: exactly the former exclusion of nhflow_amr / fnpf_amr, the B 96
    // ranges from the ends of the domain in x (zones given by B 107 / B 108 were
    // not excluded, and stay so, so that existing cases refine as before)
    const double big = 1.0e20;
    
    if(p->B98==2 && p->B96_1>0.0)
    {
        fbox.push_back(-big); fbox.push_back(p->global_xmin+p->B96_1);
        fbox.push_back(-big); fbox.push_back(big);
    }
    
    if((p->B99==1 || p->B99==2) && p->B96_2>0.0)
    {
        fbox.push_back(p->global_xmax-p->B96_2); fbox.push_back(big);
        fbox.push_back(-big); fbox.push_back(big);
    }
    
    // zones given by B 520: their box
    double bxs,bxe,bys,bye;
    
    for(const bc_zone &z : relax)
    if(z.user)
    {
        z.box(bxs,bxe,bys,bye);
        fbox.push_back(bxs); fbox.push_back(bxe);
        fbox.push_back(bys); fbox.push_back(bye);
    }
    
    for(const bc_zone &z : beach)
    if(z.user)
    {
        z.box(bxs,bxe,bys,bye);
        fbox.push_back(bxs); fbox.push_back(bxe);
        fbox.push_back(bys); fbox.push_back(bye);
    }
}
