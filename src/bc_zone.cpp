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
#include<iostream>
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
    // B 520-524 (control.h, read by read_control and sent to all ranks)
    auto fail = [&](const std::string &msg)
    {
        if(p->mpirank==0)
        std::cout<<std::endl<<"!!! bc_zone: "<<msg<<" !!!"<<std::endl<<std::endl;
        std::exit(1);
    };
    
    const int n = p->B520;
    
    if(n==0)
    {
        if(p->B521>0 || p->B524>0 || p->B523>0)
        fail("B 521 / B 523 / B 524 refer to zones, but no zone is defined by B 520");
        return;
    }
    
    std::vector<int> edge(n,0);
    std::vector<double> s0(n,0.0), s1(n,0.0), width(n,0.0);
    std::vector<std::vector<int>> src(n);
    
    auto find = [&](int k)->int
    {
        for(int q=0; q<n; ++q)
        if(p->B520_id[q]==k)
        return q;
        return -1;
    };
    
    for(int q=0; q<n; ++q)
    {
        if(find(p->B520_id[q])!=q)
        fail("zone id "+std::to_string(p->B520_id[q])+" defined twice in B 520");
        if(p->B520_method[q]<1 || p->B520_method[q]>4)
        fail("B 520 method is 1 (relaxation), 2 (beach), 3 (Riemann edge) or 4 (Flather edge)");
    }
    
    for(int m=0; m<p->B521; ++m)
    {
        int q = find(p->B521_id[m]);
        if(q<0)
        fail("B 521 refers to zone "+std::to_string(p->B521_id[m])+", which has no B 520");
        
        edge[q]  = p->B521_edge[m];
        s0[q]    = p->B521_s0[m];
        s1[q]    = p->B521_s1[m];
        width[q] = p->B521_w[m];
    }
    
    for(int m=0; m<p->B524; ++m)
    {
        int q = find(p->B524_id[m]);
        if(q<0)
        fail("B 524 refers to zone "+std::to_string(p->B524_id[m])+", which has no B 520");
        
        src[q].push_back(p->B524_src[m]);
    }
    
    std::vector<int> bgid(n,0);
    
    for(int m=0; m<p->B523; ++m)
    {
        int q = find(p->B523_id[m]);
        if(q<0)
        fail("B 523 refers to zone "+std::to_string(p->B523_id[m])+", which has no B 520");
        if(p->B523_bg[m]<1)
        fail("B 523: background ids start at 1");
        
        bgid[q] = p->B523_bg[m];
    }
    
    const double fac = (p->B99==1) ? 2.0 : 1.0;
    const double ext = 10.0*p->DXM;
    
    for(int q=0; q<n; ++q)
    {
        const int k=p->B520_id[q], m=p->B520_method[q], e=edge[q];
        
        if(e<1 || e>4)
        fail("zone "+std::to_string(k)+" needs a B 521 edge (1: x-, 2: x+, 3: y-, 4: y+)");
        if(width[q]<=0.0 && m<=2)
        fail("zone "+std::to_string(k)+": the B 521 width must be positive");
        if(m>=3 && bgid[q]==0)
        fail("zone "+std::to_string(k)+": a Riemann or Flather edge needs a background (B 523)");
        
        const double a=s0[q], b=s1[q], w = m<=2 ? width[q] : p->DXM;
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
        
        const bc_method meth = m==1 ? bc_method::relax : m==2 ? bc_method::beach : m==3 ? bc_method::riemann : bc_method::flather;
        
        bc_zone z(k, meth, xs,ys,xe,ye,w, m==2 ? fac : 1.0);
        z.priority = p->B520_prio[q];
        z.user = true;
        z.sources = src[q];
        z.bg = bgid[q];
        z.edge = e;
        
        if(m==1)
        relax.push_back(z);
        
        if(m==2)
        beach.push_back(z);
        
        if(m>=3)
        {
            if(open_edge(e)!=nullptr)
            fail("zone "+std::to_string(k)+": edge "+std::to_string(e)+" has two Riemann / Flather zones");
            
            edges.push_back(z);
        }
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

bool bc_zone_set::has_background() const
{
    for(const bc_zone &z : relax)
    if(z.bg>0)
    return true;

    for(const bc_zone &z : beach)
    if(z.bg>0)
    return true;

    return !edges.empty();
}

const bc_zone* bc_zone_set::beach_zone_at(double x0, double y0) const
{
    const bc_zone *best=nullptr;

    for(const bc_zone &z : beach)
    if(z.inside(x0,y0) && (best==nullptr || z.priority>best->priority))
    best=&z;

    return best;
}

const bc_zone* bc_zone_set::open_edge(int e) const
{
    for(const bc_zone &z : edges)
    if(z.edge==e)
    return &z;

    return nullptr;
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
