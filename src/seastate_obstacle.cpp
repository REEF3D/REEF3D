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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"seastate_obstacle.h"
#include"seastate_structure.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"seastate_dispersion.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice4.h"
#include"increment.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>

seastate_obstacle::seastate_obstacle(lexer *p, const seastate_grid *g_) : g(g_), nobs(p->A722), nsig(g_->nsig),
    coast(p->A723_kr>0.0), kr2c(p->A723_kr*p->A723_kr), radius(std::min(std::max(p->A724,0),3))
{
    for(int n=0; n<nobs; ++n)
    {
    xs.push_back(p->A722_xs[n]);
    ys.push_back(p->A722_ys[n]);
    xe.push_back(p->A722_xe[n]);
    ye.push_back(p->A722_ye[n]);
    kt.push_back(p->A722_kt[n]);
    kr.push_back(p->A722_kr[n]);
    zc.push_back(p->A722_zc[n]);
    }

    type.assign(nobs,0);
    pa.assign(nobs,0.0);
    pb.assign(nobs,0.0);
    pc.assign(nobs,0.0);
    tkt.assign(nobs,std::vector<float>());
    tkr.assign(nobs,std::vector<float>());
    wobs.assign(nobs,-1);

    // structures (A 725 n type a b c)
    for(int k=0; k<p->A725; ++k)
    {
    const int n = int(std::lround(p->A725_n[k]))-1;
    if(n<0 || n>=nobs)
    continue;

    type[n] = int(std::lround(p->A725_t[k]));
    pa[n] = p->A725_a[k];
    pb[n] = p->A725_b[k];
    pc[n] = p->A725_c[k];

        // type 2: the table of seastate-obstacle-n.dat on the frequencies of the grid (read by every rank)
        if(type[n]==2)
        {
        std::vector<double> tf, ta, tb;
        std::ifstream in("seastate-obstacle-"+std::to_string(n+1)+".dat");
        std::string line;

            while(std::getline(in,line))
            {
            const size_t c = line.find('$');
            if(c!=std::string::npos)
            line = line.substr(0,c);
            std::istringstream ls(line);
            double f, a, b;
            if(ls>>f>>a>>b)
            {
            tf.push_back(f);
            ta.push_back(a);
            tb.push_back(b);
            }
            }

        tkt[n].assign(nsig,1.0f);
        tkr[n].assign(nsig,0.0f);

            for(int l=0; l<nsig && !tf.empty(); ++l)
            {
            const double f = g->f[l];
            double a, b;
                if(f<=tf.front())
                {a = ta.front(); b = tb.front();}
                else if(f>=tf.back())
                {a = ta.back(); b = tb.back();}
                else
                {
                const size_t u = size_t(std::upper_bound(tf.begin(),tf.end(),f) - tf.begin());
                const double w = (f-tf[u-1])/std::max(tf[u]-tf[u-1],1.0e-30);
                a = (1.0-w)*ta[u-1] + w*ta[u];
                b = (1.0-w)*tb[u-1] + w*tb[u];
                }
            a = std::min(std::max(a,0.0),1.0);
            b = std::min(std::max(b,0.0),1.0);
            tkt[n][l] = float(a*a);
            tkr[n][l] = float(std::min(b*b,1.0-a*a));
            }
        }
    }

    // diffuse reflection (A 726 n pown, A 723 pown)
    for(int k=0; k<p->A726; ++k)
    {
    const int n = int(std::lround(p->A726_n[k]))-1;
    if(n>=0 && n<nobs && p->A726_p[k]>0.0)
    wobs[n] = weights(p->A726_p[k]);
    }

    if(coast && p->A723_pown>0.0)
    wcoast = weights(p->A723_pown);
}

// cos^pown filter over the directions around the specular one, cut below 0.01 and normalised (as SWAN
// REFLECT, RDIFF); index nd+k for the offset k
int seastate_obstacle::weights(double pown)
{
    const int ndir = g->ndir;
    std::vector<double> w(1,1.0);

    for(int k=1; k<=ndir/2; ++k)
    {
    const double c = std::cos(k*g->dtheta);
    const double v = (c>0.0) ? std::pow(c,pown) : 0.0;
    if(!(v>0.01))
    break;
    w.push_back(v);
    }

    const int nd = int(w.size())-1;
    double s = w[0];
    for(int k=1; k<=nd; ++k)
    s += 2.0*w[k];

    std::vector<float> r(2*nd+1);
    for(int k=-nd; k<=nd; ++k)
    r[nd+k] = float(w[std::abs(k)]/s);

    wdif.push_back(r);
    return int(wdif.size())-1;
}

// the segments a-b and c-d cross (proper intersection or touching)
static bool cross(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy)
{
    auto orient = [](double px, double py, double qx, double qy, double rx, double ry)
    {
        const double v = (qx-px)*(ry-py) - (qy-py)*(rx-px);
        return (v>0.0) - (v<0.0);
    };

    const int o1 = orient(ax,ay,bx,by,cx,cy), o2 = orient(ax,ay,bx,by,dx,dy);
    const int o3 = orient(cx,cy,dx,dy,ax,ay), o4 = orient(cx,cy,dx,dy,bx,by);

    return o1*o2<=0 && o3*o4<=0 && !(o1==0 && o2==0);
}

void seastate_obstacle::build(lexer *q, fdm_seastate *e, int oi_, int oj_, int gnx_, int gny_, ghostcell *pgc_)
{
    pgc = pgc_;
    // cell centres: XP[i + marge] (the coordinate arrays have marge ghost entries, the slices margin)
    const int m = increment::marge;

    imin = q->imin;
    jmin = q->jmin;
    ni = q->imax;
    nj = q->jmax;
    oi = oi_;
    oj = oj_;
    gnx = gnx_;
    gny = gny_;

    fe.assign(size_t(ni)*nj,-1);
    fn.assign(size_t(ni)*nj,-1);
    faces.clear();
    fi.clear();
    fj.clear();
    fdir.clear();
    fkt.clear();
    fkr.clear();
    ncoast = 0;

    auto add = [&](int i, int j, int s, face f)
    {
        const size_t c = size_t(i-imin)*nj + (j-jmin);
        (s==0 ? fe : fn)[c] = int(faces.size());
        faces.push_back(f);
        fi.push_back(i);
        fj.push_back(j);
        fdir.push_back(s);
    };

    // obstacles
    if(nobs>0)
    for(int i=imin; i<imin+ni; ++i)
    for(int j=jmin; j<jmin+nj; ++j)
    for(int s=0; s<2; ++s)
    {
    const int i2 = (s==0) ? i+1 : i, j2 = (s==0) ? j : j+1;

        if(i2>=imin+ni || j2>=jmin+nj)
        continue;

    const double x1 = q->XP[i+m], y1 = q->YP[j+m];
    const double x2 = q->XP[i2+m], y2 = q->YP[j2+m];

        for(int n=0; n<nobs; ++n)
        if(cross(x1,y1,x2,y2,xs[n],ys[n],xe[n],ye[n]))
        {
        face f;
        f.obs = n;
        f.alpha = float(std::atan2(ye[n]-ys[n],xe[n]-xs[n]));
        f.kt2 = (kt[n]>=0.0 && type[n]==0) ? float(kt[n]*kt[n]) : 1.0f;
        f.kr2 = float(std::min(kr[n]*kr[n],1.0-double(f.kt2)));
        f.wd = wobs[n];

            if(type[n]==2 || type[n]==3)
            {
            f.fq = int(fkt.size()/std::max(nsig,1));
                if(type[n]==2)
                {
                fkt.insert(fkt.end(),tkt[n].begin(),tkt[n].end());
                fkr.insert(fkr.end(),tkr[n].begin(),tkr[n].end());
                }
                else
                {
                fkt.insert(fkt.end(),nsig,1.0f);
                fkr.insert(fkr.end(),nsig,0.0f);
                }
            }

        add(i,j,s,f);
        break;
        }
    }

    // coasts: faces between an active cell and a land or dry cell inside the domain
    wmask.clear();
    wraw.clear();
    if(!coast || e==nullptr)
    return;

    wmask.resize(size_t(ni)*nj);
    wraw.resize(size_t(ni)*nj);
    for(int i=imin; i<imin+ni; ++i)
    for(int j=jmin; j<jmin+nj; ++j)
    wraw[size_t(i-imin)*nj + (j-jmin)] = wmask[size_t(i-imin)*nj + (j-jmin)] = char(e->wet(i,j)==1);

    // level 0: the active cells of all ghost layers from the neighbouring ranks (the coastline uses 3 layers)
    if(pgc!=nullptr)
    {
    slice4 W(q);
        for(int i=imin; i<imin+ni; ++i)
        for(int j=jmin; j<jmin+nj; ++j)
        W(i,j) = double(wmask[size_t(i-imin)*nj + (j-jmin)]);
    pgc->gcsl_start4(q,W,50);
        for(int i=imin; i<imin+ni; ++i)
        for(int j=jmin; j<jmin+nj; ++j)
        wmask[size_t(i-imin)*nj + (j-jmin)] = char(W(i,j)>0.5);
    }

    auto inside = [&](int i, int j) {return i+oi>=0 && i+oi<gnx && j+oj>=0 && j+oj<gny;};
    auto act = [&](int i, int j) {return wmask[size_t(i-imin)*nj + (j-jmin)]!=0;};

    // the coastline through active cell (i,j): the principal direction (total least squares) of the midpoints of the
    // faces between active and land cells within A 724 cells (s: the side of the face, for the fallback)
    auto land = [&](int a, int b) {return inside(a,b) && !act(a,b);};
    auto direction = [&](int i, int j, double &alpha)
    {
        const int a0 = std::max(i-radius,imin), a1 = std::min(i+radius,imin+ni-1);
        const int b0 = std::max(j-radius,jmin), b1 = std::min(j+radius,jmin+nj-1);
        std::vector<double> px, py;

        for(int a=a0; a<=a1; ++a)
        for(int b=b0; b<=b1; ++b)
        {
            if(a+1<=a1 && ((land(a,b) && act(a+1,b)) || (act(a,b) && land(a+1,b))))
            {
            px.push_back(0.5*(q->XP[a+m]+q->XP[a+1+m]));
            py.push_back(q->YP[b+m]);
            }
            if(b+1<=b1 && ((land(a,b) && act(a,b+1)) || (act(a,b) && land(a,b+1))))
            {
            px.push_back(q->XP[a+m]);
            py.push_back(0.5*(q->YP[b+m]+q->YP[b+1+m]));
            }
        }

        double mx = 0.0, my = 0.0;
        for(size_t k=0; k<px.size(); ++k) {mx += px[k]; my += py[k];}
        const double np_ = double(px.size());
        if(np_>0.0) {mx /= np_; my /= np_;}

        double sxx = 0.0, syy = 0.0, sxy = 0.0;
        for(size_t k=0; k<px.size(); ++k)
        {
        const double dx = px[k]-mx, dy = py[k]-my;
        sxx += dx*dx;
        syy += dy*dy;
        sxy += dx*dy;
        }

        const double h = std::min(q->DXN[i+m],q->DYN[j+m]);
        if(px.size()<2 || !(sxx+syy>1.0e-6*h*h))
        return false;

        alpha = 0.5*std::atan2(2.0*sxy,sxx-syy);
        return true;
    };

    auto idx = [&](int i, int j) {return size_t(i-imin)*nj + (j-jmin);};
    auto coastal = [&](int i, int j)
    {
        return act(i,j) && ((i>imin && land(i-1,j)) || (i<imin+ni-1 && land(i+1,j)) || (j>jmin && land(i,j-1)) || (j<jmin+nj-1 && land(i,j+1)));
    };

    // first pass: the direction of every active cell next to land, as the unit vector of 2 alpha (a line has no sense);
    // second pass: their mean over the coastal cells within A 724 cells (the ghost cells from the neighbouring ranks)
    std::vector<double> c2(size_t(ni)*nj,0.0), s2(size_t(ni)*nj,0.0), ac(size_t(ni)*nj,0.0);
    std::vector<char> av(size_t(ni)*nj,0);

    if(radius>0)
    {
        for(int i=imin; i<imin+ni; ++i)
        for(int j=jmin; j<jmin+nj; ++j)
        {
        double a;
        if(coastal(i,j) && direction(i,j,a))
        {
        c2[idx(i,j)] = std::cos(2.0*a);
        s2[idx(i,j)] = std::sin(2.0*a);
        }
        }

        if(pgc!=nullptr)
        {
        slice4 C(q), S(q);
            for(int i=imin; i<imin+ni; ++i)
            for(int j=jmin; j<jmin+nj; ++j)
            {
            C(i,j) = c2[idx(i,j)];
            S(i,j) = s2[idx(i,j)];
            }
        pgc->gcsl_start4(q,C,50);
        pgc->gcsl_start4(q,S,50);
            for(int i=imin; i<imin+ni; ++i)
            for(int j=jmin; j<jmin+nj; ++j)
            {
            const bool in = inside(i,j);
            c2[idx(i,j)] = in ? C(i,j) : 0.0;
            s2[idx(i,j)] = in ? S(i,j) : 0.0;
            }
        }

        for(int i=imin; i<imin+ni; ++i)
        for(int j=jmin; j<jmin+nj; ++j)
        if(coastal(i,j))
        {
        double sc = 0.0, ss = 0.0;
            for(int a=std::max(i-radius,imin); a<=std::min(i+radius,imin+ni-1); ++a)
            for(int b=std::max(j-radius,jmin); b<=std::min(j+radius,jmin+nj-1); ++b)
            {
            sc += c2[idx(a,b)];
            ss += s2[idx(a,b)];
            }
            if(sc*sc + ss*ss>1.0e-12)
            {
            ac[idx(i,j)] = 0.5*std::atan2(ss,sc);
            av[idx(i,j)] = 1;
            }
        }
    }

    // the coastline of the active cell of a face, else (A 724 0, no direction) the face itself
    auto coastline = [&](int i, int j, int s)
    {
        const double pi = 3.14159265358979323846;
        if(av[idx(i,j)])
        return ac[idx(i,j)];
        return (s==0) ? 0.5*pi : 0.0;
    };

    for(int i=imin; i<imin+ni; ++i)
    for(int j=jmin; j<jmin+nj; ++j)
    for(int s=0; s<2; ++s)
    {
    const int i2 = (s==0) ? i+1 : i, j2 = (s==0) ? j : j+1;

        if(i2>=imin+ni || j2>=jmin+nj)
        continue;

    const size_t c = size_t(i-imin)*nj + (j-jmin);
        if((s==0 ? fe : fn)[c]>=0)
        continue;

    const bool a1 = act(i,j), a2 = act(i2,j2);
        if(a1==a2 || (!a1 && !inside(i,j)) || (!a2 && !inside(i2,j2)))
        continue;

    face f;
    f.obs = -1;
    f.kt2 = 0.0f;
    f.kr2 = float(kr2c);
    f.alpha = float(a1 ? coastline(i,j,s) : coastline(i2,j2,s));
    f.wd = wcoast;
    add(i,j,s,f);
    ++ncoast;
    }
}

void seastate_obstacle::update(lexer *q, fdm_seastate *e)
{
    const double pi = 3.14159265358979323846;
    const double ga = 2.6, gb = 0.15;
    seastate_param sp;

    // coasts: rebuild after a change of the active cells
    if(coast && !wraw.empty())
    {
    bool changed = false;
    for(int i=imin; i<imin+ni && !changed; ++i)
    for(int j=jmin; j<jmin+nj; ++j)
    if(wraw[size_t(i-imin)*nj + (j-jmin)]!=char(e->wet(i,j)==1))
    {
    changed = true;
    break;
    }

    // level 0: all ranks rebuild together (the exchange of the coastline directions)
    if(pgc!=nullptr)
    changed = pgc->globalimax(changed ? 1 : 0)>0;

    if(changed)
    build(q,e,oi,oj,gnx,gny,pgc);
    }

    auto hs = [&](int i, int j, double &tp)
    {
        tp = 0.0;
        if(e->wet(i,j)!=1 || e->N->spec(i,j)==nullptr)
        return 0.0;
        sp.compute(*e->grid,e->N->spec(i,j));
        tp = (sp.lpeak>=0) ? 2.0*pi/g->sig[sp.lpeak] : 0.0;
        return sp.Hs;
    };

    std::vector<double> E(nsig);

    for(size_t f=0; f<faces.size(); ++f)
    {
    const int n = faces[f].obs;

        if(n<0 || (type[n]==0 && kt[n]>=0.0) || type[n]==2)
        continue;

    const int i = fi[f], j = fj[f];
    const int i2 = (fdir[f]==0) ? i+1 : i, j2 = (fdir[f]==0) ? j : j+1;
    double tp1, tp2;
    const double h1 = hs(i,j,tp1), h2 = hs(i2,j2,tp2);
    const bool first = (h1>=h2);
    const double H = first ? h1 : h2, Tp = first ? tp1 : tp2;

        // porous structure: the spectrum of the incident cell
        if(type[n]==3)
        {
        const int ia = first ? i : i2, ja = first ? j : j2;
        float *kt2 = &fkt[size_t(faces[f].fq)*nsig], *kr2 = &fkr[size_t(faces[f].fq)*nsig];

            if(H>1.0e-6)
            {
            // the energy of the incident cell travelling towards the face (not the reflected waves)
            const float *N = e->N->spec(ia,ja);
            const seastate_grid &gr = *e->grid;
            const double sg = first ? 1.0 : -1.0;
                for(int l=0; l<nsig; ++l)
                {
                double s = 0.0;
                for(int m=0; m<gr.ndir; ++m)
                if(sg*((fdir[f]==0) ? gr.costh[m] : gr.sinth[m])>0.0)
                s += double(N[gr.bin(l,m)])*gr.wth[m];
                E[l] = s*gr.sig[l]*gr.dtheta;
                }

            double hd = 0.0, nd = 0.0;
            if(e->wet(i,j)==1)   {hd += e->depth(i,j);   nd += 1.0;}
            if(e->wet(i2,j2)==1) {hd += e->depth(i2,j2); nd += 1.0;}

            seastate_porous(gr,E.data(),hd/std::max(nd,1.0),pa[n],pb[n],pc[n],kt2,kr2);
            }
        continue;
        }

    double t = 1.0;

        if(H>1.0e-6)
        {
        const double wl = q->wd + 0.5*(e->eta(i,j) + e->eta(i2,j2));
        const double rc = zc[n]-wl;

            // d'Angremond et al. (1996)
            if(type[n]==1)
            t = seastate_dangremond(rc,H,Tp,pa[n],pb[n]);
            else
            {
            // Goda
            const double r = rc/H;
            if(r>=ga-gb)
            t = 0.0;
            else if(r>-gb-ga)
            t = 0.5*(1.0 - std::sin(pi/(2.0*ga)*(r+gb)));
            }
        }

    faces[f].kt2 = float(t*t);
    faces[f].kr2 = float(std::min(kr[n]*kr[n],1.0-t*t));
    }
}

const seastate_obstacle::face *seastate_obstacle::east(int i, int j) const
{
    if(faces.empty() || i<imin || i>=imin+ni || j<jmin || j>=jmin+nj)
    return nullptr;

    const int k = fe[size_t(i-imin)*nj + (j-jmin)];
    return k<0 ? nullptr : &faces[k];
}

const seastate_obstacle::face *seastate_obstacle::north(int i, int j) const
{
    if(faces.empty() || i<imin || i>=imin+ni || j<jmin || j>=jmin+nj)
    return nullptr;

    const int k = fn[size_t(i-imin)*nj + (j-jmin)];
    return k<0 ? nullptr : &faces[k];
}
