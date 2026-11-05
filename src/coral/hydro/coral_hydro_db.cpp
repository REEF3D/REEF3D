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

#include"coral/hydro/coral_hydro_db.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<map>
#include<sstream>
#include<stdexcept>

namespace coral
{

void hydro_db::resize(int N, int n, int m)
{
    const int nd = 6*N;
    body.assign(N, body_info());
    C = MatX::Zero(nd, nd);
    w.assign(n, 0.0);
    beta.assign(m, 0.0);
    A.assign(n, MatX::Zero(nd, nd));
    B.assign(n, MatX::Zero(nd, nd));
    Ainf = MatX::Zero(nd, nd);
    has_ainf = false;
    X.assign(n, std::vector<VecXc>(m, VecXc::Zero(nd)));
}

void hydro_db::check() const
{
    const int nd = ndof(), n = nfreq(), m = nheading();
    if(nbody() < 1)
        throw std::runtime_error("coral::hydro_db: no bodies");
    if(C.rows()!=nd || C.cols()!=nd)
        throw std::runtime_error("coral::hydro_db: hydrostatics has the wrong size");
    if((int)A.size()!=n || (int)B.size()!=n || (int)X.size()!=n)
        throw std::runtime_error("coral::hydro_db: arrays and frequencies do not match");
    for(int k=0; k<n; ++k)
    {
        if(k>0 && !(w[k] > w[k-1]))
            throw std::runtime_error("coral::hydro_db: frequencies must be ascending");
        if(A[k].rows()!=nd || B[k].rows()!=nd || (int)X[k].size()!=m)
            throw std::runtime_error("coral::hydro_db: wrong size at a frequency");
    }
}

// linear interpolation in w: index k and weight f (value = (1-f) v[k] + f v[k+1]); -1 below, n-1 above
static void bracket(const std::vector<double> &w, double om, int &k, double &f)
{
    const int n = (int)w.size();
    if(n==0 || om < w.front()) {k = -1; f = 0.0; return;}
    if(om >= w.back()) {k = n-1; f = 0.0; return;}
    k = (int)(std::upper_bound(w.begin(), w.end(), om) - w.begin()) - 1;
    f = (om - w[k])/(w[k+1] - w[k]);
}

MatX hydro_db::added_mass(double om) const
{
    int k; double f;
    bracket(w, om, k, f);
    if(k < 0) return A.front();
    if(k >= nfreq()-1) return A.back();
    return (1.0 - f)*A[k] + f*A[k+1];
}

MatX hydro_db::damping(double om) const
{
    int k; double f;
    bracket(w, om, k, f);
    if(k < 0 || (k >= nfreq()-1 && om > w.back())) return MatX::Zero(ndof(), ndof());
    if(k >= nfreq()-1) return B.back();
    return (1.0 - f)*B[k] + f*B[k+1];
}

VecXc hydro_db::excitation(double om, double hd) const
{
    const int nd = ndof();
    if(nfreq()==0 || nheading()==0)
        return VecXc::Zero(nd);
    // heading: linear between the given headings, nearest outside
    int j = 0; double fb = 0.0;
    if(nheading() > 1)
    {
        if(hd <= beta.front()) {j = 0; fb = 0.0;}
        else if(hd >= beta.back()) {j = nheading()-1; fb = 0.0;}
        else
        {
            j = (int)(std::upper_bound(beta.begin(), beta.end(), hd) - beta.begin()) - 1;
            fb = (hd - beta[j])/(beta[j+1] - beta[j]);
        }
    }
    auto at = [&](int kk)
    {
        if(fb==0.0 || j+1 >= nheading()) return X[kk][j];
        return VecXc((1.0 - fb)*X[kk][j] + fb*X[kk][j+1]);
    };
    int k; double f;
    bracket(w, om, k, f);
    if(k < 0 || (k >= nfreq()-1 && om > w.back())) return VecXc::Zero(nd);
    if(k >= nfreq()-1) return at(k);
    return (1.0 - f)*at(k) + f*at(k+1);
}

// ---------------------------------------------------------------- CORAL format

static std::vector<std::string> tokens(const std::string &file)
{
    std::ifstream in(file);
    if(!in)
        throw std::runtime_error("coral::hydro_db: cannot open " + file);
    std::vector<std::string> t;
    std::string line;
    while(std::getline(in, line))
    {
        const size_t c = line.find('#');
        if(c != std::string::npos)
            line = line.substr(0, c);
        std::istringstream s(line);
        std::string x;
        while(s >> x)
            t.push_back(x);
    }
    return t;
}

void hydro_db::read(const std::string &file)
{
    const std::vector<std::string> t = tokens(file);
    size_t p = 0;
    auto need = [&](size_t n) {if(p + n > t.size()) throw std::runtime_error("coral::hydro_db: " + file + " ends early");};
    auto str = [&]() {need(1); return t[p++];};
    auto num = [&]() {need(1); return std::stod(t[p++]);};
    auto mat = [&](MatX &M) {const int nd = ndof(); M.resize(nd, nd); for(int i=0; i<nd; ++i) for(int j=0; j<nd; ++j) M(i,j) = num();};

    if(str() != "coral_hydro_db")
        throw std::runtime_error("coral::hydro_db: " + file + " is not a CORAL hydrodynamic database");
    num();    // version
    int N = -1, nf = -1, nh = -1;
    std::vector<double> wl, bl;
    std::vector<body_info> bi;
    int k = -1;
    while(p < t.size())
    {
        const std::string key = str();
        if(key=="rho") rho = num();
        else if(key=="g") g = num();
        else if(key=="depth") depth = num();
        else if(key=="bodies")
        {
            N = (int)num();
            resize(N, 0, 0);
            bi.assign(N, body_info());
        }
        else if(key=="body")
        {
            const int i = (int)num() - 1;
            if(i < 0 || i >= N) throw std::runtime_error("coral::hydro_db: body index out of range");
            while(p < t.size())
            {
                const std::string k2 = t[p];
                if(k2=="name") {++p; bi[i].name = str();}
                else if(k2=="ref") {++p; bi[i].ref(0) = num(); bi[i].ref(1) = num(); bi[i].ref(2) = num();}
                else if(k2=="volume") {++p; bi[i].volume = num();}
                else if(k2=="cob") {++p; bi[i].cob(0) = num(); bi[i].cob(1) = num(); bi[i].cob(2) = num();}
                else break;
            }
        }
        else if(key=="hydrostatics") {if(N<1) throw std::runtime_error("coral::hydro_db: bodies first"); mat(C);}
        else if(key=="frequencies") {nf = (int)num(); wl.resize(nf); for(auto &x : wl) x = num();}
        else if(key=="headings") {nh = (int)num(); bl.resize(nh); for(auto &x : bl) x = num()*M_PI/180.0;}
        else if(key=="added_mass_inf") {mat(Ainf); has_ainf = true;}
        else if(key=="frequency")
        {
            if(nf < 0 || nh < 0 || N < 1)
                throw std::runtime_error("coral::hydro_db: bodies, frequencies and headings before the first frequency");
            if(A.size() != (size_t)nf)
            {
                const MatX Cs = C, Ai = Ainf;
                const bool ha = has_ainf;
                resize(N, nf, nh);
                C = Cs; Ainf = Ai; has_ainf = ha;
                w = wl; beta = bl;
            }
            k = (int)num() - 1;
            if(k < 0 || k >= nf) throw std::runtime_error("coral::hydro_db: frequency index out of range");
        }
        else if(key=="added_mass") {if(k<0) throw std::runtime_error("coral::hydro_db: added_mass before frequency"); mat(A[k]);}
        else if(key=="damping") {if(k<0) throw std::runtime_error("coral::hydro_db: damping before frequency"); mat(B[k]);}
        else if(key=="excitation")
        {
            if(k<0) throw std::runtime_error("coral::hydro_db: excitation before frequency");
            const int j = (int)num() - 1;
            if(j < 0 || j >= nh) throw std::runtime_error("coral::hydro_db: heading index out of range");
            for(int i=0; i<ndof(); ++i)
            {
                const double re = num(), im = num();
                X[k][j](i) = std::complex<double>(re, im);
            }
        }
        else
            throw std::runtime_error("coral::hydro_db: unknown keyword " + key + " in " + file);
    }
    body = bi;
    if(A.size() != (size_t)std::max(nf, 0))
    {
        const MatX Cs = C, Ai = Ainf;
        const bool ha = has_ainf;
        resize(N, std::max(nf, 0), std::max(nh, 0));
        C = Cs; Ainf = Ai; has_ainf = ha;
        w = wl; beta = bl;
        body = bi;
    }
    check();
}

void hydro_db::write(const std::string &file) const
{
    check();
    std::ofstream out(file);
    if(!out)
        throw std::runtime_error("coral::hydro_db: cannot write " + file);
    out<<std::setprecision(17);
    const int nd = ndof();
    auto mat = [&](const MatX &M) {for(int i=0; i<nd; ++i) {for(int j=0; j<nd; ++j) out<<(j?" ":"")<<M(i,j); out<<"\n";}};
    out<<"# CORAL hydrodynamic database (REEF3D::CORAL)\ncoral_hydro_db 1\n";
    out<<"rho "<<rho<<"\ng "<<g<<"\ndepth "<<depth<<"\nbodies "<<nbody()<<"\n";
    for(int i=0; i<nbody(); ++i)
    {
        const body_info &b = body[i];
        out<<"body "<<i+1<<" name "<<(b.name.empty() ? "body"+std::to_string(i+1) : b.name)
           <<" ref "<<b.ref(0)<<" "<<b.ref(1)<<" "<<b.ref(2);
        if(b.volume >= 0.0)
            out<<" volume "<<b.volume<<" cob "<<b.cob(0)<<" "<<b.cob(1)<<" "<<b.cob(2);
        out<<"\n";
    }
    out<<"hydrostatics\n"; mat(C);
    out<<"frequencies "<<nfreq()<<"\n";
    for(int k=0; k<nfreq(); ++k) out<<(k?" ":"")<<w[k];
    out<<"\nheadings "<<nheading()<<"\n";
    for(int j=0; j<nheading(); ++j) out<<(j?" ":"")<<beta[j]*180.0/M_PI;
    out<<"\n";
    if(has_ainf) {out<<"added_mass_inf\n"; mat(Ainf);}
    for(int k=0; k<nfreq(); ++k)
    {
        out<<"frequency "<<k+1<<"   # w = "<<w[k]<<" rad/s\nadded_mass\n"; mat(A[k]);
        out<<"damping\n"; mat(B[k]);
        for(int j=0; j<nheading(); ++j)
        {
            out<<"excitation "<<j+1<<"\n";
            for(int i=0; i<nd; ++i) out<<X[k][j](i).real()<<" "<<X[k][j](i).imag()<<"\n";
        }
    }
}

// ---------------------------------------------------------------- WAMIT

static int rot(int i) {return ((i-1)%6) >= 3 ? 1 : 0;}

void hydro_db::read_wamit(const std::string &root, int N, double ulen, double rho_, double g_, double depth_)
{
    rho = rho_;
    g = g_;
    depth = depth_;
    const int nd = 6*N;
    const double L = ulen;

    // .1: added mass and damping
    std::ifstream f1(root + ".1");
    if(!f1)
        throw std::runtime_error("coral::hydro_db: cannot open " + root + ".1");
    std::map<double, std::pair<MatX,MatX>> per;    // keyed by omega
    MatX Ai = MatX::Zero(nd, nd);
    bool hai = false;
    std::string line;
    while(std::getline(f1, line))
    {
        std::istringstream s(line);
        double P; int i, j; double a, b = 0.0;
        if(!(s >> P >> i >> j >> a))
            continue;
        const bool hb = (bool)(s >> b);
        if(i < 1 || j < 1 || i > nd || j > nd)
            continue;
        const int kk = 3 + rot(i) + rot(j);
        const double Ad = a*rho*std::pow(L, kk);
        if(P < 0.0) {Ai(i-1,j-1) = Ad; hai = true; continue;}    // infinite frequency
        if(P == 0.0) continue;                                     // zero frequency
        const double om = 2.0*M_PI/P;
        auto &e = per[om];
        if(e.first.size()==0) {e.first = MatX::Zero(nd, nd); e.second = MatX::Zero(nd, nd);}
        e.first(i-1,j-1) = Ad;
        if(hb) e.second(i-1,j-1) = b*rho*om*std::pow(L, kk);
    }
    // .3: excitation (diffraction + Froude-Krylov)
    std::map<double, std::map<double, VecXc>> exc;
    std::ifstream f3(root + ".3");
    if(f3)
    {
        while(std::getline(f3, line))
        {
            std::istringstream s(line);
            double P, bdeg, mod, pha, re, im; int i;
            if(!(s >> P >> bdeg >> i >> mod >> pha >> re >> im))
                continue;
            if(P <= 0.0 || i < 1 || i > nd)
                continue;
            const double om = 2.0*M_PI/P;
            auto &v = exc[om][bdeg];
            if(v.size()==0) v = VecXc::Zero(nd);
            v(i-1) = std::complex<double>(re, im)*rho*g*std::pow(L, 2 + rot(i));
        }
    }
    std::vector<double> bl;
    if(!exc.empty())
        for(const auto &e : exc.begin()->second)
            bl.push_back(e.first);
    if(bl.empty())
        bl.push_back(0.0);

    resize(N, (int)per.size(), (int)bl.size());
    int k = 0;
    for(const auto &e : per)
    {
        w[k] = e.first;
        A[k] = e.second.first;
        B[k] = e.second.second;
        auto it = exc.find(e.first);
        for(int j=0; j<(int)bl.size(); ++j)
        {
            X[k][j] = VecXc::Zero(nd);
            if(it != exc.end())
            {
                auto jt = it->second.find(bl[j]);
                if(jt != it->second.end())
                    X[k][j] = jt->second;
            }
        }
        ++k;
    }
    for(int j=0; j<(int)bl.size(); ++j)
        beta[j] = bl[j]*M_PI/180.0;
    Ainf = Ai;
    has_ainf = hai;

    // .hst: hydrostatics
    std::ifstream fh(root + ".hst");
    if(fh)
    {
        while(std::getline(fh, line))
        {
            std::istringstream s(line);
            int i, j; double c;
            if(!(s >> i >> j >> c))
                continue;
            if(i < 1 || j < 1 || i > nd || j > nd)
                continue;
            C(i-1,j-1) = c*rho*g*std::pow(L, 2 + rot(i) + rot(j));
        }
    }
    for(int i=0; i<N; ++i)
        body[i].name = "body" + std::to_string(i+1);
    check();
}

}
