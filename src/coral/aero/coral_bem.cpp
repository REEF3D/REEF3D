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

#include"coral/aero/coral_bem.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<limits>
#include<sstream>

namespace coral
{

namespace
{
    double prandtl(const bem_element &e, double phi, const bem_options &o)
    {
        const double s = std::fabs(std::sin(phi));
        if(s < 1.0e-12)
            return 1.0;
        double F = 1.0;
        if(o.tip_loss && e.r < e.R)
            F *= 2.0/M_PI*std::acos(std::clamp(std::exp(-0.5*e.B*(e.R - e.r)/(e.r*s)), 0.0, 1.0));
        if(o.tip_loss && e.r >= e.R)
            F = 0.0;
        if(o.hub_loss && e.R_hub > 0.0 && e.r > e.R_hub)
            F *= 2.0/M_PI*std::acos(std::clamp(std::exp(-0.5*e.B*(e.r - e.R_hub)/(e.R_hub*s)), 0.0, 1.0));
        return F;
    }

    struct coef_phi {double cl, cd, cm, cn, ct;};

    coef_phi coefs(const bem_element &e, double phi, const bem_options &o)
    {
        const airfoil_coef c = e.af ? e.af->coefficients(phi - e.theta, e.re) : airfoil_coef{};
        const double cd = o.drag_in_induction ? c.cd : 0.0;
        return coef_phi{c.cl, c.cd, c.cm, c.cl*std::cos(phi) + cd*std::sin(phi), c.cl*std::sin(phi) - cd*std::cos(phi)};
    }

    // Brent's method on [a, b] with f(a) f(b) < 0
    template<class Fn> double brent(Fn f, double a, double b, double fa, double fb, double tol, bool &ok)
    {
        double c = a, fc = fa, d = b - a, e = d;
        ok = false;
        for(int it=0; it<200; ++it)
        {
            if((fb > 0.0) == (fc > 0.0)) {c = a; fc = fa; d = e = b - a;}
            if(std::fabs(fc) < std::fabs(fb)) {a = b; b = c; c = a; fa = fb; fb = fc; fc = fa;}
            const double t = 2.0*std::numeric_limits<double>::epsilon()*std::fabs(b) + 0.5*tol;
            const double m = 0.5*(c - b);
            if(std::fabs(m) <= t || fb==0.0) {ok = true; return b;}
            if(std::fabs(e) >= t && std::fabs(fa) > std::fabs(fb))
            {
                double p, q, r;
                const double s = fb/fa;
                if(a==c) {p = 2.0*m*s; q = 1.0 - s;}
                else
                {
                    q = fa/fc; r = fb/fc;
                    p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
                    q = (q - 1.0)*(r - 1.0)*(s - 1.0);
                }
                if(p > 0.0) q = -q; else p = -p;
                if(2.0*p < std::min(3.0*m*q - std::fabs(t*q), std::fabs(e*q))) {e = d; d = p/q;}
                else {d = m; e = d;}
            }
            else {d = m; e = d;}
            a = b; fa = fb;
            b += std::fabs(d) > t ? d : (m > 0.0 ? t : -t);
            fb = f(b);
        }
        return b;
    }
}

double bem_residual(const bem_element &e, double Vx, double Vy, double phi, const bem_options &o, double *pa, double *pap)
{
    const double sp = std::sin(phi), cp = std::cos(phi);
    const double sig = e.B*e.chord/(2.0*M_PI*e.r);
    double F = prandtl(e, phi, o);
    const coef_phi c = coefs(e, phi, o);
    const double lam = Vy/Vx;
    if(F <= 0.0)
        F = 1.0e-12;
    const double k = sig*c.cn/(4.0*F*sp*sp);
    const double kp = o.tangential ? sig*c.ct/(4.0*F*sp*cp) : 0.0;
    double a, ap = kp/(1.0 - kp);
    if(!o.tangential) ap = 0.0;
    if(phi > 0.0)
    {
        if(k <= 2.0/3.0 || !o.buhl)
            a = k/(1.0 + k);
        else
        {
            const double g1 = 2.0*F*k - (10.0/9.0 - F), g2 = 2.0*F*k - F*(4.0/3.0 - F), g3 = 2.0*F*k - (25.0/9.0 - 2.0*F);
            a = std::fabs(g3) < 1.0e-6 ? 1.0 - 1.0/(2.0*std::sqrt(g2)) : (g1 - std::sqrt(std::max(g2, 0.0)))/g3;
        }
        if(pa) *pa = a;
        if(pap) *pap = ap;
        return sp/(1.0 - a) - cp/lam*(1.0/(1.0 + ap));
    }
    a = k > 1.0 ? k/(k - 1.0) : 0.0;
    if(pa) *pa = a;
    if(pap) *pap = ap;
    return sp*(1.0 - k) - cp/lam*(1.0 - kp);
}

bem_result bem_solve(const bem_element &e, double Vx, double Vy, double rho, const bem_options &o, double skew)
{
    bem_result res;
    if((o.tip_loss && e.r >= e.R) || (o.hub_loss && e.R_hub > 0.0 && e.r <= e.R_hub))
    {
        res.phi = std::atan2(Vx, Vy);
        res.alpha = res.phi - e.theta;
        res.F = 0.0;
        return res;
    }
    double phi;
    if(!o.induction || Vx <= 0.0 || Vy <= 0.0 || e.r <= 0.0)
    {
        phi = std::atan2(Vx, Vy);
        res.region = 0;
        res.a = res.ap = 0.0;
    }
    else
    {
        const double eps = 1.0e-6;
        auto f = [&](double p) {return bem_residual(e, Vx, Vy, p, o);};
        double lo, hi;
        double flo, fhi;
        bool found = false;
        const double br[3][2] = {{eps, 0.5*M_PI}, {-0.25*M_PI, -eps}, {0.5*M_PI, M_PI - eps}};
        for(int q=0; q<3 && !found; ++q)
        {
            lo = br[q][0]; hi = br[q][1];
            flo = f(lo); fhi = f(hi);
            if(flo*fhi < 0.0 || flo==0.0 || fhi==0.0)
            {
                found = true;
                res.region = q + 1;
                if(flo==0.0) phi = lo;
                else if(fhi==0.0) phi = hi;
                else phi = brent(f, lo, hi, flo, fhi, o.tol, res.converged);
            }
        }
        if(!found)
        {
            // no sign change at the ends of the three intervals (a pole of a' = k'/(1 - k') inside, or a polar that
            // is not smooth): scan them in 200 parts each for a sign change that is a root, not a pole
            for(int q=0; q<3 && !found; ++q)
            {
                const int nscan = 200;
                double p0 = br[q][0], f0 = f(p0);
                for(int j=1; j<=nscan && !found; ++j)
                {
                    const double p1 = br[q][0] + (br[q][1] - br[q][0])*j/nscan, f1 = f(p1);
                    if(f0*f1 <= 0.0)
                    {
                        bool ok;
                        const double r = f0==0.0 ? p0 : f1==0.0 ? p1 : brent(f, p0, p1, f0, f1, o.tol, ok);
                        if(std::fabs(f(r)) < 1.0e-6*(1.0 + std::fabs(f0) + std::fabs(f1)))
                        {
                            found = true;
                            phi = r;
                            res.region = q + 1;
                        }
                    }
                    p0 = p1; f0 = f1;
                }
            }
        }
        if(!found)
        {
            // no root: no induction
            phi = std::atan2(Vx, Vy);
            res.region = 0;
            res.converged = false;
        }
        else
            bem_residual(e, Vx, Vy, phi, o, &res.a, &res.ap);
        res.F = prandtl(e, phi, o);
    }
    if(skew != 0.0 && res.region != 0)
    {
        res.a *= 1.0 + skew;
        phi = std::atan2(Vx*(1.0 - res.a), Vy*(1.0 + res.ap));
    }
    const double ux = Vx*(1.0 - res.a), uy = Vy*(1.0 + res.ap);
    res.W = std::hypot(ux, uy);
    res.phi = phi;
    res.alpha = phi - e.theta;
    const airfoil_coef c = e.af ? e.af->coefficients(res.alpha, e.re) : airfoil_coef{};
    res.cl = c.cl; res.cd = c.cd; res.cm = c.cm;
    res.cn = c.cl*std::cos(phi) + c.cd*std::sin(phi);
    res.ct = c.cl*std::sin(phi) - c.cd*std::cos(phi);
    const double q = 0.5*rho*res.W*res.W;
    res.fn = q*e.chord*res.cn;
    res.ft = q*e.chord*res.ct;
    res.mz = q*e.chord*e.chord*res.cm;
    return res;
}

rotor_result rotor_steady(const rotor_def &rd, double U, double Omega, double pitch, double rho, const bem_options &o)
{
    rotor_result out;
    const int n = (int)rd.stations.size();
    out.el.resize(n);
    std::vector<double> fn(n), ftr(n);
    for(int i=0; i<n; ++i)
    {
        const blade_station &s = rd.stations[i];
        bem_element e;
        e.r = s.r; e.chord = s.chord; e.theta = s.twist + pitch;
        e.R = rd.R_tip; e.R_hub = rd.R_hub; e.B = rd.B;
        e.af = s.airfoil >= 0 && s.airfoil < (int)rd.airfoils.size() ? &rd.airfoils[s.airfoil] : nullptr;
        out.el[i] = bem_solve(e, U, Omega*s.r, rho, o);
        fn[i] = out.el[i].fn;
        ftr[i] = out.el[i].ft*s.r;
    }
    for(int i=0; i+1<n; ++i)
    {
        const double dr = rd.stations[i+1].r - rd.stations[i].r;
        out.T += 0.5*dr*(fn[i] + fn[i+1]);
        out.Q += 0.5*dr*(ftr[i] + ftr[i+1]);
    }
    out.T *= rd.B;
    out.Q *= rd.B;
    out.P = out.Q*Omega;
    const double A = M_PI*rd.R_tip*rd.R_tip;
    if(U > 0.0)
    {
        out.Cp = out.P/(0.5*rho*A*U*U*U);
        out.Ct = out.T/(0.5*rho*A*U*U);
        out.lambda = Omega*rd.R_tip/U;
    }
    return out;
}

bool read_aerodyn_blade(const std::string &file, double R_hub, std::vector<blade_station> &st, std::string &msg)
{
    std::ifstream in(file);
    if(!in)
    {
        msg = "cannot open '" + file + "'";
        return false;
    }
    std::stringstream ss;
    ss << in.rdbuf();
    if(!read_aerodyn_blade_string(ss.str(), R_hub, st, msg))
    {
        msg = file + ": " + msg;
        return false;
    }
    return true;
}

bool read_aerodyn_blade_string(const std::string &text, double R_hub, std::vector<blade_station> &st, std::string &msg)
{
    std::istringstream in(text);
    std::string line;
    int n = -1;
    st.clear();
    while(std::getline(in, line))
    {
        std::istringstream ls(line);
        std::vector<std::string> w;
        std::string x;
        while(ls >> x) w.push_back(x);
        if(n < 0)
        {
            if(w.size() >= 2 && w[1]=="NumBlNds")
            {
                try {n = std::stoi(w[0]);} catch(...) {msg = "NumBlNds is not a number"; return false;}
                if(n < 2) {msg = "NumBlNds must be at least 2"; return false;}
            }
            continue;
        }
        if((int)st.size()==n)
            break;
        if(w.size() < 7)
            continue;
        double v[7];
        bool num = true;
        for(int q=0; q<7 && num; ++q)
        {
            try {size_t k; v[q] = std::stod(w[q], &k); num = k==w[q].size();} catch(...) {num = false;}
        }
        if(!num)
            continue;      // the header lines of the table
        blade_station b;
        b.r = R_hub + v[0];
        b.twist = v[4]*M_PI/180.0;
        b.chord = v[5];
        b.airfoil = (int)v[6] - 1;
        if(!st.empty() && b.r <= st.back().r)
        {
            msg = "BlSpn must increase";
            return false;
        }
        st.push_back(b);
    }
    if(n < 0) {msg = "NumBlNds not found"; return false;}
    if((int)st.size() != n) {msg = std::to_string(st.size()) + " of " + std::to_string(n) + " blade nodes found"; return false;}
    return true;
}

}
