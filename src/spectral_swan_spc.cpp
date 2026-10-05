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

#include"spectral_swan_spc.h"
#include"spectral_grid.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<numeric>
#include<sstream>

namespace
{
    // next line that is not a comment ($) or empty; first word in key
    bool next_line(std::ifstream &in, std::string &line, std::string &key)
    {
        while(std::getline(in,line))
        {
            if(!line.empty() && line.back()=='\r')
            line.pop_back();

            const size_t c = line.find('$');
            const std::string data = line.substr(0,c);
            std::istringstream s(data);
            key.clear();
            s>>key;

            if(!key.empty())
            return true;
        }
        return false;
    }

    int first_int(const std::string &line)
    {
        std::istringstream s(line);
        int n=-1;
        s>>n;
        return n;
    }
}

bool spectral_swan_spc::read(const std::string &file, std::string &error)
{
    const double pi = 3.14159265358979323846;

    std::ifstream in(file.c_str());

    if(!in.is_open())
    {
    error = "cannot open "+file;
    return false;
    }

    f.clear(); dir.clear(); E.clear();

    std::string line, key;
    bool nautical=false, energy=false;
    int nloc=0;
    double exception=-99.0;

    if(!next_line(in,line,key) || key.rfind("SWAN",0)!=0)
    {
    error = file+" is not a SWAN spectral file (first line SWAN)";
    return false;
    }

    while(next_line(in,line,key))
    {
        if(key=="TIME")
        {
        next_line(in,line,key);     // time coding option
        }
        else if(key=="LOCATIONS" || key=="LONLAT")
        {
        next_line(in,line,key);
        nloc = first_int(line);

            for(int n=0; n<nloc; ++n)
            {
            next_line(in,line,key);
                if(n==0)
                {
                std::istringstream s(line);
                s>>x>>y;
                }
            }
        }
        else if(key=="AFREQ" || key=="RFREQ")
        {
        next_line(in,line,key);
        const int nf = first_int(line);

            for(int n=0; n<nf; ++n)
            {
            next_line(in,line,key);
            f.push_back(std::stod(key));
            }
        }
        else if(key=="NDIR" || key=="CDIR")
        {
        nautical = (key=="NDIR");
        next_line(in,line,key);
        const int nd = first_int(line);

            for(int n=0; n<nd; ++n)
            {
            next_line(in,line,key);
            dir.push_back(std::stod(key));
            }
        }
        else if(key=="QUANT")
        {
        next_line(in,line,key);
        const int nq = first_int(line);

            for(int n=0; n<nq; ++n)
            {
            next_line(in,line,key);             // name
                if(n==0)
                energy = (key=="EnDens");
            next_line(in,line,key);             // unit
            next_line(in,line,key);             // exception value
                if(n==0)
                exception = std::stod(key);
            }
        }
        else if(key=="FACTOR" || key=="ZERO" || key=="NODATA" || key=="LOCATION")
        break;
    }

    if(f.empty() || dir.empty())
    {
    error = file+": no frequencies or no directions (only 2D spectra are supported)";
    return false;
    }

    const size_t nf = f.size(), nd = dir.size();
    E.assign(nf*nd,0.0);

    if(key=="LOCATION")
    {
    error = file+": 1D spectra (QUANT 3) are not supported, write a 2D spectrum (SPECOUT ... SPEC2D)";
    return false;
    }

    if(key=="FACTOR")
    {
    next_line(in,line,key);
    const double factor = std::stod(key);

        for(size_t l=0; l<nf; ++l)
        {
        if(!next_line(in,line,key))
        {
        error = file+": table ends early";
        return false;
        }

        std::istringstream s(line);

            for(size_t m=0; m<nd; ++m)
            {
            double v=0.0;
                if(!(s>>v))
                {
                // long rows may continue on the next line
                if(!next_line(in,line,key))
                {
                error = file+": table ends early";
                return false;
                }
                s.clear();
                s.str(line);
                s>>v;
                }
            E[l*nd+m] = (v==exception) ? 0.0 : std::max(v*factor,0.0);
            }
        }
    }
    // ZERO or NODATA: zero spectrum

    if(energy)
    for(double &v : E)
    v /= 1025.0*9.81;

    // directions: Cartesian direction of propagation in [0,360), sorted ascending
    for(double &d : dir)
    {
    if(nautical)
    d = 270.0 - d;
    d = std::fmod(d,360.0);
    if(d<0.0)
    d += 360.0;
    }

    std::vector<size_t> o(nd);
    std::iota(o.begin(),o.end(),size_t(0));
    std::sort(o.begin(),o.end(),[&](size_t a, size_t b){return dir[a]<dir[b];});

    std::vector<double> d2(nd), E2(nf*nd);
    for(size_t m=0; m<nd; ++m)
    {
    d2[m] = dir[o[m]];
        for(size_t l=0; l<nf; ++l)
        E2[l*nd+m] = E[l*nd+o[m]];
    }
    dir.swap(d2);
    E.swap(E2);

    (void)pi;
    return true;
}

void spectral_swan_spc::to_grid(const spectral_grid &g, std::vector<float> &N) const
{
    const double pi = 3.14159265358979323846;
    const size_t nf = f.size(), nd = dir.size();

    N.assign(g.nbin,0.0f);

    if(nf==0 || nd==0)
    return;

    for(int l=0; l<g.nsig; ++l)
    {
    const double fl = g.f[l];

        if(fl<f.front() || fl>f.back())
        continue;

    size_t a = std::upper_bound(f.begin(),f.end(),fl) - f.begin();
    a = std::min(std::max(a,size_t(1)),nf-1);
    const double wf = (nf>1) ? (fl-f[a-1])/(f[a]-f[a-1]) : 0.0;

        for(int m=0; m<g.ndir; ++m)
        {
        const double th = g.theta[m]*180.0/pi;

        // periodic bracket in direction
        size_t b = std::upper_bound(dir.begin(),dir.end(),th) - dir.begin();
        const size_t m1 = (b==0) ? nd-1 : b-1;
        const size_t m2 = (b==nd) ? 0 : b;
        double d1 = dir[m1], d2 = dir[m2];
        if(d1>th) d1 -= 360.0;
        if(d2<th) d2 += 360.0;
        const double wd = (d2>d1) ? (th-d1)/(d2-d1) : 0.0;

        const double e1 = (1.0-wd)*E[(a-1)*nd+m1] + wd*E[(a-1)*nd+m2];
        const double e2 = (1.0-wd)*E[a*nd+m1]     + wd*E[a*nd+m2];
        const double Ef = (1.0-wf)*e1 + wf*e2;          // m2/Hz/degr

        N[g.bin(l,m)] = float(Ef*180.0/pi/(2.0*pi)/g.sig[l]);
        }
    }
}
