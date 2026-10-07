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


#include"seastate_forcing.h"
#include"seastate_swan_spc.h"
#include"seastate_grid.h"
#include<algorithm>
#include<cmath>
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

    // YYYYMMDD.HHMMSS
    bool is_datetime(const std::string &key)
    {
        const size_t d = key.find('.');
        if(d!=8 || key.size()<9)
        return false;
        for(size_t c=0; c<key.size(); ++c)
        if(c!=d && !std::isdigit(static_cast<unsigned char>(key[c])))
        return false;
        return true;
    }

    // days since 1970-01-01 of a proleptic Gregorian date (H. Hinnant)
    long days_from_civil(long y, long m, long d)
    {
        y -= m<=2;
        const long era = (y>=0 ? y : y-399)/400;
        const long yoe = y - era*400;
        const long doy = (153*(m + (m>2 ? -3 : 9)) + 2)/5 + d-1;
        const long doe = yoe*365 + yoe/4 - yoe/100 + doy;
        return era*146097 + doe - 719468;
    }
}

double seastate_datetime(double v)
{
    const long date = long(std::floor(v + 1.0e-9));
    long hms = std::lround((v - double(date))*1.0e6);
    long days = 0;

    if(hms>=240000)     // rounding at the end of a day
    {
    hms -= 240000;
    days = 1;
    }

    const long y = date/10000, m = (date/100)%100, d = date%100;
    const long hh = hms/10000, mm = (hms/100)%100, ss = hms%100;

    return double(days_from_civil(y,m,d) + days)*86400.0 + double(hh*3600 + mm*60 + ss);
}

// ------------------------------------------------------------------ SWAN spectra series

bool seastate_spc_series::open(const std::string &file_, const seastate_grid &g_, std::string &error)
{
    g = &g_;
    file = file_;
    in.open(file.c_str());

    if(!in.is_open())
    {
    error = "cannot open "+file;
    return false;
    }

    std::string line, key;

    if(!next_line(in,line,key) || key.rfind("SWAN",0)!=0)
    {
    error = file+" is not a SWAN spectral file (first line SWAN)";
    return false;
    }

    bool data = false;

    while(next_line(in,line,key))
    {
        if(key=="TIME")
        {
        timed = true;
        next_line(in,line,key);
        }
        else if(key=="LONLAT")
        {
        error = file+": LONLAT locations; convert the spectra to model coordinates with tools/seastate_forcing.py";
        return false;
        }
        else if(key=="LOCATIONS")
        {
        next_line(in,line,key);
        nloc = first_int(line);

            for(int n=0; n<nloc; ++n)
            {
            next_line(in,line,key);
            std::istringstream s(line);
            double x=0.0, y=0.0;
            s>>x>>y;
            xs.push_back(x);
            ys.push_back(y);
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
            next_line(in,line,key);
                if(n==0)
                energy = (key=="EnDens");
            next_line(in,line,key);
            next_line(in,line,key);
                if(n==0)
                exception = std::stod(key);
            }
        }
        else if(is_datetime(key) || key=="FACTOR" || key=="ZERO" || key=="NODATA" || key=="LOCATION")
        {
        data = true;
        break;
        }
    }

    if(nloc<1 || f.empty() || dir.empty() || !data)
    {
    error = file+": needs LOCATIONS, frequencies, directions (2D spectra) and data";
    return false;
    }

    if(key=="LOCATION")
    {
    error = file+": 1D spectra are not supported, write 2D spectra (SPEC2D)";
    return false;
    }

    if(timed && !is_datetime(key))
    {
    error = file+": TIME given, but no date-time line YYYYMMDD.HHMMSS before the data";
    return false;
    }

    pending = line;

    // Cartesian directions of propagation in [0,360), sorted
    for(double &d : dir)
    {
    if(nautical)
    d = 270.0 - d;
    d = std::fmod(d,360.0);
    if(d<0.0)
    d += 360.0;
    }

    order.resize(dir.size());
    std::iota(order.begin(),order.end(),size_t(0));
    std::sort(order.begin(),order.end(),[&](size_t a, size_t b){return dir[a]<dir[b];});

    std::vector<double> ds(dir.size());
    for(size_t m=0; m<dir.size(); ++m)
    ds[m] = dir[order[m]];
    dir.swap(ds);

    if(!read_record(rec[0],rtime[0],error))
    return false;

    if(eof)
    {
    rec[1] = rec[0];
    rtime[1] = rtime[0];
    }
    else if(!read_record(rec[1],rtime[1],error))
    return false;

    return true;
}

bool seastate_spc_series::read_record(std::vector<std::vector<float>> &out, double &t, std::string &error)
{
    std::string line = pending, key;
    {
    std::istringstream s(line.substr(0,line.find('$')));
    s>>key;
    }

    t = 0.0;

    if(timed)
    {
        if(!is_datetime(key))
        {
        error = file+": expected a date-time line, found '"+key+"'";
        return false;
        }

    t = seastate_datetime(std::stod(key));

        if(!next_line(in,line,key))
        {
        error = file+": record without data";
        return false;
        }
    }

    const size_t nf = f.size(), nd = dir.size();
    std::vector<double> raw(nf*nd), E(nf*nd);
    seastate_swan_spc conv;
    conv.f = f;
    conv.dir = dir;

    out.assign(nloc,std::vector<float>());

    for(int n=0; n<nloc; ++n)
    {
        if(n>0 && !next_line(in,line,key))
        {
        error = file+": record ends early";
        return false;
        }

    std::fill(raw.begin(),raw.end(),0.0);

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
                        if(!next_line(in,line,key))
                        {
                        error = file+": table ends early";
                        return false;
                        }
                    s.clear();
                    s.str(line);
                    s>>v;
                    }
                raw[l*nd+m] = (v==exception) ? 0.0 : std::max(v*factor,0.0);
                }
            }
        }
        else if(key!="ZERO" && key!="NODATA")
        {
        error = file+": expected FACTOR, ZERO or NODATA, found '"+key+"'";
        return false;
        }

        for(size_t l=0; l<nf; ++l)
        for(size_t m=0; m<nd; ++m)
        E[l*nd+m] = raw[l*nd+order[m]]/(energy ? 1025.0*9.81 : 1.0);

    conv.E = E;
    conv.to_grid(*g,out[n]);
    }

    ++records;

    // first line of the next record
    if(next_line(in,line,key))
    pending = line;
    else
    eof = true;

    if(!timed)
    eof = true;

    return true;
}

bool seastate_spc_series::advance(double t, std::string &error)
{
    while(t>rtime[1] && !eof)
    {
    std::swap(rec[0],rec[1]);
    rtime[0] = rtime[1];

    double tn;
        if(!read_record(rec[1],tn,error))
        return false;

        if(!(tn>rtime[0]))
        {
        error = file+": times must increase from record to record";
        return false;
        }

    rtime[1] = tn;
    }

    return true;
}

double seastate_spc_series::weight(double t) const
{
    if(!(rtime[1]>rtime[0]))
    return t>=rtime[1] ? 1.0 : 0.0;

    return std::min(std::max((t-rtime[0])/(rtime[1]-rtime[0]),0.0),1.0);
}

// ------------------------------------------------------------------ wind field series

bool seastate_wind_series::token(std::string &tok)
{
    while(pos>=buf.size())
    {
    std::string line;
        if(!std::getline(in,line))
        return false;

    const size_t c = line.find('$');
    std::istringstream s(line.substr(0,c));
    buf.clear();
    pos = 0;
    std::string w;
        while(s>>w)
        buf.push_back(w);
    }

    tok = buf[pos++];
    return true;
}

bool seastate_wind_series::open(const std::string &file_, std::string &error)
{
    file = file_;
    in.open(file.c_str());

    if(!in.is_open())
    {
    error = "cannot open "+file;
    return false;
    }

    std::string title;
    std::getline(in,title);

    std::string a[6];
    for(int k=0; k<6; ++k)
        if(!token(a[k]))
        {
        error = file+": header needs nx ny x0 y0 dx dy";
        return false;
        }

    nx = std::stoi(a[0]);
    ny = std::stoi(a[1]);
    x0 = std::stod(a[2]);
    y0 = std::stod(a[3]);
    dx = std::stod(a[4]);
    dy = std::stod(a[5]);

    if(nx<1 || ny<1 || !(dx>0.0) || !(dy>0.0))
    {
    error = file+": nx, ny must be positive and dx, dy > 0";
    return false;
    }

    if(!read_record(u[0],v[0],rtime[0],error))
    return false;

    std::string e2;
    if(!read_record(u[1],v[1],rtime[1],e2))
    {
        if(!eof)
        {
        error = e2;
        return false;
        }
    u[1] = u[0];
    v[1] = v[0];
    rtime[1] = rtime[0];
    }

    return true;
}

bool seastate_wind_series::read_record(std::vector<double> &uu, std::vector<double> &vv, double &t, std::string &error)
{
    std::string tok;

    if(!token(tok))
    {
    eof = true;
    error = file+": no more records";
    return false;
    }

    if(tok.find('.')!=8)
    {
    error = file+": expected a date-time YYYYMMDD.HHMMSS, found '"+tok+"'";
    return false;
    }

    t = seastate_datetime(std::stod(tok));

    const size_t n = size_t(nx)*ny;
    uu.resize(n);
    vv.resize(n);

    for(int c=0; c<2; ++c)
    for(size_t k=0; k<n; ++k)
    {
        if(!token(tok))
        {
        error = file+": record ends early";
        return false;
        }
    (c==0 ? uu : vv)[k] = std::stod(tok);
    }

    ++records;
    return true;
}

bool seastate_wind_series::advance(double t, std::string &error)
{
    while(t>rtime[1] && !eof)
    {
    std::vector<double> un, vn;
    double tn;

        if(!read_record(un,vn,tn,error))
        {
            if(eof)
            {
            error.clear();
            return true;
            }
        return false;
        }

        if(!(tn>rtime[1]))
        {
        error = file+": times must increase from record to record";
        return false;
        }

    u[0].swap(u[1]); v[0].swap(v[1]); rtime[0] = rtime[1];
    u[1].swap(un);   v[1].swap(vn);   rtime[1] = tn;
    }

    return true;
}

void seastate_wind_series::space(const std::vector<double> &a, double x, double y, double &val) const
{
    const double fx = std::min(std::max((x-x0)/dx,0.0),double(nx-1));
    const double fy = std::min(std::max((y-y0)/dy,0.0),double(ny-1));
    const int i0 = std::min(int(fx),nx-1), j0 = std::min(int(fy),ny-1);
    const int i1 = std::min(i0+1,nx-1), j1 = std::min(j0+1,ny-1);
    const double wx = fx - double(i0), wy = fy - double(j0);

    val = (1.0-wy)*((1.0-wx)*a[size_t(j0)*nx+i0] + wx*a[size_t(j0)*nx+i1])
        +      wy *((1.0-wx)*a[size_t(j1)*nx+i0] + wx*a[size_t(j1)*nx+i1]);
}

void seastate_wind_series::at(double t, double x, double y, double &uu, double &vv) const
{
    double w = 0.0;
    if(rtime[1]>rtime[0])
    w = std::min(std::max((t-rtime[0])/(rtime[1]-rtime[0]),0.0),1.0);
    else if(t>=rtime[1])
    w = 1.0;

    double ua, ub, va, vb;
    space(u[0],x,y,ua);
    space(u[1],x,y,ub);
    space(v[0],x,y,va);
    space(v[1],x,y,vb);

    uu = (1.0-w)*ua + w*ub;
    vv = (1.0-w)*va + w*vb;
}
