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

#include"met_file.h"
#include<algorithm>
#include<cmath>
#include<cstdint>
#include<cstdio>
#include<cstring>

/*--------------------------------------------------------------------
Binary met file reader, see met_file.h.
--------------------------------------------------------------------*/

namespace
{
    bool host_little()
    {
        const uint16_t one = 1;
        unsigned char b;
        std::memcpy(&b,&one,1);
        return b==1;
    }

    template<class T> void swap_bytes(T &v)
    {
        unsigned char b[sizeof(T)];
        std::memcpy(b,&v,sizeof(T));
        std::reverse(b,b+sizeof(T));
        std::memcpy(&v,b,sizeof(T));
    }

    // days since 1970-01-01 of a proleptic Gregorian date (H. Hinnant, days_from_civil)
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

double met_datetime(double v)
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

std::string met_date(double t)
{
    // seconds since 1970 -> YYYY-MM-DD HH:MM:SS (H. Hinnant, civil_from_days)
    long s = long(std::floor(t+0.5));
    long z = (s>=0 ? s : s-86399)/86400;
    const long sod = s - z*86400;
    z += 719468;
    const long era = (z>=0 ? z : z-146096)/146097;
    const long doe = z - era*146097;
    const long yoe = (doe - doe/1460 + doe/36524 - doe/146096)/365;
    const long doy = doe - (365*yoe + yoe/4 - yoe/100);
    const long mp = (5*doy + 2)/153;
    const long d = doy - (153*mp+2)/5 + 1;
    const long m = mp<10 ? mp+3 : mp-9;
    const long y = yoe + era*400 + (m<=2);
    char buf[128];
    snprintf(buf,sizeof(buf),"%04ld-%02ld-%02ld %02ld:%02ld:%02ld",y,m,d,sod/3600,(sod/60)%60,sod%60);
    return buf;
}

long long met_file::offset(int n) const
{
    return (long long)header_bytes + (long long)n*((long long)record_head_bytes + (long long)nvar*nx*ny*4);
}

bool met_file::open(const std::string &file_, std::string &error)
{
    file = file_;
    in.open(file.c_str(), std::ios::binary);

    if(!in.is_open())
    {
    error = "cannot open "+file;
    return false;
    }

    char magic[8];
    in.read(magic,8);

    if(!in || std::memcmp(magic,"REEF3DMT",8)!=0)
    {
    error = file+": not a REEF3D met file (magic REEF3DMT)";
    return false;
    }

    swap = !host_little();

    int32_t iv[8];
    in.read((char*)iv,sizeof(iv));
    double dv[13];
    in.read((char*)dv,sizeof(dv));

    if(!in)
    {
    error = file+": header too short";
    return false;
    }

    if(swap)
    {
    for(auto &a : iv) swap_bytes(a);
    for(auto &a : dv) swap_bytes(a);
    }

    if(iv[0]!=1)
    {
    error = file+": format version "+std::to_string(iv[0])+", this REEF3D reads version 1";
    return false;
    }

    nx = iv[1];  ny = iv[2];  nvar = iv[3];  moving = iv[4];  nrec = iv[5];  frame = iv[6];  epsg = iv[7];
    x0 = dv[0];  y0 = dv[1];  dx = dv[2];  dy = dv[3];
    lon0 = dv[4];  lat0 = dv[5];  rot = dv[6];  X0 = dv[7];  Y0 = dv[8];  avg = dv[9];

    if(nx<1 || ny<1 || nrec<1 || (nvar!=2 && nvar!=3) || (moving!=0 && moving!=1))
    {
    error = file+": bad header (nx "+std::to_string(nx)+", ny "+std::to_string(ny)+", nvar "+std::to_string(nvar)
            +", moving "+std::to_string(moving)+", nrec "+std::to_string(nrec)+")";
    return false;
    }

    if((nx>1 && !(dx>0.0)) || (ny>1 && !(dy>0.0)))
    {
    error = file+": the grid spacing must be positive";
    return false;
    }

    in.seekg(0,std::ios::end);
    const long long size = (long long)in.tellg();

    if(size!=offset(nrec))
    {
    error = file+": "+std::to_string(size)+" bytes, the header gives "+std::to_string(offset(nrec));
    return false;
    }

    // record times and origins
    times.resize(nrec);
    xr.resize(nrec);
    yr.resize(nrec);

    for(int n=0; n<nrec; ++n)
    {
    double h[3];
    in.seekg(offset(n));
    in.read((char*)h,sizeof(h));

        if(swap)
        for(auto &a : h) swap_bytes(a);

    times[n] = h[0];
    xr[n] = moving==1 ? h[1] : x0;
    yr[n] = moving==1 ? h[2] : y0;

        if(n>0 && !(times[n]>times[n-1]))
        {
        error = file+": the record times must increase (record "+std::to_string(n)+")";
        return false;
        }
    }

    if(!in)
    {
    error = file+": read error in the record headers";
    return false;
    }

    return true;
}

void met_file::window(double xs, double xe, double ys, double ye)
{
    if(!haswin)
    {
    wxs = xs;  wxe = xe;  wys = ys;  wye = ye;
    haswin = true;
    }
    else
    {
        if(xs>=wxs && xe<=wxe && ys>=wys && ye<=wye)
        return;

    wxs = std::min(wxs,xs);  wxe = std::max(wxe,xe);
    wys = std::min(wys,ys);  wye = std::max(wye,ye);
    }

    // the loaded records are read again with the new window
    rec[0].n = -1;
    rec[1].n = -1;
}

bool met_file::read_record(int n, record &r, std::string &error)
{
    r.n = n;
    r.t = times[n];
    r.xr0 = xr[n];
    r.yr0 = yr[n];

    // window of nodes: the box and two nodes around it, clamped to the grid
    if(haswin)
    {
    r.i0 = nx>1 ? int(std::floor((wxs - r.xr0)/dx)) - 2 : 0;
    r.i1 = nx>1 ? int(std::floor((wxe - r.xr0)/dx)) + 3 : 0;
    r.j0 = ny>1 ? int(std::floor((wys - r.yr0)/dy)) - 2 : 0;
    r.j1 = ny>1 ? int(std::floor((wye - r.yr0)/dy)) + 3 : 0;

    r.i0 = std::min(std::max(r.i0,0),nx-1);
    r.i1 = std::min(std::max(r.i1,0),nx-1);
    r.j0 = std::min(std::max(r.j0,0),ny-1);
    r.j1 = std::min(std::max(r.j1,0),ny-1);
    }
    else
    {
    r.i0 = 0;  r.i1 = nx-1;
    r.j0 = 0;  r.j1 = ny-1;
    }

    const int wi = r.i1 - r.i0 + 1;
    const int wj = r.j1 - r.j0 + 1;

    std::vector<float> *f[3] = {&r.u, &r.v, &r.pa};

    for(int q=0; q<3; ++q)
    f[q]->assign(size_t(wi)*size_t(wj), 101325.0f);

    for(int q=0; q<nvar; ++q)
    for(int j=r.j0; j<=r.j1; ++j)
    {
    const long long pos = offset(n) + record_head_bytes + 4LL*((long long)(q*ny + j)*nx + r.i0);
    in.seekg(pos);
    in.read((char*)&(*f[q])[size_t(j-r.j0)*wi], 4LL*wi);
    }

    if(!in)
    {
    error = file+": read error in record "+std::to_string(n);
    return false;
    }

    values_read += (long long)nvar*wi*wj;

    if(swap)
    for(int q=0; q<nvar; ++q)
    for(auto &a : *f[q]) swap_bytes(a);

    return true;
}

bool met_file::advance(double t, std::string &error)
{
    int k = 0;

    if(nrec>1)
    {
    // the last record at or before t, clamped so that k+1 exists
    k = int(std::upper_bound(times.begin(),times.end(),t) - times.begin()) - 1;
    k = std::min(std::max(k,0),nrec-2);
    }

    const int k1 = nrec>1 ? k+1 : 0;

    if(rec[0].n==k && rec[1].n==k1)
    return true;

    if(rec[1].n==k && nrec>1)
    {
    std::swap(rec[0],rec[1]);

        if(!read_record(k1,rec[1],error))
        return false;

    return true;
    }

    if(!read_record(k,rec[0],error))
    return false;

    if(nrec>1)
    {
        if(!read_record(k1,rec[1],error))
        return false;
    }
    else
    rec[1] = rec[0];

    return true;
}

double met_file::weight(double t) const
{
    if(nrec<2 || rec[1].t<=rec[0].t)
    return 0.0;

    const double w = (t - rec[0].t)/(rec[1].t - rec[0].t);

    return std::min(std::max(w,0.0),1.0);
}

void met_file::space(const record &r, double x, double y, double &u, double &v, double &pa) const
{
    // node i and weight wx along x, clamped to the grid and to the window of the record
    int i = 0, j = 0;
    double wx = 0.0, wy = 0.0;

    if(nx>1)
    {
    double fx = (x - r.xr0)/dx;
    fx = std::min(std::max(fx, double(r.i0)), double(r.i1));
    i = std::min(int(std::floor(fx)), std::max(r.i1-1, r.i0));
    wx = r.i1>r.i0 ? fx - double(i) : 0.0;
    }

    if(ny>1)
    {
    double fy = (y - r.yr0)/dy;
    fy = std::min(std::max(fy, double(r.j0)), double(r.j1));
    j = std::min(int(std::floor(fy)), std::max(r.j1-1, r.j0));
    wy = r.j1>r.j0 ? fy - double(j) : 0.0;
    }

    const int wi = r.i1 - r.i0 + 1;
    const int a  = (j - r.j0)*wi + (i - r.i0);
    const int ai = (wx>0.0) ? a+1 : a;
    const int aj = (wy>0.0) ? a+wi : a;
    const int aij = (wx>0.0) ? aj+1 : aj;

    // bilinear as nested linear interpolation f0 + w (f1 - f0): a uniform field is reproduced bitwise
    auto lerp = [](double f0, double f1, double w) {return f0 + w*(f1 - f0);};
    auto bil = [&](const std::vector<float> &f)
    {return lerp(lerp(f[a],f[ai],wx), lerp(f[aj],f[aij],wx), wy);};

    u  = bil(r.u);
    v  = bil(r.v);
    pa = bil(r.pa);
}

void met_file::at(double t, double x, double y, double &u, double &v, double &pa) const
{
    space(rec[0],x,y,u,v,pa);

    const double w = weight(t);

    if(w>0.0)
    {
    double u1, v1, pa1;
    space(rec[1],x,y,u1,v1,pa1);

    // f0 + w (f1 - f0): a field constant in time is reproduced bitwise
    u  = u  + w*(u1 - u);
    v  = v  + w*(v1 - v);
    pa = pa + w*(pa1 - pa);
    }
}

void met_file::origin(double t, double &ox, double &oy) const
{
    ox = rec[0].xr0;
    oy = rec[0].yr0;

    const double w = weight(t);

    if(w>0.0)
    {
    ox = ox + w*(rec[1].xr0 - ox);
    oy = oy + w*(rec[1].yr0 - oy);
    }
}

bool met_file::sample(int n, int var, double x, double y, double &val, std::string &error)
{
    if(var>=nvar)
    {
    val = 101325.0;
    return true;
    }

    // nodes and weights as space(), clamped to the whole grid
    int i = 0, j = 0;
    double wx = 0.0, wy = 0.0;

    if(nx>1)
    {
    double fx = (x - xr[n])/dx;
    fx = std::min(std::max(fx,0.0),double(nx-1));
    i = std::min(int(std::floor(fx)),nx-2);
    wx = fx - double(i);
    }

    if(ny>1)
    {
    double fy = (y - yr[n])/dy;
    fy = std::min(std::max(fy,0.0),double(ny-1));
    j = std::min(int(std::floor(fy)),ny-2);
    wy = fy - double(j);
    }

    const int ni = nx>1 ? 2 : 1;
    float f[2][2] = {{0.0f,0.0f},{0.0f,0.0f}};

    for(int b=0; b<(ny>1 ? 2 : 1); ++b)
    {
    in.seekg(offset(n) + record_head_bytes + 4LL*((long long)(var*ny + j + b)*nx + i));
    in.read((char*)f[b],4LL*ni);
    }

    if(!in)
    {
    error = file+": read error in record "+std::to_string(n);
    return false;
    }

    if(swap)
    for(auto &row : f)
    for(auto &a : row)
    swap_bytes(a);

    auto lerp = [](double f0, double f1, double w) {return f0 + w*(f1 - f0);};
    const double r0 = nx>1 ? lerp(f[0][0],f[0][1],wx) : f[0][0];
    const double r1 = nx>1 ? lerp(f[1][0],f[1][1],wx) : f[1][0];

    val = ny>1 ? lerp(r0,r1,wy) : r0;

    return true;
}
