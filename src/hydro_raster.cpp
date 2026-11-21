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

#include"hydro_raster.h"
#include<fstream>
#include<sstream>
#include<cmath>
#include<algorithm>
#include<cctype>

bool hydro_raster::read(const std::string &file, std::string &error)
{
    std::ifstream in(file);

    if(!in)
    {
        error = "no raster file '" + file + "'";
        return false;
    }

    // header: keyword value pairs until the first numeric token
    bool xc = false, yc = false, hasdx = false, hasdy = false;
    double cs = -1.0;
    nc = nr = 0;
    has_nodata = false;
    std::string key;
    std::streampos pos = in.tellg();

    while(in>>key)
    {
        if(!key.empty() && (std::isdigit((unsigned char)key[0]) || key[0]=='-' || key[0]=='+' || key[0]=='.'))
        {
            in.clear();
            in.seekg(pos);
            break;
        }

        std::string k = key;
        std::transform(k.begin(),k.end(),k.begin(),[](unsigned char ch){return std::tolower(ch);});
        double val;

        if(!(in>>val))
        {
            error = "raster '" + file + "': no value for " + key;
            return false;
        }

        if(k=="ncols") nc = int(val);
        else if(k=="nrows") nr = int(val);
        else if(k=="xllcorner") {x0 = val; xc = false;}
        else if(k=="yllcorner") {y0 = val; yc = false;}
        else if(k=="xllcenter") {x0 = val; xc = true;}
        else if(k=="yllcenter") {y0 = val; yc = true;}
        else if(k=="cellsize") cs = val;
        else if(k=="dx") {dx = val; hasdx = true;}
        else if(k=="dy") {dy = val; hasdy = true;}
        else if(k=="nodata_value") {nodata = val; has_nodata = true;}
        else
        {
            error = "raster '" + file + "': unknown header keyword " + key;
            return false;
        }

        pos = in.tellg();
    }

    if(cs>0.0)
    {
        if(!hasdx) dx = cs;
        if(!hasdy) dy = cs;
    }

    if(nc<=0 || nr<=0 || dx<=0.0 || dy<=0.0)
    {
        error = "raster '" + file + "': ncols, nrows and cellsize (or dx, dy) must be positive";
        return false;
    }

    // centre-registered: the corner half a cell further out
    if(xc) x0 -= 0.5*dx;
    if(yc) y0 -= 0.5*dy;

    v.assign(size_t(nc)*size_t(nr),0.0);

    for(int r=0; r<nr; ++r)
    for(int c=0; c<nc; ++c)
    {
        double val;

        if(!(in>>val))
        {
            std::ostringstream os;
            os<<"raster '"<<file<<"': "<<size_t(nc)*size_t(nr)<<" values expected, "<<size_t(r)*size_t(nc)+size_t(c)<<" read";
            error = os.str();
            return false;
        }

        // the first row of the file is the northern one
        v[size_t(nr-1-r)*size_t(nc) + size_t(c)] = val;
    }

    return true;
}

bool hydro_raster::at(double x, double y, double &val) const
{
    const double fx = (x-x0)/dx, fy = (y-y0)/dy;

    if(!(fx>=0.0) || !(fy>=0.0) || fx>double(nc) || fy>double(nr))
    return false;

    const int ic = std::min(int(std::floor(fx)),nc-1), jr = std::min(int(std::floor(fy)),nr-1);
    val = v[size_t(jr)*size_t(nc) + size_t(ic)];

    if(has_nodata && std::fabs(val-nodata)<=1.0e-9*std::max(1.0,std::fabs(nodata)))
    return false;

    return true;
}
