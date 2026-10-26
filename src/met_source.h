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

#ifndef MET_SOURCE_H_
#define MET_SOURCE_H_

#include"met_file.h"
#include"met_frame.h"
#include"met_vortex.h"
#include<cmath>
#include<string>

/*--------------------------------------------------------------------
REEF3D atmospheric forcing (W 200): sources of the wind and the air
pressure. A source gives u10, v10 [m/s] along the model axes and the
air pressure pa [Pa] at any point (model coordinates [m]) and absolute
time t [s since 1970]. No lexer, unit test met_test.cpp.

  met_source_uniform  W 200 1: W 201 U10 dir, W 202 pa
  met_source_file     W 200 2: met-forcing.r3m (met_file)
  met_source_vortex   W 200 3: parametric vortex from met-vortex.trk (met_vortex),
                      evaluated at the longitude and latitude of the point (frame
                      W 205), the wind turned from (east, north) to the model axes

Later sources (parametric vortex from a track, blends of sources with a
radial weight, moving nests) implement the same interface.
--------------------------------------------------------------------*/

class met_source
{
public:
    virtual ~met_source() = default;

    // time t is coming: read what is needed
    virtual bool advance(double t, std::string &error) {return true;}

    // values at (x,y) and the time of the last advance
    virtual void at(double t, double x, double y, double &u, double &v, double &pa) const = 0;

    // box [m] the rank needs (file sources read only that part)
    virtual void window(double xs, double xe, double ys, double ye) {}

    // first time of the source [s since 1970], -1: none (uniform in time)
    virtual double first_time() const {return -1.0;}
};

class met_source_uniform final : public met_source
{
public:
    // U10 [m/s], dir [deg] counter-clockwise from the model x axis, the direction the wind blows to
    met_source_uniform(double U10, double dir, double pa_) : pa0(pa_)
    {
    const double a = dir*3.14159265358979323846/180.0;
    u0 = U10*std::cos(a);
    v0 = U10*std::sin(a);
    }

    void at(double t, double x, double y, double &u, double &v, double &pa) const override
    {
    u = u0;
    v = v0;
    pa = pa0;
    }

private:
    double u0, v0, pa0;
};

class met_source_file final : public met_source
{
public:
    bool open(const std::string &file, std::string &error) {return f.open(file,error);}

    bool advance(double t, std::string &error) override {return f.advance(t,error);}

    void at(double t, double x, double y, double &u, double &v, double &pa) const override
    {f.at(t,x,y,u,v,pa);}

    void window(double xs, double xe, double ys, double ye) override {f.window(xs,xe,ys,ye);}

    double first_time() const override {return f.time(0);}

    met_file f;
};

class met_source_vortex final : public met_source
{
public:
    met_source_vortex(const met_frame &fr) : frame(fr) {}

    bool open(const std::string &file, std::string &error) {return v.open(file,error);}

    void at(double t, double x, double y, double &u, double &v_, double &pa) const override
    {
    double lon, lat, ue, vn;
    frame.to_geo(x,y,lon,lat);
    v.at(t,lon,lat,ue,vn,pa);

    // true north at angle chi from the model x axis, east at chi - 90 deg
    const double chi = frame.north_geo(lon,lat);
    u  = ue*std::sin(chi) + vn*std::cos(chi);
    v_ = -ue*std::cos(chi) + vn*std::sin(chi);
    }

    double first_time() const override {return v.rec.front().t;}

    met_frame frame;
    met_vortex v;
};

#endif
