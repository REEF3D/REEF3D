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

#ifndef SIXDOF_PTO_CONTROLLER_H_
#define SIXDOF_PTO_CONTROLLER_H_

#include"6DOF_pto.h"
#include<string>
#include<vector>

// PTO controllers: called once per time step with the accepted state at t^n (before the
// PTO force of that step), they act on the composite through set_param() and its gain.
//
//   tuned passive (X 506) : B_opt = sqrt(B_rad^2 + (w (m + A) - C/w)^2) at the period T from a
//                           hydrodynamic table (pto_hydro.dat), set once at the start
//   latching      (X 507) : hold the body (latch element) for t_hold after each zero of qd
//   declutching   (X 508) : generator off for t_off after each zero of qd
// The windows start at the interpolated time of the zero and end at the RK stage time, so the
// timing error is about dt/2 instead of dt and does not accumulate.

struct pto_ctrl_in
{
    double t=0.0, q=0.0, qd=0.0;
};

class pto_controller
{
public:
    virtual ~pto_controller() = default;
    virtual void update(const pto_ctrl_in&, pto_composite&) = 0;
    virtual const char* name() const = 0;
    virtual double state() const {return 0.0;}
};

// sign change of qd between successive calls (zeros of qd are skipped, not counted); the time
// of the zero is interpolated linearly between the two calls, so that the switching times of
// the controllers do not depend on the time step
class pto_zero_crossing
{
public:
    bool detect(double t, double qd);
    void reset() {sref=0;}
    double t_cross=0.0;

private:
    int sref=0;
    double t_last=0.0, qd_last=0.0;
};

class pto_ctrl_latching : public pto_controller
{
public:
    pto_ctrl_latching(double);
    void update(const pto_ctrl_in&, pto_composite&) override;
    const char* name() const override {return "latching";}
    double state() const override {return latched ? 1.0 : 0.0;}

private:
    double t_hold, t_release=0.0;
    bool latched=false;
    pto_zero_crossing zc;
};

class pto_ctrl_declutching : public pto_controller
{
public:
    pto_ctrl_declutching(double);
    void update(const pto_ctrl_in&, pto_composite&) override;
    const char* name() const override {return "declutching";}
    double state() const override {return off ? 1.0 : 0.0;}

private:
    double t_off, t_on=0.0;
    bool off=false;
    pto_zero_crossing zc;
};

// hydrodynamic table along the joint: columns w [rad/s], A [kg], B_rad [kg/s];
// separators blank, tab or comma; lines not starting with a number are skipped
class pto_hydro_table
{
public:
    bool read(const std::string&);
    double A(double) const;
    double B(double) const;
    double wmin() const {return w_.empty() ? 0.0 : w_.front();}
    double wmax() const {return w_.empty() ? 0.0 : w_.back();}
    int size() const {return int(w_.size());}

private:
    double interp(const std::vector<double>&, double) const;
    std::vector<double> w_, A_, B_;
};

// optimal passive damping of a single-DOF absorber at the frequency w
double pto_bopt(double w, double m, double A, double Brad, double C);

#endif
