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

#ifndef SIXDOF_PTO_ELEMENTS_H_
#define SIXDOF_PTO_ELEMENTS_H_

#include"6DOF_pto.h"

// linear damper (X 502): F = -B qd
class pto_damper : public pto_base
{
public:
    pto_damper(double);

    pto_output force(const pto_state&) override;
    void set_param(const std::string&, double) override;
    const char* name() const override {return "damper";}

private:
    double B;
};

// spring-damper (X 503): F = -K q - B qd, K < 0 allowed (reactive control)
class pto_springdamper : public pto_base
{
public:
    pto_springdamper(double, double);

    pto_output force(const pto_state&) override;
    void set_param(const std::string&, double) override;
    const char* name() const override {return "spring-damper";}

private:
    double K, B;
};

// end-stops (X 504): penalty spring-damper outside [qmin,qmax], pushes only, no absorbed power
class pto_endstop : public pto_base
{
public:
    pto_endstop(double, double, double, double);

    pto_output force(const pto_state&) override;
    const char* name() const override {return "end-stop";}
    bool useful() const override {return false;}
    bool stiff() const override {return true;}

private:
    double qmin, qmax, K, C;
};

// latch (X 507): stiff spring-damper holding the position q0 while on and t < t_release, set by
// the latching controller via set_param("q0"), ("t_release") and ("latch",1/0); no absorbed power
class pto_latch : public pto_base
{
public:
    pto_latch(double, double);

    pto_output force(const pto_state&) override;
    void set_param(const std::string&, double) override;
    const char* name() const override {return "latch";}
    bool useful() const override {return false;}
    bool stiff() const override {return true;}

    bool on=false;
    double q0=0.0;
    double t_release=1.0e30;    // released at this (stage) time

private:
    double K, C;
};

#endif
