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

#ifndef SIXDOF_PTO_H_
#define SIXDOF_PTO_H_

#include<memory>
#include<string>
#include<vector>
#include<cmath>

// Power take-off (PTO) library of the 6DOF kernel.
//
//   joint     (6DOF_pto_joint.h)    : body state -> PTO coordinate q, qd and the generalised
//                                     direction J (6), qd = J^T [v; w], load on the body F J
//   element   (6DOF_pto_elements.h) : F(q,qd,t) along the joint and its Jacobians dF/dq, dF/dqd
//   composite (this file)           : sum of elements = one PTO
//
// The Jacobians of stiff elements (end-stops, latch) let the 6DOF step treat them linearly
// implicit together with the added mass (X 500 2). Controllers (6DOF_pto_controller.h) act on
// the elements via set_param() and on the composite's generator gain.

struct pto_state
{
    double q=0.0;       // joint coordinate [m] or [rad], 0 at the initial position
    double qd=0.0;      // joint velocity
    double t=0.0;
};

struct pto_output
{
    double F=0.0;       // force on the body along the joint (positive in +q)
    double dFdq=0.0;    // stiffness Jacobian
    double dFdqd=0.0;   // damping Jacobian
    double P=0.0;       // absorbed power -F qd of the useful elements
    double dFdq_s=0.0;  // part of the Jacobians from stiff elements (implicit with X 500 2)
    double dFdqd_s=0.0;
};

class pto_base
{
public:
    virtual ~pto_base() = default;

    // force at the state of the current RK stage; may be called several times per time step
    virtual pto_output force(const pto_state&) = 0;

    // once per time step with the accepted state at t^n: advance internal states (hydraulics, latching)
    virtual void accept(const pto_state&, double) {}

    // controller entry point
    virtual void set_param(const std::string&, double) {}

    virtual const char* name() const = 0;

    // contributes to the absorbed power (end-stops, friction losses do not)
    virtual bool useful() const {return true;}

    // stiff element: its Jacobians go into the linearly implicit step (X 500 2); the others
    // stay explicit, where the RK3/RK4 step is exact to its order
    virtual bool stiff() const {return false;}

    double F_last=0.0;
};

class pto_composite : public pto_base
{
public:
    // Sum of the elements. The generator part (useful elements: damper, spring-damper) is
    // scaled by gain (declutching) and limited by Fmax and Pmax (X 505); the other elements
    // (end-stops, latch) act unchanged. The element forces in F_last are as delivered.

    void add(std::unique_ptr<pto_base> e) {elem.push_back(std::move(e));}

    bool empty() const {return elem.empty();}

    pto_output force(const pto_state &s) override
    {
        pto_output out, gen;

        for(auto &e : elem)
        {
            const pto_output o = e->force(s);

            e->F_last = o.F;

            if(e->useful())
            {
            gen.F     += o.F;
            gen.dFdq  += o.dFdq;
            gen.dFdqd += o.dFdqd;
            }

            else
            {
            out.F     += o.F;
            out.dFdq  += o.dFdq;
            out.dFdqd += o.dFdqd;
            }

            if(e->stiff())
            {
            out.dFdq_s  += o.dFdq;
            out.dFdqd_s += o.dFdqd;
            }
        }

        // generator: gain (0 in the declutching window [t_off0, t_off1) at the stage time),
        // force limit, power limit
        const double g = (s.t>=t_off0 && s.t<t_off1) ? 0.0 : gain;
        double F = g*gen.F;
        double dFdq = g*gen.dFdq;
        double dFdqd = g*gen.dFdqd;

        sat_last = 0;

        if(Fmax>0.0 && fabs(F)>Fmax)
        {
            F = (F>0.0 ? Fmax : -Fmax);
            dFdq = dFdqd = 0.0;
            sat_last = 1;
        }

        if(Pmax>0.0 && fabs(F*s.qd)>Pmax)
        {
            F = (F>0.0 ? 1.0 : -1.0)*Pmax/fabs(s.qd);
            dFdq = 0.0;
            dFdqd = -F/s.qd;
            sat_last = 2;
        }

        const double r = (gen.F!=0.0) ? F/gen.F : 0.0;

        for(auto &e : elem)
        if(e->useful())
        e->F_last *= r;

        out.F     += F;
        out.dFdq  += dFdq;
        out.dFdqd += dFdqd;
        out.P      = -F*s.qd;

        F_last = out.F;

        return out;
    }

    void accept(const pto_state &s, double dt) override
    {
        for(auto &e : elem)
        e->accept(s,dt);
    }

    void set_param(const std::string &key, double val) override
    {
        for(auto &e : elem)
        e->set_param(key,val);
    }

    // set_param on the elements of one type only, e.g. ("damper","B",...)
    bool set_param(const std::string &type, const std::string &key, double val)
    {
        bool found=false;

        for(auto &e : elem)
        if(type==e->name())
        {
        e->set_param(key,val);
        found=true;
        }

        return found;
    }

    const char* name() const override {return "composite";}

    std::vector<std::unique_ptr<pto_base>> elem;

    double gain=1.0;    // generator scaling
    double t_off0=0.0, t_off1=0.0;  // generator off for t_off0 <= t < t_off1 (declutching)
    double Fmax=0.0;    // generator force limit, 0: none
    double Pmax=0.0;    // generator power limit, 0: none
    int sat_last=0;     // 0: free, 1: force limited, 2: power limited
};

#endif
