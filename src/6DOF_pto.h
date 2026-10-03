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

// Power take-off (PTO) library of the 6DOF kernel.
//
//   joint     (6DOF_pto_joint.h)    : body state -> PTO coordinate q, qd and the generalised
//                                     direction J (6), qd = J^T [v; w], load on the body F J
//   element   (6DOF_pto_elements.h) : F(q,qd,t) along the joint and its Jacobians dF/dq, dF/dqd
//   composite (this file)           : sum of elements = one PTO
//
// The Jacobians of stiff elements (end-stops) let the 6DOF step treat them linearly implicit
// together with the added mass (X 500 2). Controllers act on the elements via
// set_param() (not used yet).

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
    void add(std::unique_ptr<pto_base> e) {elem.push_back(std::move(e));}

    bool empty() const {return elem.empty();}

    pto_output force(const pto_state &s) override
    {
        pto_output out;

        for(auto &e : elem)
        {
            const pto_output o = e->force(s);

            e->F_last = o.F;

            out.F     += o.F;
            out.dFdq  += o.dFdq;
            out.dFdqd += o.dFdqd;

            if(e->useful())
            out.P += o.P;

            if(e->stiff())
            {
            out.dFdq_s  += o.dFdq;
            out.dFdqd_s += o.dFdqd;
            }
        }

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

    const char* name() const override {return "composite";}

    std::vector<std::unique_ptr<pto_base>> elem;
};

#endif
