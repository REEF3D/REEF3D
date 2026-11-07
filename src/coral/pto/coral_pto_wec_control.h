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
// REEF3D::CORAL: wave energy controllers of the PTO beyond the passive, latching and declutching
// controllers of coral_pto_controller.h. They act once per converged step on the generator of the
// PTO (its spring-damper element: stiffness K, damping B and force F0, F = F0 - K q - B qd), so the
// force limits, the end stops, the losses and the generator model of the PTO stay in place.
//
//   schedule  K, B, F0 from a table over time (t k c [F0]), linear in time: gains that follow the sea
//             state of a long run (tuned per sea state, e.g. from the power matrix)
//   user      a controller of the user's own in a shared library (C interface below), as WEC-Sim's
//             Simulink controllers or the DISCON routines of wind turbines; sampled at its period
//   mpc       model predictive control of a single-DOF absorber: a linear model of the body along the
//             PTO (mass, A_inf, hydrostatic C, linear damping, radiation memory as states), the
//             excitation force predicted over the horizon (from CORAL's linear wave excitation or from
//             a table); every period T_c it maximises the absorbed energy over the horizon minus a
//             weight on the force, with |F| <= F_max and a soft stroke limit (quadratic penalty on the
//             predicted excursions beyond q_max), and holds the first force (a horizon of about two
//             wave periods; shorter ones drain the body towards the end of the horizon)
//
// User library (C linkage, compiled e.g. with  cc -shared -fPIC -O2 my_ctrl.c -o my_ctrl.so):
//   int coral_pto_control(int n_in, const double *in, int n_out, double *out, void **user);
//     in:  t, q, qd, F, P (accepted state of the PTO, F its force, P the absorbed power), then the
//          parameters of the input line ('params p1 p2 ...')
//     out: K, B, F0 of the generator (given in with their current values; leave them to keep them)
//     user: a pointer the library may set for its own data (null at the first call)
//     returns 0, otherwise the run stops with the code
//   optional: int coral_pto_control_save(void *user, double *buf, int n_max)     (returns the count)
//             int coral_pto_control_load(void **user, const double *buf, int n)  (restart, 0 if fine; *user
//                 as for coral_pto_control, null if the library was not called yet)
//             void coral_pto_control_end(void *user)
// The first call is at the start (t = t0); later calls when t reaches the next sample time.

#ifndef CORAL_PTO_WEC_CONTROL_H_
#define CORAL_PTO_WEC_CONTROL_H_

#include"coral/pto/coral_pto_controller.h"
#include"coral/core/coral_types.h"
#include<functional>
#include<string>
#include<vector>

namespace coral
{

// --- gain schedule ---
class pto_ctrl_schedule : public pto_controller
{
public:
    // t increasing; k, c, F0 one value per time
    pto_ctrl_schedule(const std::vector<double> &t, const std::vector<double> &k, const std::vector<double> &c,
                      const std::vector<double> &F0);
    void update(const pto_ctrl_in&, pto_composite&) override;
    const char* name() const override {return "schedule";}
    // reads "t k c [F0]" lines ('#' comments); false with msg on an error
    static bool read(const std::string &file, std::vector<double> &t, std::vector<double> &k,
                     std::vector<double> &c, std::vector<double> &F0, std::string &msg);
    double value(const std::vector<double>&, double t) const;

private:
    std::vector<double> t_, k_, c_, F0_;
};

// --- user library ---
class pto_ctrl_user : public pto_controller
{
public:
    // throws std::runtime_error if the library or its function cannot be loaded
    // K, B, F0: the gains of the generator before the first call
    pto_ctrl_user(const std::string &lib, const std::vector<double> &params, double period,
                  double K=0.0, double B=0.0, double F0=0.0);
    ~pto_ctrl_user() override;
    void update(const pto_ctrl_in&, pto_composite&) override;
    const char* name() const override {return "user";}
    void save_state(std::vector<double>&) const override;
    void load_state(const std::vector<double>&, size_t&) override;
    void channel_names(const std::string&, std::vector<std::string>&) const override;
    void channel_values(std::vector<double>&) const override;

private:
    using fn_t = int (*)(int, const double*, int, double*, void**);
    using save_t = int (*)(void*, double*, int);
    using load_t = int (*)(void**, const double*, int);
    using end_t = void (*)(void*);
    void *handle_ = nullptr;
    fn_t fn_ = nullptr;
    save_t save_ = nullptr;
    load_t load_ = nullptr;
    end_t end_ = nullptr;
    void *user_ = nullptr;
    std::string lib_;
    std::vector<double> params_;
    double period_ = 0.0, t_next_ = -1.0e300;
    double out_[3] = {0.0, 0.0, 0.0};
};

// --- model predictive control ---
struct mpc_model
{
    double m = 0.0;          // mass (or inertia) along the PTO
    double A_inf = 0.0;      // infinite-frequency added mass
    double C = 0.0;          // hydrostatic (or other) stiffness
    double B = 0.0;          // linear damping (viscous, drive-train loss)
    // radiation memory: xi' = Ar xi + Br v, memory force -Cr xi (real state space, may be empty)
    MatX Ar;
    VecX Br, Cr;
    // excitation force along the PTO at time t (the prediction)
    std::function<double(double)> Fe;
};

struct mpc_settings
{
    double Tc = 0.1;         // control period [s]
    int N = 40;              // horizon (periods)
    double r = 0.0;          // weight on the force: J = sum (F v_mean + r/2 F^2) Tc (0.01 Tc/(m + A_inf) is always added)
    double F_max = 0.0;      // force limit (0: none)
    double q_max = 0.0;      // stroke limit (soft, 0: none)
    double w_q = 0.0;        // weight of the stroke penalty (0: automatic)
};

class pto_ctrl_mpc : public pto_controller
{
public:
    pto_ctrl_mpc(const mpc_model&, const mpc_settings&);
    void update(const pto_ctrl_in&, pto_composite&) override;
    const char* name() const override {return "mpc";}
    void save_state(std::vector<double>&) const override;
    void load_state(const std::vector<double>&, size_t&) override;
    void channel_names(const std::string&, std::vector<std::string>&) const override;
    void channel_values(std::vector<double>&) const override;

    // the optimal force sequence for the state x, v, radiation states xi at time t (exposed for tests)
    VecX solve(double t, double x, double v, const VecX &xi);
    double force() const {return u_;}
    int iterations() const {return iters_;}

private:
    void observe(double t, double v);

    mpc_model mod_;
    mpc_settings set_;
    int nx_ = 2;
    MatX Phi_, Phir_;             // transition over Tc (whole model), radiation over the last step
    VecX Gu_, Ge_, Gr_;
    double dt_r_ = -1.0;
    VecX xi_;                     // radiation states (observer driven by the measured velocity)
    double t_obs_ = 0.0, v_obs_ = 0.0;
    bool started_ = false;
    double t_next_ = -1.0e300, u_ = 0.0;
    VecX useq_;
    int iters_ = 0;
};

// exp(A) of a small dense matrix (scaling and squaring with a Taylor series)
MatX coral_expm(const MatX&);

}

#endif
