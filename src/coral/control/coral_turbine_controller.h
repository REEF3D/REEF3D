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
// REEF3D::CORAL: wind turbine controller of the NREL 5 MW baseline type (Jonkman et al. 2009, chapter 7, the DISCON
// routine): variable-speed generator torque and gain-scheduled collective pitch on the filtered generator speed.
//
//   speed filter   first-order low pass, alpha = exp(-dt omega_c):  w_f = (1 - alpha) w + alpha w_f
//   torque         region 1   w_f <= w_cutin                          T = 0
//                  region 1.5 w_cutin < w_f < w_2                      T = slope15 (w_f - w_cutin)
//                  region 2   w_2 <= w_f < w_trans                     T = K w_f^2
//                  region 2.5 w_trans <= w_f < w_rated                 T = slope25 (w_f - w_sync)   (induction-machine line
//                             through rated torque at w_rated, w_sync = w_rated/(1 + slip/100))
//                  region 3   w_f >= w_rated or pitch >= pitch_rgn3   T = P_rated/w_f
//                  at most T_max, rate at most T_rate
//   pitch          GK = 1/(1 + theta/theta_K) (theta: the last command), e = w_f - w_ref,
//                  I = int e dt (held so that GK KI I stays within the pitch limits),
//                  theta_cmd = clamp(GK (KP e + KI I), theta_min, theta_max), rate at most theta_rate
// Speeds on the generator (high-speed) side in rad/s, torque on the generator side. The controller is pure logic: the
// model's measurement and actuators are attached by its load (controller_load) through functions, so it does not depend
// on the rotor's implementation. Gains for another rotor: coral/aero/coral_rotor_tune.h (from CORAL's own BEM).

#ifndef CORAL_TURBINE_CONTROLLER_H_
#define CORAL_TURBINE_CONTROLLER_H_

#include"coral/core/coral_load.h"
#include<functional>
#include<string>
#include<vector>

namespace coral
{
    struct turbine_controller_settings
    {
        double omega_c = 1.570796;        // filter corner [rad/s] (0.25 Hz)
        double dt = 0.0;                  // sample time [s] (0: every step)
        // torque (generator side)
        double w_cutin = 70.16224, w_rgn2 = 91.21091, K = 2.332287;      // rad/s, rad/s, N m s^2/rad^2
        double w_rated = 121.6805, P_rated = 5296610.0, slip = 10.0;     // rad/s, W (mechanical), %
        double T_max = 47402.91, T_rate = 15000.0;                        // N m, N m/s
        double pitch_rgn3 = 0.01745329;   // pitch above which region 3 [rad]
        // pitch
        double w_ref = 122.9096;          // rad/s
        double KP = 0.01882681, KI = 0.008068634, theta_K = 0.1099965;    // s, -, rad
        double theta_min = 0.0, theta_max = 1.570796, theta_rate = 0.1396263;   // rad, rad, rad/s
        bool pitch_control = true, torque_control = true;

        // the derived points of the torque curve
        double w_sync() const {return w_rated/(1.0 + 0.01*slip);}
        double slope15() const {return w_rgn2 > w_cutin ? K*w_rgn2*w_rgn2/(w_rgn2 - w_cutin) : 0.0;}
        double slope25() const;
        double w_trans() const;

        // PI gains of the pitch loop (Jonkman 2009, 7.2): a rigid drivetrain of inertia J (low-speed shaft, generator
        // N^2 J_gen included) at the rated rotor speed Omega0 with the power sensitivity dP/dtheta at zero pitch, a
        // closed-loop response of natural frequency omega_n and damping zeta: KP = 2 J Omega0 zeta omega_n/(N (-dP/dtheta)),
        // KI = J Omega0 omega_n^2/(N (-dP/dtheta)) (on the generator speed in rad/s)
        static void pi_gains(double J, double Omega0, double N, double dPdtheta0, double omega_n, double zeta, double &KP, double &KI);
    };

    class turbine_controller
    {
    public:
        explicit turbine_controller(const turbine_controller_settings &s = turbine_controller_settings()) : set(s) {}

        // start from the measured generator speed and pitch
        void initialize(double t, double w_gen, double pitch, double torque);
        // one controller step at t with the measured generator speed: new torque and pitch commands
        void update(double t, double w_gen);

        double torque_law(double w_f, double pitch) const;      // the torque before the limits
        double torque() const {return T_;}
        double pitch() const {return th_;}
        double filtered_speed() const {return wf_;}
        double integral() const {return I_;}
        int region() const {return region_;}

        void save_state(std::vector<double>&) const;
        bool load_state(const std::vector<double>&, size_t &k);

        turbine_controller_settings set;

    private:
        double t_last_ = 0.0, wf_ = 0.0, I_ = 0.0, T_ = 0.0, th_ = 0.0;
        int region_ = 0;
        bool started_ = false;
    };

    // the controller in the model: called once per converged step (commit), the commands held over the next step
    // (zero-order hold, as a controller DLL); channels <name>_wf (rpm), _T (N m), _pitch (deg), _region
    class controller_load : public load
    {
    public:
        struct io
        {
            std::function<double(const model&)> speed;            // generator speed [rad/s]
            std::function<double(const model&)> pitch;            // measured pitch [rad]
            std::function<double(const model&)> torque;           // the generator torque at the start [N m]
            std::function<void(double)> set_torque;
            std::function<void(double)> set_pitch;                // collective
        };
        controller_load(const turbine_controller_settings&, io, const std::string &name="controller");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void eval(const model&, double, VecX&) const override {}
        void initialize(const model&, double t0) override;
        void commit(const model&, double t, double dt) override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;
        void save_state(std::vector<double> &v) const override {ctrl_.save_state(v);}
        bool load_state(const std::vector<double> &v, size_t &k) override;

        const turbine_controller &controller() const {return ctrl_;}

    private:
        void apply();
        turbine_controller ctrl_;
        io io_;
        std::string name_;
    };
}

#endif
