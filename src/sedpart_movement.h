/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#ifndef SEDPART_MOVEMENT_H_
#define SEDPART_MOVEMENT_H_

#include "increment.h"
#include "slice4.h"

#include <stdio.h>
#include <fstream>

class lexer;
class fdm;
class ghostcell;
class field4a;
class sediment_fdm;
class turbulence;

class sedpart;
class particles_obj;

namespace sediment_particle
{
    namespace movement
    {
        enum seedReturn:int
        {
            REMOVE=-1,STOP,CONTINUE
        };

        class base
        /// Base class for sediment particle movement models
        {
        public:
            virtual void setup(lexer *, fdm &, double &){};
            virtual seedReturn seeding(lexer *, particles_obj &, size_t &, double, bool=false){return seedReturn::STOP;};
            virtual void transfer(lexer *, particles_obj &, size_t &){};
            virtual void remove(lexer *, particles_obj &, size_t &){};
            virtual void move(lexer *, fdm &, ghostcell &, particles_obj &, sediment_fdm &, turbulence &){};
            virtual void make_moving(lexer *, fdm &, particles_obj &){};
            virtual void erode(lexer *, fdm &, particles_obj &, sediment_fdm &){};
            virtual void deposit(lexer *, fdm &, particles_obj &, sediment_fdm &){};
            virtual void update(lexer *, fdm &, ghostcell &, particles_obj &){};
            virtual void debug(lexer *, fdm &, ghostcell &, particles_obj &, sediment_fdm &){};
            virtual double volume(lexer *, fdm &, particles_obj &){return 0;};
            virtual void writeState(lexer *, ofstream &){};
            virtual void readState(lexer *, ifstream &){};
            virtual void setupState(lexer *, fdm &, ghostcell &, particles_obj &){};
        protected:
            double drag_coefficient(double) const;
        };
        class particleStressBased_T2021 : public base, increment
        /// Model for the movement of sediment particles following Tavouktsoglou et al. (2021)
        /// @author Alexander Hanke
        /// @date 2024
        {
        public:
            particleStressBased_T2021(lexer *);
            ~particleStressBased_T2021();

            void setup(lexer *, fdm &, double &);
            seedReturn seeding(lexer *, particles_obj &, size_t &, double, bool=false);
            void transfer(lexer *, particles_obj &, size_t &);
            void remove(lexer *, particles_obj &, size_t &);
            void move(lexer *, fdm &, ghostcell &, particles_obj &, sediment_fdm &, turbulence &);
            void make_moving(lexer *, fdm &, particles_obj &);
            void erode(lexer *, fdm &, particles_obj &, sediment_fdm &);
            void deposit(lexer *, fdm &, particles_obj &, sediment_fdm &);
            void update(lexer *, fdm &, ghostcell &, particles_obj &);
            void debug(lexer *, fdm &, ghostcell &, particles_obj &, sediment_fdm &);
            double volume(lexer *, fdm &, particles_obj &);
            void writeState(lexer *, ofstream &);
            void readState(lexer *, ifstream &);
            void setupState(lexer *, fdm &, ghostcell &, particles_obj &);
        private:
            double maxParticlesPerCell(lexer *, fdm &, double,bool=true,bool=false);
            void particleStressTensor(lexer *, fdm &, ghostcell &, particles_obj &);
            void particleStressTensorUpdateIJK(lexer *, fdm &, particles_obj &);
            void updateParticleStressTensor(lexer *, fdm &, particles_obj &, int, int, int);
            double theta_s(lexer *, fdm &, particles_obj &, int, int, int) const;
            double drag_model(lexer *, double, double, double, double, double) const;
            double sedimentation_velocity(lexer *, double, double, double, double, double) const;
            void particlePerCell(lexer *, ghostcell &, particles_obj &);
            void bedReDistribution(lexer *, fdm &, ghostcell &, particles_obj &);
            void timestep(lexer *, ghostcell &, particles_obj &);
            int activateNew(lexer *, fdm &, particles_obj &);
        private:
            /// @brief Sum of particles belonging to the stationary bed
            double *cellSumTopo;
            /// @brief Sum of particles belonging to the mobile bed
            double *cellSum;
            /// @brief Stress tensor for the particle-particle interaction
            double *stressTensor;
            /// @brief Number of particles per 2D bed column
            double *columnSum;
            /// @brief Relative density of fluid and particle
            const double drho;
            /// @brief Inverse of kinetik viscosity of the fluid
            const double invKinVis;
            /// @brief Stress tensor parameter
            const double Ps;
            /// @brief Stress tensor parameter
            const double beta;
            /// @brief Stress tensor parameter
            const double epsilon;
            /// @brief Critical solid volume fraction
            const double theta_crit;

            double dx;
            slice4 bedChange;

            double velDist=2.5;
        };

        // class doi10_1002_wrcr_20303 : public base, increment
        // {
        //     void move(lexer *, fdm *, particles_obj &, sediment_fdm &);
        //     

        //     slice4 phi_old;
        //     slice4 bedParticleNumber;
        // };
    };
};

#endif