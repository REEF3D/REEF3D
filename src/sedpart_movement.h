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

#include <stdio.h>
#include <fstream>

class lexer;
class fdm;
class ghostcell;
class field4a;

class sedpart;
class particles_obj;

namespace sediment_particle
{
    namespace movement
    {
        class base
        {
        public:
            virtual void setup(lexer *, fdm &, double &){};
            virtual bool seeding(lexer *, particles_obj &, size_t &, int &){return false;};
            virtual void transfer(lexer *, particles_obj &, size_t &){};
            virtual void remove(lexer *, particles_obj &, size_t &){};
            virtual void move(lexer *, fdm &, ghostcell &, particles_obj &){};
            virtual void update(lexer *, ghostcell &, field4a &, double &){};
            virtual void debug(lexer *, fdm &, ghostcell &, particles_obj &){};
            virtual void writeState(ofstream &){};
            virtual void readState(ifstream &){};
        };
        class Tavouktsoglou : public base, increment
        {
        public:
            Tavouktsoglou(lexer *);
            ~Tavouktsoglou();

            void setup(lexer *, fdm &, double &);
            bool seeding(lexer *, particles_obj &, size_t &, int &);
            void transfer(lexer *, particles_obj &, size_t &);
            void remove(lexer *, particles_obj &, size_t &);
            void move(lexer *, fdm &, ghostcell &, particles_obj &);
            void update(lexer *, ghostcell &, field4a &, double &);
            void debug(lexer *, fdm &, ghostcell &, particles_obj &);
            void writeState(ofstream &);
            void readState(ifstream &);
        private:
            double maxParticlesPerCell(lexer *, fdm &, double,bool=true,bool=false);
            void particleStressTensor(lexer *, fdm &, ghostcell &, particles_obj &);
            void particleStressTensorUpdateIJK(lexer *, fdm &, particles_obj &);
            void updateParticleStressTensor(lexer *, fdm &, particles_obj &, int, int, int);
            double theta_s(lexer *, fdm &, particles_obj &, int, int, int) const;
            double drag_model(lexer *, double, double, double, double, double) const;
            void particlePerCell(lexer *, ghostcell &, particles_obj &);
        private:
            double *cellSumTopo;
            double *cellSum;
            double *stressTensor;
            double *columnSum;
            const double drho;
            const double kinVis;
            const double Ps;
            const double beta;
            const double epsilon;
            const double theta_crit;
        };
    };
    class state : increment
    {
    public:
        int solid_clean(lexer *p, particles_obj &, sediment_particle::movement::base &);
    };
};

#endif