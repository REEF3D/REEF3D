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

#ifndef SIXDOF_MOTIONEXT_H_
#define SIXDOF_MOTIONEXT_H_

class lexer;
class fdm;
class fdm_nhf;
class fdm2D;
class ghostcell;
class field;
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

using namespace std;

class sixdof_motionext
{
public:
    virtual void motionext_trans(lexer*, ghostcell*, Eigen::Vector3d&, Eigen::Vector3d&)=0;
    virtual void motionext_rot(lexer*, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector4d&, Eigen::Matrix<double, 3, 4>&,  Eigen::Matrix3d&)=0;

    virtual void ini(lexer*,ghostcell*)=0;
    // checkpoints (P 48, I 45): the time index in the motion file and the prescribed velocities
    virtual void ckpt_state(std::vector<double>&) const {}
    virtual void ckpt_restore(const std::vector<double>&, size_t&) {}

    // DOF modes of the body (X 11, X 111 for single bodies): 2 is prescribed by this motion
    void set_dof(const int *d) { for(int n=0; n<6; ++n) dof[n]=d[n]; }

    // the stage of the body integrator (sixdof_obj::motion_stage_begin, CFD, NHFLOW, SFLOW and FNPF): the
    // time it ends, with the pose of the body (CoG c, quaternion e) at its start.  A motion table (the motion
    // files X 240 1 and 11) takes its velocities there; not set: at p->simtime, the start of the step
    void stage_begin(lexer*, double t, const Eigen::Vector3d &c, const Eigen::Vector4d &e, bool twoD);
    // the prescribed translations and orientation at that time into the pose of the body (motion tables)
    void stage_position(lexer*, Eigen::Vector3d &c, Eigen::Vector4d &e);

protected:
    int dof[6] = {1,1,1,1,1,1};
    double tstage = -1.0e20;

    // motion table (6DOF_motionext_table.cpp): rows t, columns; the column of each translation x, y, z and
    // each rotation angle phi, theta, psi [rad] (-1: none, the DOF keeps its pose when prescribed)
    std::vector<std::vector<double>> data;
    double ts=0.0, te=0.0;
    int tcol[3] = {-1,-1,-1};
    int rcol[3] = {-1,-1,-1};

    double table_vel(int col, double t);
    double table_pos(int col, double t);
    // the stage time is set and inside the table, the offsets of the body are taken, no velocity ramp (X 206)
    bool table_stage() const;
    // the orientation from the table: no free rotational DOF and at least one prescribed
    bool table_rotation() const;
    // the angular velocity of the body (body frame) at the stage time from the rates of the Euler angles
    void table_omega(Eigen::Vector3d &omega);

private:
    int ktab = 0;               // interval of the table at the last time asked
    double c0[3] = {0.0,0.0,0.0};   // the body at the table: c = c0 + table, angles = a0 + table
    double a0[3] = {0.0,0.0,0.0};
    bool offset_set = false;
    bool twoD_ = false;
    bool ramp_ = false;

    bool table_rot_col(int q) const;
    void table_angles(double t, double *ang);

protected:

    inline std::vector<std::vector<double>> readTable(const std::string& path, size_t M) {
        std::ifstream file(path);
        if (!file) {
            throw std::runtime_error("Could not open file: " + path);
        }

        std::vector<std::vector<double>> table_data;
        std::string line;

        while (std::getline(file, line)) {
            if (line.empty()) continue; // skip blank lines

            std::istringstream iss(line);
            std::vector<double> row;
            row.reserve(M);
            double value;
            while (iss >> value) {
                row.push_back(value);
            }

            if (row.size() != M) {
                throw std::runtime_error("Row has " + std::to_string(row.size()) +
                                          " columns, expected " + std::to_string(M));
            }
            table_data.push_back(std::move(row));
        }

        return table_data;
    }
};

#endif
