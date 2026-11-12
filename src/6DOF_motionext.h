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

protected:
    int dof[6] = {1,1,1,1,1,1};

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
