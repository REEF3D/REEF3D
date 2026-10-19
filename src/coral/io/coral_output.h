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
// REEF3D::CORAL: time-series output. One line per output time: t, per body the CoG position, the
// Euler angles (roll, pitch, yaw in degrees), the CoG velocity and the angular velocity (global), then
// the output channels of all loads. The header line names the columns.

#ifndef CORAL_OUTPUT_H_
#define CORAL_OUTPUT_H_

#include<fstream>
#include<string>
#include<vector>

namespace coral
{
    class model;

    class time_series
    {
    public:
        void open(const std::string &file, const model&, bool append=false);   // append: continue the file (restart)
        void write(const model&, double t);
        void close() {out_.close();}
        void flush() {out_.flush();}
        const std::vector<std::string> &columns() const {return cols_;}

    private:
        std::ofstream out_;
        std::vector<std::string> cols_;
    };
}

#endif
