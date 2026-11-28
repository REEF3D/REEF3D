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

#include"lagoon_nets.h"
#include"lagoon_store.h"
#include"lexer.h"
#include"runlog.h"
#include<cctype>
#include<iostream>
#include<limits>
#include<map>
#include<memory>
#include<string>

namespace
{
std::string solver_of(lexer *p)
{
    return p->A10==6 ? "CFD" : p->A10==3 ? "FNPF" : p->A10==2 ? "SFLOW" : "NHFLOW";
}

// the bar nets of the run (X 320 types 1 - 3 quasi-static, 12 - 14 dynamic; sheets are no bars)
int bar_nets(lexer *p, bool &any_dynamic)
{
    int n = 0;
    any_dynamic = false;
    for(int i=0; i<p->net_count; ++i)
    {
        const int type = p->X320_type[i];
        if(type>10)
            any_dynamic = true;
        if(type>10 || type<4)
            ++n;
    }
    return n;
}

// a particle set of lines (the nets, the mooring lines): the output being gathered,
// by the number of each net or line that printed it so far, and its writer
class line_set
{
public:
    line_set(const char *what, const char *role, const char *folder, const char *cell_field)
        : what(what), role(role), folder(folder), cell_field(cell_field) {}

    bool output(lexer *p, int number, int num, lagoon_line_points &&k)
    {
        if(p->P18<=0 || off)
            return false;
        if(num_pending>=0 && num!=num_pending)
            write(p);
        if(off)
            return false;
        num_pending = num;
        time_pending = p->simtime;
        pending[number] = std::move(k);  // printed again: its later state
        return true;
    }

    void finish(lexer *p)
    {
        if(p->mpirank==0 && p->P18>0 && !off)
            write(p);
    }

private:
    std::string what, role, folder, cell_field;
    int num_pending = -1;
    double time_pending = 0.0;
    std::map<int, lagoon_line_points> pending;
    std::unique_ptr<lagoon_particles> writer;
    bool off = false;
    bool dynamic = false;  // Velocity and Acceleration as well (a dynamic net in the run)

    void start(lexer *p)
    {
        const std::string solver = solver_of(p);
        std::string run;
        if(p->plog)
            run = "{\"type\": \"run\", \"run\": " + lagoon_store::json_string(p->plog->id()) + "}";
        std::vector<lagoon_particles::field> fields = {{"Tension", 1, false}};
        int count = int(pending.size());
        if(role=="net")
            count = bar_nets(p, dynamic);
        if(dynamic)
        {
            fields.push_back({"Velocity", 3, false});
            fields.push_back({"Acceleration", 3, false});
        }
        std::string key = solver;
        for(char &c : key)
            c = char(std::tolower((unsigned char)c));
        key += "_" + role;
        writer = std::make_unique<lagoon_particles>("./REEF3D_" + solver + ".lagoon", solver, key, role,
                                                    "REEF3D_" + solver + folder, run, fields, 1, "lines",
                                                    std::vector<lagoon_particles::field>{{cell_field, 1, true}});
        if(role=="net")
            std::cout<<"LAGOON: "<<count<<" net"<<(count==1 ? "" : "s");
        else
            std::cout<<"LAGOON: the "<<what;
        std::cout<<" in the LAGOON store (particles/"<<key<<")"<<std::endl;
    }

    // the gathered output into the store: the nets or lines side by side, by their number
    void write(lexer *p)
    {
        if(num_pending<0 || pending.empty())
        {
            reset();
            return;
        }
        if(!writer)
            start(p);
        std::vector<float> xyz, tension, velocity, acceleration;
        std::vector<int32_t> connectivity, offsets, numbers;
        const float nan = std::numeric_limits<float>::quiet_NaN();
        for(auto &[number, k] : pending)
        {
            const int32_t first = int32_t(xyz.size()/3);
            const size_t n = k.xyz.size()/3;
            xyz.insert(xyz.end(), k.xyz.begin(), k.xyz.end());
            tension.insert(tension.end(), k.tension.begin(), k.tension.end());
            tension.resize(xyz.size()/3, nan);
            if(dynamic)  // a quasi-static net among dynamic ones: no velocity, no acceleration
            {
                if(k.velocity.size()==3*n)
                    velocity.insert(velocity.end(), k.velocity.begin(), k.velocity.end());
                if(k.acceleration.size()==3*n)
                    acceleration.insert(acceleration.end(), k.acceleration.begin(), k.acceleration.end());
                velocity.resize(xyz.size(), nan);
                acceleration.resize(xyz.size(), nan);
            }
            for(size_t b=0; b+1<k.bars.size(); b+=2)
            {
                connectivity.push_back(first + k.bars[b]);
                connectivity.push_back(first + k.bars[b+1]);
                offsets.push_back(int32_t(connectivity.size()));
                numbers.push_back(int32_t(number));
            }
        }
        std::vector<const void*> values = {tension.data()};
        if(dynamic)
        {
            values.push_back(velocity.data());
            values.push_back(acceleration.data());
        }
        if(!writer->output(time_pending, num_pending, xyz.size()/3, xyz.data(), values, connectivity, offsets, {numbers.data()}))
        {
            off = true;
            std::cout<<"LAGOON: WARNING the "<<what<<" could not be written to the LAGOON store, no more outputs of them"
                     <<(p->P18==1 ? " (their .vtk files are written again)" : "")<<std::endl;
        }
        reset();
    }

    void reset()
    {
        num_pending = -1;
        time_pending = 0.0;
        pending.clear();
    }
};

line_set nets("nets", "net", "_6DOF_Net", "net");
line_set moorings("mooring lines", "mooring", "_6DOF_Mooring", "line");
}

bool lagoon_nets::output(lexer *p, int net, int num, knots &&k)
{
    return nets.output(p, net, num, std::move(k));
}

void lagoon_nets::finish(lexer *p)
{
    nets.finish(p);
}

bool lagoon_moorings::output(lexer *p, int line, int num, lagoon_line_points &&points)
{
    points.velocity.clear();
    points.acceleration.clear();
    return moorings.output(p, line, num, std::move(points));
}

void lagoon_moorings::finish(lexer *p)
{
    moorings.finish(p);
}
