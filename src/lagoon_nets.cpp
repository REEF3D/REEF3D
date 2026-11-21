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
// the output being gathered: the nets that printed it so far, by their number
struct gathered
{
    int num = -1;
    double time = 0.0;
    std::map<int, lagoon_nets::knots> nets;
};

gathered pending;
std::unique_ptr<lagoon_particles> writer;
bool off = false;
bool dynamic = false;  // a dynamic net in the run: Velocity and Acceleration as well

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

void start(lexer *p)
{
    const std::string solver = solver_of(p);
    std::string run;
    if(p->plog)
        run = "{\"type\": \"run\", \"run\": " + lagoon_store::json_string(p->plog->id()) + "}";
    std::vector<lagoon_particles::field> fields = {{"Tension", 1, false}};
    const int nets = bar_nets(p, dynamic);
    if(dynamic)
    {
        fields.push_back({"Velocity", 3, false});
        fields.push_back({"Acceleration", 3, false});
    }
    std::string key = solver;
    for(char &c : key)
        c = char(std::tolower((unsigned char)c));
    writer = std::make_unique<lagoon_particles>("./REEF3D_" + solver + ".lagoon", solver, key + "_net", "net",
                                                "REEF3D_" + solver + "_6DOF_Net", run, fields, 1, "lines",
                                                std::vector<lagoon_particles::field>{{"net", 1, true}});
    std::cout<<"LAGOON: "<<nets<<" net"<<(nets==1 ? "" : "s")<<" in the LAGOON store (particles/"<<key<<"_net)"<<std::endl;
}

// the gathered output into the store: the nets side by side, by their number
void write(lexer *p)
{
    if(pending.num<0 || pending.nets.empty())
    {
        pending = gathered();
        return;
    }
    if(!writer)
        start(p);
    std::vector<float> xyz, tension, velocity, acceleration;
    std::vector<int32_t> connectivity, offsets, net;
    const float nan = std::numeric_limits<float>::quiet_NaN();
    for(auto &[number, k] : pending.nets)
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
            net.push_back(int32_t(number));
        }
    }
    std::vector<const void*> values = {tension.data()};
    if(dynamic)
    {
        values.push_back(velocity.data());
        values.push_back(acceleration.data());
    }
    if(!writer->output(pending.time, pending.num, xyz.size()/3, xyz.data(), values, connectivity, offsets, {net.data()}))
    {
        off = true;
        std::cout<<"LAGOON: WARNING the nets could not be written to the LAGOON store, no more outputs of them"
                 <<(p->P18==1 ? " (the REEF3D-Net .vtk files are written again)" : "")<<std::endl;
    }
    pending = gathered();
}
}

bool lagoon_nets::output(lexer *p, int net, int num, knots &&k)
{
    if(p->P18<=0 || off)
        return false;
    if(pending.num>=0 && num!=pending.num)
        write(p);
    if(off)
        return false;
    pending.num = num;
    pending.time = p->simtime;
    pending.nets[net] = std::move(k);  // printed again: its later state
    return true;
}

void lagoon_nets::finish(lexer *p)
{
    if(p->mpirank==0 && p->P18>0 && !off)
        write(p);
}
