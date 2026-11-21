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

#include"lagoon_membranes.h"
#include"lagoon_store.h"
#include"lexer.h"
#include"runlog.h"
#include<iostream>
#include<map>
#include<memory>
#include<string>

namespace
{
std::map<int, lagoon_membranes::panels> latest;  // the last printed state of every membrane
int pending_count = -1;                          // the time step of the output being gathered
double pending_time = 0.0;
long long stored = 0;                            // outputs in the store
std::unique_ptr<lagoon_particles> writer;
bool off = false;

void start(lexer *p)
{
    std::string run;
    if(p->plog)
        run = "{\"type\": \"run\", \"run\": " + lagoon_store::json_string(p->plog->id()) + "}";
    using field = lagoon_particles::field;
    writer = std::make_unique<lagoon_particles>("./REEF3D_NHFLOW.lagoon", "NHFLOW", "nhflow_membrane", "membrane",
                                                "REEF3D_NHFLOW_Membrane_VTP", run,
                                                std::vector<field>{{"velocity", 3, false}, {"displacement", 3, false}, {"load", 3, false}},
                                                1, "polys",
                                                std::vector<field>{{"membrane", 1, true}, {"dp", 1, false}, {"force", 3, false},
                                                                   {"tension", 1, false}, {"tag", 1, true}});
    const size_t n = latest.size();
    std::cout<<"LAGOON: "<<n<<" membrane"<<(n==1 ? "" : "s")<<" in the LAGOON store (particles/nhflow_membrane)"<<std::endl;
}

template<typename T>
void append(std::vector<T> &to, const std::vector<T> &from)
{
    to.insert(to.end(), from.begin(), from.end());
}

// the gathered output into the store: the membranes side by side, by their number
void write(lexer *p)
{
    if(pending_count<0 || latest.empty())
        return;
    if(!writer)
        start(p);
    std::vector<float> xyz, velocity, displacement, load, dp, force, tension;
    std::vector<int32_t> connectivity, offsets, membrane, tag;
    for(const auto &[number, m] : latest)
    {
        const int32_t first = int32_t(xyz.size()/3);
        append(xyz, m.xyz);
        append(velocity, m.velocity);
        append(displacement, m.displacement);
        append(load, m.load);
        for(size_t t=0; t+2<m.triangles.size(); t+=3)
        {
            for(int q=0; q<3; ++q)
                connectivity.push_back(first + m.triangles[t+q]);
            offsets.push_back(int32_t(connectivity.size()));
            membrane.push_back(int32_t(number));
        }
        append(dp, m.dp);
        append(force, m.force);
        append(tension, m.tension);
        append(tag, m.tag);
    }
    if(!writer->output(pending_time, stored, xyz.size()/3, xyz.data(),
                       {velocity.data(), displacement.data(), load.data()}, connectivity, offsets,
                       {membrane.data(), dp.data(), force.data(), tension.data(), tag.data()}))
    {
        off = true;
        std::cout<<"LAGOON: WARNING the membranes could not be written to the LAGOON store, no more outputs of them"
                 <<(p->P18==1 ? " (the REEF3D-NHFLOW-Membrane .vtp files are written again)" : "")<<std::endl;
    }
    ++stored;
    pending_count = -1;
}
}

bool lagoon_membranes::output(lexer *p, int membrane, int num, panels &&m)
{
    (void)num;
    if(p->P18<=0 || off)
        return false;
    if(pending_count>=0 && p->count!=pending_count)
        write(p);
    if(off)
        return false;
    pending_count = p->count;
    pending_time = p->simtime;
    latest[membrane] = std::move(m);  // printed again in the same step: its later state
    return true;
}

void lagoon_membranes::finish(lexer *p)
{
    if(p->mpirank==0 && p->P18>0 && !off)
        write(p);
}
