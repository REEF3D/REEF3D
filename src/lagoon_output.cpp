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

#include"lagoon_output.h"
#include"lexer.h"
#include"ghostcell.h"
#include"runlog.h"
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<fstream>
#include<iostream>
#include<sstream>
#include<sys/stat.h>

namespace
{
std::string store_path(const char *solver)
{
    return std::string("./REEF3D_") + solver + ".lagoon";
}

// the run id in a store's root metadata, or "" if there is none
std::string stored_run(const std::string &path)
{
    std::ifstream in((path + "/zarr.json").c_str());
    if(!in)
        return "";
    std::stringstream text;
    text << in.rdbuf();
    const std::string s = text.str();
    const std::string key = "\"run\": \"";
    const size_t at = s.find(key, s.find("\"reef3d_run\""));
    if(at==std::string::npos)
        return "earlier";
    const size_t end = s.find('"', at + key.size());
    return s.substr(at + key.size(), end - at - key.size());
}
}

lagoon_output::lagoon_output(lexer *p, ghostcell *pgc, const char *solver_)
    : solver(solver_), store(store_path(solver_)), ready(false), usable(true), cartesian(false), t(0), nx(0), ny(0), nz(0), rank(p->mpirank)
{
    // a store holds one run: the store of an earlier run in this folder is moved aside
    if(p->mpirank==0)
    {
        const std::string path = store_path(solver_);
        struct stat info;
        if(stat(path.c_str(), &info)==0)
        {
            const std::string aside = std::string("./REEF3D_") + solver_ + "-" + stored_run(path) + ".lagoon-earlier";
            if(std::rename(path.c_str(), aside.c_str())==0)
                std::cout<<"LAGOON: the store of the earlier run is now "<<aside<<std::endl;
        }
        std::string run;
        if(p->plog)
            run = "{\"type\": \"run\", \"run\": " + lagoon_store::json_string(p->plog->id()) + "}";
        store.create_root(solver, run);
    }
    MPI_Barrier(pgc->mpi_comm);
}

bool lagoon_output::vtu_files(lexer *p)
{
    return p->P18!=2;
}

bool lagoon_output::start(lexer *p, ghostcell *pgc, const std::vector<lagoon_store::variable> &fields)
{
    // a structured grid: no moving or tilted grids (B 180-192 change x and z per point)
    if(p->B180>0 || p->B191>0 || p->B192>0)
    {
        if(p->mpirank==0)
            std::cout<<"LAGOON: P 18 needs a fixed grid (no B 180, B 191, B 192); no LAGOON store written"<<std::endl;
        return false;
    }
    nx = p->knox + 1;
    ny = p->knoy + 1;
    nz = p->knoz + 1;
    const int gnx = p->gknox + 1;
    const int gny = p->gknoy + 1;

    // the global grid lines: each rank gives its own, an elementwise maximum joins them
    std::vector<double> x(gnx, -1.0e300), y(gny, -1.0e300);
    int i,j,k;
    for(i=-1; i<p->knox; ++i)
        x[p->origin_i + i + 1] = p->XN[IP1];
    for(j=-1; j<p->knoy; ++j)
        y[p->origin_j + j + 1] = p->YN[JP1];
    pgc->globalmax(x.data(), gnx);
    pgc->globalmax(y.data(), gny);

    // FNPF and NHFLOW: σ-levels, the same in every rank; CFD: Cartesian levels at
    // fixed heights, the grid split in z as well
    cartesian = solver=="CFD";
    std::vector<double> levels;
    if(cartesian)
    {
        levels.assign(p->gknoz + 1, -1.0e300);
        for(k=-1; k<p->knoz; ++k)
            levels[p->origin_k + k + 1] = p->ZN[KP1];
        pgc->globalmax(levels.data(), int(levels.size()));
    }
    else
    {
        sigma.resize(nz);
        for(k=-1; k<p->knoz; ++k)
            sigma[k+1] = p->ZN[KP1];
        levels = sigma;
    }

    if(p->mpirank==0)
        store.create_output("volume", cartesian ? "cartesian" : "sigma", x, y, levels, fields,
                            p->M10>0 ? p->M10 : 1, std::string("REEF3D_") + solver + "_VTU");
    store.create_block("volume", p->mpirank, p->origin_i, p->origin_j, nx, ny, nz, fields, p->mpirank,
                       cartesian, cartesian ? p->origin_k : 0);
    MPI_Barrier(pgc->mpi_comm);
    return true;
}

void lagoon_output::write_job(job *j)
{
    try
    {
        if(!j->z.empty())  // σ-grids: the heights of the levels in every column
        {
            const size_t column = size_t(nx)*ny;
            std::vector<float> offsets;
            lagoon_store::level_offsets(j->z.data(), nz, int(column), sigma, offsets);
            store.write("volume", rank, j->t, "z_bed", &j->z[0]);
            store.write("volume", rank, j->t, "z_surface", &j->z[(nz-1)*column]);
            store.write("volume", rank, j->t, "z_offset", offsets.data());
        }
        for(const std::pair<std::string, std::vector<float> > &f : j->fields)
            store.write("volume", rank, j->t, f.first, f.second.data());
    }
    catch(std::exception &error)
    {
        std::cout<<"LAGOON rank "<<rank<<": "<<error.what()<<std::endl;
        j->ok = false;
    }
}

bool lagoon_output::settle(lexer *p, ghostcell *pgc)
{
    if(pending.t<0)
        return true;
    if(worker.joinable())
        worker.join();
    // every block has the output (also a barrier): only then is it counted
    if(pgc->globalmax(pending.ok ? 0.0 : 1.0) > 0.0)
    {
        if(p->mpirank==0)
            std::cout<<"LAGOON: a rank could not write its block; no more LAGOON output"<<std::endl;
        usable = false;
        pending = job();
        return false;
    }
    if(p->mpirank==0)
        store.commit("volume", pending.t, pending.time, pending.num);
    pending = job();
    return true;
}

void lagoon_output::finish(lexer *p, ghostcell *pgc)
{
    if(usable)
        settle(p, pgc);
}

void lagoon_output::vtu_piece(lexer *p, ghostcell *pgc, const std::vector<char> &buffer, size_t data_start, int num)
{
    if(!usable)
        return;
    std::vector<lagoon_store::vtu_array> parsed;
    long long points_offset = -1;
    const bool readable = lagoon_store::parse_vtu_header(std::string(buffer.data(), data_start), parsed, points_offset);
    std::vector<lagoon_store::variable> fields;
    for(const lagoon_store::vtu_array &a : parsed)
        fields.push_back(a.var);
    if(!ready)
    {
        usable = start(p, pgc, fields);
        ready = true;
        if(!usable)
            return;
    }
    // the previous output: written by every rank?
    if(!settle(p, pgc))
        return;

    // copy what this output needs, then write it in the background
    const size_t n = size_t(nx)*ny*nz;
    pending.t = t;
    pending.num = num;
    pending.time = p->simtime;
    pending.ok = readable && size_t(p->pointnum)==n;
    if(pending.ok)
    {
        if(!cartesian)
        {
            const float *points = reinterpret_cast<const float*>(&buffer[data_start + points_offset + sizeof(int)]);
            pending.z.resize(n);
            for(size_t q=0; q<n; ++q)
                pending.z[q] = points[3*q+2];
        }
        for(const lagoon_store::vtu_array &a : parsed)
        {
            const float *values = reinterpret_cast<const float*>(&buffer[data_start + a.offset + sizeof(int)]);
            pending.fields.emplace_back(a.var.name, std::vector<float>(values, values + n*a.var.components));
        }
        worker = std::thread(&lagoon_output::write_job, this, &pending);
    }
    ++t;
}
