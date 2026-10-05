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
#include<cctype>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<fstream>
#include<iostream>
#include<iterator>
#include<limits>
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

// a structured grid: no moving or tilted grids (B 180-192 change x and z per point)
bool fixed_grid(lexer *p)
{
    if(p->B180>0 || p->B191>0 || p->B192>0)
    {
        if(p->mpirank==0)
            std::cout<<"LAGOON: P 18 needs a fixed grid (no B 180, B 191, B 192); no LAGOON store written"<<std::endl;
        return false;
    }
    return true;
}

// the global grid lines: each rank gives its own, an elementwise maximum joins them
void grid_lines(lexer *p, ghostcell *pgc, std::vector<double> &x, std::vector<double> &y)
{
    const int gnx = p->gknox + 1;
    const int gny = p->gknoy + 1;
    x.assign(gnx, -1.0e300);
    y.assign(gny, -1.0e300);
    const int marge = increment::marge;  // for IP1, JP1
    int i,j;
    for(i=-1; i<p->knox; ++i)
        x[p->origin_i + i + 1] = p->XN[IP1];
    for(j=-1; j<p->knoy; ++j)
        y[p->origin_j + j + 1] = p->YN[JP1];
    pgc->globalmax(x.data(), gnx);
    pgc->globalmax(y.data(), gny);
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
    return p->P18!=1 || failed;
}

bool lagoon_output::start(lexer *p, ghostcell *pgc, const std::vector<lagoon_store::variable> &fields)
{
    if(!fixed_grid(p))
        return false;
    nx = p->knox + 1;
    ny = p->knoy + 1;
    nz = p->knoz + 1;
    std::vector<double> x, y;
    grid_lines(p, pgc, x, y);
    int k;

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
        failed = true;
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
        if(!usable)
            failed = true;
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

// ======================================================================== surfaces
lagoon_surface::lagoon_surface(lexer *p, const char *solver_, const char *output_, const char *source_)
    : solver(solver_), output(output_), source(source_), store(store_path(solver_)),
      ready(false), usable(true), t(0), nx(0), ny(0), rank(p->mpirank)
{
}

void lagoon_surface::piece_written(lexer *p, ghostcell *pgc, lagoon_surface *&writer, const char *solver,
                                   const char *output, const char *source, const char *file, int num)
{
    if(p->P18<=0)
        return;
    if(writer==nullptr)
        writer = new lagoon_surface(p, solver, output, source);
    std::ifstream in(file, std::ios::binary);
    std::stringstream text;
    text << in.rdbuf();
    in.close();
    writer->vtp_piece(p, pgc, text.str(), num);
    if(p->P18==1 && writer->usable)  // the store has it: no VTP file
        std::remove(file);
}

bool lagoon_surface::start(lexer *p, ghostcell *pgc, const std::vector<lagoon_store::variable> &fields)
{
    if(!fixed_grid(p))
        return false;
    nx = p->knox + 1;
    ny = p->knoy + 1;
    std::vector<double> x, y;
    grid_lines(p, pgc, x, y);
    if(p->mpirank==0)
    {
        struct stat info;
        if(stat((store_path(solver.c_str()) + "/zarr.json").c_str(), &info)!=0)  // no volume output
        {
            std::string run;
            if(p->plog)
                run = "{\"type\": \"run\", \"run\": " + lagoon_store::json_string(p->plog->id()) + "}";
            store.create_root(solver, run);
        }
        store.create_output(output, "surface", x, y, std::vector<double>(), fields, p->M10>0 ? p->M10 : 1, source);
    }
    MPI_Barrier(pgc->mpi_comm);
    store.create_block(output, p->mpirank, p->origin_i, p->origin_j, nx, ny, 1, fields, p->mpirank);
    MPI_Barrier(pgc->mpi_comm);
    return true;
}

void lagoon_surface::vtp_piece(lexer *p, ghostcell *pgc, const std::string &buffer, int num)
{
    if(!usable)
        return;
    const size_t appended = buffer.find("<AppendedData");
    const size_t data_start = appended==std::string::npos ? std::string::npos : buffer.find('_', appended);
    std::vector<lagoon_store::vtu_array> parsed;
    long long points_offset = -1;
    bool readable = data_start!=std::string::npos
                    && lagoon_store::parse_vtu_header(buffer.substr(0, data_start), parsed, points_offset);
    long long npoints = -1;
    {
        const size_t at = buffer.find("NumberOfPoints=\"");
        if(at!=std::string::npos)
            npoints = std::atoll(buffer.c_str() + at + 16);
    }
    std::vector<lagoon_store::variable> fields;
    for(const lagoon_store::vtu_array &a : parsed)
        fields.push_back(a.var);
    if(!ready)
    {
        // the first output decides the variables; a piece that cannot be read stops it on every rank
        const bool fine = pgc->globalmax(readable ? 0.0 : 1.0) <= 0.0;
        usable = fine && start(p, pgc, fields);
        if(!usable)
            lagoon_output::failed = true;
        ready = true;
        if(!usable)
        {
            if(p->mpirank==0 && !fine)
                std::cout<<"LAGOON: the "<<output<<" VTP cannot be read for the store; no LAGOON "<<output<<std::endl;
            return;
        }
    }

    // the points are the grid nodes, x outermost: z and the arrays with x fastest
    const size_t n = size_t(nx)*ny;
    bool ok = readable && npoints==(long long)n && data_start + 1 + size_t(points_offset) + 4 + 12*n <= buffer.size();
    try
    {
        if(ok)
        {
            const char *data = buffer.data() + data_start + 1;
            std::vector<float> values;
            auto arranged = [&](long long offset, int components, int component)
            {
                values.assign(n, 0.0f);
                const float *from = reinterpret_cast<const float*>(data + offset + sizeof(int));
                for(int i=0; i<nx; ++i)
                    for(int j=0; j<ny; ++j)
                        values[size_t(j)*nx + i] = from[(size_t(i)*ny + j)*components + component];
            };
            arranged(points_offset, 3, 2);
            store.write(output, rank, t, "z", values.data());
            for(const lagoon_store::vtu_array &a : parsed)
            {
                if(data_start + 1 + size_t(a.offset) + 4 + size_t(4)*a.var.components*n > buffer.size())
                    throw std::runtime_error("the VTP piece is shorter than its arrays");
                std::vector<float> all(n*a.var.components);
                for(int c=0; c<a.var.components; ++c)
                {
                    arranged(a.offset, a.var.components, c);
                    for(size_t q=0; q<n; ++q)
                        all[q*a.var.components + c] = values[q];
                }
                store.write(output, rank, t, a.var.name, all.data());
            }
        }
    }
    catch(std::exception &error)
    {
        std::cout<<"LAGOON rank "<<rank<<": "<<error.what()<<std::endl;
        ok = false;
    }
    // every block has the output (also a barrier): only then is it counted
    if(pgc->globalmax(ok ? 0.0 : 1.0) > 0.0)
    {
        if(p->mpirank==0)
            std::cout<<"LAGOON: a rank could not write its "<<output<<" block; no more LAGOON "<<output<<std::endl;
        usable = false;
        lagoon_output::failed = true;
        return;
    }
    if(p->mpirank==0)
        store.commit(output, t, p->simtime, num);
    ++t;
}

// ======================================================================== AMR surfaces
lagoon_amr_output::lagoon_amr_output(const char *solver_, const std::vector<lagoon_amr::field> &fields_)
    : solver(solver_), fields(fields_), writer(nullptr), usable(true)
{
}

bool lagoon_amr_output::files_needed(lexer *p, bool stored)
{
    return !(stored && p->P18==1);
}

bool lagoon_amr_output::write(lexer *p, ghostcell *pgc, const std::vector<lagoon_amr::grid> &grids, int printcount)
{
    if(!usable || lagoon_output::failed)  // the same on every rank
        return false;

    int ok = 1;
    std::vector<double> mine;
    try
    {
        lagoon_amr::pack(grids, int(fields.size()), mine);
    }
    catch(std::exception &problem)
    {
        std::cout<<"LAGOON: "<<problem.what()<<std::endl;
        mine.clear();
        ok = 0;
    }
    if(mine.size() > size_t(std::numeric_limits<int>::max()))
    {
        mine.clear();
        ok = 0;
    }
    int ok_all = 0;
    MPI_Allreduce(&ok,&ok_all,1,MPI_INT,MPI_MIN,pgc->mpi_comm);

    // every rank's grids to rank 0
    int count = int(mine.size());
    std::vector<int> counts(p->mpi_size,0), starts(p->mpi_size,0);
    MPI_Gather(&count,1,MPI_INT,counts.data(),1,MPI_INT,0,pgc->mpi_comm);
    std::vector<double> all;
    if(p->mpirank==0)
    {
        long long total = 0;
        for(int r=0; r<p->mpi_size; ++r)
        {
            starts[r] = int(total);
            total += counts[r];
        }
        if(total > (long long)std::numeric_limits<int>::max())
            ok_all = 0;
        else
            all.resize(size_t(total));
    }
    MPI_Bcast(&ok_all,1,MPI_INT,0,pgc->mpi_comm);
    if(ok_all==1)
        MPI_Gatherv(mine.data(),count,MPI_DOUBLE,all.data(),counts.data(),starts.data(),MPI_DOUBLE,0,pgc->mpi_comm);

    if(p->mpirank==0 && ok_all==1)
    {
        // level 0 of every rank, then the patches, rank by rank (as the .vtm lists them)
        std::vector<lagoon_amr::grid> level0, patches;
        for(int r=0; r<p->mpi_size && ok_all==1; ++r)
        {
            std::vector<lagoon_amr::grid> one;
            if(!lagoon_amr::unpack(all.data() + starts[r], size_t(counts[r]), int(fields.size()), r, one))
            {
                std::cout<<"LAGOON: the AMR grids of rank "<<r<<" did not arrive whole"<<std::endl;
                ok_all = 0;
                break;
            }
            for(lagoon_amr::grid &g : one)
                (g.level==0 ? level0 : patches).push_back(std::move(g));
        }
        if(ok_all==1)
        {
            if(writer==nullptr)
            {
                std::string run;
                if(p->plog)
                    run = "{\"type\": \"run\", \"run\": " + lagoon_store::json_string(p->plog->id()) + "}";
                std::string key = solver;
                for(char &c : key)
                    c = char(std::tolower((unsigned char)c));
                writer = new lagoon_amr(store_path(solver.c_str()), solver, key + "_amr",
                                        "REEF3D_" + solver + "_AMR", run, fields);
            }
            level0.insert(level0.end(), std::make_move_iterator(patches.begin()), std::make_move_iterator(patches.end()));
            ok_all = writer->output(p->simtime, printcount, level0) ? 1 : 0;
        }
    }
    MPI_Bcast(&ok_all,1,MPI_INT,0,pgc->mpi_comm);
    if(ok_all!=1)
    {
        usable = false;
        lagoon_output::failed = true;
        if(p->mpirank==0)
            std::cout<<"LAGOON: no more "<<solver<<" AMR output in the store; the .vtr and .vtm files are written"<<std::endl;
    }
    return ok_all==1;
}
