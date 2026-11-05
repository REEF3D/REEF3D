// Standalone verification of the LAGOON particle writer (lagoon_particles): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../src lagoon_particles_test.cpp ../../src/lagoon_store.cpp -lz -o lagoon_particles_test
// Run:    ./lagoon_particles_test     (writes lagoon_particles_test.lagoon in the current folder)
// The files follow LAGOON's lagoon-format SPEC.md 0.2, section 10; LAGOON reads them
// with lagoon_format.store (Store(...).particles["nhflow_particles"]).
// Architect: Hans Bihs
#include"lagoon_store.h"
#include<cmath>
#include<cstdint>
#include<cstring>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<sys/stat.h>
#include<vector>
#include<zlib.h>

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

static std::string slurp(const std::string &file)
{
    std::ifstream in(file.c_str(), std::ios::binary);
    std::stringstream s;
    s<<in.rdbuf();
    return s.str();
}

static uint64_t u64(const unsigned char *p)
{
    uint64_t v=0;
    for(int b=7;b>=0;--b) v=(v<<8)|p[b];
    return v;
}

// inner chunk c of a shard file: gunzip and unshuffle (itemsize 4)
static std::vector<unsigned char> inner(const std::string &file, int c, size_t values)
{
    const std::string s = slurp(file);
    const size_t entries = lagoon_particles::CHUNKS_PER_SHARD;
    const unsigned char *index = reinterpret_cast<const unsigned char*>(s.data()) + s.size() - 16*entries - 4;
    const uint64_t offset = u64(index + 16*c), nbytes = u64(index + 16*c + 8);
    std::vector<unsigned char> raw(values*4), out(values*4);
    z_stream zs;
    std::memset(&zs, 0, sizeof(zs));
    inflateInit2(&zs, 15+16);
    zs.next_in = (unsigned char*)s.data() + offset;
    zs.avail_in = nbytes;
    zs.next_out = raw.data();
    zs.avail_out = raw.size();
    inflate(&zs, Z_FINISH);
    inflateEnd(&zs);
    for(size_t e=0; e<values; ++e)
        for(int b=0; b<4; ++b)
            out[4*e+b] = raw[b*values + e];
    return out;
}

static float position(long long global, int c) { return float(0.001*global + c); }

int main()
{
    std::cout<<"LAGOON particle writer"<<std::endl;
    const std::string path = "./lagoon_particles_test.lagoon";
    std::vector<lagoon_particles::field> fields = {{"velocity", 3, false}, {"id", 1, true}, {"fluid velocity", 3, false}};
    lagoon_particles writer(path, "NHFLOW", "nhflow_particles", "particles", "REEF3D_NHFLOW_Particles",
                            "{\"type\": \"run\", \"run\": \"test\"}", fields);
    // outputs of 1000, 70000, 0, 400000 and 700000 particles: inner chunks and shards are crossed
    const size_t counts[5] = {1000, 70000, 0, 400000, 700000};
    long long global = 0;
    bool all = true;
    for(int t=0; t<5; ++t)
    {
        const size_t n = counts[t];
        std::vector<float> xyz(3*n), velocity(3*n), fluid(3*n);
        std::vector<int32_t> id(n);
        for(size_t i=0; i<n; ++i)
        {
            for(int c=0; c<3; ++c)
            {
                xyz[3*i+c] = position(global + (long long)i, c);
                velocity[3*i+c] = float(t + 0.5*c);
                fluid[3*i+c] = -float(t);
            }
            id[i] = int32_t(global + (long long)i);
        }
        all = writer.output(0.25*t, t, n, xyz.data(), {velocity.data(), id.data(), fluid.data()}) && all;
        global += (long long)n;
    }
    check(all && writer.usable(), "5 outputs (0 to 700000 particles)");
    const std::string set = path + "/particles/nhflow_particles";
    const std::string attrs = slurp(set + "/zarr.json");
    check(attrs.find("\"fluid velocity\": {\"components\": 3, \"data_type\": \"float32\", \"array\": \"fluid_velocity\"}")!=std::string::npos,
          "field names and their arrays");
    check(slurp(set + "/time/zarr.json").find("\"shape\": [5]")!=std::string::npos, "time counts 5 outputs");
    check(slurp(set + "/position/zarr.json").find("\"shape\": [1171000, 3]")!=std::string::npos, "1171000 points");
    check(slurp(set + "/point_data/id/zarr.json").find("\"data_type\": \"int32\"")!=std::string::npos, "id as int32");
    // the second shard (rows 1048576 ...) holds the end of the last output
    struct stat info;
    check(stat((set + "/position/c/1/0").c_str(), &info)==0 && stat((set + "/point_data/id/c/1").c_str(), &info)==0,
          "a second shard");
    std::vector<unsigned char> chunk = inner(set + "/point_data/id/c/1", 1, lagoon_particles::ROWS);
    int32_t first;
    std::memcpy(&first, chunk.data(), 4);
    check(first == int32_t(1048576 + 65536), "ids read back from the second shard");
    chunk = inner(set + "/position/c/0/0", 2, 3*lagoon_particles::ROWS);
    float y;
    std::memcpy(&y, chunk.data() + 4*(3*5 + 1), 4);
    check(y == position(2*65536 + 5, 1), "positions read back");

    // objects (ice floes): two quads, then one breaks off; a value per cell
    {
        std::vector<lagoon_particles::field> none;
        std::vector<lagoon_particles::field> per_cell = {{"id", 1, true}, {"speed", 1, false}};
        lagoon_particles floes(path, "FNPF", "fnpf_ice", "ice", "REEF3D_FNPF_ICE", "", none, 1, "polys", per_cell);
        bool written = true;
        for(int t=0; t<6; ++t)
        {
            const int quads = t<4 ? 2 : 1;
            std::vector<float> xyz(12*quads);
            for(size_t i=0; i<xyz.size(); ++i)
                xyz[i] = float(t + 0.1*i);
            std::vector<int32_t> connectivity(4*quads), offsets(quads), id(quads);
            std::vector<float> speed(quads, float(t));
            for(int q=0; q<quads; ++q)
            {
                for(int k=0; k<4; ++k)
                    connectivity[4*q+k] = 4*q + k;
                offsets[q] = 4*(q+1);
                id[q] = 10 + q;
            }
            written = floes.output(t, t, 4*size_t(quads), xyz.data(), {}, connectivity, offsets, {id.data(), speed.data()}) && written;
        }
        const std::string ice = path + "/particles/fnpf_ice";
        check(written, "6 outputs of ice floes (cells)");
        check(slurp(ice + "/cell_sets/cell_end/zarr.json").find("\"shape\": [2]")!=std::string::npos,
              "two cell sets for six outputs");
        check(slurp(ice + "/cell_data/speed/zarr.json").find("\"shape\": [10]")!=std::string::npos,
              "a value per cell and output (2*4 + 1*2)");
        check(slurp(ice + "/zarr.json").find("\"cells\": \"polys\"")!=std::string::npos, "cells: polys");
        std::vector<int32_t> bad = {0, 1, 2};
        check(!floes.output(6, 6, 0, nullptr, {}, bad, {4}, {nullptr, nullptr}) && !floes.usable(),
              "inconsistent cells are refused");
    }

    std::cout<<(nfail ? "FAILED" : "all passed")<<std::endl;
    return nfail ? 1 : 0;
}
