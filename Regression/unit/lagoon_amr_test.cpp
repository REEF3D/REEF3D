// Standalone verification of the LAGOON AMR writer (lagoon_amr): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../src lagoon_amr_test.cpp ../../src/lagoon_store.cpp -lz -o lagoon_amr_test
// Run:    ./lagoon_amr_test     (writes lagoon_amr_test.lagoon in the current folder)
// The files follow LAGOON's lagoon-format SPEC.md 0.3, section 11; LAGOON reads them
// with lagoon_format.store (Store(...).amr["nhflow_amr"]).
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

// inner chunk c of a shard file: gunzip and unshuffle
static std::vector<unsigned char> inner(const std::string &file, int c, size_t values, int itemsize)
{
    const std::string s = slurp(file);
    const size_t entries = lagoon_rows::CHUNKS_PER_SHARD;
    const unsigned char *index = reinterpret_cast<const unsigned char*>(s.data()) + s.size() - 16*entries - 4;
    const uint64_t offset = u64(index + 16*c), nbytes = u64(index + 16*c + 8);
    std::vector<unsigned char> raw(values*itemsize), out(values*itemsize);
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
        for(int b=0; b<itemsize; ++b)
            out[itemsize*e+b] = raw[b*values + e];
    return out;
}

// a grid of nx x ny cells of size d from (x0, y0); values: eta = level + cell number / 1000, wet = cell parity
static lagoon_amr::grid make_grid(int level, int rank, double x0, double y0, int nx, int ny, double d)
{
    lagoon_amr::grid g;
    g.level = level;
    g.rank = rank;
    g.nx = nx;
    g.ny = ny;
    for(int i=0; i<=nx; ++i) g.x.push_back(x0 + d*i);
    for(int j=0; j<=ny; ++j) g.y.push_back(y0 + d*j);
    for(int c=0; c<nx*ny; ++c) g.values.push_back(level + 0.001*c);
    for(int c=0; c<nx*ny; ++c) g.values.push_back(double(c%2));
    return g;
}

// the grids of output t: level 0 of two ranks, a patch that moves, a second one from t = 1
static std::vector<lagoon_amr::grid> grids(int t)
{
    std::vector<lagoon_amr::grid> out = {make_grid(0, 0, 0.0, 0.0, 8, 6, 0.5), make_grid(0, 1, 4.0, 0.0, 8, 6, 0.5),
                                         make_grid(1, 0, 1.0 + 0.5*t, 0.5, 8, 8, 0.25)};
    if(t>=1)
        out.push_back(make_grid(2, 1, 6.0, 0.0, 12, 12, 0.125));
    return out;
}

int main()
{
    std::cout<<"LAGOON AMR writer"<<std::endl;
    const std::string path = "./lagoon_amr_test.lagoon";
    const std::vector<lagoon_amr::field> fields = {{"eta", false}, {"wetdry", true}};
    lagoon_amr writer(path, "NHFLOW", "nhflow_amr", "REEF3D_NHFLOW_AMR", "{\"type\": \"run\", \"run\": \"test\"}", fields);
    bool all = true;
    for(int t=0; t<4; ++t)
        all = writer.output(0.5*t, t, grids(t)) && all;
    check(all && writer.usable(), "4 outputs (3 and 4 grids)");

    const std::string set = path + "/amr/nhflow_amr";
    const std::string attrs = slurp(set + "/zarr.json");
    check(attrs.find("\"kind\": \"amr\"")!=std::string::npos && attrs.find("\"refinement\": 2")!=std::string::npos,
          "kind and refinement");
    check(attrs.find("\"wetdry\": {\"components\": 1, \"data_type\": \"int32\", \"array\": \"wetdry\"}")!=std::string::npos,
          "an integer field");
    check(slurp(path + "/zarr.json").find("\"version\": \"0.3\"")!=std::string::npos, "format version 0.3");
    check(slurp(path + "/amr/zarr.json").find("\"node_type\": \"group\"")!=std::string::npos, "the amr group");
    check(slurp(set + "/time/zarr.json").find("\"shape\": [4]")!=std::string::npos, "time counts 4 outputs");
    check(slurp(set + "/grid_end/zarr.json").find("\"shape\": [4]")!=std::string::npos, "grid_end per output");
    const std::string ends = slurp(set + "/grid_end/c/0");
    const unsigned char *e = reinterpret_cast<const unsigned char*>(ends.data());
    check(u64(e)==3 && u64(e+8)==7 && u64(e+16)==11 && u64(e+24)==15, "grids end at 3, 7, 11, 15");
    check(slurp(set + "/grids/size/zarr.json").find("\"shape\": [15, 2]")!=std::string::npos, "15 grids");
    check(slurp(set + "/x/zarr.json").find("\"data_type\": \"float64\"")!=std::string::npos, "coordinates as float64");
    // cells: 3 outputs of 48+48+64+144 and one of 48+48+64
    const long long cells = 3*(48+48+64+144) + (48+48+64);
    check(slurp(set + "/cell_data/eta/zarr.json").find("\"shape\": [" + std::to_string(cells) + "]")!=std::string::npos,
          "every cell");

    // read back: the level of every grid, a patch's x, a value
    std::vector<unsigned char> level = inner(set + "/grids/level/c/0", 0, lagoon_rows::ROWS, 4);
    const int32_t *lv = reinterpret_cast<const int32_t*>(level.data());
    check(lv[0]==0 && lv[1]==0 && lv[2]==1 && lv[3]==0 && lv[6]==2 && lv[14]==2, "levels read back");
    std::vector<unsigned char> xs = inner(set + "/x/c/0", 0, lagoon_rows::ROWS, 8);
    const double *x = reinterpret_cast<const double*>(xs.data());
    check(x[9]==4.0 && x[18]==1.0 && x[19]==1.25, "node coordinates read back");
    std::vector<unsigned char> eta = inner(set + "/cell_data/eta/c/0", 0, lagoon_rows::ROWS, 4);
    const float *v = reinterpret_cast<const float*>(eta.data());
    check(v[96 + 5]==float(1 + 0.005), "cell values read back");
    std::vector<unsigned char> wet = inner(set + "/cell_data/wetdry/c/0", 0, lagoon_rows::ROWS, 4);
    check(reinterpret_cast<const int32_t*>(wet.data())[7]==1, "integer cell values read back");

    // packed for MPI and back: the grids of a rank, its rank number set
    {
        std::vector<double> buffer;
        lagoon_amr::pack(grids(2), 2, buffer);
        std::vector<lagoon_amr::grid> back;
        const bool unpacked = lagoon_amr::unpack(buffer.data(), buffer.size(), 2, 5, back);
        bool same = unpacked && back.size()==4;
        for(size_t q=0; same && q<back.size(); ++q)
        {
            const lagoon_amr::grid &a = back[q], b = grids(2)[q];
            same = a.level==b.level && a.nx==b.nx && a.ny==b.ny && a.x==b.x && a.y==b.y && a.values==b.values && a.rank==5;
        }
        check(same, "pack and unpack");
        check(!lagoon_amr::unpack(buffer.data(), buffer.size() - 1, 2, 0, back), "a cut buffer is refused");
    }

    // a grid that does not fit: refused, no more output to the store
    {
        lagoon_amr bad("./lagoon_amr_bad.lagoon", "FNPF", "fnpf_amr", "REEF3D_FNPF_AMR", "", fields);
        std::vector<lagoon_amr::grid> g = {make_grid(0, 0, 0.0, 0.0, 4, 4, 1.0)};
        g[0].values.pop_back();
        check(!bad.output(0.0, 0, g) && !bad.usable(), "a grid with too few values is refused");
        check(!bad.output(0.0, 0, grids(0)), "and the set stops");
    }

    std::cout<<(nfail==0 ? "all passed" : "FAILED")<<std::endl;
    return nfail==0 ? 0 : 1;
}
