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

#ifndef LAGOON_STORE_H_
#define LAGOON_STORE_H_

// Writes LAGOON stores (LAGOON lagoon-format SPEC.md, version 0.2): one Zarr v3
// store per run, one group per output stream, one block per MPI rank.
//
// No MPI and no REEF3D types in here: every rank writes its own block's files
// (no shared files), rank 0 also writes the group metadata and, after all ranks
// have written an output (a barrier in the caller), its time (commit).
//
// Each array of a block is a sharded Zarr array: a shard file holds `shard_time`
// outputs of the whole block, cut into inner chunks of one output, a few levels
// and a tile of at most 64 x 64 points (byte shuffle + gzip). A shard grows by one
// output at a time: the new chunks are appended and a new index (with a CRC-32C
// checksum) is written after them, so a reader never sees a half-written index as
// valid; the space of the old index stays unused.

#include <map>
#include <string>
#include <utility>
#include <vector>

class lagoon_store
{
public:
    struct variable
    {
        std::string name;
        int components;     // 1: scalar, 3: vector (x, y, z)
        std::string units;  // may be empty
    };

    // path: the store directory, e.g. "./REEF3D_FNPF.lagoon"
    lagoon_store(const std::string &path, int shard_time=16, int gzip_level=1);

    // rank 0, once: the root group (run_json: the run log's run record, or "")
    void create_root(const std::string &solver, const std::string &run_json);

    // rank 0, once per output stream: the group, its grid lines and time arrays
    // grid: "sigma" (levels: their σ, 0 at the bed, 1 at the surface), "cartesian"
    // (levels: their heights z, CFD) or "surface" (levels empty)
    void create_output(const std::string &output, const std::string &grid,
                       const std::vector<double> &x, const std::vector<double> &y,
                       const std::vector<double> &levels,
                       const std::vector<variable> &variables,
                       int blocks, const std::string &source);

    // every rank, once per output stream: its block (i0, j0: index of its first
    // point in x and y; nx, ny: its number of points; nz: levels, 1 for surfaces).
    // A Cartesian block (cartesian=true) has no height arrays and may hold the
    // levels k0 ... k0+nz-1 only (CFD splits its grid in z too).
    void create_block(const std::string &output, int block, int i0, int j0,
                      int nx, int ny, int nz, const std::vector<variable> &variables, int rank,
                      bool cartesian=false, int k0=0);

    // every rank, every output t (0, 1, 2, ...): one array of its block.
    // data: nz*ny*nx*components floats, x fastest, then y, then the level, with
    // the components of a point together (as REEF3D writes its VTU points).
    // Grid arrays: "z_bed", "z_surface" (nz = 1) and "z_offset" for σ-grids, "z"
    // for surfaces.
    void write(const std::string &output, int block, int t, const std::string &name,
               const float *data);

    // rank 0, after every rank has written output t: its time and output number
    void commit(const std::string &output, int t, double time, long long step);

    // the outputs committed so far (0 for a new output group)
    int committed(const std::string &output) const;

    // level offsets of a column set (SPEC section 3): z - (z_bed + sigma (z_surface - z_bed)),
    // 0 within two float32 steps of z. z: nz*n values, levels outermost.
    static void level_offsets(const float *z, int nz, int n, const std::vector<double> &sigma,
                              std::vector<float> &offsets);

    static std::string json_string(const std::string &text);

    // the Float32 point arrays and the points of a VTK XML header (appended data):
    // name, components and offset (from the '_' that starts the appended data)
    struct vtu_array
    {
        variable var;
        long long offset;
    };
    static bool parse_vtu_header(const std::string &header, std::vector<vtu_array> &fields, long long &points);

private:
    struct array_info
    {
        std::string dir;
        int nz, ny, nx, components;   // shape after the time dimension
        int cz, cy, cx;               // inner chunk
        bool levels;                  // has the level dimension
        float fill;
    };

    std::string path;
    int shard_time;
    int gzip_level;
    std::map<std::string, array_info> arrays;   // "output/block/name"
    std::map<std::string, std::vector<double> > times;
    std::map<std::string, std::vector<long long> > steps;
    std::map<std::string, std::vector<variable> > block_variables;  // "output/block"

    void write_array_json(const array_info &a, int nt, const std::string &units,
                          const std::string &name) const;
    void write_shard(const array_info &a, int t, const float *data) const;
};

// Floating bodies of a store (LAGOON SPEC.md section 9), written by rank 0: the mesh of
// each 6DOF body once (its vertices relative to the centre of gravity at the start,
// x0) and per output the rigid motion that moves it there, x = R x0 + c: the
// translation c and the rotation R as a unit quaternion (w, x, y, z), w >= 0. The
// set is bodies/<key> (e.g. bodies/nhflow_body); its time and step arrays count an
// output once every body has it.
class lagoon_bodies
{
public:
    lagoon_bodies(const std::string &path, const std::string &solver, const std::string &key,
                  const std::string &source, const std::string &run_json, int gzip_level=1);

    // one output of one body: points vertices (3 per triangle, as REEF3D keeps them),
    // x0 and x as xyz of each point, R (row major) and c as REEF3D moved them.
    // False when the body is not in the store (the vertices are not where the motion
    // puts them, or the store cannot be written): its VTP file is needed then; from
    // then on no body is written to the store.
    bool output(int body, int points, const double *x0, const double *x, const double R[9],
                const double c[3], double time, long long step);

    bool usable() const { return ok; }

    // the unit quaternion (w, x, y, z), w >= 0, of a rotation matrix (row major)
    static void quaternion(const double R[9], double q[4]);

private:
    struct body
    {
        int points = 0;
        long long rows = 0;
        std::vector<double> translation, rotation;  // 3 and 4 per output
        bool listed = false;
        double max_error = 0.0;
    };
    std::string path, solver, key, source, run_json, dir;
    int gzip_level;
    bool ok, started;
    std::map<int, body> bodies;
    std::map<long long, std::pair<double, long long> > pending;  // row: time, step
    long long committed;
    lagoon_store store;

    void start();
    void write_set_attributes() const;
    void write_body_attributes(int number, const body &b) const;
    void write_motion(int number, const body &b, const char *name, int width) const;
};

// Particles of a store (LAGOON SPEC.md section 10), written by rank 0 once the
// particles of all ranks are gathered: the points of every output one after another
// in growing arrays (shards of 16 inner chunks of 65536 rows), the end of each
// output's points in point_end. Each output is written at once (the newest inner
// chunk is written again until it is full), so the store can be read during a run.
class lagoon_particles
{
public:
    struct field
    {
        std::string name;
        int components;  // 1 or 3
        bool integer;    // int32, else float32
    };

    // key: the particle set, e.g. "nhflow_particles"; role: how LAGOON shows it
    lagoon_particles(const std::string &path, const std::string &solver, const std::string &key,
                     const std::string &role, const std::string &source, const std::string &run_json,
                     const std::vector<field> &fields, int gzip_level=1);

    // one output: n points (xyz), and for each field n*components values, float or
    // int32 as declared, in the order of the fields. False when it could not be
    // written (no more output to the store then: the VTP files are needed).
    bool output(double time, long long step, size_t n, const float *xyz,
                const std::vector<const void*> &values);

    bool usable() const { return ok; }

    static constexpr long long ROWS = 65536;
    static constexpr int CHUNKS_PER_SHARD = 16;

private:
    struct growing
    {
        std::string dir, dtype, fill;
        int components = 1, itemsize = 4;
        long long rows = 0;
        long long shard = 0;                  // the shard being filled
        std::vector<std::string> encoded;     // its full inner chunks
        std::vector<unsigned char> tail;      // rows of the inner chunk being filled
    };
    std::string path, solver, key, role, source, run_json, dir;
    std::vector<field> fields;
    int gzip_level;
    bool ok, started;
    std::vector<growing> arrays;  // position, then the fields
    std::vector<long long> point_end;
    lagoon_store store;

    void start();
    void append(growing &a, const unsigned char *data, size_t rows);
    void write_shard(const growing &a, long long shard, const std::vector<std::string> &chunks) const;
    void write_meta(const growing &a) const;
    void write_point_end() const;
};

#endif
