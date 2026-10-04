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

// Writes LAGOON stores (LAGOON lagoon-format SPEC.md, version 0.1): one Zarr v3
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

#endif
