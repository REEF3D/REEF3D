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

#include"lagoon_store.h"

#include<algorithm>
#include<atomic>
#include<cctype>
#include<cerrno>
#include<cmath>
#include<cstdint>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<ctime>
#include<fstream>
#include<iostream>
#include<limits>
#include<sstream>
#include<stdexcept>
#include<sys/stat.h>
#include<sys/types.h>
#include<thread>
#include<zlib.h>

namespace
{
const uint64_t EMPTY = std::numeric_limits<uint64_t>::max();
const int CHUNK_VALUES = 2048;  // values an inner chunk should hold at least

void make_dirs(const std::string &dir)
{
    std::string part;
    std::stringstream ss(dir);
    std::string item;
    if(!dir.empty() && dir[0]=='/')
        part = "/";
    while(std::getline(ss,item,'/'))
    {
        if(item.empty())
            continue;
        part += item;
        if(mkdir(part.c_str(),0777)!=0 && errno!=EEXIST)
            throw std::runtime_error("lagoon_store: cannot create "+part);
        part += "/";
    }
}

// write a file so that a reader never sees half of it
void write_file(const std::string &file, const std::string &content)
{
    const std::string temporary = file + ".tmp";
    {
        std::ofstream out(temporary.c_str(), std::ios::binary);
        if(!out)
            throw std::runtime_error("lagoon_store: cannot write "+file);
        out.write(content.data(), content.size());
    }
    if(std::rename(temporary.c_str(), file.c_str())!=0)
        throw std::runtime_error("lagoon_store: cannot rename "+temporary);
}

int tile_size(int n, int most)
{
    if(n<=most)
        return std::max(n,1);
    const int tiles = (n + most - 1)/most;
    return (n + tiles - 1)/tiles;
}

uint32_t crc32c(const unsigned char *data, size_t n)
{
    static uint32_t table[256];
    static bool ready = false;
    if(!ready)
    {
        for(uint32_t i=0; i<256; ++i)
        {
            uint32_t c = i;
            for(int k=0; k<8; ++k)
                c = (c & 1) ? (c >> 1) ^ 0x82F63B78u : c >> 1;
            table[i] = c;
        }
        ready = true;
    }
    uint32_t crc = 0xFFFFFFFFu;
    for(size_t i=0; i<n; ++i)
        crc = table[(crc ^ data[i]) & 0xFF] ^ (crc >> 8);
    return crc ^ 0xFFFFFFFFu;
}

// a byte plane whose bytes are this spread out (Shannon entropy, bits per byte) is
// noise to deflate: the low bytes of the mantissas. Huffman coding would save less
// than 5 %, and searching it for matches takes most of the time of a chunk.
const double NOISE_BITS = 7.6;

double byte_entropy(const unsigned char *data, size_t n)
{
    size_t counts[256] = {0};
    for(size_t i=0; i<n; ++i)
        ++counts[data[i]];
    double bits = 0.0;
    for(int c=0; c<256; ++c)
        if(counts[c])
        {
            const double share = double(counts[c])/double(n);
            bits -= share*std::log2(share);
        }
    return bits;
}

// one z_stream per thread, reset for each chunk (deflateInit2 allocates and clears
// a few hundred kilobytes)
struct deflater
{
    z_stream zs;
    bool ready = false;
    ~deflater()
    {
        if(ready)
            deflateEnd(&zs);
    }
};

// gzip of a byte-shuffled chunk (count values of itemsize bytes, all first bytes,
// then all second bytes, ...): one gzip stream, as the gzip codec reads it, whose
// noise planes (NOISE_BITS) are stored rather than compressed (deflate level 0)
std::string gzip_shuffled(const unsigned char *shuffled, size_t count, int itemsize, int level)
{
    thread_local deflater d;
    z_stream &zs = d.zs;
    if(!d.ready)
    {
        std::memset(&zs, 0, sizeof(zs));
        if(deflateInit2(&zs, level, Z_DEFLATED, 15+16, 8, Z_DEFAULT_STRATEGY)!=Z_OK)  // 15+16: gzip
            throw std::runtime_error("lagoon_store: deflateInit2 failed");
        d.ready = true;
    }
    else if(deflateReset(&zs)!=Z_OK)
        throw std::runtime_error("lagoon_store: deflateReset failed");
    const size_t bytes = count*size_t(itemsize);
    std::string out(deflateBound(&zs, bytes) + 64*size_t(itemsize) + 64, '\0');
    size_t written = 0;
    int current = -1;
    for(int b=0; b<itemsize; ++b)
    {
        const unsigned char *plane = shuffled + size_t(b)*count;
        if(out.size() - written < count + 1024)
            out.resize(written + 2*count + 1024);
        zs.next_out = reinterpret_cast<unsigned char*>(&out[written]);
        zs.avail_out = uInt(out.size() - written);
        const int want = level>0 && count>=256 && byte_entropy(plane, count)>NOISE_BITS ? 0 : level;
        if(want!=current && deflateParams(&zs, want, Z_DEFAULT_STRATEGY)==Z_OK)
            current = want;
        written = out.size() - zs.avail_out;  // deflateParams may flush a block
        const int flush = b==itemsize-1 ? Z_FINISH : Z_NO_FLUSH;
        zs.next_in = const_cast<unsigned char*>(plane);
        zs.avail_in = uInt(count);
        int status;
        do
        {
            if(written==out.size())
                out.resize(2*out.size());
            zs.next_out = reinterpret_cast<unsigned char*>(&out[written]);
            zs.avail_out = uInt(out.size() - written);
            status = deflate(&zs, flush);
            written = out.size() - zs.avail_out;
            if(status==Z_STREAM_ERROR)
                throw std::runtime_error("lagoon_store: deflate failed");
        }
        while(flush==Z_FINISH ? status!=Z_STREAM_END : zs.avail_in>0 || zs.avail_out==0);
    }
    out.resize(written);
    return out;
}

void put_u64(std::string &s, uint64_t v)
{
    for(int b=0; b<8; ++b)
        s.push_back(char((v >> (8*b)) & 0xFF));
}

uint64_t get_u64(const unsigned char *p)
{
    uint64_t v = 0;
    for(int b=7; b>=0; --b)
        v = (v << 8) | p[b];
    return v;
}


std::string utc_now()
{
    char text[32];
    std::time_t now = std::time(nullptr);
    std::strftime(text, sizeof(text), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
    return text;
}

std::string codecs_json(int level, int itemsize=4)
{
    std::ostringstream j;
    j << "[{\"name\": \"bytes\", \"configuration\": {\"endian\": \"little\"}}, "
      << "{\"name\": \"numcodecs.shuffle\", \"configuration\": {\"elementsize\": " << itemsize << "}}, "
      << "{\"name\": \"gzip\", \"configuration\": {\"level\": " << level << "}}]";
    return j.str();
}

std::string ints_json(const std::vector<long long> &v)
{
    std::ostringstream j;
    j << "[";
    for(size_t i=0; i<v.size(); ++i)
        j << (i ? ", " : "") << v[i];
    j << "]";
    return j.str();
}

std::string names_json(const std::vector<std::string> &v)
{
    std::string j = "[";
    for(size_t i=0; i<v.size(); ++i)
        j += (i ? ", " : "") + lagoon_store::json_string(v[i]);
    return j + "]";
}

// byte shuffle (all first bytes, then all second bytes, ...) and gzip: one chunk
std::string shuffled_gzip(const void *data, size_t count, int itemsize, int level)
{
    const unsigned char *bytes = static_cast<const unsigned char*>(data);
    std::vector<unsigned char> shuffled(count*itemsize);
    for(int b=0; b<itemsize; ++b)
    {
        unsigned char *plane = &shuffled[b*count];
        for(size_t e=0; e<count; ++e)
            plane[e] = bytes[size_t(itemsize)*e + b];
    }
    return gzip_shuffled(shuffled.data(), count, itemsize, level);
}

// zarr.json of a regular-chunked, compressed array
void write_array_meta(const std::string &dir, const std::string &dtype, int itemsize,
                      const std::vector<long long> &shape, const std::vector<long long> &chunks,
                      const std::string &fill, const std::vector<std::string> &dims,
                      const std::string &attributes, int level)
{
    make_dirs(dir);
    std::ostringstream j;
    j << "{\"zarr_format\": 3, \"node_type\": \"array\", \"shape\": " << ints_json(shape)
      << ", \"data_type\": \"" << dtype << "\", "
      << "\"chunk_grid\": {\"name\": \"regular\", \"configuration\": {\"chunk_shape\": " << ints_json(chunks) << "}}, "
      << "\"chunk_key_encoding\": {\"name\": \"default\", \"configuration\": {\"separator\": \"/\"}}, "
      << "\"fill_value\": " << fill << ", "
      << "\"codecs\": " << codecs_json(level, itemsize) << ", "
      << "\"dimension_names\": " << names_json(dims) << ", "
      << "\"attributes\": {" << attributes << "}}";
    write_file(dir + "/zarr.json", j.str());
}

bool exists(const std::string &file)
{
    struct stat info;
    return stat(file.c_str(), &info)==0;
}

// a plain (unsharded, uncompressed) float64 or int64 array in one chunk
void write_plain_array(const std::string &dir, const std::vector<double> &values,
                       const std::string &dimension, const std::string &units)
{
    make_dirs(dir + "/c");
    std::ostringstream j;
    const size_t n = values.size();
    j << "{\"zarr_format\": 3, \"node_type\": \"array\", \"shape\": [" << n << "], "
      << "\"data_type\": \"float64\", "
      << "\"chunk_grid\": {\"name\": \"regular\", \"configuration\": {\"chunk_shape\": [" << std::max<size_t>(n,1) << "]}}, "
      << "\"chunk_key_encoding\": {\"name\": \"default\", \"configuration\": {\"separator\": \"/\"}}, "
      << "\"fill_value\": \"NaN\", "
      << "\"codecs\": [{\"name\": \"bytes\", \"configuration\": {\"endian\": \"little\"}}], "
      << "\"dimension_names\": [\"" << dimension << "\"], "
      << "\"attributes\": {" << (units.empty() ? "" : "\"units\": \"" + units + "\"") << "}}";
    write_file(dir + "/zarr.json", j.str());
    std::string bytes(reinterpret_cast<const char*>(values.data()), n*sizeof(double));
    write_file(dir + "/c/0", bytes);
}
}

lagoon_store::lagoon_store(const std::string &path_, int shard_time_, int gzip_level_)
    : path(path_), shard_time(std::max(shard_time_,1)), gzip_level(gzip_level_)
{
}

std::string lagoon_store::json_string(const std::string &text)
{
    std::string out = "\"";
    for(char c : text)
    {
        if(c=='"' || c=='\\')
            out += '\\', out += c;
        else if(c=='\n')
            out += "\\n";
        else if((unsigned char)c < 0x20)
        {
            char code[8];
            std::snprintf(code, sizeof(code), "\\u%04x", c);
            out += code;
        }
        else
            out += c;
    }
    return out + "\"";
}

void lagoon_store::create_root(const std::string &solver, const std::string &run_json)
{
    make_dirs(path);
    std::ostringstream j;
    j << "{\"zarr_format\": 3, \"node_type\": \"group\", \"attributes\": {"
      << "\"lagoon\": {\"format\": \"lagoon\", \"version\": \"0.3\", \"solver\": " << json_string(solver)
      << ", \"outputs\": [], \"created_by\": \"REEF3D\", \"created\": \"" << utc_now()
      << "\", \"source\": \"REEF3D\"}";
    if(!run_json.empty())
        j << ", \"reef3d_run\": " << run_json;
    j << "}}";
    write_file(path + "/zarr.json", j.str());
}

void lagoon_store::create_output(const std::string &output, const std::string &grid,
                                 const std::vector<double> &x, const std::vector<double> &y,
                                 const std::vector<double> &levels,
                                 const std::vector<variable> &variables,
                                 int blocks, const std::string &source)
{
    const std::string dir = path + "/" + output;
    make_dirs(dir + "/blocks");
    std::ostringstream j;
    j << "{\"zarr_format\": 3, \"node_type\": \"group\", \"attributes\": {\"lagoon\": {"
      << "\"role\": " << json_string(output) << ", \"grid\": " << json_string(grid)
      << ", \"shape\": {\"x\": " << x.size() << ", \"y\": " << y.size();
    if(grid=="sigma" || grid=="cartesian")
        j << ", \"level\": " << levels.size();
    j << "}, \"variables\": {";
    for(size_t v=0; v<variables.size(); ++v)
    {
        j << (v ? ", " : "") << json_string(variables[v].name) << ": {\"components\": "
          << variables[v].components;
        if(!variables[v].units.empty())
            j << ", \"units\": " << json_string(variables[v].units);
        j << "}";
    }
    j << "}, \"blocks\": " << blocks;
    if(!source.empty())
        j << ", \"source\": " << json_string(source);
    j << "}}}";
    write_file(dir + "/zarr.json", j.str());
    write_file(dir + "/blocks/zarr.json", "{\"zarr_format\": 3, \"node_type\": \"group\", \"attributes\": {}}");
    write_plain_array(dir + "/x", x, "x", "m");
    write_plain_array(dir + "/y", y, "y", "m");
    if(grid=="sigma")
        write_plain_array(dir + "/sigma", levels, "level", "");
    if(grid=="cartesian")
        write_plain_array(dir + "/z", levels, "level", "m");
    times[output].clear();
    steps[output].clear();
    commit(output, -1, 0.0, 0);  // empty time and step arrays

    // the root lists its outputs
    std::ifstream in((path + "/zarr.json").c_str());
    std::stringstream root;
    root << in.rdbuf();
    std::string text = root.str();
    const std::string key = "\"outputs\": [";
    const size_t at = text.find(key);
    if(at!=std::string::npos && text.find(json_string(output), at)==std::string::npos)
    {
        const size_t end = at + key.size();
        const bool empty = text[end]==']';
        text.insert(end, json_string(output) + (empty ? "" : ", "));
        write_file(path + "/zarr.json", text);
    }
}

void lagoon_store::create_block(const std::string &output, int block, int i0, int j0,
                                int nx, int ny, int nz, const std::vector<variable> &variables, int rank,
                                bool cartesian, int k0)
{
    char name[16];
    std::snprintf(name, sizeof(name), "b%04d", block);
    const std::string dir = path + "/" + output + "/blocks/" + name;
    make_dirs(dir);
    std::ostringstream j;
    j << "{\"zarr_format\": 3, \"node_type\": \"group\", \"attributes\": {\"lagoon\": {"
      << "\"start\": [" << i0 << ", " << j0 << (cartesian ? ", " + std::to_string(k0) : std::string())
      << "], \"size\": [" << nx << ", " << ny << (cartesian ? ", " + std::to_string(nz) : std::string()) << "], "
      << "\"rank\": " << rank << "}}}";
    write_file(dir + "/zarr.json", j.str());

    const int cy = tile_size(ny, 64);
    const int cx = tile_size(nx, std::max(64, 2*CHUNK_VALUES/cy));
    const int cz = cy*cx >= CHUNK_VALUES ? 1 : std::min(nz, (CHUNK_VALUES + cy*cx - 1)/(cy*cx));
    std::vector<variable> all;
    if(cartesian)
        ;  // the heights are the output's z
    else if(nz>1)
    {
        all.push_back({"z_bed",1,"m"});
        all.push_back({"z_surface",1,"m"});
        all.push_back({"z_offset",1,"m"});
    }
    else
        all.push_back({"z",1,"m"});
    all.insert(all.end(), variables.begin(), variables.end());
    for(const variable &v : all)
    {
        array_info a;
        a.dir = dir + "/" + v.name;
        a.levels = (cartesian || nz>1) && v.name!="z_bed" && v.name!="z_surface";
        a.nz = a.levels ? nz : 1;
        a.ny = ny;
        a.nx = nx;
        a.components = v.components;
        a.cz = a.levels ? std::max(cz,1) : 1;
        a.cy = cy;
        a.cx = cx;
        a.fill = v.name=="z_offset" ? 0.0f : std::numeric_limits<float>::quiet_NaN();
        arrays[output + "/" + std::to_string(block) + "/" + v.name] = a;
        if(v.name!="z_offset")  // made when first written (it is 0 until then)
            write_array_json(a, 0, v.units, v.name);
    }
}

void lagoon_store::write_array_json(const array_info &a, int nt, const std::string &units,
                                    const std::string &name) const
{
    make_dirs(a.dir);
    std::vector<long long> shape, shard, inner;
    std::vector<std::string> dims;
    shape.push_back(nt); shard.push_back(shard_time); inner.push_back(1); dims.push_back("time");
    if(a.levels)
    {
        shape.push_back(a.nz);
        shard.push_back(((a.nz + a.cz - 1)/a.cz)*a.cz);
        inner.push_back(a.cz);
        dims.push_back("level");
    }
    shape.push_back(a.ny); shard.push_back(((a.ny + a.cy - 1)/a.cy)*a.cy); inner.push_back(a.cy); dims.push_back("y");
    shape.push_back(a.nx); shard.push_back(((a.nx + a.cx - 1)/a.cx)*a.cx); inner.push_back(a.cx); dims.push_back("x");
    if(a.components>1)
    {
        shape.push_back(a.components); shard.push_back(a.components); inner.push_back(a.components);
        dims.push_back("component");
    }
    std::ostringstream j;
    j << "{\"zarr_format\": 3, \"node_type\": \"array\", \"shape\": " << ints_json(shape)
      << ", \"data_type\": \"float32\", "
      << "\"chunk_grid\": {\"name\": \"regular\", \"configuration\": {\"chunk_shape\": " << ints_json(shard) << "}}, "
      << "\"chunk_key_encoding\": {\"name\": \"default\", \"configuration\": {\"separator\": \"/\"}}, "
      << "\"fill_value\": " << (name=="z_offset" ? "0.0" : "\"NaN\"") << ", "
      << "\"codecs\": [{\"name\": \"sharding_indexed\", \"configuration\": {"
      << "\"chunk_shape\": " << ints_json(inner) << ", "
      << "\"codecs\": " << codecs_json(gzip_level) << ", "
      << "\"index_codecs\": [{\"name\": \"bytes\", \"configuration\": {\"endian\": \"little\"}}, {\"name\": \"crc32c\"}], "
      << "\"index_location\": \"end\"}}], "
      << "\"dimension_names\": " << names_json(dims) << ", "
      << "\"attributes\": {" << (units.empty() ? "" : "\"units\": " + json_string(units)) << "}}";
    write_file(a.dir + "/zarr.json", j.str());
}

void lagoon_store::write(const std::string &output, int block, int t, const std::string &name,
                         const float *data)
{
    const std::string key = output + "/" + std::to_string(block) + "/" + name;
    std::map<std::string, array_info>::const_iterator it = arrays.find(key);
    if(it==arrays.end())
        throw std::runtime_error("lagoon_store: unknown array " + key);
    const array_info &a = it->second;
    if(name=="z_offset")
    {
        // all zero (the usual case): nothing to store, the fill value is 0
        const size_t n = size_t(a.nz)*a.ny*a.nx;
        bool any = false;
        for(size_t i=0; i<n && !any; ++i)
            any = data[i]!=0.0f;
        struct stat info;
        const bool exists = stat((a.dir + "/zarr.json").c_str(), &info)==0;
        if(!any && !exists)
            return;
    }
    write_shard(a, t, data);
    write_array_json(a, t+1, name.rfind("z",0)==0 ? "m" : "", name);
}

void lagoon_store::write_shard(const array_info &a, int t, const float *data) const
{
    const int nzc = (a.nz + a.cz - 1)/a.cz;
    const int nyc = (a.ny + a.cy - 1)/a.cy;
    const int nxc = (a.nx + a.cx - 1)/a.cx;
    const size_t entries = size_t(shard_time)*nzc*nyc*nxc;  // the component dimension has one chunk
    const size_t index_bytes = 16*entries + 4;
    const int shard = t/shard_time;
    const int local = t%shard_time;

    std::string file = a.dir + "/c/" + std::to_string(shard) + (a.levels ? "/0" : "") + "/0/0" + (a.components>1 ? "/0" : "");
    make_dirs(file.substr(0, file.rfind('/')));

    // the index so far (a shard that is new, or whose last index is not valid, starts empty)
    std::vector<uint64_t> index(2*entries, EMPTY);
    long long data_end = 0;
    FILE *f = std::fopen(file.c_str(), "r+b");
    if(f)
    {
        std::fseek(f, 0, SEEK_END);
        const long long size = std::ftell(f);
        data_end = size;
        if(size >= (long long)index_bytes)
        {
            std::vector<unsigned char> raw(index_bytes);
            std::fseek(f, size - (long long)index_bytes, SEEK_SET);
            if(std::fread(raw.data(), 1, index_bytes, f)==index_bytes)
            {
                uint32_t stored = 0;
                for(int b=3; b>=0; --b)
                    stored = (stored << 8) | raw[16*entries + b];
                if(stored==crc32c(raw.data(), 16*entries))
                    for(size_t e=0; e<2*entries; ++e)
                        index[e] = get_u64(&raw[8*e]);
            }
        }
    }
    else
    {
        f = std::fopen(file.c_str(), "wb");
        if(!f)
            throw std::runtime_error("lagoon_store: cannot write " + file);
    }

    // the inner chunks of output t, compressed by `threads` threads, then appended
    // after everything there is
    const int comps = a.components;
    const size_t chunks = size_t(nzc)*nyc*nxc;
    std::vector<std::string> packed(chunks);
    uint32_t fill_bits;
    std::memcpy(&fill_bits, &a.fill, 4);
    const bool fill_nan = std::isnan(a.fill);
    std::atomic<size_t> next(0);
    auto encode = [&]()
    {
        const size_t n = size_t(a.cz)*a.cy*a.cx*comps;
        std::vector<unsigned char> shuffled(4*n);
        for(size_t c=next++; c<chunks; c=next++)
        {
            const int ic = int(c%nxc), jc = int((c/nxc)%nyc), kc = int(c/(size_t(nxc)*nyc));
            // gathered straight into the byte planes (little endian)
            unsigned char *p0 = &shuffled[0], *p1 = &shuffled[n], *p2 = &shuffled[2*n], *p3 = &shuffled[3*n];
            bool all_fill = true;
            size_t m = 0;
            for(int k=kc*a.cz; k<(kc+1)*a.cz; ++k)
            for(int j=jc*a.cy; j<(jc+1)*a.cy; ++j)
            {
                const bool row = k<a.nz && j<a.ny;
                const float *line = row ? data + (size_t(k)*a.ny + j)*a.nx*comps : nullptr;
                for(int i=ic*a.cx; i<(ic+1)*a.cx; ++i)
                for(int q=0; q<comps; ++q, ++m)
                {
                    uint32_t v = fill_bits;
                    if(row && i<a.nx)
                    {
                        std::memcpy(&v, line + size_t(i)*comps + q, 4);
                        if(v!=fill_bits && !(fill_nan && (v & 0x7fffffffu) > 0x7f800000u))
                            all_fill = false;
                    }
                    p0[m] = (unsigned char)(v);
                    p1[m] = (unsigned char)(v >> 8);
                    p2[m] = (unsigned char)(v >> 16);
                    p3[m] = (unsigned char)(v >> 24);
                }
            }
            if(!all_fill)
                packed[c] = gzip_shuffled(shuffled.data(), n, 4, gzip_level);
        }
    };
    const int workers = int(std::min<size_t>(size_t(std::max(threads, 1)), chunks)) - 1;
    std::vector<std::thread> pool;
    for(int w=0; w<workers; ++w)
        pool.emplace_back(encode);
    encode();
    for(std::thread &w : pool)
        w.join();

    std::fseek(f, data_end, SEEK_SET);
    long long offset = data_end;
    for(size_t c=0; c<chunks; ++c)
    {
        if(packed[c].empty())  // all fill: not stored
            continue;
        std::fwrite(packed[c].data(), 1, packed[c].size(), f);
        const size_t entry = size_t(local)*chunks + c;
        index[2*entry] = offset;
        index[2*entry+1] = packed[c].size();
        offset += packed[c].size();
    }
    // the new index after the chunks, with its checksum
    std::string tail;
    tail.reserve(index_bytes);
    for(size_t e=0; e<2*entries; ++e)
        put_u64(tail, index[e]);
    const uint32_t crc = crc32c(reinterpret_cast<const unsigned char*>(tail.data()), tail.size());
    for(int b=0; b<4; ++b)
        tail.push_back(char((crc >> (8*b)) & 0xFF));
    std::fwrite(tail.data(), 1, tail.size(), f);
    std::fclose(f);
}

void lagoon_store::commit(const std::string &output, int t, double time, long long step)
{
    std::vector<double> &tv = times[output];
    std::vector<long long> &sv = steps[output];
    if(t>=0)
    {
        tv.resize(t+1, std::numeric_limits<double>::quiet_NaN());
        sv.resize(t+1, 0);
        tv[t] = time;
        sv[t] = step;
    }
    const int chunk = 4096;
    const std::string dir = path + "/" + output;
    const size_t nt = tv.size();
    for(int which=0; which<2; ++which)
    {
        const std::string adir = dir + (which==0 ? "/step" : "/time");  // step first, time last
        make_dirs(adir + "/c");
        if(nt>0)
        {
            const size_t first = (t>=0 ? size_t(t) : 0)/chunk*chunk;
            std::string bytes;
            for(size_t i=first; i<first+chunk; ++i)
            {
                if(which==0)
                {
                    const long long v = i<nt ? sv[i] : 0;
                    put_u64(bytes, (uint64_t)v);
                }
                else
                {
                    const double v = i<nt ? tv[i] : std::numeric_limits<double>::quiet_NaN();
                    uint64_t bits;
                    std::memcpy(&bits, &v, 8);
                    put_u64(bytes, bits);
                }
            }
            write_file(adir + "/c/" + std::to_string(first/chunk), bytes);
        }
        std::ostringstream j;
        j << "{\"zarr_format\": 3, \"node_type\": \"array\", \"shape\": [" << nt << "], "
          << "\"data_type\": \"" << (which==0 ? "int64" : "float64") << "\", "
          << "\"chunk_grid\": {\"name\": \"regular\", \"configuration\": {\"chunk_shape\": [" << chunk << "]}}, "
          << "\"chunk_key_encoding\": {\"name\": \"default\", \"configuration\": {\"separator\": \"/\"}}, "
          << "\"fill_value\": " << (which==0 ? "0" : "\"NaN\"") << ", "
          << "\"codecs\": [{\"name\": \"bytes\", \"configuration\": {\"endian\": \"little\"}}], "
          << "\"dimension_names\": [\"time\"], "
          << "\"attributes\": {" << (which==1 ? "\"units\": \"s\"" : "") << "}}";
        write_file(adir + "/zarr.json", j.str());
    }
}

int lagoon_store::committed(const std::string &output) const
{
    std::map<std::string, std::vector<double> >::const_iterator it = times.find(output);
    return it==times.end() ? 0 : int(it->second.size());
}

void lagoon_store::level_offsets(const float *z, int nz, int n, const std::vector<double> &sigma,
                                 std::vector<float> &offsets)
{
    offsets.assign(size_t(nz)*n, 0.0f);
    for(int q=0; q<n; ++q)
    {
        const double bed = z[q];
        const double top = z[size_t(nz-1)*n + q];
        for(int k=0; k<nz; ++k)
        {
            const float zk = z[size_t(k)*n + q];
            const double offset = double(zk) - (bed + sigma[k]*(top - bed));
            const float a = std::fabs(zk);
            const double step = std::nextafter(a, std::numeric_limits<float>::infinity()) - a;
            if(std::fabs(offset) > 2.0*step)
                offsets[size_t(k)*n + q] = float(offset);
        }
    }
}

namespace
{
std::string attribute(const std::string &tag, const std::string &name)
{
    const std::string key = " " + name + "=\"";
    const size_t at = tag.find(key);
    if(at==std::string::npos)
        return "";
    const size_t end = tag.find('"', at + key.size());
    return tag.substr(at + key.size(), end - at - key.size());
}

std::string units_of(const std::string &name)
{
    if(name=="velocity") return "m/s";
    if(name=="pressure") return "Pa";
    if(name=="elevation" || name=="eta" || name=="Hs") return "m";
    if(name=="Fi") return "m2/s";
    return "";
}
}

bool lagoon_store::parse_vtu_header(const std::string &header, std::vector<vtu_array> &fields, long long &points)
{
    fields.clear();
    points = -1;
    const size_t begin = header.find("<PointData>");
    const size_t end = header.find("</PointData>");
    const size_t at_points = header.find("<Points>");
    if(begin==std::string::npos || end==std::string::npos || at_points==std::string::npos)
        return false;
    size_t at = begin;
    while((at = header.find("<DataArray", at))!=std::string::npos && at<end)
    {
        const std::string tag = header.substr(at, header.find('>', at) - at);
        const std::string name = attribute(tag, "Name");
        const std::string components = attribute(tag, "NumberOfComponents");
        const std::string offset = attribute(tag, "offset");
        if(attribute(tag, "type")!="Float32" || name.empty() || offset.empty())
            return false;
        vtu_array a;
        a.var.name = name;
        a.var.components = components.empty() ? 1 : std::atoi(components.c_str());
        a.var.units = units_of(name);
        a.offset = std::atoll(offset.c_str());
        fields.push_back(a);
        at += tag.size();
    }
    const size_t tag_at = header.find("<DataArray", at_points);
    if(tag_at==std::string::npos)
        return false;
    const std::string tag = header.substr(tag_at, header.find('>', tag_at) - tag_at);
    points = std::atoll(attribute(tag, "offset").c_str());
    return !fields.empty();
}


// ======================================================================== bodies
lagoon_bodies::lagoon_bodies(const std::string &path_, const std::string &solver_, const std::string &key_,
                             const std::string &source_, const std::string &run_json_, int gzip_level_)
    : path(path_), solver(solver_), key(key_), source(source_), run_json(run_json_),
      dir(path_ + "/bodies/" + key_), gzip_level(gzip_level_), ok(true), started(false), committed(0),
      store(path_)
{
}

void lagoon_bodies::quaternion(const double R[9], double q[4])
{
    const double m00=R[0], m01=R[1], m02=R[2], m10=R[3], m11=R[4], m12=R[5], m20=R[6], m21=R[7], m22=R[8];
    const double trace = m00 + m11 + m22;
    if(trace > 0.0)
    {
        const double s = 2.0*std::sqrt(trace + 1.0);
        q[0] = 0.25*s; q[1] = (m21 - m12)/s; q[2] = (m02 - m20)/s; q[3] = (m10 - m01)/s;
    }
    else if(m00 > m11 && m00 > m22)
    {
        const double s = 2.0*std::sqrt(1.0 + m00 - m11 - m22);
        q[0] = (m21 - m12)/s; q[1] = 0.25*s; q[2] = (m01 + m10)/s; q[3] = (m02 + m20)/s;
    }
    else if(m11 > m22)
    {
        const double s = 2.0*std::sqrt(1.0 + m11 - m00 - m22);
        q[0] = (m02 - m20)/s; q[1] = (m01 + m10)/s; q[2] = 0.25*s; q[3] = (m12 + m21)/s;
    }
    else
    {
        const double s = 2.0*std::sqrt(1.0 + m22 - m00 - m11);
        q[0] = (m10 - m01)/s; q[1] = (m02 + m20)/s; q[2] = (m12 + m21)/s; q[3] = 0.25*s;
    }
    const double norm = std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    const double sign = q[0] < 0.0 ? -1.0 : 1.0;
    for(int i=0; i<4; ++i)
        q[i] *= sign/norm;
}

void lagoon_bodies::start()
{
    make_dirs(path);
    if(!exists(path + "/zarr.json"))  // no P 18 volume output: the store is made here
        store.create_root(solver, run_json);
    make_dirs(path + "/bodies");
    if(!exists(path + "/bodies/zarr.json"))
        write_file(path + "/bodies/zarr.json", "{\"zarr_format\": 3, \"node_type\": \"group\", \"attributes\": {}}");
    make_dirs(dir);
    write_set_attributes();
    store.commit("bodies/" + key, -1, 0.0, 0);  // empty time and step arrays
    started = true;
}

void lagoon_bodies::write_set_attributes() const
{
    std::vector<long long> listed;
    for(const std::pair<const int, body> &b : bodies)
        if(b.second.listed)
            listed.push_back(b.first);
    std::ostringstream j;
    j << "{\"zarr_format\": 3, \"node_type\": \"group\", \"attributes\": {\"lagoon\": {"
      << "\"kind\": \"bodies\", \"bodies\": " << ints_json(listed)
      << ", \"dataset\": " << lagoon_store::json_string(key)
      << ", \"source\": " << lagoon_store::json_string(source) << "}}}";
    write_file(dir + "/zarr.json", j.str());
}

void lagoon_bodies::write_body_attributes(int number, const body &b) const
{
    std::ostringstream j;
    j.precision(17);
    j << "{\"zarr_format\": 3, \"node_type\": \"group\", \"attributes\": {\"lagoon\": {"
      << "\"rigid\": true, \"vertices\": " << b.points << ", \"fields\": {}, \"max_error\": "
      << b.max_error << "}}}";
    write_file(dir + "/body_" + std::to_string(number) + "/zarr.json", j.str());
}

void lagoon_bodies::write_motion(int number, const body &b, const char *name, int width) const
{
    const std::vector<double> &values = width==3 ? b.translation : b.rotation;
    const long long rows = b.rows;
    const long long chunk = 4096;
    const long long first = (rows-1)/chunk*chunk;  // the chunk of the newest row, written again
    std::vector<double> part(size_t(chunk)*width, std::numeric_limits<double>::quiet_NaN());
    for(long long r=first; r<rows; ++r)
        for(int c=0; c<width; ++c)
            part[size_t(r-first)*width + c] = values[size_t(r)*width + c];
    const std::string adir = dir + "/body_" + std::to_string(number) + "/" + name;
    make_dirs(adir + "/c/" + std::to_string(first/chunk));
    write_file(adir + "/c/" + std::to_string(first/chunk) + "/0", shuffled_gzip(part.data(), part.size(), 8, gzip_level));
    write_array_meta(adir, "float64", 8, {rows, width}, {chunk, width}, "\"NaN\"",
                     {"time", width==3 ? "xyz" : "wxyz"}, width==3 ? "\"units\": \"m\"" : "", gzip_level);
}

bool lagoon_bodies::output(int number, int points, const double *x0, const double *x, const double R[9],
                           const double c[3], double time, long long step)
{
    if(!ok)
        return false;
    try
    {
        if(!started)
            start();
        body &b = bodies[number];
        if(b.points==0)  // the mesh, once: vertices relative to the centre of gravity, triangles
        {
            b.points = points;
            const std::string bdir = dir + "/body_" + std::to_string(number);
            make_dirs(bdir);
            write_body_attributes(number, b);
            std::vector<float> vertices(size_t(points)*3);
            for(size_t i=0; i<vertices.size(); ++i)
                vertices[i] = float(x0[i]);
            make_dirs(bdir + "/vertices/c/0");
            write_file(bdir + "/vertices/c/0/0", shuffled_gzip(vertices.data(), vertices.size(), 4, gzip_level));
            write_array_meta(bdir + "/vertices", "float32", 4, {points, 3}, {std::max(points,1), 3}, "\"NaN\"",
                             {"vertex", "xyz"}, "\"units\": \"m\"", gzip_level);
            const int triangles = points/3;
            std::vector<int32_t> corners(size_t(triangles)*3);
            for(size_t i=0; i<corners.size(); ++i)
                corners[i] = int32_t(i);
            make_dirs(bdir + "/triangles/c/0");
            write_file(bdir + "/triangles/c/0/0", shuffled_gzip(corners.data(), corners.size(), 4, gzip_level));
            write_array_meta(bdir + "/triangles", "int32", 4, {triangles, 3}, {std::max(triangles,1), 3}, "0",
                             {"triangle", "corner"}, "", gzip_level);
        }
        if(points!=b.points)
        {
            std::cout<<"LAGOON: body "<<number<<" has another mesh now; no more bodies in the LAGOON store"<<std::endl;
            ok = false;
            return false;
        }

        // the motion as stored, and how far it puts each vertex from where it is
        double q[4];
        quaternion(R, q);
        const double w=q[0], qx=q[1], qy=q[2], qz=q[3];
        const double M[9] = {1-2*(qy*qy+qz*qz), 2*(qx*qy-w*qz),   2*(qx*qz+w*qy),
                             2*(qx*qy+w*qz),   1-2*(qx*qx+qz*qz), 2*(qy*qz-w*qx),
                             2*(qx*qz-w*qy),   2*(qy*qz+w*qx),   1-2*(qx*qx+qy*qy)};
        double largest = 1.0, error = 0.0;
        for(int i=0; i<points; ++i)
        {
            const double v0=float(x0[3*i]), v1=float(x0[3*i+1]), v2=float(x0[3*i+2]);
            for(int r=0; r<3; ++r)
            {
                const float moved = float(M[3*r]*v0 + M[3*r+1]*v1 + M[3*r+2]*v2 + c[r]);
                const float there = float(x[3*i+r]);
                error = std::max(error, double(std::fabs(moved - there)));
                largest = std::max(largest, double(std::fabs(there)));
            }
        }
        const float big = float(largest);
        const double tolerance = 4.0*double(std::nextafter(big, 2.0f*big) - big);
        if(!(error <= tolerance))
        {
            std::cout<<"LAGOON: body "<<number<<" is "<<error<<" m off its rigid motion at output "<<step
                     <<"; no more bodies in the LAGOON store (the VTP files are written)"<<std::endl;
            ok = false;
            return false;
        }

        b.translation.insert(b.translation.end(), c, c+3);
        b.rotation.insert(b.rotation.end(), q, q+4);
        ++b.rows;
        write_motion(number, b, "translation", 3);
        write_motion(number, b, "rotation", 4);
        if(error > b.max_error || b.rows==1)
        {
            b.max_error = std::max(b.max_error, error);
            write_body_attributes(number, b);
        }
        if(!b.listed)
        {
            b.listed = true;
            write_set_attributes();
        }
        pending.insert(std::make_pair(b.rows-1, std::make_pair(time, step)));

        // an output counts once every body has it: time last
        long long rows = b.rows;
        for(const std::pair<const int, body> &other : bodies)
            rows = std::min(rows, other.second.rows);
        while(committed < rows)
        {
            const std::pair<double, long long> &when = pending[committed];
            store.commit("bodies/" + key, int(committed), when.first, when.second);
            pending.erase(committed);
            ++committed;
        }
        return true;
    }
    catch(std::exception &problem)
    {
        std::cout<<"LAGOON: "<<problem.what()<<"; no more bodies in the LAGOON store"<<std::endl;
        ok = false;
        return false;
    }
}

// ======================================================================== particles
lagoon_particles::lagoon_particles(const std::string &path_, const std::string &solver_, const std::string &key_,
                                   const std::string &role_, const std::string &source_,
                                   const std::string &run_json_, const std::vector<field> &fields_, int gzip_level_,
                                   const std::string &cells_, const std::vector<field> &cell_fields_)
    : path(path_), solver(solver_), key(key_), role(role_), source(source_), run_json(run_json_),
      dir(path_ + "/particles/" + key_), cells(cells_), fields(fields_), cell_fields(cells_.empty() ? std::vector<field>() : cell_fields_),
      gzip_level(gzip_level_), ok(true), started(false), store(path_), rw(gzip_level_)
{
}

namespace
{
// a field's array: its name with characters other than [A-Za-z0-9_.-] as '_'
std::string array_name(const std::string &name)
{
    std::string out = name;
    for(char &c : out)
        if(!(std::isalnum((unsigned char)c) || c=='_' || c=='.' || c=='-'))
            c = '_';
    return out.empty() ? "_" : out;
}

void write_group(const std::string &dir)
{
    make_dirs(dir);
    write_file(dir + "/zarr.json", "{\"zarr_format\": 3, \"node_type\": \"group\", \"attributes\": {}}");
}
}

lagoon_particles::growing lagoon_particles::make(const std::string &where, const field &f) const
{
    return lagoon_rows::make(where, f.integer ? "int32" : "float32", f.components);
}

std::string lagoon_particles::fields_json(const std::vector<field> &list) const
{
    std::ostringstream j;
    j << "{";
    for(size_t f=0; f<list.size(); ++f)
        j << (f ? ", " : "") << lagoon_store::json_string(list[f].name) << ": {\"components\": " << list[f].components
          << ", \"data_type\": \"" << (list[f].integer ? "int32" : "float32") << "\", \"array\": "
          << lagoon_store::json_string(array_name(list[f].name)) << "}";
    j << "}";
    return j.str();
}

void lagoon_particles::start()
{
    make_dirs(path);
    if(!exists(path + "/zarr.json"))
        store.create_root(solver, run_json);
    make_dirs(path + "/particles");
    if(!exists(path + "/particles/zarr.json"))
        write_group(path + "/particles");
    write_group(dir + "/point_data");

    std::ostringstream j;
    j << "{\"zarr_format\": 3, \"node_type\": \"group\", \"attributes\": {\"lagoon\": {"
      << "\"kind\": \"particles\", \"cells\": " << (cells.empty() ? "null" : lagoon_store::json_string(cells))
      << ", \"point_fields\": " << fields_json(fields) << ", \"cell_fields\": " << fields_json(cell_fields)
      << ", \"dataset\": " << lagoon_store::json_string(key)
      << ", \"role\": " << lagoon_store::json_string(role)
      << ", \"source\": " << lagoon_store::json_string(source) << "}}}";
    write_file(dir + "/zarr.json", j.str());

    arrays.clear();
    growing position = make(dir + "/position", {"position", 3, false});
    position.metre = true;
    arrays.push_back(position);
    for(const field &f : fields)
        arrays.push_back(make(dir + "/point_data/" + array_name(f.name), f));
    for(const growing &a : arrays)
        write_meta(a);
    point_end.clear();
    write_rows("point_end", point_end, "time");
    cell_arrays.clear();
    if(!cells.empty())
    {
        write_group(dir + "/cell_sets");
        write_group(dir + "/cell_data");
        cell_arrays.push_back(make(dir + "/cell_sets/offsets", {"offsets", 1, true}));
        cell_arrays.push_back(make(dir + "/cell_sets/connectivity", {"connectivity", 1, true}));
        for(const field &f : cell_fields)
            cell_arrays.push_back(make(dir + "/cell_data/" + array_name(f.name), f));
        for(const growing &a : cell_arrays)
            write_meta(a);
        write_rows("cell_set", cell_set, "time");
        write_rows("cell_sets/cell_end", cell_end, "set");
        write_rows("cell_sets/connectivity_end", connectivity_end, "set");
        if(!cell_fields.empty())
            write_rows("cell_data_end", cell_data_end, "time");
    }
    store.commit("particles/" + key, -1, 0.0, 0);  // empty time and step
    started = true;
}

// ======================================================================== growing arrays
lagoon_rows::array lagoon_rows::make(const std::string &dir, const std::string &dtype, int components,
                                     bool metre, const std::string &dimension)
{
    array a;
    a.dir = dir;
    a.dtype = dtype;
    a.fill = (dtype=="float32" || dtype=="float64") ? "\"NaN\"" : "0";
    a.components = components;
    a.itemsize = (dtype=="float64" || dtype=="int64") ? 8 : 4;
    a.metre = metre;
    a.dimension = dimension;
    return a;
}

void lagoon_rows::write_meta(const array &a) const
{
    std::ostringstream j;
    std::vector<long long> shape = {a.rows}, shard = {ROWS*CHUNKS_PER_SHARD}, inner = {ROWS};
    std::vector<std::string> dims = {a.dimension};
    if(a.components>1)
    {
        shape.push_back(a.components); shard.push_back(a.components); inner.push_back(a.components);
        dims.push_back("component");
    }
    j << "{\"zarr_format\": 3, \"node_type\": \"array\", \"shape\": " << ints_json(shape)
      << ", \"data_type\": \"" << a.dtype << "\", "
      << "\"chunk_grid\": {\"name\": \"regular\", \"configuration\": {\"chunk_shape\": " << ints_json(shard) << "}}, "
      << "\"chunk_key_encoding\": {\"name\": \"default\", \"configuration\": {\"separator\": \"/\"}}, "
      << "\"fill_value\": " << a.fill << ", "
      << "\"codecs\": [{\"name\": \"sharding_indexed\", \"configuration\": {"
      << "\"chunk_shape\": " << ints_json(inner) << ", "
      << "\"codecs\": " << codecs_json(gzip_level, a.itemsize) << ", "
      << "\"index_codecs\": [{\"name\": \"bytes\", \"configuration\": {\"endian\": \"little\"}}, {\"name\": \"crc32c\"}], "
      << "\"index_location\": \"end\"}}], "
      << "\"dimension_names\": " << names_json(dims) << ", "
      << "\"attributes\": {" << (a.metre ? "\"units\": \"m\"" : "") << "}}";
    make_dirs(a.dir);
    write_file(a.dir + "/zarr.json", j.str());
}

void lagoon_rows::write_shard(const array &a, long long shard, const std::vector<std::string> &chunks) const
{
    std::string content;
    std::vector<uint64_t> index(2*CHUNKS_PER_SHARD, EMPTY);
    uint64_t offset = 0;
    for(size_t c=0; c<chunks.size() && c<size_t(CHUNKS_PER_SHARD); ++c)
    {
        if(chunks[c].empty())
            continue;
        index[2*c] = offset;
        index[2*c+1] = chunks[c].size();
        content += chunks[c];
        offset += chunks[c].size();
    }
    std::string tail;
    for(uint64_t v : index)
        put_u64(tail, v);
    const uint32_t crc = crc32c(reinterpret_cast<const unsigned char*>(tail.data()), tail.size());
    for(int b=0; b<4; ++b)
        tail.push_back(char((crc >> (8*b)) & 0xFF));
    content += tail;
    const std::string folder = a.dir + "/c" + (a.components>1 ? "/" + std::to_string(shard) : "");
    make_dirs(folder);
    write_file(folder + "/" + (a.components>1 ? "0" : std::to_string(shard)), content);
}

void lagoon_rows::append(array &a, const unsigned char *data, size_t rows) const
{
    const size_t row_bytes = size_t(a.itemsize)*a.components;
    const size_t chunk_bytes = size_t(ROWS)*row_bytes;
    const size_t values_per_chunk = size_t(ROWS)*a.components;
    size_t done = 0;
    while(done < rows)
    {
        const size_t room = (chunk_bytes - a.tail.size())/row_bytes;
        const size_t take = std::min(room, rows - done);
        a.tail.insert(a.tail.end(), data + done*row_bytes, data + (done+take)*row_bytes);
        done += take;
        if(a.tail.size()==chunk_bytes)  // a full inner chunk: final
        {
            a.encoded.push_back(shuffled_gzip(a.tail.data(), values_per_chunk, a.itemsize, gzip_level));
            a.tail.clear();
            if(int(a.encoded.size())==CHUNKS_PER_SHARD)  // a full shard: final
            {
                write_shard(a, a.shard, a.encoded);
                a.encoded.clear();
                ++a.shard;
            }
        }
    }
    if(rows==0)
        return;
    a.rows += (long long)rows;
    // the shard being filled, with the inner chunk being filled (padded with the fill value)
    if(!a.tail.empty() || !a.encoded.empty())
    {
        std::vector<std::string> chunks = a.encoded;
        if(!a.tail.empty())
        {
            std::vector<unsigned char> padded(chunk_bytes, 0);
            std::memcpy(padded.data(), a.tail.data(), a.tail.size());
            if(a.dtype=="float32")
            {
                const float nan = std::numeric_limits<float>::quiet_NaN();
                for(size_t b=a.tail.size(); b<chunk_bytes; b+=4)
                    std::memcpy(&padded[b], &nan, 4);
            }
            else if(a.dtype=="float64")
            {
                const double nan = std::numeric_limits<double>::quiet_NaN();
                for(size_t b=a.tail.size(); b<chunk_bytes; b+=8)
                    std::memcpy(&padded[b], &nan, 8);
            }
            chunks.push_back(shuffled_gzip(padded.data(), values_per_chunk, a.itemsize, gzip_level));
        }
        write_shard(a, a.shard, chunks);
    }
    write_meta(a);
}

void lagoon_rows::write_ends(const std::string &adir, const std::vector<long long> &rows, const char *dimension) const
{
    const long long chunk = 4096;
    const size_t nt = rows.size();
    make_dirs(adir + "/c");
    if(nt>0)
    {
        const size_t first = (nt-1)/chunk*chunk;
        std::string bytes;
        for(size_t i=first; i<first+chunk; ++i)
            put_u64(bytes, (uint64_t)(i<nt ? rows[i] : 0));
        write_file(adir + "/c/" + std::to_string(first/chunk), bytes);
    }
    std::ostringstream j;
    j << "{\"zarr_format\": 3, \"node_type\": \"array\", \"shape\": [" << nt << "], "
      << "\"data_type\": \"int64\", "
      << "\"chunk_grid\": {\"name\": \"regular\", \"configuration\": {\"chunk_shape\": [" << chunk << "]}}, "
      << "\"chunk_key_encoding\": {\"name\": \"default\", \"configuration\": {\"separator\": \"/\"}}, "
      << "\"fill_value\": 0, "
      << "\"codecs\": [{\"name\": \"bytes\", \"configuration\": {\"endian\": \"little\"}}], "
      << "\"dimension_names\": [\"" << dimension << "\"], \"attributes\": {}}";
    write_file(adir + "/zarr.json", j.str());
}

bool lagoon_particles::output(double time, long long step, size_t n, const float *xyz,
                              const std::vector<const void*> &values,
                              const std::vector<int32_t> &connectivity, const std::vector<int32_t> &offsets,
                              const std::vector<const void*> &cell_values)
{
    if(!ok)
        return false;
    try
    {
        if(values.size()!=fields.size() || cell_values.size()!=cell_fields.size())
            throw std::runtime_error("lagoon_particles: the arrays given are not those declared");
        if(cells.empty() && !offsets.empty())
            throw std::runtime_error("lagoon_particles: cells given for a set of points");
        if(!offsets.empty() && size_t(offsets.back())!=connectivity.size())
            throw std::runtime_error("lagoon_particles: the last offset is not the number of corners");
        if(!started)
            start();
        append(arrays[0], reinterpret_cast<const unsigned char*>(xyz), n);
        for(size_t f=0; f<fields.size(); ++f)
            append(arrays[f+1], static_cast<const unsigned char*>(values[f]), n);
        if(!cells.empty())
        {
            const bool same = !cell_end.empty() && connectivity==last_connectivity && offsets==last_offsets;
            if(!same)  // a new cell set: its offsets and corners, then where they end
            {
                append(cell_arrays[0], reinterpret_cast<const unsigned char*>(offsets.data()), offsets.size());
                append(cell_arrays[1], reinterpret_cast<const unsigned char*>(connectivity.data()), connectivity.size());
                cell_end.push_back((cell_end.empty() ? 0 : cell_end.back()) + (long long)offsets.size());
                connectivity_end.push_back((connectivity_end.empty() ? 0 : connectivity_end.back()) + (long long)connectivity.size());
                write_rows("cell_sets/cell_end", cell_end, "set");
                write_rows("cell_sets/connectivity_end", connectivity_end, "set");
                last_connectivity = connectivity;
                last_offsets = offsets;
            }
            for(size_t f=0; f<cell_fields.size(); ++f)
                append(cell_arrays[f+2], static_cast<const unsigned char*>(cell_values[f]), offsets.size());
            cell_set.push_back((long long)cell_end.size() - 1);
            write_rows("cell_set", cell_set, "time");
            if(!cell_fields.empty())
            {
                cell_data_end.push_back((cell_data_end.empty() ? 0 : cell_data_end.back()) + (long long)offsets.size());
                write_rows("cell_data_end", cell_data_end, "time");
            }
        }
        point_end.push_back((point_end.empty() ? 0 : point_end.back()) + (long long)n);
        write_rows("point_end", point_end, "time");
        store.commit("particles/" + key, int(point_end.size()-1), time, step);  // step, then time
        return true;
    }
    catch(std::exception &problem)
    {
        std::cout<<"LAGOON: "<<problem.what()<<"; no more "<<role<<" in the LAGOON store"<<std::endl;
        ok = false;
        return false;
    }
}

// ======================================================================== AMR surfaces
lagoon_amr::lagoon_amr(const std::string &path_, const std::string &solver_, const std::string &key_,
                       const std::string &source_, const std::string &run_json_,
                       const std::vector<field> &fields_, int gzip_level_)
    : path(path_), solver(solver_), key(key_), source(source_), run_json(run_json_),
      dir(path_ + "/amr/" + key_), fields(fields_), gzip_level(gzip_level_), ok(true), started(false),
      x_end(0), y_end(0), cell_end(0), store(path_), rw(gzip_level_)
{
}

void lagoon_amr::start()
{
    make_dirs(path);
    if(!exists(path + "/zarr.json"))
        store.create_root(solver, run_json);
    make_dirs(path + "/amr");
    if(!exists(path + "/amr/zarr.json"))
        write_group(path + "/amr");
    write_group(dir + "/grids");
    write_group(dir + "/cell_data");

    std::ostringstream j;
    j << "{\"zarr_format\": 3, \"node_type\": \"group\", \"attributes\": {\"lagoon\": {"
      << "\"kind\": \"amr\", \"fields\": {";
    for(size_t f=0; f<fields.size(); ++f)
        j << (f ? ", " : "") << lagoon_store::json_string(fields[f].name) << ": {\"components\": 1, \"data_type\": \""
          << (fields[f].integer ? "int32" : "float32") << "\", \"array\": "
          << lagoon_store::json_string(array_name(fields[f].name)) << "}";
    j << "}, \"dataset\": " << lagoon_store::json_string(key)
      << ", \"role\": \"amr\", \"refinement\": 2"
      << ", \"source\": " << lagoon_store::json_string(source) << "}}}";
    write_file(dir + "/zarr.json", j.str());

    table = {lagoon_rows::make(dir + "/grids/level", "int32", 1, false, "grid"),
             lagoon_rows::make(dir + "/grids/rank", "int32", 1, false, "grid"),
             lagoon_rows::make(dir + "/grids/size", "int32", 2, false, "grid"),
             lagoon_rows::make(dir + "/grids/x_end", "int64", 1, false, "grid"),
             lagoon_rows::make(dir + "/grids/y_end", "int64", 1, false, "grid"),
             lagoon_rows::make(dir + "/grids/cell_end", "int64", 1, false, "grid")};
    nodes = {lagoon_rows::make(dir + "/x", "float64", 1, true, "node"),
             lagoon_rows::make(dir + "/y", "float64", 1, true, "node")};
    cells.clear();
    for(const field &f : fields)
        cells.push_back(lagoon_rows::make(dir + "/cell_data/" + array_name(f.name),
                                          f.integer ? "int32" : "float32", 1, false, "cell"));
    for(const lagoon_rows::array &a : table)
        rw.write_meta(a);
    for(const lagoon_rows::array &a : nodes)
        rw.write_meta(a);
    for(const lagoon_rows::array &a : cells)
        rw.write_meta(a);
    grid_end.clear();
    rw.write_ends(dir + "/grid_end", grid_end, "time");
    store.commit("amr/" + key, -1, 0.0, 0);  // empty time and step
    started = true;
}

bool lagoon_amr::output(double time, long long step, const std::vector<grid> &grids)
{
    if(!ok)
        return false;
    try
    {
        const size_t nf = fields.size();
        std::vector<int32_t> level, rank, size;
        std::vector<int64_t> xe, ye, ce;
        std::vector<double> x, y;
        std::vector<std::vector<unsigned char> > values(nf);
        for(const grid &g : grids)
        {
            const size_t n = size_t(g.nx)*size_t(g.ny);
            if(g.nx<1 || g.ny<1 || g.x.size()!=size_t(g.nx+1) || g.y.size()!=size_t(g.ny+1)
               || g.values.size()!=nf*n)
                throw std::runtime_error("lagoon_amr: a grid's sizes do not fit");
            level.push_back(g.level);
            rank.push_back(g.rank);
            size.push_back(g.nx);
            size.push_back(g.ny);
            x.insert(x.end(), g.x.begin(), g.x.end());
            y.insert(y.end(), g.y.begin(), g.y.end());
            x_end += g.nx + 1;
            y_end += g.ny + 1;
            cell_end += (long long)n;
            xe.push_back(x_end);
            ye.push_back(y_end);
            ce.push_back(cell_end);
            for(size_t f=0; f<nf; ++f)
            {
                std::vector<unsigned char> &out = values[f];
                const size_t at = out.size();
                out.resize(at + 4*n);
                for(size_t c=0; c<n; ++c)
                {
                    const double v = g.values[f*n + c];
                    if(fields[f].integer)
                    {
                        const int32_t i = int32_t(std::lround(v));
                        std::memcpy(&out[at + 4*c], &i, 4);
                    }
                    else
                    {
                        const float r = float(v);
                        std::memcpy(&out[at + 4*c], &r, 4);
                    }
                }
            }
        }
        if(!started)
            start();
        const size_t ng = grids.size();
        // the rows first, the table, then where the output ends, time last
        for(size_t f=0; f<nf; ++f)
            rw.append(cells[f], values[f].data(), values[f].size()/4);
        rw.append(nodes[0], reinterpret_cast<const unsigned char*>(x.data()), x.size());
        rw.append(nodes[1], reinterpret_cast<const unsigned char*>(y.data()), y.size());
        rw.append(table[0], reinterpret_cast<const unsigned char*>(level.data()), ng);
        rw.append(table[1], reinterpret_cast<const unsigned char*>(rank.data()), ng);
        rw.append(table[2], reinterpret_cast<const unsigned char*>(size.data()), ng);
        rw.append(table[3], reinterpret_cast<const unsigned char*>(xe.data()), ng);
        rw.append(table[4], reinterpret_cast<const unsigned char*>(ye.data()), ng);
        rw.append(table[5], reinterpret_cast<const unsigned char*>(ce.data()), ng);
        grid_end.push_back((grid_end.empty() ? 0 : grid_end.back()) + (long long)ng);
        rw.write_ends(dir + "/grid_end", grid_end, "time");
        store.commit("amr/" + key, int(grid_end.size()-1), time, step);  // step, then time
        return true;
    }
    catch(std::exception &problem)
    {
        std::cout<<"LAGOON: "<<problem.what()<<"; no more AMR output in the LAGOON store"<<std::endl;
        ok = false;
        return false;
    }
}

void lagoon_amr::pack(const std::vector<grid> &grids, int nf, std::vector<double> &out)
{
    out.push_back(double(grids.size()));
    for(const grid &g : grids)
    {
        out.push_back(g.level);
        out.push_back(g.nx);
        out.push_back(g.ny);
        out.insert(out.end(), g.x.begin(), g.x.end());
        out.insert(out.end(), g.y.begin(), g.y.end());
        const size_t n = size_t(nf)*size_t(g.nx)*size_t(g.ny);
        if(g.values.size()!=n)
            throw std::runtime_error("lagoon_amr: a grid's values do not fit its size");
        out.insert(out.end(), g.values.begin(), g.values.end());
    }
}

bool lagoon_amr::unpack(const double *data, size_t n, int nf, int rank, std::vector<grid> &out)
{
    if(n<1)
        return false;
    size_t at = 0;
    const long long count = (long long)data[at++];
    for(long long q=0; q<count; ++q)
    {
        if(at + 3 > n)
            return false;
        grid g;
        g.level = int(data[at]);
        g.nx = int(data[at+1]);
        g.ny = int(data[at+2]);
        g.rank = rank;
        at += 3;
        const size_t need = size_t(g.nx+1) + size_t(g.ny+1) + size_t(nf)*size_t(g.nx)*size_t(g.ny);
        if(g.nx<1 || g.ny<1 || at + need > n)
            return false;
        g.x.assign(data + at, data + at + g.nx + 1);
        at += g.nx + 1;
        g.y.assign(data + at, data + at + g.ny + 1);
        at += g.ny + 1;
        const size_t nv = size_t(nf)*size_t(g.nx)*size_t(g.ny);
        g.values.assign(data + at, data + at + nv);
        at += nv;
        out.push_back(std::move(g));
    }
    return at==n;
}
