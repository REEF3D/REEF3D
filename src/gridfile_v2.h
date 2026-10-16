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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef GRIDFILE_V2_H_
#define GRIDFILE_V2_H_

// DIVEMesh -> REEF3D grid format, version 2 (reader; the writer is in DIVEMesh)
//
// file     = magic[8] | int32 version | int32 endian marker | section* | "END "
// section  = char tag[4] | int64 nbytes | payload[nbytes]
//
// DIVEMesh_Grid/grid-geometry.dat (one file, all ranks)
//   GHDR  int32 nint, int32[nint], int32 ndbl, double[ndbl]
//   OBJS  int32 nobj, per object: int32 role, keyword, index, raymode, invert,
//                                 int64 tri offset, int64 tri count, int32 npar, double[npar]
//   VERT  int64 nvert, double[3*nvert]  unique vertices x y z
//   TIDX  table(ntri x 3)                vertex indices of the triangles
//   CURV  int32 nint=3, int32 layout version, ni, nj, int32 ndbl, double[ndbl] (R 2, R 3 wl,
//         R 3 margin, R 5, raster spacing), double x[n], y[n], zbed[n] with n = (ni+1)(nj+1),
//         node (i,j) at i + (ni+1)*j: river corridor grid (R 1 1), i along the river from the
//         inflow, j across from the right bank; optional (geo_mesh::curv_*)
//
// DIVEMesh_Grid/grid-%06i.dat (one file per rank)
//   HEAD  int32 nint, int32[nint], int32 ndbl, double[ndbl]
//   FLAG  run-length coded int32 cell flags, i-j-k order: int64 nrun, (int32 value, int32 length)[nrun]
//   NODE  double XN[knox+1+2*marge], YN[...], ZN[...]
//   SURF  table(nsurf x 5)               i j k side group
//   PARA  6 x table(n x 3)               i j k
//   PACO  6 x table(n x 3)               i j k
//   SLFL  run-length coded int32 column flags, i-j order
//   SLPA  4 x table(n x 2)               i j
//   SLPC  4 x table(n x 2)               i j
//
// table(n x c): int64 n, int32 c, then column by column the differences to the
//               previous row, zigzag and varint coded (1 byte for |difference| < 64)
//   GEOB  double[knox*knoy]              geodat bed level (if G 10 > 0)
//   DATA  double[knox*knoy]              interpolated data (if D 10 > 0)

#include<vector>
#include<string>
#include<cstring>
#include<cstdint>
#include<fstream>
#include<iostream>
#include<cstdlib>

namespace gridv2
{
    const char magic_grid[8] = {'D','M','G','R','I','D','0','2'};
    const char magic_geom[8] = {'D','M','G','E','O','M','0','2'};
    const int version = 2;
    const int endian = 0x01020304;

    struct section
    {
        const char *data = nullptr;
        size_t size = 0;
        size_t pos = 0;
        bool good = true;

        bool get(void *val, size_t bytes)
        {
            if(pos+bytes>size)
            {
                good = false;
                std::memset(val,0,bytes);
                return false;
            }

            std::memcpy(val,data+pos,bytes);
            pos += bytes;
            return true;
        }

        int get_int()          { int v=0; get(&v,sizeof(int)); return v; }
        long long get_i64()    { int64_t v=0; get(&v,sizeof(int64_t)); return (long long)v; }
        double get_double()    { double v=0.0; get(&v,sizeof(double)); return v; }

        void get_list(std::vector<int> &iv, std::vector<double> &dv)
        {
            const int ni = get_int();
            iv.resize(ni>0 ? ni : 0);
            for(int &v : iv)
            v = get_int();

            const int nd = get_int();
            dv.resize(nd>0 ? nd : 0);
            for(double &v : dv)
            v = get_double();
        }

        uint64_t get_varint()
        {
            uint64_t v=0;
            int shift=0;

            while(pos<size && shift<64)
            {
                const unsigned char b = (unsigned char)data[pos++];
                v |= uint64_t(b & 127u) << shift;

                if((b & 128u)==0)
                return v;

                shift += 7;
            }

            good = false;
            return 0;
        }

        // integer table of nrow x ncol, row-major result; false if the shape differs
        bool get_table(std::vector<int> &out, long long nrow, int ncol)
        {
            const long long n = get_i64();
            const int c = get_int();

            if(n!=nrow || c!=ncol)
            {
                good = false;
                return false;
            }

            out.assign(size_t(nrow)*size_t(ncol),0);

            for(int cc=0; cc<ncol; ++cc)
            {
                long long prev=0;

                for(long long r=0; r<nrow; ++r)
                {
                    const uint64_t z = get_varint();
                    const long long d = (long long)(z >> 1) ^ -(long long)(z & 1u);
                    prev += d;
                    out[size_t(r)*size_t(ncol)+cc] = int(prev);
                }
            }

            return good;
        }

        // run-length coded ints, expanded to num values
        bool get_rle(std::vector<int> &out, size_t num)
        {
            out.resize(num);

            const long long nrun = get_i64();
            size_t q=0;

            for(long long r=0; r<nrun; ++r)
            {
                const int val = get_int();
                const int len = get_int();

                if(len<0 || q+size_t(len)>num)
                {
                    good = false;
                    return false;
                }

                for(int l=0; l<len; ++l)
                out[q++] = val;
            }

            if(q!=num)
            good = false;

            return good;
        }

        bool done() const { return pos==size; }
    };

    struct reader
    {
        std::vector<char> buf;
        std::vector<std::string> tag;
        std::vector<size_t> off, len;
        std::string error;
        std::string name;
        bool good = true;

        bool open(const char *fname, const char magic[8])
        {
            name = fname;

            std::ifstream in(fname, std::ios_base::binary | std::ios_base::ate);

            if(!in)
            {
                error = "Could not open DIVEMesh grid file: " + name;
                return false;
            }

            const std::streamsize n = in.tellg();
            in.seekg(0);

            buf.resize(size_t(n>0 ? n : 0));

            if(n>0 && !in.read(buf.data(),n))
            {
                error = "Could not read DIVEMesh grid file: " + name;
                return false;
            }

            if(buf.size()<16 || std::memcmp(buf.data(),magic,8)!=0)
            {
                error = "DIVEMesh grid file " + name + " is not in grid format v2";
                return false;
            }

            int ver=0, end=0;
            std::memcpy(&ver,buf.data()+8,sizeof(int));
            std::memcpy(&end,buf.data()+12,sizeof(int));

            if(end!=endian)
            {
                error = "DIVEMesh grid file " + name + " has a different byte order";
                return false;
            }

            if(ver!=version)
            {
                error = "DIVEMesh grid file " + name + ": unsupported grid format version " + std::to_string(ver);
                return false;
            }

            // section index
            size_t pos = 16;
            bool closed = false;

            while(pos+12<=buf.size())
            {
                std::string t(buf.data()+pos,4);
                int64_t nb=0;
                std::memcpy(&nb,buf.data()+pos+4,sizeof(int64_t));
                pos += 12;

                if(nb<0 || pos+size_t(nb)>buf.size())
                {
                    error = "DIVEMesh grid file " + name + " is truncated";
                    return false;
                }

                if(t=="END ")
                {
                    closed = true;
                    break;
                }

                tag.push_back(t);
                off.push_back(pos);
                len.push_back(size_t(nb));

                pos += size_t(nb);
            }

            if(!closed)
            {
                error = "DIVEMesh grid file " + name + " is truncated";
                return false;
            }

            return true;
        }

        bool find(const char *t, section &sc)
        {
            for(size_t q=0; q<tag.size(); ++q)
            if(tag[q]==std::string(t,4))
            {
                sc.data = buf.data()+off[q];
                sc.size = len[q];
                sc.pos = 0;
                sc.good = true;
                return true;
            }

            return false;
        }

        void check(const section &sc)
        {
            if(!sc.good)
            good = false;
        }

        bool ok() const { return good; }

        template<class L>
        void fail(L *p, const char *msg)
        {
            std::cout<<std::endl<<"!!! rank "<<p->mpirank<<": "<<name<<": "<<msg<<" !!!"<<std::endl<<std::endl;
            std::exit(1);
        }
    };
}

#endif
