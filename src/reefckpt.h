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

#ifndef REEFCKPT_H_
#define REEFCKPT_H_

//  reefckpt: checkpoints that do not depend on the decomposition (P 48, I 45).
//
//  REEF3D_Checkpoint/<count>/main.dat     rank 0: the boxes of the pieces, named scalars, named blobs
//  REEF3D_Checkpoint/<count>/piece-%06i   every rank: its arrays with their ghost cells, double/int
//  REEF3D_Checkpoint/latest               the folder of the last complete checkpoint (written last)
//
//  Arrays are given in the layouts of REEF3D (margin 3 on every side):
//    S  slice: (i-imin)*jmax + (j-jmin), i in [-3,knox+3), j in [-3,knoy+3)
//    F  node field FIJK: k in [-3,knoz+4) (kmaxF)
//    C  cell field IJK:  k in [-3,knoz+3) (kmax)
//  Reading on the same decomposition takes the own piece as it is (ghost cells included, bitwise).
//  On another decomposition a cell (I,J) takes its value from the piece whose interior holds
//  (I,J) clamped to the domain: the interior from its owner, the ghost cells of a boundary from the
//  piece that set them, the parallel ghost cells from the neighbour's interior.  The decomposition
//  must not cut the vertical direction (the sigma grids of FNPF, NHFLOW).

#include<string>
#include<vector>
#include<map>
#include<fstream>

class lexer;
class ghostcell;

// an output file of the run: on a restart from a checkpoint (I 45 1) an existing file is continued
// (appended); true: a new file, its header is to be written
bool reefckpt_open(std::ofstream &o, const char *name, int I45);

class reefckpt
{
public:
    // ---- writing (all ranks)
    void write_open(lexer*, int count);
    void scalar(const std::string&, double);                     // rank 0's value
    void blob(const std::string&, const std::vector<double>&);   // rank 0's data
    void rank_blob(const std::string&, const std::vector<double>&);   // every rank its own (in its piece)
    void put_S(const std::string&, const double*);
    void put_Si(const std::string&, const int*);
    void put_F(const std::string&, const double*);
    void put_Fi(const std::string&, const int*);
    void put_C(const std::string&, const double*);
    void put_Ci(const std::string&, const int*);
    void write_close(lexer*);

    // ---- reading (all ranks); dir empty: the latest checkpoint
    bool read_open(lexer*, const std::string &dir, std::string &msg);
    bool has(const std::string&) const;
    double get_scalar(const std::string&, double def=0.0) const;
    std::vector<double> get_blob(const std::string&) const;
    // the rank blobs of all pieces (in rank order of the checkpoint)
    std::vector<std::vector<double>> get_rank_blobs(const std::string&);
    bool get_S(const std::string&, double*);
    bool get_Si(const std::string&, int*);
    bool get_F(const std::string&, double*);
    bool get_Fi(const std::string&, int*);
    bool get_C(const std::string&, double*);
    bool get_Ci(const std::string&, int*);
    bool same_decomposition() const { return same; }
    int pieces() const { return (int)box.size(); }
    std::string folder() const { return dirname; }

private:
    enum { TS=0, TSI=1, TF=2, TFI=3, TC=4, TCI=5, TB=6 };
    struct rec { std::string name; int type; long long n; };

    void put(const std::string&, int, const void*, long long, size_t);
    long long count_of(int) const;
    bool get(const std::string&, int, void*, size_t);

    lexer *p = nullptr;
    std::string dirname;
    std::ofstream out;
    std::vector<std::pair<std::string,double>> sc;
    std::vector<std::pair<std::string,std::vector<double>>> bl;

    // reading
    struct pbox { int oi, oj, ok, ni, nj, nk; };
    std::vector<pbox> box;
    std::map<std::string,double> rsc;
    std::map<std::string,std::vector<double>> rbl;
    std::vector<std::map<std::string,std::pair<long long,long long>>> pidx;   // per piece: name -> offset, n
    std::vector<int> ptype;
    std::map<std::string,int> rtype;
    bool same = false;
    int mine = -1;          // the piece of this rank on the same decomposition
    int GNX=0, GNY=0, GNZ=0;
};

#endif
