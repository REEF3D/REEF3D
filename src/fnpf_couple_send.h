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

#ifndef FNPF_COUPLE_SEND_H_
#define FNPF_COUPLE_SEND_H_

#include<mpi.h>
#include<vector>

class lexer;
class fdm_fnpf;
class ghostcell;

// REEF3D::FNPF side of the runtime coupling (reefcouple, B 92 62 in the receiver):
// handshake with each receiver model, eta and u,v,w on the sigma nodes after every time step
class fnpf_couple_send
{
public:
    fnpf_couple_send(lexer*, fdm_fnpf*, ghostcell*);
    ~fnpf_couple_send();
    
    // after the time step (p->simtime already advanced): the state at p->simtime
    void send(lexer*, fdm_fnpf*, ghostcell*);
    
    // all receivers have stopped: the FNPF model can end
    bool done() const;
    
    // end markers, outstanding sends, stop messages
    void finish();
    
private:
    struct dest
    {
        int rank;                     // rank in the receiver model
        int i0,i1,j0,j1;              // global FNPF index box sent to it
        int ncol;
    };
    
    struct link
    {
        MPI_Comm ic = MPI_COMM_NULL;
        int app = -1;
        bool open = true;
        std::vector<dest> d;
        std::vector<std::vector<double>> buf[8];
        std::vector<MPI_Request> req[8];
        std::vector<MPI_Request> endreq;
        double endval = 0.0;
    };
    
    void handshake(lexer*, fdm_fnpf*, ghostcell*, link&);
    void pack(lexer*, fdm_fnpf*, const dest&, std::vector<double>&);
    void close(link&);
    void check_stop(lexer*, ghostcell*);
    
    std::vector<link> L;
    int level = 0;
    int nz = 0;
    int i0g,i1g,j0g,j1g;              // global index box of this rank (interior)
    lexer *p;
    ghostcell *pgc;
    bool finished = false;
};

#endif
