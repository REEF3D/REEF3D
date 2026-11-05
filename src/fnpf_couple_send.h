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
class slice;
class reefcouple_bodies;

// REEF3D::FNPF side of the runtime coupling (reefcouple, B 92 62 in the receiver):
// handshake with each receiver model, eta, u, v, w (and Fi) on the sigma nodes after every time step;
// feedback from the receivers (I 243 bodies, I 244 overlap box) before every time step
class fnpf_couple_send
{
public:
    fnpf_couple_send(lexer*, fdm_fnpf*, ghostcell*);
    ~fnpf_couple_send();
    
    // before the time step: the feedback samples of the receivers (waits for them), the overlap targets
    void pre_step(lexer*, fdm_fnpf*, ghostcell*);
    
    // after the time step (p->simtime already advanced): the state at p->simtime
    void send(lexer*, fdm_fnpf*, ghostcell*);
    
    // all receivers have stopped: the FNPF model can end
    bool done() const;
    
    // end markers, outstanding sends, stop messages
    void finish();
    
    // relaxation towards the overlap targets (I 244 of a receiver), called from iowave::eta_relax and
    // iowave::fifsf_relax in every stage
    static void relax_eta(lexer*, slice&);
    static void relax_fifsf(lexer*, slice&);
    
    // the bodies driven by a receiver (I 243): CoG, velocity, angular velocity, rotation at time t
    bool body_state(int, double, double*) const;
    
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
        bool stopped = false;         // stop received (rank 0)
        double twait = -1.0e20;       // the receiver waits for a level at this time (rank 0)
        bool fi = false;              // Fi in the data (FNPF receiver)
        std::vector<dest> d;
        std::vector<std::vector<double>> buf[8];
        std::vector<MPI_Request> req[8];
        std::vector<MPI_Request> endreq;
        double endval = 0.0;
        
        // feedback
        int need = 1;                 // levels the receiver interpolates ahead (1 linear, 2 cubic)
        int nb = 0;                   // bodies
        int nc = 0;                   // overlap columns
        int fbsize = 0;               // doubles per sample
        std::vector<int> ci,cj;       // overlap columns: global FNPF indices (sorted by j, then i)
        std::vector<double> cr;       // relaxation weight of the targets
        std::vector<int> own;         // the overlap columns of this rank
        std::vector<double> s1,s2;    // the last two samples (t, dt, bodies, eta, u, v, w)
        std::vector<double> te,tf;    // targets eta, Fifsf at the overlap columns
        std::vector<double> phin;     // FNPF Fifsf at the overlap columns of the last step
        double tprev = -1.0e20;
        bool targets = false;
    };
    
    void handshake(lexer*, fdm_fnpf*, ghostcell*, link&);
    void pack(lexer*, fdm_fnpf*, const link&, const dest&, std::vector<double>&);
    void close(link&);
    bool poll(link&, double, double);
    void overlap_targets(lexer*, fdm_fnpf*, ghostcell*, link&);
    void relax(lexer*, slice&, bool);
    
    std::vector<link> L;
    std::vector<double> levt;         // times of the sent levels
    std::vector<double> Xg,Yg;        // global FNPF cell centres
    int level = 0;
    int nz = 0;
    int i0g,i1g,j0g,j1g;              // global index box of this rank (interior)
    int body_link = -1;               // the link whose receiver drives the bodies
    reefcouple_bodies *pbody = nullptr;
    lexer *p;
    ghostcell *pgc;
    bool finished = false;
};

#endif
