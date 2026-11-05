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

#ifndef REEFCOUPLE_H_
#define REEFCOUPLE_H_

#include<mpi.h>
#include<vector>

class lexer;
class ghostcell;
class fnpf_couple_send;
class wave_lib_fnpf_couple;

/*--------------------------------------------------------------------
REEF3D runtime coupling FNPF -> CFD / NHFLOW / FNPF (B 92 62), with
optional feedback to FNPF (I 243 bodies, I 244 overlap box)

The models run side by side in one MPI job (reef3d_mpmd), each with its
own working directory:

  mpirun -np 2 -wdir fnpf REEF3D : -np 8 -wdir cfd REEF3D

  source:    an FNPF model (A 10 3): the one the receiver names with
             I 246 <app>, by default the lowest app with A 10 3 that is
             no receiver; it needs no coupling key, it sends as soon as
             a receiver exists. An FNPF receiver can be a source itself
             (chain: regional FNPF -> local FNPF -> CFD, I 246 in the CFD)
  receivers: CFD (A 10 6), NHFLOW (A 10 5) or FNPF (A 10 3, nesting)
             models with the wave type B 92 62; several receivers can
             share one FNPF source

The receiver is a wave library (wave_lib_fnpf_couple): the relaxation
zones (B 98 2, B 520 zones on any edge), Dirichlet inflow and the full
initialisation (I 30 1) of the receiver use the FNPF eta, velocities
(and Fi for an FNPF receiver) instead of a wave theory.

  frame:  receiver point (x,y) -> wave frame of B 105 (angle, origin)
          -> FNPF point (xg + I 231, yg + I 232), z + I 233; the
          velocities are rotated back by the B 105 angle
  time:   FNPF time = receiver wave time + I 241
          I 242: time interpolation, 1 linear, 3 cubic (default)
  2D/3D:  FNPF 2D (j_dir 0) -> 3D receiver: long-crested in the B 105
          direction; FNPF 3D -> 2D receiver: the FNPF line y = I 232
          (rotated by B 105); 3D -> 3D, 2D -> 2D

Data exchange: at the start each receiver rank tells the FNPF ranks the
FNPF index window that covers its subdomain; each FNPF rank then sends
its part of the window (eta, u, v, w, Fi on the sigma nodes) after every
FNPF time step, MPI_Issend, at most reefcouple::depth levels ahead of
the receiver. The receiver interpolates in time between the levels.

Feedback (two-way, CFD and NHFLOW receivers): at the start of each of
its steps the receiver's rank 0 sends a sample (MPI_Isend) to FNPF
rank 0:
  I 243 1   the state of its 6DOF bodies; the FNPF bodies of the same
            number (X 10) follow them (sixdof_external)
  I 244 box (CFD) eta and surface velocity at the FNPF columns inside
            the box xs xe ys ye (receiver frame); FNPF relaxes eta and
            Fifsf there towards them, ramp width I 245; Fifsf from the
            surface velocity along FNPF x lines (dPhi/dx = u + w deta/dx)
FNPF waits at the start of its step for a sample at least as late as
the level the receiver is waiting for, and extrapolates the last two
samples to its stage times.
--------------------------------------------------------------------*/

class reefcouple
{
public:
    // collective on MPI_COMM_WORLD, after lexer_read: roles of the apps, intercommunicators
    static void roles(lexer*, ghostcell*);

    // called by ghostcell::final: closes the exchange (end markers, stop, drain)
    static void finish();

    static bool source() {return !recv_ic.empty();}
    static bool receiver() {return src_ic!=MPI_COMM_NULL;}

    static MPI_Comm src_ic;                 // receiver: intercommunicator to the FNPF source
    static std::vector<MPI_Comm> recv_ic;   // source: intercommunicators to the receivers
    static std::vector<int> recv_app;       // source: app numbers of the receivers
    static std::vector<int> recv_model;     // source: A 10 of the receivers

    static fnpf_couple_send *sender;
    static wave_lib_fnpf_couple *recv;

    // protocol
    static const int version = 2;
    static const int depth = 4;             // levels the source may run ahead
    static const int tag_bed  = 4811;
    static const int tag_data = 4812;
    static const int tag_end  = 4813;
    static const int tag_stop = 4814;
    static const int tag_fb   = 4815;
    static const int tag_fbc  = 4816;
};

#endif
