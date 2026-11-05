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

#ifndef FNPF_LAPLACE_BODY_BC_H_
#define FNPF_LAPLACE_BODY_BC_H_

// Resolved bodies (fnpf_6DOF) in the FNPF Laplace assembly.
// Body nodes (FBF > 0) are decoupled identity rows. On fluid/body faces the staircase
// Neumann condition d(f)/dx_face = (FBu,FBv,FBw).e_face is imposed by eliminating the
// ghost value, as for the A329 inflow; vertical faces use d(f)/dz = sigz*d(f)/dsig.
// Called after the KBEDBC block, which would otherwise re-couple the bed-slope terms.

struct fnpf_body_bc
{
    const bool on;
    const double *const FBF;
    const double *const FBu;
    const double *const FBv;
    const double *const FBw;
    const int *const flag7;
    
    fnpf_body_bc(bool on_, const double *F, const double *U, const double *V, const double *W, const int *fl)
        : on(on_), FBF(F), FBu(U), FBv(V), FBw(W), flag7(fl) {}
    
    inline bool body(int q) const
    {
        return on && FBF[q]>0.5;
    }
    
    inline bool face(int r) const
    {
        return FBF[r]>0.5 && flag7[r]>0;
    }
    
    // dxs/dxn/dye/dyw: node distances to the i-1/i+1/j-1/j+1 neighbours,
    // dzt/dzb: sigma distances to k+1/k-1, sz: sigz of the column
    inline void faces(int q, int sI, int sJ,
                      double dxs, double dxn, double dye, double dyw, double dzt, double dzb, double sz,
                      double &mp, double &ms, double &mn, double &me, double &mw, double &mt, double &mb,
                      double &rv) const
    {
        if(!on)
        return;
        
        if(face(q-sI))
        {
        rv += ms*FBu[q-sI]*dxs;
        mp += ms;
        ms = 0.0;
        }
        
        if(face(q+sI))
        {
        rv -= mn*FBu[q+sI]*dxn;
        mp += mn;
        mn = 0.0;
        }
        
        if(face(q-sJ))
        {
        rv += me*FBv[q-sJ]*dye;
        mp += me;
        me = 0.0;
        }
        
        if(face(q+sJ))
        {
        rv -= mw*FBv[q+sJ]*dyw;
        mp += mw;
        mw = 0.0;
        }
        
        if(face(q+1))
        {
        rv -= mt*FBw[q+1]*dzt/sz;
        mp += mt;
        mt = 0.0;
        }
        
        if(face(q-1))
        {
        rv += mb*FBw[q-1]*dzb/sz;
        mp += mb;
        mb = 0.0;
        }
    }
};

#endif
