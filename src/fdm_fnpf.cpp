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

#include"fdm_fnpf.h"
#include"lexer.h"

namespace
{
    /*!
     * @brief Rows the FNPF Laplace assembly actually writes.
     *
     * fnpf_laplace_cds2 is the only file on this path that touches M or
     * rhsvec, and both of its passes run a plain counter from 0 across LOOP.
     * Every solver's fill then reads them back through a counter built from
     * the same LOOP - reefmg's fill_matrix8, hypre_struct_fill8,
     * hypre_aij_F_fill and hypre_sstruct_fill8 all do - so nothing indexes
     * past this. lexer's veclength adds a ghost allowance on top (gcbnum plus
     * gcpara_sum, inflated further by the safety margins in
     * lexer_vectorsize), which no FNPF path ever reads: on a 133x200x10 local
     * box that was 60,360 of 326,360 rows, 18.5%.
     *
     * Safe because flag4 is fixed for the run - no fnpf_* source writes it -
     * so the count cannot grow after construction.
     */
    int fnpf_laplace_rows(lexer *p)
    {
        int i,j,k,n=0;
        LOOP ++n;
        return n;
    }

    /*!
     * @brief Rows the 2D free-surface damping system writes.
     *
     * fnpf_fsfbc::damping fills N and rvec with a counter run from 0 across
     * SLICELOOP4, and the solver behind it - sflow_bicgstab, constructed in
     * fnpf_fsfbc - indexes xvec and rhsvec the same way: for var 4 it sets
     * flagslice to p->flagslice4 with ulast and vlast both zero, which makes
     * its SLICEFLEXLOOP exactly SLICELOOP4.
     *
     * vec2D's default constructor takes lexer's veclength, the 3D row count,
     * where matrix2D correctly takes vec2Dlength; on a 133x200x10 local box
     * that made xvec and rvec 326,360 entries each for 26,600 rows of data.
     */
    int fnpf_damping_rows(lexer *p)
    {
        //  Both consumers - fnpf_fsfbc::damping and the wetting-drying
        //  variant - are wrapped in `if(p->A350>0 && damping)`, and they are
        //  the only files in the tree that touch N, xvec or rvec.  With
        //  breaking off the system is never assembled or solved, so it needs
        //  no storage at all.
        if(p->A350<=0)
        return 0;

        int i,j,n=0;
        SLICELOOP4 ++n;
        return n;
    }
}

fdm_fnpf::fdm_fnpf(lexer *p) : test(p),Fifsf(p),
                              nodeval(p),eta(p),etaloc(p),
                              wet_n(p),breaking(p),breaklog(p),bc(p),
                              eta_n(p),WL(p),bed(p),depth(p),Fz(p),K(p),
                              Fx(p),Fy(p),
                              Ex(p),Ey(p),Exx(p),Eyy(p),
                              Bx(p),By(p),Bxx(p),Byy(p),
                              coastline(p),vb(p),
                              test2D(p),Hs(p),
                              nodeval2D(p),breaking_print(p),
                              laprows(fnpf_laplace_rows(p)),
                              slicerows(fnpf_damping_rows(p)),
                              rhsvec(p,laprows),rvec(p,slicerows),xvec(p,slicerows),
                              N(p,slicerows),M(p,laprows)
{   
    p->Darray(p->sig,p->imax*p->jmax*(p->kmax+2));
    p->Darray(p->sigx,p->imax*p->jmax*(p->kmax+2));
    p->Darray(p->sigy,p->imax*p->jmax*(p->kmax+2));
    p->Darray(p->sigz,p->imax*p->jmax);
    p->Darray(p->sigxx,p->imax*p->jmax*(p->kmax+2));
    
    p->Darray(U,p->imax*p->jmax*(p->kmax+2));
    p->Darray(V,p->imax*p->jmax*(p->kmax+2));
    p->Darray(W,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fi,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Uin,p->imax*p->jmax*(p->kmax+2));
}













