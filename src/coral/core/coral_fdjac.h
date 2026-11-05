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
// REEF3D::CORAL: finite-difference Jacobian of the residual of the multibody system with respect to the
// generalized coordinates (6 per body, then the scalar coordinates), perturbing the nodes of one colour
// together (distance-2 colouring of the coupling graph, model::coloring): rows of node i only change
// with i and its neighbours, constraint rows only with the bodies of the constraint. The entries are
// added as triplets (row, column, value). Used by the time integrator, the static solver and the
// linearisation; with fd_only the loads with an analytic tangent are left out of the colouring (the
// residual must then leave them out as well).

#ifndef CORAL_FDJAC_H_
#define CORAL_FDJAC_H_

#include"coral/core/coral_model.h"
#include"coral/core/coral_types.h"
#include<vector>

namespace coral
{
    // T: triplets of the rows n + nc and the first n columns. perturb(col, eps) sets/clears a
    // perturbation of generalized coordinate col (eps = 0: clear); res(r) evaluates the residual of the
    // current perturbation; eps_of(col) gives the step; r0 the unperturbed residual. central: central
    // differences (two evaluations per colour and component).
    template<class Perturb, class Residual, class Eps>
    void fd_jacobian(const model &m, std::vector<triplet> &T, const VecX &r0, Perturb perturb, Residual res, Eps eps_of,
                     bool central, bool fd_only=false)
    {
        const int nb = m.nbody(), nn = nb + m.ncoord(), n = m.ndof();
        std::vector<int> color;
        std::vector<std::vector<int>> nbr;
        const int ncol = m.coloring(color, nbr, fd_only);
        // constraint rows of each body
        std::vector<std::vector<int>> crow(nn);
        int row = n;
        for(const auto &c : m.constraint_list)
        {
            for(int q=0; q<c->size(); ++q)
            {
                if(c->body_a >= 0) crow[c->body_a].push_back(row + q);
                if(c->body_b >= 0) crow[c->body_b].push_back(row + q);
            }
            row += c->size();
        }
        VecX rp, rm;
        std::vector<int> group;
        for(int c=0; c<ncol; ++c)
        {
            group.clear();
            for(int j=0; j<nn; ++j)
                if(color[j]==c) group.push_back(j);
            for(int k=0; k<6; ++k)
            {
                bool any = false;
                for(int j : group)
                    if(k < m.node_dofs(j)) {perturb(m.node_dof0(j)+k, eps_of(m.node_dof0(j)+k)); any = true;}
                if(!any)
                    break;
                res(rp);
                if(central)
                {
                    for(int j : group)
                        if(k < m.node_dofs(j)) perturb(m.node_dof0(j)+k, -eps_of(m.node_dof0(j)+k));
                    res(rm);
                }
                for(int j : group)
                    if(k < m.node_dofs(j)) perturb(m.node_dof0(j)+k, 0.0);
                for(int j : group)
                {
                    if(k >= m.node_dofs(j))
                        continue;
                    const int col = m.node_dof0(j)+k;
                    const double e = eps_of(col);
                    auto add = [&](int r)
                    {
                        const double v = central ? (rp(r) - rm(r))/(2.0*e) : (rp(r) - r0(r))/e;
                        if(v != 0.0) T.emplace_back(r, col, v);
                    };
                    for(int i : nbr[j])
                        for(int q=0; q<m.node_dofs(i); ++q)
                            add(m.node_dof0(i)+q);
                    for(int q : crow[j])
                        add(q);
                }
            }
        }
    }
}

#endif
