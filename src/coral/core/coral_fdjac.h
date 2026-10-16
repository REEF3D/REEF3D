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
// 6 generalized coordinates of every body, perturbing the bodies of one colour together (distance-2
// colouring of the coupling graph, model::coloring): rows of body i only change with i and its
// neighbours, constraint rows only with the bodies of the constraint. Used by the time integrator and
// the static solver.

#ifndef CORAL_FDJAC_H_
#define CORAL_FDJAC_H_

#include"coral/core/coral_model.h"
#include"coral/core/coral_types.h"
#include<vector>

namespace coral
{
    // J: rows n + nc, the first n columns are filled. perturb(col, eps) sets/clears a perturbation
    // of generalized coordinate col (eps = 0: clear); res(r) evaluates the residual of the current
    // perturbation; eps_of(col) gives the step; r0 the unperturbed residual. central: central
    // differences (two evaluations per colour and component).
    template<class Perturb, class Residual, class Eps>
    void fd_jacobian(const model &m, MatX &J, const VecX &r0, Perturb perturb, Residual res, Eps eps_of, bool central)
    {
        const int nb = m.nbody(), n = 6*nb;
        std::vector<int> color;
        std::vector<std::vector<int>> nbr;
        const int ncol = m.coloring(color, nbr);
        // constraint rows of each body
        std::vector<std::vector<int>> crow(nb);
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
            for(int j=0; j<nb; ++j)
                if(color[j]==c) group.push_back(j);
            for(int k=0; k<6; ++k)
            {
                for(int j : group) perturb(6*j+k, eps_of(6*j+k));
                res(rp);
                if(central)
                {
                    for(int j : group) perturb(6*j+k, -eps_of(6*j+k));
                    res(rm);
                }
                for(int j : group) perturb(6*j+k, 0.0);
                for(int j : group)
                {
                    const int col = 6*j+k;
                    const double e = eps_of(col);
                    for(int i : nbr[j])
                        J.block(6*i, col, 6, 1) = central ? VecX((rp.segment(6*i,6) - rm.segment(6*i,6))/(2.0*e))
                                                          : VecX((rp.segment(6*i,6) - r0.segment(6*i,6))/e);
                    for(int q : crow[j])
                        J(q, col) = central ? (rp(q) - rm(q))/(2.0*e) : (rp(q) - r0(q))/e;
                }
            }
        }
    }
}

#endif
