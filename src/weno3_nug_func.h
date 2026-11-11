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

#ifndef WENO3_NUG_FUNC_H_
#define WENO3_NUG_FUNC_H_

#include "increment.h"
#include <array>
#include <vector>

class lexer;
class fdm;
class field;

using namespace std;

class weno3_nug_func : public increment
{
public:
    weno3_nug_func(lexer*);
    virtual ~weno3_nug_func() = default;

    void precalc_qf(lexer*);
    void precalc_cf(lexer*);
    void precalc_isf(lexer*);

    void ini(lexer*);

    void is_min_x();
    void is_min_y();
    void is_min_z();

    void is_max_x();
    void is_max_y();
    void is_max_z();

    void weight_min_x();
    void weight_min_y();
    void weight_min_z();

    void weight_max_x();
    void weight_max_y();
    void weight_max_z();

    using qf_t = std::array<std::array<std::array<double, 2>, 4>, 2>;
    using cf_t = std::array<std::array<double, 4>, 2>;

    // coefficient tables: shared by the instances on the rank grid (static storage, built for the
    // first lexer); an instance on another grid (a CFD mesh refinement patch) builds its own, as
    // weno_nug_func. The pointers address the tables in use.
    qf_t *qfx = nullptr, *qfy = nullptr, *qfz = nullptr;
    cf_t *cfx = nullptr, *cfy = nullptr, *cfz = nullptr;
    cf_t *isfx = nullptr, *isfy = nullptr, *isfz = nullptr;

    static inline std::vector<qf_t> s_qfx, s_qfy, s_qfz;
    static inline std::vector<cf_t> s_cfx, s_cfy, s_cfz;
    static inline std::vector<cf_t> s_isfx, s_isfy, s_isfz;

    std::vector<qf_t> own_qfx, own_qfy, own_qfz;
    std::vector<cf_t> own_cfx, own_cfy, own_cfz;
    std::vector<cf_t> own_isfx, own_isfy, own_isfz;

    static inline lexer *s_lexer = nullptr;          // the lexer the shared tables were built for
    void own_ini(lexer*);

    double q1,q2,q3;

    double is1x,is2x,is3x;
    double is1y,is2y,is3y;
    double is1z,is2z,is3z;
    double w1x,w2x,w3x;
    double w1y,w2y,w3y;
    double w1z,w2z,w3z;

    int uf,vf,wf;

    static inline bool iniflag = false;

    static constexpr double epsilon = 0.0;
    static constexpr double psi = 1.0e-6;
};

#endif
