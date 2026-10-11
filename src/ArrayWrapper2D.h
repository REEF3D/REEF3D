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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#ifndef ARRAYWRAPPER2D_H_
#define ARRAYWRAPPER2D_H_

#include <cstddef>
#include <vector>

class lexer;

class ArrayWrapper2D final
{
public:
    explicit ArrayWrapper2D(lexer *pp);

    ArrayWrapper2D(const ArrayWrapper2D&)            = delete;
    ArrayWrapper2D &operator=(const ArrayWrapper2D&) = delete;
    ArrayWrapper2D(ArrayWrapper2D&&)                 = delete;
    ArrayWrapper2D &operator=(ArrayWrapper2D&&)      = delete;

    ~ArrayWrapper2D();

    // Deferred sizing, like ArrayWrapper3D::resize(): the lexer builds its
    // wrappers before imax/jmax are final.
    // -10 is a sentinel value for "uninitialized" (the default)
    void resize(int default_value = -10);

    void setVal(int val, bool includeGhost = false);

    inline int &operator()(int ii, int jj) noexcept;
    inline const int &operator()(int ii, int jj) const noexcept;

    // IJ-style flat access, level 0 only — the same escape hatch (and the same
    // restriction) as ArrayWrapper3D::operator[]. IJ is
    // (i-imin)*jmax + (j-jmin), so decode against jmax.
    inline int &operator[](int index) noexcept;
    inline const int &operator[](int index) const noexcept;

    // unshifted: callers index the raw pointer with IJ, not with i/j. Not an
    // implicit int* conversion, see ArrayWrapper3D::data().
    int *data() noexcept {return m_data.data();}
    const int *data() const noexcept {return m_data.data();}

    // whole array including ghost cells; 0 before resize()
    std::size_t size() const noexcept {return m_data.size();}
    int *begin() noexcept {return m_data.data();}
    const int *begin() const noexcept {return m_data.data();}
    int *end() noexcept {return m_data.data()+m_data.size();}
    const int *end() const noexcept {return m_data.data()+m_data.size();}

private:
    /// Recomputes the cached addressing below. Must be called after anything
    /// that resizes @p m_data. See ArrayWrapper3D::cache_addressing for why the
    /// origin is folded into the pointer and why the stride is a long.
    void cache_addressing() noexcept;

    lexer *const p;

    std::vector<int> m_data;

    int *m_base = nullptr; ///< origin-folded base: m_base[ii*m_js + jj]
    long m_js   = 0;       ///< i-stride (jmax); long on purpose
};

#endif
