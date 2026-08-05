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

#ifndef ARRAYWRAPPER3D_H_
#define ARRAYWRAPPER3D_H_

#include <cstddef>
#include <vector>

class lexer;

class ArrayWrapper3D final
{
public:
    explicit ArrayWrapper3D(lexer *pp, unsigned int _data_location = 4);

    ArrayWrapper3D(const ArrayWrapper3D&)            = delete;
    ArrayWrapper3D &operator=(const ArrayWrapper3D&) = delete;
    ArrayWrapper3D(ArrayWrapper3D&&)                 = delete;
    ArrayWrapper3D &operator=(ArrayWrapper3D&&)      = delete;

    ~ArrayWrapper3D();

    void resize(int default_value = 0);

    void setVal(int val, bool includeGhost = false);

    inline int &operator()(int i, int j, int k) noexcept;
    inline const int &operator()(int i, int j, int k) const noexcept;

    inline int &operator[](int index) noexcept;
    inline const int &operator[](int index) const noexcept;

    // unshifted: callers index the raw pointer with IJK/FIJK, not with i/j/k.
    // Deliberately not an implicit int* conversion, so leftover raw-array code
    // (delete[], pointer tests, int* parameters) fails to compile instead of
    // silently taking the wrapper's storage.
    int *data() noexcept {return m_data.data();}
    const int *data() const noexcept {return m_data.data();}

    // whole array including ghost cells (and the slack plane for location 7);
    // 0 before resize()
    std::size_t size() const noexcept {return m_data.size();}
    int *begin() noexcept {return m_data.data();}
    const int *begin() const noexcept {return m_data.data();}
    int *end() noexcept {return m_data.data()+m_data.size();}
    const int *end() const noexcept {return m_data.data()+m_data.size();}

private:
    /// Recomputes the cached addressing below from @p m_data and the lexer grid
    /// metrics. Must be called after anything that resizes @p m_data.
    void cache_addressing() noexcept;

    lexer *const p;

    std::vector<int> m_data;

    int *m_base = nullptr; ///< origin-folded base: m_base[i*m_js + j*m_ks + k]
    long m_js   = 0;       ///< i-stride (jmax*kz());  long on purpose, see cache_addressing
    long m_ks   = 0;       ///< j-stride (kz());       long on purpose, see cache_addressing

    /// 1/2/3/4/7.
    unsigned int data_location = 0;
};

#endif
