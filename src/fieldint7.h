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

#ifndef FIELDINT7_H_
#define FIELDINT7_H_

#include "fieldint.h"

/*!
 * @brief Integer counterpart of field7.
 *
 * Sigma-grid vertical-node layout, stride p->kmaxF, addressed with the FIJK
 * family. See field7.h for what the slack plane is for.
 */
class fieldint7 final : public fieldint
{
public:
    fieldint7(lexer *pp) : fieldint(pp, pp->kmaxF,
                                    static_cast<std::size_t>(pp->imax)*static_cast<std::size_t>(pp->jmax)) {}
    ~fieldint7() override;

    /*!
     * FBASELOOP rather than LOOP, as in field7: LOOP stops at KMAX_LOOP and
     * would leave the top node plane untouched, and its PCHECK reads the
     * IJK-strided flag4. Shadows the non-virtual field_base::setVal — hold
     * fieldint7 by concrete type, never through a fieldint&.
     */
    void setVal(int val, bool includeGhost = false)
    {
        if(includeGhost)
        {
            std::fill(V.begin(), V.end(), val);
            return;
        }

        int i,j,k;
        FBASELOOP
        operator()(i,j,k) = val;
    }
};

#endif
