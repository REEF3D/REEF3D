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

#ifndef VTP_PIECE_BUFFER_H_
#define VTP_PIECE_BUFFER_H_

#include<ios>
#include<string>

// The appended (binary) data of a VTP piece made in memory: write() appends to a string, without
// the sentry and stream state of std::ostream::write for every value. The bytes are the same.
// The piece is then the header (an std::ostringstream), these bytes and the footer, and goes to
// lagoon_surface::piece, which writes the file where it is wanted and hands it to the LAGOON store.
struct vtp_piece_buffer
{
    std::string s;

    void write(const char *c, std::streamsize n)
    {
        s.append(c,size_t(n));
    }
};

#endif
