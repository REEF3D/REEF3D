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


#ifndef GEOREF_H_
#define GEOREF_H_

#include<fstream>
#include<sstream>
#include<string>

// The georeference of a case: where its model coordinates lie on the Earth.
// Specification: LAGOON, lagoon-format/SPEC.md section 1.2.
//
// The case keeps it in ./georef.json (written by LAGOON's set-up or by hand), e.g.
//
//     {"version": 1, "crs": "EPSG:25833", "origin": [270000.0, 7040000.0],
//      "rotation_deg": 30.0, "vertical_datum": "EPSG:5941", "z_offset": 0.0}
//
// REEF3D does not compute with it: it copies the block unchanged into the run log's
// run record and into the root of the LAGOON stores it writes, so the results keep
// their place on the map. georef_json() gives the block as one line of JSON, or ""
// when there is no file (or it is not a JSON object).

inline std::string georef_json(const char *path="./georef.json")
{
    std::ifstream in(path, std::ios::binary);
    if(!in)
    return "";

    std::ostringstream text;
    text<<in.rdbuf();
    std::string s = text.str();

    if(s.size()>65536)
    return "";

    // one line: JSON allows any white space between its tokens, and a line break
    // inside a string is not JSON anyway
    for(char &c : s)
    if(c=='\n' || c=='\r' || c=='\t')
    c=' ';

    const size_t first = s.find_first_not_of(' ');
    const size_t last = s.find_last_not_of(' ');

    if(first==std::string::npos || s[first]!='{' || s[last]!='}' || s.find('\0')!=std::string::npos)
    return "";

    return s.substr(first, last-first+1);
}

#endif
