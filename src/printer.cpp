/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include "printer.h"

#include<sys/types.h>
#include<cstdio>
#include<cstring>

void printer::writeFile(const char* filename, const size_t file_size)
{
    const size_t SMALL_FILE_THRESHOLD = 10 * 1024 * 1024;   // 10MB

    if(file_size < SMALL_FILE_THRESHOLD)
    {
        // Small files: Simple fwrite with modest buffer
        FILE* file = fopen(filename, "wb");
        if(file)
        {
            setvbuf(file, nullptr, _IOFBF, 32768); // 32KB buffer
            fwrite(buffer.data(), buffer.size(), 1, file);
            fclose(file);
        }
    }
    else
    {
        // Medium files: Optimized C I/O with larger buffer
        FILE* file = fopen(filename, "wb");
        if(file)
        {
            setvbuf(file, nullptr, _IOFBF, 131072); // 128KB buffer
            fwrite(buffer.data(), buffer.size(), 1, file);
            fclose(file);
        }
    }
    // Large files: Use mmap or similar for efficient I/O, if available
}
