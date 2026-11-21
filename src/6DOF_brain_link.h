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

#ifndef SIXDOF_BRAIN_LINK_H_
#define SIXDOF_BRAIN_LINK_H_

#include<string>
#include<vector>

//  The link to an external brain (X 191 id port): a TCP connection of rank 0 to a process on this
//  machine (127.0.0.1:port, e.g. tools/fish_brain/brain_server.py), one line of JSON each way per
//  decision. REEF3D sends the senses (and the decision of its own brain), the external brain answers
//  with the gait targets. POSIX sockets; other systems: not available.

class sixdof_brain_link
{
public:

    ~sixdof_brain_link();

    // connect, retrying for up to wait seconds; false and a message if it fails
    bool connect_to(int port, double wait, std::string &msg);

    // one line out, one line back (without the newline); false if the link broke or timed out
    bool exchange(const std::string &out, std::string &in, double timeout);

    // the number after "key": in a line of JSON (def if absent)
    static double number(const std::string &json, const char *key, double def);
    static bool has(const std::string &json, const char *key);

private:

    int fd_ = -1;
    std::string buf_;
};

#endif
