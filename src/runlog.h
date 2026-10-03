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

#ifndef RUNLOG_H_
#define RUNLOG_H_

#include <fstream>
#include <set>
#include <string>

class lexer;

// The run log: which outputs belong to which run, and at which time.
// Specification: LAGOON, lagoon-format/RUN.md, version 1.
//
// One file per solver, REEF3D_Case/REEF3D_<SOLVER>_run.jsonl, one JSON object per
// line, append only: a "run" line at the start, a "stream" line when a writer
// first writes, an "output" line per output and an "end" line at the end.
// Only rank 0 writes; on the other ranks every call returns at once.
//
// Writers call it after the parallel file of an output is complete, e.g.
//
//     if(p->plog)
//     p->plog->written(p,num,"fsf","free_surface",name,p->M10);

class runlog
{
public:
    runlog(lexer*, const char *version);
    ~runlog();

    // a writer finished the file of output number step: path is the file it wrote
    // (the .pvtu/.pvtp for parallel output). Folder, name pattern and format are
    // taken from the path, the stream is registered on first use.
    void written(lexer*, int step, const char *stream, const char *role,
                 const char *path, int pieces=0);

    // the same with folder, pattern ({step:08d}) and format given explicitly
    void output(lexer*, int step, const char *stream, const char *role,
                const char *folder, const char *files, const char *format, int pieces);

    // a table that grows over the run (gauges, forces, motions), or a static file
    void table(lexer*, const char *stream, const char *role, const char *folder,
               const char *file, const char *format="dat");

    // the same, from the path of the file
    void table_file(lexer*, const char *stream, const char *role, const char *path);

    // the run ends: "finished", "stopped" or "error"
    void end(lexer*, const char *status);

    const std::string& id() const {return run_id;}

    static std::string solver_name(lexer*);
    static std::string sha256_file(const std::string &path);

private:
    void start(lexer*);
    void write(const std::string &line);
    void stream(const char *stream, const std::string &json);
    static std::string quote(const std::string&);
    static std::string number(double);
    static std::string utc(const char *format);
    static std::string format_of(const std::string &file);

    std::ofstream out;
    std::set<std::string> known;
    std::string run_id, version_text, ctrl_copy_path, started_text;
    bool active;
    bool ended;
    bool started;
};

#endif
