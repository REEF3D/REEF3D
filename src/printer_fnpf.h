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

#ifndef PRINTER_FNPF_H_
#define PRINTER_FNPF_H_

#include"printer.h"
#include"increment.h"

#include"vtks.h"

class fdm_fnpf;
class fnpf_force_ale;
class potentialfile_out;
class ioflow;
class fnpf_print_wsf;
class fnpf_print_wsf_theory;
class fnpf_print_wsfline;
class fnpf_print_wsfline_y;
class fnpf_vtp_fsf;
class fnpf_vtp_bed;
class fnpf_state;
class fnpf_breaking_log;
class fnpf_print_Hs;
class fnpf_vel_probe;
class fnpf_vel_probe_theory;
class fnpf_runup;
class fnpf_print_kinematics;
class fnpf_print_scour;
class lagoon_output;

using namespace std;

class printer_fnpf final : public printer, public increment
{

public:
    void amr_columns(lexer*, const amr_column_lookup*) override;

    printer_fnpf(lexer*,fdm_fnpf*,ghostcell*);
    virtual ~printer_fnpf() = default;
    void start(lexer*,fdm_fnpf*,ghostcell*,ioflow*) override final;
    void print_stop(lexer*,fdm_fnpf*,ghostcell*) override final;
    void print_end(lexer*,fdm_fnpf*,ghostcell*) override final;
    void ckpt_state(std::vector<double>&) override;
    void ckpt_restore(const std::vector<double>&) override;
    void flush_outputs() override;

private:
    void print(lexer*,fdm_fnpf*,ghostcell*);
    void print_fsf(lexer*,fdm_fnpf*,ghostcell*);
    int last_print = -1;    // the iteration of the last volume output (print_end)
    int last_fsf = -1;      // and of the last free surface output
    void parallel(lexer*,int);

    vtk3D *outputFormat;
    bool initial_print = true;
    size_t file_offset = 0;

    char name[200];
    int n,iin,offset[200];
    float ffn;
    int gcval_phi,gcval_phiext;
    double *printtime_wT;
    double *printfsftime_wT;

    int printcount;
    bool bed_printed = false;   // the bed of this run (also after a restart, I 45: a new LAGOON store)
    int nP35 = 0, nP185 = 0;        // sizes of printtime_wT, printfsftime_wT

    fnpf_print_wsf *pwsf = nullptr;
    fnpf_print_wsf_theory *pwsf_theory;
    fnpf_print_wsfline *pwsfline;
    fnpf_print_wsfline_y *pwsfline_y;
    potentialfile_out *ppotentialfile;
    fnpf_vtp_fsf *pfsf = nullptr;
    fnpf_vtp_bed *pbed = nullptr;
    fnpf_state *pstate;
    fnpf_breaking_log *pbreaklog;
    fnpf_force_ale **pforce_ale;
    fnpf_runup **prunup;
    fnpf_print_Hs *phs;
    fnpf_vel_probe *pvel;
    fnpf_vel_probe_theory *pveltheo;
    fnpf_print_kinematics **pkin;
    fnpf_print_scour *pscour = nullptr;  // P 83, P 84: scour gauges
    lagoon_output *plagoon = nullptr;  // P 18: LAGOON store of the volume
};

#endif
