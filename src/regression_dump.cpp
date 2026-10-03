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

#include"regression_dump.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"concentration.h"
#include"fdm_nhf.h"
#include"fdm_fnpf.h"
#include<cstdlib>
#include<cstring>
#include<cstdint>
#include<iomanip>
#include<sstream>
#include<sys/stat.h>
#include<sys/types.h>

regression_dump::regression_dump(lexer *p) : is_active(false), every(0), last_written(-1)
{
    const char *d = std::getenv("REEF3D_REGRESSION_DIR");

    if(d==nullptr || d[0]=='\0')
    return;

    is_active = true;
    dir = d;

    const char *e = std::getenv("REEF3D_REGRESSION_EVERY");
    if(e!=nullptr)
    every = std::atoi(e);

    mkdir(dir.c_str(),0777);

    std::ostringstream fn;
    fn<<dir<<"/steps_r"<<p->mpirank<<".txt";
    steplog.open(fn.str().c_str());
    steplog<<std::hexfloat;

    if(p->mpirank==0)
    std::cout<<"regression_dump: active, writing to "<<dir<<" every "<<every<<std::endl;
}

regression_dump::~regression_dump()
{
    if(steplog.is_open())
    steplog.close();
}

void regression_dump::add(const char *name, int n)
{
    names.push_back(name);
    data.emplace_back();
    data.back().reserve(n);
}

void regression_dump::cfd_collect(lexer *p, fdm *a, turbulence *pturb, concentration *pconc)
{
    names.clear();
    data.clear();

    add("u",p->cellnum);
    ULOOP
    data.back().push_back(a->u(i,j,k));

    add("v",p->cellnum);
    VLOOP
    data.back().push_back(a->v(i,j,k));

    add("w",p->cellnum);
    WLOOP
    data.back().push_back(a->w(i,j,k));

    add("press",p->cellnum);
    LOOP
    data.back().push_back(a->press(i,j,k));

    add("phi",p->cellnum);
    LOOP
    data.back().push_back(a->phi(i,j,k));

    add("eddyv",p->cellnum);
    LOOP
    data.back().push_back(a->eddyv(i,j,k));

    add("kin",p->cellnum);
    LOOP
    data.back().push_back(pturb->kinval(i,j,k));

    add("eps",p->cellnum);
    LOOP
    data.back().push_back(pturb->epsval(i,j,k));

    add("ro",p->cellnum);
    LOOP
    data.back().push_back(a->ro(i,j,k));

    add("visc",p->cellnum);
    LOOP
    data.back().push_back(a->visc(i,j,k));

    add("topo",p->cellnum);
    LOOP
    data.back().push_back(a->topo(i,j,k));

    add("conc",p->cellnum);
    LOOP
    data.back().push_back(pconc->val(i,j,k));

    // discrete divergence of the velocity field, as seen by the pressure projection
    add("div",p->cellnum);
    LOOP
    data.back().push_back((a->u(i,j,k)-a->u(i-1,j,k))/p->DXN[IP]
                         +(a->v(i,j,k)-a->v(i,j-1,k))/p->DYN[JP]*p->y_dir
                         +(a->w(i,j,k)-a->w(i,j,k-1))/p->DZN[KP]);

    add("vol",p->cellnum);
    LOOP
    data.back().push_back(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);

    add("flag4",p->cellnum);
    BASELOOP
    data.back().push_back(double(p->flag4[IJK]));
}

void regression_dump::write_state(lexer *p)
{
    std::ostringstream fn;
    fn<<dir<<"/state_"<<std::setw(8)<<std::setfill('0')<<p->count<<"_r"<<p->mpirank<<".bin";

    std::ofstream out(fn.str().c_str(), std::ios::binary);

    const char magic[8] = {'R','3','D','R','E','G','0','1'};
    out.write(magic,8);

    int32_t ival;
    double dval;

    ival=p->mpirank;  out.write((char*)&ival,sizeof(ival));
    ival=p->M10;      out.write((char*)&ival,sizeof(ival));
    ival=p->count;    out.write((char*)&ival,sizeof(ival));
    dval=p->simtime;  out.write((char*)&dval,sizeof(dval));
    ival=int32_t(names.size()); out.write((char*)&ival,sizeof(ival));

    for(size_t q=0; q<names.size(); ++q)
    {
        char name[16];
        std::memset(name,0,16);
        std::strncpy(name,names[q].c_str(),15);
        out.write(name,16);

        int64_t num = int64_t(data[q].size());
        out.write((char*)&num,sizeof(num));
        out.write((char*)data[q].data(),num*sizeof(double));
    }

    out.close();

    last_written=p->count;
}

void regression_dump::cfd_state(lexer *p, fdm *a, turbulence *pturb, concentration *pconc)
{
    if(last_written==p->count)
    return;

    cfd_collect(p,a,pturb,pconc);
    write_state(p);
}

void regression_dump::cfd_ini(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb, concentration *pconc)
{
    if(!is_active)
    return;

    steplog<<"# count simtime dt  sum(u^2) sum(v^2) sum(w^2) sum(press^2) sum(phi^2) sum(eddyv^2)  (rank-local, hexfloat)"<<std::endl;

    cfd_state(p,a,pturb,pconc);
}

void regression_dump::cfd_step(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb, concentration *pconc)
{
    if(!is_active)
    return;

    double su=0.0, sv=0.0, sw=0.0, sp=0.0, sphi=0.0, sev=0.0;

    ULOOP
    su += a->u(i,j,k)*a->u(i,j,k);

    VLOOP
    sv += a->v(i,j,k)*a->v(i,j,k);

    WLOOP
    sw += a->w(i,j,k)*a->w(i,j,k);

    LOOP
    {
    sp   += a->press(i,j,k)*a->press(i,j,k);
    sphi += a->phi(i,j,k)*a->phi(i,j,k);
    sev  += a->eddyv(i,j,k)*a->eddyv(i,j,k);
    }

    steplog<<p->count<<" "<<p->simtime<<" "<<p->dt<<"  "
           <<su<<" "<<sv<<" "<<sw<<" "<<sp<<" "<<sphi<<" "<<sev<<"\n";

    if(every>0 && p->count%every==0)
    cfd_state(p,a,pturb,pconc);
}

void regression_dump::cfd_final(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb, concentration *pconc)
{
    if(!is_active)
    return;

    steplog.flush();
    cfd_state(p,a,pturb,pconc);
}

// ---------------------------------------------------------------------
// NHFLOW
// ---------------------------------------------------------------------

void regression_dump::nhflow_collect(lexer *p, fdm_nhf *d)
{
    names.clear();
    data.clear();

    add("U",p->cellnum);
    LOOP
    data.back().push_back(d->U[IJK]);

    add("V",p->cellnum);
    LOOP
    data.back().push_back(d->V[IJK]);

    add("W",p->cellnum);
    LOOP
    data.back().push_back(d->W[IJK]);

    add("UH",p->cellnum);
    LOOP
    data.back().push_back(d->UH[IJK]);

    add("VH",p->cellnum);
    LOOP
    data.back().push_back(d->VH[IJK]);

    add("WH",p->cellnum);
    LOOP
    data.back().push_back(d->WH[IJK]);

    add("P",p->cellnum);
    LOOP
    data.back().push_back(d->P[IJK]);

    add("eta",p->cellnum);
    SLICELOOP4
    data.back().push_back(d->eta(i,j));

    add("WL",p->cellnum);
    SLICELOOP4
    data.back().push_back(d->WL(i,j));
}

void regression_dump::nhflow_state(lexer *p, fdm_nhf *d)
{
    if(last_written==p->count)
    return;

    nhflow_collect(p,d);
    write_state(p);
}

void regression_dump::nhflow_ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(!is_active)
    return;

    steplog<<"# count simtime dt  sum(U^2) sum(V^2) sum(W^2) sum(P^2) sum(eta^2) sum(WL^2)  (rank-local, hexfloat)"<<std::endl;

    nhflow_state(p,d);
}

void regression_dump::nhflow_step(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(!is_active)
    return;

    double su=0.0, sv=0.0, sw=0.0, sp=0.0, se=0.0, swl=0.0;

    LOOP
    {
    su += d->U[IJK]*d->U[IJK];
    sv += d->V[IJK]*d->V[IJK];
    sw += d->W[IJK]*d->W[IJK];
    sp += d->P[IJK]*d->P[IJK];
    }

    SLICELOOP4
    {
    se  += d->eta(i,j)*d->eta(i,j);
    swl += d->WL(i,j)*d->WL(i,j);
    }

    steplog<<p->count<<" "<<p->simtime<<" "<<p->dt<<"  "
           <<su<<" "<<sv<<" "<<sw<<" "<<sp<<" "<<se<<" "<<swl<<"\n";

    if(every>0 && p->count%every==0)
    nhflow_state(p,d);
}

void regression_dump::nhflow_final(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(!is_active)
    return;

    steplog.flush();
    nhflow_state(p,d);
}

// ---------------------------------------------------------------------
// FNPF
// ---------------------------------------------------------------------

void regression_dump::fnpf_collect(lexer *p, fdm_fnpf *c)
{
    names.clear();
    data.clear();

    add("Fi",p->cellnum);
    FLOOP
    data.back().push_back(c->Fi[FIJK]);

    add("U",p->cellnum);
    FLOOP
    data.back().push_back(c->U[FIJK]);

    add("V",p->cellnum);
    FLOOP
    data.back().push_back(c->V[FIJK]);

    add("W",p->cellnum);
    FLOOP
    data.back().push_back(c->W[FIJK]);

    add("eta",p->cellnum);
    SLICELOOP4
    data.back().push_back(c->eta(i,j));

    add("Fifsf",p->cellnum);
    SLICELOOP4
    data.back().push_back(c->Fifsf(i,j));

    add("WL",p->cellnum);
    SLICELOOP4
    data.back().push_back(c->WL(i,j));
}

void regression_dump::fnpf_state(lexer *p, fdm_fnpf *c)
{
    if(last_written==p->count)
    return;

    fnpf_collect(p,c);
    write_state(p);
}

void regression_dump::fnpf_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(!is_active)
    return;

    steplog<<"# count simtime dt  sum(Fi^2) sum(U^2) sum(W^2) sum(eta^2) sum(Fifsf^2) sum(WL^2)  (rank-local, hexfloat)"<<std::endl;

    fnpf_state(p,c);
}

void regression_dump::fnpf_step(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(!is_active)
    return;

    double sfi=0.0, su=0.0, sw=0.0, se=0.0, sff=0.0, swl=0.0;

    FLOOP
    {
    sfi += c->Fi[FIJK]*c->Fi[FIJK];
    su  += c->U[FIJK]*c->U[FIJK];
    sw  += c->W[FIJK]*c->W[FIJK];
    }

    SLICELOOP4
    {
    se  += c->eta(i,j)*c->eta(i,j);
    sff += c->Fifsf(i,j)*c->Fifsf(i,j);
    swl += c->WL(i,j)*c->WL(i,j);
    }

    steplog<<p->count<<" "<<p->simtime<<" "<<p->dt<<"  "
           <<sfi<<" "<<su<<" "<<sw<<" "<<se<<" "<<sff<<" "<<swl<<"\n";

    if(every>0 && p->count%every==0)
    fnpf_state(p,c);
}

void regression_dump::fnpf_final(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(!is_active)
    return;

    steplog.flush();
    fnpf_state(p,c);
}
