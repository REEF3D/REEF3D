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

#include "runlog.h"
#include "lexer.h"
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <random>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

runlog::runlog(lexer *p, const char *version) : active(p->mpirank==0), ended(false), started(false)
{
    if(!active)
    return;

    mkdir("./REEF3D_Case",0777);

    // run id: start time (UTC) and a random suffix, so it sorts by time and is unique
    std::random_device seed;
    std::ostringstream id;
    id<<utc("%Y%m%dT%H%M%SZ")<<"-"<<std::hex<<std::setw(4)<<std::setfill('0')<<(seed()&0xffff);
    run_id = id.str();

    std::string solver = solver_name(p);
    out.open("./REEF3D_Case/REEF3D_"+solver+"_run.jsonl", std::ios::app);
    if(!out.is_open())
    {
        active=false;
        return;
    }

    // keep the control file of this run next to the log
    std::string ctrl_copy = "REEF3D_Case/ctrl-"+run_id+".txt";
    {
        std::ifstream src("./ctrl.txt", std::ios::binary);
        if(src.is_open())
        {
            std::ofstream dst("./"+ctrl_copy, std::ios::binary);
            dst<<src.rdbuf();
        }
    }

    version_text = version;
    ctrl_copy_path = ctrl_copy;
    started_text = utc("%Y-%m-%dT%H:%M:%SZ");
}

void runlog::start(lexer *p)
{
    // written lazily, at the first record: by then the solver has set up the still
    // water level and the domain
    if(started)
    return;
    started=true;

    std::string solver = solver_name(p);
    std::string version = version_text;
    std::string ctrl_copy = ctrl_copy_path;
    char host[256] = "";
    gethostname(host,sizeof(host)-1);

    std::ostringstream line;
    line<<"{\"type\":\"run\",\"format\":\"reef3d-run\",\"version\":1"
        <<",\"run\":"<<quote(run_id)
        <<",\"solver\":"<<quote(solver)
        <<",\"reef3d\":"<<quote(version)
        <<",\"started\":"<<quote(started_text)
        <<",\"host\":"<<quote(host)
        <<",\"ranks\":"<<p->M10
        <<",\"ctrl\":{\"file\":\"ctrl.txt\",\"copy\":"<<quote(ctrl_copy)
        <<",\"sha256\":"<<quote(sha256_file("./ctrl.txt"))<<"}"
        <<",\"units\":{\"length\":\"m\",\"time\":\"s\",\"mass\":\"kg\",\"force\":\"N\"}"
        <<",\"still_water_level\":"<<number(p->phimean)
        <<",\"domain\":{\"min\":["<<number(p->global_xmin)<<","<<number(p->global_ymin)<<","<<number(p->global_zmin)
        <<"],\"max\":["<<number(p->global_xmax)<<","<<number(p->global_ymax)<<","<<number(p->global_zmax)<<"]}"
        <<",\"streams\":{}}";
    write(line.str());
}

runlog::~runlog()
{
    if(active && !ended)
    out.close();
}

std::string runlog::solver_name(lexer *p)
{
    switch(p->A10)
    {
        case 2: return "SFLOW";
        case 3: return "FNPF";
        case 5: return "NHFLOW";
        case 6: return "CFD";
        default: return "REEF3D";
    }
}

void runlog::output(lexer *p, int step, const char *stream_name, const char *role,
                    const char *folder, const char *files, const char *format, int pieces)
{
    if(!active || ended)
    return;

    start(p);

    if(known.find(stream_name)==known.end())
    {
        std::ostringstream spec;
        spec<<"{\"role\":"<<quote(role)
            <<",\"folder\":"<<quote(folder)
            <<",\"files\":"<<quote(files)
            <<",\"format\":"<<quote(format);
        if(pieces>0)
        spec<<",\"pieces\":"<<pieces;
        if(step>0)
        spec<<",\"first_step\":"<<step;
        spec<<"}";
        stream(stream_name,spec.str());
    }

    std::ostringstream line;
    line<<"{\"type\":\"output\",\"run\":"<<quote(run_id)
        <<",\"step\":"<<step
        <<",\"time\":"<<number(p->simtime)
        <<",\"iteration\":"<<p->count
        <<",\"streams\":["<<quote(stream_name)<<"]}";
    write(line.str());
}

void runlog::written(lexer *p, int step, const char *stream_name, const char *role,
                     const char *path, int pieces)
{
    if(!active || ended)
    return;

    std::string full = path;
    if(full.rfind("./",0)==0)
    full = full.substr(2);
    size_t slash = full.rfind('/');
    std::string folder = slash==std::string::npos ? "." : full.substr(0,slash);
    std::string file = slash==std::string::npos ? full : full.substr(slash+1);

    // the pattern: the last run of digits that reads as the step becomes {step:0Nd}
    std::string pattern = file;
    size_t end = file.size();
    while(end>0)
    {
        size_t last = file.find_last_of("0123456789",end-1);
        if(last==std::string::npos)
        break;
        size_t first = file.find_last_not_of("0123456789",last);
        first = first==std::string::npos ? 0 : first+1;
        if(std::stol(file.substr(first,last-first+1))==step)
        {
            int width = int(last-first+1);
            std::string field = width>1 ? "{step:0"+std::to_string(width)+"d}" : "{step}";
            pattern = file.substr(0,first)+field+file.substr(last+1);
            break;
        }
        end = first;
    }
    output(p,step,stream_name,role,folder.c_str(),pattern.c_str(),format_of(file).c_str(),pieces);
}

std::string runlog::format_of(const std::string &file)
{
    size_t dot = file.rfind('.');
    std::string ext = dot==std::string::npos ? "" : file.substr(dot+1);
    if(ext=="pvtu" || ext=="vtu") return "vtu";
    if(ext=="pvtp" || ext=="vtp") return "vtp";
    if(ext=="pvts" || ext=="vts") return "vts";
    if(ext=="pvtr" || ext=="vtr") return "vtr";
    if(ext=="vtm") return "vtm";
    if(ext=="vtk") return "vtk-legacy";
    if(ext=="stl") return "stl";
    return ext.empty() ? "dat" : ext;
}

void runlog::table(lexer *p, const char *stream_name, const char *role, const char *folder,
                   const char *file, const char *format)
{
    if(!active || ended || known.find(stream_name)!=known.end())
    return;

    start(p);

    std::ostringstream spec;
    spec<<"{\"role\":"<<quote(role)
        <<",\"folder\":"<<quote(folder)
        <<",\"file\":"<<quote(file)
        <<",\"format\":"<<quote(format)<<"}";
    stream(stream_name,spec.str());
}

void runlog::table_file(lexer *p, const char *stream_name, const char *role, const char *path)
{
    std::string full = path;
    if(full.rfind("./",0)==0)
    full = full.substr(2);
    size_t slash = full.rfind('/');
    std::string folder = slash==std::string::npos ? "." : full.substr(0,slash);
    std::string file = slash==std::string::npos ? full : full.substr(slash+1);
    table(p,stream_name,role,folder.c_str(),file.c_str(),format_of(file).c_str());
}

void runlog::end(lexer *p, const char *status)
{
    if(!active || ended)
    return;

    start(p);

    std::ostringstream line;
    line<<"{\"type\":\"end\",\"run\":"<<quote(run_id)
        <<",\"time\":"<<number(p->simtime)
        <<",\"iteration\":"<<p->count
        <<",\"status\":"<<quote(status)
        <<",\"ended\":"<<quote(utc("%Y-%m-%dT%H:%M:%SZ"))<<"}";
    write(line.str());
    out.close();
    ended=true;
}

void runlog::stream(const char *stream_name, const std::string &json)
{
    if(!started)
    return;
    known.insert(stream_name);
    write("{\"type\":\"stream\",\"run\":"+quote(run_id)+",\"name\":"+quote(stream_name)
          +",\"stream\":"+json+"}");
}

void runlog::write(const std::string &line)
{
    // one whole line, flushed: a run that is killed leaves a valid log
    out<<line<<'\n'<<std::flush;
}

std::string runlog::quote(const std::string &text)
{
    std::string result="\"";
    for(char c : text)
    {
        if(c=='"' || c=='\\')
        {
            result+='\\';
            result+=c;
        }
        else if(static_cast<unsigned char>(c)<0x20)
        {
            char buffer[8];
            std::snprintf(buffer,sizeof(buffer),"\\u%04x",c);
            result+=buffer;
        }
        else
        result+=c;
    }
    return result+"\"";
}

std::string runlog::number(double value)
{
    // 17 significant digits: doubles survive the round trip
    if(value!=value || value>1.0e300 || value<-1.0e300)
    return "null";

    std::ostringstream text;
    text<<std::setprecision(17)<<value;
    return text.str();
}

std::string runlog::utc(const char *format)
{
    std::time_t now = std::time(nullptr);
    std::tm parts{};
    gmtime_r(&now,&parts);
    char buffer[64];
    std::strftime(buffer,sizeof(buffer),format,&parts);
    return buffer;
}

// ---------------------------------------------------------------- SHA-256
// FIPS 180-4, for the control file fingerprint (same set-up = same hash)

namespace
{
    const uint32_t K[64] = {
        0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
        0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
        0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
        0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
        0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
        0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
        0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
        0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};

    inline uint32_t rotr(uint32_t x, int n) {return (x>>n)|(x<<(32-n));}

    void block(uint32_t *h, const unsigned char *data)
    {
        uint32_t w[64];
        for(int i=0; i<16; ++i)
        w[i] = (uint32_t(data[4*i])<<24)|(uint32_t(data[4*i+1])<<16)|(uint32_t(data[4*i+2])<<8)|uint32_t(data[4*i+3]);
        for(int i=16; i<64; ++i)
        {
            uint32_t s0 = rotr(w[i-15],7)^rotr(w[i-15],18)^(w[i-15]>>3);
            uint32_t s1 = rotr(w[i-2],17)^rotr(w[i-2],19)^(w[i-2]>>10);
            w[i] = w[i-16]+s0+w[i-7]+s1;
        }
        uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],k=h[7];
        for(int i=0; i<64; ++i)
        {
            uint32_t t1 = k+(rotr(e,6)^rotr(e,11)^rotr(e,25))+((e&f)^(~e&g))+K[i]+w[i];
            uint32_t t2 = (rotr(a,2)^rotr(a,13)^rotr(a,22))+((a&b)^(a&c)^(b&c));
            k=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
        }
        h[0]+=a; h[1]+=b; h[2]+=c; h[3]+=d; h[4]+=e; h[5]+=f; h[6]+=g; h[7]+=k;
    }
}

std::string runlog::sha256_file(const std::string &path)
{
    std::ifstream in(path, std::ios::binary);
    if(!in.is_open())
    return "";

    std::vector<unsigned char> data((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    uint64_t bits = uint64_t(data.size())*8;
    data.push_back(0x80);
    while(data.size()%64!=56)
    data.push_back(0x00);
    for(int i=7; i>=0; --i)
    data.push_back(static_cast<unsigned char>(bits>>(8*i)));

    uint32_t h[8] = {0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    for(size_t i=0; i<data.size(); i+=64)
    block(h,&data[i]);

    std::ostringstream hex;
    for(int i=0; i<8; ++i)
    hex<<std::hex<<std::setw(8)<<std::setfill('0')<<h[i];
    return hex.str();
}
