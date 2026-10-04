// Standalone verification of the LAGOON store writer (lagoon_store): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../src lagoon_store_test.cpp ../../src/lagoon_store.cpp -lz -o lagoon_store_test
// Run:    ./lagoon_store_test     (writes and reads back lagoon_store_test.lagoon in the current folder)
// The files follow LAGOON's lagoon-format SPEC.md 0.1; LAGOON reads them with lagoon_format.store.
#include"lagoon_store.h"
#include<cmath>
#include<cstdint>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<vector>
#include<zlib.h>

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

static std::string slurp(const std::string &file)
{
    std::ifstream in(file.c_str(), std::ios::binary);
    std::stringstream s;
    s<<in.rdbuf();
    return s.str();
}

static uint64_t u64(const unsigned char *p)
{
    uint64_t v=0;
    for(int b=7;b>=0;--b) v=(v<<8)|p[b];
    return v;
}

static uint32_t crc32c(const unsigned char *d, size_t n)
{
    uint32_t c=0xFFFFFFFFu;
    for(size_t i=0;i<n;++i)
    {
        c^=d[i];
        for(int k=0;k<8;++k) c=(c&1) ? (c>>1)^0x82F63B78u : c>>1;
    }
    return c^0xFFFFFFFFu;
}

// one inner chunk of a shard: gunzip, unshuffle (4-byte elements); empty if not written
static std::vector<float> inner_chunk(const std::string &shard, size_t entries, size_t entry, size_t values, bool &crc_ok)
{
    const unsigned char *d = reinterpret_cast<const unsigned char*>(shard.data());
    const size_t index_bytes = 16*entries;
    const unsigned char *index = d + shard.size() - index_bytes - 4;
    uint32_t stored=0;
    for(int b=3;b>=0;--b) stored=(stored<<8)|index[index_bytes+b];
    crc_ok = stored==crc32c(index,index_bytes);
    const uint64_t offset=u64(index+16*entry), size=u64(index+16*entry+8);
    if(offset==UINT64_MAX && size==UINT64_MAX)
        return {};
    std::vector<unsigned char> raw(values*4);
    z_stream zs;
    std::memset(&zs,0,sizeof(zs));
    inflateInit2(&zs,15+32);
    zs.next_in=const_cast<unsigned char*>(d+offset);
    zs.avail_in=size;
    zs.next_out=raw.data();
    zs.avail_out=raw.size();
    const int status=inflate(&zs,Z_FINISH);
    inflateEnd(&zs);
    if(status!=Z_STREAM_END || zs.total_out!=raw.size())
        return {};
    std::vector<float> out(values);
    unsigned char *o=reinterpret_cast<unsigned char*>(out.data());
    for(size_t e=0;e<values;++e)
        for(int b=0;b<4;++b)
            o[4*e+b]=raw[b*values+e];
    return out;
}

int main()
{
    char s[300];

    std::cout<<"VTU header"<<std::endl;
    {
        const std::string header =
            "<VTKFile type=\"UnstructuredGrid\">\n<UnstructuredGrid>\n<Piece NumberOfPoints=\"60\" NumberOfCells=\"24\">\n"
            "<PointData>\n"
            "<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\"/>\n"
            "<DataArray type=\"Float32\" Name=\"pressure\" format=\"appended\" offset=\"724\"/>\n"
            "</PointData>\n<Points>\n<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"968\"/>\n</Points>\n"
            "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"1692\"/>\n</Cells>\n"
            "</Piece>\n</UnstructuredGrid>\n<AppendedData encoding=\"raw\">\n_";
        std::vector<lagoon_store::vtu_array> fields;
        long long points=-1;
        const bool ok = lagoon_store::parse_vtu_header(header,fields,points);
        check(ok && fields.size()==2, "two point arrays found");
        check(ok && fields[0].var.name=="velocity" && fields[0].var.components==3 && fields[0].offset==0
                 && fields[0].var.units=="m/s", "velocity: 3 components at offset 0, in m/s");
        check(ok && fields[1].var.name=="pressure" && fields[1].var.components==1 && fields[1].offset==724, "pressure: scalar at offset 724");
        check(points==968, "points at offset 968 (not the cells)");
    }

    std::cout<<"level offsets"<<std::endl;
    {
        const std::vector<double> sigma={0.0,0.3,0.7,1.0};
        std::vector<float> z(4*2), off;
        for(int k=0;k<4;++k) { z[k*2]=float(-2.0+sigma[k]*2.1); z[k*2+1]=float(-1.0+sigma[k]*1.1); }
        lagoon_store::level_offsets(z.data(),4,2,sigma,off);
        bool zero=true; for(float v:off) zero = zero && v==0.0f;
        check(zero, "levels at their sigma: all offsets exactly 0");
        z[1*2+1] = float(-1.0+0.25*1.1);
        lagoon_store::level_offsets(z.data(),4,2,sigma,off);
        snprintf(s,300,"a level off its sigma: offset %.6f (exact %.6f), the others 0",off[1*2+1],-0.05*1.1);
        check(std::fabs(off[1*2+1]+0.055)<1e-6 && off[1*2]==0.0f && off[2*2+1]==0.0f, s);
    }

    std::cout<<"store files"<<std::endl;
    {
        const std::string path="lagoon_store_test.lagoon";
        if(std::system(("rm -rf "+path).c_str())!=0) std::cout<<"  (could not clear "<<path<<")"<<std::endl;
        const int NX=5, NY=3, NZ=4, T=2;
        std::vector<double> x={0,1,2,3,4}, y={0,1,2}, sigma={0.0,0.3,0.7,1.0};
        std::vector<lagoon_store::variable> vars={{"velocity",3,"m/s"},{"pressure",1,"Pa"}};
        lagoon_store store(path,T,1);
        store.create_root("NHFLOW","{\"type\": \"run\", \"run\": \"R1\"}");
        store.create_output("volume","sigma",x,y,sigma,vars,1,"test");
        store.create_block("volume",0,0,0,NX,NY,NZ,vars,0);
        std::vector<std::vector<float> > pressure(3);
        long long size_after_first=0;
        for(int t=0;t<3;++t)
        {
            std::vector<float> bed(NX*NY,-2.0f), top(NX*NY), u(3*NX*NY*NZ), p(NX*NY*NZ), off(NX*NY*NZ,0.0f);
            for(int q=0;q<NX*NY;++q) top[q]=float(0.1*t+0.01*q);
            for(int n=0;n<NX*NY*NZ;++n) { p[n]=float(1000*t+n); u[3*n]=float(t); u[3*n+1]=float(n); u[3*n+2]=-float(n); }
            pressure[t]=p;
            store.write("volume",0,t,"z_bed",bed.data());
            store.write("volume",0,t,"z_surface",top.data());
            store.write("volume",0,t,"z_offset",off.data());
            store.write("volume",0,t,"velocity",u.data());
            store.write("volume",0,t,"pressure",p.data());
            store.commit("volume",t,0.5*t,10*t);
            if(t==0) size_after_first=(long long)slurp(path+"/volume/blocks/b0000/pressure/c/0/0/0/0").size();
        }
        const std::string root=slurp(path+"/zarr.json");
        check(root.find("\"format\": \"lagoon\"")!=std::string::npos && root.find("\"outputs\": [\"volume\"]")!=std::string::npos
              && root.find("\"run\": \"R1\"")!=std::string::npos, "root group: format, outputs, run");
        const std::string meta=slurp(path+"/volume/blocks/b0000/pressure/zarr.json");
        check(meta.find("\"shape\": [3, 4, 3, 5]")!=std::string::npos, "pressure: shape (3 outputs, 4 levels, 3, 5)");
        check(meta.find("\"chunk_shape\": [2, 4, 3, 5]")!=std::string::npos && meta.find("\"chunk_shape\": [1, 4, 3, 5]")!=std::string::npos,
              "shards of 2 outputs, inner chunks of one output (small block: all levels)");
        check(slurp(path+"/volume/time/zarr.json").find("\"shape\": [3]")!=std::string::npos, "time: 3 outputs committed");
        const std::string time=slurp(path+"/volume/time/c/0");
        double t2; std::memcpy(&t2,time.data()+16,8);
        check(time.size()==4096*8 && t2==1.0, "time of output 2 is 1.0 s");
        check(!std::ifstream((path+"/volume/blocks/b0000/z_offset/zarr.json").c_str()), "z_offset all zero: not stored");

        const std::string shard0=slurp(path+"/volume/blocks/b0000/pressure/c/0/0/0/0");
        const std::string shard1=slurp(path+"/volume/blocks/b0000/pressure/c/1/0/0/0");
        check((long long)shard0.size()>size_after_first, "shard 0 grew when output 1 was appended");
        bool crc0=false, crc1=false;
        const std::vector<float> p1=inner_chunk(shard0,2,1,NX*NY*NZ,crc0);
        const std::vector<float> p0=inner_chunk(shard0,2,0,NX*NY*NZ,crc0);
        const std::vector<float> p2=inner_chunk(shard1,2,0,NX*NY*NZ,crc1);
        const std::vector<float> none=inner_chunk(shard1,2,1,NX*NY*NZ,crc1);
        check(crc0 && crc1, "shard indexes carry a valid CRC-32C");
        check(p0==pressure[0] && p1==pressure[1], "outputs 0 and 1 read back from shard 0 (after the append)");
        check(p2==pressure[2], "output 2 reads back from shard 1");
        check(none.empty(), "output 3 is not written: empty index entry");
        bool v_ok=false;
        const std::string vshard=slurp(path+"/volume/blocks/b0000/velocity/c/0/0/0/0/0");
        const std::vector<float> v1=inner_chunk(vshard,2,1,3*NX*NY*NZ,v_ok);
        check(v_ok && v1.size()==size_t(3*NX*NY*NZ) && v1[3*7]==1.0f && v1[3*7+1]==7.0f && v1[3*7+2]==-7.0f,
              "velocity: the components of a point together");
    }

    std::cout<<(nfail ? "FAILED: " : "all tests passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
