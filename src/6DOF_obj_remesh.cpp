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

#include"6DOF_obj_remesh.h"
#include<algorithm>
#include<cmath>
#include<cstdint>
#include<unordered_map>
#include<iomanip>

namespace
{
using vec3 = sixdof_remesh::vec3;

constexpr double PI_ = 3.14159265358979323846;

inline vec3 operator+(const vec3 &a, const vec3 &b) {return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};}
inline vec3 operator-(const vec3 &a, const vec3 &b) {return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};}
inline vec3 operator*(double s, const vec3 &a) {return {s*a[0], s*a[1], s*a[2]};}
inline double dot(const vec3 &a, const vec3 &b) {return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];}
inline vec3 cross(const vec3 &a, const vec3 &b) {return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};}
inline double norm(const vec3 &a) {return std::sqrt(dot(a,a));}
inline double dist2(const vec3 &a, const vec3 &b) {vec3 d=a-b; return dot(d,d);}
inline vec3 unit(const vec3 &a) {double n=norm(a); return n>0.0 ? (1.0/n)*a : vec3{0.0,0.0,0.0};}

inline std::uint64_t ekey(int a, int b)
{
    if(a>b) std::swap(a,b);
    return (std::uint64_t(std::uint32_t(a))<<32) | std::uint64_t(std::uint32_t(b));
}

// ------------------------------------------------------------------ closest points
vec3 closest_on_segment(const vec3 &p, const vec3 &a, const vec3 &b)
{
    vec3 ab = b-a;
    double l2 = dot(ab,ab);
    if(l2<=0.0) return a;
    double t = dot(p-a,ab)/l2;
    t = t<0.0 ? 0.0 : (t>1.0 ? 1.0 : t);
    return a + t*ab;
}

// Ericson, Real-Time Collision Detection, 5.1.5
vec3 closest_on_triangle(const vec3 &p, const vec3 &a, const vec3 &b, const vec3 &c)
{
    vec3 ab=b-a, ac=c-a, ap=p-a;
    double d1=dot(ab,ap), d2=dot(ac,ap);
    if(d1<=0.0 && d2<=0.0) return a;
    vec3 bp=p-b;
    double d3=dot(ab,bp), d4=dot(ac,bp);
    if(d3>=0.0 && d4<=d3) return b;
    double vc=d1*d4-d3*d2;
    if(vc<=0.0 && d1>=0.0 && d3<=0.0) {double v=d1/(d1-d3); return a + v*ab;}
    vec3 cp=p-c;
    double d5=dot(ab,cp), d6=dot(ac,cp);
    if(d6>=0.0 && d5<=d6) return c;
    double vb=d5*d2-d1*d6;
    if(vb<=0.0 && d2>=0.0 && d6<=0.0) {double w=d2/(d2-d6); return a + w*ac;}
    double va=d3*d6-d5*d4;
    if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0) {double w=(d4-d3)/((d4-d3)+(d5-d6)); return b + w*(c-b);}
    double den=va+vb+vc;
    if(den==0.0) return closest_on_segment(p,a,b);
    double v=vb/den, w=vc/den;
    return a + v*ab + w*ac;
}

// ------------------------------------------------------------------ BVH for closest point queries
struct bvh
{
    struct node {double lo[3], hi[3]; int left=-1, right=-1, start=0, count=0;};
    std::vector<node> nodes;
    std::vector<int> idx;
    std::vector<std::array<double,6> > box;

    void build(const std::vector<std::array<double,6> > &b)
    {
        box = b;
        int n = int(b.size());
        idx.resize(n);
        for(int i=0; i<n; ++i) idx[i]=i;
        nodes.clear();
        if(n==0) return;
        nodes.reserve(2*n/2+8);
        build_rec(0,n);
    }

    int build_rec(int s, int e)
    {
        node nd;
        for(int d=0; d<3; ++d) {nd.lo[d]=1.0e300; nd.hi[d]=-1.0e300;}
        double clo[3]={1.0e300,1.0e300,1.0e300}, chi[3]={-1.0e300,-1.0e300,-1.0e300};
        for(int i=s; i<e; ++i)
        {
            const auto &bb = box[idx[i]];
            for(int d=0; d<3; ++d)
            {
                nd.lo[d]=std::min(nd.lo[d],bb[d]);
                nd.hi[d]=std::max(nd.hi[d],bb[3+d]);
                double c=0.5*(bb[d]+bb[3+d]);
                clo[d]=std::min(clo[d],c);
                chi[d]=std::max(chi[d],c);
            }
        }
        int id = int(nodes.size());
        nodes.push_back(nd);
        if(e-s<=4)
        {
            nodes[id].start=s;
            nodes[id].count=e-s;
            return id;
        }
        int ax=0;
        for(int d=1; d<3; ++d)
        if(chi[d]-clo[d] > chi[ax]-clo[ax]) ax=d;
        int m=(s+e)/2;
        std::nth_element(idx.begin()+s, idx.begin()+m, idx.begin()+e, [&](int a, int b)
        {
            return box[a][ax]+box[a][3+ax] < box[b][ax]+box[b][3+ax];
        });
        int l = build_rec(s,m);
        int r = build_rec(m,e);
        nodes[id].left=l;
        nodes[id].right=r;
        return id;
    }

    static double boxdist2(const node &nd, const vec3 &q)
    {
        double d2=0.0;
        for(int d=0; d<3; ++d)
        {
            double v=0.0;
            if(q[d]<nd.lo[d]) v=nd.lo[d]-q[d];
            else if(q[d]>nd.hi[d]) v=q[d]-nd.hi[d];
            d2+=v*v;
        }
        return d2;
    }

    // f(prim) returns squared distance (or a value >= best to reject); best is updated
    template<class F>
    int closest(const vec3 &q, F &&f, double &best) const
    {
        int bestid=-1;
        if(nodes.empty()) return -1;
        int stack[128];
        int sp=0;
        stack[sp++]=0;
        while(sp>0)
        {
            const node &nd = nodes[stack[--sp]];
            if(boxdist2(nd,q)>=best) continue;
            if(nd.left<0)
            {
                for(int i=nd.start; i<nd.start+nd.count; ++i)
                {
                    double d2=f(idx[i]);
                    if(d2<best) {best=d2; bestid=idx[i];}
                }
                continue;
            }
            double dl=boxdist2(nodes[nd.left],q), dr=boxdist2(nodes[nd.right],q);
            if(sp+2>128) continue;
            if(dl<dr) {stack[sp++]=nd.right; stack[sp++]=nd.left;}
            else      {stack[sp++]=nd.left;  stack[sp++]=nd.right;}
        }
        return bestid;
    }
};

// ------------------------------------------------------------------ indexed triangle mesh
struct mesh
{
    std::vector<vec3> P;
    std::vector<std::array<int,3> > F;
    std::vector<char> fal, val;
    std::vector<std::vector<int> > vf;
    std::vector<int> lev;        // 0 smooth, 1 crease (on a feature curve), 2 corner (fixed)
    std::vector<int> chain;      // feature curve of a crease vertex
    std::vector<double> vtarget; // ideal valence of corner vertices
    std::unordered_map<std::uint64_t,int> fe;   // feature edge -> feature curve id
    std::vector<char> chain_bnd; // feature curve is an open boundary
    std::vector<char> fpatch;    // face belongs to a hole filling patch
    std::vector<char> vpatch;    // vertex is a free (inner) vertex of a hole filling patch

    int nfaces_alive=0;

    int add_vertex(const vec3 &p, int l, int ch)
    {
        P.push_back(p);
        val.push_back(1);
        vf.emplace_back();
        lev.push_back(l);
        chain.push_back(ch);
        vtarget.push_back(6.0);
        vpatch.push_back(0);
        return int(P.size())-1;
    }

    int add_face(int a, int b, int c)
    {
        F.push_back({a,b,c});
        fal.push_back(1);
        fpatch.push_back(0);
        int f=int(F.size())-1;
        vf[a].push_back(f);
        vf[b].push_back(f);
        vf[c].push_back(f);
        ++nfaces_alive;
        return f;
    }

    static void erase_val(std::vector<int> &v, int x)
    {
        auto it=std::find(v.begin(),v.end(),x);
        if(it!=v.end()) {*it=v.back(); v.pop_back();}
    }

    void kill_face(int f)
    {
        for(int q=0; q<3; ++q) erase_val(vf[F[f][q]],f);
        fal[f]=0;
        --nfaces_alive;
    }

    static bool has(const std::array<int,3> &t, int v) {return t[0]==v || t[1]==v || t[2]==v;}

    void edge_faces(int a, int b, std::vector<int> &out) const
    {
        out.clear();
        for(int f: vf[a]) if(has(F[f],b)) out.push_back(f);
    }

    void neighbors(int v, std::vector<int> &out) const
    {
        out.clear();
        for(int f: vf[v])
        for(int q=0; q<3; ++q)
        if(F[f][q]!=v) out.push_back(F[f][q]);
        std::sort(out.begin(),out.end());
        out.erase(std::unique(out.begin(),out.end()),out.end());
    }

    int valence(int v) const
    {
        // a smooth vertex has no boundary or non-manifold edge: closed fan, valence = #faces
        if(lev[v]==0) return int(vf[v].size());
        static thread_local std::vector<int> nb;
        neighbors(v,nb);
        return int(nb.size());
    }

    bool is_feature(int a, int b) const {return fe.count(ekey(a,b))>0;}

    vec3 fnormal(int f) const   // area weighted (|n| = 2A)
    {
        const auto &t=F[f];
        return cross(P[t[1]]-P[t[0]], P[t[2]]-P[t[0]]);
    }

    // third vertex of face f opposite to edge (a,b); dir=+1 if a->b is a directed edge of f
    int opposite(int f, int a, int b, int &dir) const
    {
        const auto &t=F[f];
        for(int q=0; q<3; ++q)
        {
            int u=t[q], w=t[(q+1)%3];
            if(u==a && w==b) {dir=+1; return t[(q+2)%3];}
            if(u==b && w==a) {dir=-1; return t[(q+2)%3];}
        }
        dir=0;
        return -1;
    }

    // unique alive edges (a<b)
    void edges(std::vector<std::pair<int,int> > &E) const
    {
        E.clear();
        E.reserve(3*F.size()/2+16);
        for(size_t f=0; f<F.size(); ++f)
        if(fal[f])
        for(int q=0; q<3; ++q)
        {
            int a=F[f][q], b=F[f][(q+1)%3];
            E.emplace_back(std::min(a,b),std::max(a,b));
        }
        std::sort(E.begin(),E.end());
        E.erase(std::unique(E.begin(),E.end()),E.end());
    }
};

double tri_quality(const vec3 &a, const vec3 &b, const vec3 &c)
{
    double A=0.5*norm(cross(b-a,c-a));
    double s=dist2(a,b)+dist2(b,c)+dist2(c,a);
    return s>0.0 ? 4.0*std::sqrt(3.0)*A/s : 0.0;
}

double angle_at(const vec3 &p, const vec3 &a, const vec3 &b)
{
    vec3 u=a-p, v=b-p;
    double nu=norm(u), nv=norm(v);
    if(nu<=0.0 || nv<=0.0) return 0.0;
    double c=dot(u,v)/(nu*nv);
    c = c>1.0 ? 1.0 : (c<-1.0 ? -1.0 : c);
    return std::acos(c);
}

double tri_minangle(const vec3 &a, const vec3 &b, const vec3 &c)
{
    return std::min(angle_at(a,b,c), std::min(angle_at(b,c,a), angle_at(c,a,b)));
}

void soup_stats(const std::vector<vec3> &T, double &area, double &vol, double &amin, double &qmean, double &qmin, double &fq05)
{
    area=vol=0.0;
    amin=PI_;
    qmean=0.0;
    qmin=1.0;
    fq05=0.0;
    size_t n=T.size()/3;
    for(size_t t=0; t<n; ++t)
    {
        const vec3 &a=T[3*t], &b=T[3*t+1], &c=T[3*t+2];
        area+=0.5*norm(cross(b-a,c-a));
        vol+=dot(a,cross(b,c))/6.0;
        amin=std::min(amin,tri_minangle(a,b,c));
        double q=tri_quality(a,b,c);
        qmean+=q;
        qmin=std::min(qmin,q);
        if(q<0.5) fq05+=1.0;
    }
    if(n>0) {qmean/=double(n); fq05/=double(n);}
    amin*=180.0/PI_;
}

// ------------------------------------------------------------------ remesher
class remesher
{
public:
    remesher(const sixdof_remesh::metric_func &hf, const sixdof_remesh::params &pr, sixdof_remesh::stats &s)
    : H(hf), prm(pr), st(s) {}

    bool run(const std::vector<vec3> &in, std::vector<vec3> &out);

private:
    const sixdof_remesh::metric_func &H;
    const sixdof_remesh::params &prm;
    sixdof_remesh::stats &st;

    mesh M;
    double diag=1.0, cos_feat=0.0, tiny_area=0.0;

    // reference surface
    std::vector<vec3> rP;
    std::vector<std::array<int,3> > rF;
    std::vector<vec3> rN;
    bvh rbvh;
    std::vector<std::array<int,3> > rS;   // feature segments (a,b,chain)
    bvh sbvh;

    std::vector<int> tmp1, tmp2, tmp3;

    // target lengths along the axes at x
    vec3 hv(const vec3 &x) const
    {
        vec3 s=H(x[0],x[1],x[2]);
        for(int d=0; d<3; ++d) s[d] = s[d]>1.0e-12*diag ? s[d] : 1.0e-12*diag;
        return s;
    }
    // isotropic equivalent target length
    double hiso(const vec3 &x) const
    {
        vec3 s=hv(x);
        return std::cbrt(s[0]*s[1]*s[2]);
    }
    static vec3 scale(const vec3 &d, const vec3 &s) {return {d[0]/s[0], d[1]/s[1], d[2]/s[2]};}
    // squared metric length of the edge a-b (metric at the midpoint), 1 = target
    double ml2(const vec3 &a, const vec3 &b) const
    {
        vec3 d=scale(b-a,hv(0.5*(a+b)));
        return dot(d,d);
    }
    // area of triangle a,b,c in the metric (at its centroid), target triangle: sqrt(3)/4
    double marea(const vec3 &a, const vec3 &b, const vec3 &c) const
    {
        vec3 s=hv((1.0/3.0)*(a+b+c));
        return 0.5*norm(cross(scale(b-a,s),scale(c-a,s)));
    }
    // min angle of the quad split (p,q,r),(q,p,s) evaluated in the metric
    void metric_quad(const vec3 &pa, const vec3 &pb, const vec3 &pc, const vec3 &pd,
                     vec3 &qa, vec3 &qb, vec3 &qc, vec3 &qd) const
    {
        vec3 s=hv(0.25*(pa+pb+pc+pd));
        qa=scale(pa,s); qb=scale(pb,s); qc=scale(pc,s); qd=scale(pd,s);
    }

    void weld(const std::vector<vec3> &in);
    int  repair_tjunctions(double tol, bool guard);
    void boundary_halfedges(std::vector<std::pair<int,int> > &H);
    void merge_vertices(int v, int u);
    int  close_gaps();
    void fill_holes();
    void fill_loop(const std::vector<int> &L);
    void refine_and_fair_patches();
    void fair_patch(int mode);
    void detect_features();
    void build_reference();
    long estimate_triangles() const;

    vec3 project_surface(const vec3 &q, const vec3 &n) const;
    vec3 project_chain(const vec3 &q, int ch) const;

    double target_valence(int v) const
    {
        if(M.lev[v]==2) return M.vtarget[v];
        if(M.lev[v]==1 && M.chain[v]>=0 && M.chain_bnd[M.chain[v]]) return 4.0;
        return 6.0;
    }

    void split(int a, int b);
    bool collapse(int a, int b);   // removes a
    bool flip_ok(int a, int b, int &f1, int &f2, int &c, int &d);
    void do_flip(int a, int b, int f1, int f2, int c, int d);
    vec3 vnormal(int v) const;
    bool faces_ok_after_move(int v, const vec3 &np) const;

    int split_long();
    int collapse_short();
    int flip_valence();
    int flip_delaunay(bool patch_only=false);
    void relax();
};

void remesher::weld(const std::vector<vec3> &in)
{
    vec3 lo{1.0e300,1.0e300,1.0e300}, hi{-1.0e300,-1.0e300,-1.0e300};
    for(const auto &p: in)
    for(int d=0; d<3; ++d) {lo[d]=std::min(lo[d],p[d]); hi[d]=std::max(hi[d],p[d]);}
    diag = norm(hi-lo);
    if(diag<=0.0) diag=1.0;
    const double tol = prm.merge_tol*diag;
    tiny_area = 1.0e-14*diag*diag;

    struct key3 {long long x,y,z; bool operator==(const key3 &o) const {return x==o.x && y==o.y && z==o.z;}};
    struct hash3 {size_t operator()(const key3 &k) const {return size_t(k.x*73856093LL ^ k.y*19349663LL ^ k.z*83492791LL);}};
    std::unordered_map<key3,std::vector<int>,hash3> grid;
    grid.reserve(in.size());

    std::vector<int> id(in.size());
    for(size_t i=0; i<in.size(); ++i)
    {
        const vec3 &p=in[i];
        long long kx=(long long)std::floor((p[0]-lo[0])/tol);
        long long ky=(long long)std::floor((p[1]-lo[1])/tol);
        long long kz=(long long)std::floor((p[2]-lo[2])/tol);
        int found=-1;
        for(int dx=-1; dx<=1 && found<0; ++dx)
        for(int dy=-1; dy<=1 && found<0; ++dy)
        for(int dz=-1; dz<=1 && found<0; ++dz)
        {
            auto it=grid.find(key3{kx+dx,ky+dy,kz+dz});
            if(it==grid.end()) continue;
            for(int v: it->second)
            if(dist2(M.P[v],p)<=tol*tol) {found=v; break;}
        }
        if(found<0)
        {
            found=M.add_vertex(p,0,-1);
            grid[key3{kx,ky,kz}].push_back(found);
        }
        id[i]=found;
    }

    // faces: drop collapsed and duplicate triangles
    std::unordered_map<std::uint64_t,std::vector<int> > seen;
    for(size_t t=0; t<in.size()/3; ++t)
    {
        int a=id[3*t], b=id[3*t+1], c=id[3*t+2];
        if(a==b || b==c || c==a) continue;
        std::array<int,3> s{a,b,c};
        std::sort(s.begin(),s.end());
        std::uint64_t k=ekey(s[0],s[1])*1000003ULL ^ std::uint64_t(s[2]);
        bool dup=false;
        for(int f: seen[k])
        {
            std::array<int,3> g=M.F[f];
            std::sort(g.begin(),g.end());
            if(g==s) {dup=true; break;}
        }
        if(dup) continue;
        seen[k].push_back(M.add_face(a,b,c));
    }
}

// Split faces whose boundary edge passes through a boundary vertex of a neighbouring face
// (non-conforming STL). Returns number of splits.
int remesher::repair_tjunctions(double tol, bool guard)
{
    int nfix=0;

    for(int pass=0; pass<20; ++pass)
    {
        std::vector<std::pair<int,int> > E;
        M.edges(E);
        std::vector<std::array<int,3> > bnd;   // a, b, face  with a->b directed in face
        std::vector<int> bv;
        for(auto &e: E)
        {
            M.edge_faces(e.first,e.second,tmp1);
            if(tmp1.size()!=1) continue;
            int dir;
            M.opposite(tmp1[0],e.first,e.second,dir);
            if(dir>0) bnd.push_back({e.first,e.second,tmp1[0]});
            else      bnd.push_back({e.second,e.first,tmp1[0]});
            bv.push_back(e.first);
            bv.push_back(e.second);
        }
        if(bnd.empty()) break;
        std::sort(bv.begin(),bv.end());
        bv.erase(std::unique(bv.begin(),bv.end()),bv.end());

        // BVH over boundary vertices
        std::vector<std::array<double,6> > bb(bv.size());
        for(size_t i=0; i<bv.size(); ++i)
        {
            const vec3 &p=M.P[bv[i]];
            bb[i]={p[0],p[1],p[2],p[0],p[1],p[2]};
        }
        bvh vb;
        vb.build(bb);

        int nsplit=0;
        std::vector<char> touched(M.F.size(),0);
        for(auto &e: bnd)
        {
            int a=e[0], b=e[1], f=e[2];
            if(touched[f] || !M.fal[f]) continue;
            const vec3 &pa=M.P[a], &pb=M.P[b];
            vec3 ab=pb-pa;
            double l2=dot(ab,ab);
            if(l2<=0.0) continue;
            int dir0;
            const int cself=M.opposite(f,a,b,dir0);
            // gap closing: the snapping distance is also limited by the edge length
            const double etol = guard ? std::min(tol,0.3*std::sqrt(l2)) : tol;
            // candidate vertex closest to the open edge: scan vertices within tol of the segment
            int best=-1;
            double bestt=2.0;
            // walk the BVH by brute-force on boxes overlapping the segment's bbox
            std::vector<int> stack{0};
            double lo[3], hi[3];
            for(int d=0; d<3; ++d) {lo[d]=std::min(pa[d],pb[d])-tol; hi[d]=std::max(pa[d],pb[d])+tol;}
            while(!stack.empty() && !vb.nodes.empty())
            {
                const auto &nd=vb.nodes[stack.back()];
                stack.pop_back();
                bool ov=true;
                for(int d=0; d<3; ++d) if(nd.hi[d]<lo[d] || nd.lo[d]>hi[d]) ov=false;
                if(!ov) continue;
                if(nd.left>=0) {stack.push_back(nd.left); stack.push_back(nd.right); continue;}
                for(int i=nd.start; i<nd.start+nd.count; ++i)
                {
                    int v=bv[vb.idx[i]];
                    if(v==a || v==b || v==cself) continue;
                    double t=dot(M.P[v]-pa,ab)/l2;
                    if(t<=1.0e-6 || t>=1.0-1.0e-6) continue;
                    if(guard && (t*t*l2<=etol*etol || (1.0-t)*(1.0-t)*l2<=etol*etol)) continue;
                    if(dist2(pa+t*ab,M.P[v])>etol*etol) continue;
                    if(t<bestt) {bestt=t; best=v;}
                }
            }
            if(best<0) continue;
            int dir;
            int c=M.opposite(f,a,b,dir);
            M.kill_face(f);
            int g1=M.add_face(a,best,c);
            int g2=M.add_face(best,b,c);
            M.fpatch[g1]=M.fpatch[g2]=M.fpatch[f];
            touched.resize(M.F.size(),0);
            touched[g1]=touched[g2]=1;
            ++nsplit;
        }
        nfix+=nsplit;
        if(nsplit==0) break;
    }
    return nfix;
}

// ------------------------------------------------------------------ open boundaries: gaps and holes
void remesher::boundary_halfedges(std::vector<std::pair<int,int> > &H)
{
    H.clear();
    for(size_t f=0; f<M.F.size(); ++f)
    if(M.fal[f])
    for(int q=0; q<3; ++q)
    {
        int a=M.F[f][q], b=M.F[f][(q+1)%3];
        M.edge_faces(a,b,tmp3);
        if(tmp3.size()==1) H.emplace_back(a,b);
    }
    std::sort(H.begin(),H.end());
}

// merge vertex v into u (at the midpoint), remove collapsed and duplicate faces
void remesher::merge_vertices(int v, int u)
{
    M.P[u]=0.5*(M.P[u]+M.P[v]);
    std::vector<int> fv=M.vf[v];
    for(int f: fv)
    {
        auto &t=M.F[f];
        for(int q=0; q<3; ++q) if(t[q]==v) t[q]=u;
        M.vf[u].push_back(f);
    }
    M.vf[v].clear();
    M.val[v]=0;

    std::vector<int> fu=M.vf[u];
    std::sort(fu.begin(),fu.end());
    fu.erase(std::unique(fu.begin(),fu.end()),fu.end());
    for(int f: fu)
    {
        if(!M.fal[f]) continue;
        const auto &t=M.F[f];
        if(t[0]==t[1] || t[1]==t[2] || t[2]==t[0]) M.kill_face(f);
    }
    // duplicate faces: same orientation -> keep one, opposite orientation -> zero volume fin, remove both
    fu=M.vf[u];
    for(size_t i=0; i<fu.size(); ++i)
    for(size_t j=i+1; j<fu.size(); ++j)
    {
        int f=fu[i], g=fu[j];
        if(!M.fal[f] || !M.fal[g]) continue;
        std::array<int,3> s=M.F[f], r=M.F[g];
        std::sort(s.begin(),s.end());
        std::sort(r.begin(),r.end());
        if(s!=r) continue;
        int d1,d2;
        M.opposite(f,M.F[f][0],M.F[f][1],d1);
        M.opposite(g,M.F[f][0],M.F[f][1],d2);
        M.kill_face(g);
        if(d1!=d2) M.kill_face(f);
    }
}

// Close cracks: T-junctions and boundary vertices within the gap tolerance
int remesher::close_gaps()
{
    const double gtol=prm.gap_tol*diag;
    int nmerge=0;
    std::vector<std::pair<int,int> > H;

    for(int pass=0; pass<10; ++pass)
    {
        boundary_halfedges(H);
        if(H.empty()) break;

        std::vector<double> Lb(M.P.size(),1.0e300);
        std::vector<int> bv;
        for(auto &e: H)
        {
            double L=std::sqrt(dist2(M.P[e.first],M.P[e.second]));
            Lb[e.first]=std::min(Lb[e.first],L);
            Lb[e.second]=std::min(Lb[e.second],L);
            bv.push_back(e.first);
            bv.push_back(e.second);
        }
        std::sort(bv.begin(),bv.end());
        bv.erase(std::unique(bv.begin(),bv.end()),bv.end());

        std::vector<std::array<double,6> > bb(bv.size());
        for(size_t i=0; i<bv.size(); ++i)
        {
            const vec3 &p=M.P[bv[i]];
            bb[i]={p[0],p[1],p[2],p[0],p[1],p[2]};
        }
        bvh vb;
        vb.build(bb);

        std::vector<char> done(M.P.size(),0);
        int nm=0;
        for(int v: bv)
        {
            if(!M.val[v] || done[v]) continue;
            double r=std::min(gtol,0.3*Lb[v]);
            double best=r*r;
            int bi=vb.closest(M.P[v],[&](int i)
            {
                int u=bv[i];
                if(u==v || !M.val[u] || done[u]) return 1.0e300;
                double d2=dist2(M.P[u],M.P[v]);
                double ru=0.3*Lb[u];
                if(d2>ru*ru) return 1.0e300;
                return d2;
            },best);
            if(bi<0) continue;
            int u=bv[bi];
            M.neighbors(v,tmp2);
            if(std::binary_search(tmp2.begin(),tmp2.end(),u)) continue;
            merge_vertices(v,u);
            done[v]=done[u]=1;
            ++nm;
        }
        nmerge+=nm;

        int nt=repair_tjunctions(gtol,true);
        st.n_tjunctions+=nt;

        if(nm==0 && nt==0) break;
    }
    return nmerge;
}

void remesher::fill_holes()
{
    std::vector<std::pair<int,int> > H;
    boundary_halfedges(H);
    if(H.empty()) return;

    std::unordered_map<int,std::vector<int> > out;
    for(auto &e: H) out[e.first].push_back(e.second);
    std::unordered_map<std::uint64_t,char> used;
    auto dkey=[](int a, int b){return (std::uint64_t(std::uint32_t(a))<<32) | std::uint64_t(std::uint32_t(b));};

    for(auto &e0: H)
    {
        if(used.count(dkey(e0.first,e0.second))) continue;
        used[dkey(e0.first,e0.second)]=1;

        const int a0=e0.first;
        std::vector<int> L{a0};
        std::unordered_map<int,int> where;
        where[a0]=0;
        int cur=e0.second;
        bool fail=false;
        size_t steps=0;

        while(cur!=a0)
        {
            if(++steps>H.size()+1) {fail=true; break;}
            auto it=where.find(cur);
            if(it!=where.end())
            {
                // the walk passed a pinch vertex: emit the closed sub-loop
                int j=it->second;
                std::vector<int> sub(L.begin()+j,L.end());
                fill_loop(sub);
                for(size_t k=j+1; k<L.size(); ++k) where.erase(L[k]);
                L.resize(j+1);
            }
            else
            {
                where[cur]=int(L.size());
                L.push_back(cur);
            }
            int nb=-1;
            for(int b: out[cur]) if(!used.count(dkey(cur,b))) {nb=b; break;}
            if(nb<0) {fail=true; break;}
            used[dkey(cur,nb)]=1;
            cur=nb;
        }
        if(fail) {++st.n_holes_open; continue;}
        fill_loop(L);
    }
}

// Triangulate one boundary loop. L follows the boundary half-edges L[i]->L[i+1] of the mesh,
// so the patch faces use the reversed orientation.
void remesher::fill_loop(const std::vector<int> &L)
{
    const int n=int(L.size());
    if(n<3) {++st.n_holes_open; return;}
    ++st.n_holes_filled;
    st.n_hole_edges+=n;

    if(n==3)
    {
        int f=M.add_face(L[2],L[1],L[0]);
        M.fpatch[f]=1;
        return;
    }

    // Newell normal of the loop; patch normals point along -N
    vec3 N{0.0,0.0,0.0};
    for(int i=0; i<n; ++i)
    {
        const vec3 &p=M.P[L[i]], &q=M.P[L[(i+1)%n]];
        N[0]+=(p[1]-q[1])*(p[2]+q[2]);
        N[1]+=(p[2]-q[2])*(p[0]+q[0]);
        N[2]+=(p[0]-q[0])*(p[1]+q[1]);
    }
    const vec3 Nt=-1.0*unit(N);

    if(n<=400)
    {
        // minimum weight triangulation (dynamic programming over the polygon)
        // weight: area, penalised for patch normals turning away from the loop normal,
        // plus the squared edge lengths (favours short diagonals and compact triangles)
        std::vector<double> W(size_t(n)*n,0.0);
        std::vector<int> Lam(size_t(n)*n,-1);
        auto tw=[&](int i, int m, int k)
        {
            const vec3 &pi=M.P[L[i]], &pm=M.P[L[m]], &pk=M.P[L[k]];
            vec3 c=cross(pm-pk,pi-pk);
            double A2=norm(c);
            double cs = A2>0.0 ? dot(c,Nt)/A2 : 0.0;
            return 0.5*A2*(2.0-cs) + 0.1*(dist2(pi,pm)+dist2(pm,pk)+dist2(pk,pi));
        };
        for(int len=2; len<n; ++len)
        for(int i=0; i+len<n; ++i)
        {
            int k=i+len;
            double best=1.0e300;
            int bm=-1;
            for(int m=i+1; m<k; ++m)
            {
                double w=W[size_t(i)*n+m]+W[size_t(m)*n+k]+tw(i,m,k);
                if(w<best) {best=w; bm=m;}
            }
            W[size_t(i)*n+k]=best;
            Lam[size_t(i)*n+k]=bm;
        }
        std::vector<std::pair<int,int> > stack{{0,n-1}};
        while(!stack.empty())
        {
            auto [i,k]=stack.back();
            stack.pop_back();
            if(k-i<2) continue;
            int m=Lam[size_t(i)*n+k];
            int f=M.add_face(L[k],L[m],L[i]);
            M.fpatch[f]=1;
            stack.push_back({i,m});
            stack.push_back({m,k});
        }
    }
    else
    {
        // very long loop: fan around the centroid, the fairing moves the centre
        vec3 c{0.0,0.0,0.0};
        for(int v: L) c=c+M.P[v];
        c=(1.0/n)*c;
        int cv=M.add_vertex(c,0,-1);
        M.vpatch[cv]=1;
        for(int i=0; i<n; ++i)
        {
            int f=M.add_face(L[(i+1)%n],L[i],cv);
            M.fpatch[f]=1;
        }
    }
}

// Refine the patches towards the target size and fair them with the rim fixed
void remesher::refine_and_fair_patches()
{
    bool any=false;
    for(size_t f=0; f<M.F.size(); ++f) if(M.fal[f] && M.fpatch[f]) {any=true; break;}
    if(!any) return;

    for(int round=0; round<4; ++round)
    {
        for(int pass=0; pass<30; ++pass)
        {
            std::vector<std::pair<int,int> > E;
            for(size_t f=0; f<M.F.size(); ++f)
            if(M.fal[f] && M.fpatch[f])
            for(int q=0; q<3; ++q)
            {
                int a=M.F[f][q], b=M.F[f][(q+1)%3];
                E.emplace_back(std::min(a,b),std::max(a,b));
            }
            std::sort(E.begin(),E.end());
            E.erase(std::unique(E.begin(),E.end()),E.end());

            std::vector<std::pair<double,int> > cand;
            for(size_t i=0; i<E.size(); ++i)
            {
                int a=E[i].first, b=E[i].second;
                // inner patch edges and rim edges (a rim midpoint lies on the original
                // edge and stays fixed, so the original surface is not changed)
                M.edge_faces(a,b,tmp1);
                if(tmp1.size()!=2) continue;
                double r2=ml2(M.P[a],M.P[b]);
                if(r2>16.0/9.0) cand.emplace_back(-r2,int(i));
            }
            if(cand.empty()) break;
            std::sort(cand.begin(),cand.end());
            for(auto &c: cand)
            {
                if(M.nfaces_alive>prm.max_tri) break;
                split(E[c.second].first,E[c.second].second);
            }
        }
        flip_delaunay(true);
        fair_patch(1);
        if(prm.hole_fill==2) fair_patch(2);
        flip_delaunay(true);
    }
}

// mode 1: harmonic (membrane, uniform Laplacian = 0), mode 2: biharmonic (Laplacian^2 = 0)
// Gauss-Seidel on the free patch vertices, all other vertices fixed.
void remesher::fair_patch(int mode)
{
    std::vector<int> V;
    for(size_t v=0; v<M.P.size(); ++v)
    if(M.val[v] && M.vpatch[v] && !M.vf[v].empty()) V.push_back(int(v));
    if(V.empty()) return;

    std::unordered_map<int,std::vector<int> > nb;
    auto nbs=[&](int v) -> const std::vector<int>&
    {
        auto it=nb.find(v);
        if(it!=nb.end()) return it->second;
        std::vector<int> t;
        M.neighbors(v,t);
        return nb.emplace(v,std::move(t)).first->second;
    };
    auto lap=[&](int v)
    {
        const auto &N=nbs(v);
        vec3 c{0.0,0.0,0.0};
        for(int j: N) c=c+M.P[j];
        return (1.0/double(N.size()))*c - M.P[v];
    };

    double hmean=0.0;
    for(int v: V) hmean+=hiso(M.P[v]);
    hmean/=double(V.size());
    const double eps=1.0e-7*hmean;

    const int maxit = mode==1 ? 5000 : 3000;
    for(int it=0; it<maxit; ++it)
    {
        double maxd=0.0;
        for(int v: V)
        {
            vec3 d;
            if(mode==1)
            {
                d=1.6*lap(v);   // SOR
            }
            else
            {
                const auto &N=nbs(v);
                const double nv=double(N.size());
                vec3 r{0.0,0.0,0.0};
                double cvv=1.0;
                for(int j: N)
                {
                    r=r+(1.0/nv)*lap(j);
                    cvv+=1.0/(nv*double(nbs(j).size()));
                }
                r=r-lap(v);
                d=(-1.0/cvv)*r;
            }
            M.P[v]=M.P[v]+d;
            maxd=std::max(maxd,norm(d));
        }
        if(maxd<eps) break;
    }
}

void remesher::detect_features()
{
    std::vector<std::pair<int,int> > E;
    M.edges(E);
    std::vector<int> ftype;               // per feature edge: 0 dihedral, 1 boundary, 2 non-manifold/inconsistent
    std::vector<std::pair<int,int> > FE;

    for(auto &e: E)
    {
        int a=e.first, b=e.second;
        M.edge_faces(a,b,tmp1);
        int type=-1;
        if(tmp1.size()==1) {type=1; ++st.n_boundary_edges;}
        else if(tmp1.size()>2) {type=2; ++st.n_nonmanifold_edges;}
        else
        {
            int d1,d2;
            M.opposite(tmp1[0],a,b,d1);
            M.opposite(tmp1[1],a,b,d2);
            if(d1==d2) {type=2; ++st.n_inconsistent_edges;}
            else
            {
                vec3 n1=M.fnormal(tmp1[0]), n2=M.fnormal(tmp1[1]);
                double l1=norm(n1), l2=norm(n2);
                if(l1>2.0*tiny_area && l2>2.0*tiny_area && dot(n1,n2)/(l1*l2) < cos_feat) type=0;
            }
        }
        if(type>=0) {FE.emplace_back(a,b); ftype.push_back(type);}
    }

    // vertex classification
    std::vector<std::vector<int> > vfe(M.P.size());
    for(size_t i=0; i<FE.size(); ++i)
    {
        vfe[FE[i].first].push_back(int(i));
        vfe[FE[i].second].push_back(int(i));
    }
    for(size_t v=0; v<M.P.size(); ++v)
    {
        int n=int(vfe[v].size());
        if(n==0) M.lev[v]=0;
        else if(n!=2) M.lev[v]=2;
        else
        {
            int e0=vfe[v][0], e1=vfe[v][1];
            if(ftype[e0]!=ftype[e1] || ftype[e0]==2) M.lev[v]=2;
            else
            {
                int u = FE[e0].first==int(v) ? FE[e0].second : FE[e0].first;
                int w = FE[e1].first==int(v) ? FE[e1].second : FE[e1].first;
                vec3 d0=unit(M.P[v]-M.P[u]), d1=unit(M.P[w]-M.P[v]);
                M.lev[v] = dot(d0,d1) < cos_feat ? 2 : 1;
            }
        }
    }

    // chain feature edges into curves between corners
    std::vector<int> echain(FE.size(),-1);
    int nch=0;
    for(size_t i0=0; i0<FE.size(); ++i0)
    {
        if(echain[i0]>=0) continue;
        int ch=nch++;
        M.chain_bnd.push_back(ftype[i0]==1);
        std::vector<int> stack{int(i0)};
        echain[i0]=ch;
        while(!stack.empty())
        {
            int e=stack.back();
            stack.pop_back();
            for(int v: {FE[e].first, FE[e].second})
            {
                if(M.lev[v]!=1) continue;
                M.chain[v]=ch;
                for(int e2: vfe[v])
                if(echain[e2]<0) {echain[e2]=ch; stack.push_back(e2);}
            }
        }
    }
    for(size_t i=0; i<FE.size(); ++i)
    M.fe[ekey(FE[i].first,FE[i].second)]=echain[i];

    st.n_feature_edges=int(FE.size());
    st.n_chains=nch;

    // ideal valence of corners from their angle sum (invariant: corners never move)
    for(size_t v=0; v<M.P.size(); ++v)
    {
        if(M.lev[v]!=2) continue;
        ++st.n_corners;
        double th=0.0;
        bool bnd=false;
        for(int f: M.vf[v])
        {
            const auto &t=M.F[f];
            int q = t[0]==int(v) ? 0 : (t[1]==int(v) ? 1 : 2);
            th+=angle_at(M.P[v],M.P[t[(q+1)%3]],M.P[t[(q+2)%3]]);
        }
        for(int e: vfe[v]) if(ftype[e]==1) bnd=true;
        M.vtarget[v] = std::max(3.0, th/(PI_/3.0)) + (bnd ? 1.0 : 0.0);
    }
}

void remesher::build_reference()
{
    rP=M.P;
    std::vector<std::array<double,6> > bb;
    for(size_t f=0; f<M.F.size(); ++f)
    {
        if(!M.fal[f]) continue;
        vec3 n=M.fnormal(int(f));
        if(norm(n)<=2.0*tiny_area) continue;
        const auto &t=M.F[f];
        rF.push_back(t);
        rN.push_back(unit(n));
        std::array<double,6> b;
        for(int d=0; d<3; ++d)
        {
            b[d]  =std::min(rP[t[0]][d],std::min(rP[t[1]][d],rP[t[2]][d]));
            b[3+d]=std::max(rP[t[0]][d],std::max(rP[t[1]][d],rP[t[2]][d]));
        }
        bb.push_back(b);
    }
    rbvh.build(bb);

    std::vector<std::array<double,6> > sb;
    for(auto &kv: M.fe)
    {
        int a=int(kv.first>>32), b=int(kv.first & 0xffffffffULL);
        rS.push_back({a,b,kv.second});
    }
    std::sort(rS.begin(),rS.end());
    for(auto &s: rS)
    {
        std::array<double,6> b;
        for(int d=0; d<3; ++d)
        {
            b[d]=std::min(rP[s[0]][d],rP[s[1]][d]);
            b[3+d]=std::max(rP[s[0]][d],rP[s[1]][d]);
        }
        sb.push_back(b);
    }
    sbvh.build(sb);
}

long remesher::estimate_triangles() const
{
    double n=0.0;
    for(size_t f=0; f<rF.size(); ++f)
    {
        const auto &t=rF[f];
        n+=marea(rP[t[0]],rP[t[1]],rP[t[2]])/(0.25*std::sqrt(3.0));
    }
    return long(n)+1;
}

vec3 remesher::project_surface(const vec3 &q, const vec3 &n) const
{
    double best=1.0e300;
    vec3 bp=q;
    bool hasn = dot(n,n)>0.0;
    rbvh.closest(q,[&](int f)
    {
        if(hasn && dot(rN[f],n)<=0.0) return 1.0e300;
        const auto &t=rF[f];
        vec3 c=closest_on_triangle(q,rP[t[0]],rP[t[1]],rP[t[2]]);
        double d2=dist2(c,q);
        if(d2<best) bp=c;
        return d2;
    },best);
    return bp;
}

vec3 remesher::project_chain(const vec3 &q, int ch) const
{
    double best=1.0e300;
    vec3 bp=q;
    sbvh.closest(q,[&](int s)
    {
        if(rS[s][2]!=ch) return 1.0e300;
        vec3 c=closest_on_segment(q,rP[rS[s][0]],rP[rS[s][1]]);
        double d2=dist2(c,q);
        if(d2<best) bp=c;
        return d2;
    },best);
    return bp;
}

void remesher::split(int a, int b)
{
    M.edge_faces(a,b,tmp1);
    if(tmp1.empty()) return;
    std::uint64_t k=ekey(a,b);
    auto it=M.fe.find(k);
    int ch = it!=M.fe.end() ? it->second : -1;
    int m=M.add_vertex(0.5*(M.P[a]+M.P[b]), ch>=0 ? 1 : 0, ch);
    if(ch>=0)
    {
        M.P[m]=project_chain(M.P[m],ch);
        M.fe.erase(k);
        M.fe[ekey(a,m)]=ch;
        M.fe[ekey(m,b)]=ch;
    }
    std::vector<int> ef=tmp1;
    bool allpatch=true;
    for(int f: ef) if(!M.fpatch[f]) allpatch=false;
    M.vpatch[m] = allpatch ? 1 : 0;
    for(int f: ef)
    {
        auto t=M.F[f];
        int q=0;
        for(; q<3; ++q)
        {
            int u=t[q], w=t[(q+1)%3];
            if((u==a && w==b) || (u==b && w==a)) break;
        }
        int u=t[q], w=t[(q+1)%3], c=t[(q+2)%3];
        M.F[f]={u,m,c};
        M.vf[m].push_back(f);
        mesh::erase_val(M.vf[w],f);
        M.F.push_back({m,w,c});
        M.fal.push_back(1);
        M.fpatch.push_back(M.fpatch[f]);
        int g=int(M.F.size())-1;
        M.vf[m].push_back(g);
        M.vf[w].push_back(g);
        M.vf[c].push_back(g);
        ++M.nfaces_alive;
    }
}

bool remesher::faces_ok_after_move(int v, const vec3 &np) const
{
    for(int f: M.vf[v])
    {
        const auto &t=M.F[f];
        vec3 p0=M.P[t[0]], p1=M.P[t[1]], p2=M.P[t[2]];
        vec3 n0=cross(p1-p0,p2-p0);
        if(t[0]==v) p0=np; else if(t[1]==v) p1=np; else p2=np;
        vec3 n1=cross(p1-p0,p2-p0);
        double l1=norm(n1), l0=norm(n0);
        if(l1<=2.0*tiny_area) return false;
        if(l0>2.0*tiny_area && dot(n0,n1)<0.2*l0*l1) return false;
    }
    return true;
}

bool remesher::collapse(int a, int b)
{
    if(!M.val[a] || !M.val[b]) return false;
    if(M.lev[a]==2) return false;
    bool fedge=M.is_feature(a,b);
    if(M.lev[a]==1 && !fedge) return false;

    std::vector<int> ef;
    M.edge_faces(a,b,ef);
    if(ef.empty() || ef.size()>2) return false;

    // link condition
    std::vector<int> Na, Nb, opp;
    M.neighbors(a,Na);
    M.neighbors(b,Nb);
    for(int f: ef) {int d; opp.push_back(M.opposite(f,a,b,d));}
    std::sort(opp.begin(),opp.end());
    if(opp.size()==2 && opp[0]==opp[1]) return false;
    std::vector<int> common;
    std::set_intersection(Na.begin(),Na.end(),Nb.begin(),Nb.end(),std::back_inserter(common));
    if(common!=opp) return false;
    for(int c: opp) if(M.valence(c)<=3) return false;
    if(ef.size()==2 && Na.size()==3 && Nb.size()==3) return false;   // tetrahedron

    // no new long edges
    const vec3 &pb=M.P[b];
    for(int x: Na)
    {
        if(x==b) continue;
        if(ml2(pb,M.P[x]) > 16.0/9.0) return false;
    }

    // no flipped or degenerate faces, limited normal rotation
    for(int f: M.vf[a])
    {
        if(std::find(ef.begin(),ef.end(),f)!=ef.end()) continue;
        const auto &t=M.F[f];
        vec3 p0=M.P[t[0]], p1=M.P[t[1]], p2=M.P[t[2]];
        vec3 n0=cross(p1-p0,p2-p0);
        if(t[0]==a) p0=pb; else if(t[1]==a) p1=pb; else p2=pb;
        vec3 n1=cross(p1-p0,p2-p0);
        double l0=norm(n0), l1=norm(n1);
        if(l1<=2.0*tiny_area) return false;
        if(l0>2.0*tiny_area && dot(n0,n1)<0.7*l0*l1) return false;
        if(tri_minangle(p0,p1,p2)<1.0*PI_/180.0) return false;
    }

    // perform
    for(int f: ef) M.kill_face(f);
    std::vector<int> fa=M.vf[a];
    for(int f: fa)
    {
        auto &t=M.F[f];
        for(int q=0; q<3; ++q) if(t[q]==a) t[q]=b;
        M.vf[b].push_back(f);
    }
    M.vf[a].clear();
    M.val[a]=0;

    for(int x: Na)
    {
        auto it=M.fe.find(ekey(a,x));
        if(it==M.fe.end()) continue;
        int ch=it->second;
        M.fe.erase(it);
        if(x!=b) M.fe[ekey(b,x)]=ch;
    }
    return true;
}

bool remesher::flip_ok(int a, int b, int &f1, int &f2, int &c, int &d)
{
    if(M.is_feature(a,b)) return false;
    M.edge_faces(a,b,tmp1);
    if(tmp1.size()!=2) return false;
    int d1,d2;
    int o1=M.opposite(tmp1[0],a,b,d1);
    int o2=M.opposite(tmp1[1],a,b,d2);
    if(d1==d2) return false;
    if(d1>0) {f1=tmp1[0]; c=o1; f2=tmp1[1]; d=o2;}
    else     {f1=tmp1[1]; c=o2; f2=tmp1[0]; d=o1;}
    if(c==d) return false;
    for(int f: M.vf[c]) if(mesh::has(M.F[f],d)) return false;   // edge (c,d) exists

    // geometry: f1=(a,b,c), f2=(b,a,d) -> g1=(c,a,d), g2=(d,b,c)
    const vec3 &pa=M.P[a], &pb=M.P[b], &pc=M.P[c], &pd=M.P[d];
    vec3 n1=cross(pb-pa,pc-pa), n2=cross(pa-pb,pd-pb);
    vec3 g1=cross(pa-pc,pd-pc), g2=cross(pb-pd,pc-pd);
    double lg1=norm(g1), lg2=norm(g2);
    if(lg1<=2.0*tiny_area || lg2<=2.0*tiny_area) return false;
    vec3 nm=unit(n1)+unit(n2);
    if(dot(g1,nm)<=0.0 || dot(g2,nm)<=0.0) return false;
    if(dot(g1,g2) < cos_feat*lg1*lg2) return false;
    return true;
}

void remesher::do_flip(int a, int b, int f1, int f2, int c, int d)
{
    M.F[f1]={c,a,d};
    M.F[f2]={d,b,c};
    mesh::erase_val(M.vf[a],f2);
    mesh::erase_val(M.vf[b],f1);
    M.vf[c].push_back(f2);
    M.vf[d].push_back(f1);
}

int remesher::split_long()
{
    int total=0;
    std::vector<std::pair<int,int> > E;
    for(int pass=0; pass<40; ++pass)
    {
        M.edges(E);
        std::vector<std::pair<double,int> > cand;
        for(size_t i=0; i<E.size(); ++i)
        {
            const vec3 &pa=M.P[E[i].first], &pb=M.P[E[i].second];
            double r2=ml2(pa,pb);
            if(r2 > 16.0/9.0) cand.emplace_back(-r2,int(i));
        }
        if(cand.empty()) break;
        std::sort(cand.begin(),cand.end());
        int n=0;
        for(auto &c: cand)
        {
            if(M.nfaces_alive>prm.max_tri) break;
            split(E[c.second].first,E[c.second].second);
            ++n;
        }
        total+=n;
        if(M.nfaces_alive>prm.max_tri) break;
    }
    return total;
}

int remesher::collapse_short()
{
    int total=0;
    std::vector<std::pair<int,int> > E;
    for(int pass=0; pass<2; ++pass)
    {
        M.edges(E);
        std::vector<std::pair<double,int> > cand;
        for(size_t i=0; i<E.size(); ++i)
        {
            const vec3 &pa=M.P[E[i].first], &pb=M.P[E[i].second];
            double r2=ml2(pa,pb);
            if(r2 < 16.0/25.0) cand.emplace_back(r2,int(i));
        }
        if(cand.empty()) break;
        std::sort(cand.begin(),cand.end());
        int n=0;
        for(auto &c: cand)
        {
            int a=E[c.second].first, b=E[c.second].second;
            if(!M.val[a] || !M.val[b]) continue;
            const vec3 &pa=M.P[a], &pb=M.P[b];
            if(ml2(pa,pb) >= 16.0/25.0) continue;
            // remove the less constrained vertex first
            bool done=false;
            if(M.lev[a]<=M.lev[b]) done = collapse(a,b) || collapse(b,a);
            else                   done = collapse(b,a) || collapse(a,b);
            if(done) ++n;
        }
        total+=n;
        if(n==0) break;
    }
    return total;
}

int remesher::flip_valence()
{
    int n=0;
    std::vector<std::pair<int,int> > E;
    M.edges(E);
    for(auto &e: E)
    {
        int a=e.first, b=e.second, f1, f2, c, d;
        if(!flip_ok(a,b,f1,f2,c,d)) continue;
        int va=M.valence(a), vb=M.valence(b), vc=M.valence(c), vd=M.valence(d);
        if(va<=3 || vb<=3) continue;
        double ta=target_valence(a), tb=target_valence(b), tc=target_valence(c), td=target_valence(d);
        auto sq=[](double x){return x*x;};
        double before=sq(va-ta)+sq(vb-tb)+sq(vc-tc)+sq(vd-td);
        double after =sq(va-1-ta)+sq(vb-1-tb)+sq(vc+1-tc)+sq(vd+1-td);
        if(after<before)
        {
            // do not trade valence for very poor angles
            vec3 pa,pb,pc,pd;
            metric_quad(M.P[a],M.P[b],M.P[c],M.P[d],pa,pb,pc,pd);
            double mb=std::min(tri_minangle(pa,pb,pc),tri_minangle(pb,pa,pd));
            double ma=std::min(tri_minangle(pc,pa,pd),tri_minangle(pd,pb,pc));
            if(ma<0.5*mb) continue;
            do_flip(a,b,f1,f2,c,d);
            ++n;
        }
    }
    return n;
}

int remesher::flip_delaunay(bool patch_only)
{
    int total=0;
    std::vector<std::pair<int,int> > E;
    for(int pass=0; pass<10; ++pass)
    {
        int n=0;
        M.edges(E);
        for(auto &e: E)
        {
            int a=e.first, b=e.second, f1, f2, c, d;
            if(!flip_ok(a,b,f1,f2,c,d)) continue;
            if(patch_only && (!M.fpatch[f1] || !M.fpatch[f2])) continue;
            if(M.valence(a)<=3 || M.valence(b)<=3) continue;
            vec3 pa,pb,pc,pd;
            metric_quad(M.P[a],M.P[b],M.P[c],M.P[d],pa,pb,pc,pd);
            if(angle_at(pc,pa,pb)+angle_at(pd,pa,pb) <= PI_+1.0e-6) continue;
            double mb=std::min(tri_minangle(pa,pb,pc),tri_minangle(pb,pa,pd));
            double ma=std::min(tri_minangle(pc,pa,pd),tri_minangle(pd,pb,pc));
            if(ma<=mb*(1.0+1.0e-6)) continue;
            do_flip(a,b,f1,f2,c,d);
            ++n;
        }
        total+=n;
        if(n==0) break;
    }
    return total;
}

vec3 remesher::vnormal(int v) const
{
    vec3 n{0.0,0.0,0.0};
    for(int f: M.vf[v]) n=n+M.fnormal(f);
    return unit(n);
}

void remesher::relax()
{
    const double lambda=0.8;
    for(size_t v=0; v<M.P.size(); ++v)
    {
        if(!M.val[v] || M.lev[v]==2 || M.vf[v].empty()) continue;
        const vec3 p=M.P[v];
        vec3 np;
        if(M.lev[v]==0)
        {
            vec3 c{0.0,0.0,0.0};
            double wsum=0.0;
            for(int f: M.vf[v])
            {
                const auto &t=M.F[f];
                vec3 g=(1.0/3.0)*(M.P[t[0]]+M.P[t[1]]+M.P[t[2]]);
                // CVT in the metric: weight = metric area^2 / area (isotropic: A/h^4)
                double A=0.5*norm(M.fnormal(f));
                if(A<=tiny_area) continue;
                double Am=marea(M.P[t[0]],M.P[t[1]],M.P[t[2]]);
                double w=Am*Am/A;
                c=c+w*g;
                wsum+=w;
            }
            if(wsum<=0.0) continue;
            c=(1.0/wsum)*c;
            vec3 n=vnormal(int(v));
            vec3 d=c-p;
            d=d-dot(d,n)*n;
            np=p+lambda*d;
            np=project_surface(np,n);
        }
        else
        {
            if(M.chain[v]<0) continue;
            M.neighbors(int(v),tmp2);
            int u=-1, w=-1, cnt=0;
            for(int x: tmp2)
            if(M.is_feature(int(v),x)) {if(cnt==0) u=x; else w=x; ++cnt;}
            if(cnt!=2) continue;
            const vec3 &pu=M.P[u], &pw=M.P[w];
            // weighted so that both feature edges approach their target lengths
            // target length along each feature edge: |d| / metric length
            double lu=std::sqrt(dist2(p,pu)), lw=std::sqrt(dist2(p,pw));
            double mu=std::sqrt(ml2(p,pu)), mw=std::sqrt(ml2(p,pw));
            if(lu<=0.0 || lw<=0.0 || mu<=0.0 || mw<=0.0) continue;
            double wu=mu/lu, ww=mw/lw;
            vec3 c=(1.0/(wu+ww))*(wu*pu+ww*pw);
            np=p+lambda*(c-p);
            np=project_chain(np,M.chain[v]);
        }
        if(faces_ok_after_move(int(v),np)) M.P[v]=np;
    }
}

bool remesher::run(const std::vector<vec3> &in, std::vector<vec3> &out)
{
    cos_feat=std::cos(prm.feature_angle*PI_/180.0);
    st.ntri_in=int(in.size()/3);
    soup_stats(in,st.area_in,st.vol_in,st.minangle_in,st.q_mean_in,st.q_min_in,st.frac_q05_in);

    if(in.size()<9) return false;

    weld(in);
    st.n_tjunctions=repair_tjunctions(prm.tjunction_tol*diag,false);
    if(prm.hole_fill>0)
    {
        st.n_gap_merges=close_gaps();
        fill_holes();
        refine_and_fair_patches();
    }
    detect_features();
    build_reference();

    st.ntri_estimate=estimate_triangles();
    if(st.ntri_estimate>prm.max_tri) return false;

    for(int it=0; it<prm.iterations; ++it)
    {
        split_long();
        collapse_short();
        flip_valence();
        relax();
    }
    for(int it=0; it<prm.smooth_final; ++it)
    {
        flip_delaunay();
        relax();
    }
    flip_delaunay();

    // output
    out.clear();
    out.reserve(3*size_t(M.nfaces_alive));
    for(size_t f=0; f<M.F.size(); ++f)
    if(M.fal[f])
    for(int q=0; q<3; ++q) out.push_back(M.P[M.F[f][q]]);

    st.ntri_out=int(out.size()/3);
    st.nvert_out=0;
    for(size_t v=0; v<M.P.size(); ++v) if(M.val[v] && !M.vf[v].empty()) ++st.nvert_out;
    soup_stats(out,st.area_out,st.vol_out,st.minangle_out,st.q_mean_out,st.q_min_out,st.frac_q05_out);

    std::vector<std::pair<int,int> > E;
    M.edges(E);
    double s=0.0;
    st.Lh_min=1.0e300;
    st.Lh_max=0.0;
    for(auto &e: E)
    {
        const vec3 &pa=M.P[e.first], &pb=M.P[e.second];
        double r=std::sqrt(ml2(pa,pb));
        s+=r;
        st.Lh_min=std::min(st.Lh_min,r);
        st.Lh_max=std::max(st.Lh_max,r);
    }
    st.Lh_mean = E.empty() ? 0.0 : s/double(E.size());

    // quality in the metric (shape relative to the grid cells)
    st.q_mean_metric=0.0;
    st.q_min_metric=1.0;
    for(size_t t=0; t<out.size()/3; ++t)
    {
        vec3 sc=hv((1.0/3.0)*(out[3*t]+out[3*t+1]+out[3*t+2]));
        double q=tri_quality(scale(out[3*t],sc),scale(out[3*t+1],sc),scale(out[3*t+2],sc));
        st.q_mean_metric+=q;
        st.q_min_metric=std::min(st.q_min_metric,q);
    }
    if(!out.empty()) st.q_mean_metric/=double(out.size()/3);
    return st.ntri_out>0;
}

} // namespace


bool sixdof_remesh::remesh(const std::vector<vec3> &in, std::vector<vec3> &out, const sizing_func &h, const params &prm, stats &st)
{
    metric_func H=[&h](double x, double y, double z)
    {
        double v=h(x,y,z);
        return vec3{v,v,v};
    };
    return remesh(in,out,H,prm,st);
}

bool sixdof_remesh::remesh(const std::vector<vec3> &in, std::vector<vec3> &out, const metric_func &H, const params &prm, stats &st)
{
    st=stats();
    remesher r(H,prm,st);
    bool ok=r.run(in,out);
    if(!ok) out=in;
    st.ok=ok;
    return ok;
}

void sixdof_remesh::print_stats(std::ostream &os, const stats &st)
{
    const std::ios_base::fmtflags fl=os.flags();
    const std::streamsize pr=os.precision();
    os<<std::setprecision(4);
    os<<"  triangles      : "<<st.ntri_in<<" -> "<<st.ntri_out<<"  (vertices "<<st.nvert_out<<")"<<std::endl;
    os<<"  features       : "<<st.n_feature_edges<<" edges, "<<st.n_chains<<" curves, "<<st.n_corners<<" corners";
    if(st.n_boundary_edges>0) os<<", "<<st.n_boundary_edges<<" open boundary edges";
    if(st.n_nonmanifold_edges>0) os<<", "<<st.n_nonmanifold_edges<<" non-manifold edges";
    if(st.n_inconsistent_edges>0) os<<", "<<st.n_inconsistent_edges<<" edges with inconsistent orientation";
    os<<std::endl;
    if(st.n_tjunctions>0 || st.n_gap_merges>0 || st.n_holes_filled>0 || st.n_holes_open>0)
    {
        os<<"  repairs        : "<<st.n_tjunctions<<" T-junctions, "<<st.n_gap_merges<<" gap vertices merged, "
          <<st.n_holes_filled<<" holes filled ("<<st.n_hole_edges<<" rim edges)";
        if(st.n_holes_open>0) os<<", "<<st.n_holes_open<<" holes could not be filled";
        os<<std::endl;
    }
    os<<"  min angle [deg]: "<<st.minangle_in<<" -> "<<st.minangle_out<<std::endl;
    os<<"  quality mean   : "<<st.q_mean_in<<" -> "<<st.q_mean_out<<"   min: "<<st.q_min_in<<" -> "<<st.q_min_out
      <<"   q<0.5: "<<100.0*st.frac_q05_in<<"% -> "<<100.0*st.frac_q05_out<<"%"<<std::endl;
    os<<"  edge/target    : mean "<<st.Lh_mean<<"  min "<<st.Lh_min<<"  max "<<st.Lh_max
      <<"   (target = X 186 cells)   quality in grid metric: mean "<<st.q_mean_metric<<"  min "<<st.q_min_metric<<std::endl;
    os<<"  area           : "<<st.area_in<<" -> "<<st.area_out<<"   volume: "<<st.vol_in<<" -> "<<st.vol_out<<std::endl;
    os.flags(fl);
    os.precision(pr);
}
