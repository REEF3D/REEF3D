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


#include"seastate_collocation.h"
#include"seastate_forcing.h"
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<fstream>
#include<iomanip>
#include<map>
#include<numeric>
#include<sstream>

double seastate_colloc_stats::rmse() const
{
    if(n<=0)
    return 0.0;
    const double d2 = (sm2 - 2.0*som + so2)/n;
    return std::sqrt(std::max(d2,0.0));
}

double seastate_colloc_stats::si() const
{
    if(n<=0 || !(so>0.0))
    return 0.0;
    const double b = bias();
    const double var = (sm2 - 2.0*som + so2)/n - b*b;
    return std::sqrt(std::max(var,0.0))/(so/n);
}

double seastate_colloc_stats::r() const
{
    if(n<2)
    return 0.0;
    const double co = so2/n - (so/n)*(so/n), cm = sm2/n - (sm/n)*(sm/n), cx = som/n - (so/n)*(sm/n);
    return (co>0.0 && cm>0.0) ? cx/std::sqrt(co*cm) : 0.0;
}

static std::string trim(const std::string &s)
{
    const size_t a = s.find_first_not_of(" \t\r\"");
    if(a==std::string::npos)
    return "";
    const size_t b = s.find_last_not_of(" \t\r\"");
    return s.substr(a,b-a+1);
}

static std::vector<std::string> split(const std::string &line)
{
    std::vector<std::string> v;
    std::string c;
    std::istringstream in(line);
    while(std::getline(in,c,','))
    v.push_back(trim(c));
    if(!line.empty() && line.back()==',')
    v.push_back("");
    return v;
}

bool seastate_collocation::parse_time(const std::string &s0, double &t)
{
    const std::string s = trim(s0);

    // YYYYMMDD.HHMMSS
    if(s.size()>=8 && s.find('-')==std::string::npos && s.find(':')==std::string::npos)
    {
    char *end = nullptr;
    const double v = std::strtod(s.c_str(),&end);
        if(end==s.c_str() || *end!='\0' || !(v>=19000101.0 && v<30000101.0))
        return false;
    t = seastate_datetime(v);
    return true;
    }

    // YYYY-MM-DD[ T]HH:MM:SS[.sss]
    int Y, M, D, h = 0, m = 0;
    double sec = 0.0;
    char sep = ' ';
    const int k = std::sscanf(s.c_str(),"%d-%d-%d%c%d:%d:%lf",&Y,&M,&D,&sep,&h,&m,&sec);
    if(k!=3 && k<6)
    return false;
    if(k==3)
    {
    h = 0;
    m = 0;
    sec = 0.0;
    }
    if(M<1 || M>12 || D<1 || D>31 || h<0 || h>23 || m<0 || m>59 || sec<0.0 || sec>=61.0)
    return false;

    const double whole = std::floor(sec);
    const double v = Y*10000.0 + M*100.0 + D + (h*10000.0 + m*100.0 + whole)*1.0e-6;
    t = seastate_datetime(v) + (sec - whole);
    return true;
}

std::string seastate_collocation::stamp(double t)
{
    // proleptic Gregorian (H. Hinnant), as the nesting output
    const long s = long(std::floor(t + 0.5));
    long z = (s>=0 ? s : s-86399)/86400;
    const long sod = s - z*86400;
    z += 719468;
    const long era = (z>=0 ? z : z-146096)/146097;
    const long doe = z - era*146097;
    const long yoe = (doe - doe/1460 + doe/36524 - doe/146096)/365;
    const long doy = doe - (365*yoe + yoe/4 - yoe/100);
    const long mp = (5*doy + 2)/153;
    const long d = doy - (153*mp + 2)/5 + 1;
    const long m = mp + (mp<10 ? 3 : -9);
    const long y = yoe + era*400 + (m<=2 ? 1 : 0);

    char buf[64];
    std::snprintf(buf,sizeof(buf),"%04ld%02ld%02ld.%02ld%02ld%02ld",y,m,d,sod/3600,(sod/60)%60,sod%60);
    return std::string(buf);
}

bool seastate_collocation::read(const std::string &file, double tref, std::string &error)
{
    obs.clear();
    sources.clear();
    ids.clear();

    std::ifstream in(file);
    if(!in.good())
    {
    error = "cannot open "+file;
    return false;
    }

    std::string line;
    std::vector<std::string> head;
    long lno = 0;

    while(std::getline(in,line))
    {
    ++lno;
    const std::string l = trim(line.substr(0,line.find('#')));
    if(l.empty())
    continue;
    head = split(l);
    break;
    }

    auto col = [&](const char *name) -> int
    {
        for(size_t k=0; k<head.size(); ++k)
        {
        std::string h = head[k];
        std::transform(h.begin(),h.end(),h.begin(),::tolower);
        if(h==name)
        return int(k);
        }
    return -1;
    };

    const int ct = col("time"), cm = col("t"), cx = col("x"), cy = col("y"), ch = col("hs"), cs = col("source"), ci = col("id");
    const int ctp = col("tp"), ct1 = col("tm01"), ct0 = col("tm10");
    model_time = (ct<0 && cm>=0);
    has_tp = (ctp>=0);
    has_tm01 = (ct1>=0);
    has_tm10 = (ct0>=0);

    if((ct<0 && cm<0) || cx<0 || cy<0 || ch<0)
    {
    error = file+": the header needs the columns time (or t), x, y, hs";
    return false;
    }

    std::map<std::string,int> smap, imap;
    const int need = std::max({ct,cm,cx,cy,ch,cs,ci});

    while(std::getline(in,line))
    {
    ++lno;
    const std::string l = trim(line.substr(0,line.find('#')));
    if(l.empty())
    continue;

    const std::vector<std::string> c = split(l);
        if(int(c.size())<=need)
        {
        error = file+", line "+std::to_string(lno)+": "+std::to_string(c.size())+" columns, the header has "+std::to_string(head.size());
        return false;
        }

    seastate_obs o;
        if(model_time)
        o.t = std::strtod(c[cm].c_str(),nullptr);
        else
        {
        double ta;
            if(!parse_time(c[ct],ta))
            {
            error = file+", line "+std::to_string(lno)+": time '"+c[ct]+"' (YYYYMMDD.HHMMSS or YYYY-MM-DD HH:MM:SS)";
            return false;
            }
        o.t = ta - tref;
        }

    char *e1 = nullptr, *e2 = nullptr, *e3 = nullptr;
    o.x = std::strtod(c[cx].c_str(),&e1);
    o.y = std::strtod(c[cy].c_str(),&e2);
    o.hs = std::strtod(c[ch].c_str(),&e3);
        if(e1==c[cx].c_str() || e2==c[cy].c_str() || e3==c[ch].c_str() || !std::isfinite(o.x) || !std::isfinite(o.y) || !std::isfinite(o.hs))
        {
        error = file+", line "+std::to_string(lno)+": x, y and hs must be numbers";
        return false;
        }

    // the observed periods: empty, nan or not positive -> none
    auto period = [&](int cc) -> double
    {
        if(cc<0 || cc>=int(c.size()) || c[cc].empty())
        return std::nan("");
    char *e = nullptr;
    const double v = std::strtod(c[cc].c_str(),&e);
    return (e!=c[cc].c_str() && std::isfinite(v) && v>0.0) ? v : std::nan("");
    };
    o.tp = period(ctp);
    o.tm01 = period(ct1);
    o.tm10 = period(ct0);

    const std::string sname = (cs>=0 && !c[cs].empty()) ? c[cs] : "all";
    const std::string iname = (ci>=0) ? c[ci] : "";
        if(!smap.count(sname))
        {
        smap[sname] = int(sources.size());
        sources.push_back(sname);
        }
        if(!imap.count(iname))
        {
        imap[iname] = int(ids.size());
        ids.push_back(iname);
        }
    o.src = smap[sname];
    o.id = imap[iname];
    obs.push_back(o);
    }

    // by time; equal times keep the order of the file
    std::stable_sort(obs.begin(),obs.end(),[](const seastate_obs &a, const seastate_obs &b){return a.t<b.t;});

    per_source.assign(sources.size(),seastate_colloc_stats());
    total = seastate_colloc_stats();
    reg.assign(size_t(reg_off(nkind))*(sources.size()+1),seastate_colloc_stats());
    per.assign(3*(sources.size()+1),seastate_colloc_stats());
    qo.assign(sources.size()+1,std::vector<float>());
    qm.assign(sources.size()+1,std::vector<float>());
    return true;
}

void seastate_collocation::count(int src, double o, double m)
{
    if(src>=0 && src<int(per_source.size()))
    per_source[src].add(o,m);
    total.add(o,m);
}

// ------------------------------------------------------------------ the regimes, the periods, the quantiles

static const double hs_edges[] = {0.0, 1.0, 2.0, 4.0, 6.0, 8.0};
static const double depth_edges[] = {0.0, 10.0, 20.0, 50.0, 100.0, 200.0};
static const double age_edges[] = {0.0, 0.8, 1.2, 2.0};

const char *seastate_collocation::kind_name(int k)
{
    return k==0 ? "observed Hs [m]" : (k==1 ? "model depth [m]" : "wave age cp/U10");
}

int seastate_collocation::nclass(int k)
{
    return k==0 ? 6 : (k==1 ? 6 : 4);
}

std::string seastate_collocation::class_name(int k, int c)
{
    const double *e = (k==0) ? hs_edges : (k==1 ? depth_edges : age_edges);
    const int n = nclass(k);
    std::ostringstream s;
    s<<e[c];
    if(c+1<n)
    s<<"-"<<e[c+1];
    else
    s<<"+";
    return s.str();
}

double seastate_collocation::phase_speed(double T, double d)
{
    const double g = 9.81, pi = 3.14159265358979323846;
    if(!(T>0.0))
    return 0.0;
    const double w = 2.0*pi/T;
    double k = w*w/g;
    if(d>0.0)
    {
    // Newton on w^2 = g k tanh(k d) from the deep-water k (or the shallow-water k if larger)
    k = std::max(k,w/std::sqrt(g*d));
        for(int it=0; it<50; ++it)
        {
        const double th = std::tanh(k*d);
        const double f = g*k*th - w*w, df = g*th + g*k*d*(1.0 - th*th);
        const double dk = f/df;
        k -= dk;
            if(std::fabs(dk)<1.0e-12*k)
            break;
        }
    }
    return w/k;
}

int seastate_collocation::classify(int k, const seastate_matchup &m)
{
    double v;
    if(k==0)
    v = m.hs_obs;
    else if(k==1)
    v = m.depth;
    else
    {
        if(!(m.u10>=1.0))
        return nclass(2)-1;
    v = phase_speed(m.tp,m.depth)/m.u10;
    }
    const double *e = (k==0) ? hs_edges : (k==1 ? depth_edges : age_edges);
    int c = 0;
    while(c+1<nclass(k) && v>=e[c+1])
    ++c;
    return c;
}

int seastate_collocation::reg_off(int k) const
{
    int o = 0;
    for(int q=0; q<k; ++q)
    o += nclass(q);
    return o;
}

const seastate_colloc_stats &seastate_collocation::regime(int k, int c, int src) const
{
    const size_t ns = sources.size()+1;
    return reg[size_t(reg_off(k)+c)*ns + size_t(src<0 ? int(sources.size()) : src)];
}

const seastate_colloc_stats &seastate_collocation::period(int q, int src) const
{
    const size_t ns = sources.size()+1;
    return per[size_t(q)*ns + size_t(src<0 ? int(sources.size()) : src)];
}

void seastate_collocation::count(const seastate_matchup &m)
{
    count(m.src,m.hs_obs,m.hs);
    const size_t ns = sources.size()+1;
    if(reg.size()!=size_t(reg_off(nkind))*ns)
    reg.assign(size_t(reg_off(nkind))*ns,seastate_colloc_stats());
    if(per.size()!=3*ns)
    per.assign(3*ns,seastate_colloc_stats());
    if(qo.size()!=ns)
    {
    qo.assign(ns,std::vector<float>());
    qm.assign(ns,std::vector<float>());
    }

    const bool sv = (m.src>=0 && m.src<int(sources.size()));
    for(int k=0; k<nkind; ++k)
    {
    const size_t r = size_t(reg_off(k)+classify(k,m))*ns;
        if(sv)
        reg[r + size_t(m.src)].add(m.hs_obs,m.hs);
    reg[r + ns-1].add(m.hs_obs,m.hs);
    }

    const double po[3] = {m.tp_obs, m.tm01_obs, m.tm10_obs}, pm[3] = {m.tp, m.tm01, m.tm10};
    for(int q=0; q<3; ++q)
        if(std::isfinite(po[q]) && po[q]>0.0)
        {
            if(sv)
            per[size_t(q)*ns + size_t(m.src)].add(po[q],pm[q]);
        per[size_t(q)*ns + ns-1].add(po[q],pm[q]);
        }

    if(sv)
    {
    qo[size_t(m.src)].push_back(float(m.hs_obs));
    qm[size_t(m.src)].push_back(float(m.hs));
    }
    qo[ns-1].push_back(float(m.hs_obs));
    qm[ns-1].push_back(float(m.hs));
}

double seastate_collocation::quantile(int src, double f, bool model) const
{
    const size_t s = size_t(src<0 ? int(sources.size()) : src);
    if(s>=qo.size())
    return 0.0;
    std::vector<float> v = model ? qm[s] : qo[s];
    if(v.empty())
    return 0.0;
    std::sort(v.begin(),v.end());
    // linear between the order statistics (the definition of numpy.quantile, 'linear')
    const double h = std::min(std::max(f,0.0),1.0)*double(v.size()-1);
    const size_t a = size_t(std::floor(h));
    const size_t b = std::min(a+1,v.size()-1);
    return double(v[a]) + (h - double(a))*(double(v[b]) - double(v[a]));
}

static void stats_row(std::ostringstream &s, const seastate_colloc_stats &a)
{
    s<<" \t "<<a.n<<" \t "<<a.mean_obs()<<" \t "<<a.bias()<<" \t "<<a.rmse()<<" \t "<<a.si()<<" \t "<<a.r();
}

std::string seastate_collocation::regimes_table() const
{
    std::ostringstream s;
    s<<"# REEF3D::SEASTATE online collocation (A 787 1): Hs [m] by regime, model - observation"<<"\n";
    s<<"# kind class source: N, mean observation, bias, RMSE, scatter index, correlation; source -1 all"<<"\n";
    for(int k=0; k<nkind; ++k)
    {
    s<<"#   kind "<<k<<": "<<kind_name(k)<<", classes";
    for(int c=0; c<nclass(k); ++c)
    s<<" "<<c<<" ("<<class_name(k,c)<<")";
    s<<"\n";
    }
    for(size_t q=0; q<sources.size(); ++q)
    s<<"#   source "<<q<<" "<<sources[q]<<"\n";
    s<<std::setprecision(8);
    for(int k=0; k<nkind; ++k)
    for(int c=0; c<nclass(k); ++c)
    for(int q=0; q<=int(sources.size()); ++q)
    {
    const int src = (q<int(sources.size())) ? q : -1;
    s<<k<<" \t "<<c<<" \t "<<src;
    stats_row(s,regime(k,c,src));
    s<<"\n";
    }
    return s.str();
}

std::string seastate_collocation::periods_table() const
{
    std::ostringstream s;
    s<<"# REEF3D::SEASTATE online collocation (A 787 1): periods [s] where observed, model - observation"<<"\n";
    s<<"# period source: N, mean observation, bias, RMSE, scatter index, correlation; period 0 Tp, 1 Tm01, 2 Tm-1,0;"
     <<" source -1 all"<<"\n";
    s<<"#   the file has: tp "<<(has_tp ? 1 : 0)<<", tm01 "<<(has_tm01 ? 1 : 0)<<", tm10 "<<(has_tm10 ? 1 : 0)<<"\n";
    for(size_t q=0; q<sources.size(); ++q)
    s<<"#   source "<<q<<" "<<sources[q]<<"\n";
    s<<std::setprecision(8);
    for(int k=0; k<3; ++k)
    for(int q=0; q<=int(sources.size()); ++q)
    {
    const int src = (q<int(sources.size())) ? q : -1;
    s<<k<<" \t "<<src;
    stats_row(s,period(k,src));
    s<<"\n";
    }
    return s.str();
}

std::string seastate_collocation::quantiles_table() const
{
    static const double fq[] = {0.5, 0.9, 0.95, 0.99, 0.999, 1.0};
    std::ostringstream s;
    s<<"# REEF3D::SEASTATE online collocation (A 787 1): quantiles of Hs [m], the observed and the model sorted each on its own"<<"\n";
    s<<"# source quantile: N, observed, model; quantile 1 the largest; source -1 all"<<"\n";
    for(size_t q=0; q<sources.size(); ++q)
    s<<"#   source "<<q<<" "<<sources[q]<<"\n";
    s<<std::setprecision(8);
    for(int q=0; q<=int(sources.size()); ++q)
    {
    const int src = (q<int(sources.size())) ? q : -1;
    const size_t n = qo.empty() ? 0 : qo[size_t(q)].size();
        for(double f : fq)
        s<<src<<" \t "<<f<<" \t "<<n<<" \t "<<quantile(src,f,false)<<" \t "<<quantile(src,f,true)<<"\n";
    }
    return s.str();
}

std::string seastate_collocation::summary() const
{
    std::ostringstream s;
    s<<std::fixed;
    auto row = [&](const std::string &name, const seastate_colloc_stats &a)
    {
    s<<"   "<<std::left<<std::setw(12)<<name<<std::right<<std::setw(9)<<a.n<<std::setprecision(3)
     <<std::setw(10)<<a.mean_obs()<<std::setw(10)<<a.bias()<<std::setw(10)<<a.rmse()<<std::setw(10)<<a.si()<<std::setw(10)<<a.r()<<"\n";
    };
    s<<"   "<<std::left<<std::setw(12)<<"source"<<std::right<<std::setw(9)<<"N"<<std::setw(10)<<"mean obs"<<std::setw(10)<<"bias"
     <<std::setw(10)<<"RMSE"<<std::setw(10)<<"SI"<<std::setw(10)<<"r"<<"\n";
    for(size_t k=0; k<sources.size(); ++k)
    row(sources[k],per_source[k]);
    if(sources.size()>1)
    row("all",total);
    return s.str();
}
