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

#include"6DOF_brain_link.h"
#include<cstdlib>
#include<cstring>
#include<chrono>
#include<thread>

#ifndef _WIN32
#include<sys/socket.h>
#include<sys/select.h>
#include<netinet/in.h>
#include<netinet/tcp.h>
#include<arpa/inet.h>
#include<unistd.h>
#endif

sixdof_brain_link::~sixdof_brain_link()
{
#ifndef _WIN32
    if(fd_>=0)
    {
        const char bye[] = "{\"bye\":1}\n";
        (void)!::send(fd_, bye, sizeof(bye)-1, 0);
        ::close(fd_);
    }
#endif
}

bool sixdof_brain_link::connect_to(int port, double wait, std::string &msg)
{
#ifdef _WIN32
    msg = "the external brain (X 191) needs POSIX sockets";
    return false;
#else
    const auto t0 = std::chrono::steady_clock::now();

    while(true)
    {
        fd_ = ::socket(AF_INET, SOCK_STREAM, 0);

        if(fd_<0)
        {
            msg = "no socket for the external brain (X 191)";
            return false;
        }

        sockaddr_in a;
        std::memset(&a, 0, sizeof(a));
        a.sin_family = AF_INET;
        a.sin_port = htons(port);
        a.sin_addr.s_addr = htonl(INADDR_LOOPBACK);

        if(::connect(fd_, (sockaddr*)&a, sizeof(a))==0)
        {
            int one = 1;
            ::setsockopt(fd_, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one));
            return true;
        }

        ::close(fd_);
        fd_ = -1;

        if(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()>wait)
        {
            msg = "no external brain (X 191) on 127.0.0.1:" + std::to_string(port);
            return false;
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(250));
    }
#endif
}

bool sixdof_brain_link::exchange(const std::string &out, std::string &in, double timeout)
{
#ifdef _WIN32
    return false;
#else
    if(fd_<0)
    return false;

    const std::string line = out + "\n";
    size_t sent = 0;

    while(sent<line.size())
    {
        const ssize_t n = ::send(fd_, line.data()+sent, line.size()-sent, 0);
        if(n<=0)
        return false;
        sent += size_t(n);
    }

    const auto t0 = std::chrono::steady_clock::now();

    while(true)
    {
        const size_t e = buf_.find('\n');

        if(e!=std::string::npos)
        {
            in = buf_.substr(0, e);
            buf_.erase(0, e+1);
            return true;
        }

        const double left = timeout - std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if(left<=0.0)
        return false;

        fd_set rs;
        FD_ZERO(&rs);
        FD_SET(fd_, &rs);
        timeval tv;
        tv.tv_sec = long(left);
        tv.tv_usec = long((left - double(tv.tv_sec))*1.0e6);

        if(::select(fd_+1, &rs, nullptr, nullptr, &tv)<=0)
        return false;

        char b[4096];
        const ssize_t n = ::recv(fd_, b, sizeof(b), 0);
        if(n<=0)
        return false;
        buf_.append(b, size_t(n));
    }
#endif
}

bool sixdof_brain_link::has(const std::string &json, const char *key)
{
    return json.find(std::string("\"") + key + "\"")!=std::string::npos;
}

double sixdof_brain_link::number(const std::string &json, const char *key, double def)
{
    const std::string k = std::string("\"") + key + "\"";
    size_t q = json.find(k);
    if(q==std::string::npos)
    return def;
    q = json.find(':', q + k.size());
    if(q==std::string::npos)
    return def;
    const char *s = json.c_str() + q + 1;
    char *end = nullptr;
    const double v = std::strtod(s, &end);
    return (end==s) ? def : v;
}
