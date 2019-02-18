// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givmemorer.h
// Copyright(c)'2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas
// Time-stamp: <21 May 12 15:35:16 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

#include "givaro/givconfig.h"
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

namespace Givaro {


    struct Memorer {
        void lookup() { this->process_mem_usage(); }
        double usage() const { return vm_usage; }
        double resident() const { return resident_set; }
        friend std::ostream& operator<< (std::ostream& out, const Memorer& M) {
            return out << M.vm_usage << " (" << M.resident_set << ')';
        }


    private:
        double vm_usage, resident_set;

        //////////////////////////////////////////////////////////////////////////////
        //
        // process_mem_usage(double &, double &) - takes two doubles by reference,
        // attempts to read the system-dependent data for a process' virtual memory
        // size and resident set size, and return the results in KB.
        //
        // On failure, returns 0.0, 0.0
        void process_mem_usage() {
            using std::ios_base;
            using std::ifstream;
            using std::string;

            vm_usage     = 0.0;
            resident_set = 0.0;

            // 'file' stat seems to give the most reliable results
            //
            ifstream stat_stream("/proc/self/stat",ios_base::in);

            // dummy vars for leading entries in stat that we don't care about
            //
            string pid, comm, state, ppid, pgrp, session, tty_nr;
            string tpgid, flags, minflt, cminflt, majflt, cmajflt;
            string utime, stime, cutime, cstime, priority, nice;
            string O, itrealvalue, starttime;

            // the two fields we want
            //
            unsigned long vsize;
            long rss;

            stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
            >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
            >> utime >> stime >> cutime >> cstime >> priority >> nice
            >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

            stat_stream.close();

            long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
            vm_usage     = vsize / 1024.0;
            resident_set = rss * page_size_kb;
        }
    };


} // namespace Givaro

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
