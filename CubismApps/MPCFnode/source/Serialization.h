/*
 *  Serialization.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 10/07/2017
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef SERIALIZATION_H_UQGXZSUK
#define SERIALIZATION_H_UQGXZSUK

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>

#include <Cubism/ArgumentParser.h>
#include "Types.h"
using namespace cubism;

class SerializationMeta
{
public:
    SerializationMeta(int& ia, int& ib, int& ic, int& id, Real& fa, Real& fb, Real& fc) :
        restart_id(ia), step_id(ib), step_lastdump(ic), dumpcount(id),
        t(fa), t_nextdump(fb), t_lastdump(fc)
    {}

    int& restart_id;
    int& step_id;
    int& step_lastdump;
    int& dumpcount;

    Real& t;
    Real& t_nextdump;
    Real& t_lastdump;
};

#define __SERIALIZER_META_DATA(__SEP__) \
    _m.restart_id __SEP__ _m.t __SEP__ _m.step_id __SEP__ _m.t_nextdump __SEP__ _m.step_lastdump __SEP__ _m.t_lastdump __SEP__ _m.dumpcount

template <typename TGrid>
class Serialization
{
public:
    Serialization(ArgumentParser& p) :
        m_parser(p),
        m_serializer(NULL),
        m_deserializer(NULL)
    {}

    // serializer / deserializer interfaces
    typedef void (*TSerializer)(const TGrid&, const int, const Real, const std::string, const std::string, const bool);
    typedef void (*TDeserializer)(TGrid&, const std::string, const std::string);

    inline void set_serializer(TSerializer s) { assert(s != NULL); m_serializer = s; }
    inline void set_deserializer(TDeserializer d) { assert(d != NULL); m_deserializer = d; }

    void serialize(const TGrid& grid, const SerializationMeta& meta, const bool verbose=true)
    {
        if (m_serializer)
        {
            const std::string path = m_parser("-fpath").asString(".");
            std::stringstream datafile;
            datafile << "restart_data_" << meta.restart_id;

            m_serializer(grid, meta.step_id, meta.t, datafile.str().c_str(), path.c_str(), false);

            if (verbose)
                std::cout << "SERIALIZATION: written file " << path << "/" << datafile.str() << std::endl;

            if (verbose)
            {
                std::stringstream statusfile;
                statusfile << path + "/restart_HEAD";
                std::ofstream status_0(statusfile.str().c_str());

                statusfile.str("");
                statusfile.clear();
                statusfile << path + "/restart_status_" << meta.restart_id;
                std::ofstream status_1(statusfile.str().c_str());

                _metadata(status_0, meta); // HEAD
                _metadata(status_1, meta); // this one

                std::cout << "SERIALIZATION: written meta data " << statusfile.str() << std::endl;
            }

            // if specified, remove previous state
            if (verbose && meta.restart_id > 1)
            {
                const bool removeold = m_parser("-removeold").asBool(false);
                const bool keep_meta = m_parser("-removeold_keep_meta").asBool(false);
                const std::string suffix = m_parser("restart_format").asString("zbin");
                if (removeold)
                {
                    std::cout << "SERIALIZATION: removing previous state... " << std::endl;

                    if (!keep_meta)
                    {
                        // remove old status file
                        std::stringstream statusfile;
                        statusfile << path + "/restart_status_" << meta.restart_id - 1;
                        std::remove(statusfile.str().c_str());
                    }

                    // remove old data file
                    std::stringstream datafile;
                    datafile << path + "/restart_data_" << meta.restart_id - 1 << "." << suffix;
                    std::remove(datafile.str().c_str());
                }
            }
        }
        else
        {
            if (verbose)
                std::cerr << "ERROR: No serializer has been assigned!" << std::endl;
            abort();
        }
    }

    void deserialize(TGrid& grid, SerializationMeta& meta, const bool verbose=true)
    {
        if (m_deserializer)
        {
            const std::string path = m_parser("fpath").asString(".");
            if (!(check_state("zbin", path, verbose) || check_state("h5", path, verbose)))
            {
                if (verbose)
                    std::cerr << "ERROR: Inconsistent serialization state!" << std::endl;
                abort();
            }

            const std::string fmeta = path + "/restart_HEAD";
            std::ifstream head(fmeta.c_str());
            _metadata(head, meta);
            assert(meta.restart_id >= 0);
            assert(meta.t >= 0);
            assert(meta.step_id >= 0);
            assert(meta.t_nextdump >= 0);
            assert(meta.step_lastdump >= 0);
            assert(meta.t_lastdump >= 0);

            std::stringstream datafile;
            datafile << "restart_data_" << meta.restart_id;
            m_deserializer(grid, datafile.str().c_str(), path.c_str());
            if (verbose)
            {
                std::cout << "DESERIALIZATION META DATA: restart_status_" << meta.restart_id << std::endl;
                std::cout << "DESERIALIZATION: Time is " << std::scientific << meta.t << " and step_id is " << meta.step_id << std::endl;
            }
        }
        else
        {
            if (verbose)
                std::cerr << "ERROR: No deserializer has been assigned!" << std::endl;
            abort();
        }
    }

    void stream_to_file(const TGrid& grid, const int step, const Real time, std::string filename, const bool verbose=true)
    {
        if (m_serializer)
        {
            const std::string path = m_parser("-fpath").asString(".");
            m_serializer(grid, step, time, filename.c_str(), path.c_str(), true);

            if (verbose)
                std::cout << "STREAM TO FILE: written file " << path << "/" << filename << std::endl;
        }
        else
        {
            if (verbose)
                std::cerr << "ERROR: No serializer has been assigned!" << std::endl;
            abort();
        }
    }

    void stream_from_file(TGrid& grid, std::string filename, const bool verbose=true)
    {
        if (m_deserializer)
        {
            const std::string path = m_parser("fpath").asString(".");
            m_deserializer(grid, filename.c_str(), path.c_str());

            if (verbose)
                std::cout << "STREAM FROM FILE: read file " << path << "/" << filename << std::endl;
        }
        else
        {
            if (verbose)
                std::cerr << "ERROR: No deserializer has been assigned!" << std::endl;
            abort();
        }
    }

    bool check_state(const std::string suffix, const std::string& path, const bool verbose=true)
    {
        // check for existence of most recent restart state for given suffix
        const std::string fmeta = path + "/restart_HEAD";
        std::ifstream head(fmeta.c_str());
        bool bstatus = head.good();
        if (bstatus)
        {
            int pair_id;
            head >> pair_id;
            head.close();
            std::ostringstream datafile;
            datafile << path << "/restart_data_" << pair_id << "." << suffix;
            std::ifstream datastatus(datafile.str().c_str());
            bstatus = datastatus.good();
            datastatus.close();
        }
        else
        {
            if (verbose)
                std::cerr << "ERROR: Can not find '" << fmeta << "'" << std::endl;
            abort();
        }
        return bstatus;
    }

private:
    ArgumentParser& m_parser;
    TSerializer m_serializer;
    TDeserializer m_deserializer;

    inline void _metadata(std::ofstream& _of, const SerializationMeta& _m)
    {
        _of.setf(std::ios::scientific, std::ios::floatfield);
        _of << __SERIALIZER_META_DATA(<< '\t' <<);
    }

    inline void _metadata(std::ifstream& _if, SerializationMeta& _m)
    {
        _if >> __SERIALIZER_META_DATA(>>);
    }
};

#undef __SERIALIZER_META_DATA

#endif /* SERIALIZATION_H_UQGXZSUK */
