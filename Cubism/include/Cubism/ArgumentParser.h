/*
 *  ArgumentParser.h
 *  Cubism
 *
 *	This argument parser assumes that all arguments are optional ie, each of the argument names is preceded by a '-'
 *		all arguments are however NOT optional to avoid a mess with default values and returned values when not found!
 *
 *	More converter could be required:
 *		add as needed
 *			TypeName as{TypeName}() in Value
 *
 *  Created by Christian Conti on 6/7/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <iosfwd>  // Forward declaration of <iostream>
#include <map>
#include <string>

#include "Common.h"

CUBISM_NAMESPACE_BEGIN

class Value
{
private:
    std::string content;

public:
    Value() = default;
    Value(const std::string &content_) : content(content_) {}
    Value(const Value& c) = default;

    Value& operator=(const Value& rhs)
    {
        if (this != &rhs)
            content = rhs.content;
        return *this;
    }
    Value& operator+=(const Value& rhs)
    {
        content += " " + rhs.content;
        return *this;
    }
    Value operator+(const Value& rhs) { return Value(content + " " + rhs.content); }

    double asDouble(double def = 0);
    int asInt(int def = 0);
    bool asBool(bool def = false);
    std::string asString(const std::string &def = std::string());
    friend std::ostream& operator<<(std::ostream& lhs, const Value& rhs);
};


class CommandlineParser
{
private:
    const int iArgC;
    char** vArgV;
    bool bStrictMode, bVerbose;

protected:
    std::map<std::string,Value> mapArguments;

public:
    CommandlineParser(int argc, char ** argv);

    Value& operator()(std::string key);
    bool check(std::string key) const;

    int getargc() const { return iArgC; }
    char** getargv() const { return vArgV; }

    void set_strict_mode()
    {
        bStrictMode = true;
    }

    void unset_strict_mode()
    {
        bStrictMode = false;
    }

    void mute()
    {
        bVerbose = false;
    }

    void loud()
    {
        bVerbose = true;
    }

    void save_options(const std::string &path = ".");
    void print_args();
};


class ArgumentParser : public CommandlineParser
{
    typedef std::map<std::string, Value> ArgMap;
    typedef std::map<std::string, Value*> pArgMap;
    typedef std::map<std::string, ArgMap* > FileMap;

    const char commentStart;

    // keep a reference from option origin
    ArgMap  from_commandline;
    FileMap from_files;
    pArgMap from_code;

    // for runtime interaction (we keep the original map)
    ArgMap mapRuntime;

    // helper
    void _ignoreComments(std::istream& stream, char commentChar);
    void _parseFile(std::ifstream& stream, ArgMap& container);

public:
    ArgumentParser(const int _argc, char ** _argv, const char cstart='#'):
        CommandlineParser(_argc, _argv), commentStart(cstart)
    {
        from_commandline = mapArguments;
    }

    virtual ~ArgumentParser()
    {
        for (FileMap::iterator it = from_files.begin(); it != from_files.end(); it++)
            delete it->second;
    }

    void readFile(const std::string &filepath);
    Value& operator()(std::string key);

    inline bool exist(const std::string &key) const { return check(key); }

    void write_runtime_environment() const;
    void read_runtime_environment();

    Value& parseRuntime(std::string key);
    void print_args(void);
};

CUBISM_NAMESPACE_END
