#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex> // C++11

#include "Cubism/ArgumentParser.h"

CUBISM_NAMESPACE_BEGIN

///////////////////////////////////////////////////////////
// Value
///////////////////////////////////////////////////////////
double Value::asDouble(double def)
{
    if (content == "")
    {
        std::ostringstream sbuf;
        sbuf << def;
        content = sbuf.str();
    }
    return (double) atof(content.c_str());
}

int Value::asInt(int def)
{
    if (content == "")
    {
        std::ostringstream sbuf;
        sbuf << def;
        content = sbuf.str();
    }
    return atoi(content.c_str());
}

bool Value::asBool(bool def)
{
    if (content == "")
    {
        if (def) content = "true";
        else     content = "false";
    }
    if (content == "0") return false;
    if (content == "false") return false;

    return true;
}

std::string Value::asString(const std::string &def)
{
    if (content == "") content = def;

    return content;
}

std::ostream& operator<<(std::ostream& lhs, const Value& rhs)
{
    lhs << rhs.content;
    return lhs;
}


///////////////////////////////////////////////////////////
// CommandlineParser
///////////////////////////////////////////////////////////
static inline void _normalizeKey(std::string& key)
{
    if (key[0] == '-') key.erase(0,1);
    if (key[0] == '+') key.erase(0,1);
}

static inline bool _existKey(const std::string& key,
                             const std::map<std::string, Value>& container)
{
    return container.find(key) != container.end();
}

Value& CommandlineParser::operator()(std::string key)
{
    _normalizeKey(key);
    if (bStrictMode)
    {
        if (!_existKey(key,mapArguments))
        {
            printf("Runtime option NOT SPECIFIED! ABORTING! name: %s\n",key.data());
            abort();
        }
    }

    if (bVerbose) printf("%s is %s\n", key.data(), mapArguments[key].asString().data());
    return mapArguments[key];
}

bool CommandlineParser::check(std::string key) const
{
    _normalizeKey(key);
    return _existKey(key,mapArguments);
}

CommandlineParser::CommandlineParser(const int argc, char **argv)
    : iArgC(argc), vArgV(argv), bStrictMode(false), bVerbose(true)
{
    // parse commandline <key> <value> pairs.  Key passed on the command
    // line must start with a leading dash (-). For example:
    // -mykey myvalue0 [myvalue1 ...]
    for (int i=1; i<argc; i++)
        if (argv[i][0] == '-')
        {
            std::string values = "";
            int itemCount = 0;

            // check if the current key i is a list of values. If yes,
            // concatenate them into a string
            for (int j=i+1; j<argc; j++)
            {
                // if the current value is numeric and (possibly) negative,
                // do not interpret it as a key
                const bool leadingDash = (argv[j][0] == '-');
                const bool isNumeric = std::regex_match(argv[j], std::regex("(\\+|-)?[0-9]*(\\.[0-9]*)?((e|E)(\\+|-)?[0-9]*)?"));
                if (leadingDash && !isNumeric)
                    break;
                else
                {
                    if (strcmp(values.c_str(), ""))
                        values += ' ';

                    values += argv[j];
                    itemCount++;
                }
            }

            if (itemCount == 0)
                values = "true";

            std::string key(argv[i]);
            key.erase(0,1); // remove leading '-'
            if (key[0] == '+') // for key concatenation
            {
                key.erase(0,1);
                if (!_existKey(key,mapArguments))
                    mapArguments[key] = Value(values); // skip leading white space
                else
                    mapArguments[key] += Value(values);
            }
            else // regular key
            {
                if (!_existKey(key,mapArguments))
                    mapArguments[key] = Value(values);
            }

            i += itemCount;
        }

    mute();
    //printf("found %ld arguments of %d\n",mapArguments.size(),argc);
}

void CommandlineParser::save_options(const std::string &path)
{
    std::string options;
    for(std::map<std::string,Value>::iterator it=mapArguments.begin(); it!=mapArguments.end(); it++)
    {
        options+= it->first + " " + it->second.asString() + " ";
    }
    std::string filepath = path + "/argumentparser.log";
    FILE * f = fopen(filepath.data(), "a");
    if (f == NULL)
    {
        fprintf(stderr, "impossible to write %s.\n", filepath.data());
        return;
    }
    fprintf(f, "%s\n", options.data());
    fclose(f);
}

void CommandlineParser::print_args()
{
    for(std::map<std::string,Value>::iterator it=mapArguments.begin(); it!=mapArguments.end(); it++)
    {
        std::cout.width(50);
        std::cout.fill('.');
        std::cout << std::left << it->first;
        std::cout << ": " << it->second.asString() << std::endl;
    }
}


///////////////////////////////////////////////////////////
// ArgumentParser
///////////////////////////////////////////////////////////
void ArgumentParser::_ignoreComments(std::istream& stream, const char commentChar)
{
    stream >> std::ws;
    int nextchar = stream.peek();
    while (nextchar == commentChar)
    {
        stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        stream >> std::ws;
        nextchar = stream.peek();
    }
}

void ArgumentParser::_parseFile(std::ifstream& stream, ArgMap& container)
{
    // read (key value) pairs from input file, ignore comments
    // beginning with commentStart
    _ignoreComments(stream, commentStart);
    while (!stream.eof())
    {
        std::string line, key, val;
        std::getline(stream, line);
        std::istringstream lineStream(line);
        lineStream >> key;
        lineStream >> val;
        _ignoreComments(lineStream, commentStart);
        while(!lineStream.eof())
        {
            std::string multiVal;
            lineStream >> multiVal;
            val += (" " + multiVal);
            _ignoreComments(lineStream, commentStart);
        }

        const Value V(val);
        if (key[0] == '-')
            key.erase(0,1);

        if (key[0] == '+')
        {
            key.erase(0,1);
            if (!_existKey(key,container)) // skip leading white space
                container[key] = V;
            else
                container[key] += V;
        }
        else if (!_existKey(key,container))
            container[key] = V;
        _ignoreComments(stream, commentStart);
    }
}

void ArgumentParser::readFile(const std::string &filepath)
{
    from_files[filepath] = new ArgMap;
    ArgMap& myFMap = *(from_files[filepath]);

    std::ifstream confFile(filepath.c_str());
    if (confFile.good())
    {
        _parseFile(confFile, mapArguments);
        confFile.clear();
        confFile.seekg(0, std::ios::beg);
        _parseFile(confFile, myFMap); // we keep a reference for each separate file read
    }
    confFile.close();
}

Value& ArgumentParser::operator()(std::string key)
{
    _normalizeKey(key);
    const bool bDefaultInCode = !_existKey(key,mapArguments);
    Value& retval = CommandlineParser::operator()(key);
    if (bDefaultInCode) from_code[key] = &retval;
    return retval;
}

void ArgumentParser::write_runtime_environment() const
{
    time_t rawtime;
    std::time(&rawtime);
    struct tm* timeinfo = std::localtime(&rawtime);
    char buf[256];
    std::strftime(buf, 256, "%A, %h %d %Y, %r", timeinfo);

    std::ofstream runtime("runtime_environment.conf");
    runtime << commentStart << " RUNTIME ENVIRONMENT SETTINGS" << std::endl;
    runtime << commentStart << " ============================" << std::endl;
    runtime << commentStart << " " << buf << std::endl;
    runtime << commentStart << " Use this file to set runtime parameter interactively." << std::endl;
    runtime << commentStart << " The parameter are read every \"refreshperiod\" steps." << std::endl;
    runtime << commentStart << " When editing this file, you may use comments and string concatenation." << std::endl;
    runtime << commentStart << " The simulation can be terminated without killing it by setting \"exit\" to true." << std::endl;
    runtime << commentStart << " (This will write a serialized restart state. Set \"exitsave\" to false if not desired.)" << std::endl;
    runtime << commentStart << std::endl;
    runtime << commentStart << " !!! WARNING !!! EDITING THIS FILE CAN POTENTIALLY CRASH YOUR SIMULATION !!! WARNING !!!" << std::endl;
    for (typename std::map<std::string,Value>::const_iterator it = mapArguments.begin(); it != mapArguments.end(); ++it)
        runtime << it->first << '\t' << it->second << std::endl;
}

void ArgumentParser::read_runtime_environment()
{
    mapRuntime.clear();
    std::ifstream runtime("runtime_environment.conf");
    if (runtime.good())
        _parseFile(runtime, mapRuntime);
    runtime.close();
}

Value& ArgumentParser::parseRuntime(std::string key)
{
    _normalizeKey(key);
    if (!_existKey(key,mapRuntime))
    {
        printf("ERROR: Runtime parsing for key %s NOT FOUND!! Check your runtime_environment.conf file\n",key.data());
        abort();
    }
    return mapRuntime[key];
}

void ArgumentParser::print_args()
{
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "* Summary:" << std::endl;
    std::cout << "*    Parameter read from command line:                " << from_commandline.size() << std::endl;
    size_t nFiles = 0;
    size_t nFileParameter = 0;
    for (FileMap::const_iterator it=from_files.begin(); it!=from_files.end(); ++it)
    {
        if (it->second->size() > 0)
        {
            ++nFiles;
            nFileParameter += it->second->size();
        }
    }
    std::cout << "*    Parameter read from " << std::setw(3) << std::right << nFiles << " file(s):                 " << nFileParameter << std::endl;
    std::cout << "*    Parameter read from defaults in code:            " << from_code.size() << std::endl;
    std::cout << "*    Total number of parameter read from all sources: " << mapArguments.size() << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    // command line given arguments
    if (!from_commandline.empty())
    {
        std::cout << "* Command Line:" << std::endl;
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        for(ArgMap::iterator it=from_commandline.begin(); it!=from_commandline.end(); it++)
        {
            std::cout.width(50);
            std::cout.fill('.');
            std::cout << std::left << it->first;
            std::cout << ": " << it->second.asString() << std::endl;
        }
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    }

    // options read from input files
    if (!from_files.empty())
    {
        for (FileMap::iterator itFile=from_files.begin(); itFile!=from_files.end(); itFile++)
        {
            if (!itFile->second->empty())
            {
                std::cout << "* File: " << itFile->first << std::endl;
                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
                ArgMap& fileArgs = *(itFile->second);
                for(ArgMap::iterator it=fileArgs.begin(); it!=fileArgs.end(); it++)
                {
                    std::cout.width(50);
                    std::cout.fill('.');
                    std::cout << std::left << it->first;
                    std::cout << ": " << it->second.asString() << std::endl;
                }
                std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            }
        }
    }

    // defaults defined in code
    if (!from_code.empty())
    {
        std::cout << "* Defaults in Code:" << std::endl;
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        for(pArgMap::iterator it=from_code.begin(); it!=from_code.end(); it++)
        {
            std::cout.width(50);
            std::cout.fill('.');
            std::cout << std::left << it->first;
            std::cout << ": " << it->second->asString() << std::endl;
        }
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    }
}

CUBISM_NAMESPACE_END
