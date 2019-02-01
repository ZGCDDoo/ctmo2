#pragma once

#include <iostream>
#include <regex>
#include <tuple>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace CMDParser
{

class CMDInfo
{

  public:
    friend CMDInfo GetProgramOptions(int argc, char **argv);

    CMDInfo(){};

    CMDInfo(const CMDInfo &cmdInfo) = default;

    CMDInfo(const std::string &prefixIn, const int &iterIn, const std::string &suffixIn, const bool &doSCIn = false,
            const bool &exitFromCMDIn = false) : fnamePrefix_(prefixIn),
                                                 iter_(iterIn),
                                                 fnameSuffix_(suffixIn),
                                                 doSC_(doSCIn),
                                                 exitFromCMD_(exitFromCMDIn) {}

    std::string fileName() const
    {
        return (fnamePrefix_ + std::to_string(iter_) + fnameSuffix_);
    }

    //setters
    std::string fnamePrefix() const { return fnamePrefix_; }
    int iter() const { return iter_; }
    std::string fnameSuffix() const { return fnameSuffix_; }
    bool doSC() const { return doSC_; }
    bool exitFromCMD() const { return exitFromCMD_; }

  private:
    std::string fnamePrefix_{""};
    int iter_{0};
    std::string fnameSuffix_{""};
    bool doSC_{true};
    bool exitFromCMD_{false};
};

CMDInfo GetProgramOptions(int argc, char **argv)
{

    // Define and parse the program options

    namespace po = boost::program_options;
    po::options_description desc("Example usage: ctmo params1.json. \n\nAllowed Options:");
    desc.add_options()("help,h", "Print help messages.")("fname,f", po::value<std::string>()->required(), "simulation filename (in json format).")("no-sc,n", "Don't perform the selfconsistency nor prepare the next iteration.");

    po::positional_options_description positional;
    positional.add("fname", -1);

    po::variables_map vm;
    CMDInfo cmdInfo;

    try
    {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(positional).run(), vm);

        if (vm.count("help"))
        {
            std::cout << desc << std::endl;
            cmdInfo.exitFromCMD_ = true;
        }

        po::notify(vm); // throws on error, so do after help in case
                        // there are any problems
        cmdInfo.doSC_ = !vm.count("no-sc");
    }
    catch (po::error &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << desc << std::endl;
        cmdInfo.exitFromCMD_ = true;
    }

    const std::string jsonFileName = vm["fname"].as<std::string>();
    std::cout << "jsonFileName = " << jsonFileName << std::endl;

    std::regex numberRegex("\\d+");

    std::smatch numberMatch;
    std::string prefix;
    int iter = -999;
    std::string suffix;
    if (std::regex_search(jsonFileName, numberMatch, numberRegex))
    {
        std::cout << "Prefix: '" << numberMatch.prefix() << "'\n";
        prefix = numberMatch.prefix().str();
        std::cout << "Suffix: '" << numberMatch.suffix() << "\'\n\n";
        iter = std::stoi(numberMatch[0]);
        std::cout << "iter = " << iter << std::endl;
        suffix = numberMatch.suffix();
    }

    CMDInfo cmdInfoResult(prefix, iter, suffix, cmdInfo.doSC(), cmdInfo.exitFromCMD());

    std::cout << cmdInfoResult.fileName() << std::endl;

    // check if the params filename exists, else exit:
    using boost::filesystem::exists;
    if (!exists(cmdInfoResult.fileName()))
    {
        std::cerr << "File: " << cmdInfoResult.fileName() << " does not exist." << std::endl;
        std::cerr << "Are there multiple numbers in your filename or bizarre filename ?" << std::endl;
        throw std::runtime_error("Aborting because of bad file name!");
    }

    return cmdInfoResult;
}

} // namespace CMDParser