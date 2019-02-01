#pragma once

#include "../../deps/nlohmann_json/json.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <cstdlib>
#include <ctime>
#include "Logging.hpp"
#include "CMDParser.hpp"
#include "Conventions.hpp"

namespace IO
{
namespace FS
{

void WriteToFile(const size_t &iter, double &value, const std::string &fname)
{
    std::ofstream fout(fname, std::ios_base::out | std::ios_base::app);
    fout << iter << " " << value << std::endl;
    fout.close();
}

void WriteToFile(const size_t &iter, const std::vector<double> &stats, const std::string &fname)
{
    assert(stats.size() == 2); //mean and stddev
    std::ofstream fout(fname, std::ios_base::out | std::ios_base::app);
    fout << iter << " " << stats[0] << " " << stats[1] << std::endl;
    fout.close();
}

size_t CalculateNextSeed()
{
    std::srand(std::time(nullptr));
    for (int i = 0; i < std::rand() % 100; i++)
    {
        std::rand();
    }

    Utilities::EngineTypeFibonacci3217_t randFib(std::rand());
    Utilities::UniformRngFibonacci3217_t urandFib(randFib, Utilities::UniformDistribution_t(0.0, 1.0));

    for (int i = 0; i < std::rand() % 100; i++)
    {
        urandFib();
    }

    double nextSeed = urandFib() * double(std::rand());
    return (static_cast<size_t>(nextSeed));
}

void PrepareNextIter(const CMDParser::CMDInfo &cmdInfo)
{
    using boost::filesystem::copy_file;
    using boost::filesystem::exists;
    using boost::filesystem::path;

    auto nameCon = Conventions::BuildFileNameConventions();
    const std::string datExt = nameCon.at("datExt");
    const size_t iter = cmdInfo.iter();
    const std::string iterStr = std::to_string(iter);

    const path path_hyb(nameCon.at("hybUpFile"));
    copy_file(nameCon.at("hybUpFile"), path_hyb.stem().string() + std::to_string(iter + 1) + path_hyb.extension().string());
    const std::vector<std::string> filePaths = {nameCon.at("greenUpFile"),
                                                nameCon.at("selfUpFile")};
    for (const auto &filePath : filePaths)
    {
        const path path_tmp(filePath);
        copy_file(filePath, path_tmp.stem().string() + iterStr + path_tmp.extension().string());
    }

    if (Logging::LevelIsTrace())
    {
        const std::vector<std::string> filePathsTrace = {nameCon.at("obsJsonFile"),
                                                         nameCon.at("gtauFile")};
        for (const auto &filePath : filePathsTrace)
        {
            if (exists(filePath))
            {
                const path path_tmp(filePath);
                copy_file(filePath, path_tmp.stem().string() + iterStr + path_tmp.extension().string());
            }
        }
    }

#ifdef AFM
    const path path_hybDown(nameCon.at("hybDownFile"));
    copy_file(nameCon.at("hybDownFile"), path_hybDown.stem().string() + std::to_string(iter + 1) + path_hybDown.extension().string());
    const std::vector<std::string> filePathsDown = {nameCon.at("greenDownFile"),
                                                    nameCon.at("selfDownFile")};
    for (const auto &filePath : filePathsDown)
    {
        const path path_tmp(filePath);
        copy_file(filePath, path_tmp.stem().string() + iterStr + path_tmp.extension().string());
    }
#endif

    const std::string fname = cmdInfo.fileName();
    std::ifstream fin(fname);
    Json params;
    fin >> params;
    fin.close();

    Json results;
    fin.open(nameCon.at("obsJsonFile"));
    fin >> results;
    fin.close();

    for (Json::iterator it = results.begin(); it != results.end(); ++it)
    {
        std::vector<double> stats = it.value();
        WriteToFile(iter, stats, it.key() + datExt);
    }

    size_t nextSeed = CalculateNextSeed();
    params["monteCarlo"]["seed"] = static_cast<size_t>(nextSeed);

    params["model"]["hybUpFile"] = std::string("hybUp") + std::to_string(iter + 1) + datExt;

#ifdef AFM
    params["model"]["hybDownFile"] = std::string("hybDown") + std::to_string(iter + 1) + datExt;
#endif

    if (params["model"].find("n") != params["model"].end())
    {
        double nParams = params["model"]["n"];
        double nResult = results["n"].at(0); //mean of n from simulation
        params["model"]["mu"] = double(params["model"]["mu"]) - double(params["solver"]["S"]) * (nResult - nParams);
    }
    std::ofstream fout(cmdInfo.fnamePrefix() + std::to_string(iter + 1) + cmdInfo.fnameSuffix());
    fout << std::setw(4) << params << std::endl;
    fout.close();

    return;
}

} // namespace FS
} // namespace IO