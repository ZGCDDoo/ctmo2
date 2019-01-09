#pragma once

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "MPITools.hpp"

namespace Logging
{

const std::string ROOT = "ROOT";

void InitFile(const Json &jjLog, const std::string &loggerName = ROOT)
{
    const std::string logLevel = jjLog["level"].get<std::string>();
    const std::string file_sink = jjLog["file"].get<std::string>();

    if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
    {
        auto logger_ = spdlog::basic_logger_st(loggerName, file_sink);

        if (logLevel == "trace")
        {
            logger_->set_level(spdlog::level::trace);
        }
        else if (logLevel == "debug")
        {
            logger_->set_level(spdlog::level::debug);
        }
        else if (logLevel == "info")
        {
            logger_->set_level(spdlog::level::info);
        }
        else if (logLevel == "warn")
        {
            logger_->set_level(spdlog::level::warn);
        }
        else
        {
            logger_->set_level(spdlog::level::critical);
        }

        const std::string msg = "Logger " + loggerName + " initialized";
        logger_->info(msg);
    }
}

void InitStdout(const Json &jjLog, const std::string loggerName = ROOT)
{
    const std::string logLevel = jjLog["level"].get<std::string>();

    if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
    {
        auto logger_ = spdlog::stdout_color_st(loggerName);

        if (logLevel == "trace")
        {
            logger_->set_level(spdlog::level::trace);
        }
        else if (logLevel == "debug")
        {
            logger_->set_level(spdlog::level::debug);
        }
        else if (logLevel == "info")
        {
            logger_->set_level(spdlog::level::info);
        }
        else if (logLevel == "warn")
        {
            logger_->set_level(spdlog::level::warn);
        }
        else
        {
            logger_->set_level(spdlog::level::critical);
        }

        const std::string msg = "Logger " + loggerName + " initialized";
        logger_->info(msg);
    }
}

void Init(const Json &jjLog, const std::string &loggerName = ROOT)
{
    const bool logToFile = jjLog["logToFile"].get<bool>();

    if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
    {
        if (logToFile)
        {
            InitFile(jjLog, loggerName);
        }
        else
        {
            InitStdout(jjLog, loggerName);
        }
    }
}

void Trace(const std::string &msg, const std::string &loggerName = ROOT)
{
    // #ifdef HAVEMPI
    //     mpi::environment env;
    //     mpi::communicator world;
    // #endif
    if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
    {
        auto logger = spdlog::get(loggerName);
        if (!logger)
        {
            logger = spdlog::stdout_color_st(loggerName);
            logger->set_level(spdlog::level::trace);
            logger = spdlog::get(loggerName);
        }
        logger->trace(msg);
    }
}

void Debug(const std::string &msg, const std::string &loggerName = ROOT)
{
    // #ifdef HAVEMPI
    //     mpi::environment env;
    //     mpi::communicator world;
    // #endif
    if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
    {
        auto logger = spdlog::get(loggerName);
        if (!logger)
        {
            logger = spdlog::stdout_color_st(loggerName);
            logger->set_level(spdlog::level::trace);
            logger = spdlog::get(loggerName);
        }
        logger->debug(msg);
    }
}

// template <typename... TArgs>
void Info(const std::string &msg, const std::string &loggerName = ROOT)
{
    // #ifdef HAVEMPI
    // mpi::environment env;
    // mpi::communicator world;
    // #endif
    if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
    {
        auto logger = spdlog::get(loggerName);

        if (!logger)
        {
            logger = spdlog::stdout_color_st(loggerName);
            logger->set_level(spdlog::level::trace);
            logger = spdlog::get(loggerName);
        }

        logger = spdlog::get(loggerName);
        logger->info(msg);
    }
}

void Warn(const std::string &msg, const std::string &loggerName = ROOT)
{
    // #ifdef HAVEMPI
    //     mpi::environment env;
    //     mpi::communicator world;
    // #endif
    if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
    {
        auto logger = spdlog::get(loggerName);

        if (!logger)
        {
            logger = spdlog::stdout_color_st(loggerName);
            logger->set_level(spdlog::level::trace);
            logger = spdlog::get(loggerName);
        }

        logger = spdlog::get(loggerName);
        logger->warn(msg);
    }
}

void Critical(const std::string &msg, const std::string &loggerName = ROOT)
{
    if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
    {
        auto logger = spdlog::get(loggerName);

        if (!logger)
        {
            logger = spdlog::stdout_color_st(loggerName);
            logger->set_level(spdlog::level::trace);
            logger = spdlog::get(loggerName);
        }

        logger = spdlog::get(loggerName);
        logger->critical(msg);
    }
}

bool LevelIsTrace(const std::string &loggerName = ROOT)
{
    auto logger = spdlog::get(loggerName);
    return (spdlog::level::level_enum::trace == logger->level());
}

} // namespace Logging
