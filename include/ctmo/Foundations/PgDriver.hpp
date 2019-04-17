//
// Created by charles-david on 14/04/19.
//
#pragma once
#ifndef CTMO_PGDRIVER_HPP
#define CTMO_PGDRIVER_HPP


#include <pqxx/pqxx>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <chrono>
#include <ctime>


#define POSTGRES_CONNECT_TIME_DELAY 60

namespace DB
{
struct Timer
{
    Timer() = default;

    void Start(double duration)
    {
        duration_ = duration;
        start_ = std::chrono::steady_clock::now();
    };

    static void PrintTime()
    {
        if (mpiUt::Tools::Rank() == mpiUt::Tools::master)
        {
            auto timeChrono = std::chrono::system_clock::now();
            std::time_t timeNow = std::chrono::system_clock::to_time_t(timeChrono);
            std::cout << "\t " << std::ctime(&timeNow) << std::endl;
        }
    }

    bool End()
    {
        return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_).count() >
               duration_;
    };

private:
    double duration_{0.0};
    std::chrono::steady_clock::time_point start_;
};

class PgDriver
{
public:
    static std::unique_ptr<PgDriver> getInstance()
    {
        if (!pInstance_)
        {
            pInstance_ = std::make_unique<PgDriver>();
        }
        return pInstance_;

    }

    void Init(const std::string &configFile = "/etc/ctmo/ctmo.cnf")
    {
        // Read settings
        const boost::property_tree::ptree pt;
        const boost::property_tree::ini_parser::read_ini(configFile, pt);

        // Read mysql information
        pgHost_ = pt.get<std::string>("POSTGRESQL.host");
        pgUser_ = pt.get<std::string>("POSTGRESQL.user");
        pgPassword_ = pt.get<std::string>("POSTGRESQL.password");
        pgDatabase_ = pt.get<std::string>("POSTGRESQL.database");

    }

    bool SaveBytea(const std::string &binaryData)
    {
        pqxx::work txn{connection_};

        pqxx::result r = txn.exec("SELECT name, salary FROM Employee");
        for (auto row: r)
            std::cout
                    // Address column by name.  Use c_str() to get C-style string.
                    << row["name"].c_str()
                    << " makes "
                    // Address column by zero-based index.  Use as<int>() to parse as int.
                    << row[1].as<int>()
                    << "."
                    << std::endl;

        // Not really needed, since we made no changes, but good habit to be
        // explicit about when the transaction is done.
        txn.commit();

        // Connection object goes out of scope here.  It closes automatically.
        return true;
    }

protected:
    PgDriver()
    {
        timer_.Start(POSTGRES_CONNECT_TIME_DELAY);
    };

    virtual ~PgDriver();

private:

    bool Connect()
    {
        // Check if we need to reconnect
        if (timer_.End())
        {
            if (connected_)
            {
                connection_.disconnect();
                connected_ = false;
            }
            timer_.Start(POSTGRES_CONNECT_TIME_DELAY);
        }


        try
        {
            connection_ = pqxx::connection(
                    "username=" + pgUser_ +
                    "host=" + pgHost_ +
                    "password=" + pgPassword_ +
                    "dbname=" + pgDatabase_
            );
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << std::endl;
        }

        if (connection_.is_open())
        {
            printf("Impossible to connect database");
            return false;
        }
        connected_ = true;
        return true;
    }


    static PgDriver pInstance_ = nullptr;
    std::string pgHost_{""};
    std::string pgUser_{""};
    std::string pgPassword_{""};
    std::string pgDatabase_{""};
    bool connected_{false};
    Timer timer_;
    pqxx::connection connection_;

};
}


#endif //CTMO_PGDRIVER_HPP
