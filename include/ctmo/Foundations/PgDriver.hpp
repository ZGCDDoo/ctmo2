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
    static PgDriver *getInstance()
    {
        if (!pInstance_)
        {
            pInstance_ = new PgDriver();
        }
        return pInstance_;

    }



    void SaveBytea(const std::string &binaryData, const int &simulation_id)
    {
        pqxx::connection conn(
                "user=" + pgUser_ + " " +
                "host=" + pgHost_ + " " +
                "password=" + pgPassword_ + " "
                "dbname=" + pgDatabase_
        );
        pqxx::work wTxn{conn};
//        conn.prepare("insertBatch", "INSERT INTO" + pgTable_ + "(simulation_id , batch) VALUES ($1, $2)");
//        pqxx::result r = wTxn.prepared("insertBatch")(simulation_id)(pqxx::binarystring(binaryData)).exec();


        wTxn.exec_params(
                "INSERT INTO " + pgTable_ + "(simulation_id , batch) VALUES ($1, $2)",
                simulation_id,
                pqxx::binarystring(binaryData)
        );
        wTxn.commit();


    }

protected:


private:
    void Init(const std::string &configFile = "/etc/ctmo/ctmo.cnf")
    {
        // Read settings
        boost::property_tree::ptree pt;
        boost::property_tree::ini_parser::read_ini(configFile, pt);

        // Read mysql information
        pgHost_ = pt.get<std::string>("POSTGRESQL.host");
        pgUser_ = pt.get<std::string>("POSTGRESQL.user");
        pgPassword_ = pt.get<std::string>("POSTGRESQL.password");
        pgDatabase_ = pt.get<std::string>("POSTGRESQL.database");
        pgTable_ = pt.get<std::string>("POSTGRESQL.table");

    }

    PgDriver()
    {
        Init();
    };

    virtual ~PgDriver();

//    bool Connect()
//    {
//        // Check if we need to reconnect
//        if (timer_.End())
//        {
//            if (connected_)
//            {
//                connection_.disconnect();
//                connected_ = false;
//            }
//            timer_.Start(POSTGRES_CONNECT_TIME_DELAY);
//        }
//
//
//        try
//        {
//            connection_ = pqxx::connection(
//                    "username=" + pgUser_ +
//                    "host=" + pgHost_ +
//                    "password=" + pgPassword_ +
//                    "dbname=" + pgDatabase_
//            );
//        }
//        catch (const std::exception &e)
//        {
//            std::cerr << e.what() << std::endl;
//        }
//
//        if (!connection_.is_open())
//        {
//            printf("Impossible to connect database");
//            return false;
//        }
//        connected_ = true;
//        return true;
//    }


    static PgDriver *pInstance_;
    std::string pgHost_{""};
    std::string pgUser_{""};
    std::string pgPassword_{""};
    std::string pgDatabase_{""};
    std::string pgTable_{""};


};

PgDriver *PgDriver::pInstance_ = nullptr;

PgDriver::~PgDriver()
{};
}


#endif //CTMO_PGDRIVER_HPP
