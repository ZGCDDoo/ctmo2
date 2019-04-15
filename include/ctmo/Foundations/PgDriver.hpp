//
// Created by charles-david on 14/04/19.
//
#pragma once
#ifndef CTMO_PGDRIVER_HPP
#define CTMO_PGDRIVER_HPP

namespace DB
{
class PgDriver
{
   explicit PgDriver(const std::string& configFile="/etc/ctmo/ctmo.cnf");

};
}


#endif //CTMO_PGDRIVER_HPP
