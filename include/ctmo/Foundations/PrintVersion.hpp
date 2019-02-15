#include "Version_ctmo.hpp"
#include <iostream>

namespace PrintVersion
{
void PrintVersion()
{

    const std::string gitBranch = GIT_BRANCH;
    const std::string gitHash = GIT_COMMIT_HASH;

    std::cout << "\n\n\n";
    std::cout << "\t\t\t\t ============= CTMO2.0 =============" << std::endl;
    std::cout << "\t\t\t\t -----------------------------------" << std::endl;
    std::cout << "\t\t\t\t -----------------------------------" << std::endl;
    std::cout << "\t\t\t\t\t\t .... \n \n ";
    std::cout << "\t\t\t\t\t Git Branch = " << gitBranch << std::endl;
    std::cout << "\t\t\t\t\t gitHash = " << gitHash << std::endl;
    std::cout << "\t\t\t\t\t\t .... \n \n";
    std::cout << "\t\t\t\t -----------------------------------" << std::endl;
    std::cout << "\t\t\t\t -----------------------------------" << std::endl;
    std::cout << "\t\t\t\t ===================================" << std::endl;
    std::cout << "\n\n\n";
}

std::string GetVersion()
{

    const std::string gitBranch = GIT_BRANCH;
    const std::string gitHash = GIT_COMMIT_HASH;
    const std::string version_message =
        std::string("\n\n\n") + std::string("\t ============= CTMO2.0 =============") +
        std::string("\n\t -----------------------------------") + std::string("\n\t -----------------------------------") +
        std::string("\n\t\t Git Branch = ") + gitBranch + std::string("\n\t\t gitHash = ") + gitHash +
        std::string("\n\t -----------------------------------") + std::string("\n\t -----------------------------------") +
        std::string("\n\t ===================================") + std::string("\n\n\n");

    return version_message;
}

} // namespace PrintVersion