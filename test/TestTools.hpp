#include <boost/filesystem.hpp>
#include <iostream>
#include "../src/Includes/Utilities/Utilities.hpp"

namespace TestTools
{

using namespace boost::filesystem;

void RemoveFilesForTests()
{
    std::vector<std::string> files = {"tloc.arma", "tktilde.arma", "hybFM.arma", "config.dat"};

    for (const std::string &ss : files)
    {
        path filepath(ss);
        if (exists(filepath) && is_regular_file(filepath))
        {
            remove(filepath);
        }
    }
}

Json BuildJson()
{
    Json tJson = R"(
    {   "NOrb": 2,
        "NKPTS": 100,
        "tParameters": 
            {"00": 
                {"tIntra": 0.0, "tx": -1.0, "ty": -1.0, "tz": 0.0, "tx=y": -0.40, "tx=-y": 0.0, "tx=z": -0.0, "tx=-z": 0.0, "ty=z": 0.0, "ty=-z": 0.0, "t2x" : 0.0, "t2y": 0.0, "t2z": 0.0},
            "01":
                {"tIntra": -0.01, "tx": -1.02, "ty": -1.02, "tz": 0.0, "tx=y": 0.230, "tx=-y": 0.230, "tx=z": -0.0, "tx=-z": 0.0, "ty=z": 0.0, "ty=-z": 0.0, "t2x" : -0.90, "t2y": -0.90, "t2z": 0.0},
            "11":
                {"tIntra": 0.0, "tx": -1.0, "ty": -1.0, "tz": 0.0, "tx=y": -0.40, "tx=-y": 0.0, "tx=z": -0.0, "tx=-z": 0.0, "ty=z": 0.0, "ty=-z": 0.0, "t2x" : 0.0, "t2y": 0.0, "t2z": 0.0}

            }
    }
    )"_json;

    return tJson;
}
} // namespace TestTools