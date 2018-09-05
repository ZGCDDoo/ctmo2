

#include <gtest/gtest.h>
#include "../src/Includes/Models/UTensorSimple.hpp"

const double DELTA = 1e-11;
const std::string fname = "../test/data/cdmft_triangle/testtriangle.json";

using namespace Models;

TEST(UTensorSimple, Init)
{
    std::ifstream fin(fname);
    Json jj;
    fin >> jj;
    fin.close();
    UTensor ut(jj);

    ASSERT_NEAR(ut.auxMu(), 0.39418507926716284, DELTA);
    ASSERT_DOUBLE_EQ(ut.U(), 3.0);
    ASSERT_DOUBLE_EQ(ut.JH(), 0.0);
    ASSERT_DOUBLE_EQ(ut.UPrime(), 3.0);
}
