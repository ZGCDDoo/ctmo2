

#include <gtest/gtest.h>
#include "../src/Includes/Models/UTensorSimple.hpp"

const double DELTA = 1e-11;
const std::string fname = "../test/data/cdmft_triangle/testtriangle.json";

using namespace Models;

TEST(UTensorSimple, Init)
{
    //====================  Test for Diagonal Orbitals, with UPrime=JH=0  =================
    std::ifstream fin(fname);
    Json jj;
    fin >> jj;
    fin.close();
    UTensor ut(jj);

    ASSERT_NEAR(ut.auxMu(), 1.8941850792671628 - 1.5, DELTA);
    ASSERT_DOUBLE_EQ(ut.U(), 3.0);
    ASSERT_DOUBLE_EQ(ut.JH(), 0.0);
    ASSERT_DOUBLE_EQ(ut.UPrime(), 0.0);

    //======================  End Test for Diagonal Orbitals, with UPrime=JH=0  =============================

    // Test with JH=0, but UPrime=3.21

    jj["solver"]["isOrbitalDiagonal"] = false;
    jj["model"]["UPrime"] = 3.21;
    jj["model"]["nOrb"] = 2;

    UTensor ut2(jj);

    ASSERT_NEAR(ut2.auxMu(), 1.8941850792671628 - 1.5 - (2.0 - 1.0) * 3.21 / 2.0, DELTA);
    ASSERT_DOUBLE_EQ(ut2.U(), 3.0);
    ASSERT_DOUBLE_EQ(ut2.JH(), 0.0);
    ASSERT_DOUBLE_EQ(ut2.UPrime(), 3.21);

    //============ End Test with JH=0, but UPrime=3.21 ==============

    //============ Test with JH=1.33, but UPrime=3.21 ========================

    jj["model"]["J_H"] = 1.33;

    UTensor ut3(jj);

    ASSERT_NEAR(ut3.auxMu(), -2.150814920732837, DELTA);
    ASSERT_DOUBLE_EQ(ut3.U(), 3.0);
    ASSERT_DOUBLE_EQ(ut3.JH(), 1.33);
    ASSERT_DOUBLE_EQ(ut3.UPrime(), 3.21);
    //======== End Test with JH=1.33, but UPrime=3.21 ===============
}
