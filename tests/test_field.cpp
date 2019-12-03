#include "gtest/gtest.h"

#include "SpaceCharge/fields.hpp"
#include "SpaceCharge/units.hpp"

class BunchTest : public ::testing::Test {

protected:
  BunchTest()
      : part("proton", 1, SpaceCharge::cst::mproton,
             SpaceCharge::cst::lfactor::beta, 0.5),
        bunch(part, 62.5 * SpaceCharge::uni::milli,
              1.0 / (352.0 * SpaceCharge::uni::mega),
              SpaceCharge::cst::dir::z) {};
  void SetUp() override {}

  // void TearDown() override {}

  SpaceCharge::Particle<float> part;
  SpaceCharge::GaussianBunch<float> bunch;
  //SpaceCharge::Field<float> field;
};

TEST_F(BunchTest, BunchTiming) {
  EXPECT_FLOAT_EQ(bunch.getFreqRF(), 352.0 * SpaceCharge::uni::mega);
  EXPECT_FLOAT_EQ(bunch.getBunchPeriod(), 2.840909091 * SpaceCharge::uni::nano);
}

TEST_F(BunchTest, BuncheCharge) {
  EXPECT_FLOAT_EQ(bunch.getCharges(),
                  0.0625 * 1.0 / (352.0 * SpaceCharge::uni::mega));
  auto ncharges = (0.0625 * 0.00286 / SpaceCharge::cst::e) /
                  (0.00286 * 352.0 * SpaceCharge::uni::mega);
  EXPECT_FLOAT_EQ(bunch.getNbParticles(), ncharges);
}

int main(int argc, char **argv) {
  using namespace testing;
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}