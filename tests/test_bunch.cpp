#include "gtest/gtest.h"

#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/core/units.hpp"

class BunchTest : public ::testing::Test {

protected:
  BunchTest()
      : part("proton", 1, SpaceCharge::cst::mproton,
             SpaceCharge::cst::lfactor::ec, 90e6*SpaceCharge::uni::p::eV),
        bunch(part, 62.5 * SpaceCharge::uni::milli,
              1.0 / (352.0 * SpaceCharge::uni::mega),
              SpaceCharge::cst::dir::z){};
  void SetUp() override {

  }

  // void TearDown() override {}

  SpaceCharge::Particle<double> part;
  SpaceCharge::GaussianBunch<double> bunch;
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

TEST_F(BunchTest, BuncheSigma) {
  Eigen::Matrix<double, 3, 1> sigma;
  sigma << 1.25e-3,1.25e-3,2.8e-3;
  bunch.setSigma(sigma);
  ASSERT_NEAR(bunch.getSigmaB()(3),0.003068579377427003,1e-9);
}

int main(int argc, char **argv) {
  using namespace testing;
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}