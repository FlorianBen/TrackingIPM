#include "gtest/gtest.h"

#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/field/fields.hpp"
#include "SpaceCharge/field/bunch.hpp"
#include "SpaceCharge/core/units.hpp"

class FieldBunchTest : public ::testing::Test {

protected:
  FieldBunchTest()
      : part("proton", 1, SpaceCharge::cst::mproton,
             SpaceCharge::cst::lfactor::beta, 0.5),
        bunch(part, 62.5 * SpaceCharge::uni::milli,
              1.0 / (352.0 * SpaceCharge::uni::mega),
              SpaceCharge::cst::dir::z){};
  void SetUp() override {
    // field.addBunch(bunch);
    field.usePeriodicity(true);
    Eigen::Matrix<double, 4, 1> pos;
    pos << 0.0, 0.0, 0.0, -1.0;
    std::cout << bunch.EfieldAt(pos) << std::endl;
  }

  // void TearDown() override {}

  SpaceCharge::Particle<double> part;
  SpaceCharge::GaussianBunch<double> bunch;
  SpaceCharge::BunchEMField<double> field;
};

TEST_F(FieldBunchTest, Potential) {
  ASSERT_NEAR(1.0, 1.0, 1e-9);
}

int main(int argc, char **argv) {
  using namespace testing;
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}