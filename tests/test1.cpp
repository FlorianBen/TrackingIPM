#include "gtest/gtest.h"

#include "SpaceCharge/particle.hpp"

class ParticleTest : public ::testing::Test {

protected:
  ParticleTest()
      : part("proton", 1, SpaceCharge::cst::mproton,
             SpaceCharge::cst::lfactor::beta, 0.5){};
  void SetUp() override {}

  // void TearDown() override {}

  SpaceCharge::Particle<float> part;
};

TEST_F(ParticleTest, ParticleIsProtonCharge) { EXPECT_EQ(part.getCharge(), 1); }

TEST_F(ParticleTest, ParticleIsProtonMass) {
  EXPECT_FLOAT_EQ(part.getMass(), SpaceCharge::cst::mproton);
}

TEST_F(ParticleTest, ParticleIsProtonSpeed) {
  EXPECT_FLOAT_EQ(part.getSpeed(), SpaceCharge::cst::sol/2.0);
  
}

int main(int argc, char **argv) {
  using namespace testing;
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}