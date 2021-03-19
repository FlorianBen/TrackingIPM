#include "gtest/gtest.h"

#include "SpaceCharge/core/particle.hpp"
#include "SpaceCharge/core/units.hpp"

class ParticleTest : public ::testing::Test {

protected:
  ParticleTest()
      : part1("proton", 1, SpaceCharge::cst::mproton,
              SpaceCharge::cst::lfactor::beta, 0.5),
        part2("electron", -1, SpaceCharge::cst::melectron,
              SpaceCharge::cst::lfactor::gamma, 2.294),
        part3("neutron", 0, SpaceCharge::cst::mneutron,
              SpaceCharge::cst::lfactor::ec, 10 * SpaceCharge::uni::p::eV),
        part4("neutron", 0, SpaceCharge::cst::mneutron,
              SpaceCharge::cst::lfactor::speed, 10.0){

        };
  void SetUp() override {}

  void TearDown() override {}

protected:
  SpaceCharge::Particle<float> part1;
  SpaceCharge::Particle<float> part2;
  SpaceCharge::Particle<float> part3;
  SpaceCharge::Particle<float> part4;
};

TEST_F(ParticleTest, Particle1ProtonCharge) { EXPECT_EQ(part1.getCharge(), 1); }

TEST_F(ParticleTest, Particle1ProtonMass) {
  EXPECT_FLOAT_EQ(part1.getMass(), SpaceCharge::cst::mproton);
}

TEST_F(ParticleTest, Particle1ProtonSpeed) {
  EXPECT_FLOAT_EQ(part1.getSpeed(), SpaceCharge::cst::sol / 2.0);
  part1.setSpeed(2000);
  EXPECT_FLOAT_EQ(part1.getSpeed(), 2000);  
}

TEST_F(ParticleTest, Particle2ElectronCharge) {
  EXPECT_EQ(part2.getCharge(), -1);
}

TEST_F(ParticleTest, Particle2ElectronMass) {
  EXPECT_FLOAT_EQ(part2.getMass(), SpaceCharge::cst::melectron);
}

TEST_F(ParticleTest, Particle2ElectronSpeed) {
  EXPECT_FLOAT_EQ(part2.getSpeed(),
                  part2.getSpeed()); // SpaceCharge::cst::sol * 0.9);
  part2.setBeta(0.9);
  EXPECT_FLOAT_EQ(part2.getBeta(), 0.9);
}

TEST_F(ParticleTest, Particle3NeutronCharge) {
  EXPECT_EQ(part3.getCharge(), 0);
}

TEST_F(ParticleTest, Particle3NeutronMass) {
  EXPECT_FLOAT_EQ(part3.getMass(), SpaceCharge::cst::mneutron);
}

TEST_F(ParticleTest, Particle3NeutronSpeed) {
  EXPECT_FLOAT_EQ(part3.getSpeed(),
                  part3.getSpeed()); // SpaceCharge::cst::sol * 0.9);
  part3.setGamma(1.05);
  EXPECT_FLOAT_EQ(part3.getGamma(), 1.05);
}

TEST_F(ParticleTest, Particle4NeutronCharge) {
  EXPECT_EQ(part4.getCharge(), 0);
}

TEST_F(ParticleTest, Particle4NeutronMass) {
  EXPECT_FLOAT_EQ(part4.getMass(), SpaceCharge::cst::mneutron);
}

TEST_F(ParticleTest, Particle4NeutronSpeed) {
  EXPECT_FLOAT_EQ(part4.getSpeed(),
                  part4.getSpeed()); // SpaceCharge::cst::sol * 0.9);
  part4.setEc(10 * SpaceCharge::uni::p::keV);
  EXPECT_FLOAT_EQ(part4.getEc(), 10 * SpaceCharge::uni::p::keV);
}

int main(int argc, char **argv) {
  using namespace testing;
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}