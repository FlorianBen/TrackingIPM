#include "gtest/gtest.h"

#include "SpaceCharge/alogger.hpp"

class LoggerTest : public ::testing::Test {

protected:
  LoggerTest() {}
  void SetUp() override {}
};

TEST_F(LoggerTest, Logger) {
  EXPECT_FALSE(SpaceCharge::Logger::GetLogger());
  SpaceCharge::Logger::Init();
  EXPECT_TRUE(SpaceCharge::Logger::GetLogger());
  EXPECT_NO_THROW(SC_INFO("Print info"));
  EXPECT_NO_THROW(SC_TRACE("Print trance"));
  EXPECT_NO_THROW(SC_WARN("Print warn"));
  EXPECT_NO_THROW(SC_CRITICAL("Print critical"));
}

int main(int argc, char **argv) {
  using namespace testing;
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}