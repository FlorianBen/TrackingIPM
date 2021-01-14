#ifndef ALOGGER_HPP
#define ALOGGER_HPP
#include "spdlog/spdlog.h"

#include <memory>

namespace SpaceCharge {
class Logger {
public:
  static void Init();
  inline static std::shared_ptr<spdlog::logger> &GetLogger() { return logger; }

private:
  static std::shared_ptr<spdlog::logger> logger;
};
} // namespace SpaceCharge

// Client log macros
#define SC_TRACE(...)         ::SpaceCharge::Logger::GetLogger()->trace(__VA_ARGS__)
#define SC_INFO(...)          ::SpaceCharge::Logger::GetLogger()->info(__VA_ARGS__)
#define SC_WARN(...)          ::SpaceCharge::Logger::GetLogger()->warn(__VA_ARGS__)
#define SC_ERROR(...)         ::SpaceCharge::Logger::GetLogger()->error(__VA_ARGS__)
#define SC_CRITICAL(...)      ::SpaceCharge::Logger::GetLogger()->critical(__VA_ARGS__)

#endif