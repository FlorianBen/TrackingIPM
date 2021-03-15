#ifndef ALOGGER_HPP
#define ALOGGER_HPP
#include "spdlog/spdlog.h"

#include <memory>

namespace SpaceCharge {
/**
 * \class Logger alogger.hpp
 * \brief Internal logger.
 **/
class Logger {
public:
  /**
   * \brief Initialize the logger.
   * Call this function before using the logger.
   **/
  static void Init();
  /**
   * \brief Get the logger.
   **/
  inline static std::shared_ptr<spdlog::logger> &GetLogger() { return logger; }

private:
  static std::shared_ptr<spdlog::logger> logger;
};
} // namespace SpaceCharge

// Log macros
#define SC_TRACE(...) ::SpaceCharge::Logger::GetLogger()->trace(__VA_ARGS__)
#define SC_INFO(...) ::SpaceCharge::Logger::GetLogger()->info(__VA_ARGS__)
#define SC_WARN(...) ::SpaceCharge::Logger::GetLogger()->warn(__VA_ARGS__)
#define SC_ERROR(...) ::SpaceCharge::Logger::GetLogger()->error(__VA_ARGS__)
#define SC_CRITICAL(...)                                                       \
  ::SpaceCharge::Logger::GetLogger()->critical(__VA_ARGS__)

#endif