#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "SpaceCharge/alogger.hpp"

namespace SpaceCharge {
std::shared_ptr<spdlog::logger> Logger::logger;

void Logger::Init() {
  std::vector<spdlog::sink_ptr> logSinks;

  logSinks.emplace_back(
      std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
  logSinks[0]->set_pattern("%+ [from thread %t]");

  logSinks.emplace_back(
      std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs.txt", true));
  logSinks[1]->set_pattern("%+ [from thread %t]");

  logger = std::make_shared<spdlog::logger>("Logger", begin(logSinks),
                                            end(logSinks));
  spdlog::register_logger(logger);
}
} // namespace SpaceCharge