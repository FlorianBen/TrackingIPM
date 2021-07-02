#include "SpaceCharge/field/field_fdm.hpp"

#include "SpaceCharge/core/alogger.hpp"

namespace SpaceCharge {

FieldFDM::FieldFDM(int nx, int ny, double dx, double dy)
    : nx(nx), ny(ny), dx(dx), dy(dy), mat(nx, ny) {
  Logger::GetLogger()->info("FieldFDM: Matrix ({},{}) created", nx, ny);
}

FieldFDM::~FieldFDM() {}

} // namespace SpaceCharge
