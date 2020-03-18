#include <iostream>

#include <tbb/parallel_for.h>

#include "field_h5.hpp"
#include "vector.hpp"

bool save_data(Field &fieldv, std::vector<float> &time);

int main(int argc, char *argv[]) {
  /* code */
  SpaceCharge::Particle<double> part("proton", 1, SpaceCharge::cst::mproton,
                                     SpaceCharge::cst::lfactor::beta, 0.5);

  SpaceCharge::Particle<double> part2("electron", -1,
                                      SpaceCharge::cst::melectron,
                                      SpaceCharge::cst::lfactor::beta, 0.8);

  constexpr int nsize = 400;
  constexpr int tsize = 12000;

  auto step = 100.0e-3 / nsize;
  auto offset = -50e-3;

  auto t_step = 5e-12;

  float z[nsize], t[tsize];
  // float Ez[tsize][nsize];
  std::vector<float> timev(tsize);

  float *z_buffer = new float[tsize * nsize];

  Field fieldv(tsize*nsize);

  tbb::parallel_for(0, tsize, [&](int k) {
    if (k % 100 == 0) {
      std::cout << k << "/" << tsize << std::endl;
    }
    auto tt = k * t_step;
    t[k] = tt;
    timev.at(k) = tt;
    for (auto i = 0; i < nsize; i++) {
      SpaceCharge::FieldBunch<double> fields;
      std::unique_ptr<SpaceCharge::Bunch<double>> bunch2(
          new SpaceCharge::GaussianBunch<double>(
              part, 62.5 * SpaceCharge::uni::milli,
              1.0 / (352.0 * SpaceCharge::uni::mega),
              SpaceCharge::cst::dir::z));

      fields.usePeriodicity(true);

      fields.addBunch(std::move(bunch2));

      Eigen::Matrix<double, 4, 1> pos;
      pos(0) = SpaceCharge::cst::sol * tt;
      pos(3) = i * step + offset;
      z[i] = i * step + offset;

      Eigen::Matrix<double, 4, 1> temp = fields.EfieldAt(pos);

      // fieldv.push_back(
      FloatVector test;
      test[0] = (float)temp(1);
      test[1] = (float)temp(2);
      test[2] = (float)temp(3);
      fieldv[k * nsize + i] = test;
      z_buffer[k * nsize + i] = (float)temp(3);
      ;
    }
  });
  save_data(fieldv, timev);

  delete[] z_buffer;
  return 0;
}

bool save_data(Field &fieldv, std::vector<float> &time) {
  using namespace hdf5;
  file::File file = file::create("write_field.h5", file::AccessFlags::TRUNCATE);
  node::Group root_group = file.root();
  node::Group my_group = root_group.create_group("Scalar");
  using data_type = std::vector<float>;
  node::Dataset dataset = my_group.create_dataset(
      "time", datatype::create<data_type>(), dataspace::create(time));
  dataset.write(time);

  VectorAppender field_appender(create_vector_dataset("Field", root_group),
                                "field");

  std::for_each(fieldv.begin(), fieldv.end(), field_appender);

  return true;
}