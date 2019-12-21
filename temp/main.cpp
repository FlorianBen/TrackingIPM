#include <iostream>

#include <tbb/parallel_for.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>

#include "SpaceCharge/bunch.hpp"
#include "SpaceCharge/fields.hpp"
#include "SpaceCharge/particle.hpp"

int main(int argc, char *argv[]) {
  /* code */

  auto app = new TApplication("TempApp", &argc, argv);

  SpaceCharge::Particle<double> part("proton", 1, SpaceCharge::cst::mproton,
                                     SpaceCharge::cst::lfactor::beta, 0.5);

  SpaceCharge::Particle<double> part2("electron", -1,
                                      SpaceCharge::cst::melectron,
                                      SpaceCharge::cst::lfactor::beta, 0.8);

  constexpr int nsize = 100;
  constexpr int tsize = 2000;

  auto step = 100.0e-3 / nsize;
  auto offset = -50e-3;

  auto t_step = 5e-12;

  float z[nsize], t[tsize];
  // float Ez[tsize][nsize];

  float *z_buffer = new float[tsize * nsize];

  tbb::parallel_for(0, tsize, [&](int k) {
    // for (auto k = 0; k < tsize; k++) {
    if (k % 100 == 0) {
      std::cout << k << "/" << tsize << std::endl;
    }
    auto tt = k * t_step;
    t[k] = tt;
    for (auto i = 0; i < nsize; i++) {
      SpaceCharge::FieldBunch<double> fields;

      SpaceCharge::GaussianBunch<double> bunch(
          part, 62.5 * SpaceCharge::uni::milli,
          1.0 / (352.0 * SpaceCharge::uni::mega), SpaceCharge::cst::dir::z);

      std::unique_ptr<SpaceCharge::Bunch<double>> bunch2(
          new SpaceCharge::GaussianBunch<double>(
              part, 62.5 * SpaceCharge::uni::milli,
              1.0 / (352.0 * SpaceCharge::uni::mega),
              SpaceCharge::cst::dir::z));

      std::unique_ptr<SpaceCharge::Bunch<double>> bunch3(
          new SpaceCharge::GaussianBunch<double>(
              part2, 30.5 * SpaceCharge::uni::milli,
              1.0 / (352.0 * SpaceCharge::uni::mega),
              SpaceCharge::cst::dir::z));

      // fields.usePeriodicity(false);

      fields.addBunch(std::move(bunch2));

      Eigen::Matrix<double, 4, 1> pos;
      pos(0) = SpaceCharge::cst::sol * tt;
      pos(3) = i * step + offset;
      z[i] = i * step + offset;

      // Ez[k][i] = (float)bunch.EfieldAt(pos)(3);
      z_buffer[k * nsize + i] = (float)fields.EfieldAt(pos)(3);
      ;
    }
  });

  auto c1 = new TCanvas("c1", "A Simple Graph Example", 200, 10, 1200, 800);

  auto gr = new TGraph(nsize); //, z, z_buffer[0]);
  gr->SetMaximum(90000);
  gr->SetMinimum(-90000);
  gr->SetTitle("Time t: 0.0");
  gr->Draw("AC");

  for (auto k = 0; k < tsize; k++) {
    for (auto i = 0; i < nsize; i++) {
      gr->SetPoint(i, z[i], z_buffer[k * nsize + i]);
    }
    std::stringstream stream;
    stream << std::scientific << t[k];
    std::string title = "Time t: ";
    title = title + stream.str();
    gr->SetTitle(title.c_str());
    gr->SetMaximum(90000);
    gr->SetMinimum(-90000);
    c1->Modified();
    c1->Update();
  }
  // app->Run();
  delete[] z_buffer;
  return 0;
}
