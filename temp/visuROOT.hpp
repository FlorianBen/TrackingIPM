///////////////////////////////////////////////////////////////////////////////
//                                Declarations                               //
///////////////////////////////////////////////////////////////////////////////
#ifndef VISUROOT_HPP
#define VISUROOT_HPP

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TProfile.h>
#include <cmath>
#include <iostream>

class VisuROOT {
public:
  VisuROOT(int argc, char **argv);
  virtual ~VisuROOT();

  TH1F *getHistogramXini();
  TH1F *getHistogramYini();
  TH1F *getHistogramZini();
  TH1F *getHistogramXfin();
  TH1F *getHistogramYfin();
  TH1F *getHistogramZfin();

  void plot();

private:
  void formating();

  TApplication *app;
  TCanvas *canvas1;
  TH1F *hist_x_ini;
  TH1F *hist_y_ini;
  TH1F *hist_z_ini;
  TH1F *hist_x_fin;
  TH1F *hist_y_fin;
  TH1F *hist_z_fin;
};

class MCPROOT {
public:
  //! Default constructor
  MCPROOT(double MCPx, double MCPz, double span)
      : size_MCPx(MCPx), size_MCPz(MCPz), size_span(span) {
    canvas1 = new TCanvas("canv1", "Canvas 2", 1920, 1280);
    canvas2 = new TCanvas("canv2", "Canvas 3", 1920, 1280);
    canvas1->Divide(1, 2);
    canvas2->Divide(2, 1);
    make_bins();
    make_hist();
  }

  //! Destructor
  virtual ~MCPROOT() noexcept {
    delete hist_mcp_init_IPM1;
    delete hist_mcp_fin_IPM1;
    delete hist_mcp_init_IPM2;
    delete hist_mcp_fin_IPM2;
    delete hist_mcp_init_IPM3;
    delete hist_mcp_fin_IPM3;

    delete canvas1;
    delete canvas2;

    delete projX_ini_IPM1;
    delete projY_ini;

    delete projX_fin;
    delete projY_fin;

    delete[] binsx;
    delete[] binsz;
    delete[] binsz;
  }

  void FillInit_IPM1(double x, double z) { hist_mcp_init_IPM1->Fill(x, z, 1); }
  void FillFinal_IPM1(double x, double z) { hist_mcp_fin_IPM1->Fill(x, z, 1); }
  void FillInit_IPM2(double x, double z) { hist_mcp_init_IPM2->Fill(x, z, 1); }
  void FillFinal_IPM2(double x, double z) { hist_mcp_fin_IPM2->Fill(x, z, 1); }
  void FillInit_IPM3(double x, double z) { hist_mcp_init_IPM3->Fill(x, z, 1); }
  void FillFinal_IPM3(double x, double z) { hist_mcp_fin_IPM3->Fill(x, z, 1); }

  void Draw() {
    canvas1->cd(1);
    hist_mcp_init_IPM1->GetXaxis()->SetTitle("Transversal direction (m)");
    hist_mcp_init_IPM1->GetYaxis()->SetTitle("Longitudinal direction (m)");
    hist_mcp_init_IPM1->Draw("COLZ");
    // canvas1->cd(2);
    // hist_mcp_fin_IPM1->GetXaxis()->SetTitle("Transversal direction (m)");
    // hist_mcp_fin_IPM1->GetYaxis()->SetTitle("Longitudinal direction (m)");
    // hist_mcp_fin_IPM1->Draw("COLZ");

    // SetRealAspectRatio(canvas1);

    canvas1->Modified();
    canvas1->Update();

    canvas2->cd(1);
    projX_ini_IPM1 = hist_mcp_init_IPM1->ProjectionX();
    projX_ini_IPM1->SetLineColor(4);
    projX_ini_IPM1->SetTitle("Transversal profile");
    projX_ini_IPM1->GetXaxis()->SetTitle("Transversal direction (m)");
    projX_ini_IPM1->GetYaxis()->SetTitle("Counts");

    // projX_fin = hist_mcp_fin_IPM1->ProjectionX();
    // projX_fin->SetLineColor(2);
    // projX_fin->SetTitle("Transversal profile");
    // projX_fin->GetXaxis()->SetTitle("Transversal direction (m)");
    // projX_fin->GetYaxis()->SetTitle("Counts");

    // // projX_ini_IPM1->Draw("CP");
    projX_ini_IPM1->Fit("gaus");
    TF1 *f = (TF1 *)projX_ini_IPM1->GetListOfFunctions()->FindObject("gaus");
    f->SetLineColor(4);
    canvas2->Update();

    // TPaveStats *stats1 = (TPaveStats *)(projX_ini_IPM1->FindObject("stats"));
    // stats1->SetY1NDC(.85);
    // stats1->SetY2NDC(.95);
    // stats1->SetX1NDC(.8);
    // stats1->SetX2NDC(.95);
    // stats1->SetTextColor(4);

    // projX_fin->Fit("gaus", "", "SAMES");
    // canvas2->Update();
    // TPaveStats *stats2 = (TPaveStats *)(projX_fin->FindObject("stats"));
    // stats2->SetY1NDC(.75);
    // stats2->SetY2NDC(.85);
    // stats2->SetX1NDC(.80);
    // stats2->SetX2NDC(.95);
    // stats2->SetTextColor(2);

    canvas2->cd(2);
    projY_ini = hist_mcp_init_IPM1->ProjectionY();
    projY_ini->SetLineColor(4);
    projY_ini->SetTitle("Longitudinal profile");
    projY_ini->GetXaxis()->SetTitle("Longitudinal direction (m)");
    projY_ini->GetXaxis()->SetTitleOffset(1);
    projY_ini->GetYaxis()->SetTitle("Counts");

    // projY_fin = hist_mcp_fin_IPM1->ProjectionY();
    // projY_fin->SetLineColor(2);
    // projY_fin->SetTitle("Longitudinal profile");
    // projY_fin->GetXaxis()->SetTitle("Longitudinal direction (m)");
    // projY_fin->GetXaxis()->SetTitleOffset(1);
    // projY_fin->GetYaxis()->SetTitle("Counts");

    projY_ini->Draw("");
    canvas2->Update();
    TPaveStats *stats3 = (TPaveStats *)(projY_ini->FindObject("stats"));
    stats3->SetY1NDC(.85);
    stats3->SetY2NDC(.95);
    stats3->SetX1NDC(.8);
    stats3->SetX2NDC(.95);
    stats3->SetTextColor(4);

    // projY_fin->Draw("SAMES");
    // canvas2->Update();
    // TPaveStats *stats4 = (TPaveStats *)(projY_fin->FindObject("stats"));
    // stats4->SetY1NDC(.75);
    // stats4->SetY2NDC(.85);
    // stats4->SetX1NDC(.80);
    // stats4->SetX2NDC(.95);
    // stats4->SetTextColor(2);

    canvas2->Modified();
    canvas2->Update();
    // canvas2->SaveAs("profile.png");
  }

protected:
private:
  TH2F *hist_mcp_init_IPM1;
  TH2F *hist_mcp_fin_IPM1;
  TH2F *hist_mcp_init_IPM2;
  TH2F *hist_mcp_fin_IPM2;
  TH2F *hist_mcp_init_IPM3;
  TH2F *hist_mcp_fin_IPM3;

  TH1D *projX_ini_IPM1;
  TH1D *projY_ini;

  TH1D *projX_fin;
  TH1D *projY_fin;

  TCanvas *canvas1;
  TCanvas *canvas2;

  double size_MCPx;
  double size_MCPz;
  double size_span;

  int sizex;
  int sizez;

  double *binsx;
  double *binsz;

  void make_bins() {
    sizex = std::round(size_MCPx / size_span);
    sizez = std::round(size_MCPz / size_span);
    binsx = new double[sizex];
    binsz = new double[sizez];
    for (auto i = 0; i < sizex; i++) {
      binsx[i] = i * size_span - size_MCPx / 2;
    }

    for (auto i = 0; i < sizez; i++) {
      binsz[i] = i * size_span - size_MCPz / 2;
    }
  }

  void make_hist() {
    hist_mcp_init_IPM1 = new TH2F("hist_init_1", "Plane detection perfect IPM1",
                                  sizex - 1, binsx, sizez - 1, binsz);
    hist_mcp_fin_IPM1 = new TH2F("hist_final_1", "Plane detection real IPM1",
                                 sizex - 1, binsx, sizez - 1, binsz);
    hist_mcp_init_IPM2 = new TH2F("hist_init_2", "Plane detection perfect IPM2",
                                  sizex - 1, binsx, sizez - 1, binsz);
    hist_mcp_fin_IPM2 = new TH2F("hist_final_2", "Plane detection real IPM2",
                                 sizex - 1, binsx, sizez - 1, binsz);
    hist_mcp_init_IPM3 = new TH2F("hist_init_3", "Plane detection perfect IPM3",
                                  sizex - 1, binsx, sizez - 1, binsz);
    hist_mcp_fin_IPM3 = new TH2F("hist_final_3", "Plane detection real IPM3",
                                 sizex - 1, binsx, sizez - 1, binsz);
  }

  bool SetRealAspectRatio(TCanvas *const c, const Int_t axis = 1) {
    if (!c) {
      std::cerr << "Error in SetRealAspectRatio: canvas is NULL" << std::endl;
      return false;
    }
    // Get the current min-max values if SetUserRange has been called
    c->Update();
    const Double_t xmin = c->GetUxmin();
    const Double_t xmax = c->GetUxmax();
    const Double_t ymin = c->GetUymin();
    const Double_t ymax = c->GetUymax();
    // Get the length of zoomed x and y axes
    const Double_t xlength = xmax - xmin;
    const Double_t ylength = ymax - ymin;
    const Double_t ratio = xlength / ylength;

    // Get how many pixels are occupied by the canvas
    const Int_t npx = c->GetWw();
    const Int_t npy = c->GetWh();

    // Get x-y coordinates at the edges of the canvas (extrapolating outside
    // the axes, NOT at the edges of the histogram)
    const Double_t x1 = c->GetX1();
    const Double_t y1 = c->GetY1();
    const Double_t x2 = c->GetX2();
    const Double_t y2 = c->GetY2();

    // Get the length of extrapolated x and y axes
    const Double_t xlength2 = x2 - x1;
    const Double_t ylength2 = y2 - y1;
    const Double_t ratio2 = xlength2 / ylength2;

    // Get now number of pixels including window borders
    const Int_t bnpx = c->GetWindowWidth();
    const Int_t bnpy = c->GetWindowHeight();
    if (axis == 1) {
      c->SetCanvasSize(TMath::Nint(npy * ratio2), npy);
      c->SetWindowSize((bnpx - npx) + TMath::Nint(npy * ratio2), bnpy);
    } else if (axis == 2) {
      c->SetCanvasSize(npx, TMath::Nint(npx / ratio2));
      c->SetWindowSize(bnpx, (bnpy - npy) + TMath::Nint(npx / ratio2));
    } else {
      std::cerr << "Error in SetRealAspectRatio: axis value " << axis
                << " is neither 1 (resize along x-axis) nor 2 (resize along "
                   "y-axis)."
                << std::endl;
      return false;
    }
    return true;
  }
};

#endif /* VISUROOT_HPP */
