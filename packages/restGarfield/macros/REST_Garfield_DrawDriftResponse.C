#ifndef REST_GARFIELD_DRAWDRIFTRESPONSE_C
#define REST_GARFIELD_DRAWDRIFTRESPONSE_C

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "TArrow.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TRestDetectorDriftResponse.h"

/// \brief Draws a drift response computed by restGarfield (TRestDetectorDriftResponse, the `drift` section):
/// what the electric field does to the electrons on their way to the readout, from the summary of each start
/// point of the grid (TRestDetectorDriftResponse::GetNodeSummary). Maps on the readout plane of the layer of
/// start points at the distance z (by default the farthest layer, where the field acts longest):
/// - survival: fraction of the electrons that reach the readout (lost at the walls, the field shapers, out of
///   the gas or attached);
/// - displacement: mean shift of the arrival from the start point projected on the readout (u, v of the
///   readout plane; x, y for a normal along z), as a map of its size with arrows (scaled, see the title); 0
///   in a uniform field, the distortion of the field otherwise;
/// - radial displacement versus the radius of the start point: positive, the electrons are pushed outwards;
/// - mean drift time;
/// - drift time versus z at the centre of the grid and at its edges along x (the drift velocity is the
///   inverse of the slope, constant in a uniform field);
/// - transverse spread (the standard deviation of the arrival, the diffusion), sqrt(sigmaU^2 + sigmaV^2).
///
/// The tree DriftResponse of the file has the same summary per start point to draw other things, e.g.
/// DriftResponse->Draw("meanU:x", "z == 90 && abs(y) < 5").
///
///    restRoot
///    REST_Garfield_DrawDriftResponse("AlphaCAMM_drift.root")
///    REST_Garfield_DrawDriftResponse("AlphaCAMM_drift.root", 45)   // the layer closest to z = 45 mm
///
/// \param fileName Response file (`~` is expanded to HOME).
/// \param z Coordinate (mm) of the layer of start points along the third axis of the grid (z); negative: the
/// last layer.
/// \return The canvas (3x2 pads), or nullptr if the file cannot be read.
TCanvas* REST_Garfield_DrawDriftResponse(const std::string& fileName, double z = -1) {
    std::string path = fileName;
    const char* home = std::getenv("HOME");
    if (home && !path.empty() && path[0] == '~') path = home + path.substr(1);
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "REST_Garfield_DrawDriftResponse: cannot open " << path << std::endl;
        return nullptr;
    }
    std::unique_ptr<TRestDetectorDriftResponse> read;
    try {
        read = std::make_unique<TRestDetectorDriftResponse>(TRestDetectorDriftResponse::Read(file.get()));
    } catch (const std::exception& error) {
        std::cerr << "REST_Garfield_DrawDriftResponse: " << error.what() << std::endl;
        return nullptr;
    }
    const TRestDetectorDriftResponse& response = *read;

    const int nx = response.fGridNodes[0], ny = response.fGridNodes[1], nz = response.fGridNodes[2];
    const double xMin = response.fXRange.first, xMax = response.fXRange.second;
    const double yMin = response.fYRange.first, yMax = response.fYRange.second;
    const double zMin = response.fZRange.first, zMax = response.fZRange.second;
    const double dx = (xMax - xMin) / (nx - 1), dy = (yMax - yMin) / (ny - 1), dz = (zMax - zMin) / (nz - 1);
    const int k = z < 0 ? nz - 1 : std::clamp(static_cast<int>(std::lround((z - zMin) / dz)), 0, nz - 1);
    const double layer = zMin + k * dz;

    // Maps of the layer: a bin per start point
    auto map = [&](const std::string& name, const std::string& title) {
        auto* histogram = new TH2D((name + "_" + std::to_string(k)).c_str(), title.c_str(), nx, xMin - dx / 2,
                                   xMax + dx / 2, ny, yMin - dy / 2, yMax + dy / 2);
        histogram->SetDirectory(nullptr);
        histogram->SetStats(false);
        return histogram;
    };
    const std::string at = " (z = " + std::to_string(std::lround(layer)) + " mm)";
    auto* survival = map("survival", "Survival" + at + ";x (mm);y (mm)");
    auto* displacement = map("displacement", "");
    auto* time = map("time", "Mean drift time (#mus)" + at + ";x (mm);y (mm)");
    auto* spread = map("spread", "Transverse spread (mm)" + at + ";x (mm);y (mm)");
    auto* radial = new TGraph();
    radial->SetTitle(("Radial displacement" + at + ";r (mm);(u, v) #upoint #hat{r} (mm)").c_str());
    double largest = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            const size_t index = response.NodeIndex(i, j, k);
            const auto summary = response.GetNodeSummary(index);
            const auto position = response.GetNodePosition(index);
            survival->SetBinContent(i + 1, j + 1, summary.fSurvival);
            if (summary.fSurvival <= 0) continue;
            const double shift = std::hypot(summary.fMeanU, summary.fMeanV);
            largest = std::max(largest, shift);
            displacement->SetBinContent(i + 1, j + 1, shift);
            time->SetBinContent(i + 1, j + 1, summary.fMeanT);
            spread->SetBinContent(i + 1, j + 1, std::hypot(summary.fSigmaU, summary.fSigmaV));
            const double r = std::hypot(position.X(), position.Y());
            if (r > 0) {
                radial->AddPoint(r, (summary.fMeanU * position.X() + summary.fMeanV * position.Y()) / r);
            }
        }
    }
    // Arrows of the mean displacement, the largest one as long as 0.9 of the spacing of the grid
    const double scale = largest > 0 ? 0.9 * std::min(dx, dy) / largest : 1;
    char title[200];
    snprintf(title, sizeof(title), "Mean displacement (mm)%s, arrows x %.3g (largest %.3g mm);x (mm);y (mm)",
             at.c_str(), scale, largest);
    displacement->SetTitle(title);

    // Drift time versus z at the centre of the grid and at its edges along x
    auto* versusZ = new TMultiGraph("", "Mean drift time;z (mm);t (#mus)");
    auto* legend = new TLegend(0.12, 0.7, 0.5, 0.88);
    const int columns[3] = {nx / 2, 0, nx - 1};
    const int colors[3] = {1, 2, 4};
    for (int c = 0; c < 3; c++) {
        auto* graph = new TGraph();
        for (int kk = 0; kk < nz; kk++) {
            const size_t index = response.NodeIndex(columns[c], ny / 2, kk);
            const auto summary = response.GetNodeSummary(index);
            if (summary.fSurvival > 0) graph->AddPoint(response.GetNodePosition(index).Z(), summary.fMeanT);
        }
        graph->SetLineColor(colors[c]);
        graph->SetMarkerColor(colors[c]);
        graph->SetMarkerStyle(20);
        versusZ->Add(graph, "LP");
        const auto position = response.GetNodePosition(response.NodeIndex(columns[c], ny / 2, 0));
        const std::string label = "x = " + std::to_string(std::lround(position.X())) +
                                  ", y = " + std::to_string(std::lround(position.Y())) + " mm";
        legend->AddEntry(graph, label.c_str(), "lp");
    }

    auto* canvas = new TCanvas(("DriftResponse_" + std::to_string(k)).c_str(), path.c_str(), 1500, 900);
    canvas->Divide(3, 2);
    canvas->cd(1)->SetRightMargin(0.15);
    survival->SetMinimum(0);
    survival->SetMaximum(1);
    survival->Draw("colz");
    canvas->cd(2)->SetRightMargin(0.15);
    displacement->Draw("colz");
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            const size_t index = response.NodeIndex(i, j, k);
            const auto summary = response.GetNodeSummary(index);
            if (summary.fSurvival <= 0 || std::hypot(summary.fMeanU, summary.fMeanV) <= 0) continue;
            const auto position = response.GetNodePosition(index);
            auto* arrow = new TArrow(position.X(), position.Y(), position.X() + scale * summary.fMeanU,
                                     position.Y() + scale * summary.fMeanV, 0.005, "|>");
            arrow->Draw();
        }
    }
    canvas->cd(3);
    radial->SetMarkerStyle(20);
    radial->SetMarkerSize(0.6);
    radial->Draw("AP");
    canvas->cd(4)->SetRightMargin(0.15);
    time->Draw("colz");
    canvas->cd(5);
    versusZ->Draw("A");
    legend->Draw();
    canvas->cd(6)->SetRightMargin(0.15);
    spread->Draw("colz");
    canvas->Update();
    return canvas;
}

#endif
