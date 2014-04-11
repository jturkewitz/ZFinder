// Standard Library
#include <algorithm>  //std::equal_range, std::sort, std::unique
#include <cmath>  // std::abs
#include <iostream>
#include <string>

// ROOT
#include <TCanvas.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>

// Plotting
#include "cross_check_plotter.h"

CrossCheckPlotter::CrossCheckPlotter(
        TFile* data_tfile,
        TFile* mc_tfile,
        std::string dir
        ) {
    // Use the same directory name for both MC and data
    setup(data_tfile, mc_tfile, dir, dir);
}

CrossCheckPlotter::CrossCheckPlotter(
        TFile* data_tfile,
        TFile* mc_tfile,
        std::string data_dir,
        std::string mc_dir
        ) {
    setup(data_tfile, mc_tfile, data_dir, mc_dir);
}

CrossCheckPlotter::~CrossCheckPlotter() {
    /* Clean up our heap variables */
    delete style_;
}

void CrossCheckPlotter::setup(
        TFile* data_tfile,
        TFile* mc_tfile,
        std::string data_dir,
        std::string mc_dir
        ) {
    // Save the pointers to the files
    data_tfile_ = data_tfile;
    mc_tfile_ = mc_tfile;

    // Save the directories
    data_dir_name_ = data_dir;
    mc_dir_name_ = mc_dir;

    // Set up the config map
    init_config_map();

    // Set up the style
    set_plot_style();
}

double CrossCheckPlotter::get_maximum(
        const TH1* const DATA_HISTO,
        const TH1* const MC_HISTO
        ) {
    /* Figure out the largest Y value */
    const double DATA_MAX = DATA_HISTO->GetMaximum();
    const double MC_MAX = MC_HISTO->GetMaximum();
    if (MC_MAX > DATA_MAX) {
        return MC_MAX;
    } else {
        return DATA_MAX;
    }
}

std::vector<double> CrossCheckPlotter::get_rebinning(
        std::vector<double> desired_bins,
        const TH1* const HISTO
        ) {
    /* Given a desired set of bin edges, and a histogram, finds the bin edges
     * in the histogram that most closely approximate the desired edges. This
     * is done by, for each desired edge, choosing the bin edge from the
     * histogram gram that is closest in terms of linear distance.
     */
    std::vector<double> out_vec;
    // Fill the old bins
    std::vector<double> old_bins;
    const int N_BINS = HISTO->GetXaxis()->GetNbins();
    for (int i_bin = 1; i_bin <= N_BINS; ++i_bin) {
        double bin_edge = HISTO->GetXaxis()->GetBinLowEdge(i_bin);
        old_bins.push_back(bin_edge);
    }
    // Add the high edge, which isn't included but is needed
    old_bins.push_back(HISTO->GetXaxis()->GetBinUpEdge(N_BINS));
    std::sort(old_bins.begin(), old_bins.end());

    // Loop through the desired bins.
    //
    // We use a binary search to find the actual
    // bins on either side of the desired bin, we then compare the distance and
    // pick the closest. At the end we remove duplicate entries. equal_range
    // will return a pointer to the first entry that is not less than our
    // desired bin, and one that is strictly greater than our test bin. This
    // means we might need to move one of the pointers back, and of course we
    // need to check that they don't run off the edge of the vector.
    std::sort(desired_bins.begin(), desired_bins.end());
    for (auto& i : desired_bins) {
        const double DESIRED_BIN = i;
        auto bounds = std::equal_range(old_bins.begin(), old_bins.end(), DESIRED_BIN);
        // Move the first pointer back if needed
        if (       bounds.first == bounds.second
                && bounds.first != old_bins.begin()
                && *bounds.first != DESIRED_BIN
           ) {
            --bounds.first;
        }
        // Check distance, pick the closest (or the first for a tie)
        const double FIRST_BIN = *bounds.first;
        const double SECOND_BIN = *bounds.second;
        const double FIRST_DIST = std::abs(FIRST_BIN - DESIRED_BIN);
        const double SECOND_DIST = std::abs(SECOND_BIN - DESIRED_BIN);
        if (SECOND_DIST >= FIRST_DIST) {
            out_vec.push_back(FIRST_BIN);
        } else {
            out_vec.push_back(SECOND_BIN);
        }
    }

    // Remove the non-unique entries. Unique moves the non-uniques to the end
    // and returns a pointer to the first non-unique item, erase then removes
    // them.
    std::sort(out_vec.begin(), out_vec.end());
    out_vec.erase(std::unique(out_vec.begin(), out_vec.end()), out_vec.end());

    // Return our vector
    return out_vec;
}

void CrossCheckPlotter::plot(
        const PlotType PLOT_TYPE,
        const std::string FILE_NAME
        ) {
    // Get our config_pair and some of the values from it
    config_map::iterator i_config_pair = conf_map_.find(PLOT_TYPE);
    if (i_config_pair == conf_map_.end()) {
        std::cout << "Missing PLOT_TYPE in conf_map_!" << std::endl;
        return;
    }
    PlotConfig plot_config = i_config_pair->second;
    const std::string HISTO_NAME = plot_config.histo_name;

    // Open the histograms
    const std::string DATA_HISTO_NAME = data_dir_name_ + "/" + HISTO_NAME;
    TH1D* data_histo;
    data_tfile_->GetObject(DATA_HISTO_NAME.c_str(), data_histo);
    if (!data_histo) {
        std::cout << "Can not open the Data Histogram!" << std::endl;
        return;
    }

    const std::string MC_HISTO_NAME = mc_dir_name_ + "/" + HISTO_NAME;
    TH1D* mc_histo;
    mc_tfile_->GetObject(MC_HISTO_NAME.c_str(), mc_histo);
    if (!mc_histo) {
        std::cout << "Can not open the MC Histogram!" << std::endl;
        return;
    }

    // Rebin if the binning is greater than 0 in size. If it is size one assume
    // we want a simple rebinning (where N bins are combined to 1), otherwise
    // the vector is the edges of the bins.
    if (plot_config.binning.size() == 1) {
        mc_histo->Rebin(static_cast<int>(plot_config.binning[0]));
        data_histo->Rebin(static_cast<int>(plot_config.binning[0]));
    } else if (plot_config.binning.size() > 1) {
        std::vector<double> new_binning = get_rebinning(
                plot_config.binning,
                data_histo
                );
        mc_histo = dynamic_cast<TH1D*>(
                mc_histo->Rebin(
                    new_binning.size() - 1,
                    "mc_rebinned_histo",
                    &new_binning[0]  // double*
                    )
                );
        data_histo = dynamic_cast<TH1D*>(
                data_histo->Rebin(
                    new_binning.size() - 1,
                    "data_rebinned_histo",
                    &new_binning[0]  // double*
                    )
                );
    }

    // Normalize areas
    const double MC_AREA = mc_histo->Integral();
    const double DATA_AREA = data_histo->Integral();
    mc_histo->Scale(DATA_AREA / MC_AREA);
    // Update uncertainties after rescaling
    data_histo->Sumw2();
    mc_histo->Sumw2();

    // Make a canvas to hold it
    TCanvas canvas("canvas", "canvas", X_VAL_, Y_VAL_);
    canvas.cd();

    // Set up the styles of the histograms
    style_->cd();
    // Title
    data_histo->SetTitle(0);  // Remove the title, we'll place it by hand
    mc_histo->SetTitle(0);
    // Axis labels
    data_histo->GetXaxis()->SetTitle(plot_config.x_label.c_str());
    mc_histo->GetXaxis()->SetTitle(plot_config.x_label.c_str());
    data_histo->GetYaxis()->SetTitle(plot_config.y_label.c_str());
    mc_histo->GetYaxis()->SetTitle(plot_config.y_label.c_str());
    // Position of axis labels
    mc_histo->GetYaxis()->SetTitleOffset(1.25);
    mc_histo->GetXaxis()->SetTitleOffset(1.1);
    // Marker, line, and fill style
    data_histo->SetMarkerStyle(kFullCircle);
    data_histo->SetMarkerColor(kBlack);
    data_histo->SetLineColor(kBlack);
    mc_histo->SetLineColor(kBlue);
    mc_histo->SetFillColor(kBlue);
    const int FORWARD_HATCH = 3004;
    //const int BACK_HATCH = 3005;
    mc_histo->SetFillStyle(FORWARD_HATCH);
    // Log
    canvas.SetLogy(plot_config.logy);


    // Set the plot range maximum based on the highest peak in either histo
    const double NEW_MAX = 1.05 * get_maximum(data_histo, mc_histo);
    data_histo->SetMaximum(NEW_MAX);

    // Set up the legend using the plot edges to set its location
    const double LEG_HEIGHT = 0.15;
    const double LEG_LENGTH = 0.15;
    TLegend legend(RIGHT_EDGE_ - LEG_LENGTH, TOP_EDGE_ - LEG_HEIGHT, RIGHT_EDGE_, TOP_EDGE_);
    legend.SetFillColor(kWhite);
    legend.AddEntry(data_histo, "Data", "p");
    legend.AddEntry(mc_histo, "MC", "f");
    legend.SetBorderSize(1);  // Remove drop shadow
    legend.SetFillStyle(0);  // Transparent

    // Add title
    TLatex *plot_title = NULL;
    if (plot_config.title != "") {
        const std::string TITLE = plot_config.title;
        plot_title = new TLatex(0.18, 0.93, TITLE.c_str());
        plot_title->SetNDC();
        plot_title->SetTextFont(42);
        plot_title->SetTextColor(1);
        plot_title->SetTextSize(0.06);
        plot_title->SetTextAlign(22);
        plot_title->SetTextAngle(0);
    }

    // Draw the histograms
    mc_histo->Draw("HIST");
    data_histo->Draw("E SAME");
    legend.Draw();
    if (plot_title != NULL) { plot_title->Draw(); }

    // Save the plot as a png
    canvas.Print(FILE_NAME.c_str(), "png");

    // Clean up
    delete plot_title;
}

void CrossCheckPlotter::set_plot_style() {
    style_ = new TStyle("style_","Style for P-TDR");

    // For the canvas:
    style_->SetCanvasBorderMode(0);
    style_->SetCanvasColor(kWhite);
    style_->SetCanvasDefX(0);  //Position on screen
    style_->SetCanvasDefY(0);

    // For the Pad:
    style_->SetPadBorderMode(0);
    style_->SetPadColor(kWhite);
    style_->SetPadGridX(false);
    style_->SetPadGridY(false);
    style_->SetGridColor(kBlack);
    style_->SetGridStyle(3);
    style_->SetGridWidth(1);

    // For the frame:
    style_->SetFrameBorderMode(0);
    style_->SetFrameBorderSize(1);
    style_->SetFrameFillColor(kWhite);
    style_->SetFrameFillStyle(0);
    style_->SetFrameLineColor(kBlack);
    style_->SetFrameLineStyle(1);
    style_->SetFrameLineWidth(1);

    // For the histo:
    // style_->SetHistFillColor(1);
    style_->SetHistFillStyle(0); //
    style_->SetHistLineColor(kBlack);
    style_->SetHistLineStyle(0);
    style_->SetHistLineWidth(1);

    style_->SetEndErrorSize(2);
    style_->SetErrorX(0.);

    style_->SetMarkerStyle(20);

    //For the fit/function:
    style_->SetOptFit(1);
    style_->SetFitFormat("5.4g");
    style_->SetFuncColor(kRed);
    style_->SetFuncStyle(1);
    style_->SetFuncWidth(1);

    //For the date:
    style_->SetOptDate(0);

    // For the statistics box:
    style_->SetOptFile(0);
    style_->SetOptStat(0);  // To display the mean and RMS: SetOptStat("mr");
    style_->SetStatColor(kWhite);
    style_->SetStatFont(42);
    style_->SetStatFontSize(0.025);
    style_->SetStatTextColor(kBlack);
    style_->SetStatFormat("6.4g");
    style_->SetStatBorderSize(1);
    style_->SetStatH(0.1);
    style_->SetStatW(0.15);

    // Margins:
    style_->SetPadTopMargin(1 - TOP_EDGE_);
    style_->SetPadBottomMargin(BOTTOM_EDGE_);
    style_->SetPadLeftMargin(LEFT_EDGE_);
    style_->SetPadRightMargin(1 - RIGHT_EDGE_);

    // For the Global title:
    style_->SetOptTitle(0);
    style_->SetTitleFont(42);
    style_->SetTitleColor(kBlack);
    style_->SetTitleTextColor(kBlack);
    style_->SetTitleFillColor(kWhite);  //10 is roughly kWhite, 10% grey?
    style_->SetTitleFontSize(0.05);

    // For the axis titles:
    style_->SetTitleColor(kBlack, "XYZ");
    style_->SetTitleFont(42, "XYZ");
    style_->SetTitleSize(0.06, "XYZ");
    style_->SetTitleXOffset(0.9);
    style_->SetTitleYOffset(1.25);

    // For the axis labels:
    style_->SetLabelColor(kBlack, "XYZ");
    style_->SetLabelFont(42, "XYZ");
    style_->SetLabelOffset(0.007, "XYZ");
    style_->SetLabelSize(0.05, "XYZ");

    // For the axis:
    style_->SetAxisColor(kBlack, "XYZ");
    style_->SetStripDecimals(true);
    style_->SetTickLength(0.03, "XYZ");
    //style_->SetNdivisions(510, "XYZ");
    style_->SetPadTickX(true);  // To get tick marks on the opposite side of the frame
    style_->SetPadTickY(true);

    // Change for log plots:
    style_->SetOptLogx(false);
    style_->SetOptLogy(false);
    style_->SetOptLogz(false);

    // Set the style
    style_->cd();
}

void CrossCheckPlotter::init_config_map() {
    /*
     * Here we fill the PlotConfigs. Unfortunately these must be set by hand,
     * and there should be one for every value in the PlotType enum.
     */
    // Z_MASS
    conf_map_.insert(
            config_pair(
                Z_MASS_ALL,
                PlotConfig(
                    "m_{ee} [GeV]",  // x_label
                    "Events",        // y_label
                    "",              // title
                    "Z0 Mass: All",  // histogram name (for reading in)
                    true,            // log Y axis
                    {}               // Desired new binning
                    )
                )
            );
    conf_map_.insert(
            config_pair(
                Z_MASS_COARSE,
                PlotConfig(
                    "m_{ee} [GeV]",
                    "Events",
                    "",
                    "Z0 Mass: Coarse",
                    true,
                    {}
                    )
                )
            );
    conf_map_.insert(
            config_pair(
                Z_MASS_FINE,
                PlotConfig(
                    "m_{ee} [GeV]",
                    "Events",
                    "",
                    "Z0 Mass: Fine",
                    true,
                    {}
                    )
                )
            );
    // Z Rapidity
    conf_map_.insert(
            config_pair(
                Z_RAPIDITY,
                PlotConfig(
                    "Y_{Z}",
                    "Events",
                    "",
                    "Z0 Rapidity",
                    true,
                    {}
                    )
                )
            );
    // PT
    conf_map_.insert(
            config_pair(
                Z_PT,
                PlotConfig(
                    "Z p_{T} [GeV]",
                    "Events",
                    "",
                    "Z0 p_{T}",
                    true,
                    {5}  // with one entry, just calls histo->Rebin(5);
                    )
                )
            );
    conf_map_.insert(
            config_pair(
                E0_PT,
                PlotConfig(
                    "e_{0} p_{T} [GeV]",
                    "Events",
                    "",
                    "p_{T,e_{0}}",
                    true,
                    {5}
                    )
                )
            );
    conf_map_.insert(
            config_pair(
                E1_PT,
                PlotConfig(
                    "e_{1} p_{T} [GeV]",
                    "Events",
                    "",
                    "p_{T,e_{1}}",
                    true,
                    {5}
                    )
                )
            );
    // Eta
    conf_map_.insert(
            config_pair(
                E0_ETA,
                PlotConfig(
                    "#eta_{e_{0}}",
                    "Events",
                    "",
                    "#eta_{e_{0}}",
                    true,
                    {}
                    )
                )
            );
    conf_map_.insert(
            config_pair(
                E1_ETA,
                PlotConfig(
                    "#eta_{e_{1}}",
                    "Events",
                    "",
                    "#eta_{e_{1}}",
                    true,
                    {}
                    )
                )
            );
    // Phi
    conf_map_.insert(
            config_pair(
                E0_PHI,
                PlotConfig(
                    "#phi_{e_{0}}",
                    "Events",
                    "",
                    "#phi_{e_{0}}",
                    false,
                    {}
                    )
                )
            );
    conf_map_.insert(
            config_pair(
                E1_PHI,
                PlotConfig(
                    "#phi_{e_{1}}",
                    "Events",
                    "",
                    "#phi_{e_{1}}",
                    false,
                    {}
                    )
                )
            );
    // Charge
    conf_map_.insert(
            config_pair(
                E0_CHARGE,
                PlotConfig(
                    "q_{e_{0}}",
                    "Events",
                    "",
                    "charge_{e_{0}}",
                    false,
                    {}
                    )
                )
            );
    conf_map_.insert(
            config_pair(
                E1_CHARGE,
                PlotConfig(
                    "q_{e_{1}}",
                    "Events",
                    "",
                    "charge_{e_{1}}",
                    false,
                    {}
                    )
                )
            );
    // Phi*
    conf_map_.insert(
            config_pair(
                PHISTAR,
                PlotConfig(
                    "#phi*",
                    "Events",
                    "",
                    "#phi*",
                    true,
                    // multiple entries means these are new bin edges
                    {0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0}
                    )
                )
            );
    // Vertexes
    conf_map_.insert(
            config_pair(
                N_VERTS,
                PlotConfig(
                    "Number of Vertexes",
                    "Events",
                    "",
                    "N_{Vertices}",
                    true,
                    {}
                    )
                )
            );
    // Electrons
    conf_map_.insert(
            config_pair(
                N_E,
                PlotConfig(
                    "Number of Electrons",
                    "Events",
                    "",
                    "N_{e}",
                    true,
                    {}
                    )
                )
            );
}