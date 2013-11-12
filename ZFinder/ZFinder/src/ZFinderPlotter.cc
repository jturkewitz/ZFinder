#include "ZFinder/ZFinder/interface/ZFinderPlotter.h"

// Standard Library
#include <string>  // string

// Root
#include <TCanvas.h>  // TCanvas

// ZFinder Code
#include "ZFinderElectron.h"  // ZFinderElectron

// Constructor
ZFinderPlotter::ZFinderPlotter(TDirectory* tdir) {
    /* 
     * Initialize a set of histograms and associate them with a given TDirectory.
     */ 

    // Set tdir to the internal directory and then set it as the current
    // working directory
    tdir_ = tdir;
    tdir_->cd();

    // Set up histograms
    // z0_mass_coarse_
    const std::string z0_mass_coarse_name = "Z0 Mass: Coarse";
    z0_mass_coarse_ = new TH1I(z0_mass_coarse_name.c_str(), z0_mass_coarse_name.c_str(), 100, 50., 150.);
    z0_mass_coarse_->SetDirectory(tdir_);
    z0_mass_coarse_->GetXaxis()->SetTitle("m_{ee} [GeV]");
    z0_mass_coarse_->GetYaxis()->SetTitle("Counts / GeV");

    // z0_mass_fine_
    const std::string z0_mass_fine_name = "Z0 Mass: Fine";
    z0_mass_fine_ = new TH1I(z0_mass_fine_name.c_str(), z0_mass_fine_name.c_str(), 80, 80., 100.);
    z0_mass_fine_->SetDirectory(tdir_);
    z0_mass_fine_->GetXaxis()->SetTitle("m_{ee} [GeV]");
    z0_mass_fine_->GetYaxis()->SetTitle("Counts / 0.25 GeV");

    // z0_rapidity_
    const std::string z0_rapidity_name = "Z0 Rapidity";
    z0_rapidity_ = new TH1I(z0_rapidity_name.c_str(), z0_rapidity_name.c_str(), 100, -5., 5.);
    z0_rapidity_->SetDirectory(tdir_);
    z0_rapidity_->GetXaxis()->SetTitle("Y_{ee}");
    z0_rapidity_->GetYaxis()->SetTitle("Counts");

    // z0_pt
    const std::string z0_pt_name = "Z0 p_{T}";
    z0_pt = new TH1I(z0_pt_name.c_str(), z0_pt_name.c_str(), 200, 0., 200.);
    z0_pt->SetDirectory(tdir_);
    z0_pt->GetXaxis()->SetTitle("p_{T,Z}");
    z0_pt->GetYaxis()->SetTitle("Counts / GeV");

    // e0_pt
    const std::string e0_pt_name = "p_{T,e_{0}}";
    e0_pt = new TH1I(e0_pt_name.c_str(), e0_pt_name.c_str(), 200, 0., 200.);
    e0_pt->SetDirectory(tdir_);
    e0_pt->GetXaxis()->SetTitle("p_{T,e_{0}}");
    e0_pt->GetYaxis()->SetTitle("Counts / GeV");

    // e1_pt
    const std::string e1_pt_name = "p_{T,e_{1}}";
    e1_pt = new TH1I(e1_pt_name.c_str(), e1_pt_name.c_str(), 200, 0., 200.);
    e1_pt->SetDirectory(tdir_);
    e1_pt->GetXaxis()->SetTitle("p_{T,e_{0}}");
    e1_pt->GetYaxis()->SetTitle("Counts / GeV");

    // e0_eta_
    const std::string e0_eta_name = "#eta_{e_{0}}";
    e0_eta_ = new TH1I(e0_eta_name.c_str(), e0_eta_name.c_str(), 50, -5., 5.);
    e0_eta_->SetDirectory(tdir_);
    e0_eta_->GetXaxis()->SetTitle("#eta_{e_{0}}");
    e0_eta_->GetYaxis()->SetTitle("Counts");

    // e1_eta_
    const std::string e1_eta_name = "#eta_{e_{1}}";
    e1_eta_ = new TH1I(e1_eta_name.c_str(), e1_eta_name.c_str(), 50, -5., 5.);
    e1_eta_->SetDirectory(tdir_);
    e1_eta_->GetXaxis()->SetTitle("#eta_{e_{1}}");
    e1_eta_->GetYaxis()->SetTitle("Counts");

    // e0_phi_
    const std::string e0_phi_name = "#phi_{e_{0}}";
    e0_phi_ = new TH1I(e0_phi_name.c_str(), e0_phi_name.c_str(), 60, -3.15, 3.15);
    e0_phi_->SetDirectory(tdir_);
    e0_phi_->GetXaxis()->SetTitle("#phi_{e_{0}}");
    e0_phi_->GetYaxis()->SetTitle("Counts");

    // e1_phi_
    const std::string e1_phi_name = "#phi_{e_{1}}";
    e1_phi_ = new TH1I(e1_phi_name.c_str(), e1_phi_name.c_str(), 50, -3.15, 3.15);
    e1_phi_->SetDirectory(tdir_);
    e1_phi_->GetXaxis()->SetTitle("#phi_{e_{1}}");
    e1_phi_->GetYaxis()->SetTitle("counts");

    // phistar
    const std::string phistar_name = "#phi^{*}";
    phistar_ = new TH1I(phistar_name.c_str(), phistar_name.c_str(), 100, 0., 1.);
    phistar_->SetDirectory(tdir_);
    phistar_->GetXaxis()->SetTitle("#phi^{*}");
    phistar_->GetYaxis()->SetTitle("Counts");

    // pileup
    const std::string pileup_name = "Pileup";
    pileup_ = new TH1I(pileup_name.c_str(), pileup_name.c_str(), 100, 0., 100.);
    pileup_->SetDirectory(tdir_);
    pileup_->GetXaxis()->SetTitle("Pileup");
    pileup_->GetYaxis()->SetTitle("Counts");

}

void ZFinderPlotter::Fill(ZFinderEvent const * const z_event, const short electron_0 = 0, const short electron_1 = 1) {
    /* 
     * Given a ZEvent, fills all the histograms. 
     *
     * electron_0 and electron_1 can be used to assign ZEvent.eN to the given
     * number in the histogram. For example, assigning electron_0 = 1 will fill
     * the e0 histograms with data from ZEvent->e1.
     */
    // Z Info
    z0_mass_coarse_->Fill(ZEvent->z.m);
    z0_mass_fine_->Fill(ZEvent->z.m);
    z0_rapidity_->Fill(ZEvent->z.y);
    z0_pt->Fill(ZEvent->z.pt);
    phistar->Fill(ZEvent->z.phistar);

    // Fill the histograms with the information from the approriate electron
    if ( electron_0 == 0 && electron_1 == 1 ) {
        e0_pt->Fill(ZEvent->e0->pt);
        e0_eta_->Fill(ZEvent->e0->eta);
        e0_phi_->Fill(ZEvent->e0->phi);
        e1_pt->Fill(ZEvent->e1->pt);
        e1_eta_->Fill(ZEvent->e1->eta);
        e1_phi_->Fill(ZEvent->e1->phi);
    } else if ( electron_0 == 1 && electron_1 == 0 ) {
        e0_pt->Fill(ZEvent->e1->pt);
        e0_eta_->Fill(ZEvent->e1->eta);
        e0_phi_->Fill(ZEvent->e1->phi);
        e1_pt->Fill(ZEvent->e0->pt);
        e1_eta_->Fill(ZEvent->e0->eta);
        e1_phi_->Fill(ZEvent->e0->phi);
    }

    // Event Info
    pileup->Fill(ZEvent->vert.num);
}

void ZFinderPlotter::Print() {
    // Get the name of the TDir
    std::string basename;
    basename.assign(tdir->GetName());

    // Set Image Sizes
    const int X_SIZE = 1280;
    const int Y_SIZE = 640;

    // Write all PNGs
    std::string z0_mass_coarse_Str = basename + "_";
    TCanvas* z0_mass_coarse_C = new TCanvas(z0_mass_coarse_.c_str(), z0_mass_coarse_.c_str(), X_SIZE, Y_SIZE);
    z0_mass_coarse_->Draw();
    z0_mass_coarse_C->Print(z0_mass_coarse_Str.c_str());

    std::string z0_mass_fine_Str = basename + "_";
    TCanvas* z0_mass_fine_C = new TCanvas(z0_mass_fine_.c_str(), z0_mass_fine_.c_str(), X_SIZE, Y_SIZE);
    z0_mass_fine_->Draw();
    z0_mass_fine_C->Print(z0_mass_fine_Str.c_str());

    std::string z0_rapidity_Str = basename + "_";
    TCanvas* z0_rapidity_C = new TCanvas(z0_rapidity_.c_str(), z0_rapidity_.c_str(), X_SIZE, Y_SIZE);
    z0_rapidity_->Draw();
    z0_rapidity_C->Print(z0_rapidity_Str.c_str());

    std::string z0_ptStr = basename + "_";
    TCanvas* z0_ptC = new TCanvas(z0_pt.c_str(), z0_pt.c_str(), X_SIZE, Y_SIZE);
    z0_pt->Draw();
    z0_ptC->Print(z0_ptStr.c_str());

    std::string e0_ptStr = basename + "_";
    TCanvas* e0_ptC = new TCanvas(e0_pt.c_str(), e0_pt.c_str(), X_SIZE, Y_SIZE);
    e0_pt->Draw();
    e0_ptC->Print(e0_ptStr.c_str());

    std::string e1_ptStr = basename + "_";
    TCanvas* e1_ptC = new TCanvas(e1_pt.c_str(), e1_pt.c_str(), X_SIZE, Y_SIZE);
    e1_pt->Draw();
    e1_ptC->Print(e1_ptStr.c_str());

    std::string e0_eta_Str = basename + "_";
    TCanvas* e0_eta_C = new TCanvas(e0_eta_.c_str(), e0_eta_.c_str(), X_SIZE, Y_SIZE);
    e0_eta_->Draw();
    e0_eta_C->Print(e0_eta_Str.c_str());

    std::string e1_eta_Str = basename + "_";
    TCanvas* e1_eta_C = new TCanvas(e1_eta_.c_str(), e1_eta_.c_str(), X_SIZE, Y_SIZE);
    e1_eta_->Draw();
    e1_eta_C->Print(e1_eta_Str.c_str());

    std::string e0_phi_Str = basename + "_";
    TCanvas* e0_phi_C = new TCanvas(e0_phi_.c_str(), e0_phi_.c_str(), X_SIZE, Y_SIZE);
    e0_phi_->Draw();
    e0_phi_C->Print(e0_phi_Str.c_str());

    std::string e1_phi_Str = basename + "_";
    TCanvas* e1_phi_C = new TCanvas(e1_phi_.c_str(), e1_phi_.c_str(), X_SIZE, Y_SIZE);
    e1_phi_->Draw();
    e1_phi_C->Print(e1_phi_Str.c_str());

    std::string phistarStr = basename + "_";
    TCanvas* phistarC = new TCanvas(phistar_.c_str(), phistar_.c_str(), X_SIZE, Y_SIZE);
    phistar_->Draw();
    phistarC->Print(phistarStr.c_str());

    std::string pileupStr = basename + "_";
    TCanvas* pileupC = new TCanvas(pileup.c_str(), pileup.c_str(), X_SIZE, Y_SIZE);
    pileup->Draw();
    pileupC->Print(pileupStr.c_str());
}