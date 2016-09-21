#include <string.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>
#include <TKey.h>
#include <TEfficiency.h>
#include <Riostream.h>
#include "tdrStyle.C"
void calculate_jpsi_efficiencies (string file_name, string output_dir = "~/public_html/ZPhysics/tmp/Test90/", int polarization = 0  )
{
  TFile *theFile0 = new TFile( file_name.c_str());
  char *out_file_name;

  //TODO put this in rootrc
  gStyle->SetLineWidth(2.);
  gStyle->SetHistLineWidth(2.5);
  gROOT->ForceStyle();
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.13);
  //TH1D *h_n_truth_matched_muons_all = (TH1D*) theFile0->Get("ZFinder/All/N_truth_matched_jpsi_muons");
  //TH1D *h_n_truth_matched_muons_dimuon = (TH1D*) theFile0->Get("ZFinder/Dimuon_Jpsi/N_truth_matched_jpsi_muons");
  //TH1D *h_n_truth_matched_muons_dimuon_soft = (TH1D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Soft/N_truth_matched_jpsi_muons");
  //TH1D *h_n_truth_matched_muons_dimuon_vtx_comp = (TH1D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/N_truth_matched_jpsi_muons");
  //TH1D *h_n_truth_matched_muons_dimuon_primary_vert = (TH1D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Primary_Vertex/N_truth_matched_jpsi_muons");
  //TH1D *h_n_truth_matched_muons_jpsi = (TH1D*) theFile0->Get("ZFinder/Jpsi/N_truth_matched_jpsi_muons");

  //TH1D *h_jpsi_mass_fine_jpsi_mc = (TH1D*) theFile0->Get("ZFinder/MC_Jpsi/jpsi Mass: Fine");

  //double mc_events = h_jpsi_mass_fine_jpsi_mc->Integral();

  //double two_truth_matched_muons_dimuon = h_n_truth_matched_muons_dimuon->GetBinContent(3) ;
  //double two_truth_matched_muons_dimuon_soft = h_n_truth_matched_muons_dimuon_soft->GetBinContent(3);
  //double two_truth_matched_muons_dimuon_vtx_comp = h_n_truth_matched_muons_dimuon_vtx_comp->GetBinContent(3);
  //double two_truth_matched_muons_dimuon_primary_vert  = h_n_truth_matched_muons_dimuon_primary_vert->GetBinContent(3);
  //double two_truth_matched_muons_jpsi = h_n_truth_matched_muons_jpsi->GetBinContent(3);

  //std::cout << "mc_events " << mc_events << std::endl;
  //std::cout << "two_truth_matched_muons_dimuon " << two_truth_matched_muons_dimuon << std::endl;
  //std::cout << "two_truth_matched_muons_dimuon_soft " << two_truth_matched_muons_dimuon_soft  << std::endl;
  //std::cout << "two_truth_matched_muons_dimuon_vtx_comp " << two_truth_matched_muons_dimuon_vtx_comp  << std::endl;
  //std::cout << "two_truth_matched_muons_dimuon_primary_vert " << two_truth_matched_muons_dimuon_primary_vert  << std::endl;
  //std::cout << "two_truth_matched_muons_jpsi " << two_truth_matched_muons_jpsi  << std::endl;


  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Jpsi/jpsi_pt_vs_rap");
  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap_polarization_long");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Jpsi/jpsi_pt_vs_rap_polarization_long");
  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap_polarization_TPlusZero");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Jpsi/jpsi_pt_vs_rap_polarization_TPlusZero");

  char *hist_name_gen;
  char *hist_name_reco;
  if (polarization == 0) {
    hist_name_gen = "ZFinder/MC_All/jpsi_pt_vs_rap_finer";
    hist_name_reco = "ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap_finer";
  }
  else if (polarization == 1) {
    //hist_name_gen = "ZFinder/MC_All/jpsi_pt_vs_rap_finer_pos";
    //hist_name_reco = "ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap_finer_pos";
    hist_name_gen = "ZFinder/MC_All/jpsi_pt_vs_rap_finer_pos_0p1";
    hist_name_reco = "ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap_finer_pos_0p1";
  }
  else if (polarization == -1) {
    //hist_name_gen = "ZFinder/MC_All/jpsi_pt_vs_rap_finer_neg";
    //hist_name_reco = "ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap_finer_neg";
    hist_name_gen = "ZFinder/MC_All/jpsi_pt_vs_rap_finer_neg_0p1";
    hist_name_reco = "ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap_finer_neg_0p1";
  }
  else {
    std::cout << "Unknown polarization" << std::endl;
  }
  TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get(hist_name_gen);
  TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get(hist_name_reco);

  TH1D *jpsi_pt_reco = (TH1D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi p_{T}");
  TH1D *jpsi_pt_mc = (TH1D*) theFile0->Get("ZFinder/MC_All/jpsi p_{T}");

  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap_finer");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap_finer");

  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap");

  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap");


  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap_polarization_long");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap_polarization_long");
  //TH2D *jpsi_pt_vs_rap_mc = (TH2D*) theFile0->Get("ZFinder/MC_All/jpsi_pt_vs_rap_polarization_TPlusZero");
  //TH2D *jpsi_pt_vs_rap_jpsi = (TH2D*) theFile0->Get("ZFinder/Dimuon_Jpsi_Vertex_Compatible/jpsi_pt_vs_rap_polarization_TPlusZero");

  //TODO uncomment following two lines as we only use 1 bin of pT
  //jpsi_pt_vs_rap_mc->Rebin2D(21,1);
  //jpsi_pt_vs_rap_jpsi->Rebin2D(21,1);

  bool do_not_draw = true;
  //bool do_not_draw = false;

  if(do_not_draw) {
    jpsi_pt_vs_rap_mc->Rebin2D(3,1);
    jpsi_pt_vs_rap_jpsi->Rebin2D(3,1);
    //jpsi_pt_vs_rap_mc->Rebin2D(21,1);
    //jpsi_pt_vs_rap_jpsi->Rebin2D(21,1);
  }
  else {
    jpsi_pt_vs_rap_mc->Rebin2D(3,2);
    jpsi_pt_vs_rap_jpsi->Rebin2D(3,2);
  }

  TH1D *jpsi_pt_clone = jpsi_pt_mc->Clone();
  jpsi_pt_clone->Divide(jpsi_pt_reco, jpsi_pt_mc, 1.0, 1.0, "B");

  TCanvas *c0 = new TCanvas("c0", "c0",800,600);
  c0->cd();
  //c1->SetLogy();
  jpsi_pt_clone->Draw();
  jpsi_pt_clone->GetXaxis()->SetRangeUser(8.5,100);
  jpsi_pt_clone->GetYaxis()->SetTitle("Efficiency / GeV");
  jpsi_pt_clone->GetXaxis()->SetTitle("p_{T,J/#psi} (GeV)");
  std::string image_name0 = "";
  image_name0.append(output_dir);
  image_name0.append("jpsi_efficiency_pt");
  image_name0.append(".png");
  c0->Print(image_name0.c_str() , "png");
    
  std::cout << jpsi_pt_vs_rap_mc->Integral() << std::endl;
  std::cout << jpsi_pt_vs_rap_jpsi->Integral() << std::endl;
  std::cout << "eff overall: " << jpsi_pt_vs_rap_jpsi->Integral() / jpsi_pt_vs_rap_mc->Integral() << std::endl;
  jpsi_pt_vs_rap_mc->Sumw2();
  jpsi_pt_vs_rap_jpsi->Sumw2();
  TH2D *acc_eff_map = jpsi_pt_vs_rap_mc->Clone();
  acc_eff_map->Divide(jpsi_pt_vs_rap_jpsi, jpsi_pt_vs_rap_mc, 1.0, 1.0, "B");

  acc_eff_map->GetXaxis()->SetTitle("J/#psi Rapidity");
  acc_eff_map->GetYaxis()->SetTitle("J/#psi p_{T} (GeV)");
  acc_eff_map->GetYaxis()->SetTitleOffset(1.1);

  if(!do_not_draw) {
    acc_eff_map->GetYaxis()->SetRangeUser(8.5,18.5);
    acc_eff_map->GetYaxis()->SetMoreLogLabels(1);
  }
  acc_eff_map->GetZaxis()->SetRangeUser(0,1);



  //acc_eff_map->Divide(jpsi_pt_vs_rap_jpsi, jpsi_pt_vs_rap_mc, 1.0, 1.0);
  //acc_eff_map->Divide(jpsi_pt_vs_rap_jpsi, jpsi_pt_vs_rap_mc, "cl=0.683 b(1,1) mode");



  TCanvas *c1 = new TCanvas("c1", "c1",800,600);
  c1->cd();
  c1->SetLogy();
  if(do_not_draw) {
    acc_eff_map->Draw("colz");
  }
  else {
    acc_eff_map->Draw("textcolz");
  }
  

  //TEfficiency * pEff = 0;
  //if (TEfficiency::CheckConsistency(*jpsi_pt_vs_rap_jpsi, *jpsi_pt_vs_rap_mc, "w")) { 
  //  pEff = new TEfficiency(*jpsi_pt_vs_rap_jpsi, *jpsi_pt_vs_rap_mc);
  //}
  //pEff->Draw("colz");

  int nxbins = acc_eff_map->GetNbinsX();
  int nybins = acc_eff_map->GetNbinsY();
  for (int x=1;x<=nxbins;++x)
  {
    for (int y=1;y<=nybins;++y)
    {
      double acc_eff = acc_eff_map->GetBinContent(x,y);
      double acc_error_low = acc_eff_map->GetBinError(x,y);
      double acc_error_high = acc_eff_map->GetBinError(x,y);
      //double acc_error_low = acc_eff_map->GetBinErrorLow(x,y);
      //double acc_error_low = acc_eff_map->GetBinErrorLow(x,y);
      //double acc_error_high = acc_eff_map->GetBinErrorHigh(x,y);
      TAxis *xaxis = acc_eff_map->GetXaxis();  
      double xlow = xaxis->GetBinLowEdge(x);
      double xwidth = xaxis->GetBinWidth(x);
      double xhigh = xlow + xwidth;
      TAxis *yaxis = acc_eff_map->GetYaxis();  
      double ylow = yaxis->GetBinLowEdge(y);
      double ywidth = yaxis->GetBinWidth(y);
      double yhigh = ylow + ywidth;
      //std::cout << "xlow " << xlow  << " xhigh " << xhigh << " ylow " << ylow << " yhigh " << yhigh << " acc_eff " << acc_eff <<
      //  " acc_error_low " << acc_error_low << " acc_error_high " << acc_error_high << std::endl;
      std::cout << "{" << xlow << ",  " << xhigh << ",  " << ylow << ",  " << yhigh << ",  " <<
        acc_eff << ",  " << acc_error_low << ",  " << acc_error_high << "   }," << std::endl;
    }
  }

  //TODO verify just dividing seems to give similar/same results instead of using TEfficiency
  //pEff->SetStatisticOption(1);
  //int nxbins = jpsi_pt_vs_rap_mc->GetNbinsX();
  //int nybins = jpsi_pt_vs_rap_mc->GetNbinsY();
  //for (int x=1;x<=nxbins;++x)
  //{
  //  for (int y=1;y<=nybins;++y)
  //  {
  //    int bin = pEff->GetGlobalBin(x,y);
  //    double acc_eff = pEff->GetEfficiency(bin);
  //    double acc_error_low = pEff->GetEfficiencyErrorLow(bin);
  //    double acc_error_high = pEff->GetEfficiencyErrorUp(bin);
  //    TAxis *xaxis = jpsi_pt_vs_rap_mc->GetXaxis();  
  //    double xlow = xaxis->GetBinLowEdge(x);
  //    double xwidth = xaxis->GetBinWidth(x);
  //    double xhigh = xlow + xwidth;
  //    TAxis *yaxis = jpsi_pt_vs_rap_mc->GetYaxis();  
  //    double ylow = yaxis->GetBinLowEdge(y);
  //    double ywidth = yaxis->GetBinWidth(y);
  //    double yhigh = ylow + ywidth;
  //    //std::cout << "xlow " << xlow  << " xhigh " << xhigh << " ylow " << ylow << " yhigh " << yhigh << " acc_eff " << acc_eff <<
  //    //  " acc_error_low " << acc_error_low << " acc_error_high " << acc_error_high << std::endl;
  //    //eta low eta high pt_low pt_high acc_eff acc_eff_err_low acc_eff_err_high
  //    std::cout << "{" << xlow << ",  " << xhigh << ",  " << ylow << ",  " << yhigh << ",  " <<
  //      acc_eff << ",  " << acc_error_low << ",  " << acc_error_high << "   }," << std::endl;
  //  }
  //}

  std::string image_name1 = "";
  std::string path = output_dir;
  image_name1.append(output_dir);
  image_name1.append("acceptance_map");
  if (polarization == 0) {
    out_file_name = "acc_eff_map_pol_0.root";
    image_name1.append("_pol0");
  }
  else if (polarization == 1) {
    out_file_name = "acc_eff_map_pol_pos.root";
    image_name1.append("_pol1");
  }
  else if (polarization == -1) {
    out_file_name = "acc_eff_map_pol_neg.root";
    image_name1.append("_polneg1");
  }
  else {
    std::cout << "Unknown polarization" << std::endl;
  }
  image_name1.append(".png");
  c1->Print(image_name1.c_str() , "png");

  TFile output(out_file_name,"new");
  acc_eff_map->Write();
}
