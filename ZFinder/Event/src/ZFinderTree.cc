#include "ZFinder/Event/interface/ZFinderTree.h"

// Standard Library
#include <algorithm>
#include <vector>  // std::min


#include <typeinfo>


// ZFinder Code


namespace zf {
  // Constructor
  //ZFinderTree::ZFinderTree(const ZDefinition& zdef, TFileDirectory& tdir, const bool IS_MC) : IS_MC_(IS_MC) {
  ZFinderTree::ZFinderTree(TFileDirectory& tdir, const bool IS_MC) : IS_MC_(IS_MC) {
    // Make the directory to save files to
    //tdir.cd();

    // Make the Tree to write to
    //tree_ = new TTree(zdef.NAME.c_str(), zdef.NAME.c_str());
    tree_ = new TTree("zfinder_tree", "zfinder_tree");
    //const std::string Z_CODE = "z_m/D:z_pt:z_y:z_phi:z_phistar:z_eta:z_vtx_prob:z_vtx_x:z_vtx_y:z_vtx_z:daughter0_pt:daughter1_pt:daughter0_charge/I:daughter1_charge";
    //const std::string Z_CODE = "z_m/D:z_pt:z_y:z_phi:z_phistar:z_eta:z_vtx_prob:z_vtx_x:z_vtx_y:z_vtx_z:daughter0_pt:daughter1_pt:daughter0_charge/I:daughter1_charge";
    const std::string Z_CODE = "z_m/D:z_pt:z_y:z_phi:z_phistar:z_eta:z_vtx_prob:z_vtx_x:z_vtx_y:z_vtx_z:daughter0_pt:daughter0_eta:daughter0_phi:daughter1_pt:daughter1_eta:daughter1_phi:daughter0_charge/I:daughter1_charge";
    //const std::string Z_CODE_INT = "daughter0_charge/I:daughter1_charge";
    const std::string JPSI_CODE = "jpsi_m/D:jpsi_pt:jpsi_y:jpsi_phi:jpsi_eta:jpsi_vtx_prob:jpsi_vtx_x:jpsi_vtx_y:jpsi_vtx_z:jpsi_tau_xy:jpsi_tau_z:jpsi_distance_xy:jpsi_distance_z:jpsi_eff:jpsi_acc_eff:jpsi_scale_factor:muon0_pt:muon0_eta:muon0_phi:muon1_pt:muon1_eta:muon1_phi:muon0_charge/I:muon1_charge:has_muons_in_eta_window:has_high_pt_muons";
    //const std::string CODE = "z_m/D:z_y:z_phistar_born:z_phistar_dressed:z_phistar_naked:z_phistar_sc:z_pt:z_eta:e_pt0:e_pt1:e_eta0:e_eta1:e_phi0:e_phi1:e_rnine0:e_rnine1:n_true_pileup:e_charge0/I:e_charge1:n_verts:t0tight/O:t1tight";
    tree_->Branch("reco_z", &reco_z_, Z_CODE.c_str());
    //tree_->Branch("reco_z", &reco_z_, Z_CODE_INT.c_str());
    tree_->Branch("reco_z_from_muons", &reco_z_from_muons_, Z_CODE.c_str());
    //tree_->Branch("reco_z_from_muons", &reco_z_from_muons_, Z_CODE_INT.c_str());
    tree_->Branch("reco_jpsi", &reco_jpsi_, JPSI_CODE.c_str());

    //TODO jpsi->ee
    //----------------
    tree_->Branch("reco_jpsi_from_electrons", &reco_jpsi_from_electrons_, JPSI_CODE.c_str());
    //------------------

    //tree_->Branch("reco", &reco_, CODE.c_str());

    //if (IS_MC_) {
    //  tree_->Branch("truth_z", &truth_z_, Z_CODE.c_str());
    //  tree_->Branch("truth_jpsi", &truth_jpsi_, JPSI_CODE.c_str());
    //}
    
    
    tree_->Branch("truth_z_muons", &truth_z_muons_, Z_CODE.c_str());
    tree_->Branch("truth_z_electrons", &truth_z_electrons_, Z_CODE.c_str());

    tree_->Branch("truth_jpsi", &truth_jpsi_, JPSI_CODE.c_str());
    //const std::string EVENT_CODE = "event_number/i:run_number:is_mc/O";

    //const std::string EVENT_CODE = "event_weight/D:event_number/I:run_number:n_verts:truth_n_verts:is_mc/O:found_high_pt_muons_from_z:found_good_muons_from_z:found_dimuon_z_compatible_vertex:found_z_to_muons_mass:found_high_pt_electrons_from_z:found_good_electrons_from_z:found_dielectron_z_compatible_vertex:found_z_to_electrons_mass:found_dimuon_jpsi_with_muons_in_eta_window:found_dimuon_jpsi_with_high_pt_muons:found_dimuon_jpsi_with_soft_id_and_high_pt_muons:found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex:found_good_dimuon_jpsi_compatible_with_primary_vertex:found_jpsi";

    //TODO jpsi->ee
    //-----------------
    const std::string EVENT_CODE = "event_weight/D:event_number/I:run_number:n_verts:truth_n_verts:is_mc/O:found_high_pt_muons_from_z:found_good_muons_from_z:found_dimuon_z_compatible_vertex:found_z_to_muons_mass:found_high_pt_electrons_from_z:found_good_electrons_from_z:found_dielectron_z_compatible_vertex:found_z_to_electrons_mass:found_dimuon_jpsi_with_muons_in_eta_window:found_dimuon_jpsi_with_high_pt_muons:found_dimuon_jpsi_with_soft_id_and_high_pt_muons:found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex:found_good_dimuon_jpsi_compatible_with_primary_vertex:found_jpsi:found_dimuon_jpsi_from_electrons_with_muons_in_eta_window:found_dimuon_jpsi_from_electrons_with_high_pt_muons:found_dimuon_jpsi_from_electrons_with_soft_id_and_high_pt_muons:found_dimuon_jpsi_from_electrons_with_good_muons_and_compatible_muon_vertex:found_good_dimuon_jpsi_from_electrons_compatible_with_primary_vertex:found_jpsi_from_electrons";
    tree_->Branch("event_info", &event_, EVENT_CODE.c_str());
    //-------------------

    //if (IS_MC_) {
    //  tree_->Branch("weight_size", &weight_size_, "weight_size/I");
    //  tree_->Branch("weights", weights_, "weights[weight_size]/D");
    //  tree_->Branch("weight_ids", weight_ids_, "weight_ids[weight_size]/I");
    //  tree_->Branch("weight_cteq_size", &weight_cteq_size_, "weight_cteq_size/I");
    //  tree_->Branch("weights_cteq", weights_cteq_, "weights_cteq[weight_cteq_size]/D");
    //  tree_->Branch("weight_mstw_size", &weight_mstw_size_, "weight_mstw_size/I");
    //  tree_->Branch("weights_mstw", weights_mstw_, "weights_mstw[weight_mstw_size]/D");
    //  tree_->Branch("weight_nnpdf_size", &weight_nnpdf_size_, "weight_nnpdf_size/I");
    //  tree_->Branch("weights_nnpdf", weights_nnpdf_, "weights_nnpdf[weight_nnpdf_size]/D");
    //  tree_->Branch("weight_fsr", &weight_fsr_, "weight_fsr/D");
    //}
  }

  ZFinderTree::~ZFinderTree() {
    tree_->Write();
    // Clean up our pointer
    delete tree_;
  }

  void ZFinderTree::Fill(const ZFinderEvent& zfe) {
    // If the ZDefinition doesn't pass, escape
    //if (!zf_event.ZDefPassed(zdef_name_)) {
    //  return;
    //}
    
    //TODO maybe want to add some similar cut here to winnow down dataset

    // Clear our branches
    reco_z_.clear_values();
    reco_z_from_muons_.clear_values();
    truth_z_electrons_.clear_values();
    truth_z_muons_.clear_values();
    
    reco_jpsi_.clear_values();


    //TODO jpsi->ee testing
    reco_jpsi_from_electrons_.clear_values();


    truth_jpsi_.clear_values();
   
    event_.clear_values();

    //weight_id_vector_.clear();

    //// Set the weights
    //if (IS_MC_) {
    //  // The weight from the generator
    //  const double GEN_WEIGHT = zf_event.weight_natural_mc;
    //  weight_id_vector_.push_back(std::make_pair(WeightID::GEN_MC, GEN_WEIGHT));

    //  // The pileup weight
    //  const double VERT_WEIGHT = zf_event.weight_vertex;
    //  weight_id_vector_.push_back(std::make_pair(WeightID::PILEUP, VERT_WEIGHT));
    //  const double VERT_WEIGHT_PLUS = zf_event.weight_vertex_plus;
    //  weight_id_vector_.push_back(std::make_pair(WeightID::PILEUP_PLUS, VERT_WEIGHT_PLUS));
    //  const double VERT_WEIGHT_MINUS = zf_event.weight_vertex_minus;
    //  weight_id_vector_.push_back(std::make_pair(WeightID::PILEUP_MINUS, VERT_WEIGHT_MINUS));

    //  // Get the scale factors weight and fill up weight_id_vector_ with them
    //  const cutlevel_vector* clv = zf_event.GetZDef(zdef_name_);
    //  FillCutWeights(clv);

    //  // Loop over out vector of weights and save them to the arrays so that
    //  // the tree can grab the values. If the vector is longer than
    //  // MAX_SIZE_, we stop there so as not to overflow our array!
    //  weight_size_ = weight_id_vector_.size();
    //  for (int i = 0; i < std::min(weight_size_, MAX_SIZE_); ++i) {
    //    const int WEIGHT_ID = weight_id_vector_.at(i).first;
    //    const double WEIGHT = weight_id_vector_.at(i).second;
    //    weights_[i] = WEIGHT;
    //    weight_ids_[i] = WEIGHT_ID;
    //  }
    //  weight_cteq_size_ = zf_event.weights_cteq.size();
    //  for (int i = 0; i < std::min(weight_cteq_size_, MAX_SIZE_PDF_); ++i) {
    //    weights_cteq_[i] = zf_event.weights_cteq[i];
    //  }
    //  weight_mstw_size_ = zf_event.weights_mstw.size();
    //  for (int i = 0; i < std::min(weight_mstw_size_, MAX_SIZE_PDF_); ++i) {
    //    weights_mstw_[i] = zf_event.weights_mstw[i];
    //  }
    //  weight_nnpdf_size_ = zf_event.weights_nnpdf.size();
    //  for (int i = 0; i < std::min(weight_nnpdf_size_, MAX_SIZE_PDF_); ++i) {
    //    weights_nnpdf_[i] = zf_event.weights_nnpdf[i];
    //  }

    //  weight_fsr_=zf_event.weight_fsr;
    //}


    // Reco
    reco_z_.z_m = zfe.reco_z.m;
    reco_z_.z_pt = zfe.reco_z.pt;
    reco_z_.z_y = zfe.reco_z.y;
    reco_z_.z_phi = zfe.reco_z.phi;
    reco_z_.z_phistar = zfe.reco_z.phistar;
    reco_z_.z_eta = zfe.reco_z.eta;
    reco_z_.z_vtx_prob = zfe.reco_z.vtx_prob;
    reco_z_.z_vtx_x = zfe.reco_z.vtx_x;
    reco_z_.z_vtx_y = zfe.reco_z.vtx_y;
    reco_z_.z_vtx_z = zfe.reco_z.vtx_z;

    if (zfe.e0 != NULL && zfe.e1 != NULL){
      reco_z_.daughter0_pt = zfe.e0->pt;
      reco_z_.daughter0_eta = zfe.e0->eta;
      reco_z_.daughter0_phi = zfe.e0->phi;
      
      reco_z_.daughter1_pt = zfe.e1->pt;
      reco_z_.daughter1_eta = zfe.e1->eta;
      reco_z_.daughter1_phi = zfe.e1->phi;

      reco_z_.daughter0_charge = zfe.e0->charge;
      reco_z_.daughter1_charge = zfe.e1->charge;
    }

    reco_z_from_muons_.z_m = zfe.reco_z_from_muons.m;
    reco_z_from_muons_.z_pt = zfe.reco_z_from_muons.pt;
    reco_z_from_muons_.z_y = zfe.reco_z_from_muons.y;
    reco_z_from_muons_.z_phi = zfe.reco_z_from_muons.phi;
    reco_z_from_muons_.z_phistar = zfe.reco_z_from_muons.phistar;
    reco_z_from_muons_.z_eta = zfe.reco_z_from_muons.eta;
    reco_z_from_muons_.z_vtx_prob = zfe.reco_z_from_muons.vtx_prob;
    reco_z_from_muons_.z_vtx_x = zfe.reco_z_from_muons.vtx_x;
    reco_z_from_muons_.z_vtx_y = zfe.reco_z_from_muons.vtx_y;
    reco_z_from_muons_.z_vtx_z = zfe.reco_z_from_muons.vtx_z;

    if (zfe.reco_z_from_muons.m > -1) {
      reco_z_from_muons_.daughter0_pt = zfe.z_muon0.pt();
      reco_z_from_muons_.daughter0_eta = zfe.z_muon0.eta();
      reco_z_from_muons_.daughter0_phi = zfe.z_muon0.phi();
      reco_z_from_muons_.daughter1_pt = zfe.z_muon1.pt();
      reco_z_from_muons_.daughter1_eta = zfe.z_muon1.eta();
      reco_z_from_muons_.daughter1_phi = zfe.z_muon1.phi();

      reco_z_from_muons_.daughter0_charge = zfe.z_muon0.charge();
      reco_z_from_muons_.daughter1_charge = zfe.z_muon1.charge();
    }
    
    
    //if (zfe.found_jpsi) {}
    if (true) {
      int n_jpsi = 0;
      for (unsigned int i = 0; i < zfe.reco_jpsi.m.size() ; ++i ) {
        //TODO should fold in eta window cut here??
        //TODO should fold in jpsi pT window cut here??
        bool APPLY_MUON_MIN_PT_;
        bool APPLY_SOFT_MUONS_;
        bool APPLY_JPSI_MASS_WINDOW_;
        bool APPLY_VERTEX_Z_POS_WINDOW_;
        bool APPLY_DIMUON_VTX_COMPATIBILITY_;

        if(zfe.found_jpsi) {
          APPLY_MUON_MIN_PT_ = true;
          APPLY_SOFT_MUONS_ = true;
          APPLY_JPSI_MASS_WINDOW_ = true;
          APPLY_VERTEX_Z_POS_WINDOW_ = true;
          APPLY_DIMUON_VTX_COMPATIBILITY_ = true;
        }
        else {
          APPLY_MUON_MIN_PT_ = false;
          APPLY_SOFT_MUONS_ = false;
          APPLY_JPSI_MASS_WINDOW_ = false;
          APPLY_VERTEX_Z_POS_WINDOW_ = false;
          APPLY_DIMUON_VTX_COMPATIBILITY_ = false;
        }


        if (APPLY_MUON_MIN_PT_ && (!zfe.reco_jpsi.has_high_pt_muons.at(i) || !zfe.reco_jpsi.has_muons_in_eta_window.at(i) || 
              !zfe.reco_jpsi.is_high_pt.at(i) ) )  {
          continue;
        }
        if (APPLY_SOFT_MUONS_ && !zfe.reco_jpsi.has_soft_id_muons.at(i) ) {
          continue;
        }
        if (APPLY_JPSI_MASS_WINDOW_ && !zfe.reco_jpsi.is_within_jpsi_mass_window.at(i) ) {
          continue;
        }
        //TODO should this cut be relative to the z position??
        if (APPLY_VERTEX_Z_POS_WINDOW_ && !zfe.reco_jpsi.has_dimuon_vertex_compatible_with_primary_vertex.at(i) ) {
          continue;
        }
        if (APPLY_DIMUON_VTX_COMPATIBILITY_ && !zfe.reco_jpsi.has_muons_with_compatible_vertex.at(i) ) {
          continue;
        }
        n_jpsi++;
        //TODO decide on how to do this for a tree (and histogram method too), for now just take first jpsi
        //this is definitely a kludge, think of proper way to do this
        if (n_jpsi > 1) {
          continue;
        }

        reco_jpsi_.jpsi_m = zfe.reco_jpsi.m.at(i);
        reco_jpsi_.jpsi_pt = zfe.reco_jpsi.pt.at(i);
        reco_jpsi_.jpsi_y = zfe.reco_jpsi.y.at(i);
        reco_jpsi_.jpsi_phi = zfe.reco_jpsi.phi.at(i);
        reco_jpsi_.jpsi_eta = zfe.reco_jpsi.eta.at(i);

        reco_jpsi_.jpsi_tau_xy = zfe.reco_jpsi.tau_xy.at(i);
        reco_jpsi_.jpsi_tau_z = zfe.reco_jpsi.tau_z.at(i);

        reco_jpsi_.jpsi_distance_xy = zfe.reco_jpsi.distance_xy.at(i);
        reco_jpsi_.jpsi_distance_z = zfe.reco_jpsi.distance_z.at(i);

        reco_jpsi_.jpsi_vtx_prob = zfe.reco_jpsi.y.at(i);
        reco_jpsi_.jpsi_vtx_x = zfe.reco_jpsi.vtx_x.at(i);
        reco_jpsi_.jpsi_vtx_y = zfe.reco_jpsi.vtx_y.at(i);
        reco_jpsi_.jpsi_vtx_z = zfe.reco_jpsi.vtx_z.at(i);

        reco_jpsi_.jpsi_eff = zfe.reco_jpsi.jpsi_efficiency.at(i);
        reco_jpsi_.jpsi_acc_eff = zfe.reco_jpsi.jpsi_acc_eff.at(i);
        reco_jpsi_.jpsi_scale_factor = zfe.reco_jpsi.jpsi_scale_factor.at(i);

        reco_jpsi_.muon0_pt = zfe.reco_jpsi.muon0.at(i).pt();
        reco_jpsi_.muon0_eta = zfe.reco_jpsi.muon0.at(i).eta();
        reco_jpsi_.muon0_phi = zfe.reco_jpsi.muon0.at(i).phi();

        reco_jpsi_.muon1_pt = zfe.reco_jpsi.muon1.at(i).pt();
        reco_jpsi_.muon1_eta = zfe.reco_jpsi.muon1.at(i).eta();
        reco_jpsi_.muon1_phi = zfe.reco_jpsi.muon1.at(i).phi();

        reco_jpsi_.muon0_charge = zfe.reco_jpsi.muon0.at(i).charge();
        reco_jpsi_.muon1_charge = zfe.reco_jpsi.muon1.at(i).charge();

        reco_jpsi_.has_muons_in_eta_window = zfe.reco_jpsi.has_muons_in_eta_window.at(i);
        reco_jpsi_.has_high_pt_muons = zfe.reco_jpsi.has_high_pt_muons.at(i);
      }
    }


    //TODO jpsi->ee
    //---------------------------------------------------------------------------------------
    if (true) {
      int n_jpsi_from_electrons = 0;
      for (unsigned int i = 0; i < zfe.reco_jpsi_from_electrons.m.size() ; ++i ) {
        //TODO should fold in eta window cut here??
        //TODO should fold in jpsi pT window cut here??
        bool APPLY_MUON_MIN_PT_from_electrons_;
        bool APPLY_SOFT_MUONS_from_electrons_;
        bool APPLY_JPSI_MASS_WINDOW_from_electrons_;
        bool APPLY_VERTEX_Z_POS_WINDOW_from_electrons_;
        bool APPLY_DIMUON_VTX_COMPATIBILITY_from_electrons_;

        if(zfe.found_jpsi_from_electrons) {
          APPLY_MUON_MIN_PT_from_electrons_ = true;
          APPLY_SOFT_MUONS_from_electrons_ = true;
          APPLY_JPSI_MASS_WINDOW_from_electrons_ = true;
          APPLY_VERTEX_Z_POS_WINDOW_from_electrons_ = true;
          APPLY_DIMUON_VTX_COMPATIBILITY_from_electrons_ = true;
        }
        else {
          APPLY_MUON_MIN_PT_from_electrons_ = false;
          APPLY_SOFT_MUONS_from_electrons_ = false;
          APPLY_JPSI_MASS_WINDOW_from_electrons_ = false;
          APPLY_VERTEX_Z_POS_WINDOW_from_electrons_ = false;
          APPLY_DIMUON_VTX_COMPATIBILITY_from_electrons_ = false;
        }

        if (APPLY_MUON_MIN_PT_from_electrons_ && (!zfe.reco_jpsi_from_electrons.has_high_pt_muons.at(i) || !zfe.reco_jpsi_from_electrons.has_muons_in_eta_window.at(i) || 
              !zfe.reco_jpsi_from_electrons.is_high_pt.at(i) ) )  {
          continue;
        }
        if (APPLY_SOFT_MUONS_from_electrons_ && !zfe.reco_jpsi_from_electrons.has_soft_id_muons.at(i) ) {
          continue;
        }
        if (APPLY_JPSI_MASS_WINDOW_from_electrons_ && !zfe.reco_jpsi_from_electrons.is_within_jpsi_mass_window.at(i) ) {
          continue;
        }
        //TODO should this cut be relative to the z position??
        if (APPLY_VERTEX_Z_POS_WINDOW_from_electrons_ && !zfe.reco_jpsi_from_electrons.has_dimuon_vertex_compatible_with_primary_vertex.at(i) ) {
          continue;
        }
        if (APPLY_DIMUON_VTX_COMPATIBILITY_from_electrons_ && !zfe.reco_jpsi_from_electrons.has_muons_with_compatible_vertex.at(i) ) {
          continue;
        }
        n_jpsi_from_electrons++;
        //TODO decide on how to do this for a tree (and histogram method too), for now just take first jpsi
        //this is definitely a kludge, think of proper way to do this
        if (n_jpsi_from_electrons > 1) {
          continue;
        }

        reco_jpsi_from_electrons_.jpsi_m = zfe.reco_jpsi_from_electrons.m.at(i);
        reco_jpsi_from_electrons_.jpsi_pt = zfe.reco_jpsi_from_electrons.pt.at(i);
        reco_jpsi_from_electrons_.jpsi_y = zfe.reco_jpsi_from_electrons.y.at(i);
        reco_jpsi_from_electrons_.jpsi_phi = zfe.reco_jpsi_from_electrons.phi.at(i);
        reco_jpsi_from_electrons_.jpsi_eta = zfe.reco_jpsi_from_electrons.eta.at(i);

        reco_jpsi_from_electrons_.jpsi_tau_xy = zfe.reco_jpsi_from_electrons.tau_xy.at(i);
        reco_jpsi_from_electrons_.jpsi_tau_z = zfe.reco_jpsi_from_electrons.tau_z.at(i);

        reco_jpsi_from_electrons_.jpsi_distance_xy = zfe.reco_jpsi_from_electrons.distance_xy.at(i);
        reco_jpsi_from_electrons_.jpsi_distance_z = zfe.reco_jpsi_from_electrons.distance_z.at(i);

        reco_jpsi_from_electrons_.jpsi_vtx_prob = zfe.reco_jpsi_from_electrons.y.at(i);
        reco_jpsi_from_electrons_.jpsi_vtx_x = zfe.reco_jpsi_from_electrons.vtx_x.at(i);
        reco_jpsi_from_electrons_.jpsi_vtx_y = zfe.reco_jpsi_from_electrons.vtx_y.at(i);
        reco_jpsi_from_electrons_.jpsi_vtx_z = zfe.reco_jpsi_from_electrons.vtx_z.at(i);

        //reco_jpsi_from_electrons_.jpsi_eff = zfe.reco_jpsi.jpsi_efficiency.at(i);
        //reco_jpsi_from_electrons_.jpsi_acc_eff = zfe.reco_jpsi.jpsi_acc_eff.at(i);
        //reco_jpsi_from_electrons_.jpsi_scale_factor = zfe.reco_jpsi.jpsi_scale_factor.at(i);

        //reco_jpsi_.muon0_pt = zfe.reco_jpsi.muon0.at(i).pt();
        //reco_jpsi_.muon0_eta = zfe.reco_jpsi.muon0.at(i).eta();
        //reco_jpsi_.muon0_phi = zfe.reco_jpsi.muon0.at(i).phi();

        //reco_jpsi_.muon1_pt = zfe.reco_jpsi.muon1.at(i).pt();
        //reco_jpsi_.muon1_eta = zfe.reco_jpsi.muon1.at(i).eta();
        //reco_jpsi_.muon1_phi = zfe.reco_jpsi.muon1.at(i).phi();

        //reco_jpsi_.muon0_charge = zfe.reco_jpsi.muon0.at(i).charge();
        //reco_jpsi_.muon1_charge = zfe.reco_jpsi.muon1.at(i).charge();

        reco_jpsi_from_electrons_.has_muons_in_eta_window = zfe.reco_jpsi_from_electrons.has_muons_in_eta_window.at(i);
        reco_jpsi_from_electrons_.has_high_pt_muons = zfe.reco_jpsi_from_electrons.has_high_pt_muons.at(i);
      }
    }
    //------------------------------------------------------------------------------------------




    // Truth
    if (!zfe.is_real_data) {
      int n_jpsi = 0;
      for (unsigned int i = 0; i < zfe.truth_jpsi.m.size() ; ++i ) {
        n_jpsi++;
        //TODO decide on how to do this for a tree (and histogram method too), for now just take first jpsi
        //this is definitely a kludge, think of proper way to do this truth should only ever have 1 jpsi
        if (n_jpsi > 1) {
          continue;
        }
        truth_jpsi_.jpsi_m = zfe.truth_jpsi.m.at(i);
        truth_jpsi_.jpsi_pt = zfe.truth_jpsi.pt.at(i);
        truth_jpsi_.jpsi_y = zfe.truth_jpsi.y.at(i);
        truth_jpsi_.jpsi_phi = zfe.truth_jpsi.phi.at(i);
        truth_jpsi_.jpsi_eta = zfe.truth_jpsi.eta.at(i);

        //truth_jpsi_.jpsi_distance_xy = zfe.truth_jpsi.distance_xy.at(i);
        //truth_jpsi_.jpsi_distance_z = zfe.truth_jpsi.distance_z.at(i);

        //truth_jpsi_.jpsi_vtx_prob = zfe.truth_jpsi.y.at(i);
        truth_jpsi_.jpsi_vtx_x = zfe.truth_jpsi.vtx_x.at(i);
        truth_jpsi_.jpsi_vtx_y = zfe.truth_jpsi.vtx_y.at(i);
        truth_jpsi_.jpsi_vtx_z = zfe.truth_jpsi.vtx_z.at(i);

        //truth_jpsi_.jpsi_eff = zfe.truth_jpsi.jpsi_efficiency.at(i);
        //truth_jpsi_.jpsi_acc_eff = zfe.truth_jpsi.jpsi_acc_eff.at(i);
        //truth_jpsi_.jpsi_scale_factor = zfe.truth_jpsi.jpsi_scale_factor.at(i);
        truth_jpsi_.muon0_pt = zfe.jpsi_muon0.at(i)->pt();
        truth_jpsi_.muon0_eta = zfe.jpsi_muon0.at(i)->eta();
        truth_jpsi_.muon0_phi = zfe.jpsi_muon0.at(i)->phi();

        truth_jpsi_.muon1_pt = zfe.jpsi_muon1.at(i)->pt();
        truth_jpsi_.muon1_eta = zfe.jpsi_muon1.at(i)->eta();
        truth_jpsi_.muon1_phi = zfe.jpsi_muon1.at(i)->phi();

        truth_jpsi_.muon0_charge = zfe.jpsi_muon0.at(i)->charge();
        truth_jpsi_.muon1_charge = zfe.jpsi_muon1.at(i)->charge();

        truth_jpsi_.has_muons_in_eta_window = zfe.truth_jpsi.has_muons_in_eta_window.at(i);
        truth_jpsi_.has_high_pt_muons = zfe.truth_jpsi.has_high_pt_muons.at(i);
      }

      truth_z_electrons_.z_m = zfe.truth_z_electrons.m;
      truth_z_electrons_.z_pt = zfe.truth_z_electrons.pt;
      truth_z_electrons_.z_y = zfe.truth_z_electrons.y;
      truth_z_electrons_.z_phistar = zfe.truth_z_electrons.phistar;
      truth_z_electrons_.z_eta = zfe.truth_z_electrons.eta;

      truth_z_muons_.z_m = zfe.truth_z_muons.m;
      truth_z_muons_.z_pt = zfe.truth_z_muons.pt;
      truth_z_muons_.z_y = zfe.truth_z_muons.y;
      truth_z_muons_.z_phistar = zfe.truth_z_muons.phistar;
      truth_z_muons_.z_eta = zfe.truth_z_muons.eta;
    }
    // Truth
    //if (IS_MC_ && !zf_event.is_real_data) {
    //  truth_.z_m = zf_event.truth_z.m;
    //  truth_.z_y = zf_event.truth_z.y;
    //  truth_.z_phistar_dressed = zf_event.truth_z.phistar;
    //  truth_.z_phistar_born = zf_event.truth_z.bornPhistar;
    //  truth_.z_phistar_naked = zf_event.truth_z.nakedPhistar;
    //  truth_.z_phistar_sc = zf_event.truth_z.scPhistar;
    //  truth_.z_pt = zf_event.truth_z.pt;
    //  truth_.z_eta = zf_event.truth_z.eta;
    //  truth_.n_verts = zf_event.truth_vert.num;
    //  truth_.n_true_pileup = zf_event.truth_vert.true_num;
    //  if (zf_event.e0_truth != nullptr) {
    //    truth_.e_pt[0] = zf_event.e0_truth->pt();
    //    truth_.e_eta[0] = zf_event.e0_truth->eta();
    //    truth_.e_phi[0] = zf_event.e0_truth->phi();
    //    truth_.e_rnine[0] = zf_event.e0_truth->r9();
    //    truth_.e_charge[0] = zf_event.e0_truth->charge();
    //  }
    //  if (zf_event.e1_truth != nullptr) {
    //    truth_.e_pt[1] = zf_event.e1_truth->pt();
    //    truth_.e_eta[1] = zf_event.e1_truth->eta();
    //    truth_.e_phi[1] = zf_event.e1_truth->phi();
    //    truth_.e_rnine[1] = zf_event.e1_truth->r9();
    //    truth_.e_charge[1] = zf_event.e1_truth->charge();
    //    if (zf_event.GetZDef(zdef_name_) != nullptr) {
    //      const cutlevel_vector* clv = zf_event.GetZDef(zdef_name_);
    //      truth_.t0tight = clv->back().second.t0p1_pass;
    //      truth_.t1tight = clv->back().second.t1p0_pass;
    //    }
    //  }
    //}

    // General Event info

    event_.event_weight = zfe.event_weight;
    event_.event_number = zfe.id.event_num;
    event_.run_number = zfe.id.run_num;

    event_.n_verts = zfe.reco_vert.num;
    event_.truth_n_verts = zfe.truth_vert.num;
    event_.is_mc = !zfe.is_real_data;

    event_.found_high_pt_muons_from_z = zfe.found_high_pt_muons_from_z;
    event_.found_good_muons_from_z =  zfe.found_good_muons_from_z;
    event_.found_dimuon_z_compatible_vertex = zfe.found_dimuon_z_compatible_vertex;
    event_.found_z_to_muons_mass = zfe.found_z_to_muons_mass;

    event_.found_high_pt_electrons_from_z = zfe.found_high_pt_electrons_from_z;
    event_.found_good_electrons_from_z = zfe.found_good_electrons_from_z;
    event_.found_dielectron_z_compatible_vertex = zfe.found_dielectron_z_compatible_vertex;
    event_.found_z_to_electrons_mass = zfe.found_z_to_electrons_mass;

    event_.found_dimuon_jpsi_with_muons_in_eta_window = zfe.found_dimuon_jpsi_with_muons_in_eta_window;
    event_.found_dimuon_jpsi_with_high_pt_muons = zfe.found_dimuon_jpsi_with_high_pt_muons;
    event_.found_dimuon_jpsi_with_soft_id_and_high_pt_muons = zfe.found_dimuon_jpsi_with_soft_id_and_high_pt_muons;
    event_.found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex = zfe.found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex;
    event_.found_good_dimuon_jpsi_compatible_with_primary_vertex = zfe.found_good_dimuon_jpsi_compatible_with_primary_vertex;
    event_.found_jpsi = zfe.found_jpsi;


    //TODO jpsi->ee
    //------------
    event_.found_dimuon_jpsi_from_electrons_with_muons_in_eta_window = zfe.found_dimuon_jpsi_from_electrons_with_muons_in_eta_window;
    event_.found_dimuon_jpsi_from_electrons_with_high_pt_muons = zfe.found_dimuon_jpsi_from_electrons_with_high_pt_muons;
    event_.found_dimuon_jpsi_from_electrons_with_soft_id_and_high_pt_muons = zfe.found_dimuon_jpsi_from_electrons_with_soft_id_and_high_pt_muons;
    event_.found_dimuon_jpsi_from_electrons_with_good_muons_and_compatible_muon_vertex = zfe.found_dimuon_jpsi_from_electrons_with_good_muons_and_compatible_muon_vertex;
    event_.found_good_dimuon_jpsi_from_electrons_compatible_with_primary_vertex = zfe.found_good_dimuon_jpsi_from_electrons_compatible_with_primary_vertex;
    event_.found_jpsi_from_electrons = zfe.found_jpsi_from_electrons;
    //------------

    // Fill if there is a good Z in either truth or reco
    //if (zf_event.truth_z.m > -1 || zf_event.reco_z.m > -1) {}
    //if ((zfe.found_z_to_muons || zfe.found_z_to_electrons) && zfe.found_jpsi) {
    //}
    //if (true) {


    //if (( (zfe.found_high_pt_muons_from_z && zfe.found_good_muons_from_z && zfe.found_dimuon_z_compatible_vertex && zfe.found_z_to_muons_mass)  || 
    //        (zfe.found_high_pt_electrons_from_z && zfe.found_good_electrons_from_z && zfe.found_dielectron_z_compatible_vertex && zfe.found_z_to_electrons_mass)) 
    //    && (zfe.found_jpsi && zfe.found_good_dimuon_jpsi_compatible_with_primary_vertex)) {
    if (true) {
      tree_->Fill();
      //tree_->Write();
    }
  }

  //void ZFinderTree::FillCutWeights(cutlevel_vector const * const CUT_LEVEL_VECTOR) {
  //  if (CUT_LEVEL_VECTOR != nullptr) {
  //    // Check if we pass via Tag 0 Probe 1, or Tag 1 Probe 1
  //    CutLevel last_cutlevel = CUT_LEVEL_VECTOR->back().second;
  //    bool t0p1 = last_cutlevel.t0p1_pass;  // Takes Priority
  //    bool t1p0 = last_cutlevel.t1p0_pass;

  //    // Loop over our vector and at each level pull out the right
  //    for (auto& i_cutlevel : *CUT_LEVEL_VECTOR) {
  //      // Get the cut names and try to find the matching WeightIDs
  //      const std::string TAG_CUT = i_cutlevel.second.tag_cut;
  //      const std::string PROBE_CUT = i_cutlevel.second.probe_cut;
  //      std::map<std::string, int>::const_iterator tag_it = STR_TO_WEIGHTID.find(TAG_CUT);
  //      std::map<std::string, int>::const_iterator probe_it = STR_TO_WEIGHTID.find(PROBE_CUT);
  //      // Get the weights for the tag and probe
  //      double tag_weight = 1.;
  //      double probe_weight = 1.;
  //      if (t0p1) {
  //        tag_weight = i_cutlevel.second.t0p1_tag_weight;
  //        probe_weight = i_cutlevel.second.t0p1_probe_weight;
  //      }
  //      else if (t1p0) {
  //        tag_weight = i_cutlevel.second.t1p0_tag_weight;
  //        probe_weight = i_cutlevel.second.t1p0_probe_weight;
  //      }
  //      // The Tag has a WeightID
  //      if (tag_it != STR_TO_WEIGHTID.end()) {
  //        const int WEIGHTID = tag_it->second;
  //        weight_id_vector_.push_back(std::make_pair(WEIGHTID, tag_weight));
  //      }
  //      if (probe_it != STR_TO_WEIGHTID.end()) {
  //        const int WEIGHTID = probe_it->second;
  //        weight_id_vector_.push_back(std::make_pair(WEIGHTID, probe_weight));
  //      }
  //    }
  //  }
  //}

  TFile* ZFinderTree::GetCurrentFile() {
    return tree_->GetCurrentFile();
  }
}  // namespace zf
