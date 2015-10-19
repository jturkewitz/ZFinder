#ifndef ZFINDER_ZFINDERTREE_H_
#define ZFINDER_ZFINDERTREE_H_

// Standard Library
#include <string>  // string
#include <utility>  // pair

// ROOT
#include "TBranch.h"  // TBranch
#include "TDirectory.h"  // TDirectory
#include "TTree.h"  // TTree

// CMSSW
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ZFinder Code
#include "ZFinderEvent.h"  // ZFinderEvent

namespace zf {
  class ZFinderTree {
    public:
      // Constructor
      ZFinderTree(
          TFileDirectory& tdir,
          const bool IS_MC = false
          );

      // destructor
      ~ZFinderTree();

      // Add event
      void Fill(const ZFinderEvent& zf_event);

      // Wrapper around TTree::GetCurrentFile()
      TFile* GetCurrentFile();

    protected:
      // Structs that map to the branches
      struct jpsi_branch {
        void clear_values() {
          jpsi_m = -1000;
          jpsi_pt = -1000;
          jpsi_y = -1000;
          jpsi_phi = -1000;
          jpsi_eta = -1000;
          jpsi_vtx_prob = -1000;
          jpsi_vtx_x = -1000;
          jpsi_vtx_y = -1000;
          jpsi_vtx_z = -1000;

          jpsi_tau_xy = -1000;
          jpsi_tau_z =  -1000;
          jpsi_distance_xy = -1000;
          jpsi_distance_z = -1000;

          jpsi_eff = -1000;
          jpsi_acc_eff = -1000;
          jpsi_scale_factor = -1000;

          muon0_pt = -1000;
          muon0_eta = -1000;
          muon0_phi = -1000;

          muon1_pt = -1000;
          muon1_eta = -1000;
          muon1_phi = -1000;

          muon0_charge = -1000;
          muon1_charge = -1000;

        }
        // Constructor
        jpsi_branch() {
          clear_values();
        }
        double jpsi_m;
        double jpsi_pt;
        double jpsi_y;
        double jpsi_phi;
        double jpsi_eta;

        double jpsi_vtx_prob;
        double jpsi_vtx_x;
        double jpsi_vtx_y;
        double jpsi_vtx_z;

        double jpsi_tau_xy;
        double jpsi_tau_z;
        double jpsi_distance_xy;
        double jpsi_distance_z;


        double jpsi_eff;
        double jpsi_acc_eff;
        double jpsi_scale_factor;

        double muon0_pt;
        double muon0_eta;
        double muon0_phi;

        double muon1_pt;
        double muon1_eta;
        double muon1_phi;

        int muon0_charge;
        int muon1_charge;

      } reco_jpsi_, truth_jpsi_;

      struct z_branch {
        void clear_values() {
          z_m = -1000;
          z_pt = -1000;
          z_y = -1000;
          z_phi = -1000;
          z_phistar = -1000;
          z_eta = -1000;
          z_vtx_prob = -1000;
          z_vtx_x = -1000;
          z_vtx_y = -1000;
          z_vtx_z = -1000;

          daughter0_pt = -1000;
          daughter0_eta = -1000;
          daughter0_phi = -1000;

          daughter1_pt = -1000;
          daughter1_eta = -1000;
          daughter1_phi = -1000;

          daughter0_charge = 0;
          daughter1_charge = 0;
        }
        // Constructor
        z_branch() {
          clear_values();
        }
        double z_m;
        double z_pt;
        double z_y;
        double z_phi;
        double z_phistar;
        double z_eta;
        double z_vtx_prob;
        double z_vtx_x;
        double z_vtx_y;
        double z_vtx_z;

        double daughter0_pt;
        double daughter0_eta;
        double daughter0_phi;

        double daughter1_pt;
        double daughter1_eta;
        double daughter1_phi;

        int daughter0_charge;
        int daughter1_charge;
      } reco_z_, reco_z_from_muons_, truth_z_;

      struct event_branch {
        void clear_values() {
          event_number = 0;
          run_number = 0;
          n_verts = 0;
          found_z_to_electrons = false;
          found_z_to_muons = false;
          found_jpsi = false;
          is_mc = false;
        }

        // Constructor
        event_branch() {
          clear_values();
        }
        unsigned int event_number;
        unsigned int run_number;
        int n_verts;

        bool is_mc;
        bool found_z_to_electrons;
        bool found_z_to_muons;
        bool found_jpsi;
      } event_;

      // File Directory to write to
      TDirectory* tdir_;

      // Use the MC or reco data
      const bool IS_MC_;

      // The tuples
      TTree* tree_;

      //// Set up a variable size branch for the weights
      //int weight_size_;
      //double weight_fsr_;
      //int weight_cteq_size_;
      //int weight_mstw_size_;
      //int weight_nnpdf_size_;
      //static constexpr int MAX_SIZE_ = 100;
      //static constexpr int MAX_SIZE_PDF_ = 110;
      //
      //// Although vectors seem like the right solution, since TTrees need
      //// the memory used for the array to be static, an array is
      //// (unfortunately) the best choice
      //double weights_[MAX_SIZE_];
      //int weight_ids_[MAX_SIZE_];
      //double weights_cteq_[MAX_SIZE_PDF_];
      //double weights_mstw_[MAX_SIZE_PDF_];
      //double weights_nnpdf_[MAX_SIZE_PDF_];
      //
      //// We insert the weights and the IDs into this vector, and then
      //// read it out into the array before filling the tree
      //std::vector<std::pair<int, double>> weight_id_vector_;
      //
      //// Get the weight of the cuts
      //void FillCutWeights(cutlevel_vector const * const CUT_LEVEL_VECTOR);
      double GetTotalWeight(cutlevel_vector const * const CUT_LEVEL_VECTOR);
  };
}  // namespace zf
#endif  // ZFINDER_ZFINDERTREE_H_
