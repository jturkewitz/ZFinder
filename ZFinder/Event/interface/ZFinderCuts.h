#ifndef ZFINDER_ZFINDERCUTS_H_
#define ZFINDER_ZFINDERCUTS_H_

namespace zf {
  //From muon pog: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs
  //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPOG
  const double MIN_JPSI_MASS = 2.85; //GeV
  const double MAX_JPSI_MASS = 3.35; //GeV
  const double MIN_Z_MASS = 60.0; //GeV
  const double MAX_Z_MASS = 120.0; //GeV
  const double MIN_JPSI_LEADING_MUON_PT = 3.5; //GeV
  const double MIN_JPSI_SUBLEADING_MUON_PT = 3.5; //GeV
  const double MIN_JPSI_PT = 8.5; //GeV
  const double MIN_PROMPT_JPSI_TIME = -0.0003; // ns
  const double MAX_PROMPT_JPSI_TIME = 0.0003; // ns
  const double MIN_Z_MUON_PT = 20.0; //GeV
  const double MIN_ELECTRON_PT = 20.0; //GeV
  const double MAX_JPSI_VERTEX_Z_DISPLACEMENT = 1.0 ; //cm
  const double MIN_VERTEX_PROB = 0.005; //
  const double MAX_JPSI_MUON_ETA = 2.1; //limited by tag and probe method to measure scale factors
  const double MAX_JPSI_RAP = 2.1; //
  const double MAX_DELTAR_TRUTH_MATCHED_JPSI_MUONS = 0.015;
  const double MIN_DELTAR_DISTINCT_Z_JPSI_MUONS = 0.015;
}
#endif  // ZFINDER_ZFINDERCUTS_H_
