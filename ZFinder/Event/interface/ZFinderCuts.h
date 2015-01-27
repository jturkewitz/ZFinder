#ifndef ZFINDER_ZFINDERCUTS_H_
#define ZFINDER_ZFINDERCUTS_H_

namespace zf {
  //TODO decide on good values for cut levels
  const double MIN_JPSI_MASS = 3.0; //GeV
  const double MAX_JPSI_MASS = 3.2; //GeV
  const double MIN_Z_MASS = 60.0; //GeV
  const double MAX_Z_MASS = 120.0; //GeV
  //const double MIN_JPSI_LEADING_MUON_PT = 3.5; //GeV
  //const double MIN_JPSI_SUBLEADING_MUON_PT = 3.5; //GeV
  const double MIN_JPSI_LEADING_MUON_PT = 5.5; //GeV
  const double MIN_JPSI_SUBLEADING_MUON_PT = 3.5; //GeV
  //TODO testing
  //const double MIN_JPSI_LEADING_MUON_PT = 2.0; //GeV
  //const double MIN_JPSI_SUBLEADING_MUON_PT = 2.0; //GeV
  //TODO testing
  //const double MIN_JPSI_LEADING_MUON_PT = 5.5; //GeV
  //const double MIN_JPSI_SUBLEADING_MUON_PT = 4.0; //GeV
  //const double MIN_JPSI_PT = 7.0; //GeV
  const double MIN_JPSI_PT = 8.0; //GeV
  const double MIN_Z_MUON_PT = 20.0; //GeV
  const double MIN_ELECTRON_PT = 20.0; //GeV
  //const double MAX_JPSI_VERTEX_Z_DISPLACEMENT = 1.0 ; //cm
  const double MAX_JPSI_VERTEX_Z_DISPLACEMENT = 1.0 ; //cm
  const double MIN_VERTEX_PROB = 0.005; //
  const double MAX_JPSI_MUON_ETA = 2.1; //
  //const double MAX_DELTAR_TRUTH_MATCHED_JPSI_MUONS = 0.01; 
  const double MAX_DELTAR_TRUTH_MATCHED_JPSI_MUONS = 0.015; 
  const double MIN_DELTAR_DISTINCT_Z_JPSI_MUONS = 0.015; 
}
#endif  // ZFINDER_ZFINDERCUTS_H_
