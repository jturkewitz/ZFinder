#ifndef ZFINDER_ZFINDERCUTS_H_
#define ZFINDER_ZFINDERCUTS_H_

namespace zf {
  const double MIN_JPSI_MASS = 3.0; //GeV
  const double MAX_JPSI_MASS = 3.2; //GeV
  const double MIN_Z_MASS = 60.0; //GeV
  const double MAX_Z_MASS = 120.0; //GeV
  const double MIN_MUON_PT = 5.0; //GeV
  const double MIN_ELECTRON_PT = 20.0; //GeV
  const double MAX_JPSI_VERTEX_Z_DISPLACEMENT = 1.0 ; //cm
  const double MIN_VERTEX_PROB = 0.005; //from double jpsi paper (I think - TODO verify this)
}
#endif  // ZFINDER_ZFINDERCUTS_H_
