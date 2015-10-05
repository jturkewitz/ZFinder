//first bin overall, then 8.5-10, 10-14, 14-18, 18-30,30-100 jpsi pt
//const double ACC_EFF[6] = {0.294805,0.206764,0.333633,0.485107,0.593919,0.721481}; //vertex_comp no primary vert requirement
//double acc_eff_weight[6] = {0.410434,0.296667,0.466977,0.647788,0.753557,0.837575}; //vertex_comp no primary vert requirement
//double acc_eff_weight[6] = {0.236977,0.161808,0.266988,0.403564,0.513895,0.66276}; //vertex_comp no primary vert requirement

//Using eta pt 2.5 for high eta, and eta up to 2.4 (questionable because of how soft muon id was defined)
const double ACC_EFF[6] = {0.374168,0.291977,0.41167,0.549881,0.646403,0.755547}; //vertex_comp no primary vert requirement
const double NUM_ZTOMUMU = 8.481e6;
const double NUM_ZTOEE = 5.274e6;
