#ifndef ZFINDER_JPSIEFFICIENCYTABLES_H_
#define ZFINDER_JPSIEFFICIENCYTABLES_H_

namespace zf {
  //From muon pog: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs
  // use pt<20 GeV from J/psi tag and probe, pt>20 GeV from Z tag and probe
  //eta min, eta max, pt min, pt max, efficiency, efficiency low error, efficiency high error
  const int SOFT_MUON_DATA_EFF_TABLE_ROWS = 61;
  const double SOFT_MUON_DATA_EFF_TABLE[SOFT_MUON_DATA_EFF_TABLE_ROWS][7] = {
    {0.0,  0.9,  3.5,   3.75,   0.811856170528,   0.0112954910528,   0.0114656571218   },
    {0.0,  0.9,  3.75,  4.0,    0.881767345824,   0.012245801121,    0.0124237690768   },
    {0.0,  0.9,  4.0,   4.5,    0.943709980278,   0.00940750516001,  0.00948700810439  },
    {0.0,  0.9,  4.5,   5.0,    0.974017575705,   0.0106858312862,   0.0107839463275   },
    {0.0,  0.9,  5.0,   6.0,    0.990266366625,   0.00873550095391,  0.0087836938685   },
    {0.0,  0.9,  6.0,   8.0,    0.98746674202,    0.00828241108268,  0.00825866283577  },
    {0.0,  0.9,  8.0,   10.0,   0.987973499684,   0.0076974955487,   0.00777505459907  },
    {0.0,  0.9,  10.0,  15.0,   0.99108118942,    0.00839860892911,  0.00847898277457  },
    {0.0,  0.9,  15.0,  20.0,   0.974052166395,   0.0198718692077,   0.0200955536036   },
    {0.0,  0.9,  20.0,  25.0,   0.938016073206,   0.00143755876021,  0.00143660099368  },
    {0.0,  0.9,  25.0,  30.0,   0.940760716847,   0.000594856884316, 0.000593995932486 },
    {0.0,  0.9,  30.0,  35.0,   0.942722733581,   0.000384301241971, 0.000383567346104 },
    {0.0,  0.9,  35.0,  40.0,   0.943334836206,   0.000284631050334, 0.00028426231427  },
    {0.0,  0.9,  40.0,  50.0,   0.943513041248,   0.000180760146635, 0.000180536081935 },
    {0.0,  0.9,  50.0,  60.0,   0.944585154364,   0.000462199163049, 0.000461019677484 },
    {0.0,  0.9,  60.0,  90.0,   0.944095588202,   0.000788330835306, 0.000786297327674 },
    {0.0,  0.9,  90.0,  140.0,  0.95690888425,    0.00256660044354,  0.00255966353756  },
    {0.0,  0.9,  140.0, 300.0,  0.949071618981,   0.0146616490715,   0.0149181273786   },
    {0.9,  1.2,  2.5,   2.75,   0.148668208996,   0.0132791072103,   0.0142194750875   },
    {0.9,  1.2,  2.75,  3.0,    0.248426179686,   0.013721595176,    0.0144357130192   },
    {0.9,  1.2,  3.0,   3.25,   0.392190415602,   0.0149953845681,   0.0155629019732   },
    {0.9,  1.2,  3.25,  3.5,    0.61906542648,    0.0195478735288,   0.020271596778    },
    {0.9,  1.2,  3.5,   3.75,   0.840082298909,   0.0240619492506,   0.024980450526    },
    {0.9,  1.2,  3.75,  4.0,    0.915805818544,   0.0266319550521,   0.0266319550521   },
    {0.9,  1.2,  4.0,   4.5,    0.953071316808,   0.0209764049986,   0.0215989785467   },
    {0.9,  1.2,  4.5,   5.0,    0.952562425234,   0.0231799379694,   0.0239014512613   },
    {0.9,  1.2,  5.0,   6.0,    0.983103773768,   0.0185193501533,   0.0168962262317   },
    {0.9,  1.2,  6.0,   8.0,    0.999999984657,   0.0158243522204,   1.5343295634e-08  },
    {0.9,  1.2,  8.0,   10.0,   0.972759293054,   0.0180116461328,   0.0183363698581   },
    {0.9,  1.2,  10.0,  20.0,   0.958498959696,   0.0173158866519,   0.0176101034266   },
    {0.9,  1.2,  20.0,  25.0,   0.945671706945,   0.00220550596817,  0.00220204027635  },
    {0.9,  1.2,  25.0,  30.0,   0.94579830735,    0.00106706034241,  0.00106348077172  },
    {0.9,  1.2,  30.0,  35.0,   0.947399710663,   0.000758298668893, 0.000755579196844 },
    {0.9,  1.2,  35.0,  40.0,   0.946744088016,   0.000542887311769, 0.00054115216541  },
    {0.9,  1.2,  40.0,  50.0,   0.947401516825,   0.000329129938245, 0.000328305625735 },
    {0.9,  1.2,  50.0,  60.0,   0.948048324823,   0.000862036425886, 0.000858061044621 },
    {0.9,  1.2,  60.0,  90.0,   0.950906751625,   0.00148646024291,  0.0014795571817   },
    {0.9,  1.2,  90.0,  140.0,  0.961582520164,   0.00502285078418,  0.00500305435252  },
    {0.9,  1.2,  140.0, 300.0,  0.97556223873,    0.0309430792706,   0.0244377612698   },
    {1.2,  2.1,  2.0,   2.5,    0.788063673924,   0.0147043424413,   0.0151590107844   },
    {1.2,  2.1,  2.5,   2.75,   0.896084278543,   0.0208839221933,   0.0215931491884   },
    {1.2,  2.1,  2.75,  3.0,    0.904630140119,   0.0202011352972,   0.0209742764143   },
    {1.2,  2.1,  3.0,   3.25,   0.97686992351,    0.0229127893885,   0.02313007649     },
    {1.2,  2.1,  3.25,  3.5,    0.958662030673,   0.0214885545447,   0.0224660977611   },
    {1.2,  2.1,  3.5,   3.75,   0.962333142748,   0.0226447092518,   0.0233808224058   },
    {1.2,  2.1,  3.75,  4.0,    0.987563311969,   0.0237310784311,   0.0124366880308   },
    {1.2,  2.1,  4.0,   4.5,    0.983945630692,   0.0171539500361,   0.0160543693078   },
    {1.2,  2.1,  4.5,   5.0,    0.973083375624,   0.0178090956938,   0.0182947889492   },
    {1.2,  2.1,  5.0,   6.0,    0.994616576898,   0.015680557913,    0.00538342310241  },
    {1.2,  2.1,  6.0,   8.0,    0.963029661266,   0.0138960839461,   0.0141209048164   },
    {1.2,  2.1,  8.0,   10.0,   0.999999999999,   0.00475923518369,  5.51780843239e-13 },
    {1.2,  2.1,  10.0,  20.0,   0.985272164617,   0.0189464304749,   0.0147278353833   },
    {1.2,  2.1,  20.0,  25.0,   0.970947235338,   0.00113936581869,  0.00113577755535  },
    {1.2,  2.1,  25.0,  30.0,   0.968019925997,   0.000555180821957, 0.000553714820424 },
    {1.2,  2.1,  30.0,  35.0,   0.96866112088,    0.000412550992355, 0.000411355507302 },
    {1.2,  2.1,  35.0,  40.0,   0.968353136107,   0.00031619116437,  0.000315318329187 },
    {1.2,  2.1,  40.0,  50.0,   0.968408636868,   0.00018878956112,  0.000188419587359 },
    {1.2,  2.1,  50.0,  60.0,   0.968635160396,   0.000529253890032, 0.000527646297193 },
    {1.2,  2.1,  60.0,  90.0,   0.965968559039,   0.00104465490546,  0.00104154694487  },
    {1.2,  2.1,  90.0,  140.0,  0.993703616522,   0.00421815405176,  0.00419250724001  },
    {1.2,  2.1,  140.0, 300.0,  0.96114824326,    0.026100922782,    0.0271238072942   }
  };

  const double SOFT_MUON_SCALE_FACTOR_TABLE[SOFT_MUON_DATA_EFF_TABLE_ROWS][7] = {
    {0.0,  0.9,  3.5,   3.75,   1.06175988556,    0.0148411237678,   0.0150626552697   },
    {0.0,  0.9,  3.75,  4.0,    1.03839181307,    0.0144666299957,   0.014675277874    },
    {0.0,  0.9,  4.0,   4.5,    1.04257114671,    0.0104182229996,   0.0105058437829   },
    {0.0,  0.9,  4.5,   5.0,    1.03492614746,    0.0113734388698,   0.0114773213542   },
    {0.0,  0.9,  5.0,   6.0,    1.03162166163,    0.00911176523849,  0.00916175895683  },
    {0.0,  0.9,  6.0,   8.0,    1.01993291255,    0.00856369336426,  0.00853908992198  },
    {0.0,  0.9,  8.0,   10.0,   1.01674439921,    0.00797614131407,  0.00805412148744  },
    {0.0,  0.9,  10.0,  15.0,   1.01775152966,    0.00868752966077,  0.0087668140643   },
    {0.0,  0.9,  15.0,  20.0,   0.996885928808,   0.0204816109998,   0.0206931324255   },
    {0.0,  0.9,  20.0,  25.0,   0.982074130575,   0.00170612864269,  0.00170220849432  },
    {0.0,  0.9,  25.0,  30.0,   0.985531372426,   0.000835805814475, 0.000832743363302 },
    {0.0,  0.9,  30.0,  35.0,   0.986203005944,   0.000582892613207, 0.000579320834879 },
    {0.0,  0.9,  35.0,  40.0,   0.98722816383,    0.00045103156208,  0.000448892158719 },
    {0.0,  0.9,  40.0,  50.0,   0.985879162208,   0.00028729842992,  0.000287489603719 },
    {0.0,  0.9,  50.0,  60.0,   0.981736707188,   0.000681866570282, 0.000678923568891 },
    {0.0,  0.9,  60.0,  90.0,   0.973947058198,   0.00103288297786,  0.00102760754367  },
    {0.0,  0.9,  90.0,  140.0,  0.97727539521,    0.00303302682596,  0.00300895375529  },
    {0.0,  0.9,  140.0, 300.0,  0.984921437542,   0.0165778523062,   0.0167264819187   },
    {0.9,  1.2,  2.5,   2.75,   1.28127488444,    0.116155748854,    0.124175962259    },
    {0.9,  1.2,  2.75,  3.0,    1.49431507538,    0.0838829895094,   0.088125700488    },
    {0.9,  1.2,  3.0,   3.25,   1.58141055953,    0.0614996215937,   0.0637563655318   },
    {0.9,  1.2,  3.25,  3.5,    1.18594878529,    0.0377507305157,   0.0391264338333   },
    {0.9,  1.2,  3.5,   3.75,   1.11611378642,    0.0321065571164,   0.0333191996937   },
    {0.9,  1.2,  3.75,  4.0,    1.06403310568,    0.0310210353962,   0.0310200845383   },
    {0.9,  1.2,  4.0,   4.5,    1.04460205883,    0.0230323186219,   0.0237129780957   },
    {0.9,  1.2,  4.5,   5.0,    1.01388427876,    0.0247061049087,   0.0254723909564   },
    {0.9,  1.2,  5.0,   6.0,    1.03018737949,    0.0194294608593,   0.0177303241869   },
    {0.9,  1.2,  6.0,   8.0,    1.03616447631,    0.0164168131938,   0.000804508055488 },
    {0.9,  1.2,  8.0,   10.0,   1.00589645852,    0.0187278037325,   0.019056172089    },
    {0.9,  1.2,  10.0,  20.0,   0.990065233559,   0.0180080791509,   0.018302039577    },
    {0.9,  1.2,  20.0,  25.0,   0.993219848431,   0.00265265918313,  0.00264224906867  },
    {0.9,  1.2,  25.0,  30.0,   0.992860860266,   0.0015339802173,   0.00152380856829  },
    {0.9,  1.2,  30.0,  35.0,   0.994521208456,   0.00118163259904,  0.00117358848581  },
    {0.9,  1.2,  35.0,  40.0,   0.993348050584,   0.000875704816882, 0.000871302667468 },
    {0.9,  1.2,  40.0,  50.0,   0.993530673156,   0.000542316335752, 0.000540421926313 },
    {0.9,  1.2,  50.0,  60.0,   0.992882711196,   0.00135032345436,  0.00134068067333  },
    {0.9,  1.2,  60.0,  90.0,   0.984319473215,   0.00203377543914,  0.00201529405734  },
    {0.9,  1.2,  90.0,  140.0,  0.992636933575,   0.00626539790305,  0.00618103742429  },
    {0.9,  1.2,  140.0, 300.0,  0.999005948461,   0.0336665550366,   0.0271553611489   },
    {1.2,  2.1,  2.0,   2.5,    1.06358769813,    0.0198816464497,   0.0204941296075   },
    {1.2,  2.1,  2.5,   2.75,   1.08257564176,    0.0252688181352,   0.0261242115962   },
    {1.2,  2.1,  2.75,  3.0,    1.03359650812,    0.0231129623349,   0.023994938073    },
    {1.2,  2.1,  3.0,   3.25,   1.06809590016,    0.0250777170426,   0.0253148218038   },
    {1.2,  2.1,  3.25,  3.5,    1.02748075322,    0.023054336101,    0.0241007934019   },
    {1.2,  2.1,  3.5,   3.75,   1.02038251672,    0.0240328345046,   0.0248123354305   },
    {1.2,  2.1,  3.75,  4.0,    1.04120365,       0.0250433003993,   0.0131562721628   },
    {1.2,  2.1,  4.0,   4.5,    1.02999048709,    0.0179736524377,   0.0168235590171   },
    {1.2,  2.1,  4.5,   5.0,    1.01399826285,    0.0185772110819,   0.0190824723584   },
    {1.2,  2.1,  5.0,   6.0,    1.03040406334,    0.0162578053458,   0.00561437509595  },
    {1.2,  2.1,  6.0,   8.0,    0.992502787142,   0.0143325770374,   0.0145639175714   },
    {1.2,  2.1,  8.0,   10.0,   1.02888571576,    0.00514981310476,  0.00155859801463  },
    {1.2,  2.1,  10.0,  20.0,   1.01081753829,    0.019530607445,    0.0152220822886   },
    {1.2,  2.1,  20.0,  25.0,   1.00640890256,    0.00133729534218,  0.00133180603791  },
    {1.2,  2.1,  25.0,  30.0,   1.00087363207,    0.000765754549018, 0.000762243286366 },
    {1.2,  2.1,  30.0,  35.0,   0.999850597269,   0.000604223180681, 0.000601414967594 },
    {1.2,  2.1,  35.0,  40.0,   0.9979503498,     0.000483779942964, 0.000481698409755 },
    {1.2,  2.1,  40.0,  50.0,   0.997282501473,   0.000294047630611, 0.000293165705946 },
    {1.2,  2.1,  50.0,  60.0,   0.994577834542,   0.000752496883996, 0.000748345021647 },
    {1.2,  2.1,  60.0,  90.0,   0.986885875814,   0.00130663307114,  0.00129886088897  },
    {1.2,  2.1,  90.0,  140.0,  1.0095174065,     0.0048054135234,   0.00475828615565  },
    {1.2,  2.1,  140.0, 300.0,  0.969208290895,   0.0281130616917,   0.0285333159178   }
  };

  const int SOFT_MUON_DATA_ACC_EFF_TABLE_ROWS = 49;
  const double SOFT_MUON_DATA_ACC_EFF_TABLE[SOFT_MUON_DATA_EFF_TABLE_ROWS][7] = {
    {-2.1,  -1.5,  8,  8.5,  0.0524149,  0.00122028,  0.00122028   },
    {-2.1,  -1.5,  8.5,  9,  0.0813407,  0.00174599,  0.00174599   },
    {-2.1,  -1.5,  9,  9.5,  0.108138,  0.00227629,  0.00227629   },
    {-2.1,  -1.5,  9.5,  10,  0.13054,  0.00282069,  0.00282069   },
    {-2.1,  -1.5,  10,  12,  0.190683,  0.00223761,  0.00223761   },
    {-2.1,  -1.5,  12,  15,  0.296755,  0.00369543,  0.00369543   },
    {-2.1,  -1.5,  15,  100,  0.458527,  0.00543686,  0.00543686   },
    {-1.5,  -0.9,  8,  8.5,  0.0673722,  0.00123457,  0.00123457   },
    {-1.5,  -0.9,  8.5,  9,  0.10291,  0.00174055,  0.00174055   },
    {-1.5,  -0.9,  9,  9.5,  0.126088,  0.00219989,  0.00219989   },
    {-1.5,  -0.9,  9.5,  10,  0.164135,  0.00280498,  0.00280498   },
    {-1.5,  -0.9,  10,  12,  0.234229,  0.00218765,  0.00218765   },
    {-1.5,  -0.9,  12,  15,  0.376742,  0.00354744,  0.00354744   },
    {-1.5,  -0.9,  15,  100,  0.52832,  0.00499396,  0.00499396   },
    {-0.9,  -0.3,  8,  8.5,  0.0652169,  0.00117112,  0.00117112   },
    {-0.9,  -0.3,  8.5,  9,  0.0993263,  0.00166509,  0.00166509   },
    {-0.9,  -0.3,  9,  9.5,  0.129577,  0.00212884,  0.00212884   },
    {-0.9,  -0.3,  9.5,  10,  0.163615,  0.00275307,  0.00275307   },
    {-0.9,  -0.3,  10,  12,  0.241517,  0.0021267,  0.0021267   },
    {-0.9,  -0.3,  12,  15,  0.382701,  0.00342107,  0.00342107   },
    {-0.9,  -0.3,  15,  100,  0.545537,  0.00474708,  0.00474708   },
    {-0.3,  0.3,  8,  8.5,  0.0645992,  0.0011577,  0.0011577   },
    {-0.3,  0.3,  8.5,  9,  0.100385,  0.00163271,  0.00163271   },
    {-0.3,  0.3,  9,  9.5,  0.134553,  0.00216176,  0.00216176   },
    {-0.3,  0.3,  9.5,  10,  0.158135,  0.00265052,  0.00265052   },
    {-0.3,  0.3,  10,  12,  0.242834,  0.0021135,  0.0021135   },
    {-0.3,  0.3,  12,  15,  0.378013,  0.00337519,  0.00337519   },
    {-0.3,  0.3,  15,  100,  0.534651,  0.0046368,  0.0046368   },
    {0.3,  0.9,  8,  8.5,  0.0657064,  0.00117573,  0.00117573   },
    {0.3,  0.9,  8.5,  9,  0.0943671,  0.00161157,  0.00161157   },
    {0.3,  0.9,  9,  9.5,  0.132424,  0.00216862,  0.00216862   },
    {0.3,  0.9,  9.5,  10,  0.166409,  0.00273751,  0.00273751   },
    {0.3,  0.9,  10,  12,  0.241519,  0.00212504,  0.00212504   },
    {0.3,  0.9,  12,  15,  0.381022,  0.00341388,  0.00341388   },
    {0.3,  0.9,  15,  100,  0.538559,  0.00475284,  0.00475284   },
    {0.9,  1.5,  8,  8.5,  0.066637,  0.00122199,  0.00122199   },
    {0.9,  1.5,  8.5,  9,  0.101144,  0.00172595,  0.00172595   },
    {0.9,  1.5,  9,  9.5,  0.127711,  0.0022096,  0.0022096   },
    {0.9,  1.5,  9.5,  10,  0.175082,  0.00289825,  0.00289825   },
    {0.9,  1.5,  10,  12,  0.230901,  0.00216556,  0.00216556   },
    {0.9,  1.5,  12,  15,  0.3705,  0.00352466,  0.00352466   },
    {0.9,  1.5,  15,  100,  0.508935,  0.0049203,  0.0049203   },
    {1.5,  2.1,  8,  8.5,  0.0547558,  0.0012484,  0.0012484   },
    {1.5,  2.1,  8.5,  9,  0.0815384,  0.00172765,  0.00172765   },
    {1.5,  2.1,  9,  9.5,  0.108196,  0.00227465,  0.00227465   },
    {1.5,  2.1,  9.5,  10,  0.126153,  0.00280314,  0.00280314   },
    {1.5,  2.1,  10,  12,  0.190532,  0.00224661,  0.00224661   },
    {1.5,  2.1,  12,  15,  0.300408,  0.00370019,  0.00370019   },
    {1.5,  2.1,  15,  100,  0.460558,  0.00549688,  0.00549688   }
  };

}
#endif //ZFINDER_JPSIEFFICIENCYTABLES_H_
