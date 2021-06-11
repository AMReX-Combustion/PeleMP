#include <PeleLM.H>
#include <pelelm_prob.H>

std::string
read_data_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

void
read_data(const std::string& myfile)
{
  std::string firstline;
  std::string remaininglines;
  std::string colname;
  unsigned int pos1;
  unsigned int pos2;
  int variable_count;
  int line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_data_file(infile);
  if (!amrex::FileSystem::Exists(myfile))
    amrex::Abort("Data file does not exist");
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  std::stringstream ss(firstline);
  variable_count = 0;
  std::vector<std::string> var_names;
  while (std::getline(ss, colname, ' ')) {
    var_names.push_back(colname);
    variable_count++;
  }

  amrex::Vector<std::string> spec_names;
  pele::physics::eos::speciesNames(spec_names);
  int Ycol = -1;
  int Tcol = -1;
  int Zcol = -1;
  for (int i = 0; i < variable_count; i++) {
    std::string cur_name = var_names[i];
    if (cur_name == "Z" || cur_name == "z")
      Zcol = i;
    else if (cur_name == "T" || cur_name == "temp")
      Tcol = i;
    else if (cur_name != spec_names[i - 2]) {
      Print() << cur_name << " " << spec_names[i-2] << std::endl;
      amrex::Abort("Variables do not match");
    }
  }
  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  if (Zcol != 0 || Tcol != 1) amrex::Abort("Columns must be ordered Z, T, Y");
  amrex::Print() << line_count << " lines found in data file" << std::endl;
  PeleLM::prob_parm->N_equil = line_count;
  std::vector<amrex::Real> T_vals(line_count);
  std::vector<amrex::Real> Z_vals(line_count);
  std::vector<amrex::Real> Y_vals(line_count * NUM_SPECIES);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  for (unsigned int i = 0; i < line_count; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> Z_vals[i];
    sinput >> T_vals[i];
    for (unsigned int j = 0; j < NUM_SPECIES; j++) {
      sinput >> Y_vals[j * line_count + i];
    }
  }

  // Renormalize the mass fractions so they sum to 1
  for (unsigned int i = 0; i < line_count; i++) {
    amrex::Real sumY = 0.;
    for (int n = 0; n < NUM_SPECIES; n++) {
      amrex::Real Yval = Y_vals[n * line_count + i];
      Yval = std::min(1., std::max(Yval, 0.));
      sumY += Yval;
    }
    sumY = 1. / sumY;
    for (int n = 0; n < NUM_SPECIES; n++) {
      amrex::Real Yval = Y_vals[n * line_count + i];
      Yval = std::min(1., std::max(Yval, 0.));
      Y_vals[n * line_count + i] = Yval * sumY;
    }
  }

  PeleLM::prob_parm->Z_equil = (amrex::Real *) The_Arena()->alloc(line_count * sizeof(amrex::Real));
  PeleLM::prob_parm->T_equil = (amrex::Real *) The_Arena()->alloc(line_count * sizeof(amrex::Real));
  PeleLM::prob_parm->Y_equil = (amrex::Real *) The_Arena()->alloc(line_count * NUM_SPECIES * sizeof(amrex::Real));
  for (unsigned int i = 0; i < line_count; i++) {
    PeleLM::prob_parm->Z_equil[i] = Z_vals[i];
    PeleLM::prob_parm->T_equil[i] = T_vals[i];
    for (unsigned int j = 0; j < NUM_SPECIES; ++j) {
      const int k = j * line_count + i;
      PeleLM::prob_parm->Y_equil[k] = Y_vals[k];
    }
  }
}

extern "C" {
void
amrex_probinit(
  const int* init,
  const int* name,
  const int* namelen,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("ref_p", PeleLM::prob_parm->P_mean);
  pp.query("fuel_temp", PeleLM::prob_parm->fuel_T);
  pp.query("oxid_temp", PeleLM::prob_parm->oxid_T);
  pp.query("ref_temp", PeleLM::prob_parm->T_ref);
  pp.query("init_vel", PeleLM::prob_parm->init_vel);
  std::string fuel_name = "C2H4";
  pp.query("fuel_name", fuel_name);
  pp.query("fuel_dia", PeleLM::prob_parm->fuel_dia);
  pp.query("oxid_dia", PeleLM::prob_parm->oxid_dia);
  pp.query("fuel_vel", PeleLM::prob_parm->fuel_vel);
  pp.query("oxid_vel", PeleLM::prob_parm->oxid_vel);
  pp.query("oxid_N2", PeleLM::prob_parm->oxid_N2);
  pp.query("oxid_O2", PeleLM::prob_parm->oxid_O2);
  pp.query("smooth_factor", PeleLM::prob_parm->smooth_b);

  pp.query("jet_center", PeleLM::prob_parm->jet_center);
  amrex::Vector<std::string> spec_names(NUM_SPECIES);
  pele::physics::eos::speciesNames(spec_names);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    if (spec_names[n] == fuel_name) {
      PeleLM::prob_parm->fuel_indx = n;
    }
  }
  if (PeleLM::prob_parm->fuel_indx < 0)
    amrex::Abort("Fuel species not found from reaction data");
  auto eos = pele::physics::PhysicsType::eos();
  {
    int ecompCHON[NUM_SPECIES * 4];
    pele::physics::eos::element_compositionCHON(ecompCHON);
    int numc = ecompCHON[PeleLM::prob_parm->fuel_indx * 4];
    int numh = ecompCHON[PeleLM::prob_parm->fuel_indx * 4 + 1];
    amrex::Real m = (Real)numc;
    amrex::Real n = (Real)numh;
    amrex::Real s = 32.*(m + n / 4.)/(12. * m + n);
    PeleLM::prob_parm->S_equil = s / PeleLM::prob_parm->oxid_O2;
  }

  std::string filename;
  pp.get("data_file", filename);
  read_data(filename);

  amrex::Real moments[NUM_SOOT_MOMENTS + 1];
  SootData* const sd = PeleLM::soot_model->getSootData();
  sd->initialSmallMomVals(moments);
  for (int n = 0; n < NUM_SOOT_MOMENTS + 1; ++n) {
    PeleLM::prob_parm->soot_vals[n] = moments[n];
  }
}
}
