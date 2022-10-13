#include "prob.H"

std::string
read_inject_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

void
read_inject(const std::string myfile)
{
  std::string firstline, remaininglines;
  int line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_inject_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in inject file"
                 << std::endl;

  PeleC::prob_parm_host->inject_N = line_count;
  PeleC::prob_parm_host->inject_time.resize(PeleC::prob_parm_host->inject_N);
  PeleC::prob_parm_host->inject_mass.resize(PeleC::prob_parm_host->inject_N);
  PeleC::prob_parm_host->inject_vel.resize(PeleC::prob_parm_host->inject_N);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  for (int i = 0; i < PeleC::prob_parm_host->inject_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> PeleC::prob_parm_host->inject_time[i];
    sinput >> PeleC::prob_parm_host->inject_mass[i];
    sinput >> PeleC::prob_parm_host->inject_vel[i];
  }
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("inject_file", PeleC::prob_parm_host->input_file);
  pp.query("init_v", PeleC::h_prob_parm_device->v0);
  pp.get("ref_p", PeleC::h_prob_parm_device->p0);
  pp.get("ref_T", PeleC::h_prob_parm_device->T0);
  pp.query("init_N2", PeleC::h_prob_parm_device->Y_N2);
  pp.query("init_O2", PeleC::h_prob_parm_device->Y_O2);

  // Initial density, velocity, and material properties
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[N2_ID] = PeleC::h_prob_parm_device->Y_N2;
  massfrac[O2_ID] = PeleC::h_prob_parm_device->Y_O2;
  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p0, massfrac, PeleC::h_prob_parm_device->T0,
    PeleC::h_prob_parm_device->rho0, eint);
  eos.RTY2Cs(
    PeleC::h_prob_parm_device->rho0, PeleC::h_prob_parm_device->T0, massfrac,
    cs);
  eos.TY2Cp(PeleC::h_prob_parm_device->T0, massfrac, cp);
  // Read injection ROI file if specified
  if (PeleC::prob_parm_host->input_file != "") {
    read_inject(PeleC::prob_parm_host->input_file);
  }
}
}

void
pc_prob_close()
{
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
