#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 300.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real v0 = 0.;
AMREX_GPU_DEVICE_MANAGED amrex::Real jet_vel = 5000.;
AMREX_GPU_DEVICE_MANAGED amrex::Real jet_dia = 1.E-2;
AMREX_GPU_DEVICE_MANAGED amrex::Real part_mean_dia = 1.E-3;
AMREX_GPU_DEVICE_MANAGED amrex::Real part_stdev_dia = 1.E-4;
AMREX_GPU_DEVICE_MANAGED amrex::Real mass_flow_rate = 2.3;
AMREX_GPU_DEVICE_MANAGED unsigned int inject_N = 0;
amrex::Gpu::ManagedVector<amrex::Real>* inject_time = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* inject_mass = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* inject_vel = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* d_inject_time = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* d_inject_mass = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* d_inject_vel = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real part_rho = 0.693;
AMREX_GPU_DEVICE_MANAGED amrex::Real part_temp = 300.;
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_O2 = 0.233;
AMREX_GPU_DEVICE_MANAGED amrex::Real Y_N2 = 0.767;
AMREX_GPU_DEVICE_MANAGED amrex::Real jet_start_time = 0.;
AMREX_GPU_DEVICE_MANAGED amrex::Real jet_end_time = 10000.;
AMREX_GPU_DEVICE_MANAGED amrex::Real spray_angle = 20.;

std::string input_file = "";
} // namespace ProbParm

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
  amrex::Print() << line_count << " data lines found in inject file" << std::endl;

  ProbParm::inject_N = line_count;
  ProbParm::inject_time->resize(ProbParm::inject_N);
  ProbParm::inject_mass->resize(ProbParm::inject_N);
  ProbParm::inject_vel->resize(ProbParm::inject_N);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  int pos1 = 0;
  int pos2 = 0;
  for (int i = 0; i < ProbParm::inject_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    amrex::Real test, test1, test2;
    sinput >> (*ProbParm::inject_time)[i];
    sinput >> (*ProbParm::inject_mass)[i];
    sinput >> (*ProbParm::inject_vel)[i];
  }
  ProbParm::d_inject_time = ProbParm::inject_time->dataPtr();
  ProbParm::d_inject_mass = ProbParm::inject_mass->dataPtr();
  ProbParm::d_inject_vel = ProbParm::inject_vel->dataPtr();
  ProbParm::jet_start_time = ProbParm::d_inject_time[0];
  ProbParm::jet_end_time = ProbParm::d_inject_time[line_count - 1];
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
  pp.query("inject_file", ProbParm::input_file);
  pp.query("init_v", ProbParm::v0);
  pp.query("ref_p", ProbParm::p0);
  pp.query("ref_T", ProbParm::T0);
  pp.query("init_N2", ProbParm::Y_N2);
  pp.query("init_O2", ProbParm::Y_O2);
  pp.query("jet_vel", ProbParm::jet_vel);
  pp.query("jet_dia", ProbParm::jet_dia);
  pp.query("part_mean_dia", ProbParm::part_mean_dia);
  pp.query("part_stdev_dia", ProbParm::part_stdev_dia);
  pp.query("part_rho", ProbParm::part_rho);
  pp.query("part_temp", ProbParm::part_temp);
  pp.query("mass_flow_rate", ProbParm::mass_flow_rate);
  pp.query("spray_angle_deg", ProbParm::spray_angle);
  // Convert to radians
  ProbParm::spray_angle *= M_PI/180.;

  // Initial density, velocity, and material properties
  amrex::Real eint, cs, cp;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[N2_ID] = ProbParm::Y_N2;
  massfrac[O2_ID] = ProbParm::Y_O2;
  EOS::PYT2RE(ProbParm::p0, massfrac, ProbParm::T0, ProbParm::rho0, eint);
  EOS::RTY2Cs(ProbParm::rho0, ProbParm::T0, massfrac, cs);
  EOS::TY2Cp(ProbParm::T0, massfrac, cp);
  if (ProbParm::input_file != "") {
    ProbParm::inject_time = new amrex::Gpu::ManagedVector<amrex::Real>;
    ProbParm::inject_mass = new amrex::Gpu::ManagedVector<amrex::Real>;
    ProbParm::inject_vel = new amrex::Gpu::ManagedVector<amrex::Real>;
    read_inject(ProbParm::input_file);
  }
}
}

void pc_prob_close()
{
  delete ProbParm::inject_time;
  delete ProbParm::inject_mass;
  delete ProbParm::inject_vel;

  ProbParm::inject_time = nullptr;
  ProbParm::inject_mass = nullptr;
  ProbParm::inject_vel = nullptr;
}
