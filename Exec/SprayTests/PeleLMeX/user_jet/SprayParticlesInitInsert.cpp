
#include "SprayParticles.H"
#include "SprayInjection.H"
#include "pelelm_prob.H"

// This demonstrates a way to override the get_new_particle function in the
// SprayJet class
class ThisJet : public SprayJet
{
public:
  ThisJet(const std::string& jet_name, const amrex::Geometry& geom);

  ~ThisJet() override {}

  bool get_new_particle(
    const amrex::Real time,
    const amrex::Real& phi_radial,
    const amrex::Real& cur_radius,
    amrex::Real& umag,
    amrex::Real& theta_spread,
    amrex::Real& phi_swirl,
    amrex::Real& dia_part,
    amrex::Real& T_part,
    amrex::Real* Y_part) override;

protected:
  int smd_col = 0;
  int temp_col = 1;
  int vel_col = 2;
  int data_len;
  std::unique_ptr<DistBase> normDist;
  amrex::Vector<amrex::Real> jet_radius_vec;
  std::array<amrex::Vector<amrex::Real>, 3> mean_vals;
  std::array<amrex::Vector<amrex::Real>, 3> std_vals;
};

ThisJet::ThisJet(const std::string& jet_name, const amrex::Geometry& geom)
{
  std::string ppname = "spray." + jet_name;
  amrex::ParmParse pp(ppname);
  std::vector<amrex::Real> jcent(AMREX_SPACEDIM);
  pp.getarr("jet_cent", jcent);
  std::vector<amrex::Real> jnorm(AMREX_SPACEDIM);
  pp.getarr("jet_norm", jnorm);
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    m_cent[dir] = jcent[dir];
    m_norm[dir] = jnorm[dir];
  }
  amrex::Real mag = m_norm.vectorLength();
  m_norm /= mag;
  check_jet_cent(geom);
  pp.get("spread_angle", m_spreadAngle);
  m_spreadAngle *= M_PI / 180.; // Assumes spread angle is in degrees
  pp.query("start_time", m_startTime);
  pp.query("end_time", m_endTime);
  pp.get("mass_flow_rate", m_massFlow);
  std::string dist_type = "Normal";
  normDist = DistBase::create(dist_type);
  normDist->init(1., 1.);
  std::string infile;
  pp.get("injection_file", infile);
  std::ifstream iffile(infile);
  const std::string memfile =
    static_cast<std::stringstream const&>(std::stringstream() << iffile.rdbuf())
      .str();
  if (!amrex::FileSystem::Exists(infile)) {
    amrex::Abort("injection_file does not exist");
  }
  iffile.close();
  std::istringstream iss(memfile);
  std::string firstline, remaininglines;
  std::getline(iss, firstline);
  int line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  jet_radius_vec.resize(line_count);
  for (int col = 0; col < 3; ++col) {
    mean_vals[col].resize(line_count);
    std_vals[col].resize(line_count);
  }
  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  m_avgDia = 0.;
  m_jetT = 0.;
  m_jetVel = 0.;
  amrex::Real max_jet_dia = 0.;
  for (unsigned int i = 0; i < line_count; ++i) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> jet_radius_vec[i];
    jet_radius_vec[i] /= 1000.; // Locations are in mm
    max_jet_dia = amrex::max<amrex::Real>(max_jet_dia, 2. * jet_radius_vec[i]);
    sinput >> mean_vals[temp_col][i];
    sinput >> std_vals[temp_col][i];
    sinput >> mean_vals[smd_col][i];
    sinput >> std_vals[smd_col][i];
    mean_vals[smd_col][i] *= 1.E-6; // Diameters are in um
    std_vals[smd_col][i] *= 1.E-6;
    sinput >> mean_vals[vel_col][i];
    sinput >> std_vals[vel_col][i];
    m_avgDia += mean_vals[smd_col][i];
    m_jetT += mean_vals[temp_col][i];
    m_jetVel = amrex::max(m_jetVel, mean_vals[vel_col][i]);
  }
  // Must provide values for m_avgDia, m_jetT, m_jetVel, and m_jetY
  // The average values are used during the jet calculations
  m_avgDia = m_avgDia / static_cast<amrex::Real>(line_count);
  m_jetT = m_jetT / static_cast<amrex::Real>(line_count);
  if (SPRAY_FUEL_NUM > 1) {
    std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
    pp.getarr("Y", in_Y_jet);
    amrex::Real sumY = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      m_jetY[spf] = in_Y_jet[spf];
      sumY += in_Y_jet[spf];
    }
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      m_jetY[spf] /= sumY;
    }
  } else {
    m_jetY[0] = 1.;
  }
  data_len = line_count;
  if (pp.contains("jet_dia")) {
    pp.get("jet_dia", m_jetDia);
  } else {
    m_jetDia = max_jet_dia;
  }
}

bool
ThisJet::get_new_particle(
  const amrex::Real time,
  const amrex::Real& phi_radial,
  const amrex::Real& cur_radius,
  amrex::Real& umag,
  amrex::Real& theta_spread,
  amrex::Real& phi_swirl,
  amrex::Real& dia_part,
  amrex::Real& T_part,
  amrex::Real* Y_part)
{
  // Interpolate values from data tables
  amrex::Real Tmean, Tstd, SMDmean, SMDstd, Umean, Ustd, dxdx;
  int iloc = 0;
  if (cur_radius < jet_radius_vec[0]) {
    iloc = 0;
  }
  for (int i = 0; i < data_len - 1; ++i) {
    amrex::Real x1 = jet_radius_vec[i];
    amrex::Real x2 = jet_radius_vec[i + 1];
    if (cur_radius > x1 && cur_radius < x2) {
      dxdx = (cur_radius - x1) / (x2 - x1);
      iloc = i;
    }
  }
  for (int jc = 0; jc < 3; ++jc) {
    amrex::Real y1 = mean_vals[jc][iloc];
    amrex::Real y2 = mean_vals[jc][iloc + 1];
    amrex::Real yval = y1 + (y2 - y1) * dxdx;
    amrex::Real ys1 = std_vals[jc][iloc];
    amrex::Real ys2 = std_vals[jc][iloc + 1];
    amrex::Real ysval = ys1 + (ys2 - ys1) * dxdx;
    if (jc == temp_col) {
      Tmean = yval;
      Tstd = ysval;
    } else if (jc == smd_col) {
      SMDmean = yval;
      SMDstd = ysval;
    } else {
      Umean = yval;
      Ustd = ysval;
    }
  }
  // Change normalized distribution to match each variable
  T_part = Tmean + Tstd * (normDist->get_dia() - 1.);
  umag = Umean + Ustd * (normDist->get_dia() - 1.);
  dia_part = SMDmean + SMDstd * (normDist->get_dia() - 1.);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_part[spf] = m_jetY[spf];
  }
  phi_swirl = 0.;
  theta_spread = cur_radius / m_jetDia * m_spreadAngle;
  return true;
}

bool
SprayParticleContainer::injectParticles(
  amrex::Real time,
  amrex::Real dt,
  int nstep,
  int lev,
  int finest_level,
  ProbParm const& prob_parm)
{
  amrex::ignore_unused(nstep, finest_level);
  SprayJet* js = m_sprayJets[0].get();
  if (lev != 0) {
    return false;
  }
  if (!js->jet_active(time)) {
    return false;
  }

  sprayInjection(time, js, dt, lev);

  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(
  const bool init_parts, ProbParm const& prob_parm)
{
  amrex::ignore_unused(prob_parm);
  m_sprayJets.resize(1);
  std::string jet_name = "jet1";
  m_sprayJets[0] = std::make_unique<ThisJet>(jet_name, Geom(0));
  // Start without any particles
  m_injectVel = m_sprayJets[0]->jet_vel();
  return;
}
