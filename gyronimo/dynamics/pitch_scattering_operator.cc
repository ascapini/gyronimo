// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// ::gyronimo:: is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ::gyronimo:: is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ::gyronimo::.  If not, see <https://www.gnu.org/licenses/>.

// @pitch_scattering_operator.cc, this file is part of ::gyronimo::

#include <gyronimo/dynamics/pitch_scattering_operator.hh>
#include <gyronimo/core/codata.hh>
#include <random>
#include <numbers>
#include <ctime>
#include <vector>


namespace gyronimo {

pitch_scattering_operator::pitch_scattering_operator(
	double V_ref,
	double tp_mass, double tp_temperature, double tp_charge,
	double freq)
        //double bckg_mass, double bckg_temperature, double bckg_charge, double bckg_density,
        //double coulomb_log)//, double dt)
	: generator_r_(std::time(0)), generator_phi_(std::time(0) + 42),
	  uniform_(0.0, 2.0*std::numbers::pi), gaussian_(0.0, 1.0),
	  V_ref_(V_ref),
	  tp_mass_(tp_mass), tp_temperature_(tp_temperature * 1.60218e-19), tp_charge_(tp_charge),
	  const_freq_(freq ? freq : 0.0){ 
	  
}

/*pitch_scattering_operator::pitch_state pitch_scattering_operator::operator()(
	pitch_state ps, const double dt) {
  double pitch_in = ps.vpp/ps.vel;
  double r = gaussian_(generator_r_);
  double dteta = r * std::sqrt(2 * this->frequency(ps.vel * V_ref_) * dt); //dt in constructor, freq in ?
  double phi = uniform_(generator_phi_);
  double pitch_out = std::sin(dteta)*std::sin(phi)*std::sqrt(1-pitch_in*pitch_in) + pitch_in*std::cos(dteta);
  double vpar_out = ps.vel * pitch_out;
  pitch_scattering_operator::pitch_state final_state;
  final_state.vpp = vpar_out;
  final_state.vel = ps.vel;
  return final_state;
}

pitch_scattering_operator::pitch_state pitch_scattering_operator::generate_state(double vpp, double v) const {
  pitch_scattering_operator::pitch_state state;
  state.vpp = vpp;
  state.vel = v;
  return state;
}*/
pitch_scattering_operator::pitch_state pitch_scattering_operator::operator()(
	pitch_state ps, const double dt) {
  if(const_freq_ == 0.0){	
    std::size_t num_species = background_.size();
    double tot_freq = 0;
    for(int i = 0; i < num_species; i++) {
      tot_freq += this->frequency(ps.second * V_ref_, background_[i]);
      };
    double pitch_in = ps.first/ps.second;
    double r = gaussian_(generator_r_);
    double dteta = r * std::sqrt(2 * tot_freq * dt); //dt in constructor, freq in ?
    double phi = uniform_(generator_phi_);
    double pitch_out = std::sin(dteta)*std::sin(phi)*std::sqrt(1-pitch_in*pitch_in) + pitch_in*std::cos(dteta);
    double vpar_out = ps.second * pitch_out;
    return {vpar_out, ps.second};
  }
  else {
  double pitch_in = ps.first/ps.second;
    double r = gaussian_(generator_r_);
    double dteta = r * std::sqrt(2 * const_freq_ * dt); //dt in constructor, freq in ?
    double phi = uniform_(generator_phi_);
    double pitch_out = std::sin(dteta)*std::sin(phi)*std::sqrt(1-pitch_in*pitch_in) + pitch_in*std::cos(dteta);
    double vpar_out = ps.second * pitch_out;
    return {vpar_out, ps.second};
  }
}

pitch_scattering_operator::pitch_state pitch_scattering_operator::generate_state(double vpp, double v) const {
  return {vpp, v};
}

double pitch_scattering_operator::frequency(const double v, pitch_scattering_operator::background_pop specie) const {  //takes as input particle velocity in SI units, non normalized
  double x_bckg = v/std::sqrt(2*specie.temperature/specie.mass);
  double x_tp = v/std::sqrt(2*tp_temperature_/tp_mass_);
  return specie.frequency_coeff * (std::erf(x_bckg) * (1 - 1 / (2 * std::pow(x_bckg,2))) + std::exp(-std::pow(x_bckg,2))/(x_bckg * std::sqrt(std::numbers::pi))) / std::pow(x_tp,3);
}

void pitch_scattering_operator::add_pop(double bckg_mass, double bckg_charge, double bckg_temperature, double bckg_density, double coulomb_log) {
  pitch_scattering_operator::background_pop pop;
  pop.mass = bckg_mass;
  pop.charge = bckg_charge;
  pop.temperature = bckg_temperature * 1.60218e-19;
  pop.density = bckg_density;
  pop.coulomb_log = coulomb_log;
  pop.frequency_coeff = pop.density * std::pow(tp_charge_*pop.charge,2) * pop.coulomb_log * std::pow(codata::mu0,2) * std::pow(codata::c,4) / (4 * std::numbers::pi *  std::pow(tp_mass_,2) * std::pow(2*tp_temperature_/tp_mass_,1.5));
  background_.push_back(pop);
}



}
