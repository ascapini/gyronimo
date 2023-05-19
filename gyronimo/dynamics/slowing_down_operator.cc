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

// @slowing_down_operator.cc, this file is part of ::gyronimo::

#include <gyronimo/dynamics/slowing_down_operator.hh>
#include <gyronimo/core/codata.hh>
#include <random>
#include <vector>
#include <numbers>
#include <ctime>
#include <iostream>
#include <string>

namespace gyronimo {

slowing_down_operator::slowing_down_operator(
        double V_ref,
	double tp_mass, double tp_temperature, double tp_charge,
	double freq)
        //double bckg_mass, double bckg_temperature, double bckg_charge, double bckg_density,
        //double coulomb_log)//, double dt)
	: generator_r_(std::time(0)), gaussian_(0.0, 1.0),
	  V_ref_(V_ref),
	  tp_mass_(tp_mass), tp_temperature_(tp_temperature * 1.60218e-19), tp_charge_(tp_charge),
	  const_freq_(freq ? freq : 0.0){
	   	
}

double slowing_down_operator::operator()(
	double v, const double dt) {
  std::size_t num_species = background_.size();
  double square_in_velocity = v*v; //normalized
  for(int i = 0; i < num_species; i++) {
    square_in_velocity += this->delta_v2(v, dt, background_[i])/std::pow(V_ref_,2);
  };
  return std::sqrt(v*v); 
}



//pitch_scattering_operator::coll_state pitch_scattering_operator::generate_state(double vpp, double v) const {
//  return {vpp, v};
//}

double slowing_down_operator::delta_v2(double v, const double dt, slowing_down_operator::background_pop specie) {
  if(const_freq_ == 0.0){
    double vel = v * V_ref_;
    double x_bckg = std::pow(vel,2)/(2*specie.temperature/specie.mass);
    double r = gaussian_(generator_r_);
    double freq = this->frequency(vel, specie);
    return -2 * freq * dt / (1 + specie.mass/tp_mass_) * (std::pow(vel,2) - 1 / (std::erf(std::sqrt(x_bckg)) - 2 * std::sqrt(x_bckg/std::numbers::pi) * std::exp(-x_bckg)) * 2 * std::pow(x_bckg,1.5) / std::sqrt(std::numbers::pi) * specie.mass / tp_mass_ * (2*specie.temperature/specie.mass) * std::exp(-x_bckg)) + r * 2 * std::sqrt((2*specie.temperature/specie.mass) * std::pow(vel,2) * freq * dt / (1 + tp_mass_/specie.mass));
  }
  else{
    double vel = v * V_ref_;
    double x_bckg = std::pow(vel,2)/(2*specie.temperature/specie.mass);
    double r = gaussian_(generator_r_);
    double freq = const_freq_;
    return -2 * freq * dt / (1 + specie.mass/tp_mass_) * (std::pow(vel,2) - 1 / (std::erf(std::sqrt(x_bckg)) - 2 * std::sqrt(x_bckg/std::numbers::pi) * std::exp(-x_bckg)) * 2 * std::pow(x_bckg,1.5) / std::sqrt(std::numbers::pi) * specie.mass / tp_mass_ * (2*specie.temperature/specie.mass) * std::exp(-x_bckg)) + r * 2 * std::sqrt((2*specie.temperature/specie.mass) * std::pow(vel,2) * freq * dt / (1 + tp_mass_/specie.mass));
  };
}
    

double slowing_down_operator::frequency(const double v, slowing_down_operator::background_pop specie) const {  //takes as input particle velocity in SI units, non normalized
  double x_bckg = std::pow(v,2) / (2*specie.temperature/specie.mass);
  double freq = specie.frequency_coeff * std::pow(2*tp_temperature_ /tp_mass_,1.5) / std::pow(v,3);
  //double x_tp = v/std::sqrt(2*tp_temperature_/tp_mass_;
  return freq * (1 + tp_mass_/specie.mass) * (std::erf(std::sqrt(x_bckg)) - 2 * std::sqrt(x_bckg/std::numbers::pi) * std::exp(-x_bckg));
}

void slowing_down_operator::add_pop(double bckg_mass, double bckg_charge, double bckg_temperature, double bckg_density, double coulomb_log) {
  slowing_down_operator::background_pop pop;
  pop.mass = bckg_mass;
  pop.charge = bckg_charge;
  pop.temperature = bckg_temperature * 1.60218e-19;
  pop.density = bckg_density;
  pop.coulomb_log = coulomb_log;
  pop.frequency_coeff = pop.density * std::pow(tp_charge_*pop.charge,2) * pop.coulomb_log * std::pow(codata::mu0,2) * std::pow(codata::c,4) / (4 * std::numbers::pi *  std::pow(tp_mass_,2) * std::pow(2*tp_temperature_/tp_mass_,1.5));
  background_.push_back(pop);
}

}
