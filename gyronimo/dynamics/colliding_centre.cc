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

// @colliding_centre.cc, this file is part of ::gyronimo::

/* *************************************************************
   *************                                   *************
   *************          IDEAS / TODOs            *************     
   *************                                   *************
   *************************************************************
   
   COLLISIONS
   - Add slowing down due to ions for slow particles
   - Add pitch scattering with other species
   
   
   - Follow same normalization as guiding centre!
   - !!!! before going on make code cleaner, it is a real mess
   - generalize observer, so that one can choose its own when 
     using the code
   - initializer from custom distribution function OK
   - build distribution function from particles
   - momentum cons
   - express everything in terms of energy, mu and pitch? 
   - time normalization????
   - time difference between step advance and step collide
   - dt inside or outside of cc initialization?
   - understand units for freq calculation
   - add seed for rng
   - advance t->t+dt after advancement, before collisions (?)
   - remove dt_ or choose something different?
   - better way to get vpp sign from colliding centre state
   - add stepper into the class so that it is created only once
     per particle (or even outside of the class)
   - definition of an operator instead of do_step?
   - check normalization of mu_tilde and other quantities
   - delete unnecessary included packages
  (- from kauffmann 2012, substitute rnd number from gaussian
     with the other option)exponen
*/



/* *************************************************************
   *************                                   *************
   *************               DONE                *************     
   *************                                   *************
   *************************************************************
   
   - time step vs char time of polynomial
*/

  
#include <iostream>
#include <cmath>
#include <gyronimo/core/error.hh>
#include <gyronimo/core/codata.hh>
#include <gyronimo/dynamics/colliding_centre.hh>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <random>

namespace gyronimo {

//! Class constructor.
/*!
   
*/
colliding_centre::colliding_centre(
    double Lref, double Vref, double qom, double dt,
    pitch_scattering_operator ps_operator, slowing_down_operator sd_operator,
    const IR3field_c1* B, const IR3field* E)//t_in, dt)
    : //generator_r(3), generator_phi(42),
      //uniform_(0.0, 2.0*std::numbers::pi), gaussian_(0.0, 1.0),
      ps_operator_(ps_operator), sd_operator_(sd_operator),
      qom_tilde_(qom),
      electric_field_(E), magnetic_field_(B),
      Lref_(Lref), Vref_(Vref), Tref_(Lref/Vref),
      Bfield_time_factor_(Lref/(Vref*B->t_factor())),
      Efield_time_factor_(E ? Lref/(Vref*E->t_factor()) : 0.0),
      dt_(dt)
      {
  Oref_tilde_ = qom*codata::e/codata::m_proton*B->m_factor()*Tref_;
  iOref_tilde_ = 1.0/Oref_tilde_;
  if (E ?
      std::fabs(1.0 - E->m_factor()/(Vref_*B->m_factor())) > 1.0e-14 : false)
          error(__func__, __FILE__, __LINE__, "inconsistent E field.", 1);   
}
   
void colliding_centre::step_advance(coll_state* cs, double time) {
  guiding_centre gc(Lref_, Vref_, qom_tilde_, cs->at(4), magnetic_field_, electric_field_);
  IR3 position = this->get_position(*cs);
  double Btime = time*Bfield_time_factor_;
  double en_tilde = cs->at(3)*cs->at(3) + cs->at(4)*magnetic_field_->magnitude(position, Btime);
  guiding_centre::vpp_sign sign;
  if (this->get_vpp(*cs) >= 0) {
    sign = guiding_centre::vpp_sign::plus;}
  else { 
    sign = guiding_centre::vpp_sign::minus;};
  guiding_centre::state instate = gc.generate_state(position, en_tilde, sign, time);
  odeint_adapter adapter(&gc);
  boost::numeric::odeint::runge_kutta4<guiding_centre::state> rk;
  rk.do_step(adapter, instate, time, dt_);
  cs->at(0)=instate[0]; cs->at(1)=instate[1]; cs->at(2)=instate[2]; cs->at(3)=instate[3];
}
  
/*void colliding_centre::step_advance(coll_state* cs, double time) { //w/ time normalized
  double tau = time/Tref_;
  double dtau = dt_/Tref_;
  guiding_centre gc(Lref_, Vref_, qom_tilde_, cs->at(4), magnetic_field_, electric_field_);
  IR3 position = this->get_position(*cs);
  double Btime = tau*Bfield_time_factor_;
  double en_tilde = cs->at(3)*cs->at(3) + cs->at(4)*magnetic_field_->magnitude(position, Btime);
  guiding_centre::vpp_sign sign;
  if (this->get_vpp(*cs) >= 0) {
    sign = guiding_centre::vpp_sign::plus;}
  else { 
    sign = guiding_centre::vpp_sign::minus;};
  guiding_centre::state instate = gc.generate_state(position, en_tilde, sign, tau);
  odeint_adapter adapter(&gc);
  boost::numeric::odeint::runge_kutta4<guiding_centre::state> rk;
  rk.do_step(adapter, instate, tau, dtau);
  cs->at(0)=instate[0]; cs->at(1)=instate[1]; cs->at(2)=instate[2]; cs->at(3)=instate[3];
}*/


void colliding_centre::step_pitch_scattering(coll_state* cs, double time) { 
  
  //advancing pitch 
  double vpp = cs->at(3);
  double v = std::sqrt(this->energy_parallel(*cs) + this->energy_perpendicular(*cs, time));
  pitch_scattering_operator::pitch_state in_state = ps_operator_.generate_state(vpp, v);
  pitch_scattering_operator::pitch_state fin_state = ps_operator_(in_state, dt_ * Tref_);
  //sub vpp
  cs->at(3) = fin_state.first;
  //sub mu
  double v_perp_out = std::sqrt(v*v - std::pow(fin_state.first,2));
  cs->at(4) = std::pow(v_perp_out,2)/magnetic_field_->magnitude(this->get_position(*cs), time*Bfield_time_factor_);
}
  
void colliding_centre::step_slowing_down(coll_state* cs, double time) { //advancing energy
  double initial_v = std::sqrt(this->energy_parallel(*cs) + this->energy_perpendicular(*cs, time));  
  double new_v = sd_operator_(initial_v, dt_ * Tref_);
  double ratio = new_v/initial_v;
  cs->at(3) = cs->at(3) * ratio;
  //double vperp_out = std::sqrt(new_v*new_v - std::pow(cs->at(3),2));
  //cs->at(4) = std::pow(vperp_out,2)/magnetic_field_->magnitude(this->get_position(*cs), time*Bfield_time_factor_);//cs->at(4) * std::pow(ratio,2);
  cs->at(4) = cs->at(4) * std::pow(ratio,2);
  /*double vel = std::sqrt(this->energy_parallel(*cs) + this->energy_perpendicular(*cs, time));
  double pitch_in = this->get_vpp(*cs)/vel;
  double r = gaussian_(generator_r);
  double dteta = r * std::sqrt(2 * freq * dt_);
  double phi = uniform_(generator_phi);
  double pitch_out = std::sin(dteta)*std::sin(phi)*std::sqrt(1-pitch_in*pitch_in) + pitch_in*std::cos(dteta);
  double vpar_out = vel * pitch_out;
  double vperp_out = std::sqrt(vel*vel - vpar_out*vpar_out);
  cs->at(3)=vpar_out; cs->at(4)=vperp_out*vperp_out/magnetic_field_->magnitude(this->get_position(*cs), time*Bfield_time_factor_);*/ 
}

void colliding_centre::observer(coll_state cs, double time) const { //const coll_state& cs)
  IR3 x = this->get_position(cs);
  std::cout << time * Tref_ << " " 
  << x[IR3::u] << " " << x[IR3::v] << " "
  << this->get_vpp(cs) << " " 
  << this->get_pitch(cs, time) << " "
  << this->energy_parallel(cs) << " " 
  << this->energy_perpendicular(cs, time) << " "
  << this->get_mu(cs) << "\n";
} 
  
double colliding_centre::get_vpp(const coll_state& s) const {
  return s[3];
}

IR3 colliding_centre::get_position(const coll_state& s) const {
  return {Lref_*s[0], Lref_*s[1], Lref_*s[2]};
}

double colliding_centre::get_mu(const coll_state& s) const {
  return s[4];
}

double colliding_centre::get_pitch(const coll_state& s, double time) const{
  double v = std::sqrt(this->energy_parallel(s) + this->energy_perpendicular(s, time));
  double vpp = this->get_vpp(s);
  return vpp/v;
}

//double colliding_centre::coll_freq(const double v) const {
  

colliding_centre::coll_state colliding_centre::generate_state(
    const IR3& position, double energy_tilde,
    vpp_sign sign, double mu_tilde_, double time) const {
  double iLref = 1.0/Lref_;
  double Btime = time*Bfield_time_factor_;
  double B = magnetic_field_->magnitude(position, Btime);
  double vpp = (double)(sign)*std::sqrt(energy_tilde - mu_tilde_*B);
  //std::cout << vpp <<std::endl;
  return {iLref*position[0], iLref*position[1], iLref*position[2], vpp, mu_tilde_};
}

double colliding_centre::energy_parallel(const coll_state &s) const {
  return s[3]*s[3];
}

double colliding_centre::energy_perpendicular(const coll_state &s, double time) const {
  double Btime = time*Bfield_time_factor_;
  double B = magnetic_field_->magnitude(this->get_position(s), Btime);
  return colliding_centre::get_mu(s)*B;
}

void colliding_centre::do_step(coll_state* cs, double time) {
  this->step_advance(cs, time);
  this->step_pitch_scattering(cs, time + dt_);
  this->step_slowing_down(cs, time + dt_);
  this->observer(*cs, time + dt_);
}


//initializer from distribution function
//void colliding_centre::initializer(void function)

} // end namespace gyronimo
