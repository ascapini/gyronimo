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

// @colliding_centre.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_COLLIDING_CENTRE
#define GYRONIMO_COLLIDING_CENTRE

#include <array>
#include <numbers>
#include <gyronimo/dynamics/slowing_down_operator.hh>
#include <gyronimo/dynamics/pitch_scattering_operator.hh>
#include <gyronimo/dynamics/guiding_centre.hh>
#include <gyronimo/dynamics/odeint_adapter.hh>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

#include <random> // for std::mt19937

namespace gyronimo {

//! Advancement of particle including collision effect.
/*!
    Takes colliding centre state as an input, in particular the initial guiding 
    centre and its equations. Then creates a guiding centre object and advances 
    it one time-step using te boost library.                                       could insert integration algorithm an input, toghether with timestep,...
    After this first step, computes new values of v_par and magnetic moment.       how do I mantain a generic type for the integration algorithm and method, so that we 
*/
class colliding_centre {
 public:
  typedef std::array<double, 5> coll_state; //array containing position, v_par and mu, or their derivatives
  typedef boost::numeric::odeint::runge_kutta4<guiding_centre::state> runge;
  enum vpp_sign : int {minus = -1, plus = 1};
  
 // boost::numeric::odeint::runge_kutta4<guiding_centre::state> rk_;
    
  colliding_centre(
      double Lref, double Vref, double qom, double dt, pitch_scattering_operator ps_operator, slowing_down_operator sd_operator, const IR3field_c1* B, const IR3field* E = nullptr);
  ~colliding_centre() {};
    
    
  double get_vpp(const coll_state& s) const;
  IR3 get_position(const coll_state& s) const;
  double get_mu(const coll_state& s) const;
  double energy_parallel(const coll_state& s) const;
  double energy_perpendicular(const coll_state& s, double time) const;
  guiding_centre::vpp_sign get_sign(const coll_state cs) const;
  coll_state generate_state(
      const IR3& position,
      double energy_tilde, vpp_sign sign, double mu_tilde, double time = 0) const;
  double get_pitch(const coll_state& s, double time) const;   
//  double coll_freq(const double v) const; //calculate collision frequency from particle velocity
        
// Advances state using gc equationn with runge-kutta4; updates position and vpp in the coll_state.   
  void step_advance(coll_state* cs, double time); //const

  void step_pitch_scattering(coll_state* cs, double time); //const
  void step_slowing_down(coll_state* cs, double time);
  void observer(coll_state cs, double time) const;
  
  void do_step(coll_state* cs, double time); //const
    
  double Lref() const {return Lref_;}; 
  double Tref() const {return Tref_;}; 
  double Vref() const {return Vref_;}; 
  double dt() const {return dt_;};
  double qom_tilde() const {return qom_tilde_;};
  double Oref_tilde() const {return Oref_tilde_;}; 
  const IR3field* electric_field() const {return electric_field_;};
  const IR3field_c1* magnetic_field() const {return magnetic_field_;}; 
    
 private:
 // std::mt19937 generator_r, generator_phi;
 // std::uniform_real_distribution<double> uniform_;
 // std::normal_distribution<double> gaussian_;
  pitch_scattering_operator ps_operator_;
  slowing_down_operator sd_operator_;
  const IR3field* electric_field_;
  const IR3field_c1* magnetic_field_;
  const double qom_tilde_;
  const double Lref_, Vref_, Tref_;
  const double Bfield_time_factor_, Efield_time_factor_;
  const double dt_;
  //const boost::numeric::odeint::runge_kutta4<guiding_centre::state> rk_;
  double Oref_tilde_, iOref_tilde_;
  
  //
  
};

} // end namespace gyronimo.

#endif // GYRONIMO_COLLIDING_CENTRE.
        
        
        
        
         
