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

// @pitch_scattering_operator.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_PITCH_SCATTERING_OPERATOR
#define GYRONIMO_PITCH_SCATTERING_OPERATOR

#include <random> //for rng and distributions
#include <utility> //for pairs
#include <cmath>
#include <vector>

/*

***********************************************************************
*******************************   TODO   ******************************
***********************************************************************
- Add time normalization when we'll express dt using a reference time 
  or frequency. In the calculations dt_ has to be expressed in s
- Add dependencies on plasma parameters

*/


namespace gyronimo {

//
/*

*/
class pitch_scattering_operator {
  public:
    typedef std::pair<double, double> pitch_state; // pair containing vpp, v of a particle
    typedef struct background_pop {
                      double mass;
                      double charge;
                      double temperature;
                      double density;
                      double coulomb_log;
                      double frequency_coeff;
                      } background_pop;
    //typedef struct {double vpp; double vel;} pitch_state;
    
    pitch_scattering_operator(
    	double V_ref,
    	double tp_mass, double tp_temperature, double tp_charge,
    	double freq = 0.0);
        //double bckg_mass, double bckg_temperature, double bckg_charge, double bckg_density,
        //double coulomb_log);//, double dt);
    ~pitch_scattering_operator() {};
    
    pitch_state operator()(pitch_state ps, const double dt);
  
    pitch_state generate_state(double vpp, double v) const;
    
    double frequency(const double v, background_pop specie) const;
    
    void add_pop(double bckg_mass, double bckg_charge, double bckg_temperature, double bckg_density, double coulomb_log);

  private:
    std::mt19937 generator_r_, generator_phi_;
    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> gaussian_;
    std::vector<background_pop> background_;
    const double V_ref_;
    const double tp_mass_; //all in SI units except for T in ev
    const double tp_temperature_;
    const double tp_charge_;
    const double const_freq_;
    //const double bckg_mass_;
    //const double bckg_temperature_;
    //const double bckg_charge_;
    //const double bckg_density_; 
    //const double coulomb_log_;
  //  const double dt_;
    //double frequency_const_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_PITCH_SCATTERING_OPERATOR. 
