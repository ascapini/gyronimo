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

// @slowing_down_operator.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_SLOWING_DOWN_OPERATOR
#define GYRONIMO_SLOWING_DOWN_OPERATOR

#include <random> //for rng and distributions
#include <utility> //for pairs
#include <vector>
#include <cmath>

/*

***********************************************************************
*******************************   TODO   ******************************
***********************************************************************

- How to deal with giga standard deviation
- Substitute incomplete gamma function in the algorithm
- Add time normalization when we'll express dt using a reference time 
  or frequency. In the calculations dt_ has to be expressed in s
- Add dependencies on plasma parameters
- Organize particle properties in struct?
*/


namespace gyronimo {

//
/*

*/
class slowing_down_operator {
  public:
    
    typedef struct background_pop {
                       double mass;
                       double charge;
                       double temperature;
                       double density;
                       double coulomb_log;
                       double frequency_coeff;
                       } background_pop;
                       
    //typedef std::pair<double, double> coll_state; // pair containing vpp, v of a particle
    
    slowing_down_operator(
        double V_ref,
    	double tp_mass, double tp_temperature, double tp_charge,
    	double freq = 0.0);
        ////double bckg_mass, double bckg_temperature, double bckg_charge, double bckg_density,
       // double coulomb_log);//, double dt);
    ~slowing_down_operator() {};
    
    double operator()(double v, const double dt);
  
 //   double generate_state(double vpp, double v) const;
    
    double delta_v2(double v, const double dt, background_pop specie);
    
    double frequency(const double v, background_pop specie) const;
    
    void add_pop(double bckg_mass, double bckg_charge, double bckg_temperature, double bckg_density, double coulomb_log);
    
  private:
    std::mt19937 generator_r_;
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
//    const double dt_;
    //double frequency_coeff_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_SLOWING_DOWN_OPERATOR. 
