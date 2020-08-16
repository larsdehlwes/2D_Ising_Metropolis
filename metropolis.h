#ifndef METROPOLIS_H
#define METROPOLIS_H

#ifdef DEBUG 
#define D(x) (x)
#else
#define D(x) do{}while(0)
#endif

#include <fmt/core.h>
#include <string>
#include <iostream>
#include <math.h>
#include <chrono>
#include "configuration.h"
#include "avg_stdev.h"
#include <vector>

template <uint16_t ARRAY_LEN>
class metropolis: public configuration<ARRAY_LEN>
{
  public:
    metropolis(float _beta, float bias = 1);                                                                                 // constructor 
    metropolis(std::string _filename, float _beta, float bias = 1);                                                          // constructor with filename argument
    uint16_t* wiggle_random_spin();                                                                                          // choose a random spin and flip it if the condition is met
    int8_t energy_change_upon_flip(uint16_t i, uint16_t j);                                                                  // return the energy change upon flipping the spin at (i,j)
    void draw_information();                                                                                                 // display the number of cycles and magnetization
    void datawrite();                                                                                                        // append the current magnetization and energy to the datafile
    double run(uint32_t mincycles = 4000, uint32_t cycles = 10000, uint32_t eval_cycles = 1, uint32_t frame_cycles = 1);     // runs the Monte-Carlo simulation
    double mean_magnetization;                                                                                               // average abolute value of the magnetization per spin
    double mean_magnetization_squared;                                                                                       // average square of the magnetization per spin
    double mean_magnetization_fourth;                                                                                        // average fourth power of the magnetization per spin
    double mean_energy;                                                                                                      // average energy per spin
    double mean_energy_squared;                                                                                              // the square of the energy per spin
  private:
    float beta;                                                                                                              // beta (-> temperature)
    int64_t iter;                                                                                                            // iterations carried out
};

template <uint16_t ARRAY_LEN>
metropolis<ARRAY_LEN>::metropolis(float _beta, float bias) : configuration<ARRAY_LEN>(fmt::format("results/beta={:.4f}_N={:d}_bias={:.2f}",_beta,ARRAY_LEN,bias),bias)
{
  beta = _beta;
  iter = 0;
}

template <uint16_t ARRAY_LEN>
metropolis<ARRAY_LEN>::metropolis(std::string _filename, float _beta, float bias) : configuration<ARRAY_LEN>(_filename,bias)
{
  beta = _beta;
  iter = 0;
}

template <uint16_t ARRAY_LEN> void
metropolis<ARRAY_LEN>::draw_information()
{
  if(ARRAY_LEN >= 200){
    cv::putText(this->bgr, "iter = ", cv::Point(ARRAY_LEN/100,ARRAY_LEN+ARRAY_LEN/17), cv::FONT_HERSHEY_DUPLEX, (float) ARRAY_LEN/500., cv::Scalar(255,0,0), ARRAY_LEN/200);
    cv::putText(this->bgr, fmt::format("{:.0f}", (float) iter/(((uint32_t) ARRAY_LEN) * ((uint32_t) ARRAY_LEN))), cv::Point(ARRAY_LEN/4.54,ARRAY_LEN+ARRAY_LEN/17), cv::FONT_HERSHEY_DUPLEX, (float) ARRAY_LEN/500., cv::Scalar(255,0,0), ARRAY_LEN/200);
    cv::putText(this->bgr, "m = ", cv::Point(ARRAY_LEN/2+5,ARRAY_LEN+ARRAY_LEN/17), cv::FONT_HERSHEY_DUPLEX, (float) ARRAY_LEN/500., cv::Scalar(255,0,0), ARRAY_LEN/200);
    cv::putText(this->bgr, fmt::format("{:.4f}", this->get_magnetization()), cv::Point(ARRAY_LEN/2+ARRAY_LEN/5.5,ARRAY_LEN+ARRAY_LEN/17), cv::FONT_HERSHEY_DUPLEX, (float) ARRAY_LEN/500., cv::Scalar(255,0,0), ARRAY_LEN/200);
  }
}

template <uint16_t ARRAY_LEN> int8_t
metropolis<ARRAY_LEN>::energy_change_upon_flip(uint16_t i, uint16_t j){
  return 2 * (-!this->get_spin(i,j) + this->get_spin(i,j)) * (-4. + 2. * (this->get_spin(this->idx(i-1),j) + this->get_spin(this->idx(i+1),j) + this->get_spin(i,this->idx(j-1)) + this->get_spin(i,this->idx(j+1)) ) );
}

template <uint16_t ARRAY_LEN> void
metropolis<ARRAY_LEN>::datawrite()
{
  this->datafile << fmt::format("{:.2f}",(float) this->iter/(((uint32_t) ARRAY_LEN) * ((uint32_t) ARRAY_LEN))) << "\t" << fmt::format("{:.6f}",this->get_magnetization()) <<  "\t" << fmt::format("{:.6f}",this->get_energy()) << std::endl;
}

template < uint16_t ARRAY_LEN > uint16_t*
metropolis<ARRAY_LEN>::wiggle_random_spin(){
  std::chrono::steady_clock::time_point begin;
  std::chrono::steady_clock::time_point end;
  D(begin = std::chrono::steady_clock::now());
  uint16_t i = this->int_distribution(this->rng);
  uint16_t j = this->int_distribution(this->rng);
  int8_t energy_change = energy_change_upon_flip(i,j);
  if(energy_change <= 0){
    this->invert_spin(i,j);
  } 
  else{
    double rnd = this->real_distribution(this->rng);
    if(rnd < exp(-beta*energy_change)) this->invert_spin(i,j);
  } 
  uint16_t* out = new uint16_t[2];
  out[0] = i;
  out[1] = j;
  iter++;
  D(end = std::chrono::steady_clock::now());
  D(std::cout << "Wiggle took " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "ns." << std::endl);
  return out;
}

template <uint16_t ARRAY_LEN> double 
metropolis<ARRAY_LEN>::run(uint32_t mincycles, uint32_t cycles, uint32_t eval_cycles, uint32_t frame_cycles){
  uint32_t k = 0;
  uint32_t cycle = 0;
  int32_t counter = 0;
  this->gray2bgr();
  this->draw_information();
  this->imshow();
  this->vidwrite();
  this->datawrite();
  bool start_averaging = false;
  uint32_t initial_cycle = 0;
  uint16_t averaging_over = (eval_cycles > 1)? 1000/eval_cycles : 1000;
  std::vector<double> indices;
  std::vector<double> last_magnetization_values;
  char key = 0;
  uint32_t last_frame = 0;
  while(key != 27 && cycle*start_averaging < initial_cycle + cycles)
  {
    for(uint32_t i = 0; i < eval_cycles*((uint32_t) ARRAY_LEN)*((uint32_t) ARRAY_LEN); i++){
      uint16_t* ptr = this->wiggle_random_spin();
      delete[] ptr;
    }
    cycle = iter/(((uint32_t) ARRAY_LEN) * ((uint32_t) ARRAY_LEN));
    double magnetization = this->get_magnetization();
    double energy = this->get_energy();
    last_magnetization_values.push_back(magnetization);
    indices.push_back(((double) iter)/(((uint32_t) ARRAY_LEN) * ((uint32_t) ARRAY_LEN)));
    if(k > averaging_over)
    {
      last_magnetization_values.erase(last_magnetization_values.begin());
      indices.erase(indices.begin());
      if(!start_averaging && std::abs(slope(indices,last_magnetization_values)) < 0.000001 && cycle > mincycles)
      {
        std::cout << "Target slope " << std::abs(slope(indices,last_magnetization_values)) << " reached at " << cycle << "." << std::endl;
        start_averaging = true;
        initial_cycle = cycle;
      }
    }
    if(start_averaging){
      mean_magnetization = (counter == 0) ? std::abs(magnetization) : (mean_magnetization*counter + std::abs(magnetization))/(counter+1);
      mean_magnetization_squared = (counter == 0) ? magnetization*magnetization : (mean_magnetization_squared*counter + magnetization*magnetization)/(counter+1);
      mean_magnetization_fourth = (counter == 0) ? magnetization*magnetization*magnetization*magnetization : (mean_magnetization_fourth*counter + magnetization*magnetization*magnetization*magnetization)/(counter+1);
      mean_energy = (counter == 0) ? energy : (mean_energy*counter + energy)/(counter+1);
      mean_energy_squared = (counter == 0) ? energy*energy : (mean_energy_squared*counter + energy*energy)/(counter+1);
      counter++;
    }
    if(cycle >= last_frame + frame_cycles){
      this->gray2bgr();
      this->draw_information();
      this->imshow();
      this->datawrite();
      this->vidwrite();
      last_frame = cycle;
      key = cv::waitKey(1);
    }
    k++;
  }
  this->destroyWindow();
  this->vidrelease();
  this->datafile.close();
  return mean_magnetization;
}
#endif
