#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

#define DISPLAY                      // if defined, a window will open and display the current configuration

#include "configuration.h"
#include "metropolis.h"

#define L 256                        // system length

int main(int argc, char *argv[]){
  if(argc < 3){
    std::cout << "Usage:\n\tOption 1: " << argv[0] << " basename temperature\n\tOption 2: " << argv[0] << " basename temperature_start temperature_end temperature_step\n\tOption 3: " << argv[0] << " basename temp1 temp2 temp3 temp4 ..." <<     std::endl;
    return 0;
  }
  std::string results_base_filename = argv[1];
  std::vector<float> temperature_list;
  if(argc == 5){
    for(uint8_t k = 0; k < (atof(argv[3])-atof(argv[2]))/(atof(argv[4]))+1; k++)
    {
      temperature_list.push_back(atof(argv[2])+atof(argv[4])*k);
    }
  }
  else{
    for(uint8_t k = 2; k < argc; k++)
    {
      temperature_list.push_back(atof(argv[k]));
    }
  }
  std::cout << "Will use the following temperatures: ";
  for(uint16_t i = 0; i < temperature_list.size(); i++){
    std::cout << temperature_list[i] << " ";
  }
  std::cout << std::endl;
  std::ofstream results_dist(results_base_filename+"_dist.dat",std::ofstream::out);
  std::ofstream results_stdev(results_base_filename+"_stdev.dat",std::ofstream::out);
  results_dist << "L\tT\tbias\tmag\tmag2\tmag4\te\te2\tx\tc\tU_L\n";
  results_stdev << "L\tT\tavg_mag\tstdev_mag\tavg_mag2\tstdev_mag2\tavg_mag4\tstdev_mag4\tavg_e\tstdev_e\tavg_e2\tstdev_e2\tavg_x\tstdev_x\tavg_c\tstdev_c\tavg_U_L\tstdev_U_L\n";
  // for each temperature in the list do...
  for(uint16_t i = 0; i < temperature_list.size(); i++){
    float T = temperature_list[i];
    float beta = 1./T;
    std::vector<double> mag_list;
    std::vector<double> mag2_list;
    std::vector<double> mag4_list;
    std::vector<double> e_list;
    std::vector<double> e2_list;
    std::vector<double> x_list;
    std::vector<double> c_list;
    std::vector<double> U_L_list;
    // for each temperature, use 10 different initial conditions with differnt bias
    for(uint8_t k = 1; k < 11; k++){
      float bias = std::exp(0.2*k);
      std::chrono::steady_clock::time_point begin;
      std::chrono::steady_clock::time_point end;
      metropolis<L> metrop(beta,bias);
      begin = std::chrono::steady_clock::now();
      //uint32_t frame_cycles = (L < 256)? 2*512/L*512/L : 10*L/256;
      uint32_t frame_cycles = (L < 256)? 2*512/L*512/L : 20*L/256;
      //uint32_t total_cycles = (L < 32)? 50000*128/L*128/L : 12500*512/L;
      uint32_t total_cycles = (L < 32)? 50000*128/L*128/L : 1000*512/L;
      //double magnetization = metrop.run(5000,total_cycles,1,frame_cycles);
      double magnetization = metrop.run(2000,total_cycles,1,frame_cycles);
      end = std::chrono::steady_clock::now();
      std::cout << "run() took " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << " seconds:" << std::endl;
      double susceptibility = (metrop.mean_magnetization_squared-metrop.mean_magnetization*metrop.mean_magnetization)/T*L*L;
      double heat_capacity = (metrop.mean_energy_squared-metrop.mean_energy*metrop.mean_energy)/(T*T)*L*L;
      double binder_cumulant = 1-metrop.mean_magnetization_fourth/(3.*metrop.mean_magnetization_squared*metrop.mean_magnetization_squared);
      std::cout << "L = " << L << ", T = " << T << ", bias = " << bias << ": m = " << magnetization << ", m^2 = " << metrop.mean_magnetization_squared << ", e = " << metrop.mean_energy << ", e^2 = " << metrop.mean_energy_squared << ", x = " << susceptibility << ", c = " << heat_capacity << ", U_L = " << binder_cumulant << std::endl;
      mag_list.push_back(magnetization);
      mag2_list.push_back(metrop.mean_magnetization_squared);
      mag4_list.push_back(metrop.mean_magnetization_fourth);
      e_list.push_back(metrop.mean_energy);
      e2_list.push_back(metrop.mean_energy_squared);
      x_list.push_back(susceptibility);
      c_list.push_back(heat_capacity);
      U_L_list.push_back(binder_cumulant);
      results_dist << L << "\t" << T << "\t" << bias << "\t" << magnetization << "\t" << metrop.mean_magnetization_squared << "\t" << metrop.mean_magnetization_fourth << "\t" << metrop.mean_energy << "\t" << metrop.mean_energy_squared << "\t" << susceptibility << "\t" << heat_capacity << "\t" << binder_cumulant << std::endl;
    }
    std::cout << "Summary: L = " << L << ", T = " << T << ": m = " << avg(mag_list) << " +- " << stdev(mag_list) << ", e = " << avg(e_list) << " +- " << stdev(e_list) << ", x = " << avg(x_list) << " +- " << stdev(x_list) << ", c = " << avg    (c_list) << " +- " << stdev(c_list) << ", U_L = " << avg(U_L_list) << " +- " << stdev(U_L_list) << std::endl;
    results_stdev << L << "\t" << T << "\t" << avg(mag_list) << "\t" << stdev(mag_list) << "\t" << avg(mag2_list) << "\t" << stdev(mag2_list) << "\t" << avg(mag4_list) << "\t" << stdev(mag4_list) << "\t" << avg(e_list) << "\t" << stdev(e_list) << "\t" << avg(e2_list) << "\t" << stdev(e2_list) << "\t" << avg(x_list) << "\t" << stdev(x_list) << "\t" << avg(c_list) << "\t" << stdev(c_list) << "\t" << avg(U_L_list) << "\t" << stdev(U_L_list) << std::endl;
  }
  results_dist.close();
  results_stdev.close();
}

