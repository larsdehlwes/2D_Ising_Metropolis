#ifndef CONF_H
#define CONF_H

#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do{}while(0)
#endif

#ifdef DISPLAY
#define Display(x) (x)
#else
#define Display(x) do{}while(0)
#endif

#include <string>
#include <fmt/core.h>
#include <fstream>
#include <random>
#include <opencv2/opencv.hpp>

template <uint16_t ARRAY_LEN>
class configuration
{
  public:
    configuration(std::string _filename, float bias);           // constructor
    bool get_spin(uint16_t i, uint16_t j);                      // returns the state of the spin at position (i,j)
    void gray2bgr();                                            // converts the grayscale image img to the blue-green-red image bgr
    void imshow();                                              // shows the blue-green-red image bgr in the window created in the constructor
    void destroyWindow();                                       // closes the window created in the constructor
    void vidwrite();                                            // appends the frame bgr to the videofile
    void vidrelease();                                          // saves and closes the videofile
    float get_magnetization();                                  // returns the magnetization of the current state
    float get_energy();                                         // returns the energy of the current state
  protected:
    uint16_t idx(int32_t x);                                    // index helper function for periodic boundary conditions
    void invert_spin(uint16_t i, uint16_t j);                   // inverts the spin at (i,j)
    void set_spin(uint16_t i, uint16_t j, bool newspin);        // sets the spin at (i,j)
    std::mt19937 rng;                                           // 32-bit Mersenne Twister pseudo-random generator
    std::uniform_int_distribution<uint16_t> int_distribution;   // converts the 32-bit random numbers to integer range
    std::uniform_real_distribution<double> real_distribution;   // converts the 32-bit random numbers to real interval
    cv::Mat bgr;                                                // blue-green-red image used to display the spinsystem and information
    std::ofstream datafile;                                     // datafile used to log the evolution of the configuration
  private:
    const uint16_t length = ARRAY_LEN;                          // length of the system
    bool spin[ARRAY_LEN][ARRAY_LEN];                            // state of the spinsystem
    uint8_t spinimg[ARRAY_LEN][ARRAY_LEN];                      // uint8_t representation of the spinsystem for the grayscale image img
    cv::Mat img;                                                // grayscale image based directly on the above array spinimg
    std::string videofilename;                                  // name of the datafile
    std::string datafilename;                                   // name of the videofile
    cv::VideoWriter video;                                      // tool to append frames to a video
};

template <uint16_t ARRAY_LEN>
configuration<ARRAY_LEN>::configuration(std::string _filename, float bias) : rng(std::random_device{}()) , int_distribution{0,ARRAY_LEN-1} , real_distribution{0.0,1.0} , datafile(_filename+".dat",std::ofstream::out) , img(ARRAY_LEN,ARRAY_LEN,CV_8U,spinimg) , video(_filename+".mkv",cv::VideoWriter::fourcc('X','2','6','4'),30, cv::Size(ARRAY_LEN,ARRAY_LEN+(ARRAY_LEN >= 200)*ARRAY_LEN/15))
{
  videofilename = _filename+".mkv";
  datafilename = _filename+".dat";
  std::uniform_real_distribution<float> biased_distribution(0.0,1.+1./bias);
  for(uint16_t i = 0; i < ARRAY_LEN; i++)
  {
    for(uint16_t j = 0; j < ARRAY_LEN; j++)
    {
      this->spin[i][j] = (int) biased_distribution(rng);
      this->spinimg[i][j] = this->spin[i][j]*255;
    }
  }
  Display(cv::namedWindow(videofilename,cv::WINDOW_NORMAL));
  cv::cvtColor(img,bgr,cv::COLOR_GRAY2BGR);
}

template <uint16_t ARRAY_LEN> void
configuration<ARRAY_LEN>::gray2bgr()
{
  cv::cvtColor(img,bgr,cv::COLOR_GRAY2BGR);
  if(ARRAY_LEN >= 200) bgr.push_back(cv::Mat(cv::Size(ARRAY_LEN,ARRAY_LEN/15), CV_8UC3, cv::Scalar(0,0,0)));
}

template <uint16_t ARRAY_LEN> void
configuration<ARRAY_LEN>::imshow()
{
  Display(cv::imshow(videofilename,bgr));
}

template <uint16_t ARRAY_LEN> void
configuration<ARRAY_LEN>::destroyWindow()
{
  Display(cv::destroyWindow(videofilename));
}

template <uint16_t ARRAY_LEN> void
configuration<ARRAY_LEN>::vidwrite()
{
  video.write(bgr);
}

template <uint16_t ARRAY_LEN> void
configuration<ARRAY_LEN>::vidrelease()
{
  video.release();
}

template <uint16_t ARRAY_LEN> bool
configuration<ARRAY_LEN>::get_spin(uint16_t i, uint16_t j){
  return spin[i][j];
}

template <uint16_t ARRAY_LEN> float
configuration<ARRAY_LEN>::get_magnetization(){
  uint64_t sum = 0;
  for(uint16_t i = 0; i < ARRAY_LEN; i++){
    for(uint16_t j = 0; j < ARRAY_LEN; j++){
      sum += spin[i][j];
    }
  }
  return -1.+2.*((float) sum)/(((uint32_t) ARRAY_LEN) * ((uint32_t) ARRAY_LEN));
}

template <uint16_t ARRAY_LEN> float
configuration<ARRAY_LEN>::get_energy(){
  int64_t sum = 0;
  for(uint16_t i = 0; i < ARRAY_LEN; i++){
    for(uint16_t j = 0; j < ARRAY_LEN; j++){
        sum += (this->get_spin(i,j)-!this->get_spin(i,j))*(this->get_spin(this->idx(i+1),j)-!this->get_spin(this->idx(i+1),j)+this->get_spin(i,this->idx(j+1))-!this->get_spin(i,this->idx(j+1)));
    }
  }
  return -((float) sum)/(((uint32_t) ARRAY_LEN) * ((uint32_t) ARRAY_LEN));
}

template <uint16_t ARRAY_LEN> void
configuration<ARRAY_LEN>::set_spin(uint16_t i, uint16_t j, bool newspin){
  spin[i][j] = newspin;
  img.at<uchar>(i,j) = 255;
}

template <uint16_t ARRAY_LEN> void
configuration<ARRAY_LEN>::invert_spin(uint16_t i, uint16_t j){
  spin[i][j] = !spin[i][j];
  spinimg[i][j] = spin[i][j]*255;
}

template <uint16_t ARRAY_LEN> uint16_t
configuration<ARRAY_LEN>::idx(int32_t x)
{
  return (ARRAY_LEN + x % ARRAY_LEN) % ARRAY_LEN;
}

#endif
