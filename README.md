## Installation

### Dependencies

#### Install the dependencies (Archlinux)
	sudo pacman -S cmake fmt ffmpeg vtk qt5-base opencv
#### Install the dependencies (Ubuntu)
	sudo apt-get install cmake tcl-vtk qt5-default libgtk3.0-dev ffmpeg
From within Ubuntu it is best to compile OpenCV 4.4.0 from source as well as fmtlib.

## Usage
1. Setup cmake: `cmake .`
2. Choose the desired system length in line 11 of main.cpp.
3. Compile the program: `make`
4. Run the program using the following syntax:
   * Option 1: `./main basename temperature`
   * Option 2: `./main basename temperature_start temperature_end temperature_step`
   * Option 3: `./main basename temp1 temp2 temp3 temp4 ...`

## Wiki
An in-depth discussion of the code and results that can be achieved with it can be found [here](https://theoreticalphysics.info/index.php/2D_Ising_Model:_Monte_Carlo_Simulations_using_the_Metropolis_Algorithm).
