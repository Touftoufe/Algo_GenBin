# Algo_GenBin

The hole project is here, working well, i will comment it later.

# ALL THE FILES HAVE TO BE IN THE SAME FOLDER
the file result is not used by the program

# Parameters
  - The cross rate = 0.2
  - The mutation rate = 0.015
  - Number of chromosomes = 100
  - The first generation is randomly created
  - Error < 0.204
  
# Our constants:
  - lambda_0 E [6555, 6570] coded in 19 bit integer ~ [0, 524287] <=> [6555, 6607,4287]
  - lambda_L E [7, 15.191] coded in 13 bit integer ~ [0, 8191] <=> [7, 15.191]
  - Y_0 E [0.2, 0.25] coded in 16 bit integer ~ [0, 65535] <=> [0.2, 0.265535]
  - Y_max E [0.72, 0.785535] coded in 16 bit integer ~ [0, 65535]

# Ploting the curves
  The programm uses Matlab or Octave to plot the curves, you have the possibility to choose which one to use during the execution.
  
  NOTE: Make sure that you put the right path to the octave bin folder on the nongenfun.h, by default it is "C:\\Octave\\Octave-4.4.1\\bin"
