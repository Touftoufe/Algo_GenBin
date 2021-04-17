# Algo_GenBin

This is a school project in C language programming
The hole project is here, working well, with all the comments.

# ALL THE FILES HAVE TO BE IN THE SAME FOLDER
the file "result" is not used by the program

# Parameters
  - The cross rate = 0.2
  - The mutation rate = 0.015
  - Number of chromosomes = 100
  - The first generation is randomly created
  - Error < 0.204
  - max number of itterations 5000
  
  NOTE: the program is not stopped until it reaches an error < 0.204. If the number of itterations reaches 5000 (~ 20 s), it means that the program is stuck and is not converging to the solution, so the chromosomes are reinitialised using the function "init()"
  
# Our constants with their conversion rules:
  - lambda_0 E [6555, 6570] coded in 19 bit integer ~ [0, 524287] <=> [6555, 6607,4287]
  - lambda_L E [7, 15.191] coded in 13 bit integer ~ [0, 8191] <=> [7, 15.191]
  - Y_0 E [0.2, 0.25] coded in 16 bit integer ~ [0, 65535] <=> [0.2, 0.265535]
  - Y_max E [0.72, 0.785535] coded in 16 bit integer ~ [0, 65535]

# Ploting the curves
  The program uses Matlab or Octave to plot the curves, you have the possibility to choose which one to use during the execution.
  
  NOTE: Make sure that you put the right path to the octave bin folder on the nongenfun.h, by default it is "C:\\Octave\\Octave-4.4.1\\bin"
