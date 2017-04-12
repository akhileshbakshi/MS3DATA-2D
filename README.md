# MS3DATA-2D
Multiphase-flow Statistics using 3D Detection and Tracking Algorithm (2D Version) 

Reference Publication: 
Bakshi, A., Altantzis, C., Bates, R.B., and Ghoniem, A. F.,"Multiphase-flow Statistics using 3D Detection and Tracking Algorithm (MS3DATA): Methodology and application to large-scale fluidized beds.", Chemical Engineering Journal 293 (2016): 355-364 

Description: 
MS3DATA is a scalable interphase detection and tracking algorithm. The tool is tailored for bubble detection in solid-gas flows using open-source software MFiX, but can be easily extended for detection of air bubbles in liquids, droplets, clusters and so on in other areas of multiphase flow, as long as a critical threshold can be defined for isolating interphases/areas of interest. The algorithm is agnostic to the source of data (i.e. experimental or simulation data) as long as data can be output in the format [frame#, cell#, field-var]. This version of the code is for 2D statistics only i.e. statistics in 2D simulations and/or slice-based statistics in 3D simulations. The 3D version of the code will be released soon!

Schematic: 
Step 1. Generate formatted output [frame#, cell#, field-var]. sol_flux.f provided is tailored for MFiX simulations
Step 2. Edit geometry.xlsx with geometric parameters to replicate grid of simulation data 
Step 3. Define postprocessing parameters in bubbles2D.m and run simulation! 

Sample formatted data for 300 time frames is provided in bubblestats2D.txt for better understanding. 
Whether you have fun with the code or not (debugging!), feel free to reach out to me at abakshi@mit.edu with commends, questions and feedback!

