# MS3DATA-2D
Multiphase-flow Statistics using 3D Detection and Tracking Algorithm (2D Version) 

Reference Publication: 
Bakshi, A., Altantzis, C., Bates, R.B., and Ghoniem, A. F.,"Multiphase-flow Statistics using 3D Detection and Tracking Algorithm (MS3DATA): Methodology and application to large-scale fluidized beds.", Chemical Engineering Journal 293 (2016): 355-364 

I would also like to acknowledge the contribution of MIT '19 Anastasiya Bershanska towards Lagrangian tracking of bubbles

Description: 
MS3DATA is a scalable interphase detection and tracking algorithm. The tool is developed for bubble detection in solid-gas flows using open-source software MFiX, but can be easily extended to other applications of multiphase flow (e.g. detection of air bubbles in liquids, droplets, clusters and so on) as long as a threshold criterion for demarcating interphases is possible. The algorithm is indifferent to the data-source (i.e. experimental or simulation data) as long as data can be output in the format [frame#, cell#, field-var]. This version of the code is for 2D statistics only i.e. statistics in 2D simulations and/or slice-based statistics in 3D simulations. The 3D version of the code will be released soon!

Schematic:
In solid-gas flows, bubbles can be defined as areas with local void fraction > 0.7. The basic schematic of the code is as follows: 
Step 1. generate formatted output  [frame#, cell#, voidage] from CFD simulations (see sol_flux.f for more details).
Step 2. input geometric parameters in geometry.xlsx (to replicate simulation grid) 
Step 3. define bubble detection parameters in bubbles2D.m and run script 
Step 4. post-process matrix "bubblepropertiestotal" for suitable statistics 

Sample data for 300 time frames is provided in bubblestats2D.txt for step-by-step implementation and better understanding. 
Whether you have fun with the code or not (debugging!), feel free to reach out to me at abakshi@mit.edu with comments, questions and feedback!

