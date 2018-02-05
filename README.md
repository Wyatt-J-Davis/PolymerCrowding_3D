# PolymerCrowding_3D
## Description
Monte Carlo simulation of polymer-nanoparticle mixture. Written primarily by Wei Kang Lim, this code has been
heavily modified by Wyatt Davis in Java utilizing the Open Source Physics library. This version collects data on conformation of self-avoiding 
walk polymer. Coarse-grained model approximates chain polymers as ellipsoids. Results of simulation to be published in Journal of Chemical Physics.

## Instructions
1. Place Source code in osp_project/src/org/opensourcephysics/sip/CPM/
2. Compile and run code from osp_project/ using bash script
3. Specify run parameters using GUI
4. The simulation will store raw polymer eigenvalue data in osp_project/data/ for five independant runs
5. It is up to the user to analyze raw geometric data in meaningful way. I recommend writing a custom data analysis script in python.

## Requirements
- JDK 1.8 (http://www.oracle.com/technetwork/java/javase/downloads/index.html)
- Open Source Physics Library (https://www.compadre.org/osp/)
