# MC-protein-folding

This repo is for the "Monte Carlo Simulations of Protein Folding" project. This project is my final project in the computational physics course at Brown University, taken in spring 2018.

Protein folding is a hot topic both in biology and biophysics. There are many scientists working on the simulation of protein folding system. Although Molecular Dynamics method is a good technique to probe this problem, the calculation is huge when the size of the system becomes large. In this case, Monte Carlo method can give us a rough result quickly and help us understand this system. In this project, I simulate the process of the protein folding with the Monte Carlo method and monitor two key quantities: the energy of the system and the end-to-end distance which is the distance between the first amino acid and the last one. Moreover, I study the effect of temperature, interaction coefficient, and the number of amino acids. 
I find the following conclusions: 

  - The extent of folding tends to be high in low temperature.
  - Phase transition is observed in the protein system and the transition temperature is independent of N.
  - The protein system approaches to self-avoiding random walk system when temperature is high enough.

## Code List
Below is the list explaining the goal of each code file. 

- simulating protein folding: protein.py
- analyzing the simulation data: analysis_equilibrium.py
- plotting the evolution of the energy and length: plot_equilibrium.py
- plotting the phase diagram: plot-phase.py
- simulation & analysis & plotting when temperature is low: low_temperature.py

## Report & Presentations
- The proposal for this project:  "proposal_YangruiHu.pdf"
- The final report for this project: "FinalReport_Yangrui.pdf"
- The slides for the presentation: "MC_protein_folding.pdf"
