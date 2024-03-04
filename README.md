# MC-protein-folding

This repo is for the "Monte Carlo Simulations of Protein Folding" project. This project is my final project in the computational physics course at Brown University, taken in Spring 2018.

Protein folding is a prominent subject in both biology and biophysics, captivating the attention of numerous researchers. While Molecular Dynamics offers a powerful tool for probing this problem, its computational demands escalate significantly with larger systems. In such cases, the Monte Carlo method provides a swift albeit approximate solution, giving hints for the system's dynamics. In this project, I employ the Monte Carlo method to simulate protein folding processes, focusing on two crucial parameters: the system's energy and the end-to-end distance, which measures the distance between the first and last amino acids.

The main conclusions are as follows:

  - Protein folding tends to be more pronounced at lower temperatures.
  - Phase transition phenomena occur within the protein system, with the transition temperature independent of the number of amino acids (N).
  - At sufficiently high temperatures, the protein system converges towards a self-avoiding random walk configuration. 

## Code List
Below is the list explaining the purpose of each code file:

- simulating protein folding: protein.py
- analyzing the simulation data: analysis_equilibrium.py
- plotting the evolution of the energy and length: plot_equilibrium.py
- plotting the phase diagram: plot-phase.py
- simulation & analysis & plotting when temperature is low: low_temperature.py

## Report & Presentations
- The proposal for this project:  "proposal_YangruiHu.pdf"
- The final report for this project: "FinalReport_Yangrui.pdf"
- The slides for the presentation: "MC_protein_folding.pdf"
