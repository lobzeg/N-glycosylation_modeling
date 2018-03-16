# N-glycosylation_modeling

This repository contains scripts for a stochastic (Gillespie algorithm based) chemical-kinetic model of the N-glycosylation process in the Golgi Apparatus.

The file simulation.m sets up a simulation, given the parameters of the system. It has two parts: the adjustment of the system (determining the concentrations of glycans available for each reaction in each compartment) and the actual simulation, including multiple runs of a single glycan through the system and recording the data about the final structure of the glycan. In the provided script, the glycan type (complex, high-mannose or hybrid) is the structural feature being tracked; however, any other structural feature of N-linked glycans can be tracked instead.

Golgi_sim_full.m runs a single glycan through the system; adjustment.m does the same, but with adjustment of concentrations of glycans available for each reaction.

Both golgi_sim_full.m and adjustment.m call available_reactions.m and perform_reaction.m. Available_reactions.m determines which reactions are available at a certain time point, depending on the Golgi compartment and structure of the glycan. Perform_reaction.m randomly selects one of the available reactions and performs the corresponding structural changes to the glycan.
