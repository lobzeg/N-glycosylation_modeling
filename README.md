# N-glycosylation_modeling

This repository contains scripts for a stochastic (Gillespie algorithm based) chemical-kinetic model of N-glycosylation process 
in the Golgi Apparatus.

File simulation.m sets up a simulation for given parameters of the system. It has two parts: adjustment of the system 
(determining concentrations of glycans available for each reaction in each compartment) and actual symulations which includes
multiple runs of a single glycan through the system and recording of the data about final structure of a glycan. In the provided 
script type of glycan (complex, high-mannose or hybrid is the structural feature being tracked, however any other structural 
features of N-linked glycans can be tracked instead.

Golgi_sim_full.m runs a single glycan through the system, adjustment.m does the same, but with adjustment of concentrations
of glycans available for each reaction.

Both golgi_sim_full.m and adjustment.m call available_reactions.m and perform_reaction.m.
Available_reactions.m determines which reactions are available at the certain time point, depending on the Golgi compartment 
and structure of a glycan. Perform_reaction.m randomly selects one of available reactions and performs corresponding 
structural changes to a glycan.
