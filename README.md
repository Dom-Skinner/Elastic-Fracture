# Elastic-Fracture
This repository contains the code for the paper "Viscous control of shallow
elastic fracture: peeling without precursors" [1]

You give the code the value of lambda (and sometimes L) and you get out the
stationary solution with values of KI,KII and the field h.

/full_fluid_problem/ contains the code that solves for a fracture without a
dry tip. Scaled_K_of_c_march.m is the main script
/full_fluid_problem_h_advance/ contains the code that solves for a fracture
with a dry tip, and so makes use of L as a parameter. Scaled_K_of_c_march.m 
is the main script
/linear_perturbation_problem/ contains code which solves the linear 
perturbation problem. The main script is h_find.m

/Tex/ contains some early version of the write up.
/Graphs/ contains plotting scripts for the graphs used in the paper and others
/Data/ contains results from old saved runs

[1] https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/viscous-control-of-shallow-elastic-fracture-peeling-without-precursors/A24C1E0EE18F78A1FFC00BAC7E411F51