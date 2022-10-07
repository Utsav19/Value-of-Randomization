## The value of randomized strategies in distributionally robust risk averse network interdiction problems
This repository by [Utsav Sadana](https://utsav19.github.io/) and [Erick Delage](http://tintin.hec.ca/pages/erick.delage/)  contains the code to replicate the numerical experiments in the paper "[The value of randomized strategies in distributionally robust risk averse network interdiction problems](https://arxiv.org/abs/2003.07915)".


# Dependencies

* [MATLAB](https://www.mathworks.com/products/matlab.html)
* [YALMIP](https://yalmip.github.io/) 
* [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) 
* [Gurobi](https://www.gurobi.com/) v9.0 or above


# Description of functions


A description of the folders containing the functions used in the implementations is given below:

- sbb_our_model: Functions to implement spatial branch and bound algorithm coupled with columnn generation procedure (sbb-CG)
- Heuristic: Functions to implement the heuristic combined with CG
- Consgen: Functions to implement constraint generation algorithm
- Gurobi: Contains the code that uses Gurobi's bilinear solver and functions implementation of spatial branch and bound without CG
- Loizou_model: Functions to implement the model in [Distributionally Robust Game Theory](https://arxiv.org/abs/1512.03253) for the network interdiction problem using our CG algorithm
- network_instances: Contains csv files to generate the Sioux Falls network and Nobel-US network. Also, it contains the functions to randomly generate network instances and capacities of the networ


To get the tables and figures in the paper, run the following files:

- ConvergenceGridNetwork.m: Convergence of sbb-CG with random instances as well as with instances for which there is value in randomization
- Compare_sbb_heuristic.m: Comparison of computation times for heuristic and sbb-CG
- ConvergenceGridNetwork.m Comparison of Gurobi's bilinear solver, constraint generation algorithm, and our sbb algorithm without CG and with CG
- ConvergenceGridNetwork.m Comparison of Gurobi's bilinear solver, constraint generation algorithm, and our sbb without CG and with CG
- VRAM.m: value of risk-averse model
- 
