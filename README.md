## The value of randomized strategies in distributionally robust risk averse network interdiction problems
This repository by [Utsav Sadana](https://utsav19.github.io/) and [Erick Delage](http://tintin.hec.ca/pages/erick.delage/)  contains the code to replicate the numerical experiments in the paper "[The value of randomized strategies in distributionally robust risk averse network interdiction problems](https://arxiv.org/abs/2003.07915)".


# Dependencies

* [MATLAB](https://www.mathworks.com/products/matlab.html)
* [YALMIP](https://yalmip.github.io/) 
* [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) 
* [Gurobi](https://www.gurobi.com/) v9.0 or above


Description  


A description of the folders containing the functions used in the implementations is given below:

- sbb_our_model: Functions to implement spatial branch and bound algorithm coupled with columnn generation procedure
- Consgen: Functions to implement constraint generation algorithm
- Gurobi: Contains the code that uses Gurobi's bilinear solver and implementation of spatial branch and bound without column generation
- Loizou_model: Functions to implement (https://arxiv.org/abs/1512.032530 for our network interdiction problem using our column generation algorithm
- network_instances: Contains csv files to generate the Sioux Falls network and Nobel-US network. Also, it contains the functions to randomly generate network instances and capacities of the network
