# CS205 Final project
# Implementation and parallelization of Redfield equations

Florian Hase, Hannah Sim, and Teresa Tamayo

In this project, we are going to implement and parallelize a 
method for computing the time evolution of the density matrix in an open quantum, 
the Redfield method. This method applies to some photosynthetically active protein complexes[1,2]. 

<center>
<img src="files/FMO.png" width="200">
</center>

**Figure:** The Fenna-Matthews-Olson (FMO) complex.
Protein from green sulfur bacteria and involved in
the excitation energy transfer from light-harvesting chlorosomes.

The Theoretical and Physical Chemistry group of Prof. Alan Aspuru-Guzik, 
has wide experience in the field, where they have studied the
exciton energy transfer in diverse systems.


## <i class="fa fa-check-square" aria-hidden="true"></i> Redfield method

The evolution of small systems is usually influenced by the interaction of the surroundings, 
given that in general is impossible to isolate it. Hence, the dynamics of a quantum system depends 
substantially on the interaction of the external environment. 
However, we are usually unable to keep track the evolution of the complete systems and their surroundings. 
In this situation, we use equations that account the influence of the surroundings on the systems, 
but not keeping track the environment evolutions. 
These ubiquitous phenomena determine the physics and chemistry in diverse fields, one of the most important is the exciton transfer dynamics[1],
One way to describe these dynamics is to utilize the Redfield master equations, where we assume that the transport is a dissipative dynamics for the reduced excitonic density matrix[3,4]. This model considers the interactions between the environment and the system weak and that the system depends only on its present state. The density matrix of the excitonic system can be time evolved according to the Liouville-von Neumann equation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{d}{d&space;t}\rho(t)&space;=&space;\frac{-i}{\hbar}[H_s,\rho(t)]&space;&plus;&space;L&space;(\rho(t))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{d}{d&space;t}\rho(t)&space;=&space;\frac{-i}{\hbar}[H_s,\rho(t)]&space;&plus;&space;L&space;(\rho(t))" title="\frac{d}{d t}\rho(t) = \frac{-i}{\hbar}[H_s,\rho(t)] + L (\rho(t))," /></a>

where *L* denotes the Lindblad operator in the secular Redfield approximation. 

The first term takes into account the evolution of the system without the presence of the environment, while the second keeps track on the environmental influence of the system[3,4]. 

Solving this equation involves some matrices operations in large matrices, and parallelization in needed to have a better
representation of the system.








## <i class="fa fa-check-square" aria-hidden="true"></i>  Previous implementations
There are a wide variety of methods that can obtain the exciton dynamics, some of them
departs on the Redfield equations.
One straightforward implementation is build and diagonalize the Liouville superoperator,
and if we use a time-independent Hamiltonian, the compuational effor scales 
<a href="https://www.codecogs.com/eqnedit.php?latex=N^6" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N^6" title="N^6" /></a>
 where <a href="https://www.codecogs.com/eqnedit.php?latex=N^6" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N" title="N" /></a>, 
is the number of the basis of vibrational states.
This strategy can be numericall unstable in some case.
Finally, there are other approaches where by rewritting the Redfield equation and making futher
approximation, there is only needed matrix-matrix multiplications, hence the computational time scales
is 
<a href="https://www.codecogs.com/eqnedit.php?latex=N^6" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N^3" title="N^3" /></a>[5].



## <i class="fa fa-check-square" aria-hidden="true"></i>  Methodology
Our first approach will be implement a serial version of the
solver with a variety of numerical implementations to
perform the time propagations, including different numbers 
of vibrational states in the bath.
Later, we explore different parallelization schemes
with hybrid architectures.
Finally, we will compare our results with an state-of-the-art 
implementation of hierarchical equations of motion [6,7].




## <i class="fa fa-check-square" aria-hidden="true"></i>  Python implementations

To understand the secular Redfield approximation for propagating excitons under a given Hamiltonian in more detail we implemented a Python version of the Redfield method. With this implementation we were able to identify the scaling of the method and determine parts of the algorithm which are time consuming but suited for parallelization. 

In particular, we found that the secular Redfield approximation can be devided into two major parts. The Lindblad operator *L* contains a set of rates which model the transition of a particular excitonic state into another excitonic state or the decay of the excitonic state into the ground state or the target state. These rates are calculated based on the excitonic eigenstates of the Hamiltonian and the spectral density, for which parameters are expected from the user. 

Although computing the rates involves a number of matrix operations (diagonalizing the Hamiltonian, matrix multiplications, etc.) this is not the time limiting step for computing the exciton dynamics because rates can be calculated once at the beginning of the simulation. 

Instead, we found that we should focus our parallelization strategies on the propagation of the density matrix for precomputed Hamiltonians and Lindblad operators. Propagating the system in the excitonic eigenbasis yields the additional advantage of the Hamiltonian being diagonalized, which simplifies the computation of the commutator of the Hamiltonian and the density matrix (we save one matrix multiplication). 

However, the action of the Lindblad operator on the density matrix needs to be calculated at every iteration step. Even with precomputed Lindblad operator matrices, this operation requires at least four matrix multiplications, which is reflected in the scaling plots we obtained for the Python implementation of the code. We expect the algorithm to greatly benefit from a parallelization of this operation. 

<center>
<img src="files/runtimes.png" width="200"><img src="files/runtimes_loglog.png" width="200">
</center>

We found that both implementations of the secular Redfield method show a scaling of roughly 
<a href="https://www.codecogs.com/eqnedit.php?latex=N^6" target="_blank"><img src="https://latex.codecogs.com/gif.latex?N^3.5" title="N^3.5" /></a> with *N* the number of excitonic sites in the system.



## <i class="fa fa-check-square" aria-hidden="true"></i>  References
[1] P. Rebentrost, R. Chakraborty, and A. Aspuru-Guzik, *J. Chem. Phys.*, **131**, 184102 (2009).

[2] P., M. Mohseni, I. Kassal, S. Lloyd, and A. Aspuru-Guzik, *New Journal of Physics*, **11**, 033003 (2009).

[3] V. May and O. Kühn, *Charge and Energy Transfer Dynamics in Molecular Systems*, Wiley, New York, 2004. 

[4] H. P. Breuer, and F. Petruccione, *The theory of open quantum systems*, Oxford, New York, 2010.

[5] I. Kondov, U. Kleinekathöfer, and M. Schreiber, *J. Chem. Phys.*,  **114**, 1497 (2001).

[6] C. Kreisbeck and T. Kramer, *J. Phys. Chem.* **3**, 2828, (2011).

[7] C. Kreisbeck and T. Kramer "Exciton Dynamics Lab for Light-Harvesting Complexes (GPU-HEOM)," https://nanohub.org/resources/gpuheompop, (DOI: 10.4231/D3QB9V630), (2014).
