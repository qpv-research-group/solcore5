===============
QW unit creator
===============

This module defines a class derived from *solcore3.Structure* that allows to solve the Schrodinger equation and the kp model. It also prepares the properties of the structure (bandedges, efective density of states (DOS), etc) in order to have a meaningful set of properties for the PDD. Without this preparation, the structure is just a colection of layers with bulk-like properties, as it is illustrated in the figure:

.. _figure-QWunit: 

.. figure:: Figures/QWunit.png
	:align: center

	Figure 1

The QWunit class
----------------

.. py:class:: QWunit (*args, **kwargs)

	Asembles a group of layers as a quantum well, calculating its properties as a whole. As with any **solcore3.Structure**, a QWunit usually gets as input a list of **solcore3.Layers** and a **solcore3.material** as substrate (defined with the keyword 'substrate'). In this case, the class constructor is disigned to work specifically for the PDD solver so it requires the following:

	- It needs a minimum of 3 layers.
	- There must be exactly one layer with the role = "well".
	- Top and bottom layers are assumed to be the barriers.
	- Anything not a barrier nor a well will be considered as "interlayer".
	- A substrate.

	Since the PDD solver can not work with superlatices, there is no point of considering the solution of more than one QW. All QWs entering into the solver are independent, regardless of the width of the barriers. To account for the possible presence of nearby QWs, perdiodic boundary conditions can be used to solve the Schrodinger equation. If barriers + interlayers are thick enought, this make no difference but if they are thin, it affects the number and energy of the confined levels. It does not 'create' minibands, though, so the results with thin barriers must be used with caution. 
	
	
	.. py:method:: solve([V=0, WLsteps=(300e-9, 1100e-9, 201), wavelengths=None, periodic=True, filter_strength=0.0, blur=None, blurmode="left", offset=0,use_kp=True, use_Adachi=False, calculate_absorption=True, alpha_params=None, T=293])
	
		Solves the structure, calculating the energy levels, the absorption, etc. First, it calls the Schrodinger solver (**solcore3.qm.schrodinger**) and then it uses the result to obtain an efective band profile for the conduction and valence band, an efective density of states considering both, the bulk and quantum levels, and the absorption coefficient. This is done for each layer of the structure (see the following methods for details). The input parameters of this method are the same that for the Schrodinger solver except for the wavelengths definition that can be used as default (WLsteps) or be added with 'wavelengths'.
		
		The method returns as output the output dictionary of the Schrodinger solver.
		
	.. py:method:: RecalculateBandEdges(use_kp, SR)
	
		From the perspective of the PDD solver, the actual bandgap and electron affinity of each layer depend on the energy levels, i.e. the minimum energy for electrons is not the band edge of the conduction band, but the ground confined level. The same applies to holes, being the actual band edge the maximum between the ground states of light holes and heavy holes (see :ref:`figure-QWunit`). 
		
		In this method we modify that by creating an effective electron affinity and band gaps for all layers in the QW. For the barriers, the electron affinity and band gap are the same than in bulk, modified by the kp calculation, if necesary. For interlayers, it depends on what is higher, the bandedges of the interlayer or the confined carrier levels.
		
		It requires as input kp should be used to recalculate the band positons and effective mases ('use_kp=True') and the output of the Schrodinger solver.  
		
	.. py:method:: RecalculateDensityOfStates()
	
		Calculates the effective density of states for each layer in the QW. The general rule is:
		- Barriers have the bulk density of states
		- QW have ALL the density of states asociated with the confined states + bulk density of states above the barrier
		- Interlayers have only the bulk density of states above the barrier

		This simplification is similar to that in Nelson et al. [#Ref2]_ and allow us to keep the bulk-like form of the carrier densities in the drift diffusion equations under the Boltzmann aproximation. From a physical porint of view, it probably can be done better.
		
	.. py:method:: CalculateAbsorption(use_Adachi, SR):
	
		If required, this function calculates the absorption of the QW, putting together the absorption of the confined levels and the absorption of the bulk. As with the density of states, the rules are:
		
		- Barriers have the bulk absorption
		- Interlayers have the bulk absorption from the barrier energy and zero below that
		- Wells have the absorption of the confined levels below the barrier energy and of the bulk above it. 
		
		The calculation is similar to that in [#Ref3]_ but, as with the DOS, it can probably be improved. 
		
		The required input is if Adachi method [#Ref1]_ for calculating the absorption should be used ('use_Adachi=True') and the output of the Schrodinger solver. 
		
All functions description
-------------------------

.. automodule:: solcore.poisson_drift_diffusion.QWunit
    :members:
    :undoc-members:
		
References
----------

.. [#Ref2] J. Nelson, M. Paxman, K. W. J. Barnham, J. S. Roberts, and C. Button, “Steady-state carrier escape from single quantum wells,” IEEE J. Quantum Electron., vol. 29, no. 6, pp. 1460–1468, 1993.

.. [#Ref3] C. I. Cabrera, J. C. Rimada, J. P. Connolly, and L. Hernandez, “Modelling of GaAsP/InGaAs/GaAs strain-balanced multiple-quantum well solar cells,” J. Appl. Phys., 113, 024512, (2013).

.. [#Ref1] S. Adachi “Optical dispersion relations for GaP, GaAs, GaSb, InP, InAs, InSb, AlxGa1−xAs, and In1−xGaxAsyP1−y,” J. Appl. Phys.,66, 6030 (1989).