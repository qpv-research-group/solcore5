================
Device Structure
================

This module contains the functions necesary to build a sample structure that can be read by the DD solver, as well functions for saving and loading that structure to/from a file. Typically, you will need to use only the :ref:`main-functions`. 

.. _main-functions:

Main functions
--------------

.. py:function:: CreateDeviceStructure(name [, role='device', T=293, layers=[], comments='', repeat=1, substrate=DefaultMaterial, reflection=None])

	Creates a dictionary with subdctionaries and lists storing all the parameters and material information necesary to solve the Poisson - Drift-Diffusion equations.
	
	*Layers* is a list of objects of class *solcore3.Layer* or other dictionary created previously with this function. If a file containing the absorption coeficient of the layer is to be used (2 columns, 1st the wavelengths 2nd the absorption), it must be included here as an option (with keyword 'absorption_file', relative path) when creating the *solcore3.material* of that layer. Any property not defined explicitly in the definition of the materials for the layers and that can not be calculated by :literal:`solcore3` is taken from the default material (see section :ref:`default-material`). *substrate* must be a *solcore3.material*, otherwise the default material is used. 

	The path to an external file with the reflection of the sample (2 columns, 1st the wavelengths 2nd the reflection) can also be included. **Note**: the reflection is not saved in the *json* file nor externally, so this file should always go with the json file to be used. 

   	**Output**: dictionary with all the information of the structure.

	**Example 1:**
	::
		
		# Example 1
		# ---------
		# This example illustrates the creation of a p-i-n solar cell structure including AlGaAs window and back surface field layers. 
		# A custom absorption coefficient for the GaAs of the intrinsic layer is used. 
		
		import solcore3.PDD as PDD
		from solcore3 import Layer, material
		
		# First, we create the materials, overriding any default property we want, such as the doping or the absorption coefficient
		window      = material('AlGaAs')    (T=300, Na = 1e24, Al=0.8)
		p_GaAs      = material('GaAs')      (T=300, Na = 1e24)
		i_GaAs      = material('GaAs')      (T=300, absorption_file = 'MyGaAs_Absorption.dat')
		n_GaAs      = material('GaAs')      (T=300, Nd = 1e23)
		bsf         = material('AlGaAs')    (T=300, Nd = 1e24, Al=0.4)
		
		# Then, we create the structure, adding the layers as a list
		MyDevice = PDD.CreateDeviceStructure( 'TestDevice', T = 300, layers = [
		                            		Layer(width = 30e-9,         material = window,           role="Window"),
		                            		Layer(width = 400e-9,        material = p_GaAs,           role="Emitter"),
		                            		Layer(width = 400e-9,        material = i_GaAs,           role="Intrinsic"),
		                            		Layer(width = 2000e-9,       material = n_GaAs,           role="Base"),
		                            		Layer(width = 200e-9,        material = bsf,              role="BSF")
		                            		])
											

.. py:function:: SolveQWproperties(device[, calculate_absorption=True, wavelengths=None, periodic=True, filter_strength=0.0, blur=None, blurmode="left",use_kp=True, use_Adachi=False, alpha_params=None])
   
	Considers the device as a QW and solves its properties, including the modification of the bandeges due to strain, the efective mases and the absorption coefficient. Without calling this function, the structure is just a colection of layers with bulk-like properties. 
	
	**Output**: A dictionary with the output of the Schrodinger solver.  
	
	**Example 2:**
	::
		
		# Example 2
		# ---------
		# This example illustrates the creation of a QW structure with 40 wells, solve its quantum properties and 
		# add it to the intrinsic region of a p-i-n solar cell
		
		import solcore3.PDD as PDD
		from solcore3 import Layer, material
		
		# First, we create the materials of the QW
		QWmat       = material('InGaAs')    (T=300, In=0.2)
		Bmat        = material('GaAsP')     (T=300, P=0.1)
		i_GaAs      = material('GaAs')      (T=300)
		
		# The QW is 7 nm wide, with GaAs interlayers 2 nm thick at each side and GaAsP barriers 10 nm thick. 
		# The final device will have 40 of these QWs.
		QW = PDD.CreateDeviceStructure( 'QW', T = 300, repeat = 40, layers = [                            
		                            Layer(width = 10e-9,          material = Bmat,           role="barrier"),
		                            Layer(width = 2e-9,           material = i_GaAs,         role="interlayer"),
		                            Layer(width = 7e-9,           material = QWmat,          role="well"),
		                            Layer(width = 2e-9,           material = i_GaAs,         role="interlayer"),
		                            Layer(width = 10e-9,          material = Bmat,           role="barrier"),
		                                 ])

		# We solve the quantum properties of the QW, leaving the default values of all parameters
		PDD.SolveQWproperties(QW) 
		
		# We create the other materials we need for the device
		window      = material('AlGaAs')    (T=300, Na = 1e24, Al=0.8)
		p_GaAs      = material('GaAs')      (T=300, Na = 1e24)
		n_GaAs      = material('GaAs')      (T=300, Nd = 1e23)
		bsf         = material('AlGaAs')    (T=300, Nd = 1e24, Al=0.4)
		
		# And finally we create another p-i-n structure incorporating the QWs in the intrinsic region.
		MyDevice = PDD.CreateDeviceStructure( 'TestDevice', T = 300, layers = [
		                            Layer(width = 30e-9,         material = window,           role="Window"),
		                            Layer(width = 400e-9,        material = p_GaAs,           role="Emitter"),
		                            Layer(width = 10e-9,        material = i_GaAs,           role="Intrinsic"),
		                            QW,
		                            Layer(width = 10e-9,        material = i_GaAs,           role="Intrinsic"),
		                            Layer(width = 2000e-9,       material = n_GaAs,           role="Base"),
		                            Layer(width = 200e-9,        material = bsf,              role="BSF")
		                            ])
									

.. py:function:: Save(device, filename [, save_absorptions_individually=False, remove_absorption_from_json=False, override_absorption=False, directory='default', yaml=False])

	Saves the structure in a text file. The default file format is *json*, although *yaml* can also be selected. The file extension is added automatically to *filename*. The layout of the file is not very nice but it is human readable, cross platform and independent of the programing language. In Python, *json* is part of the standard library, contrary to *yaml* that has to be installed separately.	
	
	An option is given to save the absorption coefficients in individual files inside the default subdirectory (= 'filename_inputs'), as well as removing them from the *json/yaml* file, greatly improving the readability. 
	
	**Example 3:**
	::
	
		# Example 3
		# ---------
		# Following the previous example, we can save the properties of the QW and the final device. 
		# We choose to keep the absorption coeficients in separate files but leave the default directory.
		# The directory containing the absorption coefficients will be 'myQW_inputs' and the name of the files will have the following pattern:
		#
		#	<deviceName>_<layerNumber>_<layerComposition>_<layerRole>.dat
		#
		# The absorption of the well itself will be: 'QW_2_In0.2GaAs_well.dat'		
		PDD.Save(QW, 'myQW', save_absorptions_individually=True, remove_absorption_from_json=True)
		
		# And the same for the complete device:
		PDD.Save(MyDevice, 'MyDevice', save_absorptions_individually=True, remove_absorption_from_json=True)
		
		
.. py:function:: Load(filename[, yaml=False])

	Reads a structure from a *json/yaml* file created with *Save*. Other files created with a valid *json/yaml* syntax are accepted without any checking of whether they are valid structures or not. If they are not valid, they will produce an error later in the program . The safest option is to create the structure with *Save* or, at least, use an existing file with a structure as template in case you prefere doing it 'by hand'. Absorption coeficients from external files will be loaded, if present.
		
	**Output**: dictionary with all the information of the structure, as in *CreateDeviceStructure* above.
	
	**Example 4:**
	::
	
		# Example 4
		# ---------
		# Once saved, the structures can be used in other Python scripts. In this final example, 
		# we load 'MySample.json' of the previous example and change the thickness and doping of the emitter.
		
		import solcore3.PDD as PDD
		
		# NewDevice now contains the structure of MyDevice
		NewDevice = PDD.Load('MyDevice')
		
		NewDevice['layers'][1]['properties']['width'] = 200e-9
		NewDevice['layers'][1]['properties']['Na'] = 5e24
		
		# We save the structure of the new device
		PDD.Save(NewDevice, 'NewDevice_ThinEmitter', save_absorptions_individually=True, remove_absorption_from_json=True)
		

Other functions
---------------   

.. py:function:: AddLayers(device, layers)

	Adds a list of objects of class *solcore3.Layer* or a full device structure created with *CreateDeviceStructure* to the structure. They are appended to the end of the current tree of layers. *Layers* must be a list, even if it is only one layer. 

.. py:function:: RemoveLayer(device, i)

	Removes layer 'i' from the structure.
   
.. py:function:: GetLayerProperties(layer, T)

	Reads the properties of a *solcore3.Layer*. The output is a dictionary with all the properties of the layer - or the default properties, if not available. Is called internally by *AddLayers* above. It calls the *mobility* module to calculate the carrier mobilities if they are not included as an input. 
   
.. py:function:: LoadAbsorption(layer, T[, wavelengths=None])

	Loads the absorption of the layer, either from the corresponding file or calculating it from the material properties. It is called automatically when running a drift difusion experiments that requires the absorption. 
   
.. py:function:: InLineComposition(layer)

	Produces a string with the composition of the layer, e.g.: :code:`In0.2GaAs`
   
.. py:function:: SolcoreMaterialToStr(material_input)

	Given a *solcore3.material*, produces a dictionary with its composition. e.g.: :code:`{ 'material' = 'InGaAs', 'element' = 'In', 'fraction' = '0.2' }`. Currently it only works with binaries and ternaies alloys. 
   
.. py:function:: ToSolcoreMaterial(comp, T[, execute=False])
   
	Produces either a string with the defeinition of a *solcore3.material* (default) or a *solcore3.material* itself based on a dictionary with the composition of the material as produced by *SolcoreMaterialToStr*. 

.. py:function:: ToLayer(width, material, role)

	Produces a *solcore3.Layer* based on the temperature, composition and thicknesses (and nothing else) of a material, as produced by *ToSolcoreMaterial*. 
   
.. py:function:: ToStructure(device)

	Takes a list of *solcore3.Layers* and create an opbject of class *solcore3.Structure*. 

.. _default-material:

Default material
----------------

The default material is **GaAs at 293K**. In general, all the intrinsic properties of the materials are taken from the literature and calculated by the :literal:`solcore3` material system. If not found there (for example if the material or alloy is not in the database) or if they are extrinsic, the default values below are used. The extrinsic properties, that depend on the quality of the material or its doping, are just asigned a 'reasonable' value. The user must make sure that this 'reasonable' value is indeed resonable for the intended application and override it with his/her own, otherwise.

The total list of parameters and the default values are:

============================ ============================ ====================== =============================================
Parameter                    Default value                Units                  Notes
============================ ============================ ====================== =============================================
band_gap                     Calculated for GaAs at 293K  J                      \-
electron_affinity            Calculated for GaAs at 293K  J                      \-
eff_mass_electron_Gamma	     Calculated for GaAs at 293K  relative to m0         \-
eff_mass_hh_z                Calculated for GaAs at 293K  relative to m0         \-
eff_mass_lh_z                Calculated for GaAs at 293K  relative to m0         \-
permittivity                 12.9                         relative to epsilon0   \-
electron_mobility            Calculated for GaAs at 293K  m2 V-1 s-1             \-
hole_mobility                Calculated for GaAs at 293K  m2 V-1 s-1             \-
electron_minority_lifetime   3e-6                         s                      SRH [#f1]_ recombination time for electrons
hole_minority_lifetime       2.5e-7                       s                      SRH recombination time for holes
electron_auger_recombination 1e-42                        m6 s-1                 \-
hole_auger_recombination     1e-42                        m6 s-1                 \-
radiative_recombination	     7.2e-16                      m3 s-1                 \-
Nd                           1                            m-3                    Density of donors
Na                           1                            m-3                    Density of acceptors
sn                           1e6                          m s-1                  Surface recombination velocity for electrons
sp                           1e6                          m s-1                  Surface recombination velocity for holes
absorption_file              None                         \-                     \-
============================ ============================ ====================== =============================================

.. [#f1] Shockley-Read-Hall

All functions description
-------------------------

.. automodule:: solcore.poisson_drift_diffusion.DeviceStructure
    :members:
    :undoc-members:
