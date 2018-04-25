# S4 - Rigorous Coupled Wave Analysis

Solcore’s capacity for modelling periodic nanophotonic structures is provided through an interface with the S4 Python extension, an implementation of RCWA (Rigorous Coupled Wave Analysis) developed in the Fan Group in the Stanford Electrical Engineering Department. The documentation for S4 can be found [here](http://web.stanford.edu/group/fan/S4/). The basic mechanics of the RCWA module are:


* Solcore translates the relevant inputs (list) into the appropriate notation for S4
*	Solcore calls S4 and feeds in the relevant inputs
*	Solcore translates the outputs from S4 into the relevant outputs, which match the outputs from Solcore’s other optics methods (absorption, reflection, transmission, absorption profile with depth).

To use S4 from within Solcore, first [make sure S4 is installed from the custom branch compatible with Python3](https://dalonsoa.github.io/solcore5/html/installation.html#installing-solcore). `Layers` in the `SolarCell` object are defined in the usual Solcore way, but now have an additional `geometry` attribute, e.g.:

```python
# define materials
Air = material('Air')(T=T)
TiO2 = material('TiO2', sopra=True)(T=T)  # for the nanoparticles
# define a layer with circular discs
NP_layer = [Layer(si('50nm'), Air, geometry=[{'type': 'circle', 'mat': TiO2, 'center': (200, 200), 'radius': 50}])]
```

The geometry attribute is a **list of dictionaries** containing relevant entries. You can add more than one shape per layer by simply adding more dictionaries to the list; each dictionary defines one shape, which is assumed to be periodic in two directions. The necessary information to define shapes is:

For all shapes:

- ‘type’: 'circle', 'ellipse', 'rectangle' or 'polygon'
- ‘mat’: the material the shape is made of; a Solcore material object.
- ‘center’: a tuple giving x and y coordinates (in nm) of the centre of the shape: `(x, y)`

Shape-dependent:

- Circle:
  - 'radius': a number in nm


- Ellipse: 
  - ‘angle’: a number in degrees, defining the angle by which the x-axis of the shape should be rotated (counter-clockwise).
  - ‘halfwidths’: a tuple of halfwidths in the *x* and *y* directions: `(hw_x, hw_y)`
  
-	Rectangle: ‘angle’ and ‘halfwidths’, as before

- Polygon: 
  - ‘angle’ as before
  - ‘vertices’: a tuple of tuples; each entry in the outer tuple are the x- and y-coordinates of the vertices of the (unrotated) polygon, one after another, in counter-clockwise order, e.g. `((x1, y1), (x2, y2), (x3, y3))`. Coordinates are in nm.

Additionally, you should set:

- the size (periodicity) of the unit cell in the *x* and *y* directions (i.e. the in-plane directions) in the options:

```python
opts.size = [400, 600] # [x-periodicity, y-periodicity]
```

- the number of Fourier orders to keep in the calculations in `opts.orders`. The calculation should converge for a higher number of orders, but computational time increases dramatically with the number of orders.

where `opts` is the Options object you will pass to e.g. `solar_cell_solver`. Note that **all dimensional information for the size and geometries should be in nm**!

## Description of functions

### optics/rcwa.py


`solve_rcwa(solar_cell, options)`

- very similar to solve_tmm, but calls relevant RCWA functions instead. Note that the RAT variable calculated using calculate_rat_rcwa is an input into calculate_absorption_profile (which, obviously, calculates the depth-dependent absorption profile). 

- Parameters:
  - solar_cell - a SolarCell object
  - options - options for the solver, e.g. wavelengths, number of Fourier orders, unit cell size.
  
- Returns: none - the function will modify the solar_cell object passed to it, updating each layer with the attributes diff_absorption (fraction of photons absorbed per unit length as a function of the position and the wavelength) and absorbed (total absorption in the layer). The total reflected, transmitted and absorbed intensity of the whole solar cell are also added to the solar_cell object.


`absorbed(self, z)`

-  same as TMM

`calculate_absorption_tmm(tmm_out)`

- `diff_absorption(z)`

- same as TMM


### absorption_calculator/rigorous_coupled_wave.py:

`calculate_rat_rcwa(structure, size, orders, wavelength, theta, phi, pol)`

- Calculates the total reflected, absorbed and transmitted intensity of the whole solar cell structure for the wavelengths and angles defined using an RCWA method implemented using the S4 package. It can handle different incident angles theta and phi and polarizations ('s', 'p' or 'u', the latter meaning unpolarised). It first creates an S4 simulation object using `initialise_S` and sets the incident angle and polarization using S4's SetExcitationPlaneWave, then loops through the wavelengths (updated in S4 using SetFrequency), using update_epsilon to update the optical constants of all the materials. It calls `rcwa_rat` to actually calculate R, A and T at each wavelength.

- Parameters:
  - structure: A solcore SolarCell object with layers and materials or a OptiStack object.
  - options: options for the calculation 
  - wavelength: Wavelengths (in nm) in which calculate the data.
  - angle: Angle (in degrees) of the incident light. Default: 0 (normal incidence).
  - pol: Polarisation of the light: 's', 'p' or 'u'. Default: 'u' (unpolarised).
  
- Returns:
  - output: A dictionary with the R, A and T at the specified wavelengths and angle.

`rcwa_rat(S, n_layers)`

- This is called from `calculate_rat_rcwa` and calculates R and T using S4's `GetPowerFlux` function which gives the forward and backward Poynting vector; the sum gives the effective real & imaginary power flux. T is given by the power flux entering the transmission medium, while R is given by 1 - (power flux entering top of solar cell).

- Parameters:
  - S: an S4 Simulation object
  - n_layers: number of layers in the stack, including incidence and transmission medium.

- Returns: a dictionary with R and T at a specific wavelength (A can be calculated from this since R + A + T = 1)

`initialise_S(stack, size, orders)`

- This is called from `calculate_rat_rcwa` and `calculate_absorption_profile_rcwa` and creates the S4 simulation object (class S4.Simulation) from the Solcore information, as well as OptiStack objects which are used to get the optical constants and thicknesses of the base layers, and the optical constants of the shapes defined in the geometry attribute of the layers. It does the following:
  - Initialises the S4 simulation object with basis vectors of the correct length and the desired number of Fourier orders using S4.New.
  - Makes a list of all the geometry objects contained in the Layers of the SolarCell object stack, and passes that to `necessary_materials ` to identify the materials used in the shapes for which we need optical constant information. The OptiStack function ignores the geometry attribute, so if we ignored this step we would not have optical constant data for materials which are only used in the geometry attribute and not in a base layer.
  - Makes an OptiStack object called stack_OS from the SolarCell stack itself. 
  - S4 uses named materials, so Solcore creates the required number of materials (S.SetMaterial) and creates the base layers (S.AddLayer), then loops through each Layer's geometry list to add the shapes to the layer.

- Parameters:
  - stack: a Solcore SolarCell object
  - size: a list of length 2 with the periodicity of the structure in the x and y directions
  - orders: the number of Fourier orders to retain in the RCWA calculation.
    
- Returns: 
  - S: a new S4 simulation object
  - stack_OS: an OptiStack object (of the Solcore stack)
  - shape_mats_OS: an OptiStack object, containing only the materials used in geometry objects. This is just a workaround to easily get the optical constants for those materials – the order of the layers doesn’t represent an actual physical structure.

`necessary_materials(geom_list)`

- Called from initialise_S; this simply extracts all the materials used in the geometries so that they can be turned into an OptiStack object, allowing the optical constants to be accessed easily.

- Parameters: 
  - geom_list: a list of all the geometry attributes for all the layers in the stack.

- Return: a list of the names of materials used in the layer geometries

`update_epsilon(S, stack_OS, shape_mats_OS, wl)`

- This updates the materials in the S4 simulations structure S using the S4 function S.SetMaterial using the information contained in stack_OS (for the base layers) and shape_mats_OS (for the shapes within those layers).

- Paramaters:
  - S: S4 Simulation object
  - stack_OS: OptiStack object of the base layers
  - shape_mats_OS: OptiStack object of the materials used in the geometry shapes
  - wl: the wavelength for which the optical constants should be set.

- Return: the S4 simulation object S with updated optical materials (i.e. optical constants) in order to loop through the wavelengths

`calculate_absorption_profile_rcwa(structure, size, orders, wavelength, rat_output, z_limit=None, steps_size=2, dist=None, theta=0, phi=0, pol='u')`

-	Mostly copied from TMM code, but calls `rcwa_position_resolved` to calculate the absorption profile (depth-dependent) in the stack. Like `calculate_rat_rcwa`, it creates an S4 simulation object and gets the optical constants from OptiStack objects using `initialise_S`, sets the parameters for the incident light using the S4 function SetExcitationPlaneWave and SetFrequency. It loops through the wavelengths, and at each wavelength it loops through the z-coordinates in dist. It then has to identify which Layer in the structure we are looking at, and how deep inside that layer the z-coordinate is, which is done using the function find_in_structure_with_inf from the tmm module. Once this has been identified, `rcwa_position_resolved` is called. 

- Parameters:
  - structure: SolarCell object
  - size: [x, y] size of unit cell (x and y periodicity)
  - orders: Fourier orders to retain in calculation
  - wavelength
  - rat_output: output from the calculate_rat_rcwa function. This is used to normalise the differential absorption.
  - z_limit: depth up to which absorption profile should be calculated (default None, if None, set to thickness of cell inside function)
  - steps_size: step size in nm at which to calculate depth-dependent absorption (default 2)
  - dist: array of z-locations at which to calculate depth-dependent absorption (default None, if None, calculated inside function from z_limit and steps_size.
  - theta: incident polar angle in degrees
  - phi: incident azimuthal angle in degrees
  
-	Return: output, an array describing depth-dependent absorption. First index is the wavelength, second index is the depth. 

`rcwa_position_resolved(S, layer, depth, A)`

- This function manually differentiates the absorption with respect to depth. The absorbed power across a small interval is calculated using S4's GetPowerFlux function (like in `rcwa_rat`).-delta is an arbitrarily small distance, in order to do very simple numerical differentiation. This was chosen because it gives good results (consistent with TMM for planar structures, integrated A sums correctly).

- Parameters:
  - S: S4 Simulation object
  - layer: which layer of the structure we are looking at
  - depth: how deep in that layer we are
  - A: total absorption

- Return: the absorbed energy density in layer at depth.

