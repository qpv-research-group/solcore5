Defining new materials
=======================

.. toctree::
    :maxdepth: 1

    Materials

There are two main ways to use a material which is not built in
to the Solcore database, which can also be connected. They are:

1. Downloading and using the database from `refractiveindex.info <https://refractiveindex.info/>`_
2. Providing n and k data, and other parameters, to ``create_new_material``

In order to control where custom materials get saved, you need to tell Solcore where to
create and look for new materials by adding some entries to your user configuration file
(by default, a hidden folder called .solcore_config.txt in your home directory):

- Path the refractiveindex.info database will be downloaded to is set under [Others] with flag ``nk``.
- Path where n and k data will be saved is set under [Others] with flag ``custom_mats``.
- Path where the file containing parameters of custom materials will be created is set under [Parameters] with flag ``custom``. 

The following code snippet sets the location for each of these within a folder called Solcore, which
is a sub-directory of your home folder (you could also manually add the correct paths to the config
file).

.. code-block:: Python

    import os
    from solcore.config_tools import add_source
	
    home_folder = os.path.expanduser('~')
    custom_nk_path = os.path.join(home_folder, 'Solcore/custommats')
    nk_db_path = os.path.join(home_folder, 'Solcore/NK.db')
    param_path = os.path.join(home_folder, 'Solcore/custom_params.txt')

    add_source('Others', 'custom_mats', custom_nk_path)
    add_source('Others', 'nk', nk_db_path)
    add_source('Parameters', 'custom', param_path)
	

Adding new materials to the database
-------------------------------------

.. automodule:: solcore.material_system.create_new_material
    :members:
    :undoc-members:
	
The correct format for the n and k files are tab or space separated text files,
with the first column the wavelength **in metres** and the second column n or k.
The file containing other parameters should take 

When you add a new material to the database, a new folder will be created for it
as a subfolder of the path specified in your (user or default) configuration under
the ``custom_mats`` path. The n and k data files you provide will be copied into that folder
(they are renamed automatically). Any other parameters you supply will be copied into
the file specified under the ``custom`` path.

Using refractiveindex.info
---------------------------

.. automodule:: solcore.absorption_calculator.nk_db
    :members:
    :undoc-members:

Before the first use, you will need to download the database:

.. code-block:: Python

	from solcore.absorption_calculator download_db
	download_db()

``download_db`` takes two (optional) arguments: the URL of the database to be downloaded,
(default is the most recent, hardcoded into the function), and how many interpolation points 
to use when saving the database (default is 200).

The code which is used to download and get data from the refractiveindex.info
database is based on code from `Hugo Guillen
<https://github.com/HugoGuillen/refractiveindex.info-sqlite>`_.

You can now directly use materials from the .db file created this way by 
referencing them via their **pageid**. To locate which database entry you want to use,
you can search the database; the code below searches the database for entries matching
'Diamond' and then uses the pageid of the first result to create an instance of this
new Diamond material. In general, though, it is a good idea to check explicitly which
of the database entries is appropriate (e.g. in terms of the wavelength range, and which type
of data is available) rather than simply using the first result.

.. code-block:: Python

	results = search_db('Diamond')
	Diamond = material(name = str(results[0][0]), nk_db = True)()


Adding materials from refractiveindex.info to the database
---------------------------------------------------------------

There is a convenient function, ``create_nk_txt``, to generate the n and k data files needed to
add a new material to the Solcore database directly from the downloaded refractiveindex.info database:

.. code-block:: Python

	results = search_db('Diamond')
	create_nk_txt(pageid=results[0][0], file='C_Diamond')
	create_new_material(mat_name = 'Diamond', n_source='C_Diamond_n.txt', k_source='C_Diamond_k.txt')

This searches the refractiveindex.info database for entries matching 'Diamond', and then creates
files with the n and k data from the first matching database entry in the format required by 
``create_new_material``.
