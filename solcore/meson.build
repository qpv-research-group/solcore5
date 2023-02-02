install_subdir('absorption_calculator', install_dir: py.get_install_dir() / 'solcore')
install_subdir('analytic_solar_cells', install_dir: py.get_install_dir() / 'solcore')
install_subdir('data_analysis_tools', install_dir: py.get_install_dir() / 'solcore')
install_subdir('graphing', install_dir: py.get_install_dir() / 'solcore')
install_subdir('light_source', install_dir: py.get_install_dir() / 'solcore')
install_subdir('material_data', install_dir: py.get_install_dir() / 'solcore')
install_subdir('material_system', install_dir: py.get_install_dir() / 'solcore')
install_subdir('optics', install_dir: py.get_install_dir() / 'solcore')
install_subdir('optimization', install_dir: py.get_install_dir() / 'solcore')
install_subdir('parameter_system', install_dir: py.get_install_dir() / 'solcore')
subdir('poisson_drift_diffusion')
install_subdir('quantum_mechanics', install_dir: py.get_install_dir() / 'solcore')
install_subdir('spice', install_dir: py.get_install_dir() / 'solcore')
install_subdir('units_system', install_dir: py.get_install_dir() / 'solcore')


python_sources = [
  '__init__.py',
  'solcore_config.txt',
  'config_tools.py',
  'constants.py',
  'crystals.py',
  'interpolate.py',
  'registries.py',
  'science_tracker.py',
  'singleton.py',
  'smooth.py',
  'solar_cell.py',
  'solar_cell_solver.py',
  'source_managed_class.py',
  'state.py',
  'strain_balancing.py',
  'structure.py'
]

py.install_sources(
  python_sources,
  subdir: 'solcore'
)