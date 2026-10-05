Python API
==========

Everything is importable from the top-level package (``import physim as ps``). The classes are implemented in
Rust (``physim._core``); the submodules are pure Python.

World and trajectories
----------------------

.. autoclass:: physim.World
   :members:

.. autoclass:: physim.Trajectory
   :members:

.. autoclass:: physim.Event
   :members:

.. autodata:: physim.INTEGRATORS
   :annotation: list of integrator names accepted by World(integrator=...)

Forces
------

.. autoclass:: physim.UniformField
   :members:
.. autoclass:: physim.NewtonianGravity
   :members:
.. autoclass:: physim.TreeGravity
   :members:
.. autoclass:: physim.ParticleMesh
   :members:
.. autoclass:: physim.Spring
   :members:
.. autoclass:: physim.AnchorSpring
   :members:
.. autoclass:: physim.DampedSpring
   :members:
.. autoclass:: physim.ModulatedSpring
   :members:
.. autoclass:: physim.SpringNetwork
   :members:
.. autoclass:: physim.LinearDrag
   :members:
.. autoclass:: physim.QuadraticDrag
   :members:
.. autoclass:: physim.PowerLaw
   :members:
.. autoclass:: physim.Yukawa
   :members:
.. autoclass:: physim.PlummerPotential
   :members:
.. autoclass:: physim.HernquistPotential
   :members:
.. autoclass:: physim.HarmonicTrap
   :members:
.. autoclass:: physim.HenonHeiles
   :members:
.. autoclass:: physim.PeriodicForce
   :members:
.. autoclass:: physim.PostNewtonian
   :members:
.. autoclass:: physim.J2Oblateness
   :members:
.. autoclass:: physim.ElectricField
   :members:
.. autoclass:: physim.MagneticField
   :members:
.. autoclass:: physim.Coulomb
   :members:
.. autoclass:: physim.FieldForce
   :members:
.. autoclass:: physim.LennardJones
   :members:
.. autoclass:: physim.Morse
   :members:
.. autoclass:: physim.TabulatedPair
   :members:
.. autoclass:: physim.SoftContact
   :members:
.. autoclass:: physim.CustomForce
   :members:

Rigid bodies
------------

.. autoclass:: physim.RigidSystem
   :members:
.. autoclass:: physim.RigidTrajectory
   :members:
.. autoclass:: physim.BodyGravity
   :members:
.. autoclass:: physim.BodySpring
   :members:

Fields on grids
---------------

.. autoclass:: physim.WaveEquation
   :members:
.. autoclass:: physim.HeatEquation
   :members:
.. autofunction:: physim.solve_poisson

Quantum mechanics
-----------------

.. autoclass:: physim.Schrodinger
   :members:

.. automodule:: physim.quantum
   :members:

Nuclear experiment planner
--------------------------

The setup file and its Python form; see :doc:`../nuclear-setup` for every field.

.. autoclass:: physim.nuclear.Experiment
   :members:
.. autoclass:: physim.nuclear.Beam
   :members:
.. autoclass:: physim.nuclear.Target
   :members:
.. autoclass:: physim.nuclear.Layer
   :members:
.. autoclass:: physim.nuclear.Detector
   :members:
.. autoclass:: physim.nuclear.Run
   :members:
.. autoclass:: physim.nuclear.SetupError
.. autoclass:: physim.nuclear.Quantity
   :members:
.. autofunction:: physim.nuclear.example_names
.. automodule:: physim.nuclear.data
   :members:

.. automodule:: physim.nuclear.kinematics
   :members:

.. automodule:: physim.nuclear.stopping
   :members:

.. automodule:: physim.nuclear.rutherford
   :members:

.. automodule:: physim.nuclear.detectors
   :members:

.. automodule:: physim.nuclear.coulex
   :members:

.. automodule:: physim.nuclear.gamma
   :members:

.. automodule:: physim.nuclear.levels
   :members:

.. automodule:: physim.nuclear.catalogue
   :members:

.. automodule:: physim.nuclear.scene
   :members:

.. automodule:: physim.nuclear.ensdf
   :members:

.. automodule:: physim.nuclear.rates
   :members:

.. automodule:: physim.nuclear.events
   :members:

.. automodule:: physim.nuclear.planner
   :members:

.. automodule:: physim.nuclear.guide
   :members: Help, Step, help_for, steps_for, reading

.. automodule:: physim.nuclear.report
   :members:

.. automodule:: physim.nuclear.rootio
   :members:

.. automodule:: physim.nuclear.paper
   :members:

.. automodule:: physim.nuclear.app
   :members: figure_geometry, figure_kinematics, figure_strips, figure_energy_loss, figure_spectra,
      figure_trajectories, figure_sweep, themed, report_zip, main

.. automodule:: physim.nuclear.validation
   :members:

.. automodule:: physim.nuclear.export
   :members:

.. automodule:: physim.nuclear.plot
   :members:

.. autofunction:: physim.nuclear.parse_nuclide
.. autofunction:: physim.nuclear.parse_material

Random numbers
--------------

.. automodule:: physim.random
   :members:

Scenarios
---------

.. automodule:: physim.scenarios
   :members:

Units
-----

.. automodule:: physim.units
   :members:

Plotting
--------

.. automodule:: physim.plot
   :members:

Analysis
--------

.. automodule:: physim.analysis
   :members:

Files and checkpoints
---------------------

.. automodule:: physim.io
   :members:

Geometry and inertia helpers
----------------------------

.. automodule:: physim.geometry
   :members:
