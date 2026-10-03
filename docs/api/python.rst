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
