Welcome to py_gearworks's documentation! 
============================================================

(version |release|)

Introduction
------------
`Python Gearworks <https://github.com/GarryBGoode/py_gearworks>`_ is a Python module for accurate and feature-rich gear geometry generation.
It is built upon the `build123d <https://github.com/gumyr/build123d>`_ Python package, but can be used with any CAD workflow via STEP export.

Py_gearworks is built on geometric calculations via numpy and scipy, and it is designed to circumvent the common limitations of CAD programs when designing gear surfaces.
It is also designed to be easily extensible via object oriented design, so that new gear types can be added with minimal effort.

Example
^^^^^^^

Create a simple :py:class:`SpurGear <py_gearworks.wrapper.SpurGear>` and export it to a STEP file.

::

    from py_gearworks import *
    from build123d import export_step

    gear1 = SpurGear(number_of_teeth=12, module=1.2, height=5)
    gear_part_1 = gear1.build_part()
    export_step(gear_part_1, "gear1.step")


Features
---------

Py_gearworks is currently focused on geometry generation of gears, and geometric alignment functions.
Gear strength calculations, efficiency calculations are out of scope for this project, but may be added in the future.

+-------------------+-------------------+-------------------------------------+
| Gear Types        | Profile Mods      | Position & Alignment                |
+===================+===================+=====================================+
| Spur              | Undercut          | Std. position                       |
+-------------------+-------------------+-------------------------------------+
| Helical           | Profile shift     | Backlash-controlled position        |
+-------------------+-------------------+-------------------------------------+
| Bevel             | Root/tip fillet   | Axis alignment (bevels and helicals)|
+-------------------+-------------------+-------------------------------------+
| Cycloid           | Crowning          | Rotation to align teeth to mesh     |
+-------------------+-------------------+-------------------------------------+
| Inside-ring       |                   | Axial offset (only spur and helical)|
+-------------------+-------------------+-------------------------------------+

Not implemented
^^^^^^^^^^^^^^^^

- Worm gears
- Hypoid gears
- Face / Crown gears


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   examples_lander
   development
   api
   



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
