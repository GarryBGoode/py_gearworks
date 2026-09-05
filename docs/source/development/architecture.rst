Architecture
============

Files
-----------------

Py_gearworks file structure was originally organized into layers as follows:

* Wrapper
    * This is meant for user interaction
* build123d conversion
    * This converts from custom / internal py gearworks objects to build123d objects
* Core layer
    * core.py, base_classes.py, gearteeth.py
    * This contains the gear-creator engine, key classes, basic gear tooth geometry calculations
* Curve
    * Classes to describe and manipulate curves
* Function defs
    * Just math

Python Modules, dependencies
----------------------------

Py_gearworks was developed with the goal of minimizing dependencies.

The main dependencies are:

* build123d - for CAD geometry generation and export. Ocp_vscode kind of comes with build123d.
* numpy: math
* scipy: more math

For testing, pytest, matplotlib, shapely are used. For documentation, sphinx and sphinx_rtd_theme are used.
