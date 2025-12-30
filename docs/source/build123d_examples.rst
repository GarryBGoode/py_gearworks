Build123d Workflow Examples
===========================

Simple gears on a backplate
----------------------------

Various properties and methods are made available in the class :py:class:`GearInfoMixin <py_gearworks.wrapper.GearInfoMixin>`. 
The following example demonstrates the creation of a gear-pair and attaching them to a base-plate: 

.. image:: ./assets/gears_on_plate.png
  :align: center
  :alt: token gears_on_plate

Highlights:

* You can use \
  :py:attr:`center_location_bottom <py_gearworks.wrapper.GearInfoMixin.center_location_bottom>`, \
  :py:attr:`center_location_top <py_gearworks.wrapper.GearInfoMixin.center_location_top>`, \
  :py:attr:`face_location_bottom <py_gearworks.wrapper.GearInfoMixin.face_location_bottom>`, \
  :py:attr:`face_location_bottom <py_gearworks.wrapper.GearInfoMixin.face_location_bottom>` to align parts with gear centers.
* Note that **center** refers to the pitch circle center and **face** refers to the face (surface) of the gear. These are different for bevel gears.
* :py:attr:`gear.center <py_gearworks.wrapper.GearInfoMixin.center>` is a numpy array, often needs to be converted to a `Vector` for build123d via :py:func:`np2v() <py_gearworks.conv_build123d.np2v>`.
* Ideal center distance can be retrieved after calling the :py:meth:`mesh_to() <py_gearworks.wrapper.InvoluteGear.mesh_to>` method, and calculating the difference of :py:attr:`gear.center <py_gearworks.wrapper.GearInfoMixin.center>` values.


.. literalinclude:: ../../examples/build123d_example.py



Crescent Gear Pump
-------------------

This example demonstrates building a gear pump. The design is missing fasteners and seals,
but showcases the gear generator and its helper functions for build-123d workflow.

Highlights:

* You can use `center_location_bottom` and `center_location_top` to align parts with gear centers.
* The :py:attr:`radii_data_top <py_gearworks.wrapper.GearInfoMixin.radii_data_top>` method generates reference curves for the gear.
* The :py:class:`LineOfAction <py_gearworks.wrapper.GearInfoMixin.LineOfAction>` class is available for generating the line of action between gears.
* Sometimes converter functions are needed such as :py:func:`arc_to_b123d() <py_gearworks.conv_build123d.arc_to_b123d>` and :py:func:`line_to_b123d() <py_gearworks.conv_build123d.line_to_b123d>`. These convert between py_gearworks' own geometry classes and build123d geometry.
* Animation is used from `ocp_vscode` to visualize the gear meshing. Animation can be non-intuitive, but explaining it is beyond this example.


.. image:: ./assets/gearpump_1.png
  :align: center
  :alt: token gearpump_1

.. literalinclude:: ../../examples/build123d_crescent_gearpump.py

Here you can find the sketches that helped the construction of the crescent and the fluid channels.

.. image:: ./assets/gearpump_2.png
  :align: center
  :alt: token gearpump_with_circles
