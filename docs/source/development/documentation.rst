Documentation
=================

Py_gearworks documentation uses Sphinx and the Read the Docs theme.
Building the documentation
--------------------------

Run the following commands from the repository root. A virtual environment is
recommended:

.. code-block:: powershell

	python -m venv venv
	.\venv\Scripts\Activate.ps1
	python -m pip install --upgrade pip
	python -m pip install -e ".[docs]"

Build the HTML documentation from the ``docs`` directory:

.. code-block:: powershell

	cd docs
	.\make.bat html

On Linux and macOS, activate the environment with
``source venv/bin/activate`` and run ``make html`` instead.

The generated site is written to ``docs/build/html/index.html``. Open that
file in a browser to review the result.

Cleaning a build
----------------

Remove generated documentation before a fresh build with:

.. code-block:: powershell

	cd docs
	.\make.bat clean
	.\make.bat html

Generating API documentation
------------------------

This command scans the source code and generates API documentation in the ``docs/source/gendocs`` directory. The generated files are included in the Sphinx build.

.. code-block:: powershell

    cd docs
    python -m sphinx.ext.apidoc -f -e -o source ../src/py_gearworks
