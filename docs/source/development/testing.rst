Testing
==================

Tests
-----

Tests are executed using ``pytest`` with extension ``pytest-xdist`` for parallel execution.
The project contains regular pytest tests and doctests as well, embedded in docstrings.

The command to execute tests locally:

.. code-block:: console

   pytest -doctest -n auto

Gear tests of core classes target accurate geometry, checking interference with the help of shapely.

Gear tests of build123d conversion attempt generating build123d gear objects, and check volume and simple things like that.

These 2 types of tests are parameterized to a large number of gear types and features, thus it takes a long time to run them even with parallelism.

Other miscellaneous tests check a few other things.

Github actions: check.yml
-------------------------

Tests are also executed on github actions, the definition is under ``.github/workflows/check.yml``. Execution time can be about 2 hours. (Github actions don't allow or use parallelism).
Local test execution is recommended for development, although you might get 100% cpu load for 20 minutes, might be better than waiting 2 hours for results.

Github actions checking executes on pull requests and pushes to main branch.