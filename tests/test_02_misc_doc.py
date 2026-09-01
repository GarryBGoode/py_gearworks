# Copyright 2026 Gergely Bencsik
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import py_gearworks.wrapper as pgw
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import doctest
import sys
import os
import glob
import runpy
import inspect

import pytest


def test_example():
    """Test all the examples in the examples.py file. Can't make a passing assertion,
    except that there should not be any errors when running the examples."""
    current_dir = os.path.dirname(__file__)
    relative_path = os.path.join(current_dir, "..", "examples")
    sys.path.append(relative_path)

    import examples

    # Get all functions from the examples module
    example_functions = [
        func
        for name, func in inspect.getmembers(examples, inspect.isfunction)
        if func.__module__ == "examples"
    ]

    # Call each example function
    for func in example_functions:
        print(f"Running {func.__name__}...")
        func()

def _example_scripts():
    """Return all standalone example scripts, excluding examples.py which is
    exercised by test_example above."""
    examples_dir = os.path.join(os.path.dirname(__file__), "..", "examples")
    scripts = sorted(glob.glob(os.path.join(examples_dir, "*.py")))
    return [s for s in scripts if os.path.basename(s) != "examples.py"]


@pytest.mark.parametrize(
    "script", _example_scripts(), ids=lambda s: os.path.basename(s)
)
def test_example_2(script):
    """Every script under examples/ (except examples.py) must run without error."""
    cwd = os.getcwd()
    try:
        os.chdir(os.path.dirname(script))
        runpy.run_path(script, run_name="__main__")
    finally:
        os.chdir(cwd)
        plt.close("all")



def test_doctest():
    """
    Run doctests on the py_gearworks module.
    """
    # Run doctests and check for failures
    doctest_results = doctest.testmod(m=pgw)
    assert doctest_results.failed == 0


if __name__ == "__main__":
    test_example()
    test_doctest()
