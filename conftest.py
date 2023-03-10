"""Configuration file for pytest.

Directs the test directory to have access to files as if it were in src.
"""

import pytest
import sys
import os


def pytest_configure(config):
    """Allows plugins and conftest files to perform initial configuration.

    This hook is called for every plugin and initial conftest file after command line
    options have been parsed.
    """

    sys.path.append(os.path.abspath("./src/"))
    print(sys.path)
    return
