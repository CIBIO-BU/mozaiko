"""
Unit tests for the setup.py file.
"""

import os
import unittest

from setuptools import find_packages


class TestSetup(unittest.TestCase):
    """
    Class to test the setup.py file.
    """

    def test_setup(self):
        """
        Test if the setup.py file exists.
        """
        self.assertTrue(os.path.exists("setup.py"))

    def test_package_metadata(self):
        """
        Test if the package metadata is correctly defined.
        """
        with open("setup.py", encoding="UTF-8") as file:
            metadata = file.read()
        self.assertIn('name="dnaquaimg"', metadata)
        self.assertIn('version="0.1.0"', metadata)
        self.assertIn('url="https://github.com/CIBIO-BU/DNAquaIMG"', metadata)

    def test_find_packages_presence(self):
        """
        Test if the find_packages function is called in the setup.py file.
        """
        with open("setup.py", encoding="UTF-8") as file:
            metadata = file.read()
        self.assertIn("packages=find_packages()", metadata)

    def test_find_packages_call(self):
        """
        Test the find_packages function call mathces the actual packages in the repository.
        """
        packages = find_packages()
        expected_packages = ["src", "tests"]
        self.assertTrue(all(pkg in packages for pkg in expected_packages))
