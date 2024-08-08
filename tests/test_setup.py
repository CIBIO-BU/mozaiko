import os
import unittest


class TestSetup(unittest.TestCase):
    def test_setup(self):
        self.assertTrue(os.path.exists("setup.py"))

    def test_package_metadata(self):
        with open("setup.py") as file:
            metadata = file.read()
        self.assertIn('name="dnaquaimg"', metadata)
        self.assertIn('version="0.1.0"', metadata)
        self.assertIn('url="https://github.com/CIBIO-BU/DNAquaIMG"', metadata)

    def test_install_requires(self):
        with open("setup.py") as file:
            metadata = file.read()
        self.assertIn('install_requires=["numpy==2.0.0", "pandas==2.2.2"]', metadata)

    def test_find_packages(self):
        with open("setup.py") as file:
            metadata = file.read()
        self.assertIn("packages=find_packages()", metadata)
