"""
This module contains the methods needed to perform the in-silico amplification analysis.
"""

import subprocess
import sys


class InSilicoAmplification:
    """
    This class contains the methods needed to perform the in-silico amplification analysis.
    """

    def __init__(self, data):
        self.data = data

    def _check_if_cutadapt_installed(self):
        """
        Function to check if Cutadapt is installed.
        """
        print("mosaiko INFO: Checking if cutadapt is installed...")
        try:
            subprocess.run(["cutadapt", "--version"], check=True)

            print("mosaiko INFO: Cutadapt is installed.")

        except FileNotFoundError:
            print(
                "mosaiko INFO: Cutadapt is not installed. Please install Cutadapt before running this \
                    script."
            )
            print(
                "Cutadapt can be found at "
                + "https://cutadapt.readthedocs.io/en/stable/installation.html"
            )
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)

    def _filter_sequences_by_prcnt_ambiguous_bases():
        pass
