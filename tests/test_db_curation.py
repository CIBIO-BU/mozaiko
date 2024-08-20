import subprocess
import unittest
from unittest.mock import mock_open, patch

from src.reference_database.db_curation import CrabsScriptGenerator


class TestCrabsScriptGenerator(unittest.TestCase):

    def setUp(self):
        self.generator = CrabsScriptGenerator()

    @patch("subprocess.run")
    def test_check_if_crabs_installed(self, mock_subprocess):
        self.generator._check_if_crabs_installed()
        mock_subprocess.assert_called_with(
            ["crabs", "-h"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

    patch("subprocess.run", side_effect=FileNotFoundError)

    def test_check_if_crabs_installed_not_installed(self, mock_run):
        with self.assertRaises(SystemExit):
            self.generator._check_if_crabs_installed()

        mock_run.assert_called_once_with(
            ["crabs", "-h"], check=True, stdout=-3, stderr=-3
        )

    @patch(
        "builtins.open", new_callable=mock_open, read_data='{"input": "input.fasta"}'
    )
    @patch("json.load")
    def test_load_parameters(self, mock_json_load, mock_file):
        mock_json_load.return_value = {"input": "input.fasta"}
        self.generator._load_parameters("dummy.json")
        mock_file.assert_called_with("dummy.json", encoding="UTF-8")
        mock_json_load.assert_called_once()

    @patch("os.path.exists", return_value=True)
    def test_download_taxonomy_files_exists(self, mock_path_exists):
        self.generator._download_taxonomy_files()
        mock_path_exists.assert_called_with("taxonomy_files")

    @patch("os.path.exists", return_value=False)
    @patch("os.makedirs")
    @patch("os.chdir")
    @patch("subprocess.run")
    def test_download_taxonomy_files_not_exists(
        self, mock_run, mock_chdir, mock_mkdir, mock_path_exists
    ):
        self.generator._download_taxonomy_files()
        mock_mkdir.assert_called_with("taxonomy_files")
        mock_chdir.assert_called_with("taxonomy_files")
        mock_run.assert_called_with(
            "crabs db_download --source taxonomy", shell=True, check=True
        )

    @patch("subprocess.run")
    @patch(
        "src.reference_database.db_curation.CrabsScriptGenerator._download_taxonomy_files"
    )
    @patch("src.reference_database.db_curation.CrabsScriptGenerator._load_parameters")
    def test_run_assign_tax_command(
        self, mock_load_params, mock_download_files, mock_run
    ):
        self.generator.params = {
            "input": "input_file",
            "output": "output_file",
            "acc2tax": "acc2tax_file",
            "taxid": "taxid_file",
            "name": "name_file",
            "web": "web",
            "ranks": "rank",
            "missing": "missing",
        }

        self.generator.run_assign_tax_command("dummy.json")
        mock_load_params.assert_called_with("dummy.json")
        mock_download_files.assert_called_once()
        mock_run.assert_called_with(
            "crabs assign_tax --input input_file --output output_file --acc2tax acc2tax_file \
                  --taxid taxid_file --name name_file --web web --rank rank --missing missing",
            shell=True,
            check=True,
        )
