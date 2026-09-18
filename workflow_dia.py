# Imports
import time
from datetime import datetime
from src.mozaiko.reference_database.sequence_import import *
from src.mozaiko.in_silico_analysis.amplification import InSilicoAmplification
from src.mozaiko.marker_scoring.metrics_system import *

# Input Data
database_file = "data/input_data/DIA/diat_barcode_hrm_dreped.fasta"
primer_table = "data/input_data/DIA/diat-barcode-primers.tsv"
run_name = 'DIA-test-7Aug'
otl_folder = "/home/camilababo/Documents/DNAquaIMG/countries-otls/harmonized/dia"
output_folder = 'data/output_data/' + run_name

# Logging setup
log_dir = Path("data/output_data") / run_name / "logs"
log_dir.mkdir(parents=True, exist_ok=True)

log_file = log_dir / f"{run_name}_workflow.log"


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


log_handle = open(log_file, "w", buffering=1)

# Everything printed to stdout/stderr goes to both terminal + log
sys.stdout = Tee(sys.__stdout__, log_handle)
sys.stderr = Tee(sys.__stderr__, log_handle)


# Print starting date/time
start_info = datetime.now()
start_time = time.time()
print(f"mozaiko INFO: Starting mozaiko workflow for {run_name} at {start_info.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"mozaiko INFO: Log file: {log_file}")

# Data Import
custom_fasta_import = CustomFastaImport(database_file)
custom_fasta_import.read_fasta(database_file, check_taxid=False)
custom_fasta_import.pre_process_harmonized_fasta_database()
data = custom_fasta_import.data

# In Sillico Analysis
insil = InSilicoAmplification(custom_fasta_import.database_fasta_file, run_name=run_name)
insil.run_in_silico_analysis(primer_table,  minimum_percentage_identity=0.50)

# Primer Evaluation
MetricsSystemExecutor.evaluate_several_OTLs(otl_folder=otl_folder,
                      output_folder=output_folder,
                      primer_table=primer_table,
                      save_intermediate_ranks=True,
                      run_catnip=True,
                      thresholds = [10.0, 5.0, 2.0],
                      ranking_mode='flat',
                      min_barcode=5)

# Log time
end_time = time.time()
elapsed_time = end_time - start_time
print(f"mozaiko INFO: Total execution time: {elapsed_time/60:.2f} minutes.")
