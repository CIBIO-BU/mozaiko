# Log time
import time
start_time = time.time()

# Imports
from src.reference_database.sequence_import import *
from src.in_silico_analysis.amplification import InSilicoAmplification
from src.marker_scoring.metrics_system import *

# Files
database_file = "data/input_data/BMI/bmi_database.fasta"
primer_table = "data/input_data/BMI/bmi-primers.tsv"
run_name = 'BMI-NEW-CATNIP'
otl_folder = "/home/camilababo/Documents/DNAquaIMG/countries-otls/harmonized/bmi/filtered_phyl"
output_folder = 'data/output_data/' + run_name

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
                      thresholds = [10.0, 5.0, 2.0])

# Log time
end_time = time.time()
elapsed_time = end_time - start_time
print(f"mozaiko INFO: Total execution time: {elapsed_time/60:.2f} minutes.")