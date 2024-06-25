ambiguous_bases = set("RYWSMKHBVDN")

def calculate_ambiguous_percentage(sequence):
    """
    Calculate the percentage of ambiguous bases in a DNA sequence.

    Parameters:
    sequence (str): The DNA sequence.

    Returns:
    float: The percentage of ambiguous bases in the sequence.
    """
    return sum(base in ambiguous_bases for base in sequence) / len(sequence)

def write_filtered_sequence(output_handle, record): #TODO: Check if it makes sense to move this function to the data_transformer.py file
    """
    Write a filtered sequence to the output file.
    """
    sequence = str(record.seq)
    output_handle.write(f">{record.description}\n{sequence}\n")