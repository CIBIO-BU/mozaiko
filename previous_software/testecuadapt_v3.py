#%%
import os
import subprocess as sp
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from typing import DefaultDict
from Bio.SeqUtils import MeltingTemp
from Bio.SeqUtils import gc_fraction
import re

#### Functions
def filter_sequences(infile, out_file, max_ambiguous_percentage=0.05):
    """
    Filters sequences based on the maximum allowed percentage of ambiguous bases.

    Parameters:
    infile (str): Path to the input file containing DNA sequences in FASTA format.
    out_file (str): Path to the output file to write the filtered sequences.
    max_ambiguous_percentage (float): Maximum allowed percentage of ambiguous bases in a sequence.

    Returns:
    None
    """
    ambiguous_bases = set("RYWSMKHBVDN")
    with open(out_file, 'w') as output_handle:
        for record in SeqIO.parse(infile, 'fasta'):
            sequence = str(record.seq)
            ambiguous_percentage = sum(base in ambiguous_bases for base in sequence) / len(sequence)
            if ambiguous_percentage <= max_ambiguous_percentage:
                # Write the header and sequence on separate lines
                output_handle.write(f">{record.description}\n{sequence}\n")

def filter_sequences_getprimers(infile, out_file, FWDlen, REVlen, max_ambiguous_percentage=0.05):
    """
    Filters sequences based on the maximum allowed percentage of ambiguous bases and extracts primers.

    Parameters:
    infile (str): Path to the input file containing DNA sequences in FASTA format.
    out_file (str): Path to the output file to write the filtered sequences.
    FWDlen (int): Length of the forward primer.
    REVlen (int): Length of the reverse primer.
    max_ambiguous_percentage (float): Maximum allowed percentage of ambiguous bases in a sequence.

    Returns:
    DefaultDict: Dictionary of extracted primers.
    
    """
    ambiguous_bases = set("RYWSMKHBVDN")
    primers = DefaultDict(lambda:[None, None, None]) # Initialize the primers dictionary
    with open(out_file, 'w') as output_handle:
        for record in SeqIO.parse(infile, 'fasta'):
            sequence = str(record.seq)
            ambiguous_percentage = sum(base in ambiguous_bases for base in sequence) / len(sequence)
            if ambiguous_percentage <= max_ambiguous_percentage:
                ac = record.id # Accession number
                primers[ac][0] = sequence[0:FWDlen] # Forward primer
                primers[ac][1] = sequence[len(sequence)-REVlen:] # Reverse primer
                primers[ac][2] = len(sequence) - REVlen - FWDlen # Amplicon length
                # Write the header and sequence on separate lines
                output_handle.write(f">{record.description}\n{sequence}\n")
    return primers


def fastatodataframe(inFile):
    """
    Transform a fasta file into a pandas DataFrame.

    Parameters:
    inFile (str): Path to the input fasta file.

    Returns:
    pd.DataFrame: DataFrame containing the sequences.
    
    """
    seq_object = SeqIO.parse(inFile, "fasta")
    sequences = [] # This will store the entries (header and sequence)
    for seq in seq_object: # Store the entries hits in a list
        sequences.append(seq)
    samples = [] 
    tax = []
    seq_lenghts = []
    sequences2 = [] # 
    for record in sequences: # Iterate over the entries
        seq = str(record.seq) # Get the sequence
        length = len(seq)
        sample = record.id # Accession number
        taxid = re.search(r'(?<=taxid=)([0-9]+)(?=;)',record.description).group(1) # Finds a sequence of digits between 'taxid=' and ';'
        samples.append(sample)
        tax.append(taxid)
        seq_lenghts.append(length)
        sequences2.append(seq)
    data =  pd.DataFrame(list(zip(samples, tax, seq_lenghts, sequences2)), columns=["AC", "TaxID", "Length", "Sequence"]) # Create a DataFrame from a list of tuples
    # Each tuple contains the accession number, taxid, sequence length, and sequence
    # Elements are then assigned to columns
    return(data)


def defineresolution(df):
    """
    Defines the resolution rank of the sequences in the DataFrame.

    Parameters:
    df (pd.DataFrame): DataFrame containing the sequences.

    Algorithm:
    1. Iterate over the unique sequences in the DataFrame.
    2. Create a temporary DataFrame with the attributes of the sequence.
    3. Count the number of unique values in each column.
    4. Identify the resolution rank.
    5. Store the resolution rank in a dictionary.

    * Taxonomic resolution is the level to which a taxon is identified.

    Returns:
    DefaultDict: Dictionary containing the sequences and their resolution rank.

    """
    seqsRanks = DefaultDict(lambda:[None]) # Initialize the dictionary of sequences and their resolution rank
    sequences = set(df['Sequence'].tolist()) 
    for s in sequences:
        temp=df.loc[df['Sequence']==s, ] # Get the rows of the DataFrame that match the sequence -> Creates a temporary DataFrame with the sequence attributes
        a = temp[['Species', 'Genus', 'Family', 'Order', 'Class', 'phylum', 'Kigdom']].nunique().to_frame() # Count the number of unique values in each column
        a = a.transpose()
        rank = (a == 1).idxmax(axis=1).iloc[0] # Defines the resolution rank as the highest taxonomic rank with a single unique value
        # e.g. If there are 3 unique values in the 'Species' column, 1 in the 'Genus', and 1 in the 'Phylum': the resolution rank is 'Genus'
        # idxmax(axis=1) returns the index of first occurrence
        seqsRanks[s][0]=rank # Why the use [0]?
    return(seqsRanks)

def frominsilico(fasta): # Function name is not descriptive of its actions
    """
    Transform a fasta file into a list of sequence identifiers.
    
    Parameters:
    fasta (str): Path to the input fasta file.

    Returns:
    list: List of sequence identifiers.
    """
    Ids = []
    for record in SeqIO.parse(fasta, 'fasta'):
        Ids.append(record.id)
    return(Ids)


def runcut(fasta, dataframe, e, finalname):
    """
    Run cutadapt on the input fasta file.
    
    * cutadapt is a tool that finds and removes adapter sequences, primers,
    poly-A tails and other types of unwanted sequence from your high-throughput 
    sequencing reads.


    Parameters:
    fasta (str): Path to the input fasta file.
    dataframe (pd.DataFrame): DataFrame containing the primer sequences.
    e (int): Maximum allowed error rate.
    finalname (str): Suffix to add to the output file name.

    Returns:
    None
    
    """
    for index, row in dataframe.iterrows(): 
        AssayName = row['AssayName']
        print(AssayName)
        FwSeq = row['FwSeq'].replace('I', 'N') # Replace 'I' with 'N' in the forward primer sequence -> Why?
        RvSeq = row['RvSeq'].replace('I', 'N')
        br = row['BarcodeRegion']
        ERROR = int(e)
        RvSeq_CORRECT = str(Seq(RvSeq).reverse_complement()) # Reverse complement of the reverse primer sequence
        OVERLAP = str(min([len(FwSeq), len(RvSeq_CORRECT)])) # Find the minimum length between the forward and reverse primer sequences
        # Region of overlap between the forward and reverse primer sequences shouldn't exceed the length of the shortest primer
        # so that all every position of the longer primer is covered by the shorter primer
        # Avoids mismatches and gaps
        ADAPTER = FwSeq + '...' + RvSeq_CORRECT # Concatenate the forward and reverse primer sequences with '...'
        OUTNAME = './cutadapt/'+ br +'/'+ AssayName+finalname # Output file name based on the barcode region and assay name
        action = 'retain'
        # "With --action=retain, the read is trimmed, but the adapter sequence itself is not removed. Up- or downstream sequences 
        # are removed in the same way as for the trim action. For linked adapters, both adapter sequences are kept."
        #a = sp.Popen(['cutadapt', '-g', ADAPTER, '-o', OUTNAME, fasta, '--no-indels', '-e', str(ERROR), '--overlap', OVERLAP, '--action', 
                      #str(action),'--discard-untrimmed', '--revcom', '-M', '600', '--quiet'])
        #a.wait()

        # Improved implementation -> separate the command:

        command = [
            'cutadapt', '-g', ADAPTER, '-o', OUTNAME, fasta, '--no-indels',
            '-e', str(ERROR), '--overlap', OVERLAP, '--action', str(action),
            '--discard-untrimmed', '--revcom', '-M', '600', '--quiet'
        ] 
        process = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)

        process.wait()


def checkcutadapt(cutoutput, iberiaTaxonomy):
    """
    This functions creates a summary for the cutadapt output file.
    
    ALgorithm:
    1. Read and tranform the .cutadapt file.
    2. Merge with taxonomic infomation.
    3. Calculate taxonomic resolution.
    4. Merges resolution with cutadapt and taxonomic information.
    5. Creates a summary file with the results.

    Parameters:
    cutoutput (str): Path to the cutadapt output file.
    iberiaTaxonomy (pd.DataFrame): DataFrame containing the taxonomy information.

    Returns:
    None
    """
    cutoutput = str(cutoutput)
    inserts = fastatodataframe(cutoutput)

    inserts['AC'] = inserts['AC'].apply(lambda x: x.split(';')[0]) 
    inserts['TaxID'] = inserts['TaxID'].astype(int)
    inserts['Sequence'] = inserts['Sequence'].apply(lambda x: x.upper()) 
    inserts = pd.merge(inserts, iberiaTaxonomy[['Kigdom', 'phylum', 'Class', 'Order', 'Family', 'Genus','Species', 'TaxID']], how= 'outer', on='TaxID')
    
    # Calculates taxonomic resolution (to which level a taxon is identified)
    ranks = defineresolution(inserts)
    ranks = pd.DataFrame([{'Sequence': key, 'Resolution': value[0]} for key, value in ranks.items()])
    ranks = ranks.loc[~ranks['Sequence'].isna()] # Remove sequences with no information
    inserts = pd.merge(inserts, ranks, on='Sequence', how='left')
    
    nan_species_df = inserts.loc[inserts['AC'].isna()]
    grouped_df = inserts.groupby(['Species', 'Resolution']).agg({'Sequence': 'count', 'Length': 'mean'}).reset_index()
    all = pd.merge(grouped_df, nan_species_df['Species'], on = 'Species', how= 'outer')
    all.to_csv(cutoutput.replace('.cutadapt', '.summary'), sep = '\t', index = None)
    
    # Count species with a non-zero sequence
    #species_with_sequence = all[all['Sequence'] != 0]['Species'].nunique()

    # Calculate the sum of sequences for species with a non-zero sequence
    #total_sequence_sum = all[all['Sequence'] != 0]['Sequence'].sum()

    # Find the minimum and maximum length
    #min_length = all['Length'].min()
    #max_length = all['Length'].max()

    # Count species with no information (NaN) in the Sequence column
    #species_no_sequence_info = all['Sequence'].isna().sum()

    # Create a new DataFrame to store the results
    #result_df = pd.DataFrame({
    #    'Metric': ['Species Input', 'TotalSpeciesAmplified', 'MinLength', 'MaxLength', 'SpeciesNOTAmplified'],
    #    'Value': [species_with_sequence, total_sequence_sum, min_length, max_length, species_no_sequence_info]
    #})
    #result_df.to_csv(cutoutput.replace('.cutadapt', '.Bc'), sep = '\t', index = None)

def primersettings(primers, FWD, REV_CORRECT, outname):
    """
    Processes initial data for the primer sequences.

    Parameters:
    primers (DefaultDict): Dictionary containing the extracted primers.
    FWD (str): Forward primer sequence.
    REV_CORRECT (str): Reverse primer sequence.
    outname (str): Output file name.

    Returns:
    None
    
    """
    ## calculate metling temperature primers through GC content w/ BioPython
    primers = pd.DataFrame([{'seqid': key, 'PrimerF': value[0], 'PrimerR': value[1], 'InsertLength': value[2]} for key, value in primers.items()])
    primers['PrimerF_Mt'] = round(primers['PrimerF'].apply(MeltingTemp.Tm_GC, strict=False),2)
    primers['PrimerR_Mt'] = round(primers['PrimerR'].apply(MeltingTemp.Tm_GC, strict=False),2)
    
    # check mismatches
    # pairs nucleotides of primer seq. x (c1) with FWD and REV_CORRECT (c2) to count mismatches
    primers['PrimerF_Mis'] = primers['PrimerF'].apply(lambda x: sum(c1 != c2 for c1, c2 in zip(x, FWD))) 
    primers['PrimerR_Mis'] = primers['PrimerR'].apply(lambda x: sum(c1 != c2 for c1, c2 in zip(x, REV_CORRECT)))
    
    ### last 3 bps of each primers
    ### Important in evaluating binding efficency
    primers['PrimerF_3'] = primers['PrimerF'].apply(lambda x: x[len(x)-3:])
    primers['PrimerR_3'] = primers['PrimerR'].apply(lambda x: x[:3])
    
    ### mismatch 3 bps of each primers (Fwd is last 3 bps and Rvs )
    FWD3 = FWD[len(FWD)-3:]
    REV_CORRECT3 = REV_CORRECT[:3]
    primers['PrimerF_3Mis'] = primers['PrimerF_3'].apply(lambda x: sum(c1 != c2 for c1, c2 in zip(x, FWD3)))
    primers['PrimerR_3Mis'] = primers['PrimerR_3'].apply(lambda x: sum(c1 != c2 for c1, c2 in zip(x, REV_CORRECT3)))
    
    ### Calculates GC content per primer
    primers['GC_F'] = primers['PrimerF'].apply(lambda x: round(gc_fraction(x)*100, 2))
    primers['GC_R'] = primers['PrimerR'].apply(lambda x: round(gc_fraction(x)*100, 2))

    primers.to_csv(outname, sep='\t', index=None)
    print('Primers initial data is processed!')

def runcutnoprimers(dataframe):
    """
    Run cutadapt to trim primer sequences and discard them.

    *few differences compared to runcut function.
    
    Parameters:
    dataframe (pd.DataFrame): DataFrame containing the primer sequences.

    Returns:
    None
    """
    for index, row in dataframe.iterrows():
        AssayName = row['AssayName']
        print(AssayName)
        FwSeq = row['FwSeq'].replace('I', 'N')
        RvSeq = row['RvSeq'].replace('I', 'N')
        br = row['BarcodeRegion']
        ERROR = 3
        RvSeq_CORRECT = str(Seq(RvSeq).reverse_complement())
        OVERLAP = str(min([len(FwSeq), len(RvSeq_CORRECT)]))
        ADAPTER = FwSeq + '...' + RvSeq_CORRECT
        fasta = './cutadapt/'+ br +'/'+ AssayName+'.clean.cutadapt'
        OUTNAME = './cutadapt/'+ br +'/'+ AssayName+'.noPrimer.cutadapt'
        if os.path.isfile(fasta):
            
            command = ['cutadapt', '-g', ADAPTER, '-o', OUTNAME, fasta, '--no-indels', '-e', 
                        str(ERROR), '--overlap', OVERLAP,'--discard-untrimmed',
                          '--revcom', '--quiet']
            
            # No action is mentioned, meaning it goes to deafault -> 'trim'
            
            process = sp.Popen(command)
            
            process.wait()


def count_seqsTaxId(fasta):
    """
    Count the number of sequences by TaxID in the input fasta file.

    *The goal of this function is to keep tabs on the number of sequences that get
    processed and removed at each step of the pipeline (AKA taxonomic distribution). 
    At the end, it allows us to compute statistics on primer efficiency 
    and sequence recovery rates across different taxonomic groups.

    Parameters:
    fasta (str): Path to the input fasta file.

    Returns:
    pd.DataFrame: DataFrame containing the TaxID and the number of sequences.
    """
    taxa = DefaultDict(lambda:0) # Initialize an dictionary of taxa and their counts to 0
    
    for record in SeqIO.parse(fasta, 'fasta'):
        # Retrieves a specific taxa ID and counts the number of sequences found 
        species = re.search(r'(?<=taxid=)([0-9]+)(?=;)', record.description).group(1)
        taxa[species]+=1
    
    taxa = pd.DataFrame([{'TaxID': key, 'Nseqs': value} for key, value in taxa.items()])
    
    return(taxa)

#%%
TAXONOMY = 'taxonomy.txt'
FASTA = 'reference.fasta'

PRIMERLIST = pd.read_csv('previous_software/AllPrimers_gene.txt', sep='\t', usecols=['BarcodeRegion','AssayName', 'FwSeq', 'RvSeq', 'AvSize','minL','maxL'])
BR = set(PRIMERLIST['BarcodeRegion'].tolist())

# -------------------------------------------------------------------------------
print("############# Running Step 1: Checking barcode regions and running cutadapt #############")

for B in BR: #for barcode in barcode region?
    print(B) # Print the barcode  
    
    if not os.path.exists('./cutadapt/'+B):
        os.makedirs('./cutadapt/'+B)
    
    temp = PRIMERLIST.loc[PRIMERLIST['BarcodeRegion']==B] # Subsets primer list by barcode region so that
    # only the primers for the current barcode region are considered to runcut
    #runcut(FASTA, temp, B)
    runcut(FASTA, temp, 3, '.cutadapt')

print('############# Step 1 Done ##############')

# -------------------------------------------------------------------------------
print("############# Running Step 2: Removing empty files #############")
for B in BR:
    cutadapt_outs = [x for x in os.listdir('./cutadapt/'+B) if os.path.getsize('./cutadapt/'+B+'/'+x)==0]  
    for f in cutadapt_outs:
        os.remove('./cutadapt/'+B+'/'+f)

print('############# Step 2 Done ##############')

# -------------------------------------------------------------------------------
print('############# Step 3: Filtering sequences and retrieve primer data #############')

for B in BR:
    # print(B)
    temp = PRIMERLIST.loc[PRIMERLIST['BarcodeRegion']==B, ['AssayName', 'FwSeq', 'RvSeq']] # Subset and select cols. 'AssayName', 'FwSeq', and 'RvSeq'
    for index, row in temp.iterrows():
        t = str(row['AssayName']) # Get the assay name
        FWD = row['FwSeq']
        FWDlen = len(FWD)
        REV = row['RvSeq']
        REVlen = len(REV)
        REV_CORRECT = str(Seq(REV).reverse_complement())
        if os.path.isfile('./cutadapt/'+B+'/'+ t+'.cutadapt'): # check if the file exists for the current barcode region and assay
            # filter sequences by ambiguous bases and extract primers
            primers = filter_sequences_getprimers('./cutadapt/'+B+'/'+ t+'.cutadapt', './cutadapt/'+B+'/'+ t+'.clean.cutadapt', FWDlen, REVlen)
            if primers:
                #primersettings(primers, FWD, REV_CORRECT, './cutadapt/'+B+'/'+ t+'.settings')
                primersettings(primers, FWD, REV_CORRECT, './cutadapt/'+B+'/'+ t+'.primerstatistics')

print('############# Step 3 Done ##############')
# -------------------------------------------------------------------------------


print('\n ############# Step 4: Removing primer sequences from in-silico data ############# \n')

for B in BR:
    print(B)
    temp = PRIMERLIST.loc[PRIMERLIST['BarcodeRegion']==B]
    runcutnoprimers(temp, B)

print('############# Step 4 Done ##############')
# -------------------------------------------------------------------------------

print('\n ############# Step 5: Running CRABS pga and cleaning it ############# \n')

##removing empty files
for B in BR:
    cutadaptouts = [x for x in os.listdir('./cutadapt/'+B) if os.path.getsize('./cutadapt/'+B+'/'+x)==0]  
    for f in cutadaptouts:
        os.remove('./cutadapt/'+B+'/'+f)

# run pga

""" *Info on PGA: pairwise global alignment. 
Amplicon regions retrieved through cutadapt are used as seed sequences. (--database)
Includes alignments that have the lenght of the foward or reserve primer-binding region (--strict).
    This is to minimize the risk of including erroneous sequences.
"""


for B in BR:
    print(" Running PGA for barcode region: ", B)

    temp = PRIMERLIST.loc[PRIMERLIST['BarcodeRegion']==B, ['AssayName', 'FwSeq', 'RvSeq']] # susbset primer list by barcode region
    for index, row in temp.iterrows():
        t = str(row['AssayName'])
        print(f'Checking {t}')
        FWD = row['FwSeq'].replace('I', 'N')
        REV = row['RvSeq'].replace('I','N')
        REV_CORRECT = str(Seq(REV).reverse_complement())
        if os.path.isfile('./cutadapt/'+B+'/'+ t+'.clean.cutadapt'):
            try:
                command = ['crabs', 'pga', '--input', FASTA, '--output', './cutadapt/'+B+'/'+ t+'.pga', 
                           '--database', './cutadapt/'+B+'/'+ t +'.clean.cutadapt', '--fwd', FWD,'--rev', 
                           REV_CORRECT, '--speed', 'slow', '--percid', str('0.95'), 
                           '--coverage', str('0.99'), '--filter_method', 'strict']

                process = sp.Popen(command)

                process.wait()
                
                # filter sequences based on ambiguous bases
                filter_sequences('./cutadapt/'+B+'/'+ t+'.pga', './cutadapt/'+B+'/'+ t+'.clean.pga')
            
            except:
                print(f'Primer pair{t} from {B} has found inserts in for all sequences in initial fasta')
                # Doesn't perform PGA for sequences that have inserts in all sequences
                # Inserts can create mismatches and gaps in the alignment due to length differences

print('############# Step 5 Done ##############')
# -------------------------------------------------------------------------------

print('\n ############# Step 6: Running cutadapt with max. error of 5 ############# \n')


##### runc cutadapt with error 5; add here previuosly to noPrimer (add? part of an update?)
for B in BR:
    # print(B)
    temp = PRIMERLIST.loc[PRIMERLIST['BarcodeRegion']==B]
    #runcut(FASTA, temp, B)
    runcut(FASTA, temp, 5, '.maxerror')
    print(f"Running cutadapt with max error = {5}")


#### remove empty files:
for B in BR:
    cutadaptouts = [x for x in os.listdir('./cutadapt/'+B) if x.endswith('.maxerror') and os.path.getsize('./cutadapt/'+B+'/'+x)==0]  
    for f in cutadaptouts:
      os.remove('./cutadapt/'+B+'/'+f)


### confirm number of primers mismatches:
for B in BR:
    print(B)
    temp = PRIMERLIST.loc[PRIMERLIST['BarcodeRegion']==B, ['AssayName']]
    for index, row in temp.iterrows():
        t = str(row['AssayName'])
        
        #Counting number of sequences by Taxid on inital input fasta file
        temp = count_seqsTaxId(FASTA)
        temp.rename(columns={'Nseqs':'OriginalFasta'}, inplace=True)
        
        #Counting number of sequences by Taxid after insilico
        if os.path.exists(f'./cutadapt/{B}/{t}.noPrimer.cutadapt'):
            temp2 = count_seqsTaxId(f'./cutadapt/{B}/{t}.noPrimer.cutadapt')
            temp = pd.merge(temp, temp2, how='left', on='TaxID')
            temp.rename(columns= {'Nseqs':'1stPCR'}, inplace=True)
            del temp2
        
        #Counting number of sequences by Taxid after inslilico w/ max error of 5
        if os.path.exists(f'./cutadapt/{B}/{t}.noPrimer.cutadapt'):
            temp2 = count_seqsTaxId(f'./cutadapt/{B}/{t}.maxerror')
            temp = pd.merge(temp, temp2, how='left', on='TaxID')
            temp.rename(columns= {'Nseqs':'insilicoMaxError'}, inplace=True)
            del temp2  
        
        #Counting number of sequences by Taxid after cleaning insilico (filter by ambiguous bases)
        if os.path.exists(f'./cutadapt/{B}/{t}.clean.cutadapt'):
            temp2 = count_seqsTaxId(f'./cutadapt/{B}/{t}.clean.cutadapt')
            temp = pd.merge(temp, temp2, how='left', on='TaxID')
            temp.rename(columns= {'Nseqs':'1stPCRclean'}, inplace=True)
            del temp2
        
        # Counting sequences after pga
        if os.path.exists(f'./cutadapt/{B}/{t}.pga'):
            temp2 = count_seqsTaxId(f'./cutadapt/{B}/{t}.pga')
            temp = pd.merge(temp, temp2, how='left', on='TaxID')
            temp.rename(columns= {'Nseqs':'PGA'}, inplace=True)
            temp['PGA'] = temp['PGA']-temp['1stPCRclean']
            del temp2
        
        # Couting sequences after cleaning pga clean (filter by ambiguous bases)
        if os.path.exists(f'./cutadapt/{B}/{t}.clean.pga'):
            temp2 = count_seqsTaxId(f'./cutadapt/{B}/{t}.clean.pga')
            temp = pd.merge(temp, temp2, how='left', on='TaxID')
            temp.rename(columns= {'Nseqs':'PGAclean'}, inplace=True)
            del temp2
            temp['PGAclean'] = temp['PGAclean']-temp['1stPCRclean']
            temp['TotalSeqs'] = temp['1stPCRclean'] + temp['PGAclean']
            temp.to_csv(f'./cutadapt/{B}/{t}.NumberSeqs.tsv', index=None, sep = '\t')

print('############# Step 6 Done ##############')
# -------------------------------------------------------------------------------


#%%
#######end of update

#read resolution
iberiaTaxonomy = pd.read_csv(TAXONOMY, sep = '\t')
iberiaTaxonomy['TaxID'] = iberiaTaxonomy['TaxID'].astype(int)
#iberiaTaxonomy = iberiaTaxonomy.loc[~iberiaTaxonomy['Sample'].isna()].reset_index(drop=True)


#Check resolution
for B in BR:
    print(B)
    temp = PRIMERLIST.loc[PRIMERLIST['BarcodeRegion']==B, 'AssayName'].tolist()
    for t in temp:
        if os.path.exists('./cutadapt/'+B+'/'+ t+'.cutadapt'):
            checkcutadapt('./cutadapt/'+B+'/'+ t+'.cutadapt', iberiaTaxonomy)



for B in BR:
    temp = PRIMERLIST.loc[PRIMERLIST['BarcodeRegion']==B, 'AssayName'].tolist()
    for t in temp:
        if os.path.exists('./cutadapt/'+B+'/'+ t+'.summary'):
            resolution = pd.read_csv('./cutadapt/'+B+'/'+ t+'.summary', sep = '\t')
            resolution.groupby(['Resolution'])['Species'].count().to_csv('./cutadapt/'+B+'/'+ t+'.counts', sep = '\t')





# %%
