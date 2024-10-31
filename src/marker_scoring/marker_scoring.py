class ReferenceDatabaseQuality:
    def __init__(self, all_inserts_file, otl):
        self.all_inserts_file = all_inserts_file
        self.otl = otl

    def calculate_number_of_barcodes_per_taxon(self,
                                               all_inserts_file):
        """
        This method calculates the total number of barcodes that exists per taxon.

        Parameters:
        - all_inserts_file: File containing all inserts present in the reference database, either
        successfully amplified or not.
        """
        pass

    def calculate_percentage_of_taxa_w_x_barcodes(self,
                                                  barcodes_per_species_dict,
                                                  otl,
                                                  barcode_threshold=1):
        """
        This method calculates the percentage of taxa with more than X barcodes.

        Parameters:
        - barcodes_per_species_dict: Dictionary containing the exisiting number of barcodes per
        taxon. Species name are considered as keys, and number of barcodes are values.
        - otl: Operational Taxonomy List. List of taxons considered in a country for biomonitoring
        purposes.
        - barcode_threshold: Number of barcodes to consider as a threshold. Only taxons with more
        than this value will be considered.

        Output:
        - percentage_of_taxa: float
            Decimal percentage of taxa with more than X barcodes
        """
        pass

    def barcode_coverage_score(self, barcoded_taxa_one, barcoded_taxa_two):
        """
        This method calculates the Barcode Coverage Score (BCS).
        * Will be used as acceptance criteria for ranking primer pair performance (>0.65)

        Parameters:
        - barcoded_taxa_one: float
            Percetntage of taxa with more than five barcodes.
        - barcoded_taxa_two: float
            Percentage of taxa with more than one barcode.

        Output:
        - bcs: float
            Barcode Coverage Score
        """
        bcs = 0.75 * barcoded_taxa_one + 0.25 * barcoded_taxa_two

        return bcs