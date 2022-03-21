# 1. finding all fragments in fasta or multifasta file
# 2. extracting translations for each fragments
# 3. extracting ORF for each translation

import Utils
from msr import MatrixSynthReaction


WORK_PATH = "C:\\bio\\orf-finder"


class OrfFinder:
    def __init__(self, src_path, snk_path):
        self.data_src = src_path
        self.data_snk = snk_path

        self.data_fragments = []
        self.potential_genes = []

    def get_fragments_from_data(self):
        pass

    def get_translation_from_fragment(self, fragment, shift):
        """
        :type shift: int
        might be in [1, 2, 3, -1, -2, -3]
        """
        pass

    def get_ORFs(self, frames):
        pass


