# 1. finding all fragments in fasta or multifasta file
# 2. extracting translations for each fragments
# 3. extracting ORF for each translation

from Utils import IOUtils
from msr import MatrixSynthReaction

WORK_PATH = "C:\\bio\\orf-finder\\"


class OrfFinder:
    def __init__(self, src_path, snk_path, data):
        self.data_src = src_path
        self.data_snk = snk_path
        self.data = data

        self.data_fragments = []
        self.potential_genes = []

    def set_data_fragments(self, fragments):
        self.data_fragments = fragments

    def get_data_fragments(self):
        return self.data_fragments

    def parse_fragments_from_data(self):
        d = len(self.data)
        header = None
        fragments = []
        seq = ''
        for i in range(d):
            line = self.data[i][:-1]
            if line[0] == '>':
                if header:
                    fragments.append((header, seq))
                header = str(line[1:])
                seq = ""
            else:
                seq += line
        fragments.append((header, seq))
        self.set_data_fragments(fragments)

    def get_translation_from_fragment(self, fragment, shift):
        """
        :type shift: int
        might be in [1, 2, 3, -1, -2, -3]
        """
        pass

    def get_ORFs(self, frames):
        pass


if __name__ == '__main__':
    src_path = WORK_PATH + "dna-seq.fna"
    snk_path = WORK_PATH + "orfs.txt"
    io = IOUtils()
    data = io.open_file(src_path)
    orf_finder = OrfFinder(src_path, snk_path, data)

    orf_finder.parse_fragments_from_data()
    fragments = orf_finder.get_data_fragments()
    # print(*fragments, sep='\\n')
    print(fragments[0][1])
