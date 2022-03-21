# 1. finding all fragments in fasta or multifasta file
# 2. extracting translations for each fragments
# 3. extracting ORF for each translation
from typing import Any

from Utils import IOUtils
from msr import MatrixSynthReaction

WORK_PATH = "C:\\bio\\orf-finder\\"


class OrfFinder:
    def __init__(self, src_path, snk_path, data):

        self.human_sartcodons = ["ATG"]
        self.human_stopcodons = ["TAA", "TAG", "TGA"]

        self.mitochondrion_startcodons = ["ATA", "ATT", "GTG", "TTG"]

        self.data_src = src_path
        self.data_snk = snk_path
        self.data = data

        self.data_fragments = []
        self.frames = []
        self.orfs = []

    def set_data_fragments(self, fragments):
        self.data_fragments = fragments

    def get_data_fragments(self):
        return self.data_fragments

    def set_frames(self, frames):
        self.frames = frames

    def get_frames(self):
        return self.frames

    def set_orfs(self, orfs):
        self.orfs = orfs

    def get_orfs(self):
        return self.orfs

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

    def parse_translation_from_fragments(self, fragments):
        """
        :type shift: int
        might be in [1, 2, 3, -1, -2, -3]
        """
        frames = []  # storing the six frame translation that it zould be extacted from the fragments
        for index, value in enumerate(fragments):  # looping over the fragments extracted
            dna = value[1]  # extract the fragment
            description = value[0]  # extact the desciption even were not use it, just for learning purpose
            reverseCdna = []  # storing the reverse compliments
            # create the positive frames
            # split the frames into codons for better performance
            frames.append([dna[i:i + 3] for i in range(0, len(dna), 3)])
            frames.append([dna[i:i + 3] for i in range(1, len(dna), 3)])
            frames.append([dna[i:i + 3] for i in range(2, len(dna), 3)])
            # reverse compliment of the fragment
            reverse = {"A": "T", "C": "G", "T": "A", "G": "C"}
            for i in range(len(dna)):
                reverseCdna.append(reverse[dna[-i - 1]]) if dna[-i - 1] in reverse.keys() else reverseCdna.append(
                    dna[-i - 1])  # if any contamination found we keep it for further more check
            reverseCdna = ''.join(reverseCdna)  # joining
            # create the negative frames
            frames.append([reverseCdna[i:i + 3] for i in range(0, len(reverseCdna), 3)])
            frames.append([reverseCdna[i:i + 3] for i in range(1, len(reverseCdna), 3)])
            frames.append([reverseCdna[i:i + 3] for i in range(2, len(reverseCdna), 3)])
        self.set_frames(frames)

    def parse_orfs(self, frames, startcodons, stopcodons):
        orfs = []
        frames = self.get_frames()
        for i in range(0, len(frames), 1):
            start = 0
            while start < len(frames[i]):
                if frames[i][start] in startcodons:
                    for stop in range(start + 1, len(frames[i]), 1):
                        if frames[i][stop] in stopcodons:
                            orfs.append(frames[i][start:stop])  # retrieve the orf
                            start = stop + 1  # avoiding multiple start codons
                            break
                start += 1
        self.set_orfs(orfs)


if __name__ == '__main__':
    src_path = WORK_PATH + "dna-seq.fna"
    snk_path = WORK_PATH + "orfs.txt"
    io = IOUtils()
    data = io.open_file(src_path)

    orf_finder = OrfFinder(src_path, snk_path, data)
    orf_finder.parse_fragments_from_data()
    orf_finder.parse_translation_from_fragments(orf_finder.data_fragments)
    orf_finder.parse_orfs(orf_finder.get_frames(),
                          startcodons=orf_finder.human_sartcodons,
                          stopcodons=orf_finder.human_stopcodons)
    print(*orf_finder.get_orfs(), sep='\n')
