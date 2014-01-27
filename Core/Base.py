from collections import defaultdict


class MetadataLogger(object):
    # This class provides a simple interface for assigning additional metadata to any object in our
    # data model. Examples: storing ANNOVAR columns like depth, base count, dbSNP id, quality
    # information for variants, prediction scores for peptides etc.
    # Normally used by custom toolbox functions and importers.

    def __init__(self):
        self.metadata = defaultdict(list)

    def log_metadata(self, label, value):
        self.metadata[label].append(value)

    def get_metadata(self, label, only_first=False):
        # although defaultdict *would* return [] if it didn't find label in self.metadata, it would
        # come with the side effect of adding label as a key to the defaultdict, so a getter is
        # justified in this case.
        if not only_first:
            return self.metadata[label] if label in self.metadata else []
        else:
            return self.metadata[label][0] if self.metadata[label] else None


class FrameshiftNode(object):
    # This class helps us generate all possible frameshift combinations on a transcript.
    # We represent them as a binary tree, where a path from the root to a leaves is a valid
    # frameshift combination. The first level of the tree corresponds to the first frameshift
    # appearing on the transcript, the second level to the second frameshift and so on.
    # The path from the root to a leaf is a sequence of steps each taken either to the left or
    # the right, where a left step means frameshift is "on" and a right step means "off".
    # As homozygous frameshifts can't be "off" their right child is always a None leaf.
    # TODO: explain how the path enumerator function gives you the legal combinations in
    # such an order that you can use it to your advantage when pruning equivalent combinations
    # thanks to early stop codons.

    def __init__(self, fs, next=None, genotype='het'):
        assert genotype in ('hom', 'het'), "Frameshift genotype must be either 'hom' or 'het'"
        self.fs = fs
        self.on = next  # next FrameshiftNode on transcript, None if this is the last one
        self.off = None if genotype=='hom' else next  # if FS is 'hom' it can't be "off"

    def traverse(self):
        tron = self.on.traverse() if self.on is not None else []
        troff = self.off.traverse() if self.off is not None else []

        if self.on is None and self.off is None:
            return [[]]

        myresult = []

        for a in tron:
            myresult.append([self] + a)
        for a in troff:
            myresult.append(a)
        #print 'traversed ', self.fs, myresult
        return myresult

    def __repr__(self):
        return str(self.fs)
