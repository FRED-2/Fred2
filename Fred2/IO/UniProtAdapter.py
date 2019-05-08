# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: IO.UniProtAdapter
   :synopsis: DB-Adapter class for UniProt
.. moduleauthor:: walzer
.. deprecated:: 1.0
"""

import warnings
import bisect

from Bio import SeqIO
from Fred2.Core.Base import deprecated


class UniProtDB:
    @deprecated  # TODO: refactor ... function based on old code
    def __init__(self, name='fdb'):
        """
        UniProtDB class to give quick access to entries (fast exact match searches) and convenient ways to produce
        combined fasta files. Search is done with python's fast search  based on a mix between boyer-moore and horspool
        (http://svn.python.org/view/python/trunk/Objects/stringlib/fastsearch.h?revision=68811&view=markup)

        :param name: a name for the UniProtDB object

        Usage tutorials:
            import UniProtDB
            db = UniProtDB.UniProtDB('uniprot') #give it a name
            db.read_seqs('/path/to/file.fasta')
            l = list(SeqIO.parse(f, "fasta"))
            d = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
            db = UniProtDB.UniProtDB('something')
            db.read_seqs(l)
            db.read_seqs(d)
        """
        self.name = name
        self.collection = {}  # all the biopython seq records in a dict keyed by the id of the record
        self.searchstring = ''  # all sequences concatenated with a '#'
        self.accs = list()  # all accessions in respective order to searchstring
        self.idx = list()  # all indices of starting strings in the searchstring in respective order

    def read_seqs(self, sequence_file):
        """
        read sequences from uniprot files (.dat or .fasta) or from lists or dicts of BioPython SeqRecords
        and make them available for fast search. Appending also with this function.

        :param sequence_file: uniprot files (.dat or .fasta)
        :return:
        """
        recs = sequence_file
        if not isinstance(sequence_file, dict) and not isinstance(sequence_file, list):
            try:
                with open(sequence_file, 'rb') as f:
                    if sequence_file.endswith('.fa') or sequence_file.endswith('.fasta'):
                        recs = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
                    else:  # assume it is a dat file
                        recs = SeqIO.to_dict(SeqIO.parse(open(sequence_file), 'swiss'))
            except:
                warnings.warn("Could not read file", UserWarning)
                return
        if isinstance(sequence_file, list):
            recs = SeqIO.to_dict(sequence_file)
        if recs:
            self.collection.update(recs)
            self.searchstring = '#'.join([str(x.seq) for x in list(self.collection.values())]).decode('ascii')
            self.accs = list(self.collection.keys())
            self.idx = list()
            self.idx.append(0)
            for i, v in enumerate(self.collection.values()):
                self.idx.append(1 + self.idx[-1] + len(list(self.collection.values())[i].seq))
        return

    def write_seqs(self, name):
        """
        writes all fasta entries in the current object into one fasta file

        :param name: the complete path with file name where the fasta is going to be written
        """
        with open(name, "w") as output:
            SeqIO.write(list(self.collection.values()), output, "fasta")

    def exists(self, seq):
        """
        fast check if given sequence exists (as subsequence) in one of the UniProtDB objects collection of sequences.

        :param seq: the subsequence to be searched for
        :return: True, if it is found somewhere, False otherwise
        """
        if isinstance(seq, str):
            index = self.searchstring.find(seq)
            if index >= 0:
                return True
            else:
                return False
        return None

    def search(self, seq):
        """
        search for first occurrence of given sequence(s) in the UniProtDB objects collection returning (each) the fasta
        header front part of the first occurrence.

        :param seq: a string interpreted as a single sequence or a list (of str) interpreted as a coll. of sequences
        :return: a dictionary of sequences to lists (of ids, 'null' if n/a)
        """
        if isinstance(seq, str):
            ids = 'null'
            index = self.searchstring.find(seq)
            if index >= 0:
                j = bisect.bisect(self.idx, index) - 1
                ids = self.accs[j]
            return {seq: ids}
        if isinstance(seq, list):
            ids = list()
            for i in seq:
                ids.append('null')
            for i, v in enumerate(seq):
                index = self.searchstring.find(v)
                if index >= 0:
                    j = bisect.bisect(self.idx, index) - 1
                    ids[i] = self.accs[j]
            return dict(list(zip(seq, ids)))
        return None

    def search_all(self, seq):
        """
        search for all occurrences of given sequence(s) in the UniProtDB objects collection returning (each) the
        fasta header front part of all occurrences.

        :param seq: a string interpreted as a single sequence or a list (of str) interpreted as a coll. of sequences
        :return: a dictionary of the given sequences to lists (of ids, 'null' if n/a)
        """
        if isinstance(seq, str):
            ids = 'null'
            index = 0
            searchstring_length = len(seq)
            while index < len(self.searchstring):
                index = self.searchstring.find(seq, index)
                if index == -1:
                    break
                j = bisect.bisect(self.idx, index) - 1
                if ids == 'null':
                    ids = self.accs[j]
                else:
                    ids = ids + ',' + self.accs[j]
                index += searchstring_length
            return {seq: ids}
        if isinstance(seq, list):
            ids = list()
            for i in seq:
                ids.append('null')
            for i, v in enumerate(seq):
                index = 0
                searchstring_length = len(v)
                while index < len(self.searchstring):
                    index = self.searchstring.find(v, index)
                    if index == -1:
                        break
                    j = bisect.bisect(self.idx, index) - 1
                    if ids[i] == 'null':
                        ids[i] = self.accs[j]
                    else:
                        ids[i] = ids[i] + ',' + self.accs[j]
                    index += searchstring_length
            return dict(list(zip(seq, ids)))
        return None