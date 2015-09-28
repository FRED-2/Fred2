# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

__author__ = 'walzer'

import warnings
import bisect

from Bio import SeqIO
from Fred2.IO.ADBAdapter import ADBAdapter, EAdapterFields


class EnsemblDB(ADBAdapter):
    def __init__(self, name='fdb'):
        """
        EnsembleDB class to give quick access to entries (fast exact match searches) and convenient ways to produce
        combined fasta files. Search is done with python's fast search  based on a mix between boyer-moore and horspool
        (http://svn.python.org/view/python/trunk/Objects/stringlib/fastsearch.h?revision=68811&view=markup)
        :param name: a name for the EnsembleDB object
        Usage tutorials:
            import EnsembleDB
            db = EnsembleDB.EnsembleDB('Ensemble') #give it a name
            db.read_seqs('/path/to/file.fasta')
            l = list(SeqIO.parse(f, "fasta"))
            d = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
            db = EnsembleDB.EnsembleDB('something')
            db.read_seqs(l)
            db.read_seqs(d)
        """
        self.name = name
        self.collection = {}  # all the biopython seq records in a dict keyed by the id of the record
        self.searchstring = ''  # all sequences concatenated with a '#'
        self.accs = list()  # all accessions in respective order to searchstring
        self.idx = list()  # all indices of starting strings in the searchstring in respective order
        self.ensg2enst = dict()
        self.ensg2ensp = dict()
        self.enst2ensg = dict()
        self.enst2ensp = dict()
        self.ensp2ensg = dict()
        self.ensp2enst = dict()

    def read_seqs(self, sequence_file):
        """
        read sequences from Ensemble protein files (.fasta) or from lists or dicts of BioPython SeqRecords
        and make them available for fast search. Appending also with this function.
        :param sequence_file: Ensemble files (.dat or .fasta)
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
            self.searchstring = '#'.join([str(x.seq) for x in self.collection.values()]).decode('ascii')
            self.accs = self.collection.keys()
            self.idx = list()
            self.idx.append(0)
            for i, v in enumerate(self.collection.values()):
                self.idx.append(1 + self.idx[-1] + len(self.collection.values()[i].seq))

            for i in recs.items():
                ensg = None
                enst = None
                ensp = None
                if i[0].startswith('ENSG'):
                    ensg = i[0]
                elif i[0].startswith('ENST'):
                    enst = i[0]
                elif i[0].startswith('ENSP'):
                    ensp = i[0]
                ks = i[1].description.split(' ')
                for j in ks:
                    if j.startswith('transcript:'):
                        enst = j.strip('transcript:')
                    elif j.startswith('gene:'):
                        ensg = j.strip('gene:')
                if ensg not in self.ensg2enst:
                    self.ensg2enst[ensg] = list()
                    self.ensg2ensp[ensg] = list()
                self.ensg2enst[ensg].append(enst)
                self.ensg2ensp[ensg].append(ensp)
                if enst:
                    self.enst2ensp[enst] = ensp
                    self.enst2ensg[enst] = ensg
                else:
                    warnings.warn("Unparsable filecontents", UserWarning)
                if ensp:
                    self.ensp2ensg[ensp] = ensg
                    self.ensp2enst[ensp] = enst
                else:
                    warnings.warn("Unparsable filecontents", UserWarning)
        return

    def map_enst(self, enst):
        """
        looks up enst from the mapping and returns a ensg and a ensp
        :param enst: the transcript to map
        :return: tuple of gene and protein (might be None)
        """
        ensg = None
        ensp = None
        if enst in self.enst2ensg:
            ensg = self.enst2ensg[enst]
        if enst in self.enst2ensp:
            ensp = self.enst2ensp[enst]
        return {'Ensembl Gene ID': ensg, 'Ensembl Transcript ID': enst, 'Ensembl Protein ID': ensp}

    def map_ensp(self, ensp):
        """
        looks up ensp from the mapping and returns a ensg and a enst
        :param ensp: the protein to map
        :return: tuple of gene and protein (should not be None)
        """
        ensg = None
        enst = None
        if ensp in self.ensp2ensg:
            ensg = self.ensp2ensg[ensp]
        if enst in self.ensp2enst:
            enst = self.ensp2enst[ensp]
        return {'Ensembl Gene ID': ensg, 'Ensembl Transcript ID': enst, 'Ensembl Protein ID': ensp}

    def map_ensg(self, ensg):
        warnings.warn('mapping ensg not implemented', NotImplementedError)

    def get_transcript_sequence(self, transcript_id):
        if transcript_id in self.collection:
            return str(self.collection[transcript_id].seq)
        else:
            return None

    def get_product_sequence(self, product_id):
        if product_id in self.collection:
            return self.collection[product_id]
        else:
            return None

    def get_transcript_information(self, transcript_id):
        if transcript_id in self.collection:
            return {EAdapterFields.SEQ: str(self.collection[transcript_id].seq),
                    EAdapterFields.GENE: self.collection[transcript_id].description.split('gene:')[1].split(' ')[0],
                    EAdapterFields.STRAND: "-" if int(self.collection[transcript_id].description.split('chromosome:')[1].split(' ')[0].split(':')[-1]) < 0 else "+"}
        else:
            return None

    def write_seqs(self, name):
        """
            writes all fasta entries in the current object into one fasta file
            :param name: the complete path with file name where the fasta is going to be written
            """
        with open(name, "w") as output:
            SeqIO.write(self.collection.values(), output, "fasta")

    def exists(self, seq):
        """
            fast check if given sequence exists (as subsequence) in one of the EnsembleDB objects collection of sequences.
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
            search for first occurrence of given sequence(s) in the EnsembleDB objects collection returning (each) the fasta
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
            return dict(zip(seq, ids))
        return None

    def search_all(self, seq):
        """
            search for all occurrences of given sequence(s) in the EnsembleDB objects collection returning (each) the
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
            return dict(zip(seq, ids))
        return None