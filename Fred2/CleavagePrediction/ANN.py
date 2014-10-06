# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
__author__ = 'schubert'

import collections
import warnings
import sys
import subprocess
import pandas

from collections import defaultdict
from tempfile import NamedTemporaryFile

from Fred2.Core.Base import ACleavageSitePrediction, AExternal
from Fred2.Core.Protein import Protein
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Result import CleavageSitePredictionResult


class NetChop(ACleavageSitePrediction, AExternal):
    """
    Implements NetChop Cleavage Site Prediction (v. 3.1)
    """
    __supported_length = [sys.maxint]
    __name = "netchop"
    __cleavage_pos = 0
    __command = "netChop %s | grep -v '#' > %s"

    @property
    def command(self):
        return self.__command

    @property
    def supportedLength(self):
        return self.__supported_length

    @property
    def name(self):
        return self.__name

    @property
    def cleavagePos(self):
        return self.__cleavage_pos

    def predict(self, _aa_seq, **kwargs):

        if isinstance(_aa_seq, collections.Iterable):
            if any((not isinstance(p, Peptide)) and (not isinstance(p, Protein)) for p in _aa_seq):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(p):p for p in _aa_seq}
        else:
            if (not isinstance(_aa_seq, Peptide)) or (not isinstance(_aa_seq, Protein)):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(_aa_seq):_aa_seq}

        tmp_out = NamedTemporaryFile(delete=False)
        tmp_file = NamedTemporaryFile(delete=False)
        tmp_file.write("\n".join(">pep_%i\n%s"%(i,str(p)) for i, p in enumerate(pep_seqs.iterkeys())))
        tmp_file.close()

        r = subprocess.call(self.command%(tmp_file.name, tmp_out.name), shell=True)

        if r != 0:
            warnings.warn("An unknown error occurred for method %s"%self.name)
            sys.exit(-1)

        return self.parse_external_result(tmp_out)


    def parse_external_result(self, _file):
        result = {"Seq":{}, self.name:{}}
        count = 0
        is_new_seq = 0

        for l in _file:
            l = l.strip()
            if not len(l):
                continue
            #print l
            if not is_new_seq % 4 and is_new_seq:
                #print "New seq starts", l
                count +=1
                is_new_seq = 0
            elif l[0] == "-":
                #print "in counter", l
                is_new_seq += 1
            elif l[0].isdigit():
                pos,aa,_,s,_ = l.split()
                pos = int(pos)
                seq_id = "seq_%i"%count
                result["Seq"][(seq_id, pos)] = aa
                result[self.name][(seq_id, pos)] = float(s)
        print result
        df_result = CleavageSitePredictionResult.from_dict(result)
        print df_result.index
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i,j)) for i,j in df_result.index],
                                                        names=['ID','Pos'])
        return df_result
