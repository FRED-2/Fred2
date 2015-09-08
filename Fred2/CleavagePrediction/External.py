# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: CleavagePrediction.ANN
   :synopsis: ANN-based cleavage prediction methods.
.. moduleauthor:: schubert

"""

import sys
import os
import subprocess
import pandas

from tempfile import NamedTemporaryFile

from Fred2.Core.Base import ACleavageSitePrediction, AExternal
from Fred2.Core.Protein import Protein
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Result import CleavageSitePredictionResult


class NetChop_3_1(ACleavageSitePrediction, AExternal):
    """
    Implements NetChop Cleavage Site Prediction (v. 3.1)
    """
    __supported_length = [sys.maxint]
    __name = "netchop"
    __cleavage_pos = 0
    __command = "netChop {input} {options} | grep -v '#' > {out}"
    __version = "3.1"

    @property
    def version(self):
        return self.__version

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

        if isinstance(_aa_seq, Peptide) or isinstance(_aa_seq, Protein):
            pep_seqs = {str(_aa_seq): _aa_seq}
        else:
            if any((not isinstance(p, Peptide)) and (not isinstance(p, Protein)) for p in _aa_seq):
                raise ValueError("Input is not of type Protein or Peptide")
            pep_seqs = {str(p): p for p in _aa_seq}

        tmp_out = NamedTemporaryFile(delete=False)
        tmp_file = NamedTemporaryFile(delete=False)
        tmp_file.write("\n".join(">pep_%i\n%s"%(i, str(p)) for i, p in enumerate(pep_seqs.iterkeys())))
        tmp_file.close()

        #allowe customary executable specification
        if "path" in kwargs:
            exe = self.command.split()[0]
            _command = self.command.replace(exe, kwargs["path"])
        else:
            _command = self.command

        try:
            stdo = None
            stde = None
            cmd = _command.format(peptides=tmp_file.name, options=kwargs.get("options", ""), out=tmp_out.name)
            p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait() #block the rest
            stdo, stde = p.communicate()
            stdr = p.returncode
            if stdr > 0:
                raise RuntimeError("Unsuccessful execution of " + cmd + " (EXIT!=0) with error: " + stde)
        except Exception as e:
            raise RuntimeError(e)

        result = self.parse_external_result(tmp_out)

        df_result = CleavageSitePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i,j)) for i, j in df_result.index],
                                                        names=['ID', 'Pos'])
        os.remove(tmp_file.name)
        tmp_out.close()
        os.remove(tmp_out.name)

        return df_result

    def parse_external_result(self, _file):
        result = {"Seq": {}, self.name: {}}
        count = 0
        is_new_seq = 0

        for l in _file:
            l = l.strip()
            if not len(l):
                continue
            #print l
            if not is_new_seq % 4 and is_new_seq:
                #print "New seq starts", l
                count += 1
                is_new_seq = 0
            elif l[0] == "-":
                #print "in counter", l
                is_new_seq += 1
            elif l[0].isdigit():
                pos, aa, _, s, _ = l.split()
                pos = int(pos) - 1
                seq_id = "seq_%i"%count
                result["Seq"][(seq_id, pos)] = aa
                result[self.name][(seq_id, pos)] = float(s)

        return result

    def get_external_version(self):
        #cannot be determined method does not support --version or something similar
        return None