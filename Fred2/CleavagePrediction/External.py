# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: CleavagePrediction.ANN
   :synopsis: This module contains all classes for external cleavage prediction methods.
.. moduleauthor:: schubert

"""

import sys
import os
import subprocess
import pandas
import abc

from tempfile import NamedTemporaryFile

from Fred2.Core.Base import ACleavageSitePrediction, AExternal
from Fred2.Core.Protein import Protein
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Result import CleavageSitePredictionResult


class AExternalCleavageSitePrediction(ACleavageSitePrediction, AExternal):
    """
        Abstract base class for external cleavage site prediction methods.
        Implements predict functionality.
    """
    @abc.abstractmethod
    def prepare_input(self, input, file):
        """
        Prepares the data :attr:_input and writes them to :attr:_file in the special format used by the external tool

        :param list(str) input: The input data (here peptide sequences)
        :param File file: A file handler with which the data are written to file
        """
        raise NotImplementedError

    def predict(self, aa_seq, command=None, options=None, **kwargs):
        """
        Overwrites ACleavageSitePrediction.predict

        :param aa_seq: A list of or a single :class:`~Fred2.Core.Peptide.Peptide` or :class:`~Fred2.Core.Protein.Protein` object
        :type aa_seq: list(:class:`~Fred2.Core.Peptide.Peptide`/:class:`~Fred2.Core.Protein.Protein`) or :class:`~Fred2.Core.Peptide.Peptide`/:class:`~Fred2.Core.Protein.Protein`
        :param str command: The path to a alternative binary (can be used if binary is not globally executable)
        :param str options: A string of additional options directly past to the external tool
        :return: A :class:`~Fred2.Core.CleavageSitePredictionResult` object
        :rtype: :class:`~Fred2.Core.CleavageSitePredictionResult`
        """
        if not self.is_in_path() and "path" not in kwargs and command is None:
            raise RuntimeError("{name} {version} could not be found in PATH".format(name=self.name,
                                                                                    version=self.version))
        external_version = self.get_external_version(path=command)
        if self.version != external_version and external_version is not None:
            raise RuntimeError("Internal version {internal_version} does "
                               "not match external version {external_version}".format(internal_version=self.version,
                                                                                      external_version=external_version))

        if isinstance(aa_seq, Peptide) or isinstance(aa_seq, Protein):
            pep_seqs = {str(aa_seq): aa_seq}
        else:
            pep_seqs = {}
            for p in aa_seq:
                if not isinstance(p, Peptide) and not isinstance(p, Protein):
                    raise ValueError("Input is not of type Protein or Peptide")
                pep_seqs[str(p)] = p

        chunksize = len(pep_seqs)
        if 'chunks' in kwargs:
            chunksize = kwargs['chunks']

        result = {}
        peps = list(pep_seqs.values())

        for i in range(0, len(peps), chunksize):
            tmp_out = NamedTemporaryFile(delete=False)
            tmp_file = NamedTemporaryFile(delete=False)
            self.prepare_input(peps[i:i+chunksize], tmp_file)
            tmp_file.close()

            #allowe customary executable specification
            if command is not None:
                exe = self.command.split()[0]
                _command = self.command.replace(exe, command)
            else:
                _command = self.command

            try:
                stdo = None
                stde = None
                cmd = _command.format(input=tmp_file.name, options="" if options is None else options, out=tmp_out.name)
                p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdo, stde = p.communicate()
                stdr = p.returncode
                if stdr > 0:
                    raise RuntimeError("Unsuccessful execution of " + cmd + " (EXIT!=0) with error: " + stde)
            except Exception as e:
                raise RuntimeError(e)

            result.update(self.parse_external_result(tmp_out))
            os.remove(tmp_file.name)
            tmp_out.close()
            os.remove(tmp_out.name)

        df_result = CleavageSitePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i,j)) for i, j in df_result.index],
                                                        names=['ID', 'Pos'])

        return df_result


class NetChop_3_1(AExternalCleavageSitePrediction, AExternal):
    """
    Implements NetChop Cleavage Site Prediction (v. 3.1).

    .. note::

        Nielsen, M., Lundegaard, C., Lund, O., & Kesmir, C. (2005). The role of the proteasome in generating cytotoxic
        T-cell epitopes: insights obtained from improved predictions of proteasomal cleavage.
        Immunogenetics, 57(1-2), 33-41.
    """
    __supported_length = [sys.maxsize]
    __name = "netchop"
    __cleavage_pos = 0
    __command = "netchop {input} {options} | grep -v '#' > {out}"
    __version = "3.1"

    @property
    def version(self):
        """
        The version of the Method
        """
        return self.__version

    @property
    def command(self):
        """
        Defines the commandline call for external tool
        """
        return self.__command

    @property
    def supportedLength(self):
        """The supported lengths of the predictor"""
        return self.__supported_length

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    @property
    def cleavagePos(self):
        """Parameter specifying the position of aa (within the prediction window) after which the sequence is cleaved"""
        return self.__cleavage_pos

    def parse_external_result(self, file, hash_dict):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :param dict(int,Protein) hash_dict: A hash dict mapping from internal ID to Protein object
                                            du to ID length restriction of NetChop
        :return: Returns a dictionary with the prediction results
        :rtype: dict(str,dict((str,int),float))
        """
        result = {"Seq": {}, self.name: {}}
        count = 0
        is_new_seq = 0

        for l in file:
            l = l.strip()
            if not len(l):
                continue
            if not is_new_seq % 4 and is_new_seq:
                count += 1
                is_new_seq = 0
            elif l[0] == "-":
                is_new_seq += 1
            elif l[0].isdigit():
                pos, aa, _, s, ident = l.split()
                pos = int(pos) - 1
                try:
                    seq_id = hash_dict[int(ident)].transcript_id
                except:
                    seq_id = int(ident)
                result["Seq"][(seq_id, pos)] = aa
                result[self.name][(seq_id, pos)] = float(s)

        return result

    def get_external_version(self, path=None):
        """
        Returns the external version of the tool by executing
        >{command} --version

        might be dependent on the method and has to be overwritten
        therefore it is declared abstract to enforce the user to
        overwrite the method. The function in the base class can be called
        with super()

        :param str path: - Optional specification of executable path if deviant from :attr:self.__command
        :return: The external version of the tool or None if tool does not support versioning
        :rtype: str
        """
        #cannot be determined method does not support --version or something similar
        return None

    def prepare_input(self, input, file):
        """
        Prepares the data and writes them to _file in the special format used by the external tool

        :param dict(int,Protein/Peptide) input: The input data (here peptide sequences)
        :param File file: A file handler with which the data are written to file
        """
        file.write("\n".join(">%i\n%s"%(i , str(p)) for i, p in input))

    def predict(self, aa_seq, command=None, options=None, **kwargs):
        """
        Overwrites ACleavageSitePrediction.predict

        :param aa_seq: A list of or a single :class:`~Fred2.Core.Peptide.Peptide` or :class:`~Fred2.Core.Protein.Protein` object
        :type aa_seq: list(:class:`~Fred2.Core.Peptide.Peptide`/:class:`~Fred2.Core.Protein.Protein`) or :class:`~Fred2.Core.Peptide.Peptide`/:class:`~Fred2.Core.Protein.Protein`
        :param str command: The path to a alternative binary (can be used if binary is not globally executable)
        :param str options: A string of additional options directly past to the external tool
        :return: A :class:`~Fred2.Core.CleavageSitePredictionResult` object
        :rtype: :class:`~Fred2.Core.CleavageSitePredictionResult`
        """
        if not self.is_in_path() and "path" not in kwargs and command is None:
            raise RuntimeError("{name} {version} could not be found in PATH".format(name=self.name,
                                                                                    version=self.version))
        external_version = self.get_external_version(path=command)
        if self.version != external_version and external_version is not None:
            raise RuntimeError("Internal version {internal_version} does "
                               "not match external version {external_version}".format(internal_version=self.version,
                                                                                      external_version=external_version))

        #Since NetChop 3.1 cuts identifiers to 10-digits we use
        #An integer hashing to generate unique ids for
        #With is we can predict 99999999999 protein sequences simultaniously
        #After reaching the limit one could revers the counting in the negative
        #direction
        if isinstance(aa_seq, Peptide) or isinstance(aa_seq, Protein):
            pep_seqs = {0: str(aa_seq)}
        else:
            pep_seqs = {}
            for i,p in enumerate(aa_seq):
                if not isinstance(p, Peptide) and not isinstance(p, Protein):
                    raise ValueError("Input is not of type Protein or Peptide")
                if i < 99999999999:
                    pep_seqs[i] = p
                else:
                    pep_seqs[i-99999999999] = p

        chunksize = len(pep_seqs)
        if 'chunks' in kwargs:
            chunksize = kwargs['chunks']

        result = {}
        peps = list(pep_seqs.items())

        for i in range(0, len(peps), chunksize):

            tmp_out = NamedTemporaryFile(delete=False)
            tmp_file = NamedTemporaryFile(delete=False)
            self.prepare_input(peps[i:i+chunksize], tmp_file)
            tmp_file.close()

            #allowe customary executable specification
            if command is not None:
                exe = self.command.split()[0]
                _command = self.command.replace(exe, command)
            else:
                _command = self.command

            try:
                stdo = None
                stde = None
                cmd = _command.format(input=tmp_file.name, options="" if options is None else options, out=tmp_out.name)
                p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdo, stde = p.communicate()
                stdr = p.returncode
                if stdr > 0:
                    raise RuntimeError("Unsuccessful execution of " + cmd + " (EXIT!=0) with error: " + stde)
            except Exception as e:
                raise RuntimeError(e)

            result.update(self.parse_external_result(tmp_out, pep_seqs))

            os.remove(tmp_file.name)
            tmp_out.close()
            os.remove(tmp_out.name)

        df_result = CleavageSitePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i, j)) for i, j in df_result.index],
                                                        names=['ID', 'Pos'])

        return df_result
