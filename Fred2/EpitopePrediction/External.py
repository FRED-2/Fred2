# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: EpitopePrediction.ANN
   :synopsis: This module contains all classes for ANN-based epitope prediction methods.
.. moduleauthor:: schubert, walzer

"""
import abc

import itertools
import warnings
import logging
import pandas
import subprocess
import csv
import os
import math

from collections import defaultdict

from Fred2.Core.Allele import Allele, CombinedAllele, MouseAllele
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Result import EpitopePredictionResult
from Fred2.Core.Base import AEpitopePrediction, AExternal
from tempfile import NamedTemporaryFile


class AExternalEpitopePrediction(AEpitopePrediction, AExternal):
    """
        Abstract class representing an external prediction function. Implementations shall wrap external binaries by
        following the given abstraction.
    """

    @abc.abstractmethod
    def prepare_input(self, input, file):
        """
        Prepares input for external tools
        and writes them to _file in the specific format

        NO return value!

        :param: list(str) _input: The :class:`~Fred2.Core.Peptide.Peptide` sequences to write into file
        :param File file: File-handler to input file for external tool
        """
        return NotImplementedError

    def predict(self, peptides, alleles=None, command=None, options=None, **kwargs):
        """
        Overwrites AEpitopePrediction.predict

        :param peptides: A list of or a single :class:`~Fred2.Core.Peptide.Peptide` object
        :type peptides: list(:class:`~Fred2.Core.Peptide.Peptide`) or :class:`~Fred2.Core.Peptide.Peptide`
        :param alleles: A list of or a single :class:`~Fred2.Core.Allele.Allele` object. If no
                        :class:`~Fred2.Core.Allele.Allele` are provided, predictions are made for all
                        :class:`~Fred2.Core.Allele.Allele` supported by the prediction method
        :type alleles: list(:class:`~Fred2.Core.Allele.Allele`)/:class:`~Fred2.Core.Allele.Allele`
        :param str command: The path to a alternative binary (can be used if binary is not globally executable)
        :param str options: A string of additional options directly past to the external tool.
        :keyword chunksize: denotes the chunksize in which the number of peptides are bulk processed
        :return: A :class:`~Fred2.Core.Result.EpitopePredictionResult` object
        :rtype: :class:`~Fred2.Core.Result.EpitopePredictionResult`
        """
        if not self.is_in_path() and command is None:
            raise RuntimeError("{name} {version} could not be found in PATH".format(name=self.name,
                                                                                    version=self.version))
        external_version = self.get_external_version(path=command)
        if self.version != external_version and external_version is not None:
            raise RuntimeError("Internal version {internal_version} does "
                               "not match external version {external_version}".format(internal_version=self.version,
                                                                                      external_version=external_version))

        if isinstance(peptides, Peptide):
            pep_seqs = {str(peptides): peptides}
        else:
            pep_seqs = {}
            for p in peptides:
                if not isinstance(p, Peptide):
                    raise ValueError("Input is not of type Protein or Peptide")
                pep_seqs[str(p)] = p

        chunksize = len(pep_seqs)
        if 'chunks' in kwargs:
            chunksize = kwargs['chunks']

        if alleles is None:
            al = [Allele(a) for a in self.supportedAlleles]
            allales_string = {conv_a: a for conv_a, a in zip(self.convert_alleles(al), al)}
        else:
            if isinstance(alleles, Allele):
                alleles = [alleles]
            if any(not isinstance(p, Allele) for p in alleles):
                raise ValueError("Input is not of type Allele")
            allales_string = {conv_a: a for conv_a, a in zip(self.convert_alleles(alleles), alleles)}

        result = defaultdict(defaultdict)

        # group alleles in blocks of 80 alleles (NetMHC can't deal with more)
        _MAX_ALLELES = 50

        # allow custom executable specification
        if command is not None:
            exe = self.command.split()[0]
            _command = self.command.replace(exe, command)
        else:
            _command = self.command

        allele_groups = []
        c_a = 0
        allele_group = []
        for a in allales_string.keys():
            if c_a >= _MAX_ALLELES:
                c_a = 0
                allele_groups.append(allele_group)
                if str(allales_string[a]) not in self.supportedAlleles:
                    logging.warn("Allele %s is not supported by %s" % (str(allales_string[a]), self.name))
                    allele_group = []
                    continue
                allele_group = [a]
            else:
                if str(allales_string[a]) not in self.supportedAlleles:
                    logging.warn("Allele %s is not supported by %s" % (str(allales_string[a]), self.name))
                    continue
                allele_group.append(a)
                c_a += 1

        if len(allele_group) > 0:
            allele_groups.append(allele_group)
        # export peptides to peptide list

        pep_groups = list(pep_seqs.keys())
        pep_groups.sort(key=len)
        for length, peps in itertools.groupby(pep_groups, key=len):
            if length not in self.supportedLength:
                logging.warn("Peptide length must be at least %i or at most %i for %s but is %i" % (min(self.supportedLength), max(self.supportedLength),
                                                                                       self.name, length))
                continue
            peps = list(peps)
            for i in range(0, len(peps), chunksize):
                tmp_out = NamedTemporaryFile(delete=False)
                tmp_file = NamedTemporaryFile(delete=False)
                self.prepare_input(peps[i:i+chunksize], tmp_file)
                #            tmp_file.write("\n".join(">pepe_%i\n%s"%(i, p) for i, p in enumerate(peps))
                #                           if self.name.lower() in ["netmhcii","netctlpan"] else "\n".join(peps))
                tmp_file.close()

                # generate cmd command
                for allele_group in allele_groups:
                    try:
                        stdo = None
                        stde = None
                        cmd = _command.format(peptides=tmp_file.name, alleles=",".join(allele_group),
                                              options="" if options is None else options, out=tmp_out.name, length=str(length))
                        p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE)
                        # p.wait() communicate already waits for the process https://docs.python.org/2.7/library/subprocess.html#subprocess.Popen.communicate
                        stdo, stde = p.communicate()
                        stdr = p.returncode
                        if stdr > 0:
                            raise RuntimeError("Unsuccessful execution of " + cmd + " (EXIT!=0) with error: " + stde)
                    except Exception as e:
                        raise RuntimeError(e)

                    res_tmp = self.parse_external_result(tmp_out.name)
                    for al, ep_dict in res_tmp.items():
                        for p, v in ep_dict.items():
                            result[allales_string[al]][pep_seqs[p]] = v
                os.remove(tmp_file.name)
                tmp_out.close()
                os.remove(tmp_out.name)

        if not result:
            raise ValueError("No predictions could be made with " + self.name +
                             " for given input. Check your epitope length and HLA allele combination.")


        df_result = EpitopePredictionResult.from_dict(result)
        df_result.index = pandas.MultiIndex.from_tuples([tuple((i, self.name)) for i in df_result.index],
                                                        names=['Seq', 'Method'])
        return df_result


class NetMHC_3_4(AExternalEpitopePrediction):
    """
    Implements the NetMHC binding (in current form for netMHC3.4).

    .. note::

        NetMHC-3.0: accurate web accessible predictions of human, mouse and monkey MHC class I affinities for peptides
        of length 8-11. Lundegaard C, Lamberth K, Harndahl M, Buus S, Lund O, Nielsen M.
        Nucleic Acids Res. 1;36(Web Server issue):W509-12. 2008

        Accurate approximation method for prediction of class I MHC affinities for peptides of length 8, 10 and 11 using
        prediction tools trained on 9mers. Lundegaard C, Lund O, Nielsen M. Bioinformatics, 24(11):1397-98, 2008.

    """

    __alleles = frozenset(['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:02', 'HLA-A*02:03', 'HLA-A*02:06', 'HLA-A*02:11', 'HLA-A*02:12', 'HLA-A*02:16',
                           'HLA-A*02:17', 'HLA-A*02:19', 'HLA-A*02:50', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*24:03',
                           'HLA-A*25:01', 'HLA-A*26:01', 'HLA-A*26:02', 'HLA-A*26:03', 'HLA-A*29:02', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01',
                           'HLA-A*32:01', 'HLA-A*32:07', 'HLA-A*32:15', 'HLA-A*33:01', 'HLA-A*66:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-A*68:23',
                           'HLA-A*69:01', 'HLA-A*80:01', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*08:02', 'HLA-B*08:03', 'HLA-B*14:02', 'HLA-B*15:01',
                           'HLA-B*15:02', 'HLA-B*15:03', 'HLA-B*15:09', 'HLA-B*15:17', 'HLA-B*18:01', 'HLA-B*27:05', 'HLA-B*27:20', 'HLA-B*35:01',
                           'HLA-B*35:03', 'HLA-B*38:01', 'HLA-B*39:01', 'HLA-B*40:01', 'HLA-B*40:02', 'HLA-B*40:13', 'HLA-B*42:01', 'HLA-B*44:02',
                           'HLA-B*44:03', 'HLA-B*45:01', 'HLA-B*46:01', 'HLA-B*48:01', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*54:01', 'HLA-B*57:01',
                           'HLA-B*58:01', 'HLA-B*73:01', 'HLA-B*83:01', 'HLA-C*03:03', 'HLA-C*04:01', 'HLA-C*05:01', 'HLA-C*06:02', 'HLA-C*07:01',
                           'HLA-C*07:02', 'HLA-C*08:02', 'HLA-C*12:03', 'HLA-C*14:02', 'HLA-C*15:02', 'HLA-E*01:01',
                           'H2-Db', 'H2-Dd', 'H2-Kb', 'H2-Kd', 'H2-Kk', 'H2-Ld'])
    __supported_length = frozenset([8, 9, 10, 11])
    __name = "netmhc"
    __command = "netMHC -p {peptides} -a {alleles} -x {out} {options}"
    __version = "3.4"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    def _represent(self, allele):
        """
        Internal function transforming an allele object into its representative string
        :param allele: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: str
        """
        if isinstance(allele, MouseAllele):
            return "H-2-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)
        else:
            return "HLA-%s%s:%s" % (allele.locus, allele.supertype, allele.subtype)

    def convert_alleles(self, alleles):
        """
        Converts :class:`~Fred2.Core.Allele.Allele` into the internal :class:`~Fred2.Core.Allele.Allele` representation
        of the predictor and returns a string representation

        :param alleles: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`
        :rtype: list(str)
        """
        return [self._represent(a) for a in alleles]

    @property
    def supportedAlleles(self):
        """
        A list of valid allele models
        """
        return self.__alleles

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    @property
    def command(self):
        """
        Defines the commandline call for external tool
        """
        return self.__command

    @property
    def supportedLength(self):
        """
        A list of supported :class:`~Fred2.Core.Peptide.Peptide` lengths
        """
        return self.__supported_length

    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        result = defaultdict(defaultdict)
        f = csv.reader(open(file, "r"), delimiter='\t')
        next(f)
        next(f)
        alleles = [x.split()[0] for x in f.next()[3:]]
        for l in f:
            if not l:
                continue
            pep_seq = l[2]
            for ic_50, a in zip(l[3:], alleles):
                sc = 1.0 - math.log(float(ic_50), 50000)
                result[a][pep_seq] = sc if sc > 0.0 else 0.0
        return dict(result)

    def get_external_version(self, path=None):
        """
        Returns the external version of the tool by executing
        >{command} --version

        might be dependent on the method and has to be overwritten
        therefore it is declared abstract to enforce the user to
        overwrite the method. The function in the base class can be called
        with super()

        :param str path: Optional specification of executable path if deviant from :attr:`self.__command`
        :return: The external version of the tool or None if tool does not support versioning
        :rtype: dict
        """
        return super(NetMHC_3_4, self).get_external_version()

    def prepare_input(self, input, file):
        """
        Prepares input for external tools
        and writes them to file in the specific format

        NO return value!

        :param: list(str) input: The : sequences to write into _file
        :param File file: File-handler to input file for external tool
        """
        file.write("\n".join(input))


class NetMHC_3_0(NetMHC_3_4):
    """
    Implements the NetMHC binding (for netMHC3.0)::


    .. note::

        NetMHC-3.0: accurate web accessible predictions of human, mouse and monkey MHC class I affinities for peptides
        of length 8-11. Lundegaard C, Lamberth K, Harndahl M, Buus S, Lund O, Nielsen M.
        Nucleic Acids Res. 1;36(Web Server issue):W509-12. 2008

        Accurate approximation method for prediction of class I MHC affinities for peptides of length 8, 10 and 11
        using prediction tools trained on 9mers. Lundegaard C, Lund O, Nielsen M. Bioinformatics, 24(11):1397-98, 2008.
    """

    __alleles = frozenset(['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:02', 'HLA-A*02:03', 'HLA-A*02:04', 'HLA-A*02:06', 'HLA-A*02:11', 'HLA-A*02:12',
                           'HLA-A*02:16', 'HLA-A*02:19', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*24:03', 'HLA-A*26:01',
                           'HLA-A*26:02', 'HLA-A*29:02', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*33:01', 'HLA-A*68:01', 'HLA-A*68:02',
                           'HLA-A*69:01', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*08:02', 'HLA-B*15:01', 'HLA-B*18:01', 'HLA-B*27:05', 'HLA-B*35:01',
                           'HLA-B*39:01', 'HLA-B*40:01', 'HLA-B*40:02', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*45:01', 'HLA-B*51:01', 'HLA-B*53:01',
                           'HLA-B*54:01', 'HLA-B*57:01', 'HLA-B*58:01',
                           'H2-Db', 'H2-Dd', 'H2-Kb', 'H2-Kd', 'H2-Kk', 'H2-Ld'])  # no PSSM predictors

    __supported_length = frozenset([8, 9, 10, 11])
    __name = "netmhc"
    __version = "3.0a"
    __command = "netMHC-3.0 -p {peptides} -a {alleles} -x {out} -l {length} {options}"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    @property
    def command(self):
        """
        Defines the commandline call for external tool
        """
        return self.__command

    @property
    def supportedAlleles(self):
        """
        A list of valid :class:`~Fred2.Core.Allele.Allele` models
        """
        return self.__alleles

    @property
    def supportedLength(self):
        """
        A list of supported :class:`~Fred2.Core.Peptide.Peptide` lengths
        """
        return self.__supported_length

    def _represent(self, allele):
        """
        Internal function transforming an allele object into its representative string
        :param allele: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: str
        """
        if isinstance(allele, MouseAllele):
            return "H-2-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)
        else:
            return "HLA-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)

    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        result = defaultdict(dict)
        with open(file, 'r') as f:
            next(f, None)  # skip first line with logging stuff
            next(f, None)  # skip first line with nothing
            csvr = csv.reader(f, delimiter='\t')
            alleles = [x.split()[0] for x in csvr.next()[3:]]
            for l in csvr:
                if not l:
                    continue
                pep_seq = l[2]
                for ic_50, a in zip(l[3:], alleles):
                    sc = 1.0 - math.log(float(ic_50), 50000)
                    result[a][pep_seq] = sc if sc > 0.0 else 0.0
        if 'Average' in result:
            result.pop('Average')
        return dict(result)


class NetMHC_4_0(NetMHC_3_4):
    """
    Implements the NetMHC 4.0 binding

    .. note::
        Andreatta M, Nielsen M. Gapped sequence alignment using artificial neural networks:
        application to the MHC class I system. Bioinformatics (2016) Feb 15;32(4):511-7
    """
    __command = "netMHC -p {peptides} -a {alleles} -xls -xlsfile {out} {options}"
    __version = "4.0"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def command(self):
        """
        Defines the commandline call for external tool
        """
        return self.__command

    def _represent(self, allele):
        """
        Internal function transforming an allele object into its representative string
        :param allele: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: str
        """
        if isinstance(allele, MouseAllele):
            return "H-2-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)
        else:
            return "HLA-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)

    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        result = defaultdict(defaultdict)
        f = csv.reader(open(file, "r"), delimiter='\t')
        pos_factor = 3
        alleles = [x.split()[0] for x in [x for x in next(f) if x.strip() != ""]]
        next(f)
        for l in f:
            if not l:
                continue
            pep_seq = l[1]
            for i, a in enumerate(alleles):
                ic_50 = l[(i+1)*pos_factor]
                sc = 1.0 - math.log(float(ic_50), 50000)
                result[a][pep_seq] = sc if sc > 0.0 else 0.0
        return dict(result)

    def get_external_version(self, path=None):
        """
        Returns the external version of the tool by executing
        >{command} --version

        might be dependent on the method and has to be overwritten
        therefore it is declared abstract to enforce the user to
        overwrite the method. The function in the base class can be called
        with super()

        :param str path: Optional specification of executable path if deviant from :attr:`self.__command`
        :return: The external version of the tool or None if tool does not support versioning
        :rtype: str
        """
        # can not be determined netmhcpan does not support --version or similar
        return None


class NetMHCpan_2_4(AExternalEpitopePrediction):
    """
    Implements the NetMHC binding (in current form for netMHCpan 2.4).
    Supported  MHC alleles currently only restricted to HLA alleles.

    .. note::

        Nielsen, Morten, et al. "NetMHCpan, a method for quantitative predictions of peptide binding to any HLA-A and-B
        locus protein of known sequence." PloS one 2.8 (2007): e796.
    """
    __supported_length = frozenset([8, 9, 10, 11])
    __name = "netmhcpan"
    __command = "netMHCpan-2.4 -p {peptides} -a {alleles} {options} -ic50 -xls -xlsfile {out}"
    __alleles = frozenset(
        ['HLA-A*01:01', 'HLA-A*01:02', 'HLA-A*01:03', 'HLA-A*01:06', 'HLA-A*01:07', 'HLA-A*01:08', 'HLA-A*01:09', 'HLA-A*01:10', 'HLA-A*01:12',
         'HLA-A*01:13', 'HLA-A*01:14', 'HLA-A*01:17', 'HLA-A*01:19', 'HLA-A*01:20', 'HLA-A*01:21', 'HLA-A*01:23', 'HLA-A*01:24', 'HLA-A*01:25',
         'HLA-A*01:26', 'HLA-A*01:28', 'HLA-A*01:29', 'HLA-A*01:30', 'HLA-A*01:32', 'HLA-A*01:33', 'HLA-A*01:35', 'HLA-A*01:36', 'HLA-A*01:37',
         'HLA-A*01:38', 'HLA-A*01:39', 'HLA-A*01:40', 'HLA-A*01:41', 'HLA-A*01:42', 'HLA-A*01:43', 'HLA-A*01:44', 'HLA-A*01:45', 'HLA-A*01:46',
         'HLA-A*01:47', 'HLA-A*01:48', 'HLA-A*01:49', 'HLA-A*01:50', 'HLA-A*01:51', 'HLA-A*01:54', 'HLA-A*01:55', 'HLA-A*01:58', 'HLA-A*01:59',
         'HLA-A*01:60', 'HLA-A*01:61', 'HLA-A*01:62', 'HLA-A*01:63', 'HLA-A*01:64', 'HLA-A*01:65', 'HLA-A*01:66', 'HLA-A*02:01', 'HLA-A*02:02',
         'HLA-A*02:03', 'HLA-A*02:04', 'HLA-A*02:05', 'HLA-A*02:06', 'HLA-A*02:07', 'HLA-A*02:08', 'HLA-A*02:09', 'HLA-A*02:10', 'HLA-A*02:101',
         'HLA-A*02:102', 'HLA-A*02:103', 'HLA-A*02:104', 'HLA-A*02:105', 'HLA-A*02:106', 'HLA-A*02:107', 'HLA-A*02:108', 'HLA-A*02:109',
         'HLA-A*02:11', 'HLA-A*02:110', 'HLA-A*02:111', 'HLA-A*02:112', 'HLA-A*02:114', 'HLA-A*02:115', 'HLA-A*02:116', 'HLA-A*02:117',
         'HLA-A*02:118', 'HLA-A*02:119', 'HLA-A*02:12', 'HLA-A*02:120', 'HLA-A*02:121', 'HLA-A*02:122', 'HLA-A*02:123', 'HLA-A*02:124',
         'HLA-A*02:126', 'HLA-A*02:127', 'HLA-A*02:128', 'HLA-A*02:129', 'HLA-A*02:13', 'HLA-A*02:130', 'HLA-A*02:131', 'HLA-A*02:132',
         'HLA-A*02:133', 'HLA-A*02:134', 'HLA-A*02:135', 'HLA-A*02:136', 'HLA-A*02:137', 'HLA-A*02:138', 'HLA-A*02:139', 'HLA-A*02:14',
         'HLA-A*02:140', 'HLA-A*02:141', 'HLA-A*02:142', 'HLA-A*02:143', 'HLA-A*02:144', 'HLA-A*02:145', 'HLA-A*02:146', 'HLA-A*02:147',
         'HLA-A*02:148', 'HLA-A*02:149', 'HLA-A*02:150', 'HLA-A*02:151', 'HLA-A*02:152', 'HLA-A*02:153', 'HLA-A*02:154', 'HLA-A*02:155',
         'HLA-A*02:156', 'HLA-A*02:157', 'HLA-A*02:158', 'HLA-A*02:159', 'HLA-A*02:16', 'HLA-A*02:160', 'HLA-A*02:161', 'HLA-A*02:162',
         'HLA-A*02:163', 'HLA-A*02:164', 'HLA-A*02:165', 'HLA-A*02:166', 'HLA-A*02:167', 'HLA-A*02:168', 'HLA-A*02:169', 'HLA-A*02:17',
         'HLA-A*02:170', 'HLA-A*02:171', 'HLA-A*02:172', 'HLA-A*02:173', 'HLA-A*02:174', 'HLA-A*02:175', 'HLA-A*02:176', 'HLA-A*02:177',
         'HLA-A*02:178', 'HLA-A*02:179', 'HLA-A*02:18', 'HLA-A*02:180', 'HLA-A*02:181', 'HLA-A*02:182', 'HLA-A*02:183', 'HLA-A*02:184',
         'HLA-A*02:185', 'HLA-A*02:186', 'HLA-A*02:187', 'HLA-A*02:188', 'HLA-A*02:189', 'HLA-A*02:19', 'HLA-A*02:190', 'HLA-A*02:191',
         'HLA-A*02:192', 'HLA-A*02:193', 'HLA-A*02:194', 'HLA-A*02:195', 'HLA-A*02:196', 'HLA-A*02:197', 'HLA-A*02:198', 'HLA-A*02:199',
         'HLA-A*02:20', 'HLA-A*02:200', 'HLA-A*02:201', 'HLA-A*02:202', 'HLA-A*02:203', 'HLA-A*02:204', 'HLA-A*02:205', 'HLA-A*02:206',
         'HLA-A*02:207', 'HLA-A*02:208', 'HLA-A*02:209', 'HLA-A*02:21', 'HLA-A*02:210', 'HLA-A*02:211', 'HLA-A*02:212', 'HLA-A*02:213',
         'HLA-A*02:214', 'HLA-A*02:215', 'HLA-A*02:216', 'HLA-A*02:217', 'HLA-A*02:218', 'HLA-A*02:219', 'HLA-A*02:22', 'HLA-A*02:220',
         'HLA-A*02:221', 'HLA-A*02:224', 'HLA-A*02:228', 'HLA-A*02:229', 'HLA-A*02:230', 'HLA-A*02:231', 'HLA-A*02:232', 'HLA-A*02:233',
         'HLA-A*02:234', 'HLA-A*02:235', 'HLA-A*02:236', 'HLA-A*02:237', 'HLA-A*02:238', 'HLA-A*02:239', 'HLA-A*02:24', 'HLA-A*02:240',
         'HLA-A*02:241', 'HLA-A*02:242', 'HLA-A*02:243', 'HLA-A*02:244', 'HLA-A*02:245', 'HLA-A*02:246', 'HLA-A*02:247', 'HLA-A*02:248',
         'HLA-A*02:249', 'HLA-A*02:25', 'HLA-A*02:251', 'HLA-A*02:252', 'HLA-A*02:253', 'HLA-A*02:254', 'HLA-A*02:255', 'HLA-A*02:256',
         'HLA-A*02:257', 'HLA-A*02:258', 'HLA-A*02:259', 'HLA-A*02:26', 'HLA-A*02:260', 'HLA-A*02:261', 'HLA-A*02:262', 'HLA-A*02:263',
         'HLA-A*02:264', 'HLA-A*02:265', 'HLA-A*02:266', 'HLA-A*02:27', 'HLA-A*02:28', 'HLA-A*02:29', 'HLA-A*02:30', 'HLA-A*02:31', 'HLA-A*02:33',
         'HLA-A*02:34', 'HLA-A*02:35', 'HLA-A*02:36', 'HLA-A*02:37', 'HLA-A*02:38', 'HLA-A*02:39', 'HLA-A*02:40', 'HLA-A*02:41', 'HLA-A*02:42',
         'HLA-A*02:44', 'HLA-A*02:45', 'HLA-A*02:46', 'HLA-A*02:47', 'HLA-A*02:48', 'HLA-A*02:49', 'HLA-A*02:50', 'HLA-A*02:51', 'HLA-A*02:52',
         'HLA-A*02:54', 'HLA-A*02:55', 'HLA-A*02:56', 'HLA-A*02:57', 'HLA-A*02:58', 'HLA-A*02:59', 'HLA-A*02:60', 'HLA-A*02:61', 'HLA-A*02:62',
         'HLA-A*02:63', 'HLA-A*02:64', 'HLA-A*02:65', 'HLA-A*02:66', 'HLA-A*02:67', 'HLA-A*02:68', 'HLA-A*02:69', 'HLA-A*02:70', 'HLA-A*02:71',
         'HLA-A*02:72', 'HLA-A*02:73', 'HLA-A*02:74', 'HLA-A*02:75', 'HLA-A*02:76', 'HLA-A*02:77', 'HLA-A*02:78', 'HLA-A*02:79', 'HLA-A*02:80',
         'HLA-A*02:81', 'HLA-A*02:84', 'HLA-A*02:85', 'HLA-A*02:86', 'HLA-A*02:87', 'HLA-A*02:89', 'HLA-A*02:90', 'HLA-A*02:91', 'HLA-A*02:92',
         'HLA-A*02:93', 'HLA-A*02:95', 'HLA-A*02:96', 'HLA-A*02:97', 'HLA-A*02:99', 'HLA-A*03:01', 'HLA-A*03:02', 'HLA-A*03:04', 'HLA-A*03:05',
         'HLA-A*03:06', 'HLA-A*03:07', 'HLA-A*03:08', 'HLA-A*03:09', 'HLA-A*03:10', 'HLA-A*03:12', 'HLA-A*03:13', 'HLA-A*03:14', 'HLA-A*03:15',
         'HLA-A*03:16', 'HLA-A*03:17', 'HLA-A*03:18', 'HLA-A*03:19', 'HLA-A*03:20', 'HLA-A*03:22', 'HLA-A*03:23', 'HLA-A*03:24', 'HLA-A*03:25',
         'HLA-A*03:26', 'HLA-A*03:27', 'HLA-A*03:28', 'HLA-A*03:29', 'HLA-A*03:30', 'HLA-A*03:31', 'HLA-A*03:32', 'HLA-A*03:33', 'HLA-A*03:34',
         'HLA-A*03:35', 'HLA-A*03:37', 'HLA-A*03:38', 'HLA-A*03:39', 'HLA-A*03:40', 'HLA-A*03:41', 'HLA-A*03:42', 'HLA-A*03:43', 'HLA-A*03:44',
         'HLA-A*03:45', 'HLA-A*03:46', 'HLA-A*03:47', 'HLA-A*03:48', 'HLA-A*03:49', 'HLA-A*03:50', 'HLA-A*03:51', 'HLA-A*03:52', 'HLA-A*03:53',
         'HLA-A*03:54', 'HLA-A*03:55', 'HLA-A*03:56', 'HLA-A*03:57', 'HLA-A*03:58', 'HLA-A*03:59', 'HLA-A*03:60', 'HLA-A*03:61', 'HLA-A*03:62',
         'HLA-A*03:63', 'HLA-A*03:64', 'HLA-A*03:65', 'HLA-A*03:66', 'HLA-A*03:67', 'HLA-A*03:70', 'HLA-A*03:71', 'HLA-A*03:72', 'HLA-A*03:73',
         'HLA-A*03:74', 'HLA-A*03:75', 'HLA-A*03:76', 'HLA-A*03:77', 'HLA-A*03:78', 'HLA-A*03:79', 'HLA-A*03:80', 'HLA-A*03:81', 'HLA-A*03:82',
         'HLA-A*11:01', 'HLA-A*11:02', 'HLA-A*11:03', 'HLA-A*11:04', 'HLA-A*11:05', 'HLA-A*11:06', 'HLA-A*11:07', 'HLA-A*11:08', 'HLA-A*11:09',
         'HLA-A*11:10', 'HLA-A*11:11', 'HLA-A*11:12', 'HLA-A*11:13', 'HLA-A*11:14', 'HLA-A*11:15', 'HLA-A*11:16', 'HLA-A*11:17', 'HLA-A*11:18',
         'HLA-A*11:19', 'HLA-A*11:20', 'HLA-A*11:22', 'HLA-A*11:23', 'HLA-A*11:24', 'HLA-A*11:25', 'HLA-A*11:26', 'HLA-A*11:27', 'HLA-A*11:29',
         'HLA-A*11:30', 'HLA-A*11:31', 'HLA-A*11:32', 'HLA-A*11:33', 'HLA-A*11:34', 'HLA-A*11:35', 'HLA-A*11:36', 'HLA-A*11:37', 'HLA-A*11:38',
         'HLA-A*11:39', 'HLA-A*11:40', 'HLA-A*11:41', 'HLA-A*11:42', 'HLA-A*11:43', 'HLA-A*11:44', 'HLA-A*11:45', 'HLA-A*11:46', 'HLA-A*11:47',
         'HLA-A*11:48', 'HLA-A*11:49', 'HLA-A*11:51', 'HLA-A*11:53', 'HLA-A*11:54', 'HLA-A*11:55', 'HLA-A*11:56', 'HLA-A*11:57', 'HLA-A*11:58',
         'HLA-A*11:59', 'HLA-A*11:60', 'HLA-A*11:61', 'HLA-A*11:62', 'HLA-A*11:63', 'HLA-A*11:64', 'HLA-A*23:01', 'HLA-A*23:02', 'HLA-A*23:03',
         'HLA-A*23:04', 'HLA-A*23:05', 'HLA-A*23:06', 'HLA-A*23:09', 'HLA-A*23:10', 'HLA-A*23:12', 'HLA-A*23:13', 'HLA-A*23:14', 'HLA-A*23:15',
         'HLA-A*23:16', 'HLA-A*23:17', 'HLA-A*23:18', 'HLA-A*23:20', 'HLA-A*23:21', 'HLA-A*23:22', 'HLA-A*23:23', 'HLA-A*23:24', 'HLA-A*23:25',
         'HLA-A*23:26', 'HLA-A*24:02', 'HLA-A*24:03', 'HLA-A*24:04', 'HLA-A*24:05', 'HLA-A*24:06', 'HLA-A*24:07', 'HLA-A*24:08', 'HLA-A*24:10',
         'HLA-A*24:100', 'HLA-A*24:101', 'HLA-A*24:102', 'HLA-A*24:103', 'HLA-A*24:104', 'HLA-A*24:105', 'HLA-A*24:106', 'HLA-A*24:107',
         'HLA-A*24:108', 'HLA-A*24:109', 'HLA-A*24:110', 'HLA-A*24:111', 'HLA-A*24:112', 'HLA-A*24:113', 'HLA-A*24:114', 'HLA-A*24:115',
         'HLA-A*24:116', 'HLA-A*24:117', 'HLA-A*24:118', 'HLA-A*24:119', 'HLA-A*24:120', 'HLA-A*24:121', 'HLA-A*24:122', 'HLA-A*24:123',
         'HLA-A*24:124', 'HLA-A*24:125', 'HLA-A*24:126', 'HLA-A*24:127', 'HLA-A*24:128', 'HLA-A*24:129', 'HLA-A*24:13', 'HLA-A*24:130',
         'HLA-A*24:131', 'HLA-A*24:133', 'HLA-A*24:134', 'HLA-A*24:135', 'HLA-A*24:136', 'HLA-A*24:137', 'HLA-A*24:138', 'HLA-A*24:139',
         'HLA-A*24:14', 'HLA-A*24:140', 'HLA-A*24:141', 'HLA-A*24:142', 'HLA-A*24:143', 'HLA-A*24:144', 'HLA-A*24:15', 'HLA-A*24:17', 'HLA-A*24:18',
         'HLA-A*24:19', 'HLA-A*24:20', 'HLA-A*24:21', 'HLA-A*24:22', 'HLA-A*24:23', 'HLA-A*24:24', 'HLA-A*24:25', 'HLA-A*24:26', 'HLA-A*24:27',
         'HLA-A*24:28', 'HLA-A*24:29', 'HLA-A*24:30', 'HLA-A*24:31', 'HLA-A*24:32', 'HLA-A*24:33', 'HLA-A*24:34', 'HLA-A*24:35', 'HLA-A*24:37',
         'HLA-A*24:38', 'HLA-A*24:39', 'HLA-A*24:41', 'HLA-A*24:42', 'HLA-A*24:43', 'HLA-A*24:44', 'HLA-A*24:46', 'HLA-A*24:47', 'HLA-A*24:49',
         'HLA-A*24:50', 'HLA-A*24:51', 'HLA-A*24:52', 'HLA-A*24:53', 'HLA-A*24:54', 'HLA-A*24:55', 'HLA-A*24:56', 'HLA-A*24:57', 'HLA-A*24:58',
         'HLA-A*24:59', 'HLA-A*24:61', 'HLA-A*24:62', 'HLA-A*24:63', 'HLA-A*24:64', 'HLA-A*24:66', 'HLA-A*24:67', 'HLA-A*24:68', 'HLA-A*24:69',
         'HLA-A*24:70', 'HLA-A*24:71', 'HLA-A*24:72', 'HLA-A*24:73', 'HLA-A*24:74', 'HLA-A*24:75', 'HLA-A*24:76', 'HLA-A*24:77', 'HLA-A*24:78',
         'HLA-A*24:79', 'HLA-A*24:80', 'HLA-A*24:81', 'HLA-A*24:82', 'HLA-A*24:85', 'HLA-A*24:87', 'HLA-A*24:88', 'HLA-A*24:89', 'HLA-A*24:91',
         'HLA-A*24:92', 'HLA-A*24:93', 'HLA-A*24:94', 'HLA-A*24:95', 'HLA-A*24:96', 'HLA-A*24:97', 'HLA-A*24:98', 'HLA-A*24:99', 'HLA-A*25:01',
         'HLA-A*25:02', 'HLA-A*25:03', 'HLA-A*25:04', 'HLA-A*25:05', 'HLA-A*25:06', 'HLA-A*25:07', 'HLA-A*25:08', 'HLA-A*25:09', 'HLA-A*25:10',
         'HLA-A*25:11', 'HLA-A*25:13', 'HLA-A*26:01', 'HLA-A*26:02', 'HLA-A*26:03', 'HLA-A*26:04', 'HLA-A*26:05', 'HLA-A*26:06', 'HLA-A*26:07',
         'HLA-A*26:08', 'HLA-A*26:09', 'HLA-A*26:10', 'HLA-A*26:12', 'HLA-A*26:13', 'HLA-A*26:14', 'HLA-A*26:15', 'HLA-A*26:16', 'HLA-A*26:17',
         'HLA-A*26:18', 'HLA-A*26:19', 'HLA-A*26:20', 'HLA-A*26:21', 'HLA-A*26:22', 'HLA-A*26:23', 'HLA-A*26:24', 'HLA-A*26:26', 'HLA-A*26:27',
         'HLA-A*26:28', 'HLA-A*26:29', 'HLA-A*26:30', 'HLA-A*26:31', 'HLA-A*26:32', 'HLA-A*26:33', 'HLA-A*26:34', 'HLA-A*26:35', 'HLA-A*26:36',
         'HLA-A*26:37', 'HLA-A*26:38', 'HLA-A*26:39', 'HLA-A*26:40', 'HLA-A*26:41', 'HLA-A*26:42', 'HLA-A*26:43', 'HLA-A*26:45', 'HLA-A*26:46',
         'HLA-A*26:47', 'HLA-A*26:48', 'HLA-A*26:49', 'HLA-A*26:50', 'HLA-A*29:01', 'HLA-A*29:02', 'HLA-A*29:03', 'HLA-A*29:04', 'HLA-A*29:05',
         'HLA-A*29:06', 'HLA-A*29:07', 'HLA-A*29:09', 'HLA-A*29:10', 'HLA-A*29:11', 'HLA-A*29:12', 'HLA-A*29:13', 'HLA-A*29:14', 'HLA-A*29:15',
         'HLA-A*29:16', 'HLA-A*29:17', 'HLA-A*29:18', 'HLA-A*29:19', 'HLA-A*29:20', 'HLA-A*29:21', 'HLA-A*29:22', 'HLA-A*30:01', 'HLA-A*30:02',
         'HLA-A*30:03', 'HLA-A*30:04', 'HLA-A*30:06', 'HLA-A*30:07', 'HLA-A*30:08', 'HLA-A*30:09', 'HLA-A*30:10', 'HLA-A*30:11', 'HLA-A*30:12',
         'HLA-A*30:13', 'HLA-A*30:15', 'HLA-A*30:16', 'HLA-A*30:17', 'HLA-A*30:18', 'HLA-A*30:19', 'HLA-A*30:20', 'HLA-A*30:22', 'HLA-A*30:23',
         'HLA-A*30:24', 'HLA-A*30:25', 'HLA-A*30:26', 'HLA-A*30:28', 'HLA-A*30:29', 'HLA-A*30:30', 'HLA-A*30:31', 'HLA-A*30:32', 'HLA-A*30:33',
         'HLA-A*30:34', 'HLA-A*30:35', 'HLA-A*30:36', 'HLA-A*30:37', 'HLA-A*30:38', 'HLA-A*30:39', 'HLA-A*30:40', 'HLA-A*30:41', 'HLA-A*31:01',
         'HLA-A*31:02', 'HLA-A*31:03', 'HLA-A*31:04', 'HLA-A*31:05', 'HLA-A*31:06', 'HLA-A*31:07', 'HLA-A*31:08', 'HLA-A*31:09', 'HLA-A*31:10',
         'HLA-A*31:11', 'HLA-A*31:12', 'HLA-A*31:13', 'HLA-A*31:15', 'HLA-A*31:16', 'HLA-A*31:17', 'HLA-A*31:18', 'HLA-A*31:19', 'HLA-A*31:20',
         'HLA-A*31:21', 'HLA-A*31:22', 'HLA-A*31:23', 'HLA-A*31:24', 'HLA-A*31:25', 'HLA-A*31:26', 'HLA-A*31:27', 'HLA-A*31:28', 'HLA-A*31:29',
         'HLA-A*31:30', 'HLA-A*31:31', 'HLA-A*31:32', 'HLA-A*31:33', 'HLA-A*31:34', 'HLA-A*31:35', 'HLA-A*31:36', 'HLA-A*31:37', 'HLA-A*32:01',
         'HLA-A*32:02', 'HLA-A*32:03', 'HLA-A*32:04', 'HLA-A*32:05', 'HLA-A*32:06', 'HLA-A*32:07', 'HLA-A*32:08', 'HLA-A*32:09', 'HLA-A*32:10',
         'HLA-A*32:12', 'HLA-A*32:13', 'HLA-A*32:14', 'HLA-A*32:15', 'HLA-A*32:16', 'HLA-A*32:17', 'HLA-A*32:18', 'HLA-A*32:20', 'HLA-A*32:21',
         'HLA-A*32:22', 'HLA-A*32:23', 'HLA-A*32:24', 'HLA-A*32:25', 'HLA-A*33:01', 'HLA-A*33:03', 'HLA-A*33:04', 'HLA-A*33:05', 'HLA-A*33:06',
         'HLA-A*33:07', 'HLA-A*33:08', 'HLA-A*33:09', 'HLA-A*33:10', 'HLA-A*33:11', 'HLA-A*33:12', 'HLA-A*33:13', 'HLA-A*33:14', 'HLA-A*33:15',
         'HLA-A*33:16', 'HLA-A*33:17', 'HLA-A*33:18', 'HLA-A*33:19', 'HLA-A*33:20', 'HLA-A*33:21', 'HLA-A*33:22', 'HLA-A*33:23', 'HLA-A*33:24',
         'HLA-A*33:25', 'HLA-A*33:26', 'HLA-A*33:27', 'HLA-A*33:28', 'HLA-A*33:29', 'HLA-A*33:30', 'HLA-A*33:31', 'HLA-A*34:01', 'HLA-A*34:02',
         'HLA-A*34:03', 'HLA-A*34:04', 'HLA-A*34:05', 'HLA-A*34:06', 'HLA-A*34:07', 'HLA-A*34:08', 'HLA-A*36:01', 'HLA-A*36:02', 'HLA-A*36:03',
         'HLA-A*36:04', 'HLA-A*36:05', 'HLA-A*43:01', 'HLA-A*66:01', 'HLA-A*66:02', 'HLA-A*66:03', 'HLA-A*66:04', 'HLA-A*66:05', 'HLA-A*66:06',
         'HLA-A*66:07', 'HLA-A*66:08', 'HLA-A*66:09', 'HLA-A*66:10', 'HLA-A*66:11', 'HLA-A*66:12', 'HLA-A*66:13', 'HLA-A*66:14', 'HLA-A*66:15',
         'HLA-A*68:01', 'HLA-A*68:02', 'HLA-A*68:03', 'HLA-A*68:04', 'HLA-A*68:05', 'HLA-A*68:06', 'HLA-A*68:07', 'HLA-A*68:08', 'HLA-A*68:09',
         'HLA-A*68:10', 'HLA-A*68:12', 'HLA-A*68:13', 'HLA-A*68:14', 'HLA-A*68:15', 'HLA-A*68:16', 'HLA-A*68:17', 'HLA-A*68:19', 'HLA-A*68:20',
         'HLA-A*68:21', 'HLA-A*68:22', 'HLA-A*68:23', 'HLA-A*68:24', 'HLA-A*68:25', 'HLA-A*68:26', 'HLA-A*68:27', 'HLA-A*68:28', 'HLA-A*68:29',
         'HLA-A*68:30', 'HLA-A*68:31', 'HLA-A*68:32', 'HLA-A*68:33', 'HLA-A*68:34', 'HLA-A*68:35', 'HLA-A*68:36', 'HLA-A*68:37', 'HLA-A*68:38',
         'HLA-A*68:39', 'HLA-A*68:40', 'HLA-A*68:41', 'HLA-A*68:42', 'HLA-A*68:43', 'HLA-A*68:44', 'HLA-A*68:45', 'HLA-A*68:46', 'HLA-A*68:47',
         'HLA-A*68:48', 'HLA-A*68:50', 'HLA-A*68:51', 'HLA-A*68:52', 'HLA-A*68:53', 'HLA-A*68:54', 'HLA-A*69:01', 'HLA-A*74:01', 'HLA-A*74:02',
         'HLA-A*74:03', 'HLA-A*74:04', 'HLA-A*74:05', 'HLA-A*74:06', 'HLA-A*74:07', 'HLA-A*74:08', 'HLA-A*74:09', 'HLA-A*74:10', 'HLA-A*74:11',
         'HLA-A*74:13', 'HLA-A*80:01', 'HLA-A*80:02', 'HLA-B*07:02', 'HLA-B*07:03', 'HLA-B*07:04', 'HLA-B*07:05', 'HLA-B*07:06', 'HLA-B*07:07',
         'HLA-B*07:08', 'HLA-B*07:09', 'HLA-B*07:10', 'HLA-B*07:100', 'HLA-B*07:101', 'HLA-B*07:102', 'HLA-B*07:103', 'HLA-B*07:104',
         'HLA-B*07:105', 'HLA-B*07:106', 'HLA-B*07:107', 'HLA-B*07:108', 'HLA-B*07:109', 'HLA-B*07:11', 'HLA-B*07:110', 'HLA-B*07:112',
         'HLA-B*07:113', 'HLA-B*07:114', 'HLA-B*07:115', 'HLA-B*07:12', 'HLA-B*07:13', 'HLA-B*07:14', 'HLA-B*07:15', 'HLA-B*07:16', 'HLA-B*07:17',
         'HLA-B*07:18', 'HLA-B*07:19', 'HLA-B*07:20', 'HLA-B*07:21', 'HLA-B*07:22', 'HLA-B*07:23', 'HLA-B*07:24', 'HLA-B*07:25', 'HLA-B*07:26',
         'HLA-B*07:27', 'HLA-B*07:28', 'HLA-B*07:29', 'HLA-B*07:30', 'HLA-B*07:31', 'HLA-B*07:32', 'HLA-B*07:33', 'HLA-B*07:34', 'HLA-B*07:35',
         'HLA-B*07:36', 'HLA-B*07:37', 'HLA-B*07:38', 'HLA-B*07:39', 'HLA-B*07:40', 'HLA-B*07:41', 'HLA-B*07:42', 'HLA-B*07:43', 'HLA-B*07:44',
         'HLA-B*07:45', 'HLA-B*07:46', 'HLA-B*07:47', 'HLA-B*07:48', 'HLA-B*07:50', 'HLA-B*07:51', 'HLA-B*07:52', 'HLA-B*07:53', 'HLA-B*07:54',
         'HLA-B*07:55', 'HLA-B*07:56', 'HLA-B*07:57', 'HLA-B*07:58', 'HLA-B*07:59', 'HLA-B*07:60', 'HLA-B*07:61', 'HLA-B*07:62', 'HLA-B*07:63',
         'HLA-B*07:64', 'HLA-B*07:65', 'HLA-B*07:66', 'HLA-B*07:68', 'HLA-B*07:69', 'HLA-B*07:70', 'HLA-B*07:71', 'HLA-B*07:72', 'HLA-B*07:73',
         'HLA-B*07:74', 'HLA-B*07:75', 'HLA-B*07:76', 'HLA-B*07:77', 'HLA-B*07:78', 'HLA-B*07:79', 'HLA-B*07:80', 'HLA-B*07:81', 'HLA-B*07:82',
         'HLA-B*07:83', 'HLA-B*07:84', 'HLA-B*07:85', 'HLA-B*07:86', 'HLA-B*07:87', 'HLA-B*07:88', 'HLA-B*07:89', 'HLA-B*07:90', 'HLA-B*07:91',
         'HLA-B*07:92', 'HLA-B*07:93', 'HLA-B*07:94', 'HLA-B*07:95', 'HLA-B*07:96', 'HLA-B*07:97', 'HLA-B*07:98', 'HLA-B*07:99', 'HLA-B*08:01',
         'HLA-B*08:02', 'HLA-B*08:03', 'HLA-B*08:04', 'HLA-B*08:05', 'HLA-B*08:07', 'HLA-B*08:09', 'HLA-B*08:10', 'HLA-B*08:11', 'HLA-B*08:12',
         'HLA-B*08:13', 'HLA-B*08:14', 'HLA-B*08:15', 'HLA-B*08:16', 'HLA-B*08:17', 'HLA-B*08:18', 'HLA-B*08:20', 'HLA-B*08:21', 'HLA-B*08:22',
         'HLA-B*08:23', 'HLA-B*08:24', 'HLA-B*08:25', 'HLA-B*08:26', 'HLA-B*08:27', 'HLA-B*08:28', 'HLA-B*08:29', 'HLA-B*08:31', 'HLA-B*08:32',
         'HLA-B*08:33', 'HLA-B*08:34', 'HLA-B*08:35', 'HLA-B*08:36', 'HLA-B*08:37', 'HLA-B*08:38', 'HLA-B*08:39', 'HLA-B*08:40', 'HLA-B*08:41',
         'HLA-B*08:42', 'HLA-B*08:43', 'HLA-B*08:44', 'HLA-B*08:45', 'HLA-B*08:46', 'HLA-B*08:47', 'HLA-B*08:48', 'HLA-B*08:49', 'HLA-B*08:50',
         'HLA-B*08:51', 'HLA-B*08:52', 'HLA-B*08:53', 'HLA-B*08:54', 'HLA-B*08:55', 'HLA-B*08:56', 'HLA-B*08:57', 'HLA-B*08:58', 'HLA-B*08:59',
         'HLA-B*08:60', 'HLA-B*08:61', 'HLA-B*08:62', 'HLA-B*13:01', 'HLA-B*13:02', 'HLA-B*13:03', 'HLA-B*13:04', 'HLA-B*13:06', 'HLA-B*13:09',
         'HLA-B*13:10', 'HLA-B*13:11', 'HLA-B*13:12', 'HLA-B*13:13', 'HLA-B*13:14', 'HLA-B*13:15', 'HLA-B*13:16', 'HLA-B*13:17', 'HLA-B*13:18',
         'HLA-B*13:19', 'HLA-B*13:20', 'HLA-B*13:21', 'HLA-B*13:22', 'HLA-B*13:23', 'HLA-B*13:25', 'HLA-B*13:26', 'HLA-B*13:27', 'HLA-B*13:28',
         'HLA-B*13:29', 'HLA-B*13:30', 'HLA-B*13:31', 'HLA-B*13:32', 'HLA-B*13:33', 'HLA-B*13:34', 'HLA-B*13:35', 'HLA-B*13:36', 'HLA-B*13:37',
         'HLA-B*13:38', 'HLA-B*13:39', 'HLA-B*14:01', 'HLA-B*14:02', 'HLA-B*14:03', 'HLA-B*14:04', 'HLA-B*14:05', 'HLA-B*14:06', 'HLA-B*14:08',
         'HLA-B*14:09', 'HLA-B*14:10', 'HLA-B*14:11', 'HLA-B*14:12', 'HLA-B*14:13', 'HLA-B*14:14', 'HLA-B*14:15', 'HLA-B*14:16', 'HLA-B*14:17',
         'HLA-B*14:18', 'HLA-B*15:01', 'HLA-B*15:02', 'HLA-B*15:03', 'HLA-B*15:04', 'HLA-B*15:05', 'HLA-B*15:06', 'HLA-B*15:07', 'HLA-B*15:08',
         'HLA-B*15:09', 'HLA-B*15:10', 'HLA-B*15:101', 'HLA-B*15:102', 'HLA-B*15:103', 'HLA-B*15:104', 'HLA-B*15:105', 'HLA-B*15:106',
         'HLA-B*15:107', 'HLA-B*15:108', 'HLA-B*15:109', 'HLA-B*15:11', 'HLA-B*15:110', 'HLA-B*15:112', 'HLA-B*15:113', 'HLA-B*15:114',
         'HLA-B*15:115', 'HLA-B*15:116', 'HLA-B*15:117', 'HLA-B*15:118', 'HLA-B*15:119', 'HLA-B*15:12', 'HLA-B*15:120', 'HLA-B*15:121',
         'HLA-B*15:122', 'HLA-B*15:123', 'HLA-B*15:124', 'HLA-B*15:125', 'HLA-B*15:126', 'HLA-B*15:127', 'HLA-B*15:128', 'HLA-B*15:129',
         'HLA-B*15:13', 'HLA-B*15:131', 'HLA-B*15:132', 'HLA-B*15:133', 'HLA-B*15:134', 'HLA-B*15:135', 'HLA-B*15:136', 'HLA-B*15:137',
         'HLA-B*15:138', 'HLA-B*15:139', 'HLA-B*15:14', 'HLA-B*15:140', 'HLA-B*15:141', 'HLA-B*15:142', 'HLA-B*15:143', 'HLA-B*15:144',
         'HLA-B*15:145', 'HLA-B*15:146', 'HLA-B*15:147', 'HLA-B*15:148', 'HLA-B*15:15', 'HLA-B*15:150', 'HLA-B*15:151', 'HLA-B*15:152',
         'HLA-B*15:153', 'HLA-B*15:154', 'HLA-B*15:155', 'HLA-B*15:156', 'HLA-B*15:157', 'HLA-B*15:158', 'HLA-B*15:159', 'HLA-B*15:16',
         'HLA-B*15:160', 'HLA-B*15:161', 'HLA-B*15:162', 'HLA-B*15:163', 'HLA-B*15:164', 'HLA-B*15:165', 'HLA-B*15:166', 'HLA-B*15:167',
         'HLA-B*15:168', 'HLA-B*15:169', 'HLA-B*15:17', 'HLA-B*15:170', 'HLA-B*15:171', 'HLA-B*15:172', 'HLA-B*15:173', 'HLA-B*15:174',
         'HLA-B*15:175', 'HLA-B*15:176', 'HLA-B*15:177', 'HLA-B*15:178', 'HLA-B*15:179', 'HLA-B*15:18', 'HLA-B*15:180', 'HLA-B*15:183',
         'HLA-B*15:184', 'HLA-B*15:185', 'HLA-B*15:186', 'HLA-B*15:187', 'HLA-B*15:188', 'HLA-B*15:189', 'HLA-B*15:19', 'HLA-B*15:191',
         'HLA-B*15:192', 'HLA-B*15:193', 'HLA-B*15:194', 'HLA-B*15:195', 'HLA-B*15:196', 'HLA-B*15:197', 'HLA-B*15:198', 'HLA-B*15:199',
         'HLA-B*15:20', 'HLA-B*15:200', 'HLA-B*15:201', 'HLA-B*15:202', 'HLA-B*15:21', 'HLA-B*15:23', 'HLA-B*15:24', 'HLA-B*15:25', 'HLA-B*15:27',
         'HLA-B*15:28', 'HLA-B*15:29', 'HLA-B*15:30', 'HLA-B*15:31', 'HLA-B*15:32', 'HLA-B*15:33', 'HLA-B*15:34', 'HLA-B*15:35', 'HLA-B*15:36',
         'HLA-B*15:37', 'HLA-B*15:38', 'HLA-B*15:39', 'HLA-B*15:40', 'HLA-B*15:42', 'HLA-B*15:43', 'HLA-B*15:44', 'HLA-B*15:45', 'HLA-B*15:46',
         'HLA-B*15:47', 'HLA-B*15:48', 'HLA-B*15:49', 'HLA-B*15:50', 'HLA-B*15:51', 'HLA-B*15:52', 'HLA-B*15:53', 'HLA-B*15:54', 'HLA-B*15:55',
         'HLA-B*15:56', 'HLA-B*15:57', 'HLA-B*15:58', 'HLA-B*15:60', 'HLA-B*15:61', 'HLA-B*15:62', 'HLA-B*15:63', 'HLA-B*15:64', 'HLA-B*15:65',
         'HLA-B*15:66', 'HLA-B*15:67', 'HLA-B*15:68', 'HLA-B*15:69', 'HLA-B*15:70', 'HLA-B*15:71', 'HLA-B*15:72', 'HLA-B*15:73', 'HLA-B*15:74',
         'HLA-B*15:75', 'HLA-B*15:76', 'HLA-B*15:77', 'HLA-B*15:78', 'HLA-B*15:80', 'HLA-B*15:81', 'HLA-B*15:82', 'HLA-B*15:83', 'HLA-B*15:84',
         'HLA-B*15:85', 'HLA-B*15:86', 'HLA-B*15:87', 'HLA-B*15:88', 'HLA-B*15:89', 'HLA-B*15:90', 'HLA-B*15:91', 'HLA-B*15:92', 'HLA-B*15:93',
         'HLA-B*15:95', 'HLA-B*15:96', 'HLA-B*15:97', 'HLA-B*15:98', 'HLA-B*15:99', 'HLA-B*18:01', 'HLA-B*18:02', 'HLA-B*18:03', 'HLA-B*18:04',
         'HLA-B*18:05', 'HLA-B*18:06', 'HLA-B*18:07', 'HLA-B*18:08', 'HLA-B*18:09', 'HLA-B*18:10', 'HLA-B*18:11', 'HLA-B*18:12', 'HLA-B*18:13',
         'HLA-B*18:14', 'HLA-B*18:15', 'HLA-B*18:18', 'HLA-B*18:19', 'HLA-B*18:20', 'HLA-B*18:21', 'HLA-B*18:22', 'HLA-B*18:24', 'HLA-B*18:25',
         'HLA-B*18:26', 'HLA-B*18:27', 'HLA-B*18:28', 'HLA-B*18:29', 'HLA-B*18:30', 'HLA-B*18:31', 'HLA-B*18:32', 'HLA-B*18:33', 'HLA-B*18:34',
         'HLA-B*18:35', 'HLA-B*18:36', 'HLA-B*18:37', 'HLA-B*18:38', 'HLA-B*18:39', 'HLA-B*18:40', 'HLA-B*18:41', 'HLA-B*18:42', 'HLA-B*18:43',
         'HLA-B*18:44', 'HLA-B*18:45', 'HLA-B*18:46', 'HLA-B*18:47', 'HLA-B*18:48', 'HLA-B*18:49', 'HLA-B*18:50', 'HLA-B*27:01', 'HLA-B*27:02',
         'HLA-B*27:03', 'HLA-B*27:04', 'HLA-B*27:05', 'HLA-B*27:06', 'HLA-B*27:07', 'HLA-B*27:08', 'HLA-B*27:09', 'HLA-B*27:10', 'HLA-B*27:11',
         'HLA-B*27:12', 'HLA-B*27:13', 'HLA-B*27:14', 'HLA-B*27:15', 'HLA-B*27:16', 'HLA-B*27:17', 'HLA-B*27:18', 'HLA-B*27:19', 'HLA-B*27:20',
         'HLA-B*27:21', 'HLA-B*27:23', 'HLA-B*27:24', 'HLA-B*27:25', 'HLA-B*27:26', 'HLA-B*27:27', 'HLA-B*27:28', 'HLA-B*27:29', 'HLA-B*27:30',
         'HLA-B*27:31', 'HLA-B*27:32', 'HLA-B*27:33', 'HLA-B*27:34', 'HLA-B*27:35', 'HLA-B*27:36', 'HLA-B*27:37', 'HLA-B*27:38', 'HLA-B*27:39',
         'HLA-B*27:40', 'HLA-B*27:41', 'HLA-B*27:42', 'HLA-B*27:43', 'HLA-B*27:44', 'HLA-B*27:45', 'HLA-B*27:46', 'HLA-B*27:47', 'HLA-B*27:48',
         'HLA-B*27:49', 'HLA-B*27:50', 'HLA-B*27:51', 'HLA-B*27:52', 'HLA-B*27:53', 'HLA-B*27:54', 'HLA-B*27:55', 'HLA-B*27:56', 'HLA-B*27:57',
         'HLA-B*27:58', 'HLA-B*27:60', 'HLA-B*27:61', 'HLA-B*27:62', 'HLA-B*27:63', 'HLA-B*27:67', 'HLA-B*27:68', 'HLA-B*27:69', 'HLA-B*35:01',
         'HLA-B*35:02', 'HLA-B*35:03', 'HLA-B*35:04', 'HLA-B*35:05', 'HLA-B*35:06', 'HLA-B*35:07', 'HLA-B*35:08', 'HLA-B*35:09', 'HLA-B*35:10',
         'HLA-B*35:100', 'HLA-B*35:101', 'HLA-B*35:102', 'HLA-B*35:103', 'HLA-B*35:104', 'HLA-B*35:105', 'HLA-B*35:106', 'HLA-B*35:107',
         'HLA-B*35:108', 'HLA-B*35:109', 'HLA-B*35:11', 'HLA-B*35:110', 'HLA-B*35:111', 'HLA-B*35:112', 'HLA-B*35:113', 'HLA-B*35:114',
         'HLA-B*35:115', 'HLA-B*35:116', 'HLA-B*35:117', 'HLA-B*35:118', 'HLA-B*35:119', 'HLA-B*35:12', 'HLA-B*35:120', 'HLA-B*35:121',
         'HLA-B*35:122', 'HLA-B*35:123', 'HLA-B*35:124', 'HLA-B*35:125', 'HLA-B*35:126', 'HLA-B*35:127', 'HLA-B*35:128', 'HLA-B*35:13',
         'HLA-B*35:131', 'HLA-B*35:132', 'HLA-B*35:133', 'HLA-B*35:135', 'HLA-B*35:136', 'HLA-B*35:137', 'HLA-B*35:138', 'HLA-B*35:139',
         'HLA-B*35:14', 'HLA-B*35:140', 'HLA-B*35:141', 'HLA-B*35:142', 'HLA-B*35:143', 'HLA-B*35:144', 'HLA-B*35:15', 'HLA-B*35:16', 'HLA-B*35:17',
         'HLA-B*35:18', 'HLA-B*35:19', 'HLA-B*35:20', 'HLA-B*35:21', 'HLA-B*35:22', 'HLA-B*35:23', 'HLA-B*35:24', 'HLA-B*35:25', 'HLA-B*35:26',
         'HLA-B*35:27', 'HLA-B*35:28', 'HLA-B*35:29', 'HLA-B*35:30', 'HLA-B*35:31', 'HLA-B*35:32', 'HLA-B*35:33', 'HLA-B*35:34', 'HLA-B*35:35',
         'HLA-B*35:36', 'HLA-B*35:37', 'HLA-B*35:38', 'HLA-B*35:39', 'HLA-B*35:41', 'HLA-B*35:42', 'HLA-B*35:43', 'HLA-B*35:44', 'HLA-B*35:45',
         'HLA-B*35:46', 'HLA-B*35:47', 'HLA-B*35:48', 'HLA-B*35:49', 'HLA-B*35:50', 'HLA-B*35:51', 'HLA-B*35:52', 'HLA-B*35:54', 'HLA-B*35:55',
         'HLA-B*35:56', 'HLA-B*35:57', 'HLA-B*35:58', 'HLA-B*35:59', 'HLA-B*35:60', 'HLA-B*35:61', 'HLA-B*35:62', 'HLA-B*35:63', 'HLA-B*35:64',
         'HLA-B*35:66', 'HLA-B*35:67', 'HLA-B*35:68', 'HLA-B*35:69', 'HLA-B*35:70', 'HLA-B*35:71', 'HLA-B*35:72', 'HLA-B*35:74', 'HLA-B*35:75',
         'HLA-B*35:76', 'HLA-B*35:77', 'HLA-B*35:78', 'HLA-B*35:79', 'HLA-B*35:80', 'HLA-B*35:81', 'HLA-B*35:82', 'HLA-B*35:83', 'HLA-B*35:84',
         'HLA-B*35:85', 'HLA-B*35:86', 'HLA-B*35:87', 'HLA-B*35:88', 'HLA-B*35:89', 'HLA-B*35:90', 'HLA-B*35:91', 'HLA-B*35:92', 'HLA-B*35:93',
         'HLA-B*35:94', 'HLA-B*35:95', 'HLA-B*35:96', 'HLA-B*35:97', 'HLA-B*35:98', 'HLA-B*35:99', 'HLA-B*37:01', 'HLA-B*37:02', 'HLA-B*37:04',
         'HLA-B*37:05', 'HLA-B*37:06', 'HLA-B*37:07', 'HLA-B*37:08', 'HLA-B*37:09', 'HLA-B*37:10', 'HLA-B*37:11', 'HLA-B*37:12', 'HLA-B*37:13',
         'HLA-B*37:14', 'HLA-B*37:15', 'HLA-B*37:17', 'HLA-B*37:18', 'HLA-B*37:19', 'HLA-B*37:20', 'HLA-B*37:21', 'HLA-B*37:22', 'HLA-B*37:23',
         'HLA-B*38:01', 'HLA-B*38:02', 'HLA-B*38:03', 'HLA-B*38:04', 'HLA-B*38:05', 'HLA-B*38:06', 'HLA-B*38:07', 'HLA-B*38:08', 'HLA-B*38:09',
         'HLA-B*38:10', 'HLA-B*38:11', 'HLA-B*38:12', 'HLA-B*38:13', 'HLA-B*38:14', 'HLA-B*38:15', 'HLA-B*38:16', 'HLA-B*38:17', 'HLA-B*38:18',
         'HLA-B*38:19', 'HLA-B*38:20', 'HLA-B*38:21', 'HLA-B*38:22', 'HLA-B*38:23', 'HLA-B*39:01', 'HLA-B*39:02', 'HLA-B*39:03', 'HLA-B*39:04',
         'HLA-B*39:05', 'HLA-B*39:06', 'HLA-B*39:07', 'HLA-B*39:08', 'HLA-B*39:09', 'HLA-B*39:10', 'HLA-B*39:11', 'HLA-B*39:12', 'HLA-B*39:13',
         'HLA-B*39:14', 'HLA-B*39:15', 'HLA-B*39:16', 'HLA-B*39:17', 'HLA-B*39:18', 'HLA-B*39:19', 'HLA-B*39:20', 'HLA-B*39:22', 'HLA-B*39:23',
         'HLA-B*39:24', 'HLA-B*39:26', 'HLA-B*39:27', 'HLA-B*39:28', 'HLA-B*39:29', 'HLA-B*39:30', 'HLA-B*39:31', 'HLA-B*39:32', 'HLA-B*39:33',
         'HLA-B*39:34', 'HLA-B*39:35', 'HLA-B*39:36', 'HLA-B*39:37', 'HLA-B*39:39', 'HLA-B*39:41', 'HLA-B*39:42', 'HLA-B*39:43', 'HLA-B*39:44',
         'HLA-B*39:45', 'HLA-B*39:46', 'HLA-B*39:47', 'HLA-B*39:48', 'HLA-B*39:49', 'HLA-B*39:50', 'HLA-B*39:51', 'HLA-B*39:52', 'HLA-B*39:53',
         'HLA-B*39:54', 'HLA-B*39:55', 'HLA-B*39:56', 'HLA-B*39:57', 'HLA-B*39:58', 'HLA-B*39:59', 'HLA-B*39:60', 'HLA-B*40:01', 'HLA-B*40:02',
         'HLA-B*40:03', 'HLA-B*40:04', 'HLA-B*40:05', 'HLA-B*40:06', 'HLA-B*40:07', 'HLA-B*40:08', 'HLA-B*40:09', 'HLA-B*40:10', 'HLA-B*40:100',
         'HLA-B*40:101', 'HLA-B*40:102', 'HLA-B*40:103', 'HLA-B*40:104', 'HLA-B*40:105', 'HLA-B*40:106', 'HLA-B*40:107', 'HLA-B*40:108',
         'HLA-B*40:109', 'HLA-B*40:11', 'HLA-B*40:110', 'HLA-B*40:111', 'HLA-B*40:112', 'HLA-B*40:113', 'HLA-B*40:114', 'HLA-B*40:115',
         'HLA-B*40:116', 'HLA-B*40:117', 'HLA-B*40:119', 'HLA-B*40:12', 'HLA-B*40:120', 'HLA-B*40:121', 'HLA-B*40:122', 'HLA-B*40:123',
         'HLA-B*40:124', 'HLA-B*40:125', 'HLA-B*40:126', 'HLA-B*40:127', 'HLA-B*40:128', 'HLA-B*40:129', 'HLA-B*40:13', 'HLA-B*40:130',
         'HLA-B*40:131', 'HLA-B*40:132', 'HLA-B*40:134', 'HLA-B*40:135', 'HLA-B*40:136', 'HLA-B*40:137', 'HLA-B*40:138', 'HLA-B*40:139',
         'HLA-B*40:14', 'HLA-B*40:140', 'HLA-B*40:141', 'HLA-B*40:143', 'HLA-B*40:145', 'HLA-B*40:146', 'HLA-B*40:147', 'HLA-B*40:15',
         'HLA-B*40:16', 'HLA-B*40:18', 'HLA-B*40:19', 'HLA-B*40:20', 'HLA-B*40:21', 'HLA-B*40:23', 'HLA-B*40:24', 'HLA-B*40:25', 'HLA-B*40:26',
         'HLA-B*40:27', 'HLA-B*40:28', 'HLA-B*40:29', 'HLA-B*40:30', 'HLA-B*40:31', 'HLA-B*40:32', 'HLA-B*40:33', 'HLA-B*40:34', 'HLA-B*40:35',
         'HLA-B*40:36', 'HLA-B*40:37', 'HLA-B*40:38', 'HLA-B*40:39', 'HLA-B*40:40', 'HLA-B*40:42', 'HLA-B*40:43', 'HLA-B*40:44', 'HLA-B*40:45',
         'HLA-B*40:46', 'HLA-B*40:47', 'HLA-B*40:48', 'HLA-B*40:49', 'HLA-B*40:50', 'HLA-B*40:51', 'HLA-B*40:52', 'HLA-B*40:53', 'HLA-B*40:54',
         'HLA-B*40:55', 'HLA-B*40:56', 'HLA-B*40:57', 'HLA-B*40:58', 'HLA-B*40:59', 'HLA-B*40:60', 'HLA-B*40:61', 'HLA-B*40:62', 'HLA-B*40:63',
         'HLA-B*40:64', 'HLA-B*40:65', 'HLA-B*40:66', 'HLA-B*40:67', 'HLA-B*40:68', 'HLA-B*40:69', 'HLA-B*40:70', 'HLA-B*40:71', 'HLA-B*40:72',
         'HLA-B*40:73', 'HLA-B*40:74', 'HLA-B*40:75', 'HLA-B*40:76', 'HLA-B*40:77', 'HLA-B*40:78', 'HLA-B*40:79', 'HLA-B*40:80', 'HLA-B*40:81',
         'HLA-B*40:82', 'HLA-B*40:83', 'HLA-B*40:84', 'HLA-B*40:85', 'HLA-B*40:86', 'HLA-B*40:87', 'HLA-B*40:88', 'HLA-B*40:89', 'HLA-B*40:90',
         'HLA-B*40:91', 'HLA-B*40:92', 'HLA-B*40:93', 'HLA-B*40:94', 'HLA-B*40:95', 'HLA-B*40:96', 'HLA-B*40:97', 'HLA-B*40:98', 'HLA-B*40:99',
         'HLA-B*41:01', 'HLA-B*41:02', 'HLA-B*41:03', 'HLA-B*41:04', 'HLA-B*41:05', 'HLA-B*41:06', 'HLA-B*41:07', 'HLA-B*41:08', 'HLA-B*41:09',
         'HLA-B*41:10', 'HLA-B*41:11', 'HLA-B*41:12', 'HLA-B*42:01', 'HLA-B*42:02', 'HLA-B*42:04', 'HLA-B*42:05', 'HLA-B*42:06', 'HLA-B*42:07',
         'HLA-B*42:08', 'HLA-B*42:09', 'HLA-B*42:10', 'HLA-B*42:11', 'HLA-B*42:12', 'HLA-B*42:13', 'HLA-B*42:14', 'HLA-B*44:02', 'HLA-B*44:03',
         'HLA-B*44:04', 'HLA-B*44:05', 'HLA-B*44:06', 'HLA-B*44:07', 'HLA-B*44:08', 'HLA-B*44:09', 'HLA-B*44:10', 'HLA-B*44:100', 'HLA-B*44:101',
         'HLA-B*44:102', 'HLA-B*44:103', 'HLA-B*44:104', 'HLA-B*44:105', 'HLA-B*44:106', 'HLA-B*44:107', 'HLA-B*44:109', 'HLA-B*44:11',
         'HLA-B*44:110', 'HLA-B*44:12', 'HLA-B*44:13', 'HLA-B*44:14', 'HLA-B*44:15', 'HLA-B*44:16', 'HLA-B*44:17', 'HLA-B*44:18', 'HLA-B*44:20',
         'HLA-B*44:21', 'HLA-B*44:22', 'HLA-B*44:24', 'HLA-B*44:25', 'HLA-B*44:26', 'HLA-B*44:27', 'HLA-B*44:28', 'HLA-B*44:29', 'HLA-B*44:30',
         'HLA-B*44:31', 'HLA-B*44:32', 'HLA-B*44:33', 'HLA-B*44:34', 'HLA-B*44:35', 'HLA-B*44:36', 'HLA-B*44:37', 'HLA-B*44:38', 'HLA-B*44:39',
         'HLA-B*44:40', 'HLA-B*44:41', 'HLA-B*44:42', 'HLA-B*44:43', 'HLA-B*44:44', 'HLA-B*44:45', 'HLA-B*44:46', 'HLA-B*44:47', 'HLA-B*44:48',
         'HLA-B*44:49', 'HLA-B*44:50', 'HLA-B*44:51', 'HLA-B*44:53', 'HLA-B*44:54', 'HLA-B*44:55', 'HLA-B*44:57', 'HLA-B*44:59', 'HLA-B*44:60',
         'HLA-B*44:62', 'HLA-B*44:63', 'HLA-B*44:64', 'HLA-B*44:65', 'HLA-B*44:66', 'HLA-B*44:67', 'HLA-B*44:68', 'HLA-B*44:69', 'HLA-B*44:70',
         'HLA-B*44:71', 'HLA-B*44:72', 'HLA-B*44:73', 'HLA-B*44:74', 'HLA-B*44:75', 'HLA-B*44:76', 'HLA-B*44:77', 'HLA-B*44:78', 'HLA-B*44:79',
         'HLA-B*44:80', 'HLA-B*44:81', 'HLA-B*44:82', 'HLA-B*44:83', 'HLA-B*44:84', 'HLA-B*44:85', 'HLA-B*44:86', 'HLA-B*44:87', 'HLA-B*44:88',
         'HLA-B*44:89', 'HLA-B*44:90', 'HLA-B*44:91', 'HLA-B*44:92', 'HLA-B*44:93', 'HLA-B*44:94', 'HLA-B*44:95', 'HLA-B*44:96', 'HLA-B*44:97',
         'HLA-B*44:98', 'HLA-B*44:99', 'HLA-B*45:01', 'HLA-B*45:02', 'HLA-B*45:03', 'HLA-B*45:04', 'HLA-B*45:05', 'HLA-B*45:06', 'HLA-B*45:07',
         'HLA-B*45:08', 'HLA-B*45:09', 'HLA-B*45:10', 'HLA-B*45:11', 'HLA-B*45:12', 'HLA-B*46:01', 'HLA-B*46:02', 'HLA-B*46:03', 'HLA-B*46:04',
         'HLA-B*46:05', 'HLA-B*46:06', 'HLA-B*46:08', 'HLA-B*46:09', 'HLA-B*46:10', 'HLA-B*46:11', 'HLA-B*46:12', 'HLA-B*46:13', 'HLA-B*46:14',
         'HLA-B*46:16', 'HLA-B*46:17', 'HLA-B*46:18', 'HLA-B*46:19', 'HLA-B*46:20', 'HLA-B*46:21', 'HLA-B*46:22', 'HLA-B*46:23', 'HLA-B*46:24',
         'HLA-B*47:01', 'HLA-B*47:02', 'HLA-B*47:03', 'HLA-B*47:04', 'HLA-B*47:05', 'HLA-B*47:06', 'HLA-B*47:07', 'HLA-B*48:01', 'HLA-B*48:02',
         'HLA-B*48:03', 'HLA-B*48:04', 'HLA-B*48:05', 'HLA-B*48:06', 'HLA-B*48:07', 'HLA-B*48:08', 'HLA-B*48:09', 'HLA-B*48:10', 'HLA-B*48:11',
         'HLA-B*48:12', 'HLA-B*48:13', 'HLA-B*48:14', 'HLA-B*48:15', 'HLA-B*48:16', 'HLA-B*48:17', 'HLA-B*48:18', 'HLA-B*48:19', 'HLA-B*48:20',
         'HLA-B*48:21', 'HLA-B*48:22', 'HLA-B*48:23', 'HLA-B*49:01', 'HLA-B*49:02', 'HLA-B*49:03', 'HLA-B*49:04', 'HLA-B*49:05', 'HLA-B*49:06',
         'HLA-B*49:07', 'HLA-B*49:08', 'HLA-B*49:09', 'HLA-B*49:10', 'HLA-B*50:01', 'HLA-B*50:02', 'HLA-B*50:04', 'HLA-B*50:05', 'HLA-B*50:06',
         'HLA-B*50:07', 'HLA-B*50:08', 'HLA-B*50:09', 'HLA-B*51:01', 'HLA-B*51:02', 'HLA-B*51:03', 'HLA-B*51:04', 'HLA-B*51:05', 'HLA-B*51:06',
         'HLA-B*51:07', 'HLA-B*51:08', 'HLA-B*51:09', 'HLA-B*51:12', 'HLA-B*51:13', 'HLA-B*51:14', 'HLA-B*51:15', 'HLA-B*51:16', 'HLA-B*51:17',
         'HLA-B*51:18', 'HLA-B*51:19', 'HLA-B*51:20', 'HLA-B*51:21', 'HLA-B*51:22', 'HLA-B*51:23', 'HLA-B*51:24', 'HLA-B*51:26', 'HLA-B*51:28',
         'HLA-B*51:29', 'HLA-B*51:30', 'HLA-B*51:31', 'HLA-B*51:32', 'HLA-B*51:33', 'HLA-B*51:34', 'HLA-B*51:35', 'HLA-B*51:36', 'HLA-B*51:37',
         'HLA-B*51:38', 'HLA-B*51:39', 'HLA-B*51:40', 'HLA-B*51:42', 'HLA-B*51:43', 'HLA-B*51:45', 'HLA-B*51:46', 'HLA-B*51:48', 'HLA-B*51:49',
         'HLA-B*51:50', 'HLA-B*51:51', 'HLA-B*51:52', 'HLA-B*51:53', 'HLA-B*51:54', 'HLA-B*51:55', 'HLA-B*51:56', 'HLA-B*51:57', 'HLA-B*51:58',
         'HLA-B*51:59', 'HLA-B*51:60', 'HLA-B*51:61', 'HLA-B*51:62', 'HLA-B*51:63', 'HLA-B*51:64', 'HLA-B*51:65', 'HLA-B*51:66', 'HLA-B*51:67',
         'HLA-B*51:68', 'HLA-B*51:69', 'HLA-B*51:70', 'HLA-B*51:71', 'HLA-B*51:72', 'HLA-B*51:73', 'HLA-B*51:74', 'HLA-B*51:75', 'HLA-B*51:76',
         'HLA-B*51:77', 'HLA-B*51:78', 'HLA-B*51:79', 'HLA-B*51:80', 'HLA-B*51:81', 'HLA-B*51:82', 'HLA-B*51:83', 'HLA-B*51:84', 'HLA-B*51:85',
         'HLA-B*51:86', 'HLA-B*51:87', 'HLA-B*51:88', 'HLA-B*51:89', 'HLA-B*51:90', 'HLA-B*51:91', 'HLA-B*51:92', 'HLA-B*51:93', 'HLA-B*51:94',
         'HLA-B*51:95', 'HLA-B*51:96', 'HLA-B*52:01', 'HLA-B*52:02', 'HLA-B*52:03', 'HLA-B*52:04', 'HLA-B*52:05', 'HLA-B*52:06', 'HLA-B*52:07',
         'HLA-B*52:08', 'HLA-B*52:09', 'HLA-B*52:10', 'HLA-B*52:11', 'HLA-B*52:12', 'HLA-B*52:13', 'HLA-B*52:14', 'HLA-B*52:15', 'HLA-B*52:16',
         'HLA-B*52:17', 'HLA-B*52:18', 'HLA-B*52:19', 'HLA-B*52:20', 'HLA-B*52:21', 'HLA-B*53:01', 'HLA-B*53:02', 'HLA-B*53:03', 'HLA-B*53:04',
         'HLA-B*53:05', 'HLA-B*53:06', 'HLA-B*53:07', 'HLA-B*53:08', 'HLA-B*53:09', 'HLA-B*53:10', 'HLA-B*53:11', 'HLA-B*53:12', 'HLA-B*53:13',
         'HLA-B*53:14', 'HLA-B*53:15', 'HLA-B*53:16', 'HLA-B*53:17', 'HLA-B*53:18', 'HLA-B*53:19', 'HLA-B*53:20', 'HLA-B*53:21', 'HLA-B*53:22',
         'HLA-B*53:23', 'HLA-B*54:01', 'HLA-B*54:02', 'HLA-B*54:03', 'HLA-B*54:04', 'HLA-B*54:06', 'HLA-B*54:07', 'HLA-B*54:09', 'HLA-B*54:10',
         'HLA-B*54:11', 'HLA-B*54:12', 'HLA-B*54:13', 'HLA-B*54:14', 'HLA-B*54:15', 'HLA-B*54:16', 'HLA-B*54:17', 'HLA-B*54:18', 'HLA-B*54:19',
         'HLA-B*54:20', 'HLA-B*54:21', 'HLA-B*54:22', 'HLA-B*54:23', 'HLA-B*55:01', 'HLA-B*55:02', 'HLA-B*55:03', 'HLA-B*55:04', 'HLA-B*55:05',
         'HLA-B*55:07', 'HLA-B*55:08', 'HLA-B*55:09', 'HLA-B*55:10', 'HLA-B*55:11', 'HLA-B*55:12', 'HLA-B*55:13', 'HLA-B*55:14', 'HLA-B*55:15',
         'HLA-B*55:16', 'HLA-B*55:17', 'HLA-B*55:18', 'HLA-B*55:19', 'HLA-B*55:20', 'HLA-B*55:21', 'HLA-B*55:22', 'HLA-B*55:23', 'HLA-B*55:24',
         'HLA-B*55:25', 'HLA-B*55:26', 'HLA-B*55:27', 'HLA-B*55:28', 'HLA-B*55:29', 'HLA-B*55:30', 'HLA-B*55:31', 'HLA-B*55:32', 'HLA-B*55:33',
         'HLA-B*55:34', 'HLA-B*55:35', 'HLA-B*55:36', 'HLA-B*55:37', 'HLA-B*55:38', 'HLA-B*55:39', 'HLA-B*55:40', 'HLA-B*55:41', 'HLA-B*55:42',
         'HLA-B*55:43', 'HLA-B*56:01', 'HLA-B*56:02', 'HLA-B*56:03', 'HLA-B*56:04', 'HLA-B*56:05', 'HLA-B*56:06', 'HLA-B*56:07', 'HLA-B*56:08',
         'HLA-B*56:09', 'HLA-B*56:10', 'HLA-B*56:11', 'HLA-B*56:12', 'HLA-B*56:13', 'HLA-B*56:14', 'HLA-B*56:15', 'HLA-B*56:16', 'HLA-B*56:17',
         'HLA-B*56:18', 'HLA-B*56:20', 'HLA-B*56:21', 'HLA-B*56:22', 'HLA-B*56:23', 'HLA-B*56:24', 'HLA-B*56:25', 'HLA-B*56:26', 'HLA-B*56:27',
         'HLA-B*56:29', 'HLA-B*57:01', 'HLA-B*57:02', 'HLA-B*57:03', 'HLA-B*57:04', 'HLA-B*57:05', 'HLA-B*57:06', 'HLA-B*57:07', 'HLA-B*57:08',
         'HLA-B*57:09', 'HLA-B*57:10', 'HLA-B*57:11', 'HLA-B*57:12', 'HLA-B*57:13', 'HLA-B*57:14', 'HLA-B*57:15', 'HLA-B*57:16', 'HLA-B*57:17',
         'HLA-B*57:18', 'HLA-B*57:19', 'HLA-B*57:20', 'HLA-B*57:21', 'HLA-B*57:22', 'HLA-B*57:23', 'HLA-B*57:24', 'HLA-B*57:25', 'HLA-B*57:26',
         'HLA-B*57:27', 'HLA-B*57:29', 'HLA-B*57:30', 'HLA-B*57:31', 'HLA-B*57:32', 'HLA-B*58:01', 'HLA-B*58:02', 'HLA-B*58:04', 'HLA-B*58:05',
         'HLA-B*58:06', 'HLA-B*58:07', 'HLA-B*58:08', 'HLA-B*58:09', 'HLA-B*58:11', 'HLA-B*58:12', 'HLA-B*58:13', 'HLA-B*58:14', 'HLA-B*58:15',
         'HLA-B*58:16', 'HLA-B*58:18', 'HLA-B*58:19', 'HLA-B*58:20', 'HLA-B*58:21', 'HLA-B*58:22', 'HLA-B*58:23', 'HLA-B*58:24', 'HLA-B*58:25',
         'HLA-B*58:26', 'HLA-B*58:27', 'HLA-B*58:28', 'HLA-B*58:29', 'HLA-B*58:30', 'HLA-B*59:01', 'HLA-B*59:02', 'HLA-B*59:03', 'HLA-B*59:04',
         'HLA-B*59:05', 'HLA-B*67:01', 'HLA-B*67:02', 'HLA-B*73:01', 'HLA-B*73:02', 'HLA-B*78:01', 'HLA-B*78:02', 'HLA-B*78:03', 'HLA-B*78:04',
         'HLA-B*78:05', 'HLA-B*78:06', 'HLA-B*78:07', 'HLA-B*81:01', 'HLA-B*81:02', 'HLA-B*81:03', 'HLA-B*81:05', 'HLA-B*82:01', 'HLA-B*82:02',
         'HLA-B*82:03', 'HLA-B*83:01', 'HLA-C*01:02', 'HLA-C*01:03', 'HLA-C*01:04', 'HLA-C*01:05', 'HLA-C*01:06', 'HLA-C*01:07', 'HLA-C*01:08',
         'HLA-C*01:09', 'HLA-C*01:10', 'HLA-C*01:11', 'HLA-C*01:12', 'HLA-C*01:13', 'HLA-C*01:14', 'HLA-C*01:15', 'HLA-C*01:16', 'HLA-C*01:17',
         'HLA-C*01:18', 'HLA-C*01:19', 'HLA-C*01:20', 'HLA-C*01:21', 'HLA-C*01:22', 'HLA-C*01:23', 'HLA-C*01:24', 'HLA-C*01:25', 'HLA-C*01:26',
         'HLA-C*01:27', 'HLA-C*01:28', 'HLA-C*01:29', 'HLA-C*01:30', 'HLA-C*01:31', 'HLA-C*01:32', 'HLA-C*01:33', 'HLA-C*01:34', 'HLA-C*01:35',
         'HLA-C*01:36', 'HLA-C*01:38', 'HLA-C*01:39', 'HLA-C*01:40', 'HLA-C*02:02', 'HLA-C*02:03', 'HLA-C*02:04', 'HLA-C*02:05', 'HLA-C*02:06',
         'HLA-C*02:07', 'HLA-C*02:08', 'HLA-C*02:09', 'HLA-C*02:10', 'HLA-C*02:11', 'HLA-C*02:12', 'HLA-C*02:13', 'HLA-C*02:14', 'HLA-C*02:15',
         'HLA-C*02:16', 'HLA-C*02:17', 'HLA-C*02:18', 'HLA-C*02:19', 'HLA-C*02:20', 'HLA-C*02:21', 'HLA-C*02:22', 'HLA-C*02:23', 'HLA-C*02:24',
         'HLA-C*02:26', 'HLA-C*02:27', 'HLA-C*02:28', 'HLA-C*02:29', 'HLA-C*02:30', 'HLA-C*02:31', 'HLA-C*02:32', 'HLA-C*02:33', 'HLA-C*02:34',
         'HLA-C*02:35', 'HLA-C*02:36', 'HLA-C*02:37', 'HLA-C*02:39', 'HLA-C*02:40', 'HLA-C*03:01', 'HLA-C*03:02', 'HLA-C*03:03', 'HLA-C*03:04',
         'HLA-C*03:05', 'HLA-C*03:06', 'HLA-C*03:07', 'HLA-C*03:08', 'HLA-C*03:09', 'HLA-C*03:10', 'HLA-C*03:11', 'HLA-C*03:12', 'HLA-C*03:13',
         'HLA-C*03:14', 'HLA-C*03:15', 'HLA-C*03:16', 'HLA-C*03:17', 'HLA-C*03:18', 'HLA-C*03:19', 'HLA-C*03:21', 'HLA-C*03:23', 'HLA-C*03:24',
         'HLA-C*03:25', 'HLA-C*03:26', 'HLA-C*03:27', 'HLA-C*03:28', 'HLA-C*03:29', 'HLA-C*03:30', 'HLA-C*03:31', 'HLA-C*03:32', 'HLA-C*03:33',
         'HLA-C*03:34', 'HLA-C*03:35', 'HLA-C*03:36', 'HLA-C*03:37', 'HLA-C*03:38', 'HLA-C*03:39', 'HLA-C*03:40', 'HLA-C*03:41', 'HLA-C*03:42',
         'HLA-C*03:43', 'HLA-C*03:44', 'HLA-C*03:45', 'HLA-C*03:46', 'HLA-C*03:47', 'HLA-C*03:48', 'HLA-C*03:49', 'HLA-C*03:50', 'HLA-C*03:51',
         'HLA-C*03:52', 'HLA-C*03:53', 'HLA-C*03:54', 'HLA-C*03:55', 'HLA-C*03:56', 'HLA-C*03:57', 'HLA-C*03:58', 'HLA-C*03:59', 'HLA-C*03:60',
         'HLA-C*03:61', 'HLA-C*03:62', 'HLA-C*03:63', 'HLA-C*03:64', 'HLA-C*03:65', 'HLA-C*03:66', 'HLA-C*03:67', 'HLA-C*03:68', 'HLA-C*03:69',
         'HLA-C*03:70', 'HLA-C*03:71', 'HLA-C*03:72', 'HLA-C*03:73', 'HLA-C*03:74', 'HLA-C*03:75', 'HLA-C*03:76', 'HLA-C*03:77', 'HLA-C*03:78',
         'HLA-C*03:79', 'HLA-C*03:80', 'HLA-C*03:81', 'HLA-C*03:82', 'HLA-C*03:83', 'HLA-C*03:84', 'HLA-C*03:85', 'HLA-C*03:86', 'HLA-C*03:87',
         'HLA-C*03:88', 'HLA-C*03:89', 'HLA-C*03:90', 'HLA-C*03:91', 'HLA-C*03:92', 'HLA-C*03:93', 'HLA-C*03:94', 'HLA-C*04:01', 'HLA-C*04:03',
         'HLA-C*04:04', 'HLA-C*04:05', 'HLA-C*04:06', 'HLA-C*04:07', 'HLA-C*04:08', 'HLA-C*04:10', 'HLA-C*04:11', 'HLA-C*04:12', 'HLA-C*04:13',
         'HLA-C*04:14', 'HLA-C*04:15', 'HLA-C*04:16', 'HLA-C*04:17', 'HLA-C*04:18', 'HLA-C*04:19', 'HLA-C*04:20', 'HLA-C*04:23', 'HLA-C*04:24',
         'HLA-C*04:25', 'HLA-C*04:26', 'HLA-C*04:27', 'HLA-C*04:28', 'HLA-C*04:29', 'HLA-C*04:30', 'HLA-C*04:31', 'HLA-C*04:32', 'HLA-C*04:33',
         'HLA-C*04:34', 'HLA-C*04:35', 'HLA-C*04:36', 'HLA-C*04:37', 'HLA-C*04:38', 'HLA-C*04:39', 'HLA-C*04:40', 'HLA-C*04:41', 'HLA-C*04:42',
         'HLA-C*04:43', 'HLA-C*04:44', 'HLA-C*04:45', 'HLA-C*04:46', 'HLA-C*04:47', 'HLA-C*04:48', 'HLA-C*04:49', 'HLA-C*04:50', 'HLA-C*04:51',
         'HLA-C*04:52', 'HLA-C*04:53', 'HLA-C*04:54', 'HLA-C*04:55', 'HLA-C*04:56', 'HLA-C*04:57', 'HLA-C*04:58', 'HLA-C*04:60', 'HLA-C*04:61',
         'HLA-C*04:62', 'HLA-C*04:63', 'HLA-C*04:64', 'HLA-C*04:65', 'HLA-C*04:66', 'HLA-C*04:67', 'HLA-C*04:68', 'HLA-C*04:69', 'HLA-C*04:70',
         'HLA-C*05:01', 'HLA-C*05:03', 'HLA-C*05:04', 'HLA-C*05:05', 'HLA-C*05:06', 'HLA-C*05:08', 'HLA-C*05:09', 'HLA-C*05:10', 'HLA-C*05:11',
         'HLA-C*05:12', 'HLA-C*05:13', 'HLA-C*05:14', 'HLA-C*05:15', 'HLA-C*05:16', 'HLA-C*05:17', 'HLA-C*05:18', 'HLA-C*05:19', 'HLA-C*05:20',
         'HLA-C*05:21', 'HLA-C*05:22', 'HLA-C*05:23', 'HLA-C*05:24', 'HLA-C*05:25', 'HLA-C*05:26', 'HLA-C*05:27', 'HLA-C*05:28', 'HLA-C*05:29',
         'HLA-C*05:30', 'HLA-C*05:31', 'HLA-C*05:32', 'HLA-C*05:33', 'HLA-C*05:34', 'HLA-C*05:35', 'HLA-C*05:36', 'HLA-C*05:37', 'HLA-C*05:38',
         'HLA-C*05:39', 'HLA-C*05:40', 'HLA-C*05:41', 'HLA-C*05:42', 'HLA-C*05:43', 'HLA-C*05:44', 'HLA-C*05:45', 'HLA-C*06:02', 'HLA-C*06:03',
         'HLA-C*06:04', 'HLA-C*06:05', 'HLA-C*06:06', 'HLA-C*06:07', 'HLA-C*06:08', 'HLA-C*06:09', 'HLA-C*06:10', 'HLA-C*06:11', 'HLA-C*06:12',
         'HLA-C*06:13', 'HLA-C*06:14', 'HLA-C*06:15', 'HLA-C*06:17', 'HLA-C*06:18', 'HLA-C*06:19', 'HLA-C*06:20', 'HLA-C*06:21', 'HLA-C*06:22',
         'HLA-C*06:23', 'HLA-C*06:24', 'HLA-C*06:25', 'HLA-C*06:26', 'HLA-C*06:27', 'HLA-C*06:28', 'HLA-C*06:29', 'HLA-C*06:30', 'HLA-C*06:31',
         'HLA-C*06:32', 'HLA-C*06:33', 'HLA-C*06:34', 'HLA-C*06:35', 'HLA-C*06:36', 'HLA-C*06:37', 'HLA-C*06:38', 'HLA-C*06:39', 'HLA-C*06:40',
         'HLA-C*06:41', 'HLA-C*06:42', 'HLA-C*06:43', 'HLA-C*06:44', 'HLA-C*06:45', 'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*07:03', 'HLA-C*07:04',
         'HLA-C*07:05', 'HLA-C*07:06', 'HLA-C*07:07', 'HLA-C*07:08', 'HLA-C*07:09', 'HLA-C*07:10', 'HLA-C*07:100', 'HLA-C*07:101', 'HLA-C*07:102',
         'HLA-C*07:103', 'HLA-C*07:105', 'HLA-C*07:106', 'HLA-C*07:107', 'HLA-C*07:108', 'HLA-C*07:109', 'HLA-C*07:11', 'HLA-C*07:110',
         'HLA-C*07:111', 'HLA-C*07:112', 'HLA-C*07:113', 'HLA-C*07:114', 'HLA-C*07:115', 'HLA-C*07:116', 'HLA-C*07:117', 'HLA-C*07:118',
         'HLA-C*07:119', 'HLA-C*07:12', 'HLA-C*07:120', 'HLA-C*07:122', 'HLA-C*07:123', 'HLA-C*07:124', 'HLA-C*07:125', 'HLA-C*07:126',
         'HLA-C*07:127', 'HLA-C*07:128', 'HLA-C*07:129', 'HLA-C*07:13', 'HLA-C*07:130', 'HLA-C*07:131', 'HLA-C*07:132', 'HLA-C*07:133',
         'HLA-C*07:134', 'HLA-C*07:135', 'HLA-C*07:136', 'HLA-C*07:137', 'HLA-C*07:138', 'HLA-C*07:139', 'HLA-C*07:14', 'HLA-C*07:140',
         'HLA-C*07:141', 'HLA-C*07:142', 'HLA-C*07:143', 'HLA-C*07:144', 'HLA-C*07:145', 'HLA-C*07:146', 'HLA-C*07:147', 'HLA-C*07:148',
         'HLA-C*07:149', 'HLA-C*07:15', 'HLA-C*07:16', 'HLA-C*07:17', 'HLA-C*07:18', 'HLA-C*07:19', 'HLA-C*07:20', 'HLA-C*07:21', 'HLA-C*07:22',
         'HLA-C*07:23', 'HLA-C*07:24', 'HLA-C*07:25', 'HLA-C*07:26', 'HLA-C*07:27', 'HLA-C*07:28', 'HLA-C*07:29', 'HLA-C*07:30', 'HLA-C*07:31',
         'HLA-C*07:35', 'HLA-C*07:36', 'HLA-C*07:37', 'HLA-C*07:38', 'HLA-C*07:39', 'HLA-C*07:40', 'HLA-C*07:41', 'HLA-C*07:42', 'HLA-C*07:43',
         'HLA-C*07:44', 'HLA-C*07:45', 'HLA-C*07:46', 'HLA-C*07:47', 'HLA-C*07:48', 'HLA-C*07:49', 'HLA-C*07:50', 'HLA-C*07:51', 'HLA-C*07:52',
         'HLA-C*07:53', 'HLA-C*07:54', 'HLA-C*07:56', 'HLA-C*07:57', 'HLA-C*07:58', 'HLA-C*07:59', 'HLA-C*07:60', 'HLA-C*07:62', 'HLA-C*07:63',
         'HLA-C*07:64', 'HLA-C*07:65', 'HLA-C*07:66', 'HLA-C*07:67', 'HLA-C*07:68', 'HLA-C*07:69', 'HLA-C*07:70', 'HLA-C*07:71', 'HLA-C*07:72',
         'HLA-C*07:73', 'HLA-C*07:74', 'HLA-C*07:75', 'HLA-C*07:76', 'HLA-C*07:77', 'HLA-C*07:78', 'HLA-C*07:79', 'HLA-C*07:80', 'HLA-C*07:81',
         'HLA-C*07:82', 'HLA-C*07:83', 'HLA-C*07:84', 'HLA-C*07:85', 'HLA-C*07:86', 'HLA-C*07:87', 'HLA-C*07:88', 'HLA-C*07:89', 'HLA-C*07:90',
         'HLA-C*07:91', 'HLA-C*07:92', 'HLA-C*07:93', 'HLA-C*07:94', 'HLA-C*07:95', 'HLA-C*07:96', 'HLA-C*07:97', 'HLA-C*07:99', 'HLA-C*08:01',
         'HLA-C*08:02', 'HLA-C*08:03', 'HLA-C*08:04', 'HLA-C*08:05', 'HLA-C*08:06', 'HLA-C*08:07', 'HLA-C*08:08', 'HLA-C*08:09', 'HLA-C*08:10',
         'HLA-C*08:11', 'HLA-C*08:12', 'HLA-C*08:13', 'HLA-C*08:14', 'HLA-C*08:15', 'HLA-C*08:16', 'HLA-C*08:17', 'HLA-C*08:18', 'HLA-C*08:19',
         'HLA-C*08:20', 'HLA-C*08:21', 'HLA-C*08:22', 'HLA-C*08:23', 'HLA-C*08:24', 'HLA-C*08:25', 'HLA-C*08:27', 'HLA-C*08:28', 'HLA-C*08:29',
         'HLA-C*08:30', 'HLA-C*08:31', 'HLA-C*08:32', 'HLA-C*08:33', 'HLA-C*08:34', 'HLA-C*08:35', 'HLA-C*12:02', 'HLA-C*12:03', 'HLA-C*12:04',
         'HLA-C*12:05', 'HLA-C*12:06', 'HLA-C*12:07', 'HLA-C*12:08', 'HLA-C*12:09', 'HLA-C*12:10', 'HLA-C*12:11', 'HLA-C*12:12', 'HLA-C*12:13',
         'HLA-C*12:14', 'HLA-C*12:15', 'HLA-C*12:16', 'HLA-C*12:17', 'HLA-C*12:18', 'HLA-C*12:19', 'HLA-C*12:20', 'HLA-C*12:21', 'HLA-C*12:22',
         'HLA-C*12:23', 'HLA-C*12:24', 'HLA-C*12:25', 'HLA-C*12:26', 'HLA-C*12:27', 'HLA-C*12:28', 'HLA-C*12:29', 'HLA-C*12:30', 'HLA-C*12:31',
         'HLA-C*12:32', 'HLA-C*12:33', 'HLA-C*12:34', 'HLA-C*12:35', 'HLA-C*12:36', 'HLA-C*12:37', 'HLA-C*12:38', 'HLA-C*12:40', 'HLA-C*12:41',
         'HLA-C*12:43', 'HLA-C*12:44', 'HLA-C*14:02', 'HLA-C*14:03', 'HLA-C*14:04', 'HLA-C*14:05', 'HLA-C*14:06', 'HLA-C*14:08', 'HLA-C*14:09',
         'HLA-C*14:10', 'HLA-C*14:11', 'HLA-C*14:12', 'HLA-C*14:13', 'HLA-C*14:14', 'HLA-C*14:15', 'HLA-C*14:16', 'HLA-C*14:17', 'HLA-C*14:18',
         'HLA-C*14:19', 'HLA-C*14:20', 'HLA-C*15:02', 'HLA-C*15:03', 'HLA-C*15:04', 'HLA-C*15:05', 'HLA-C*15:06', 'HLA-C*15:07', 'HLA-C*15:08',
         'HLA-C*15:09', 'HLA-C*15:10', 'HLA-C*15:11', 'HLA-C*15:12', 'HLA-C*15:13', 'HLA-C*15:15', 'HLA-C*15:16', 'HLA-C*15:17', 'HLA-C*15:18',
         'HLA-C*15:19', 'HLA-C*15:20', 'HLA-C*15:21', 'HLA-C*15:22', 'HLA-C*15:23', 'HLA-C*15:24', 'HLA-C*15:25', 'HLA-C*15:26', 'HLA-C*15:27',
         'HLA-C*15:28', 'HLA-C*15:29', 'HLA-C*15:30', 'HLA-C*15:31', 'HLA-C*15:33', 'HLA-C*15:34', 'HLA-C*15:35', 'HLA-C*16:01', 'HLA-C*16:02',
         'HLA-C*16:04', 'HLA-C*16:06', 'HLA-C*16:07', 'HLA-C*16:08', 'HLA-C*16:09', 'HLA-C*16:10', 'HLA-C*16:11', 'HLA-C*16:12', 'HLA-C*16:13',
         'HLA-C*16:14', 'HLA-C*16:15', 'HLA-C*16:17', 'HLA-C*16:18', 'HLA-C*16:19', 'HLA-C*16:20', 'HLA-C*16:21', 'HLA-C*16:22', 'HLA-C*16:23',
         'HLA-C*16:24', 'HLA-C*16:25', 'HLA-C*16:26', 'HLA-C*17:01', 'HLA-C*17:02', 'HLA-C*17:03', 'HLA-C*17:04', 'HLA-C*17:05', 'HLA-C*17:06',
         'HLA-C*17:07', 'HLA-C*18:01', 'HLA-C*18:02', 'HLA-C*18:03', 'HLA-E*01:01', 'HLA-G*01:01', 'HLA-G*01:02', 'HLA-G*01:03', 'HLA-G*01:04',
         'HLA-G*01:06', 'HLA-G*01:07', 'HLA-G*01:08', 'HLA-G*01:09',
         'H2-Db', 'H2-Dd', 'H2-Kb', 'H2-Kd', 'H2-Kk', 'H2-Ld'])
    __version = "2.4"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def supportedAlleles(self):
        """
        A list of valid :class:`~Fred2.Core.Allele.Allele` models
        """
        return self.__alleles

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    @property
    def command(self):
        """
        Defines the commandline call for external tool
        """
        return self.__command

    @property
    def supportedLength(self):
        """
        A list of supported :class:`~Fred2.Core.Peptide.Peptide` lengths
        """
        return self.__supported_length

    def _represent(self, allele):
        """
        Internal function transforming an allele object into its representative string
        :param allele: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: str
        """
        if isinstance(allele, MouseAllele):
            return "H-2-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)
        else:
            return "HLA-%s%s:%s" % (allele.locus, allele.supertype, allele.subtype)

    def convert_alleles(self, alleles):
        """
        Converts :class:`~Fred2.Core.Allele.Allele` into the internal :class:`~Fred2.Core.Allele.Allele` representation
        of the predictor and returns a string representation

        :param alleles: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`
        :rtype: list(str)
        """
        return [self._represent(a) for a in alleles]

    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        result = defaultdict(dict)
        with open(file, "r") as f:
            f = csv.reader(f, delimiter='\t')
            alleles = f.next()[3:-1]
            ic_pos = 3
            for row in f:
                pep_seq = row[1]
                for i, a in enumerate(alleles):
                    result[a][pep_seq] = float(row[ic_pos + i])
        return result

    def get_external_version(self, path=None):
        """
        Returns the external version of the tool by executing
        >{command} --version

        might be dependent on the method and has to be overwritten
        therefore it is declared abstract to enforce the user to
        overwrite the method. The function in the base class can be called
        with super()

        :param str path: Optional specification of executable path if deviant from :attr:`self.__command`
        :return: The external version of the tool or None if tool does not support versioning
        :rtype: str
        """
        # can not be determined netmhcpan does not support --version or similar
        return None

    def prepare_input(self, input, file):
        """
        Prepares input for external tools and writes them to file in the specific format

        NO return value!

        :param: list(str) input: The :class:`~Fred2.Core.Peptide.Peptide` sequences to write into file
        :param File file: File-handler to input file for external tool
        """
        file.write("\n".join(input))


class NetMHCpan_2_8(AExternalEpitopePrediction):
    """
    Implements the NetMHC binding (in current form for netMHCpan 2.8).
    Supported  MHC alleles currently only restricted to HLA alleles.

    .. note::

        Nielsen, Morten, et al. "NetMHCpan, a method for quantitative predictions of peptide binding to any HLA-A and-B
        locus protein of known sequence." PloS one 2.8 (2007): e796.
    """
    __version = "2.8"
    __supported_length = frozenset([8, 9, 10, 11, 12, 13, 14])
    __name = "netmhcpan"
    __command = "netMHCpan -p {peptides} -a {alleles} {options} -ic50 -xls -xlsfile {out}"
    __alleles = frozenset(
        ['HLA-A*01:01', 'HLA-A*01:02', 'HLA-A*01:03', 'HLA-A*01:06', 'HLA-A*01:07', 'HLA-A*01:08', 'HLA-A*01:09', 'HLA-A*01:10', 'HLA-A*01:12',
         'HLA-A*01:13', 'HLA-A*01:14', 'HLA-A*01:17', 'HLA-A*01:19', 'HLA-A*01:20', 'HLA-A*01:21', 'HLA-A*01:23', 'HLA-A*01:24', 'HLA-A*01:25',
         'HLA-A*01:26', 'HLA-A*01:28', 'HLA-A*01:29', 'HLA-A*01:30', 'HLA-A*01:32', 'HLA-A*01:33', 'HLA-A*01:35', 'HLA-A*01:36', 'HLA-A*01:37',
         'HLA-A*01:38', 'HLA-A*01:39', 'HLA-A*01:40', 'HLA-A*01:41', 'HLA-A*01:42', 'HLA-A*01:43', 'HLA-A*01:44', 'HLA-A*01:45', 'HLA-A*01:46',
         'HLA-A*01:47', 'HLA-A*01:48', 'HLA-A*01:49', 'HLA-A*01:50', 'HLA-A*01:51', 'HLA-A*01:54', 'HLA-A*01:55', 'HLA-A*01:58', 'HLA-A*01:59',
         'HLA-A*01:60', 'HLA-A*01:61', 'HLA-A*01:62', 'HLA-A*01:63', 'HLA-A*01:64', 'HLA-A*01:65', 'HLA-A*01:66', 'HLA-A*02:01', 'HLA-A*02:02',
         'HLA-A*02:03', 'HLA-A*02:04', 'HLA-A*02:05', 'HLA-A*02:06', 'HLA-A*02:07', 'HLA-A*02:08', 'HLA-A*02:09', 'HLA-A*02:10', 'HLA-A*02:101',
         'HLA-A*02:102', 'HLA-A*02:103', 'HLA-A*02:104', 'HLA-A*02:105', 'HLA-A*02:106', 'HLA-A*02:107', 'HLA-A*02:108', 'HLA-A*02:109',
         'HLA-A*02:11', 'HLA-A*02:110', 'HLA-A*02:111', 'HLA-A*02:112', 'HLA-A*02:114', 'HLA-A*02:115', 'HLA-A*02:116', 'HLA-A*02:117',
         'HLA-A*02:118', 'HLA-A*02:119', 'HLA-A*02:12', 'HLA-A*02:120', 'HLA-A*02:121', 'HLA-A*02:122', 'HLA-A*02:123', 'HLA-A*02:124',
         'HLA-A*02:126', 'HLA-A*02:127', 'HLA-A*02:128', 'HLA-A*02:129', 'HLA-A*02:13', 'HLA-A*02:130', 'HLA-A*02:131', 'HLA-A*02:132',
         'HLA-A*02:133', 'HLA-A*02:134', 'HLA-A*02:135', 'HLA-A*02:136', 'HLA-A*02:137', 'HLA-A*02:138', 'HLA-A*02:139', 'HLA-A*02:14',
         'HLA-A*02:140', 'HLA-A*02:141', 'HLA-A*02:142', 'HLA-A*02:143', 'HLA-A*02:144', 'HLA-A*02:145', 'HLA-A*02:146', 'HLA-A*02:147',
         'HLA-A*02:148', 'HLA-A*02:149', 'HLA-A*02:150', 'HLA-A*02:151', 'HLA-A*02:152', 'HLA-A*02:153', 'HLA-A*02:154', 'HLA-A*02:155',
         'HLA-A*02:156', 'HLA-A*02:157', 'HLA-A*02:158', 'HLA-A*02:159', 'HLA-A*02:16', 'HLA-A*02:160', 'HLA-A*02:161', 'HLA-A*02:162',
         'HLA-A*02:163', 'HLA-A*02:164', 'HLA-A*02:165', 'HLA-A*02:166', 'HLA-A*02:167', 'HLA-A*02:168', 'HLA-A*02:169', 'HLA-A*02:17',
         'HLA-A*02:170', 'HLA-A*02:171', 'HLA-A*02:172', 'HLA-A*02:173', 'HLA-A*02:174', 'HLA-A*02:175', 'HLA-A*02:176', 'HLA-A*02:177',
         'HLA-A*02:178', 'HLA-A*02:179', 'HLA-A*02:18', 'HLA-A*02:180', 'HLA-A*02:181', 'HLA-A*02:182', 'HLA-A*02:183', 'HLA-A*02:184',
         'HLA-A*02:185', 'HLA-A*02:186', 'HLA-A*02:187', 'HLA-A*02:188', 'HLA-A*02:189', 'HLA-A*02:19', 'HLA-A*02:190', 'HLA-A*02:191',
         'HLA-A*02:192', 'HLA-A*02:193', 'HLA-A*02:194', 'HLA-A*02:195', 'HLA-A*02:196', 'HLA-A*02:197', 'HLA-A*02:198', 'HLA-A*02:199',
         'HLA-A*02:20', 'HLA-A*02:200', 'HLA-A*02:201', 'HLA-A*02:202', 'HLA-A*02:203', 'HLA-A*02:204', 'HLA-A*02:205', 'HLA-A*02:206',
         'HLA-A*02:207', 'HLA-A*02:208', 'HLA-A*02:209', 'HLA-A*02:21', 'HLA-A*02:210', 'HLA-A*02:211', 'HLA-A*02:212', 'HLA-A*02:213',
         'HLA-A*02:214', 'HLA-A*02:215', 'HLA-A*02:216', 'HLA-A*02:217', 'HLA-A*02:218', 'HLA-A*02:219', 'HLA-A*02:22', 'HLA-A*02:220',
         'HLA-A*02:221', 'HLA-A*02:224', 'HLA-A*02:228', 'HLA-A*02:229', 'HLA-A*02:230', 'HLA-A*02:231', 'HLA-A*02:232', 'HLA-A*02:233',
         'HLA-A*02:234', 'HLA-A*02:235', 'HLA-A*02:236', 'HLA-A*02:237', 'HLA-A*02:238', 'HLA-A*02:239', 'HLA-A*02:24', 'HLA-A*02:240',
         'HLA-A*02:241', 'HLA-A*02:242', 'HLA-A*02:243', 'HLA-A*02:244', 'HLA-A*02:245', 'HLA-A*02:246', 'HLA-A*02:247', 'HLA-A*02:248',
         'HLA-A*02:249', 'HLA-A*02:25', 'HLA-A*02:251', 'HLA-A*02:252', 'HLA-A*02:253', 'HLA-A*02:254', 'HLA-A*02:255', 'HLA-A*02:256',
         'HLA-A*02:257', 'HLA-A*02:258', 'HLA-A*02:259', 'HLA-A*02:26', 'HLA-A*02:260', 'HLA-A*02:261', 'HLA-A*02:262', 'HLA-A*02:263',
         'HLA-A*02:264', 'HLA-A*02:265', 'HLA-A*02:266', 'HLA-A*02:27', 'HLA-A*02:28', 'HLA-A*02:29', 'HLA-A*02:30', 'HLA-A*02:31', 'HLA-A*02:33',
         'HLA-A*02:34', 'HLA-A*02:35', 'HLA-A*02:36', 'HLA-A*02:37', 'HLA-A*02:38', 'HLA-A*02:39', 'HLA-A*02:40', 'HLA-A*02:41', 'HLA-A*02:42',
         'HLA-A*02:44', 'HLA-A*02:45', 'HLA-A*02:46', 'HLA-A*02:47', 'HLA-A*02:48', 'HLA-A*02:49', 'HLA-A*02:50', 'HLA-A*02:51', 'HLA-A*02:52',
         'HLA-A*02:54', 'HLA-A*02:55', 'HLA-A*02:56', 'HLA-A*02:57', 'HLA-A*02:58', 'HLA-A*02:59', 'HLA-A*02:60', 'HLA-A*02:61', 'HLA-A*02:62',
         'HLA-A*02:63', 'HLA-A*02:64', 'HLA-A*02:65', 'HLA-A*02:66', 'HLA-A*02:67', 'HLA-A*02:68', 'HLA-A*02:69', 'HLA-A*02:70', 'HLA-A*02:71',
         'HLA-A*02:72', 'HLA-A*02:73', 'HLA-A*02:74', 'HLA-A*02:75', 'HLA-A*02:76', 'HLA-A*02:77', 'HLA-A*02:78', 'HLA-A*02:79', 'HLA-A*02:80',
         'HLA-A*02:81', 'HLA-A*02:84', 'HLA-A*02:85', 'HLA-A*02:86', 'HLA-A*02:87', 'HLA-A*02:89', 'HLA-A*02:90', 'HLA-A*02:91', 'HLA-A*02:92',
         'HLA-A*02:93', 'HLA-A*02:95', 'HLA-A*02:96', 'HLA-A*02:97', 'HLA-A*02:99', 'HLA-A*03:01', 'HLA-A*03:02', 'HLA-A*03:04', 'HLA-A*03:05',
         'HLA-A*03:06', 'HLA-A*03:07', 'HLA-A*03:08', 'HLA-A*03:09', 'HLA-A*03:10', 'HLA-A*03:12', 'HLA-A*03:13', 'HLA-A*03:14', 'HLA-A*03:15',
         'HLA-A*03:16', 'HLA-A*03:17', 'HLA-A*03:18', 'HLA-A*03:19', 'HLA-A*03:20', 'HLA-A*03:22', 'HLA-A*03:23', 'HLA-A*03:24', 'HLA-A*03:25',
         'HLA-A*03:26', 'HLA-A*03:27', 'HLA-A*03:28', 'HLA-A*03:29', 'HLA-A*03:30', 'HLA-A*03:31', 'HLA-A*03:32', 'HLA-A*03:33', 'HLA-A*03:34',
         'HLA-A*03:35', 'HLA-A*03:37', 'HLA-A*03:38', 'HLA-A*03:39', 'HLA-A*03:40', 'HLA-A*03:41', 'HLA-A*03:42', 'HLA-A*03:43', 'HLA-A*03:44',
         'HLA-A*03:45', 'HLA-A*03:46', 'HLA-A*03:47', 'HLA-A*03:48', 'HLA-A*03:49', 'HLA-A*03:50', 'HLA-A*03:51', 'HLA-A*03:52', 'HLA-A*03:53',
         'HLA-A*03:54', 'HLA-A*03:55', 'HLA-A*03:56', 'HLA-A*03:57', 'HLA-A*03:58', 'HLA-A*03:59', 'HLA-A*03:60', 'HLA-A*03:61', 'HLA-A*03:62',
         'HLA-A*03:63', 'HLA-A*03:64', 'HLA-A*03:65', 'HLA-A*03:66', 'HLA-A*03:67', 'HLA-A*03:70', 'HLA-A*03:71', 'HLA-A*03:72', 'HLA-A*03:73',
         'HLA-A*03:74', 'HLA-A*03:75', 'HLA-A*03:76', 'HLA-A*03:77', 'HLA-A*03:78', 'HLA-A*03:79', 'HLA-A*03:80', 'HLA-A*03:81', 'HLA-A*03:82',
         'HLA-A*11:01', 'HLA-A*11:02', 'HLA-A*11:03', 'HLA-A*11:04', 'HLA-A*11:05', 'HLA-A*11:06', 'HLA-A*11:07', 'HLA-A*11:08', 'HLA-A*11:09',
         'HLA-A*11:10', 'HLA-A*11:11', 'HLA-A*11:12', 'HLA-A*11:13', 'HLA-A*11:14', 'HLA-A*11:15', 'HLA-A*11:16', 'HLA-A*11:17', 'HLA-A*11:18',
         'HLA-A*11:19', 'HLA-A*11:20', 'HLA-A*11:22', 'HLA-A*11:23', 'HLA-A*11:24', 'HLA-A*11:25', 'HLA-A*11:26', 'HLA-A*11:27', 'HLA-A*11:29',
         'HLA-A*11:30', 'HLA-A*11:31', 'HLA-A*11:32', 'HLA-A*11:33', 'HLA-A*11:34', 'HLA-A*11:35', 'HLA-A*11:36', 'HLA-A*11:37', 'HLA-A*11:38',
         'HLA-A*11:39', 'HLA-A*11:40', 'HLA-A*11:41', 'HLA-A*11:42', 'HLA-A*11:43', 'HLA-A*11:44', 'HLA-A*11:45', 'HLA-A*11:46', 'HLA-A*11:47',
         'HLA-A*11:48', 'HLA-A*11:49', 'HLA-A*11:51', 'HLA-A*11:53', 'HLA-A*11:54', 'HLA-A*11:55', 'HLA-A*11:56', 'HLA-A*11:57', 'HLA-A*11:58',
         'HLA-A*11:59', 'HLA-A*11:60', 'HLA-A*11:61', 'HLA-A*11:62', 'HLA-A*11:63', 'HLA-A*11:64', 'HLA-A*23:01', 'HLA-A*23:02', 'HLA-A*23:03',
         'HLA-A*23:04', 'HLA-A*23:05', 'HLA-A*23:06', 'HLA-A*23:09', 'HLA-A*23:10', 'HLA-A*23:12', 'HLA-A*23:13', 'HLA-A*23:14', 'HLA-A*23:15',
         'HLA-A*23:16', 'HLA-A*23:17', 'HLA-A*23:18', 'HLA-A*23:20', 'HLA-A*23:21', 'HLA-A*23:22', 'HLA-A*23:23', 'HLA-A*23:24', 'HLA-A*23:25',
         'HLA-A*23:26', 'HLA-A*24:02', 'HLA-A*24:03', 'HLA-A*24:04', 'HLA-A*24:05', 'HLA-A*24:06', 'HLA-A*24:07', 'HLA-A*24:08', 'HLA-A*24:10',
         'HLA-A*24:100', 'HLA-A*24:101', 'HLA-A*24:102', 'HLA-A*24:103', 'HLA-A*24:104', 'HLA-A*24:105', 'HLA-A*24:106', 'HLA-A*24:107',
         'HLA-A*24:108', 'HLA-A*24:109', 'HLA-A*24:110', 'HLA-A*24:111', 'HLA-A*24:112', 'HLA-A*24:113', 'HLA-A*24:114', 'HLA-A*24:115',
         'HLA-A*24:116', 'HLA-A*24:117', 'HLA-A*24:118', 'HLA-A*24:119', 'HLA-A*24:120', 'HLA-A*24:121', 'HLA-A*24:122', 'HLA-A*24:123',
         'HLA-A*24:124', 'HLA-A*24:125', 'HLA-A*24:126', 'HLA-A*24:127', 'HLA-A*24:128', 'HLA-A*24:129', 'HLA-A*24:13', 'HLA-A*24:130',
         'HLA-A*24:131', 'HLA-A*24:133', 'HLA-A*24:134', 'HLA-A*24:135', 'HLA-A*24:136', 'HLA-A*24:137', 'HLA-A*24:138', 'HLA-A*24:139',
         'HLA-A*24:14', 'HLA-A*24:140', 'HLA-A*24:141', 'HLA-A*24:142', 'HLA-A*24:143', 'HLA-A*24:144', 'HLA-A*24:15', 'HLA-A*24:17', 'HLA-A*24:18',
         'HLA-A*24:19', 'HLA-A*24:20', 'HLA-A*24:21', 'HLA-A*24:22', 'HLA-A*24:23', 'HLA-A*24:24', 'HLA-A*24:25', 'HLA-A*24:26', 'HLA-A*24:27',
         'HLA-A*24:28', 'HLA-A*24:29', 'HLA-A*24:30', 'HLA-A*24:31', 'HLA-A*24:32', 'HLA-A*24:33', 'HLA-A*24:34', 'HLA-A*24:35', 'HLA-A*24:37',
         'HLA-A*24:38', 'HLA-A*24:39', 'HLA-A*24:41', 'HLA-A*24:42', 'HLA-A*24:43', 'HLA-A*24:44', 'HLA-A*24:46', 'HLA-A*24:47', 'HLA-A*24:49',
         'HLA-A*24:50', 'HLA-A*24:51', 'HLA-A*24:52', 'HLA-A*24:53', 'HLA-A*24:54', 'HLA-A*24:55', 'HLA-A*24:56', 'HLA-A*24:57', 'HLA-A*24:58',
         'HLA-A*24:59', 'HLA-A*24:61', 'HLA-A*24:62', 'HLA-A*24:63', 'HLA-A*24:64', 'HLA-A*24:66', 'HLA-A*24:67', 'HLA-A*24:68', 'HLA-A*24:69',
         'HLA-A*24:70', 'HLA-A*24:71', 'HLA-A*24:72', 'HLA-A*24:73', 'HLA-A*24:74', 'HLA-A*24:75', 'HLA-A*24:76', 'HLA-A*24:77', 'HLA-A*24:78',
         'HLA-A*24:79', 'HLA-A*24:80', 'HLA-A*24:81', 'HLA-A*24:82', 'HLA-A*24:85', 'HLA-A*24:87', 'HLA-A*24:88', 'HLA-A*24:89', 'HLA-A*24:91',
         'HLA-A*24:92', 'HLA-A*24:93', 'HLA-A*24:94', 'HLA-A*24:95', 'HLA-A*24:96', 'HLA-A*24:97', 'HLA-A*24:98', 'HLA-A*24:99', 'HLA-A*25:01',
         'HLA-A*25:02', 'HLA-A*25:03', 'HLA-A*25:04', 'HLA-A*25:05', 'HLA-A*25:06', 'HLA-A*25:07', 'HLA-A*25:08', 'HLA-A*25:09', 'HLA-A*25:10',
         'HLA-A*25:11', 'HLA-A*25:13', 'HLA-A*26:01', 'HLA-A*26:02', 'HLA-A*26:03', 'HLA-A*26:04', 'HLA-A*26:05', 'HLA-A*26:06', 'HLA-A*26:07',
         'HLA-A*26:08', 'HLA-A*26:09', 'HLA-A*26:10', 'HLA-A*26:12', 'HLA-A*26:13', 'HLA-A*26:14', 'HLA-A*26:15', 'HLA-A*26:16', 'HLA-A*26:17',
         'HLA-A*26:18', 'HLA-A*26:19', 'HLA-A*26:20', 'HLA-A*26:21', 'HLA-A*26:22', 'HLA-A*26:23', 'HLA-A*26:24', 'HLA-A*26:26', 'HLA-A*26:27',
         'HLA-A*26:28', 'HLA-A*26:29', 'HLA-A*26:30', 'HLA-A*26:31', 'HLA-A*26:32', 'HLA-A*26:33', 'HLA-A*26:34', 'HLA-A*26:35', 'HLA-A*26:36',
         'HLA-A*26:37', 'HLA-A*26:38', 'HLA-A*26:39', 'HLA-A*26:40', 'HLA-A*26:41', 'HLA-A*26:42', 'HLA-A*26:43', 'HLA-A*26:45', 'HLA-A*26:46',
         'HLA-A*26:47', 'HLA-A*26:48', 'HLA-A*26:49', 'HLA-A*26:50', 'HLA-A*29:01', 'HLA-A*29:02', 'HLA-A*29:03', 'HLA-A*29:04', 'HLA-A*29:05',
         'HLA-A*29:06', 'HLA-A*29:07', 'HLA-A*29:09', 'HLA-A*29:10', 'HLA-A*29:11', 'HLA-A*29:12', 'HLA-A*29:13', 'HLA-A*29:14', 'HLA-A*29:15',
         'HLA-A*29:16', 'HLA-A*29:17', 'HLA-A*29:18', 'HLA-A*29:19', 'HLA-A*29:20', 'HLA-A*29:21', 'HLA-A*29:22', 'HLA-A*30:01', 'HLA-A*30:02',
         'HLA-A*30:03', 'HLA-A*30:04', 'HLA-A*30:06', 'HLA-A*30:07', 'HLA-A*30:08', 'HLA-A*30:09', 'HLA-A*30:10', 'HLA-A*30:11', 'HLA-A*30:12',
         'HLA-A*30:13', 'HLA-A*30:15', 'HLA-A*30:16', 'HLA-A*30:17', 'HLA-A*30:18', 'HLA-A*30:19', 'HLA-A*30:20', 'HLA-A*30:22', 'HLA-A*30:23',
         'HLA-A*30:24', 'HLA-A*30:25', 'HLA-A*30:26', 'HLA-A*30:28', 'HLA-A*30:29', 'HLA-A*30:30', 'HLA-A*30:31', 'HLA-A*30:32', 'HLA-A*30:33',
         'HLA-A*30:34', 'HLA-A*30:35', 'HLA-A*30:36', 'HLA-A*30:37', 'HLA-A*30:38', 'HLA-A*30:39', 'HLA-A*30:40', 'HLA-A*30:41', 'HLA-A*31:01',
         'HLA-A*31:02', 'HLA-A*31:03', 'HLA-A*31:04', 'HLA-A*31:05', 'HLA-A*31:06', 'HLA-A*31:07', 'HLA-A*31:08', 'HLA-A*31:09', 'HLA-A*31:10',
         'HLA-A*31:11', 'HLA-A*31:12', 'HLA-A*31:13', 'HLA-A*31:15', 'HLA-A*31:16', 'HLA-A*31:17', 'HLA-A*31:18', 'HLA-A*31:19', 'HLA-A*31:20',
         'HLA-A*31:21', 'HLA-A*31:22', 'HLA-A*31:23', 'HLA-A*31:24', 'HLA-A*31:25', 'HLA-A*31:26', 'HLA-A*31:27', 'HLA-A*31:28', 'HLA-A*31:29',
         'HLA-A*31:30', 'HLA-A*31:31', 'HLA-A*31:32', 'HLA-A*31:33', 'HLA-A*31:34', 'HLA-A*31:35', 'HLA-A*31:36', 'HLA-A*31:37', 'HLA-A*32:01',
         'HLA-A*32:02', 'HLA-A*32:03', 'HLA-A*32:04', 'HLA-A*32:05', 'HLA-A*32:06', 'HLA-A*32:07', 'HLA-A*32:08', 'HLA-A*32:09', 'HLA-A*32:10',
         'HLA-A*32:12', 'HLA-A*32:13', 'HLA-A*32:14', 'HLA-A*32:15', 'HLA-A*32:16', 'HLA-A*32:17', 'HLA-A*32:18', 'HLA-A*32:20', 'HLA-A*32:21',
         'HLA-A*32:22', 'HLA-A*32:23', 'HLA-A*32:24', 'HLA-A*32:25', 'HLA-A*33:01', 'HLA-A*33:03', 'HLA-A*33:04', 'HLA-A*33:05', 'HLA-A*33:06',
         'HLA-A*33:07', 'HLA-A*33:08', 'HLA-A*33:09', 'HLA-A*33:10', 'HLA-A*33:11', 'HLA-A*33:12', 'HLA-A*33:13', 'HLA-A*33:14', 'HLA-A*33:15',
         'HLA-A*33:16', 'HLA-A*33:17', 'HLA-A*33:18', 'HLA-A*33:19', 'HLA-A*33:20', 'HLA-A*33:21', 'HLA-A*33:22', 'HLA-A*33:23', 'HLA-A*33:24',
         'HLA-A*33:25', 'HLA-A*33:26', 'HLA-A*33:27', 'HLA-A*33:28', 'HLA-A*33:29', 'HLA-A*33:30', 'HLA-A*33:31', 'HLA-A*34:01', 'HLA-A*34:02',
         'HLA-A*34:03', 'HLA-A*34:04', 'HLA-A*34:05', 'HLA-A*34:06', 'HLA-A*34:07', 'HLA-A*34:08', 'HLA-A*36:01', 'HLA-A*36:02', 'HLA-A*36:03',
         'HLA-A*36:04', 'HLA-A*36:05', 'HLA-A*43:01', 'HLA-A*66:01', 'HLA-A*66:02', 'HLA-A*66:03', 'HLA-A*66:04', 'HLA-A*66:05', 'HLA-A*66:06',
         'HLA-A*66:07', 'HLA-A*66:08', 'HLA-A*66:09', 'HLA-A*66:10', 'HLA-A*66:11', 'HLA-A*66:12', 'HLA-A*66:13', 'HLA-A*66:14', 'HLA-A*66:15',
         'HLA-A*68:01', 'HLA-A*68:02', 'HLA-A*68:03', 'HLA-A*68:04', 'HLA-A*68:05', 'HLA-A*68:06', 'HLA-A*68:07', 'HLA-A*68:08', 'HLA-A*68:09',
         'HLA-A*68:10', 'HLA-A*68:12', 'HLA-A*68:13', 'HLA-A*68:14', 'HLA-A*68:15', 'HLA-A*68:16', 'HLA-A*68:17', 'HLA-A*68:19', 'HLA-A*68:20',
         'HLA-A*68:21', 'HLA-A*68:22', 'HLA-A*68:23', 'HLA-A*68:24', 'HLA-A*68:25', 'HLA-A*68:26', 'HLA-A*68:27', 'HLA-A*68:28', 'HLA-A*68:29',
         'HLA-A*68:30', 'HLA-A*68:31', 'HLA-A*68:32', 'HLA-A*68:33', 'HLA-A*68:34', 'HLA-A*68:35', 'HLA-A*68:36', 'HLA-A*68:37', 'HLA-A*68:38',
         'HLA-A*68:39', 'HLA-A*68:40', 'HLA-A*68:41', 'HLA-A*68:42', 'HLA-A*68:43', 'HLA-A*68:44', 'HLA-A*68:45', 'HLA-A*68:46', 'HLA-A*68:47',
         'HLA-A*68:48', 'HLA-A*68:50', 'HLA-A*68:51', 'HLA-A*68:52', 'HLA-A*68:53', 'HLA-A*68:54', 'HLA-A*69:01', 'HLA-A*74:01', 'HLA-A*74:02',
         'HLA-A*74:03', 'HLA-A*74:04', 'HLA-A*74:05', 'HLA-A*74:06', 'HLA-A*74:07', 'HLA-A*74:08', 'HLA-A*74:09', 'HLA-A*74:10', 'HLA-A*74:11',
         'HLA-A*74:13', 'HLA-A*80:01', 'HLA-A*80:02', 'HLA-B*07:02', 'HLA-B*07:03', 'HLA-B*07:04', 'HLA-B*07:05', 'HLA-B*07:06', 'HLA-B*07:07',
         'HLA-B*07:08', 'HLA-B*07:09', 'HLA-B*07:10', 'HLA-B*07:100', 'HLA-B*07:101', 'HLA-B*07:102', 'HLA-B*07:103', 'HLA-B*07:104',
         'HLA-B*07:105', 'HLA-B*07:106', 'HLA-B*07:107', 'HLA-B*07:108', 'HLA-B*07:109', 'HLA-B*07:11', 'HLA-B*07:110', 'HLA-B*07:112',
         'HLA-B*07:113', 'HLA-B*07:114', 'HLA-B*07:115', 'HLA-B*07:12', 'HLA-B*07:13', 'HLA-B*07:14', 'HLA-B*07:15', 'HLA-B*07:16', 'HLA-B*07:17',
         'HLA-B*07:18', 'HLA-B*07:19', 'HLA-B*07:20', 'HLA-B*07:21', 'HLA-B*07:22', 'HLA-B*07:23', 'HLA-B*07:24', 'HLA-B*07:25', 'HLA-B*07:26',
         'HLA-B*07:27', 'HLA-B*07:28', 'HLA-B*07:29', 'HLA-B*07:30', 'HLA-B*07:31', 'HLA-B*07:32', 'HLA-B*07:33', 'HLA-B*07:34', 'HLA-B*07:35',
         'HLA-B*07:36', 'HLA-B*07:37', 'HLA-B*07:38', 'HLA-B*07:39', 'HLA-B*07:40', 'HLA-B*07:41', 'HLA-B*07:42', 'HLA-B*07:43', 'HLA-B*07:44',
         'HLA-B*07:45', 'HLA-B*07:46', 'HLA-B*07:47', 'HLA-B*07:48', 'HLA-B*07:50', 'HLA-B*07:51', 'HLA-B*07:52', 'HLA-B*07:53', 'HLA-B*07:54',
         'HLA-B*07:55', 'HLA-B*07:56', 'HLA-B*07:57', 'HLA-B*07:58', 'HLA-B*07:59', 'HLA-B*07:60', 'HLA-B*07:61', 'HLA-B*07:62', 'HLA-B*07:63',
         'HLA-B*07:64', 'HLA-B*07:65', 'HLA-B*07:66', 'HLA-B*07:68', 'HLA-B*07:69', 'HLA-B*07:70', 'HLA-B*07:71', 'HLA-B*07:72', 'HLA-B*07:73',
         'HLA-B*07:74', 'HLA-B*07:75', 'HLA-B*07:76', 'HLA-B*07:77', 'HLA-B*07:78', 'HLA-B*07:79', 'HLA-B*07:80', 'HLA-B*07:81', 'HLA-B*07:82',
         'HLA-B*07:83', 'HLA-B*07:84', 'HLA-B*07:85', 'HLA-B*07:86', 'HLA-B*07:87', 'HLA-B*07:88', 'HLA-B*07:89', 'HLA-B*07:90', 'HLA-B*07:91',
         'HLA-B*07:92', 'HLA-B*07:93', 'HLA-B*07:94', 'HLA-B*07:95', 'HLA-B*07:96', 'HLA-B*07:97', 'HLA-B*07:98', 'HLA-B*07:99', 'HLA-B*08:01',
         'HLA-B*08:02', 'HLA-B*08:03', 'HLA-B*08:04', 'HLA-B*08:05', 'HLA-B*08:07', 'HLA-B*08:09', 'HLA-B*08:10', 'HLA-B*08:11', 'HLA-B*08:12',
         'HLA-B*08:13', 'HLA-B*08:14', 'HLA-B*08:15', 'HLA-B*08:16', 'HLA-B*08:17', 'HLA-B*08:18', 'HLA-B*08:20', 'HLA-B*08:21', 'HLA-B*08:22',
         'HLA-B*08:23', 'HLA-B*08:24', 'HLA-B*08:25', 'HLA-B*08:26', 'HLA-B*08:27', 'HLA-B*08:28', 'HLA-B*08:29', 'HLA-B*08:31', 'HLA-B*08:32',
         'HLA-B*08:33', 'HLA-B*08:34', 'HLA-B*08:35', 'HLA-B*08:36', 'HLA-B*08:37', 'HLA-B*08:38', 'HLA-B*08:39', 'HLA-B*08:40', 'HLA-B*08:41',
         'HLA-B*08:42', 'HLA-B*08:43', 'HLA-B*08:44', 'HLA-B*08:45', 'HLA-B*08:46', 'HLA-B*08:47', 'HLA-B*08:48', 'HLA-B*08:49', 'HLA-B*08:50',
         'HLA-B*08:51', 'HLA-B*08:52', 'HLA-B*08:53', 'HLA-B*08:54', 'HLA-B*08:55', 'HLA-B*08:56', 'HLA-B*08:57', 'HLA-B*08:58', 'HLA-B*08:59',
         'HLA-B*08:60', 'HLA-B*08:61', 'HLA-B*08:62', 'HLA-B*13:01', 'HLA-B*13:02', 'HLA-B*13:03', 'HLA-B*13:04', 'HLA-B*13:06', 'HLA-B*13:09',
         'HLA-B*13:10', 'HLA-B*13:11', 'HLA-B*13:12', 'HLA-B*13:13', 'HLA-B*13:14', 'HLA-B*13:15', 'HLA-B*13:16', 'HLA-B*13:17', 'HLA-B*13:18',
         'HLA-B*13:19', 'HLA-B*13:20', 'HLA-B*13:21', 'HLA-B*13:22', 'HLA-B*13:23', 'HLA-B*13:25', 'HLA-B*13:26', 'HLA-B*13:27', 'HLA-B*13:28',
         'HLA-B*13:29', 'HLA-B*13:30', 'HLA-B*13:31', 'HLA-B*13:32', 'HLA-B*13:33', 'HLA-B*13:34', 'HLA-B*13:35', 'HLA-B*13:36', 'HLA-B*13:37',
         'HLA-B*13:38', 'HLA-B*13:39', 'HLA-B*14:01', 'HLA-B*14:02', 'HLA-B*14:03', 'HLA-B*14:04', 'HLA-B*14:05', 'HLA-B*14:06', 'HLA-B*14:08',
         'HLA-B*14:09', 'HLA-B*14:10', 'HLA-B*14:11', 'HLA-B*14:12', 'HLA-B*14:13', 'HLA-B*14:14', 'HLA-B*14:15', 'HLA-B*14:16', 'HLA-B*14:17',
         'HLA-B*14:18', 'HLA-B*15:01', 'HLA-B*15:02', 'HLA-B*15:03', 'HLA-B*15:04', 'HLA-B*15:05', 'HLA-B*15:06', 'HLA-B*15:07', 'HLA-B*15:08',
         'HLA-B*15:09', 'HLA-B*15:10', 'HLA-B*15:101', 'HLA-B*15:102', 'HLA-B*15:103', 'HLA-B*15:104', 'HLA-B*15:105', 'HLA-B*15:106',
         'HLA-B*15:107', 'HLA-B*15:108', 'HLA-B*15:109', 'HLA-B*15:11', 'HLA-B*15:110', 'HLA-B*15:112', 'HLA-B*15:113', 'HLA-B*15:114',
         'HLA-B*15:115', 'HLA-B*15:116', 'HLA-B*15:117', 'HLA-B*15:118', 'HLA-B*15:119', 'HLA-B*15:12', 'HLA-B*15:120', 'HLA-B*15:121',
         'HLA-B*15:122', 'HLA-B*15:123', 'HLA-B*15:124', 'HLA-B*15:125', 'HLA-B*15:126', 'HLA-B*15:127', 'HLA-B*15:128', 'HLA-B*15:129',
         'HLA-B*15:13', 'HLA-B*15:131', 'HLA-B*15:132', 'HLA-B*15:133', 'HLA-B*15:134', 'HLA-B*15:135', 'HLA-B*15:136', 'HLA-B*15:137',
         'HLA-B*15:138', 'HLA-B*15:139', 'HLA-B*15:14', 'HLA-B*15:140', 'HLA-B*15:141', 'HLA-B*15:142', 'HLA-B*15:143', 'HLA-B*15:144',
         'HLA-B*15:145', 'HLA-B*15:146', 'HLA-B*15:147', 'HLA-B*15:148', 'HLA-B*15:15', 'HLA-B*15:150', 'HLA-B*15:151', 'HLA-B*15:152',
         'HLA-B*15:153', 'HLA-B*15:154', 'HLA-B*15:155', 'HLA-B*15:156', 'HLA-B*15:157', 'HLA-B*15:158', 'HLA-B*15:159', 'HLA-B*15:16',
         'HLA-B*15:160', 'HLA-B*15:161', 'HLA-B*15:162', 'HLA-B*15:163', 'HLA-B*15:164', 'HLA-B*15:165', 'HLA-B*15:166', 'HLA-B*15:167',
         'HLA-B*15:168', 'HLA-B*15:169', 'HLA-B*15:17', 'HLA-B*15:170', 'HLA-B*15:171', 'HLA-B*15:172', 'HLA-B*15:173', 'HLA-B*15:174',
         'HLA-B*15:175', 'HLA-B*15:176', 'HLA-B*15:177', 'HLA-B*15:178', 'HLA-B*15:179', 'HLA-B*15:18', 'HLA-B*15:180', 'HLA-B*15:183',
         'HLA-B*15:184', 'HLA-B*15:185', 'HLA-B*15:186', 'HLA-B*15:187', 'HLA-B*15:188', 'HLA-B*15:189', 'HLA-B*15:19', 'HLA-B*15:191',
         'HLA-B*15:192', 'HLA-B*15:193', 'HLA-B*15:194', 'HLA-B*15:195', 'HLA-B*15:196', 'HLA-B*15:197', 'HLA-B*15:198', 'HLA-B*15:199',
         'HLA-B*15:20', 'HLA-B*15:200', 'HLA-B*15:201', 'HLA-B*15:202', 'HLA-B*15:21', 'HLA-B*15:23', 'HLA-B*15:24', 'HLA-B*15:25', 'HLA-B*15:27',
         'HLA-B*15:28', 'HLA-B*15:29', 'HLA-B*15:30', 'HLA-B*15:31', 'HLA-B*15:32', 'HLA-B*15:33', 'HLA-B*15:34', 'HLA-B*15:35', 'HLA-B*15:36',
         'HLA-B*15:37', 'HLA-B*15:38', 'HLA-B*15:39', 'HLA-B*15:40', 'HLA-B*15:42', 'HLA-B*15:43', 'HLA-B*15:44', 'HLA-B*15:45', 'HLA-B*15:46',
         'HLA-B*15:47', 'HLA-B*15:48', 'HLA-B*15:49', 'HLA-B*15:50', 'HLA-B*15:51', 'HLA-B*15:52', 'HLA-B*15:53', 'HLA-B*15:54', 'HLA-B*15:55',
         'HLA-B*15:56', 'HLA-B*15:57', 'HLA-B*15:58', 'HLA-B*15:60', 'HLA-B*15:61', 'HLA-B*15:62', 'HLA-B*15:63', 'HLA-B*15:64', 'HLA-B*15:65',
         'HLA-B*15:66', 'HLA-B*15:67', 'HLA-B*15:68', 'HLA-B*15:69', 'HLA-B*15:70', 'HLA-B*15:71', 'HLA-B*15:72', 'HLA-B*15:73', 'HLA-B*15:74',
         'HLA-B*15:75', 'HLA-B*15:76', 'HLA-B*15:77', 'HLA-B*15:78', 'HLA-B*15:80', 'HLA-B*15:81', 'HLA-B*15:82', 'HLA-B*15:83', 'HLA-B*15:84',
         'HLA-B*15:85', 'HLA-B*15:86', 'HLA-B*15:87', 'HLA-B*15:88', 'HLA-B*15:89', 'HLA-B*15:90', 'HLA-B*15:91', 'HLA-B*15:92', 'HLA-B*15:93',
         'HLA-B*15:95', 'HLA-B*15:96', 'HLA-B*15:97', 'HLA-B*15:98', 'HLA-B*15:99', 'HLA-B*18:01', 'HLA-B*18:02', 'HLA-B*18:03', 'HLA-B*18:04',
         'HLA-B*18:05', 'HLA-B*18:06', 'HLA-B*18:07', 'HLA-B*18:08', 'HLA-B*18:09', 'HLA-B*18:10', 'HLA-B*18:11', 'HLA-B*18:12', 'HLA-B*18:13',
         'HLA-B*18:14', 'HLA-B*18:15', 'HLA-B*18:18', 'HLA-B*18:19', 'HLA-B*18:20', 'HLA-B*18:21', 'HLA-B*18:22', 'HLA-B*18:24', 'HLA-B*18:25',
         'HLA-B*18:26', 'HLA-B*18:27', 'HLA-B*18:28', 'HLA-B*18:29', 'HLA-B*18:30', 'HLA-B*18:31', 'HLA-B*18:32', 'HLA-B*18:33', 'HLA-B*18:34',
         'HLA-B*18:35', 'HLA-B*18:36', 'HLA-B*18:37', 'HLA-B*18:38', 'HLA-B*18:39', 'HLA-B*18:40', 'HLA-B*18:41', 'HLA-B*18:42', 'HLA-B*18:43',
         'HLA-B*18:44', 'HLA-B*18:45', 'HLA-B*18:46', 'HLA-B*18:47', 'HLA-B*18:48', 'HLA-B*18:49', 'HLA-B*18:50', 'HLA-B*27:01', 'HLA-B*27:02',
         'HLA-B*27:03', 'HLA-B*27:04', 'HLA-B*27:05', 'HLA-B*27:06', 'HLA-B*27:07', 'HLA-B*27:08', 'HLA-B*27:09', 'HLA-B*27:10', 'HLA-B*27:11',
         'HLA-B*27:12', 'HLA-B*27:13', 'HLA-B*27:14', 'HLA-B*27:15', 'HLA-B*27:16', 'HLA-B*27:17', 'HLA-B*27:18', 'HLA-B*27:19', 'HLA-B*27:20',
         'HLA-B*27:21', 'HLA-B*27:23', 'HLA-B*27:24', 'HLA-B*27:25', 'HLA-B*27:26', 'HLA-B*27:27', 'HLA-B*27:28', 'HLA-B*27:29', 'HLA-B*27:30',
         'HLA-B*27:31', 'HLA-B*27:32', 'HLA-B*27:33', 'HLA-B*27:34', 'HLA-B*27:35', 'HLA-B*27:36', 'HLA-B*27:37', 'HLA-B*27:38', 'HLA-B*27:39',
         'HLA-B*27:40', 'HLA-B*27:41', 'HLA-B*27:42', 'HLA-B*27:43', 'HLA-B*27:44', 'HLA-B*27:45', 'HLA-B*27:46', 'HLA-B*27:47', 'HLA-B*27:48',
         'HLA-B*27:49', 'HLA-B*27:50', 'HLA-B*27:51', 'HLA-B*27:52', 'HLA-B*27:53', 'HLA-B*27:54', 'HLA-B*27:55', 'HLA-B*27:56', 'HLA-B*27:57',
         'HLA-B*27:58', 'HLA-B*27:60', 'HLA-B*27:61', 'HLA-B*27:62', 'HLA-B*27:63', 'HLA-B*27:67', 'HLA-B*27:68', 'HLA-B*27:69', 'HLA-B*35:01',
         'HLA-B*35:02', 'HLA-B*35:03', 'HLA-B*35:04', 'HLA-B*35:05', 'HLA-B*35:06', 'HLA-B*35:07', 'HLA-B*35:08', 'HLA-B*35:09', 'HLA-B*35:10',
         'HLA-B*35:100', 'HLA-B*35:101', 'HLA-B*35:102', 'HLA-B*35:103', 'HLA-B*35:104', 'HLA-B*35:105', 'HLA-B*35:106', 'HLA-B*35:107',
         'HLA-B*35:108', 'HLA-B*35:109', 'HLA-B*35:11', 'HLA-B*35:110', 'HLA-B*35:111', 'HLA-B*35:112', 'HLA-B*35:113', 'HLA-B*35:114',
         'HLA-B*35:115', 'HLA-B*35:116', 'HLA-B*35:117', 'HLA-B*35:118', 'HLA-B*35:119', 'HLA-B*35:12', 'HLA-B*35:120', 'HLA-B*35:121',
         'HLA-B*35:122', 'HLA-B*35:123', 'HLA-B*35:124', 'HLA-B*35:125', 'HLA-B*35:126', 'HLA-B*35:127', 'HLA-B*35:128', 'HLA-B*35:13',
         'HLA-B*35:131', 'HLA-B*35:132', 'HLA-B*35:133', 'HLA-B*35:135', 'HLA-B*35:136', 'HLA-B*35:137', 'HLA-B*35:138', 'HLA-B*35:139',
         'HLA-B*35:14', 'HLA-B*35:140', 'HLA-B*35:141', 'HLA-B*35:142', 'HLA-B*35:143', 'HLA-B*35:144', 'HLA-B*35:15', 'HLA-B*35:16', 'HLA-B*35:17',
         'HLA-B*35:18', 'HLA-B*35:19', 'HLA-B*35:20', 'HLA-B*35:21', 'HLA-B*35:22', 'HLA-B*35:23', 'HLA-B*35:24', 'HLA-B*35:25', 'HLA-B*35:26',
         'HLA-B*35:27', 'HLA-B*35:28', 'HLA-B*35:29', 'HLA-B*35:30', 'HLA-B*35:31', 'HLA-B*35:32', 'HLA-B*35:33', 'HLA-B*35:34', 'HLA-B*35:35',
         'HLA-B*35:36', 'HLA-B*35:37', 'HLA-B*35:38', 'HLA-B*35:39', 'HLA-B*35:41', 'HLA-B*35:42', 'HLA-B*35:43', 'HLA-B*35:44', 'HLA-B*35:45',
         'HLA-B*35:46', 'HLA-B*35:47', 'HLA-B*35:48', 'HLA-B*35:49', 'HLA-B*35:50', 'HLA-B*35:51', 'HLA-B*35:52', 'HLA-B*35:54', 'HLA-B*35:55',
         'HLA-B*35:56', 'HLA-B*35:57', 'HLA-B*35:58', 'HLA-B*35:59', 'HLA-B*35:60', 'HLA-B*35:61', 'HLA-B*35:62', 'HLA-B*35:63', 'HLA-B*35:64',
         'HLA-B*35:66', 'HLA-B*35:67', 'HLA-B*35:68', 'HLA-B*35:69', 'HLA-B*35:70', 'HLA-B*35:71', 'HLA-B*35:72', 'HLA-B*35:74', 'HLA-B*35:75',
         'HLA-B*35:76', 'HLA-B*35:77', 'HLA-B*35:78', 'HLA-B*35:79', 'HLA-B*35:80', 'HLA-B*35:81', 'HLA-B*35:82', 'HLA-B*35:83', 'HLA-B*35:84',
         'HLA-B*35:85', 'HLA-B*35:86', 'HLA-B*35:87', 'HLA-B*35:88', 'HLA-B*35:89', 'HLA-B*35:90', 'HLA-B*35:91', 'HLA-B*35:92', 'HLA-B*35:93',
         'HLA-B*35:94', 'HLA-B*35:95', 'HLA-B*35:96', 'HLA-B*35:97', 'HLA-B*35:98', 'HLA-B*35:99', 'HLA-B*37:01', 'HLA-B*37:02', 'HLA-B*37:04',
         'HLA-B*37:05', 'HLA-B*37:06', 'HLA-B*37:07', 'HLA-B*37:08', 'HLA-B*37:09', 'HLA-B*37:10', 'HLA-B*37:11', 'HLA-B*37:12', 'HLA-B*37:13',
         'HLA-B*37:14', 'HLA-B*37:15', 'HLA-B*37:17', 'HLA-B*37:18', 'HLA-B*37:19', 'HLA-B*37:20', 'HLA-B*37:21', 'HLA-B*37:22', 'HLA-B*37:23',
         'HLA-B*38:01', 'HLA-B*38:02', 'HLA-B*38:03', 'HLA-B*38:04', 'HLA-B*38:05', 'HLA-B*38:06', 'HLA-B*38:07', 'HLA-B*38:08', 'HLA-B*38:09',
         'HLA-B*38:10', 'HLA-B*38:11', 'HLA-B*38:12', 'HLA-B*38:13', 'HLA-B*38:14', 'HLA-B*38:15', 'HLA-B*38:16', 'HLA-B*38:17', 'HLA-B*38:18',
         'HLA-B*38:19', 'HLA-B*38:20', 'HLA-B*38:21', 'HLA-B*38:22', 'HLA-B*38:23', 'HLA-B*39:01', 'HLA-B*39:02', 'HLA-B*39:03', 'HLA-B*39:04',
         'HLA-B*39:05', 'HLA-B*39:06', 'HLA-B*39:07', 'HLA-B*39:08', 'HLA-B*39:09', 'HLA-B*39:10', 'HLA-B*39:11', 'HLA-B*39:12', 'HLA-B*39:13',
         'HLA-B*39:14', 'HLA-B*39:15', 'HLA-B*39:16', 'HLA-B*39:17', 'HLA-B*39:18', 'HLA-B*39:19', 'HLA-B*39:20', 'HLA-B*39:22', 'HLA-B*39:23',
         'HLA-B*39:24', 'HLA-B*39:26', 'HLA-B*39:27', 'HLA-B*39:28', 'HLA-B*39:29', 'HLA-B*39:30', 'HLA-B*39:31', 'HLA-B*39:32', 'HLA-B*39:33',
         'HLA-B*39:34', 'HLA-B*39:35', 'HLA-B*39:36', 'HLA-B*39:37', 'HLA-B*39:39', 'HLA-B*39:41', 'HLA-B*39:42', 'HLA-B*39:43', 'HLA-B*39:44',
         'HLA-B*39:45', 'HLA-B*39:46', 'HLA-B*39:47', 'HLA-B*39:48', 'HLA-B*39:49', 'HLA-B*39:50', 'HLA-B*39:51', 'HLA-B*39:52', 'HLA-B*39:53',
         'HLA-B*39:54', 'HLA-B*39:55', 'HLA-B*39:56', 'HLA-B*39:57', 'HLA-B*39:58', 'HLA-B*39:59', 'HLA-B*39:60', 'HLA-B*40:01', 'HLA-B*40:02',
         'HLA-B*40:03', 'HLA-B*40:04', 'HLA-B*40:05', 'HLA-B*40:06', 'HLA-B*40:07', 'HLA-B*40:08', 'HLA-B*40:09', 'HLA-B*40:10', 'HLA-B*40:100',
         'HLA-B*40:101', 'HLA-B*40:102', 'HLA-B*40:103', 'HLA-B*40:104', 'HLA-B*40:105', 'HLA-B*40:106', 'HLA-B*40:107', 'HLA-B*40:108',
         'HLA-B*40:109', 'HLA-B*40:11', 'HLA-B*40:110', 'HLA-B*40:111', 'HLA-B*40:112', 'HLA-B*40:113', 'HLA-B*40:114', 'HLA-B*40:115',
         'HLA-B*40:116', 'HLA-B*40:117', 'HLA-B*40:119', 'HLA-B*40:12', 'HLA-B*40:120', 'HLA-B*40:121', 'HLA-B*40:122', 'HLA-B*40:123',
         'HLA-B*40:124', 'HLA-B*40:125', 'HLA-B*40:126', 'HLA-B*40:127', 'HLA-B*40:128', 'HLA-B*40:129', 'HLA-B*40:13', 'HLA-B*40:130',
         'HLA-B*40:131', 'HLA-B*40:132', 'HLA-B*40:134', 'HLA-B*40:135', 'HLA-B*40:136', 'HLA-B*40:137', 'HLA-B*40:138', 'HLA-B*40:139',
         'HLA-B*40:14', 'HLA-B*40:140', 'HLA-B*40:141', 'HLA-B*40:143', 'HLA-B*40:145', 'HLA-B*40:146', 'HLA-B*40:147', 'HLA-B*40:15',
         'HLA-B*40:16', 'HLA-B*40:18', 'HLA-B*40:19', 'HLA-B*40:20', 'HLA-B*40:21', 'HLA-B*40:23', 'HLA-B*40:24', 'HLA-B*40:25', 'HLA-B*40:26',
         'HLA-B*40:27', 'HLA-B*40:28', 'HLA-B*40:29', 'HLA-B*40:30', 'HLA-B*40:31', 'HLA-B*40:32', 'HLA-B*40:33', 'HLA-B*40:34', 'HLA-B*40:35',
         'HLA-B*40:36', 'HLA-B*40:37', 'HLA-B*40:38', 'HLA-B*40:39', 'HLA-B*40:40', 'HLA-B*40:42', 'HLA-B*40:43', 'HLA-B*40:44', 'HLA-B*40:45',
         'HLA-B*40:46', 'HLA-B*40:47', 'HLA-B*40:48', 'HLA-B*40:49', 'HLA-B*40:50', 'HLA-B*40:51', 'HLA-B*40:52', 'HLA-B*40:53', 'HLA-B*40:54',
         'HLA-B*40:55', 'HLA-B*40:56', 'HLA-B*40:57', 'HLA-B*40:58', 'HLA-B*40:59', 'HLA-B*40:60', 'HLA-B*40:61', 'HLA-B*40:62', 'HLA-B*40:63',
         'HLA-B*40:64', 'HLA-B*40:65', 'HLA-B*40:66', 'HLA-B*40:67', 'HLA-B*40:68', 'HLA-B*40:69', 'HLA-B*40:70', 'HLA-B*40:71', 'HLA-B*40:72',
         'HLA-B*40:73', 'HLA-B*40:74', 'HLA-B*40:75', 'HLA-B*40:76', 'HLA-B*40:77', 'HLA-B*40:78', 'HLA-B*40:79', 'HLA-B*40:80', 'HLA-B*40:81',
         'HLA-B*40:82', 'HLA-B*40:83', 'HLA-B*40:84', 'HLA-B*40:85', 'HLA-B*40:86', 'HLA-B*40:87', 'HLA-B*40:88', 'HLA-B*40:89', 'HLA-B*40:90',
         'HLA-B*40:91', 'HLA-B*40:92', 'HLA-B*40:93', 'HLA-B*40:94', 'HLA-B*40:95', 'HLA-B*40:96', 'HLA-B*40:97', 'HLA-B*40:98', 'HLA-B*40:99',
         'HLA-B*41:01', 'HLA-B*41:02', 'HLA-B*41:03', 'HLA-B*41:04', 'HLA-B*41:05', 'HLA-B*41:06', 'HLA-B*41:07', 'HLA-B*41:08', 'HLA-B*41:09',
         'HLA-B*41:10', 'HLA-B*41:11', 'HLA-B*41:12', 'HLA-B*42:01', 'HLA-B*42:02', 'HLA-B*42:04', 'HLA-B*42:05', 'HLA-B*42:06', 'HLA-B*42:07',
         'HLA-B*42:08', 'HLA-B*42:09', 'HLA-B*42:10', 'HLA-B*42:11', 'HLA-B*42:12', 'HLA-B*42:13', 'HLA-B*42:14', 'HLA-B*44:02', 'HLA-B*44:03',
         'HLA-B*44:04', 'HLA-B*44:05', 'HLA-B*44:06', 'HLA-B*44:07', 'HLA-B*44:08', 'HLA-B*44:09', 'HLA-B*44:10', 'HLA-B*44:100', 'HLA-B*44:101',
         'HLA-B*44:102', 'HLA-B*44:103', 'HLA-B*44:104', 'HLA-B*44:105', 'HLA-B*44:106', 'HLA-B*44:107', 'HLA-B*44:109', 'HLA-B*44:11',
         'HLA-B*44:110', 'HLA-B*44:12', 'HLA-B*44:13', 'HLA-B*44:14', 'HLA-B*44:15', 'HLA-B*44:16', 'HLA-B*44:17', 'HLA-B*44:18', 'HLA-B*44:20',
         'HLA-B*44:21', 'HLA-B*44:22', 'HLA-B*44:24', 'HLA-B*44:25', 'HLA-B*44:26', 'HLA-B*44:27', 'HLA-B*44:28', 'HLA-B*44:29', 'HLA-B*44:30',
         'HLA-B*44:31', 'HLA-B*44:32', 'HLA-B*44:33', 'HLA-B*44:34', 'HLA-B*44:35', 'HLA-B*44:36', 'HLA-B*44:37', 'HLA-B*44:38', 'HLA-B*44:39',
         'HLA-B*44:40', 'HLA-B*44:41', 'HLA-B*44:42', 'HLA-B*44:43', 'HLA-B*44:44', 'HLA-B*44:45', 'HLA-B*44:46', 'HLA-B*44:47', 'HLA-B*44:48',
         'HLA-B*44:49', 'HLA-B*44:50', 'HLA-B*44:51', 'HLA-B*44:53', 'HLA-B*44:54', 'HLA-B*44:55', 'HLA-B*44:57', 'HLA-B*44:59', 'HLA-B*44:60',
         'HLA-B*44:62', 'HLA-B*44:63', 'HLA-B*44:64', 'HLA-B*44:65', 'HLA-B*44:66', 'HLA-B*44:67', 'HLA-B*44:68', 'HLA-B*44:69', 'HLA-B*44:70',
         'HLA-B*44:71', 'HLA-B*44:72', 'HLA-B*44:73', 'HLA-B*44:74', 'HLA-B*44:75', 'HLA-B*44:76', 'HLA-B*44:77', 'HLA-B*44:78', 'HLA-B*44:79',
         'HLA-B*44:80', 'HLA-B*44:81', 'HLA-B*44:82', 'HLA-B*44:83', 'HLA-B*44:84', 'HLA-B*44:85', 'HLA-B*44:86', 'HLA-B*44:87', 'HLA-B*44:88',
         'HLA-B*44:89', 'HLA-B*44:90', 'HLA-B*44:91', 'HLA-B*44:92', 'HLA-B*44:93', 'HLA-B*44:94', 'HLA-B*44:95', 'HLA-B*44:96', 'HLA-B*44:97',
         'HLA-B*44:98', 'HLA-B*44:99', 'HLA-B*45:01', 'HLA-B*45:02', 'HLA-B*45:03', 'HLA-B*45:04', 'HLA-B*45:05', 'HLA-B*45:06', 'HLA-B*45:07',
         'HLA-B*45:08', 'HLA-B*45:09', 'HLA-B*45:10', 'HLA-B*45:11', 'HLA-B*45:12', 'HLA-B*46:01', 'HLA-B*46:02', 'HLA-B*46:03', 'HLA-B*46:04',
         'HLA-B*46:05', 'HLA-B*46:06', 'HLA-B*46:08', 'HLA-B*46:09', 'HLA-B*46:10', 'HLA-B*46:11', 'HLA-B*46:12', 'HLA-B*46:13', 'HLA-B*46:14',
         'HLA-B*46:16', 'HLA-B*46:17', 'HLA-B*46:18', 'HLA-B*46:19', 'HLA-B*46:20', 'HLA-B*46:21', 'HLA-B*46:22', 'HLA-B*46:23', 'HLA-B*46:24',
         'HLA-B*47:01', 'HLA-B*47:02', 'HLA-B*47:03', 'HLA-B*47:04', 'HLA-B*47:05', 'HLA-B*47:06', 'HLA-B*47:07', 'HLA-B*48:01', 'HLA-B*48:02',
         'HLA-B*48:03', 'HLA-B*48:04', 'HLA-B*48:05', 'HLA-B*48:06', 'HLA-B*48:07', 'HLA-B*48:08', 'HLA-B*48:09', 'HLA-B*48:10', 'HLA-B*48:11',
         'HLA-B*48:12', 'HLA-B*48:13', 'HLA-B*48:14', 'HLA-B*48:15', 'HLA-B*48:16', 'HLA-B*48:17', 'HLA-B*48:18', 'HLA-B*48:19', 'HLA-B*48:20',
         'HLA-B*48:21', 'HLA-B*48:22', 'HLA-B*48:23', 'HLA-B*49:01', 'HLA-B*49:02', 'HLA-B*49:03', 'HLA-B*49:04', 'HLA-B*49:05', 'HLA-B*49:06',
         'HLA-B*49:07', 'HLA-B*49:08', 'HLA-B*49:09', 'HLA-B*49:10', 'HLA-B*50:01', 'HLA-B*50:02', 'HLA-B*50:04', 'HLA-B*50:05', 'HLA-B*50:06',
         'HLA-B*50:07', 'HLA-B*50:08', 'HLA-B*50:09', 'HLA-B*51:01', 'HLA-B*51:02', 'HLA-B*51:03', 'HLA-B*51:04', 'HLA-B*51:05', 'HLA-B*51:06',
         'HLA-B*51:07', 'HLA-B*51:08', 'HLA-B*51:09', 'HLA-B*51:12', 'HLA-B*51:13', 'HLA-B*51:14', 'HLA-B*51:15', 'HLA-B*51:16', 'HLA-B*51:17',
         'HLA-B*51:18', 'HLA-B*51:19', 'HLA-B*51:20', 'HLA-B*51:21', 'HLA-B*51:22', 'HLA-B*51:23', 'HLA-B*51:24', 'HLA-B*51:26', 'HLA-B*51:28',
         'HLA-B*51:29', 'HLA-B*51:30', 'HLA-B*51:31', 'HLA-B*51:32', 'HLA-B*51:33', 'HLA-B*51:34', 'HLA-B*51:35', 'HLA-B*51:36', 'HLA-B*51:37',
         'HLA-B*51:38', 'HLA-B*51:39', 'HLA-B*51:40', 'HLA-B*51:42', 'HLA-B*51:43', 'HLA-B*51:45', 'HLA-B*51:46', 'HLA-B*51:48', 'HLA-B*51:49',
         'HLA-B*51:50', 'HLA-B*51:51', 'HLA-B*51:52', 'HLA-B*51:53', 'HLA-B*51:54', 'HLA-B*51:55', 'HLA-B*51:56', 'HLA-B*51:57', 'HLA-B*51:58',
         'HLA-B*51:59', 'HLA-B*51:60', 'HLA-B*51:61', 'HLA-B*51:62', 'HLA-B*51:63', 'HLA-B*51:64', 'HLA-B*51:65', 'HLA-B*51:66', 'HLA-B*51:67',
         'HLA-B*51:68', 'HLA-B*51:69', 'HLA-B*51:70', 'HLA-B*51:71', 'HLA-B*51:72', 'HLA-B*51:73', 'HLA-B*51:74', 'HLA-B*51:75', 'HLA-B*51:76',
         'HLA-B*51:77', 'HLA-B*51:78', 'HLA-B*51:79', 'HLA-B*51:80', 'HLA-B*51:81', 'HLA-B*51:82', 'HLA-B*51:83', 'HLA-B*51:84', 'HLA-B*51:85',
         'HLA-B*51:86', 'HLA-B*51:87', 'HLA-B*51:88', 'HLA-B*51:89', 'HLA-B*51:90', 'HLA-B*51:91', 'HLA-B*51:92', 'HLA-B*51:93', 'HLA-B*51:94',
         'HLA-B*51:95', 'HLA-B*51:96', 'HLA-B*52:01', 'HLA-B*52:02', 'HLA-B*52:03', 'HLA-B*52:04', 'HLA-B*52:05', 'HLA-B*52:06', 'HLA-B*52:07',
         'HLA-B*52:08', 'HLA-B*52:09', 'HLA-B*52:10', 'HLA-B*52:11', 'HLA-B*52:12', 'HLA-B*52:13', 'HLA-B*52:14', 'HLA-B*52:15', 'HLA-B*52:16',
         'HLA-B*52:17', 'HLA-B*52:18', 'HLA-B*52:19', 'HLA-B*52:20', 'HLA-B*52:21', 'HLA-B*53:01', 'HLA-B*53:02', 'HLA-B*53:03', 'HLA-B*53:04',
         'HLA-B*53:05', 'HLA-B*53:06', 'HLA-B*53:07', 'HLA-B*53:08', 'HLA-B*53:09', 'HLA-B*53:10', 'HLA-B*53:11', 'HLA-B*53:12', 'HLA-B*53:13',
         'HLA-B*53:14', 'HLA-B*53:15', 'HLA-B*53:16', 'HLA-B*53:17', 'HLA-B*53:18', 'HLA-B*53:19', 'HLA-B*53:20', 'HLA-B*53:21', 'HLA-B*53:22',
         'HLA-B*53:23', 'HLA-B*54:01', 'HLA-B*54:02', 'HLA-B*54:03', 'HLA-B*54:04', 'HLA-B*54:06', 'HLA-B*54:07', 'HLA-B*54:09', 'HLA-B*54:10',
         'HLA-B*54:11', 'HLA-B*54:12', 'HLA-B*54:13', 'HLA-B*54:14', 'HLA-B*54:15', 'HLA-B*54:16', 'HLA-B*54:17', 'HLA-B*54:18', 'HLA-B*54:19',
         'HLA-B*54:20', 'HLA-B*54:21', 'HLA-B*54:22', 'HLA-B*54:23', 'HLA-B*55:01', 'HLA-B*55:02', 'HLA-B*55:03', 'HLA-B*55:04', 'HLA-B*55:05',
         'HLA-B*55:07', 'HLA-B*55:08', 'HLA-B*55:09', 'HLA-B*55:10', 'HLA-B*55:11', 'HLA-B*55:12', 'HLA-B*55:13', 'HLA-B*55:14', 'HLA-B*55:15',
         'HLA-B*55:16', 'HLA-B*55:17', 'HLA-B*55:18', 'HLA-B*55:19', 'HLA-B*55:20', 'HLA-B*55:21', 'HLA-B*55:22', 'HLA-B*55:23', 'HLA-B*55:24',
         'HLA-B*55:25', 'HLA-B*55:26', 'HLA-B*55:27', 'HLA-B*55:28', 'HLA-B*55:29', 'HLA-B*55:30', 'HLA-B*55:31', 'HLA-B*55:32', 'HLA-B*55:33',
         'HLA-B*55:34', 'HLA-B*55:35', 'HLA-B*55:36', 'HLA-B*55:37', 'HLA-B*55:38', 'HLA-B*55:39', 'HLA-B*55:40', 'HLA-B*55:41', 'HLA-B*55:42',
         'HLA-B*55:43', 'HLA-B*56:01', 'HLA-B*56:02', 'HLA-B*56:03', 'HLA-B*56:04', 'HLA-B*56:05', 'HLA-B*56:06', 'HLA-B*56:07', 'HLA-B*56:08',
         'HLA-B*56:09', 'HLA-B*56:10', 'HLA-B*56:11', 'HLA-B*56:12', 'HLA-B*56:13', 'HLA-B*56:14', 'HLA-B*56:15', 'HLA-B*56:16', 'HLA-B*56:17',
         'HLA-B*56:18', 'HLA-B*56:20', 'HLA-B*56:21', 'HLA-B*56:22', 'HLA-B*56:23', 'HLA-B*56:24', 'HLA-B*56:25', 'HLA-B*56:26', 'HLA-B*56:27',
         'HLA-B*56:29', 'HLA-B*57:01', 'HLA-B*57:02', 'HLA-B*57:03', 'HLA-B*57:04', 'HLA-B*57:05', 'HLA-B*57:06', 'HLA-B*57:07', 'HLA-B*57:08',
         'HLA-B*57:09', 'HLA-B*57:10', 'HLA-B*57:11', 'HLA-B*57:12', 'HLA-B*57:13', 'HLA-B*57:14', 'HLA-B*57:15', 'HLA-B*57:16', 'HLA-B*57:17',
         'HLA-B*57:18', 'HLA-B*57:19', 'HLA-B*57:20', 'HLA-B*57:21', 'HLA-B*57:22', 'HLA-B*57:23', 'HLA-B*57:24', 'HLA-B*57:25', 'HLA-B*57:26',
         'HLA-B*57:27', 'HLA-B*57:29', 'HLA-B*57:30', 'HLA-B*57:31', 'HLA-B*57:32', 'HLA-B*58:01', 'HLA-B*58:02', 'HLA-B*58:04', 'HLA-B*58:05',
         'HLA-B*58:06', 'HLA-B*58:07', 'HLA-B*58:08', 'HLA-B*58:09', 'HLA-B*58:11', 'HLA-B*58:12', 'HLA-B*58:13', 'HLA-B*58:14', 'HLA-B*58:15',
         'HLA-B*58:16', 'HLA-B*58:18', 'HLA-B*58:19', 'HLA-B*58:20', 'HLA-B*58:21', 'HLA-B*58:22', 'HLA-B*58:23', 'HLA-B*58:24', 'HLA-B*58:25',
         'HLA-B*58:26', 'HLA-B*58:27', 'HLA-B*58:28', 'HLA-B*58:29', 'HLA-B*58:30', 'HLA-B*59:01', 'HLA-B*59:02', 'HLA-B*59:03', 'HLA-B*59:04',
         'HLA-B*59:05', 'HLA-B*67:01', 'HLA-B*67:02', 'HLA-B*73:01', 'HLA-B*73:02', 'HLA-B*78:01', 'HLA-B*78:02', 'HLA-B*78:03', 'HLA-B*78:04',
         'HLA-B*78:05', 'HLA-B*78:06', 'HLA-B*78:07', 'HLA-B*81:01', 'HLA-B*81:02', 'HLA-B*81:03', 'HLA-B*81:05', 'HLA-B*82:01', 'HLA-B*82:02',
         'HLA-B*82:03', 'HLA-B*83:01', 'HLA-C*01:02', 'HLA-C*01:03', 'HLA-C*01:04', 'HLA-C*01:05', 'HLA-C*01:06', 'HLA-C*01:07', 'HLA-C*01:08',
         'HLA-C*01:09', 'HLA-C*01:10', 'HLA-C*01:11', 'HLA-C*01:12', 'HLA-C*01:13', 'HLA-C*01:14', 'HLA-C*01:15', 'HLA-C*01:16', 'HLA-C*01:17',
         'HLA-C*01:18', 'HLA-C*01:19', 'HLA-C*01:20', 'HLA-C*01:21', 'HLA-C*01:22', 'HLA-C*01:23', 'HLA-C*01:24', 'HLA-C*01:25', 'HLA-C*01:26',
         'HLA-C*01:27', 'HLA-C*01:28', 'HLA-C*01:29', 'HLA-C*01:30', 'HLA-C*01:31', 'HLA-C*01:32', 'HLA-C*01:33', 'HLA-C*01:34', 'HLA-C*01:35',
         'HLA-C*01:36', 'HLA-C*01:38', 'HLA-C*01:39', 'HLA-C*01:40', 'HLA-C*02:02', 'HLA-C*02:03', 'HLA-C*02:04', 'HLA-C*02:05', 'HLA-C*02:06',
         'HLA-C*02:07', 'HLA-C*02:08', 'HLA-C*02:09', 'HLA-C*02:10', 'HLA-C*02:11', 'HLA-C*02:12', 'HLA-C*02:13', 'HLA-C*02:14', 'HLA-C*02:15',
         'HLA-C*02:16', 'HLA-C*02:17', 'HLA-C*02:18', 'HLA-C*02:19', 'HLA-C*02:20', 'HLA-C*02:21', 'HLA-C*02:22', 'HLA-C*02:23', 'HLA-C*02:24',
         'HLA-C*02:26', 'HLA-C*02:27', 'HLA-C*02:28', 'HLA-C*02:29', 'HLA-C*02:30', 'HLA-C*02:31', 'HLA-C*02:32', 'HLA-C*02:33', 'HLA-C*02:34',
         'HLA-C*02:35', 'HLA-C*02:36', 'HLA-C*02:37', 'HLA-C*02:39', 'HLA-C*02:40', 'HLA-C*03:01', 'HLA-C*03:02', 'HLA-C*03:03', 'HLA-C*03:04',
         'HLA-C*03:05', 'HLA-C*03:06', 'HLA-C*03:07', 'HLA-C*03:08', 'HLA-C*03:09', 'HLA-C*03:10', 'HLA-C*03:11', 'HLA-C*03:12', 'HLA-C*03:13',
         'HLA-C*03:14', 'HLA-C*03:15', 'HLA-C*03:16', 'HLA-C*03:17', 'HLA-C*03:18', 'HLA-C*03:19', 'HLA-C*03:21', 'HLA-C*03:23', 'HLA-C*03:24',
         'HLA-C*03:25', 'HLA-C*03:26', 'HLA-C*03:27', 'HLA-C*03:28', 'HLA-C*03:29', 'HLA-C*03:30', 'HLA-C*03:31', 'HLA-C*03:32', 'HLA-C*03:33',
         'HLA-C*03:34', 'HLA-C*03:35', 'HLA-C*03:36', 'HLA-C*03:37', 'HLA-C*03:38', 'HLA-C*03:39', 'HLA-C*03:40', 'HLA-C*03:41', 'HLA-C*03:42',
         'HLA-C*03:43', 'HLA-C*03:44', 'HLA-C*03:45', 'HLA-C*03:46', 'HLA-C*03:47', 'HLA-C*03:48', 'HLA-C*03:49', 'HLA-C*03:50', 'HLA-C*03:51',
         'HLA-C*03:52', 'HLA-C*03:53', 'HLA-C*03:54', 'HLA-C*03:55', 'HLA-C*03:56', 'HLA-C*03:57', 'HLA-C*03:58', 'HLA-C*03:59', 'HLA-C*03:60',
         'HLA-C*03:61', 'HLA-C*03:62', 'HLA-C*03:63', 'HLA-C*03:64', 'HLA-C*03:65', 'HLA-C*03:66', 'HLA-C*03:67', 'HLA-C*03:68', 'HLA-C*03:69',
         'HLA-C*03:70', 'HLA-C*03:71', 'HLA-C*03:72', 'HLA-C*03:73', 'HLA-C*03:74', 'HLA-C*03:75', 'HLA-C*03:76', 'HLA-C*03:77', 'HLA-C*03:78',
         'HLA-C*03:79', 'HLA-C*03:80', 'HLA-C*03:81', 'HLA-C*03:82', 'HLA-C*03:83', 'HLA-C*03:84', 'HLA-C*03:85', 'HLA-C*03:86', 'HLA-C*03:87',
         'HLA-C*03:88', 'HLA-C*03:89', 'HLA-C*03:90', 'HLA-C*03:91', 'HLA-C*03:92', 'HLA-C*03:93', 'HLA-C*03:94', 'HLA-C*04:01', 'HLA-C*04:03',
         'HLA-C*04:04', 'HLA-C*04:05', 'HLA-C*04:06', 'HLA-C*04:07', 'HLA-C*04:08', 'HLA-C*04:10', 'HLA-C*04:11', 'HLA-C*04:12', 'HLA-C*04:13',
         'HLA-C*04:14', 'HLA-C*04:15', 'HLA-C*04:16', 'HLA-C*04:17', 'HLA-C*04:18', 'HLA-C*04:19', 'HLA-C*04:20', 'HLA-C*04:23', 'HLA-C*04:24',
         'HLA-C*04:25', 'HLA-C*04:26', 'HLA-C*04:27', 'HLA-C*04:28', 'HLA-C*04:29', 'HLA-C*04:30', 'HLA-C*04:31', 'HLA-C*04:32', 'HLA-C*04:33',
         'HLA-C*04:34', 'HLA-C*04:35', 'HLA-C*04:36', 'HLA-C*04:37', 'HLA-C*04:38', 'HLA-C*04:39', 'HLA-C*04:40', 'HLA-C*04:41', 'HLA-C*04:42',
         'HLA-C*04:43', 'HLA-C*04:44', 'HLA-C*04:45', 'HLA-C*04:46', 'HLA-C*04:47', 'HLA-C*04:48', 'HLA-C*04:49', 'HLA-C*04:50', 'HLA-C*04:51',
         'HLA-C*04:52', 'HLA-C*04:53', 'HLA-C*04:54', 'HLA-C*04:55', 'HLA-C*04:56', 'HLA-C*04:57', 'HLA-C*04:58', 'HLA-C*04:60', 'HLA-C*04:61',
         'HLA-C*04:62', 'HLA-C*04:63', 'HLA-C*04:64', 'HLA-C*04:65', 'HLA-C*04:66', 'HLA-C*04:67', 'HLA-C*04:68', 'HLA-C*04:69', 'HLA-C*04:70',
         'HLA-C*05:01', 'HLA-C*05:03', 'HLA-C*05:04', 'HLA-C*05:05', 'HLA-C*05:06', 'HLA-C*05:08', 'HLA-C*05:09', 'HLA-C*05:10', 'HLA-C*05:11',
         'HLA-C*05:12', 'HLA-C*05:13', 'HLA-C*05:14', 'HLA-C*05:15', 'HLA-C*05:16', 'HLA-C*05:17', 'HLA-C*05:18', 'HLA-C*05:19', 'HLA-C*05:20',
         'HLA-C*05:21', 'HLA-C*05:22', 'HLA-C*05:23', 'HLA-C*05:24', 'HLA-C*05:25', 'HLA-C*05:26', 'HLA-C*05:27', 'HLA-C*05:28', 'HLA-C*05:29',
         'HLA-C*05:30', 'HLA-C*05:31', 'HLA-C*05:32', 'HLA-C*05:33', 'HLA-C*05:34', 'HLA-C*05:35', 'HLA-C*05:36', 'HLA-C*05:37', 'HLA-C*05:38',
         'HLA-C*05:39', 'HLA-C*05:40', 'HLA-C*05:41', 'HLA-C*05:42', 'HLA-C*05:43', 'HLA-C*05:44', 'HLA-C*05:45', 'HLA-C*06:02', 'HLA-C*06:03',
         'HLA-C*06:04', 'HLA-C*06:05', 'HLA-C*06:06', 'HLA-C*06:07', 'HLA-C*06:08', 'HLA-C*06:09', 'HLA-C*06:10', 'HLA-C*06:11', 'HLA-C*06:12',
         'HLA-C*06:13', 'HLA-C*06:14', 'HLA-C*06:15', 'HLA-C*06:17', 'HLA-C*06:18', 'HLA-C*06:19', 'HLA-C*06:20', 'HLA-C*06:21', 'HLA-C*06:22',
         'HLA-C*06:23', 'HLA-C*06:24', 'HLA-C*06:25', 'HLA-C*06:26', 'HLA-C*06:27', 'HLA-C*06:28', 'HLA-C*06:29', 'HLA-C*06:30', 'HLA-C*06:31',
         'HLA-C*06:32', 'HLA-C*06:33', 'HLA-C*06:34', 'HLA-C*06:35', 'HLA-C*06:36', 'HLA-C*06:37', 'HLA-C*06:38', 'HLA-C*06:39', 'HLA-C*06:40',
         'HLA-C*06:41', 'HLA-C*06:42', 'HLA-C*06:43', 'HLA-C*06:44', 'HLA-C*06:45', 'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*07:03', 'HLA-C*07:04',
         'HLA-C*07:05', 'HLA-C*07:06', 'HLA-C*07:07', 'HLA-C*07:08', 'HLA-C*07:09', 'HLA-C*07:10', 'HLA-C*07:100', 'HLA-C*07:101', 'HLA-C*07:102',
         'HLA-C*07:103', 'HLA-C*07:105', 'HLA-C*07:106', 'HLA-C*07:107', 'HLA-C*07:108', 'HLA-C*07:109', 'HLA-C*07:11', 'HLA-C*07:110',
         'HLA-C*07:111', 'HLA-C*07:112', 'HLA-C*07:113', 'HLA-C*07:114', 'HLA-C*07:115', 'HLA-C*07:116', 'HLA-C*07:117', 'HLA-C*07:118',
         'HLA-C*07:119', 'HLA-C*07:12', 'HLA-C*07:120', 'HLA-C*07:122', 'HLA-C*07:123', 'HLA-C*07:124', 'HLA-C*07:125', 'HLA-C*07:126',
         'HLA-C*07:127', 'HLA-C*07:128', 'HLA-C*07:129', 'HLA-C*07:13', 'HLA-C*07:130', 'HLA-C*07:131', 'HLA-C*07:132', 'HLA-C*07:133',
         'HLA-C*07:134', 'HLA-C*07:135', 'HLA-C*07:136', 'HLA-C*07:137', 'HLA-C*07:138', 'HLA-C*07:139', 'HLA-C*07:14', 'HLA-C*07:140',
         'HLA-C*07:141', 'HLA-C*07:142', 'HLA-C*07:143', 'HLA-C*07:144', 'HLA-C*07:145', 'HLA-C*07:146', 'HLA-C*07:147', 'HLA-C*07:148',
         'HLA-C*07:149', 'HLA-C*07:15', 'HLA-C*07:16', 'HLA-C*07:17', 'HLA-C*07:18', 'HLA-C*07:19', 'HLA-C*07:20', 'HLA-C*07:21', 'HLA-C*07:22',
         'HLA-C*07:23', 'HLA-C*07:24', 'HLA-C*07:25', 'HLA-C*07:26', 'HLA-C*07:27', 'HLA-C*07:28', 'HLA-C*07:29', 'HLA-C*07:30', 'HLA-C*07:31',
         'HLA-C*07:35', 'HLA-C*07:36', 'HLA-C*07:37', 'HLA-C*07:38', 'HLA-C*07:39', 'HLA-C*07:40', 'HLA-C*07:41', 'HLA-C*07:42', 'HLA-C*07:43',
         'HLA-C*07:44', 'HLA-C*07:45', 'HLA-C*07:46', 'HLA-C*07:47', 'HLA-C*07:48', 'HLA-C*07:49', 'HLA-C*07:50', 'HLA-C*07:51', 'HLA-C*07:52',
         'HLA-C*07:53', 'HLA-C*07:54', 'HLA-C*07:56', 'HLA-C*07:57', 'HLA-C*07:58', 'HLA-C*07:59', 'HLA-C*07:60', 'HLA-C*07:62', 'HLA-C*07:63',
         'HLA-C*07:64', 'HLA-C*07:65', 'HLA-C*07:66', 'HLA-C*07:67', 'HLA-C*07:68', 'HLA-C*07:69', 'HLA-C*07:70', 'HLA-C*07:71', 'HLA-C*07:72',
         'HLA-C*07:73', 'HLA-C*07:74', 'HLA-C*07:75', 'HLA-C*07:76', 'HLA-C*07:77', 'HLA-C*07:78', 'HLA-C*07:79', 'HLA-C*07:80', 'HLA-C*07:81',
         'HLA-C*07:82', 'HLA-C*07:83', 'HLA-C*07:84', 'HLA-C*07:85', 'HLA-C*07:86', 'HLA-C*07:87', 'HLA-C*07:88', 'HLA-C*07:89', 'HLA-C*07:90',
         'HLA-C*07:91', 'HLA-C*07:92', 'HLA-C*07:93', 'HLA-C*07:94', 'HLA-C*07:95', 'HLA-C*07:96', 'HLA-C*07:97', 'HLA-C*07:99', 'HLA-C*08:01',
         'HLA-C*08:02', 'HLA-C*08:03', 'HLA-C*08:04', 'HLA-C*08:05', 'HLA-C*08:06', 'HLA-C*08:07', 'HLA-C*08:08', 'HLA-C*08:09', 'HLA-C*08:10',
         'HLA-C*08:11', 'HLA-C*08:12', 'HLA-C*08:13', 'HLA-C*08:14', 'HLA-C*08:15', 'HLA-C*08:16', 'HLA-C*08:17', 'HLA-C*08:18', 'HLA-C*08:19',
         'HLA-C*08:20', 'HLA-C*08:21', 'HLA-C*08:22', 'HLA-C*08:23', 'HLA-C*08:24', 'HLA-C*08:25', 'HLA-C*08:27', 'HLA-C*08:28', 'HLA-C*08:29',
         'HLA-C*08:30', 'HLA-C*08:31', 'HLA-C*08:32', 'HLA-C*08:33', 'HLA-C*08:34', 'HLA-C*08:35', 'HLA-C*12:02', 'HLA-C*12:03', 'HLA-C*12:04',
         'HLA-C*12:05', 'HLA-C*12:06', 'HLA-C*12:07', 'HLA-C*12:08', 'HLA-C*12:09', 'HLA-C*12:10', 'HLA-C*12:11', 'HLA-C*12:12', 'HLA-C*12:13',
         'HLA-C*12:14', 'HLA-C*12:15', 'HLA-C*12:16', 'HLA-C*12:17', 'HLA-C*12:18', 'HLA-C*12:19', 'HLA-C*12:20', 'HLA-C*12:21', 'HLA-C*12:22',
         'HLA-C*12:23', 'HLA-C*12:24', 'HLA-C*12:25', 'HLA-C*12:26', 'HLA-C*12:27', 'HLA-C*12:28', 'HLA-C*12:29', 'HLA-C*12:30', 'HLA-C*12:31',
         'HLA-C*12:32', 'HLA-C*12:33', 'HLA-C*12:34', 'HLA-C*12:35', 'HLA-C*12:36', 'HLA-C*12:37', 'HLA-C*12:38', 'HLA-C*12:40', 'HLA-C*12:41',
         'HLA-C*12:43', 'HLA-C*12:44', 'HLA-C*14:02', 'HLA-C*14:03', 'HLA-C*14:04', 'HLA-C*14:05', 'HLA-C*14:06', 'HLA-C*14:08', 'HLA-C*14:09',
         'HLA-C*14:10', 'HLA-C*14:11', 'HLA-C*14:12', 'HLA-C*14:13', 'HLA-C*14:14', 'HLA-C*14:15', 'HLA-C*14:16', 'HLA-C*14:17', 'HLA-C*14:18',
         'HLA-C*14:19', 'HLA-C*14:20', 'HLA-C*15:02', 'HLA-C*15:03', 'HLA-C*15:04', 'HLA-C*15:05', 'HLA-C*15:06', 'HLA-C*15:07', 'HLA-C*15:08',
         'HLA-C*15:09', 'HLA-C*15:10', 'HLA-C*15:11', 'HLA-C*15:12', 'HLA-C*15:13', 'HLA-C*15:15', 'HLA-C*15:16', 'HLA-C*15:17', 'HLA-C*15:18',
         'HLA-C*15:19', 'HLA-C*15:20', 'HLA-C*15:21', 'HLA-C*15:22', 'HLA-C*15:23', 'HLA-C*15:24', 'HLA-C*15:25', 'HLA-C*15:26', 'HLA-C*15:27',
         'HLA-C*15:28', 'HLA-C*15:29', 'HLA-C*15:30', 'HLA-C*15:31', 'HLA-C*15:33', 'HLA-C*15:34', 'HLA-C*15:35', 'HLA-C*16:01', 'HLA-C*16:02',
         'HLA-C*16:04', 'HLA-C*16:06', 'HLA-C*16:07', 'HLA-C*16:08', 'HLA-C*16:09', 'HLA-C*16:10', 'HLA-C*16:11', 'HLA-C*16:12', 'HLA-C*16:13',
         'HLA-C*16:14', 'HLA-C*16:15', 'HLA-C*16:17', 'HLA-C*16:18', 'HLA-C*16:19', 'HLA-C*16:20', 'HLA-C*16:21', 'HLA-C*16:22', 'HLA-C*16:23',
         'HLA-C*16:24', 'HLA-C*16:25', 'HLA-C*16:26', 'HLA-C*17:01', 'HLA-C*17:02', 'HLA-C*17:03', 'HLA-C*17:04', 'HLA-C*17:05', 'HLA-C*17:06',
         'HLA-C*17:07', 'HLA-C*18:01', 'HLA-C*18:02', 'HLA-C*18:03', 'HLA-E*01:01', 'HLA-G*01:01', 'HLA-G*01:02', 'HLA-G*01:03', 'HLA-G*01:04',
         'HLA-G*01:06', 'HLA-G*01:07', 'HLA-G*01:08', 'HLA-G*01:09',
         'H2-Db', 'H2-Dd', 'H2-Kb', 'H2-Kd', 'H2-Kk', 'H2-Ld', "H2-Qa1", "H2-Qa2"])

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    def _represent(self, allele):
        """
        Internal function transforming an allele object into its representative string
        :param allele: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: str
        """
        if isinstance(allele, MouseAllele):
            return "H-2-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)
        else:
            return "HLA-%s%s:%s" % (allele.locus, allele.supertype, allele.subtype)

    def convert_alleles(self, alleles):
        """
        Converts :class:`~Fred2.Core.Allele.Allele` into the internal :class:`~Fred2.Core.Allele.Allele` representation
        of the predictor and returns a string representation

        :param alleles: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`
        :rtype: list(str)
        """
        return [self._represent(a) for a in alleles]

    @property
    def supportedAlleles(self):
        """A list of valid :class:`~Fred2.Core.Allele.Allele` models"""
        return self.__alleles

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    @property
    def command(self):
        """
        Defines the commandline call for external tool
        """
        return self.__command

    @property
    def supportedLength(self):
        """
        A list of supported :class:`~Fred2.Core.Peptide.Peptide` lengths
        """
        return self.__supported_length

    def get_external_version(self, path=None):
        """
        Returns the external version of the tool by executing
        >{command} --version

        might be dependent on the method and has to be overwritten
        therefore it is declared abstract to enforce the user to
        overwrite the method. The function in the base class can be called
        with super()

        :param str path: Optional specification of executable path if deviant from :attr:`self.__command`
        :return: The external version of the tool or None if tool does not support versioning
        :rtype: str
        """
        return None

    def prepare_input(self, input, file):
        """
        Prepares input for external tools
        and writes them to file in the specific format

        No return value!

        :param: list(str) input: The :class:`~Fred2.Core.Peptide.Peptide` sequences to write into file
        :param File file: File-handler to input file for external tool
        """
        file.write("\n".join(input))

    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        result = defaultdict(defaultdict)
        f = csv.reader(open(file, "r"), delimiter='\t')
        alleles = list([x for x in next(f) if x != ""])
        next(f)
        ic_pos = 3
        for row in f:
            pep_seq = row[1]
            for i, a in enumerate(alleles):
                result[a][pep_seq] = float(row[ic_pos + i * 3])
        return result


class NetMHCpan_3_0(NetMHCpan_2_8):
    """
        Implements the NetMHC binding version 3.0
        Supported  MHC alleles currently only restricted to HLA alleles.

    .. note::

        Nielsen, M., & Andreatta, M. (2016).
        NetMHCpan-3.0; improved prediction of binding to MHC class I molecules integrating information from multiple
        receptor and peptide length datasets. Genome Medicine, 8(1), 1.
    """

    __version = "3.0"
    __command = "netMHCpan -p {peptides} -a {alleles} {options} -xls -xlsfile {out}"


    @property
    def version(self):
        return self.__version

    @property
    def command(self):
        return self.__command

    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        result = defaultdict(defaultdict)
        f = csv.reader(open(file, "r"), delimiter='\t')
        alleles = list([x for x in next(f) if x != ""])
        next(f)
        ic_pos = 4
        for row in f:
            pep_seq = row[1]
            for i, a in enumerate(alleles):
                result[a][pep_seq] = float(row[ic_pos + i * 4])
        return result


class NetMHCstabpan_1_0(AExternalEpitopePrediction):
    """
    Implements a wrapper to NetMHCstabpan 1.0

    .. note:

    Pan-specific prediction of peptide-MHC-I complex stability; a correlate of T cell immunogenicity
    M Rasmussen, E Fenoy, M Nielsen, Buus S, Accepted JI June, 2016
    """
    __name = "netMHCstabpan"
    __length = frozenset([8, 9, 10, 11])
    __version = "1.0"
    __command = "netMHCstabpan -p {peptides} -a {alleles} {options} -xls -xlsfile {out}"
    __alleles = frozenset(['HLA-A*01:01', 'HLA-A*01:02', 'HLA-A*01:03', 'HLA-A*01:06', 'HLA-A*01:07', 'HLA-A*01:08', 'HLA-A*01:09', 'HLA-A*01:10', 'HLA-A*01:12',
                            'HLA-A*01:13', 'HLA-A*01:14', 'HLA-A*01:17', 'HLA-A*01:19', 'HLA-A*01:20', 'HLA-A*01:21', 'HLA-A*01:23', 'HLA-A*01:24',
                            'HLA-A*01:25', 'HLA-A*01:26', 'HLA-A*01:28', 'HLA-A*01:29', 'HLA-A*01:30', 'HLA-A*01:32', 'HLA-A*01:33', 'HLA-A*01:35',
                            'HLA-A*01:36', 'HLA-A*01:37', 'HLA-A*01:38', 'HLA-A*01:39', 'HLA-A*01:40', 'HLA-A*01:41', 'HLA-A*01:42', 'HLA-A*01:43',
                            'HLA-A*01:44', 'HLA-A*01:45', 'HLA-A*01:46', 'HLA-A*01:47', 'HLA-A*01:48', 'HLA-A*01:49', 'HLA-A*01:50', 'HLA-A*01:51',
                            'HLA-A*01:54', 'HLA-A*01:55', 'HLA-A*01:58', 'HLA-A*01:59', 'HLA-A*01:60', 'HLA-A*01:61', 'HLA-A*01:62', 'HLA-A*01:63',
                            'HLA-A*01:64', 'HLA-A*01:65', 'HLA-A*01:66', 'HLA-A*02:01', 'HLA-A*02:02', 'HLA-A*02:03', 'HLA-A*02:04', 'HLA-A*02:05',
                            'HLA-A*02:06', 'HLA-A*02:07', 'HLA-A*02:08', 'HLA-A*02:09', 'HLA-A*02:10', 'HLA-A*02:11', 'HLA-A*02:12', 'HLA-A*02:13',
                            'HLA-A*02:14', 'HLA-A*02:16', 'HLA-A*02:17', 'HLA-A*02:18', 'HLA-A*02:19', 'HLA-A*02:20', 'HLA-A*02:21', 'HLA-A*02:22',
                            'HLA-A*02:24', 'HLA-A*02:25', 'HLA-A*02:26', 'HLA-A*02:27', 'HLA-A*02:28', 'HLA-A*02:29', 'HLA-A*02:30', 'HLA-A*02:31',
                            'HLA-A*02:33', 'HLA-A*02:34', 'HLA-A*02:35', 'HLA-A*02:36', 'HLA-A*02:37', 'HLA-A*02:38', 'HLA-A*02:39', 'HLA-A*02:40',
                            'HLA-A*02:41', 'HLA-A*02:42', 'HLA-A*02:44', 'HLA-A*02:45', 'HLA-A*02:46', 'HLA-A*02:47', 'HLA-A*02:48', 'HLA-A*02:49',
                            'HLA-A*02:50', 'HLA-A*02:51', 'HLA-A*02:52', 'HLA-A*02:54', 'HLA-A*02:55', 'HLA-A*02:56', 'HLA-A*02:57', 'HLA-A*02:58',
                            'HLA-A*02:59', 'HLA-A*02:60', 'HLA-A*02:61', 'HLA-A*02:62', 'HLA-A*02:63', 'HLA-A*02:64', 'HLA-A*02:65', 'HLA-A*02:66',
                            'HLA-A*02:67', 'HLA-A*02:68', 'HLA-A*02:69', 'HLA-A*02:70', 'HLA-A*02:71', 'HLA-A*02:72', 'HLA-A*02:73', 'HLA-A*02:74',
                            'HLA-A*02:75', 'HLA-A*02:76', 'HLA-A*02:77', 'HLA-A*02:78', 'HLA-A*02:79', 'HLA-A*02:80', 'HLA-A*02:81', 'HLA-A*02:84',
                            'HLA-A*02:85', 'HLA-A*02:86', 'HLA-A*02:87', 'HLA-A*02:89', 'HLA-A*02:90', 'HLA-A*02:91', 'HLA-A*02:92', 'HLA-A*02:93',
                            'HLA-A*02:95', 'HLA-A*02:96', 'HLA-A*02:97', 'HLA-A*02:99', 'HLA-A*02:101', 'HLA-A*02:102', 'HLA-A*02:103', 'HLA-A*02:104',
                            'HLA-A*02:105', 'HLA-A*02:106', 'HLA-A*02:107', 'HLA-A*02:108', 'HLA-A*02:109', 'HLA-A*02:110', 'HLA-A*02:111',
                            'HLA-A*02:112', 'HLA-A*02:114', 'HLA-A*02:115', 'HLA-A*02:116', 'HLA-A*02:117', 'HLA-A*02:118', 'HLA-A*02:119',
                            'HLA-A*02:120', 'HLA-A*02:121', 'HLA-A*02:122', 'HLA-A*02:123', 'HLA-A*02:124', 'HLA-A*02:126', 'HLA-A*02:127',
                            'HLA-A*02:128', 'HLA-A*02:129', 'HLA-A*02:130', 'HLA-A*02:131', 'HLA-A*02:132', 'HLA-A*02:133', 'HLA-A*02:134',
                            'HLA-A*02:135', 'HLA-A*02:136', 'HLA-A*02:137', 'HLA-A*02:138', 'HLA-A*02:139', 'HLA-A*02:140', 'HLA-A*02:141',
                            'HLA-A*02:142', 'HLA-A*02:143', 'HLA-A*02:144', 'HLA-A*02:145', 'HLA-A*02:146', 'HLA-A*02:147', 'HLA-A*02:148',
                            'HLA-A*02:149', 'HLA-A*02:150', 'HLA-A*02:151', 'HLA-A*02:152', 'HLA-A*02:153', 'HLA-A*02:154', 'HLA-A*02:155',
                            'HLA-A*02:156', 'HLA-A*02:157', 'HLA-A*02:158', 'HLA-A*02:159', 'HLA-A*02:160', 'HLA-A*02:161', 'HLA-A*02:162',
                            'HLA-A*02:163', 'HLA-A*02:164', 'HLA-A*02:165', 'HLA-A*02:166', 'HLA-A*02:167', 'HLA-A*02:168', 'HLA-A*02:169',
                            'HLA-A*02:170', 'HLA-A*02:171', 'HLA-A*02:172', 'HLA-A*02:173', 'HLA-A*02:174', 'HLA-A*02:175', 'HLA-A*02:176',
                            'HLA-A*02:177', 'HLA-A*02:178', 'HLA-A*02:179', 'HLA-A*02:180', 'HLA-A*02:181', 'HLA-A*02:182', 'HLA-A*02:183',
                            'HLA-A*02:184', 'HLA-A*02:185', 'HLA-A*02:186', 'HLA-A*02:187', 'HLA-A*02:188', 'HLA-A*02:189', 'HLA-A*02:190',
                            'HLA-A*02:191', 'HLA-A*02:192', 'HLA-A*02:193', 'HLA-A*02:194', 'HLA-A*02:195', 'HLA-A*02:196', 'HLA-A*02:197',
                            'HLA-A*02:198', 'HLA-A*02:199', 'HLA-A*02:200', 'HLA-A*02:201', 'HLA-A*02:202', 'HLA-A*02:203', 'HLA-A*02:204',
                            'HLA-A*02:205', 'HLA-A*02:206', 'HLA-A*02:207', 'HLA-A*02:208', 'HLA-A*02:209', 'HLA-A*02:210', 'HLA-A*02:211',
                            'HLA-A*02:212', 'HLA-A*02:213', 'HLA-A*02:214', 'HLA-A*02:215', 'HLA-A*02:216', 'HLA-A*02:217', 'HLA-A*02:218',
                            'HLA-A*02:219', 'HLA-A*02:220', 'HLA-A*02:221', 'HLA-A*02:224', 'HLA-A*02:228', 'HLA-A*02:229', 'HLA-A*02:230',
                            'HLA-A*02:231', 'HLA-A*02:232', 'HLA-A*02:233', 'HLA-A*02:234', 'HLA-A*02:235', 'HLA-A*02:236', 'HLA-A*02:237',
                            'HLA-A*02:238', 'HLA-A*02:239', 'HLA-A*02:240', 'HLA-A*02:241', 'HLA-A*02:242', 'HLA-A*02:243', 'HLA-A*02:244',
                            'HLA-A*02:245', 'HLA-A*02:246', 'HLA-A*02:247', 'HLA-A*02:248', 'HLA-A*02:249', 'HLA-A*02:251', 'HLA-A*02:252',
                            'HLA-A*02:253', 'HLA-A*02:254', 'HLA-A*02:255', 'HLA-A*02:256', 'HLA-A*02:257', 'HLA-A*02:258', 'HLA-A*02:259',
                            'HLA-A*02:260', 'HLA-A*02:261', 'HLA-A*02:262', 'HLA-A*02:263', 'HLA-A*02:264', 'HLA-A*02:265', 'HLA-A*02:266',
                            'HLA-A*03:01', 'HLA-A*03:02', 'HLA-A*03:04', 'HLA-A*03:05', 'HLA-A*03:06', 'HLA-A*03:07', 'HLA-A*03:08', 'HLA-A*03:09',
                            'HLA-A*03:10', 'HLA-A*03:12', 'HLA-A*03:13', 'HLA-A*03:14', 'HLA-A*03:15', 'HLA-A*03:16', 'HLA-A*03:17', 'HLA-A*03:18',
                            'HLA-A*03:19', 'HLA-A*03:20', 'HLA-A*03:22', 'HLA-A*03:23', 'HLA-A*03:24', 'HLA-A*03:25', 'HLA-A*03:26', 'HLA-A*03:27',
                            'HLA-A*03:28', 'HLA-A*03:29', 'HLA-A*03:30', 'HLA-A*03:31', 'HLA-A*03:32', 'HLA-A*03:33', 'HLA-A*03:34', 'HLA-A*03:35',
                            'HLA-A*03:37', 'HLA-A*03:38', 'HLA-A*03:39', 'HLA-A*03:40', 'HLA-A*03:41', 'HLA-A*03:42', 'HLA-A*03:43', 'HLA-A*03:44',
                            'HLA-A*03:45', 'HLA-A*03:46', 'HLA-A*03:47', 'HLA-A*03:48', 'HLA-A*03:49', 'HLA-A*03:50', 'HLA-A*03:51', 'HLA-A*03:52',
                            'HLA-A*03:53', 'HLA-A*03:54', 'HLA-A*03:55', 'HLA-A*03:56', 'HLA-A*03:57', 'HLA-A*03:58', 'HLA-A*03:59', 'HLA-A*03:60',
                            'HLA-A*03:61', 'HLA-A*03:62', 'HLA-A*03:63', 'HLA-A*03:64', 'HLA-A*03:65', 'HLA-A*03:66', 'HLA-A*03:67', 'HLA-A*03:70',
                            'HLA-A*03:71', 'HLA-A*03:72', 'HLA-A*03:73', 'HLA-A*03:74', 'HLA-A*03:75', 'HLA-A*03:76', 'HLA-A*03:77', 'HLA-A*03:78',
                            'HLA-A*03:79', 'HLA-A*03:80', 'HLA-A*03:81', 'HLA-A*03:82', 'HLA-A*11:01', 'HLA-A*11:02', 'HLA-A*11:03', 'HLA-A*11:04',
                            'HLA-A*11:05', 'HLA-A*11:06', 'HLA-A*11:07', 'HLA-A*11:08', 'HLA-A*11:09', 'HLA-A*11:10', 'HLA-A*11:11', 'HLA-A*11:12',
                            'HLA-A*11:13', 'HLA-A*11:14', 'HLA-A*11:15', 'HLA-A*11:16', 'HLA-A*11:17', 'HLA-A*11:18', 'HLA-A*11:19', 'HLA-A*11:20',
                            'HLA-A*11:22', 'HLA-A*11:23', 'HLA-A*11:24', 'HLA-A*11:25', 'HLA-A*11:26', 'HLA-A*11:27', 'HLA-A*11:29', 'HLA-A*11:30',
                            'HLA-A*11:31', 'HLA-A*11:32', 'HLA-A*11:33', 'HLA-A*11:34', 'HLA-A*11:35', 'HLA-A*11:36', 'HLA-A*11:37', 'HLA-A*11:38',
                            'HLA-A*11:39', 'HLA-A*11:40', 'HLA-A*11:41', 'HLA-A*11:42', 'HLA-A*11:43', 'HLA-A*11:44', 'HLA-A*11:45', 'HLA-A*11:46',
                            'HLA-A*11:47', 'HLA-A*11:48', 'HLA-A*11:49', 'HLA-A*11:51', 'HLA-A*11:53', 'HLA-A*11:54', 'HLA-A*11:55', 'HLA-A*11:56',
                            'HLA-A*11:57', 'HLA-A*11:58', 'HLA-A*11:59', 'HLA-A*11:60', 'HLA-A*11:61', 'HLA-A*11:62', 'HLA-A*11:63', 'HLA-A*11:64',
                            'HLA-A*23:01', 'HLA-A*23:02', 'HLA-A*23:03', 'HLA-A*23:04', 'HLA-A*23:05', 'HLA-A*23:06', 'HLA-A*23:09', 'HLA-A*23:10',
                            'HLA-A*23:12', 'HLA-A*23:13', 'HLA-A*23:14', 'HLA-A*23:15', 'HLA-A*23:16', 'HLA-A*23:17', 'HLA-A*23:18', 'HLA-A*23:20',
                            'HLA-A*23:21', 'HLA-A*23:22', 'HLA-A*23:23', 'HLA-A*23:24', 'HLA-A*23:25', 'HLA-A*23:26', 'HLA-A*24:02', 'HLA-A*24:03',
                            'HLA-A*24:04', 'HLA-A*24:05', 'HLA-A*24:06', 'HLA-A*24:07', 'HLA-A*24:08', 'HLA-A*24:10', 'HLA-A*24:13', 'HLA-A*24:14',
                            'HLA-A*24:15', 'HLA-A*24:17', 'HLA-A*24:18', 'HLA-A*24:19', 'HLA-A*24:20', 'HLA-A*24:21', 'HLA-A*24:22', 'HLA-A*24:23',
                            'HLA-A*24:24', 'HLA-A*24:25', 'HLA-A*24:26', 'HLA-A*24:27', 'HLA-A*24:28', 'HLA-A*24:29', 'HLA-A*24:30', 'HLA-A*24:31',
                            'HLA-A*24:32', 'HLA-A*24:33', 'HLA-A*24:34', 'HLA-A*24:35', 'HLA-A*24:37', 'HLA-A*24:38', 'HLA-A*24:39', 'HLA-A*24:41',
                            'HLA-A*24:42', 'HLA-A*24:43', 'HLA-A*24:44', 'HLA-A*24:46', 'HLA-A*24:47', 'HLA-A*24:49', 'HLA-A*24:50', 'HLA-A*24:51',
                            'HLA-A*24:52', 'HLA-A*24:53', 'HLA-A*24:54', 'HLA-A*24:55', 'HLA-A*24:56', 'HLA-A*24:57', 'HLA-A*24:58', 'HLA-A*24:59',
                            'HLA-A*24:61', 'HLA-A*24:62', 'HLA-A*24:63', 'HLA-A*24:64', 'HLA-A*24:66', 'HLA-A*24:67', 'HLA-A*24:68', 'HLA-A*24:69',
                            'HLA-A*24:70', 'HLA-A*24:71', 'HLA-A*24:72', 'HLA-A*24:73', 'HLA-A*24:74', 'HLA-A*24:75', 'HLA-A*24:76', 'HLA-A*24:77',
                            'HLA-A*24:78', 'HLA-A*24:79', 'HLA-A*24:80', 'HLA-A*24:81', 'HLA-A*24:82', 'HLA-A*24:85', 'HLA-A*24:87', 'HLA-A*24:88',
                            'HLA-A*24:89', 'HLA-A*24:91', 'HLA-A*24:92', 'HLA-A*24:93', 'HLA-A*24:94', 'HLA-A*24:95', 'HLA-A*24:96', 'HLA-A*24:97',
                            'HLA-A*24:98', 'HLA-A*24:99', 'HLA-A*24:100', 'HLA-A*24:101', 'HLA-A*24:102', 'HLA-A*24:103', 'HLA-A*24:104',
                            'HLA-A*24:105', 'HLA-A*24:106', 'HLA-A*24:107', 'HLA-A*24:108', 'HLA-A*24:109', 'HLA-A*24:110', 'HLA-A*24:111',
                            'HLA-A*24:112', 'HLA-A*24:113', 'HLA-A*24:114', 'HLA-A*24:115', 'HLA-A*24:116', 'HLA-A*24:117', 'HLA-A*24:118',
                            'HLA-A*24:119', 'HLA-A*24:120', 'HLA-A*24:121', 'HLA-A*24:122', 'HLA-A*24:123', 'HLA-A*24:124', 'HLA-A*24:125',
                            'HLA-A*24:126', 'HLA-A*24:127', 'HLA-A*24:128', 'HLA-A*24:129', 'HLA-A*24:130', 'HLA-A*24:131', 'HLA-A*24:133',
                            'HLA-A*24:134', 'HLA-A*24:135', 'HLA-A*24:136', 'HLA-A*24:137', 'HLA-A*24:138', 'HLA-A*24:139', 'HLA-A*24:140',
                            'HLA-A*24:141', 'HLA-A*24:142', 'HLA-A*24:143', 'HLA-A*24:144', 'HLA-A*25:01', 'HLA-A*25:02', 'HLA-A*25:03', 'HLA-A*25:04',
                            'HLA-A*25:05', 'HLA-A*25:06', 'HLA-A*25:07', 'HLA-A*25:08', 'HLA-A*25:09', 'HLA-A*25:10', 'HLA-A*25:11', 'HLA-A*25:13',
                            'HLA-A*26:01', 'HLA-A*26:02', 'HLA-A*26:03', 'HLA-A*26:04', 'HLA-A*26:05', 'HLA-A*26:06', 'HLA-A*26:07', 'HLA-A*26:08',
                            'HLA-A*26:09', 'HLA-A*26:10', 'HLA-A*26:12', 'HLA-A*26:13', 'HLA-A*26:14', 'HLA-A*26:15', 'HLA-A*26:16', 'HLA-A*26:17',
                            'HLA-A*26:18', 'HLA-A*26:19', 'HLA-A*26:20', 'HLA-A*26:21', 'HLA-A*26:22', 'HLA-A*26:23', 'HLA-A*26:24', 'HLA-A*26:26',
                            'HLA-A*26:27', 'HLA-A*26:28', 'HLA-A*26:29', 'HLA-A*26:30', 'HLA-A*26:31', 'HLA-A*26:32', 'HLA-A*26:33', 'HLA-A*26:34',
                            'HLA-A*26:35', 'HLA-A*26:36', 'HLA-A*26:37', 'HLA-A*26:38', 'HLA-A*26:39', 'HLA-A*26:40', 'HLA-A*26:41', 'HLA-A*26:42',
                            'HLA-A*26:43', 'HLA-A*26:45', 'HLA-A*26:46', 'HLA-A*26:47', 'HLA-A*26:48', 'HLA-A*26:49', 'HLA-A*26:50', 'HLA-A*29:01',
                            'HLA-A*29:02', 'HLA-A*29:03', 'HLA-A*29:04', 'HLA-A*29:05', 'HLA-A*29:06', 'HLA-A*29:07', 'HLA-A*29:09', 'HLA-A*29:10',
                            'HLA-A*29:11', 'HLA-A*29:12', 'HLA-A*29:13', 'HLA-A*29:14', 'HLA-A*29:15', 'HLA-A*29:16', 'HLA-A*29:17', 'HLA-A*29:18',
                            'HLA-A*29:19', 'HLA-A*29:20', 'HLA-A*29:21', 'HLA-A*29:22', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*30:03', 'HLA-A*30:04',
                            'HLA-A*30:06', 'HLA-A*30:07', 'HLA-A*30:08', 'HLA-A*30:09', 'HLA-A*30:10', 'HLA-A*30:11', 'HLA-A*30:12', 'HLA-A*30:13',
                            'HLA-A*30:15', 'HLA-A*30:16', 'HLA-A*30:17', 'HLA-A*30:18', 'HLA-A*30:19', 'HLA-A*30:20', 'HLA-A*30:22', 'HLA-A*30:23',
                            'HLA-A*30:24', 'HLA-A*30:25', 'HLA-A*30:26', 'HLA-A*30:28', 'HLA-A*30:29', 'HLA-A*30:30', 'HLA-A*30:31', 'HLA-A*30:32',
                            'HLA-A*30:33', 'HLA-A*30:34', 'HLA-A*30:35', 'HLA-A*30:36', 'HLA-A*30:37', 'HLA-A*30:38', 'HLA-A*30:39', 'HLA-A*30:40',
                            'HLA-A*30:41', 'HLA-A*31:01', 'HLA-A*31:02', 'HLA-A*31:03', 'HLA-A*31:04', 'HLA-A*31:05', 'HLA-A*31:06', 'HLA-A*31:07',
                            'HLA-A*31:08', 'HLA-A*31:09', 'HLA-A*31:10', 'HLA-A*31:11', 'HLA-A*31:12', 'HLA-A*31:13', 'HLA-A*31:15', 'HLA-A*31:16',
                            'HLA-A*31:17', 'HLA-A*31:18', 'HLA-A*31:19', 'HLA-A*31:20', 'HLA-A*31:21', 'HLA-A*31:22', 'HLA-A*31:23', 'HLA-A*31:24',
                            'HLA-A*31:25', 'HLA-A*31:26', 'HLA-A*31:27', 'HLA-A*31:28', 'HLA-A*31:29', 'HLA-A*31:30', 'HLA-A*31:31', 'HLA-A*31:32',
                            'HLA-A*31:33', 'HLA-A*31:34', 'HLA-A*31:35', 'HLA-A*31:36', 'HLA-A*31:37', 'HLA-A*32:01', 'HLA-A*32:02', 'HLA-A*32:03',
                            'HLA-A*32:04', 'HLA-A*32:05', 'HLA-A*32:06', 'HLA-A*32:07', 'HLA-A*32:08', 'HLA-A*32:09', 'HLA-A*32:10', 'HLA-A*32:12',
                            'HLA-A*32:13', 'HLA-A*32:14', 'HLA-A*32:15', 'HLA-A*32:16', 'HLA-A*32:17', 'HLA-A*32:18', 'HLA-A*32:20', 'HLA-A*32:21',
                            'HLA-A*32:22', 'HLA-A*32:23', 'HLA-A*32:24', 'HLA-A*32:25', 'HLA-A*33:01', 'HLA-A*33:03', 'HLA-A*33:04', 'HLA-A*33:05',
                            'HLA-A*33:06', 'HLA-A*33:07', 'HLA-A*33:08', 'HLA-A*33:09', 'HLA-A*33:10', 'HLA-A*33:11', 'HLA-A*33:12', 'HLA-A*33:13',
                            'HLA-A*33:14', 'HLA-A*33:15', 'HLA-A*33:16', 'HLA-A*33:17', 'HLA-A*33:18', 'HLA-A*33:19', 'HLA-A*33:20', 'HLA-A*33:21',
                            'HLA-A*33:22', 'HLA-A*33:23', 'HLA-A*33:24', 'HLA-A*33:25', 'HLA-A*33:26', 'HLA-A*33:27', 'HLA-A*33:28', 'HLA-A*33:29',
                            'HLA-A*33:30', 'HLA-A*33:31', 'HLA-A*34:01', 'HLA-A*34:02', 'HLA-A*34:03', 'HLA-A*34:04', 'HLA-A*34:05', 'HLA-A*34:06',
                            'HLA-A*34:07', 'HLA-A*34:08', 'HLA-A*36:01', 'HLA-A*36:02', 'HLA-A*36:03', 'HLA-A*36:04', 'HLA-A*36:05', 'HLA-A*43:01',
                            'HLA-A*66:01', 'HLA-A*66:02', 'HLA-A*66:03', 'HLA-A*66:04', 'HLA-A*66:05', 'HLA-A*66:06', 'HLA-A*66:07', 'HLA-A*66:08',
                            'HLA-A*66:09', 'HLA-A*66:10', 'HLA-A*66:11', 'HLA-A*66:12', 'HLA-A*66:13', 'HLA-A*66:14', 'HLA-A*66:15', 'HLA-A*68:01',
                            'HLA-A*68:02', 'HLA-A*68:03', 'HLA-A*68:04', 'HLA-A*68:05', 'HLA-A*68:06', 'HLA-A*68:07', 'HLA-A*68:08', 'HLA-A*68:09',
                            'HLA-A*68:10', 'HLA-A*68:12', 'HLA-A*68:13', 'HLA-A*68:14', 'HLA-A*68:15', 'HLA-A*68:16', 'HLA-A*68:17', 'HLA-A*68:19',
                            'HLA-A*68:20', 'HLA-A*68:21', 'HLA-A*68:22', 'HLA-A*68:23', 'HLA-A*68:24', 'HLA-A*68:25', 'HLA-A*68:26', 'HLA-A*68:27',
                            'HLA-A*68:28', 'HLA-A*68:29', 'HLA-A*68:30', 'HLA-A*68:31', 'HLA-A*68:32', 'HLA-A*68:33', 'HLA-A*68:34', 'HLA-A*68:35',
                            'HLA-A*68:36', 'HLA-A*68:37', 'HLA-A*68:38', 'HLA-A*68:39', 'HLA-A*68:40', 'HLA-A*68:41', 'HLA-A*68:42', 'HLA-A*68:43',
                            'HLA-A*68:44', 'HLA-A*68:45', 'HLA-A*68:46', 'HLA-A*68:47', 'HLA-A*68:48', 'HLA-A*68:50', 'HLA-A*68:51', 'HLA-A*68:52',
                            'HLA-A*68:53', 'HLA-A*68:54', 'HLA-A*69:01', 'HLA-A*74:01', 'HLA-A*74:02', 'HLA-A*74:03', 'HLA-A*74:04', 'HLA-A*74:05',
                            'HLA-A*74:06', 'HLA-A*74:07', 'HLA-A*74:08', 'HLA-A*74:09', 'HLA-A*74:10', 'HLA-A*74:11', 'HLA-A*74:13', 'HLA-A*80:01',
                            'HLA-A*80:02', 'HLA-B*07:02', 'HLA-B*07:03', 'HLA-B*07:04', 'HLA-B*07:05', 'HLA-B*07:06', 'HLA-B*07:07', 'HLA-B*07:08',
                            'HLA-B*07:09', 'HLA-B*07:10', 'HLA-B*07:11', 'HLA-B*07:12', 'HLA-B*07:13', 'HLA-B*07:14', 'HLA-B*07:15', 'HLA-B*07:16',
                            'HLA-B*07:17', 'HLA-B*07:18', 'HLA-B*07:19', 'HLA-B*07:20', 'HLA-B*07:21', 'HLA-B*07:22', 'HLA-B*07:23', 'HLA-B*07:24',
                            'HLA-B*07:25', 'HLA-B*07:26', 'HLA-B*07:27', 'HLA-B*07:28', 'HLA-B*07:29', 'HLA-B*07:30', 'HLA-B*07:31', 'HLA-B*07:32',
                            'HLA-B*07:33', 'HLA-B*07:34', 'HLA-B*07:35', 'HLA-B*07:36', 'HLA-B*07:37', 'HLA-B*07:38', 'HLA-B*07:39', 'HLA-B*07:40',
                            'HLA-B*07:41', 'HLA-B*07:42', 'HLA-B*07:43', 'HLA-B*07:44', 'HLA-B*07:45', 'HLA-B*07:46', 'HLA-B*07:47', 'HLA-B*07:48',
                            'HLA-B*07:50', 'HLA-B*07:51', 'HLA-B*07:52', 'HLA-B*07:53', 'HLA-B*07:54', 'HLA-B*07:55', 'HLA-B*07:56', 'HLA-B*07:57',
                            'HLA-B*07:58', 'HLA-B*07:59', 'HLA-B*07:60', 'HLA-B*07:61', 'HLA-B*07:62', 'HLA-B*07:63', 'HLA-B*07:64', 'HLA-B*07:65',
                            'HLA-B*07:66', 'HLA-B*07:68', 'HLA-B*07:69', 'HLA-B*07:70', 'HLA-B*07:71', 'HLA-B*07:72', 'HLA-B*07:73', 'HLA-B*07:74',
                            'HLA-B*07:75', 'HLA-B*07:76', 'HLA-B*07:77', 'HLA-B*07:78', 'HLA-B*07:79', 'HLA-B*07:80', 'HLA-B*07:81', 'HLA-B*07:82',
                            'HLA-B*07:83', 'HLA-B*07:84', 'HLA-B*07:85', 'HLA-B*07:86', 'HLA-B*07:87', 'HLA-B*07:88', 'HLA-B*07:89', 'HLA-B*07:90',
                            'HLA-B*07:91', 'HLA-B*07:92', 'HLA-B*07:93', 'HLA-B*07:94', 'HLA-B*07:95', 'HLA-B*07:96', 'HLA-B*07:97', 'HLA-B*07:98',
                            'HLA-B*07:99', 'HLA-B*07:100', 'HLA-B*07:101', 'HLA-B*07:102', 'HLA-B*07:103', 'HLA-B*07:104', 'HLA-B*07:105',
                            'HLA-B*07:106', 'HLA-B*07:107', 'HLA-B*07:108', 'HLA-B*07:109', 'HLA-B*07:110', 'HLA-B*07:112', 'HLA-B*07:113',
                            'HLA-B*07:114', 'HLA-B*07:115', 'HLA-B*08:01', 'HLA-B*08:02', 'HLA-B*08:03', 'HLA-B*08:04', 'HLA-B*08:05', 'HLA-B*08:07',
                            'HLA-B*08:09', 'HLA-B*08:10', 'HLA-B*08:11', 'HLA-B*08:12', 'HLA-B*08:13', 'HLA-B*08:14', 'HLA-B*08:15', 'HLA-B*08:16',
                            'HLA-B*08:17', 'HLA-B*08:18', 'HLA-B*08:20', 'HLA-B*08:21', 'HLA-B*08:22', 'HLA-B*08:23', 'HLA-B*08:24', 'HLA-B*08:25',
                            'HLA-B*08:26', 'HLA-B*08:27', 'HLA-B*08:28', 'HLA-B*08:29', 'HLA-B*08:31', 'HLA-B*08:32', 'HLA-B*08:33', 'HLA-B*08:34',
                            'HLA-B*08:35', 'HLA-B*08:36', 'HLA-B*08:37', 'HLA-B*08:38', 'HLA-B*08:39', 'HLA-B*08:40', 'HLA-B*08:41', 'HLA-B*08:42',
                            'HLA-B*08:43', 'HLA-B*08:44', 'HLA-B*08:45', 'HLA-B*08:46', 'HLA-B*08:47', 'HLA-B*08:48', 'HLA-B*08:49', 'HLA-B*08:50',
                            'HLA-B*08:51', 'HLA-B*08:52', 'HLA-B*08:53', 'HLA-B*08:54', 'HLA-B*08:55', 'HLA-B*08:56', 'HLA-B*08:57', 'HLA-B*08:58',
                            'HLA-B*08:59', 'HLA-B*08:60', 'HLA-B*08:61', 'HLA-B*08:62', 'HLA-B*13:01', 'HLA-B*13:02', 'HLA-B*13:03', 'HLA-B*13:04',
                            'HLA-B*13:06', 'HLA-B*13:09', 'HLA-B*13:10', 'HLA-B*13:11', 'HLA-B*13:12', 'HLA-B*13:13', 'HLA-B*13:14', 'HLA-B*13:15',
                            'HLA-B*13:16', 'HLA-B*13:17', 'HLA-B*13:18', 'HLA-B*13:19', 'HLA-B*13:20', 'HLA-B*13:21', 'HLA-B*13:22', 'HLA-B*13:23',
                            'HLA-B*13:25', 'HLA-B*13:26', 'HLA-B*13:27', 'HLA-B*13:28', 'HLA-B*13:29', 'HLA-B*13:30', 'HLA-B*13:31', 'HLA-B*13:32',
                            'HLA-B*13:33', 'HLA-B*13:34', 'HLA-B*13:35', 'HLA-B*13:36', 'HLA-B*13:37', 'HLA-B*13:38', 'HLA-B*13:39', 'HLA-B*14:01',
                            'HLA-B*14:02', 'HLA-B*14:03', 'HLA-B*14:04', 'HLA-B*14:05', 'HLA-B*14:06', 'HLA-B*14:08', 'HLA-B*14:09', 'HLA-B*14:10',
                            'HLA-B*14:11', 'HLA-B*14:12', 'HLA-B*14:13', 'HLA-B*14:14', 'HLA-B*14:15', 'HLA-B*14:16', 'HLA-B*14:17', 'HLA-B*14:18',
                            'HLA-B*15:01', 'HLA-B*15:02', 'HLA-B*15:03', 'HLA-B*15:04', 'HLA-B*15:05', 'HLA-B*15:06', 'HLA-B*15:07', 'HLA-B*15:08',
                            'HLA-B*15:09', 'HLA-B*15:10', 'HLA-B*15:11', 'HLA-B*15:12', 'HLA-B*15:13', 'HLA-B*15:14', 'HLA-B*15:15', 'HLA-B*15:16',
                            'HLA-B*15:17', 'HLA-B*15:18', 'HLA-B*15:19', 'HLA-B*15:20', 'HLA-B*15:21', 'HLA-B*15:23', 'HLA-B*15:24', 'HLA-B*15:25',
                            'HLA-B*15:27', 'HLA-B*15:28', 'HLA-B*15:29', 'HLA-B*15:30', 'HLA-B*15:31', 'HLA-B*15:32', 'HLA-B*15:33', 'HLA-B*15:34',
                            'HLA-B*15:35', 'HLA-B*15:36', 'HLA-B*15:37', 'HLA-B*15:38', 'HLA-B*15:39', 'HLA-B*15:40', 'HLA-B*15:42', 'HLA-B*15:43',
                            'HLA-B*15:44', 'HLA-B*15:45', 'HLA-B*15:46', 'HLA-B*15:47', 'HLA-B*15:48', 'HLA-B*15:49', 'HLA-B*15:50', 'HLA-B*15:51',
                            'HLA-B*15:52', 'HLA-B*15:53', 'HLA-B*15:54', 'HLA-B*15:55', 'HLA-B*15:56', 'HLA-B*15:57', 'HLA-B*15:58', 'HLA-B*15:60',
                            'HLA-B*15:61', 'HLA-B*15:62', 'HLA-B*15:63', 'HLA-B*15:64', 'HLA-B*15:65', 'HLA-B*15:66', 'HLA-B*15:67', 'HLA-B*15:68',
                            'HLA-B*15:69', 'HLA-B*15:70', 'HLA-B*15:71', 'HLA-B*15:72', 'HLA-B*15:73', 'HLA-B*15:74', 'HLA-B*15:75', 'HLA-B*15:76',
                            'HLA-B*15:77', 'HLA-B*15:78', 'HLA-B*15:80', 'HLA-B*15:81', 'HLA-B*15:82', 'HLA-B*15:83', 'HLA-B*15:84', 'HLA-B*15:85',
                            'HLA-B*15:86', 'HLA-B*15:87', 'HLA-B*15:88', 'HLA-B*15:89', 'HLA-B*15:90', 'HLA-B*15:91', 'HLA-B*15:92', 'HLA-B*15:93',
                            'HLA-B*15:95', 'HLA-B*15:96', 'HLA-B*15:97', 'HLA-B*15:98', 'HLA-B*15:99', 'HLA-B*15:101', 'HLA-B*15:102', 'HLA-B*15:103',
                            'HLA-B*15:104', 'HLA-B*15:105', 'HLA-B*15:106', 'HLA-B*15:107', 'HLA-B*15:108', 'HLA-B*15:109', 'HLA-B*15:110',
                            'HLA-B*15:112', 'HLA-B*15:113', 'HLA-B*15:114', 'HLA-B*15:115', 'HLA-B*15:116', 'HLA-B*15:117', 'HLA-B*15:118',
                            'HLA-B*15:119', 'HLA-B*15:120', 'HLA-B*15:121', 'HLA-B*15:122', 'HLA-B*15:123', 'HLA-B*15:124', 'HLA-B*15:125',
                            'HLA-B*15:126', 'HLA-B*15:127', 'HLA-B*15:128', 'HLA-B*15:129', 'HLA-B*15:131', 'HLA-B*15:132', 'HLA-B*15:133',
                            'HLA-B*15:134', 'HLA-B*15:135', 'HLA-B*15:136', 'HLA-B*15:137', 'HLA-B*15:138', 'HLA-B*15:139', 'HLA-B*15:140',
                            'HLA-B*15:141', 'HLA-B*15:142', 'HLA-B*15:143', 'HLA-B*15:144', 'HLA-B*15:145', 'HLA-B*15:146', 'HLA-B*15:147',
                            'HLA-B*15:148', 'HLA-B*15:150', 'HLA-B*15:151', 'HLA-B*15:152', 'HLA-B*15:153', 'HLA-B*15:154', 'HLA-B*15:155',
                            'HLA-B*15:156', 'HLA-B*15:157', 'HLA-B*15:158', 'HLA-B*15:159', 'HLA-B*15:160', 'HLA-B*15:161', 'HLA-B*15:162',
                            'HLA-B*15:163', 'HLA-B*15:164', 'HLA-B*15:165', 'HLA-B*15:166', 'HLA-B*15:167', 'HLA-B*15:168', 'HLA-B*15:169',
                            'HLA-B*15:170', 'HLA-B*15:171', 'HLA-B*15:172', 'HLA-B*15:173', 'HLA-B*15:174', 'HLA-B*15:175', 'HLA-B*15:176',
                            'HLA-B*15:177', 'HLA-B*15:178', 'HLA-B*15:179', 'HLA-B*15:180', 'HLA-B*15:183', 'HLA-B*15:184', 'HLA-B*15:185',
                            'HLA-B*15:186', 'HLA-B*15:187', 'HLA-B*15:188', 'HLA-B*15:189', 'HLA-B*15:191', 'HLA-B*15:192', 'HLA-B*15:193',
                            'HLA-B*15:194', 'HLA-B*15:195', 'HLA-B*15:196', 'HLA-B*15:197', 'HLA-B*15:198', 'HLA-B*15:199', 'HLA-B*15:200',
                            'HLA-B*15:201', 'HLA-B*15:202', 'HLA-B*18:01', 'HLA-B*18:02', 'HLA-B*18:03', 'HLA-B*18:04', 'HLA-B*18:05', 'HLA-B*18:06',
                            'HLA-B*18:07', 'HLA-B*18:08', 'HLA-B*18:09', 'HLA-B*18:10', 'HLA-B*18:11', 'HLA-B*18:12', 'HLA-B*18:13', 'HLA-B*18:14',
                            'HLA-B*18:15', 'HLA-B*18:18', 'HLA-B*18:19', 'HLA-B*18:20', 'HLA-B*18:21', 'HLA-B*18:22', 'HLA-B*18:24', 'HLA-B*18:25',
                            'HLA-B*18:26', 'HLA-B*18:27', 'HLA-B*18:28', 'HLA-B*18:29', 'HLA-B*18:30', 'HLA-B*18:31', 'HLA-B*18:32', 'HLA-B*18:33',
                            'HLA-B*18:34', 'HLA-B*18:35', 'HLA-B*18:36', 'HLA-B*18:37', 'HLA-B*18:38', 'HLA-B*18:39', 'HLA-B*18:40', 'HLA-B*18:41',
                            'HLA-B*18:42', 'HLA-B*18:43', 'HLA-B*18:44', 'HLA-B*18:45', 'HLA-B*18:46', 'HLA-B*18:47', 'HLA-B*18:48', 'HLA-B*18:49',
                            'HLA-B*18:50', 'HLA-B*27:01', 'HLA-B*27:02', 'HLA-B*27:03', 'HLA-B*27:04', 'HLA-B*27:05', 'HLA-B*27:06', 'HLA-B*27:07',
                            'HLA-B*27:08', 'HLA-B*27:09', 'HLA-B*27:10', 'HLA-B*27:11', 'HLA-B*27:12', 'HLA-B*27:13', 'HLA-B*27:14', 'HLA-B*27:15',
                            'HLA-B*27:16', 'HLA-B*27:17', 'HLA-B*27:18', 'HLA-B*27:19', 'HLA-B*27:20', 'HLA-B*27:21', 'HLA-B*27:23', 'HLA-B*27:24',
                            'HLA-B*27:25', 'HLA-B*27:26', 'HLA-B*27:27', 'HLA-B*27:28', 'HLA-B*27:29', 'HLA-B*27:30', 'HLA-B*27:31', 'HLA-B*27:32',
                            'HLA-B*27:33', 'HLA-B*27:34', 'HLA-B*27:35', 'HLA-B*27:36', 'HLA-B*27:37', 'HLA-B*27:38', 'HLA-B*27:39', 'HLA-B*27:40',
                            'HLA-B*27:41', 'HLA-B*27:42', 'HLA-B*27:43', 'HLA-B*27:44', 'HLA-B*27:45', 'HLA-B*27:46', 'HLA-B*27:47', 'HLA-B*27:48',
                            'HLA-B*27:49', 'HLA-B*27:50', 'HLA-B*27:51', 'HLA-B*27:52', 'HLA-B*27:53', 'HLA-B*27:54', 'HLA-B*27:55', 'HLA-B*27:56',
                            'HLA-B*27:57', 'HLA-B*27:58', 'HLA-B*27:60', 'HLA-B*27:61', 'HLA-B*27:62', 'HLA-B*27:63', 'HLA-B*27:67', 'HLA-B*27:68',
                            'HLA-B*27:69', 'HLA-B*35:01', 'HLA-B*35:02', 'HLA-B*35:03', 'HLA-B*35:04', 'HLA-B*35:05', 'HLA-B*35:06', 'HLA-B*35:07',
                            'HLA-B*35:08', 'HLA-B*35:09', 'HLA-B*35:10', 'HLA-B*35:11', 'HLA-B*35:12', 'HLA-B*35:13', 'HLA-B*35:14', 'HLA-B*35:15',
                            'HLA-B*35:16', 'HLA-B*35:17', 'HLA-B*35:18', 'HLA-B*35:19', 'HLA-B*35:20', 'HLA-B*35:21', 'HLA-B*35:22', 'HLA-B*35:23',
                            'HLA-B*35:24', 'HLA-B*35:25', 'HLA-B*35:26', 'HLA-B*35:27', 'HLA-B*35:28', 'HLA-B*35:29', 'HLA-B*35:30', 'HLA-B*35:31',
                            'HLA-B*35:32', 'HLA-B*35:33', 'HLA-B*35:34', 'HLA-B*35:35', 'HLA-B*35:36', 'HLA-B*35:37', 'HLA-B*35:38', 'HLA-B*35:39',
                            'HLA-B*35:41', 'HLA-B*35:42', 'HLA-B*35:43', 'HLA-B*35:44', 'HLA-B*35:45', 'HLA-B*35:46', 'HLA-B*35:47', 'HLA-B*35:48',
                            'HLA-B*35:49', 'HLA-B*35:50', 'HLA-B*35:51', 'HLA-B*35:52', 'HLA-B*35:54', 'HLA-B*35:55', 'HLA-B*35:56', 'HLA-B*35:57',
                            'HLA-B*35:58', 'HLA-B*35:59', 'HLA-B*35:60', 'HLA-B*35:61', 'HLA-B*35:62', 'HLA-B*35:63', 'HLA-B*35:64', 'HLA-B*35:66',
                            'HLA-B*35:67', 'HLA-B*35:68', 'HLA-B*35:69', 'HLA-B*35:70', 'HLA-B*35:71', 'HLA-B*35:72', 'HLA-B*35:74', 'HLA-B*35:75',
                            'HLA-B*35:76', 'HLA-B*35:77', 'HLA-B*35:78', 'HLA-B*35:79', 'HLA-B*35:80', 'HLA-B*35:81', 'HLA-B*35:82', 'HLA-B*35:83',
                            'HLA-B*35:84', 'HLA-B*35:85', 'HLA-B*35:86', 'HLA-B*35:87', 'HLA-B*35:88', 'HLA-B*35:89', 'HLA-B*35:90', 'HLA-B*35:91',
                            'HLA-B*35:92', 'HLA-B*35:93', 'HLA-B*35:94', 'HLA-B*35:95', 'HLA-B*35:96', 'HLA-B*35:97', 'HLA-B*35:98', 'HLA-B*35:99',
                            'HLA-B*35:100', 'HLA-B*35:101', 'HLA-B*35:102', 'HLA-B*35:103', 'HLA-B*35:104', 'HLA-B*35:105', 'HLA-B*35:106',
                            'HLA-B*35:107', 'HLA-B*35:108', 'HLA-B*35:109', 'HLA-B*35:110', 'HLA-B*35:111', 'HLA-B*35:112', 'HLA-B*35:113',
                            'HLA-B*35:114', 'HLA-B*35:115', 'HLA-B*35:116', 'HLA-B*35:117', 'HLA-B*35:118', 'HLA-B*35:119', 'HLA-B*35:120',
                            'HLA-B*35:121', 'HLA-B*35:122', 'HLA-B*35:123', 'HLA-B*35:124', 'HLA-B*35:125', 'HLA-B*35:126', 'HLA-B*35:127',
                            'HLA-B*35:128', 'HLA-B*35:131', 'HLA-B*35:132', 'HLA-B*35:133', 'HLA-B*35:135', 'HLA-B*35:136', 'HLA-B*35:137',
                            'HLA-B*35:138', 'HLA-B*35:139', 'HLA-B*35:140', 'HLA-B*35:141', 'HLA-B*35:142', 'HLA-B*35:143', 'HLA-B*35:144',
                            'HLA-B*37:01', 'HLA-B*37:02', 'HLA-B*37:04', 'HLA-B*37:05', 'HLA-B*37:06', 'HLA-B*37:07', 'HLA-B*37:08', 'HLA-B*37:09',
                            'HLA-B*37:10', 'HLA-B*37:11', 'HLA-B*37:12', 'HLA-B*37:13', 'HLA-B*37:14', 'HLA-B*37:15', 'HLA-B*37:17', 'HLA-B*37:18',
                            'HLA-B*37:19', 'HLA-B*37:20', 'HLA-B*37:21', 'HLA-B*37:22', 'HLA-B*37:23', 'HLA-B*38:01', 'HLA-B*38:02', 'HLA-B*38:03',
                            'HLA-B*38:04', 'HLA-B*38:05', 'HLA-B*38:06', 'HLA-B*38:07', 'HLA-B*38:08', 'HLA-B*38:09', 'HLA-B*38:10', 'HLA-B*38:11',
                            'HLA-B*38:12', 'HLA-B*38:13', 'HLA-B*38:14', 'HLA-B*38:15', 'HLA-B*38:16', 'HLA-B*38:17', 'HLA-B*38:18', 'HLA-B*38:19',
                            'HLA-B*38:20', 'HLA-B*38:21', 'HLA-B*38:22', 'HLA-B*38:23', 'HLA-B*39:01', 'HLA-B*39:02', 'HLA-B*39:03', 'HLA-B*39:04',
                            'HLA-B*39:05', 'HLA-B*39:06', 'HLA-B*39:07', 'HLA-B*39:08', 'HLA-B*39:09', 'HLA-B*39:10', 'HLA-B*39:11', 'HLA-B*39:12',
                            'HLA-B*39:13', 'HLA-B*39:14', 'HLA-B*39:15', 'HLA-B*39:16', 'HLA-B*39:17', 'HLA-B*39:18', 'HLA-B*39:19', 'HLA-B*39:20',
                            'HLA-B*39:22', 'HLA-B*39:23', 'HLA-B*39:24', 'HLA-B*39:26', 'HLA-B*39:27', 'HLA-B*39:28', 'HLA-B*39:29', 'HLA-B*39:30',
                            'HLA-B*39:31', 'HLA-B*39:32', 'HLA-B*39:33', 'HLA-B*39:34', 'HLA-B*39:35', 'HLA-B*39:36', 'HLA-B*39:37', 'HLA-B*39:39',
                            'HLA-B*39:41', 'HLA-B*39:42', 'HLA-B*39:43', 'HLA-B*39:44', 'HLA-B*39:45', 'HLA-B*39:46', 'HLA-B*39:47', 'HLA-B*39:48',
                            'HLA-B*39:49', 'HLA-B*39:50', 'HLA-B*39:51', 'HLA-B*39:52', 'HLA-B*39:53', 'HLA-B*39:54', 'HLA-B*39:55', 'HLA-B*39:56',
                            'HLA-B*39:57', 'HLA-B*39:58', 'HLA-B*39:59', 'HLA-B*39:60', 'HLA-B*40:01', 'HLA-B*40:02', 'HLA-B*40:03', 'HLA-B*40:04',
                            'HLA-B*40:05', 'HLA-B*40:06', 'HLA-B*40:07', 'HLA-B*40:08', 'HLA-B*40:09', 'HLA-B*40:10', 'HLA-B*40:11', 'HLA-B*40:12',
                            'HLA-B*40:13', 'HLA-B*40:14', 'HLA-B*40:15', 'HLA-B*40:16', 'HLA-B*40:18', 'HLA-B*40:19', 'HLA-B*40:20', 'HLA-B*40:21',
                            'HLA-B*40:23', 'HLA-B*40:24', 'HLA-B*40:25', 'HLA-B*40:26', 'HLA-B*40:27', 'HLA-B*40:28', 'HLA-B*40:29', 'HLA-B*40:30',
                            'HLA-B*40:31', 'HLA-B*40:32', 'HLA-B*40:33', 'HLA-B*40:34', 'HLA-B*40:35', 'HLA-B*40:36', 'HLA-B*40:37', 'HLA-B*40:38',
                            'HLA-B*40:39', 'HLA-B*40:40', 'HLA-B*40:42', 'HLA-B*40:43', 'HLA-B*40:44', 'HLA-B*40:45', 'HLA-B*40:46', 'HLA-B*40:47',
                            'HLA-B*40:48', 'HLA-B*40:49', 'HLA-B*40:50', 'HLA-B*40:51', 'HLA-B*40:52', 'HLA-B*40:53', 'HLA-B*40:54', 'HLA-B*40:55',
                            'HLA-B*40:56', 'HLA-B*40:57', 'HLA-B*40:58', 'HLA-B*40:59', 'HLA-B*40:60', 'HLA-B*40:61', 'HLA-B*40:62', 'HLA-B*40:63',
                            'HLA-B*40:64', 'HLA-B*40:65', 'HLA-B*40:66', 'HLA-B*40:67', 'HLA-B*40:68', 'HLA-B*40:69', 'HLA-B*40:70', 'HLA-B*40:71',
                            'HLA-B*40:72', 'HLA-B*40:73', 'HLA-B*40:74', 'HLA-B*40:75', 'HLA-B*40:76', 'HLA-B*40:77', 'HLA-B*40:78', 'HLA-B*40:79',
                            'HLA-B*40:80', 'HLA-B*40:81', 'HLA-B*40:82', 'HLA-B*40:83', 'HLA-B*40:84', 'HLA-B*40:85', 'HLA-B*40:86', 'HLA-B*40:87',
                            'HLA-B*40:88', 'HLA-B*40:89', 'HLA-B*40:90', 'HLA-B*40:91', 'HLA-B*40:92', 'HLA-B*40:93', 'HLA-B*40:94', 'HLA-B*40:95',
                            'HLA-B*40:96', 'HLA-B*40:97', 'HLA-B*40:98', 'HLA-B*40:99', 'HLA-B*40:100', 'HLA-B*40:101', 'HLA-B*40:102', 'HLA-B*40:103',
                            'HLA-B*40:104', 'HLA-B*40:105', 'HLA-B*40:106', 'HLA-B*40:107', 'HLA-B*40:108', 'HLA-B*40:109', 'HLA-B*40:110',
                            'HLA-B*40:111', 'HLA-B*40:112', 'HLA-B*40:113', 'HLA-B*40:114', 'HLA-B*40:115', 'HLA-B*40:116', 'HLA-B*40:117',
                            'HLA-B*40:119', 'HLA-B*40:120', 'HLA-B*40:121', 'HLA-B*40:122', 'HLA-B*40:123', 'HLA-B*40:124', 'HLA-B*40:125',
                            'HLA-B*40:126', 'HLA-B*40:127', 'HLA-B*40:128', 'HLA-B*40:129', 'HLA-B*40:130', 'HLA-B*40:131', 'HLA-B*40:132',
                            'HLA-B*40:134', 'HLA-B*40:135', 'HLA-B*40:136', 'HLA-B*40:137', 'HLA-B*40:138', 'HLA-B*40:139', 'HLA-B*40:140',
                            'HLA-B*40:141', 'HLA-B*40:143', 'HLA-B*40:145', 'HLA-B*40:146', 'HLA-B*40:147', 'HLA-B*41:01', 'HLA-B*41:02', 'HLA-B*41:03',
                            'HLA-B*41:04', 'HLA-B*41:05', 'HLA-B*41:06', 'HLA-B*41:07', 'HLA-B*41:08', 'HLA-B*41:09', 'HLA-B*41:10', 'HLA-B*41:11',
                            'HLA-B*41:12', 'HLA-B*42:01', 'HLA-B*42:02', 'HLA-B*42:04', 'HLA-B*42:05', 'HLA-B*42:06', 'HLA-B*42:07', 'HLA-B*42:08',
                            'HLA-B*42:09', 'HLA-B*42:10', 'HLA-B*42:11', 'HLA-B*42:12', 'HLA-B*42:13', 'HLA-B*42:14', 'HLA-B*44:02', 'HLA-B*44:03',
                            'HLA-B*44:04', 'HLA-B*44:05', 'HLA-B*44:06', 'HLA-B*44:07', 'HLA-B*44:08', 'HLA-B*44:09', 'HLA-B*44:10', 'HLA-B*44:11',
                            'HLA-B*44:12', 'HLA-B*44:13', 'HLA-B*44:14', 'HLA-B*44:15', 'HLA-B*44:16', 'HLA-B*44:17', 'HLA-B*44:18', 'HLA-B*44:20',
                            'HLA-B*44:21', 'HLA-B*44:22', 'HLA-B*44:24', 'HLA-B*44:25', 'HLA-B*44:26', 'HLA-B*44:27', 'HLA-B*44:28', 'HLA-B*44:29',
                            'HLA-B*44:30', 'HLA-B*44:31', 'HLA-B*44:32', 'HLA-B*44:33', 'HLA-B*44:34', 'HLA-B*44:35', 'HLA-B*44:36', 'HLA-B*44:37',
                            'HLA-B*44:38', 'HLA-B*44:39', 'HLA-B*44:40', 'HLA-B*44:41', 'HLA-B*44:42', 'HLA-B*44:43', 'HLA-B*44:44', 'HLA-B*44:45',
                            'HLA-B*44:46', 'HLA-B*44:47', 'HLA-B*44:48', 'HLA-B*44:49', 'HLA-B*44:50', 'HLA-B*44:51', 'HLA-B*44:53', 'HLA-B*44:54',
                            'HLA-B*44:55', 'HLA-B*44:57', 'HLA-B*44:59', 'HLA-B*44:60', 'HLA-B*44:62', 'HLA-B*44:63', 'HLA-B*44:64', 'HLA-B*44:65',
                            'HLA-B*44:66', 'HLA-B*44:67', 'HLA-B*44:68', 'HLA-B*44:69', 'HLA-B*44:70', 'HLA-B*44:71', 'HLA-B*44:72', 'HLA-B*44:73',
                            'HLA-B*44:74', 'HLA-B*44:75', 'HLA-B*44:76', 'HLA-B*44:77', 'HLA-B*44:78', 'HLA-B*44:79', 'HLA-B*44:80', 'HLA-B*44:81',
                            'HLA-B*44:82', 'HLA-B*44:83', 'HLA-B*44:84', 'HLA-B*44:85', 'HLA-B*44:86', 'HLA-B*44:87', 'HLA-B*44:88', 'HLA-B*44:89',
                            'HLA-B*44:90', 'HLA-B*44:91', 'HLA-B*44:92', 'HLA-B*44:93', 'HLA-B*44:94', 'HLA-B*44:95', 'HLA-B*44:96', 'HLA-B*44:97',
                            'HLA-B*44:98', 'HLA-B*44:99', 'HLA-B*44:100', 'HLA-B*44:101', 'HLA-B*44:102', 'HLA-B*44:103', 'HLA-B*44:104',
                            'HLA-B*44:105', 'HLA-B*44:106', 'HLA-B*44:107', 'HLA-B*44:109', 'HLA-B*44:110', 'HLA-B*45:01', 'HLA-B*45:02', 'HLA-B*45:03',
                            'HLA-B*45:04', 'HLA-B*45:05', 'HLA-B*45:06', 'HLA-B*45:07', 'HLA-B*45:08', 'HLA-B*45:09', 'HLA-B*45:10', 'HLA-B*45:11',
                            'HLA-B*45:12', 'HLA-B*46:01', 'HLA-B*46:02', 'HLA-B*46:03', 'HLA-B*46:04', 'HLA-B*46:05', 'HLA-B*46:06', 'HLA-B*46:08',
                            'HLA-B*46:09', 'HLA-B*46:10', 'HLA-B*46:11', 'HLA-B*46:12', 'HLA-B*46:13', 'HLA-B*46:14', 'HLA-B*46:16', 'HLA-B*46:17',
                            'HLA-B*46:18', 'HLA-B*46:19', 'HLA-B*46:20', 'HLA-B*46:21', 'HLA-B*46:22', 'HLA-B*46:23', 'HLA-B*46:24', 'HLA-B*47:01',
                            'HLA-B*47:02', 'HLA-B*47:03', 'HLA-B*47:04', 'HLA-B*47:05', 'HLA-B*47:06', 'HLA-B*47:07', 'HLA-B*48:01', 'HLA-B*48:02',
                            'HLA-B*48:03', 'HLA-B*48:04', 'HLA-B*48:05', 'HLA-B*48:06', 'HLA-B*48:07', 'HLA-B*48:08', 'HLA-B*48:09', 'HLA-B*48:10',
                            'HLA-B*48:11', 'HLA-B*48:12', 'HLA-B*48:13', 'HLA-B*48:14', 'HLA-B*48:15', 'HLA-B*48:16', 'HLA-B*48:17', 'HLA-B*48:18',
                            'HLA-B*48:19', 'HLA-B*48:20', 'HLA-B*48:21', 'HLA-B*48:22', 'HLA-B*48:23', 'HLA-B*49:01', 'HLA-B*49:02', 'HLA-B*49:03',
                            'HLA-B*49:04', 'HLA-B*49:05', 'HLA-B*49:06', 'HLA-B*49:07', 'HLA-B*49:08', 'HLA-B*49:09', 'HLA-B*49:10', 'HLA-B*50:01',
                            'HLA-B*50:02', 'HLA-B*50:04', 'HLA-B*50:05', 'HLA-B*50:06', 'HLA-B*50:07', 'HLA-B*50:08', 'HLA-B*50:09', 'HLA-B*51:01',
                            'HLA-B*51:02', 'HLA-B*51:03', 'HLA-B*51:04', 'HLA-B*51:05', 'HLA-B*51:06', 'HLA-B*51:07', 'HLA-B*51:08', 'HLA-B*51:09',
                            'HLA-B*51:12', 'HLA-B*51:13', 'HLA-B*51:14', 'HLA-B*51:15', 'HLA-B*51:16', 'HLA-B*51:17', 'HLA-B*51:18', 'HLA-B*51:19',
                            'HLA-B*51:20', 'HLA-B*51:21', 'HLA-B*51:22', 'HLA-B*51:23', 'HLA-B*51:24', 'HLA-B*51:26', 'HLA-B*51:28', 'HLA-B*51:29',
                            'HLA-B*51:30', 'HLA-B*51:31', 'HLA-B*51:32', 'HLA-B*51:33', 'HLA-B*51:34', 'HLA-B*51:35', 'HLA-B*51:36', 'HLA-B*51:37',
                            'HLA-B*51:38', 'HLA-B*51:39', 'HLA-B*51:40', 'HLA-B*51:42', 'HLA-B*51:43', 'HLA-B*51:45', 'HLA-B*51:46', 'HLA-B*51:48',
                            'HLA-B*51:49', 'HLA-B*51:50', 'HLA-B*51:51', 'HLA-B*51:52', 'HLA-B*51:53', 'HLA-B*51:54', 'HLA-B*51:55', 'HLA-B*51:56',
                            'HLA-B*51:57', 'HLA-B*51:58', 'HLA-B*51:59', 'HLA-B*51:60', 'HLA-B*51:61', 'HLA-B*51:62', 'HLA-B*51:63', 'HLA-B*51:64',
                            'HLA-B*51:65', 'HLA-B*51:66', 'HLA-B*51:67', 'HLA-B*51:68', 'HLA-B*51:69', 'HLA-B*51:70', 'HLA-B*51:71', 'HLA-B*51:72',
                            'HLA-B*51:73', 'HLA-B*51:74', 'HLA-B*51:75', 'HLA-B*51:76', 'HLA-B*51:77', 'HLA-B*51:78', 'HLA-B*51:79', 'HLA-B*51:80',
                            'HLA-B*51:81', 'HLA-B*51:82', 'HLA-B*51:83', 'HLA-B*51:84', 'HLA-B*51:85', 'HLA-B*51:86', 'HLA-B*51:87', 'HLA-B*51:88',
                            'HLA-B*51:89', 'HLA-B*51:90', 'HLA-B*51:91', 'HLA-B*51:92', 'HLA-B*51:93', 'HLA-B*51:94', 'HLA-B*51:95', 'HLA-B*51:96',
                            'HLA-B*52:01', 'HLA-B*52:02', 'HLA-B*52:03', 'HLA-B*52:04', 'HLA-B*52:05', 'HLA-B*52:06', 'HLA-B*52:07', 'HLA-B*52:08',
                            'HLA-B*52:09', 'HLA-B*52:10', 'HLA-B*52:11', 'HLA-B*52:12', 'HLA-B*52:13', 'HLA-B*52:14', 'HLA-B*52:15', 'HLA-B*52:16',
                            'HLA-B*52:17', 'HLA-B*52:18', 'HLA-B*52:19', 'HLA-B*52:20', 'HLA-B*52:21', 'HLA-B*53:01', 'HLA-B*53:02', 'HLA-B*53:03',
                            'HLA-B*53:04', 'HLA-B*53:05', 'HLA-B*53:06', 'HLA-B*53:07', 'HLA-B*53:08', 'HLA-B*53:09', 'HLA-B*53:10', 'HLA-B*53:11',
                            'HLA-B*53:12', 'HLA-B*53:13', 'HLA-B*53:14', 'HLA-B*53:15', 'HLA-B*53:16', 'HLA-B*53:17', 'HLA-B*53:18', 'HLA-B*53:19',
                            'HLA-B*53:20', 'HLA-B*53:21', 'HLA-B*53:22', 'HLA-B*53:23', 'HLA-B*54:01', 'HLA-B*54:02', 'HLA-B*54:03', 'HLA-B*54:04',
                            'HLA-B*54:06', 'HLA-B*54:07', 'HLA-B*54:09', 'HLA-B*54:10', 'HLA-B*54:11', 'HLA-B*54:12', 'HLA-B*54:13', 'HLA-B*54:14',
                            'HLA-B*54:15', 'HLA-B*54:16', 'HLA-B*54:17', 'HLA-B*54:18', 'HLA-B*54:19', 'HLA-B*54:20', 'HLA-B*54:21', 'HLA-B*54:22',
                            'HLA-B*54:23', 'HLA-B*55:01', 'HLA-B*55:02', 'HLA-B*55:03', 'HLA-B*55:04', 'HLA-B*55:05', 'HLA-B*55:07', 'HLA-B*55:08',
                            'HLA-B*55:09', 'HLA-B*55:10', 'HLA-B*55:11', 'HLA-B*55:12', 'HLA-B*55:13', 'HLA-B*55:14', 'HLA-B*55:15', 'HLA-B*55:16',
                            'HLA-B*55:17', 'HLA-B*55:18', 'HLA-B*55:19', 'HLA-B*55:20', 'HLA-B*55:21', 'HLA-B*55:22', 'HLA-B*55:23', 'HLA-B*55:24',
                            'HLA-B*55:25', 'HLA-B*55:26', 'HLA-B*55:27', 'HLA-B*55:28', 'HLA-B*55:29', 'HLA-B*55:30', 'HLA-B*55:31', 'HLA-B*55:32',
                            'HLA-B*55:33', 'HLA-B*55:34', 'HLA-B*55:35', 'HLA-B*55:36', 'HLA-B*55:37', 'HLA-B*55:38', 'HLA-B*55:39', 'HLA-B*55:40',
                            'HLA-B*55:41', 'HLA-B*55:42', 'HLA-B*55:43', 'HLA-B*56:01', 'HLA-B*56:02', 'HLA-B*56:03', 'HLA-B*56:04', 'HLA-B*56:05',
                            'HLA-B*56:06', 'HLA-B*56:07', 'HLA-B*56:08', 'HLA-B*56:09', 'HLA-B*56:10', 'HLA-B*56:11', 'HLA-B*56:12', 'HLA-B*56:13',
                            'HLA-B*56:14', 'HLA-B*56:15', 'HLA-B*56:16', 'HLA-B*56:17', 'HLA-B*56:18', 'HLA-B*56:20', 'HLA-B*56:21', 'HLA-B*56:22',
                            'HLA-B*56:23', 'HLA-B*56:24', 'HLA-B*56:25', 'HLA-B*56:26', 'HLA-B*56:27', 'HLA-B*56:29', 'HLA-B*57:01', 'HLA-B*57:02',
                            'HLA-B*57:03', 'HLA-B*57:04', 'HLA-B*57:05', 'HLA-B*57:06', 'HLA-B*57:07', 'HLA-B*57:08', 'HLA-B*57:09', 'HLA-B*57:10',
                            'HLA-B*57:11', 'HLA-B*57:12', 'HLA-B*57:13', 'HLA-B*57:14', 'HLA-B*57:15', 'HLA-B*57:16', 'HLA-B*57:17', 'HLA-B*57:18',
                            'HLA-B*57:19', 'HLA-B*57:20', 'HLA-B*57:21', 'HLA-B*57:22', 'HLA-B*57:23', 'HLA-B*57:24', 'HLA-B*57:25', 'HLA-B*57:26',
                            'HLA-B*57:27', 'HLA-B*57:29', 'HLA-B*57:30', 'HLA-B*57:31', 'HLA-B*57:32', 'HLA-B*58:01', 'HLA-B*58:02', 'HLA-B*58:04',
                            'HLA-B*58:05', 'HLA-B*58:06', 'HLA-B*58:07', 'HLA-B*58:08', 'HLA-B*58:09', 'HLA-B*58:11', 'HLA-B*58:12', 'HLA-B*58:13',
                            'HLA-B*58:14', 'HLA-B*58:15', 'HLA-B*58:16', 'HLA-B*58:18', 'HLA-B*58:19', 'HLA-B*58:20', 'HLA-B*58:21', 'HLA-B*58:22',
                            'HLA-B*58:23', 'HLA-B*58:24', 'HLA-B*58:25', 'HLA-B*58:26', 'HLA-B*58:27', 'HLA-B*58:28', 'HLA-B*58:29', 'HLA-B*58:30',
                            'HLA-B*59:01', 'HLA-B*59:02', 'HLA-B*59:03', 'HLA-B*59:04', 'HLA-B*59:05', 'HLA-B*67:01', 'HLA-B*67:02', 'HLA-B*73:01',
                            'HLA-B*73:02', 'HLA-B*78:01', 'HLA-B*78:02', 'HLA-B*78:03', 'HLA-B*78:04', 'HLA-B*78:05', 'HLA-B*78:06', 'HLA-B*78:07',
                            'HLA-B*81:01', 'HLA-B*81:02', 'HLA-B*81:03', 'HLA-B*81:05', 'HLA-B*82:01', 'HLA-B*82:02', 'HLA-B*82:03', 'HLA-B*83:01',
                            'HLA-C*01:02', 'HLA-C*01:03', 'HLA-C*01:04', 'HLA-C*01:05', 'HLA-C*01:06', 'HLA-C*01:07', 'HLA-C*01:08', 'HLA-C*01:09',
                            'HLA-C*01:10', 'HLA-C*01:11', 'HLA-C*01:12', 'HLA-C*01:13', 'HLA-C*01:14', 'HLA-C*01:15', 'HLA-C*01:16', 'HLA-C*01:17',
                            'HLA-C*01:18', 'HLA-C*01:19', 'HLA-C*01:20', 'HLA-C*01:21', 'HLA-C*01:22', 'HLA-C*01:23', 'HLA-C*01:24', 'HLA-C*01:25',
                            'HLA-C*01:26', 'HLA-C*01:27', 'HLA-C*01:28', 'HLA-C*01:29', 'HLA-C*01:30', 'HLA-C*01:31', 'HLA-C*01:32', 'HLA-C*01:33',
                            'HLA-C*01:34', 'HLA-C*01:35', 'HLA-C*01:36', 'HLA-C*01:38', 'HLA-C*01:39', 'HLA-C*01:40', 'HLA-C*02:02', 'HLA-C*02:03',
                            'HLA-C*02:04', 'HLA-C*02:05', 'HLA-C*02:06', 'HLA-C*02:07', 'HLA-C*02:08', 'HLA-C*02:09', 'HLA-C*02:10', 'HLA-C*02:11',
                            'HLA-C*02:12', 'HLA-C*02:13', 'HLA-C*02:14', 'HLA-C*02:15', 'HLA-C*02:16', 'HLA-C*02:17', 'HLA-C*02:18', 'HLA-C*02:19',
                            'HLA-C*02:20', 'HLA-C*02:21', 'HLA-C*02:22', 'HLA-C*02:23', 'HLA-C*02:24', 'HLA-C*02:26', 'HLA-C*02:27', 'HLA-C*02:28',
                            'HLA-C*02:29', 'HLA-C*02:30', 'HLA-C*02:31', 'HLA-C*02:32', 'HLA-C*02:33', 'HLA-C*02:34', 'HLA-C*02:35', 'HLA-C*02:36',
                            'HLA-C*02:37', 'HLA-C*02:39', 'HLA-C*02:40', 'HLA-C*03:01', 'HLA-C*03:02', 'HLA-C*03:03', 'HLA-C*03:04', 'HLA-C*03:05',
                            'HLA-C*03:06', 'HLA-C*03:07', 'HLA-C*03:08', 'HLA-C*03:09', 'HLA-C*03:10', 'HLA-C*03:11', 'HLA-C*03:12', 'HLA-C*03:13',
                            'HLA-C*03:14', 'HLA-C*03:15', 'HLA-C*03:16', 'HLA-C*03:17', 'HLA-C*03:18', 'HLA-C*03:19', 'HLA-C*03:21', 'HLA-C*03:23',
                            'HLA-C*03:24', 'HLA-C*03:25', 'HLA-C*03:26', 'HLA-C*03:27', 'HLA-C*03:28', 'HLA-C*03:29', 'HLA-C*03:30', 'HLA-C*03:31',
                            'HLA-C*03:32', 'HLA-C*03:33', 'HLA-C*03:34', 'HLA-C*03:35', 'HLA-C*03:36', 'HLA-C*03:37', 'HLA-C*03:38', 'HLA-C*03:39',
                            'HLA-C*03:40', 'HLA-C*03:41', 'HLA-C*03:42', 'HLA-C*03:43', 'HLA-C*03:44', 'HLA-C*03:45', 'HLA-C*03:46', 'HLA-C*03:47',
                            'HLA-C*03:48', 'HLA-C*03:49', 'HLA-C*03:50', 'HLA-C*03:51', 'HLA-C*03:52', 'HLA-C*03:53', 'HLA-C*03:54', 'HLA-C*03:55',
                            'HLA-C*03:56', 'HLA-C*03:57', 'HLA-C*03:58', 'HLA-C*03:59', 'HLA-C*03:60', 'HLA-C*03:61', 'HLA-C*03:62', 'HLA-C*03:63',
                            'HLA-C*03:64', 'HLA-C*03:65', 'HLA-C*03:66', 'HLA-C*03:67', 'HLA-C*03:68', 'HLA-C*03:69', 'HLA-C*03:70', 'HLA-C*03:71',
                            'HLA-C*03:72', 'HLA-C*03:73', 'HLA-C*03:74', 'HLA-C*03:75', 'HLA-C*03:76', 'HLA-C*03:77', 'HLA-C*03:78', 'HLA-C*03:79',
                            'HLA-C*03:80', 'HLA-C*03:81', 'HLA-C*03:82', 'HLA-C*03:83', 'HLA-C*03:84', 'HLA-C*03:85', 'HLA-C*03:86', 'HLA-C*03:87',
                            'HLA-C*03:88', 'HLA-C*03:89', 'HLA-C*03:90', 'HLA-C*03:91', 'HLA-C*03:92', 'HLA-C*03:93', 'HLA-C*03:94', 'HLA-C*04:01',
                            'HLA-C*04:03', 'HLA-C*04:04', 'HLA-C*04:05', 'HLA-C*04:06', 'HLA-C*04:07', 'HLA-C*04:08', 'HLA-C*04:10', 'HLA-C*04:11',
                            'HLA-C*04:12', 'HLA-C*04:13', 'HLA-C*04:14', 'HLA-C*04:15', 'HLA-C*04:16', 'HLA-C*04:17', 'HLA-C*04:18', 'HLA-C*04:19',
                            'HLA-C*04:20', 'HLA-C*04:23', 'HLA-C*04:24', 'HLA-C*04:25', 'HLA-C*04:26', 'HLA-C*04:27', 'HLA-C*04:28', 'HLA-C*04:29',
                            'HLA-C*04:30', 'HLA-C*04:31', 'HLA-C*04:32', 'HLA-C*04:33', 'HLA-C*04:34', 'HLA-C*04:35', 'HLA-C*04:36', 'HLA-C*04:37',
                            'HLA-C*04:38', 'HLA-C*04:39', 'HLA-C*04:40', 'HLA-C*04:41', 'HLA-C*04:42', 'HLA-C*04:43', 'HLA-C*04:44', 'HLA-C*04:45',
                            'HLA-C*04:46', 'HLA-C*04:47', 'HLA-C*04:48', 'HLA-C*04:49', 'HLA-C*04:50', 'HLA-C*04:51', 'HLA-C*04:52', 'HLA-C*04:53',
                            'HLA-C*04:54', 'HLA-C*04:55', 'HLA-C*04:56', 'HLA-C*04:57', 'HLA-C*04:58', 'HLA-C*04:60', 'HLA-C*04:61', 'HLA-C*04:62',
                            'HLA-C*04:63', 'HLA-C*04:64', 'HLA-C*04:65', 'HLA-C*04:66', 'HLA-C*04:67', 'HLA-C*04:68', 'HLA-C*04:69', 'HLA-C*04:70',
                            'HLA-C*05:01', 'HLA-C*05:03', 'HLA-C*05:04', 'HLA-C*05:05', 'HLA-C*05:06', 'HLA-C*05:08', 'HLA-C*05:09', 'HLA-C*05:10',
                            'HLA-C*05:11', 'HLA-C*05:12', 'HLA-C*05:13', 'HLA-C*05:14', 'HLA-C*05:15', 'HLA-C*05:16', 'HLA-C*05:17', 'HLA-C*05:18',
                            'HLA-C*05:19', 'HLA-C*05:20', 'HLA-C*05:21', 'HLA-C*05:22', 'HLA-C*05:23', 'HLA-C*05:24', 'HLA-C*05:25', 'HLA-C*05:26',
                            'HLA-C*05:27', 'HLA-C*05:28', 'HLA-C*05:29', 'HLA-C*05:30', 'HLA-C*05:31', 'HLA-C*05:32', 'HLA-C*05:33', 'HLA-C*05:34',
                            'HLA-C*05:35', 'HLA-C*05:36', 'HLA-C*05:37', 'HLA-C*05:38', 'HLA-C*05:39', 'HLA-C*05:40', 'HLA-C*05:41', 'HLA-C*05:42',
                            'HLA-C*05:43', 'HLA-C*05:44', 'HLA-C*05:45', 'HLA-C*06:02', 'HLA-C*06:03', 'HLA-C*06:04', 'HLA-C*06:05', 'HLA-C*06:06',
                            'HLA-C*06:07', 'HLA-C*06:08', 'HLA-C*06:09', 'HLA-C*06:10', 'HLA-C*06:11', 'HLA-C*06:12', 'HLA-C*06:13', 'HLA-C*06:14',
                            'HLA-C*06:15', 'HLA-C*06:17', 'HLA-C*06:18', 'HLA-C*06:19', 'HLA-C*06:20', 'HLA-C*06:21', 'HLA-C*06:22', 'HLA-C*06:23',
                            'HLA-C*06:24', 'HLA-C*06:25', 'HLA-C*06:26', 'HLA-C*06:27', 'HLA-C*06:28', 'HLA-C*06:29', 'HLA-C*06:30', 'HLA-C*06:31',
                            'HLA-C*06:32', 'HLA-C*06:33', 'HLA-C*06:34', 'HLA-C*06:35', 'HLA-C*06:36', 'HLA-C*06:37', 'HLA-C*06:38', 'HLA-C*06:39',
                            'HLA-C*06:40', 'HLA-C*06:41', 'HLA-C*06:42', 'HLA-C*06:43', 'HLA-C*06:44', 'HLA-C*06:45', 'HLA-C*07:01', 'HLA-C*07:02',
                            'HLA-C*07:03', 'HLA-C*07:04', 'HLA-C*07:05', 'HLA-C*07:06', 'HLA-C*07:07', 'HLA-C*07:08', 'HLA-C*07:09', 'HLA-C*07:10',
                            'HLA-C*07:11', 'HLA-C*07:12', 'HLA-C*07:13', 'HLA-C*07:14', 'HLA-C*07:15', 'HLA-C*07:16', 'HLA-C*07:17', 'HLA-C*07:18',
                            'HLA-C*07:19', 'HLA-C*07:20', 'HLA-C*07:21', 'HLA-C*07:22', 'HLA-C*07:23', 'HLA-C*07:24', 'HLA-C*07:25', 'HLA-C*07:26',
                            'HLA-C*07:27', 'HLA-C*07:28', 'HLA-C*07:29', 'HLA-C*07:30', 'HLA-C*07:31', 'HLA-C*07:35', 'HLA-C*07:36', 'HLA-C*07:37',
                            'HLA-C*07:38', 'HLA-C*07:39', 'HLA-C*07:40', 'HLA-C*07:41', 'HLA-C*07:42', 'HLA-C*07:43', 'HLA-C*07:44', 'HLA-C*07:45',
                            'HLA-C*07:46', 'HLA-C*07:47', 'HLA-C*07:48', 'HLA-C*07:49', 'HLA-C*07:50', 'HLA-C*07:51', 'HLA-C*07:52', 'HLA-C*07:53',
                            'HLA-C*07:54', 'HLA-C*07:56', 'HLA-C*07:57', 'HLA-C*07:58', 'HLA-C*07:59', 'HLA-C*07:60', 'HLA-C*07:62', 'HLA-C*07:63',
                            'HLA-C*07:64', 'HLA-C*07:65', 'HLA-C*07:66', 'HLA-C*07:67', 'HLA-C*07:68', 'HLA-C*07:69', 'HLA-C*07:70', 'HLA-C*07:71',
                            'HLA-C*07:72', 'HLA-C*07:73', 'HLA-C*07:74', 'HLA-C*07:75', 'HLA-C*07:76', 'HLA-C*07:77', 'HLA-C*07:78', 'HLA-C*07:79',
                            'HLA-C*07:80', 'HLA-C*07:81', 'HLA-C*07:82', 'HLA-C*07:83', 'HLA-C*07:84', 'HLA-C*07:85', 'HLA-C*07:86', 'HLA-C*07:87',
                            'HLA-C*07:88', 'HLA-C*07:89', 'HLA-C*07:90', 'HLA-C*07:91', 'HLA-C*07:92', 'HLA-C*07:93', 'HLA-C*07:94', 'HLA-C*07:95',
                            'HLA-C*07:96', 'HLA-C*07:97', 'HLA-C*07:99', 'HLA-C*07:100', 'HLA-C*07:101', 'HLA-C*07:102', 'HLA-C*07:103', 'HLA-C*07:105',
                            'HLA-C*07:106', 'HLA-C*07:107', 'HLA-C*07:108', 'HLA-C*07:109', 'HLA-C*07:110', 'HLA-C*07:111', 'HLA-C*07:112',
                            'HLA-C*07:113', 'HLA-C*07:114', 'HLA-C*07:115', 'HLA-C*07:116', 'HLA-C*07:117', 'HLA-C*07:118', 'HLA-C*07:119',
                            'HLA-C*07:120', 'HLA-C*07:122', 'HLA-C*07:123', 'HLA-C*07:124', 'HLA-C*07:125', 'HLA-C*07:126', 'HLA-C*07:127',
                            'HLA-C*07:128', 'HLA-C*07:129', 'HLA-C*07:130', 'HLA-C*07:131', 'HLA-C*07:132', 'HLA-C*07:133', 'HLA-C*07:134',
                            'HLA-C*07:135', 'HLA-C*07:136', 'HLA-C*07:137', 'HLA-C*07:138', 'HLA-C*07:139', 'HLA-C*07:140', 'HLA-C*07:141',
                            'HLA-C*07:142', 'HLA-C*07:143', 'HLA-C*07:144', 'HLA-C*07:145', 'HLA-C*07:146', 'HLA-C*07:147', 'HLA-C*07:148',
                            'HLA-C*07:149', 'HLA-C*08:01', 'HLA-C*08:02', 'HLA-C*08:03', 'HLA-C*08:04', 'HLA-C*08:05', 'HLA-C*08:06', 'HLA-C*08:07',
                            'HLA-C*08:08', 'HLA-C*08:09', 'HLA-C*08:10', 'HLA-C*08:11', 'HLA-C*08:12', 'HLA-C*08:13', 'HLA-C*08:14', 'HLA-C*08:15',
                            'HLA-C*08:16', 'HLA-C*08:17', 'HLA-C*08:18', 'HLA-C*08:19', 'HLA-C*08:20', 'HLA-C*08:21', 'HLA-C*08:22', 'HLA-C*08:23',
                            'HLA-C*08:24', 'HLA-C*08:25', 'HLA-C*08:27', 'HLA-C*08:28', 'HLA-C*08:29', 'HLA-C*08:30', 'HLA-C*08:31', 'HLA-C*08:32',
                            'HLA-C*08:33', 'HLA-C*08:34', 'HLA-C*08:35', 'HLA-C*12:02', 'HLA-C*12:03', 'HLA-C*12:04', 'HLA-C*12:05', 'HLA-C*12:06',
                            'HLA-C*12:07', 'HLA-C*12:08', 'HLA-C*12:09', 'HLA-C*12:10', 'HLA-C*12:11', 'HLA-C*12:12', 'HLA-C*12:13', 'HLA-C*12:14',
                            'HLA-C*12:15', 'HLA-C*12:16', 'HLA-C*12:17', 'HLA-C*12:18', 'HLA-C*12:19', 'HLA-C*12:20', 'HLA-C*12:21', 'HLA-C*12:22',
                            'HLA-C*12:23', 'HLA-C*12:24', 'HLA-C*12:25', 'HLA-C*12:26', 'HLA-C*12:27', 'HLA-C*12:28', 'HLA-C*12:29', 'HLA-C*12:30',
                            'HLA-C*12:31', 'HLA-C*12:32', 'HLA-C*12:33', 'HLA-C*12:34', 'HLA-C*12:35', 'HLA-C*12:36', 'HLA-C*12:37', 'HLA-C*12:38',
                            'HLA-C*12:40', 'HLA-C*12:41', 'HLA-C*12:43', 'HLA-C*12:44', 'HLA-C*14:02', 'HLA-C*14:03', 'HLA-C*14:04', 'HLA-C*14:05',
                            'HLA-C*14:06', 'HLA-C*14:08', 'HLA-C*14:09', 'HLA-C*14:10', 'HLA-C*14:11', 'HLA-C*14:12', 'HLA-C*14:13', 'HLA-C*14:14',
                            'HLA-C*14:15', 'HLA-C*14:16', 'HLA-C*14:17', 'HLA-C*14:18', 'HLA-C*14:19', 'HLA-C*14:20', 'HLA-C*15:02', 'HLA-C*15:03',
                            'HLA-C*15:04', 'HLA-C*15:05', 'HLA-C*15:06', 'HLA-C*15:07', 'HLA-C*15:08', 'HLA-C*15:09', 'HLA-C*15:10', 'HLA-C*15:11',
                            'HLA-C*15:12', 'HLA-C*15:13', 'HLA-C*15:15', 'HLA-C*15:16', 'HLA-C*15:17', 'HLA-C*15:18', 'HLA-C*15:19', 'HLA-C*15:20',
                            'HLA-C*15:21', 'HLA-C*15:22', 'HLA-C*15:23', 'HLA-C*15:24', 'HLA-C*15:25', 'HLA-C*15:26', 'HLA-C*15:27', 'HLA-C*15:28',
                            'HLA-C*15:29', 'HLA-C*15:30', 'HLA-C*15:31', 'HLA-C*15:33', 'HLA-C*15:34', 'HLA-C*15:35', 'HLA-C*16:01', 'HLA-C*16:02',
                            'HLA-C*16:04', 'HLA-C*16:06', 'HLA-C*16:07', 'HLA-C*16:08', 'HLA-C*16:09', 'HLA-C*16:10', 'HLA-C*16:11', 'HLA-C*16:12',
                            'HLA-C*16:13', 'HLA-C*16:14', 'HLA-C*16:15', 'HLA-C*16:17', 'HLA-C*16:18', 'HLA-C*16:19', 'HLA-C*16:20', 'HLA-C*16:21',
                            'HLA-C*16:22', 'HLA-C*16:23', 'HLA-C*16:24', 'HLA-C*16:25', 'HLA-C*16:26', 'HLA-C*17:01', 'HLA-C*17:02', 'HLA-C*17:03',
                            'HLA-C*17:04', 'HLA-C*17:05', 'HLA-C*17:06', 'HLA-C*17:07', 'HLA-C*18:01', 'HLA-C*18:02', 'HLA-C*18:03', 'HLA-G*01:01',
                            'HLA-G*01:02', 'HLA-G*01:03', 'HLA-G*01:04', 'HLA-G*01:06', 'HLA-G*01:07', 'HLA-G*01:08', 'HLA-G*01:09', 'HLA-E*01:01'])

    @property
    def command(self):
        return self.__command

    @property
    def name(self):
        return self.__name

    @property
    def version(self):
        return self.__version

    @property
    def supportedAlleles(self):
        return self.__alleles

    @property
    def supportedLength(self):
        return self.__length

    def convert_alleles(self, alleles):
        """
        Converts :class:`~Fred2.Core.Allele.Allele` into the internal allele representation of the predictor
        and returns a string representation

        :param  alleles: The :class:`~Fred2.Core.Allele.Allele` for which the
                         internal predictor representation is needed
        :type alleles: list(:class:`~Fred2.Core.Allele.Allele`)
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`
        :rtype: list(str)
        """
        return ["HLA-%s%s:%s" % (a.locus, a.supertype, a.subtype) for a in alleles]

    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        result = defaultdict(dict)
        with open(file, "r") as f:
            f = csv.reader(f, delimiter='\t')
            alleles = [x for x in next(f) if x.strip() != ""]
            print(alleles)
            ic_pos = 4
            offset = 3
            header = next(f)
            print("\t".join(header))
            if "Aff(nM)" in header:
                ic_pos = 9
                offset = 8
            for row in f:
                print("\t".join(row))
                pep_seq = row[1]
                for i, a in enumerate(alleles):
                    result[a][pep_seq] = float(row[ic_pos +i*offset])
        return result

    def get_external_version(self, path=None):
        """
        Returns the external version of the tool by executing
        >{command} --version

        might be dependent on the method and has to be overwritten
        therefore it is declared abstract to enforce the user to
        overwrite the method. The function in the base class can be called
        with super()

        :param str path: Optional specification of executable path if deviant from :attr:`self.__command`
        :return: The external version of the tool or None if tool does not support versioning
        :rtype: str
        """
        # can not be determined netmhcpan does not support --version or similar
        return None

    def prepare_input(self, input, file):
        """
        Prepares input for external tools and writes them to file in the specific format

        NO return value!

        :param: list(str) input: The :class:`~Fred2.Core.Peptide.Peptide` sequences to write into file
        :param File file: File-handler to input file for external tool
        """
        file.write("\n".join(input))



class NetMHCII_2_2(AExternalEpitopePrediction):
    """
    Implements a wrapper for NetMHCII

    .. note::

        Nielsen, M., & Lund, O. (2009). NN-align. An artificial neural network-based alignment algorithm for MHC class
        II peptide binding prediction. BMC Bioinformatics, 10(1), 296.

        Nielsen, M., Lundegaard, C., & Lund, O. (2007). Prediction of MHC class II binding affinity using SMM-align,
        a novel stabilization matrix alignment method. BMC Bioinformatics, 8(1), 238.
    """
    __supported_length = frozenset([15])
    __name = "netmhcII"
    __command = 'netMHCII {peptides} -a {alleles} {options} | grep -v "#" > {out}'
    __alleles = frozenset(
        ['HLA-DRB1*01:01', 'HLA-DRB1*03:01', 'HLA-DRB1*04:01', 'HLA-DRB1*04:04', 'HLA-DRB1*04:05', 'HLA-DRB1*07:01', 'HLA-DRB1*08:02', 'HLA-DRB1*09:01',
         'HLA-DRB1*11:01', 'HLA-DRB1*13:02', 'HLA-DRB1*15:01', 'HLA-DRB3*01:01', 'HLA-DRB4*01:01', 'HLA-DRB5*01:01',
         'H2-IAb', 'H2-IAd'])
    __version = "2.2"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def command(self):
        """
        Defines the commandline call for external tool
        """
        return self.__command

    @property
    def supportedLength(self):
        """
        A list of supported :class:`~Fred2.Core.Peptide.Peptide` lengths
        """
        return self.__supported_length

    @property
    def supportedAlleles(self):
        """A list of valid :class:`~Fred2.Core.Allele.Allele` models"""
        return self.__alleles

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    def _represent(self, allele):
        """
        Internal function transforming an allele object into its representative string
        :param allele: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: str
        """
        if isinstance(allele, MouseAllele):
            return "H-2-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)
        else:
            return "HLA-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)

    def convert_alleles(self, alleles):
        """
        Converts :class:`~Fred2.Core.Allele.Allele` into the internal :class:`~Fred2.Core.Allele.Allele` representation
        of the predictor and returns a string representation

        :param alleles: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`
        :rtype: list(str)
        """
        return [self._represent(a) for a in alleles]

    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        result = defaultdict(defaultdict)
        f = csv.reader(open(file, "r"), delimiter='\t')
        for r in f:
            if not r:
                continue

            row = r[0].split()
            if not len(row):
                continue

            if "HLA-" not in row[0]:
                continue
            result[row[0]][row[2]] = float(row[4])
        return result

    def get_external_version(self, path=None):
        """
        Returns the external version of the tool by executing
        >{command} --version

        might be dependent on the method and has to be overwritten
        therefore it is declared abstract to enforce the user to
        overwrite the method. The function in the base class can be called
        with super()

        :param str path: Optional specification of executable path if deviant from :attr:`self.__command`
        :return: The external version of the tool or None if tool does not support versioning
        :rtype: str
        """
        return None

    def prepare_input(self, input, file):
        """
        Prepares input for external tools
        and writes them to _file in the specific format

        No return value!

        :param: list(str) input: The :class:`~Fred2.Core.Peptide.Peptide` sequences to write into file
        :param File file: File-handler to input file for external tool
        """
        file.write("\n".join(">pepe_%i\n%s" % (i, p) for i, p in enumerate(input)))


class NetMHCIIpan_3_0(AExternalEpitopePrediction):
    """
    Implements a wrapper for NetMHCIIpan.

    .. note::

        Andreatta, M., Karosiene, E., Rasmussen, M., Stryhn, A., Buus, S., & Nielsen, M. (2015).
        Accurate pan-specific prediction of peptide-MHC class II binding affinity with improved binding
        core identification. Immunogenetics, 1-10.
    """

    __supported_length = frozenset([9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
    __name = "netmhcIIpan"
    __command = "netMHCIIpan -f {peptides} -inptype 1 -a {alleles} {options} -xls -xlsfile {out}"
    __alleles = frozenset(
        ['HLA-DRB1*01:01', 'HLA-DRB1*01:02', 'HLA-DRB1*01:03', 'HLA-DRB1*01:04', 'HLA-DRB1*01:05', 'HLA-DRB1*01:06',
         'HLA-DRB1*01:07', 'HLA-DRB1*01:08', 'HLA-DRB1*01:09', 'HLA-DRB1*01:10', 'HLA-DRB1*01:11', 'HLA-DRB1*01:12',
         'HLA-DRB1*01:13', 'HLA-DRB1*01:14', 'HLA-DRB1*01:15', 'HLA-DRB1*01:16', 'HLA-DRB1*01:17', 'HLA-DRB1*01:18',
         'HLA-DRB1*01:19', 'HLA-DRB1*01:20', 'HLA-DRB1*01:21', 'HLA-DRB1*01:22', 'HLA-DRB1*01:23', 'HLA-DRB1*01:24',
         'HLA-DRB1*01:25', 'HLA-DRB1*01:26', 'HLA-DRB1*01:27', 'HLA-DRB1*01:28', 'HLA-DRB1*01:29', 'HLA-DRB1*01:30',
         'HLA-DRB1*01:31', 'HLA-DRB1*01:32', 'HLA-DRB1*03:01', 'HLA-DRB1*03:02', 'HLA-DRB1*03:03', 'HLA-DRB1*03:04',
         'HLA-DRB1*03:05', 'HLA-DRB1*03:06', 'HLA-DRB1*03:07', 'HLA-DRB1*03:08', 'HLA-DRB1*03:10', 'HLA-DRB1*03:11',
         'HLA-DRB1*03:13', 'HLA-DRB1*03:14', 'HLA-DRB1*03:15', 'HLA-DRB1*03:17', 'HLA-DRB1*03:18', 'HLA-DRB1*03:19',
         'HLA-DRB1*03:20', 'HLA-DRB1*03:21', 'HLA-DRB1*03:22', 'HLA-DRB1*03:23', 'HLA-DRB1*03:24', 'HLA-DRB1*03:25',
         'HLA-DRB1*03:26', 'HLA-DRB1*03:27', 'HLA-DRB1*03:28', 'HLA-DRB1*03:29', 'HLA-DRB1*03:30', 'HLA-DRB1*03:31',
         'HLA-DRB1*03:32', 'HLA-DRB1*03:33', 'HLA-DRB1*03:34', 'HLA-DRB1*03:35', 'HLA-DRB1*03:36', 'HLA-DRB1*03:37',
         'HLA-DRB1*03:38', 'HLA-DRB1*03:39', 'HLA-DRB1*03:40', 'HLA-DRB1*03:41', 'HLA-DRB1*03:42', 'HLA-DRB1*03:43',
         'HLA-DRB1*03:44', 'HLA-DRB1*03:45', 'HLA-DRB1*03:46', 'HLA-DRB1*03:47', 'HLA-DRB1*03:48', 'HLA-DRB1*03:49',
         'HLA-DRB1*03:50', 'HLA-DRB1*03:51', 'HLA-DRB1*03:52', 'HLA-DRB1*03:53', 'HLA-DRB1*03:54', 'HLA-DRB1*03:55',
         'HLA-DRB1*04:01', 'HLA-DRB1*04:02', 'HLA-DRB1*04:03', 'HLA-DRB1*04:04', 'HLA-DRB1*04:05', 'HLA-DRB1*04:06',
         'HLA-DRB1*04:07', 'HLA-DRB1*04:08', 'HLA-DRB1*04:09', 'HLA-DRB1*04:10', 'HLA-DRB1*04:11', 'HLA-DRB1*04:12',
         'HLA-DRB1*04:13', 'HLA-DRB1*04:14', 'HLA-DRB1*04:15', 'HLA-DRB1*04:16', 'HLA-DRB1*04:17', 'HLA-DRB1*04:18',
         'HLA-DRB1*04:19', 'HLA-DRB1*04:21', 'HLA-DRB1*04:22', 'HLA-DRB1*04:23', 'HLA-DRB1*04:24', 'HLA-DRB1*04:26',
         'HLA-DRB1*04:27', 'HLA-DRB1*04:28', 'HLA-DRB1*04:29', 'HLA-DRB1*04:30', 'HLA-DRB1*04:31', 'HLA-DRB1*04:33',
         'HLA-DRB1*04:34', 'HLA-DRB1*04:35', 'HLA-DRB1*04:36', 'HLA-DRB1*04:37', 'HLA-DRB1*04:38', 'HLA-DRB1*04:39',
         'HLA-DRB1*04:40', 'HLA-DRB1*04:41', 'HLA-DRB1*04:42', 'HLA-DRB1*04:43', 'HLA-DRB1*04:44', 'HLA-DRB1*04:45',
         'HLA-DRB1*04:46', 'HLA-DRB1*04:47', 'HLA-DRB1*04:48', 'HLA-DRB1*04:49', 'HLA-DRB1*04:50', 'HLA-DRB1*04:51',
         'HLA-DRB1*04:52', 'HLA-DRB1*04:53', 'HLA-DRB1*04:54', 'HLA-DRB1*04:55', 'HLA-DRB1*04:56', 'HLA-DRB1*04:57',
         'HLA-DRB1*04:58', 'HLA-DRB1*04:59', 'HLA-DRB1*04:60', 'HLA-DRB1*04:61', 'HLA-DRB1*04:62', 'HLA-DRB1*04:63',
         'HLA-DRB1*04:64', 'HLA-DRB1*04:65', 'HLA-DRB1*04:66', 'HLA-DRB1*04:67', 'HLA-DRB1*04:68', 'HLA-DRB1*04:69',
         'HLA-DRB1*04:70', 'HLA-DRB1*04:71', 'HLA-DRB1*04:72', 'HLA-DRB1*04:73', 'HLA-DRB1*04:74', 'HLA-DRB1*04:75',
         'HLA-DRB1*04:76', 'HLA-DRB1*04:77', 'HLA-DRB1*04:78', 'HLA-DRB1*04:79', 'HLA-DRB1*04:80', 'HLA-DRB1*04:82',
         'HLA-DRB1*04:83', 'HLA-DRB1*04:84', 'HLA-DRB1*04:85', 'HLA-DRB1*04:86', 'HLA-DRB1*04:87', 'HLA-DRB1*04:88',
         'HLA-DRB1*04:89', 'HLA-DRB1*04:91', 'HLA-DRB1*07:01', 'HLA-DRB1*07:03', 'HLA-DRB1*07:04', 'HLA-DRB1*07:05',
         'HLA-DRB1*07:06', 'HLA-DRB1*07:07', 'HLA-DRB1*07:08', 'HLA-DRB1*07:09', 'HLA-DRB1*07:11', 'HLA-DRB1*07:12',
         'HLA-DRB1*07:13', 'HLA-DRB1*07:14', 'HLA-DRB1*07:15', 'HLA-DRB1*07:16', 'HLA-DRB1*07:17', 'HLA-DRB1*07:19',
         'HLA-DRB1*08:01', 'HLA-DRB1*08:02', 'HLA-DRB1*08:03', 'HLA-DRB1*08:04', 'HLA-DRB1*08:05', 'HLA-DRB1*08:06',
         'HLA-DRB1*08:07', 'HLA-DRB1*08:08', 'HLA-DRB1*08:09', 'HLA-DRB1*08:10', 'HLA-DRB1*08:11', 'HLA-DRB1*08:12',
         'HLA-DRB1*08:13', 'HLA-DRB1*08:14', 'HLA-DRB1*08:15', 'HLA-DRB1*08:16', 'HLA-DRB1*08:18', 'HLA-DRB1*08:19',
         'HLA-DRB1*08:20', 'HLA-DRB1*08:21', 'HLA-DRB1*08:22', 'HLA-DRB1*08:23', 'HLA-DRB1*08:24', 'HLA-DRB1*08:25',
         'HLA-DRB1*08:26', 'HLA-DRB1*08:27', 'HLA-DRB1*08:28', 'HLA-DRB1*08:29', 'HLA-DRB1*08:30', 'HLA-DRB1*08:31',
         'HLA-DRB1*08:32', 'HLA-DRB1*08:33', 'HLA-DRB1*08:34', 'HLA-DRB1*08:35', 'HLA-DRB1*08:36', 'HLA-DRB1*08:37',
         'HLA-DRB1*08:38', 'HLA-DRB1*08:39', 'HLA-DRB1*08:40', 'HLA-DRB1*09:01', 'HLA-DRB1*09:02', 'HLA-DRB1*09:03',
         'HLA-DRB1*09:04', 'HLA-DRB1*09:05', 'HLA-DRB1*09:06', 'HLA-DRB1*09:07', 'HLA-DRB1*09:08', 'HLA-DRB1*09:09',
         'HLA-DRB1*10:01', 'HLA-DRB1*10:02', 'HLA-DRB1*10:03', 'HLA-DRB1*11:01', 'HLA-DRB1*11:02', 'HLA-DRB1*11:03',
         'HLA-DRB1*11:04', 'HLA-DRB1*11:05', 'HLA-DRB1*11:06', 'HLA-DRB1*11:07', 'HLA-DRB1*11:08', 'HLA-DRB1*11:09',
         'HLA-DRB1*11:10', 'HLA-DRB1*11:11', 'HLA-DRB1*11:12', 'HLA-DRB1*11:13', 'HLA-DRB1*11:14', 'HLA-DRB1*11:15',
         'HLA-DRB1*11:16', 'HLA-DRB1*11:17', 'HLA-DRB1*11:18', 'HLA-DRB1*11:19', 'HLA-DRB1*11:20', 'HLA-DRB1*11:21',
         'HLA-DRB1*11:24', 'HLA-DRB1*11:25', 'HLA-DRB1*11:27', 'HLA-DRB1*11:28', 'HLA-DRB1*11:29', 'HLA-DRB1*11:30',
         'HLA-DRB1*11:31', 'HLA-DRB1*11:32', 'HLA-DRB1*11:33', 'HLA-DRB1*11:34', 'HLA-DRB1*11:35', 'HLA-DRB1*11:36',
         'HLA-DRB1*11:37', 'HLA-DRB1*11:38', 'HLA-DRB1*11:39', 'HLA-DRB1*11:41', 'HLA-DRB1*11:42', 'HLA-DRB1*11:43',
         'HLA-DRB1*11:44', 'HLA-DRB1*11:45', 'HLA-DRB1*11:46', 'HLA-DRB1*11:47', 'HLA-DRB1*11:48', 'HLA-DRB1*11:49',
         'HLA-DRB1*11:50', 'HLA-DRB1*11:51', 'HLA-DRB1*11:52', 'HLA-DRB1*11:53', 'HLA-DRB1*11:54', 'HLA-DRB1*11:55',
         'HLA-DRB1*11:56', 'HLA-DRB1*11:57', 'HLA-DRB1*11:58', 'HLA-DRB1*11:59', 'HLA-DRB1*11:60', 'HLA-DRB1*11:61',
         'HLA-DRB1*11:62', 'HLA-DRB1*11:63', 'HLA-DRB1*11:64', 'HLA-DRB1*11:65', 'HLA-DRB1*11:66', 'HLA-DRB1*11:67',
         'HLA-DRB1*11:68', 'HLA-DRB1*11:69', 'HLA-DRB1*11:70', 'HLA-DRB1*11:72', 'HLA-DRB1*11:73', 'HLA-DRB1*11:74',
         'HLA-DRB1*11:75', 'HLA-DRB1*11:76', 'HLA-DRB1*11:77', 'HLA-DRB1*11:78', 'HLA-DRB1*11:79', 'HLA-DRB1*11:80',
         'HLA-DRB1*11:81', 'HLA-DRB1*11:82', 'HLA-DRB1*11:83', 'HLA-DRB1*11:84', 'HLA-DRB1*11:85', 'HLA-DRB1*11:86',
         'HLA-DRB1*11:87', 'HLA-DRB1*11:88', 'HLA-DRB1*11:89', 'HLA-DRB1*11:90', 'HLA-DRB1*11:91', 'HLA-DRB1*11:92',
         'HLA-DRB1*11:93', 'HLA-DRB1*11:94', 'HLA-DRB1*11:95', 'HLA-DRB1*11:96', 'HLA-DRB1*12:01', 'HLA-DRB1*12:02',
         'HLA-DRB1*12:03', 'HLA-DRB1*12:04', 'HLA-DRB1*12:05', 'HLA-DRB1*12:06', 'HLA-DRB1*12:07', 'HLA-DRB1*12:08',
         'HLA-DRB1*12:09', 'HLA-DRB1*12:10', 'HLA-DRB1*12:11', 'HLA-DRB1*12:12', 'HLA-DRB1*12:13', 'HLA-DRB1*12:14',
         'HLA-DRB1*12:15', 'HLA-DRB1*12:16', 'HLA-DRB1*12:17', 'HLA-DRB1*12:18', 'HLA-DRB1*12:19', 'HLA-DRB1*12:20',
         'HLA-DRB1*12:21', 'HLA-DRB1*12:22', 'HLA-DRB1*12:23', 'HLA-DRB1*13:01', 'HLA-DRB1*13:02', 'HLA-DRB1*13:03',
         'HLA-DRB1*13:04', 'HLA-DRB1*13:05', 'HLA-DRB1*13:06', 'HLA-DRB1*13:07', 'HLA-DRB1*13:08', 'HLA-DRB1*13:09',
         'HLA-DRB1*13:10', 'HLA-DRB1*13:100', 'HLA-DRB1*13:101', 'HLA-DRB1*13:11', 'HLA-DRB1*13:12', 'HLA-DRB1*13:13',
         'HLA-DRB1*13:14', 'HLA-DRB1*13:15', 'HLA-DRB1*13:16', 'HLA-DRB1*13:17', 'HLA-DRB1*13:18', 'HLA-DRB1*13:19',
         'HLA-DRB1*13:20', 'HLA-DRB1*13:21', 'HLA-DRB1*13:22', 'HLA-DRB1*13:23', 'HLA-DRB1*13:24', 'HLA-DRB1*13:26',
         'HLA-DRB1*13:27', 'HLA-DRB1*13:29', 'HLA-DRB1*13:30', 'HLA-DRB1*13:31', 'HLA-DRB1*13:32', 'HLA-DRB1*13:33',
         'HLA-DRB1*13:34', 'HLA-DRB1*13:35', 'HLA-DRB1*13:36', 'HLA-DRB1*13:37', 'HLA-DRB1*13:38', 'HLA-DRB1*13:39',
         'HLA-DRB1*13:41', 'HLA-DRB1*13:42', 'HLA-DRB1*13:43', 'HLA-DRB1*13:44', 'HLA-DRB1*13:46', 'HLA-DRB1*13:47',
         'HLA-DRB1*13:48', 'HLA-DRB1*13:49', 'HLA-DRB1*13:50', 'HLA-DRB1*13:51', 'HLA-DRB1*13:52', 'HLA-DRB1*13:53',
         'HLA-DRB1*13:54', 'HLA-DRB1*13:55', 'HLA-DRB1*13:56', 'HLA-DRB1*13:57', 'HLA-DRB1*13:58', 'HLA-DRB1*13:59',
         'HLA-DRB1*13:60', 'HLA-DRB1*13:61', 'HLA-DRB1*13:62', 'HLA-DRB1*13:63', 'HLA-DRB1*13:64', 'HLA-DRB1*13:65',
         'HLA-DRB1*13:66', 'HLA-DRB1*13:67', 'HLA-DRB1*13:68', 'HLA-DRB1*13:69', 'HLA-DRB1*13:70', 'HLA-DRB1*13:71',
         'HLA-DRB1*13:72', 'HLA-DRB1*13:73', 'HLA-DRB1*13:74', 'HLA-DRB1*13:75', 'HLA-DRB1*13:76', 'HLA-DRB1*13:77',
         'HLA-DRB1*13:78', 'HLA-DRB1*13:79', 'HLA-DRB1*13:80', 'HLA-DRB1*13:81', 'HLA-DRB1*13:82', 'HLA-DRB1*13:83',
         'HLA-DRB1*13:84', 'HLA-DRB1*13:85', 'HLA-DRB1*13:86', 'HLA-DRB1*13:87', 'HLA-DRB1*13:88', 'HLA-DRB1*13:89',
         'HLA-DRB1*13:90', 'HLA-DRB1*13:91', 'HLA-DRB1*13:92', 'HLA-DRB1*13:93', 'HLA-DRB1*13:94', 'HLA-DRB1*13:95',
         'HLA-DRB1*13:96', 'HLA-DRB1*13:97', 'HLA-DRB1*13:98', 'HLA-DRB1*13:99', 'HLA-DRB1*14:01', 'HLA-DRB1*14:02',
         'HLA-DRB1*14:03', 'HLA-DRB1*14:04', 'HLA-DRB1*14:05', 'HLA-DRB1*14:06', 'HLA-DRB1*14:07', 'HLA-DRB1*14:08',
         'HLA-DRB1*14:09', 'HLA-DRB1*14:10', 'HLA-DRB1*14:11', 'HLA-DRB1*14:12', 'HLA-DRB1*14:13', 'HLA-DRB1*14:14',
         'HLA-DRB1*14:15', 'HLA-DRB1*14:16', 'HLA-DRB1*14:17', 'HLA-DRB1*14:18', 'HLA-DRB1*14:19', 'HLA-DRB1*14:20',
         'HLA-DRB1*14:21', 'HLA-DRB1*14:22', 'HLA-DRB1*14:23', 'HLA-DRB1*14:24', 'HLA-DRB1*14:25', 'HLA-DRB1*14:26',
         'HLA-DRB1*14:27', 'HLA-DRB1*14:28', 'HLA-DRB1*14:29', 'HLA-DRB1*14:30', 'HLA-DRB1*14:31', 'HLA-DRB1*14:32',
         'HLA-DRB1*14:33', 'HLA-DRB1*14:34', 'HLA-DRB1*14:35', 'HLA-DRB1*14:36', 'HLA-DRB1*14:37', 'HLA-DRB1*14:38',
         'HLA-DRB1*14:39', 'HLA-DRB1*14:40', 'HLA-DRB1*14:41', 'HLA-DRB1*14:42', 'HLA-DRB1*14:43', 'HLA-DRB1*14:44',
         'HLA-DRB1*14:45', 'HLA-DRB1*14:46', 'HLA-DRB1*14:47', 'HLA-DRB1*14:48', 'HLA-DRB1*14:49', 'HLA-DRB1*14:50',
         'HLA-DRB1*14:51', 'HLA-DRB1*14:52', 'HLA-DRB1*14:53', 'HLA-DRB1*14:54', 'HLA-DRB1*14:55', 'HLA-DRB1*14:56',
         'HLA-DRB1*14:57', 'HLA-DRB1*14:58', 'HLA-DRB1*14:59', 'HLA-DRB1*14:60', 'HLA-DRB1*14:61', 'HLA-DRB1*14:62',
         'HLA-DRB1*14:63', 'HLA-DRB1*14:64', 'HLA-DRB1*14:65', 'HLA-DRB1*14:67', 'HLA-DRB1*14:68', 'HLA-DRB1*14:69',
         'HLA-DRB1*14:70', 'HLA-DRB1*14:71', 'HLA-DRB1*14:72', 'HLA-DRB1*14:73', 'HLA-DRB1*14:74', 'HLA-DRB1*14:75',
         'HLA-DRB1*14:76', 'HLA-DRB1*14:77', 'HLA-DRB1*14:78', 'HLA-DRB1*14:79', 'HLA-DRB1*14:80', 'HLA-DRB1*14:81',
         'HLA-DRB1*14:82', 'HLA-DRB1*14:83', 'HLA-DRB1*14:84', 'HLA-DRB1*14:85', 'HLA-DRB1*14:86', 'HLA-DRB1*14:87',
         'HLA-DRB1*14:88', 'HLA-DRB1*14:89', 'HLA-DRB1*14:90', 'HLA-DRB1*14:91', 'HLA-DRB1*14:93', 'HLA-DRB1*14:94',
         'HLA-DRB1*14:95', 'HLA-DRB1*14:96', 'HLA-DRB1*14:97', 'HLA-DRB1*14:98', 'HLA-DRB1*14:99', 'HLA-DRB1*15:01',
         'HLA-DRB1*15:02', 'HLA-DRB1*15:03', 'HLA-DRB1*15:04', 'HLA-DRB1*15:05', 'HLA-DRB1*15:06', 'HLA-DRB1*15:07',
         'HLA-DRB1*15:08', 'HLA-DRB1*15:09', 'HLA-DRB1*15:10', 'HLA-DRB1*15:11', 'HLA-DRB1*15:12', 'HLA-DRB1*15:13',
         'HLA-DRB1*15:14', 'HLA-DRB1*15:15', 'HLA-DRB1*15:16', 'HLA-DRB1*15:18', 'HLA-DRB1*15:19', 'HLA-DRB1*15:20',
         'HLA-DRB1*15:21', 'HLA-DRB1*15:22', 'HLA-DRB1*15:23', 'HLA-DRB1*15:24', 'HLA-DRB1*15:25', 'HLA-DRB1*15:26',
         'HLA-DRB1*15:27', 'HLA-DRB1*15:28', 'HLA-DRB1*15:29', 'HLA-DRB1*15:30', 'HLA-DRB1*15:31', 'HLA-DRB1*15:32',
         'HLA-DRB1*15:33', 'HLA-DRB1*15:34', 'HLA-DRB1*15:35', 'HLA-DRB1*15:36', 'HLA-DRB1*15:37', 'HLA-DRB1*15:38',
         'HLA-DRB1*15:39', 'HLA-DRB1*15:40', 'HLA-DRB1*15:41', 'HLA-DRB1*15:42', 'HLA-DRB1*15:43', 'HLA-DRB1*15:44',
         'HLA-DRB1*15:45', 'HLA-DRB1*15:46', 'HLA-DRB1*15:47', 'HLA-DRB1*15:48', 'HLA-DRB1*15:49', 'HLA-DRB1*16:01',
         'HLA-DRB1*16:02', 'HLA-DRB1*16:03', 'HLA-DRB1*16:04', 'HLA-DRB1*16:05', 'HLA-DRB1*16:07', 'HLA-DRB1*16:08',
         'HLA-DRB1*16:09', 'HLA-DRB1*16:10', 'HLA-DRB1*16:11', 'HLA-DRB1*16:12', 'HLA-DRB1*16:14', 'HLA-DRB1*16:15',
         'HLA-DRB1*16:16', 'HLA-DRB3*01:01', 'HLA-DRB3*01:04', 'HLA-DRB3*01:05', 'HLA-DRB3*01:08', 'HLA-DRB3*01:09',
         'HLA-DRB3*01:11', 'HLA-DRB3*01:12', 'HLA-DRB3*01:13', 'HLA-DRB3*01:14', 'HLA-DRB3*02:01', 'HLA-DRB3*02:02',
         'HLA-DRB3*02:04', 'HLA-DRB3*02:05', 'HLA-DRB3*02:09', 'HLA-DRB3*02:10', 'HLA-DRB3*02:11', 'HLA-DRB3*02:12',
         'HLA-DRB3*02:13', 'HLA-DRB3*02:14', 'HLA-DRB3*02:15', 'HLA-DRB3*02:16', 'HLA-DRB3*02:17', 'HLA-DRB3*02:18',
         'HLA-DRB3*02:19', 'HLA-DRB3*02:20', 'HLA-DRB3*02:21', 'HLA-DRB3*02:22', 'HLA-DRB3*02:23', 'HLA-DRB3*02:24',
         'HLA-DRB3*02:25', 'HLA-DRB3*03:01', 'HLA-DRB3*03:03', 'HLA-DRB4*01:01', 'HLA-DRB4*01:03', 'HLA-DRB4*01:04',
         'HLA-DRB4*01:06', 'HLA-DRB4*01:07', 'HLA-DRB4*01:08', 'HLA-DRB5*01:01', 'HLA-DRB5*01:02', 'HLA-DRB5*01:03',
         'HLA-DRB5*01:04', 'HLA-DRB5*01:05', 'HLA-DRB5*01:06', 'HLA-DRB5*01:08N', 'HLA-DRB5*01:11', 'HLA-DRB5*01:12',
         'HLA-DRB5*01:13', 'HLA-DRB5*01:14', 'HLA-DRB5*02:02', 'HLA-DRB5*02:03', 'HLA-DRB5*02:04', 'HLA-DRB5*02:05',
         'HLA-DPA1*01:03-HLA-DPB1*01:01', 'HLA-DPA1*01:03-HLA-DPB1*02:01', 'HLA-DPA1*01:03-HLA-DPB1*02:02', 'HLA-DPA1*01:03-HLA-DPB1*03:01',
         'HLA-DPA1*01:03-HLA-DPB1*04:01', 'HLA-DPA1*01:03-HLA-DPB1*04:02', 'HLA-DPA1*01:03-HLA-DPB1*05:01', 'HLA-DPA1*01:03-HLA-DPB1*06:01',
         'HLA-DPA1*01:03-HLA-DPB1*08:01', 'HLA-DPA1*01:03-HLA-DPB1*09:01', 'HLA-DPA1*01:03-HLA-DPB1*10:001', 'HLA-DPA1*01:03-HLA-DPB1*10:01',
         'HLA-DPA1*01:03-HLA-DPB1*10:101', 'HLA-DPA1*01:03-HLA-DPB1*10:201',
         'HLA-DPA1*01:03-HLA-DPB1*10:301', 'HLA-DPA1*01:03-HLA-DPB1*10:401',
         'HLA-DPA1*01:03-HLA-DPB1*10:501', 'HLA-DPA1*01:03-HLA-DPB1*10:601', 'HLA-DPA1*01:03-HLA-DPB1*10:701', 'HLA-DPA1*01:03-HLA-DPB1*10:801',
         'HLA-DPA1*01:03-HLA-DPB1*10:901', 'HLA-DPA1*01:03-HLA-DPB1*11:001',
         'HLA-DPA1*01:03-HLA-DPB1*11:01', 'HLA-DPA1*01:03-HLA-DPB1*11:101', 'HLA-DPA1*01:03-HLA-DPB1*11:201', 'HLA-DPA1*01:03-HLA-DPB1*11:301',
         'HLA-DPA1*01:03-HLA-DPB1*11:401', 'HLA-DPA1*01:03-HLA-DPB1*11:501',
         'HLA-DPA1*01:03-HLA-DPB1*11:601', 'HLA-DPA1*01:03-HLA-DPB1*11:701', 'HLA-DPA1*01:03-HLA-DPB1*11:801', 'HLA-DPA1*01:03-HLA-DPB1*11:901',
         'HLA-DPA1*01:03-HLA-DPB1*12:101', 'HLA-DPA1*01:03-HLA-DPB1*12:201',
         'HLA-DPA1*01:03-HLA-DPB1*12:301', 'HLA-DPA1*01:03-HLA-DPB1*12:401', 'HLA-DPA1*01:03-HLA-DPB1*12:501', 'HLA-DPA1*01:03-HLA-DPB1*12:601',
         'HLA-DPA1*01:03-HLA-DPB1*12:701', 'HLA-DPA1*01:03-HLA-DPB1*12:801',
         'HLA-DPA1*01:03-HLA-DPB1*12:901', 'HLA-DPA1*01:03-HLA-DPB1*13:001', 'HLA-DPA1*01:03-HLA-DPB1*13:01', 'HLA-DPA1*01:03-HLA-DPB1*13:101',
         'HLA-DPA1*01:03-HLA-DPB1*13:201', 'HLA-DPA1*01:03-HLA-DPB1*13:301',
         'HLA-DPA1*01:03-HLA-DPB1*13:401', 'HLA-DPA1*01:03-HLA-DPB1*14:01', 'HLA-DPA1*01:03-HLA-DPB1*15:01', 'HLA-DPA1*01:03-HLA-DPB1*16:01',
         'HLA-DPA1*01:03-HLA-DPB1*17:01', 'HLA-DPA1*01:03-HLA-DPB1*18:01',
         'HLA-DPA1*01:03-HLA-DPB1*19:01', 'HLA-DPA1*01:03-HLA-DPB1*20:01', 'HLA-DPA1*01:03-HLA-DPB1*21:01', 'HLA-DPA1*01:03-HLA-DPB1*22:01',
         'HLA-DPA1*01:03-HLA-DPB1*23:01', 'HLA-DPA1*01:03-HLA-DPB1*24:01',
         'HLA-DPA1*01:03-HLA-DPB1*25:01', 'HLA-DPA1*01:03-HLA-DPB1*26:01', 'HLA-DPA1*01:03-HLA-DPB1*27:01', 'HLA-DPA1*01:03-HLA-DPB1*28:01',
         'HLA-DPA1*01:03-HLA-DPB1*29:01', 'HLA-DPA1*01:03-HLA-DPB1*30:01',
         'HLA-DPA1*01:03-HLA-DPB1*31:01', 'HLA-DPA1*01:03-HLA-DPB1*32:01', 'HLA-DPA1*01:03-HLA-DPB1*33:01', 'HLA-DPA1*01:03-HLA-DPB1*34:01',
         'HLA-DPA1*01:03-HLA-DPB1*35:01', 'HLA-DPA1*01:03-HLA-DPB1*36:01',
         'HLA-DPA1*01:03-HLA-DPB1*37:01', 'HLA-DPA1*01:03-HLA-DPB1*38:01', 'HLA-DPA1*01:03-HLA-DPB1*39:01', 'HLA-DPA1*01:03-HLA-DPB1*40:01',
         'HLA-DPA1*01:03-HLA-DPB1*41:01', 'HLA-DPA1*01:03-HLA-DPB1*44:01',
         'HLA-DPA1*01:03-HLA-DPB1*45:01', 'HLA-DPA1*01:03-HLA-DPB1*46:01', 'HLA-DPA1*01:03-HLA-DPB1*47:01', 'HLA-DPA1*01:03-HLA-DPB1*48:01',
         'HLA-DPA1*01:03-HLA-DPB1*49:01', 'HLA-DPA1*01:03-HLA-DPB1*50:01',
         'HLA-DPA1*01:03-HLA-DPB1*51:01', 'HLA-DPA1*01:03-HLA-DPB1*52:01', 'HLA-DPA1*01:03-HLA-DPB1*53:01', 'HLA-DPA1*01:03-HLA-DPB1*54:01',
         'HLA-DPA1*01:03-HLA-DPB1*55:01', 'HLA-DPA1*01:03-HLA-DPB1*56:01',
         'HLA-DPA1*01:03-HLA-DPB1*58:01', 'HLA-DPA1*01:03-HLA-DPB1*59:01', 'HLA-DPA1*01:03-HLA-DPB1*60:01', 'HLA-DPA1*01:03-HLA-DPB1*62:01',
         'HLA-DPA1*01:03-HLA-DPB1*63:01', 'HLA-DPA1*01:03-HLA-DPB1*65:01',
         'HLA-DPA1*01:03-HLA-DPB1*66:01', 'HLA-DPA1*01:03-HLA-DPB1*67:01', 'HLA-DPA1*01:03-HLA-DPB1*68:01', 'HLA-DPA1*01:03-HLA-DPB1*69:01',
         'HLA-DPA1*01:03-HLA-DPB1*70:01', 'HLA-DPA1*01:03-HLA-DPB1*71:01',
         'HLA-DPA1*01:03-HLA-DPB1*72:01', 'HLA-DPA1*01:03-HLA-DPB1*73:01', 'HLA-DPA1*01:03-HLA-DPB1*74:01', 'HLA-DPA1*01:03-HLA-DPB1*75:01',
         'HLA-DPA1*01:03-HLA-DPB1*76:01', 'HLA-DPA1*01:03-HLA-DPB1*77:01',
         'HLA-DPA1*01:03-HLA-DPB1*78:01', 'HLA-DPA1*01:03-HLA-DPB1*79:01', 'HLA-DPA1*01:03-HLA-DPB1*80:01', 'HLA-DPA1*01:03-HLA-DPB1*81:01',
         'HLA-DPA1*01:03-HLA-DPB1*82:01', 'HLA-DPA1*01:03-HLA-DPB1*83:01',
         'HLA-DPA1*01:03-HLA-DPB1*84:01', 'HLA-DPA1*01:03-HLA-DPB1*85:01', 'HLA-DPA1*01:03-HLA-DPB1*86:01', 'HLA-DPA1*01:03-HLA-DPB1*87:01',
         'HLA-DPA1*01:03-HLA-DPB1*88:01', 'HLA-DPA1*01:03-HLA-DPB1*89:01',
         'HLA-DPA1*01:03-HLA-DPB1*90:01', 'HLA-DPA1*01:03-HLA-DPB1*91:01', 'HLA-DPA1*01:03-HLA-DPB1*92:01', 'HLA-DPA1*01:03-HLA-DPB1*93:01',
         'HLA-DPA1*01:03-HLA-DPB1*94:01', 'HLA-DPA1*01:03-HLA-DPB1*95:01',
         'HLA-DPA1*01:03-HLA-DPB1*96:01', 'HLA-DPA1*01:03-HLA-DPB1*97:01', 'HLA-DPA1*01:03-HLA-DPB1*98:01', 'HLA-DPA1*01:03-HLA-DPB1*99:01',
         'HLA-DPA1*01:04-HLA-DPB1*01:01', 'HLA-DPA1*01:04-HLA-DPB1*02:01',
         'HLA-DPA1*01:04-HLA-DPB1*02:02', 'HLA-DPA1*01:04-HLA-DPB1*03:01', 'HLA-DPA1*01:04-HLA-DPB1*04:01', 'HLA-DPA1*01:04-HLA-DPB1*04:02',
         'HLA-DPA1*01:04-HLA-DPB1*05:01', 'HLA-DPA1*01:04-HLA-DPB1*06:01',
         'HLA-DPA1*01:04-HLA-DPB1*08:01', 'HLA-DPA1*01:04-HLA-DPB1*09:01', 'HLA-DPA1*01:04-HLA-DPB1*10:001', 'HLA-DPA1*01:04-HLA-DPB1*10:01',
         'HLA-DPA1*01:04-HLA-DPB1*10:101', 'HLA-DPA1*01:04-HLA-DPB1*10:201',
         'HLA-DPA1*01:04-HLA-DPB1*10:301', 'HLA-DPA1*01:04-HLA-DPB1*10:401', 'HLA-DPA1*01:04-HLA-DPB1*10:501', 'HLA-DPA1*01:04-HLA-DPB1*10:601',
         'HLA-DPA1*01:04-HLA-DPB1*10:701', 'HLA-DPA1*01:04-HLA-DPB1*10:801',
         'HLA-DPA1*01:04-HLA-DPB1*10:901', 'HLA-DPA1*01:04-HLA-DPB1*11:001', 'HLA-DPA1*01:04-HLA-DPB1*11:01', 'HLA-DPA1*01:04-HLA-DPB1*11:101',
         'HLA-DPA1*01:04-HLA-DPB1*11:201', 'HLA-DPA1*01:04-HLA-DPB1*11:301',
         'HLA-DPA1*01:04-HLA-DPB1*11:401', 'HLA-DPA1*01:04-HLA-DPB1*11:501', 'HLA-DPA1*01:04-HLA-DPB1*11:601', 'HLA-DPA1*01:04-HLA-DPB1*11:701',
         'HLA-DPA1*01:04-HLA-DPB1*11:801', 'HLA-DPA1*01:04-HLA-DPB1*11:901',
         'HLA-DPA1*01:04-HLA-DPB1*12:101', 'HLA-DPA1*01:04-HLA-DPB1*12:201', 'HLA-DPA1*01:04-HLA-DPB1*12:301', 'HLA-DPA1*01:04-HLA-DPB1*12:401',
         'HLA-DPA1*01:04-HLA-DPB1*12:501', 'HLA-DPA1*01:04-HLA-DPB1*12:601',
         'HLA-DPA1*01:04-HLA-DPB1*12:701', 'HLA-DPA1*01:04-HLA-DPB1*12:801', 'HLA-DPA1*01:04-HLA-DPB1*12:901', 'HLA-DPA1*01:04-HLA-DPB1*13:001',
         'HLA-DPA1*01:04-HLA-DPB1*13:01', 'HLA-DPA1*01:04-HLA-DPB1*13:101',
         'HLA-DPA1*01:04-HLA-DPB1*13:201', 'HLA-DPA1*01:04-HLA-DPB1*13:301', 'HLA-DPA1*01:04-HLA-DPB1*13:401', 'HLA-DPA1*01:04-HLA-DPB1*14:01',
         'HLA-DPA1*01:04-HLA-DPB1*15:01', 'HLA-DPA1*01:04-HLA-DPB1*16:01',
         'HLA-DPA1*01:04-HLA-DPB1*17:01', 'HLA-DPA1*01:04-HLA-DPB1*18:01', 'HLA-DPA1*01:04-HLA-DPB1*19:01', 'HLA-DPA1*01:04-HLA-DPB1*20:01',
         'HLA-DPA1*01:04-HLA-DPB1*21:01', 'HLA-DPA1*01:04-HLA-DPB1*22:01',
         'HLA-DPA1*01:04-HLA-DPB1*23:01', 'HLA-DPA1*01:04-HLA-DPB1*24:01', 'HLA-DPA1*01:04-HLA-DPB1*25:01', 'HLA-DPA1*01:04-HLA-DPB1*26:01',
         'HLA-DPA1*01:04-HLA-DPB1*27:01', 'HLA-DPA1*01:04-HLA-DPB1*28:01',
         'HLA-DPA1*01:04-HLA-DPB1*29:01', 'HLA-DPA1*01:04-HLA-DPB1*30:01', 'HLA-DPA1*01:04-HLA-DPB1*31:01', 'HLA-DPA1*01:04-HLA-DPB1*32:01',
         'HLA-DPA1*01:04-HLA-DPB1*33:01', 'HLA-DPA1*01:04-HLA-DPB1*34:01',
         'HLA-DPA1*01:04-HLA-DPB1*35:01', 'HLA-DPA1*01:04-HLA-DPB1*36:01', 'HLA-DPA1*01:04-HLA-DPB1*37:01', 'HLA-DPA1*01:04-HLA-DPB1*38:01',
         'HLA-DPA1*01:04-HLA-DPB1*39:01', 'HLA-DPA1*01:04-HLA-DPB1*40:01',
         'HLA-DPA1*01:04-HLA-DPB1*41:01', 'HLA-DPA1*01:04-HLA-DPB1*44:01', 'HLA-DPA1*01:04-HLA-DPB1*45:01', 'HLA-DPA1*01:04-HLA-DPB1*46:01',
         'HLA-DPA1*01:04-HLA-DPB1*47:01', 'HLA-DPA1*01:04-HLA-DPB1*48:01',
         'HLA-DPA1*01:04-HLA-DPB1*49:01', 'HLA-DPA1*01:04-HLA-DPB1*50:01', 'HLA-DPA1*01:04-HLA-DPB1*51:01', 'HLA-DPA1*01:04-HLA-DPB1*52:01',
         'HLA-DPA1*01:04-HLA-DPB1*53:01', 'HLA-DPA1*01:04-HLA-DPB1*54:01',
         'HLA-DPA1*01:04-HLA-DPB1*55:01', 'HLA-DPA1*01:04-HLA-DPB1*56:01', 'HLA-DPA1*01:04-HLA-DPB1*58:01', 'HLA-DPA1*01:04-HLA-DPB1*59:01',
         'HLA-DPA1*01:04-HLA-DPB1*60:01', 'HLA-DPA1*01:04-HLA-DPB1*62:01',
         'HLA-DPA1*01:04-HLA-DPB1*63:01', 'HLA-DPA1*01:04-HLA-DPB1*65:01', 'HLA-DPA1*01:04-HLA-DPB1*66:01', 'HLA-DPA1*01:04-HLA-DPB1*67:01',
         'HLA-DPA1*01:04-HLA-DPB1*68:01', 'HLA-DPA1*01:04-HLA-DPB1*69:01',
         'HLA-DPA1*01:04-HLA-DPB1*70:01', 'HLA-DPA1*01:04-HLA-DPB1*71:01', 'HLA-DPA1*01:04-HLA-DPB1*72:01', 'HLA-DPA1*01:04-HLA-DPB1*73:01',
         'HLA-DPA1*01:04-HLA-DPB1*74:01', 'HLA-DPA1*01:04-HLA-DPB1*75:01',
         'HLA-DPA1*01:04-HLA-DPB1*76:01', 'HLA-DPA1*01:04-HLA-DPB1*77:01', 'HLA-DPA1*01:04-HLA-DPB1*78:01', 'HLA-DPA1*01:04-HLA-DPB1*79:01',
         'HLA-DPA1*01:04-HLA-DPB1*80:01', 'HLA-DPA1*01:04-HLA-DPB1*81:01',
         'HLA-DPA1*01:04-HLA-DPB1*82:01', 'HLA-DPA1*01:04-HLA-DPB1*83:01', 'HLA-DPA1*01:04-HLA-DPB1*84:01', 'HLA-DPA1*01:04-HLA-DPB1*85:01',
         'HLA-DPA1*01:04-HLA-DPB1*86:01', 'HLA-DPA1*01:04-HLA-DPB1*87:01',
         'HLA-DPA1*01:04-HLA-DPB1*88:01', 'HLA-DPA1*01:04-HLA-DPB1*89:01', 'HLA-DPA1*01:04-HLA-DPB1*90:01', 'HLA-DPA1*01:04-HLA-DPB1*91:01',
         'HLA-DPA1*01:04-HLA-DPB1*92:01', 'HLA-DPA1*01:04-HLA-DPB1*93:01',
         'HLA-DPA1*01:04-HLA-DPB1*94:01', 'HLA-DPA1*01:04-HLA-DPB1*95:01', 'HLA-DPA1*01:04-HLA-DPB1*96:01', 'HLA-DPA1*01:04-HLA-DPB1*97:01',
         'HLA-DPA1*01:04-HLA-DPB1*98:01', 'HLA-DPA1*01:04-HLA-DPB1*99:01',
         'HLA-DPA1*01:05-HLA-DPB1*01:01', 'HLA-DPA1*01:05-HLA-DPB1*02:01', 'HLA-DPA1*01:05-HLA-DPB1*02:02', 'HLA-DPA1*01:05-HLA-DPB1*03:01',
         'HLA-DPA1*01:05-HLA-DPB1*04:01', 'HLA-DPA1*01:05-HLA-DPB1*04:02',
         'HLA-DPA1*01:05-HLA-DPB1*05:01', 'HLA-DPA1*01:05-HLA-DPB1*06:01', 'HLA-DPA1*01:05-HLA-DPB1*08:01', 'HLA-DPA1*01:05-HLA-DPB1*09:01',
         'HLA-DPA1*01:05-HLA-DPB1*10:001', 'HLA-DPA1*01:05-HLA-DPB1*10:01',
         'HLA-DPA1*01:05-HLA-DPB1*10:101', 'HLA-DPA1*01:05-HLA-DPB1*10:201', 'HLA-DPA1*01:05-HLA-DPB1*10:301', 'HLA-DPA1*01:05-HLA-DPB1*10:401',
         'HLA-DPA1*01:05-HLA-DPB1*10:501', 'HLA-DPA1*01:05-HLA-DPB1*10:601',
         'HLA-DPA1*01:05-HLA-DPB1*10:701', 'HLA-DPA1*01:05-HLA-DPB1*10:801', 'HLA-DPA1*01:05-HLA-DPB1*10:901', 'HLA-DPA1*01:05-HLA-DPB1*11:001',
         'HLA-DPA1*01:05-HLA-DPB1*11:01', 'HLA-DPA1*01:05-HLA-DPB1*11:101',
         'HLA-DPA1*01:05-HLA-DPB1*11:201', 'HLA-DPA1*01:05-HLA-DPB1*11:301', 'HLA-DPA1*01:05-HLA-DPB1*11:401', 'HLA-DPA1*01:05-HLA-DPB1*11:501',
         'HLA-DPA1*01:05-HLA-DPB1*11:601', 'HLA-DPA1*01:05-HLA-DPB1*11:701',
         'HLA-DPA1*01:05-HLA-DPB1*11:801', 'HLA-DPA1*01:05-HLA-DPB1*11:901', 'HLA-DPA1*01:05-HLA-DPB1*12:101', 'HLA-DPA1*01:05-HLA-DPB1*12:201',
         'HLA-DPA1*01:05-HLA-DPB1*12:301', 'HLA-DPA1*01:05-HLA-DPB1*12:401',
         'HLA-DPA1*01:05-HLA-DPB1*12:501', 'HLA-DPA1*01:05-HLA-DPB1*12:601', 'HLA-DPA1*01:05-HLA-DPB1*12:701', 'HLA-DPA1*01:05-HLA-DPB1*12:801',
         'HLA-DPA1*01:05-HLA-DPB1*12:901', 'HLA-DPA1*01:05-HLA-DPB1*13:001',
         'HLA-DPA1*01:05-HLA-DPB1*13:01', 'HLA-DPA1*01:05-HLA-DPB1*13:101', 'HLA-DPA1*01:05-HLA-DPB1*13:201', 'HLA-DPA1*01:05-HLA-DPB1*13:301',
         'HLA-DPA1*01:05-HLA-DPB1*13:401', 'HLA-DPA1*01:05-HLA-DPB1*14:01',
         'HLA-DPA1*01:05-HLA-DPB1*15:01', 'HLA-DPA1*01:05-HLA-DPB1*16:01', 'HLA-DPA1*01:05-HLA-DPB1*17:01', 'HLA-DPA1*01:05-HLA-DPB1*18:01',
         'HLA-DPA1*01:05-HLA-DPB1*19:01', 'HLA-DPA1*01:05-HLA-DPB1*20:01',
         'HLA-DPA1*01:05-HLA-DPB1*21:01', 'HLA-DPA1*01:05-HLA-DPB1*22:01', 'HLA-DPA1*01:05-HLA-DPB1*23:01', 'HLA-DPA1*01:05-HLA-DPB1*24:01',
         'HLA-DPA1*01:05-HLA-DPB1*25:01', 'HLA-DPA1*01:05-HLA-DPB1*26:01',
         'HLA-DPA1*01:05-HLA-DPB1*27:01', 'HLA-DPA1*01:05-HLA-DPB1*28:01', 'HLA-DPA1*01:05-HLA-DPB1*29:01', 'HLA-DPA1*01:05-HLA-DPB1*30:01',
         'HLA-DPA1*01:05-HLA-DPB1*31:01', 'HLA-DPA1*01:05-HLA-DPB1*32:01',
         'HLA-DPA1*01:05-HLA-DPB1*33:01', 'HLA-DPA1*01:05-HLA-DPB1*34:01', 'HLA-DPA1*01:05-HLA-DPB1*35:01', 'HLA-DPA1*01:05-HLA-DPB1*36:01',
         'HLA-DPA1*01:05-HLA-DPB1*37:01', 'HLA-DPA1*01:05-HLA-DPB1*38:01',
         'HLA-DPA1*01:05-HLA-DPB1*39:01', 'HLA-DPA1*01:05-HLA-DPB1*40:01', 'HLA-DPA1*01:05-HLA-DPB1*41:01', 'HLA-DPA1*01:05-HLA-DPB1*44:01',
         'HLA-DPA1*01:05-HLA-DPB1*45:01', 'HLA-DPA1*01:05-HLA-DPB1*46:01',
         'HLA-DPA1*01:05-HLA-DPB1*47:01', 'HLA-DPA1*01:05-HLA-DPB1*48:01', 'HLA-DPA1*01:05-HLA-DPB1*49:01', 'HLA-DPA1*01:05-HLA-DPB1*50:01',
         'HLA-DPA1*01:05-HLA-DPB1*51:01', 'HLA-DPA1*01:05-HLA-DPB1*52:01',
         'HLA-DPA1*01:05-HLA-DPB1*53:01', 'HLA-DPA1*01:05-HLA-DPB1*54:01', 'HLA-DPA1*01:05-HLA-DPB1*55:01', 'HLA-DPA1*01:05-HLA-DPB1*56:01',
         'HLA-DPA1*01:05-HLA-DPB1*58:01', 'HLA-DPA1*01:05-HLA-DPB1*59:01',
         'HLA-DPA1*01:05-HLA-DPB1*60:01', 'HLA-DPA1*01:05-HLA-DPB1*62:01', 'HLA-DPA1*01:05-HLA-DPB1*63:01', 'HLA-DPA1*01:05-HLA-DPB1*65:01',
         'HLA-DPA1*01:05-HLA-DPB1*66:01', 'HLA-DPA1*01:05-HLA-DPB1*67:01',
         'HLA-DPA1*01:05-HLA-DPB1*68:01', 'HLA-DPA1*01:05-HLA-DPB1*69:01', 'HLA-DPA1*01:05-HLA-DPB1*70:01', 'HLA-DPA1*01:05-HLA-DPB1*71:01',
         'HLA-DPA1*01:05-HLA-DPB1*72:01', 'HLA-DPA1*01:05-HLA-DPB1*73:01',
         'HLA-DPA1*01:05-HLA-DPB1*74:01', 'HLA-DPA1*01:05-HLA-DPB1*75:01', 'HLA-DPA1*01:05-HLA-DPB1*76:01', 'HLA-DPA1*01:05-HLA-DPB1*77:01',
         'HLA-DPA1*01:05-HLA-DPB1*78:01', 'HLA-DPA1*01:05-HLA-DPB1*79:01',
         'HLA-DPA1*01:05-HLA-DPB1*80:01', 'HLA-DPA1*01:05-HLA-DPB1*81:01', 'HLA-DPA1*01:05-HLA-DPB1*82:01', 'HLA-DPA1*01:05-HLA-DPB1*83:01',
         'HLA-DPA1*01:05-HLA-DPB1*84:01', 'HLA-DPA1*01:05-HLA-DPB1*85:01',
         'HLA-DPA1*01:05-HLA-DPB1*86:01', 'HLA-DPA1*01:05-HLA-DPB1*87:01', 'HLA-DPA1*01:05-HLA-DPB1*88:01', 'HLA-DPA1*01:05-HLA-DPB1*89:01',
         'HLA-DPA1*01:05-HLA-DPB1*90:01', 'HLA-DPA1*01:05-HLA-DPB1*91:01',
         'HLA-DPA1*01:05-HLA-DPB1*92:01', 'HLA-DPA1*01:05-HLA-DPB1*93:01', 'HLA-DPA1*01:05-HLA-DPB1*94:01', 'HLA-DPA1*01:05-HLA-DPB1*95:01',
         'HLA-DPA1*01:05-HLA-DPB1*96:01', 'HLA-DPA1*01:05-HLA-DPB1*97:01',
         'HLA-DPA1*01:05-HLA-DPB1*98:01', 'HLA-DPA1*01:05-HLA-DPB1*99:01', 'HLA-DPA1*01:06-HLA-DPB1*01:01', 'HLA-DPA1*01:06-HLA-DPB1*02:01',
         'HLA-DPA1*01:06-HLA-DPB1*02:02', 'HLA-DPA1*01:06-HLA-DPB1*03:01',
         'HLA-DPA1*01:06-HLA-DPB1*04:01', 'HLA-DPA1*01:06-HLA-DPB1*04:02', 'HLA-DPA1*01:06-HLA-DPB1*05:01', 'HLA-DPA1*01:06-HLA-DPB1*06:01',
         'HLA-DPA1*01:06-HLA-DPB1*08:01', 'HLA-DPA1*01:06-HLA-DPB1*09:01',
         'HLA-DPA1*01:06-HLA-DPB1*10:001', 'HLA-DPA1*01:06-HLA-DPB1*10:01', 'HLA-DPA1*01:06-HLA-DPB1*10:101', 'HLA-DPA1*01:06-HLA-DPB1*10:201',
         'HLA-DPA1*01:06-HLA-DPB1*10:301', 'HLA-DPA1*01:06-HLA-DPB1*10:401',
         'HLA-DPA1*01:06-HLA-DPB1*10:501', 'HLA-DPA1*01:06-HLA-DPB1*10:601', 'HLA-DPA1*01:06-HLA-DPB1*10:701', 'HLA-DPA1*01:06-HLA-DPB1*10:801',
         'HLA-DPA1*01:06-HLA-DPB1*10:901', 'HLA-DPA1*01:06-HLA-DPB1*11:001',
         'HLA-DPA1*01:06-HLA-DPB1*11:01', 'HLA-DPA1*01:06-HLA-DPB1*11:101', 'HLA-DPA1*01:06-HLA-DPB1*11:201', 'HLA-DPA1*01:06-HLA-DPB1*11:301',
         'HLA-DPA1*01:06-HLA-DPB1*11:401', 'HLA-DPA1*01:06-HLA-DPB1*11:501',
         'HLA-DPA1*01:06-HLA-DPB1*11:601', 'HLA-DPA1*01:06-HLA-DPB1*11:701', 'HLA-DPA1*01:06-HLA-DPB1*11:801', 'HLA-DPA1*01:06-HLA-DPB1*11:901',
         'HLA-DPA1*01:06-HLA-DPB1*12:101', 'HLA-DPA1*01:06-HLA-DPB1*12:201',
         'HLA-DPA1*01:06-HLA-DPB1*12:301', 'HLA-DPA1*01:06-HLA-DPB1*12:401', 'HLA-DPA1*01:06-HLA-DPB1*12:501', 'HLA-DPA1*01:06-HLA-DPB1*12:601',
         'HLA-DPA1*01:06-HLA-DPB1*12:701', 'HLA-DPA1*01:06-HLA-DPB1*12:801',
         'HLA-DPA1*01:06-HLA-DPB1*12:901', 'HLA-DPA1*01:06-HLA-DPB1*13:001', 'HLA-DPA1*01:06-HLA-DPB1*13:01', 'HLA-DPA1*01:06-HLA-DPB1*13:101',
         'HLA-DPA1*01:06-HLA-DPB1*13:201', 'HLA-DPA1*01:06-HLA-DPB1*13:301',
         'HLA-DPA1*01:06-HLA-DPB1*13:401', 'HLA-DPA1*01:06-HLA-DPB1*14:01', 'HLA-DPA1*01:06-HLA-DPB1*15:01', 'HLA-DPA1*01:06-HLA-DPB1*16:01',
         'HLA-DPA1*01:06-HLA-DPB1*17:01', 'HLA-DPA1*01:06-HLA-DPB1*18:01',
         'HLA-DPA1*01:06-HLA-DPB1*19:01', 'HLA-DPA1*01:06-HLA-DPB1*20:01', 'HLA-DPA1*01:06-HLA-DPB1*21:01', 'HLA-DPA1*01:06-HLA-DPB1*22:01',
         'HLA-DPA1*01:06-HLA-DPB1*23:01', 'HLA-DPA1*01:06-HLA-DPB1*24:01',
         'HLA-DPA1*01:06-HLA-DPB1*25:01', 'HLA-DPA1*01:06-HLA-DPB1*26:01', 'HLA-DPA1*01:06-HLA-DPB1*27:01', 'HLA-DPA1*01:06-HLA-DPB1*28:01',
         'HLA-DPA1*01:06-HLA-DPB1*29:01', 'HLA-DPA1*01:06-HLA-DPB1*30:01',
         'HLA-DPA1*01:06-HLA-DPB1*31:01', 'HLA-DPA1*01:06-HLA-DPB1*32:01', 'HLA-DPA1*01:06-HLA-DPB1*33:01', 'HLA-DPA1*01:06-HLA-DPB1*34:01',
         'HLA-DPA1*01:06-HLA-DPB1*35:01', 'HLA-DPA1*01:06-HLA-DPB1*36:01',
         'HLA-DPA1*01:06-HLA-DPB1*37:01', 'HLA-DPA1*01:06-HLA-DPB1*38:01', 'HLA-DPA1*01:06-HLA-DPB1*39:01', 'HLA-DPA1*01:06-HLA-DPB1*40:01',
         'HLA-DPA1*01:06-HLA-DPB1*41:01', 'HLA-DPA1*01:06-HLA-DPB1*44:01',
         'HLA-DPA1*01:06-HLA-DPB1*45:01', 'HLA-DPA1*01:06-HLA-DPB1*46:01', 'HLA-DPA1*01:06-HLA-DPB1*47:01', 'HLA-DPA1*01:06-HLA-DPB1*48:01',
         'HLA-DPA1*01:06-HLA-DPB1*49:01', 'HLA-DPA1*01:06-HLA-DPB1*50:01',
         'HLA-DPA1*01:06-HLA-DPB1*51:01', 'HLA-DPA1*01:06-HLA-DPB1*52:01', 'HLA-DPA1*01:06-HLA-DPB1*53:01', 'HLA-DPA1*01:06-HLA-DPB1*54:01',
         'HLA-DPA1*01:06-HLA-DPB1*55:01', 'HLA-DPA1*01:06-HLA-DPB1*56:01',
         'HLA-DPA1*01:06-HLA-DPB1*58:01', 'HLA-DPA1*01:06-HLA-DPB1*59:01', 'HLA-DPA1*01:06-HLA-DPB1*60:01', 'HLA-DPA1*01:06-HLA-DPB1*62:01',
         'HLA-DPA1*01:06-HLA-DPB1*63:01', 'HLA-DPA1*01:06-HLA-DPB1*65:01',
         'HLA-DPA1*01:06-HLA-DPB1*66:01', 'HLA-DPA1*01:06-HLA-DPB1*67:01', 'HLA-DPA1*01:06-HLA-DPB1*68:01', 'HLA-DPA1*01:06-HLA-DPB1*69:01',
         'HLA-DPA1*01:06-HLA-DPB1*70:01', 'HLA-DPA1*01:06-HLA-DPB1*71:01',
         'HLA-DPA1*01:06-HLA-DPB1*72:01', 'HLA-DPA1*01:06-HLA-DPB1*73:01', 'HLA-DPA1*01:06-HLA-DPB1*74:01', 'HLA-DPA1*01:06-HLA-DPB1*75:01',
         'HLA-DPA1*01:06-HLA-DPB1*76:01', 'HLA-DPA1*01:06-HLA-DPB1*77:01',
         'HLA-DPA1*01:06-HLA-DPB1*78:01', 'HLA-DPA1*01:06-HLA-DPB1*79:01', 'HLA-DPA1*01:06-HLA-DPB1*80:01', 'HLA-DPA1*01:06-HLA-DPB1*81:01',
         'HLA-DPA1*01:06-HLA-DPB1*82:01', 'HLA-DPA1*01:06-HLA-DPB1*83:01',
         'HLA-DPA1*01:06-HLA-DPB1*84:01', 'HLA-DPA1*01:06-HLA-DPB1*85:01', 'HLA-DPA1*01:06-HLA-DPB1*86:01', 'HLA-DPA1*01:06-HLA-DPB1*87:01',
         'HLA-DPA1*01:06-HLA-DPB1*88:01', 'HLA-DPA1*01:06-HLA-DPB1*89:01',
         'HLA-DPA1*01:06-HLA-DPB1*90:01', 'HLA-DPA1*01:06-HLA-DPB1*91:01', 'HLA-DPA1*01:06-HLA-DPB1*92:01', 'HLA-DPA1*01:06-HLA-DPB1*93:01',
         'HLA-DPA1*01:06-HLA-DPB1*94:01', 'HLA-DPA1*01:06-HLA-DPB1*95:01',
         'HLA-DPA1*01:06-HLA-DPB1*96:01', 'HLA-DPA1*01:06-HLA-DPB1*97:01', 'HLA-DPA1*01:06-HLA-DPB1*98:01', 'HLA-DPA1*01:06-HLA-DPB1*99:01',
         'HLA-DPA1*01:07-HLA-DPB1*01:01', 'HLA-DPA1*01:07-HLA-DPB1*02:01',
         'HLA-DPA1*01:07-HLA-DPB1*02:02', 'HLA-DPA1*01:07-HLA-DPB1*03:01', 'HLA-DPA1*01:07-HLA-DPB1*04:01', 'HLA-DPA1*01:07-HLA-DPB1*04:02',
         'HLA-DPA1*01:07-HLA-DPB1*05:01', 'HLA-DPA1*01:07-HLA-DPB1*06:01',
         'HLA-DPA1*01:07-HLA-DPB1*08:01', 'HLA-DPA1*01:07-HLA-DPB1*09:01', 'HLA-DPA1*01:07-HLA-DPB1*10:001', 'HLA-DPA1*01:07-HLA-DPB1*10:01',
         'HLA-DPA1*01:07-HLA-DPB1*10:101', 'HLA-DPA1*01:07-HLA-DPB1*10:201',
         'HLA-DPA1*01:07-HLA-DPB1*10:301', 'HLA-DPA1*01:07-HLA-DPB1*10:401', 'HLA-DPA1*01:07-HLA-DPB1*10:501', 'HLA-DPA1*01:07-HLA-DPB1*10:601',
         'HLA-DPA1*01:07-HLA-DPB1*10:701', 'HLA-DPA1*01:07-HLA-DPB1*10:801',
         'HLA-DPA1*01:07-HLA-DPB1*10:901', 'HLA-DPA1*01:07-HLA-DPB1*11:001', 'HLA-DPA1*01:07-HLA-DPB1*11:01', 'HLA-DPA1*01:07-HLA-DPB1*11:101',
         'HLA-DPA1*01:07-HLA-DPB1*11:201', 'HLA-DPA1*01:07-HLA-DPB1*11:301',
         'HLA-DPA1*01:07-HLA-DPB1*11:401', 'HLA-DPA1*01:07-HLA-DPB1*11:501', 'HLA-DPA1*01:07-HLA-DPB1*11:601', 'HLA-DPA1*01:07-HLA-DPB1*11:701',
         'HLA-DPA1*01:07-HLA-DPB1*11:801', 'HLA-DPA1*01:07-HLA-DPB1*11:901',
         'HLA-DPA1*01:07-HLA-DPB1*12:101', 'HLA-DPA1*01:07-HLA-DPB1*12:201', 'HLA-DPA1*01:07-HLA-DPB1*12:301', 'HLA-DPA1*01:07-HLA-DPB1*12:401',
         'HLA-DPA1*01:07-HLA-DPB1*12:501', 'HLA-DPA1*01:07-HLA-DPB1*12:601',
         'HLA-DPA1*01:07-HLA-DPB1*12:701', 'HLA-DPA1*01:07-HLA-DPB1*12:801', 'HLA-DPA1*01:07-HLA-DPB1*12:901', 'HLA-DPA1*01:07-HLA-DPB1*13:001',
         'HLA-DPA1*01:07-HLA-DPB1*13:01', 'HLA-DPA1*01:07-HLA-DPB1*13:101',
         'HLA-DPA1*01:07-HLA-DPB1*13:201', 'HLA-DPA1*01:07-HLA-DPB1*13:301', 'HLA-DPA1*01:07-HLA-DPB1*13:401', 'HLA-DPA1*01:07-HLA-DPB1*14:01',
         'HLA-DPA1*01:07-HLA-DPB1*15:01', 'HLA-DPA1*01:07-HLA-DPB1*16:01',
         'HLA-DPA1*01:07-HLA-DPB1*17:01', 'HLA-DPA1*01:07-HLA-DPB1*18:01', 'HLA-DPA1*01:07-HLA-DPB1*19:01', 'HLA-DPA1*01:07-HLA-DPB1*20:01',
         'HLA-DPA1*01:07-HLA-DPB1*21:01', 'HLA-DPA1*01:07-HLA-DPB1*22:01',
         'HLA-DPA1*01:07-HLA-DPB1*23:01', 'HLA-DPA1*01:07-HLA-DPB1*24:01', 'HLA-DPA1*01:07-HLA-DPB1*25:01', 'HLA-DPA1*01:07-HLA-DPB1*26:01',
         'HLA-DPA1*01:07-HLA-DPB1*27:01', 'HLA-DPA1*01:07-HLA-DPB1*28:01',
         'HLA-DPA1*01:07-HLA-DPB1*29:01', 'HLA-DPA1*01:07-HLA-DPB1*30:01', 'HLA-DPA1*01:07-HLA-DPB1*31:01', 'HLA-DPA1*01:07-HLA-DPB1*32:01',
         'HLA-DPA1*01:07-HLA-DPB1*33:01', 'HLA-DPA1*01:07-HLA-DPB1*34:01',
         'HLA-DPA1*01:07-HLA-DPB1*35:01', 'HLA-DPA1*01:07-HLA-DPB1*36:01', 'HLA-DPA1*01:07-HLA-DPB1*37:01', 'HLA-DPA1*01:07-HLA-DPB1*38:01',
         'HLA-DPA1*01:07-HLA-DPB1*39:01', 'HLA-DPA1*01:07-HLA-DPB1*40:01',
         'HLA-DPA1*01:07-HLA-DPB1*41:01', 'HLA-DPA1*01:07-HLA-DPB1*44:01', 'HLA-DPA1*01:07-HLA-DPB1*45:01', 'HLA-DPA1*01:07-HLA-DPB1*46:01',
         'HLA-DPA1*01:07-HLA-DPB1*47:01', 'HLA-DPA1*01:07-HLA-DPB1*48:01',
         'HLA-DPA1*01:07-HLA-DPB1*49:01', 'HLA-DPA1*01:07-HLA-DPB1*50:01', 'HLA-DPA1*01:07-HLA-DPB1*51:01', 'HLA-DPA1*01:07-HLA-DPB1*52:01',
         'HLA-DPA1*01:07-HLA-DPB1*53:01', 'HLA-DPA1*01:07-HLA-DPB1*54:01',
         'HLA-DPA1*01:07-HLA-DPB1*55:01', 'HLA-DPA1*01:07-HLA-DPB1*56:01', 'HLA-DPA1*01:07-HLA-DPB1*58:01', 'HLA-DPA1*01:07-HLA-DPB1*59:01',
         'HLA-DPA1*01:07-HLA-DPB1*60:01', 'HLA-DPA1*01:07-HLA-DPB1*62:01',
         'HLA-DPA1*01:07-HLA-DPB1*63:01', 'HLA-DPA1*01:07-HLA-DPB1*65:01', 'HLA-DPA1*01:07-HLA-DPB1*66:01', 'HLA-DPA1*01:07-HLA-DPB1*67:01',
         'HLA-DPA1*01:07-HLA-DPB1*68:01', 'HLA-DPA1*01:07-HLA-DPB1*69:01',
         'HLA-DPA1*01:07-HLA-DPB1*70:01', 'HLA-DPA1*01:07-HLA-DPB1*71:01', 'HLA-DPA1*01:07-HLA-DPB1*72:01', 'HLA-DPA1*01:07-HLA-DPB1*73:01',
         'HLA-DPA1*01:07-HLA-DPB1*74:01', 'HLA-DPA1*01:07-HLA-DPB1*75:01',
         'HLA-DPA1*01:07-HLA-DPB1*76:01', 'HLA-DPA1*01:07-HLA-DPB1*77:01', 'HLA-DPA1*01:07-HLA-DPB1*78:01', 'HLA-DPA1*01:07-HLA-DPB1*79:01',
         'HLA-DPA1*01:07-HLA-DPB1*80:01', 'HLA-DPA1*01:07-HLA-DPB1*81:01',
         'HLA-DPA1*01:07-HLA-DPB1*82:01', 'HLA-DPA1*01:07-HLA-DPB1*83:01', 'HLA-DPA1*01:07-HLA-DPB1*84:01', 'HLA-DPA1*01:07-HLA-DPB1*85:01',
         'HLA-DPA1*01:07-HLA-DPB1*86:01', 'HLA-DPA1*01:07-HLA-DPB1*87:01',
         'HLA-DPA1*01:07-HLA-DPB1*88:01', 'HLA-DPA1*01:07-HLA-DPB1*89:01', 'HLA-DPA1*01:07-HLA-DPB1*90:01', 'HLA-DPA1*01:07-HLA-DPB1*91:01',
         'HLA-DPA1*01:07-HLA-DPB1*92:01', 'HLA-DPA1*01:07-HLA-DPB1*93:01',
         'HLA-DPA1*01:07-HLA-DPB1*94:01', 'HLA-DPA1*01:07-HLA-DPB1*95:01', 'HLA-DPA1*01:07-HLA-DPB1*96:01', 'HLA-DPA1*01:07-HLA-DPB1*97:01',
         'HLA-DPA1*01:07-HLA-DPB1*98:01', 'HLA-DPA1*01:07-HLA-DPB1*99:01',
         'HLA-DPA1*01:08-HLA-DPB1*01:01', 'HLA-DPA1*01:08-HLA-DPB1*02:01', 'HLA-DPA1*01:08-HLA-DPB1*02:02', 'HLA-DPA1*01:08-HLA-DPB1*03:01',
         'HLA-DPA1*01:08-HLA-DPB1*04:01', 'HLA-DPA1*01:08-HLA-DPB1*04:02',
         'HLA-DPA1*01:08-HLA-DPB1*05:01', 'HLA-DPA1*01:08-HLA-DPB1*06:01', 'HLA-DPA1*01:08-HLA-DPB1*08:01', 'HLA-DPA1*01:08-HLA-DPB1*09:01',
         'HLA-DPA1*01:08-HLA-DPB1*10:001', 'HLA-DPA1*01:08-HLA-DPB1*10:01',
         'HLA-DPA1*01:08-HLA-DPB1*10:101', 'HLA-DPA1*01:08-HLA-DPB1*10:201', 'HLA-DPA1*01:08-HLA-DPB1*10:301', 'HLA-DPA1*01:08-HLA-DPB1*10:401',
         'HLA-DPA1*01:08-HLA-DPB1*10:501', 'HLA-DPA1*01:08-HLA-DPB1*10:601',
         'HLA-DPA1*01:08-HLA-DPB1*10:701', 'HLA-DPA1*01:08-HLA-DPB1*10:801', 'HLA-DPA1*01:08-HLA-DPB1*10:901', 'HLA-DPA1*01:08-HLA-DPB1*11:001',
         'HLA-DPA1*01:08-HLA-DPB1*11:01', 'HLA-DPA1*01:08-HLA-DPB1*11:101',
         'HLA-DPA1*01:08-HLA-DPB1*11:201', 'HLA-DPA1*01:08-HLA-DPB1*11:301', 'HLA-DPA1*01:08-HLA-DPB1*11:401', 'HLA-DPA1*01:08-HLA-DPB1*11:501',
         'HLA-DPA1*01:08-HLA-DPB1*11:601', 'HLA-DPA1*01:08-HLA-DPB1*11:701',
         'HLA-DPA1*01:08-HLA-DPB1*11:801', 'HLA-DPA1*01:08-HLA-DPB1*11:901', 'HLA-DPA1*01:08-HLA-DPB1*12:101', 'HLA-DPA1*01:08-HLA-DPB1*12:201',
         'HLA-DPA1*01:08-HLA-DPB1*12:301', 'HLA-DPA1*01:08-HLA-DPB1*12:401',
         'HLA-DPA1*01:08-HLA-DPB1*12:501', 'HLA-DPA1*01:08-HLA-DPB1*12:601', 'HLA-DPA1*01:08-HLA-DPB1*12:701', 'HLA-DPA1*01:08-HLA-DPB1*12:801',
         'HLA-DPA1*01:08-HLA-DPB1*12:901', 'HLA-DPA1*01:08-HLA-DPB1*13:001',
         'HLA-DPA1*01:08-HLA-DPB1*13:01', 'HLA-DPA1*01:08-HLA-DPB1*13:101', 'HLA-DPA1*01:08-HLA-DPB1*13:201', 'HLA-DPA1*01:08-HLA-DPB1*13:301',
         'HLA-DPA1*01:08-HLA-DPB1*13:401', 'HLA-DPA1*01:08-HLA-DPB1*14:01',
         'HLA-DPA1*01:08-HLA-DPB1*15:01', 'HLA-DPA1*01:08-HLA-DPB1*16:01', 'HLA-DPA1*01:08-HLA-DPB1*17:01', 'HLA-DPA1*01:08-HLA-DPB1*18:01',
         'HLA-DPA1*01:08-HLA-DPB1*19:01', 'HLA-DPA1*01:08-HLA-DPB1*20:01',
         'HLA-DPA1*01:08-HLA-DPB1*21:01', 'HLA-DPA1*01:08-HLA-DPB1*22:01', 'HLA-DPA1*01:08-HLA-DPB1*23:01', 'HLA-DPA1*01:08-HLA-DPB1*24:01',
         'HLA-DPA1*01:08-HLA-DPB1*25:01', 'HLA-DPA1*01:08-HLA-DPB1*26:01',
         'HLA-DPA1*01:08-HLA-DPB1*27:01', 'HLA-DPA1*01:08-HLA-DPB1*28:01', 'HLA-DPA1*01:08-HLA-DPB1*29:01', 'HLA-DPA1*01:08-HLA-DPB1*30:01',
         'HLA-DPA1*01:08-HLA-DPB1*31:01', 'HLA-DPA1*01:08-HLA-DPB1*32:01',
         'HLA-DPA1*01:08-HLA-DPB1*33:01', 'HLA-DPA1*01:08-HLA-DPB1*34:01', 'HLA-DPA1*01:08-HLA-DPB1*35:01', 'HLA-DPA1*01:08-HLA-DPB1*36:01',
         'HLA-DPA1*01:08-HLA-DPB1*37:01', 'HLA-DPA1*01:08-HLA-DPB1*38:01',
         'HLA-DPA1*01:08-HLA-DPB1*39:01', 'HLA-DPA1*01:08-HLA-DPB1*40:01', 'HLA-DPA1*01:08-HLA-DPB1*41:01', 'HLA-DPA1*01:08-HLA-DPB1*44:01',
         'HLA-DPA1*01:08-HLA-DPB1*45:01', 'HLA-DPA1*01:08-HLA-DPB1*46:01',
         'HLA-DPA1*01:08-HLA-DPB1*47:01', 'HLA-DPA1*01:08-HLA-DPB1*48:01', 'HLA-DPA1*01:08-HLA-DPB1*49:01', 'HLA-DPA1*01:08-HLA-DPB1*50:01',
         'HLA-DPA1*01:08-HLA-DPB1*51:01', 'HLA-DPA1*01:08-HLA-DPB1*52:01',
         'HLA-DPA1*01:08-HLA-DPB1*53:01', 'HLA-DPA1*01:08-HLA-DPB1*54:01', 'HLA-DPA1*01:08-HLA-DPB1*55:01', 'HLA-DPA1*01:08-HLA-DPB1*56:01',
         'HLA-DPA1*01:08-HLA-DPB1*58:01', 'HLA-DPA1*01:08-HLA-DPB1*59:01',
         'HLA-DPA1*01:08-HLA-DPB1*60:01', 'HLA-DPA1*01:08-HLA-DPB1*62:01', 'HLA-DPA1*01:08-HLA-DPB1*63:01', 'HLA-DPA1*01:08-HLA-DPB1*65:01',
         'HLA-DPA1*01:08-HLA-DPB1*66:01', 'HLA-DPA1*01:08-HLA-DPB1*67:01',
         'HLA-DPA1*01:08-HLA-DPB1*68:01', 'HLA-DPA1*01:08-HLA-DPB1*69:01', 'HLA-DPA1*01:08-HLA-DPB1*70:01', 'HLA-DPA1*01:08-HLA-DPB1*71:01',
         'HLA-DPA1*01:08-HLA-DPB1*72:01', 'HLA-DPA1*01:08-HLA-DPB1*73:01',
         'HLA-DPA1*01:08-HLA-DPB1*74:01', 'HLA-DPA1*01:08-HLA-DPB1*75:01', 'HLA-DPA1*01:08-HLA-DPB1*76:01', 'HLA-DPA1*01:08-HLA-DPB1*77:01',
         'HLA-DPA1*01:08-HLA-DPB1*78:01', 'HLA-DPA1*01:08-HLA-DPB1*79:01',
         'HLA-DPA1*01:08-HLA-DPB1*80:01', 'HLA-DPA1*01:08-HLA-DPB1*81:01', 'HLA-DPA1*01:08-HLA-DPB1*82:01', 'HLA-DPA1*01:08-HLA-DPB1*83:01',
         'HLA-DPA1*01:08-HLA-DPB1*84:01', 'HLA-DPA1*01:08-HLA-DPB1*85:01',
         'HLA-DPA1*01:08-HLA-DPB1*86:01', 'HLA-DPA1*01:08-HLA-DPB1*87:01', 'HLA-DPA1*01:08-HLA-DPB1*88:01', 'HLA-DPA1*01:08-HLA-DPB1*89:01',
         'HLA-DPA1*01:08-HLA-DPB1*90:01', 'HLA-DPA1*01:08-HLA-DPB1*91:01',
         'HLA-DPA1*01:08-HLA-DPB1*92:01', 'HLA-DPA1*01:08-HLA-DPB1*93:01', 'HLA-DPA1*01:08-HLA-DPB1*94:01', 'HLA-DPA1*01:08-HLA-DPB1*95:01',
         'HLA-DPA1*01:08-HLA-DPB1*96:01', 'HLA-DPA1*01:08-HLA-DPB1*97:01',
         'HLA-DPA1*01:08-HLA-DPB1*98:01', 'HLA-DPA1*01:08-HLA-DPB1*99:01', 'HLA-DPA1*01:09-HLA-DPB1*01:01', 'HLA-DPA1*01:09-HLA-DPB1*02:01',
         'HLA-DPA1*01:09-HLA-DPB1*02:02', 'HLA-DPA1*01:09-HLA-DPB1*03:01',
         'HLA-DPA1*01:09-HLA-DPB1*04:01', 'HLA-DPA1*01:09-HLA-DPB1*04:02', 'HLA-DPA1*01:09-HLA-DPB1*05:01', 'HLA-DPA1*01:09-HLA-DPB1*06:01',
         'HLA-DPA1*01:09-HLA-DPB1*08:01', 'HLA-DPA1*01:09-HLA-DPB1*09:01',
         'HLA-DPA1*01:09-HLA-DPB1*10:001', 'HLA-DPA1*01:09-HLA-DPB1*10:01', 'HLA-DPA1*01:09-HLA-DPB1*10:101', 'HLA-DPA1*01:09-HLA-DPB1*10:201',
         'HLA-DPA1*01:09-HLA-DPB1*10:301', 'HLA-DPA1*01:09-HLA-DPB1*10:401',
         'HLA-DPA1*01:09-HLA-DPB1*10:501', 'HLA-DPA1*01:09-HLA-DPB1*10:601', 'HLA-DPA1*01:09-HLA-DPB1*10:701', 'HLA-DPA1*01:09-HLA-DPB1*10:801',
         'HLA-DPA1*01:09-HLA-DPB1*10:901', 'HLA-DPA1*01:09-HLA-DPB1*11:001',
         'HLA-DPA1*01:09-HLA-DPB1*11:01', 'HLA-DPA1*01:09-HLA-DPB1*11:101', 'HLA-DPA1*01:09-HLA-DPB1*11:201', 'HLA-DPA1*01:09-HLA-DPB1*11:301',
         'HLA-DPA1*01:09-HLA-DPB1*11:401', 'HLA-DPA1*01:09-HLA-DPB1*11:501',
         'HLA-DPA1*01:09-HLA-DPB1*11:601', 'HLA-DPA1*01:09-HLA-DPB1*11:701', 'HLA-DPA1*01:09-HLA-DPB1*11:801', 'HLA-DPA1*01:09-HLA-DPB1*11:901',
         'HLA-DPA1*01:09-HLA-DPB1*12:101', 'HLA-DPA1*01:09-HLA-DPB1*12:201',
         'HLA-DPA1*01:09-HLA-DPB1*12:301', 'HLA-DPA1*01:09-HLA-DPB1*12:401', 'HLA-DPA1*01:09-HLA-DPB1*12:501', 'HLA-DPA1*01:09-HLA-DPB1*12:601',
         'HLA-DPA1*01:09-HLA-DPB1*12:701', 'HLA-DPA1*01:09-HLA-DPB1*12:801',
         'HLA-DPA1*01:09-HLA-DPB1*12:901', 'HLA-DPA1*01:09-HLA-DPB1*13:001', 'HLA-DPA1*01:09-HLA-DPB1*13:01', 'HLA-DPA1*01:09-HLA-DPB1*13:101',
         'HLA-DPA1*01:09-HLA-DPB1*13:201', 'HLA-DPA1*01:09-HLA-DPB1*13:301',
         'HLA-DPA1*01:09-HLA-DPB1*13:401', 'HLA-DPA1*01:09-HLA-DPB1*14:01', 'HLA-DPA1*01:09-HLA-DPB1*15:01', 'HLA-DPA1*01:09-HLA-DPB1*16:01',
         'HLA-DPA1*01:09-HLA-DPB1*17:01', 'HLA-DPA1*01:09-HLA-DPB1*18:01',
         'HLA-DPA1*01:09-HLA-DPB1*19:01', 'HLA-DPA1*01:09-HLA-DPB1*20:01', 'HLA-DPA1*01:09-HLA-DPB1*21:01', 'HLA-DPA1*01:09-HLA-DPB1*22:01',
         'HLA-DPA1*01:09-HLA-DPB1*23:01', 'HLA-DPA1*01:09-HLA-DPB1*24:01',
         'HLA-DPA1*01:09-HLA-DPB1*25:01', 'HLA-DPA1*01:09-HLA-DPB1*26:01', 'HLA-DPA1*01:09-HLA-DPB1*27:01', 'HLA-DPA1*01:09-HLA-DPB1*28:01',
         'HLA-DPA1*01:09-HLA-DPB1*29:01', 'HLA-DPA1*01:09-HLA-DPB1*30:01',
         'HLA-DPA1*01:09-HLA-DPB1*31:01', 'HLA-DPA1*01:09-HLA-DPB1*32:01', 'HLA-DPA1*01:09-HLA-DPB1*33:01', 'HLA-DPA1*01:09-HLA-DPB1*34:01',
         'HLA-DPA1*01:09-HLA-DPB1*35:01', 'HLA-DPA1*01:09-HLA-DPB1*36:01',
         'HLA-DPA1*01:09-HLA-DPB1*37:01', 'HLA-DPA1*01:09-HLA-DPB1*38:01', 'HLA-DPA1*01:09-HLA-DPB1*39:01', 'HLA-DPA1*01:09-HLA-DPB1*40:01',
         'HLA-DPA1*01:09-HLA-DPB1*41:01', 'HLA-DPA1*01:09-HLA-DPB1*44:01',
         'HLA-DPA1*01:09-HLA-DPB1*45:01', 'HLA-DPA1*01:09-HLA-DPB1*46:01', 'HLA-DPA1*01:09-HLA-DPB1*47:01', 'HLA-DPA1*01:09-HLA-DPB1*48:01',
         'HLA-DPA1*01:09-HLA-DPB1*49:01', 'HLA-DPA1*01:09-HLA-DPB1*50:01',
         'HLA-DPA1*01:09-HLA-DPB1*51:01', 'HLA-DPA1*01:09-HLA-DPB1*52:01', 'HLA-DPA1*01:09-HLA-DPB1*53:01', 'HLA-DPA1*01:09-HLA-DPB1*54:01',
         'HLA-DPA1*01:09-HLA-DPB1*55:01', 'HLA-DPA1*01:09-HLA-DPB1*56:01',
         'HLA-DPA1*01:09-HLA-DPB1*58:01', 'HLA-DPA1*01:09-HLA-DPB1*59:01', 'HLA-DPA1*01:09-HLA-DPB1*60:01', 'HLA-DPA1*01:09-HLA-DPB1*62:01',
         'HLA-DPA1*01:09-HLA-DPB1*63:01', 'HLA-DPA1*01:09-HLA-DPB1*65:01',
         'HLA-DPA1*01:09-HLA-DPB1*66:01', 'HLA-DPA1*01:09-HLA-DPB1*67:01', 'HLA-DPA1*01:09-HLA-DPB1*68:01', 'HLA-DPA1*01:09-HLA-DPB1*69:01',
         'HLA-DPA1*01:09-HLA-DPB1*70:01', 'HLA-DPA1*01:09-HLA-DPB1*71:01',
         'HLA-DPA1*01:09-HLA-DPB1*72:01', 'HLA-DPA1*01:09-HLA-DPB1*73:01', 'HLA-DPA1*01:09-HLA-DPB1*74:01', 'HLA-DPA1*01:09-HLA-DPB1*75:01',
         'HLA-DPA1*01:09-HLA-DPB1*76:01', 'HLA-DPA1*01:09-HLA-DPB1*77:01',
         'HLA-DPA1*01:09-HLA-DPB1*78:01', 'HLA-DPA1*01:09-HLA-DPB1*79:01', 'HLA-DPA1*01:09-HLA-DPB1*80:01', 'HLA-DPA1*01:09-HLA-DPB1*81:01',
         'HLA-DPA1*01:09-HLA-DPB1*82:01', 'HLA-DPA1*01:09-HLA-DPB1*83:01',
         'HLA-DPA1*01:09-HLA-DPB1*84:01', 'HLA-DPA1*01:09-HLA-DPB1*85:01', 'HLA-DPA1*01:09-HLA-DPB1*86:01', 'HLA-DPA1*01:09-HLA-DPB1*87:01',
         'HLA-DPA1*01:09-HLA-DPB1*88:01', 'HLA-DPA1*01:09-HLA-DPB1*89:01',
         'HLA-DPA1*01:09-HLA-DPB1*90:01', 'HLA-DPA1*01:09-HLA-DPB1*91:01', 'HLA-DPA1*01:09-HLA-DPB1*92:01', 'HLA-DPA1*01:09-HLA-DPB1*93:01',
         'HLA-DPA1*01:09-HLA-DPB1*94:01', 'HLA-DPA1*01:09-HLA-DPB1*95:01',
         'HLA-DPA1*01:09-HLA-DPB1*96:01', 'HLA-DPA1*01:09-HLA-DPB1*97:01', 'HLA-DPA1*01:09-HLA-DPB1*98:01', 'HLA-DPA1*01:09-HLA-DPB1*99:01',
         'HLA-DPA1*01:10-HLA-DPB1*01:01', 'HLA-DPA1*01:10-HLA-DPB1*02:01',
         'HLA-DPA1*01:10-HLA-DPB1*02:02', 'HLA-DPA1*01:10-HLA-DPB1*03:01', 'HLA-DPA1*01:10-HLA-DPB1*04:01', 'HLA-DPA1*01:10-HLA-DPB1*04:02',
         'HLA-DPA1*01:10-HLA-DPB1*05:01', 'HLA-DPA1*01:10-HLA-DPB1*06:01',
         'HLA-DPA1*01:10-HLA-DPB1*08:01', 'HLA-DPA1*01:10-HLA-DPB1*09:01', 'HLA-DPA1*01:10-HLA-DPB1*10:001', 'HLA-DPA1*01:10-HLA-DPB1*10:01',
         'HLA-DPA1*01:10-HLA-DPB1*10:101', 'HLA-DPA1*01:10-HLA-DPB1*10:201',
         'HLA-DPA1*01:10-HLA-DPB1*10:301', 'HLA-DPA1*01:10-HLA-DPB1*10:401', 'HLA-DPA1*01:10-HLA-DPB1*10:501', 'HLA-DPA1*01:10-HLA-DPB1*10:601',
         'HLA-DPA1*01:10-HLA-DPB1*10:701', 'HLA-DPA1*01:10-HLA-DPB1*10:801',
         'HLA-DPA1*01:10-HLA-DPB1*10:901', 'HLA-DPA1*01:10-HLA-DPB1*11:001', 'HLA-DPA1*01:10-HLA-DPB1*11:01', 'HLA-DPA1*01:10-HLA-DPB1*11:101',
         'HLA-DPA1*01:10-HLA-DPB1*11:201', 'HLA-DPA1*01:10-HLA-DPB1*11:301',
         'HLA-DPA1*01:10-HLA-DPB1*11:401', 'HLA-DPA1*01:10-HLA-DPB1*11:501', 'HLA-DPA1*01:10-HLA-DPB1*11:601', 'HLA-DPA1*01:10-HLA-DPB1*11:701',
         'HLA-DPA1*01:10-HLA-DPB1*11:801', 'HLA-DPA1*01:10-HLA-DPB1*11:901',
         'HLA-DPA1*01:10-HLA-DPB1*12:101', 'HLA-DPA1*01:10-HLA-DPB1*12:201', 'HLA-DPA1*01:10-HLA-DPB1*12:301', 'HLA-DPA1*01:10-HLA-DPB1*12:401',
         'HLA-DPA1*01:10-HLA-DPB1*12:501', 'HLA-DPA1*01:10-HLA-DPB1*12:601',
         'HLA-DPA1*01:10-HLA-DPB1*12:701', 'HLA-DPA1*01:10-HLA-DPB1*12:801', 'HLA-DPA1*01:10-HLA-DPB1*12:901', 'HLA-DPA1*01:10-HLA-DPB1*13:001',
         'HLA-DPA1*01:10-HLA-DPB1*13:01', 'HLA-DPA1*01:10-HLA-DPB1*13:101',
         'HLA-DPA1*01:10-HLA-DPB1*13:201', 'HLA-DPA1*01:10-HLA-DPB1*13:301', 'HLA-DPA1*01:10-HLA-DPB1*13:401', 'HLA-DPA1*01:10-HLA-DPB1*14:01',
         'HLA-DPA1*01:10-HLA-DPB1*15:01', 'HLA-DPA1*01:10-HLA-DPB1*16:01',
         'HLA-DPA1*01:10-HLA-DPB1*17:01', 'HLA-DPA1*01:10-HLA-DPB1*18:01', 'HLA-DPA1*01:10-HLA-DPB1*19:01', 'HLA-DPA1*01:10-HLA-DPB1*20:01',
         'HLA-DPA1*01:10-HLA-DPB1*21:01', 'HLA-DPA1*01:10-HLA-DPB1*22:01',
         'HLA-DPA1*01:10-HLA-DPB1*23:01', 'HLA-DPA1*01:10-HLA-DPB1*24:01', 'HLA-DPA1*01:10-HLA-DPB1*25:01', 'HLA-DPA1*01:10-HLA-DPB1*26:01',
         'HLA-DPA1*01:10-HLA-DPB1*27:01', 'HLA-DPA1*01:10-HLA-DPB1*28:01',
         'HLA-DPA1*01:10-HLA-DPB1*29:01', 'HLA-DPA1*01:10-HLA-DPB1*30:01', 'HLA-DPA1*01:10-HLA-DPB1*31:01', 'HLA-DPA1*01:10-HLA-DPB1*32:01',
         'HLA-DPA1*01:10-HLA-DPB1*33:01', 'HLA-DPA1*01:10-HLA-DPB1*34:01',
         'HLA-DPA1*01:10-HLA-DPB1*35:01', 'HLA-DPA1*01:10-HLA-DPB1*36:01', 'HLA-DPA1*01:10-HLA-DPB1*37:01', 'HLA-DPA1*01:10-HLA-DPB1*38:01',
         'HLA-DPA1*01:10-HLA-DPB1*39:01', 'HLA-DPA1*01:10-HLA-DPB1*40:01',
         'HLA-DPA1*01:10-HLA-DPB1*41:01', 'HLA-DPA1*01:10-HLA-DPB1*44:01', 'HLA-DPA1*01:10-HLA-DPB1*45:01', 'HLA-DPA1*01:10-HLA-DPB1*46:01',
         'HLA-DPA1*01:10-HLA-DPB1*47:01', 'HLA-DPA1*01:10-HLA-DPB1*48:01',
         'HLA-DPA1*01:10-HLA-DPB1*49:01', 'HLA-DPA1*01:10-HLA-DPB1*50:01', 'HLA-DPA1*01:10-HLA-DPB1*51:01', 'HLA-DPA1*01:10-HLA-DPB1*52:01',
         'HLA-DPA1*01:10-HLA-DPB1*53:01', 'HLA-DPA1*01:10-HLA-DPB1*54:01',
         'HLA-DPA1*01:10-HLA-DPB1*55:01', 'HLA-DPA1*01:10-HLA-DPB1*56:01', 'HLA-DPA1*01:10-HLA-DPB1*58:01', 'HLA-DPA1*01:10-HLA-DPB1*59:01',
         'HLA-DPA1*01:10-HLA-DPB1*60:01', 'HLA-DPA1*01:10-HLA-DPB1*62:01',
         'HLA-DPA1*01:10-HLA-DPB1*63:01', 'HLA-DPA1*01:10-HLA-DPB1*65:01', 'HLA-DPA1*01:10-HLA-DPB1*66:01', 'HLA-DPA1*01:10-HLA-DPB1*67:01',
         'HLA-DPA1*01:10-HLA-DPB1*68:01', 'HLA-DPA1*01:10-HLA-DPB1*69:01',
         'HLA-DPA1*01:10-HLA-DPB1*70:01', 'HLA-DPA1*01:10-HLA-DPB1*71:01', 'HLA-DPA1*01:10-HLA-DPB1*72:01', 'HLA-DPA1*01:10-HLA-DPB1*73:01',
         'HLA-DPA1*01:10-HLA-DPB1*74:01', 'HLA-DPA1*01:10-HLA-DPB1*75:01',
         'HLA-DPA1*01:10-HLA-DPB1*76:01', 'HLA-DPA1*01:10-HLA-DPB1*77:01', 'HLA-DPA1*01:10-HLA-DPB1*78:01', 'HLA-DPA1*01:10-HLA-DPB1*79:01',
         'HLA-DPA1*01:10-HLA-DPB1*80:01', 'HLA-DPA1*01:10-HLA-DPB1*81:01',
         'HLA-DPA1*01:10-HLA-DPB1*82:01', 'HLA-DPA1*01:10-HLA-DPB1*83:01', 'HLA-DPA1*01:10-HLA-DPB1*84:01', 'HLA-DPA1*01:10-HLA-DPB1*85:01',
         'HLA-DPA1*01:10-HLA-DPB1*86:01', 'HLA-DPA1*01:10-HLA-DPB1*87:01',
         'HLA-DPA1*01:10-HLA-DPB1*88:01', 'HLA-DPA1*01:10-HLA-DPB1*89:01', 'HLA-DPA1*01:10-HLA-DPB1*90:01', 'HLA-DPA1*01:10-HLA-DPB1*91:01',
         'HLA-DPA1*01:10-HLA-DPB1*92:01', 'HLA-DPA1*01:10-HLA-DPB1*93:01',
         'HLA-DPA1*01:10-HLA-DPB1*94:01', 'HLA-DPA1*01:10-HLA-DPB1*95:01', 'HLA-DPA1*01:10-HLA-DPB1*96:01', 'HLA-DPA1*01:10-HLA-DPB1*97:01',
         'HLA-DPA1*01:10-HLA-DPB1*98:01', 'HLA-DPA1*01:10-HLA-DPB1*99:01',
         'HLA-DPA1*02:01-HLA-DPB1*01:01', 'HLA-DPA1*02:01-HLA-DPB1*02:01', 'HLA-DPA1*02:01-HLA-DPB1*02:02', 'HLA-DPA1*02:01-HLA-DPB1*03:01',
         'HLA-DPA1*02:01-HLA-DPB1*04:01', 'HLA-DPA1*02:01-HLA-DPB1*04:02',
         'HLA-DPA1*02:01-HLA-DPB1*05:01', 'HLA-DPA1*02:01-HLA-DPB1*06:01', 'HLA-DPA1*02:01-HLA-DPB1*08:01', 'HLA-DPA1*02:01-HLA-DPB1*09:01',
         'HLA-DPA1*02:01-HLA-DPB1*10:001', 'HLA-DPA1*02:01-HLA-DPB1*10:01',
         'HLA-DPA1*02:01-HLA-DPB1*10:101', 'HLA-DPA1*02:01-HLA-DPB1*10:201', 'HLA-DPA1*02:01-HLA-DPB1*10:301', 'HLA-DPA1*02:01-HLA-DPB1*10:401',
         'HLA-DPA1*02:01-HLA-DPB1*10:501', 'HLA-DPA1*02:01-HLA-DPB1*10:601',
         'HLA-DPA1*02:01-HLA-DPB1*10:701', 'HLA-DPA1*02:01-HLA-DPB1*10:801', 'HLA-DPA1*02:01-HLA-DPB1*10:901', 'HLA-DPA1*02:01-HLA-DPB1*11:001',
         'HLA-DPA1*02:01-HLA-DPB1*11:01', 'HLA-DPA1*02:01-HLA-DPB1*11:101',
         'HLA-DPA1*02:01-HLA-DPB1*11:201', 'HLA-DPA1*02:01-HLA-DPB1*11:301', 'HLA-DPA1*02:01-HLA-DPB1*11:401', 'HLA-DPA1*02:01-HLA-DPB1*11:501',
         'HLA-DPA1*02:01-HLA-DPB1*11:601', 'HLA-DPA1*02:01-HLA-DPB1*11:701',
         'HLA-DPA1*02:01-HLA-DPB1*11:801', 'HLA-DPA1*02:01-HLA-DPB1*11:901', 'HLA-DPA1*02:01-HLA-DPB1*12:101', 'HLA-DPA1*02:01-HLA-DPB1*12:201',
         'HLA-DPA1*02:01-HLA-DPB1*12:301', 'HLA-DPA1*02:01-HLA-DPB1*12:401',
         'HLA-DPA1*02:01-HLA-DPB1*12:501', 'HLA-DPA1*02:01-HLA-DPB1*12:601', 'HLA-DPA1*02:01-HLA-DPB1*12:701', 'HLA-DPA1*02:01-HLA-DPB1*12:801',
         'HLA-DPA1*02:01-HLA-DPB1*12:901', 'HLA-DPA1*02:01-HLA-DPB1*13:001',
         'HLA-DPA1*02:01-HLA-DPB1*13:01', 'HLA-DPA1*02:01-HLA-DPB1*13:101', 'HLA-DPA1*02:01-HLA-DPB1*13:201', 'HLA-DPA1*02:01-HLA-DPB1*13:301',
         'HLA-DPA1*02:01-HLA-DPB1*13:401', 'HLA-DPA1*02:01-HLA-DPB1*14:01',
         'HLA-DPA1*02:01-HLA-DPB1*15:01', 'HLA-DPA1*02:01-HLA-DPB1*16:01', 'HLA-DPA1*02:01-HLA-DPB1*17:01', 'HLA-DPA1*02:01-HLA-DPB1*18:01',
         'HLA-DPA1*02:01-HLA-DPB1*19:01', 'HLA-DPA1*02:01-HLA-DPB1*20:01',
         'HLA-DPA1*02:01-HLA-DPB1*21:01', 'HLA-DPA1*02:01-HLA-DPB1*22:01', 'HLA-DPA1*02:01-HLA-DPB1*23:01', 'HLA-DPA1*02:01-HLA-DPB1*24:01',
         'HLA-DPA1*02:01-HLA-DPB1*25:01', 'HLA-DPA1*02:01-HLA-DPB1*26:01',
         'HLA-DPA1*02:01-HLA-DPB1*27:01', 'HLA-DPA1*02:01-HLA-DPB1*28:01', 'HLA-DPA1*02:01-HLA-DPB1*29:01', 'HLA-DPA1*02:01-HLA-DPB1*30:01',
         'HLA-DPA1*02:01-HLA-DPB1*31:01', 'HLA-DPA1*02:01-HLA-DPB1*32:01',
         'HLA-DPA1*02:01-HLA-DPB1*33:01', 'HLA-DPA1*02:01-HLA-DPB1*34:01', 'HLA-DPA1*02:01-HLA-DPB1*35:01', 'HLA-DPA1*02:01-HLA-DPB1*36:01',
         'HLA-DPA1*02:01-HLA-DPB1*37:01', 'HLA-DPA1*02:01-HLA-DPB1*38:01',
         'HLA-DPA1*02:01-HLA-DPB1*39:01', 'HLA-DPA1*02:01-HLA-DPB1*40:01', 'HLA-DPA1*02:01-HLA-DPB1*41:01', 'HLA-DPA1*02:01-HLA-DPB1*44:01',
         'HLA-DPA1*02:01-HLA-DPB1*45:01', 'HLA-DPA1*02:01-HLA-DPB1*46:01',
         'HLA-DPA1*02:01-HLA-DPB1*47:01', 'HLA-DPA1*02:01-HLA-DPB1*48:01', 'HLA-DPA1*02:01-HLA-DPB1*49:01', 'HLA-DPA1*02:01-HLA-DPB1*50:01',
         'HLA-DPA1*02:01-HLA-DPB1*51:01', 'HLA-DPA1*02:01-HLA-DPB1*52:01',
         'HLA-DPA1*02:01-HLA-DPB1*53:01', 'HLA-DPA1*02:01-HLA-DPB1*54:01', 'HLA-DPA1*02:01-HLA-DPB1*55:01', 'HLA-DPA1*02:01-HLA-DPB1*56:01',
         'HLA-DPA1*02:01-HLA-DPB1*58:01', 'HLA-DPA1*02:01-HLA-DPB1*59:01',
         'HLA-DPA1*02:01-HLA-DPB1*60:01', 'HLA-DPA1*02:01-HLA-DPB1*62:01', 'HLA-DPA1*02:01-HLA-DPB1*63:01', 'HLA-DPA1*02:01-HLA-DPB1*65:01',
         'HLA-DPA1*02:01-HLA-DPB1*66:01', 'HLA-DPA1*02:01-HLA-DPB1*67:01',
         'HLA-DPA1*02:01-HLA-DPB1*68:01', 'HLA-DPA1*02:01-HLA-DPB1*69:01', 'HLA-DPA1*02:01-HLA-DPB1*70:01', 'HLA-DPA1*02:01-HLA-DPB1*71:01',
         'HLA-DPA1*02:01-HLA-DPB1*72:01', 'HLA-DPA1*02:01-HLA-DPB1*73:01',
         'HLA-DPA1*02:01-HLA-DPB1*74:01', 'HLA-DPA1*02:01-HLA-DPB1*75:01', 'HLA-DPA1*02:01-HLA-DPB1*76:01', 'HLA-DPA1*02:01-HLA-DPB1*77:01',
         'HLA-DPA1*02:01-HLA-DPB1*78:01', 'HLA-DPA1*02:01-HLA-DPB1*79:01',
         'HLA-DPA1*02:01-HLA-DPB1*80:01', 'HLA-DPA1*02:01-HLA-DPB1*81:01', 'HLA-DPA1*02:01-HLA-DPB1*82:01', 'HLA-DPA1*02:01-HLA-DPB1*83:01',
         'HLA-DPA1*02:01-HLA-DPB1*84:01', 'HLA-DPA1*02:01-HLA-DPB1*85:01',
         'HLA-DPA1*02:01-HLA-DPB1*86:01', 'HLA-DPA1*02:01-HLA-DPB1*87:01', 'HLA-DPA1*02:01-HLA-DPB1*88:01', 'HLA-DPA1*02:01-HLA-DPB1*89:01',
         'HLA-DPA1*02:01-HLA-DPB1*90:01', 'HLA-DPA1*02:01-HLA-DPB1*91:01',
         'HLA-DPA1*02:01-HLA-DPB1*92:01', 'HLA-DPA1*02:01-HLA-DPB1*93:01', 'HLA-DPA1*02:01-HLA-DPB1*94:01', 'HLA-DPA1*02:01-HLA-DPB1*95:01',
         'HLA-DPA1*02:01-HLA-DPB1*96:01', 'HLA-DPA1*02:01-HLA-DPB1*97:01',
         'HLA-DPA1*02:01-HLA-DPB1*98:01', 'HLA-DPA1*02:01-HLA-DPB1*99:01', 'HLA-DPA1*02:02-HLA-DPB1*01:01', 'HLA-DPA1*02:02-HLA-DPB1*02:01',
         'HLA-DPA1*02:02-HLA-DPB1*02:02', 'HLA-DPA1*02:02-HLA-DPB1*03:01',
         'HLA-DPA1*02:02-HLA-DPB1*04:01', 'HLA-DPA1*02:02-HLA-DPB1*04:02', 'HLA-DPA1*02:02-HLA-DPB1*05:01', 'HLA-DPA1*02:02-HLA-DPB1*06:01',
         'HLA-DPA1*02:02-HLA-DPB1*08:01', 'HLA-DPA1*02:02-HLA-DPB1*09:01',
         'HLA-DPA1*02:02-HLA-DPB1*10:001', 'HLA-DPA1*02:02-HLA-DPB1*10:01', 'HLA-DPA1*02:02-HLA-DPB1*10:101', 'HLA-DPA1*02:02-HLA-DPB1*10:201',
         'HLA-DPA1*02:02-HLA-DPB1*10:301', 'HLA-DPA1*02:02-HLA-DPB1*10:401',
         'HLA-DPA1*02:02-HLA-DPB1*10:501', 'HLA-DPA1*02:02-HLA-DPB1*10:601', 'HLA-DPA1*02:02-HLA-DPB1*10:701', 'HLA-DPA1*02:02-HLA-DPB1*10:801',
         'HLA-DPA1*02:02-HLA-DPB1*10:901', 'HLA-DPA1*02:02-HLA-DPB1*11:001',
         'HLA-DPA1*02:02-HLA-DPB1*11:01', 'HLA-DPA1*02:02-HLA-DPB1*11:101', 'HLA-DPA1*02:02-HLA-DPB1*11:201', 'HLA-DPA1*02:02-HLA-DPB1*11:301',
         'HLA-DPA1*02:02-HLA-DPB1*11:401', 'HLA-DPA1*02:02-HLA-DPB1*11:501',
         'HLA-DPA1*02:02-HLA-DPB1*11:601', 'HLA-DPA1*02:02-HLA-DPB1*11:701', 'HLA-DPA1*02:02-HLA-DPB1*11:801', 'HLA-DPA1*02:02-HLA-DPB1*11:901',
         'HLA-DPA1*02:02-HLA-DPB1*12:101', 'HLA-DPA1*02:02-HLA-DPB1*12:201',
         'HLA-DPA1*02:02-HLA-DPB1*12:301', 'HLA-DPA1*02:02-HLA-DPB1*12:401', 'HLA-DPA1*02:02-HLA-DPB1*12:501', 'HLA-DPA1*02:02-HLA-DPB1*12:601',
         'HLA-DPA1*02:02-HLA-DPB1*12:701', 'HLA-DPA1*02:02-HLA-DPB1*12:801',
         'HLA-DPA1*02:02-HLA-DPB1*12:901', 'HLA-DPA1*02:02-HLA-DPB1*13:001', 'HLA-DPA1*02:02-HLA-DPB1*13:01', 'HLA-DPA1*02:02-HLA-DPB1*13:101',
         'HLA-DPA1*02:02-HLA-DPB1*13:201', 'HLA-DPA1*02:02-HLA-DPB1*13:301',
         'HLA-DPA1*02:02-HLA-DPB1*13:401', 'HLA-DPA1*02:02-HLA-DPB1*14:01', 'HLA-DPA1*02:02-HLA-DPB1*15:01', 'HLA-DPA1*02:02-HLA-DPB1*16:01',
         'HLA-DPA1*02:02-HLA-DPB1*17:01', 'HLA-DPA1*02:02-HLA-DPB1*18:01',
         'HLA-DPA1*02:02-HLA-DPB1*19:01', 'HLA-DPA1*02:02-HLA-DPB1*20:01', 'HLA-DPA1*02:02-HLA-DPB1*21:01', 'HLA-DPA1*02:02-HLA-DPB1*22:01',
         'HLA-DPA1*02:02-HLA-DPB1*23:01', 'HLA-DPA1*02:02-HLA-DPB1*24:01',
         'HLA-DPA1*02:02-HLA-DPB1*25:01', 'HLA-DPA1*02:02-HLA-DPB1*26:01', 'HLA-DPA1*02:02-HLA-DPB1*27:01', 'HLA-DPA1*02:02-HLA-DPB1*28:01',
         'HLA-DPA1*02:02-HLA-DPB1*29:01', 'HLA-DPA1*02:02-HLA-DPB1*30:01',
         'HLA-DPA1*02:02-HLA-DPB1*31:01', 'HLA-DPA1*02:02-HLA-DPB1*32:01', 'HLA-DPA1*02:02-HLA-DPB1*33:01', 'HLA-DPA1*02:02-HLA-DPB1*34:01',
         'HLA-DPA1*02:02-HLA-DPB1*35:01', 'HLA-DPA1*02:02-HLA-DPB1*36:01',
         'HLA-DPA1*02:02-HLA-DPB1*37:01', 'HLA-DPA1*02:02-HLA-DPB1*38:01', 'HLA-DPA1*02:02-HLA-DPB1*39:01', 'HLA-DPA1*02:02-HLA-DPB1*40:01',
         'HLA-DPA1*02:02-HLA-DPB1*41:01', 'HLA-DPA1*02:02-HLA-DPB1*44:01',
         'HLA-DPA1*02:02-HLA-DPB1*45:01', 'HLA-DPA1*02:02-HLA-DPB1*46:01', 'HLA-DPA1*02:02-HLA-DPB1*47:01', 'HLA-DPA1*02:02-HLA-DPB1*48:01',
         'HLA-DPA1*02:02-HLA-DPB1*49:01', 'HLA-DPA1*02:02-HLA-DPB1*50:01',
         'HLA-DPA1*02:02-HLA-DPB1*51:01', 'HLA-DPA1*02:02-HLA-DPB1*52:01', 'HLA-DPA1*02:02-HLA-DPB1*53:01', 'HLA-DPA1*02:02-HLA-DPB1*54:01',
         'HLA-DPA1*02:02-HLA-DPB1*55:01', 'HLA-DPA1*02:02-HLA-DPB1*56:01',
         'HLA-DPA1*02:02-HLA-DPB1*58:01', 'HLA-DPA1*02:02-HLA-DPB1*59:01', 'HLA-DPA1*02:02-HLA-DPB1*60:01', 'HLA-DPA1*02:02-HLA-DPB1*62:01',
         'HLA-DPA1*02:02-HLA-DPB1*63:01', 'HLA-DPA1*02:02-HLA-DPB1*65:01',
         'HLA-DPA1*02:02-HLA-DPB1*66:01', 'HLA-DPA1*02:02-HLA-DPB1*67:01', 'HLA-DPA1*02:02-HLA-DPB1*68:01', 'HLA-DPA1*02:02-HLA-DPB1*69:01',
         'HLA-DPA1*02:02-HLA-DPB1*70:01', 'HLA-DPA1*02:02-HLA-DPB1*71:01',
         'HLA-DPA1*02:02-HLA-DPB1*72:01', 'HLA-DPA1*02:02-HLA-DPB1*73:01', 'HLA-DPA1*02:02-HLA-DPB1*74:01', 'HLA-DPA1*02:02-HLA-DPB1*75:01',
         'HLA-DPA1*02:02-HLA-DPB1*76:01', 'HLA-DPA1*02:02-HLA-DPB1*77:01',
         'HLA-DPA1*02:02-HLA-DPB1*78:01', 'HLA-DPA1*02:02-HLA-DPB1*79:01', 'HLA-DPA1*02:02-HLA-DPB1*80:01', 'HLA-DPA1*02:02-HLA-DPB1*81:01',
         'HLA-DPA1*02:02-HLA-DPB1*82:01', 'HLA-DPA1*02:02-HLA-DPB1*83:01',
         'HLA-DPA1*02:02-HLA-DPB1*84:01', 'HLA-DPA1*02:02-HLA-DPB1*85:01', 'HLA-DPA1*02:02-HLA-DPB1*86:01', 'HLA-DPA1*02:02-HLA-DPB1*87:01',
         'HLA-DPA1*02:02-HLA-DPB1*88:01', 'HLA-DPA1*02:02-HLA-DPB1*89:01',
         'HLA-DPA1*02:02-HLA-DPB1*90:01', 'HLA-DPA1*02:02-HLA-DPB1*91:01', 'HLA-DPA1*02:02-HLA-DPB1*92:01', 'HLA-DPA1*02:02-HLA-DPB1*93:01',
         'HLA-DPA1*02:02-HLA-DPB1*94:01', 'HLA-DPA1*02:02-HLA-DPB1*95:01',
         'HLA-DPA1*02:02-HLA-DPB1*96:01', 'HLA-DPA1*02:02-HLA-DPB1*97:01', 'HLA-DPA1*02:02-HLA-DPB1*98:01', 'HLA-DPA1*02:02-HLA-DPB1*99:01',
         'HLA-DPA1*02:03-HLA-DPB1*01:01', 'HLA-DPA1*02:03-HLA-DPB1*02:01',
         'HLA-DPA1*02:03-HLA-DPB1*02:02', 'HLA-DPA1*02:03-HLA-DPB1*03:01', 'HLA-DPA1*02:03-HLA-DPB1*04:01', 'HLA-DPA1*02:03-HLA-DPB1*04:02',
         'HLA-DPA1*02:03-HLA-DPB1*05:01', 'HLA-DPA1*02:03-HLA-DPB1*06:01',
         'HLA-DPA1*02:03-HLA-DPB1*08:01', 'HLA-DPA1*02:03-HLA-DPB1*09:01', 'HLA-DPA1*02:03-HLA-DPB1*10:001', 'HLA-DPA1*02:03-HLA-DPB1*10:01',
         'HLA-DPA1*02:03-HLA-DPB1*10:101', 'HLA-DPA1*02:03-HLA-DPB1*10:201',
         'HLA-DPA1*02:03-HLA-DPB1*10:301', 'HLA-DPA1*02:03-HLA-DPB1*10:401', 'HLA-DPA1*02:03-HLA-DPB1*10:501', 'HLA-DPA1*02:03-HLA-DPB1*10:601',
         'HLA-DPA1*02:03-HLA-DPB1*10:701', 'HLA-DPA1*02:03-HLA-DPB1*10:801',
         'HLA-DPA1*02:03-HLA-DPB1*10:901', 'HLA-DPA1*02:03-HLA-DPB1*11:001', 'HLA-DPA1*02:03-HLA-DPB1*11:01', 'HLA-DPA1*02:03-HLA-DPB1*11:101',
         'HLA-DPA1*02:03-HLA-DPB1*11:201', 'HLA-DPA1*02:03-HLA-DPB1*11:301',
         'HLA-DPA1*02:03-HLA-DPB1*11:401', 'HLA-DPA1*02:03-HLA-DPB1*11:501', 'HLA-DPA1*02:03-HLA-DPB1*11:601', 'HLA-DPA1*02:03-HLA-DPB1*11:701',
         'HLA-DPA1*02:03-HLA-DPB1*11:801', 'HLA-DPA1*02:03-HLA-DPB1*11:901',
         'HLA-DPA1*02:03-HLA-DPB1*12:101', 'HLA-DPA1*02:03-HLA-DPB1*12:201', 'HLA-DPA1*02:03-HLA-DPB1*12:301', 'HLA-DPA1*02:03-HLA-DPB1*12:401',
         'HLA-DPA1*02:03-HLA-DPB1*12:501', 'HLA-DPA1*02:03-HLA-DPB1*12:601',
         'HLA-DPA1*02:03-HLA-DPB1*12:701', 'HLA-DPA1*02:03-HLA-DPB1*12:801', 'HLA-DPA1*02:03-HLA-DPB1*12:901', 'HLA-DPA1*02:03-HLA-DPB1*13:001',
         'HLA-DPA1*02:03-HLA-DPB1*13:01', 'HLA-DPA1*02:03-HLA-DPB1*13:101',
         'HLA-DPA1*02:03-HLA-DPB1*13:201', 'HLA-DPA1*02:03-HLA-DPB1*13:301', 'HLA-DPA1*02:03-HLA-DPB1*13:401', 'HLA-DPA1*02:03-HLA-DPB1*14:01',
         'HLA-DPA1*02:03-HLA-DPB1*15:01', 'HLA-DPA1*02:03-HLA-DPB1*16:01',
         'HLA-DPA1*02:03-HLA-DPB1*17:01', 'HLA-DPA1*02:03-HLA-DPB1*18:01', 'HLA-DPA1*02:03-HLA-DPB1*19:01', 'HLA-DPA1*02:03-HLA-DPB1*20:01',
         'HLA-DPA1*02:03-HLA-DPB1*21:01', 'HLA-DPA1*02:03-HLA-DPB1*22:01',
         'HLA-DPA1*02:03-HLA-DPB1*23:01', 'HLA-DPA1*02:03-HLA-DPB1*24:01', 'HLA-DPA1*02:03-HLA-DPB1*25:01', 'HLA-DPA1*02:03-HLA-DPB1*26:01',
         'HLA-DPA1*02:03-HLA-DPB1*27:01', 'HLA-DPA1*02:03-HLA-DPB1*28:01',
         'HLA-DPA1*02:03-HLA-DPB1*29:01', 'HLA-DPA1*02:03-HLA-DPB1*30:01', 'HLA-DPA1*02:03-HLA-DPB1*31:01', 'HLA-DPA1*02:03-HLA-DPB1*32:01',
         'HLA-DPA1*02:03-HLA-DPB1*33:01', 'HLA-DPA1*02:03-HLA-DPB1*34:01',
         'HLA-DPA1*02:03-HLA-DPB1*35:01', 'HLA-DPA1*02:03-HLA-DPB1*36:01', 'HLA-DPA1*02:03-HLA-DPB1*37:01', 'HLA-DPA1*02:03-HLA-DPB1*38:01',
         'HLA-DPA1*02:03-HLA-DPB1*39:01', 'HLA-DPA1*02:03-HLA-DPB1*40:01',
         'HLA-DPA1*02:03-HLA-DPB1*41:01', 'HLA-DPA1*02:03-HLA-DPB1*44:01', 'HLA-DPA1*02:03-HLA-DPB1*45:01', 'HLA-DPA1*02:03-HLA-DPB1*46:01',
         'HLA-DPA1*02:03-HLA-DPB1*47:01', 'HLA-DPA1*02:03-HLA-DPB1*48:01',
         'HLA-DPA1*02:03-HLA-DPB1*49:01', 'HLA-DPA1*02:03-HLA-DPB1*50:01', 'HLA-DPA1*02:03-HLA-DPB1*51:01', 'HLA-DPA1*02:03-HLA-DPB1*52:01',
         'HLA-DPA1*02:03-HLA-DPB1*53:01', 'HLA-DPA1*02:03-HLA-DPB1*54:01',
         'HLA-DPA1*02:03-HLA-DPB1*55:01', 'HLA-DPA1*02:03-HLA-DPB1*56:01', 'HLA-DPA1*02:03-HLA-DPB1*58:01', 'HLA-DPA1*02:03-HLA-DPB1*59:01',
         'HLA-DPA1*02:03-HLA-DPB1*60:01', 'HLA-DPA1*02:03-HLA-DPB1*62:01',
         'HLA-DPA1*02:03-HLA-DPB1*63:01', 'HLA-DPA1*02:03-HLA-DPB1*65:01', 'HLA-DPA1*02:03-HLA-DPB1*66:01', 'HLA-DPA1*02:03-HLA-DPB1*67:01',
         'HLA-DPA1*02:03-HLA-DPB1*68:01', 'HLA-DPA1*02:03-HLA-DPB1*69:01',
         'HLA-DPA1*02:03-HLA-DPB1*70:01', 'HLA-DPA1*02:03-HLA-DPB1*71:01', 'HLA-DPA1*02:03-HLA-DPB1*72:01', 'HLA-DPA1*02:03-HLA-DPB1*73:01',
         'HLA-DPA1*02:03-HLA-DPB1*74:01', 'HLA-DPA1*02:03-HLA-DPB1*75:01',
         'HLA-DPA1*02:03-HLA-DPB1*76:01', 'HLA-DPA1*02:03-HLA-DPB1*77:01', 'HLA-DPA1*02:03-HLA-DPB1*78:01', 'HLA-DPA1*02:03-HLA-DPB1*79:01',
         'HLA-DPA1*02:03-HLA-DPB1*80:01', 'HLA-DPA1*02:03-HLA-DPB1*81:01',
         'HLA-DPA1*02:03-HLA-DPB1*82:01', 'HLA-DPA1*02:03-HLA-DPB1*83:01', 'HLA-DPA1*02:03-HLA-DPB1*84:01', 'HLA-DPA1*02:03-HLA-DPB1*85:01',
         'HLA-DPA1*02:03-HLA-DPB1*86:01', 'HLA-DPA1*02:03-HLA-DPB1*87:01',
         'HLA-DPA1*02:03-HLA-DPB1*88:01', 'HLA-DPA1*02:03-HLA-DPB1*89:01', 'HLA-DPA1*02:03-HLA-DPB1*90:01', 'HLA-DPA1*02:03-HLA-DPB1*91:01',
         'HLA-DPA1*02:03-HLA-DPB1*92:01', 'HLA-DPA1*02:03-HLA-DPB1*93:01',
         'HLA-DPA1*02:03-HLA-DPB1*94:01', 'HLA-DPA1*02:03-HLA-DPB1*95:01', 'HLA-DPA1*02:03-HLA-DPB1*96:01', 'HLA-DPA1*02:03-HLA-DPB1*97:01',
         'HLA-DPA1*02:03-HLA-DPB1*98:01', 'HLA-DPA1*02:03-HLA-DPB1*99:01',
         'HLA-DPA1*02:04-HLA-DPB1*01:01', 'HLA-DPA1*02:04-HLA-DPB1*02:01', 'HLA-DPA1*02:04-HLA-DPB1*02:02', 'HLA-DPA1*02:04-HLA-DPB1*03:01',
         'HLA-DPA1*02:04-HLA-DPB1*04:01', 'HLA-DPA1*02:04-HLA-DPB1*04:02',
         'HLA-DPA1*02:04-HLA-DPB1*05:01', 'HLA-DPA1*02:04-HLA-DPB1*06:01', 'HLA-DPA1*02:04-HLA-DPB1*08:01', 'HLA-DPA1*02:04-HLA-DPB1*09:01',
         'HLA-DPA1*02:04-HLA-DPB1*10:001', 'HLA-DPA1*02:04-HLA-DPB1*10:01',
         'HLA-DPA1*02:04-HLA-DPB1*10:101', 'HLA-DPA1*02:04-HLA-DPB1*10:201', 'HLA-DPA1*02:04-HLA-DPB1*10:301', 'HLA-DPA1*02:04-HLA-DPB1*10:401',
         'HLA-DPA1*02:04-HLA-DPB1*10:501', 'HLA-DPA1*02:04-HLA-DPB1*10:601',
         'HLA-DPA1*02:04-HLA-DPB1*10:701', 'HLA-DPA1*02:04-HLA-DPB1*10:801', 'HLA-DPA1*02:04-HLA-DPB1*10:901', 'HLA-DPA1*02:04-HLA-DPB1*11:001',
         'HLA-DPA1*02:04-HLA-DPB1*11:01', 'HLA-DPA1*02:04-HLA-DPB1*11:101',
         'HLA-DPA1*02:04-HLA-DPB1*11:201', 'HLA-DPA1*02:04-HLA-DPB1*11:301', 'HLA-DPA1*02:04-HLA-DPB1*11:401', 'HLA-DPA1*02:04-HLA-DPB1*11:501',
         'HLA-DPA1*02:04-HLA-DPB1*11:601', 'HLA-DPA1*02:04-HLA-DPB1*11:701',
         'HLA-DPA1*02:04-HLA-DPB1*11:801', 'HLA-DPA1*02:04-HLA-DPB1*11:901', 'HLA-DPA1*02:04-HLA-DPB1*12:101', 'HLA-DPA1*02:04-HLA-DPB1*12:201',
         'HLA-DPA1*02:04-HLA-DPB1*12:301', 'HLA-DPA1*02:04-HLA-DPB1*12:401',
         'HLA-DPA1*02:04-HLA-DPB1*12:501', 'HLA-DPA1*02:04-HLA-DPB1*12:601', 'HLA-DPA1*02:04-HLA-DPB1*12:701', 'HLA-DPA1*02:04-HLA-DPB1*12:801',
         'HLA-DPA1*02:04-HLA-DPB1*12:901', 'HLA-DPA1*02:04-HLA-DPB1*13:001',
         'HLA-DPA1*02:04-HLA-DPB1*13:01', 'HLA-DPA1*02:04-HLA-DPB1*13:101', 'HLA-DPA1*02:04-HLA-DPB1*13:201', 'HLA-DPA1*02:04-HLA-DPB1*13:301',
         'HLA-DPA1*02:04-HLA-DPB1*13:401', 'HLA-DPA1*02:04-HLA-DPB1*14:01',
         'HLA-DPA1*02:04-HLA-DPB1*15:01', 'HLA-DPA1*02:04-HLA-DPB1*16:01', 'HLA-DPA1*02:04-HLA-DPB1*17:01', 'HLA-DPA1*02:04-HLA-DPB1*18:01',
         'HLA-DPA1*02:04-HLA-DPB1*19:01', 'HLA-DPA1*02:04-HLA-DPB1*20:01',
         'HLA-DPA1*02:04-HLA-DPB1*21:01', 'HLA-DPA1*02:04-HLA-DPB1*22:01', 'HLA-DPA1*02:04-HLA-DPB1*23:01', 'HLA-DPA1*02:04-HLA-DPB1*24:01',
         'HLA-DPA1*02:04-HLA-DPB1*25:01', 'HLA-DPA1*02:04-HLA-DPB1*26:01',
         'HLA-DPA1*02:04-HLA-DPB1*27:01', 'HLA-DPA1*02:04-HLA-DPB1*28:01', 'HLA-DPA1*02:04-HLA-DPB1*29:01', 'HLA-DPA1*02:04-HLA-DPB1*30:01',
         'HLA-DPA1*02:04-HLA-DPB1*31:01', 'HLA-DPA1*02:04-HLA-DPB1*32:01',
         'HLA-DPA1*02:04-HLA-DPB1*33:01', 'HLA-DPA1*02:04-HLA-DPB1*34:01', 'HLA-DPA1*02:04-HLA-DPB1*35:01', 'HLA-DPA1*02:04-HLA-DPB1*36:01',
         'HLA-DPA1*02:04-HLA-DPB1*37:01', 'HLA-DPA1*02:04-HLA-DPB1*38:01',
         'HLA-DPA1*02:04-HLA-DPB1*39:01', 'HLA-DPA1*02:04-HLA-DPB1*40:01', 'HLA-DPA1*02:04-HLA-DPB1*41:01', 'HLA-DPA1*02:04-HLA-DPB1*44:01',
         'HLA-DPA1*02:04-HLA-DPB1*45:01', 'HLA-DPA1*02:04-HLA-DPB1*46:01',
         'HLA-DPA1*02:04-HLA-DPB1*47:01', 'HLA-DPA1*02:04-HLA-DPB1*48:01', 'HLA-DPA1*02:04-HLA-DPB1*49:01', 'HLA-DPA1*02:04-HLA-DPB1*50:01',
         'HLA-DPA1*02:04-HLA-DPB1*51:01', 'HLA-DPA1*02:04-HLA-DPB1*52:01',
         'HLA-DPA1*02:04-HLA-DPB1*53:01', 'HLA-DPA1*02:04-HLA-DPB1*54:01', 'HLA-DPA1*02:04-HLA-DPB1*55:01', 'HLA-DPA1*02:04-HLA-DPB1*56:01',
         'HLA-DPA1*02:04-HLA-DPB1*58:01', 'HLA-DPA1*02:04-HLA-DPB1*59:01',
         'HLA-DPA1*02:04-HLA-DPB1*60:01', 'HLA-DPA1*02:04-HLA-DPB1*62:01', 'HLA-DPA1*02:04-HLA-DPB1*63:01', 'HLA-DPA1*02:04-HLA-DPB1*65:01',
         'HLA-DPA1*02:04-HLA-DPB1*66:01', 'HLA-DPA1*02:04-HLA-DPB1*67:01',
         'HLA-DPA1*02:04-HLA-DPB1*68:01', 'HLA-DPA1*02:04-HLA-DPB1*69:01', 'HLA-DPA1*02:04-HLA-DPB1*70:01', 'HLA-DPA1*02:04-HLA-DPB1*71:01',
         'HLA-DPA1*02:04-HLA-DPB1*72:01', 'HLA-DPA1*02:04-HLA-DPB1*73:01',
         'HLA-DPA1*02:04-HLA-DPB1*74:01', 'HLA-DPA1*02:04-HLA-DPB1*75:01', 'HLA-DPA1*02:04-HLA-DPB1*76:01', 'HLA-DPA1*02:04-HLA-DPB1*77:01',
         'HLA-DPA1*02:04-HLA-DPB1*78:01', 'HLA-DPA1*02:04-HLA-DPB1*79:01',
         'HLA-DPA1*02:04-HLA-DPB1*80:01', 'HLA-DPA1*02:04-HLA-DPB1*81:01', 'HLA-DPA1*02:04-HLA-DPB1*82:01', 'HLA-DPA1*02:04-HLA-DPB1*83:01',
         'HLA-DPA1*02:04-HLA-DPB1*84:01', 'HLA-DPA1*02:04-HLA-DPB1*85:01',
         'HLA-DPA1*02:04-HLA-DPB1*86:01', 'HLA-DPA1*02:04-HLA-DPB1*87:01', 'HLA-DPA1*02:04-HLA-DPB1*88:01', 'HLA-DPA1*02:04-HLA-DPB1*89:01',
         'HLA-DPA1*02:04-HLA-DPB1*90:01', 'HLA-DPA1*02:04-HLA-DPB1*91:01',
         'HLA-DPA1*02:04-HLA-DPB1*92:01', 'HLA-DPA1*02:04-HLA-DPB1*93:01', 'HLA-DPA1*02:04-HLA-DPB1*94:01', 'HLA-DPA1*02:04-HLA-DPB1*95:01',
         'HLA-DPA1*02:04-HLA-DPB1*96:01', 'HLA-DPA1*02:04-HLA-DPB1*97:01',
         'HLA-DPA1*02:04-HLA-DPB1*98:01', 'HLA-DPA1*02:04-HLA-DPB1*99:01', 'HLA-DPA1*03:01-HLA-DPB1*01:01', 'HLA-DPA1*03:01-HLA-DPB1*02:01',
         'HLA-DPA1*03:01-HLA-DPB1*02:02', 'HLA-DPA1*03:01-HLA-DPB1*03:01',
         'HLA-DPA1*03:01-HLA-DPB1*04:01', 'HLA-DPA1*03:01-HLA-DPB1*04:02', 'HLA-DPA1*03:01-HLA-DPB1*05:01', 'HLA-DPA1*03:01-HLA-DPB1*06:01',
         'HLA-DPA1*03:01-HLA-DPB1*08:01', 'HLA-DPA1*03:01-HLA-DPB1*09:01',
         'HLA-DPA1*03:01-HLA-DPB1*10:001', 'HLA-DPA1*03:01-HLA-DPB1*10:01', 'HLA-DPA1*03:01-HLA-DPB1*10:101', 'HLA-DPA1*03:01-HLA-DPB1*10:201',
         'HLA-DPA1*03:01-HLA-DPB1*10:301', 'HLA-DPA1*03:01-HLA-DPB1*10:401',
         'HLA-DPA1*03:01-HLA-DPB1*10:501', 'HLA-DPA1*03:01-HLA-DPB1*10:601', 'HLA-DPA1*03:01-HLA-DPB1*10:701', 'HLA-DPA1*03:01-HLA-DPB1*10:801',
         'HLA-DPA1*03:01-HLA-DPB1*10:901', 'HLA-DPA1*03:01-HLA-DPB1*11:001',
         'HLA-DPA1*03:01-HLA-DPB1*11:01', 'HLA-DPA1*03:01-HLA-DPB1*11:101', 'HLA-DPA1*03:01-HLA-DPB1*11:201', 'HLA-DPA1*03:01-HLA-DPB1*11:301',
         'HLA-DPA1*03:01-HLA-DPB1*11:401', 'HLA-DPA1*03:01-HLA-DPB1*11:501',
         'HLA-DPA1*03:01-HLA-DPB1*11:601', 'HLA-DPA1*03:01-HLA-DPB1*11:701', 'HLA-DPA1*03:01-HLA-DPB1*11:801', 'HLA-DPA1*03:01-HLA-DPB1*11:901',
         'HLA-DPA1*03:01-HLA-DPB1*12:101', 'HLA-DPA1*03:01-HLA-DPB1*12:201',
         'HLA-DPA1*03:01-HLA-DPB1*12:301', 'HLA-DPA1*03:01-HLA-DPB1*12:401', 'HLA-DPA1*03:01-HLA-DPB1*12:501', 'HLA-DPA1*03:01-HLA-DPB1*12:601',
         'HLA-DPA1*03:01-HLA-DPB1*12:701', 'HLA-DPA1*03:01-HLA-DPB1*12:801',
         'HLA-DPA1*03:01-HLA-DPB1*12:901', 'HLA-DPA1*03:01-HLA-DPB1*13:001', 'HLA-DPA1*03:01-HLA-DPB1*13:01', 'HLA-DPA1*03:01-HLA-DPB1*13:101',
         'HLA-DPA1*03:01-HLA-DPB1*13:201', 'HLA-DPA1*03:01-HLA-DPB1*13:301',
         'HLA-DPA1*03:01-HLA-DPB1*13:401', 'HLA-DPA1*03:01-HLA-DPB1*14:01', 'HLA-DPA1*03:01-HLA-DPB1*15:01', 'HLA-DPA1*03:01-HLA-DPB1*16:01',
         'HLA-DPA1*03:01-HLA-DPB1*17:01', 'HLA-DPA1*03:01-HLA-DPB1*18:01',
         'HLA-DPA1*03:01-HLA-DPB1*19:01', 'HLA-DPA1*03:01-HLA-DPB1*20:01', 'HLA-DPA1*03:01-HLA-DPB1*21:01', 'HLA-DPA1*03:01-HLA-DPB1*22:01',
         'HLA-DPA1*03:01-HLA-DPB1*23:01', 'HLA-DPA1*03:01-HLA-DPB1*24:01',
         'HLA-DPA1*03:01-HLA-DPB1*25:01', 'HLA-DPA1*03:01-HLA-DPB1*26:01', 'HLA-DPA1*03:01-HLA-DPB1*27:01', 'HLA-DPA1*03:01-HLA-DPB1*28:01',
         'HLA-DPA1*03:01-HLA-DPB1*29:01', 'HLA-DPA1*03:01-HLA-DPB1*30:01',
         'HLA-DPA1*03:01-HLA-DPB1*31:01', 'HLA-DPA1*03:01-HLA-DPB1*32:01', 'HLA-DPA1*03:01-HLA-DPB1*33:01', 'HLA-DPA1*03:01-HLA-DPB1*34:01',
         'HLA-DPA1*03:01-HLA-DPB1*35:01', 'HLA-DPA1*03:01-HLA-DPB1*36:01',
         'HLA-DPA1*03:01-HLA-DPB1*37:01', 'HLA-DPA1*03:01-HLA-DPB1*38:01', 'HLA-DPA1*03:01-HLA-DPB1*39:01', 'HLA-DPA1*03:01-HLA-DPB1*40:01',
         'HLA-DPA1*03:01-HLA-DPB1*41:01', 'HLA-DPA1*03:01-HLA-DPB1*44:01',
         'HLA-DPA1*03:01-HLA-DPB1*45:01', 'HLA-DPA1*03:01-HLA-DPB1*46:01', 'HLA-DPA1*03:01-HLA-DPB1*47:01', 'HLA-DPA1*03:01-HLA-DPB1*48:01',
         'HLA-DPA1*03:01-HLA-DPB1*49:01', 'HLA-DPA1*03:01-HLA-DPB1*50:01',
         'HLA-DPA1*03:01-HLA-DPB1*51:01', 'HLA-DPA1*03:01-HLA-DPB1*52:01', 'HLA-DPA1*03:01-HLA-DPB1*53:01', 'HLA-DPA1*03:01-HLA-DPB1*54:01',
         'HLA-DPA1*03:01-HLA-DPB1*55:01', 'HLA-DPA1*03:01-HLA-DPB1*56:01',
         'HLA-DPA1*03:01-HLA-DPB1*58:01', 'HLA-DPA1*03:01-HLA-DPB1*59:01', 'HLA-DPA1*03:01-HLA-DPB1*60:01', 'HLA-DPA1*03:01-HLA-DPB1*62:01',
         'HLA-DPA1*03:01-HLA-DPB1*63:01', 'HLA-DPA1*03:01-HLA-DPB1*65:01',
         'HLA-DPA1*03:01-HLA-DPB1*66:01', 'HLA-DPA1*03:01-HLA-DPB1*67:01', 'HLA-DPA1*03:01-HLA-DPB1*68:01', 'HLA-DPA1*03:01-HLA-DPB1*69:01',
         'HLA-DPA1*03:01-HLA-DPB1*70:01', 'HLA-DPA1*03:01-HLA-DPB1*71:01',
         'HLA-DPA1*03:01-HLA-DPB1*72:01', 'HLA-DPA1*03:01-HLA-DPB1*73:01', 'HLA-DPA1*03:01-HLA-DPB1*74:01', 'HLA-DPA1*03:01-HLA-DPB1*75:01',
         'HLA-DPA1*03:01-HLA-DPB1*76:01', 'HLA-DPA1*03:01-HLA-DPB1*77:01',
         'HLA-DPA1*03:01-HLA-DPB1*78:01', 'HLA-DPA1*03:01-HLA-DPB1*79:01', 'HLA-DPA1*03:01-HLA-DPB1*80:01', 'HLA-DPA1*03:01-HLA-DPB1*81:01',
         'HLA-DPA1*03:01-HLA-DPB1*82:01', 'HLA-DPA1*03:01-HLA-DPB1*83:01',
         'HLA-DPA1*03:01-HLA-DPB1*84:01', 'HLA-DPA1*03:01-HLA-DPB1*85:01', 'HLA-DPA1*03:01-HLA-DPB1*86:01', 'HLA-DPA1*03:01-HLA-DPB1*87:01',
         'HLA-DPA1*03:01-HLA-DPB1*88:01', 'HLA-DPA1*03:01-HLA-DPB1*89:01',
         'HLA-DPA1*03:01-HLA-DPB1*90:01', 'HLA-DPA1*03:01-HLA-DPB1*91:01', 'HLA-DPA1*03:01-HLA-DPB1*92:01', 'HLA-DPA1*03:01-HLA-DPB1*93:01',
         'HLA-DPA1*03:01-HLA-DPB1*94:01', 'HLA-DPA1*03:01-HLA-DPB1*95:01',
         'HLA-DPA1*03:01-HLA-DPB1*96:01', 'HLA-DPA1*03:01-HLA-DPB1*97:01', 'HLA-DPA1*03:01-HLA-DPB1*98:01', 'HLA-DPA1*03:01-HLA-DPB1*99:01',
         'HLA-DPA1*03:02-HLA-DPB1*01:01', 'HLA-DPA1*03:02-HLA-DPB1*02:01',
         'HLA-DPA1*03:02-HLA-DPB1*02:02', 'HLA-DPA1*03:02-HLA-DPB1*03:01', 'HLA-DPA1*03:02-HLA-DPB1*04:01', 'HLA-DPA1*03:02-HLA-DPB1*04:02',
         'HLA-DPA1*03:02-HLA-DPB1*05:01', 'HLA-DPA1*03:02-HLA-DPB1*06:01',
         'HLA-DPA1*03:02-HLA-DPB1*08:01', 'HLA-DPA1*03:02-HLA-DPB1*09:01', 'HLA-DPA1*03:02-HLA-DPB1*10:001', 'HLA-DPA1*03:02-HLA-DPB1*10:01',
         'HLA-DPA1*03:02-HLA-DPB1*10:101', 'HLA-DPA1*03:02-HLA-DPB1*10:201',
         'HLA-DPA1*03:02-HLA-DPB1*10:301', 'HLA-DPA1*03:02-HLA-DPB1*10:401', 'HLA-DPA1*03:02-HLA-DPB1*10:501', 'HLA-DPA1*03:02-HLA-DPB1*10:601',
         'HLA-DPA1*03:02-HLA-DPB1*10:701', 'HLA-DPA1*03:02-HLA-DPB1*10:801',
         'HLA-DPA1*03:02-HLA-DPB1*10:901', 'HLA-DPA1*03:02-HLA-DPB1*11:001', 'HLA-DPA1*03:02-HLA-DPB1*11:01', 'HLA-DPA1*03:02-HLA-DPB1*11:101',
         'HLA-DPA1*03:02-HLA-DPB1*11:201', 'HLA-DPA1*03:02-HLA-DPB1*11:301',
         'HLA-DPA1*03:02-HLA-DPB1*11:401', 'HLA-DPA1*03:02-HLA-DPB1*11:501', 'HLA-DPA1*03:02-HLA-DPB1*11:601', 'HLA-DPA1*03:02-HLA-DPB1*11:701',
         'HLA-DPA1*03:02-HLA-DPB1*11:801', 'HLA-DPA1*03:02-HLA-DPB1*11:901',
         'HLA-DPA1*03:02-HLA-DPB1*12:101', 'HLA-DPA1*03:02-HLA-DPB1*12:201', 'HLA-DPA1*03:02-HLA-DPB1*12:301', 'HLA-DPA1*03:02-HLA-DPB1*12:401',
         'HLA-DPA1*03:02-HLA-DPB1*12:501', 'HLA-DPA1*03:02-HLA-DPB1*12:601',
         'HLA-DPA1*03:02-HLA-DPB1*12:701', 'HLA-DPA1*03:02-HLA-DPB1*12:801', 'HLA-DPA1*03:02-HLA-DPB1*12:901', 'HLA-DPA1*03:02-HLA-DPB1*13:001',
         'HLA-DPA1*03:02-HLA-DPB1*13:01', 'HLA-DPA1*03:02-HLA-DPB1*13:101',
         'HLA-DPA1*03:02-HLA-DPB1*13:201', 'HLA-DPA1*03:02-HLA-DPB1*13:301', 'HLA-DPA1*03:02-HLA-DPB1*13:401', 'HLA-DPA1*03:02-HLA-DPB1*14:01',
         'HLA-DPA1*03:02-HLA-DPB1*15:01', 'HLA-DPA1*03:02-HLA-DPB1*16:01',
         'HLA-DPA1*03:02-HLA-DPB1*17:01', 'HLA-DPA1*03:02-HLA-DPB1*18:01', 'HLA-DPA1*03:02-HLA-DPB1*19:01', 'HLA-DPA1*03:02-HLA-DPB1*20:01',
         'HLA-DPA1*03:02-HLA-DPB1*21:01', 'HLA-DPA1*03:02-HLA-DPB1*22:01',
         'HLA-DPA1*03:02-HLA-DPB1*23:01', 'HLA-DPA1*03:02-HLA-DPB1*24:01', 'HLA-DPA1*03:02-HLA-DPB1*25:01', 'HLA-DPA1*03:02-HLA-DPB1*26:01',
         'HLA-DPA1*03:02-HLA-DPB1*27:01', 'HLA-DPA1*03:02-HLA-DPB1*28:01',
         'HLA-DPA1*03:02-HLA-DPB1*29:01', 'HLA-DPA1*03:02-HLA-DPB1*30:01', 'HLA-DPA1*03:02-HLA-DPB1*31:01', 'HLA-DPA1*03:02-HLA-DPB1*32:01',
         'HLA-DPA1*03:02-HLA-DPB1*33:01', 'HLA-DPA1*03:02-HLA-DPB1*34:01',
         'HLA-DPA1*03:02-HLA-DPB1*35:01', 'HLA-DPA1*03:02-HLA-DPB1*36:01', 'HLA-DPA1*03:02-HLA-DPB1*37:01', 'HLA-DPA1*03:02-HLA-DPB1*38:01',
         'HLA-DPA1*03:02-HLA-DPB1*39:01', 'HLA-DPA1*03:02-HLA-DPB1*40:01',
         'HLA-DPA1*03:02-HLA-DPB1*41:01', 'HLA-DPA1*03:02-HLA-DPB1*44:01', 'HLA-DPA1*03:02-HLA-DPB1*45:01', 'HLA-DPA1*03:02-HLA-DPB1*46:01',
         'HLA-DPA1*03:02-HLA-DPB1*47:01', 'HLA-DPA1*03:02-HLA-DPB1*48:01',
         'HLA-DPA1*03:02-HLA-DPB1*49:01', 'HLA-DPA1*03:02-HLA-DPB1*50:01', 'HLA-DPA1*03:02-HLA-DPB1*51:01', 'HLA-DPA1*03:02-HLA-DPB1*52:01',
         'HLA-DPA1*03:02-HLA-DPB1*53:01', 'HLA-DPA1*03:02-HLA-DPB1*54:01',
         'HLA-DPA1*03:02-HLA-DPB1*55:01', 'HLA-DPA1*03:02-HLA-DPB1*56:01', 'HLA-DPA1*03:02-HLA-DPB1*58:01', 'HLA-DPA1*03:02-HLA-DPB1*59:01',
         'HLA-DPA1*03:02-HLA-DPB1*60:01', 'HLA-DPA1*03:02-HLA-DPB1*62:01',
         'HLA-DPA1*03:02-HLA-DPB1*63:01', 'HLA-DPA1*03:02-HLA-DPB1*65:01', 'HLA-DPA1*03:02-HLA-DPB1*66:01', 'HLA-DPA1*03:02-HLA-DPB1*67:01',
         'HLA-DPA1*03:02-HLA-DPB1*68:01', 'HLA-DPA1*03:02-HLA-DPB1*69:01',
         'HLA-DPA1*03:02-HLA-DPB1*70:01', 'HLA-DPA1*03:02-HLA-DPB1*71:01', 'HLA-DPA1*03:02-HLA-DPB1*72:01', 'HLA-DPA1*03:02-HLA-DPB1*73:01',
         'HLA-DPA1*03:02-HLA-DPB1*74:01', 'HLA-DPA1*03:02-HLA-DPB1*75:01',
         'HLA-DPA1*03:02-HLA-DPB1*76:01', 'HLA-DPA1*03:02-HLA-DPB1*77:01', 'HLA-DPA1*03:02-HLA-DPB1*78:01', 'HLA-DPA1*03:02-HLA-DPB1*79:01',
         'HLA-DPA1*03:02-HLA-DPB1*80:01', 'HLA-DPA1*03:02-HLA-DPB1*81:01',
         'HLA-DPA1*03:02-HLA-DPB1*82:01', 'HLA-DPA1*03:02-HLA-DPB1*83:01', 'HLA-DPA1*03:02-HLA-DPB1*84:01', 'HLA-DPA1*03:02-HLA-DPB1*85:01',
         'HLA-DPA1*03:02-HLA-DPB1*86:01', 'HLA-DPA1*03:02-HLA-DPB1*87:01',
         'HLA-DPA1*03:02-HLA-DPB1*88:01', 'HLA-DPA1*03:02-HLA-DPB1*89:01', 'HLA-DPA1*03:02-HLA-DPB1*90:01', 'HLA-DPA1*03:02-HLA-DPB1*91:01',
         'HLA-DPA1*03:02-HLA-DPB1*92:01', 'HLA-DPA1*03:02-HLA-DPB1*93:01',
         'HLA-DPA1*03:02-HLA-DPB1*94:01', 'HLA-DPA1*03:02-HLA-DPB1*95:01', 'HLA-DPA1*03:02-HLA-DPB1*96:01', 'HLA-DPA1*03:02-HLA-DPB1*97:01',
         'HLA-DPA1*03:02-HLA-DPB1*98:01', 'HLA-DPA1*03:02-HLA-DPB1*99:01',
         'HLA-DPA1*03:03-HLA-DPB1*01:01', 'HLA-DPA1*03:03-HLA-DPB1*02:01', 'HLA-DPA1*03:03-HLA-DPB1*02:02', 'HLA-DPA1*03:03-HLA-DPB1*03:01',
         'HLA-DPA1*03:03-HLA-DPB1*04:01', 'HLA-DPA1*03:03-HLA-DPB1*04:02',
         'HLA-DPA1*03:03-HLA-DPB1*05:01', 'HLA-DPA1*03:03-HLA-DPB1*06:01', 'HLA-DPA1*03:03-HLA-DPB1*08:01', 'HLA-DPA1*03:03-HLA-DPB1*09:01',
         'HLA-DPA1*03:03-HLA-DPB1*10:001', 'HLA-DPA1*03:03-HLA-DPB1*10:01',
         'HLA-DPA1*03:03-HLA-DPB1*10:101', 'HLA-DPA1*03:03-HLA-DPB1*10:201', 'HLA-DPA1*03:03-HLA-DPB1*10:301', 'HLA-DPA1*03:03-HLA-DPB1*10:401',
         'HLA-DPA1*03:03-HLA-DPB1*10:501', 'HLA-DPA1*03:03-HLA-DPB1*10:601',
         'HLA-DPA1*03:03-HLA-DPB1*10:701', 'HLA-DPA1*03:03-HLA-DPB1*10:801', 'HLA-DPA1*03:03-HLA-DPB1*10:901', 'HLA-DPA1*03:03-HLA-DPB1*11:001',
         'HLA-DPA1*03:03-HLA-DPB1*11:01', 'HLA-DPA1*03:03-HLA-DPB1*11:101',
         'HLA-DPA1*03:03-HLA-DPB1*11:201', 'HLA-DPA1*03:03-HLA-DPB1*11:301', 'HLA-DPA1*03:03-HLA-DPB1*11:401', 'HLA-DPA1*03:03-HLA-DPB1*11:501',
         'HLA-DPA1*03:03-HLA-DPB1*11:601', 'HLA-DPA1*03:03-HLA-DPB1*11:701',
         'HLA-DPA1*03:03-HLA-DPB1*11:801', 'HLA-DPA1*03:03-HLA-DPB1*11:901', 'HLA-DPA1*03:03-HLA-DPB1*12:101', 'HLA-DPA1*03:03-HLA-DPB1*12:201',
         'HLA-DPA1*03:03-HLA-DPB1*12:301', 'HLA-DPA1*03:03-HLA-DPB1*12:401',
         'HLA-DPA1*03:03-HLA-DPB1*12:501', 'HLA-DPA1*03:03-HLA-DPB1*12:601', 'HLA-DPA1*03:03-HLA-DPB1*12:701', 'HLA-DPA1*03:03-HLA-DPB1*12:801',
         'HLA-DPA1*03:03-HLA-DPB1*12:901', 'HLA-DPA1*03:03-HLA-DPB1*13:001',
         'HLA-DPA1*03:03-HLA-DPB1*13:01', 'HLA-DPA1*03:03-HLA-DPB1*13:101', 'HLA-DPA1*03:03-HLA-DPB1*13:201', 'HLA-DPA1*03:03-HLA-DPB1*13:301',
         'HLA-DPA1*03:03-HLA-DPB1*13:401', 'HLA-DPA1*03:03-HLA-DPB1*14:01',
         'HLA-DPA1*03:03-HLA-DPB1*15:01', 'HLA-DPA1*03:03-HLA-DPB1*16:01', 'HLA-DPA1*03:03-HLA-DPB1*17:01', 'HLA-DPA1*03:03-HLA-DPB1*18:01',
         'HLA-DPA1*03:03-HLA-DPB1*19:01', 'HLA-DPA1*03:03-HLA-DPB1*20:01',
         'HLA-DPA1*03:03-HLA-DPB1*21:01', 'HLA-DPA1*03:03-HLA-DPB1*22:01', 'HLA-DPA1*03:03-HLA-DPB1*23:01', 'HLA-DPA1*03:03-HLA-DPB1*24:01',
         'HLA-DPA1*03:03-HLA-DPB1*25:01', 'HLA-DPA1*03:03-HLA-DPB1*26:01',
         'HLA-DPA1*03:03-HLA-DPB1*27:01', 'HLA-DPA1*03:03-HLA-DPB1*28:01', 'HLA-DPA1*03:03-HLA-DPB1*29:01', 'HLA-DPA1*03:03-HLA-DPB1*30:01',
         'HLA-DPA1*03:03-HLA-DPB1*31:01', 'HLA-DPA1*03:03-HLA-DPB1*32:01',
         'HLA-DPA1*03:03-HLA-DPB1*33:01', 'HLA-DPA1*03:03-HLA-DPB1*34:01', 'HLA-DPA1*03:03-HLA-DPB1*35:01', 'HLA-DPA1*03:03-HLA-DPB1*36:01',
         'HLA-DPA1*03:03-HLA-DPB1*37:01', 'HLA-DPA1*03:03-HLA-DPB1*38:01',
         'HLA-DPA1*03:03-HLA-DPB1*39:01', 'HLA-DPA1*03:03-HLA-DPB1*40:01', 'HLA-DPA1*03:03-HLA-DPB1*41:01', 'HLA-DPA1*03:03-HLA-DPB1*44:01',
         'HLA-DPA1*03:03-HLA-DPB1*45:01', 'HLA-DPA1*03:03-HLA-DPB1*46:01',
         'HLA-DPA1*03:03-HLA-DPB1*47:01', 'HLA-DPA1*03:03-HLA-DPB1*48:01', 'HLA-DPA1*03:03-HLA-DPB1*49:01', 'HLA-DPA1*03:03-HLA-DPB1*50:01',
         'HLA-DPA1*03:03-HLA-DPB1*51:01', 'HLA-DPA1*03:03-HLA-DPB1*52:01',
         'HLA-DPA1*03:03-HLA-DPB1*53:01', 'HLA-DPA1*03:03-HLA-DPB1*54:01', 'HLA-DPA1*03:03-HLA-DPB1*55:01', 'HLA-DPA1*03:03-HLA-DPB1*56:01',
         'HLA-DPA1*03:03-HLA-DPB1*58:01', 'HLA-DPA1*03:03-HLA-DPB1*59:01',
         'HLA-DPA1*03:03-HLA-DPB1*60:01', 'HLA-DPA1*03:03-HLA-DPB1*62:01', 'HLA-DPA1*03:03-HLA-DPB1*63:01', 'HLA-DPA1*03:03-HLA-DPB1*65:01',
         'HLA-DPA1*03:03-HLA-DPB1*66:01', 'HLA-DPA1*03:03-HLA-DPB1*67:01',
         'HLA-DPA1*03:03-HLA-DPB1*68:01', 'HLA-DPA1*03:03-HLA-DPB1*69:01', 'HLA-DPA1*03:03-HLA-DPB1*70:01', 'HLA-DPA1*03:03-HLA-DPB1*71:01',
         'HLA-DPA1*03:03-HLA-DPB1*72:01', 'HLA-DPA1*03:03-HLA-DPB1*73:01',
         'HLA-DPA1*03:03-HLA-DPB1*74:01', 'HLA-DPA1*03:03-HLA-DPB1*75:01', 'HLA-DPA1*03:03-HLA-DPB1*76:01', 'HLA-DPA1*03:03-HLA-DPB1*77:01',
         'HLA-DPA1*03:03-HLA-DPB1*78:01', 'HLA-DPA1*03:03-HLA-DPB1*79:01',
         'HLA-DPA1*03:03-HLA-DPB1*80:01', 'HLA-DPA1*03:03-HLA-DPB1*81:01', 'HLA-DPA1*03:03-HLA-DPB1*82:01', 'HLA-DPA1*03:03-HLA-DPB1*83:01',
         'HLA-DPA1*03:03-HLA-DPB1*84:01', 'HLA-DPA1*03:03-HLA-DPB1*85:01',
         'HLA-DPA1*03:03-HLA-DPB1*86:01', 'HLA-DPA1*03:03-HLA-DPB1*87:01', 'HLA-DPA1*03:03-HLA-DPB1*88:01', 'HLA-DPA1*03:03-HLA-DPB1*89:01',
         'HLA-DPA1*03:03-HLA-DPB1*90:01', 'HLA-DPA1*03:03-HLA-DPB1*91:01',
         'HLA-DPA1*03:03-HLA-DPB1*92:01', 'HLA-DPA1*03:03-HLA-DPB1*93:01', 'HLA-DPA1*03:03-HLA-DPB1*94:01', 'HLA-DPA1*03:03-HLA-DPB1*95:01',
         'HLA-DPA1*03:03-HLA-DPB1*96:01', 'HLA-DPA1*03:03-HLA-DPB1*97:01',
         'HLA-DPA1*03:03-HLA-DPB1*98:01', 'HLA-DPA1*03:03-HLA-DPB1*99:01', 'HLA-DPA1*04:01-HLA-DPB1*01:01', 'HLA-DPA1*04:01-HLA-DPB1*02:01',
         'HLA-DPA1*04:01-HLA-DPB1*02:02', 'HLA-DPA1*04:01-HLA-DPB1*03:01',
         'HLA-DPA1*04:01-HLA-DPB1*04:01', 'HLA-DPA1*04:01-HLA-DPB1*04:02', 'HLA-DPA1*04:01-HLA-DPB1*05:01', 'HLA-DPA1*04:01-HLA-DPB1*06:01',
         'HLA-DPA1*04:01-HLA-DPB1*08:01', 'HLA-DPA1*04:01-HLA-DPB1*09:01',
         'HLA-DPA1*04:01-HLA-DPB1*10:001', 'HLA-DPA1*04:01-HLA-DPB1*10:01', 'HLA-DPA1*04:01-HLA-DPB1*10:101', 'HLA-DPA1*04:01-HLA-DPB1*10:201',
         'HLA-DPA1*04:01-HLA-DPB1*10:301', 'HLA-DPA1*04:01-HLA-DPB1*10:401',
         'HLA-DPA1*04:01-HLA-DPB1*10:501', 'HLA-DPA1*04:01-HLA-DPB1*10:601', 'HLA-DPA1*04:01-HLA-DPB1*10:701', 'HLA-DPA1*04:01-HLA-DPB1*10:801',
         'HLA-DPA1*04:01-HLA-DPB1*10:901', 'HLA-DPA1*04:01-HLA-DPB1*11:001',
         'HLA-DPA1*04:01-HLA-DPB1*11:01', 'HLA-DPA1*04:01-HLA-DPB1*11:101', 'HLA-DPA1*04:01-HLA-DPB1*11:201', 'HLA-DPA1*04:01-HLA-DPB1*11:301',
         'HLA-DPA1*04:01-HLA-DPB1*11:401', 'HLA-DPA1*04:01-HLA-DPB1*11:501',
         'HLA-DPA1*04:01-HLA-DPB1*11:601', 'HLA-DPA1*04:01-HLA-DPB1*11:701', 'HLA-DPA1*04:01-HLA-DPB1*11:801', 'HLA-DPA1*04:01-HLA-DPB1*11:901',
         'HLA-DPA1*04:01-HLA-DPB1*12:101', 'HLA-DPA1*04:01-HLA-DPB1*12:201',
         'HLA-DPA1*04:01-HLA-DPB1*12:301', 'HLA-DPA1*04:01-HLA-DPB1*12:401', 'HLA-DPA1*04:01-HLA-DPB1*12:501', 'HLA-DPA1*04:01-HLA-DPB1*12:601',
         'HLA-DPA1*04:01-HLA-DPB1*12:701', 'HLA-DPA1*04:01-HLA-DPB1*12:801',
         'HLA-DPA1*04:01-HLA-DPB1*12:901', 'HLA-DPA1*04:01-HLA-DPB1*13:001', 'HLA-DPA1*04:01-HLA-DPB1*13:01', 'HLA-DPA1*04:01-HLA-DPB1*13:101',
         'HLA-DPA1*04:01-HLA-DPB1*13:201', 'HLA-DPA1*04:01-HLA-DPB1*13:301',
         'HLA-DPA1*04:01-HLA-DPB1*13:401', 'HLA-DPA1*04:01-HLA-DPB1*14:01', 'HLA-DPA1*04:01-HLA-DPB1*15:01', 'HLA-DPA1*04:01-HLA-DPB1*16:01',
         'HLA-DPA1*04:01-HLA-DPB1*17:01', 'HLA-DPA1*04:01-HLA-DPB1*18:01',
         'HLA-DPA1*04:01-HLA-DPB1*19:01', 'HLA-DPA1*04:01-HLA-DPB1*20:01', 'HLA-DPA1*04:01-HLA-DPB1*21:01', 'HLA-DPA1*04:01-HLA-DPB1*22:01',
         'HLA-DPA1*04:01-HLA-DPB1*23:01', 'HLA-DPA1*04:01-HLA-DPB1*24:01',
         'HLA-DPA1*04:01-HLA-DPB1*25:01', 'HLA-DPA1*04:01-HLA-DPB1*26:01', 'HLA-DPA1*04:01-HLA-DPB1*27:01', 'HLA-DPA1*04:01-HLA-DPB1*28:01',
         'HLA-DPA1*04:01-HLA-DPB1*29:01', 'HLA-DPA1*04:01-HLA-DPB1*30:01',
         'HLA-DPA1*04:01-HLA-DPB1*31:01', 'HLA-DPA1*04:01-HLA-DPB1*32:01', 'HLA-DPA1*04:01-HLA-DPB1*33:01', 'HLA-DPA1*04:01-HLA-DPB1*34:01',
         'HLA-DPA1*04:01-HLA-DPB1*35:01', 'HLA-DPA1*04:01-HLA-DPB1*36:01',
         'HLA-DPA1*04:01-HLA-DPB1*37:01', 'HLA-DPA1*04:01-HLA-DPB1*38:01', 'HLA-DPA1*04:01-HLA-DPB1*39:01', 'HLA-DPA1*04:01-HLA-DPB1*40:01',
         'HLA-DPA1*04:01-HLA-DPB1*41:01', 'HLA-DPA1*04:01-HLA-DPB1*44:01',
         'HLA-DPA1*04:01-HLA-DPB1*45:01', 'HLA-DPA1*04:01-HLA-DPB1*46:01', 'HLA-DPA1*04:01-HLA-DPB1*47:01', 'HLA-DPA1*04:01-HLA-DPB1*48:01',
         'HLA-DPA1*04:01-HLA-DPB1*49:01', 'HLA-DPA1*04:01-HLA-DPB1*50:01',
         'HLA-DPA1*04:01-HLA-DPB1*51:01', 'HLA-DPA1*04:01-HLA-DPB1*52:01', 'HLA-DPA1*04:01-HLA-DPB1*53:01', 'HLA-DPA1*04:01-HLA-DPB1*54:01',
         'HLA-DPA1*04:01-HLA-DPB1*55:01', 'HLA-DPA1*04:01-HLA-DPB1*56:01',
         'HLA-DPA1*04:01-HLA-DPB1*58:01', 'HLA-DPA1*04:01-HLA-DPB1*59:01', 'HLA-DPA1*04:01-HLA-DPB1*60:01', 'HLA-DPA1*04:01-HLA-DPB1*62:01',
         'HLA-DPA1*04:01-HLA-DPB1*63:01', 'HLA-DPA1*04:01-HLA-DPB1*65:01',
         'HLA-DPA1*04:01-HLA-DPB1*66:01', 'HLA-DPA1*04:01-HLA-DPB1*67:01', 'HLA-DPA1*04:01-HLA-DPB1*68:01', 'HLA-DPA1*04:01-HLA-DPB1*69:01',
         'HLA-DPA1*04:01-HLA-DPB1*70:01', 'HLA-DPA1*04:01-HLA-DPB1*71:01',
         'HLA-DPA1*04:01-HLA-DPB1*72:01', 'HLA-DPA1*04:01-HLA-DPB1*73:01', 'HLA-DPA1*04:01-HLA-DPB1*74:01', 'HLA-DPA1*04:01-HLA-DPB1*75:01',
         'HLA-DPA1*04:01-HLA-DPB1*76:01', 'HLA-DPA1*04:01-HLA-DPB1*77:01',
         'HLA-DPA1*04:01-HLA-DPB1*78:01', 'HLA-DPA1*04:01-HLA-DPB1*79:01', 'HLA-DPA1*04:01-HLA-DPB1*80:01', 'HLA-DPA1*04:01-HLA-DPB1*81:01',
         'HLA-DPA1*04:01-HLA-DPB1*82:01', 'HLA-DPA1*04:01-HLA-DPB1*83:01',
         'HLA-DPA1*04:01-HLA-DPB1*84:01', 'HLA-DPA1*04:01-HLA-DPB1*85:01', 'HLA-DPA1*04:01-HLA-DPB1*86:01', 'HLA-DPA1*04:01-HLA-DPB1*87:01',
         'HLA-DPA1*04:01-HLA-DPB1*88:01', 'HLA-DPA1*04:01-HLA-DPB1*89:01',
         'HLA-DPA1*04:01-HLA-DPB1*90:01', 'HLA-DPA1*04:01-HLA-DPB1*91:01', 'HLA-DPA1*04:01-HLA-DPB1*92:01', 'HLA-DPA1*04:01-HLA-DPB1*93:01',
         'HLA-DPA1*04:01-HLA-DPB1*94:01', 'HLA-DPA1*04:01-HLA-DPB1*95:01',
         'HLA-DPA1*04:01-HLA-DPB1*96:01', 'HLA-DPA1*04:01-HLA-DPB1*97:01', 'HLA-DPA1*04:01-HLA-DPB1*98:01', 'HLA-DPA1*04:01-HLA-DPB1*99:01',
         'HLA-DQA1*01:01-HLA-DQB1*02:01', 'HLA-DQA1*01:01-HLA-DQB1*02:02',
         'HLA-DQA1*01:01-HLA-DQB1*02:03', 'HLA-DQA1*01:01-HLA-DQB1*02:04', 'HLA-DQA1*01:01-HLA-DQB1*02:05', 'HLA-DQA1*01:01-HLA-DQB1*02:06',
         'HLA-DQA1*01:01-HLA-DQB1*03:01', 'HLA-DQA1*01:01-HLA-DQB1*03:02',
         'HLA-DQA1*01:01-HLA-DQB1*03:03', 'HLA-DQA1*01:01-HLA-DQB1*03:04', 'HLA-DQA1*01:01-HLA-DQB1*03:05', 'HLA-DQA1*01:01-HLA-DQB1*03:06',
         'HLA-DQA1*01:01-HLA-DQB1*03:07', 'HLA-DQA1*01:01-HLA-DQB1*03:08',
         'HLA-DQA1*01:01-HLA-DQB1*03:09', 'HLA-DQA1*01:01-HLA-DQB1*03:10', 'HLA-DQA1*01:01-HLA-DQB1*03:11', 'HLA-DQA1*01:01-HLA-DQB1*03:12',
         'HLA-DQA1*01:01-HLA-DQB1*03:13', 'HLA-DQA1*01:01-HLA-DQB1*03:14',
         'HLA-DQA1*01:01-HLA-DQB1*03:15', 'HLA-DQA1*01:01-HLA-DQB1*03:16', 'HLA-DQA1*01:01-HLA-DQB1*03:17', 'HLA-DQA1*01:01-HLA-DQB1*03:18',
         'HLA-DQA1*01:01-HLA-DQB1*03:19', 'HLA-DQA1*01:01-HLA-DQB1*03:20',
         'HLA-DQA1*01:01-HLA-DQB1*03:21', 'HLA-DQA1*01:01-HLA-DQB1*03:22', 'HLA-DQA1*01:01-HLA-DQB1*03:23', 'HLA-DQA1*01:01-HLA-DQB1*03:24',
         'HLA-DQA1*01:01-HLA-DQB1*03:25', 'HLA-DQA1*01:01-HLA-DQB1*03:26',
         'HLA-DQA1*01:01-HLA-DQB1*03:27', 'HLA-DQA1*01:01-HLA-DQB1*03:28', 'HLA-DQA1*01:01-HLA-DQB1*03:29', 'HLA-DQA1*01:01-HLA-DQB1*03:30',
         'HLA-DQA1*01:01-HLA-DQB1*03:31', 'HLA-DQA1*01:01-HLA-DQB1*03:32',
         'HLA-DQA1*01:01-HLA-DQB1*03:33', 'HLA-DQA1*01:01-HLA-DQB1*03:34', 'HLA-DQA1*01:01-HLA-DQB1*03:35', 'HLA-DQA1*01:01-HLA-DQB1*03:36',
         'HLA-DQA1*01:01-HLA-DQB1*03:37', 'HLA-DQA1*01:01-HLA-DQB1*03:38',
         'HLA-DQA1*01:01-HLA-DQB1*04:01', 'HLA-DQA1*01:01-HLA-DQB1*04:02', 'HLA-DQA1*01:01-HLA-DQB1*04:03', 'HLA-DQA1*01:01-HLA-DQB1*04:04',
         'HLA-DQA1*01:01-HLA-DQB1*04:05', 'HLA-DQA1*01:01-HLA-DQB1*04:06',
         'HLA-DQA1*01:01-HLA-DQB1*04:07', 'HLA-DQA1*01:01-HLA-DQB1*04:08', 'HLA-DQA1*01:01-HLA-DQB1*05:01', 'HLA-DQA1*01:01-HLA-DQB1*05:02',
         'HLA-DQA1*01:01-HLA-DQB1*05:03', 'HLA-DQA1*01:01-HLA-DQB1*05:05',
         'HLA-DQA1*01:01-HLA-DQB1*05:06', 'HLA-DQA1*01:01-HLA-DQB1*05:07', 'HLA-DQA1*01:01-HLA-DQB1*05:08', 'HLA-DQA1*01:01-HLA-DQB1*05:09',
         'HLA-DQA1*01:01-HLA-DQB1*05:10', 'HLA-DQA1*01:01-HLA-DQB1*05:11',
         'HLA-DQA1*01:01-HLA-DQB1*05:12', 'HLA-DQA1*01:01-HLA-DQB1*05:13', 'HLA-DQA1*01:01-HLA-DQB1*05:14', 'HLA-DQA1*01:01-HLA-DQB1*06:01',
         'HLA-DQA1*01:01-HLA-DQB1*06:02', 'HLA-DQA1*01:01-HLA-DQB1*06:03',
         'HLA-DQA1*01:01-HLA-DQB1*06:04', 'HLA-DQA1*01:01-HLA-DQB1*06:07', 'HLA-DQA1*01:01-HLA-DQB1*06:08', 'HLA-DQA1*01:01-HLA-DQB1*06:09',
         'HLA-DQA1*01:01-HLA-DQB1*06:10', 'HLA-DQA1*01:01-HLA-DQB1*06:11',
         'HLA-DQA1*01:01-HLA-DQB1*06:12', 'HLA-DQA1*01:01-HLA-DQB1*06:14', 'HLA-DQA1*01:01-HLA-DQB1*06:15', 'HLA-DQA1*01:01-HLA-DQB1*06:16',
         'HLA-DQA1*01:01-HLA-DQB1*06:17', 'HLA-DQA1*01:01-HLA-DQB1*06:18',
         'HLA-DQA1*01:01-HLA-DQB1*06:19', 'HLA-DQA1*01:01-HLA-DQB1*06:21', 'HLA-DQA1*01:01-HLA-DQB1*06:22', 'HLA-DQA1*01:01-HLA-DQB1*06:23',
         'HLA-DQA1*01:01-HLA-DQB1*06:24', 'HLA-DQA1*01:01-HLA-DQB1*06:25',
         'HLA-DQA1*01:01-HLA-DQB1*06:27', 'HLA-DQA1*01:01-HLA-DQB1*06:28', 'HLA-DQA1*01:01-HLA-DQB1*06:29', 'HLA-DQA1*01:01-HLA-DQB1*06:30',
         'HLA-DQA1*01:01-HLA-DQB1*06:31', 'HLA-DQA1*01:01-HLA-DQB1*06:32',
         'HLA-DQA1*01:01-HLA-DQB1*06:33', 'HLA-DQA1*01:01-HLA-DQB1*06:34', 'HLA-DQA1*01:01-HLA-DQB1*06:35', 'HLA-DQA1*01:01-HLA-DQB1*06:36',
         'HLA-DQA1*01:01-HLA-DQB1*06:37', 'HLA-DQA1*01:01-HLA-DQB1*06:38',
         'HLA-DQA1*01:01-HLA-DQB1*06:39', 'HLA-DQA1*01:01-HLA-DQB1*06:40', 'HLA-DQA1*01:01-HLA-DQB1*06:41', 'HLA-DQA1*01:01-HLA-DQB1*06:42',
         'HLA-DQA1*01:01-HLA-DQB1*06:43', 'HLA-DQA1*01:01-HLA-DQB1*06:44',
         'HLA-DQA1*01:02-HLA-DQB1*02:01', 'HLA-DQA1*01:02-HLA-DQB1*02:02', 'HLA-DQA1*01:02-HLA-DQB1*02:03', 'HLA-DQA1*01:02-HLA-DQB1*02:04',
         'HLA-DQA1*01:02-HLA-DQB1*02:05', 'HLA-DQA1*01:02-HLA-DQB1*02:06',
         'HLA-DQA1*01:02-HLA-DQB1*03:01', 'HLA-DQA1*01:02-HLA-DQB1*03:02', 'HLA-DQA1*01:02-HLA-DQB1*03:03', 'HLA-DQA1*01:02-HLA-DQB1*03:04',
         'HLA-DQA1*01:02-HLA-DQB1*03:05', 'HLA-DQA1*01:02-HLA-DQB1*03:06',
         'HLA-DQA1*01:02-HLA-DQB1*03:07', 'HLA-DQA1*01:02-HLA-DQB1*03:08', 'HLA-DQA1*01:02-HLA-DQB1*03:09', 'HLA-DQA1*01:02-HLA-DQB1*03:10',
         'HLA-DQA1*01:02-HLA-DQB1*03:11', 'HLA-DQA1*01:02-HLA-DQB1*03:12',
         'HLA-DQA1*01:02-HLA-DQB1*03:13', 'HLA-DQA1*01:02-HLA-DQB1*03:14', 'HLA-DQA1*01:02-HLA-DQB1*03:15', 'HLA-DQA1*01:02-HLA-DQB1*03:16',
         'HLA-DQA1*01:02-HLA-DQB1*03:17', 'HLA-DQA1*01:02-HLA-DQB1*03:18',
         'HLA-DQA1*01:02-HLA-DQB1*03:19', 'HLA-DQA1*01:02-HLA-DQB1*03:20', 'HLA-DQA1*01:02-HLA-DQB1*03:21', 'HLA-DQA1*01:02-HLA-DQB1*03:22',
         'HLA-DQA1*01:02-HLA-DQB1*03:23', 'HLA-DQA1*01:02-HLA-DQB1*03:24',
         'HLA-DQA1*01:02-HLA-DQB1*03:25', 'HLA-DQA1*01:02-HLA-DQB1*03:26', 'HLA-DQA1*01:02-HLA-DQB1*03:27', 'HLA-DQA1*01:02-HLA-DQB1*03:28',
         'HLA-DQA1*01:02-HLA-DQB1*03:29', 'HLA-DQA1*01:02-HLA-DQB1*03:30',
         'HLA-DQA1*01:02-HLA-DQB1*03:31', 'HLA-DQA1*01:02-HLA-DQB1*03:32', 'HLA-DQA1*01:02-HLA-DQB1*03:33', 'HLA-DQA1*01:02-HLA-DQB1*03:34',
         'HLA-DQA1*01:02-HLA-DQB1*03:35', 'HLA-DQA1*01:02-HLA-DQB1*03:36',
         'HLA-DQA1*01:02-HLA-DQB1*03:37', 'HLA-DQA1*01:02-HLA-DQB1*03:38', 'HLA-DQA1*01:02-HLA-DQB1*04:01', 'HLA-DQA1*01:02-HLA-DQB1*04:02',
         'HLA-DQA1*01:02-HLA-DQB1*04:03', 'HLA-DQA1*01:02-HLA-DQB1*04:04',
         'HLA-DQA1*01:02-HLA-DQB1*04:05', 'HLA-DQA1*01:02-HLA-DQB1*04:06', 'HLA-DQA1*01:02-HLA-DQB1*04:07', 'HLA-DQA1*01:02-HLA-DQB1*04:08',
         'HLA-DQA1*01:02-HLA-DQB1*05:01', 'HLA-DQA1*01:02-HLA-DQB1*05:02',
         'HLA-DQA1*01:02-HLA-DQB1*05:03', 'HLA-DQA1*01:02-HLA-DQB1*05:05', 'HLA-DQA1*01:02-HLA-DQB1*05:06', 'HLA-DQA1*01:02-HLA-DQB1*05:07',
         'HLA-DQA1*01:02-HLA-DQB1*05:08', 'HLA-DQA1*01:02-HLA-DQB1*05:09',
         'HLA-DQA1*01:02-HLA-DQB1*05:10', 'HLA-DQA1*01:02-HLA-DQB1*05:11', 'HLA-DQA1*01:02-HLA-DQB1*05:12', 'HLA-DQA1*01:02-HLA-DQB1*05:13',
         'HLA-DQA1*01:02-HLA-DQB1*05:14', 'HLA-DQA1*01:02-HLA-DQB1*06:01',
         'HLA-DQA1*01:02-HLA-DQB1*06:02', 'HLA-DQA1*01:02-HLA-DQB1*06:03', 'HLA-DQA1*01:02-HLA-DQB1*06:04', 'HLA-DQA1*01:02-HLA-DQB1*06:07',
         'HLA-DQA1*01:02-HLA-DQB1*06:08', 'HLA-DQA1*01:02-HLA-DQB1*06:09',
         'HLA-DQA1*01:02-HLA-DQB1*06:10', 'HLA-DQA1*01:02-HLA-DQB1*06:11', 'HLA-DQA1*01:02-HLA-DQB1*06:12', 'HLA-DQA1*01:02-HLA-DQB1*06:14',
         'HLA-DQA1*01:02-HLA-DQB1*06:15', 'HLA-DQA1*01:02-HLA-DQB1*06:16',
         'HLA-DQA1*01:02-HLA-DQB1*06:17', 'HLA-DQA1*01:02-HLA-DQB1*06:18', 'HLA-DQA1*01:02-HLA-DQB1*06:19', 'HLA-DQA1*01:02-HLA-DQB1*06:21',
         'HLA-DQA1*01:02-HLA-DQB1*06:22', 'HLA-DQA1*01:02-HLA-DQB1*06:23',
         'HLA-DQA1*01:02-HLA-DQB1*06:24', 'HLA-DQA1*01:02-HLA-DQB1*06:25', 'HLA-DQA1*01:02-HLA-DQB1*06:27', 'HLA-DQA1*01:02-HLA-DQB1*06:28',
         'HLA-DQA1*01:02-HLA-DQB1*06:29', 'HLA-DQA1*01:02-HLA-DQB1*06:30',
         'HLA-DQA1*01:02-HLA-DQB1*06:31', 'HLA-DQA1*01:02-HLA-DQB1*06:32', 'HLA-DQA1*01:02-HLA-DQB1*06:33', 'HLA-DQA1*01:02-HLA-DQB1*06:34',
         'HLA-DQA1*01:02-HLA-DQB1*06:35', 'HLA-DQA1*01:02-HLA-DQB1*06:36',
         'HLA-DQA1*01:02-HLA-DQB1*06:37', 'HLA-DQA1*01:02-HLA-DQB1*06:38', 'HLA-DQA1*01:02-HLA-DQB1*06:39', 'HLA-DQA1*01:02-HLA-DQB1*06:40',
         'HLA-DQA1*01:02-HLA-DQB1*06:41', 'HLA-DQA1*01:02-HLA-DQB1*06:42',
         'HLA-DQA1*01:02-HLA-DQB1*06:43', 'HLA-DQA1*01:02-HLA-DQB1*06:44', 'HLA-DQA1*01:03-HLA-DQB1*02:01', 'HLA-DQA1*01:03-HLA-DQB1*02:02',
         'HLA-DQA1*01:03-HLA-DQB1*02:03', 'HLA-DQA1*01:03-HLA-DQB1*02:04',
         'HLA-DQA1*01:03-HLA-DQB1*02:05', 'HLA-DQA1*01:03-HLA-DQB1*02:06', 'HLA-DQA1*01:03-HLA-DQB1*03:01', 'HLA-DQA1*01:03-HLA-DQB1*03:02',
         'HLA-DQA1*01:03-HLA-DQB1*03:03', 'HLA-DQA1*01:03-HLA-DQB1*03:04',
         'HLA-DQA1*01:03-HLA-DQB1*03:05', 'HLA-DQA1*01:03-HLA-DQB1*03:06', 'HLA-DQA1*01:03-HLA-DQB1*03:07', 'HLA-DQA1*01:03-HLA-DQB1*03:08',
         'HLA-DQA1*01:03-HLA-DQB1*03:09', 'HLA-DQA1*01:03-HLA-DQB1*03:10',
         'HLA-DQA1*01:03-HLA-DQB1*03:11', 'HLA-DQA1*01:03-HLA-DQB1*03:12', 'HLA-DQA1*01:03-HLA-DQB1*03:13', 'HLA-DQA1*01:03-HLA-DQB1*03:14',
         'HLA-DQA1*01:03-HLA-DQB1*03:15', 'HLA-DQA1*01:03-HLA-DQB1*03:16',
         'HLA-DQA1*01:03-HLA-DQB1*03:17', 'HLA-DQA1*01:03-HLA-DQB1*03:18', 'HLA-DQA1*01:03-HLA-DQB1*03:19', 'HLA-DQA1*01:03-HLA-DQB1*03:20',
         'HLA-DQA1*01:03-HLA-DQB1*03:21', 'HLA-DQA1*01:03-HLA-DQB1*03:22',
         'HLA-DQA1*01:03-HLA-DQB1*03:23', 'HLA-DQA1*01:03-HLA-DQB1*03:24', 'HLA-DQA1*01:03-HLA-DQB1*03:25', 'HLA-DQA1*01:03-HLA-DQB1*03:26',
         'HLA-DQA1*01:03-HLA-DQB1*03:27', 'HLA-DQA1*01:03-HLA-DQB1*03:28',
         'HLA-DQA1*01:03-HLA-DQB1*03:29', 'HLA-DQA1*01:03-HLA-DQB1*03:30', 'HLA-DQA1*01:03-HLA-DQB1*03:31', 'HLA-DQA1*01:03-HLA-DQB1*03:32',
         'HLA-DQA1*01:03-HLA-DQB1*03:33', 'HLA-DQA1*01:03-HLA-DQB1*03:34',
         'HLA-DQA1*01:03-HLA-DQB1*03:35', 'HLA-DQA1*01:03-HLA-DQB1*03:36', 'HLA-DQA1*01:03-HLA-DQB1*03:37', 'HLA-DQA1*01:03-HLA-DQB1*03:38',
         'HLA-DQA1*01:03-HLA-DQB1*04:01', 'HLA-DQA1*01:03-HLA-DQB1*04:02',
         'HLA-DQA1*01:03-HLA-DQB1*04:03', 'HLA-DQA1*01:03-HLA-DQB1*04:04', 'HLA-DQA1*01:03-HLA-DQB1*04:05', 'HLA-DQA1*01:03-HLA-DQB1*04:06',
         'HLA-DQA1*01:03-HLA-DQB1*04:07', 'HLA-DQA1*01:03-HLA-DQB1*04:08',
         'HLA-DQA1*01:03-HLA-DQB1*05:01', 'HLA-DQA1*01:03-HLA-DQB1*05:02', 'HLA-DQA1*01:03-HLA-DQB1*05:03', 'HLA-DQA1*01:03-HLA-DQB1*05:05',
         'HLA-DQA1*01:03-HLA-DQB1*05:06', 'HLA-DQA1*01:03-HLA-DQB1*05:07',
         'HLA-DQA1*01:03-HLA-DQB1*05:08', 'HLA-DQA1*01:03-HLA-DQB1*05:09', 'HLA-DQA1*01:03-HLA-DQB1*05:10', 'HLA-DQA1*01:03-HLA-DQB1*05:11',
         'HLA-DQA1*01:03-HLA-DQB1*05:12', 'HLA-DQA1*01:03-HLA-DQB1*05:13',
         'HLA-DQA1*01:03-HLA-DQB1*05:14', 'HLA-DQA1*01:03-HLA-DQB1*06:01', 'HLA-DQA1*01:03-HLA-DQB1*06:02', 'HLA-DQA1*01:03-HLA-DQB1*06:03',
         'HLA-DQA1*01:03-HLA-DQB1*06:04', 'HLA-DQA1*01:03-HLA-DQB1*06:07',
         'HLA-DQA1*01:03-HLA-DQB1*06:08', 'HLA-DQA1*01:03-HLA-DQB1*06:09', 'HLA-DQA1*01:03-HLA-DQB1*06:10', 'HLA-DQA1*01:03-HLA-DQB1*06:11',
         'HLA-DQA1*01:03-HLA-DQB1*06:12', 'HLA-DQA1*01:03-HLA-DQB1*06:14',
         'HLA-DQA1*01:03-HLA-DQB1*06:15', 'HLA-DQA1*01:03-HLA-DQB1*06:16', 'HLA-DQA1*01:03-HLA-DQB1*06:17', 'HLA-DQA1*01:03-HLA-DQB1*06:18',
         'HLA-DQA1*01:03-HLA-DQB1*06:19', 'HLA-DQA1*01:03-HLA-DQB1*06:21',
         'HLA-DQA1*01:03-HLA-DQB1*06:22', 'HLA-DQA1*01:03-HLA-DQB1*06:23', 'HLA-DQA1*01:03-HLA-DQB1*06:24', 'HLA-DQA1*01:03-HLA-DQB1*06:25',
         'HLA-DQA1*01:03-HLA-DQB1*06:27', 'HLA-DQA1*01:03-HLA-DQB1*06:28',
         'HLA-DQA1*01:03-HLA-DQB1*06:29', 'HLA-DQA1*01:03-HLA-DQB1*06:30', 'HLA-DQA1*01:03-HLA-DQB1*06:31', 'HLA-DQA1*01:03-HLA-DQB1*06:32',
         'HLA-DQA1*01:03-HLA-DQB1*06:33', 'HLA-DQA1*01:03-HLA-DQB1*06:34',
         'HLA-DQA1*01:03-HLA-DQB1*06:35', 'HLA-DQA1*01:03-HLA-DQB1*06:36', 'HLA-DQA1*01:03-HLA-DQB1*06:37', 'HLA-DQA1*01:03-HLA-DQB1*06:38',
         'HLA-DQA1*01:03-HLA-DQB1*06:39', 'HLA-DQA1*01:03-HLA-DQB1*06:40',
         'HLA-DQA1*01:03-HLA-DQB1*06:41', 'HLA-DQA1*01:03-HLA-DQB1*06:42', 'HLA-DQA1*01:03-HLA-DQB1*06:43', 'HLA-DQA1*01:03-HLA-DQB1*06:44',
         'HLA-DQA1*01:04-HLA-DQB1*02:01', 'HLA-DQA1*01:04-HLA-DQB1*02:02',
         'HLA-DQA1*01:04-HLA-DQB1*02:03', 'HLA-DQA1*01:04-HLA-DQB1*02:04', 'HLA-DQA1*01:04-HLA-DQB1*02:05', 'HLA-DQA1*01:04-HLA-DQB1*02:06',
         'HLA-DQA1*01:04-HLA-DQB1*03:01', 'HLA-DQA1*01:04-HLA-DQB1*03:02',
         'HLA-DQA1*01:04-HLA-DQB1*03:03', 'HLA-DQA1*01:04-HLA-DQB1*03:04', 'HLA-DQA1*01:04-HLA-DQB1*03:05', 'HLA-DQA1*01:04-HLA-DQB1*03:06',
         'HLA-DQA1*01:04-HLA-DQB1*03:07', 'HLA-DQA1*01:04-HLA-DQB1*03:08',
         'HLA-DQA1*01:04-HLA-DQB1*03:09', 'HLA-DQA1*01:04-HLA-DQB1*03:10', 'HLA-DQA1*01:04-HLA-DQB1*03:11', 'HLA-DQA1*01:04-HLA-DQB1*03:12',
         'HLA-DQA1*01:04-HLA-DQB1*03:13', 'HLA-DQA1*01:04-HLA-DQB1*03:14',
         'HLA-DQA1*01:04-HLA-DQB1*03:15', 'HLA-DQA1*01:04-HLA-DQB1*03:16', 'HLA-DQA1*01:04-HLA-DQB1*03:17', 'HLA-DQA1*01:04-HLA-DQB1*03:18',
         'HLA-DQA1*01:04-HLA-DQB1*03:19', 'HLA-DQA1*01:04-HLA-DQB1*03:20',
         'HLA-DQA1*01:04-HLA-DQB1*03:21', 'HLA-DQA1*01:04-HLA-DQB1*03:22', 'HLA-DQA1*01:04-HLA-DQB1*03:23', 'HLA-DQA1*01:04-HLA-DQB1*03:24',
         'HLA-DQA1*01:04-HLA-DQB1*03:25', 'HLA-DQA1*01:04-HLA-DQB1*03:26',
         'HLA-DQA1*01:04-HLA-DQB1*03:27', 'HLA-DQA1*01:04-HLA-DQB1*03:28', 'HLA-DQA1*01:04-HLA-DQB1*03:29', 'HLA-DQA1*01:04-HLA-DQB1*03:30',
         'HLA-DQA1*01:04-HLA-DQB1*03:31', 'HLA-DQA1*01:04-HLA-DQB1*03:32',
         'HLA-DQA1*01:04-HLA-DQB1*03:33', 'HLA-DQA1*01:04-HLA-DQB1*03:34', 'HLA-DQA1*01:04-HLA-DQB1*03:35', 'HLA-DQA1*01:04-HLA-DQB1*03:36',
         'HLA-DQA1*01:04-HLA-DQB1*03:37', 'HLA-DQA1*01:04-HLA-DQB1*03:38',
         'HLA-DQA1*01:04-HLA-DQB1*04:01', 'HLA-DQA1*01:04-HLA-DQB1*04:02', 'HLA-DQA1*01:04-HLA-DQB1*04:03', 'HLA-DQA1*01:04-HLA-DQB1*04:04',
         'HLA-DQA1*01:04-HLA-DQB1*04:05', 'HLA-DQA1*01:04-HLA-DQB1*04:06',
         'HLA-DQA1*01:04-HLA-DQB1*04:07', 'HLA-DQA1*01:04-HLA-DQB1*04:08', 'HLA-DQA1*01:04-HLA-DQB1*05:01', 'HLA-DQA1*01:04-HLA-DQB1*05:02',
         'HLA-DQA1*01:04-HLA-DQB1*05:03', 'HLA-DQA1*01:04-HLA-DQB1*05:05',
         'HLA-DQA1*01:04-HLA-DQB1*05:06', 'HLA-DQA1*01:04-HLA-DQB1*05:07', 'HLA-DQA1*01:04-HLA-DQB1*05:08', 'HLA-DQA1*01:04-HLA-DQB1*05:09',
         'HLA-DQA1*01:04-HLA-DQB1*05:10', 'HLA-DQA1*01:04-HLA-DQB1*05:11',
         'HLA-DQA1*01:04-HLA-DQB1*05:12', 'HLA-DQA1*01:04-HLA-DQB1*05:13', 'HLA-DQA1*01:04-HLA-DQB1*05:14', 'HLA-DQA1*01:04-HLA-DQB1*06:01',
         'HLA-DQA1*01:04-HLA-DQB1*06:02', 'HLA-DQA1*01:04-HLA-DQB1*06:03',
         'HLA-DQA1*01:04-HLA-DQB1*06:04', 'HLA-DQA1*01:04-HLA-DQB1*06:07', 'HLA-DQA1*01:04-HLA-DQB1*06:08', 'HLA-DQA1*01:04-HLA-DQB1*06:09',
         'HLA-DQA1*01:04-HLA-DQB1*06:10', 'HLA-DQA1*01:04-HLA-DQB1*06:11',
         'HLA-DQA1*01:04-HLA-DQB1*06:12', 'HLA-DQA1*01:04-HLA-DQB1*06:14', 'HLA-DQA1*01:04-HLA-DQB1*06:15', 'HLA-DQA1*01:04-HLA-DQB1*06:16',
         'HLA-DQA1*01:04-HLA-DQB1*06:17', 'HLA-DQA1*01:04-HLA-DQB1*06:18',
         'HLA-DQA1*01:04-HLA-DQB1*06:19', 'HLA-DQA1*01:04-HLA-DQB1*06:21', 'HLA-DQA1*01:04-HLA-DQB1*06:22', 'HLA-DQA1*01:04-HLA-DQB1*06:23',
         'HLA-DQA1*01:04-HLA-DQB1*06:24', 'HLA-DQA1*01:04-HLA-DQB1*06:25',
         'HLA-DQA1*01:04-HLA-DQB1*06:27', 'HLA-DQA1*01:04-HLA-DQB1*06:28', 'HLA-DQA1*01:04-HLA-DQB1*06:29', 'HLA-DQA1*01:04-HLA-DQB1*06:30',
         'HLA-DQA1*01:04-HLA-DQB1*06:31', 'HLA-DQA1*01:04-HLA-DQB1*06:32',
         'HLA-DQA1*01:04-HLA-DQB1*06:33', 'HLA-DQA1*01:04-HLA-DQB1*06:34', 'HLA-DQA1*01:04-HLA-DQB1*06:35', 'HLA-DQA1*01:04-HLA-DQB1*06:36',
         'HLA-DQA1*01:04-HLA-DQB1*06:37', 'HLA-DQA1*01:04-HLA-DQB1*06:38',
         'HLA-DQA1*01:04-HLA-DQB1*06:39', 'HLA-DQA1*01:04-HLA-DQB1*06:40', 'HLA-DQA1*01:04-HLA-DQB1*06:41', 'HLA-DQA1*01:04-HLA-DQB1*06:42',
         'HLA-DQA1*01:04-HLA-DQB1*06:43', 'HLA-DQA1*01:04-HLA-DQB1*06:44',
         'HLA-DQA1*01:05-HLA-DQB1*02:01', 'HLA-DQA1*01:05-HLA-DQB1*02:02', 'HLA-DQA1*01:05-HLA-DQB1*02:03', 'HLA-DQA1*01:05-HLA-DQB1*02:04',
         'HLA-DQA1*01:05-HLA-DQB1*02:05', 'HLA-DQA1*01:05-HLA-DQB1*02:06',
         'HLA-DQA1*01:05-HLA-DQB1*03:01', 'HLA-DQA1*01:05-HLA-DQB1*03:02', 'HLA-DQA1*01:05-HLA-DQB1*03:03', 'HLA-DQA1*01:05-HLA-DQB1*03:04',
         'HLA-DQA1*01:05-HLA-DQB1*03:05', 'HLA-DQA1*01:05-HLA-DQB1*03:06',
         'HLA-DQA1*01:05-HLA-DQB1*03:07', 'HLA-DQA1*01:05-HLA-DQB1*03:08', 'HLA-DQA1*01:05-HLA-DQB1*03:09', 'HLA-DQA1*01:05-HLA-DQB1*03:10',
         'HLA-DQA1*01:05-HLA-DQB1*03:11', 'HLA-DQA1*01:05-HLA-DQB1*03:12',
         'HLA-DQA1*01:05-HLA-DQB1*03:13', 'HLA-DQA1*01:05-HLA-DQB1*03:14', 'HLA-DQA1*01:05-HLA-DQB1*03:15', 'HLA-DQA1*01:05-HLA-DQB1*03:16',
         'HLA-DQA1*01:05-HLA-DQB1*03:17', 'HLA-DQA1*01:05-HLA-DQB1*03:18',
         'HLA-DQA1*01:05-HLA-DQB1*03:19', 'HLA-DQA1*01:05-HLA-DQB1*03:20', 'HLA-DQA1*01:05-HLA-DQB1*03:21', 'HLA-DQA1*01:05-HLA-DQB1*03:22',
         'HLA-DQA1*01:05-HLA-DQB1*03:23', 'HLA-DQA1*01:05-HLA-DQB1*03:24',
         'HLA-DQA1*01:05-HLA-DQB1*03:25', 'HLA-DQA1*01:05-HLA-DQB1*03:26', 'HLA-DQA1*01:05-HLA-DQB1*03:27', 'HLA-DQA1*01:05-HLA-DQB1*03:28',
         'HLA-DQA1*01:05-HLA-DQB1*03:29', 'HLA-DQA1*01:05-HLA-DQB1*03:30',
         'HLA-DQA1*01:05-HLA-DQB1*03:31', 'HLA-DQA1*01:05-HLA-DQB1*03:32', 'HLA-DQA1*01:05-HLA-DQB1*03:33', 'HLA-DQA1*01:05-HLA-DQB1*03:34',
         'HLA-DQA1*01:05-HLA-DQB1*03:35', 'HLA-DQA1*01:05-HLA-DQB1*03:36',
         'HLA-DQA1*01:05-HLA-DQB1*03:37', 'HLA-DQA1*01:05-HLA-DQB1*03:38', 'HLA-DQA1*01:05-HLA-DQB1*04:01', 'HLA-DQA1*01:05-HLA-DQB1*04:02',
         'HLA-DQA1*01:05-HLA-DQB1*04:03', 'HLA-DQA1*01:05-HLA-DQB1*04:04',
         'HLA-DQA1*01:05-HLA-DQB1*04:05', 'HLA-DQA1*01:05-HLA-DQB1*04:06', 'HLA-DQA1*01:05-HLA-DQB1*04:07', 'HLA-DQA1*01:05-HLA-DQB1*04:08',
         'HLA-DQA1*01:05-HLA-DQB1*05:01', 'HLA-DQA1*01:05-HLA-DQB1*05:02',
         'HLA-DQA1*01:05-HLA-DQB1*05:03', 'HLA-DQA1*01:05-HLA-DQB1*05:05', 'HLA-DQA1*01:05-HLA-DQB1*05:06', 'HLA-DQA1*01:05-HLA-DQB1*05:07',
         'HLA-DQA1*01:05-HLA-DQB1*05:08', 'HLA-DQA1*01:05-HLA-DQB1*05:09',
         'HLA-DQA1*01:05-HLA-DQB1*05:10', 'HLA-DQA1*01:05-HLA-DQB1*05:11', 'HLA-DQA1*01:05-HLA-DQB1*05:12', 'HLA-DQA1*01:05-HLA-DQB1*05:13',
         'HLA-DQA1*01:05-HLA-DQB1*05:14', 'HLA-DQA1*01:05-HLA-DQB1*06:01',
         'HLA-DQA1*01:05-HLA-DQB1*06:02', 'HLA-DQA1*01:05-HLA-DQB1*06:03', 'HLA-DQA1*01:05-HLA-DQB1*06:04', 'HLA-DQA1*01:05-HLA-DQB1*06:07',
         'HLA-DQA1*01:05-HLA-DQB1*06:08', 'HLA-DQA1*01:05-HLA-DQB1*06:09',
         'HLA-DQA1*01:05-HLA-DQB1*06:10', 'HLA-DQA1*01:05-HLA-DQB1*06:11', 'HLA-DQA1*01:05-HLA-DQB1*06:12', 'HLA-DQA1*01:05-HLA-DQB1*06:14',
         'HLA-DQA1*01:05-HLA-DQB1*06:15', 'HLA-DQA1*01:05-HLA-DQB1*06:16',
         'HLA-DQA1*01:05-HLA-DQB1*06:17', 'HLA-DQA1*01:05-HLA-DQB1*06:18', 'HLA-DQA1*01:05-HLA-DQB1*06:19', 'HLA-DQA1*01:05-HLA-DQB1*06:21',
         'HLA-DQA1*01:05-HLA-DQB1*06:22', 'HLA-DQA1*01:05-HLA-DQB1*06:23',
         'HLA-DQA1*01:05-HLA-DQB1*06:24', 'HLA-DQA1*01:05-HLA-DQB1*06:25', 'HLA-DQA1*01:05-HLA-DQB1*06:27', 'HLA-DQA1*01:05-HLA-DQB1*06:28',
         'HLA-DQA1*01:05-HLA-DQB1*06:29', 'HLA-DQA1*01:05-HLA-DQB1*06:30',
         'HLA-DQA1*01:05-HLA-DQB1*06:31', 'HLA-DQA1*01:05-HLA-DQB1*06:32', 'HLA-DQA1*01:05-HLA-DQB1*06:33', 'HLA-DQA1*01:05-HLA-DQB1*06:34',
         'HLA-DQA1*01:05-HLA-DQB1*06:35', 'HLA-DQA1*01:05-HLA-DQB1*06:36',
         'HLA-DQA1*01:05-HLA-DQB1*06:37', 'HLA-DQA1*01:05-HLA-DQB1*06:38', 'HLA-DQA1*01:05-HLA-DQB1*06:39', 'HLA-DQA1*01:05-HLA-DQB1*06:40',
         'HLA-DQA1*01:05-HLA-DQB1*06:41', 'HLA-DQA1*01:05-HLA-DQB1*06:42',
         'HLA-DQA1*01:05-HLA-DQB1*06:43', 'HLA-DQA1*01:05-HLA-DQB1*06:44', 'HLA-DQA1*01:06-HLA-DQB1*02:01', 'HLA-DQA1*01:06-HLA-DQB1*02:02',
         'HLA-DQA1*01:06-HLA-DQB1*02:03', 'HLA-DQA1*01:06-HLA-DQB1*02:04',
         'HLA-DQA1*01:06-HLA-DQB1*02:05', 'HLA-DQA1*01:06-HLA-DQB1*02:06', 'HLA-DQA1*01:06-HLA-DQB1*03:01', 'HLA-DQA1*01:06-HLA-DQB1*03:02',
         'HLA-DQA1*01:06-HLA-DQB1*03:03', 'HLA-DQA1*01:06-HLA-DQB1*03:04',
         'HLA-DQA1*01:06-HLA-DQB1*03:05', 'HLA-DQA1*01:06-HLA-DQB1*03:06', 'HLA-DQA1*01:06-HLA-DQB1*03:07', 'HLA-DQA1*01:06-HLA-DQB1*03:08',
         'HLA-DQA1*01:06-HLA-DQB1*03:09', 'HLA-DQA1*01:06-HLA-DQB1*03:10',
         'HLA-DQA1*01:06-HLA-DQB1*03:11', 'HLA-DQA1*01:06-HLA-DQB1*03:12', 'HLA-DQA1*01:06-HLA-DQB1*03:13', 'HLA-DQA1*01:06-HLA-DQB1*03:14',
         'HLA-DQA1*01:06-HLA-DQB1*03:15', 'HLA-DQA1*01:06-HLA-DQB1*03:16',
         'HLA-DQA1*01:06-HLA-DQB1*03:17', 'HLA-DQA1*01:06-HLA-DQB1*03:18', 'HLA-DQA1*01:06-HLA-DQB1*03:19', 'HLA-DQA1*01:06-HLA-DQB1*03:20',
         'HLA-DQA1*01:06-HLA-DQB1*03:21', 'HLA-DQA1*01:06-HLA-DQB1*03:22',
         'HLA-DQA1*01:06-HLA-DQB1*03:23', 'HLA-DQA1*01:06-HLA-DQB1*03:24', 'HLA-DQA1*01:06-HLA-DQB1*03:25', 'HLA-DQA1*01:06-HLA-DQB1*03:26',
         'HLA-DQA1*01:06-HLA-DQB1*03:27', 'HLA-DQA1*01:06-HLA-DQB1*03:28',
         'HLA-DQA1*01:06-HLA-DQB1*03:29', 'HLA-DQA1*01:06-HLA-DQB1*03:30', 'HLA-DQA1*01:06-HLA-DQB1*03:31', 'HLA-DQA1*01:06-HLA-DQB1*03:32',
         'HLA-DQA1*01:06-HLA-DQB1*03:33', 'HLA-DQA1*01:06-HLA-DQB1*03:34',
         'HLA-DQA1*01:06-HLA-DQB1*03:35', 'HLA-DQA1*01:06-HLA-DQB1*03:36', 'HLA-DQA1*01:06-HLA-DQB1*03:37', 'HLA-DQA1*01:06-HLA-DQB1*03:38',
         'HLA-DQA1*01:06-HLA-DQB1*04:01', 'HLA-DQA1*01:06-HLA-DQB1*04:02',
         'HLA-DQA1*01:06-HLA-DQB1*04:03', 'HLA-DQA1*01:06-HLA-DQB1*04:04', 'HLA-DQA1*01:06-HLA-DQB1*04:05', 'HLA-DQA1*01:06-HLA-DQB1*04:06',
         'HLA-DQA1*01:06-HLA-DQB1*04:07', 'HLA-DQA1*01:06-HLA-DQB1*04:08',
         'HLA-DQA1*01:06-HLA-DQB1*05:01', 'HLA-DQA1*01:06-HLA-DQB1*05:02', 'HLA-DQA1*01:06-HLA-DQB1*05:03', 'HLA-DQA1*01:06-HLA-DQB1*05:05',
         'HLA-DQA1*01:06-HLA-DQB1*05:06', 'HLA-DQA1*01:06-HLA-DQB1*05:07',
         'HLA-DQA1*01:06-HLA-DQB1*05:08', 'HLA-DQA1*01:06-HLA-DQB1*05:09', 'HLA-DQA1*01:06-HLA-DQB1*05:10', 'HLA-DQA1*01:06-HLA-DQB1*05:11',
         'HLA-DQA1*01:06-HLA-DQB1*05:12', 'HLA-DQA1*01:06-HLA-DQB1*05:13',
         'HLA-DQA1*01:06-HLA-DQB1*05:14', 'HLA-DQA1*01:06-HLA-DQB1*06:01', 'HLA-DQA1*01:06-HLA-DQB1*06:02', 'HLA-DQA1*01:06-HLA-DQB1*06:03',
         'HLA-DQA1*01:06-HLA-DQB1*06:04', 'HLA-DQA1*01:06-HLA-DQB1*06:07',
         'HLA-DQA1*01:06-HLA-DQB1*06:08', 'HLA-DQA1*01:06-HLA-DQB1*06:09', 'HLA-DQA1*01:06-HLA-DQB1*06:10', 'HLA-DQA1*01:06-HLA-DQB1*06:11',
         'HLA-DQA1*01:06-HLA-DQB1*06:12', 'HLA-DQA1*01:06-HLA-DQB1*06:14',
         'HLA-DQA1*01:06-HLA-DQB1*06:15', 'HLA-DQA1*01:06-HLA-DQB1*06:16', 'HLA-DQA1*01:06-HLA-DQB1*06:17', 'HLA-DQA1*01:06-HLA-DQB1*06:18',
         'HLA-DQA1*01:06-HLA-DQB1*06:19', 'HLA-DQA1*01:06-HLA-DQB1*06:21',
         'HLA-DQA1*01:06-HLA-DQB1*06:22', 'HLA-DQA1*01:06-HLA-DQB1*06:23', 'HLA-DQA1*01:06-HLA-DQB1*06:24', 'HLA-DQA1*01:06-HLA-DQB1*06:25',
         'HLA-DQA1*01:06-HLA-DQB1*06:27', 'HLA-DQA1*01:06-HLA-DQB1*06:28',
         'HLA-DQA1*01:06-HLA-DQB1*06:29', 'HLA-DQA1*01:06-HLA-DQB1*06:30', 'HLA-DQA1*01:06-HLA-DQB1*06:31', 'HLA-DQA1*01:06-HLA-DQB1*06:32',
         'HLA-DQA1*01:06-HLA-DQB1*06:33', 'HLA-DQA1*01:06-HLA-DQB1*06:34',
         'HLA-DQA1*01:06-HLA-DQB1*06:35', 'HLA-DQA1*01:06-HLA-DQB1*06:36', 'HLA-DQA1*01:06-HLA-DQB1*06:37', 'HLA-DQA1*01:06-HLA-DQB1*06:38',
         'HLA-DQA1*01:06-HLA-DQB1*06:39', 'HLA-DQA1*01:06-HLA-DQB1*06:40',
         'HLA-DQA1*01:06-HLA-DQB1*06:41', 'HLA-DQA1*01:06-HLA-DQB1*06:42', 'HLA-DQA1*01:06-HLA-DQB1*06:43', 'HLA-DQA1*01:06-HLA-DQB1*06:44',
         'HLA-DQA1*01:07-HLA-DQB1*02:01', 'HLA-DQA1*01:07-HLA-DQB1*02:02',
         'HLA-DQA1*01:07-HLA-DQB1*02:03', 'HLA-DQA1*01:07-HLA-DQB1*02:04', 'HLA-DQA1*01:07-HLA-DQB1*02:05', 'HLA-DQA1*01:07-HLA-DQB1*02:06',
         'HLA-DQA1*01:07-HLA-DQB1*03:01', 'HLA-DQA1*01:07-HLA-DQB1*03:02',
         'HLA-DQA1*01:07-HLA-DQB1*03:03', 'HLA-DQA1*01:07-HLA-DQB1*03:04', 'HLA-DQA1*01:07-HLA-DQB1*03:05', 'HLA-DQA1*01:07-HLA-DQB1*03:06',
         'HLA-DQA1*01:07-HLA-DQB1*03:07', 'HLA-DQA1*01:07-HLA-DQB1*03:08',
         'HLA-DQA1*01:07-HLA-DQB1*03:09', 'HLA-DQA1*01:07-HLA-DQB1*03:10', 'HLA-DQA1*01:07-HLA-DQB1*03:11', 'HLA-DQA1*01:07-HLA-DQB1*03:12',
         'HLA-DQA1*01:07-HLA-DQB1*03:13', 'HLA-DQA1*01:07-HLA-DQB1*03:14',
         'HLA-DQA1*01:07-HLA-DQB1*03:15', 'HLA-DQA1*01:07-HLA-DQB1*03:16', 'HLA-DQA1*01:07-HLA-DQB1*03:17', 'HLA-DQA1*01:07-HLA-DQB1*03:18',
         'HLA-DQA1*01:07-HLA-DQB1*03:19', 'HLA-DQA1*01:07-HLA-DQB1*03:20',
         'HLA-DQA1*01:07-HLA-DQB1*03:21', 'HLA-DQA1*01:07-HLA-DQB1*03:22', 'HLA-DQA1*01:07-HLA-DQB1*03:23', 'HLA-DQA1*01:07-HLA-DQB1*03:24',
         'HLA-DQA1*01:07-HLA-DQB1*03:25', 'HLA-DQA1*01:07-HLA-DQB1*03:26',
         'HLA-DQA1*01:07-HLA-DQB1*03:27', 'HLA-DQA1*01:07-HLA-DQB1*03:28', 'HLA-DQA1*01:07-HLA-DQB1*03:29', 'HLA-DQA1*01:07-HLA-DQB1*03:30',
         'HLA-DQA1*01:07-HLA-DQB1*03:31', 'HLA-DQA1*01:07-HLA-DQB1*03:32',
         'HLA-DQA1*01:07-HLA-DQB1*03:33', 'HLA-DQA1*01:07-HLA-DQB1*03:34', 'HLA-DQA1*01:07-HLA-DQB1*03:35', 'HLA-DQA1*01:07-HLA-DQB1*03:36',
         'HLA-DQA1*01:07-HLA-DQB1*03:37', 'HLA-DQA1*01:07-HLA-DQB1*03:38',
         'HLA-DQA1*01:07-HLA-DQB1*04:01', 'HLA-DQA1*01:07-HLA-DQB1*04:02', 'HLA-DQA1*01:07-HLA-DQB1*04:03', 'HLA-DQA1*01:07-HLA-DQB1*04:04',
         'HLA-DQA1*01:07-HLA-DQB1*04:05', 'HLA-DQA1*01:07-HLA-DQB1*04:06',
         'HLA-DQA1*01:07-HLA-DQB1*04:07', 'HLA-DQA1*01:07-HLA-DQB1*04:08', 'HLA-DQA1*01:07-HLA-DQB1*05:01', 'HLA-DQA1*01:07-HLA-DQB1*05:02',
         'HLA-DQA1*01:07-HLA-DQB1*05:03', 'HLA-DQA1*01:07-HLA-DQB1*05:05',
         'HLA-DQA1*01:07-HLA-DQB1*05:06', 'HLA-DQA1*01:07-HLA-DQB1*05:07', 'HLA-DQA1*01:07-HLA-DQB1*05:08', 'HLA-DQA1*01:07-HLA-DQB1*05:09',
         'HLA-DQA1*01:07-HLA-DQB1*05:10', 'HLA-DQA1*01:07-HLA-DQB1*05:11',
         'HLA-DQA1*01:07-HLA-DQB1*05:12', 'HLA-DQA1*01:07-HLA-DQB1*05:13', 'HLA-DQA1*01:07-HLA-DQB1*05:14', 'HLA-DQA1*01:07-HLA-DQB1*06:01',
         'HLA-DQA1*01:07-HLA-DQB1*06:02', 'HLA-DQA1*01:07-HLA-DQB1*06:03',
         'HLA-DQA1*01:07-HLA-DQB1*06:04', 'HLA-DQA1*01:07-HLA-DQB1*06:07', 'HLA-DQA1*01:07-HLA-DQB1*06:08', 'HLA-DQA1*01:07-HLA-DQB1*06:09',
         'HLA-DQA1*01:07-HLA-DQB1*06:10', 'HLA-DQA1*01:07-HLA-DQB1*06:11',
         'HLA-DQA1*01:07-HLA-DQB1*06:12', 'HLA-DQA1*01:07-HLA-DQB1*06:14', 'HLA-DQA1*01:07-HLA-DQB1*06:15', 'HLA-DQA1*01:07-HLA-DQB1*06:16',
         'HLA-DQA1*01:07-HLA-DQB1*06:17', 'HLA-DQA1*01:07-HLA-DQB1*06:18',
         'HLA-DQA1*01:07-HLA-DQB1*06:19', 'HLA-DQA1*01:07-HLA-DQB1*06:21', 'HLA-DQA1*01:07-HLA-DQB1*06:22', 'HLA-DQA1*01:07-HLA-DQB1*06:23',
         'HLA-DQA1*01:07-HLA-DQB1*06:24', 'HLA-DQA1*01:07-HLA-DQB1*06:25',
         'HLA-DQA1*01:07-HLA-DQB1*06:27', 'HLA-DQA1*01:07-HLA-DQB1*06:28', 'HLA-DQA1*01:07-HLA-DQB1*06:29', 'HLA-DQA1*01:07-HLA-DQB1*06:30',
         'HLA-DQA1*01:07-HLA-DQB1*06:31', 'HLA-DQA1*01:07-HLA-DQB1*06:32',
         'HLA-DQA1*01:07-HLA-DQB1*06:33', 'HLA-DQA1*01:07-HLA-DQB1*06:34', 'HLA-DQA1*01:07-HLA-DQB1*06:35', 'HLA-DQA1*01:07-HLA-DQB1*06:36',
         'HLA-DQA1*01:07-HLA-DQB1*06:37', 'HLA-DQA1*01:07-HLA-DQB1*06:38',
         'HLA-DQA1*01:07-HLA-DQB1*06:39', 'HLA-DQA1*01:07-HLA-DQB1*06:40', 'HLA-DQA1*01:07-HLA-DQB1*06:41', 'HLA-DQA1*01:07-HLA-DQB1*06:42',
         'HLA-DQA1*01:07-HLA-DQB1*06:43', 'HLA-DQA1*01:07-HLA-DQB1*06:44',
         'HLA-DQA1*01:08-HLA-DQB1*02:01', 'HLA-DQA1*01:08-HLA-DQB1*02:02', 'HLA-DQA1*01:08-HLA-DQB1*02:03', 'HLA-DQA1*01:08-HLA-DQB1*02:04',
         'HLA-DQA1*01:08-HLA-DQB1*02:05', 'HLA-DQA1*01:08-HLA-DQB1*02:06',
         'HLA-DQA1*01:08-HLA-DQB1*03:01', 'HLA-DQA1*01:08-HLA-DQB1*03:02', 'HLA-DQA1*01:08-HLA-DQB1*03:03', 'HLA-DQA1*01:08-HLA-DQB1*03:04',
         'HLA-DQA1*01:08-HLA-DQB1*03:05', 'HLA-DQA1*01:08-HLA-DQB1*03:06',
         'HLA-DQA1*01:08-HLA-DQB1*03:07', 'HLA-DQA1*01:08-HLA-DQB1*03:08', 'HLA-DQA1*01:08-HLA-DQB1*03:09', 'HLA-DQA1*01:08-HLA-DQB1*03:10',
         'HLA-DQA1*01:08-HLA-DQB1*03:11', 'HLA-DQA1*01:08-HLA-DQB1*03:12',
         'HLA-DQA1*01:08-HLA-DQB1*03:13', 'HLA-DQA1*01:08-HLA-DQB1*03:14', 'HLA-DQA1*01:08-HLA-DQB1*03:15', 'HLA-DQA1*01:08-HLA-DQB1*03:16',
         'HLA-DQA1*01:08-HLA-DQB1*03:17', 'HLA-DQA1*01:08-HLA-DQB1*03:18',
         'HLA-DQA1*01:08-HLA-DQB1*03:19', 'HLA-DQA1*01:08-HLA-DQB1*03:20', 'HLA-DQA1*01:08-HLA-DQB1*03:21', 'HLA-DQA1*01:08-HLA-DQB1*03:22',
         'HLA-DQA1*01:08-HLA-DQB1*03:23', 'HLA-DQA1*01:08-HLA-DQB1*03:24',
         'HLA-DQA1*01:08-HLA-DQB1*03:25', 'HLA-DQA1*01:08-HLA-DQB1*03:26', 'HLA-DQA1*01:08-HLA-DQB1*03:27', 'HLA-DQA1*01:08-HLA-DQB1*03:28',
         'HLA-DQA1*01:08-HLA-DQB1*03:29', 'HLA-DQA1*01:08-HLA-DQB1*03:30',
         'HLA-DQA1*01:08-HLA-DQB1*03:31', 'HLA-DQA1*01:08-HLA-DQB1*03:32', 'HLA-DQA1*01:08-HLA-DQB1*03:33', 'HLA-DQA1*01:08-HLA-DQB1*03:34',
         'HLA-DQA1*01:08-HLA-DQB1*03:35', 'HLA-DQA1*01:08-HLA-DQB1*03:36',
         'HLA-DQA1*01:08-HLA-DQB1*03:37', 'HLA-DQA1*01:08-HLA-DQB1*03:38', 'HLA-DQA1*01:08-HLA-DQB1*04:01', 'HLA-DQA1*01:08-HLA-DQB1*04:02',
         'HLA-DQA1*01:08-HLA-DQB1*04:03', 'HLA-DQA1*01:08-HLA-DQB1*04:04',
         'HLA-DQA1*01:08-HLA-DQB1*04:05', 'HLA-DQA1*01:08-HLA-DQB1*04:06', 'HLA-DQA1*01:08-HLA-DQB1*04:07', 'HLA-DQA1*01:08-HLA-DQB1*04:08',
         'HLA-DQA1*01:08-HLA-DQB1*05:01', 'HLA-DQA1*01:08-HLA-DQB1*05:02',
         'HLA-DQA1*01:08-HLA-DQB1*05:03', 'HLA-DQA1*01:08-HLA-DQB1*05:05', 'HLA-DQA1*01:08-HLA-DQB1*05:06', 'HLA-DQA1*01:08-HLA-DQB1*05:07',
         'HLA-DQA1*01:08-HLA-DQB1*05:08', 'HLA-DQA1*01:08-HLA-DQB1*05:09',
         'HLA-DQA1*01:08-HLA-DQB1*05:10', 'HLA-DQA1*01:08-HLA-DQB1*05:11', 'HLA-DQA1*01:08-HLA-DQB1*05:12', 'HLA-DQA1*01:08-HLA-DQB1*05:13',
         'HLA-DQA1*01:08-HLA-DQB1*05:14', 'HLA-DQA1*01:08-HLA-DQB1*06:01',
         'HLA-DQA1*01:08-HLA-DQB1*06:02', 'HLA-DQA1*01:08-HLA-DQB1*06:03', 'HLA-DQA1*01:08-HLA-DQB1*06:04', 'HLA-DQA1*01:08-HLA-DQB1*06:07',
         'HLA-DQA1*01:08-HLA-DQB1*06:08', 'HLA-DQA1*01:08-HLA-DQB1*06:09',
         'HLA-DQA1*01:08-HLA-DQB1*06:10', 'HLA-DQA1*01:08-HLA-DQB1*06:11', 'HLA-DQA1*01:08-HLA-DQB1*06:12', 'HLA-DQA1*01:08-HLA-DQB1*06:14',
         'HLA-DQA1*01:08-HLA-DQB1*06:15', 'HLA-DQA1*01:08-HLA-DQB1*06:16',
         'HLA-DQA1*01:08-HLA-DQB1*06:17', 'HLA-DQA1*01:08-HLA-DQB1*06:18', 'HLA-DQA1*01:08-HLA-DQB1*06:19', 'HLA-DQA1*01:08-HLA-DQB1*06:21',
         'HLA-DQA1*01:08-HLA-DQB1*06:22', 'HLA-DQA1*01:08-HLA-DQB1*06:23',
         'HLA-DQA1*01:08-HLA-DQB1*06:24', 'HLA-DQA1*01:08-HLA-DQB1*06:25', 'HLA-DQA1*01:08-HLA-DQB1*06:27', 'HLA-DQA1*01:08-HLA-DQB1*06:28',
         'HLA-DQA1*01:08-HLA-DQB1*06:29', 'HLA-DQA1*01:08-HLA-DQB1*06:30',
         'HLA-DQA1*01:08-HLA-DQB1*06:31', 'HLA-DQA1*01:08-HLA-DQB1*06:32', 'HLA-DQA1*01:08-HLA-DQB1*06:33', 'HLA-DQA1*01:08-HLA-DQB1*06:34',
         'HLA-DQA1*01:08-HLA-DQB1*06:35', 'HLA-DQA1*01:08-HLA-DQB1*06:36',
         'HLA-DQA1*01:08-HLA-DQB1*06:37', 'HLA-DQA1*01:08-HLA-DQB1*06:38', 'HLA-DQA1*01:08-HLA-DQB1*06:39', 'HLA-DQA1*01:08-HLA-DQB1*06:40',
         'HLA-DQA1*01:08-HLA-DQB1*06:41', 'HLA-DQA1*01:08-HLA-DQB1*06:42',
         'HLA-DQA1*01:08-HLA-DQB1*06:43', 'HLA-DQA1*01:08-HLA-DQB1*06:44', 'HLA-DQA1*01:09-HLA-DQB1*02:01', 'HLA-DQA1*01:09-HLA-DQB1*02:02',
         'HLA-DQA1*01:09-HLA-DQB1*02:03', 'HLA-DQA1*01:09-HLA-DQB1*02:04',
         'HLA-DQA1*01:09-HLA-DQB1*02:05', 'HLA-DQA1*01:09-HLA-DQB1*02:06', 'HLA-DQA1*01:09-HLA-DQB1*03:01', 'HLA-DQA1*01:09-HLA-DQB1*03:02',
         'HLA-DQA1*01:09-HLA-DQB1*03:03', 'HLA-DQA1*01:09-HLA-DQB1*03:04',
         'HLA-DQA1*01:09-HLA-DQB1*03:05', 'HLA-DQA1*01:09-HLA-DQB1*03:06', 'HLA-DQA1*01:09-HLA-DQB1*03:07', 'HLA-DQA1*01:09-HLA-DQB1*03:08',
         'HLA-DQA1*01:09-HLA-DQB1*03:09', 'HLA-DQA1*01:09-HLA-DQB1*03:10',
         'HLA-DQA1*01:09-HLA-DQB1*03:11', 'HLA-DQA1*01:09-HLA-DQB1*03:12', 'HLA-DQA1*01:09-HLA-DQB1*03:13', 'HLA-DQA1*01:09-HLA-DQB1*03:14',
         'HLA-DQA1*01:09-HLA-DQB1*03:15', 'HLA-DQA1*01:09-HLA-DQB1*03:16',
         'HLA-DQA1*01:09-HLA-DQB1*03:17', 'HLA-DQA1*01:09-HLA-DQB1*03:18', 'HLA-DQA1*01:09-HLA-DQB1*03:19', 'HLA-DQA1*01:09-HLA-DQB1*03:20',
         'HLA-DQA1*01:09-HLA-DQB1*03:21', 'HLA-DQA1*01:09-HLA-DQB1*03:22',
         'HLA-DQA1*01:09-HLA-DQB1*03:23', 'HLA-DQA1*01:09-HLA-DQB1*03:24', 'HLA-DQA1*01:09-HLA-DQB1*03:25', 'HLA-DQA1*01:09-HLA-DQB1*03:26',
         'HLA-DQA1*01:09-HLA-DQB1*03:27', 'HLA-DQA1*01:09-HLA-DQB1*03:28',
         'HLA-DQA1*01:09-HLA-DQB1*03:29', 'HLA-DQA1*01:09-HLA-DQB1*03:30', 'HLA-DQA1*01:09-HLA-DQB1*03:31', 'HLA-DQA1*01:09-HLA-DQB1*03:32',
         'HLA-DQA1*01:09-HLA-DQB1*03:33', 'HLA-DQA1*01:09-HLA-DQB1*03:34',
         'HLA-DQA1*01:09-HLA-DQB1*03:35', 'HLA-DQA1*01:09-HLA-DQB1*03:36', 'HLA-DQA1*01:09-HLA-DQB1*03:37', 'HLA-DQA1*01:09-HLA-DQB1*03:38',
         'HLA-DQA1*01:09-HLA-DQB1*04:01', 'HLA-DQA1*01:09-HLA-DQB1*04:02',
         'HLA-DQA1*01:09-HLA-DQB1*04:03', 'HLA-DQA1*01:09-HLA-DQB1*04:04', 'HLA-DQA1*01:09-HLA-DQB1*04:05', 'HLA-DQA1*01:09-HLA-DQB1*04:06',
         'HLA-DQA1*01:09-HLA-DQB1*04:07', 'HLA-DQA1*01:09-HLA-DQB1*04:08',
         'HLA-DQA1*01:09-HLA-DQB1*05:01', 'HLA-DQA1*01:09-HLA-DQB1*05:02', 'HLA-DQA1*01:09-HLA-DQB1*05:03', 'HLA-DQA1*01:09-HLA-DQB1*05:05',
         'HLA-DQA1*01:09-HLA-DQB1*05:06', 'HLA-DQA1*01:09-HLA-DQB1*05:07',
         'HLA-DQA1*01:09-HLA-DQB1*05:08', 'HLA-DQA1*01:09-HLA-DQB1*05:09', 'HLA-DQA1*01:09-HLA-DQB1*05:10', 'HLA-DQA1*01:09-HLA-DQB1*05:11',
         'HLA-DQA1*01:09-HLA-DQB1*05:12', 'HLA-DQA1*01:09-HLA-DQB1*05:13',
         'HLA-DQA1*01:09-HLA-DQB1*05:14', 'HLA-DQA1*01:09-HLA-DQB1*06:01', 'HLA-DQA1*01:09-HLA-DQB1*06:02', 'HLA-DQA1*01:09-HLA-DQB1*06:03',
         'HLA-DQA1*01:09-HLA-DQB1*06:04', 'HLA-DQA1*01:09-HLA-DQB1*06:07',
         'HLA-DQA1*01:09-HLA-DQB1*06:08', 'HLA-DQA1*01:09-HLA-DQB1*06:09', 'HLA-DQA1*01:09-HLA-DQB1*06:10', 'HLA-DQA1*01:09-HLA-DQB1*06:11',
         'HLA-DQA1*01:09-HLA-DQB1*06:12', 'HLA-DQA1*01:09-HLA-DQB1*06:14',
         'HLA-DQA1*01:09-HLA-DQB1*06:15', 'HLA-DQA1*01:09-HLA-DQB1*06:16', 'HLA-DQA1*01:09-HLA-DQB1*06:17', 'HLA-DQA1*01:09-HLA-DQB1*06:18',
         'HLA-DQA1*01:09-HLA-DQB1*06:19', 'HLA-DQA1*01:09-HLA-DQB1*06:21',
         'HLA-DQA1*01:09-HLA-DQB1*06:22', 'HLA-DQA1*01:09-HLA-DQB1*06:23', 'HLA-DQA1*01:09-HLA-DQB1*06:24', 'HLA-DQA1*01:09-HLA-DQB1*06:25',
         'HLA-DQA1*01:09-HLA-DQB1*06:27', 'HLA-DQA1*01:09-HLA-DQB1*06:28',
         'HLA-DQA1*01:09-HLA-DQB1*06:29', 'HLA-DQA1*01:09-HLA-DQB1*06:30', 'HLA-DQA1*01:09-HLA-DQB1*06:31', 'HLA-DQA1*01:09-HLA-DQB1*06:32',
         'HLA-DQA1*01:09-HLA-DQB1*06:33', 'HLA-DQA1*01:09-HLA-DQB1*06:34',
         'HLA-DQA1*01:09-HLA-DQB1*06:35', 'HLA-DQA1*01:09-HLA-DQB1*06:36', 'HLA-DQA1*01:09-HLA-DQB1*06:37', 'HLA-DQA1*01:09-HLA-DQB1*06:38',
         'HLA-DQA1*01:09-HLA-DQB1*06:39', 'HLA-DQA1*01:09-HLA-DQB1*06:40',
         'HLA-DQA1*01:09-HLA-DQB1*06:41', 'HLA-DQA1*01:09-HLA-DQB1*06:42', 'HLA-DQA1*01:09-HLA-DQB1*06:43', 'HLA-DQA1*01:09-HLA-DQB1*06:44',
         'HLA-DQA1*02:01-HLA-DQB1*02:01', 'HLA-DQA1*02:01-HLA-DQB1*02:02',
         'HLA-DQA1*02:01-HLA-DQB1*02:03', 'HLA-DQA1*02:01-HLA-DQB1*02:04', 'HLA-DQA1*02:01-HLA-DQB1*02:05', 'HLA-DQA1*02:01-HLA-DQB1*02:06',
         'HLA-DQA1*02:01-HLA-DQB1*03:01', 'HLA-DQA1*02:01-HLA-DQB1*03:02',
         'HLA-DQA1*02:01-HLA-DQB1*03:03', 'HLA-DQA1*02:01-HLA-DQB1*03:04', 'HLA-DQA1*02:01-HLA-DQB1*03:05', 'HLA-DQA1*02:01-HLA-DQB1*03:06',
         'HLA-DQA1*02:01-HLA-DQB1*03:07', 'HLA-DQA1*02:01-HLA-DQB1*03:08',
         'HLA-DQA1*02:01-HLA-DQB1*03:09', 'HLA-DQA1*02:01-HLA-DQB1*03:10', 'HLA-DQA1*02:01-HLA-DQB1*03:11', 'HLA-DQA1*02:01-HLA-DQB1*03:12',
         'HLA-DQA1*02:01-HLA-DQB1*03:13', 'HLA-DQA1*02:01-HLA-DQB1*03:14',
         'HLA-DQA1*02:01-HLA-DQB1*03:15', 'HLA-DQA1*02:01-HLA-DQB1*03:16', 'HLA-DQA1*02:01-HLA-DQB1*03:17', 'HLA-DQA1*02:01-HLA-DQB1*03:18',
         'HLA-DQA1*02:01-HLA-DQB1*03:19', 'HLA-DQA1*02:01-HLA-DQB1*03:20',
         'HLA-DQA1*02:01-HLA-DQB1*03:21', 'HLA-DQA1*02:01-HLA-DQB1*03:22', 'HLA-DQA1*02:01-HLA-DQB1*03:23', 'HLA-DQA1*02:01-HLA-DQB1*03:24',
         'HLA-DQA1*02:01-HLA-DQB1*03:25', 'HLA-DQA1*02:01-HLA-DQB1*03:26',
         'HLA-DQA1*02:01-HLA-DQB1*03:27', 'HLA-DQA1*02:01-HLA-DQB1*03:28', 'HLA-DQA1*02:01-HLA-DQB1*03:29', 'HLA-DQA1*02:01-HLA-DQB1*03:30',
         'HLA-DQA1*02:01-HLA-DQB1*03:31', 'HLA-DQA1*02:01-HLA-DQB1*03:32',
         'HLA-DQA1*02:01-HLA-DQB1*03:33', 'HLA-DQA1*02:01-HLA-DQB1*03:34', 'HLA-DQA1*02:01-HLA-DQB1*03:35', 'HLA-DQA1*02:01-HLA-DQB1*03:36',
         'HLA-DQA1*02:01-HLA-DQB1*03:37', 'HLA-DQA1*02:01-HLA-DQB1*03:38',
         'HLA-DQA1*02:01-HLA-DQB1*04:01', 'HLA-DQA1*02:01-HLA-DQB1*04:02', 'HLA-DQA1*02:01-HLA-DQB1*04:03', 'HLA-DQA1*02:01-HLA-DQB1*04:04',
         'HLA-DQA1*02:01-HLA-DQB1*04:05', 'HLA-DQA1*02:01-HLA-DQB1*04:06',
         'HLA-DQA1*02:01-HLA-DQB1*04:07', 'HLA-DQA1*02:01-HLA-DQB1*04:08', 'HLA-DQA1*02:01-HLA-DQB1*05:01', 'HLA-DQA1*02:01-HLA-DQB1*05:02',
         'HLA-DQA1*02:01-HLA-DQB1*05:03', 'HLA-DQA1*02:01-HLA-DQB1*05:05',
         'HLA-DQA1*02:01-HLA-DQB1*05:06', 'HLA-DQA1*02:01-HLA-DQB1*05:07', 'HLA-DQA1*02:01-HLA-DQB1*05:08', 'HLA-DQA1*02:01-HLA-DQB1*05:09',
         'HLA-DQA1*02:01-HLA-DQB1*05:10', 'HLA-DQA1*02:01-HLA-DQB1*05:11',
         'HLA-DQA1*02:01-HLA-DQB1*05:12', 'HLA-DQA1*02:01-HLA-DQB1*05:13', 'HLA-DQA1*02:01-HLA-DQB1*05:14', 'HLA-DQA1*02:01-HLA-DQB1*06:01',
         'HLA-DQA1*02:01-HLA-DQB1*06:02', 'HLA-DQA1*02:01-HLA-DQB1*06:03',
         'HLA-DQA1*02:01-HLA-DQB1*06:04', 'HLA-DQA1*02:01-HLA-DQB1*06:07', 'HLA-DQA1*02:01-HLA-DQB1*06:08', 'HLA-DQA1*02:01-HLA-DQB1*06:09',
         'HLA-DQA1*02:01-HLA-DQB1*06:10', 'HLA-DQA1*02:01-HLA-DQB1*06:11',
         'HLA-DQA1*02:01-HLA-DQB1*06:12', 'HLA-DQA1*02:01-HLA-DQB1*06:14', 'HLA-DQA1*02:01-HLA-DQB1*06:15', 'HLA-DQA1*02:01-HLA-DQB1*06:16',
         'HLA-DQA1*02:01-HLA-DQB1*06:17', 'HLA-DQA1*02:01-HLA-DQB1*06:18',
         'HLA-DQA1*02:01-HLA-DQB1*06:19', 'HLA-DQA1*02:01-HLA-DQB1*06:21', 'HLA-DQA1*02:01-HLA-DQB1*06:22', 'HLA-DQA1*02:01-HLA-DQB1*06:23',
         'HLA-DQA1*02:01-HLA-DQB1*06:24', 'HLA-DQA1*02:01-HLA-DQB1*06:25',
         'HLA-DQA1*02:01-HLA-DQB1*06:27', 'HLA-DQA1*02:01-HLA-DQB1*06:28', 'HLA-DQA1*02:01-HLA-DQB1*06:29', 'HLA-DQA1*02:01-HLA-DQB1*06:30',
         'HLA-DQA1*02:01-HLA-DQB1*06:31', 'HLA-DQA1*02:01-HLA-DQB1*06:32',
         'HLA-DQA1*02:01-HLA-DQB1*06:33', 'HLA-DQA1*02:01-HLA-DQB1*06:34', 'HLA-DQA1*02:01-HLA-DQB1*06:35', 'HLA-DQA1*02:01-HLA-DQB1*06:36',
         'HLA-DQA1*02:01-HLA-DQB1*06:37', 'HLA-DQA1*02:01-HLA-DQB1*06:38',
         'HLA-DQA1*02:01-HLA-DQB1*06:39', 'HLA-DQA1*02:01-HLA-DQB1*06:40', 'HLA-DQA1*02:01-HLA-DQB1*06:41', 'HLA-DQA1*02:01-HLA-DQB1*06:42',
         'HLA-DQA1*02:01-HLA-DQB1*06:43', 'HLA-DQA1*02:01-HLA-DQB1*06:44',
         'HLA-DQA1*03:01-HLA-DQB1*02:01', 'HLA-DQA1*03:01-HLA-DQB1*02:02', 'HLA-DQA1*03:01-HLA-DQB1*02:03', 'HLA-DQA1*03:01-HLA-DQB1*02:04',
         'HLA-DQA1*03:01-HLA-DQB1*02:05', 'HLA-DQA1*03:01-HLA-DQB1*02:06',
         'HLA-DQA1*03:01-HLA-DQB1*03:01', 'HLA-DQA1*03:01-HLA-DQB1*03:02', 'HLA-DQA1*03:01-HLA-DQB1*03:03', 'HLA-DQA1*03:01-HLA-DQB1*03:04',
         'HLA-DQA1*03:01-HLA-DQB1*03:05', 'HLA-DQA1*03:01-HLA-DQB1*03:06',
         'HLA-DQA1*03:01-HLA-DQB1*03:07', 'HLA-DQA1*03:01-HLA-DQB1*03:08', 'HLA-DQA1*03:01-HLA-DQB1*03:09', 'HLA-DQA1*03:01-HLA-DQB1*03:10',
         'HLA-DQA1*03:01-HLA-DQB1*03:11', 'HLA-DQA1*03:01-HLA-DQB1*03:12',
         'HLA-DQA1*03:01-HLA-DQB1*03:13', 'HLA-DQA1*03:01-HLA-DQB1*03:14', 'HLA-DQA1*03:01-HLA-DQB1*03:15', 'HLA-DQA1*03:01-HLA-DQB1*03:16',
         'HLA-DQA1*03:01-HLA-DQB1*03:17', 'HLA-DQA1*03:01-HLA-DQB1*03:18',
         'HLA-DQA1*03:01-HLA-DQB1*03:19', 'HLA-DQA1*03:01-HLA-DQB1*03:20', 'HLA-DQA1*03:01-HLA-DQB1*03:21', 'HLA-DQA1*03:01-HLA-DQB1*03:22',
         'HLA-DQA1*03:01-HLA-DQB1*03:23', 'HLA-DQA1*03:01-HLA-DQB1*03:24',
         'HLA-DQA1*03:01-HLA-DQB1*03:25', 'HLA-DQA1*03:01-HLA-DQB1*03:26', 'HLA-DQA1*03:01-HLA-DQB1*03:27', 'HLA-DQA1*03:01-HLA-DQB1*03:28',
         'HLA-DQA1*03:01-HLA-DQB1*03:29', 'HLA-DQA1*03:01-HLA-DQB1*03:30',
         'HLA-DQA1*03:01-HLA-DQB1*03:31', 'HLA-DQA1*03:01-HLA-DQB1*03:32', 'HLA-DQA1*03:01-HLA-DQB1*03:33', 'HLA-DQA1*03:01-HLA-DQB1*03:34',
         'HLA-DQA1*03:01-HLA-DQB1*03:35', 'HLA-DQA1*03:01-HLA-DQB1*03:36',
         'HLA-DQA1*03:01-HLA-DQB1*03:37', 'HLA-DQA1*03:01-HLA-DQB1*03:38', 'HLA-DQA1*03:01-HLA-DQB1*04:01', 'HLA-DQA1*03:01-HLA-DQB1*04:02',
         'HLA-DQA1*03:01-HLA-DQB1*04:03', 'HLA-DQA1*03:01-HLA-DQB1*04:04',
         'HLA-DQA1*03:01-HLA-DQB1*04:05', 'HLA-DQA1*03:01-HLA-DQB1*04:06', 'HLA-DQA1*03:01-HLA-DQB1*04:07', 'HLA-DQA1*03:01-HLA-DQB1*04:08',
         'HLA-DQA1*03:01-HLA-DQB1*05:01', 'HLA-DQA1*03:01-HLA-DQB1*05:02',
         'HLA-DQA1*03:01-HLA-DQB1*05:03', 'HLA-DQA1*03:01-HLA-DQB1*05:05', 'HLA-DQA1*03:01-HLA-DQB1*05:06', 'HLA-DQA1*03:01-HLA-DQB1*05:07',
         'HLA-DQA1*03:01-HLA-DQB1*05:08', 'HLA-DQA1*03:01-HLA-DQB1*05:09',
         'HLA-DQA1*03:01-HLA-DQB1*05:10', 'HLA-DQA1*03:01-HLA-DQB1*05:11', 'HLA-DQA1*03:01-HLA-DQB1*05:12', 'HLA-DQA1*03:01-HLA-DQB1*05:13',
         'HLA-DQA1*03:01-HLA-DQB1*05:14', 'HLA-DQA1*03:01-HLA-DQB1*06:01',
         'HLA-DQA1*03:01-HLA-DQB1*06:02', 'HLA-DQA1*03:01-HLA-DQB1*06:03', 'HLA-DQA1*03:01-HLA-DQB1*06:04', 'HLA-DQA1*03:01-HLA-DQB1*06:07',
         'HLA-DQA1*03:01-HLA-DQB1*06:08', 'HLA-DQA1*03:01-HLA-DQB1*06:09',
         'HLA-DQA1*03:01-HLA-DQB1*06:10', 'HLA-DQA1*03:01-HLA-DQB1*06:11', 'HLA-DQA1*03:01-HLA-DQB1*06:12', 'HLA-DQA1*03:01-HLA-DQB1*06:14',
         'HLA-DQA1*03:01-HLA-DQB1*06:15', 'HLA-DQA1*03:01-HLA-DQB1*06:16',
         'HLA-DQA1*03:01-HLA-DQB1*06:17', 'HLA-DQA1*03:01-HLA-DQB1*06:18', 'HLA-DQA1*03:01-HLA-DQB1*06:19', 'HLA-DQA1*03:01-HLA-DQB1*06:21',
         'HLA-DQA1*03:01-HLA-DQB1*06:22', 'HLA-DQA1*03:01-HLA-DQB1*06:23',
         'HLA-DQA1*03:01-HLA-DQB1*06:24', 'HLA-DQA1*03:01-HLA-DQB1*06:25', 'HLA-DQA1*03:01-HLA-DQB1*06:27', 'HLA-DQA1*03:01-HLA-DQB1*06:28',
         'HLA-DQA1*03:01-HLA-DQB1*06:29', 'HLA-DQA1*03:01-HLA-DQB1*06:30',
         'HLA-DQA1*03:01-HLA-DQB1*06:31', 'HLA-DQA1*03:01-HLA-DQB1*06:32', 'HLA-DQA1*03:01-HLA-DQB1*06:33', 'HLA-DQA1*03:01-HLA-DQB1*06:34',
         'HLA-DQA1*03:01-HLA-DQB1*06:35', 'HLA-DQA1*03:01-HLA-DQB1*06:36',
         'HLA-DQA1*03:01-HLA-DQB1*06:37', 'HLA-DQA1*03:01-HLA-DQB1*06:38', 'HLA-DQA1*03:01-HLA-DQB1*06:39', 'HLA-DQA1*03:01-HLA-DQB1*06:40',
         'HLA-DQA1*03:01-HLA-DQB1*06:41', 'HLA-DQA1*03:01-HLA-DQB1*06:42',
         'HLA-DQA1*03:01-HLA-DQB1*06:43', 'HLA-DQA1*03:01-HLA-DQB1*06:44', 'HLA-DQA1*03:02-HLA-DQB1*02:01', 'HLA-DQA1*03:02-HLA-DQB1*02:02',
         'HLA-DQA1*03:02-HLA-DQB1*02:03', 'HLA-DQA1*03:02-HLA-DQB1*02:04',
         'HLA-DQA1*03:02-HLA-DQB1*02:05', 'HLA-DQA1*03:02-HLA-DQB1*02:06', 'HLA-DQA1*03:02-HLA-DQB1*03:01', 'HLA-DQA1*03:02-HLA-DQB1*03:02',
         'HLA-DQA1*03:02-HLA-DQB1*03:03', 'HLA-DQA1*03:02-HLA-DQB1*03:04',
         'HLA-DQA1*03:02-HLA-DQB1*03:05', 'HLA-DQA1*03:02-HLA-DQB1*03:06', 'HLA-DQA1*03:02-HLA-DQB1*03:07', 'HLA-DQA1*03:02-HLA-DQB1*03:08',
         'HLA-DQA1*03:02-HLA-DQB1*03:09', 'HLA-DQA1*03:02-HLA-DQB1*03:10',
         'HLA-DQA1*03:02-HLA-DQB1*03:11', 'HLA-DQA1*03:02-HLA-DQB1*03:12', 'HLA-DQA1*03:02-HLA-DQB1*03:13', 'HLA-DQA1*03:02-HLA-DQB1*03:14',
         'HLA-DQA1*03:02-HLA-DQB1*03:15', 'HLA-DQA1*03:02-HLA-DQB1*03:16',
         'HLA-DQA1*03:02-HLA-DQB1*03:17', 'HLA-DQA1*03:02-HLA-DQB1*03:18', 'HLA-DQA1*03:02-HLA-DQB1*03:19', 'HLA-DQA1*03:02-HLA-DQB1*03:20',
         'HLA-DQA1*03:02-HLA-DQB1*03:21', 'HLA-DQA1*03:02-HLA-DQB1*03:22',
         'HLA-DQA1*03:02-HLA-DQB1*03:23', 'HLA-DQA1*03:02-HLA-DQB1*03:24', 'HLA-DQA1*03:02-HLA-DQB1*03:25', 'HLA-DQA1*03:02-HLA-DQB1*03:26',
         'HLA-DQA1*03:02-HLA-DQB1*03:27', 'HLA-DQA1*03:02-HLA-DQB1*03:28',
         'HLA-DQA1*03:02-HLA-DQB1*03:29', 'HLA-DQA1*03:02-HLA-DQB1*03:30', 'HLA-DQA1*03:02-HLA-DQB1*03:31', 'HLA-DQA1*03:02-HLA-DQB1*03:32',
         'HLA-DQA1*03:02-HLA-DQB1*03:33', 'HLA-DQA1*03:02-HLA-DQB1*03:34',
         'HLA-DQA1*03:02-HLA-DQB1*03:35', 'HLA-DQA1*03:02-HLA-DQB1*03:36', 'HLA-DQA1*03:02-HLA-DQB1*03:37', 'HLA-DQA1*03:02-HLA-DQB1*03:38',
         'HLA-DQA1*03:02-HLA-DQB1*04:01', 'HLA-DQA1*03:02-HLA-DQB1*04:02',
         'HLA-DQA1*03:02-HLA-DQB1*04:03', 'HLA-DQA1*03:02-HLA-DQB1*04:04', 'HLA-DQA1*03:02-HLA-DQB1*04:05', 'HLA-DQA1*03:02-HLA-DQB1*04:06',
         'HLA-DQA1*03:02-HLA-DQB1*04:07', 'HLA-DQA1*03:02-HLA-DQB1*04:08',
         'HLA-DQA1*03:02-HLA-DQB1*05:01', 'HLA-DQA1*03:02-HLA-DQB1*05:02', 'HLA-DQA1*03:02-HLA-DQB1*05:03', 'HLA-DQA1*03:02-HLA-DQB1*05:05',
         'HLA-DQA1*03:02-HLA-DQB1*05:06', 'HLA-DQA1*03:02-HLA-DQB1*05:07',
         'HLA-DQA1*03:02-HLA-DQB1*05:08', 'HLA-DQA1*03:02-HLA-DQB1*05:09', 'HLA-DQA1*03:02-HLA-DQB1*05:10', 'HLA-DQA1*03:02-HLA-DQB1*05:11',
         'HLA-DQA1*03:02-HLA-DQB1*05:12', 'HLA-DQA1*03:02-HLA-DQB1*05:13',
         'HLA-DQA1*03:02-HLA-DQB1*05:14', 'HLA-DQA1*03:02-HLA-DQB1*06:01', 'HLA-DQA1*03:02-HLA-DQB1*06:02', 'HLA-DQA1*03:02-HLA-DQB1*06:03',
         'HLA-DQA1*03:02-HLA-DQB1*06:04', 'HLA-DQA1*03:02-HLA-DQB1*06:07',
         'HLA-DQA1*03:02-HLA-DQB1*06:08', 'HLA-DQA1*03:02-HLA-DQB1*06:09', 'HLA-DQA1*03:02-HLA-DQB1*06:10', 'HLA-DQA1*03:02-HLA-DQB1*06:11',
         'HLA-DQA1*03:02-HLA-DQB1*06:12', 'HLA-DQA1*03:02-HLA-DQB1*06:14',
         'HLA-DQA1*03:02-HLA-DQB1*06:15', 'HLA-DQA1*03:02-HLA-DQB1*06:16', 'HLA-DQA1*03:02-HLA-DQB1*06:17', 'HLA-DQA1*03:02-HLA-DQB1*06:18',
         'HLA-DQA1*03:02-HLA-DQB1*06:19', 'HLA-DQA1*03:02-HLA-DQB1*06:21',
         'HLA-DQA1*03:02-HLA-DQB1*06:22', 'HLA-DQA1*03:02-HLA-DQB1*06:23', 'HLA-DQA1*03:02-HLA-DQB1*06:24', 'HLA-DQA1*03:02-HLA-DQB1*06:25',
         'HLA-DQA1*03:02-HLA-DQB1*06:27', 'HLA-DQA1*03:02-HLA-DQB1*06:28',
         'HLA-DQA1*03:02-HLA-DQB1*06:29', 'HLA-DQA1*03:02-HLA-DQB1*06:30', 'HLA-DQA1*03:02-HLA-DQB1*06:31', 'HLA-DQA1*03:02-HLA-DQB1*06:32',
         'HLA-DQA1*03:02-HLA-DQB1*06:33', 'HLA-DQA1*03:02-HLA-DQB1*06:34',
         'HLA-DQA1*03:02-HLA-DQB1*06:35', 'HLA-DQA1*03:02-HLA-DQB1*06:36', 'HLA-DQA1*03:02-HLA-DQB1*06:37', 'HLA-DQA1*03:02-HLA-DQB1*06:38',
         'HLA-DQA1*03:02-HLA-DQB1*06:39', 'HLA-DQA1*03:02-HLA-DQB1*06:40',
         'HLA-DQA1*03:02-HLA-DQB1*06:41', 'HLA-DQA1*03:02-HLA-DQB1*06:42', 'HLA-DQA1*03:02-HLA-DQB1*06:43', 'HLA-DQA1*03:02-HLA-DQB1*06:44',
         'HLA-DQA1*03:03-HLA-DQB1*02:01', 'HLA-DQA1*03:03-HLA-DQB1*02:02',
         'HLA-DQA1*03:03-HLA-DQB1*02:03', 'HLA-DQA1*03:03-HLA-DQB1*02:04', 'HLA-DQA1*03:03-HLA-DQB1*02:05', 'HLA-DQA1*03:03-HLA-DQB1*02:06',
         'HLA-DQA1*03:03-HLA-DQB1*03:01', 'HLA-DQA1*03:03-HLA-DQB1*03:02',
         'HLA-DQA1*03:03-HLA-DQB1*03:03', 'HLA-DQA1*03:03-HLA-DQB1*03:04', 'HLA-DQA1*03:03-HLA-DQB1*03:05', 'HLA-DQA1*03:03-HLA-DQB1*03:06',
         'HLA-DQA1*03:03-HLA-DQB1*03:07', 'HLA-DQA1*03:03-HLA-DQB1*03:08',
         'HLA-DQA1*03:03-HLA-DQB1*03:09', 'HLA-DQA1*03:03-HLA-DQB1*03:10', 'HLA-DQA1*03:03-HLA-DQB1*03:11', 'HLA-DQA1*03:03-HLA-DQB1*03:12',
         'HLA-DQA1*03:03-HLA-DQB1*03:13', 'HLA-DQA1*03:03-HLA-DQB1*03:14',
         'HLA-DQA1*03:03-HLA-DQB1*03:15', 'HLA-DQA1*03:03-HLA-DQB1*03:16', 'HLA-DQA1*03:03-HLA-DQB1*03:17', 'HLA-DQA1*03:03-HLA-DQB1*03:18',
         'HLA-DQA1*03:03-HLA-DQB1*03:19', 'HLA-DQA1*03:03-HLA-DQB1*03:20',
         'HLA-DQA1*03:03-HLA-DQB1*03:21', 'HLA-DQA1*03:03-HLA-DQB1*03:22', 'HLA-DQA1*03:03-HLA-DQB1*03:23', 'HLA-DQA1*03:03-HLA-DQB1*03:24',
         'HLA-DQA1*03:03-HLA-DQB1*03:25', 'HLA-DQA1*03:03-HLA-DQB1*03:26',
         'HLA-DQA1*03:03-HLA-DQB1*03:27', 'HLA-DQA1*03:03-HLA-DQB1*03:28', 'HLA-DQA1*03:03-HLA-DQB1*03:29', 'HLA-DQA1*03:03-HLA-DQB1*03:30',
         'HLA-DQA1*03:03-HLA-DQB1*03:31', 'HLA-DQA1*03:03-HLA-DQB1*03:32',
         'HLA-DQA1*03:03-HLA-DQB1*03:33', 'HLA-DQA1*03:03-HLA-DQB1*03:34', 'HLA-DQA1*03:03-HLA-DQB1*03:35', 'HLA-DQA1*03:03-HLA-DQB1*03:36',
         'HLA-DQA1*03:03-HLA-DQB1*03:37', 'HLA-DQA1*03:03-HLA-DQB1*03:38',
         'HLA-DQA1*03:03-HLA-DQB1*04:01', 'HLA-DQA1*03:03-HLA-DQB1*04:02', 'HLA-DQA1*03:03-HLA-DQB1*04:03', 'HLA-DQA1*03:03-HLA-DQB1*04:04',
         'HLA-DQA1*03:03-HLA-DQB1*04:05', 'HLA-DQA1*03:03-HLA-DQB1*04:06',
         'HLA-DQA1*03:03-HLA-DQB1*04:07', 'HLA-DQA1*03:03-HLA-DQB1*04:08', 'HLA-DQA1*03:03-HLA-DQB1*05:01', 'HLA-DQA1*03:03-HLA-DQB1*05:02',
         'HLA-DQA1*03:03-HLA-DQB1*05:03', 'HLA-DQA1*03:03-HLA-DQB1*05:05',
         'HLA-DQA1*03:03-HLA-DQB1*05:06', 'HLA-DQA1*03:03-HLA-DQB1*05:07', 'HLA-DQA1*03:03-HLA-DQB1*05:08', 'HLA-DQA1*03:03-HLA-DQB1*05:09',
         'HLA-DQA1*03:03-HLA-DQB1*05:10', 'HLA-DQA1*03:03-HLA-DQB1*05:11',
         'HLA-DQA1*03:03-HLA-DQB1*05:12', 'HLA-DQA1*03:03-HLA-DQB1*05:13', 'HLA-DQA1*03:03-HLA-DQB1*05:14', 'HLA-DQA1*03:03-HLA-DQB1*06:01',
         'HLA-DQA1*03:03-HLA-DQB1*06:02', 'HLA-DQA1*03:03-HLA-DQB1*06:03',
         'HLA-DQA1*03:03-HLA-DQB1*06:04', 'HLA-DQA1*03:03-HLA-DQB1*06:07', 'HLA-DQA1*03:03-HLA-DQB1*06:08', 'HLA-DQA1*03:03-HLA-DQB1*06:09',
         'HLA-DQA1*03:03-HLA-DQB1*06:10', 'HLA-DQA1*03:03-HLA-DQB1*06:11',
         'HLA-DQA1*03:03-HLA-DQB1*06:12', 'HLA-DQA1*03:03-HLA-DQB1*06:14', 'HLA-DQA1*03:03-HLA-DQB1*06:15', 'HLA-DQA1*03:03-HLA-DQB1*06:16',
         'HLA-DQA1*03:03-HLA-DQB1*06:17', 'HLA-DQA1*03:03-HLA-DQB1*06:18',
         'HLA-DQA1*03:03-HLA-DQB1*06:19', 'HLA-DQA1*03:03-HLA-DQB1*06:21', 'HLA-DQA1*03:03-HLA-DQB1*06:22', 'HLA-DQA1*03:03-HLA-DQB1*06:23',
         'HLA-DQA1*03:03-HLA-DQB1*06:24', 'HLA-DQA1*03:03-HLA-DQB1*06:25',
         'HLA-DQA1*03:03-HLA-DQB1*06:27', 'HLA-DQA1*03:03-HLA-DQB1*06:28', 'HLA-DQA1*03:03-HLA-DQB1*06:29', 'HLA-DQA1*03:03-HLA-DQB1*06:30',
         'HLA-DQA1*03:03-HLA-DQB1*06:31', 'HLA-DQA1*03:03-HLA-DQB1*06:32',
         'HLA-DQA1*03:03-HLA-DQB1*06:33', 'HLA-DQA1*03:03-HLA-DQB1*06:34', 'HLA-DQA1*03:03-HLA-DQB1*06:35', 'HLA-DQA1*03:03-HLA-DQB1*06:36',
         'HLA-DQA1*03:03-HLA-DQB1*06:37', 'HLA-DQA1*03:03-HLA-DQB1*06:38',
         'HLA-DQA1*03:03-HLA-DQB1*06:39', 'HLA-DQA1*03:03-HLA-DQB1*06:40', 'HLA-DQA1*03:03-HLA-DQB1*06:41', 'HLA-DQA1*03:03-HLA-DQB1*06:42',
         'HLA-DQA1*03:03-HLA-DQB1*06:43', 'HLA-DQA1*03:03-HLA-DQB1*06:44',
         'HLA-DQA1*04:01-HLA-DQB1*02:01', 'HLA-DQA1*04:01-HLA-DQB1*02:02', 'HLA-DQA1*04:01-HLA-DQB1*02:03', 'HLA-DQA1*04:01-HLA-DQB1*02:04',
         'HLA-DQA1*04:01-HLA-DQB1*02:05', 'HLA-DQA1*04:01-HLA-DQB1*02:06',
         'HLA-DQA1*04:01-HLA-DQB1*03:01', 'HLA-DQA1*04:01-HLA-DQB1*03:02', 'HLA-DQA1*04:01-HLA-DQB1*03:03', 'HLA-DQA1*04:01-HLA-DQB1*03:04',
         'HLA-DQA1*04:01-HLA-DQB1*03:05', 'HLA-DQA1*04:01-HLA-DQB1*03:06',
         'HLA-DQA1*04:01-HLA-DQB1*03:07', 'HLA-DQA1*04:01-HLA-DQB1*03:08', 'HLA-DQA1*04:01-HLA-DQB1*03:09', 'HLA-DQA1*04:01-HLA-DQB1*03:10',
         'HLA-DQA1*04:01-HLA-DQB1*03:11', 'HLA-DQA1*04:01-HLA-DQB1*03:12',
         'HLA-DQA1*04:01-HLA-DQB1*03:13', 'HLA-DQA1*04:01-HLA-DQB1*03:14', 'HLA-DQA1*04:01-HLA-DQB1*03:15', 'HLA-DQA1*04:01-HLA-DQB1*03:16',
         'HLA-DQA1*04:01-HLA-DQB1*03:17', 'HLA-DQA1*04:01-HLA-DQB1*03:18',
         'HLA-DQA1*04:01-HLA-DQB1*03:19', 'HLA-DQA1*04:01-HLA-DQB1*03:20', 'HLA-DQA1*04:01-HLA-DQB1*03:21', 'HLA-DQA1*04:01-HLA-DQB1*03:22',
         'HLA-DQA1*04:01-HLA-DQB1*03:23', 'HLA-DQA1*04:01-HLA-DQB1*03:24',
         'HLA-DQA1*04:01-HLA-DQB1*03:25', 'HLA-DQA1*04:01-HLA-DQB1*03:26', 'HLA-DQA1*04:01-HLA-DQB1*03:27', 'HLA-DQA1*04:01-HLA-DQB1*03:28',
         'HLA-DQA1*04:01-HLA-DQB1*03:29', 'HLA-DQA1*04:01-HLA-DQB1*03:30',
         'HLA-DQA1*04:01-HLA-DQB1*03:31', 'HLA-DQA1*04:01-HLA-DQB1*03:32', 'HLA-DQA1*04:01-HLA-DQB1*03:33', 'HLA-DQA1*04:01-HLA-DQB1*03:34',
         'HLA-DQA1*04:01-HLA-DQB1*03:35', 'HLA-DQA1*04:01-HLA-DQB1*03:36',
         'HLA-DQA1*04:01-HLA-DQB1*03:37', 'HLA-DQA1*04:01-HLA-DQB1*03:38', 'HLA-DQA1*04:01-HLA-DQB1*04:01', 'HLA-DQA1*04:01-HLA-DQB1*04:02',
         'HLA-DQA1*04:01-HLA-DQB1*04:03', 'HLA-DQA1*04:01-HLA-DQB1*04:04',
         'HLA-DQA1*04:01-HLA-DQB1*04:05', 'HLA-DQA1*04:01-HLA-DQB1*04:06', 'HLA-DQA1*04:01-HLA-DQB1*04:07', 'HLA-DQA1*04:01-HLA-DQB1*04:08',
         'HLA-DQA1*04:01-HLA-DQB1*05:01', 'HLA-DQA1*04:01-HLA-DQB1*05:02',
         'HLA-DQA1*04:01-HLA-DQB1*05:03', 'HLA-DQA1*04:01-HLA-DQB1*05:05', 'HLA-DQA1*04:01-HLA-DQB1*05:06', 'HLA-DQA1*04:01-HLA-DQB1*05:07',
         'HLA-DQA1*04:01-HLA-DQB1*05:08', 'HLA-DQA1*04:01-HLA-DQB1*05:09',
         'HLA-DQA1*04:01-HLA-DQB1*05:10', 'HLA-DQA1*04:01-HLA-DQB1*05:11', 'HLA-DQA1*04:01-HLA-DQB1*05:12', 'HLA-DQA1*04:01-HLA-DQB1*05:13',
         'HLA-DQA1*04:01-HLA-DQB1*05:14', 'HLA-DQA1*04:01-HLA-DQB1*06:01',
         'HLA-DQA1*04:01-HLA-DQB1*06:02', 'HLA-DQA1*04:01-HLA-DQB1*06:03', 'HLA-DQA1*04:01-HLA-DQB1*06:04', 'HLA-DQA1*04:01-HLA-DQB1*06:07',
         'HLA-DQA1*04:01-HLA-DQB1*06:08', 'HLA-DQA1*04:01-HLA-DQB1*06:09',
         'HLA-DQA1*04:01-HLA-DQB1*06:10', 'HLA-DQA1*04:01-HLA-DQB1*06:11', 'HLA-DQA1*04:01-HLA-DQB1*06:12', 'HLA-DQA1*04:01-HLA-DQB1*06:14',
         'HLA-DQA1*04:01-HLA-DQB1*06:15', 'HLA-DQA1*04:01-HLA-DQB1*06:16',
         'HLA-DQA1*04:01-HLA-DQB1*06:17', 'HLA-DQA1*04:01-HLA-DQB1*06:18', 'HLA-DQA1*04:01-HLA-DQB1*06:19', 'HLA-DQA1*04:01-HLA-DQB1*06:21',
         'HLA-DQA1*04:01-HLA-DQB1*06:22', 'HLA-DQA1*04:01-HLA-DQB1*06:23',
         'HLA-DQA1*04:01-HLA-DQB1*06:24', 'HLA-DQA1*04:01-HLA-DQB1*06:25', 'HLA-DQA1*04:01-HLA-DQB1*06:27', 'HLA-DQA1*04:01-HLA-DQB1*06:28',
         'HLA-DQA1*04:01-HLA-DQB1*06:29', 'HLA-DQA1*04:01-HLA-DQB1*06:30',
         'HLA-DQA1*04:01-HLA-DQB1*06:31', 'HLA-DQA1*04:01-HLA-DQB1*06:32', 'HLA-DQA1*04:01-HLA-DQB1*06:33', 'HLA-DQA1*04:01-HLA-DQB1*06:34',
         'HLA-DQA1*04:01-HLA-DQB1*06:35', 'HLA-DQA1*04:01-HLA-DQB1*06:36',
         'HLA-DQA1*04:01-HLA-DQB1*06:37', 'HLA-DQA1*04:01-HLA-DQB1*06:38', 'HLA-DQA1*04:01-HLA-DQB1*06:39', 'HLA-DQA1*04:01-HLA-DQB1*06:40',
         'HLA-DQA1*04:01-HLA-DQB1*06:41', 'HLA-DQA1*04:01-HLA-DQB1*06:42',
         'HLA-DQA1*04:01-HLA-DQB1*06:43', 'HLA-DQA1*04:01-HLA-DQB1*06:44', 'HLA-DQA1*04:02-HLA-DQB1*02:01', 'HLA-DQA1*04:02-HLA-DQB1*02:02',
         'HLA-DQA1*04:02-HLA-DQB1*02:03', 'HLA-DQA1*04:02-HLA-DQB1*02:04',
         'HLA-DQA1*04:02-HLA-DQB1*02:05', 'HLA-DQA1*04:02-HLA-DQB1*02:06', 'HLA-DQA1*04:02-HLA-DQB1*03:01', 'HLA-DQA1*04:02-HLA-DQB1*03:02',
         'HLA-DQA1*04:02-HLA-DQB1*03:03', 'HLA-DQA1*04:02-HLA-DQB1*03:04',
         'HLA-DQA1*04:02-HLA-DQB1*03:05', 'HLA-DQA1*04:02-HLA-DQB1*03:06', 'HLA-DQA1*04:02-HLA-DQB1*03:07', 'HLA-DQA1*04:02-HLA-DQB1*03:08',
         'HLA-DQA1*04:02-HLA-DQB1*03:09', 'HLA-DQA1*04:02-HLA-DQB1*03:10',
         'HLA-DQA1*04:02-HLA-DQB1*03:11', 'HLA-DQA1*04:02-HLA-DQB1*03:12', 'HLA-DQA1*04:02-HLA-DQB1*03:13', 'HLA-DQA1*04:02-HLA-DQB1*03:14',
         'HLA-DQA1*04:02-HLA-DQB1*03:15', 'HLA-DQA1*04:02-HLA-DQB1*03:16',
         'HLA-DQA1*04:02-HLA-DQB1*03:17', 'HLA-DQA1*04:02-HLA-DQB1*03:18', 'HLA-DQA1*04:02-HLA-DQB1*03:19', 'HLA-DQA1*04:02-HLA-DQB1*03:20',
         'HLA-DQA1*04:02-HLA-DQB1*03:21', 'HLA-DQA1*04:02-HLA-DQB1*03:22',
         'HLA-DQA1*04:02-HLA-DQB1*03:23', 'HLA-DQA1*04:02-HLA-DQB1*03:24', 'HLA-DQA1*04:02-HLA-DQB1*03:25', 'HLA-DQA1*04:02-HLA-DQB1*03:26',
         'HLA-DQA1*04:02-HLA-DQB1*03:27', 'HLA-DQA1*04:02-HLA-DQB1*03:28',
         'HLA-DQA1*04:02-HLA-DQB1*03:29', 'HLA-DQA1*04:02-HLA-DQB1*03:30', 'HLA-DQA1*04:02-HLA-DQB1*03:31', 'HLA-DQA1*04:02-HLA-DQB1*03:32',
         'HLA-DQA1*04:02-HLA-DQB1*03:33', 'HLA-DQA1*04:02-HLA-DQB1*03:34',
         'HLA-DQA1*04:02-HLA-DQB1*03:35', 'HLA-DQA1*04:02-HLA-DQB1*03:36', 'HLA-DQA1*04:02-HLA-DQB1*03:37', 'HLA-DQA1*04:02-HLA-DQB1*03:38',
         'HLA-DQA1*04:02-HLA-DQB1*04:01', 'HLA-DQA1*04:02-HLA-DQB1*04:02',
         'HLA-DQA1*04:02-HLA-DQB1*04:03', 'HLA-DQA1*04:02-HLA-DQB1*04:04', 'HLA-DQA1*04:02-HLA-DQB1*04:05', 'HLA-DQA1*04:02-HLA-DQB1*04:06',
         'HLA-DQA1*04:02-HLA-DQB1*04:07', 'HLA-DQA1*04:02-HLA-DQB1*04:08',
         'HLA-DQA1*04:02-HLA-DQB1*05:01', 'HLA-DQA1*04:02-HLA-DQB1*05:02', 'HLA-DQA1*04:02-HLA-DQB1*05:03', 'HLA-DQA1*04:02-HLA-DQB1*05:05',
         'HLA-DQA1*04:02-HLA-DQB1*05:06', 'HLA-DQA1*04:02-HLA-DQB1*05:07',
         'HLA-DQA1*04:02-HLA-DQB1*05:08', 'HLA-DQA1*04:02-HLA-DQB1*05:09', 'HLA-DQA1*04:02-HLA-DQB1*05:10', 'HLA-DQA1*04:02-HLA-DQB1*05:11',
         'HLA-DQA1*04:02-HLA-DQB1*05:12', 'HLA-DQA1*04:02-HLA-DQB1*05:13',
         'HLA-DQA1*04:02-HLA-DQB1*05:14', 'HLA-DQA1*04:02-HLA-DQB1*06:01', 'HLA-DQA1*04:02-HLA-DQB1*06:02', 'HLA-DQA1*04:02-HLA-DQB1*06:03',
         'HLA-DQA1*04:02-HLA-DQB1*06:04', 'HLA-DQA1*04:02-HLA-DQB1*06:07',
         'HLA-DQA1*04:02-HLA-DQB1*06:08', 'HLA-DQA1*04:02-HLA-DQB1*06:09', 'HLA-DQA1*04:02-HLA-DQB1*06:10', 'HLA-DQA1*04:02-HLA-DQB1*06:11',
         'HLA-DQA1*04:02-HLA-DQB1*06:12', 'HLA-DQA1*04:02-HLA-DQB1*06:14',
         'HLA-DQA1*04:02-HLA-DQB1*06:15', 'HLA-DQA1*04:02-HLA-DQB1*06:16', 'HLA-DQA1*04:02-HLA-DQB1*06:17', 'HLA-DQA1*04:02-HLA-DQB1*06:18',
         'HLA-DQA1*04:02-HLA-DQB1*06:19', 'HLA-DQA1*04:02-HLA-DQB1*06:21',
         'HLA-DQA1*04:02-HLA-DQB1*06:22', 'HLA-DQA1*04:02-HLA-DQB1*06:23', 'HLA-DQA1*04:02-HLA-DQB1*06:24', 'HLA-DQA1*04:02-HLA-DQB1*06:25',
         'HLA-DQA1*04:02-HLA-DQB1*06:27', 'HLA-DQA1*04:02-HLA-DQB1*06:28',
         'HLA-DQA1*04:02-HLA-DQB1*06:29', 'HLA-DQA1*04:02-HLA-DQB1*06:30', 'HLA-DQA1*04:02-HLA-DQB1*06:31', 'HLA-DQA1*04:02-HLA-DQB1*06:32',
         'HLA-DQA1*04:02-HLA-DQB1*06:33', 'HLA-DQA1*04:02-HLA-DQB1*06:34',
         'HLA-DQA1*04:02-HLA-DQB1*06:35', 'HLA-DQA1*04:02-HLA-DQB1*06:36', 'HLA-DQA1*04:02-HLA-DQB1*06:37', 'HLA-DQA1*04:02-HLA-DQB1*06:38',
         'HLA-DQA1*04:02-HLA-DQB1*06:39', 'HLA-DQA1*04:02-HLA-DQB1*06:40',
         'HLA-DQA1*04:02-HLA-DQB1*06:41', 'HLA-DQA1*04:02-HLA-DQB1*06:42', 'HLA-DQA1*04:02-HLA-DQB1*06:43', 'HLA-DQA1*04:02-HLA-DQB1*06:44',
         'HLA-DQA1*04:04-HLA-DQB1*02:01', 'HLA-DQA1*04:04-HLA-DQB1*02:02',
         'HLA-DQA1*04:04-HLA-DQB1*02:03', 'HLA-DQA1*04:04-HLA-DQB1*02:04', 'HLA-DQA1*04:04-HLA-DQB1*02:05', 'HLA-DQA1*04:04-HLA-DQB1*02:06',
         'HLA-DQA1*04:04-HLA-DQB1*03:01', 'HLA-DQA1*04:04-HLA-DQB1*03:02',
         'HLA-DQA1*04:04-HLA-DQB1*03:03', 'HLA-DQA1*04:04-HLA-DQB1*03:04', 'HLA-DQA1*04:04-HLA-DQB1*03:05', 'HLA-DQA1*04:04-HLA-DQB1*03:06',
         'HLA-DQA1*04:04-HLA-DQB1*03:07', 'HLA-DQA1*04:04-HLA-DQB1*03:08',
         'HLA-DQA1*04:04-HLA-DQB1*03:09', 'HLA-DQA1*04:04-HLA-DQB1*03:10', 'HLA-DQA1*04:04-HLA-DQB1*03:11', 'HLA-DQA1*04:04-HLA-DQB1*03:12',
         'HLA-DQA1*04:04-HLA-DQB1*03:13', 'HLA-DQA1*04:04-HLA-DQB1*03:14',
         'HLA-DQA1*04:04-HLA-DQB1*03:15', 'HLA-DQA1*04:04-HLA-DQB1*03:16', 'HLA-DQA1*04:04-HLA-DQB1*03:17', 'HLA-DQA1*04:04-HLA-DQB1*03:18',
         'HLA-DQA1*04:04-HLA-DQB1*03:19', 'HLA-DQA1*04:04-HLA-DQB1*03:20',
         'HLA-DQA1*04:04-HLA-DQB1*03:21', 'HLA-DQA1*04:04-HLA-DQB1*03:22', 'HLA-DQA1*04:04-HLA-DQB1*03:23', 'HLA-DQA1*04:04-HLA-DQB1*03:24',
         'HLA-DQA1*04:04-HLA-DQB1*03:25', 'HLA-DQA1*04:04-HLA-DQB1*03:26',
         'HLA-DQA1*04:04-HLA-DQB1*03:27', 'HLA-DQA1*04:04-HLA-DQB1*03:28', 'HLA-DQA1*04:04-HLA-DQB1*03:29', 'HLA-DQA1*04:04-HLA-DQB1*03:30',
         'HLA-DQA1*04:04-HLA-DQB1*03:31', 'HLA-DQA1*04:04-HLA-DQB1*03:32',
         'HLA-DQA1*04:04-HLA-DQB1*03:33', 'HLA-DQA1*04:04-HLA-DQB1*03:34', 'HLA-DQA1*04:04-HLA-DQB1*03:35', 'HLA-DQA1*04:04-HLA-DQB1*03:36',
         'HLA-DQA1*04:04-HLA-DQB1*03:37', 'HLA-DQA1*04:04-HLA-DQB1*03:38',
         'HLA-DQA1*04:04-HLA-DQB1*04:01', 'HLA-DQA1*04:04-HLA-DQB1*04:02', 'HLA-DQA1*04:04-HLA-DQB1*04:03', 'HLA-DQA1*04:04-HLA-DQB1*04:04',
         'HLA-DQA1*04:04-HLA-DQB1*04:05', 'HLA-DQA1*04:04-HLA-DQB1*04:06',
         'HLA-DQA1*04:04-HLA-DQB1*04:07', 'HLA-DQA1*04:04-HLA-DQB1*04:08', 'HLA-DQA1*04:04-HLA-DQB1*05:01', 'HLA-DQA1*04:04-HLA-DQB1*05:02',
         'HLA-DQA1*04:04-HLA-DQB1*05:03', 'HLA-DQA1*04:04-HLA-DQB1*05:05',
         'HLA-DQA1*04:04-HLA-DQB1*05:06', 'HLA-DQA1*04:04-HLA-DQB1*05:07', 'HLA-DQA1*04:04-HLA-DQB1*05:08', 'HLA-DQA1*04:04-HLA-DQB1*05:09',
         'HLA-DQA1*04:04-HLA-DQB1*05:10', 'HLA-DQA1*04:04-HLA-DQB1*05:11',
         'HLA-DQA1*04:04-HLA-DQB1*05:12', 'HLA-DQA1*04:04-HLA-DQB1*05:13', 'HLA-DQA1*04:04-HLA-DQB1*05:14', 'HLA-DQA1*04:04-HLA-DQB1*06:01',
         'HLA-DQA1*04:04-HLA-DQB1*06:02', 'HLA-DQA1*04:04-HLA-DQB1*06:03',
         'HLA-DQA1*04:04-HLA-DQB1*06:04', 'HLA-DQA1*04:04-HLA-DQB1*06:07', 'HLA-DQA1*04:04-HLA-DQB1*06:08', 'HLA-DQA1*04:04-HLA-DQB1*06:09',
         'HLA-DQA1*04:04-HLA-DQB1*06:10', 'HLA-DQA1*04:04-HLA-DQB1*06:11',
         'HLA-DQA1*04:04-HLA-DQB1*06:12', 'HLA-DQA1*04:04-HLA-DQB1*06:14', 'HLA-DQA1*04:04-HLA-DQB1*06:15', 'HLA-DQA1*04:04-HLA-DQB1*06:16',
         'HLA-DQA1*04:04-HLA-DQB1*06:17', 'HLA-DQA1*04:04-HLA-DQB1*06:18',
         'HLA-DQA1*04:04-HLA-DQB1*06:19', 'HLA-DQA1*04:04-HLA-DQB1*06:21', 'HLA-DQA1*04:04-HLA-DQB1*06:22', 'HLA-DQA1*04:04-HLA-DQB1*06:23',
         'HLA-DQA1*04:04-HLA-DQB1*06:24', 'HLA-DQA1*04:04-HLA-DQB1*06:25',
         'HLA-DQA1*04:04-HLA-DQB1*06:27', 'HLA-DQA1*04:04-HLA-DQB1*06:28', 'HLA-DQA1*04:04-HLA-DQB1*06:29', 'HLA-DQA1*04:04-HLA-DQB1*06:30',
         'HLA-DQA1*04:04-HLA-DQB1*06:31', 'HLA-DQA1*04:04-HLA-DQB1*06:32',
         'HLA-DQA1*04:04-HLA-DQB1*06:33', 'HLA-DQA1*04:04-HLA-DQB1*06:34', 'HLA-DQA1*04:04-HLA-DQB1*06:35', 'HLA-DQA1*04:04-HLA-DQB1*06:36',
         'HLA-DQA1*04:04-HLA-DQB1*06:37', 'HLA-DQA1*04:04-HLA-DQB1*06:38',
         'HLA-DQA1*04:04-HLA-DQB1*06:39', 'HLA-DQA1*04:04-HLA-DQB1*06:40', 'HLA-DQA1*04:04-HLA-DQB1*06:41', 'HLA-DQA1*04:04-HLA-DQB1*06:42',
         'HLA-DQA1*04:04-HLA-DQB1*06:43', 'HLA-DQA1*04:04-HLA-DQB1*06:44',
         'HLA-DQA1*05:01-HLA-DQB1*02:01', 'HLA-DQA1*05:01-HLA-DQB1*02:02', 'HLA-DQA1*05:01-HLA-DQB1*02:03', 'HLA-DQA1*05:01-HLA-DQB1*02:04',
         'HLA-DQA1*05:01-HLA-DQB1*02:05', 'HLA-DQA1*05:01-HLA-DQB1*02:06',
         'HLA-DQA1*05:01-HLA-DQB1*03:01', 'HLA-DQA1*05:01-HLA-DQB1*03:02', 'HLA-DQA1*05:01-HLA-DQB1*03:03', 'HLA-DQA1*05:01-HLA-DQB1*03:04',
         'HLA-DQA1*05:01-HLA-DQB1*03:05', 'HLA-DQA1*05:01-HLA-DQB1*03:06',
         'HLA-DQA1*05:01-HLA-DQB1*03:07', 'HLA-DQA1*05:01-HLA-DQB1*03:08', 'HLA-DQA1*05:01-HLA-DQB1*03:09', 'HLA-DQA1*05:01-HLA-DQB1*03:10',
         'HLA-DQA1*05:01-HLA-DQB1*03:11', 'HLA-DQA1*05:01-HLA-DQB1*03:12',
         'HLA-DQA1*05:01-HLA-DQB1*03:13', 'HLA-DQA1*05:01-HLA-DQB1*03:14', 'HLA-DQA1*05:01-HLA-DQB1*03:15', 'HLA-DQA1*05:01-HLA-DQB1*03:16',
         'HLA-DQA1*05:01-HLA-DQB1*03:17', 'HLA-DQA1*05:01-HLA-DQB1*03:18',
         'HLA-DQA1*05:01-HLA-DQB1*03:19', 'HLA-DQA1*05:01-HLA-DQB1*03:20', 'HLA-DQA1*05:01-HLA-DQB1*03:21', 'HLA-DQA1*05:01-HLA-DQB1*03:22',
         'HLA-DQA1*05:01-HLA-DQB1*03:23', 'HLA-DQA1*05:01-HLA-DQB1*03:24',
         'HLA-DQA1*05:01-HLA-DQB1*03:25', 'HLA-DQA1*05:01-HLA-DQB1*03:26', 'HLA-DQA1*05:01-HLA-DQB1*03:27', 'HLA-DQA1*05:01-HLA-DQB1*03:28',
         'HLA-DQA1*05:01-HLA-DQB1*03:29', 'HLA-DQA1*05:01-HLA-DQB1*03:30',
         'HLA-DQA1*05:01-HLA-DQB1*03:31', 'HLA-DQA1*05:01-HLA-DQB1*03:32', 'HLA-DQA1*05:01-HLA-DQB1*03:33', 'HLA-DQA1*05:01-HLA-DQB1*03:34',
         'HLA-DQA1*05:01-HLA-DQB1*03:35', 'HLA-DQA1*05:01-HLA-DQB1*03:36',
         'HLA-DQA1*05:01-HLA-DQB1*03:37', 'HLA-DQA1*05:01-HLA-DQB1*03:38', 'HLA-DQA1*05:01-HLA-DQB1*04:01', 'HLA-DQA1*05:01-HLA-DQB1*04:02',
         'HLA-DQA1*05:01-HLA-DQB1*04:03', 'HLA-DQA1*05:01-HLA-DQB1*04:04',
         'HLA-DQA1*05:01-HLA-DQB1*04:05', 'HLA-DQA1*05:01-HLA-DQB1*04:06', 'HLA-DQA1*05:01-HLA-DQB1*04:07', 'HLA-DQA1*05:01-HLA-DQB1*04:08',
         'HLA-DQA1*05:01-HLA-DQB1*05:01', 'HLA-DQA1*05:01-HLA-DQB1*05:02',
         'HLA-DQA1*05:01-HLA-DQB1*05:03', 'HLA-DQA1*05:01-HLA-DQB1*05:05', 'HLA-DQA1*05:01-HLA-DQB1*05:06', 'HLA-DQA1*05:01-HLA-DQB1*05:07',
         'HLA-DQA1*05:01-HLA-DQB1*05:08', 'HLA-DQA1*05:01-HLA-DQB1*05:09',
         'HLA-DQA1*05:01-HLA-DQB1*05:10', 'HLA-DQA1*05:01-HLA-DQB1*05:11', 'HLA-DQA1*05:01-HLA-DQB1*05:12', 'HLA-DQA1*05:01-HLA-DQB1*05:13',
         'HLA-DQA1*05:01-HLA-DQB1*05:14', 'HLA-DQA1*05:01-HLA-DQB1*06:01',
         'HLA-DQA1*05:01-HLA-DQB1*06:02', 'HLA-DQA1*05:01-HLA-DQB1*06:03', 'HLA-DQA1*05:01-HLA-DQB1*06:04', 'HLA-DQA1*05:01-HLA-DQB1*06:07',
         'HLA-DQA1*05:01-HLA-DQB1*06:08', 'HLA-DQA1*05:01-HLA-DQB1*06:09',
         'HLA-DQA1*05:01-HLA-DQB1*06:10', 'HLA-DQA1*05:01-HLA-DQB1*06:11', 'HLA-DQA1*05:01-HLA-DQB1*06:12', 'HLA-DQA1*05:01-HLA-DQB1*06:14',
         'HLA-DQA1*05:01-HLA-DQB1*06:15', 'HLA-DQA1*05:01-HLA-DQB1*06:16',
         'HLA-DQA1*05:01-HLA-DQB1*06:17', 'HLA-DQA1*05:01-HLA-DQB1*06:18', 'HLA-DQA1*05:01-HLA-DQB1*06:19', 'HLA-DQA1*05:01-HLA-DQB1*06:21',
         'HLA-DQA1*05:01-HLA-DQB1*06:22', 'HLA-DQA1*05:01-HLA-DQB1*06:23',
         'HLA-DQA1*05:01-HLA-DQB1*06:24', 'HLA-DQA1*05:01-HLA-DQB1*06:25', 'HLA-DQA1*05:01-HLA-DQB1*06:27', 'HLA-DQA1*05:01-HLA-DQB1*06:28',
         'HLA-DQA1*05:01-HLA-DQB1*06:29', 'HLA-DQA1*05:01-HLA-DQB1*06:30',
         'HLA-DQA1*05:01-HLA-DQB1*06:31', 'HLA-DQA1*05:01-HLA-DQB1*06:32', 'HLA-DQA1*05:01-HLA-DQB1*06:33', 'HLA-DQA1*05:01-HLA-DQB1*06:34',
         'HLA-DQA1*05:01-HLA-DQB1*06:35', 'HLA-DQA1*05:01-HLA-DQB1*06:36',
         'HLA-DQA1*05:01-HLA-DQB1*06:37', 'HLA-DQA1*05:01-HLA-DQB1*06:38', 'HLA-DQA1*05:01-HLA-DQB1*06:39', 'HLA-DQA1*05:01-HLA-DQB1*06:40',
         'HLA-DQA1*05:01-HLA-DQB1*06:41', 'HLA-DQA1*05:01-HLA-DQB1*06:42',
         'HLA-DQA1*05:01-HLA-DQB1*06:43', 'HLA-DQA1*05:01-HLA-DQB1*06:44', 'HLA-DQA1*05:03-HLA-DQB1*02:01', 'HLA-DQA1*05:03-HLA-DQB1*02:02',
         'HLA-DQA1*05:03-HLA-DQB1*02:03', 'HLA-DQA1*05:03-HLA-DQB1*02:04',
         'HLA-DQA1*05:03-HLA-DQB1*02:05', 'HLA-DQA1*05:03-HLA-DQB1*02:06', 'HLA-DQA1*05:03-HLA-DQB1*03:01', 'HLA-DQA1*05:03-HLA-DQB1*03:02',
         'HLA-DQA1*05:03-HLA-DQB1*03:03', 'HLA-DQA1*05:03-HLA-DQB1*03:04',
         'HLA-DQA1*05:03-HLA-DQB1*03:05', 'HLA-DQA1*05:03-HLA-DQB1*03:06', 'HLA-DQA1*05:03-HLA-DQB1*03:07', 'HLA-DQA1*05:03-HLA-DQB1*03:08',
         'HLA-DQA1*05:03-HLA-DQB1*03:09', 'HLA-DQA1*05:03-HLA-DQB1*03:10',
         'HLA-DQA1*05:03-HLA-DQB1*03:11', 'HLA-DQA1*05:03-HLA-DQB1*03:12', 'HLA-DQA1*05:03-HLA-DQB1*03:13', 'HLA-DQA1*05:03-HLA-DQB1*03:14',
         'HLA-DQA1*05:03-HLA-DQB1*03:15', 'HLA-DQA1*05:03-HLA-DQB1*03:16',
         'HLA-DQA1*05:03-HLA-DQB1*03:17', 'HLA-DQA1*05:03-HLA-DQB1*03:18', 'HLA-DQA1*05:03-HLA-DQB1*03:19', 'HLA-DQA1*05:03-HLA-DQB1*03:20',
         'HLA-DQA1*05:03-HLA-DQB1*03:21', 'HLA-DQA1*05:03-HLA-DQB1*03:22',
         'HLA-DQA1*05:03-HLA-DQB1*03:23', 'HLA-DQA1*05:03-HLA-DQB1*03:24', 'HLA-DQA1*05:03-HLA-DQB1*03:25', 'HLA-DQA1*05:03-HLA-DQB1*03:26',
         'HLA-DQA1*05:03-HLA-DQB1*03:27', 'HLA-DQA1*05:03-HLA-DQB1*03:28',
         'HLA-DQA1*05:03-HLA-DQB1*03:29', 'HLA-DQA1*05:03-HLA-DQB1*03:30', 'HLA-DQA1*05:03-HLA-DQB1*03:31', 'HLA-DQA1*05:03-HLA-DQB1*03:32',
         'HLA-DQA1*05:03-HLA-DQB1*03:33', 'HLA-DQA1*05:03-HLA-DQB1*03:34',
         'HLA-DQA1*05:03-HLA-DQB1*03:35', 'HLA-DQA1*05:03-HLA-DQB1*03:36', 'HLA-DQA1*05:03-HLA-DQB1*03:37', 'HLA-DQA1*05:03-HLA-DQB1*03:38',
         'HLA-DQA1*05:03-HLA-DQB1*04:01', 'HLA-DQA1*05:03-HLA-DQB1*04:02',
         'HLA-DQA1*05:03-HLA-DQB1*04:03', 'HLA-DQA1*05:03-HLA-DQB1*04:04', 'HLA-DQA1*05:03-HLA-DQB1*04:05', 'HLA-DQA1*05:03-HLA-DQB1*04:06',
         'HLA-DQA1*05:03-HLA-DQB1*04:07', 'HLA-DQA1*05:03-HLA-DQB1*04:08',
         'HLA-DQA1*05:03-HLA-DQB1*05:01', 'HLA-DQA1*05:03-HLA-DQB1*05:02', 'HLA-DQA1*05:03-HLA-DQB1*05:03', 'HLA-DQA1*05:03-HLA-DQB1*05:05',
         'HLA-DQA1*05:03-HLA-DQB1*05:06', 'HLA-DQA1*05:03-HLA-DQB1*05:07',
         'HLA-DQA1*05:03-HLA-DQB1*05:08', 'HLA-DQA1*05:03-HLA-DQB1*05:09', 'HLA-DQA1*05:03-HLA-DQB1*05:10', 'HLA-DQA1*05:03-HLA-DQB1*05:11',
         'HLA-DQA1*05:03-HLA-DQB1*05:12', 'HLA-DQA1*05:03-HLA-DQB1*05:13',
         'HLA-DQA1*05:03-HLA-DQB1*05:14', 'HLA-DQA1*05:03-HLA-DQB1*06:01', 'HLA-DQA1*05:03-HLA-DQB1*06:02', 'HLA-DQA1*05:03-HLA-DQB1*06:03',
         'HLA-DQA1*05:03-HLA-DQB1*06:04', 'HLA-DQA1*05:03-HLA-DQB1*06:07',
         'HLA-DQA1*05:03-HLA-DQB1*06:08', 'HLA-DQA1*05:03-HLA-DQB1*06:09', 'HLA-DQA1*05:03-HLA-DQB1*06:10', 'HLA-DQA1*05:03-HLA-DQB1*06:11',
         'HLA-DQA1*05:03-HLA-DQB1*06:12', 'HLA-DQA1*05:03-HLA-DQB1*06:14',
         'HLA-DQA1*05:03-HLA-DQB1*06:15', 'HLA-DQA1*05:03-HLA-DQB1*06:16', 'HLA-DQA1*05:03-HLA-DQB1*06:17', 'HLA-DQA1*05:03-HLA-DQB1*06:18',
         'HLA-DQA1*05:03-HLA-DQB1*06:19', 'HLA-DQA1*05:03-HLA-DQB1*06:21',
         'HLA-DQA1*05:03-HLA-DQB1*06:22', 'HLA-DQA1*05:03-HLA-DQB1*06:23', 'HLA-DQA1*05:03-HLA-DQB1*06:24', 'HLA-DQA1*05:03-HLA-DQB1*06:25',
         'HLA-DQA1*05:03-HLA-DQB1*06:27', 'HLA-DQA1*05:03-HLA-DQB1*06:28',
         'HLA-DQA1*05:03-HLA-DQB1*06:29', 'HLA-DQA1*05:03-HLA-DQB1*06:30', 'HLA-DQA1*05:03-HLA-DQB1*06:31', 'HLA-DQA1*05:03-HLA-DQB1*06:32',
         'HLA-DQA1*05:03-HLA-DQB1*06:33', 'HLA-DQA1*05:03-HLA-DQB1*06:34',
         'HLA-DQA1*05:03-HLA-DQB1*06:35', 'HLA-DQA1*05:03-HLA-DQB1*06:36', 'HLA-DQA1*05:03-HLA-DQB1*06:37', 'HLA-DQA1*05:03-HLA-DQB1*06:38',
         'HLA-DQA1*05:03-HLA-DQB1*06:39', 'HLA-DQA1*05:03-HLA-DQB1*06:40',
         'HLA-DQA1*05:03-HLA-DQB1*06:41', 'HLA-DQA1*05:03-HLA-DQB1*06:42', 'HLA-DQA1*05:03-HLA-DQB1*06:43', 'HLA-DQA1*05:03-HLA-DQB1*06:44',
         'HLA-DQA1*05:04-HLA-DQB1*02:01', 'HLA-DQA1*05:04-HLA-DQB1*02:02',
         'HLA-DQA1*05:04-HLA-DQB1*02:03', 'HLA-DQA1*05:04-HLA-DQB1*02:04', 'HLA-DQA1*05:04-HLA-DQB1*02:05', 'HLA-DQA1*05:04-HLA-DQB1*02:06',
         'HLA-DQA1*05:04-HLA-DQB1*03:01', 'HLA-DQA1*05:04-HLA-DQB1*03:02',
         'HLA-DQA1*05:04-HLA-DQB1*03:03', 'HLA-DQA1*05:04-HLA-DQB1*03:04', 'HLA-DQA1*05:04-HLA-DQB1*03:05', 'HLA-DQA1*05:04-HLA-DQB1*03:06',
         'HLA-DQA1*05:04-HLA-DQB1*03:07', 'HLA-DQA1*05:04-HLA-DQB1*03:08',
         'HLA-DQA1*05:04-HLA-DQB1*03:09', 'HLA-DQA1*05:04-HLA-DQB1*03:10', 'HLA-DQA1*05:04-HLA-DQB1*03:11', 'HLA-DQA1*05:04-HLA-DQB1*03:12',
         'HLA-DQA1*05:04-HLA-DQB1*03:13', 'HLA-DQA1*05:04-HLA-DQB1*03:14',
         'HLA-DQA1*05:04-HLA-DQB1*03:15', 'HLA-DQA1*05:04-HLA-DQB1*03:16', 'HLA-DQA1*05:04-HLA-DQB1*03:17', 'HLA-DQA1*05:04-HLA-DQB1*03:18',
         'HLA-DQA1*05:04-HLA-DQB1*03:19', 'HLA-DQA1*05:04-HLA-DQB1*03:20',
         'HLA-DQA1*05:04-HLA-DQB1*03:21', 'HLA-DQA1*05:04-HLA-DQB1*03:22', 'HLA-DQA1*05:04-HLA-DQB1*03:23', 'HLA-DQA1*05:04-HLA-DQB1*03:24',
         'HLA-DQA1*05:04-HLA-DQB1*03:25', 'HLA-DQA1*05:04-HLA-DQB1*03:26',
         'HLA-DQA1*05:04-HLA-DQB1*03:27', 'HLA-DQA1*05:04-HLA-DQB1*03:28', 'HLA-DQA1*05:04-HLA-DQB1*03:29', 'HLA-DQA1*05:04-HLA-DQB1*03:30',
         'HLA-DQA1*05:04-HLA-DQB1*03:31', 'HLA-DQA1*05:04-HLA-DQB1*03:32',
         'HLA-DQA1*05:04-HLA-DQB1*03:33', 'HLA-DQA1*05:04-HLA-DQB1*03:34', 'HLA-DQA1*05:04-HLA-DQB1*03:35', 'HLA-DQA1*05:04-HLA-DQB1*03:36',
         'HLA-DQA1*05:04-HLA-DQB1*03:37', 'HLA-DQA1*05:04-HLA-DQB1*03:38',
         'HLA-DQA1*05:04-HLA-DQB1*04:01', 'HLA-DQA1*05:04-HLA-DQB1*04:02', 'HLA-DQA1*05:04-HLA-DQB1*04:03', 'HLA-DQA1*05:04-HLA-DQB1*04:04',
         'HLA-DQA1*05:04-HLA-DQB1*04:05', 'HLA-DQA1*05:04-HLA-DQB1*04:06',
         'HLA-DQA1*05:04-HLA-DQB1*04:07', 'HLA-DQA1*05:04-HLA-DQB1*04:08', 'HLA-DQA1*05:04-HLA-DQB1*05:01', 'HLA-DQA1*05:04-HLA-DQB1*05:02',
         'HLA-DQA1*05:04-HLA-DQB1*05:03', 'HLA-DQA1*05:04-HLA-DQB1*05:05',
         'HLA-DQA1*05:04-HLA-DQB1*05:06', 'HLA-DQA1*05:04-HLA-DQB1*05:07', 'HLA-DQA1*05:04-HLA-DQB1*05:08', 'HLA-DQA1*05:04-HLA-DQB1*05:09',
         'HLA-DQA1*05:04-HLA-DQB1*05:10', 'HLA-DQA1*05:04-HLA-DQB1*05:11',
         'HLA-DQA1*05:04-HLA-DQB1*05:12', 'HLA-DQA1*05:04-HLA-DQB1*05:13', 'HLA-DQA1*05:04-HLA-DQB1*05:14', 'HLA-DQA1*05:04-HLA-DQB1*06:01',
         'HLA-DQA1*05:04-HLA-DQB1*06:02', 'HLA-DQA1*05:04-HLA-DQB1*06:03',
         'HLA-DQA1*05:04-HLA-DQB1*06:04', 'HLA-DQA1*05:04-HLA-DQB1*06:07', 'HLA-DQA1*05:04-HLA-DQB1*06:08', 'HLA-DQA1*05:04-HLA-DQB1*06:09',
         'HLA-DQA1*05:04-HLA-DQB1*06:10', 'HLA-DQA1*05:04-HLA-DQB1*06:11',
         'HLA-DQA1*05:04-HLA-DQB1*06:12', 'HLA-DQA1*05:04-HLA-DQB1*06:14', 'HLA-DQA1*05:04-HLA-DQB1*06:15', 'HLA-DQA1*05:04-HLA-DQB1*06:16',
         'HLA-DQA1*05:04-HLA-DQB1*06:17', 'HLA-DQA1*05:04-HLA-DQB1*06:18',
         'HLA-DQA1*05:04-HLA-DQB1*06:19', 'HLA-DQA1*05:04-HLA-DQB1*06:21', 'HLA-DQA1*05:04-HLA-DQB1*06:22', 'HLA-DQA1*05:04-HLA-DQB1*06:23',
         'HLA-DQA1*05:04-HLA-DQB1*06:24', 'HLA-DQA1*05:04-HLA-DQB1*06:25',
         'HLA-DQA1*05:04-HLA-DQB1*06:27', 'HLA-DQA1*05:04-HLA-DQB1*06:28', 'HLA-DQA1*05:04-HLA-DQB1*06:29', 'HLA-DQA1*05:04-HLA-DQB1*06:30',
         'HLA-DQA1*05:04-HLA-DQB1*06:31', 'HLA-DQA1*05:04-HLA-DQB1*06:32',
         'HLA-DQA1*05:04-HLA-DQB1*06:33', 'HLA-DQA1*05:04-HLA-DQB1*06:34', 'HLA-DQA1*05:04-HLA-DQB1*06:35', 'HLA-DQA1*05:04-HLA-DQB1*06:36',
         'HLA-DQA1*05:04-HLA-DQB1*06:37', 'HLA-DQA1*05:04-HLA-DQB1*06:38',
         'HLA-DQA1*05:04-HLA-DQB1*06:39', 'HLA-DQA1*05:04-HLA-DQB1*06:40', 'HLA-DQA1*05:04-HLA-DQB1*06:41', 'HLA-DQA1*05:04-HLA-DQB1*06:42',
         'HLA-DQA1*05:04-HLA-DQB1*06:43', 'HLA-DQA1*05:04-HLA-DQB1*06:44',
         'HLA-DQA1*05:05-HLA-DQB1*02:01', 'HLA-DQA1*05:05-HLA-DQB1*02:02', 'HLA-DQA1*05:05-HLA-DQB1*02:03', 'HLA-DQA1*05:05-HLA-DQB1*02:04',
         'HLA-DQA1*05:05-HLA-DQB1*02:05', 'HLA-DQA1*05:05-HLA-DQB1*02:06',
         'HLA-DQA1*05:05-HLA-DQB1*03:01', 'HLA-DQA1*05:05-HLA-DQB1*03:02', 'HLA-DQA1*05:05-HLA-DQB1*03:03', 'HLA-DQA1*05:05-HLA-DQB1*03:04',
         'HLA-DQA1*05:05-HLA-DQB1*03:05', 'HLA-DQA1*05:05-HLA-DQB1*03:06',
         'HLA-DQA1*05:05-HLA-DQB1*03:07', 'HLA-DQA1*05:05-HLA-DQB1*03:08', 'HLA-DQA1*05:05-HLA-DQB1*03:09', 'HLA-DQA1*05:05-HLA-DQB1*03:10',
         'HLA-DQA1*05:05-HLA-DQB1*03:11', 'HLA-DQA1*05:05-HLA-DQB1*03:12',
         'HLA-DQA1*05:05-HLA-DQB1*03:13', 'HLA-DQA1*05:05-HLA-DQB1*03:14', 'HLA-DQA1*05:05-HLA-DQB1*03:15', 'HLA-DQA1*05:05-HLA-DQB1*03:16',
         'HLA-DQA1*05:05-HLA-DQB1*03:17', 'HLA-DQA1*05:05-HLA-DQB1*03:18',
         'HLA-DQA1*05:05-HLA-DQB1*03:19', 'HLA-DQA1*05:05-HLA-DQB1*03:20', 'HLA-DQA1*05:05-HLA-DQB1*03:21', 'HLA-DQA1*05:05-HLA-DQB1*03:22',
         'HLA-DQA1*05:05-HLA-DQB1*03:23', 'HLA-DQA1*05:05-HLA-DQB1*03:24',
         'HLA-DQA1*05:05-HLA-DQB1*03:25', 'HLA-DQA1*05:05-HLA-DQB1*03:26', 'HLA-DQA1*05:05-HLA-DQB1*03:27', 'HLA-DQA1*05:05-HLA-DQB1*03:28',
         'HLA-DQA1*05:05-HLA-DQB1*03:29', 'HLA-DQA1*05:05-HLA-DQB1*03:30',
         'HLA-DQA1*05:05-HLA-DQB1*03:31', 'HLA-DQA1*05:05-HLA-DQB1*03:32', 'HLA-DQA1*05:05-HLA-DQB1*03:33', 'HLA-DQA1*05:05-HLA-DQB1*03:34',
         'HLA-DQA1*05:05-HLA-DQB1*03:35', 'HLA-DQA1*05:05-HLA-DQB1*03:36',
         'HLA-DQA1*05:05-HLA-DQB1*03:37', 'HLA-DQA1*05:05-HLA-DQB1*03:38', 'HLA-DQA1*05:05-HLA-DQB1*04:01', 'HLA-DQA1*05:05-HLA-DQB1*04:02',
         'HLA-DQA1*05:05-HLA-DQB1*04:03', 'HLA-DQA1*05:05-HLA-DQB1*04:04',
         'HLA-DQA1*05:05-HLA-DQB1*04:05', 'HLA-DQA1*05:05-HLA-DQB1*04:06', 'HLA-DQA1*05:05-HLA-DQB1*04:07', 'HLA-DQA1*05:05-HLA-DQB1*04:08',
         'HLA-DQA1*05:05-HLA-DQB1*05:01', 'HLA-DQA1*05:05-HLA-DQB1*05:02',
         'HLA-DQA1*05:05-HLA-DQB1*05:03', 'HLA-DQA1*05:05-HLA-DQB1*05:05', 'HLA-DQA1*05:05-HLA-DQB1*05:06', 'HLA-DQA1*05:05-HLA-DQB1*05:07',
         'HLA-DQA1*05:05-HLA-DQB1*05:08', 'HLA-DQA1*05:05-HLA-DQB1*05:09',
         'HLA-DQA1*05:05-HLA-DQB1*05:10', 'HLA-DQA1*05:05-HLA-DQB1*05:11', 'HLA-DQA1*05:05-HLA-DQB1*05:12', 'HLA-DQA1*05:05-HLA-DQB1*05:13',
         'HLA-DQA1*05:05-HLA-DQB1*05:14', 'HLA-DQA1*05:05-HLA-DQB1*06:01',
         'HLA-DQA1*05:05-HLA-DQB1*06:02', 'HLA-DQA1*05:05-HLA-DQB1*06:03', 'HLA-DQA1*05:05-HLA-DQB1*06:04', 'HLA-DQA1*05:05-HLA-DQB1*06:07',
         'HLA-DQA1*05:05-HLA-DQB1*06:08', 'HLA-DQA1*05:05-HLA-DQB1*06:09',
         'HLA-DQA1*05:05-HLA-DQB1*06:10', 'HLA-DQA1*05:05-HLA-DQB1*06:11', 'HLA-DQA1*05:05-HLA-DQB1*06:12', 'HLA-DQA1*05:05-HLA-DQB1*06:14',
         'HLA-DQA1*05:05-HLA-DQB1*06:15', 'HLA-DQA1*05:05-HLA-DQB1*06:16',
         'HLA-DQA1*05:05-HLA-DQB1*06:17', 'HLA-DQA1*05:05-HLA-DQB1*06:18', 'HLA-DQA1*05:05-HLA-DQB1*06:19', 'HLA-DQA1*05:05-HLA-DQB1*06:21',
         'HLA-DQA1*05:05-HLA-DQB1*06:22', 'HLA-DQA1*05:05-HLA-DQB1*06:23',
         'HLA-DQA1*05:05-HLA-DQB1*06:24', 'HLA-DQA1*05:05-HLA-DQB1*06:25', 'HLA-DQA1*05:05-HLA-DQB1*06:27', 'HLA-DQA1*05:05-HLA-DQB1*06:28',
         'HLA-DQA1*05:05-HLA-DQB1*06:29', 'HLA-DQA1*05:05-HLA-DQB1*06:30',
         'HLA-DQA1*05:05-HLA-DQB1*06:31', 'HLA-DQA1*05:05-HLA-DQB1*06:32', 'HLA-DQA1*05:05-HLA-DQB1*06:33', 'HLA-DQA1*05:05-HLA-DQB1*06:34',
         'HLA-DQA1*05:05-HLA-DQB1*06:35', 'HLA-DQA1*05:05-HLA-DQB1*06:36',
         'HLA-DQA1*05:05-HLA-DQB1*06:37', 'HLA-DQA1*05:05-HLA-DQB1*06:38', 'HLA-DQA1*05:05-HLA-DQB1*06:39', 'HLA-DQA1*05:05-HLA-DQB1*06:40',
         'HLA-DQA1*05:05-HLA-DQB1*06:41', 'HLA-DQA1*05:05-HLA-DQB1*06:42',
         'HLA-DQA1*05:05-HLA-DQB1*06:43', 'HLA-DQA1*05:05-HLA-DQB1*06:44', 'HLA-DQA1*05:06-HLA-DQB1*02:01', 'HLA-DQA1*05:06-HLA-DQB1*02:02',
         'HLA-DQA1*05:06-HLA-DQB1*02:03', 'HLA-DQA1*05:06-HLA-DQB1*02:04',
         'HLA-DQA1*05:06-HLA-DQB1*02:05', 'HLA-DQA1*05:06-HLA-DQB1*02:06', 'HLA-DQA1*05:06-HLA-DQB1*03:01', 'HLA-DQA1*05:06-HLA-DQB1*03:02',
         'HLA-DQA1*05:06-HLA-DQB1*03:03', 'HLA-DQA1*05:06-HLA-DQB1*03:04',
         'HLA-DQA1*05:06-HLA-DQB1*03:05', 'HLA-DQA1*05:06-HLA-DQB1*03:06', 'HLA-DQA1*05:06-HLA-DQB1*03:07', 'HLA-DQA1*05:06-HLA-DQB1*03:08',
         'HLA-DQA1*05:06-HLA-DQB1*03:09', 'HLA-DQA1*05:06-HLA-DQB1*03:10',
         'HLA-DQA1*05:06-HLA-DQB1*03:11', 'HLA-DQA1*05:06-HLA-DQB1*03:12', 'HLA-DQA1*05:06-HLA-DQB1*03:13', 'HLA-DQA1*05:06-HLA-DQB1*03:14',
         'HLA-DQA1*05:06-HLA-DQB1*03:15', 'HLA-DQA1*05:06-HLA-DQB1*03:16',
         'HLA-DQA1*05:06-HLA-DQB1*03:17', 'HLA-DQA1*05:06-HLA-DQB1*03:18', 'HLA-DQA1*05:06-HLA-DQB1*03:19', 'HLA-DQA1*05:06-HLA-DQB1*03:20',
         'HLA-DQA1*05:06-HLA-DQB1*03:21', 'HLA-DQA1*05:06-HLA-DQB1*03:22',
         'HLA-DQA1*05:06-HLA-DQB1*03:23', 'HLA-DQA1*05:06-HLA-DQB1*03:24', 'HLA-DQA1*05:06-HLA-DQB1*03:25', 'HLA-DQA1*05:06-HLA-DQB1*03:26',
         'HLA-DQA1*05:06-HLA-DQB1*03:27', 'HLA-DQA1*05:06-HLA-DQB1*03:28',
         'HLA-DQA1*05:06-HLA-DQB1*03:29', 'HLA-DQA1*05:06-HLA-DQB1*03:30', 'HLA-DQA1*05:06-HLA-DQB1*03:31', 'HLA-DQA1*05:06-HLA-DQB1*03:32',
         'HLA-DQA1*05:06-HLA-DQB1*03:33', 'HLA-DQA1*05:06-HLA-DQB1*03:34',
         'HLA-DQA1*05:06-HLA-DQB1*03:35', 'HLA-DQA1*05:06-HLA-DQB1*03:36', 'HLA-DQA1*05:06-HLA-DQB1*03:37', 'HLA-DQA1*05:06-HLA-DQB1*03:38',
         'HLA-DQA1*05:06-HLA-DQB1*04:01', 'HLA-DQA1*05:06-HLA-DQB1*04:02',
         'HLA-DQA1*05:06-HLA-DQB1*04:03', 'HLA-DQA1*05:06-HLA-DQB1*04:04', 'HLA-DQA1*05:06-HLA-DQB1*04:05', 'HLA-DQA1*05:06-HLA-DQB1*04:06',
         'HLA-DQA1*05:06-HLA-DQB1*04:07', 'HLA-DQA1*05:06-HLA-DQB1*04:08',
         'HLA-DQA1*05:06-HLA-DQB1*05:01', 'HLA-DQA1*05:06-HLA-DQB1*05:02', 'HLA-DQA1*05:06-HLA-DQB1*05:03', 'HLA-DQA1*05:06-HLA-DQB1*05:05',
         'HLA-DQA1*05:06-HLA-DQB1*05:06', 'HLA-DQA1*05:06-HLA-DQB1*05:07',
         'HLA-DQA1*05:06-HLA-DQB1*05:08', 'HLA-DQA1*05:06-HLA-DQB1*05:09', 'HLA-DQA1*05:06-HLA-DQB1*05:10', 'HLA-DQA1*05:06-HLA-DQB1*05:11',
         'HLA-DQA1*05:06-HLA-DQB1*05:12', 'HLA-DQA1*05:06-HLA-DQB1*05:13',
         'HLA-DQA1*05:06-HLA-DQB1*05:14', 'HLA-DQA1*05:06-HLA-DQB1*06:01', 'HLA-DQA1*05:06-HLA-DQB1*06:02', 'HLA-DQA1*05:06-HLA-DQB1*06:03',
         'HLA-DQA1*05:06-HLA-DQB1*06:04', 'HLA-DQA1*05:06-HLA-DQB1*06:07',
         'HLA-DQA1*05:06-HLA-DQB1*06:08', 'HLA-DQA1*05:06-HLA-DQB1*06:09', 'HLA-DQA1*05:06-HLA-DQB1*06:10', 'HLA-DQA1*05:06-HLA-DQB1*06:11',
         'HLA-DQA1*05:06-HLA-DQB1*06:12', 'HLA-DQA1*05:06-HLA-DQB1*06:14',
         'HLA-DQA1*05:06-HLA-DQB1*06:15', 'HLA-DQA1*05:06-HLA-DQB1*06:16', 'HLA-DQA1*05:06-HLA-DQB1*06:17', 'HLA-DQA1*05:06-HLA-DQB1*06:18',
         'HLA-DQA1*05:06-HLA-DQB1*06:19', 'HLA-DQA1*05:06-HLA-DQB1*06:21',
         'HLA-DQA1*05:06-HLA-DQB1*06:22', 'HLA-DQA1*05:06-HLA-DQB1*06:23', 'HLA-DQA1*05:06-HLA-DQB1*06:24', 'HLA-DQA1*05:06-HLA-DQB1*06:25',
         'HLA-DQA1*05:06-HLA-DQB1*06:27', 'HLA-DQA1*05:06-HLA-DQB1*06:28',
         'HLA-DQA1*05:06-HLA-DQB1*06:29', 'HLA-DQA1*05:06-HLA-DQB1*06:30', 'HLA-DQA1*05:06-HLA-DQB1*06:31', 'HLA-DQA1*05:06-HLA-DQB1*06:32',
         'HLA-DQA1*05:06-HLA-DQB1*06:33', 'HLA-DQA1*05:06-HLA-DQB1*06:34',
         'HLA-DQA1*05:06-HLA-DQB1*06:35', 'HLA-DQA1*05:06-HLA-DQB1*06:36', 'HLA-DQA1*05:06-HLA-DQB1*06:37', 'HLA-DQA1*05:06-HLA-DQB1*06:38',
         'HLA-DQA1*05:06-HLA-DQB1*06:39', 'HLA-DQA1*05:06-HLA-DQB1*06:40',
         'HLA-DQA1*05:06-HLA-DQB1*06:41', 'HLA-DQA1*05:06-HLA-DQB1*06:42', 'HLA-DQA1*05:06-HLA-DQB1*06:43', 'HLA-DQA1*05:06-HLA-DQB1*06:44',
         'HLA-DQA1*05:07-HLA-DQB1*02:01', 'HLA-DQA1*05:07-HLA-DQB1*02:02',
         'HLA-DQA1*05:07-HLA-DQB1*02:03', 'HLA-DQA1*05:07-HLA-DQB1*02:04', 'HLA-DQA1*05:07-HLA-DQB1*02:05', 'HLA-DQA1*05:07-HLA-DQB1*02:06',
         'HLA-DQA1*05:07-HLA-DQB1*03:01', 'HLA-DQA1*05:07-HLA-DQB1*03:02',
         'HLA-DQA1*05:07-HLA-DQB1*03:03', 'HLA-DQA1*05:07-HLA-DQB1*03:04', 'HLA-DQA1*05:07-HLA-DQB1*03:05', 'HLA-DQA1*05:07-HLA-DQB1*03:06',
         'HLA-DQA1*05:07-HLA-DQB1*03:07', 'HLA-DQA1*05:07-HLA-DQB1*03:08',
         'HLA-DQA1*05:07-HLA-DQB1*03:09', 'HLA-DQA1*05:07-HLA-DQB1*03:10', 'HLA-DQA1*05:07-HLA-DQB1*03:11', 'HLA-DQA1*05:07-HLA-DQB1*03:12',
         'HLA-DQA1*05:07-HLA-DQB1*03:13', 'HLA-DQA1*05:07-HLA-DQB1*03:14',
         'HLA-DQA1*05:07-HLA-DQB1*03:15', 'HLA-DQA1*05:07-HLA-DQB1*03:16', 'HLA-DQA1*05:07-HLA-DQB1*03:17', 'HLA-DQA1*05:07-HLA-DQB1*03:18',
         'HLA-DQA1*05:07-HLA-DQB1*03:19', 'HLA-DQA1*05:07-HLA-DQB1*03:20',
         'HLA-DQA1*05:07-HLA-DQB1*03:21', 'HLA-DQA1*05:07-HLA-DQB1*03:22', 'HLA-DQA1*05:07-HLA-DQB1*03:23', 'HLA-DQA1*05:07-HLA-DQB1*03:24',
         'HLA-DQA1*05:07-HLA-DQB1*03:25', 'HLA-DQA1*05:07-HLA-DQB1*03:26',
         'HLA-DQA1*05:07-HLA-DQB1*03:27', 'HLA-DQA1*05:07-HLA-DQB1*03:28', 'HLA-DQA1*05:07-HLA-DQB1*03:29', 'HLA-DQA1*05:07-HLA-DQB1*03:30',
         'HLA-DQA1*05:07-HLA-DQB1*03:31', 'HLA-DQA1*05:07-HLA-DQB1*03:32',
         'HLA-DQA1*05:07-HLA-DQB1*03:33', 'HLA-DQA1*05:07-HLA-DQB1*03:34', 'HLA-DQA1*05:07-HLA-DQB1*03:35', 'HLA-DQA1*05:07-HLA-DQB1*03:36',
         'HLA-DQA1*05:07-HLA-DQB1*03:37', 'HLA-DQA1*05:07-HLA-DQB1*03:38',
         'HLA-DQA1*05:07-HLA-DQB1*04:01', 'HLA-DQA1*05:07-HLA-DQB1*04:02', 'HLA-DQA1*05:07-HLA-DQB1*04:03', 'HLA-DQA1*05:07-HLA-DQB1*04:04',
         'HLA-DQA1*05:07-HLA-DQB1*04:05', 'HLA-DQA1*05:07-HLA-DQB1*04:06',
         'HLA-DQA1*05:07-HLA-DQB1*04:07', 'HLA-DQA1*05:07-HLA-DQB1*04:08', 'HLA-DQA1*05:07-HLA-DQB1*05:01', 'HLA-DQA1*05:07-HLA-DQB1*05:02',
         'HLA-DQA1*05:07-HLA-DQB1*05:03', 'HLA-DQA1*05:07-HLA-DQB1*05:05',
         'HLA-DQA1*05:07-HLA-DQB1*05:06', 'HLA-DQA1*05:07-HLA-DQB1*05:07', 'HLA-DQA1*05:07-HLA-DQB1*05:08', 'HLA-DQA1*05:07-HLA-DQB1*05:09',
         'HLA-DQA1*05:07-HLA-DQB1*05:10', 'HLA-DQA1*05:07-HLA-DQB1*05:11',
         'HLA-DQA1*05:07-HLA-DQB1*05:12', 'HLA-DQA1*05:07-HLA-DQB1*05:13', 'HLA-DQA1*05:07-HLA-DQB1*05:14', 'HLA-DQA1*05:07-HLA-DQB1*06:01',
         'HLA-DQA1*05:07-HLA-DQB1*06:02', 'HLA-DQA1*05:07-HLA-DQB1*06:03',
         'HLA-DQA1*05:07-HLA-DQB1*06:04', 'HLA-DQA1*05:07-HLA-DQB1*06:07', 'HLA-DQA1*05:07-HLA-DQB1*06:08', 'HLA-DQA1*05:07-HLA-DQB1*06:09',
         'HLA-DQA1*05:07-HLA-DQB1*06:10', 'HLA-DQA1*05:07-HLA-DQB1*06:11',
         'HLA-DQA1*05:07-HLA-DQB1*06:12', 'HLA-DQA1*05:07-HLA-DQB1*06:14', 'HLA-DQA1*05:07-HLA-DQB1*06:15', 'HLA-DQA1*05:07-HLA-DQB1*06:16',
         'HLA-DQA1*05:07-HLA-DQB1*06:17', 'HLA-DQA1*05:07-HLA-DQB1*06:18',
         'HLA-DQA1*05:07-HLA-DQB1*06:19', 'HLA-DQA1*05:07-HLA-DQB1*06:21', 'HLA-DQA1*05:07-HLA-DQB1*06:22', 'HLA-DQA1*05:07-HLA-DQB1*06:23',
         'HLA-DQA1*05:07-HLA-DQB1*06:24', 'HLA-DQA1*05:07-HLA-DQB1*06:25',
         'HLA-DQA1*05:07-HLA-DQB1*06:27', 'HLA-DQA1*05:07-HLA-DQB1*06:28', 'HLA-DQA1*05:07-HLA-DQB1*06:29', 'HLA-DQA1*05:07-HLA-DQB1*06:30',
         'HLA-DQA1*05:07-HLA-DQB1*06:31', 'HLA-DQA1*05:07-HLA-DQB1*06:32',
         'HLA-DQA1*05:07-HLA-DQB1*06:33', 'HLA-DQA1*05:07-HLA-DQB1*06:34', 'HLA-DQA1*05:07-HLA-DQB1*06:35', 'HLA-DQA1*05:07-HLA-DQB1*06:36',
         'HLA-DQA1*05:07-HLA-DQB1*06:37', 'HLA-DQA1*05:07-HLA-DQB1*06:38',
         'HLA-DQA1*05:07-HLA-DQB1*06:39', 'HLA-DQA1*05:07-HLA-DQB1*06:40', 'HLA-DQA1*05:07-HLA-DQB1*06:41', 'HLA-DQA1*05:07-HLA-DQB1*06:42',
         'HLA-DQA1*05:07-HLA-DQB1*06:43', 'HLA-DQA1*05:07-HLA-DQB1*06:44',
         'HLA-DQA1*05:08-HLA-DQB1*02:01', 'HLA-DQA1*05:08-HLA-DQB1*02:02', 'HLA-DQA1*05:08-HLA-DQB1*02:03', 'HLA-DQA1*05:08-HLA-DQB1*02:04',
         'HLA-DQA1*05:08-HLA-DQB1*02:05', 'HLA-DQA1*05:08-HLA-DQB1*02:06',
         'HLA-DQA1*05:08-HLA-DQB1*03:01', 'HLA-DQA1*05:08-HLA-DQB1*03:02', 'HLA-DQA1*05:08-HLA-DQB1*03:03', 'HLA-DQA1*05:08-HLA-DQB1*03:04',
         'HLA-DQA1*05:08-HLA-DQB1*03:05', 'HLA-DQA1*05:08-HLA-DQB1*03:06',
         'HLA-DQA1*05:08-HLA-DQB1*03:07', 'HLA-DQA1*05:08-HLA-DQB1*03:08', 'HLA-DQA1*05:08-HLA-DQB1*03:09', 'HLA-DQA1*05:08-HLA-DQB1*03:10',
         'HLA-DQA1*05:08-HLA-DQB1*03:11', 'HLA-DQA1*05:08-HLA-DQB1*03:12',
         'HLA-DQA1*05:08-HLA-DQB1*03:13', 'HLA-DQA1*05:08-HLA-DQB1*03:14', 'HLA-DQA1*05:08-HLA-DQB1*03:15', 'HLA-DQA1*05:08-HLA-DQB1*03:16',
         'HLA-DQA1*05:08-HLA-DQB1*03:17', 'HLA-DQA1*05:08-HLA-DQB1*03:18',
         'HLA-DQA1*05:08-HLA-DQB1*03:19', 'HLA-DQA1*05:08-HLA-DQB1*03:20', 'HLA-DQA1*05:08-HLA-DQB1*03:21', 'HLA-DQA1*05:08-HLA-DQB1*03:22',
         'HLA-DQA1*05:08-HLA-DQB1*03:23', 'HLA-DQA1*05:08-HLA-DQB1*03:24',
         'HLA-DQA1*05:08-HLA-DQB1*03:25', 'HLA-DQA1*05:08-HLA-DQB1*03:26', 'HLA-DQA1*05:08-HLA-DQB1*03:27', 'HLA-DQA1*05:08-HLA-DQB1*03:28',
         'HLA-DQA1*05:08-HLA-DQB1*03:29', 'HLA-DQA1*05:08-HLA-DQB1*03:30',
         'HLA-DQA1*05:08-HLA-DQB1*03:31', 'HLA-DQA1*05:08-HLA-DQB1*03:32', 'HLA-DQA1*05:08-HLA-DQB1*03:33', 'HLA-DQA1*05:08-HLA-DQB1*03:34',
         'HLA-DQA1*05:08-HLA-DQB1*03:35', 'HLA-DQA1*05:08-HLA-DQB1*03:36',
         'HLA-DQA1*05:08-HLA-DQB1*03:37', 'HLA-DQA1*05:08-HLA-DQB1*03:38', 'HLA-DQA1*05:08-HLA-DQB1*04:01', 'HLA-DQA1*05:08-HLA-DQB1*04:02',
         'HLA-DQA1*05:08-HLA-DQB1*04:03', 'HLA-DQA1*05:08-HLA-DQB1*04:04',
         'HLA-DQA1*05:08-HLA-DQB1*04:05', 'HLA-DQA1*05:08-HLA-DQB1*04:06', 'HLA-DQA1*05:08-HLA-DQB1*04:07', 'HLA-DQA1*05:08-HLA-DQB1*04:08',
         'HLA-DQA1*05:08-HLA-DQB1*05:01', 'HLA-DQA1*05:08-HLA-DQB1*05:02',
         'HLA-DQA1*05:08-HLA-DQB1*05:03', 'HLA-DQA1*05:08-HLA-DQB1*05:05', 'HLA-DQA1*05:08-HLA-DQB1*05:06', 'HLA-DQA1*05:08-HLA-DQB1*05:07',
         'HLA-DQA1*05:08-HLA-DQB1*05:08', 'HLA-DQA1*05:08-HLA-DQB1*05:09',
         'HLA-DQA1*05:08-HLA-DQB1*05:10', 'HLA-DQA1*05:08-HLA-DQB1*05:11', 'HLA-DQA1*05:08-HLA-DQB1*05:12', 'HLA-DQA1*05:08-HLA-DQB1*05:13',
         'HLA-DQA1*05:08-HLA-DQB1*05:14', 'HLA-DQA1*05:08-HLA-DQB1*06:01',
         'HLA-DQA1*05:08-HLA-DQB1*06:02', 'HLA-DQA1*05:08-HLA-DQB1*06:03', 'HLA-DQA1*05:08-HLA-DQB1*06:04', 'HLA-DQA1*05:08-HLA-DQB1*06:07',
         'HLA-DQA1*05:08-HLA-DQB1*06:08', 'HLA-DQA1*05:08-HLA-DQB1*06:09',
         'HLA-DQA1*05:08-HLA-DQB1*06:10', 'HLA-DQA1*05:08-HLA-DQB1*06:11', 'HLA-DQA1*05:08-HLA-DQB1*06:12', 'HLA-DQA1*05:08-HLA-DQB1*06:14',
         'HLA-DQA1*05:08-HLA-DQB1*06:15', 'HLA-DQA1*05:08-HLA-DQB1*06:16',
         'HLA-DQA1*05:08-HLA-DQB1*06:17', 'HLA-DQA1*05:08-HLA-DQB1*06:18', 'HLA-DQA1*05:08-HLA-DQB1*06:19', 'HLA-DQA1*05:08-HLA-DQB1*06:21',
         'HLA-DQA1*05:08-HLA-DQB1*06:22', 'HLA-DQA1*05:08-HLA-DQB1*06:23',
         'HLA-DQA1*05:08-HLA-DQB1*06:24', 'HLA-DQA1*05:08-HLA-DQB1*06:25', 'HLA-DQA1*05:08-HLA-DQB1*06:27', 'HLA-DQA1*05:08-HLA-DQB1*06:28',
         'HLA-DQA1*05:08-HLA-DQB1*06:29', 'HLA-DQA1*05:08-HLA-DQB1*06:30',
         'HLA-DQA1*05:08-HLA-DQB1*06:31', 'HLA-DQA1*05:08-HLA-DQB1*06:32', 'HLA-DQA1*05:08-HLA-DQB1*06:33', 'HLA-DQA1*05:08-HLA-DQB1*06:34',
         'HLA-DQA1*05:08-HLA-DQB1*06:35', 'HLA-DQA1*05:08-HLA-DQB1*06:36',
         'HLA-DQA1*05:08-HLA-DQB1*06:37', 'HLA-DQA1*05:08-HLA-DQB1*06:38', 'HLA-DQA1*05:08-HLA-DQB1*06:39', 'HLA-DQA1*05:08-HLA-DQB1*06:40',
         'HLA-DQA1*05:08-HLA-DQB1*06:41', 'HLA-DQA1*05:08-HLA-DQB1*06:42',
         'HLA-DQA1*05:08-HLA-DQB1*06:43', 'HLA-DQA1*05:08-HLA-DQB1*06:44', 'HLA-DQA1*05:09-HLA-DQB1*02:01', 'HLA-DQA1*05:09-HLA-DQB1*02:02',
         'HLA-DQA1*05:09-HLA-DQB1*02:03', 'HLA-DQA1*05:09-HLA-DQB1*02:04',
         'HLA-DQA1*05:09-HLA-DQB1*02:05', 'HLA-DQA1*05:09-HLA-DQB1*02:06', 'HLA-DQA1*05:09-HLA-DQB1*03:01', 'HLA-DQA1*05:09-HLA-DQB1*03:02',
         'HLA-DQA1*05:09-HLA-DQB1*03:03', 'HLA-DQA1*05:09-HLA-DQB1*03:04',
         'HLA-DQA1*05:09-HLA-DQB1*03:05', 'HLA-DQA1*05:09-HLA-DQB1*03:06', 'HLA-DQA1*05:09-HLA-DQB1*03:07', 'HLA-DQA1*05:09-HLA-DQB1*03:08',
         'HLA-DQA1*05:09-HLA-DQB1*03:09', 'HLA-DQA1*05:09-HLA-DQB1*03:10',
         'HLA-DQA1*05:09-HLA-DQB1*03:11', 'HLA-DQA1*05:09-HLA-DQB1*03:12', 'HLA-DQA1*05:09-HLA-DQB1*03:13', 'HLA-DQA1*05:09-HLA-DQB1*03:14',
         'HLA-DQA1*05:09-HLA-DQB1*03:15', 'HLA-DQA1*05:09-HLA-DQB1*03:16',
         'HLA-DQA1*05:09-HLA-DQB1*03:17', 'HLA-DQA1*05:09-HLA-DQB1*03:18', 'HLA-DQA1*05:09-HLA-DQB1*03:19', 'HLA-DQA1*05:09-HLA-DQB1*03:20',
         'HLA-DQA1*05:09-HLA-DQB1*03:21', 'HLA-DQA1*05:09-HLA-DQB1*03:22',
         'HLA-DQA1*05:09-HLA-DQB1*03:23', 'HLA-DQA1*05:09-HLA-DQB1*03:24', 'HLA-DQA1*05:09-HLA-DQB1*03:25', 'HLA-DQA1*05:09-HLA-DQB1*03:26',
         'HLA-DQA1*05:09-HLA-DQB1*03:27', 'HLA-DQA1*05:09-HLA-DQB1*03:28',
         'HLA-DQA1*05:09-HLA-DQB1*03:29', 'HLA-DQA1*05:09-HLA-DQB1*03:30', 'HLA-DQA1*05:09-HLA-DQB1*03:31', 'HLA-DQA1*05:09-HLA-DQB1*03:32',
         'HLA-DQA1*05:09-HLA-DQB1*03:33', 'HLA-DQA1*05:09-HLA-DQB1*03:34',
         'HLA-DQA1*05:09-HLA-DQB1*03:35', 'HLA-DQA1*05:09-HLA-DQB1*03:36', 'HLA-DQA1*05:09-HLA-DQB1*03:37', 'HLA-DQA1*05:09-HLA-DQB1*03:38',
         'HLA-DQA1*05:09-HLA-DQB1*04:01', 'HLA-DQA1*05:09-HLA-DQB1*04:02',
         'HLA-DQA1*05:09-HLA-DQB1*04:03', 'HLA-DQA1*05:09-HLA-DQB1*04:04', 'HLA-DQA1*05:09-HLA-DQB1*04:05', 'HLA-DQA1*05:09-HLA-DQB1*04:06',
         'HLA-DQA1*05:09-HLA-DQB1*04:07', 'HLA-DQA1*05:09-HLA-DQB1*04:08',
         'HLA-DQA1*05:09-HLA-DQB1*05:01', 'HLA-DQA1*05:09-HLA-DQB1*05:02', 'HLA-DQA1*05:09-HLA-DQB1*05:03', 'HLA-DQA1*05:09-HLA-DQB1*05:05',
         'HLA-DQA1*05:09-HLA-DQB1*05:06', 'HLA-DQA1*05:09-HLA-DQB1*05:07',
         'HLA-DQA1*05:09-HLA-DQB1*05:08', 'HLA-DQA1*05:09-HLA-DQB1*05:09', 'HLA-DQA1*05:09-HLA-DQB1*05:10', 'HLA-DQA1*05:09-HLA-DQB1*05:11',
         'HLA-DQA1*05:09-HLA-DQB1*05:12', 'HLA-DQA1*05:09-HLA-DQB1*05:13',
         'HLA-DQA1*05:09-HLA-DQB1*05:14', 'HLA-DQA1*05:09-HLA-DQB1*06:01', 'HLA-DQA1*05:09-HLA-DQB1*06:02', 'HLA-DQA1*05:09-HLA-DQB1*06:03',
         'HLA-DQA1*05:09-HLA-DQB1*06:04', 'HLA-DQA1*05:09-HLA-DQB1*06:07',
         'HLA-DQA1*05:09-HLA-DQB1*06:08', 'HLA-DQA1*05:09-HLA-DQB1*06:09', 'HLA-DQA1*05:09-HLA-DQB1*06:10', 'HLA-DQA1*05:09-HLA-DQB1*06:11',
         'HLA-DQA1*05:09-HLA-DQB1*06:12', 'HLA-DQA1*05:09-HLA-DQB1*06:14',
         'HLA-DQA1*05:09-HLA-DQB1*06:15', 'HLA-DQA1*05:09-HLA-DQB1*06:16', 'HLA-DQA1*05:09-HLA-DQB1*06:17', 'HLA-DQA1*05:09-HLA-DQB1*06:18',
         'HLA-DQA1*05:09-HLA-DQB1*06:19', 'HLA-DQA1*05:09-HLA-DQB1*06:21',
         'HLA-DQA1*05:09-HLA-DQB1*06:22', 'HLA-DQA1*05:09-HLA-DQB1*06:23', 'HLA-DQA1*05:09-HLA-DQB1*06:24', 'HLA-DQA1*05:09-HLA-DQB1*06:25',
         'HLA-DQA1*05:09-HLA-DQB1*06:27', 'HLA-DQA1*05:09-HLA-DQB1*06:28',
         'HLA-DQA1*05:09-HLA-DQB1*06:29', 'HLA-DQA1*05:09-HLA-DQB1*06:30', 'HLA-DQA1*05:09-HLA-DQB1*06:31', 'HLA-DQA1*05:09-HLA-DQB1*06:32',
         'HLA-DQA1*05:09-HLA-DQB1*06:33', 'HLA-DQA1*05:09-HLA-DQB1*06:34',
         'HLA-DQA1*05:09-HLA-DQB1*06:35', 'HLA-DQA1*05:09-HLA-DQB1*06:36', 'HLA-DQA1*05:09-HLA-DQB1*06:37', 'HLA-DQA1*05:09-HLA-DQB1*06:38',
         'HLA-DQA1*05:09-HLA-DQB1*06:39', 'HLA-DQA1*05:09-HLA-DQB1*06:40',
         'HLA-DQA1*05:09-HLA-DQB1*06:41', 'HLA-DQA1*05:09-HLA-DQB1*06:42', 'HLA-DQA1*05:09-HLA-DQB1*06:43', 'HLA-DQA1*05:09-HLA-DQB1*06:44',
         'HLA-DQA1*05:10-HLA-DQB1*02:01', 'HLA-DQA1*05:10-HLA-DQB1*02:02',
         'HLA-DQA1*05:10-HLA-DQB1*02:03', 'HLA-DQA1*05:10-HLA-DQB1*02:04', 'HLA-DQA1*05:10-HLA-DQB1*02:05', 'HLA-DQA1*05:10-HLA-DQB1*02:06',
         'HLA-DQA1*05:10-HLA-DQB1*03:01', 'HLA-DQA1*05:10-HLA-DQB1*03:02',
         'HLA-DQA1*05:10-HLA-DQB1*03:03', 'HLA-DQA1*05:10-HLA-DQB1*03:04', 'HLA-DQA1*05:10-HLA-DQB1*03:05', 'HLA-DQA1*05:10-HLA-DQB1*03:06',
         'HLA-DQA1*05:10-HLA-DQB1*03:07', 'HLA-DQA1*05:10-HLA-DQB1*03:08',
         'HLA-DQA1*05:10-HLA-DQB1*03:09', 'HLA-DQA1*05:10-HLA-DQB1*03:10', 'HLA-DQA1*05:10-HLA-DQB1*03:11', 'HLA-DQA1*05:10-HLA-DQB1*03:12',
         'HLA-DQA1*05:10-HLA-DQB1*03:13', 'HLA-DQA1*05:10-HLA-DQB1*03:14',
         'HLA-DQA1*05:10-HLA-DQB1*03:15', 'HLA-DQA1*05:10-HLA-DQB1*03:16', 'HLA-DQA1*05:10-HLA-DQB1*03:17', 'HLA-DQA1*05:10-HLA-DQB1*03:18',
         'HLA-DQA1*05:10-HLA-DQB1*03:19', 'HLA-DQA1*05:10-HLA-DQB1*03:20',
         'HLA-DQA1*05:10-HLA-DQB1*03:21', 'HLA-DQA1*05:10-HLA-DQB1*03:22', 'HLA-DQA1*05:10-HLA-DQB1*03:23', 'HLA-DQA1*05:10-HLA-DQB1*03:24',
         'HLA-DQA1*05:10-HLA-DQB1*03:25', 'HLA-DQA1*05:10-HLA-DQB1*03:26',
         'HLA-DQA1*05:10-HLA-DQB1*03:27', 'HLA-DQA1*05:10-HLA-DQB1*03:28', 'HLA-DQA1*05:10-HLA-DQB1*03:29', 'HLA-DQA1*05:10-HLA-DQB1*03:30',
         'HLA-DQA1*05:10-HLA-DQB1*03:31', 'HLA-DQA1*05:10-HLA-DQB1*03:32',
         'HLA-DQA1*05:10-HLA-DQB1*03:33', 'HLA-DQA1*05:10-HLA-DQB1*03:34', 'HLA-DQA1*05:10-HLA-DQB1*03:35', 'HLA-DQA1*05:10-HLA-DQB1*03:36',
         'HLA-DQA1*05:10-HLA-DQB1*03:37', 'HLA-DQA1*05:10-HLA-DQB1*03:38',
         'HLA-DQA1*05:10-HLA-DQB1*04:01', 'HLA-DQA1*05:10-HLA-DQB1*04:02', 'HLA-DQA1*05:10-HLA-DQB1*04:03', 'HLA-DQA1*05:10-HLA-DQB1*04:04',
         'HLA-DQA1*05:10-HLA-DQB1*04:05', 'HLA-DQA1*05:10-HLA-DQB1*04:06',
         'HLA-DQA1*05:10-HLA-DQB1*04:07', 'HLA-DQA1*05:10-HLA-DQB1*04:08', 'HLA-DQA1*05:10-HLA-DQB1*05:01', 'HLA-DQA1*05:10-HLA-DQB1*05:02',
         'HLA-DQA1*05:10-HLA-DQB1*05:03', 'HLA-DQA1*05:10-HLA-DQB1*05:05',
         'HLA-DQA1*05:10-HLA-DQB1*05:06', 'HLA-DQA1*05:10-HLA-DQB1*05:07', 'HLA-DQA1*05:10-HLA-DQB1*05:08', 'HLA-DQA1*05:10-HLA-DQB1*05:09',
         'HLA-DQA1*05:10-HLA-DQB1*05:10', 'HLA-DQA1*05:10-HLA-DQB1*05:11',
         'HLA-DQA1*05:10-HLA-DQB1*05:12', 'HLA-DQA1*05:10-HLA-DQB1*05:13', 'HLA-DQA1*05:10-HLA-DQB1*05:14', 'HLA-DQA1*05:10-HLA-DQB1*06:01',
         'HLA-DQA1*05:10-HLA-DQB1*06:02', 'HLA-DQA1*05:10-HLA-DQB1*06:03',
         'HLA-DQA1*05:10-HLA-DQB1*06:04', 'HLA-DQA1*05:10-HLA-DQB1*06:07', 'HLA-DQA1*05:10-HLA-DQB1*06:08', 'HLA-DQA1*05:10-HLA-DQB1*06:09',
         'HLA-DQA1*05:10-HLA-DQB1*06:10', 'HLA-DQA1*05:10-HLA-DQB1*06:11',
         'HLA-DQA1*05:10-HLA-DQB1*06:12', 'HLA-DQA1*05:10-HLA-DQB1*06:14', 'HLA-DQA1*05:10-HLA-DQB1*06:15', 'HLA-DQA1*05:10-HLA-DQB1*06:16',
         'HLA-DQA1*05:10-HLA-DQB1*06:17', 'HLA-DQA1*05:10-HLA-DQB1*06:18',
         'HLA-DQA1*05:10-HLA-DQB1*06:19', 'HLA-DQA1*05:10-HLA-DQB1*06:21', 'HLA-DQA1*05:10-HLA-DQB1*06:22', 'HLA-DQA1*05:10-HLA-DQB1*06:23',
         'HLA-DQA1*05:10-HLA-DQB1*06:24', 'HLA-DQA1*05:10-HLA-DQB1*06:25',
         'HLA-DQA1*05:10-HLA-DQB1*06:27', 'HLA-DQA1*05:10-HLA-DQB1*06:28', 'HLA-DQA1*05:10-HLA-DQB1*06:29', 'HLA-DQA1*05:10-HLA-DQB1*06:30',
         'HLA-DQA1*05:10-HLA-DQB1*06:31', 'HLA-DQA1*05:10-HLA-DQB1*06:32',
         'HLA-DQA1*05:10-HLA-DQB1*06:33', 'HLA-DQA1*05:10-HLA-DQB1*06:34', 'HLA-DQA1*05:10-HLA-DQB1*06:35', 'HLA-DQA1*05:10-HLA-DQB1*06:36',
         'HLA-DQA1*05:10-HLA-DQB1*06:37', 'HLA-DQA1*05:10-HLA-DQB1*06:38',
         'HLA-DQA1*05:10-HLA-DQB1*06:39', 'HLA-DQA1*05:10-HLA-DQB1*06:40', 'HLA-DQA1*05:10-HLA-DQB1*06:41', 'HLA-DQA1*05:10-HLA-DQB1*06:42',
         'HLA-DQA1*05:10-HLA-DQB1*06:43', 'HLA-DQA1*05:10-HLA-DQB1*06:44',
         'HLA-DQA1*05:11-HLA-DQB1*02:01', 'HLA-DQA1*05:11-HLA-DQB1*02:02', 'HLA-DQA1*05:11-HLA-DQB1*02:03', 'HLA-DQA1*05:11-HLA-DQB1*02:04',
         'HLA-DQA1*05:11-HLA-DQB1*02:05', 'HLA-DQA1*05:11-HLA-DQB1*02:06',
         'HLA-DQA1*05:11-HLA-DQB1*03:01', 'HLA-DQA1*05:11-HLA-DQB1*03:02', 'HLA-DQA1*05:11-HLA-DQB1*03:03', 'HLA-DQA1*05:11-HLA-DQB1*03:04',
         'HLA-DQA1*05:11-HLA-DQB1*03:05', 'HLA-DQA1*05:11-HLA-DQB1*03:06',
         'HLA-DQA1*05:11-HLA-DQB1*03:07', 'HLA-DQA1*05:11-HLA-DQB1*03:08', 'HLA-DQA1*05:11-HLA-DQB1*03:09', 'HLA-DQA1*05:11-HLA-DQB1*03:10',
         'HLA-DQA1*05:11-HLA-DQB1*03:11', 'HLA-DQA1*05:11-HLA-DQB1*03:12',
         'HLA-DQA1*05:11-HLA-DQB1*03:13', 'HLA-DQA1*05:11-HLA-DQB1*03:14', 'HLA-DQA1*05:11-HLA-DQB1*03:15', 'HLA-DQA1*05:11-HLA-DQB1*03:16',
         'HLA-DQA1*05:11-HLA-DQB1*03:17', 'HLA-DQA1*05:11-HLA-DQB1*03:18',
         'HLA-DQA1*05:11-HLA-DQB1*03:19', 'HLA-DQA1*05:11-HLA-DQB1*03:20', 'HLA-DQA1*05:11-HLA-DQB1*03:21', 'HLA-DQA1*05:11-HLA-DQB1*03:22',
         'HLA-DQA1*05:11-HLA-DQB1*03:23', 'HLA-DQA1*05:11-HLA-DQB1*03:24',
         'HLA-DQA1*05:11-HLA-DQB1*03:25', 'HLA-DQA1*05:11-HLA-DQB1*03:26', 'HLA-DQA1*05:11-HLA-DQB1*03:27', 'HLA-DQA1*05:11-HLA-DQB1*03:28',
         'HLA-DQA1*05:11-HLA-DQB1*03:29', 'HLA-DQA1*05:11-HLA-DQB1*03:30',
         'HLA-DQA1*05:11-HLA-DQB1*03:31', 'HLA-DQA1*05:11-HLA-DQB1*03:32', 'HLA-DQA1*05:11-HLA-DQB1*03:33', 'HLA-DQA1*05:11-HLA-DQB1*03:34',
         'HLA-DQA1*05:11-HLA-DQB1*03:35', 'HLA-DQA1*05:11-HLA-DQB1*03:36',
         'HLA-DQA1*05:11-HLA-DQB1*03:37', 'HLA-DQA1*05:11-HLA-DQB1*03:38', 'HLA-DQA1*05:11-HLA-DQB1*04:01', 'HLA-DQA1*05:11-HLA-DQB1*04:02',
         'HLA-DQA1*05:11-HLA-DQB1*04:03', 'HLA-DQA1*05:11-HLA-DQB1*04:04',
         'HLA-DQA1*05:11-HLA-DQB1*04:05', 'HLA-DQA1*05:11-HLA-DQB1*04:06', 'HLA-DQA1*05:11-HLA-DQB1*04:07', 'HLA-DQA1*05:11-HLA-DQB1*04:08',
         'HLA-DQA1*05:11-HLA-DQB1*05:01', 'HLA-DQA1*05:11-HLA-DQB1*05:02',
         'HLA-DQA1*05:11-HLA-DQB1*05:03', 'HLA-DQA1*05:11-HLA-DQB1*05:05', 'HLA-DQA1*05:11-HLA-DQB1*05:06', 'HLA-DQA1*05:11-HLA-DQB1*05:07',
         'HLA-DQA1*05:11-HLA-DQB1*05:08', 'HLA-DQA1*05:11-HLA-DQB1*05:09',
         'HLA-DQA1*05:11-HLA-DQB1*05:10', 'HLA-DQA1*05:11-HLA-DQB1*05:11', 'HLA-DQA1*05:11-HLA-DQB1*05:12', 'HLA-DQA1*05:11-HLA-DQB1*05:13',
         'HLA-DQA1*05:11-HLA-DQB1*05:14', 'HLA-DQA1*05:11-HLA-DQB1*06:01',
         'HLA-DQA1*05:11-HLA-DQB1*06:02', 'HLA-DQA1*05:11-HLA-DQB1*06:03', 'HLA-DQA1*05:11-HLA-DQB1*06:04', 'HLA-DQA1*05:11-HLA-DQB1*06:07',
         'HLA-DQA1*05:11-HLA-DQB1*06:08', 'HLA-DQA1*05:11-HLA-DQB1*06:09',
         'HLA-DQA1*05:11-HLA-DQB1*06:10', 'HLA-DQA1*05:11-HLA-DQB1*06:11', 'HLA-DQA1*05:11-HLA-DQB1*06:12', 'HLA-DQA1*05:11-HLA-DQB1*06:14',
         'HLA-DQA1*05:11-HLA-DQB1*06:15', 'HLA-DQA1*05:11-HLA-DQB1*06:16',
         'HLA-DQA1*05:11-HLA-DQB1*06:17', 'HLA-DQA1*05:11-HLA-DQB1*06:18', 'HLA-DQA1*05:11-HLA-DQB1*06:19', 'HLA-DQA1*05:11-HLA-DQB1*06:21',
         'HLA-DQA1*05:11-HLA-DQB1*06:22', 'HLA-DQA1*05:11-HLA-DQB1*06:23',
         'HLA-DQA1*05:11-HLA-DQB1*06:24', 'HLA-DQA1*05:11-HLA-DQB1*06:25', 'HLA-DQA1*05:11-HLA-DQB1*06:27', 'HLA-DQA1*05:11-HLA-DQB1*06:28',
         'HLA-DQA1*05:11-HLA-DQB1*06:29', 'HLA-DQA1*05:11-HLA-DQB1*06:30',
         'HLA-DQA1*05:11-HLA-DQB1*06:31', 'HLA-DQA1*05:11-HLA-DQB1*06:32', 'HLA-DQA1*05:11-HLA-DQB1*06:33', 'HLA-DQA1*05:11-HLA-DQB1*06:34',
         'HLA-DQA1*05:11-HLA-DQB1*06:35', 'HLA-DQA1*05:11-HLA-DQB1*06:36',
         'HLA-DQA1*05:11-HLA-DQB1*06:37', 'HLA-DQA1*05:11-HLA-DQB1*06:38', 'HLA-DQA1*05:11-HLA-DQB1*06:39', 'HLA-DQA1*05:11-HLA-DQB1*06:40',
         'HLA-DQA1*05:11-HLA-DQB1*06:41', 'HLA-DQA1*05:11-HLA-DQB1*06:42',
         'HLA-DQA1*05:11-HLA-DQB1*06:43', 'HLA-DQA1*05:11-HLA-DQB1*06:44', 'HLA-DQA1*06:01-HLA-DQB1*02:01', 'HLA-DQA1*06:01-HLA-DQB1*02:02',
         'HLA-DQA1*06:01-HLA-DQB1*02:03', 'HLA-DQA1*06:01-HLA-DQB1*02:04',
         'HLA-DQA1*06:01-HLA-DQB1*02:05', 'HLA-DQA1*06:01-HLA-DQB1*02:06', 'HLA-DQA1*06:01-HLA-DQB1*03:01', 'HLA-DQA1*06:01-HLA-DQB1*03:02',
         'HLA-DQA1*06:01-HLA-DQB1*03:03', 'HLA-DQA1*06:01-HLA-DQB1*03:04',
         'HLA-DQA1*06:01-HLA-DQB1*03:05', 'HLA-DQA1*06:01-HLA-DQB1*03:06', 'HLA-DQA1*06:01-HLA-DQB1*03:07', 'HLA-DQA1*06:01-HLA-DQB1*03:08',
         'HLA-DQA1*06:01-HLA-DQB1*03:09', 'HLA-DQA1*06:01-HLA-DQB1*03:10',
         'HLA-DQA1*06:01-HLA-DQB1*03:11', 'HLA-DQA1*06:01-HLA-DQB1*03:12', 'HLA-DQA1*06:01-HLA-DQB1*03:13', 'HLA-DQA1*06:01-HLA-DQB1*03:14',
         'HLA-DQA1*06:01-HLA-DQB1*03:15', 'HLA-DQA1*06:01-HLA-DQB1*03:16',
         'HLA-DQA1*06:01-HLA-DQB1*03:17', 'HLA-DQA1*06:01-HLA-DQB1*03:18', 'HLA-DQA1*06:01-HLA-DQB1*03:19', 'HLA-DQA1*06:01-HLA-DQB1*03:20',
         'HLA-DQA1*06:01-HLA-DQB1*03:21', 'HLA-DQA1*06:01-HLA-DQB1*03:22',
         'HLA-DQA1*06:01-HLA-DQB1*03:23', 'HLA-DQA1*06:01-HLA-DQB1*03:24', 'HLA-DQA1*06:01-HLA-DQB1*03:25', 'HLA-DQA1*06:01-HLA-DQB1*03:26',
         'HLA-DQA1*06:01-HLA-DQB1*03:27', 'HLA-DQA1*06:01-HLA-DQB1*03:28',
         'HLA-DQA1*06:01-HLA-DQB1*03:29', 'HLA-DQA1*06:01-HLA-DQB1*03:30', 'HLA-DQA1*06:01-HLA-DQB1*03:31', 'HLA-DQA1*06:01-HLA-DQB1*03:32',
         'HLA-DQA1*06:01-HLA-DQB1*03:33', 'HLA-DQA1*06:01-HLA-DQB1*03:34',
         'HLA-DQA1*06:01-HLA-DQB1*03:35', 'HLA-DQA1*06:01-HLA-DQB1*03:36', 'HLA-DQA1*06:01-HLA-DQB1*03:37', 'HLA-DQA1*06:01-HLA-DQB1*03:38',
         'HLA-DQA1*06:01-HLA-DQB1*04:01', 'HLA-DQA1*06:01-HLA-DQB1*04:02',
         'HLA-DQA1*06:01-HLA-DQB1*04:03', 'HLA-DQA1*06:01-HLA-DQB1*04:04', 'HLA-DQA1*06:01-HLA-DQB1*04:05', 'HLA-DQA1*06:01-HLA-DQB1*04:06',
         'HLA-DQA1*06:01-HLA-DQB1*04:07', 'HLA-DQA1*06:01-HLA-DQB1*04:08',
         'HLA-DQA1*06:01-HLA-DQB1*05:01', 'HLA-DQA1*06:01-HLA-DQB1*05:02', 'HLA-DQA1*06:01-HLA-DQB1*05:03', 'HLA-DQA1*06:01-HLA-DQB1*05:05',
         'HLA-DQA1*06:01-HLA-DQB1*05:06', 'HLA-DQA1*06:01-HLA-DQB1*05:07',
         'HLA-DQA1*06:01-HLA-DQB1*05:08', 'HLA-DQA1*06:01-HLA-DQB1*05:09', 'HLA-DQA1*06:01-HLA-DQB1*05:10', 'HLA-DQA1*06:01-HLA-DQB1*05:11',
         'HLA-DQA1*06:01-HLA-DQB1*05:12', 'HLA-DQA1*06:01-HLA-DQB1*05:13',
         'HLA-DQA1*06:01-HLA-DQB1*05:14', 'HLA-DQA1*06:01-HLA-DQB1*06:01', 'HLA-DQA1*06:01-HLA-DQB1*06:02', 'HLA-DQA1*06:01-HLA-DQB1*06:03',
         'HLA-DQA1*06:01-HLA-DQB1*06:04', 'HLA-DQA1*06:01-HLA-DQB1*06:07',
         'HLA-DQA1*06:01-HLA-DQB1*06:08', 'HLA-DQA1*06:01-HLA-DQB1*06:09', 'HLA-DQA1*06:01-HLA-DQB1*06:10', 'HLA-DQA1*06:01-HLA-DQB1*06:11',
         'HLA-DQA1*06:01-HLA-DQB1*06:12', 'HLA-DQA1*06:01-HLA-DQB1*06:14',
         'HLA-DQA1*06:01-HLA-DQB1*06:15', 'HLA-DQA1*06:01-HLA-DQB1*06:16', 'HLA-DQA1*06:01-HLA-DQB1*06:17', 'HLA-DQA1*06:01-HLA-DQB1*06:18',
         'HLA-DQA1*06:01-HLA-DQB1*06:19', 'HLA-DQA1*06:01-HLA-DQB1*06:21',
         'HLA-DQA1*06:01-HLA-DQB1*06:22', 'HLA-DQA1*06:01-HLA-DQB1*06:23', 'HLA-DQA1*06:01-HLA-DQB1*06:24', 'HLA-DQA1*06:01-HLA-DQB1*06:25',
         'HLA-DQA1*06:01-HLA-DQB1*06:27', 'HLA-DQA1*06:01-HLA-DQB1*06:28',
         'HLA-DQA1*06:01-HLA-DQB1*06:29', 'HLA-DQA1*06:01-HLA-DQB1*06:30', 'HLA-DQA1*06:01-HLA-DQB1*06:31', 'HLA-DQA1*06:01-HLA-DQB1*06:32',
         'HLA-DQA1*06:01-HLA-DQB1*06:33', 'HLA-DQA1*06:01-HLA-DQB1*06:34',
         'HLA-DQA1*06:01-HLA-DQB1*06:35', 'HLA-DQA1*06:01-HLA-DQB1*06:36', 'HLA-DQA1*06:01-HLA-DQB1*06:37', 'HLA-DQA1*06:01-HLA-DQB1*06:38',
         'HLA-DQA1*06:01-HLA-DQB1*06:39', 'HLA-DQA1*06:01-HLA-DQB1*06:40',
         'HLA-DQA1*06:01-HLA-DQB1*06:41', 'HLA-DQA1*06:01-HLA-DQB1*06:42', 'HLA-DQA1*06:01-HLA-DQB1*06:43', 'HLA-DQA1*06:01-HLA-DQB1*06:44',
         'HLA-DQA1*06:02-HLA-DQB1*02:01', 'HLA-DQA1*06:02-HLA-DQB1*02:02',
         'HLA-DQA1*06:02-HLA-DQB1*02:03', 'HLA-DQA1*06:02-HLA-DQB1*02:04', 'HLA-DQA1*06:02-HLA-DQB1*02:05', 'HLA-DQA1*06:02-HLA-DQB1*02:06',
         'HLA-DQA1*06:02-HLA-DQB1*03:01', 'HLA-DQA1*06:02-HLA-DQB1*03:02',
         'HLA-DQA1*06:02-HLA-DQB1*03:03', 'HLA-DQA1*06:02-HLA-DQB1*03:04', 'HLA-DQA1*06:02-HLA-DQB1*03:05', 'HLA-DQA1*06:02-HLA-DQB1*03:06',
         'HLA-DQA1*06:02-HLA-DQB1*03:07', 'HLA-DQA1*06:02-HLA-DQB1*03:08',
         'HLA-DQA1*06:02-HLA-DQB1*03:09', 'HLA-DQA1*06:02-HLA-DQB1*03:10', 'HLA-DQA1*06:02-HLA-DQB1*03:11', 'HLA-DQA1*06:02-HLA-DQB1*03:12',
         'HLA-DQA1*06:02-HLA-DQB1*03:13', 'HLA-DQA1*06:02-HLA-DQB1*03:14',
         'HLA-DQA1*06:02-HLA-DQB1*03:15', 'HLA-DQA1*06:02-HLA-DQB1*03:16', 'HLA-DQA1*06:02-HLA-DQB1*03:17', 'HLA-DQA1*06:02-HLA-DQB1*03:18',
         'HLA-DQA1*06:02-HLA-DQB1*03:19', 'HLA-DQA1*06:02-HLA-DQB1*03:20',
         'HLA-DQA1*06:02-HLA-DQB1*03:21', 'HLA-DQA1*06:02-HLA-DQB1*03:22', 'HLA-DQA1*06:02-HLA-DQB1*03:23', 'HLA-DQA1*06:02-HLA-DQB1*03:24',
         'HLA-DQA1*06:02-HLA-DQB1*03:25', 'HLA-DQA1*06:02-HLA-DQB1*03:26',
         'HLA-DQA1*06:02-HLA-DQB1*03:27', 'HLA-DQA1*06:02-HLA-DQB1*03:28', 'HLA-DQA1*06:02-HLA-DQB1*03:29', 'HLA-DQA1*06:02-HLA-DQB1*03:30',
         'HLA-DQA1*06:02-HLA-DQB1*03:31', 'HLA-DQA1*06:02-HLA-DQB1*03:32',
         'HLA-DQA1*06:02-HLA-DQB1*03:33', 'HLA-DQA1*06:02-HLA-DQB1*03:34', 'HLA-DQA1*06:02-HLA-DQB1*03:35', 'HLA-DQA1*06:02-HLA-DQB1*03:36',
         'HLA-DQA1*06:02-HLA-DQB1*03:37', 'HLA-DQA1*06:02-HLA-DQB1*03:38',
         'HLA-DQA1*06:02-HLA-DQB1*04:01', 'HLA-DQA1*06:02-HLA-DQB1*04:02', 'HLA-DQA1*06:02-HLA-DQB1*04:03', 'HLA-DQA1*06:02-HLA-DQB1*04:04',
         'HLA-DQA1*06:02-HLA-DQB1*04:05', 'HLA-DQA1*06:02-HLA-DQB1*04:06',
         'HLA-DQA1*06:02-HLA-DQB1*04:07', 'HLA-DQA1*06:02-HLA-DQB1*04:08', 'HLA-DQA1*06:02-HLA-DQB1*05:01', 'HLA-DQA1*06:02-HLA-DQB1*05:02',
         'HLA-DQA1*06:02-HLA-DQB1*05:03', 'HLA-DQA1*06:02-HLA-DQB1*05:05',
         'HLA-DQA1*06:02-HLA-DQB1*05:06', 'HLA-DQA1*06:02-HLA-DQB1*05:07', 'HLA-DQA1*06:02-HLA-DQB1*05:08', 'HLA-DQA1*06:02-HLA-DQB1*05:09',
         'HLA-DQA1*06:02-HLA-DQB1*05:10', 'HLA-DQA1*06:02-HLA-DQB1*05:11',
         'HLA-DQA1*06:02-HLA-DQB1*05:12', 'HLA-DQA1*06:02-HLA-DQB1*05:13', 'HLA-DQA1*06:02-HLA-DQB1*05:14', 'HLA-DQA1*06:02-HLA-DQB1*06:01',
         'HLA-DQA1*06:02-HLA-DQB1*06:02', 'HLA-DQA1*06:02-HLA-DQB1*06:03',
         'HLA-DQA1*06:02-HLA-DQB1*06:04', 'HLA-DQA1*06:02-HLA-DQB1*06:07', 'HLA-DQA1*06:02-HLA-DQB1*06:08', 'HLA-DQA1*06:02-HLA-DQB1*06:09',
         'HLA-DQA1*06:02-HLA-DQB1*06:10', 'HLA-DQA1*06:02-HLA-DQB1*06:11',
         'HLA-DQA1*06:02-HLA-DQB1*06:12', 'HLA-DQA1*06:02-HLA-DQB1*06:14', 'HLA-DQA1*06:02-HLA-DQB1*06:15', 'HLA-DQA1*06:02-HLA-DQB1*06:16',
         'HLA-DQA1*06:02-HLA-DQB1*06:17', 'HLA-DQA1*06:02-HLA-DQB1*06:18',
         'HLA-DQA1*06:02-HLA-DQB1*06:19', 'HLA-DQA1*06:02-HLA-DQB1*06:21', 'HLA-DQA1*06:02-HLA-DQB1*06:22', 'HLA-DQA1*06:02-HLA-DQB1*06:23',
         'HLA-DQA1*06:02-HLA-DQB1*06:24', 'HLA-DQA1*06:02-HLA-DQB1*06:25',
         'HLA-DQA1*06:02-HLA-DQB1*06:27', 'HLA-DQA1*06:02-HLA-DQB1*06:28', 'HLA-DQA1*06:02-HLA-DQB1*06:29', 'HLA-DQA1*06:02-HLA-DQB1*06:30',
         'HLA-DQA1*06:02-HLA-DQB1*06:31', 'HLA-DQA1*06:02-HLA-DQB1*06:32',
         'HLA-DQA1*06:02-HLA-DQB1*06:33', 'HLA-DQA1*06:02-HLA-DQB1*06:34', 'HLA-DQA1*06:02-HLA-DQB1*06:35', 'HLA-DQA1*06:02-HLA-DQB1*06:36',
         'HLA-DQA1*06:02-HLA-DQB1*06:37', 'HLA-DQA1*06:02-HLA-DQB1*06:38',
         'HLA-DQA1*06:02-HLA-DQB1*06:39', 'HLA-DQA1*06:02-HLA-DQB1*06:40', 'HLA-DQA1*06:02-HLA-DQB1*06:41', 'HLA-DQA1*06:02-HLA-DQB1*06:42',
         'HLA-DQA1*06:02-HLA-DQB1*06:43', 'HLA-DQA1*06:02-HLA-DQB1*06:44',
         'H2-IAb', 'H2-IAd'])
    __version = "3.0"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def command(self):
        """
        Defines the commandline call for external tool
        """
        return self.__command

    @property
    def supportedLength(self):
        """
        A list of supported :class:`~Fred2.Core.Peptide.Peptide` lengths
        """
        return self.__supported_length

    @property
    def supportedAlleles(self):
        """A list of valid :class:`~Fred2.Core.Allele.Allele` models"""
        return self.__alleles

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    def _represent(self, allele):
        """
        Internal function transforming an allele object into its representative string
        :param allele: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: str
        """
        if isinstance(allele, MouseAllele):
            return "H-2-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)
        elif isinstance(allele, CombinedAllele):
            return "HLA-%s%s%s-%s%s%s" % (allele.alpha_locus, allele.alpha_supertype, allele.alpha_subtype,
                                          allele.beta_locus, allele.beta_supertype, allele.beta_subtype)
        else:
            return "%s_%s%s" % (allele.locus, allele.supertype, allele.subtype)

    def convert_alleles(self, alleles):
        """
        Converts :class:`~Fred2.Core.Allele.Allele` into the internal :class:`~Fred2.Core.Allele.Allele` representation
        of the predictor and returns a string representation

        :param alleles: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`
        :rtype: list(str)
        """
        return [self._represent(a) for a in alleles]

    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        result = defaultdict(defaultdict)
        f = csv.reader(open(file, "r"), delimiter='\t')
        alleles = [x.replace("*", "_").replace(":", "") for x in set([x for x in next(f) if x != ""])]
        next(f)
        ic_pos = 3
        for row in f:
            pep_seq = row[1]
            for i, a in enumerate(alleles):
                result[a][pep_seq] = float(row[ic_pos + i * 3])
        return result

    def get_external_version(self, path=None):
        """
        Returns the external version of the tool by executing
        >{command} --version

        might be dependent on the method and has to be overwritten
        therefore it is declared abstract to enforce the user to
        overwrite the method. The function in the base class can be called
        with super()

        :param str path: Optional specification of executable path if deviant from :attr:`self.__command`
        :return: The external version of the tool or None if tool does not support versioning
        :rtype: str
        """
        return None

    def prepare_input(self, input, file):
        """
        Prepares input for external tools
        and writes them to _file in the specific format

        No return value!

        :param: list(str) input: The :class:`~Fred2.Core.Peptide.Peptide` sequences to write into file
        :param File file: File-handler to input file for external tool
        """
        file.write("\n".join(input))


class NetMHCIIpan_3_1(NetMHCIIpan_3_0):
    """
    Implementation of NetMHCIIpan 3.1 adapter.

    .. note::

        Andreatta, M., Karosiene, E., Rasmussen, M., Stryhn, A., Buus, S., & Nielsen, M. (2015). Accurate pan-specific
        prediction of peptide-MHC class II binding affinity with improved binding core identification.
        Immunogenetics, 1-10.
    """

    __supported_length = frozenset([9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
    __name = "netmhcIIpan"
    __command = "netMHCIIpan -f {peptides} -inptype 1 -a {alleles} {options} -xls -xlsfile {out}"
    __alleles = frozenset(
        ['HLA-DRB1*01:01', 'HLA-DRB1*01:02', 'HLA-DRB1*01:03', 'HLA-DRB1*01:04', 'HLA-DRB1*01:05', 'HLA-DRB1*01:06',
         'HLA-DRB1*01:07', 'HLA-DRB1*01:08', 'HLA-DRB1*01:09', 'HLA-DRB1*01:10', 'HLA-DRB1*01:11', 'HLA-DRB1*01:12',
         'HLA-DRB1*01:13', 'HLA-DRB1*01:14', 'HLA-DRB1*01:15', 'HLA-DRB1*01:16', 'HLA-DRB1*01:17', 'HLA-DRB1*01:18',
         'HLA-DRB1*01:19', 'HLA-DRB1*01:20', 'HLA-DRB1*01:21', 'HLA-DRB1*01:22', 'HLA-DRB1*01:23', 'HLA-DRB1*01:24',
         'HLA-DRB1*01:25', 'HLA-DRB1*01:26', 'HLA-DRB1*01:27', 'HLA-DRB1*01:28', 'HLA-DRB1*01:29', 'HLA-DRB1*01:30',
         'HLA-DRB1*01:31', 'HLA-DRB1*01:32', 'HLA-DRB1*03:01', 'HLA-DRB1*03:02', 'HLA-DRB1*03:03', 'HLA-DRB1*03:04',
         'HLA-DRB1*03:05', 'HLA-DRB1*03:06', 'HLA-DRB1*03:07', 'HLA-DRB1*03:08', 'HLA-DRB1*03:10', 'HLA-DRB1*03:11',
         'HLA-DRB1*03:13', 'HLA-DRB1*03:14', 'HLA-DRB1*03:15', 'HLA-DRB1*03:17', 'HLA-DRB1*03:18', 'HLA-DRB1*03:19',
         'HLA-DRB1*03:20', 'HLA-DRB1*03:21', 'HLA-DRB1*03:22', 'HLA-DRB1*03:23', 'HLA-DRB1*03:24', 'HLA-DRB1*03:25',
         'HLA-DRB1*03:26', 'HLA-DRB1*03:27', 'HLA-DRB1*03:28', 'HLA-DRB1*03:29', 'HLA-DRB1*03:30', 'HLA-DRB1*03:31',
         'HLA-DRB1*03:32', 'HLA-DRB1*03:33', 'HLA-DRB1*03:34', 'HLA-DRB1*03:35', 'HLA-DRB1*03:36', 'HLA-DRB1*03:37',
         'HLA-DRB1*03:38', 'HLA-DRB1*03:39', 'HLA-DRB1*03:40', 'HLA-DRB1*03:41', 'HLA-DRB1*03:42', 'HLA-DRB1*03:43',
         'HLA-DRB1*03:44', 'HLA-DRB1*03:45', 'HLA-DRB1*03:46', 'HLA-DRB1*03:47', 'HLA-DRB1*03:48', 'HLA-DRB1*03:49',
         'HLA-DRB1*03:50', 'HLA-DRB1*03:51', 'HLA-DRB1*03:52', 'HLA-DRB1*03:53', 'HLA-DRB1*03:54', 'HLA-DRB1*03:55',
         'HLA-DRB1*04:01', 'HLA-DRB1*04:02', 'HLA-DRB1*04:03', 'HLA-DRB1*04:04', 'HLA-DRB1*04:05', 'HLA-DRB1*04:06',
         'HLA-DRB1*04:07', 'HLA-DRB1*04:08', 'HLA-DRB1*04:09', 'HLA-DRB1*04:10', 'HLA-DRB1*04:11', 'HLA-DRB1*04:12',
         'HLA-DRB1*04:13', 'HLA-DRB1*04:14', 'HLA-DRB1*04:15', 'HLA-DRB1*04:16', 'HLA-DRB1*04:17', 'HLA-DRB1*04:18',
         'HLA-DRB1*04:19', 'HLA-DRB1*04:21', 'HLA-DRB1*04:22', 'HLA-DRB1*04:23', 'HLA-DRB1*04:24', 'HLA-DRB1*04:26',
         'HLA-DRB1*04:27', 'HLA-DRB1*04:28', 'HLA-DRB1*04:29', 'HLA-DRB1*04:30', 'HLA-DRB1*04:31', 'HLA-DRB1*04:33',
         'HLA-DRB1*04:34', 'HLA-DRB1*04:35', 'HLA-DRB1*04:36', 'HLA-DRB1*04:37', 'HLA-DRB1*04:38', 'HLA-DRB1*04:39',
         'HLA-DRB1*04:40', 'HLA-DRB1*04:41', 'HLA-DRB1*04:42', 'HLA-DRB1*04:43', 'HLA-DRB1*04:44', 'HLA-DRB1*04:45',
         'HLA-DRB1*04:46', 'HLA-DRB1*04:47', 'HLA-DRB1*04:48', 'HLA-DRB1*04:49', 'HLA-DRB1*04:50', 'HLA-DRB1*04:51',
         'HLA-DRB1*04:52', 'HLA-DRB1*04:53', 'HLA-DRB1*04:54', 'HLA-DRB1*04:55', 'HLA-DRB1*04:56', 'HLA-DRB1*04:57',
         'HLA-DRB1*04:58', 'HLA-DRB1*04:59', 'HLA-DRB1*04:60', 'HLA-DRB1*04:61', 'HLA-DRB1*04:62', 'HLA-DRB1*04:63',
         'HLA-DRB1*04:64', 'HLA-DRB1*04:65', 'HLA-DRB1*04:66', 'HLA-DRB1*04:67', 'HLA-DRB1*04:68', 'HLA-DRB1*04:69',
         'HLA-DRB1*04:70', 'HLA-DRB1*04:71', 'HLA-DRB1*04:72', 'HLA-DRB1*04:73', 'HLA-DRB1*04:74', 'HLA-DRB1*04:75',
         'HLA-DRB1*04:76', 'HLA-DRB1*04:77', 'HLA-DRB1*04:78', 'HLA-DRB1*04:79', 'HLA-DRB1*04:80', 'HLA-DRB1*04:82',
         'HLA-DRB1*04:83', 'HLA-DRB1*04:84', 'HLA-DRB1*04:85', 'HLA-DRB1*04:86', 'HLA-DRB1*04:87', 'HLA-DRB1*04:88',
         'HLA-DRB1*04:89', 'HLA-DRB1*04:91', 'HLA-DRB1*07:01', 'HLA-DRB1*07:03', 'HLA-DRB1*07:04', 'HLA-DRB1*07:05',
         'HLA-DRB1*07:06', 'HLA-DRB1*07:07', 'HLA-DRB1*07:08', 'HLA-DRB1*07:09', 'HLA-DRB1*07:11', 'HLA-DRB1*07:12',
         'HLA-DRB1*07:13', 'HLA-DRB1*07:14', 'HLA-DRB1*07:15', 'HLA-DRB1*07:16', 'HLA-DRB1*07:17', 'HLA-DRB1*07:19',
         'HLA-DRB1*08:01', 'HLA-DRB1*08:02', 'HLA-DRB1*08:03', 'HLA-DRB1*08:04', 'HLA-DRB1*08:05', 'HLA-DRB1*08:06',
         'HLA-DRB1*08:07', 'HLA-DRB1*08:08', 'HLA-DRB1*08:09', 'HLA-DRB1*08:10', 'HLA-DRB1*08:11', 'HLA-DRB1*08:12',
         'HLA-DRB1*08:13', 'HLA-DRB1*08:14', 'HLA-DRB1*08:15', 'HLA-DRB1*08:16', 'HLA-DRB1*08:18', 'HLA-DRB1*08:19',
         'HLA-DRB1*08:20', 'HLA-DRB1*08:21', 'HLA-DRB1*08:22', 'HLA-DRB1*08:23', 'HLA-DRB1*08:24', 'HLA-DRB1*08:25',
         'HLA-DRB1*08:26', 'HLA-DRB1*08:27', 'HLA-DRB1*08:28', 'HLA-DRB1*08:29', 'HLA-DRB1*08:30', 'HLA-DRB1*08:31',
         'HLA-DRB1*08:32', 'HLA-DRB1*08:33', 'HLA-DRB1*08:34', 'HLA-DRB1*08:35', 'HLA-DRB1*08:36', 'HLA-DRB1*08:37',
         'HLA-DRB1*08:38', 'HLA-DRB1*08:39', 'HLA-DRB1*08:40', 'HLA-DRB1*09:01', 'HLA-DRB1*09:02', 'HLA-DRB1*09:03',
         'HLA-DRB1*09:04', 'HLA-DRB1*09:05', 'HLA-DRB1*09:06', 'HLA-DRB1*09:07', 'HLA-DRB1*09:08', 'HLA-DRB1*09:09',
         'HLA-DRB1*10:01', 'HLA-DRB1*10:02', 'HLA-DRB1*10:03', 'HLA-DRB1*11:01', 'HLA-DRB1*11:02', 'HLA-DRB1*11:03',
         'HLA-DRB1*11:04', 'HLA-DRB1*11:05', 'HLA-DRB1*11:06', 'HLA-DRB1*11:07', 'HLA-DRB1*11:08', 'HLA-DRB1*11:09',
         'HLA-DRB1*11:10', 'HLA-DRB1*11:11', 'HLA-DRB1*11:12', 'HLA-DRB1*11:13', 'HLA-DRB1*11:14', 'HLA-DRB1*11:15',
         'HLA-DRB1*11:16', 'HLA-DRB1*11:17', 'HLA-DRB1*11:18', 'HLA-DRB1*11:19', 'HLA-DRB1*11:20', 'HLA-DRB1*11:21',
         'HLA-DRB1*11:24', 'HLA-DRB1*11:25', 'HLA-DRB1*11:27', 'HLA-DRB1*11:28', 'HLA-DRB1*11:29', 'HLA-DRB1*11:30',
         'HLA-DRB1*11:31', 'HLA-DRB1*11:32', 'HLA-DRB1*11:33', 'HLA-DRB1*11:34', 'HLA-DRB1*11:35', 'HLA-DRB1*11:36',
         'HLA-DRB1*11:37', 'HLA-DRB1*11:38', 'HLA-DRB1*11:39', 'HLA-DRB1*11:41', 'HLA-DRB1*11:42', 'HLA-DRB1*11:43',
         'HLA-DRB1*11:44', 'HLA-DRB1*11:45', 'HLA-DRB1*11:46', 'HLA-DRB1*11:47', 'HLA-DRB1*11:48', 'HLA-DRB1*11:49',
         'HLA-DRB1*11:50', 'HLA-DRB1*11:51', 'HLA-DRB1*11:52', 'HLA-DRB1*11:53', 'HLA-DRB1*11:54', 'HLA-DRB1*11:55',
         'HLA-DRB1*11:56', 'HLA-DRB1*11:57', 'HLA-DRB1*11:58', 'HLA-DRB1*11:59', 'HLA-DRB1*11:60', 'HLA-DRB1*11:61',
         'HLA-DRB1*11:62', 'HLA-DRB1*11:63', 'HLA-DRB1*11:64', 'HLA-DRB1*11:65', 'HLA-DRB1*11:66', 'HLA-DRB1*11:67',
         'HLA-DRB1*11:68', 'HLA-DRB1*11:69', 'HLA-DRB1*11:70', 'HLA-DRB1*11:72', 'HLA-DRB1*11:73', 'HLA-DRB1*11:74',
         'HLA-DRB1*11:75', 'HLA-DRB1*11:76', 'HLA-DRB1*11:77', 'HLA-DRB1*11:78', 'HLA-DRB1*11:79', 'HLA-DRB1*11:80',
         'HLA-DRB1*11:81', 'HLA-DRB1*11:82', 'HLA-DRB1*11:83', 'HLA-DRB1*11:84', 'HLA-DRB1*11:85', 'HLA-DRB1*11:86',
         'HLA-DRB1*11:87', 'HLA-DRB1*11:88', 'HLA-DRB1*11:89', 'HLA-DRB1*11:90', 'HLA-DRB1*11:91', 'HLA-DRB1*11:92',
         'HLA-DRB1*11:93', 'HLA-DRB1*11:94', 'HLA-DRB1*11:95', 'HLA-DRB1*11:96', 'HLA-DRB1*12:01', 'HLA-DRB1*12:02',
         'HLA-DRB1*12:03', 'HLA-DRB1*12:04', 'HLA-DRB1*12:05', 'HLA-DRB1*12:06', 'HLA-DRB1*12:07', 'HLA-DRB1*12:08',
         'HLA-DRB1*12:09', 'HLA-DRB1*12:10', 'HLA-DRB1*12:11', 'HLA-DRB1*12:12', 'HLA-DRB1*12:13', 'HLA-DRB1*12:14',
         'HLA-DRB1*12:15', 'HLA-DRB1*12:16', 'HLA-DRB1*12:17', 'HLA-DRB1*12:18', 'HLA-DRB1*12:19', 'HLA-DRB1*12:20',
         'HLA-DRB1*12:21', 'HLA-DRB1*12:22', 'HLA-DRB1*12:23', 'HLA-DRB1*13:01', 'HLA-DRB1*13:02', 'HLA-DRB1*13:03',
         'HLA-DRB1*13:04', 'HLA-DRB1*13:05', 'HLA-DRB1*13:06', 'HLA-DRB1*13:07', 'HLA-DRB1*13:08', 'HLA-DRB1*13:09',
         'HLA-DRB1*13:10', 'HLA-DRB1*13:100', 'HLA-DRB1*13:101', 'HLA-DRB1*13:11', 'HLA-DRB1*13:12', 'HLA-DRB1*13:13',
         'HLA-DRB1*13:14', 'HLA-DRB1*13:15', 'HLA-DRB1*13:16', 'HLA-DRB1*13:17', 'HLA-DRB1*13:18', 'HLA-DRB1*13:19',
         'HLA-DRB1*13:20', 'HLA-DRB1*13:21', 'HLA-DRB1*13:22', 'HLA-DRB1*13:23', 'HLA-DRB1*13:24', 'HLA-DRB1*13:26',
         'HLA-DRB1*13:27', 'HLA-DRB1*13:29', 'HLA-DRB1*13:30', 'HLA-DRB1*13:31', 'HLA-DRB1*13:32', 'HLA-DRB1*13:33',
         'HLA-DRB1*13:34', 'HLA-DRB1*13:35', 'HLA-DRB1*13:36', 'HLA-DRB1*13:37', 'HLA-DRB1*13:38', 'HLA-DRB1*13:39',
         'HLA-DRB1*13:41', 'HLA-DRB1*13:42', 'HLA-DRB1*13:43', 'HLA-DRB1*13:44', 'HLA-DRB1*13:46', 'HLA-DRB1*13:47',
         'HLA-DRB1*13:48', 'HLA-DRB1*13:49', 'HLA-DRB1*13:50', 'HLA-DRB1*13:51', 'HLA-DRB1*13:52', 'HLA-DRB1*13:53',
         'HLA-DRB1*13:54', 'HLA-DRB1*13:55', 'HLA-DRB1*13:56', 'HLA-DRB1*13:57', 'HLA-DRB1*13:58', 'HLA-DRB1*13:59',
         'HLA-DRB1*13:60', 'HLA-DRB1*13:61', 'HLA-DRB1*13:62', 'HLA-DRB1*13:63', 'HLA-DRB1*13:64', 'HLA-DRB1*13:65',
         'HLA-DRB1*13:66', 'HLA-DRB1*13:67', 'HLA-DRB1*13:68', 'HLA-DRB1*13:69', 'HLA-DRB1*13:70', 'HLA-DRB1*13:71',
         'HLA-DRB1*13:72', 'HLA-DRB1*13:73', 'HLA-DRB1*13:74', 'HLA-DRB1*13:75', 'HLA-DRB1*13:76', 'HLA-DRB1*13:77',
         'HLA-DRB1*13:78', 'HLA-DRB1*13:79', 'HLA-DRB1*13:80', 'HLA-DRB1*13:81', 'HLA-DRB1*13:82', 'HLA-DRB1*13:83',
         'HLA-DRB1*13:84', 'HLA-DRB1*13:85', 'HLA-DRB1*13:86', 'HLA-DRB1*13:87', 'HLA-DRB1*13:88', 'HLA-DRB1*13:89',
         'HLA-DRB1*13:90', 'HLA-DRB1*13:91', 'HLA-DRB1*13:92', 'HLA-DRB1*13:93', 'HLA-DRB1*13:94', 'HLA-DRB1*13:95',
         'HLA-DRB1*13:96', 'HLA-DRB1*13:97', 'HLA-DRB1*13:98', 'HLA-DRB1*13:99', 'HLA-DRB1*14:01', 'HLA-DRB1*14:02',
         'HLA-DRB1*14:03', 'HLA-DRB1*14:04', 'HLA-DRB1*14:05', 'HLA-DRB1*14:06', 'HLA-DRB1*14:07', 'HLA-DRB1*14:08',
         'HLA-DRB1*14:09', 'HLA-DRB1*14:10', 'HLA-DRB1*14:11', 'HLA-DRB1*14:12', 'HLA-DRB1*14:13', 'HLA-DRB1*14:14',
         'HLA-DRB1*14:15', 'HLA-DRB1*14:16', 'HLA-DRB1*14:17', 'HLA-DRB1*14:18', 'HLA-DRB1*14:19', 'HLA-DRB1*14:20',
         'HLA-DRB1*14:21', 'HLA-DRB1*14:22', 'HLA-DRB1*14:23', 'HLA-DRB1*14:24', 'HLA-DRB1*14:25', 'HLA-DRB1*14:26',
         'HLA-DRB1*14:27', 'HLA-DRB1*14:28', 'HLA-DRB1*14:29', 'HLA-DRB1*14:30', 'HLA-DRB1*14:31', 'HLA-DRB1*14:32',
         'HLA-DRB1*14:33', 'HLA-DRB1*14:34', 'HLA-DRB1*14:35', 'HLA-DRB1*14:36', 'HLA-DRB1*14:37', 'HLA-DRB1*14:38',
         'HLA-DRB1*14:39', 'HLA-DRB1*14:40', 'HLA-DRB1*14:41', 'HLA-DRB1*14:42', 'HLA-DRB1*14:43', 'HLA-DRB1*14:44',
         'HLA-DRB1*14:45', 'HLA-DRB1*14:46', 'HLA-DRB1*14:47', 'HLA-DRB1*14:48', 'HLA-DRB1*14:49', 'HLA-DRB1*14:50',
         'HLA-DRB1*14:51', 'HLA-DRB1*14:52', 'HLA-DRB1*14:53', 'HLA-DRB1*14:54', 'HLA-DRB1*14:55', 'HLA-DRB1*14:56',
         'HLA-DRB1*14:57', 'HLA-DRB1*14:58', 'HLA-DRB1*14:59', 'HLA-DRB1*14:60', 'HLA-DRB1*14:61', 'HLA-DRB1*14:62',
         'HLA-DRB1*14:63', 'HLA-DRB1*14:64', 'HLA-DRB1*14:65', 'HLA-DRB1*14:67', 'HLA-DRB1*14:68', 'HLA-DRB1*14:69',
         'HLA-DRB1*14:70', 'HLA-DRB1*14:71', 'HLA-DRB1*14:72', 'HLA-DRB1*14:73', 'HLA-DRB1*14:74', 'HLA-DRB1*14:75',
         'HLA-DRB1*14:76', 'HLA-DRB1*14:77', 'HLA-DRB1*14:78', 'HLA-DRB1*14:79', 'HLA-DRB1*14:80', 'HLA-DRB1*14:81',
         'HLA-DRB1*14:82', 'HLA-DRB1*14:83', 'HLA-DRB1*14:84', 'HLA-DRB1*14:85', 'HLA-DRB1*14:86', 'HLA-DRB1*14:87',
         'HLA-DRB1*14:88', 'HLA-DRB1*14:89', 'HLA-DRB1*14:90', 'HLA-DRB1*14:91', 'HLA-DRB1*14:93', 'HLA-DRB1*14:94',
         'HLA-DRB1*14:95', 'HLA-DRB1*14:96', 'HLA-DRB1*14:97', 'HLA-DRB1*14:98', 'HLA-DRB1*14:99', 'HLA-DRB1*15:01',
         'HLA-DRB1*15:02', 'HLA-DRB1*15:03', 'HLA-DRB1*15:04', 'HLA-DRB1*15:05', 'HLA-DRB1*15:06', 'HLA-DRB1*15:07',
         'HLA-DRB1*15:08', 'HLA-DRB1*15:09', 'HLA-DRB1*15:10', 'HLA-DRB1*15:11', 'HLA-DRB1*15:12', 'HLA-DRB1*15:13',
         'HLA-DRB1*15:14', 'HLA-DRB1*15:15', 'HLA-DRB1*15:16', 'HLA-DRB1*15:18', 'HLA-DRB1*15:19', 'HLA-DRB1*15:20',
         'HLA-DRB1*15:21', 'HLA-DRB1*15:22', 'HLA-DRB1*15:23', 'HLA-DRB1*15:24', 'HLA-DRB1*15:25', 'HLA-DRB1*15:26',
         'HLA-DRB1*15:27', 'HLA-DRB1*15:28', 'HLA-DRB1*15:29', 'HLA-DRB1*15:30', 'HLA-DRB1*15:31', 'HLA-DRB1*15:32',
         'HLA-DRB1*15:33', 'HLA-DRB1*15:34', 'HLA-DRB1*15:35', 'HLA-DRB1*15:36', 'HLA-DRB1*15:37', 'HLA-DRB1*15:38',
         'HLA-DRB1*15:39', 'HLA-DRB1*15:40', 'HLA-DRB1*15:41', 'HLA-DRB1*15:42', 'HLA-DRB1*15:43', 'HLA-DRB1*15:44',
         'HLA-DRB1*15:45', 'HLA-DRB1*15:46', 'HLA-DRB1*15:47', 'HLA-DRB1*15:48', 'HLA-DRB1*15:49', 'HLA-DRB1*16:01',
         'HLA-DRB1*16:02', 'HLA-DRB1*16:03', 'HLA-DRB1*16:04', 'HLA-DRB1*16:05', 'HLA-DRB1*16:07', 'HLA-DRB1*16:08',
         'HLA-DRB1*16:09', 'HLA-DRB1*16:10', 'HLA-DRB1*16:11', 'HLA-DRB1*16:12', 'HLA-DRB1*16:14', 'HLA-DRB1*16:15',
         'HLA-DRB1*16:16', 'HLA-DRB3*01:01', 'HLA-DRB3*01:04', 'HLA-DRB3*01:05', 'HLA-DRB3*01:08', 'HLA-DRB3*01:09',
         'HLA-DRB3*01:11', 'HLA-DRB3*01:12', 'HLA-DRB3*01:13', 'HLA-DRB3*01:14', 'HLA-DRB3*02:01', 'HLA-DRB3*02:02',
         'HLA-DRB3*02:04', 'HLA-DRB3*02:05', 'HLA-DRB3*02:09', 'HLA-DRB3*02:10', 'HLA-DRB3*02:11', 'HLA-DRB3*02:12',
         'HLA-DRB3*02:13', 'HLA-DRB3*02:14', 'HLA-DRB3*02:15', 'HLA-DRB3*02:16', 'HLA-DRB3*02:17', 'HLA-DRB3*02:18',
         'HLA-DRB3*02:19', 'HLA-DRB3*02:20', 'HLA-DRB3*02:21', 'HLA-DRB3*02:22', 'HLA-DRB3*02:23', 'HLA-DRB3*02:24',
         'HLA-DRB3*02:25', 'HLA-DRB3*03:01', 'HLA-DRB3*03:03', 'HLA-DRB4*01:01', 'HLA-DRB4*01:03', 'HLA-DRB4*01:04',
         'HLA-DRB4*01:06', 'HLA-DRB4*01:07', 'HLA-DRB4*01:08', 'HLA-DRB5*01:01', 'HLA-DRB5*01:02', 'HLA-DRB5*01:03',
         'HLA-DRB5*01:04', 'HLA-DRB5*01:05', 'HLA-DRB5*01:06', 'HLA-DRB5*01:08N', 'HLA-DRB5*01:11', 'HLA-DRB5*01:12',
         'HLA-DRB5*01:13', 'HLA-DRB5*01:14', 'HLA-DRB5*02:02', 'HLA-DRB5*02:03', 'HLA-DRB5*02:04', 'HLA-DRB5*02:05',
         'HLA-DPA1*01:03-HLA-DPB1*01:01', 'HLA-DPA1*01:03-HLA-DPB1*02:01', 'HLA-DPA1*01:03-HLA-DPB1*02:02', 'HLA-DPA1*01:03-HLA-DPB1*03:01',
         'HLA-DPA1*01:03-HLA-DPB1*04:01', 'HLA-DPA1*01:03-HLA-DPB1*04:02', 'HLA-DPA1*01:03-HLA-DPB1*05:01', 'HLA-DPA1*01:03-HLA-DPB1*06:01',
         'HLA-DPA1*01:03-HLA-DPB1*08:01', 'HLA-DPA1*01:03-HLA-DPB1*09:01', 'HLA-DPA1*01:03-HLA-DPB1*10:001', 'HLA-DPA1*01:03-HLA-DPB1*10:01',
         'HLA-DPA1*01:03-HLA-DPB1*10:101', 'HLA-DPA1*01:03-HLA-DPB1*10:201',
         'HLA-DPA1*01:03-HLA-DPB1*10:301', 'HLA-DPA1*01:03-HLA-DPB1*10:401',
         'HLA-DPA1*01:03-HLA-DPB1*10:501', 'HLA-DPA1*01:03-HLA-DPB1*10:601', 'HLA-DPA1*01:03-HLA-DPB1*10:701', 'HLA-DPA1*01:03-HLA-DPB1*10:801',
         'HLA-DPA1*01:03-HLA-DPB1*10:901', 'HLA-DPA1*01:03-HLA-DPB1*11:001',
         'HLA-DPA1*01:03-HLA-DPB1*11:01', 'HLA-DPA1*01:03-HLA-DPB1*11:101', 'HLA-DPA1*01:03-HLA-DPB1*11:201', 'HLA-DPA1*01:03-HLA-DPB1*11:301',
         'HLA-DPA1*01:03-HLA-DPB1*11:401', 'HLA-DPA1*01:03-HLA-DPB1*11:501',
         'HLA-DPA1*01:03-HLA-DPB1*11:601', 'HLA-DPA1*01:03-HLA-DPB1*11:701', 'HLA-DPA1*01:03-HLA-DPB1*11:801', 'HLA-DPA1*01:03-HLA-DPB1*11:901',
         'HLA-DPA1*01:03-HLA-DPB1*12:101', 'HLA-DPA1*01:03-HLA-DPB1*12:201',
         'HLA-DPA1*01:03-HLA-DPB1*12:301', 'HLA-DPA1*01:03-HLA-DPB1*12:401', 'HLA-DPA1*01:03-HLA-DPB1*12:501', 'HLA-DPA1*01:03-HLA-DPB1*12:601',
         'HLA-DPA1*01:03-HLA-DPB1*12:701', 'HLA-DPA1*01:03-HLA-DPB1*12:801',
         'HLA-DPA1*01:03-HLA-DPB1*12:901', 'HLA-DPA1*01:03-HLA-DPB1*13:001', 'HLA-DPA1*01:03-HLA-DPB1*13:01', 'HLA-DPA1*01:03-HLA-DPB1*13:101',
         'HLA-DPA1*01:03-HLA-DPB1*13:201', 'HLA-DPA1*01:03-HLA-DPB1*13:301',
         'HLA-DPA1*01:03-HLA-DPB1*13:401', 'HLA-DPA1*01:03-HLA-DPB1*14:01', 'HLA-DPA1*01:03-HLA-DPB1*15:01', 'HLA-DPA1*01:03-HLA-DPB1*16:01',
         'HLA-DPA1*01:03-HLA-DPB1*17:01', 'HLA-DPA1*01:03-HLA-DPB1*18:01',
         'HLA-DPA1*01:03-HLA-DPB1*19:01', 'HLA-DPA1*01:03-HLA-DPB1*20:01', 'HLA-DPA1*01:03-HLA-DPB1*21:01', 'HLA-DPA1*01:03-HLA-DPB1*22:01',
         'HLA-DPA1*01:03-HLA-DPB1*23:01', 'HLA-DPA1*01:03-HLA-DPB1*24:01',
         'HLA-DPA1*01:03-HLA-DPB1*25:01', 'HLA-DPA1*01:03-HLA-DPB1*26:01', 'HLA-DPA1*01:03-HLA-DPB1*27:01', 'HLA-DPA1*01:03-HLA-DPB1*28:01',
         'HLA-DPA1*01:03-HLA-DPB1*29:01', 'HLA-DPA1*01:03-HLA-DPB1*30:01',
         'HLA-DPA1*01:03-HLA-DPB1*31:01', 'HLA-DPA1*01:03-HLA-DPB1*32:01', 'HLA-DPA1*01:03-HLA-DPB1*33:01', 'HLA-DPA1*01:03-HLA-DPB1*34:01',
         'HLA-DPA1*01:03-HLA-DPB1*35:01', 'HLA-DPA1*01:03-HLA-DPB1*36:01',
         'HLA-DPA1*01:03-HLA-DPB1*37:01', 'HLA-DPA1*01:03-HLA-DPB1*38:01', 'HLA-DPA1*01:03-HLA-DPB1*39:01', 'HLA-DPA1*01:03-HLA-DPB1*40:01',
         'HLA-DPA1*01:03-HLA-DPB1*41:01', 'HLA-DPA1*01:03-HLA-DPB1*44:01',
         'HLA-DPA1*01:03-HLA-DPB1*45:01', 'HLA-DPA1*01:03-HLA-DPB1*46:01', 'HLA-DPA1*01:03-HLA-DPB1*47:01', 'HLA-DPA1*01:03-HLA-DPB1*48:01',
         'HLA-DPA1*01:03-HLA-DPB1*49:01', 'HLA-DPA1*01:03-HLA-DPB1*50:01',
         'HLA-DPA1*01:03-HLA-DPB1*51:01', 'HLA-DPA1*01:03-HLA-DPB1*52:01', 'HLA-DPA1*01:03-HLA-DPB1*53:01', 'HLA-DPA1*01:03-HLA-DPB1*54:01',
         'HLA-DPA1*01:03-HLA-DPB1*55:01', 'HLA-DPA1*01:03-HLA-DPB1*56:01',
         'HLA-DPA1*01:03-HLA-DPB1*58:01', 'HLA-DPA1*01:03-HLA-DPB1*59:01', 'HLA-DPA1*01:03-HLA-DPB1*60:01', 'HLA-DPA1*01:03-HLA-DPB1*62:01',
         'HLA-DPA1*01:03-HLA-DPB1*63:01', 'HLA-DPA1*01:03-HLA-DPB1*65:01',
         'HLA-DPA1*01:03-HLA-DPB1*66:01', 'HLA-DPA1*01:03-HLA-DPB1*67:01', 'HLA-DPA1*01:03-HLA-DPB1*68:01', 'HLA-DPA1*01:03-HLA-DPB1*69:01',
         'HLA-DPA1*01:03-HLA-DPB1*70:01', 'HLA-DPA1*01:03-HLA-DPB1*71:01',
         'HLA-DPA1*01:03-HLA-DPB1*72:01', 'HLA-DPA1*01:03-HLA-DPB1*73:01', 'HLA-DPA1*01:03-HLA-DPB1*74:01', 'HLA-DPA1*01:03-HLA-DPB1*75:01',
         'HLA-DPA1*01:03-HLA-DPB1*76:01', 'HLA-DPA1*01:03-HLA-DPB1*77:01',
         'HLA-DPA1*01:03-HLA-DPB1*78:01', 'HLA-DPA1*01:03-HLA-DPB1*79:01', 'HLA-DPA1*01:03-HLA-DPB1*80:01', 'HLA-DPA1*01:03-HLA-DPB1*81:01',
         'HLA-DPA1*01:03-HLA-DPB1*82:01', 'HLA-DPA1*01:03-HLA-DPB1*83:01',
         'HLA-DPA1*01:03-HLA-DPB1*84:01', 'HLA-DPA1*01:03-HLA-DPB1*85:01', 'HLA-DPA1*01:03-HLA-DPB1*86:01', 'HLA-DPA1*01:03-HLA-DPB1*87:01',
         'HLA-DPA1*01:03-HLA-DPB1*88:01', 'HLA-DPA1*01:03-HLA-DPB1*89:01',
         'HLA-DPA1*01:03-HLA-DPB1*90:01', 'HLA-DPA1*01:03-HLA-DPB1*91:01', 'HLA-DPA1*01:03-HLA-DPB1*92:01', 'HLA-DPA1*01:03-HLA-DPB1*93:01',
         'HLA-DPA1*01:03-HLA-DPB1*94:01', 'HLA-DPA1*01:03-HLA-DPB1*95:01',
         'HLA-DPA1*01:03-HLA-DPB1*96:01', 'HLA-DPA1*01:03-HLA-DPB1*97:01', 'HLA-DPA1*01:03-HLA-DPB1*98:01', 'HLA-DPA1*01:03-HLA-DPB1*99:01',
         'HLA-DPA1*01:04-HLA-DPB1*01:01', 'HLA-DPA1*01:04-HLA-DPB1*02:01',
         'HLA-DPA1*01:04-HLA-DPB1*02:02', 'HLA-DPA1*01:04-HLA-DPB1*03:01', 'HLA-DPA1*01:04-HLA-DPB1*04:01', 'HLA-DPA1*01:04-HLA-DPB1*04:02',
         'HLA-DPA1*01:04-HLA-DPB1*05:01', 'HLA-DPA1*01:04-HLA-DPB1*06:01',
         'HLA-DPA1*01:04-HLA-DPB1*08:01', 'HLA-DPA1*01:04-HLA-DPB1*09:01', 'HLA-DPA1*01:04-HLA-DPB1*10:001', 'HLA-DPA1*01:04-HLA-DPB1*10:01',
         'HLA-DPA1*01:04-HLA-DPB1*10:101', 'HLA-DPA1*01:04-HLA-DPB1*10:201',
         'HLA-DPA1*01:04-HLA-DPB1*10:301', 'HLA-DPA1*01:04-HLA-DPB1*10:401', 'HLA-DPA1*01:04-HLA-DPB1*10:501', 'HLA-DPA1*01:04-HLA-DPB1*10:601',
         'HLA-DPA1*01:04-HLA-DPB1*10:701', 'HLA-DPA1*01:04-HLA-DPB1*10:801',
         'HLA-DPA1*01:04-HLA-DPB1*10:901', 'HLA-DPA1*01:04-HLA-DPB1*11:001', 'HLA-DPA1*01:04-HLA-DPB1*11:01', 'HLA-DPA1*01:04-HLA-DPB1*11:101',
         'HLA-DPA1*01:04-HLA-DPB1*11:201', 'HLA-DPA1*01:04-HLA-DPB1*11:301',
         'HLA-DPA1*01:04-HLA-DPB1*11:401', 'HLA-DPA1*01:04-HLA-DPB1*11:501', 'HLA-DPA1*01:04-HLA-DPB1*11:601', 'HLA-DPA1*01:04-HLA-DPB1*11:701',
         'HLA-DPA1*01:04-HLA-DPB1*11:801', 'HLA-DPA1*01:04-HLA-DPB1*11:901',
         'HLA-DPA1*01:04-HLA-DPB1*12:101', 'HLA-DPA1*01:04-HLA-DPB1*12:201', 'HLA-DPA1*01:04-HLA-DPB1*12:301', 'HLA-DPA1*01:04-HLA-DPB1*12:401',
         'HLA-DPA1*01:04-HLA-DPB1*12:501', 'HLA-DPA1*01:04-HLA-DPB1*12:601',
         'HLA-DPA1*01:04-HLA-DPB1*12:701', 'HLA-DPA1*01:04-HLA-DPB1*12:801', 'HLA-DPA1*01:04-HLA-DPB1*12:901', 'HLA-DPA1*01:04-HLA-DPB1*13:001',
         'HLA-DPA1*01:04-HLA-DPB1*13:01', 'HLA-DPA1*01:04-HLA-DPB1*13:101',
         'HLA-DPA1*01:04-HLA-DPB1*13:201', 'HLA-DPA1*01:04-HLA-DPB1*13:301', 'HLA-DPA1*01:04-HLA-DPB1*13:401', 'HLA-DPA1*01:04-HLA-DPB1*14:01',
         'HLA-DPA1*01:04-HLA-DPB1*15:01', 'HLA-DPA1*01:04-HLA-DPB1*16:01',
         'HLA-DPA1*01:04-HLA-DPB1*17:01', 'HLA-DPA1*01:04-HLA-DPB1*18:01', 'HLA-DPA1*01:04-HLA-DPB1*19:01', 'HLA-DPA1*01:04-HLA-DPB1*20:01',
         'HLA-DPA1*01:04-HLA-DPB1*21:01', 'HLA-DPA1*01:04-HLA-DPB1*22:01',
         'HLA-DPA1*01:04-HLA-DPB1*23:01', 'HLA-DPA1*01:04-HLA-DPB1*24:01', 'HLA-DPA1*01:04-HLA-DPB1*25:01', 'HLA-DPA1*01:04-HLA-DPB1*26:01',
         'HLA-DPA1*01:04-HLA-DPB1*27:01', 'HLA-DPA1*01:04-HLA-DPB1*28:01',
         'HLA-DPA1*01:04-HLA-DPB1*29:01', 'HLA-DPA1*01:04-HLA-DPB1*30:01', 'HLA-DPA1*01:04-HLA-DPB1*31:01', 'HLA-DPA1*01:04-HLA-DPB1*32:01',
         'HLA-DPA1*01:04-HLA-DPB1*33:01', 'HLA-DPA1*01:04-HLA-DPB1*34:01',
         'HLA-DPA1*01:04-HLA-DPB1*35:01', 'HLA-DPA1*01:04-HLA-DPB1*36:01', 'HLA-DPA1*01:04-HLA-DPB1*37:01', 'HLA-DPA1*01:04-HLA-DPB1*38:01',
         'HLA-DPA1*01:04-HLA-DPB1*39:01', 'HLA-DPA1*01:04-HLA-DPB1*40:01',
         'HLA-DPA1*01:04-HLA-DPB1*41:01', 'HLA-DPA1*01:04-HLA-DPB1*44:01', 'HLA-DPA1*01:04-HLA-DPB1*45:01', 'HLA-DPA1*01:04-HLA-DPB1*46:01',
         'HLA-DPA1*01:04-HLA-DPB1*47:01', 'HLA-DPA1*01:04-HLA-DPB1*48:01',
         'HLA-DPA1*01:04-HLA-DPB1*49:01', 'HLA-DPA1*01:04-HLA-DPB1*50:01', 'HLA-DPA1*01:04-HLA-DPB1*51:01', 'HLA-DPA1*01:04-HLA-DPB1*52:01',
         'HLA-DPA1*01:04-HLA-DPB1*53:01', 'HLA-DPA1*01:04-HLA-DPB1*54:01',
         'HLA-DPA1*01:04-HLA-DPB1*55:01', 'HLA-DPA1*01:04-HLA-DPB1*56:01', 'HLA-DPA1*01:04-HLA-DPB1*58:01', 'HLA-DPA1*01:04-HLA-DPB1*59:01',
         'HLA-DPA1*01:04-HLA-DPB1*60:01', 'HLA-DPA1*01:04-HLA-DPB1*62:01',
         'HLA-DPA1*01:04-HLA-DPB1*63:01', 'HLA-DPA1*01:04-HLA-DPB1*65:01', 'HLA-DPA1*01:04-HLA-DPB1*66:01', 'HLA-DPA1*01:04-HLA-DPB1*67:01',
         'HLA-DPA1*01:04-HLA-DPB1*68:01', 'HLA-DPA1*01:04-HLA-DPB1*69:01',
         'HLA-DPA1*01:04-HLA-DPB1*70:01', 'HLA-DPA1*01:04-HLA-DPB1*71:01', 'HLA-DPA1*01:04-HLA-DPB1*72:01', 'HLA-DPA1*01:04-HLA-DPB1*73:01',
         'HLA-DPA1*01:04-HLA-DPB1*74:01', 'HLA-DPA1*01:04-HLA-DPB1*75:01',
         'HLA-DPA1*01:04-HLA-DPB1*76:01', 'HLA-DPA1*01:04-HLA-DPB1*77:01', 'HLA-DPA1*01:04-HLA-DPB1*78:01', 'HLA-DPA1*01:04-HLA-DPB1*79:01',
         'HLA-DPA1*01:04-HLA-DPB1*80:01', 'HLA-DPA1*01:04-HLA-DPB1*81:01',
         'HLA-DPA1*01:04-HLA-DPB1*82:01', 'HLA-DPA1*01:04-HLA-DPB1*83:01', 'HLA-DPA1*01:04-HLA-DPB1*84:01', 'HLA-DPA1*01:04-HLA-DPB1*85:01',
         'HLA-DPA1*01:04-HLA-DPB1*86:01', 'HLA-DPA1*01:04-HLA-DPB1*87:01',
         'HLA-DPA1*01:04-HLA-DPB1*88:01', 'HLA-DPA1*01:04-HLA-DPB1*89:01', 'HLA-DPA1*01:04-HLA-DPB1*90:01', 'HLA-DPA1*01:04-HLA-DPB1*91:01',
         'HLA-DPA1*01:04-HLA-DPB1*92:01', 'HLA-DPA1*01:04-HLA-DPB1*93:01',
         'HLA-DPA1*01:04-HLA-DPB1*94:01', 'HLA-DPA1*01:04-HLA-DPB1*95:01', 'HLA-DPA1*01:04-HLA-DPB1*96:01', 'HLA-DPA1*01:04-HLA-DPB1*97:01',
         'HLA-DPA1*01:04-HLA-DPB1*98:01', 'HLA-DPA1*01:04-HLA-DPB1*99:01',
         'HLA-DPA1*01:05-HLA-DPB1*01:01', 'HLA-DPA1*01:05-HLA-DPB1*02:01', 'HLA-DPA1*01:05-HLA-DPB1*02:02', 'HLA-DPA1*01:05-HLA-DPB1*03:01',
         'HLA-DPA1*01:05-HLA-DPB1*04:01', 'HLA-DPA1*01:05-HLA-DPB1*04:02',
         'HLA-DPA1*01:05-HLA-DPB1*05:01', 'HLA-DPA1*01:05-HLA-DPB1*06:01', 'HLA-DPA1*01:05-HLA-DPB1*08:01', 'HLA-DPA1*01:05-HLA-DPB1*09:01',
         'HLA-DPA1*01:05-HLA-DPB1*10:001', 'HLA-DPA1*01:05-HLA-DPB1*10:01',
         'HLA-DPA1*01:05-HLA-DPB1*10:101', 'HLA-DPA1*01:05-HLA-DPB1*10:201', 'HLA-DPA1*01:05-HLA-DPB1*10:301', 'HLA-DPA1*01:05-HLA-DPB1*10:401',
         'HLA-DPA1*01:05-HLA-DPB1*10:501', 'HLA-DPA1*01:05-HLA-DPB1*10:601',
         'HLA-DPA1*01:05-HLA-DPB1*10:701', 'HLA-DPA1*01:05-HLA-DPB1*10:801', 'HLA-DPA1*01:05-HLA-DPB1*10:901', 'HLA-DPA1*01:05-HLA-DPB1*11:001',
         'HLA-DPA1*01:05-HLA-DPB1*11:01', 'HLA-DPA1*01:05-HLA-DPB1*11:101',
         'HLA-DPA1*01:05-HLA-DPB1*11:201', 'HLA-DPA1*01:05-HLA-DPB1*11:301', 'HLA-DPA1*01:05-HLA-DPB1*11:401', 'HLA-DPA1*01:05-HLA-DPB1*11:501',
         'HLA-DPA1*01:05-HLA-DPB1*11:601', 'HLA-DPA1*01:05-HLA-DPB1*11:701',
         'HLA-DPA1*01:05-HLA-DPB1*11:801', 'HLA-DPA1*01:05-HLA-DPB1*11:901', 'HLA-DPA1*01:05-HLA-DPB1*12:101', 'HLA-DPA1*01:05-HLA-DPB1*12:201',
         'HLA-DPA1*01:05-HLA-DPB1*12:301', 'HLA-DPA1*01:05-HLA-DPB1*12:401',
         'HLA-DPA1*01:05-HLA-DPB1*12:501', 'HLA-DPA1*01:05-HLA-DPB1*12:601', 'HLA-DPA1*01:05-HLA-DPB1*12:701', 'HLA-DPA1*01:05-HLA-DPB1*12:801',
         'HLA-DPA1*01:05-HLA-DPB1*12:901', 'HLA-DPA1*01:05-HLA-DPB1*13:001',
         'HLA-DPA1*01:05-HLA-DPB1*13:01', 'HLA-DPA1*01:05-HLA-DPB1*13:101', 'HLA-DPA1*01:05-HLA-DPB1*13:201', 'HLA-DPA1*01:05-HLA-DPB1*13:301',
         'HLA-DPA1*01:05-HLA-DPB1*13:401', 'HLA-DPA1*01:05-HLA-DPB1*14:01',
         'HLA-DPA1*01:05-HLA-DPB1*15:01', 'HLA-DPA1*01:05-HLA-DPB1*16:01', 'HLA-DPA1*01:05-HLA-DPB1*17:01', 'HLA-DPA1*01:05-HLA-DPB1*18:01',
         'HLA-DPA1*01:05-HLA-DPB1*19:01', 'HLA-DPA1*01:05-HLA-DPB1*20:01',
         'HLA-DPA1*01:05-HLA-DPB1*21:01', 'HLA-DPA1*01:05-HLA-DPB1*22:01', 'HLA-DPA1*01:05-HLA-DPB1*23:01', 'HLA-DPA1*01:05-HLA-DPB1*24:01',
         'HLA-DPA1*01:05-HLA-DPB1*25:01', 'HLA-DPA1*01:05-HLA-DPB1*26:01',
         'HLA-DPA1*01:05-HLA-DPB1*27:01', 'HLA-DPA1*01:05-HLA-DPB1*28:01', 'HLA-DPA1*01:05-HLA-DPB1*29:01', 'HLA-DPA1*01:05-HLA-DPB1*30:01',
         'HLA-DPA1*01:05-HLA-DPB1*31:01', 'HLA-DPA1*01:05-HLA-DPB1*32:01',
         'HLA-DPA1*01:05-HLA-DPB1*33:01', 'HLA-DPA1*01:05-HLA-DPB1*34:01', 'HLA-DPA1*01:05-HLA-DPB1*35:01', 'HLA-DPA1*01:05-HLA-DPB1*36:01',
         'HLA-DPA1*01:05-HLA-DPB1*37:01', 'HLA-DPA1*01:05-HLA-DPB1*38:01',
         'HLA-DPA1*01:05-HLA-DPB1*39:01', 'HLA-DPA1*01:05-HLA-DPB1*40:01', 'HLA-DPA1*01:05-HLA-DPB1*41:01', 'HLA-DPA1*01:05-HLA-DPB1*44:01',
         'HLA-DPA1*01:05-HLA-DPB1*45:01', 'HLA-DPA1*01:05-HLA-DPB1*46:01',
         'HLA-DPA1*01:05-HLA-DPB1*47:01', 'HLA-DPA1*01:05-HLA-DPB1*48:01', 'HLA-DPA1*01:05-HLA-DPB1*49:01', 'HLA-DPA1*01:05-HLA-DPB1*50:01',
         'HLA-DPA1*01:05-HLA-DPB1*51:01', 'HLA-DPA1*01:05-HLA-DPB1*52:01',
         'HLA-DPA1*01:05-HLA-DPB1*53:01', 'HLA-DPA1*01:05-HLA-DPB1*54:01', 'HLA-DPA1*01:05-HLA-DPB1*55:01', 'HLA-DPA1*01:05-HLA-DPB1*56:01',
         'HLA-DPA1*01:05-HLA-DPB1*58:01', 'HLA-DPA1*01:05-HLA-DPB1*59:01',
         'HLA-DPA1*01:05-HLA-DPB1*60:01', 'HLA-DPA1*01:05-HLA-DPB1*62:01', 'HLA-DPA1*01:05-HLA-DPB1*63:01', 'HLA-DPA1*01:05-HLA-DPB1*65:01',
         'HLA-DPA1*01:05-HLA-DPB1*66:01', 'HLA-DPA1*01:05-HLA-DPB1*67:01',
         'HLA-DPA1*01:05-HLA-DPB1*68:01', 'HLA-DPA1*01:05-HLA-DPB1*69:01', 'HLA-DPA1*01:05-HLA-DPB1*70:01', 'HLA-DPA1*01:05-HLA-DPB1*71:01',
         'HLA-DPA1*01:05-HLA-DPB1*72:01', 'HLA-DPA1*01:05-HLA-DPB1*73:01',
         'HLA-DPA1*01:05-HLA-DPB1*74:01', 'HLA-DPA1*01:05-HLA-DPB1*75:01', 'HLA-DPA1*01:05-HLA-DPB1*76:01', 'HLA-DPA1*01:05-HLA-DPB1*77:01',
         'HLA-DPA1*01:05-HLA-DPB1*78:01', 'HLA-DPA1*01:05-HLA-DPB1*79:01',
         'HLA-DPA1*01:05-HLA-DPB1*80:01', 'HLA-DPA1*01:05-HLA-DPB1*81:01', 'HLA-DPA1*01:05-HLA-DPB1*82:01', 'HLA-DPA1*01:05-HLA-DPB1*83:01',
         'HLA-DPA1*01:05-HLA-DPB1*84:01', 'HLA-DPA1*01:05-HLA-DPB1*85:01',
         'HLA-DPA1*01:05-HLA-DPB1*86:01', 'HLA-DPA1*01:05-HLA-DPB1*87:01', 'HLA-DPA1*01:05-HLA-DPB1*88:01', 'HLA-DPA1*01:05-HLA-DPB1*89:01',
         'HLA-DPA1*01:05-HLA-DPB1*90:01', 'HLA-DPA1*01:05-HLA-DPB1*91:01',
         'HLA-DPA1*01:05-HLA-DPB1*92:01', 'HLA-DPA1*01:05-HLA-DPB1*93:01', 'HLA-DPA1*01:05-HLA-DPB1*94:01', 'HLA-DPA1*01:05-HLA-DPB1*95:01',
         'HLA-DPA1*01:05-HLA-DPB1*96:01', 'HLA-DPA1*01:05-HLA-DPB1*97:01',
         'HLA-DPA1*01:05-HLA-DPB1*98:01', 'HLA-DPA1*01:05-HLA-DPB1*99:01', 'HLA-DPA1*01:06-HLA-DPB1*01:01', 'HLA-DPA1*01:06-HLA-DPB1*02:01',
         'HLA-DPA1*01:06-HLA-DPB1*02:02', 'HLA-DPA1*01:06-HLA-DPB1*03:01',
         'HLA-DPA1*01:06-HLA-DPB1*04:01', 'HLA-DPA1*01:06-HLA-DPB1*04:02', 'HLA-DPA1*01:06-HLA-DPB1*05:01', 'HLA-DPA1*01:06-HLA-DPB1*06:01',
         'HLA-DPA1*01:06-HLA-DPB1*08:01', 'HLA-DPA1*01:06-HLA-DPB1*09:01',
         'HLA-DPA1*01:06-HLA-DPB1*10:001', 'HLA-DPA1*01:06-HLA-DPB1*10:01', 'HLA-DPA1*01:06-HLA-DPB1*10:101', 'HLA-DPA1*01:06-HLA-DPB1*10:201',
         'HLA-DPA1*01:06-HLA-DPB1*10:301', 'HLA-DPA1*01:06-HLA-DPB1*10:401',
         'HLA-DPA1*01:06-HLA-DPB1*10:501', 'HLA-DPA1*01:06-HLA-DPB1*10:601', 'HLA-DPA1*01:06-HLA-DPB1*10:701', 'HLA-DPA1*01:06-HLA-DPB1*10:801',
         'HLA-DPA1*01:06-HLA-DPB1*10:901', 'HLA-DPA1*01:06-HLA-DPB1*11:001',
         'HLA-DPA1*01:06-HLA-DPB1*11:01', 'HLA-DPA1*01:06-HLA-DPB1*11:101', 'HLA-DPA1*01:06-HLA-DPB1*11:201', 'HLA-DPA1*01:06-HLA-DPB1*11:301',
         'HLA-DPA1*01:06-HLA-DPB1*11:401', 'HLA-DPA1*01:06-HLA-DPB1*11:501',
         'HLA-DPA1*01:06-HLA-DPB1*11:601', 'HLA-DPA1*01:06-HLA-DPB1*11:701', 'HLA-DPA1*01:06-HLA-DPB1*11:801', 'HLA-DPA1*01:06-HLA-DPB1*11:901',
         'HLA-DPA1*01:06-HLA-DPB1*12:101', 'HLA-DPA1*01:06-HLA-DPB1*12:201',
         'HLA-DPA1*01:06-HLA-DPB1*12:301', 'HLA-DPA1*01:06-HLA-DPB1*12:401', 'HLA-DPA1*01:06-HLA-DPB1*12:501', 'HLA-DPA1*01:06-HLA-DPB1*12:601',
         'HLA-DPA1*01:06-HLA-DPB1*12:701', 'HLA-DPA1*01:06-HLA-DPB1*12:801',
         'HLA-DPA1*01:06-HLA-DPB1*12:901', 'HLA-DPA1*01:06-HLA-DPB1*13:001', 'HLA-DPA1*01:06-HLA-DPB1*13:01', 'HLA-DPA1*01:06-HLA-DPB1*13:101',
         'HLA-DPA1*01:06-HLA-DPB1*13:201', 'HLA-DPA1*01:06-HLA-DPB1*13:301',
         'HLA-DPA1*01:06-HLA-DPB1*13:401', 'HLA-DPA1*01:06-HLA-DPB1*14:01', 'HLA-DPA1*01:06-HLA-DPB1*15:01', 'HLA-DPA1*01:06-HLA-DPB1*16:01',
         'HLA-DPA1*01:06-HLA-DPB1*17:01', 'HLA-DPA1*01:06-HLA-DPB1*18:01',
         'HLA-DPA1*01:06-HLA-DPB1*19:01', 'HLA-DPA1*01:06-HLA-DPB1*20:01', 'HLA-DPA1*01:06-HLA-DPB1*21:01', 'HLA-DPA1*01:06-HLA-DPB1*22:01',
         'HLA-DPA1*01:06-HLA-DPB1*23:01', 'HLA-DPA1*01:06-HLA-DPB1*24:01',
         'HLA-DPA1*01:06-HLA-DPB1*25:01', 'HLA-DPA1*01:06-HLA-DPB1*26:01', 'HLA-DPA1*01:06-HLA-DPB1*27:01', 'HLA-DPA1*01:06-HLA-DPB1*28:01',
         'HLA-DPA1*01:06-HLA-DPB1*29:01', 'HLA-DPA1*01:06-HLA-DPB1*30:01',
         'HLA-DPA1*01:06-HLA-DPB1*31:01', 'HLA-DPA1*01:06-HLA-DPB1*32:01', 'HLA-DPA1*01:06-HLA-DPB1*33:01', 'HLA-DPA1*01:06-HLA-DPB1*34:01',
         'HLA-DPA1*01:06-HLA-DPB1*35:01', 'HLA-DPA1*01:06-HLA-DPB1*36:01',
         'HLA-DPA1*01:06-HLA-DPB1*37:01', 'HLA-DPA1*01:06-HLA-DPB1*38:01', 'HLA-DPA1*01:06-HLA-DPB1*39:01', 'HLA-DPA1*01:06-HLA-DPB1*40:01',
         'HLA-DPA1*01:06-HLA-DPB1*41:01', 'HLA-DPA1*01:06-HLA-DPB1*44:01',
         'HLA-DPA1*01:06-HLA-DPB1*45:01', 'HLA-DPA1*01:06-HLA-DPB1*46:01', 'HLA-DPA1*01:06-HLA-DPB1*47:01', 'HLA-DPA1*01:06-HLA-DPB1*48:01',
         'HLA-DPA1*01:06-HLA-DPB1*49:01', 'HLA-DPA1*01:06-HLA-DPB1*50:01',
         'HLA-DPA1*01:06-HLA-DPB1*51:01', 'HLA-DPA1*01:06-HLA-DPB1*52:01', 'HLA-DPA1*01:06-HLA-DPB1*53:01', 'HLA-DPA1*01:06-HLA-DPB1*54:01',
         'HLA-DPA1*01:06-HLA-DPB1*55:01', 'HLA-DPA1*01:06-HLA-DPB1*56:01',
         'HLA-DPA1*01:06-HLA-DPB1*58:01', 'HLA-DPA1*01:06-HLA-DPB1*59:01', 'HLA-DPA1*01:06-HLA-DPB1*60:01', 'HLA-DPA1*01:06-HLA-DPB1*62:01',
         'HLA-DPA1*01:06-HLA-DPB1*63:01', 'HLA-DPA1*01:06-HLA-DPB1*65:01',
         'HLA-DPA1*01:06-HLA-DPB1*66:01', 'HLA-DPA1*01:06-HLA-DPB1*67:01', 'HLA-DPA1*01:06-HLA-DPB1*68:01', 'HLA-DPA1*01:06-HLA-DPB1*69:01',
         'HLA-DPA1*01:06-HLA-DPB1*70:01', 'HLA-DPA1*01:06-HLA-DPB1*71:01',
         'HLA-DPA1*01:06-HLA-DPB1*72:01', 'HLA-DPA1*01:06-HLA-DPB1*73:01', 'HLA-DPA1*01:06-HLA-DPB1*74:01', 'HLA-DPA1*01:06-HLA-DPB1*75:01',
         'HLA-DPA1*01:06-HLA-DPB1*76:01', 'HLA-DPA1*01:06-HLA-DPB1*77:01',
         'HLA-DPA1*01:06-HLA-DPB1*78:01', 'HLA-DPA1*01:06-HLA-DPB1*79:01', 'HLA-DPA1*01:06-HLA-DPB1*80:01', 'HLA-DPA1*01:06-HLA-DPB1*81:01',
         'HLA-DPA1*01:06-HLA-DPB1*82:01', 'HLA-DPA1*01:06-HLA-DPB1*83:01',
         'HLA-DPA1*01:06-HLA-DPB1*84:01', 'HLA-DPA1*01:06-HLA-DPB1*85:01', 'HLA-DPA1*01:06-HLA-DPB1*86:01', 'HLA-DPA1*01:06-HLA-DPB1*87:01',
         'HLA-DPA1*01:06-HLA-DPB1*88:01', 'HLA-DPA1*01:06-HLA-DPB1*89:01',
         'HLA-DPA1*01:06-HLA-DPB1*90:01', 'HLA-DPA1*01:06-HLA-DPB1*91:01', 'HLA-DPA1*01:06-HLA-DPB1*92:01', 'HLA-DPA1*01:06-HLA-DPB1*93:01',
         'HLA-DPA1*01:06-HLA-DPB1*94:01', 'HLA-DPA1*01:06-HLA-DPB1*95:01',
         'HLA-DPA1*01:06-HLA-DPB1*96:01', 'HLA-DPA1*01:06-HLA-DPB1*97:01', 'HLA-DPA1*01:06-HLA-DPB1*98:01', 'HLA-DPA1*01:06-HLA-DPB1*99:01',
         'HLA-DPA1*01:07-HLA-DPB1*01:01', 'HLA-DPA1*01:07-HLA-DPB1*02:01',
         'HLA-DPA1*01:07-HLA-DPB1*02:02', 'HLA-DPA1*01:07-HLA-DPB1*03:01', 'HLA-DPA1*01:07-HLA-DPB1*04:01', 'HLA-DPA1*01:07-HLA-DPB1*04:02',
         'HLA-DPA1*01:07-HLA-DPB1*05:01', 'HLA-DPA1*01:07-HLA-DPB1*06:01',
         'HLA-DPA1*01:07-HLA-DPB1*08:01', 'HLA-DPA1*01:07-HLA-DPB1*09:01', 'HLA-DPA1*01:07-HLA-DPB1*10:001', 'HLA-DPA1*01:07-HLA-DPB1*10:01',
         'HLA-DPA1*01:07-HLA-DPB1*10:101', 'HLA-DPA1*01:07-HLA-DPB1*10:201',
         'HLA-DPA1*01:07-HLA-DPB1*10:301', 'HLA-DPA1*01:07-HLA-DPB1*10:401', 'HLA-DPA1*01:07-HLA-DPB1*10:501', 'HLA-DPA1*01:07-HLA-DPB1*10:601',
         'HLA-DPA1*01:07-HLA-DPB1*10:701', 'HLA-DPA1*01:07-HLA-DPB1*10:801',
         'HLA-DPA1*01:07-HLA-DPB1*10:901', 'HLA-DPA1*01:07-HLA-DPB1*11:001', 'HLA-DPA1*01:07-HLA-DPB1*11:01', 'HLA-DPA1*01:07-HLA-DPB1*11:101',
         'HLA-DPA1*01:07-HLA-DPB1*11:201', 'HLA-DPA1*01:07-HLA-DPB1*11:301',
         'HLA-DPA1*01:07-HLA-DPB1*11:401', 'HLA-DPA1*01:07-HLA-DPB1*11:501', 'HLA-DPA1*01:07-HLA-DPB1*11:601', 'HLA-DPA1*01:07-HLA-DPB1*11:701',
         'HLA-DPA1*01:07-HLA-DPB1*11:801', 'HLA-DPA1*01:07-HLA-DPB1*11:901',
         'HLA-DPA1*01:07-HLA-DPB1*12:101', 'HLA-DPA1*01:07-HLA-DPB1*12:201', 'HLA-DPA1*01:07-HLA-DPB1*12:301', 'HLA-DPA1*01:07-HLA-DPB1*12:401',
         'HLA-DPA1*01:07-HLA-DPB1*12:501', 'HLA-DPA1*01:07-HLA-DPB1*12:601',
         'HLA-DPA1*01:07-HLA-DPB1*12:701', 'HLA-DPA1*01:07-HLA-DPB1*12:801', 'HLA-DPA1*01:07-HLA-DPB1*12:901', 'HLA-DPA1*01:07-HLA-DPB1*13:001',
         'HLA-DPA1*01:07-HLA-DPB1*13:01', 'HLA-DPA1*01:07-HLA-DPB1*13:101',
         'HLA-DPA1*01:07-HLA-DPB1*13:201', 'HLA-DPA1*01:07-HLA-DPB1*13:301', 'HLA-DPA1*01:07-HLA-DPB1*13:401', 'HLA-DPA1*01:07-HLA-DPB1*14:01',
         'HLA-DPA1*01:07-HLA-DPB1*15:01', 'HLA-DPA1*01:07-HLA-DPB1*16:01',
         'HLA-DPA1*01:07-HLA-DPB1*17:01', 'HLA-DPA1*01:07-HLA-DPB1*18:01', 'HLA-DPA1*01:07-HLA-DPB1*19:01', 'HLA-DPA1*01:07-HLA-DPB1*20:01',
         'HLA-DPA1*01:07-HLA-DPB1*21:01', 'HLA-DPA1*01:07-HLA-DPB1*22:01',
         'HLA-DPA1*01:07-HLA-DPB1*23:01', 'HLA-DPA1*01:07-HLA-DPB1*24:01', 'HLA-DPA1*01:07-HLA-DPB1*25:01', 'HLA-DPA1*01:07-HLA-DPB1*26:01',
         'HLA-DPA1*01:07-HLA-DPB1*27:01', 'HLA-DPA1*01:07-HLA-DPB1*28:01',
         'HLA-DPA1*01:07-HLA-DPB1*29:01', 'HLA-DPA1*01:07-HLA-DPB1*30:01', 'HLA-DPA1*01:07-HLA-DPB1*31:01', 'HLA-DPA1*01:07-HLA-DPB1*32:01',
         'HLA-DPA1*01:07-HLA-DPB1*33:01', 'HLA-DPA1*01:07-HLA-DPB1*34:01',
         'HLA-DPA1*01:07-HLA-DPB1*35:01', 'HLA-DPA1*01:07-HLA-DPB1*36:01', 'HLA-DPA1*01:07-HLA-DPB1*37:01', 'HLA-DPA1*01:07-HLA-DPB1*38:01',
         'HLA-DPA1*01:07-HLA-DPB1*39:01', 'HLA-DPA1*01:07-HLA-DPB1*40:01',
         'HLA-DPA1*01:07-HLA-DPB1*41:01', 'HLA-DPA1*01:07-HLA-DPB1*44:01', 'HLA-DPA1*01:07-HLA-DPB1*45:01', 'HLA-DPA1*01:07-HLA-DPB1*46:01',
         'HLA-DPA1*01:07-HLA-DPB1*47:01', 'HLA-DPA1*01:07-HLA-DPB1*48:01',
         'HLA-DPA1*01:07-HLA-DPB1*49:01', 'HLA-DPA1*01:07-HLA-DPB1*50:01', 'HLA-DPA1*01:07-HLA-DPB1*51:01', 'HLA-DPA1*01:07-HLA-DPB1*52:01',
         'HLA-DPA1*01:07-HLA-DPB1*53:01', 'HLA-DPA1*01:07-HLA-DPB1*54:01',
         'HLA-DPA1*01:07-HLA-DPB1*55:01', 'HLA-DPA1*01:07-HLA-DPB1*56:01', 'HLA-DPA1*01:07-HLA-DPB1*58:01', 'HLA-DPA1*01:07-HLA-DPB1*59:01',
         'HLA-DPA1*01:07-HLA-DPB1*60:01', 'HLA-DPA1*01:07-HLA-DPB1*62:01',
         'HLA-DPA1*01:07-HLA-DPB1*63:01', 'HLA-DPA1*01:07-HLA-DPB1*65:01', 'HLA-DPA1*01:07-HLA-DPB1*66:01', 'HLA-DPA1*01:07-HLA-DPB1*67:01',
         'HLA-DPA1*01:07-HLA-DPB1*68:01', 'HLA-DPA1*01:07-HLA-DPB1*69:01',
         'HLA-DPA1*01:07-HLA-DPB1*70:01', 'HLA-DPA1*01:07-HLA-DPB1*71:01', 'HLA-DPA1*01:07-HLA-DPB1*72:01', 'HLA-DPA1*01:07-HLA-DPB1*73:01',
         'HLA-DPA1*01:07-HLA-DPB1*74:01', 'HLA-DPA1*01:07-HLA-DPB1*75:01',
         'HLA-DPA1*01:07-HLA-DPB1*76:01', 'HLA-DPA1*01:07-HLA-DPB1*77:01', 'HLA-DPA1*01:07-HLA-DPB1*78:01', 'HLA-DPA1*01:07-HLA-DPB1*79:01',
         'HLA-DPA1*01:07-HLA-DPB1*80:01', 'HLA-DPA1*01:07-HLA-DPB1*81:01',
         'HLA-DPA1*01:07-HLA-DPB1*82:01', 'HLA-DPA1*01:07-HLA-DPB1*83:01', 'HLA-DPA1*01:07-HLA-DPB1*84:01', 'HLA-DPA1*01:07-HLA-DPB1*85:01',
         'HLA-DPA1*01:07-HLA-DPB1*86:01', 'HLA-DPA1*01:07-HLA-DPB1*87:01',
         'HLA-DPA1*01:07-HLA-DPB1*88:01', 'HLA-DPA1*01:07-HLA-DPB1*89:01', 'HLA-DPA1*01:07-HLA-DPB1*90:01', 'HLA-DPA1*01:07-HLA-DPB1*91:01',
         'HLA-DPA1*01:07-HLA-DPB1*92:01', 'HLA-DPA1*01:07-HLA-DPB1*93:01',
         'HLA-DPA1*01:07-HLA-DPB1*94:01', 'HLA-DPA1*01:07-HLA-DPB1*95:01', 'HLA-DPA1*01:07-HLA-DPB1*96:01', 'HLA-DPA1*01:07-HLA-DPB1*97:01',
         'HLA-DPA1*01:07-HLA-DPB1*98:01', 'HLA-DPA1*01:07-HLA-DPB1*99:01',
         'HLA-DPA1*01:08-HLA-DPB1*01:01', 'HLA-DPA1*01:08-HLA-DPB1*02:01', 'HLA-DPA1*01:08-HLA-DPB1*02:02', 'HLA-DPA1*01:08-HLA-DPB1*03:01',
         'HLA-DPA1*01:08-HLA-DPB1*04:01', 'HLA-DPA1*01:08-HLA-DPB1*04:02',
         'HLA-DPA1*01:08-HLA-DPB1*05:01', 'HLA-DPA1*01:08-HLA-DPB1*06:01', 'HLA-DPA1*01:08-HLA-DPB1*08:01', 'HLA-DPA1*01:08-HLA-DPB1*09:01',
         'HLA-DPA1*01:08-HLA-DPB1*10:001', 'HLA-DPA1*01:08-HLA-DPB1*10:01',
         'HLA-DPA1*01:08-HLA-DPB1*10:101', 'HLA-DPA1*01:08-HLA-DPB1*10:201', 'HLA-DPA1*01:08-HLA-DPB1*10:301', 'HLA-DPA1*01:08-HLA-DPB1*10:401',
         'HLA-DPA1*01:08-HLA-DPB1*10:501', 'HLA-DPA1*01:08-HLA-DPB1*10:601',
         'HLA-DPA1*01:08-HLA-DPB1*10:701', 'HLA-DPA1*01:08-HLA-DPB1*10:801', 'HLA-DPA1*01:08-HLA-DPB1*10:901', 'HLA-DPA1*01:08-HLA-DPB1*11:001',
         'HLA-DPA1*01:08-HLA-DPB1*11:01', 'HLA-DPA1*01:08-HLA-DPB1*11:101',
         'HLA-DPA1*01:08-HLA-DPB1*11:201', 'HLA-DPA1*01:08-HLA-DPB1*11:301', 'HLA-DPA1*01:08-HLA-DPB1*11:401', 'HLA-DPA1*01:08-HLA-DPB1*11:501',
         'HLA-DPA1*01:08-HLA-DPB1*11:601', 'HLA-DPA1*01:08-HLA-DPB1*11:701',
         'HLA-DPA1*01:08-HLA-DPB1*11:801', 'HLA-DPA1*01:08-HLA-DPB1*11:901', 'HLA-DPA1*01:08-HLA-DPB1*12:101', 'HLA-DPA1*01:08-HLA-DPB1*12:201',
         'HLA-DPA1*01:08-HLA-DPB1*12:301', 'HLA-DPA1*01:08-HLA-DPB1*12:401',
         'HLA-DPA1*01:08-HLA-DPB1*12:501', 'HLA-DPA1*01:08-HLA-DPB1*12:601', 'HLA-DPA1*01:08-HLA-DPB1*12:701', 'HLA-DPA1*01:08-HLA-DPB1*12:801',
         'HLA-DPA1*01:08-HLA-DPB1*12:901', 'HLA-DPA1*01:08-HLA-DPB1*13:001',
         'HLA-DPA1*01:08-HLA-DPB1*13:01', 'HLA-DPA1*01:08-HLA-DPB1*13:101', 'HLA-DPA1*01:08-HLA-DPB1*13:201', 'HLA-DPA1*01:08-HLA-DPB1*13:301',
         'HLA-DPA1*01:08-HLA-DPB1*13:401', 'HLA-DPA1*01:08-HLA-DPB1*14:01',
         'HLA-DPA1*01:08-HLA-DPB1*15:01', 'HLA-DPA1*01:08-HLA-DPB1*16:01', 'HLA-DPA1*01:08-HLA-DPB1*17:01', 'HLA-DPA1*01:08-HLA-DPB1*18:01',
         'HLA-DPA1*01:08-HLA-DPB1*19:01', 'HLA-DPA1*01:08-HLA-DPB1*20:01',
         'HLA-DPA1*01:08-HLA-DPB1*21:01', 'HLA-DPA1*01:08-HLA-DPB1*22:01', 'HLA-DPA1*01:08-HLA-DPB1*23:01', 'HLA-DPA1*01:08-HLA-DPB1*24:01',
         'HLA-DPA1*01:08-HLA-DPB1*25:01', 'HLA-DPA1*01:08-HLA-DPB1*26:01',
         'HLA-DPA1*01:08-HLA-DPB1*27:01', 'HLA-DPA1*01:08-HLA-DPB1*28:01', 'HLA-DPA1*01:08-HLA-DPB1*29:01', 'HLA-DPA1*01:08-HLA-DPB1*30:01',
         'HLA-DPA1*01:08-HLA-DPB1*31:01', 'HLA-DPA1*01:08-HLA-DPB1*32:01',
         'HLA-DPA1*01:08-HLA-DPB1*33:01', 'HLA-DPA1*01:08-HLA-DPB1*34:01', 'HLA-DPA1*01:08-HLA-DPB1*35:01', 'HLA-DPA1*01:08-HLA-DPB1*36:01',
         'HLA-DPA1*01:08-HLA-DPB1*37:01', 'HLA-DPA1*01:08-HLA-DPB1*38:01',
         'HLA-DPA1*01:08-HLA-DPB1*39:01', 'HLA-DPA1*01:08-HLA-DPB1*40:01', 'HLA-DPA1*01:08-HLA-DPB1*41:01', 'HLA-DPA1*01:08-HLA-DPB1*44:01',
         'HLA-DPA1*01:08-HLA-DPB1*45:01', 'HLA-DPA1*01:08-HLA-DPB1*46:01',
         'HLA-DPA1*01:08-HLA-DPB1*47:01', 'HLA-DPA1*01:08-HLA-DPB1*48:01', 'HLA-DPA1*01:08-HLA-DPB1*49:01', 'HLA-DPA1*01:08-HLA-DPB1*50:01',
         'HLA-DPA1*01:08-HLA-DPB1*51:01', 'HLA-DPA1*01:08-HLA-DPB1*52:01',
         'HLA-DPA1*01:08-HLA-DPB1*53:01', 'HLA-DPA1*01:08-HLA-DPB1*54:01', 'HLA-DPA1*01:08-HLA-DPB1*55:01', 'HLA-DPA1*01:08-HLA-DPB1*56:01',
         'HLA-DPA1*01:08-HLA-DPB1*58:01', 'HLA-DPA1*01:08-HLA-DPB1*59:01',
         'HLA-DPA1*01:08-HLA-DPB1*60:01', 'HLA-DPA1*01:08-HLA-DPB1*62:01', 'HLA-DPA1*01:08-HLA-DPB1*63:01', 'HLA-DPA1*01:08-HLA-DPB1*65:01',
         'HLA-DPA1*01:08-HLA-DPB1*66:01', 'HLA-DPA1*01:08-HLA-DPB1*67:01',
         'HLA-DPA1*01:08-HLA-DPB1*68:01', 'HLA-DPA1*01:08-HLA-DPB1*69:01', 'HLA-DPA1*01:08-HLA-DPB1*70:01', 'HLA-DPA1*01:08-HLA-DPB1*71:01',
         'HLA-DPA1*01:08-HLA-DPB1*72:01', 'HLA-DPA1*01:08-HLA-DPB1*73:01',
         'HLA-DPA1*01:08-HLA-DPB1*74:01', 'HLA-DPA1*01:08-HLA-DPB1*75:01', 'HLA-DPA1*01:08-HLA-DPB1*76:01', 'HLA-DPA1*01:08-HLA-DPB1*77:01',
         'HLA-DPA1*01:08-HLA-DPB1*78:01', 'HLA-DPA1*01:08-HLA-DPB1*79:01',
         'HLA-DPA1*01:08-HLA-DPB1*80:01', 'HLA-DPA1*01:08-HLA-DPB1*81:01', 'HLA-DPA1*01:08-HLA-DPB1*82:01', 'HLA-DPA1*01:08-HLA-DPB1*83:01',
         'HLA-DPA1*01:08-HLA-DPB1*84:01', 'HLA-DPA1*01:08-HLA-DPB1*85:01',
         'HLA-DPA1*01:08-HLA-DPB1*86:01', 'HLA-DPA1*01:08-HLA-DPB1*87:01', 'HLA-DPA1*01:08-HLA-DPB1*88:01', 'HLA-DPA1*01:08-HLA-DPB1*89:01',
         'HLA-DPA1*01:08-HLA-DPB1*90:01', 'HLA-DPA1*01:08-HLA-DPB1*91:01',
         'HLA-DPA1*01:08-HLA-DPB1*92:01', 'HLA-DPA1*01:08-HLA-DPB1*93:01', 'HLA-DPA1*01:08-HLA-DPB1*94:01', 'HLA-DPA1*01:08-HLA-DPB1*95:01',
         'HLA-DPA1*01:08-HLA-DPB1*96:01', 'HLA-DPA1*01:08-HLA-DPB1*97:01',
         'HLA-DPA1*01:08-HLA-DPB1*98:01', 'HLA-DPA1*01:08-HLA-DPB1*99:01', 'HLA-DPA1*01:09-HLA-DPB1*01:01', 'HLA-DPA1*01:09-HLA-DPB1*02:01',
         'HLA-DPA1*01:09-HLA-DPB1*02:02', 'HLA-DPA1*01:09-HLA-DPB1*03:01',
         'HLA-DPA1*01:09-HLA-DPB1*04:01', 'HLA-DPA1*01:09-HLA-DPB1*04:02', 'HLA-DPA1*01:09-HLA-DPB1*05:01', 'HLA-DPA1*01:09-HLA-DPB1*06:01',
         'HLA-DPA1*01:09-HLA-DPB1*08:01', 'HLA-DPA1*01:09-HLA-DPB1*09:01',
         'HLA-DPA1*01:09-HLA-DPB1*10:001', 'HLA-DPA1*01:09-HLA-DPB1*10:01', 'HLA-DPA1*01:09-HLA-DPB1*10:101', 'HLA-DPA1*01:09-HLA-DPB1*10:201',
         'HLA-DPA1*01:09-HLA-DPB1*10:301', 'HLA-DPA1*01:09-HLA-DPB1*10:401',
         'HLA-DPA1*01:09-HLA-DPB1*10:501', 'HLA-DPA1*01:09-HLA-DPB1*10:601', 'HLA-DPA1*01:09-HLA-DPB1*10:701', 'HLA-DPA1*01:09-HLA-DPB1*10:801',
         'HLA-DPA1*01:09-HLA-DPB1*10:901', 'HLA-DPA1*01:09-HLA-DPB1*11:001',
         'HLA-DPA1*01:09-HLA-DPB1*11:01', 'HLA-DPA1*01:09-HLA-DPB1*11:101', 'HLA-DPA1*01:09-HLA-DPB1*11:201', 'HLA-DPA1*01:09-HLA-DPB1*11:301',
         'HLA-DPA1*01:09-HLA-DPB1*11:401', 'HLA-DPA1*01:09-HLA-DPB1*11:501',
         'HLA-DPA1*01:09-HLA-DPB1*11:601', 'HLA-DPA1*01:09-HLA-DPB1*11:701', 'HLA-DPA1*01:09-HLA-DPB1*11:801', 'HLA-DPA1*01:09-HLA-DPB1*11:901',
         'HLA-DPA1*01:09-HLA-DPB1*12:101', 'HLA-DPA1*01:09-HLA-DPB1*12:201',
         'HLA-DPA1*01:09-HLA-DPB1*12:301', 'HLA-DPA1*01:09-HLA-DPB1*12:401', 'HLA-DPA1*01:09-HLA-DPB1*12:501', 'HLA-DPA1*01:09-HLA-DPB1*12:601',
         'HLA-DPA1*01:09-HLA-DPB1*12:701', 'HLA-DPA1*01:09-HLA-DPB1*12:801',
         'HLA-DPA1*01:09-HLA-DPB1*12:901', 'HLA-DPA1*01:09-HLA-DPB1*13:001', 'HLA-DPA1*01:09-HLA-DPB1*13:01', 'HLA-DPA1*01:09-HLA-DPB1*13:101',
         'HLA-DPA1*01:09-HLA-DPB1*13:201', 'HLA-DPA1*01:09-HLA-DPB1*13:301',
         'HLA-DPA1*01:09-HLA-DPB1*13:401', 'HLA-DPA1*01:09-HLA-DPB1*14:01', 'HLA-DPA1*01:09-HLA-DPB1*15:01', 'HLA-DPA1*01:09-HLA-DPB1*16:01',
         'HLA-DPA1*01:09-HLA-DPB1*17:01', 'HLA-DPA1*01:09-HLA-DPB1*18:01',
         'HLA-DPA1*01:09-HLA-DPB1*19:01', 'HLA-DPA1*01:09-HLA-DPB1*20:01', 'HLA-DPA1*01:09-HLA-DPB1*21:01', 'HLA-DPA1*01:09-HLA-DPB1*22:01',
         'HLA-DPA1*01:09-HLA-DPB1*23:01', 'HLA-DPA1*01:09-HLA-DPB1*24:01',
         'HLA-DPA1*01:09-HLA-DPB1*25:01', 'HLA-DPA1*01:09-HLA-DPB1*26:01', 'HLA-DPA1*01:09-HLA-DPB1*27:01', 'HLA-DPA1*01:09-HLA-DPB1*28:01',
         'HLA-DPA1*01:09-HLA-DPB1*29:01', 'HLA-DPA1*01:09-HLA-DPB1*30:01',
         'HLA-DPA1*01:09-HLA-DPB1*31:01', 'HLA-DPA1*01:09-HLA-DPB1*32:01', 'HLA-DPA1*01:09-HLA-DPB1*33:01', 'HLA-DPA1*01:09-HLA-DPB1*34:01',
         'HLA-DPA1*01:09-HLA-DPB1*35:01', 'HLA-DPA1*01:09-HLA-DPB1*36:01',
         'HLA-DPA1*01:09-HLA-DPB1*37:01', 'HLA-DPA1*01:09-HLA-DPB1*38:01', 'HLA-DPA1*01:09-HLA-DPB1*39:01', 'HLA-DPA1*01:09-HLA-DPB1*40:01',
         'HLA-DPA1*01:09-HLA-DPB1*41:01', 'HLA-DPA1*01:09-HLA-DPB1*44:01',
         'HLA-DPA1*01:09-HLA-DPB1*45:01', 'HLA-DPA1*01:09-HLA-DPB1*46:01', 'HLA-DPA1*01:09-HLA-DPB1*47:01', 'HLA-DPA1*01:09-HLA-DPB1*48:01',
         'HLA-DPA1*01:09-HLA-DPB1*49:01', 'HLA-DPA1*01:09-HLA-DPB1*50:01',
         'HLA-DPA1*01:09-HLA-DPB1*51:01', 'HLA-DPA1*01:09-HLA-DPB1*52:01', 'HLA-DPA1*01:09-HLA-DPB1*53:01', 'HLA-DPA1*01:09-HLA-DPB1*54:01',
         'HLA-DPA1*01:09-HLA-DPB1*55:01', 'HLA-DPA1*01:09-HLA-DPB1*56:01',
         'HLA-DPA1*01:09-HLA-DPB1*58:01', 'HLA-DPA1*01:09-HLA-DPB1*59:01', 'HLA-DPA1*01:09-HLA-DPB1*60:01', 'HLA-DPA1*01:09-HLA-DPB1*62:01',
         'HLA-DPA1*01:09-HLA-DPB1*63:01', 'HLA-DPA1*01:09-HLA-DPB1*65:01',
         'HLA-DPA1*01:09-HLA-DPB1*66:01', 'HLA-DPA1*01:09-HLA-DPB1*67:01', 'HLA-DPA1*01:09-HLA-DPB1*68:01', 'HLA-DPA1*01:09-HLA-DPB1*69:01',
         'HLA-DPA1*01:09-HLA-DPB1*70:01', 'HLA-DPA1*01:09-HLA-DPB1*71:01',
         'HLA-DPA1*01:09-HLA-DPB1*72:01', 'HLA-DPA1*01:09-HLA-DPB1*73:01', 'HLA-DPA1*01:09-HLA-DPB1*74:01', 'HLA-DPA1*01:09-HLA-DPB1*75:01',
         'HLA-DPA1*01:09-HLA-DPB1*76:01', 'HLA-DPA1*01:09-HLA-DPB1*77:01',
         'HLA-DPA1*01:09-HLA-DPB1*78:01', 'HLA-DPA1*01:09-HLA-DPB1*79:01', 'HLA-DPA1*01:09-HLA-DPB1*80:01', 'HLA-DPA1*01:09-HLA-DPB1*81:01',
         'HLA-DPA1*01:09-HLA-DPB1*82:01', 'HLA-DPA1*01:09-HLA-DPB1*83:01',
         'HLA-DPA1*01:09-HLA-DPB1*84:01', 'HLA-DPA1*01:09-HLA-DPB1*85:01', 'HLA-DPA1*01:09-HLA-DPB1*86:01', 'HLA-DPA1*01:09-HLA-DPB1*87:01',
         'HLA-DPA1*01:09-HLA-DPB1*88:01', 'HLA-DPA1*01:09-HLA-DPB1*89:01',
         'HLA-DPA1*01:09-HLA-DPB1*90:01', 'HLA-DPA1*01:09-HLA-DPB1*91:01', 'HLA-DPA1*01:09-HLA-DPB1*92:01', 'HLA-DPA1*01:09-HLA-DPB1*93:01',
         'HLA-DPA1*01:09-HLA-DPB1*94:01', 'HLA-DPA1*01:09-HLA-DPB1*95:01',
         'HLA-DPA1*01:09-HLA-DPB1*96:01', 'HLA-DPA1*01:09-HLA-DPB1*97:01', 'HLA-DPA1*01:09-HLA-DPB1*98:01', 'HLA-DPA1*01:09-HLA-DPB1*99:01',
         'HLA-DPA1*01:10-HLA-DPB1*01:01', 'HLA-DPA1*01:10-HLA-DPB1*02:01',
         'HLA-DPA1*01:10-HLA-DPB1*02:02', 'HLA-DPA1*01:10-HLA-DPB1*03:01', 'HLA-DPA1*01:10-HLA-DPB1*04:01', 'HLA-DPA1*01:10-HLA-DPB1*04:02',
         'HLA-DPA1*01:10-HLA-DPB1*05:01', 'HLA-DPA1*01:10-HLA-DPB1*06:01',
         'HLA-DPA1*01:10-HLA-DPB1*08:01', 'HLA-DPA1*01:10-HLA-DPB1*09:01', 'HLA-DPA1*01:10-HLA-DPB1*10:001', 'HLA-DPA1*01:10-HLA-DPB1*10:01',
         'HLA-DPA1*01:10-HLA-DPB1*10:101', 'HLA-DPA1*01:10-HLA-DPB1*10:201',
         'HLA-DPA1*01:10-HLA-DPB1*10:301', 'HLA-DPA1*01:10-HLA-DPB1*10:401', 'HLA-DPA1*01:10-HLA-DPB1*10:501', 'HLA-DPA1*01:10-HLA-DPB1*10:601',
         'HLA-DPA1*01:10-HLA-DPB1*10:701', 'HLA-DPA1*01:10-HLA-DPB1*10:801',
         'HLA-DPA1*01:10-HLA-DPB1*10:901', 'HLA-DPA1*01:10-HLA-DPB1*11:001', 'HLA-DPA1*01:10-HLA-DPB1*11:01', 'HLA-DPA1*01:10-HLA-DPB1*11:101',
         'HLA-DPA1*01:10-HLA-DPB1*11:201', 'HLA-DPA1*01:10-HLA-DPB1*11:301',
         'HLA-DPA1*01:10-HLA-DPB1*11:401', 'HLA-DPA1*01:10-HLA-DPB1*11:501', 'HLA-DPA1*01:10-HLA-DPB1*11:601', 'HLA-DPA1*01:10-HLA-DPB1*11:701',
         'HLA-DPA1*01:10-HLA-DPB1*11:801', 'HLA-DPA1*01:10-HLA-DPB1*11:901',
         'HLA-DPA1*01:10-HLA-DPB1*12:101', 'HLA-DPA1*01:10-HLA-DPB1*12:201', 'HLA-DPA1*01:10-HLA-DPB1*12:301', 'HLA-DPA1*01:10-HLA-DPB1*12:401',
         'HLA-DPA1*01:10-HLA-DPB1*12:501', 'HLA-DPA1*01:10-HLA-DPB1*12:601',
         'HLA-DPA1*01:10-HLA-DPB1*12:701', 'HLA-DPA1*01:10-HLA-DPB1*12:801', 'HLA-DPA1*01:10-HLA-DPB1*12:901', 'HLA-DPA1*01:10-HLA-DPB1*13:001',
         'HLA-DPA1*01:10-HLA-DPB1*13:01', 'HLA-DPA1*01:10-HLA-DPB1*13:101',
         'HLA-DPA1*01:10-HLA-DPB1*13:201', 'HLA-DPA1*01:10-HLA-DPB1*13:301', 'HLA-DPA1*01:10-HLA-DPB1*13:401', 'HLA-DPA1*01:10-HLA-DPB1*14:01',
         'HLA-DPA1*01:10-HLA-DPB1*15:01', 'HLA-DPA1*01:10-HLA-DPB1*16:01',
         'HLA-DPA1*01:10-HLA-DPB1*17:01', 'HLA-DPA1*01:10-HLA-DPB1*18:01', 'HLA-DPA1*01:10-HLA-DPB1*19:01', 'HLA-DPA1*01:10-HLA-DPB1*20:01',
         'HLA-DPA1*01:10-HLA-DPB1*21:01', 'HLA-DPA1*01:10-HLA-DPB1*22:01',
         'HLA-DPA1*01:10-HLA-DPB1*23:01', 'HLA-DPA1*01:10-HLA-DPB1*24:01', 'HLA-DPA1*01:10-HLA-DPB1*25:01', 'HLA-DPA1*01:10-HLA-DPB1*26:01',
         'HLA-DPA1*01:10-HLA-DPB1*27:01', 'HLA-DPA1*01:10-HLA-DPB1*28:01',
         'HLA-DPA1*01:10-HLA-DPB1*29:01', 'HLA-DPA1*01:10-HLA-DPB1*30:01', 'HLA-DPA1*01:10-HLA-DPB1*31:01', 'HLA-DPA1*01:10-HLA-DPB1*32:01',
         'HLA-DPA1*01:10-HLA-DPB1*33:01', 'HLA-DPA1*01:10-HLA-DPB1*34:01',
         'HLA-DPA1*01:10-HLA-DPB1*35:01', 'HLA-DPA1*01:10-HLA-DPB1*36:01', 'HLA-DPA1*01:10-HLA-DPB1*37:01', 'HLA-DPA1*01:10-HLA-DPB1*38:01',
         'HLA-DPA1*01:10-HLA-DPB1*39:01', 'HLA-DPA1*01:10-HLA-DPB1*40:01',
         'HLA-DPA1*01:10-HLA-DPB1*41:01', 'HLA-DPA1*01:10-HLA-DPB1*44:01', 'HLA-DPA1*01:10-HLA-DPB1*45:01', 'HLA-DPA1*01:10-HLA-DPB1*46:01',
         'HLA-DPA1*01:10-HLA-DPB1*47:01', 'HLA-DPA1*01:10-HLA-DPB1*48:01',
         'HLA-DPA1*01:10-HLA-DPB1*49:01', 'HLA-DPA1*01:10-HLA-DPB1*50:01', 'HLA-DPA1*01:10-HLA-DPB1*51:01', 'HLA-DPA1*01:10-HLA-DPB1*52:01',
         'HLA-DPA1*01:10-HLA-DPB1*53:01', 'HLA-DPA1*01:10-HLA-DPB1*54:01',
         'HLA-DPA1*01:10-HLA-DPB1*55:01', 'HLA-DPA1*01:10-HLA-DPB1*56:01', 'HLA-DPA1*01:10-HLA-DPB1*58:01', 'HLA-DPA1*01:10-HLA-DPB1*59:01',
         'HLA-DPA1*01:10-HLA-DPB1*60:01', 'HLA-DPA1*01:10-HLA-DPB1*62:01',
         'HLA-DPA1*01:10-HLA-DPB1*63:01', 'HLA-DPA1*01:10-HLA-DPB1*65:01', 'HLA-DPA1*01:10-HLA-DPB1*66:01', 'HLA-DPA1*01:10-HLA-DPB1*67:01',
         'HLA-DPA1*01:10-HLA-DPB1*68:01', 'HLA-DPA1*01:10-HLA-DPB1*69:01',
         'HLA-DPA1*01:10-HLA-DPB1*70:01', 'HLA-DPA1*01:10-HLA-DPB1*71:01', 'HLA-DPA1*01:10-HLA-DPB1*72:01', 'HLA-DPA1*01:10-HLA-DPB1*73:01',
         'HLA-DPA1*01:10-HLA-DPB1*74:01', 'HLA-DPA1*01:10-HLA-DPB1*75:01',
         'HLA-DPA1*01:10-HLA-DPB1*76:01', 'HLA-DPA1*01:10-HLA-DPB1*77:01', 'HLA-DPA1*01:10-HLA-DPB1*78:01', 'HLA-DPA1*01:10-HLA-DPB1*79:01',
         'HLA-DPA1*01:10-HLA-DPB1*80:01', 'HLA-DPA1*01:10-HLA-DPB1*81:01',
         'HLA-DPA1*01:10-HLA-DPB1*82:01', 'HLA-DPA1*01:10-HLA-DPB1*83:01', 'HLA-DPA1*01:10-HLA-DPB1*84:01', 'HLA-DPA1*01:10-HLA-DPB1*85:01',
         'HLA-DPA1*01:10-HLA-DPB1*86:01', 'HLA-DPA1*01:10-HLA-DPB1*87:01',
         'HLA-DPA1*01:10-HLA-DPB1*88:01', 'HLA-DPA1*01:10-HLA-DPB1*89:01', 'HLA-DPA1*01:10-HLA-DPB1*90:01', 'HLA-DPA1*01:10-HLA-DPB1*91:01',
         'HLA-DPA1*01:10-HLA-DPB1*92:01', 'HLA-DPA1*01:10-HLA-DPB1*93:01',
         'HLA-DPA1*01:10-HLA-DPB1*94:01', 'HLA-DPA1*01:10-HLA-DPB1*95:01', 'HLA-DPA1*01:10-HLA-DPB1*96:01', 'HLA-DPA1*01:10-HLA-DPB1*97:01',
         'HLA-DPA1*01:10-HLA-DPB1*98:01', 'HLA-DPA1*01:10-HLA-DPB1*99:01',
         'HLA-DPA1*02:01-HLA-DPB1*01:01', 'HLA-DPA1*02:01-HLA-DPB1*02:01', 'HLA-DPA1*02:01-HLA-DPB1*02:02', 'HLA-DPA1*02:01-HLA-DPB1*03:01',
         'HLA-DPA1*02:01-HLA-DPB1*04:01', 'HLA-DPA1*02:01-HLA-DPB1*04:02',
         'HLA-DPA1*02:01-HLA-DPB1*05:01', 'HLA-DPA1*02:01-HLA-DPB1*06:01', 'HLA-DPA1*02:01-HLA-DPB1*08:01', 'HLA-DPA1*02:01-HLA-DPB1*09:01',
         'HLA-DPA1*02:01-HLA-DPB1*10:001', 'HLA-DPA1*02:01-HLA-DPB1*10:01',
         'HLA-DPA1*02:01-HLA-DPB1*10:101', 'HLA-DPA1*02:01-HLA-DPB1*10:201', 'HLA-DPA1*02:01-HLA-DPB1*10:301', 'HLA-DPA1*02:01-HLA-DPB1*10:401',
         'HLA-DPA1*02:01-HLA-DPB1*10:501', 'HLA-DPA1*02:01-HLA-DPB1*10:601',
         'HLA-DPA1*02:01-HLA-DPB1*10:701', 'HLA-DPA1*02:01-HLA-DPB1*10:801', 'HLA-DPA1*02:01-HLA-DPB1*10:901', 'HLA-DPA1*02:01-HLA-DPB1*11:001',
         'HLA-DPA1*02:01-HLA-DPB1*11:01', 'HLA-DPA1*02:01-HLA-DPB1*11:101',
         'HLA-DPA1*02:01-HLA-DPB1*11:201', 'HLA-DPA1*02:01-HLA-DPB1*11:301', 'HLA-DPA1*02:01-HLA-DPB1*11:401', 'HLA-DPA1*02:01-HLA-DPB1*11:501',
         'HLA-DPA1*02:01-HLA-DPB1*11:601', 'HLA-DPA1*02:01-HLA-DPB1*11:701',
         'HLA-DPA1*02:01-HLA-DPB1*11:801', 'HLA-DPA1*02:01-HLA-DPB1*11:901', 'HLA-DPA1*02:01-HLA-DPB1*12:101', 'HLA-DPA1*02:01-HLA-DPB1*12:201',
         'HLA-DPA1*02:01-HLA-DPB1*12:301', 'HLA-DPA1*02:01-HLA-DPB1*12:401',
         'HLA-DPA1*02:01-HLA-DPB1*12:501', 'HLA-DPA1*02:01-HLA-DPB1*12:601', 'HLA-DPA1*02:01-HLA-DPB1*12:701', 'HLA-DPA1*02:01-HLA-DPB1*12:801',
         'HLA-DPA1*02:01-HLA-DPB1*12:901', 'HLA-DPA1*02:01-HLA-DPB1*13:001',
         'HLA-DPA1*02:01-HLA-DPB1*13:01', 'HLA-DPA1*02:01-HLA-DPB1*13:101', 'HLA-DPA1*02:01-HLA-DPB1*13:201', 'HLA-DPA1*02:01-HLA-DPB1*13:301',
         'HLA-DPA1*02:01-HLA-DPB1*13:401', 'HLA-DPA1*02:01-HLA-DPB1*14:01',
         'HLA-DPA1*02:01-HLA-DPB1*15:01', 'HLA-DPA1*02:01-HLA-DPB1*16:01', 'HLA-DPA1*02:01-HLA-DPB1*17:01', 'HLA-DPA1*02:01-HLA-DPB1*18:01',
         'HLA-DPA1*02:01-HLA-DPB1*19:01', 'HLA-DPA1*02:01-HLA-DPB1*20:01',
         'HLA-DPA1*02:01-HLA-DPB1*21:01', 'HLA-DPA1*02:01-HLA-DPB1*22:01', 'HLA-DPA1*02:01-HLA-DPB1*23:01', 'HLA-DPA1*02:01-HLA-DPB1*24:01',
         'HLA-DPA1*02:01-HLA-DPB1*25:01', 'HLA-DPA1*02:01-HLA-DPB1*26:01',
         'HLA-DPA1*02:01-HLA-DPB1*27:01', 'HLA-DPA1*02:01-HLA-DPB1*28:01', 'HLA-DPA1*02:01-HLA-DPB1*29:01', 'HLA-DPA1*02:01-HLA-DPB1*30:01',
         'HLA-DPA1*02:01-HLA-DPB1*31:01', 'HLA-DPA1*02:01-HLA-DPB1*32:01',
         'HLA-DPA1*02:01-HLA-DPB1*33:01', 'HLA-DPA1*02:01-HLA-DPB1*34:01', 'HLA-DPA1*02:01-HLA-DPB1*35:01', 'HLA-DPA1*02:01-HLA-DPB1*36:01',
         'HLA-DPA1*02:01-HLA-DPB1*37:01', 'HLA-DPA1*02:01-HLA-DPB1*38:01',
         'HLA-DPA1*02:01-HLA-DPB1*39:01', 'HLA-DPA1*02:01-HLA-DPB1*40:01', 'HLA-DPA1*02:01-HLA-DPB1*41:01', 'HLA-DPA1*02:01-HLA-DPB1*44:01',
         'HLA-DPA1*02:01-HLA-DPB1*45:01', 'HLA-DPA1*02:01-HLA-DPB1*46:01',
         'HLA-DPA1*02:01-HLA-DPB1*47:01', 'HLA-DPA1*02:01-HLA-DPB1*48:01', 'HLA-DPA1*02:01-HLA-DPB1*49:01', 'HLA-DPA1*02:01-HLA-DPB1*50:01',
         'HLA-DPA1*02:01-HLA-DPB1*51:01', 'HLA-DPA1*02:01-HLA-DPB1*52:01',
         'HLA-DPA1*02:01-HLA-DPB1*53:01', 'HLA-DPA1*02:01-HLA-DPB1*54:01', 'HLA-DPA1*02:01-HLA-DPB1*55:01', 'HLA-DPA1*02:01-HLA-DPB1*56:01',
         'HLA-DPA1*02:01-HLA-DPB1*58:01', 'HLA-DPA1*02:01-HLA-DPB1*59:01',
         'HLA-DPA1*02:01-HLA-DPB1*60:01', 'HLA-DPA1*02:01-HLA-DPB1*62:01', 'HLA-DPA1*02:01-HLA-DPB1*63:01', 'HLA-DPA1*02:01-HLA-DPB1*65:01',
         'HLA-DPA1*02:01-HLA-DPB1*66:01', 'HLA-DPA1*02:01-HLA-DPB1*67:01',
         'HLA-DPA1*02:01-HLA-DPB1*68:01', 'HLA-DPA1*02:01-HLA-DPB1*69:01', 'HLA-DPA1*02:01-HLA-DPB1*70:01', 'HLA-DPA1*02:01-HLA-DPB1*71:01',
         'HLA-DPA1*02:01-HLA-DPB1*72:01', 'HLA-DPA1*02:01-HLA-DPB1*73:01',
         'HLA-DPA1*02:01-HLA-DPB1*74:01', 'HLA-DPA1*02:01-HLA-DPB1*75:01', 'HLA-DPA1*02:01-HLA-DPB1*76:01', 'HLA-DPA1*02:01-HLA-DPB1*77:01',
         'HLA-DPA1*02:01-HLA-DPB1*78:01', 'HLA-DPA1*02:01-HLA-DPB1*79:01',
         'HLA-DPA1*02:01-HLA-DPB1*80:01', 'HLA-DPA1*02:01-HLA-DPB1*81:01', 'HLA-DPA1*02:01-HLA-DPB1*82:01', 'HLA-DPA1*02:01-HLA-DPB1*83:01',
         'HLA-DPA1*02:01-HLA-DPB1*84:01', 'HLA-DPA1*02:01-HLA-DPB1*85:01',
         'HLA-DPA1*02:01-HLA-DPB1*86:01', 'HLA-DPA1*02:01-HLA-DPB1*87:01', 'HLA-DPA1*02:01-HLA-DPB1*88:01', 'HLA-DPA1*02:01-HLA-DPB1*89:01',
         'HLA-DPA1*02:01-HLA-DPB1*90:01', 'HLA-DPA1*02:01-HLA-DPB1*91:01',
         'HLA-DPA1*02:01-HLA-DPB1*92:01', 'HLA-DPA1*02:01-HLA-DPB1*93:01', 'HLA-DPA1*02:01-HLA-DPB1*94:01', 'HLA-DPA1*02:01-HLA-DPB1*95:01',
         'HLA-DPA1*02:01-HLA-DPB1*96:01', 'HLA-DPA1*02:01-HLA-DPB1*97:01',
         'HLA-DPA1*02:01-HLA-DPB1*98:01', 'HLA-DPA1*02:01-HLA-DPB1*99:01', 'HLA-DPA1*02:02-HLA-DPB1*01:01', 'HLA-DPA1*02:02-HLA-DPB1*02:01',
         'HLA-DPA1*02:02-HLA-DPB1*02:02', 'HLA-DPA1*02:02-HLA-DPB1*03:01',
         'HLA-DPA1*02:02-HLA-DPB1*04:01', 'HLA-DPA1*02:02-HLA-DPB1*04:02', 'HLA-DPA1*02:02-HLA-DPB1*05:01', 'HLA-DPA1*02:02-HLA-DPB1*06:01',
         'HLA-DPA1*02:02-HLA-DPB1*08:01', 'HLA-DPA1*02:02-HLA-DPB1*09:01',
         'HLA-DPA1*02:02-HLA-DPB1*10:001', 'HLA-DPA1*02:02-HLA-DPB1*10:01', 'HLA-DPA1*02:02-HLA-DPB1*10:101', 'HLA-DPA1*02:02-HLA-DPB1*10:201',
         'HLA-DPA1*02:02-HLA-DPB1*10:301', 'HLA-DPA1*02:02-HLA-DPB1*10:401',
         'HLA-DPA1*02:02-HLA-DPB1*10:501', 'HLA-DPA1*02:02-HLA-DPB1*10:601', 'HLA-DPA1*02:02-HLA-DPB1*10:701', 'HLA-DPA1*02:02-HLA-DPB1*10:801',
         'HLA-DPA1*02:02-HLA-DPB1*10:901', 'HLA-DPA1*02:02-HLA-DPB1*11:001',
         'HLA-DPA1*02:02-HLA-DPB1*11:01', 'HLA-DPA1*02:02-HLA-DPB1*11:101', 'HLA-DPA1*02:02-HLA-DPB1*11:201', 'HLA-DPA1*02:02-HLA-DPB1*11:301',
         'HLA-DPA1*02:02-HLA-DPB1*11:401', 'HLA-DPA1*02:02-HLA-DPB1*11:501',
         'HLA-DPA1*02:02-HLA-DPB1*11:601', 'HLA-DPA1*02:02-HLA-DPB1*11:701', 'HLA-DPA1*02:02-HLA-DPB1*11:801', 'HLA-DPA1*02:02-HLA-DPB1*11:901',
         'HLA-DPA1*02:02-HLA-DPB1*12:101', 'HLA-DPA1*02:02-HLA-DPB1*12:201',
         'HLA-DPA1*02:02-HLA-DPB1*12:301', 'HLA-DPA1*02:02-HLA-DPB1*12:401', 'HLA-DPA1*02:02-HLA-DPB1*12:501', 'HLA-DPA1*02:02-HLA-DPB1*12:601',
         'HLA-DPA1*02:02-HLA-DPB1*12:701', 'HLA-DPA1*02:02-HLA-DPB1*12:801',
         'HLA-DPA1*02:02-HLA-DPB1*12:901', 'HLA-DPA1*02:02-HLA-DPB1*13:001', 'HLA-DPA1*02:02-HLA-DPB1*13:01', 'HLA-DPA1*02:02-HLA-DPB1*13:101',
         'HLA-DPA1*02:02-HLA-DPB1*13:201', 'HLA-DPA1*02:02-HLA-DPB1*13:301',
         'HLA-DPA1*02:02-HLA-DPB1*13:401', 'HLA-DPA1*02:02-HLA-DPB1*14:01', 'HLA-DPA1*02:02-HLA-DPB1*15:01', 'HLA-DPA1*02:02-HLA-DPB1*16:01',
         'HLA-DPA1*02:02-HLA-DPB1*17:01', 'HLA-DPA1*02:02-HLA-DPB1*18:01',
         'HLA-DPA1*02:02-HLA-DPB1*19:01', 'HLA-DPA1*02:02-HLA-DPB1*20:01', 'HLA-DPA1*02:02-HLA-DPB1*21:01', 'HLA-DPA1*02:02-HLA-DPB1*22:01',
         'HLA-DPA1*02:02-HLA-DPB1*23:01', 'HLA-DPA1*02:02-HLA-DPB1*24:01',
         'HLA-DPA1*02:02-HLA-DPB1*25:01', 'HLA-DPA1*02:02-HLA-DPB1*26:01', 'HLA-DPA1*02:02-HLA-DPB1*27:01', 'HLA-DPA1*02:02-HLA-DPB1*28:01',
         'HLA-DPA1*02:02-HLA-DPB1*29:01', 'HLA-DPA1*02:02-HLA-DPB1*30:01',
         'HLA-DPA1*02:02-HLA-DPB1*31:01', 'HLA-DPA1*02:02-HLA-DPB1*32:01', 'HLA-DPA1*02:02-HLA-DPB1*33:01', 'HLA-DPA1*02:02-HLA-DPB1*34:01',
         'HLA-DPA1*02:02-HLA-DPB1*35:01', 'HLA-DPA1*02:02-HLA-DPB1*36:01',
         'HLA-DPA1*02:02-HLA-DPB1*37:01', 'HLA-DPA1*02:02-HLA-DPB1*38:01', 'HLA-DPA1*02:02-HLA-DPB1*39:01', 'HLA-DPA1*02:02-HLA-DPB1*40:01',
         'HLA-DPA1*02:02-HLA-DPB1*41:01', 'HLA-DPA1*02:02-HLA-DPB1*44:01',
         'HLA-DPA1*02:02-HLA-DPB1*45:01', 'HLA-DPA1*02:02-HLA-DPB1*46:01', 'HLA-DPA1*02:02-HLA-DPB1*47:01', 'HLA-DPA1*02:02-HLA-DPB1*48:01',
         'HLA-DPA1*02:02-HLA-DPB1*49:01', 'HLA-DPA1*02:02-HLA-DPB1*50:01',
         'HLA-DPA1*02:02-HLA-DPB1*51:01', 'HLA-DPA1*02:02-HLA-DPB1*52:01', 'HLA-DPA1*02:02-HLA-DPB1*53:01', 'HLA-DPA1*02:02-HLA-DPB1*54:01',
         'HLA-DPA1*02:02-HLA-DPB1*55:01', 'HLA-DPA1*02:02-HLA-DPB1*56:01',
         'HLA-DPA1*02:02-HLA-DPB1*58:01', 'HLA-DPA1*02:02-HLA-DPB1*59:01', 'HLA-DPA1*02:02-HLA-DPB1*60:01', 'HLA-DPA1*02:02-HLA-DPB1*62:01',
         'HLA-DPA1*02:02-HLA-DPB1*63:01', 'HLA-DPA1*02:02-HLA-DPB1*65:01',
         'HLA-DPA1*02:02-HLA-DPB1*66:01', 'HLA-DPA1*02:02-HLA-DPB1*67:01', 'HLA-DPA1*02:02-HLA-DPB1*68:01', 'HLA-DPA1*02:02-HLA-DPB1*69:01',
         'HLA-DPA1*02:02-HLA-DPB1*70:01', 'HLA-DPA1*02:02-HLA-DPB1*71:01',
         'HLA-DPA1*02:02-HLA-DPB1*72:01', 'HLA-DPA1*02:02-HLA-DPB1*73:01', 'HLA-DPA1*02:02-HLA-DPB1*74:01', 'HLA-DPA1*02:02-HLA-DPB1*75:01',
         'HLA-DPA1*02:02-HLA-DPB1*76:01', 'HLA-DPA1*02:02-HLA-DPB1*77:01',
         'HLA-DPA1*02:02-HLA-DPB1*78:01', 'HLA-DPA1*02:02-HLA-DPB1*79:01', 'HLA-DPA1*02:02-HLA-DPB1*80:01', 'HLA-DPA1*02:02-HLA-DPB1*81:01',
         'HLA-DPA1*02:02-HLA-DPB1*82:01', 'HLA-DPA1*02:02-HLA-DPB1*83:01',
         'HLA-DPA1*02:02-HLA-DPB1*84:01', 'HLA-DPA1*02:02-HLA-DPB1*85:01', 'HLA-DPA1*02:02-HLA-DPB1*86:01', 'HLA-DPA1*02:02-HLA-DPB1*87:01',
         'HLA-DPA1*02:02-HLA-DPB1*88:01', 'HLA-DPA1*02:02-HLA-DPB1*89:01',
         'HLA-DPA1*02:02-HLA-DPB1*90:01', 'HLA-DPA1*02:02-HLA-DPB1*91:01', 'HLA-DPA1*02:02-HLA-DPB1*92:01', 'HLA-DPA1*02:02-HLA-DPB1*93:01',
         'HLA-DPA1*02:02-HLA-DPB1*94:01', 'HLA-DPA1*02:02-HLA-DPB1*95:01',
         'HLA-DPA1*02:02-HLA-DPB1*96:01', 'HLA-DPA1*02:02-HLA-DPB1*97:01', 'HLA-DPA1*02:02-HLA-DPB1*98:01', 'HLA-DPA1*02:02-HLA-DPB1*99:01',
         'HLA-DPA1*02:03-HLA-DPB1*01:01', 'HLA-DPA1*02:03-HLA-DPB1*02:01',
         'HLA-DPA1*02:03-HLA-DPB1*02:02', 'HLA-DPA1*02:03-HLA-DPB1*03:01', 'HLA-DPA1*02:03-HLA-DPB1*04:01', 'HLA-DPA1*02:03-HLA-DPB1*04:02',
         'HLA-DPA1*02:03-HLA-DPB1*05:01', 'HLA-DPA1*02:03-HLA-DPB1*06:01',
         'HLA-DPA1*02:03-HLA-DPB1*08:01', 'HLA-DPA1*02:03-HLA-DPB1*09:01', 'HLA-DPA1*02:03-HLA-DPB1*10:001', 'HLA-DPA1*02:03-HLA-DPB1*10:01',
         'HLA-DPA1*02:03-HLA-DPB1*10:101', 'HLA-DPA1*02:03-HLA-DPB1*10:201',
         'HLA-DPA1*02:03-HLA-DPB1*10:301', 'HLA-DPA1*02:03-HLA-DPB1*10:401', 'HLA-DPA1*02:03-HLA-DPB1*10:501', 'HLA-DPA1*02:03-HLA-DPB1*10:601',
         'HLA-DPA1*02:03-HLA-DPB1*10:701', 'HLA-DPA1*02:03-HLA-DPB1*10:801',
         'HLA-DPA1*02:03-HLA-DPB1*10:901', 'HLA-DPA1*02:03-HLA-DPB1*11:001', 'HLA-DPA1*02:03-HLA-DPB1*11:01', 'HLA-DPA1*02:03-HLA-DPB1*11:101',
         'HLA-DPA1*02:03-HLA-DPB1*11:201', 'HLA-DPA1*02:03-HLA-DPB1*11:301',
         'HLA-DPA1*02:03-HLA-DPB1*11:401', 'HLA-DPA1*02:03-HLA-DPB1*11:501', 'HLA-DPA1*02:03-HLA-DPB1*11:601', 'HLA-DPA1*02:03-HLA-DPB1*11:701',
         'HLA-DPA1*02:03-HLA-DPB1*11:801', 'HLA-DPA1*02:03-HLA-DPB1*11:901',
         'HLA-DPA1*02:03-HLA-DPB1*12:101', 'HLA-DPA1*02:03-HLA-DPB1*12:201', 'HLA-DPA1*02:03-HLA-DPB1*12:301', 'HLA-DPA1*02:03-HLA-DPB1*12:401',
         'HLA-DPA1*02:03-HLA-DPB1*12:501', 'HLA-DPA1*02:03-HLA-DPB1*12:601',
         'HLA-DPA1*02:03-HLA-DPB1*12:701', 'HLA-DPA1*02:03-HLA-DPB1*12:801', 'HLA-DPA1*02:03-HLA-DPB1*12:901', 'HLA-DPA1*02:03-HLA-DPB1*13:001',
         'HLA-DPA1*02:03-HLA-DPB1*13:01', 'HLA-DPA1*02:03-HLA-DPB1*13:101',
         'HLA-DPA1*02:03-HLA-DPB1*13:201', 'HLA-DPA1*02:03-HLA-DPB1*13:301', 'HLA-DPA1*02:03-HLA-DPB1*13:401', 'HLA-DPA1*02:03-HLA-DPB1*14:01',
         'HLA-DPA1*02:03-HLA-DPB1*15:01', 'HLA-DPA1*02:03-HLA-DPB1*16:01',
         'HLA-DPA1*02:03-HLA-DPB1*17:01', 'HLA-DPA1*02:03-HLA-DPB1*18:01', 'HLA-DPA1*02:03-HLA-DPB1*19:01', 'HLA-DPA1*02:03-HLA-DPB1*20:01',
         'HLA-DPA1*02:03-HLA-DPB1*21:01', 'HLA-DPA1*02:03-HLA-DPB1*22:01',
         'HLA-DPA1*02:03-HLA-DPB1*23:01', 'HLA-DPA1*02:03-HLA-DPB1*24:01', 'HLA-DPA1*02:03-HLA-DPB1*25:01', 'HLA-DPA1*02:03-HLA-DPB1*26:01',
         'HLA-DPA1*02:03-HLA-DPB1*27:01', 'HLA-DPA1*02:03-HLA-DPB1*28:01',
         'HLA-DPA1*02:03-HLA-DPB1*29:01', 'HLA-DPA1*02:03-HLA-DPB1*30:01', 'HLA-DPA1*02:03-HLA-DPB1*31:01', 'HLA-DPA1*02:03-HLA-DPB1*32:01',
         'HLA-DPA1*02:03-HLA-DPB1*33:01', 'HLA-DPA1*02:03-HLA-DPB1*34:01',
         'HLA-DPA1*02:03-HLA-DPB1*35:01', 'HLA-DPA1*02:03-HLA-DPB1*36:01', 'HLA-DPA1*02:03-HLA-DPB1*37:01', 'HLA-DPA1*02:03-HLA-DPB1*38:01',
         'HLA-DPA1*02:03-HLA-DPB1*39:01', 'HLA-DPA1*02:03-HLA-DPB1*40:01',
         'HLA-DPA1*02:03-HLA-DPB1*41:01', 'HLA-DPA1*02:03-HLA-DPB1*44:01', 'HLA-DPA1*02:03-HLA-DPB1*45:01', 'HLA-DPA1*02:03-HLA-DPB1*46:01',
         'HLA-DPA1*02:03-HLA-DPB1*47:01', 'HLA-DPA1*02:03-HLA-DPB1*48:01',
         'HLA-DPA1*02:03-HLA-DPB1*49:01', 'HLA-DPA1*02:03-HLA-DPB1*50:01', 'HLA-DPA1*02:03-HLA-DPB1*51:01', 'HLA-DPA1*02:03-HLA-DPB1*52:01',
         'HLA-DPA1*02:03-HLA-DPB1*53:01', 'HLA-DPA1*02:03-HLA-DPB1*54:01',
         'HLA-DPA1*02:03-HLA-DPB1*55:01', 'HLA-DPA1*02:03-HLA-DPB1*56:01', 'HLA-DPA1*02:03-HLA-DPB1*58:01', 'HLA-DPA1*02:03-HLA-DPB1*59:01',
         'HLA-DPA1*02:03-HLA-DPB1*60:01', 'HLA-DPA1*02:03-HLA-DPB1*62:01',
         'HLA-DPA1*02:03-HLA-DPB1*63:01', 'HLA-DPA1*02:03-HLA-DPB1*65:01', 'HLA-DPA1*02:03-HLA-DPB1*66:01', 'HLA-DPA1*02:03-HLA-DPB1*67:01',
         'HLA-DPA1*02:03-HLA-DPB1*68:01', 'HLA-DPA1*02:03-HLA-DPB1*69:01',
         'HLA-DPA1*02:03-HLA-DPB1*70:01', 'HLA-DPA1*02:03-HLA-DPB1*71:01', 'HLA-DPA1*02:03-HLA-DPB1*72:01', 'HLA-DPA1*02:03-HLA-DPB1*73:01',
         'HLA-DPA1*02:03-HLA-DPB1*74:01', 'HLA-DPA1*02:03-HLA-DPB1*75:01',
         'HLA-DPA1*02:03-HLA-DPB1*76:01', 'HLA-DPA1*02:03-HLA-DPB1*77:01', 'HLA-DPA1*02:03-HLA-DPB1*78:01', 'HLA-DPA1*02:03-HLA-DPB1*79:01',
         'HLA-DPA1*02:03-HLA-DPB1*80:01', 'HLA-DPA1*02:03-HLA-DPB1*81:01',
         'HLA-DPA1*02:03-HLA-DPB1*82:01', 'HLA-DPA1*02:03-HLA-DPB1*83:01', 'HLA-DPA1*02:03-HLA-DPB1*84:01', 'HLA-DPA1*02:03-HLA-DPB1*85:01',
         'HLA-DPA1*02:03-HLA-DPB1*86:01', 'HLA-DPA1*02:03-HLA-DPB1*87:01',
         'HLA-DPA1*02:03-HLA-DPB1*88:01', 'HLA-DPA1*02:03-HLA-DPB1*89:01', 'HLA-DPA1*02:03-HLA-DPB1*90:01', 'HLA-DPA1*02:03-HLA-DPB1*91:01',
         'HLA-DPA1*02:03-HLA-DPB1*92:01', 'HLA-DPA1*02:03-HLA-DPB1*93:01',
         'HLA-DPA1*02:03-HLA-DPB1*94:01', 'HLA-DPA1*02:03-HLA-DPB1*95:01', 'HLA-DPA1*02:03-HLA-DPB1*96:01', 'HLA-DPA1*02:03-HLA-DPB1*97:01',
         'HLA-DPA1*02:03-HLA-DPB1*98:01', 'HLA-DPA1*02:03-HLA-DPB1*99:01',
         'HLA-DPA1*02:04-HLA-DPB1*01:01', 'HLA-DPA1*02:04-HLA-DPB1*02:01', 'HLA-DPA1*02:04-HLA-DPB1*02:02', 'HLA-DPA1*02:04-HLA-DPB1*03:01',
         'HLA-DPA1*02:04-HLA-DPB1*04:01', 'HLA-DPA1*02:04-HLA-DPB1*04:02',
         'HLA-DPA1*02:04-HLA-DPB1*05:01', 'HLA-DPA1*02:04-HLA-DPB1*06:01', 'HLA-DPA1*02:04-HLA-DPB1*08:01', 'HLA-DPA1*02:04-HLA-DPB1*09:01',
         'HLA-DPA1*02:04-HLA-DPB1*10:001', 'HLA-DPA1*02:04-HLA-DPB1*10:01',
         'HLA-DPA1*02:04-HLA-DPB1*10:101', 'HLA-DPA1*02:04-HLA-DPB1*10:201', 'HLA-DPA1*02:04-HLA-DPB1*10:301', 'HLA-DPA1*02:04-HLA-DPB1*10:401',
         'HLA-DPA1*02:04-HLA-DPB1*10:501', 'HLA-DPA1*02:04-HLA-DPB1*10:601',
         'HLA-DPA1*02:04-HLA-DPB1*10:701', 'HLA-DPA1*02:04-HLA-DPB1*10:801', 'HLA-DPA1*02:04-HLA-DPB1*10:901', 'HLA-DPA1*02:04-HLA-DPB1*11:001',
         'HLA-DPA1*02:04-HLA-DPB1*11:01', 'HLA-DPA1*02:04-HLA-DPB1*11:101',
         'HLA-DPA1*02:04-HLA-DPB1*11:201', 'HLA-DPA1*02:04-HLA-DPB1*11:301', 'HLA-DPA1*02:04-HLA-DPB1*11:401', 'HLA-DPA1*02:04-HLA-DPB1*11:501',
         'HLA-DPA1*02:04-HLA-DPB1*11:601', 'HLA-DPA1*02:04-HLA-DPB1*11:701',
         'HLA-DPA1*02:04-HLA-DPB1*11:801', 'HLA-DPA1*02:04-HLA-DPB1*11:901', 'HLA-DPA1*02:04-HLA-DPB1*12:101', 'HLA-DPA1*02:04-HLA-DPB1*12:201',
         'HLA-DPA1*02:04-HLA-DPB1*12:301', 'HLA-DPA1*02:04-HLA-DPB1*12:401',
         'HLA-DPA1*02:04-HLA-DPB1*12:501', 'HLA-DPA1*02:04-HLA-DPB1*12:601', 'HLA-DPA1*02:04-HLA-DPB1*12:701', 'HLA-DPA1*02:04-HLA-DPB1*12:801',
         'HLA-DPA1*02:04-HLA-DPB1*12:901', 'HLA-DPA1*02:04-HLA-DPB1*13:001',
         'HLA-DPA1*02:04-HLA-DPB1*13:01', 'HLA-DPA1*02:04-HLA-DPB1*13:101', 'HLA-DPA1*02:04-HLA-DPB1*13:201', 'HLA-DPA1*02:04-HLA-DPB1*13:301',
         'HLA-DPA1*02:04-HLA-DPB1*13:401', 'HLA-DPA1*02:04-HLA-DPB1*14:01',
         'HLA-DPA1*02:04-HLA-DPB1*15:01', 'HLA-DPA1*02:04-HLA-DPB1*16:01', 'HLA-DPA1*02:04-HLA-DPB1*17:01', 'HLA-DPA1*02:04-HLA-DPB1*18:01',
         'HLA-DPA1*02:04-HLA-DPB1*19:01', 'HLA-DPA1*02:04-HLA-DPB1*20:01',
         'HLA-DPA1*02:04-HLA-DPB1*21:01', 'HLA-DPA1*02:04-HLA-DPB1*22:01', 'HLA-DPA1*02:04-HLA-DPB1*23:01', 'HLA-DPA1*02:04-HLA-DPB1*24:01',
         'HLA-DPA1*02:04-HLA-DPB1*25:01', 'HLA-DPA1*02:04-HLA-DPB1*26:01',
         'HLA-DPA1*02:04-HLA-DPB1*27:01', 'HLA-DPA1*02:04-HLA-DPB1*28:01', 'HLA-DPA1*02:04-HLA-DPB1*29:01', 'HLA-DPA1*02:04-HLA-DPB1*30:01',
         'HLA-DPA1*02:04-HLA-DPB1*31:01', 'HLA-DPA1*02:04-HLA-DPB1*32:01',
         'HLA-DPA1*02:04-HLA-DPB1*33:01', 'HLA-DPA1*02:04-HLA-DPB1*34:01', 'HLA-DPA1*02:04-HLA-DPB1*35:01', 'HLA-DPA1*02:04-HLA-DPB1*36:01',
         'HLA-DPA1*02:04-HLA-DPB1*37:01', 'HLA-DPA1*02:04-HLA-DPB1*38:01',
         'HLA-DPA1*02:04-HLA-DPB1*39:01', 'HLA-DPA1*02:04-HLA-DPB1*40:01', 'HLA-DPA1*02:04-HLA-DPB1*41:01', 'HLA-DPA1*02:04-HLA-DPB1*44:01',
         'HLA-DPA1*02:04-HLA-DPB1*45:01', 'HLA-DPA1*02:04-HLA-DPB1*46:01',
         'HLA-DPA1*02:04-HLA-DPB1*47:01', 'HLA-DPA1*02:04-HLA-DPB1*48:01', 'HLA-DPA1*02:04-HLA-DPB1*49:01', 'HLA-DPA1*02:04-HLA-DPB1*50:01',
         'HLA-DPA1*02:04-HLA-DPB1*51:01', 'HLA-DPA1*02:04-HLA-DPB1*52:01',
         'HLA-DPA1*02:04-HLA-DPB1*53:01', 'HLA-DPA1*02:04-HLA-DPB1*54:01', 'HLA-DPA1*02:04-HLA-DPB1*55:01', 'HLA-DPA1*02:04-HLA-DPB1*56:01',
         'HLA-DPA1*02:04-HLA-DPB1*58:01', 'HLA-DPA1*02:04-HLA-DPB1*59:01',
         'HLA-DPA1*02:04-HLA-DPB1*60:01', 'HLA-DPA1*02:04-HLA-DPB1*62:01', 'HLA-DPA1*02:04-HLA-DPB1*63:01', 'HLA-DPA1*02:04-HLA-DPB1*65:01',
         'HLA-DPA1*02:04-HLA-DPB1*66:01', 'HLA-DPA1*02:04-HLA-DPB1*67:01',
         'HLA-DPA1*02:04-HLA-DPB1*68:01', 'HLA-DPA1*02:04-HLA-DPB1*69:01', 'HLA-DPA1*02:04-HLA-DPB1*70:01', 'HLA-DPA1*02:04-HLA-DPB1*71:01',
         'HLA-DPA1*02:04-HLA-DPB1*72:01', 'HLA-DPA1*02:04-HLA-DPB1*73:01',
         'HLA-DPA1*02:04-HLA-DPB1*74:01', 'HLA-DPA1*02:04-HLA-DPB1*75:01', 'HLA-DPA1*02:04-HLA-DPB1*76:01', 'HLA-DPA1*02:04-HLA-DPB1*77:01',
         'HLA-DPA1*02:04-HLA-DPB1*78:01', 'HLA-DPA1*02:04-HLA-DPB1*79:01',
         'HLA-DPA1*02:04-HLA-DPB1*80:01', 'HLA-DPA1*02:04-HLA-DPB1*81:01', 'HLA-DPA1*02:04-HLA-DPB1*82:01', 'HLA-DPA1*02:04-HLA-DPB1*83:01',
         'HLA-DPA1*02:04-HLA-DPB1*84:01', 'HLA-DPA1*02:04-HLA-DPB1*85:01',
         'HLA-DPA1*02:04-HLA-DPB1*86:01', 'HLA-DPA1*02:04-HLA-DPB1*87:01', 'HLA-DPA1*02:04-HLA-DPB1*88:01', 'HLA-DPA1*02:04-HLA-DPB1*89:01',
         'HLA-DPA1*02:04-HLA-DPB1*90:01', 'HLA-DPA1*02:04-HLA-DPB1*91:01',
         'HLA-DPA1*02:04-HLA-DPB1*92:01', 'HLA-DPA1*02:04-HLA-DPB1*93:01', 'HLA-DPA1*02:04-HLA-DPB1*94:01', 'HLA-DPA1*02:04-HLA-DPB1*95:01',
         'HLA-DPA1*02:04-HLA-DPB1*96:01', 'HLA-DPA1*02:04-HLA-DPB1*97:01',
         'HLA-DPA1*02:04-HLA-DPB1*98:01', 'HLA-DPA1*02:04-HLA-DPB1*99:01', 'HLA-DPA1*03:01-HLA-DPB1*01:01', 'HLA-DPA1*03:01-HLA-DPB1*02:01',
         'HLA-DPA1*03:01-HLA-DPB1*02:02', 'HLA-DPA1*03:01-HLA-DPB1*03:01',
         'HLA-DPA1*03:01-HLA-DPB1*04:01', 'HLA-DPA1*03:01-HLA-DPB1*04:02', 'HLA-DPA1*03:01-HLA-DPB1*05:01', 'HLA-DPA1*03:01-HLA-DPB1*06:01',
         'HLA-DPA1*03:01-HLA-DPB1*08:01', 'HLA-DPA1*03:01-HLA-DPB1*09:01',
         'HLA-DPA1*03:01-HLA-DPB1*10:001', 'HLA-DPA1*03:01-HLA-DPB1*10:01', 'HLA-DPA1*03:01-HLA-DPB1*10:101', 'HLA-DPA1*03:01-HLA-DPB1*10:201',
         'HLA-DPA1*03:01-HLA-DPB1*10:301', 'HLA-DPA1*03:01-HLA-DPB1*10:401',
         'HLA-DPA1*03:01-HLA-DPB1*10:501', 'HLA-DPA1*03:01-HLA-DPB1*10:601', 'HLA-DPA1*03:01-HLA-DPB1*10:701', 'HLA-DPA1*03:01-HLA-DPB1*10:801',
         'HLA-DPA1*03:01-HLA-DPB1*10:901', 'HLA-DPA1*03:01-HLA-DPB1*11:001',
         'HLA-DPA1*03:01-HLA-DPB1*11:01', 'HLA-DPA1*03:01-HLA-DPB1*11:101', 'HLA-DPA1*03:01-HLA-DPB1*11:201', 'HLA-DPA1*03:01-HLA-DPB1*11:301',
         'HLA-DPA1*03:01-HLA-DPB1*11:401', 'HLA-DPA1*03:01-HLA-DPB1*11:501',
         'HLA-DPA1*03:01-HLA-DPB1*11:601', 'HLA-DPA1*03:01-HLA-DPB1*11:701', 'HLA-DPA1*03:01-HLA-DPB1*11:801', 'HLA-DPA1*03:01-HLA-DPB1*11:901',
         'HLA-DPA1*03:01-HLA-DPB1*12:101', 'HLA-DPA1*03:01-HLA-DPB1*12:201',
         'HLA-DPA1*03:01-HLA-DPB1*12:301', 'HLA-DPA1*03:01-HLA-DPB1*12:401', 'HLA-DPA1*03:01-HLA-DPB1*12:501', 'HLA-DPA1*03:01-HLA-DPB1*12:601',
         'HLA-DPA1*03:01-HLA-DPB1*12:701', 'HLA-DPA1*03:01-HLA-DPB1*12:801',
         'HLA-DPA1*03:01-HLA-DPB1*12:901', 'HLA-DPA1*03:01-HLA-DPB1*13:001', 'HLA-DPA1*03:01-HLA-DPB1*13:01', 'HLA-DPA1*03:01-HLA-DPB1*13:101',
         'HLA-DPA1*03:01-HLA-DPB1*13:201', 'HLA-DPA1*03:01-HLA-DPB1*13:301',
         'HLA-DPA1*03:01-HLA-DPB1*13:401', 'HLA-DPA1*03:01-HLA-DPB1*14:01', 'HLA-DPA1*03:01-HLA-DPB1*15:01', 'HLA-DPA1*03:01-HLA-DPB1*16:01',
         'HLA-DPA1*03:01-HLA-DPB1*17:01', 'HLA-DPA1*03:01-HLA-DPB1*18:01',
         'HLA-DPA1*03:01-HLA-DPB1*19:01', 'HLA-DPA1*03:01-HLA-DPB1*20:01', 'HLA-DPA1*03:01-HLA-DPB1*21:01', 'HLA-DPA1*03:01-HLA-DPB1*22:01',
         'HLA-DPA1*03:01-HLA-DPB1*23:01', 'HLA-DPA1*03:01-HLA-DPB1*24:01',
         'HLA-DPA1*03:01-HLA-DPB1*25:01', 'HLA-DPA1*03:01-HLA-DPB1*26:01', 'HLA-DPA1*03:01-HLA-DPB1*27:01', 'HLA-DPA1*03:01-HLA-DPB1*28:01',
         'HLA-DPA1*03:01-HLA-DPB1*29:01', 'HLA-DPA1*03:01-HLA-DPB1*30:01',
         'HLA-DPA1*03:01-HLA-DPB1*31:01', 'HLA-DPA1*03:01-HLA-DPB1*32:01', 'HLA-DPA1*03:01-HLA-DPB1*33:01', 'HLA-DPA1*03:01-HLA-DPB1*34:01',
         'HLA-DPA1*03:01-HLA-DPB1*35:01', 'HLA-DPA1*03:01-HLA-DPB1*36:01',
         'HLA-DPA1*03:01-HLA-DPB1*37:01', 'HLA-DPA1*03:01-HLA-DPB1*38:01', 'HLA-DPA1*03:01-HLA-DPB1*39:01', 'HLA-DPA1*03:01-HLA-DPB1*40:01',
         'HLA-DPA1*03:01-HLA-DPB1*41:01', 'HLA-DPA1*03:01-HLA-DPB1*44:01',
         'HLA-DPA1*03:01-HLA-DPB1*45:01', 'HLA-DPA1*03:01-HLA-DPB1*46:01', 'HLA-DPA1*03:01-HLA-DPB1*47:01', 'HLA-DPA1*03:01-HLA-DPB1*48:01',
         'HLA-DPA1*03:01-HLA-DPB1*49:01', 'HLA-DPA1*03:01-HLA-DPB1*50:01',
         'HLA-DPA1*03:01-HLA-DPB1*51:01', 'HLA-DPA1*03:01-HLA-DPB1*52:01', 'HLA-DPA1*03:01-HLA-DPB1*53:01', 'HLA-DPA1*03:01-HLA-DPB1*54:01',
         'HLA-DPA1*03:01-HLA-DPB1*55:01', 'HLA-DPA1*03:01-HLA-DPB1*56:01',
         'HLA-DPA1*03:01-HLA-DPB1*58:01', 'HLA-DPA1*03:01-HLA-DPB1*59:01', 'HLA-DPA1*03:01-HLA-DPB1*60:01', 'HLA-DPA1*03:01-HLA-DPB1*62:01',
         'HLA-DPA1*03:01-HLA-DPB1*63:01', 'HLA-DPA1*03:01-HLA-DPB1*65:01',
         'HLA-DPA1*03:01-HLA-DPB1*66:01', 'HLA-DPA1*03:01-HLA-DPB1*67:01', 'HLA-DPA1*03:01-HLA-DPB1*68:01', 'HLA-DPA1*03:01-HLA-DPB1*69:01',
         'HLA-DPA1*03:01-HLA-DPB1*70:01', 'HLA-DPA1*03:01-HLA-DPB1*71:01',
         'HLA-DPA1*03:01-HLA-DPB1*72:01', 'HLA-DPA1*03:01-HLA-DPB1*73:01', 'HLA-DPA1*03:01-HLA-DPB1*74:01', 'HLA-DPA1*03:01-HLA-DPB1*75:01',
         'HLA-DPA1*03:01-HLA-DPB1*76:01', 'HLA-DPA1*03:01-HLA-DPB1*77:01',
         'HLA-DPA1*03:01-HLA-DPB1*78:01', 'HLA-DPA1*03:01-HLA-DPB1*79:01', 'HLA-DPA1*03:01-HLA-DPB1*80:01', 'HLA-DPA1*03:01-HLA-DPB1*81:01',
         'HLA-DPA1*03:01-HLA-DPB1*82:01', 'HLA-DPA1*03:01-HLA-DPB1*83:01',
         'HLA-DPA1*03:01-HLA-DPB1*84:01', 'HLA-DPA1*03:01-HLA-DPB1*85:01', 'HLA-DPA1*03:01-HLA-DPB1*86:01', 'HLA-DPA1*03:01-HLA-DPB1*87:01',
         'HLA-DPA1*03:01-HLA-DPB1*88:01', 'HLA-DPA1*03:01-HLA-DPB1*89:01',
         'HLA-DPA1*03:01-HLA-DPB1*90:01', 'HLA-DPA1*03:01-HLA-DPB1*91:01', 'HLA-DPA1*03:01-HLA-DPB1*92:01', 'HLA-DPA1*03:01-HLA-DPB1*93:01',
         'HLA-DPA1*03:01-HLA-DPB1*94:01', 'HLA-DPA1*03:01-HLA-DPB1*95:01',
         'HLA-DPA1*03:01-HLA-DPB1*96:01', 'HLA-DPA1*03:01-HLA-DPB1*97:01', 'HLA-DPA1*03:01-HLA-DPB1*98:01', 'HLA-DPA1*03:01-HLA-DPB1*99:01',
         'HLA-DPA1*03:02-HLA-DPB1*01:01', 'HLA-DPA1*03:02-HLA-DPB1*02:01',
         'HLA-DPA1*03:02-HLA-DPB1*02:02', 'HLA-DPA1*03:02-HLA-DPB1*03:01', 'HLA-DPA1*03:02-HLA-DPB1*04:01', 'HLA-DPA1*03:02-HLA-DPB1*04:02',
         'HLA-DPA1*03:02-HLA-DPB1*05:01', 'HLA-DPA1*03:02-HLA-DPB1*06:01',
         'HLA-DPA1*03:02-HLA-DPB1*08:01', 'HLA-DPA1*03:02-HLA-DPB1*09:01', 'HLA-DPA1*03:02-HLA-DPB1*10:001', 'HLA-DPA1*03:02-HLA-DPB1*10:01',
         'HLA-DPA1*03:02-HLA-DPB1*10:101', 'HLA-DPA1*03:02-HLA-DPB1*10:201',
         'HLA-DPA1*03:02-HLA-DPB1*10:301', 'HLA-DPA1*03:02-HLA-DPB1*10:401', 'HLA-DPA1*03:02-HLA-DPB1*10:501', 'HLA-DPA1*03:02-HLA-DPB1*10:601',
         'HLA-DPA1*03:02-HLA-DPB1*10:701', 'HLA-DPA1*03:02-HLA-DPB1*10:801',
         'HLA-DPA1*03:02-HLA-DPB1*10:901', 'HLA-DPA1*03:02-HLA-DPB1*11:001', 'HLA-DPA1*03:02-HLA-DPB1*11:01', 'HLA-DPA1*03:02-HLA-DPB1*11:101',
         'HLA-DPA1*03:02-HLA-DPB1*11:201', 'HLA-DPA1*03:02-HLA-DPB1*11:301',
         'HLA-DPA1*03:02-HLA-DPB1*11:401', 'HLA-DPA1*03:02-HLA-DPB1*11:501', 'HLA-DPA1*03:02-HLA-DPB1*11:601', 'HLA-DPA1*03:02-HLA-DPB1*11:701',
         'HLA-DPA1*03:02-HLA-DPB1*11:801', 'HLA-DPA1*03:02-HLA-DPB1*11:901',
         'HLA-DPA1*03:02-HLA-DPB1*12:101', 'HLA-DPA1*03:02-HLA-DPB1*12:201', 'HLA-DPA1*03:02-HLA-DPB1*12:301', 'HLA-DPA1*03:02-HLA-DPB1*12:401',
         'HLA-DPA1*03:02-HLA-DPB1*12:501', 'HLA-DPA1*03:02-HLA-DPB1*12:601',
         'HLA-DPA1*03:02-HLA-DPB1*12:701', 'HLA-DPA1*03:02-HLA-DPB1*12:801', 'HLA-DPA1*03:02-HLA-DPB1*12:901', 'HLA-DPA1*03:02-HLA-DPB1*13:001',
         'HLA-DPA1*03:02-HLA-DPB1*13:01', 'HLA-DPA1*03:02-HLA-DPB1*13:101',
         'HLA-DPA1*03:02-HLA-DPB1*13:201', 'HLA-DPA1*03:02-HLA-DPB1*13:301', 'HLA-DPA1*03:02-HLA-DPB1*13:401', 'HLA-DPA1*03:02-HLA-DPB1*14:01',
         'HLA-DPA1*03:02-HLA-DPB1*15:01', 'HLA-DPA1*03:02-HLA-DPB1*16:01',
         'HLA-DPA1*03:02-HLA-DPB1*17:01', 'HLA-DPA1*03:02-HLA-DPB1*18:01', 'HLA-DPA1*03:02-HLA-DPB1*19:01', 'HLA-DPA1*03:02-HLA-DPB1*20:01',
         'HLA-DPA1*03:02-HLA-DPB1*21:01', 'HLA-DPA1*03:02-HLA-DPB1*22:01',
         'HLA-DPA1*03:02-HLA-DPB1*23:01', 'HLA-DPA1*03:02-HLA-DPB1*24:01', 'HLA-DPA1*03:02-HLA-DPB1*25:01', 'HLA-DPA1*03:02-HLA-DPB1*26:01',
         'HLA-DPA1*03:02-HLA-DPB1*27:01', 'HLA-DPA1*03:02-HLA-DPB1*28:01',
         'HLA-DPA1*03:02-HLA-DPB1*29:01', 'HLA-DPA1*03:02-HLA-DPB1*30:01', 'HLA-DPA1*03:02-HLA-DPB1*31:01', 'HLA-DPA1*03:02-HLA-DPB1*32:01',
         'HLA-DPA1*03:02-HLA-DPB1*33:01', 'HLA-DPA1*03:02-HLA-DPB1*34:01',
         'HLA-DPA1*03:02-HLA-DPB1*35:01', 'HLA-DPA1*03:02-HLA-DPB1*36:01', 'HLA-DPA1*03:02-HLA-DPB1*37:01', 'HLA-DPA1*03:02-HLA-DPB1*38:01',
         'HLA-DPA1*03:02-HLA-DPB1*39:01', 'HLA-DPA1*03:02-HLA-DPB1*40:01',
         'HLA-DPA1*03:02-HLA-DPB1*41:01', 'HLA-DPA1*03:02-HLA-DPB1*44:01', 'HLA-DPA1*03:02-HLA-DPB1*45:01', 'HLA-DPA1*03:02-HLA-DPB1*46:01',
         'HLA-DPA1*03:02-HLA-DPB1*47:01', 'HLA-DPA1*03:02-HLA-DPB1*48:01',
         'HLA-DPA1*03:02-HLA-DPB1*49:01', 'HLA-DPA1*03:02-HLA-DPB1*50:01', 'HLA-DPA1*03:02-HLA-DPB1*51:01', 'HLA-DPA1*03:02-HLA-DPB1*52:01',
         'HLA-DPA1*03:02-HLA-DPB1*53:01', 'HLA-DPA1*03:02-HLA-DPB1*54:01',
         'HLA-DPA1*03:02-HLA-DPB1*55:01', 'HLA-DPA1*03:02-HLA-DPB1*56:01', 'HLA-DPA1*03:02-HLA-DPB1*58:01', 'HLA-DPA1*03:02-HLA-DPB1*59:01',
         'HLA-DPA1*03:02-HLA-DPB1*60:01', 'HLA-DPA1*03:02-HLA-DPB1*62:01',
         'HLA-DPA1*03:02-HLA-DPB1*63:01', 'HLA-DPA1*03:02-HLA-DPB1*65:01', 'HLA-DPA1*03:02-HLA-DPB1*66:01', 'HLA-DPA1*03:02-HLA-DPB1*67:01',
         'HLA-DPA1*03:02-HLA-DPB1*68:01', 'HLA-DPA1*03:02-HLA-DPB1*69:01',
         'HLA-DPA1*03:02-HLA-DPB1*70:01', 'HLA-DPA1*03:02-HLA-DPB1*71:01', 'HLA-DPA1*03:02-HLA-DPB1*72:01', 'HLA-DPA1*03:02-HLA-DPB1*73:01',
         'HLA-DPA1*03:02-HLA-DPB1*74:01', 'HLA-DPA1*03:02-HLA-DPB1*75:01',
         'HLA-DPA1*03:02-HLA-DPB1*76:01', 'HLA-DPA1*03:02-HLA-DPB1*77:01', 'HLA-DPA1*03:02-HLA-DPB1*78:01', 'HLA-DPA1*03:02-HLA-DPB1*79:01',
         'HLA-DPA1*03:02-HLA-DPB1*80:01', 'HLA-DPA1*03:02-HLA-DPB1*81:01',
         'HLA-DPA1*03:02-HLA-DPB1*82:01', 'HLA-DPA1*03:02-HLA-DPB1*83:01', 'HLA-DPA1*03:02-HLA-DPB1*84:01', 'HLA-DPA1*03:02-HLA-DPB1*85:01',
         'HLA-DPA1*03:02-HLA-DPB1*86:01', 'HLA-DPA1*03:02-HLA-DPB1*87:01',
         'HLA-DPA1*03:02-HLA-DPB1*88:01', 'HLA-DPA1*03:02-HLA-DPB1*89:01', 'HLA-DPA1*03:02-HLA-DPB1*90:01', 'HLA-DPA1*03:02-HLA-DPB1*91:01',
         'HLA-DPA1*03:02-HLA-DPB1*92:01', 'HLA-DPA1*03:02-HLA-DPB1*93:01',
         'HLA-DPA1*03:02-HLA-DPB1*94:01', 'HLA-DPA1*03:02-HLA-DPB1*95:01', 'HLA-DPA1*03:02-HLA-DPB1*96:01', 'HLA-DPA1*03:02-HLA-DPB1*97:01',
         'HLA-DPA1*03:02-HLA-DPB1*98:01', 'HLA-DPA1*03:02-HLA-DPB1*99:01',
         'HLA-DPA1*03:03-HLA-DPB1*01:01', 'HLA-DPA1*03:03-HLA-DPB1*02:01', 'HLA-DPA1*03:03-HLA-DPB1*02:02', 'HLA-DPA1*03:03-HLA-DPB1*03:01',
         'HLA-DPA1*03:03-HLA-DPB1*04:01', 'HLA-DPA1*03:03-HLA-DPB1*04:02',
         'HLA-DPA1*03:03-HLA-DPB1*05:01', 'HLA-DPA1*03:03-HLA-DPB1*06:01', 'HLA-DPA1*03:03-HLA-DPB1*08:01', 'HLA-DPA1*03:03-HLA-DPB1*09:01',
         'HLA-DPA1*03:03-HLA-DPB1*10:001', 'HLA-DPA1*03:03-HLA-DPB1*10:01',
         'HLA-DPA1*03:03-HLA-DPB1*10:101', 'HLA-DPA1*03:03-HLA-DPB1*10:201', 'HLA-DPA1*03:03-HLA-DPB1*10:301', 'HLA-DPA1*03:03-HLA-DPB1*10:401',
         'HLA-DPA1*03:03-HLA-DPB1*10:501', 'HLA-DPA1*03:03-HLA-DPB1*10:601',
         'HLA-DPA1*03:03-HLA-DPB1*10:701', 'HLA-DPA1*03:03-HLA-DPB1*10:801', 'HLA-DPA1*03:03-HLA-DPB1*10:901', 'HLA-DPA1*03:03-HLA-DPB1*11:001',
         'HLA-DPA1*03:03-HLA-DPB1*11:01', 'HLA-DPA1*03:03-HLA-DPB1*11:101',
         'HLA-DPA1*03:03-HLA-DPB1*11:201', 'HLA-DPA1*03:03-HLA-DPB1*11:301', 'HLA-DPA1*03:03-HLA-DPB1*11:401', 'HLA-DPA1*03:03-HLA-DPB1*11:501',
         'HLA-DPA1*03:03-HLA-DPB1*11:601', 'HLA-DPA1*03:03-HLA-DPB1*11:701',
         'HLA-DPA1*03:03-HLA-DPB1*11:801', 'HLA-DPA1*03:03-HLA-DPB1*11:901', 'HLA-DPA1*03:03-HLA-DPB1*12:101', 'HLA-DPA1*03:03-HLA-DPB1*12:201',
         'HLA-DPA1*03:03-HLA-DPB1*12:301', 'HLA-DPA1*03:03-HLA-DPB1*12:401',
         'HLA-DPA1*03:03-HLA-DPB1*12:501', 'HLA-DPA1*03:03-HLA-DPB1*12:601', 'HLA-DPA1*03:03-HLA-DPB1*12:701', 'HLA-DPA1*03:03-HLA-DPB1*12:801',
         'HLA-DPA1*03:03-HLA-DPB1*12:901', 'HLA-DPA1*03:03-HLA-DPB1*13:001',
         'HLA-DPA1*03:03-HLA-DPB1*13:01', 'HLA-DPA1*03:03-HLA-DPB1*13:101', 'HLA-DPA1*03:03-HLA-DPB1*13:201', 'HLA-DPA1*03:03-HLA-DPB1*13:301',
         'HLA-DPA1*03:03-HLA-DPB1*13:401', 'HLA-DPA1*03:03-HLA-DPB1*14:01',
         'HLA-DPA1*03:03-HLA-DPB1*15:01', 'HLA-DPA1*03:03-HLA-DPB1*16:01', 'HLA-DPA1*03:03-HLA-DPB1*17:01', 'HLA-DPA1*03:03-HLA-DPB1*18:01',
         'HLA-DPA1*03:03-HLA-DPB1*19:01', 'HLA-DPA1*03:03-HLA-DPB1*20:01',
         'HLA-DPA1*03:03-HLA-DPB1*21:01', 'HLA-DPA1*03:03-HLA-DPB1*22:01', 'HLA-DPA1*03:03-HLA-DPB1*23:01', 'HLA-DPA1*03:03-HLA-DPB1*24:01',
         'HLA-DPA1*03:03-HLA-DPB1*25:01', 'HLA-DPA1*03:03-HLA-DPB1*26:01',
         'HLA-DPA1*03:03-HLA-DPB1*27:01', 'HLA-DPA1*03:03-HLA-DPB1*28:01', 'HLA-DPA1*03:03-HLA-DPB1*29:01', 'HLA-DPA1*03:03-HLA-DPB1*30:01',
         'HLA-DPA1*03:03-HLA-DPB1*31:01', 'HLA-DPA1*03:03-HLA-DPB1*32:01',
         'HLA-DPA1*03:03-HLA-DPB1*33:01', 'HLA-DPA1*03:03-HLA-DPB1*34:01', 'HLA-DPA1*03:03-HLA-DPB1*35:01', 'HLA-DPA1*03:03-HLA-DPB1*36:01',
         'HLA-DPA1*03:03-HLA-DPB1*37:01', 'HLA-DPA1*03:03-HLA-DPB1*38:01',
         'HLA-DPA1*03:03-HLA-DPB1*39:01', 'HLA-DPA1*03:03-HLA-DPB1*40:01', 'HLA-DPA1*03:03-HLA-DPB1*41:01', 'HLA-DPA1*03:03-HLA-DPB1*44:01',
         'HLA-DPA1*03:03-HLA-DPB1*45:01', 'HLA-DPA1*03:03-HLA-DPB1*46:01',
         'HLA-DPA1*03:03-HLA-DPB1*47:01', 'HLA-DPA1*03:03-HLA-DPB1*48:01', 'HLA-DPA1*03:03-HLA-DPB1*49:01', 'HLA-DPA1*03:03-HLA-DPB1*50:01',
         'HLA-DPA1*03:03-HLA-DPB1*51:01', 'HLA-DPA1*03:03-HLA-DPB1*52:01',
         'HLA-DPA1*03:03-HLA-DPB1*53:01', 'HLA-DPA1*03:03-HLA-DPB1*54:01', 'HLA-DPA1*03:03-HLA-DPB1*55:01', 'HLA-DPA1*03:03-HLA-DPB1*56:01',
         'HLA-DPA1*03:03-HLA-DPB1*58:01', 'HLA-DPA1*03:03-HLA-DPB1*59:01',
         'HLA-DPA1*03:03-HLA-DPB1*60:01', 'HLA-DPA1*03:03-HLA-DPB1*62:01', 'HLA-DPA1*03:03-HLA-DPB1*63:01', 'HLA-DPA1*03:03-HLA-DPB1*65:01',
         'HLA-DPA1*03:03-HLA-DPB1*66:01', 'HLA-DPA1*03:03-HLA-DPB1*67:01',
         'HLA-DPA1*03:03-HLA-DPB1*68:01', 'HLA-DPA1*03:03-HLA-DPB1*69:01', 'HLA-DPA1*03:03-HLA-DPB1*70:01', 'HLA-DPA1*03:03-HLA-DPB1*71:01',
         'HLA-DPA1*03:03-HLA-DPB1*72:01', 'HLA-DPA1*03:03-HLA-DPB1*73:01',
         'HLA-DPA1*03:03-HLA-DPB1*74:01', 'HLA-DPA1*03:03-HLA-DPB1*75:01', 'HLA-DPA1*03:03-HLA-DPB1*76:01', 'HLA-DPA1*03:03-HLA-DPB1*77:01',
         'HLA-DPA1*03:03-HLA-DPB1*78:01', 'HLA-DPA1*03:03-HLA-DPB1*79:01',
         'HLA-DPA1*03:03-HLA-DPB1*80:01', 'HLA-DPA1*03:03-HLA-DPB1*81:01', 'HLA-DPA1*03:03-HLA-DPB1*82:01', 'HLA-DPA1*03:03-HLA-DPB1*83:01',
         'HLA-DPA1*03:03-HLA-DPB1*84:01', 'HLA-DPA1*03:03-HLA-DPB1*85:01',
         'HLA-DPA1*03:03-HLA-DPB1*86:01', 'HLA-DPA1*03:03-HLA-DPB1*87:01', 'HLA-DPA1*03:03-HLA-DPB1*88:01', 'HLA-DPA1*03:03-HLA-DPB1*89:01',
         'HLA-DPA1*03:03-HLA-DPB1*90:01', 'HLA-DPA1*03:03-HLA-DPB1*91:01',
         'HLA-DPA1*03:03-HLA-DPB1*92:01', 'HLA-DPA1*03:03-HLA-DPB1*93:01', 'HLA-DPA1*03:03-HLA-DPB1*94:01', 'HLA-DPA1*03:03-HLA-DPB1*95:01',
         'HLA-DPA1*03:03-HLA-DPB1*96:01', 'HLA-DPA1*03:03-HLA-DPB1*97:01',
         'HLA-DPA1*03:03-HLA-DPB1*98:01', 'HLA-DPA1*03:03-HLA-DPB1*99:01', 'HLA-DPA1*04:01-HLA-DPB1*01:01', 'HLA-DPA1*04:01-HLA-DPB1*02:01',
         'HLA-DPA1*04:01-HLA-DPB1*02:02', 'HLA-DPA1*04:01-HLA-DPB1*03:01',
         'HLA-DPA1*04:01-HLA-DPB1*04:01', 'HLA-DPA1*04:01-HLA-DPB1*04:02', 'HLA-DPA1*04:01-HLA-DPB1*05:01', 'HLA-DPA1*04:01-HLA-DPB1*06:01',
         'HLA-DPA1*04:01-HLA-DPB1*08:01', 'HLA-DPA1*04:01-HLA-DPB1*09:01',
         'HLA-DPA1*04:01-HLA-DPB1*10:001', 'HLA-DPA1*04:01-HLA-DPB1*10:01', 'HLA-DPA1*04:01-HLA-DPB1*10:101', 'HLA-DPA1*04:01-HLA-DPB1*10:201',
         'HLA-DPA1*04:01-HLA-DPB1*10:301', 'HLA-DPA1*04:01-HLA-DPB1*10:401',
         'HLA-DPA1*04:01-HLA-DPB1*10:501', 'HLA-DPA1*04:01-HLA-DPB1*10:601', 'HLA-DPA1*04:01-HLA-DPB1*10:701', 'HLA-DPA1*04:01-HLA-DPB1*10:801',
         'HLA-DPA1*04:01-HLA-DPB1*10:901', 'HLA-DPA1*04:01-HLA-DPB1*11:001',
         'HLA-DPA1*04:01-HLA-DPB1*11:01', 'HLA-DPA1*04:01-HLA-DPB1*11:101', 'HLA-DPA1*04:01-HLA-DPB1*11:201', 'HLA-DPA1*04:01-HLA-DPB1*11:301',
         'HLA-DPA1*04:01-HLA-DPB1*11:401', 'HLA-DPA1*04:01-HLA-DPB1*11:501',
         'HLA-DPA1*04:01-HLA-DPB1*11:601', 'HLA-DPA1*04:01-HLA-DPB1*11:701', 'HLA-DPA1*04:01-HLA-DPB1*11:801', 'HLA-DPA1*04:01-HLA-DPB1*11:901',
         'HLA-DPA1*04:01-HLA-DPB1*12:101', 'HLA-DPA1*04:01-HLA-DPB1*12:201',
         'HLA-DPA1*04:01-HLA-DPB1*12:301', 'HLA-DPA1*04:01-HLA-DPB1*12:401', 'HLA-DPA1*04:01-HLA-DPB1*12:501', 'HLA-DPA1*04:01-HLA-DPB1*12:601',
         'HLA-DPA1*04:01-HLA-DPB1*12:701', 'HLA-DPA1*04:01-HLA-DPB1*12:801',
         'HLA-DPA1*04:01-HLA-DPB1*12:901', 'HLA-DPA1*04:01-HLA-DPB1*13:001', 'HLA-DPA1*04:01-HLA-DPB1*13:01', 'HLA-DPA1*04:01-HLA-DPB1*13:101',
         'HLA-DPA1*04:01-HLA-DPB1*13:201', 'HLA-DPA1*04:01-HLA-DPB1*13:301',
         'HLA-DPA1*04:01-HLA-DPB1*13:401', 'HLA-DPA1*04:01-HLA-DPB1*14:01', 'HLA-DPA1*04:01-HLA-DPB1*15:01', 'HLA-DPA1*04:01-HLA-DPB1*16:01',
         'HLA-DPA1*04:01-HLA-DPB1*17:01', 'HLA-DPA1*04:01-HLA-DPB1*18:01',
         'HLA-DPA1*04:01-HLA-DPB1*19:01', 'HLA-DPA1*04:01-HLA-DPB1*20:01', 'HLA-DPA1*04:01-HLA-DPB1*21:01', 'HLA-DPA1*04:01-HLA-DPB1*22:01',
         'HLA-DPA1*04:01-HLA-DPB1*23:01', 'HLA-DPA1*04:01-HLA-DPB1*24:01',
         'HLA-DPA1*04:01-HLA-DPB1*25:01', 'HLA-DPA1*04:01-HLA-DPB1*26:01', 'HLA-DPA1*04:01-HLA-DPB1*27:01', 'HLA-DPA1*04:01-HLA-DPB1*28:01',
         'HLA-DPA1*04:01-HLA-DPB1*29:01', 'HLA-DPA1*04:01-HLA-DPB1*30:01',
         'HLA-DPA1*04:01-HLA-DPB1*31:01', 'HLA-DPA1*04:01-HLA-DPB1*32:01', 'HLA-DPA1*04:01-HLA-DPB1*33:01', 'HLA-DPA1*04:01-HLA-DPB1*34:01',
         'HLA-DPA1*04:01-HLA-DPB1*35:01', 'HLA-DPA1*04:01-HLA-DPB1*36:01',
         'HLA-DPA1*04:01-HLA-DPB1*37:01', 'HLA-DPA1*04:01-HLA-DPB1*38:01', 'HLA-DPA1*04:01-HLA-DPB1*39:01', 'HLA-DPA1*04:01-HLA-DPB1*40:01',
         'HLA-DPA1*04:01-HLA-DPB1*41:01', 'HLA-DPA1*04:01-HLA-DPB1*44:01',
         'HLA-DPA1*04:01-HLA-DPB1*45:01', 'HLA-DPA1*04:01-HLA-DPB1*46:01', 'HLA-DPA1*04:01-HLA-DPB1*47:01', 'HLA-DPA1*04:01-HLA-DPB1*48:01',
         'HLA-DPA1*04:01-HLA-DPB1*49:01', 'HLA-DPA1*04:01-HLA-DPB1*50:01',
         'HLA-DPA1*04:01-HLA-DPB1*51:01', 'HLA-DPA1*04:01-HLA-DPB1*52:01', 'HLA-DPA1*04:01-HLA-DPB1*53:01', 'HLA-DPA1*04:01-HLA-DPB1*54:01',
         'HLA-DPA1*04:01-HLA-DPB1*55:01', 'HLA-DPA1*04:01-HLA-DPB1*56:01',
         'HLA-DPA1*04:01-HLA-DPB1*58:01', 'HLA-DPA1*04:01-HLA-DPB1*59:01', 'HLA-DPA1*04:01-HLA-DPB1*60:01', 'HLA-DPA1*04:01-HLA-DPB1*62:01',
         'HLA-DPA1*04:01-HLA-DPB1*63:01', 'HLA-DPA1*04:01-HLA-DPB1*65:01',
         'HLA-DPA1*04:01-HLA-DPB1*66:01', 'HLA-DPA1*04:01-HLA-DPB1*67:01', 'HLA-DPA1*04:01-HLA-DPB1*68:01', 'HLA-DPA1*04:01-HLA-DPB1*69:01',
         'HLA-DPA1*04:01-HLA-DPB1*70:01', 'HLA-DPA1*04:01-HLA-DPB1*71:01',
         'HLA-DPA1*04:01-HLA-DPB1*72:01', 'HLA-DPA1*04:01-HLA-DPB1*73:01', 'HLA-DPA1*04:01-HLA-DPB1*74:01', 'HLA-DPA1*04:01-HLA-DPB1*75:01',
         'HLA-DPA1*04:01-HLA-DPB1*76:01', 'HLA-DPA1*04:01-HLA-DPB1*77:01',
         'HLA-DPA1*04:01-HLA-DPB1*78:01', 'HLA-DPA1*04:01-HLA-DPB1*79:01', 'HLA-DPA1*04:01-HLA-DPB1*80:01', 'HLA-DPA1*04:01-HLA-DPB1*81:01',
         'HLA-DPA1*04:01-HLA-DPB1*82:01', 'HLA-DPA1*04:01-HLA-DPB1*83:01',
         'HLA-DPA1*04:01-HLA-DPB1*84:01', 'HLA-DPA1*04:01-HLA-DPB1*85:01', 'HLA-DPA1*04:01-HLA-DPB1*86:01', 'HLA-DPA1*04:01-HLA-DPB1*87:01',
         'HLA-DPA1*04:01-HLA-DPB1*88:01', 'HLA-DPA1*04:01-HLA-DPB1*89:01',
         'HLA-DPA1*04:01-HLA-DPB1*90:01', 'HLA-DPA1*04:01-HLA-DPB1*91:01', 'HLA-DPA1*04:01-HLA-DPB1*92:01', 'HLA-DPA1*04:01-HLA-DPB1*93:01',
         'HLA-DPA1*04:01-HLA-DPB1*94:01', 'HLA-DPA1*04:01-HLA-DPB1*95:01',
         'HLA-DPA1*04:01-HLA-DPB1*96:01', 'HLA-DPA1*04:01-HLA-DPB1*97:01', 'HLA-DPA1*04:01-HLA-DPB1*98:01', 'HLA-DPA1*04:01-HLA-DPB1*99:01',
         'HLA-DQA1*01:01-HLA-DQB1*02:01', 'HLA-DQA1*01:01-HLA-DQB1*02:02',
         'HLA-DQA1*01:01-HLA-DQB1*02:03', 'HLA-DQA1*01:01-HLA-DQB1*02:04', 'HLA-DQA1*01:01-HLA-DQB1*02:05', 'HLA-DQA1*01:01-HLA-DQB1*02:06',
         'HLA-DQA1*01:01-HLA-DQB1*03:01', 'HLA-DQA1*01:01-HLA-DQB1*03:02',
         'HLA-DQA1*01:01-HLA-DQB1*03:03', 'HLA-DQA1*01:01-HLA-DQB1*03:04', 'HLA-DQA1*01:01-HLA-DQB1*03:05', 'HLA-DQA1*01:01-HLA-DQB1*03:06',
         'HLA-DQA1*01:01-HLA-DQB1*03:07', 'HLA-DQA1*01:01-HLA-DQB1*03:08',
         'HLA-DQA1*01:01-HLA-DQB1*03:09', 'HLA-DQA1*01:01-HLA-DQB1*03:10', 'HLA-DQA1*01:01-HLA-DQB1*03:11', 'HLA-DQA1*01:01-HLA-DQB1*03:12',
         'HLA-DQA1*01:01-HLA-DQB1*03:13', 'HLA-DQA1*01:01-HLA-DQB1*03:14',
         'HLA-DQA1*01:01-HLA-DQB1*03:15', 'HLA-DQA1*01:01-HLA-DQB1*03:16', 'HLA-DQA1*01:01-HLA-DQB1*03:17', 'HLA-DQA1*01:01-HLA-DQB1*03:18',
         'HLA-DQA1*01:01-HLA-DQB1*03:19', 'HLA-DQA1*01:01-HLA-DQB1*03:20',
         'HLA-DQA1*01:01-HLA-DQB1*03:21', 'HLA-DQA1*01:01-HLA-DQB1*03:22', 'HLA-DQA1*01:01-HLA-DQB1*03:23', 'HLA-DQA1*01:01-HLA-DQB1*03:24',
         'HLA-DQA1*01:01-HLA-DQB1*03:25', 'HLA-DQA1*01:01-HLA-DQB1*03:26',
         'HLA-DQA1*01:01-HLA-DQB1*03:27', 'HLA-DQA1*01:01-HLA-DQB1*03:28', 'HLA-DQA1*01:01-HLA-DQB1*03:29', 'HLA-DQA1*01:01-HLA-DQB1*03:30',
         'HLA-DQA1*01:01-HLA-DQB1*03:31', 'HLA-DQA1*01:01-HLA-DQB1*03:32',
         'HLA-DQA1*01:01-HLA-DQB1*03:33', 'HLA-DQA1*01:01-HLA-DQB1*03:34', 'HLA-DQA1*01:01-HLA-DQB1*03:35', 'HLA-DQA1*01:01-HLA-DQB1*03:36',
         'HLA-DQA1*01:01-HLA-DQB1*03:37', 'HLA-DQA1*01:01-HLA-DQB1*03:38',
         'HLA-DQA1*01:01-HLA-DQB1*04:01', 'HLA-DQA1*01:01-HLA-DQB1*04:02', 'HLA-DQA1*01:01-HLA-DQB1*04:03', 'HLA-DQA1*01:01-HLA-DQB1*04:04',
         'HLA-DQA1*01:01-HLA-DQB1*04:05', 'HLA-DQA1*01:01-HLA-DQB1*04:06',
         'HLA-DQA1*01:01-HLA-DQB1*04:07', 'HLA-DQA1*01:01-HLA-DQB1*04:08', 'HLA-DQA1*01:01-HLA-DQB1*05:01', 'HLA-DQA1*01:01-HLA-DQB1*05:02',
         'HLA-DQA1*01:01-HLA-DQB1*05:03', 'HLA-DQA1*01:01-HLA-DQB1*05:05',
         'HLA-DQA1*01:01-HLA-DQB1*05:06', 'HLA-DQA1*01:01-HLA-DQB1*05:07', 'HLA-DQA1*01:01-HLA-DQB1*05:08', 'HLA-DQA1*01:01-HLA-DQB1*05:09',
         'HLA-DQA1*01:01-HLA-DQB1*05:10', 'HLA-DQA1*01:01-HLA-DQB1*05:11',
         'HLA-DQA1*01:01-HLA-DQB1*05:12', 'HLA-DQA1*01:01-HLA-DQB1*05:13', 'HLA-DQA1*01:01-HLA-DQB1*05:14', 'HLA-DQA1*01:01-HLA-DQB1*06:01',
         'HLA-DQA1*01:01-HLA-DQB1*06:02', 'HLA-DQA1*01:01-HLA-DQB1*06:03',
         'HLA-DQA1*01:01-HLA-DQB1*06:04', 'HLA-DQA1*01:01-HLA-DQB1*06:07', 'HLA-DQA1*01:01-HLA-DQB1*06:08', 'HLA-DQA1*01:01-HLA-DQB1*06:09',
         'HLA-DQA1*01:01-HLA-DQB1*06:10', 'HLA-DQA1*01:01-HLA-DQB1*06:11',
         'HLA-DQA1*01:01-HLA-DQB1*06:12', 'HLA-DQA1*01:01-HLA-DQB1*06:14', 'HLA-DQA1*01:01-HLA-DQB1*06:15', 'HLA-DQA1*01:01-HLA-DQB1*06:16',
         'HLA-DQA1*01:01-HLA-DQB1*06:17', 'HLA-DQA1*01:01-HLA-DQB1*06:18',
         'HLA-DQA1*01:01-HLA-DQB1*06:19', 'HLA-DQA1*01:01-HLA-DQB1*06:21', 'HLA-DQA1*01:01-HLA-DQB1*06:22', 'HLA-DQA1*01:01-HLA-DQB1*06:23',
         'HLA-DQA1*01:01-HLA-DQB1*06:24', 'HLA-DQA1*01:01-HLA-DQB1*06:25',
         'HLA-DQA1*01:01-HLA-DQB1*06:27', 'HLA-DQA1*01:01-HLA-DQB1*06:28', 'HLA-DQA1*01:01-HLA-DQB1*06:29', 'HLA-DQA1*01:01-HLA-DQB1*06:30',
         'HLA-DQA1*01:01-HLA-DQB1*06:31', 'HLA-DQA1*01:01-HLA-DQB1*06:32',
         'HLA-DQA1*01:01-HLA-DQB1*06:33', 'HLA-DQA1*01:01-HLA-DQB1*06:34', 'HLA-DQA1*01:01-HLA-DQB1*06:35', 'HLA-DQA1*01:01-HLA-DQB1*06:36',
         'HLA-DQA1*01:01-HLA-DQB1*06:37', 'HLA-DQA1*01:01-HLA-DQB1*06:38',
         'HLA-DQA1*01:01-HLA-DQB1*06:39', 'HLA-DQA1*01:01-HLA-DQB1*06:40', 'HLA-DQA1*01:01-HLA-DQB1*06:41', 'HLA-DQA1*01:01-HLA-DQB1*06:42',
         'HLA-DQA1*01:01-HLA-DQB1*06:43', 'HLA-DQA1*01:01-HLA-DQB1*06:44',
         'HLA-DQA1*01:02-HLA-DQB1*02:01', 'HLA-DQA1*01:02-HLA-DQB1*02:02', 'HLA-DQA1*01:02-HLA-DQB1*02:03', 'HLA-DQA1*01:02-HLA-DQB1*02:04',
         'HLA-DQA1*01:02-HLA-DQB1*02:05', 'HLA-DQA1*01:02-HLA-DQB1*02:06',
         'HLA-DQA1*01:02-HLA-DQB1*03:01', 'HLA-DQA1*01:02-HLA-DQB1*03:02', 'HLA-DQA1*01:02-HLA-DQB1*03:03', 'HLA-DQA1*01:02-HLA-DQB1*03:04',
         'HLA-DQA1*01:02-HLA-DQB1*03:05', 'HLA-DQA1*01:02-HLA-DQB1*03:06',
         'HLA-DQA1*01:02-HLA-DQB1*03:07', 'HLA-DQA1*01:02-HLA-DQB1*03:08', 'HLA-DQA1*01:02-HLA-DQB1*03:09', 'HLA-DQA1*01:02-HLA-DQB1*03:10',
         'HLA-DQA1*01:02-HLA-DQB1*03:11', 'HLA-DQA1*01:02-HLA-DQB1*03:12',
         'HLA-DQA1*01:02-HLA-DQB1*03:13', 'HLA-DQA1*01:02-HLA-DQB1*03:14', 'HLA-DQA1*01:02-HLA-DQB1*03:15', 'HLA-DQA1*01:02-HLA-DQB1*03:16',
         'HLA-DQA1*01:02-HLA-DQB1*03:17', 'HLA-DQA1*01:02-HLA-DQB1*03:18',
         'HLA-DQA1*01:02-HLA-DQB1*03:19', 'HLA-DQA1*01:02-HLA-DQB1*03:20', 'HLA-DQA1*01:02-HLA-DQB1*03:21', 'HLA-DQA1*01:02-HLA-DQB1*03:22',
         'HLA-DQA1*01:02-HLA-DQB1*03:23', 'HLA-DQA1*01:02-HLA-DQB1*03:24',
         'HLA-DQA1*01:02-HLA-DQB1*03:25', 'HLA-DQA1*01:02-HLA-DQB1*03:26', 'HLA-DQA1*01:02-HLA-DQB1*03:27', 'HLA-DQA1*01:02-HLA-DQB1*03:28',
         'HLA-DQA1*01:02-HLA-DQB1*03:29', 'HLA-DQA1*01:02-HLA-DQB1*03:30',
         'HLA-DQA1*01:02-HLA-DQB1*03:31', 'HLA-DQA1*01:02-HLA-DQB1*03:32', 'HLA-DQA1*01:02-HLA-DQB1*03:33', 'HLA-DQA1*01:02-HLA-DQB1*03:34',
         'HLA-DQA1*01:02-HLA-DQB1*03:35', 'HLA-DQA1*01:02-HLA-DQB1*03:36',
         'HLA-DQA1*01:02-HLA-DQB1*03:37', 'HLA-DQA1*01:02-HLA-DQB1*03:38', 'HLA-DQA1*01:02-HLA-DQB1*04:01', 'HLA-DQA1*01:02-HLA-DQB1*04:02',
         'HLA-DQA1*01:02-HLA-DQB1*04:03', 'HLA-DQA1*01:02-HLA-DQB1*04:04',
         'HLA-DQA1*01:02-HLA-DQB1*04:05', 'HLA-DQA1*01:02-HLA-DQB1*04:06', 'HLA-DQA1*01:02-HLA-DQB1*04:07', 'HLA-DQA1*01:02-HLA-DQB1*04:08',
         'HLA-DQA1*01:02-HLA-DQB1*05:01', 'HLA-DQA1*01:02-HLA-DQB1*05:02',
         'HLA-DQA1*01:02-HLA-DQB1*05:03', 'HLA-DQA1*01:02-HLA-DQB1*05:05', 'HLA-DQA1*01:02-HLA-DQB1*05:06', 'HLA-DQA1*01:02-HLA-DQB1*05:07',
         'HLA-DQA1*01:02-HLA-DQB1*05:08', 'HLA-DQA1*01:02-HLA-DQB1*05:09',
         'HLA-DQA1*01:02-HLA-DQB1*05:10', 'HLA-DQA1*01:02-HLA-DQB1*05:11', 'HLA-DQA1*01:02-HLA-DQB1*05:12', 'HLA-DQA1*01:02-HLA-DQB1*05:13',
         'HLA-DQA1*01:02-HLA-DQB1*05:14', 'HLA-DQA1*01:02-HLA-DQB1*06:01',
         'HLA-DQA1*01:02-HLA-DQB1*06:02', 'HLA-DQA1*01:02-HLA-DQB1*06:03', 'HLA-DQA1*01:02-HLA-DQB1*06:04', 'HLA-DQA1*01:02-HLA-DQB1*06:07',
         'HLA-DQA1*01:02-HLA-DQB1*06:08', 'HLA-DQA1*01:02-HLA-DQB1*06:09',
         'HLA-DQA1*01:02-HLA-DQB1*06:10', 'HLA-DQA1*01:02-HLA-DQB1*06:11', 'HLA-DQA1*01:02-HLA-DQB1*06:12', 'HLA-DQA1*01:02-HLA-DQB1*06:14',
         'HLA-DQA1*01:02-HLA-DQB1*06:15', 'HLA-DQA1*01:02-HLA-DQB1*06:16',
         'HLA-DQA1*01:02-HLA-DQB1*06:17', 'HLA-DQA1*01:02-HLA-DQB1*06:18', 'HLA-DQA1*01:02-HLA-DQB1*06:19', 'HLA-DQA1*01:02-HLA-DQB1*06:21',
         'HLA-DQA1*01:02-HLA-DQB1*06:22', 'HLA-DQA1*01:02-HLA-DQB1*06:23',
         'HLA-DQA1*01:02-HLA-DQB1*06:24', 'HLA-DQA1*01:02-HLA-DQB1*06:25', 'HLA-DQA1*01:02-HLA-DQB1*06:27', 'HLA-DQA1*01:02-HLA-DQB1*06:28',
         'HLA-DQA1*01:02-HLA-DQB1*06:29', 'HLA-DQA1*01:02-HLA-DQB1*06:30',
         'HLA-DQA1*01:02-HLA-DQB1*06:31', 'HLA-DQA1*01:02-HLA-DQB1*06:32', 'HLA-DQA1*01:02-HLA-DQB1*06:33', 'HLA-DQA1*01:02-HLA-DQB1*06:34',
         'HLA-DQA1*01:02-HLA-DQB1*06:35', 'HLA-DQA1*01:02-HLA-DQB1*06:36',
         'HLA-DQA1*01:02-HLA-DQB1*06:37', 'HLA-DQA1*01:02-HLA-DQB1*06:38', 'HLA-DQA1*01:02-HLA-DQB1*06:39', 'HLA-DQA1*01:02-HLA-DQB1*06:40',
         'HLA-DQA1*01:02-HLA-DQB1*06:41', 'HLA-DQA1*01:02-HLA-DQB1*06:42',
         'HLA-DQA1*01:02-HLA-DQB1*06:43', 'HLA-DQA1*01:02-HLA-DQB1*06:44', 'HLA-DQA1*01:03-HLA-DQB1*02:01', 'HLA-DQA1*01:03-HLA-DQB1*02:02',
         'HLA-DQA1*01:03-HLA-DQB1*02:03', 'HLA-DQA1*01:03-HLA-DQB1*02:04',
         'HLA-DQA1*01:03-HLA-DQB1*02:05', 'HLA-DQA1*01:03-HLA-DQB1*02:06', 'HLA-DQA1*01:03-HLA-DQB1*03:01', 'HLA-DQA1*01:03-HLA-DQB1*03:02',
         'HLA-DQA1*01:03-HLA-DQB1*03:03', 'HLA-DQA1*01:03-HLA-DQB1*03:04',
         'HLA-DQA1*01:03-HLA-DQB1*03:05', 'HLA-DQA1*01:03-HLA-DQB1*03:06', 'HLA-DQA1*01:03-HLA-DQB1*03:07', 'HLA-DQA1*01:03-HLA-DQB1*03:08',
         'HLA-DQA1*01:03-HLA-DQB1*03:09', 'HLA-DQA1*01:03-HLA-DQB1*03:10',
         'HLA-DQA1*01:03-HLA-DQB1*03:11', 'HLA-DQA1*01:03-HLA-DQB1*03:12', 'HLA-DQA1*01:03-HLA-DQB1*03:13', 'HLA-DQA1*01:03-HLA-DQB1*03:14',
         'HLA-DQA1*01:03-HLA-DQB1*03:15', 'HLA-DQA1*01:03-HLA-DQB1*03:16',
         'HLA-DQA1*01:03-HLA-DQB1*03:17', 'HLA-DQA1*01:03-HLA-DQB1*03:18', 'HLA-DQA1*01:03-HLA-DQB1*03:19', 'HLA-DQA1*01:03-HLA-DQB1*03:20',
         'HLA-DQA1*01:03-HLA-DQB1*03:21', 'HLA-DQA1*01:03-HLA-DQB1*03:22',
         'HLA-DQA1*01:03-HLA-DQB1*03:23', 'HLA-DQA1*01:03-HLA-DQB1*03:24', 'HLA-DQA1*01:03-HLA-DQB1*03:25', 'HLA-DQA1*01:03-HLA-DQB1*03:26',
         'HLA-DQA1*01:03-HLA-DQB1*03:27', 'HLA-DQA1*01:03-HLA-DQB1*03:28',
         'HLA-DQA1*01:03-HLA-DQB1*03:29', 'HLA-DQA1*01:03-HLA-DQB1*03:30', 'HLA-DQA1*01:03-HLA-DQB1*03:31', 'HLA-DQA1*01:03-HLA-DQB1*03:32',
         'HLA-DQA1*01:03-HLA-DQB1*03:33', 'HLA-DQA1*01:03-HLA-DQB1*03:34',
         'HLA-DQA1*01:03-HLA-DQB1*03:35', 'HLA-DQA1*01:03-HLA-DQB1*03:36', 'HLA-DQA1*01:03-HLA-DQB1*03:37', 'HLA-DQA1*01:03-HLA-DQB1*03:38',
         'HLA-DQA1*01:03-HLA-DQB1*04:01', 'HLA-DQA1*01:03-HLA-DQB1*04:02',
         'HLA-DQA1*01:03-HLA-DQB1*04:03', 'HLA-DQA1*01:03-HLA-DQB1*04:04', 'HLA-DQA1*01:03-HLA-DQB1*04:05', 'HLA-DQA1*01:03-HLA-DQB1*04:06',
         'HLA-DQA1*01:03-HLA-DQB1*04:07', 'HLA-DQA1*01:03-HLA-DQB1*04:08',
         'HLA-DQA1*01:03-HLA-DQB1*05:01', 'HLA-DQA1*01:03-HLA-DQB1*05:02', 'HLA-DQA1*01:03-HLA-DQB1*05:03', 'HLA-DQA1*01:03-HLA-DQB1*05:05',
         'HLA-DQA1*01:03-HLA-DQB1*05:06', 'HLA-DQA1*01:03-HLA-DQB1*05:07',
         'HLA-DQA1*01:03-HLA-DQB1*05:08', 'HLA-DQA1*01:03-HLA-DQB1*05:09', 'HLA-DQA1*01:03-HLA-DQB1*05:10', 'HLA-DQA1*01:03-HLA-DQB1*05:11',
         'HLA-DQA1*01:03-HLA-DQB1*05:12', 'HLA-DQA1*01:03-HLA-DQB1*05:13',
         'HLA-DQA1*01:03-HLA-DQB1*05:14', 'HLA-DQA1*01:03-HLA-DQB1*06:01', 'HLA-DQA1*01:03-HLA-DQB1*06:02', 'HLA-DQA1*01:03-HLA-DQB1*06:03',
         'HLA-DQA1*01:03-HLA-DQB1*06:04', 'HLA-DQA1*01:03-HLA-DQB1*06:07',
         'HLA-DQA1*01:03-HLA-DQB1*06:08', 'HLA-DQA1*01:03-HLA-DQB1*06:09', 'HLA-DQA1*01:03-HLA-DQB1*06:10', 'HLA-DQA1*01:03-HLA-DQB1*06:11',
         'HLA-DQA1*01:03-HLA-DQB1*06:12', 'HLA-DQA1*01:03-HLA-DQB1*06:14',
         'HLA-DQA1*01:03-HLA-DQB1*06:15', 'HLA-DQA1*01:03-HLA-DQB1*06:16', 'HLA-DQA1*01:03-HLA-DQB1*06:17', 'HLA-DQA1*01:03-HLA-DQB1*06:18',
         'HLA-DQA1*01:03-HLA-DQB1*06:19', 'HLA-DQA1*01:03-HLA-DQB1*06:21',
         'HLA-DQA1*01:03-HLA-DQB1*06:22', 'HLA-DQA1*01:03-HLA-DQB1*06:23', 'HLA-DQA1*01:03-HLA-DQB1*06:24', 'HLA-DQA1*01:03-HLA-DQB1*06:25',
         'HLA-DQA1*01:03-HLA-DQB1*06:27', 'HLA-DQA1*01:03-HLA-DQB1*06:28',
         'HLA-DQA1*01:03-HLA-DQB1*06:29', 'HLA-DQA1*01:03-HLA-DQB1*06:30', 'HLA-DQA1*01:03-HLA-DQB1*06:31', 'HLA-DQA1*01:03-HLA-DQB1*06:32',
         'HLA-DQA1*01:03-HLA-DQB1*06:33', 'HLA-DQA1*01:03-HLA-DQB1*06:34',
         'HLA-DQA1*01:03-HLA-DQB1*06:35', 'HLA-DQA1*01:03-HLA-DQB1*06:36', 'HLA-DQA1*01:03-HLA-DQB1*06:37', 'HLA-DQA1*01:03-HLA-DQB1*06:38',
         'HLA-DQA1*01:03-HLA-DQB1*06:39', 'HLA-DQA1*01:03-HLA-DQB1*06:40',
         'HLA-DQA1*01:03-HLA-DQB1*06:41', 'HLA-DQA1*01:03-HLA-DQB1*06:42', 'HLA-DQA1*01:03-HLA-DQB1*06:43', 'HLA-DQA1*01:03-HLA-DQB1*06:44',
         'HLA-DQA1*01:04-HLA-DQB1*02:01', 'HLA-DQA1*01:04-HLA-DQB1*02:02',
         'HLA-DQA1*01:04-HLA-DQB1*02:03', 'HLA-DQA1*01:04-HLA-DQB1*02:04', 'HLA-DQA1*01:04-HLA-DQB1*02:05', 'HLA-DQA1*01:04-HLA-DQB1*02:06',
         'HLA-DQA1*01:04-HLA-DQB1*03:01', 'HLA-DQA1*01:04-HLA-DQB1*03:02',
         'HLA-DQA1*01:04-HLA-DQB1*03:03', 'HLA-DQA1*01:04-HLA-DQB1*03:04', 'HLA-DQA1*01:04-HLA-DQB1*03:05', 'HLA-DQA1*01:04-HLA-DQB1*03:06',
         'HLA-DQA1*01:04-HLA-DQB1*03:07', 'HLA-DQA1*01:04-HLA-DQB1*03:08',
         'HLA-DQA1*01:04-HLA-DQB1*03:09', 'HLA-DQA1*01:04-HLA-DQB1*03:10', 'HLA-DQA1*01:04-HLA-DQB1*03:11', 'HLA-DQA1*01:04-HLA-DQB1*03:12',
         'HLA-DQA1*01:04-HLA-DQB1*03:13', 'HLA-DQA1*01:04-HLA-DQB1*03:14',
         'HLA-DQA1*01:04-HLA-DQB1*03:15', 'HLA-DQA1*01:04-HLA-DQB1*03:16', 'HLA-DQA1*01:04-HLA-DQB1*03:17', 'HLA-DQA1*01:04-HLA-DQB1*03:18',
         'HLA-DQA1*01:04-HLA-DQB1*03:19', 'HLA-DQA1*01:04-HLA-DQB1*03:20',
         'HLA-DQA1*01:04-HLA-DQB1*03:21', 'HLA-DQA1*01:04-HLA-DQB1*03:22', 'HLA-DQA1*01:04-HLA-DQB1*03:23', 'HLA-DQA1*01:04-HLA-DQB1*03:24',
         'HLA-DQA1*01:04-HLA-DQB1*03:25', 'HLA-DQA1*01:04-HLA-DQB1*03:26',
         'HLA-DQA1*01:04-HLA-DQB1*03:27', 'HLA-DQA1*01:04-HLA-DQB1*03:28', 'HLA-DQA1*01:04-HLA-DQB1*03:29', 'HLA-DQA1*01:04-HLA-DQB1*03:30',
         'HLA-DQA1*01:04-HLA-DQB1*03:31', 'HLA-DQA1*01:04-HLA-DQB1*03:32',
         'HLA-DQA1*01:04-HLA-DQB1*03:33', 'HLA-DQA1*01:04-HLA-DQB1*03:34', 'HLA-DQA1*01:04-HLA-DQB1*03:35', 'HLA-DQA1*01:04-HLA-DQB1*03:36',
         'HLA-DQA1*01:04-HLA-DQB1*03:37', 'HLA-DQA1*01:04-HLA-DQB1*03:38',
         'HLA-DQA1*01:04-HLA-DQB1*04:01', 'HLA-DQA1*01:04-HLA-DQB1*04:02', 'HLA-DQA1*01:04-HLA-DQB1*04:03', 'HLA-DQA1*01:04-HLA-DQB1*04:04',
         'HLA-DQA1*01:04-HLA-DQB1*04:05', 'HLA-DQA1*01:04-HLA-DQB1*04:06',
         'HLA-DQA1*01:04-HLA-DQB1*04:07', 'HLA-DQA1*01:04-HLA-DQB1*04:08', 'HLA-DQA1*01:04-HLA-DQB1*05:01', 'HLA-DQA1*01:04-HLA-DQB1*05:02',
         'HLA-DQA1*01:04-HLA-DQB1*05:03', 'HLA-DQA1*01:04-HLA-DQB1*05:05',
         'HLA-DQA1*01:04-HLA-DQB1*05:06', 'HLA-DQA1*01:04-HLA-DQB1*05:07', 'HLA-DQA1*01:04-HLA-DQB1*05:08', 'HLA-DQA1*01:04-HLA-DQB1*05:09',
         'HLA-DQA1*01:04-HLA-DQB1*05:10', 'HLA-DQA1*01:04-HLA-DQB1*05:11',
         'HLA-DQA1*01:04-HLA-DQB1*05:12', 'HLA-DQA1*01:04-HLA-DQB1*05:13', 'HLA-DQA1*01:04-HLA-DQB1*05:14', 'HLA-DQA1*01:04-HLA-DQB1*06:01',
         'HLA-DQA1*01:04-HLA-DQB1*06:02', 'HLA-DQA1*01:04-HLA-DQB1*06:03',
         'HLA-DQA1*01:04-HLA-DQB1*06:04', 'HLA-DQA1*01:04-HLA-DQB1*06:07', 'HLA-DQA1*01:04-HLA-DQB1*06:08', 'HLA-DQA1*01:04-HLA-DQB1*06:09',
         'HLA-DQA1*01:04-HLA-DQB1*06:10', 'HLA-DQA1*01:04-HLA-DQB1*06:11',
         'HLA-DQA1*01:04-HLA-DQB1*06:12', 'HLA-DQA1*01:04-HLA-DQB1*06:14', 'HLA-DQA1*01:04-HLA-DQB1*06:15', 'HLA-DQA1*01:04-HLA-DQB1*06:16',
         'HLA-DQA1*01:04-HLA-DQB1*06:17', 'HLA-DQA1*01:04-HLA-DQB1*06:18',
         'HLA-DQA1*01:04-HLA-DQB1*06:19', 'HLA-DQA1*01:04-HLA-DQB1*06:21', 'HLA-DQA1*01:04-HLA-DQB1*06:22', 'HLA-DQA1*01:04-HLA-DQB1*06:23',
         'HLA-DQA1*01:04-HLA-DQB1*06:24', 'HLA-DQA1*01:04-HLA-DQB1*06:25',
         'HLA-DQA1*01:04-HLA-DQB1*06:27', 'HLA-DQA1*01:04-HLA-DQB1*06:28', 'HLA-DQA1*01:04-HLA-DQB1*06:29', 'HLA-DQA1*01:04-HLA-DQB1*06:30',
         'HLA-DQA1*01:04-HLA-DQB1*06:31', 'HLA-DQA1*01:04-HLA-DQB1*06:32',
         'HLA-DQA1*01:04-HLA-DQB1*06:33', 'HLA-DQA1*01:04-HLA-DQB1*06:34', 'HLA-DQA1*01:04-HLA-DQB1*06:35', 'HLA-DQA1*01:04-HLA-DQB1*06:36',
         'HLA-DQA1*01:04-HLA-DQB1*06:37', 'HLA-DQA1*01:04-HLA-DQB1*06:38',
         'HLA-DQA1*01:04-HLA-DQB1*06:39', 'HLA-DQA1*01:04-HLA-DQB1*06:40', 'HLA-DQA1*01:04-HLA-DQB1*06:41', 'HLA-DQA1*01:04-HLA-DQB1*06:42',
         'HLA-DQA1*01:04-HLA-DQB1*06:43', 'HLA-DQA1*01:04-HLA-DQB1*06:44',
         'HLA-DQA1*01:05-HLA-DQB1*02:01', 'HLA-DQA1*01:05-HLA-DQB1*02:02', 'HLA-DQA1*01:05-HLA-DQB1*02:03', 'HLA-DQA1*01:05-HLA-DQB1*02:04',
         'HLA-DQA1*01:05-HLA-DQB1*02:05', 'HLA-DQA1*01:05-HLA-DQB1*02:06',
         'HLA-DQA1*01:05-HLA-DQB1*03:01', 'HLA-DQA1*01:05-HLA-DQB1*03:02', 'HLA-DQA1*01:05-HLA-DQB1*03:03', 'HLA-DQA1*01:05-HLA-DQB1*03:04',
         'HLA-DQA1*01:05-HLA-DQB1*03:05', 'HLA-DQA1*01:05-HLA-DQB1*03:06',
         'HLA-DQA1*01:05-HLA-DQB1*03:07', 'HLA-DQA1*01:05-HLA-DQB1*03:08', 'HLA-DQA1*01:05-HLA-DQB1*03:09', 'HLA-DQA1*01:05-HLA-DQB1*03:10',
         'HLA-DQA1*01:05-HLA-DQB1*03:11', 'HLA-DQA1*01:05-HLA-DQB1*03:12',
         'HLA-DQA1*01:05-HLA-DQB1*03:13', 'HLA-DQA1*01:05-HLA-DQB1*03:14', 'HLA-DQA1*01:05-HLA-DQB1*03:15', 'HLA-DQA1*01:05-HLA-DQB1*03:16',
         'HLA-DQA1*01:05-HLA-DQB1*03:17', 'HLA-DQA1*01:05-HLA-DQB1*03:18',
         'HLA-DQA1*01:05-HLA-DQB1*03:19', 'HLA-DQA1*01:05-HLA-DQB1*03:20', 'HLA-DQA1*01:05-HLA-DQB1*03:21', 'HLA-DQA1*01:05-HLA-DQB1*03:22',
         'HLA-DQA1*01:05-HLA-DQB1*03:23', 'HLA-DQA1*01:05-HLA-DQB1*03:24',
         'HLA-DQA1*01:05-HLA-DQB1*03:25', 'HLA-DQA1*01:05-HLA-DQB1*03:26', 'HLA-DQA1*01:05-HLA-DQB1*03:27', 'HLA-DQA1*01:05-HLA-DQB1*03:28',
         'HLA-DQA1*01:05-HLA-DQB1*03:29', 'HLA-DQA1*01:05-HLA-DQB1*03:30',
         'HLA-DQA1*01:05-HLA-DQB1*03:31', 'HLA-DQA1*01:05-HLA-DQB1*03:32', 'HLA-DQA1*01:05-HLA-DQB1*03:33', 'HLA-DQA1*01:05-HLA-DQB1*03:34',
         'HLA-DQA1*01:05-HLA-DQB1*03:35', 'HLA-DQA1*01:05-HLA-DQB1*03:36',
         'HLA-DQA1*01:05-HLA-DQB1*03:37', 'HLA-DQA1*01:05-HLA-DQB1*03:38', 'HLA-DQA1*01:05-HLA-DQB1*04:01', 'HLA-DQA1*01:05-HLA-DQB1*04:02',
         'HLA-DQA1*01:05-HLA-DQB1*04:03', 'HLA-DQA1*01:05-HLA-DQB1*04:04',
         'HLA-DQA1*01:05-HLA-DQB1*04:05', 'HLA-DQA1*01:05-HLA-DQB1*04:06', 'HLA-DQA1*01:05-HLA-DQB1*04:07', 'HLA-DQA1*01:05-HLA-DQB1*04:08',
         'HLA-DQA1*01:05-HLA-DQB1*05:01', 'HLA-DQA1*01:05-HLA-DQB1*05:02',
         'HLA-DQA1*01:05-HLA-DQB1*05:03', 'HLA-DQA1*01:05-HLA-DQB1*05:05', 'HLA-DQA1*01:05-HLA-DQB1*05:06', 'HLA-DQA1*01:05-HLA-DQB1*05:07',
         'HLA-DQA1*01:05-HLA-DQB1*05:08', 'HLA-DQA1*01:05-HLA-DQB1*05:09',
         'HLA-DQA1*01:05-HLA-DQB1*05:10', 'HLA-DQA1*01:05-HLA-DQB1*05:11', 'HLA-DQA1*01:05-HLA-DQB1*05:12', 'HLA-DQA1*01:05-HLA-DQB1*05:13',
         'HLA-DQA1*01:05-HLA-DQB1*05:14', 'HLA-DQA1*01:05-HLA-DQB1*06:01',
         'HLA-DQA1*01:05-HLA-DQB1*06:02', 'HLA-DQA1*01:05-HLA-DQB1*06:03', 'HLA-DQA1*01:05-HLA-DQB1*06:04', 'HLA-DQA1*01:05-HLA-DQB1*06:07',
         'HLA-DQA1*01:05-HLA-DQB1*06:08', 'HLA-DQA1*01:05-HLA-DQB1*06:09',
         'HLA-DQA1*01:05-HLA-DQB1*06:10', 'HLA-DQA1*01:05-HLA-DQB1*06:11', 'HLA-DQA1*01:05-HLA-DQB1*06:12', 'HLA-DQA1*01:05-HLA-DQB1*06:14',
         'HLA-DQA1*01:05-HLA-DQB1*06:15', 'HLA-DQA1*01:05-HLA-DQB1*06:16',
         'HLA-DQA1*01:05-HLA-DQB1*06:17', 'HLA-DQA1*01:05-HLA-DQB1*06:18', 'HLA-DQA1*01:05-HLA-DQB1*06:19', 'HLA-DQA1*01:05-HLA-DQB1*06:21',
         'HLA-DQA1*01:05-HLA-DQB1*06:22', 'HLA-DQA1*01:05-HLA-DQB1*06:23',
         'HLA-DQA1*01:05-HLA-DQB1*06:24', 'HLA-DQA1*01:05-HLA-DQB1*06:25', 'HLA-DQA1*01:05-HLA-DQB1*06:27', 'HLA-DQA1*01:05-HLA-DQB1*06:28',
         'HLA-DQA1*01:05-HLA-DQB1*06:29', 'HLA-DQA1*01:05-HLA-DQB1*06:30',
         'HLA-DQA1*01:05-HLA-DQB1*06:31', 'HLA-DQA1*01:05-HLA-DQB1*06:32', 'HLA-DQA1*01:05-HLA-DQB1*06:33', 'HLA-DQA1*01:05-HLA-DQB1*06:34',
         'HLA-DQA1*01:05-HLA-DQB1*06:35', 'HLA-DQA1*01:05-HLA-DQB1*06:36',
         'HLA-DQA1*01:05-HLA-DQB1*06:37', 'HLA-DQA1*01:05-HLA-DQB1*06:38', 'HLA-DQA1*01:05-HLA-DQB1*06:39', 'HLA-DQA1*01:05-HLA-DQB1*06:40',
         'HLA-DQA1*01:05-HLA-DQB1*06:41', 'HLA-DQA1*01:05-HLA-DQB1*06:42',
         'HLA-DQA1*01:05-HLA-DQB1*06:43', 'HLA-DQA1*01:05-HLA-DQB1*06:44', 'HLA-DQA1*01:06-HLA-DQB1*02:01', 'HLA-DQA1*01:06-HLA-DQB1*02:02',
         'HLA-DQA1*01:06-HLA-DQB1*02:03', 'HLA-DQA1*01:06-HLA-DQB1*02:04',
         'HLA-DQA1*01:06-HLA-DQB1*02:05', 'HLA-DQA1*01:06-HLA-DQB1*02:06', 'HLA-DQA1*01:06-HLA-DQB1*03:01', 'HLA-DQA1*01:06-HLA-DQB1*03:02',
         'HLA-DQA1*01:06-HLA-DQB1*03:03', 'HLA-DQA1*01:06-HLA-DQB1*03:04',
         'HLA-DQA1*01:06-HLA-DQB1*03:05', 'HLA-DQA1*01:06-HLA-DQB1*03:06', 'HLA-DQA1*01:06-HLA-DQB1*03:07', 'HLA-DQA1*01:06-HLA-DQB1*03:08',
         'HLA-DQA1*01:06-HLA-DQB1*03:09', 'HLA-DQA1*01:06-HLA-DQB1*03:10',
         'HLA-DQA1*01:06-HLA-DQB1*03:11', 'HLA-DQA1*01:06-HLA-DQB1*03:12', 'HLA-DQA1*01:06-HLA-DQB1*03:13', 'HLA-DQA1*01:06-HLA-DQB1*03:14',
         'HLA-DQA1*01:06-HLA-DQB1*03:15', 'HLA-DQA1*01:06-HLA-DQB1*03:16',
         'HLA-DQA1*01:06-HLA-DQB1*03:17', 'HLA-DQA1*01:06-HLA-DQB1*03:18', 'HLA-DQA1*01:06-HLA-DQB1*03:19', 'HLA-DQA1*01:06-HLA-DQB1*03:20',
         'HLA-DQA1*01:06-HLA-DQB1*03:21', 'HLA-DQA1*01:06-HLA-DQB1*03:22',
         'HLA-DQA1*01:06-HLA-DQB1*03:23', 'HLA-DQA1*01:06-HLA-DQB1*03:24', 'HLA-DQA1*01:06-HLA-DQB1*03:25', 'HLA-DQA1*01:06-HLA-DQB1*03:26',
         'HLA-DQA1*01:06-HLA-DQB1*03:27', 'HLA-DQA1*01:06-HLA-DQB1*03:28',
         'HLA-DQA1*01:06-HLA-DQB1*03:29', 'HLA-DQA1*01:06-HLA-DQB1*03:30', 'HLA-DQA1*01:06-HLA-DQB1*03:31', 'HLA-DQA1*01:06-HLA-DQB1*03:32',
         'HLA-DQA1*01:06-HLA-DQB1*03:33', 'HLA-DQA1*01:06-HLA-DQB1*03:34',
         'HLA-DQA1*01:06-HLA-DQB1*03:35', 'HLA-DQA1*01:06-HLA-DQB1*03:36', 'HLA-DQA1*01:06-HLA-DQB1*03:37', 'HLA-DQA1*01:06-HLA-DQB1*03:38',
         'HLA-DQA1*01:06-HLA-DQB1*04:01', 'HLA-DQA1*01:06-HLA-DQB1*04:02',
         'HLA-DQA1*01:06-HLA-DQB1*04:03', 'HLA-DQA1*01:06-HLA-DQB1*04:04', 'HLA-DQA1*01:06-HLA-DQB1*04:05', 'HLA-DQA1*01:06-HLA-DQB1*04:06',
         'HLA-DQA1*01:06-HLA-DQB1*04:07', 'HLA-DQA1*01:06-HLA-DQB1*04:08',
         'HLA-DQA1*01:06-HLA-DQB1*05:01', 'HLA-DQA1*01:06-HLA-DQB1*05:02', 'HLA-DQA1*01:06-HLA-DQB1*05:03', 'HLA-DQA1*01:06-HLA-DQB1*05:05',
         'HLA-DQA1*01:06-HLA-DQB1*05:06', 'HLA-DQA1*01:06-HLA-DQB1*05:07',
         'HLA-DQA1*01:06-HLA-DQB1*05:08', 'HLA-DQA1*01:06-HLA-DQB1*05:09', 'HLA-DQA1*01:06-HLA-DQB1*05:10', 'HLA-DQA1*01:06-HLA-DQB1*05:11',
         'HLA-DQA1*01:06-HLA-DQB1*05:12', 'HLA-DQA1*01:06-HLA-DQB1*05:13',
         'HLA-DQA1*01:06-HLA-DQB1*05:14', 'HLA-DQA1*01:06-HLA-DQB1*06:01', 'HLA-DQA1*01:06-HLA-DQB1*06:02', 'HLA-DQA1*01:06-HLA-DQB1*06:03',
         'HLA-DQA1*01:06-HLA-DQB1*06:04', 'HLA-DQA1*01:06-HLA-DQB1*06:07',
         'HLA-DQA1*01:06-HLA-DQB1*06:08', 'HLA-DQA1*01:06-HLA-DQB1*06:09', 'HLA-DQA1*01:06-HLA-DQB1*06:10', 'HLA-DQA1*01:06-HLA-DQB1*06:11',
         'HLA-DQA1*01:06-HLA-DQB1*06:12', 'HLA-DQA1*01:06-HLA-DQB1*06:14',
         'HLA-DQA1*01:06-HLA-DQB1*06:15', 'HLA-DQA1*01:06-HLA-DQB1*06:16', 'HLA-DQA1*01:06-HLA-DQB1*06:17', 'HLA-DQA1*01:06-HLA-DQB1*06:18',
         'HLA-DQA1*01:06-HLA-DQB1*06:19', 'HLA-DQA1*01:06-HLA-DQB1*06:21',
         'HLA-DQA1*01:06-HLA-DQB1*06:22', 'HLA-DQA1*01:06-HLA-DQB1*06:23', 'HLA-DQA1*01:06-HLA-DQB1*06:24', 'HLA-DQA1*01:06-HLA-DQB1*06:25',
         'HLA-DQA1*01:06-HLA-DQB1*06:27', 'HLA-DQA1*01:06-HLA-DQB1*06:28',
         'HLA-DQA1*01:06-HLA-DQB1*06:29', 'HLA-DQA1*01:06-HLA-DQB1*06:30', 'HLA-DQA1*01:06-HLA-DQB1*06:31', 'HLA-DQA1*01:06-HLA-DQB1*06:32',
         'HLA-DQA1*01:06-HLA-DQB1*06:33', 'HLA-DQA1*01:06-HLA-DQB1*06:34',
         'HLA-DQA1*01:06-HLA-DQB1*06:35', 'HLA-DQA1*01:06-HLA-DQB1*06:36', 'HLA-DQA1*01:06-HLA-DQB1*06:37', 'HLA-DQA1*01:06-HLA-DQB1*06:38',
         'HLA-DQA1*01:06-HLA-DQB1*06:39', 'HLA-DQA1*01:06-HLA-DQB1*06:40',
         'HLA-DQA1*01:06-HLA-DQB1*06:41', 'HLA-DQA1*01:06-HLA-DQB1*06:42', 'HLA-DQA1*01:06-HLA-DQB1*06:43', 'HLA-DQA1*01:06-HLA-DQB1*06:44',
         'HLA-DQA1*01:07-HLA-DQB1*02:01', 'HLA-DQA1*01:07-HLA-DQB1*02:02',
         'HLA-DQA1*01:07-HLA-DQB1*02:03', 'HLA-DQA1*01:07-HLA-DQB1*02:04', 'HLA-DQA1*01:07-HLA-DQB1*02:05', 'HLA-DQA1*01:07-HLA-DQB1*02:06',
         'HLA-DQA1*01:07-HLA-DQB1*03:01', 'HLA-DQA1*01:07-HLA-DQB1*03:02',
         'HLA-DQA1*01:07-HLA-DQB1*03:03', 'HLA-DQA1*01:07-HLA-DQB1*03:04', 'HLA-DQA1*01:07-HLA-DQB1*03:05', 'HLA-DQA1*01:07-HLA-DQB1*03:06',
         'HLA-DQA1*01:07-HLA-DQB1*03:07', 'HLA-DQA1*01:07-HLA-DQB1*03:08',
         'HLA-DQA1*01:07-HLA-DQB1*03:09', 'HLA-DQA1*01:07-HLA-DQB1*03:10', 'HLA-DQA1*01:07-HLA-DQB1*03:11', 'HLA-DQA1*01:07-HLA-DQB1*03:12',
         'HLA-DQA1*01:07-HLA-DQB1*03:13', 'HLA-DQA1*01:07-HLA-DQB1*03:14',
         'HLA-DQA1*01:07-HLA-DQB1*03:15', 'HLA-DQA1*01:07-HLA-DQB1*03:16', 'HLA-DQA1*01:07-HLA-DQB1*03:17', 'HLA-DQA1*01:07-HLA-DQB1*03:18',
         'HLA-DQA1*01:07-HLA-DQB1*03:19', 'HLA-DQA1*01:07-HLA-DQB1*03:20',
         'HLA-DQA1*01:07-HLA-DQB1*03:21', 'HLA-DQA1*01:07-HLA-DQB1*03:22', 'HLA-DQA1*01:07-HLA-DQB1*03:23', 'HLA-DQA1*01:07-HLA-DQB1*03:24',
         'HLA-DQA1*01:07-HLA-DQB1*03:25', 'HLA-DQA1*01:07-HLA-DQB1*03:26',
         'HLA-DQA1*01:07-HLA-DQB1*03:27', 'HLA-DQA1*01:07-HLA-DQB1*03:28', 'HLA-DQA1*01:07-HLA-DQB1*03:29', 'HLA-DQA1*01:07-HLA-DQB1*03:30',
         'HLA-DQA1*01:07-HLA-DQB1*03:31', 'HLA-DQA1*01:07-HLA-DQB1*03:32',
         'HLA-DQA1*01:07-HLA-DQB1*03:33', 'HLA-DQA1*01:07-HLA-DQB1*03:34', 'HLA-DQA1*01:07-HLA-DQB1*03:35', 'HLA-DQA1*01:07-HLA-DQB1*03:36',
         'HLA-DQA1*01:07-HLA-DQB1*03:37', 'HLA-DQA1*01:07-HLA-DQB1*03:38',
         'HLA-DQA1*01:07-HLA-DQB1*04:01', 'HLA-DQA1*01:07-HLA-DQB1*04:02', 'HLA-DQA1*01:07-HLA-DQB1*04:03', 'HLA-DQA1*01:07-HLA-DQB1*04:04',
         'HLA-DQA1*01:07-HLA-DQB1*04:05', 'HLA-DQA1*01:07-HLA-DQB1*04:06',
         'HLA-DQA1*01:07-HLA-DQB1*04:07', 'HLA-DQA1*01:07-HLA-DQB1*04:08', 'HLA-DQA1*01:07-HLA-DQB1*05:01', 'HLA-DQA1*01:07-HLA-DQB1*05:02',
         'HLA-DQA1*01:07-HLA-DQB1*05:03', 'HLA-DQA1*01:07-HLA-DQB1*05:05',
         'HLA-DQA1*01:07-HLA-DQB1*05:06', 'HLA-DQA1*01:07-HLA-DQB1*05:07', 'HLA-DQA1*01:07-HLA-DQB1*05:08', 'HLA-DQA1*01:07-HLA-DQB1*05:09',
         'HLA-DQA1*01:07-HLA-DQB1*05:10', 'HLA-DQA1*01:07-HLA-DQB1*05:11',
         'HLA-DQA1*01:07-HLA-DQB1*05:12', 'HLA-DQA1*01:07-HLA-DQB1*05:13', 'HLA-DQA1*01:07-HLA-DQB1*05:14', 'HLA-DQA1*01:07-HLA-DQB1*06:01',
         'HLA-DQA1*01:07-HLA-DQB1*06:02', 'HLA-DQA1*01:07-HLA-DQB1*06:03',
         'HLA-DQA1*01:07-HLA-DQB1*06:04', 'HLA-DQA1*01:07-HLA-DQB1*06:07', 'HLA-DQA1*01:07-HLA-DQB1*06:08', 'HLA-DQA1*01:07-HLA-DQB1*06:09',
         'HLA-DQA1*01:07-HLA-DQB1*06:10', 'HLA-DQA1*01:07-HLA-DQB1*06:11',
         'HLA-DQA1*01:07-HLA-DQB1*06:12', 'HLA-DQA1*01:07-HLA-DQB1*06:14', 'HLA-DQA1*01:07-HLA-DQB1*06:15', 'HLA-DQA1*01:07-HLA-DQB1*06:16',
         'HLA-DQA1*01:07-HLA-DQB1*06:17', 'HLA-DQA1*01:07-HLA-DQB1*06:18',
         'HLA-DQA1*01:07-HLA-DQB1*06:19', 'HLA-DQA1*01:07-HLA-DQB1*06:21', 'HLA-DQA1*01:07-HLA-DQB1*06:22', 'HLA-DQA1*01:07-HLA-DQB1*06:23',
         'HLA-DQA1*01:07-HLA-DQB1*06:24', 'HLA-DQA1*01:07-HLA-DQB1*06:25',
         'HLA-DQA1*01:07-HLA-DQB1*06:27', 'HLA-DQA1*01:07-HLA-DQB1*06:28', 'HLA-DQA1*01:07-HLA-DQB1*06:29', 'HLA-DQA1*01:07-HLA-DQB1*06:30',
         'HLA-DQA1*01:07-HLA-DQB1*06:31', 'HLA-DQA1*01:07-HLA-DQB1*06:32',
         'HLA-DQA1*01:07-HLA-DQB1*06:33', 'HLA-DQA1*01:07-HLA-DQB1*06:34', 'HLA-DQA1*01:07-HLA-DQB1*06:35', 'HLA-DQA1*01:07-HLA-DQB1*06:36',
         'HLA-DQA1*01:07-HLA-DQB1*06:37', 'HLA-DQA1*01:07-HLA-DQB1*06:38',
         'HLA-DQA1*01:07-HLA-DQB1*06:39', 'HLA-DQA1*01:07-HLA-DQB1*06:40', 'HLA-DQA1*01:07-HLA-DQB1*06:41', 'HLA-DQA1*01:07-HLA-DQB1*06:42',
         'HLA-DQA1*01:07-HLA-DQB1*06:43', 'HLA-DQA1*01:07-HLA-DQB1*06:44',
         'HLA-DQA1*01:08-HLA-DQB1*02:01', 'HLA-DQA1*01:08-HLA-DQB1*02:02', 'HLA-DQA1*01:08-HLA-DQB1*02:03', 'HLA-DQA1*01:08-HLA-DQB1*02:04',
         'HLA-DQA1*01:08-HLA-DQB1*02:05', 'HLA-DQA1*01:08-HLA-DQB1*02:06',
         'HLA-DQA1*01:08-HLA-DQB1*03:01', 'HLA-DQA1*01:08-HLA-DQB1*03:02', 'HLA-DQA1*01:08-HLA-DQB1*03:03', 'HLA-DQA1*01:08-HLA-DQB1*03:04',
         'HLA-DQA1*01:08-HLA-DQB1*03:05', 'HLA-DQA1*01:08-HLA-DQB1*03:06',
         'HLA-DQA1*01:08-HLA-DQB1*03:07', 'HLA-DQA1*01:08-HLA-DQB1*03:08', 'HLA-DQA1*01:08-HLA-DQB1*03:09', 'HLA-DQA1*01:08-HLA-DQB1*03:10',
         'HLA-DQA1*01:08-HLA-DQB1*03:11', 'HLA-DQA1*01:08-HLA-DQB1*03:12',
         'HLA-DQA1*01:08-HLA-DQB1*03:13', 'HLA-DQA1*01:08-HLA-DQB1*03:14', 'HLA-DQA1*01:08-HLA-DQB1*03:15', 'HLA-DQA1*01:08-HLA-DQB1*03:16',
         'HLA-DQA1*01:08-HLA-DQB1*03:17', 'HLA-DQA1*01:08-HLA-DQB1*03:18',
         'HLA-DQA1*01:08-HLA-DQB1*03:19', 'HLA-DQA1*01:08-HLA-DQB1*03:20', 'HLA-DQA1*01:08-HLA-DQB1*03:21', 'HLA-DQA1*01:08-HLA-DQB1*03:22',
         'HLA-DQA1*01:08-HLA-DQB1*03:23', 'HLA-DQA1*01:08-HLA-DQB1*03:24',
         'HLA-DQA1*01:08-HLA-DQB1*03:25', 'HLA-DQA1*01:08-HLA-DQB1*03:26', 'HLA-DQA1*01:08-HLA-DQB1*03:27', 'HLA-DQA1*01:08-HLA-DQB1*03:28',
         'HLA-DQA1*01:08-HLA-DQB1*03:29', 'HLA-DQA1*01:08-HLA-DQB1*03:30',
         'HLA-DQA1*01:08-HLA-DQB1*03:31', 'HLA-DQA1*01:08-HLA-DQB1*03:32', 'HLA-DQA1*01:08-HLA-DQB1*03:33', 'HLA-DQA1*01:08-HLA-DQB1*03:34',
         'HLA-DQA1*01:08-HLA-DQB1*03:35', 'HLA-DQA1*01:08-HLA-DQB1*03:36',
         'HLA-DQA1*01:08-HLA-DQB1*03:37', 'HLA-DQA1*01:08-HLA-DQB1*03:38', 'HLA-DQA1*01:08-HLA-DQB1*04:01', 'HLA-DQA1*01:08-HLA-DQB1*04:02',
         'HLA-DQA1*01:08-HLA-DQB1*04:03', 'HLA-DQA1*01:08-HLA-DQB1*04:04',
         'HLA-DQA1*01:08-HLA-DQB1*04:05', 'HLA-DQA1*01:08-HLA-DQB1*04:06', 'HLA-DQA1*01:08-HLA-DQB1*04:07', 'HLA-DQA1*01:08-HLA-DQB1*04:08',
         'HLA-DQA1*01:08-HLA-DQB1*05:01', 'HLA-DQA1*01:08-HLA-DQB1*05:02',
         'HLA-DQA1*01:08-HLA-DQB1*05:03', 'HLA-DQA1*01:08-HLA-DQB1*05:05', 'HLA-DQA1*01:08-HLA-DQB1*05:06', 'HLA-DQA1*01:08-HLA-DQB1*05:07',
         'HLA-DQA1*01:08-HLA-DQB1*05:08', 'HLA-DQA1*01:08-HLA-DQB1*05:09',
         'HLA-DQA1*01:08-HLA-DQB1*05:10', 'HLA-DQA1*01:08-HLA-DQB1*05:11', 'HLA-DQA1*01:08-HLA-DQB1*05:12', 'HLA-DQA1*01:08-HLA-DQB1*05:13',
         'HLA-DQA1*01:08-HLA-DQB1*05:14', 'HLA-DQA1*01:08-HLA-DQB1*06:01',
         'HLA-DQA1*01:08-HLA-DQB1*06:02', 'HLA-DQA1*01:08-HLA-DQB1*06:03', 'HLA-DQA1*01:08-HLA-DQB1*06:04', 'HLA-DQA1*01:08-HLA-DQB1*06:07',
         'HLA-DQA1*01:08-HLA-DQB1*06:08', 'HLA-DQA1*01:08-HLA-DQB1*06:09',
         'HLA-DQA1*01:08-HLA-DQB1*06:10', 'HLA-DQA1*01:08-HLA-DQB1*06:11', 'HLA-DQA1*01:08-HLA-DQB1*06:12', 'HLA-DQA1*01:08-HLA-DQB1*06:14',
         'HLA-DQA1*01:08-HLA-DQB1*06:15', 'HLA-DQA1*01:08-HLA-DQB1*06:16',
         'HLA-DQA1*01:08-HLA-DQB1*06:17', 'HLA-DQA1*01:08-HLA-DQB1*06:18', 'HLA-DQA1*01:08-HLA-DQB1*06:19', 'HLA-DQA1*01:08-HLA-DQB1*06:21',
         'HLA-DQA1*01:08-HLA-DQB1*06:22', 'HLA-DQA1*01:08-HLA-DQB1*06:23',
         'HLA-DQA1*01:08-HLA-DQB1*06:24', 'HLA-DQA1*01:08-HLA-DQB1*06:25', 'HLA-DQA1*01:08-HLA-DQB1*06:27', 'HLA-DQA1*01:08-HLA-DQB1*06:28',
         'HLA-DQA1*01:08-HLA-DQB1*06:29', 'HLA-DQA1*01:08-HLA-DQB1*06:30',
         'HLA-DQA1*01:08-HLA-DQB1*06:31', 'HLA-DQA1*01:08-HLA-DQB1*06:32', 'HLA-DQA1*01:08-HLA-DQB1*06:33', 'HLA-DQA1*01:08-HLA-DQB1*06:34',
         'HLA-DQA1*01:08-HLA-DQB1*06:35', 'HLA-DQA1*01:08-HLA-DQB1*06:36',
         'HLA-DQA1*01:08-HLA-DQB1*06:37', 'HLA-DQA1*01:08-HLA-DQB1*06:38', 'HLA-DQA1*01:08-HLA-DQB1*06:39', 'HLA-DQA1*01:08-HLA-DQB1*06:40',
         'HLA-DQA1*01:08-HLA-DQB1*06:41', 'HLA-DQA1*01:08-HLA-DQB1*06:42',
         'HLA-DQA1*01:08-HLA-DQB1*06:43', 'HLA-DQA1*01:08-HLA-DQB1*06:44', 'HLA-DQA1*01:09-HLA-DQB1*02:01', 'HLA-DQA1*01:09-HLA-DQB1*02:02',
         'HLA-DQA1*01:09-HLA-DQB1*02:03', 'HLA-DQA1*01:09-HLA-DQB1*02:04',
         'HLA-DQA1*01:09-HLA-DQB1*02:05', 'HLA-DQA1*01:09-HLA-DQB1*02:06', 'HLA-DQA1*01:09-HLA-DQB1*03:01', 'HLA-DQA1*01:09-HLA-DQB1*03:02',
         'HLA-DQA1*01:09-HLA-DQB1*03:03', 'HLA-DQA1*01:09-HLA-DQB1*03:04',
         'HLA-DQA1*01:09-HLA-DQB1*03:05', 'HLA-DQA1*01:09-HLA-DQB1*03:06', 'HLA-DQA1*01:09-HLA-DQB1*03:07', 'HLA-DQA1*01:09-HLA-DQB1*03:08',
         'HLA-DQA1*01:09-HLA-DQB1*03:09', 'HLA-DQA1*01:09-HLA-DQB1*03:10',
         'HLA-DQA1*01:09-HLA-DQB1*03:11', 'HLA-DQA1*01:09-HLA-DQB1*03:12', 'HLA-DQA1*01:09-HLA-DQB1*03:13', 'HLA-DQA1*01:09-HLA-DQB1*03:14',
         'HLA-DQA1*01:09-HLA-DQB1*03:15', 'HLA-DQA1*01:09-HLA-DQB1*03:16',
         'HLA-DQA1*01:09-HLA-DQB1*03:17', 'HLA-DQA1*01:09-HLA-DQB1*03:18', 'HLA-DQA1*01:09-HLA-DQB1*03:19', 'HLA-DQA1*01:09-HLA-DQB1*03:20',
         'HLA-DQA1*01:09-HLA-DQB1*03:21', 'HLA-DQA1*01:09-HLA-DQB1*03:22',
         'HLA-DQA1*01:09-HLA-DQB1*03:23', 'HLA-DQA1*01:09-HLA-DQB1*03:24', 'HLA-DQA1*01:09-HLA-DQB1*03:25', 'HLA-DQA1*01:09-HLA-DQB1*03:26',
         'HLA-DQA1*01:09-HLA-DQB1*03:27', 'HLA-DQA1*01:09-HLA-DQB1*03:28',
         'HLA-DQA1*01:09-HLA-DQB1*03:29', 'HLA-DQA1*01:09-HLA-DQB1*03:30', 'HLA-DQA1*01:09-HLA-DQB1*03:31', 'HLA-DQA1*01:09-HLA-DQB1*03:32',
         'HLA-DQA1*01:09-HLA-DQB1*03:33', 'HLA-DQA1*01:09-HLA-DQB1*03:34',
         'HLA-DQA1*01:09-HLA-DQB1*03:35', 'HLA-DQA1*01:09-HLA-DQB1*03:36', 'HLA-DQA1*01:09-HLA-DQB1*03:37', 'HLA-DQA1*01:09-HLA-DQB1*03:38',
         'HLA-DQA1*01:09-HLA-DQB1*04:01', 'HLA-DQA1*01:09-HLA-DQB1*04:02',
         'HLA-DQA1*01:09-HLA-DQB1*04:03', 'HLA-DQA1*01:09-HLA-DQB1*04:04', 'HLA-DQA1*01:09-HLA-DQB1*04:05', 'HLA-DQA1*01:09-HLA-DQB1*04:06',
         'HLA-DQA1*01:09-HLA-DQB1*04:07', 'HLA-DQA1*01:09-HLA-DQB1*04:08',
         'HLA-DQA1*01:09-HLA-DQB1*05:01', 'HLA-DQA1*01:09-HLA-DQB1*05:02', 'HLA-DQA1*01:09-HLA-DQB1*05:03', 'HLA-DQA1*01:09-HLA-DQB1*05:05',
         'HLA-DQA1*01:09-HLA-DQB1*05:06', 'HLA-DQA1*01:09-HLA-DQB1*05:07',
         'HLA-DQA1*01:09-HLA-DQB1*05:08', 'HLA-DQA1*01:09-HLA-DQB1*05:09', 'HLA-DQA1*01:09-HLA-DQB1*05:10', 'HLA-DQA1*01:09-HLA-DQB1*05:11',
         'HLA-DQA1*01:09-HLA-DQB1*05:12', 'HLA-DQA1*01:09-HLA-DQB1*05:13',
         'HLA-DQA1*01:09-HLA-DQB1*05:14', 'HLA-DQA1*01:09-HLA-DQB1*06:01', 'HLA-DQA1*01:09-HLA-DQB1*06:02', 'HLA-DQA1*01:09-HLA-DQB1*06:03',
         'HLA-DQA1*01:09-HLA-DQB1*06:04', 'HLA-DQA1*01:09-HLA-DQB1*06:07',
         'HLA-DQA1*01:09-HLA-DQB1*06:08', 'HLA-DQA1*01:09-HLA-DQB1*06:09', 'HLA-DQA1*01:09-HLA-DQB1*06:10', 'HLA-DQA1*01:09-HLA-DQB1*06:11',
         'HLA-DQA1*01:09-HLA-DQB1*06:12', 'HLA-DQA1*01:09-HLA-DQB1*06:14',
         'HLA-DQA1*01:09-HLA-DQB1*06:15', 'HLA-DQA1*01:09-HLA-DQB1*06:16', 'HLA-DQA1*01:09-HLA-DQB1*06:17', 'HLA-DQA1*01:09-HLA-DQB1*06:18',
         'HLA-DQA1*01:09-HLA-DQB1*06:19', 'HLA-DQA1*01:09-HLA-DQB1*06:21',
         'HLA-DQA1*01:09-HLA-DQB1*06:22', 'HLA-DQA1*01:09-HLA-DQB1*06:23', 'HLA-DQA1*01:09-HLA-DQB1*06:24', 'HLA-DQA1*01:09-HLA-DQB1*06:25',
         'HLA-DQA1*01:09-HLA-DQB1*06:27', 'HLA-DQA1*01:09-HLA-DQB1*06:28',
         'HLA-DQA1*01:09-HLA-DQB1*06:29', 'HLA-DQA1*01:09-HLA-DQB1*06:30', 'HLA-DQA1*01:09-HLA-DQB1*06:31', 'HLA-DQA1*01:09-HLA-DQB1*06:32',
         'HLA-DQA1*01:09-HLA-DQB1*06:33', 'HLA-DQA1*01:09-HLA-DQB1*06:34',
         'HLA-DQA1*01:09-HLA-DQB1*06:35', 'HLA-DQA1*01:09-HLA-DQB1*06:36', 'HLA-DQA1*01:09-HLA-DQB1*06:37', 'HLA-DQA1*01:09-HLA-DQB1*06:38',
         'HLA-DQA1*01:09-HLA-DQB1*06:39', 'HLA-DQA1*01:09-HLA-DQB1*06:40',
         'HLA-DQA1*01:09-HLA-DQB1*06:41', 'HLA-DQA1*01:09-HLA-DQB1*06:42', 'HLA-DQA1*01:09-HLA-DQB1*06:43', 'HLA-DQA1*01:09-HLA-DQB1*06:44',
         'HLA-DQA1*02:01-HLA-DQB1*02:01', 'HLA-DQA1*02:01-HLA-DQB1*02:02',
         'HLA-DQA1*02:01-HLA-DQB1*02:03', 'HLA-DQA1*02:01-HLA-DQB1*02:04', 'HLA-DQA1*02:01-HLA-DQB1*02:05', 'HLA-DQA1*02:01-HLA-DQB1*02:06',
         'HLA-DQA1*02:01-HLA-DQB1*03:01', 'HLA-DQA1*02:01-HLA-DQB1*03:02',
         'HLA-DQA1*02:01-HLA-DQB1*03:03', 'HLA-DQA1*02:01-HLA-DQB1*03:04', 'HLA-DQA1*02:01-HLA-DQB1*03:05', 'HLA-DQA1*02:01-HLA-DQB1*03:06',
         'HLA-DQA1*02:01-HLA-DQB1*03:07', 'HLA-DQA1*02:01-HLA-DQB1*03:08',
         'HLA-DQA1*02:01-HLA-DQB1*03:09', 'HLA-DQA1*02:01-HLA-DQB1*03:10', 'HLA-DQA1*02:01-HLA-DQB1*03:11', 'HLA-DQA1*02:01-HLA-DQB1*03:12',
         'HLA-DQA1*02:01-HLA-DQB1*03:13', 'HLA-DQA1*02:01-HLA-DQB1*03:14',
         'HLA-DQA1*02:01-HLA-DQB1*03:15', 'HLA-DQA1*02:01-HLA-DQB1*03:16', 'HLA-DQA1*02:01-HLA-DQB1*03:17', 'HLA-DQA1*02:01-HLA-DQB1*03:18',
         'HLA-DQA1*02:01-HLA-DQB1*03:19', 'HLA-DQA1*02:01-HLA-DQB1*03:20',
         'HLA-DQA1*02:01-HLA-DQB1*03:21', 'HLA-DQA1*02:01-HLA-DQB1*03:22', 'HLA-DQA1*02:01-HLA-DQB1*03:23', 'HLA-DQA1*02:01-HLA-DQB1*03:24',
         'HLA-DQA1*02:01-HLA-DQB1*03:25', 'HLA-DQA1*02:01-HLA-DQB1*03:26',
         'HLA-DQA1*02:01-HLA-DQB1*03:27', 'HLA-DQA1*02:01-HLA-DQB1*03:28', 'HLA-DQA1*02:01-HLA-DQB1*03:29', 'HLA-DQA1*02:01-HLA-DQB1*03:30',
         'HLA-DQA1*02:01-HLA-DQB1*03:31', 'HLA-DQA1*02:01-HLA-DQB1*03:32',
         'HLA-DQA1*02:01-HLA-DQB1*03:33', 'HLA-DQA1*02:01-HLA-DQB1*03:34', 'HLA-DQA1*02:01-HLA-DQB1*03:35', 'HLA-DQA1*02:01-HLA-DQB1*03:36',
         'HLA-DQA1*02:01-HLA-DQB1*03:37', 'HLA-DQA1*02:01-HLA-DQB1*03:38',
         'HLA-DQA1*02:01-HLA-DQB1*04:01', 'HLA-DQA1*02:01-HLA-DQB1*04:02', 'HLA-DQA1*02:01-HLA-DQB1*04:03', 'HLA-DQA1*02:01-HLA-DQB1*04:04',
         'HLA-DQA1*02:01-HLA-DQB1*04:05', 'HLA-DQA1*02:01-HLA-DQB1*04:06',
         'HLA-DQA1*02:01-HLA-DQB1*04:07', 'HLA-DQA1*02:01-HLA-DQB1*04:08', 'HLA-DQA1*02:01-HLA-DQB1*05:01', 'HLA-DQA1*02:01-HLA-DQB1*05:02',
         'HLA-DQA1*02:01-HLA-DQB1*05:03', 'HLA-DQA1*02:01-HLA-DQB1*05:05',
         'HLA-DQA1*02:01-HLA-DQB1*05:06', 'HLA-DQA1*02:01-HLA-DQB1*05:07', 'HLA-DQA1*02:01-HLA-DQB1*05:08', 'HLA-DQA1*02:01-HLA-DQB1*05:09',
         'HLA-DQA1*02:01-HLA-DQB1*05:10', 'HLA-DQA1*02:01-HLA-DQB1*05:11',
         'HLA-DQA1*02:01-HLA-DQB1*05:12', 'HLA-DQA1*02:01-HLA-DQB1*05:13', 'HLA-DQA1*02:01-HLA-DQB1*05:14', 'HLA-DQA1*02:01-HLA-DQB1*06:01',
         'HLA-DQA1*02:01-HLA-DQB1*06:02', 'HLA-DQA1*02:01-HLA-DQB1*06:03',
         'HLA-DQA1*02:01-HLA-DQB1*06:04', 'HLA-DQA1*02:01-HLA-DQB1*06:07', 'HLA-DQA1*02:01-HLA-DQB1*06:08', 'HLA-DQA1*02:01-HLA-DQB1*06:09',
         'HLA-DQA1*02:01-HLA-DQB1*06:10', 'HLA-DQA1*02:01-HLA-DQB1*06:11',
         'HLA-DQA1*02:01-HLA-DQB1*06:12', 'HLA-DQA1*02:01-HLA-DQB1*06:14', 'HLA-DQA1*02:01-HLA-DQB1*06:15', 'HLA-DQA1*02:01-HLA-DQB1*06:16',
         'HLA-DQA1*02:01-HLA-DQB1*06:17', 'HLA-DQA1*02:01-HLA-DQB1*06:18',
         'HLA-DQA1*02:01-HLA-DQB1*06:19', 'HLA-DQA1*02:01-HLA-DQB1*06:21', 'HLA-DQA1*02:01-HLA-DQB1*06:22', 'HLA-DQA1*02:01-HLA-DQB1*06:23',
         'HLA-DQA1*02:01-HLA-DQB1*06:24', 'HLA-DQA1*02:01-HLA-DQB1*06:25',
         'HLA-DQA1*02:01-HLA-DQB1*06:27', 'HLA-DQA1*02:01-HLA-DQB1*06:28', 'HLA-DQA1*02:01-HLA-DQB1*06:29', 'HLA-DQA1*02:01-HLA-DQB1*06:30',
         'HLA-DQA1*02:01-HLA-DQB1*06:31', 'HLA-DQA1*02:01-HLA-DQB1*06:32',
         'HLA-DQA1*02:01-HLA-DQB1*06:33', 'HLA-DQA1*02:01-HLA-DQB1*06:34', 'HLA-DQA1*02:01-HLA-DQB1*06:35', 'HLA-DQA1*02:01-HLA-DQB1*06:36',
         'HLA-DQA1*02:01-HLA-DQB1*06:37', 'HLA-DQA1*02:01-HLA-DQB1*06:38',
         'HLA-DQA1*02:01-HLA-DQB1*06:39', 'HLA-DQA1*02:01-HLA-DQB1*06:40', 'HLA-DQA1*02:01-HLA-DQB1*06:41', 'HLA-DQA1*02:01-HLA-DQB1*06:42',
         'HLA-DQA1*02:01-HLA-DQB1*06:43', 'HLA-DQA1*02:01-HLA-DQB1*06:44',
         'HLA-DQA1*03:01-HLA-DQB1*02:01', 'HLA-DQA1*03:01-HLA-DQB1*02:02', 'HLA-DQA1*03:01-HLA-DQB1*02:03', 'HLA-DQA1*03:01-HLA-DQB1*02:04',
         'HLA-DQA1*03:01-HLA-DQB1*02:05', 'HLA-DQA1*03:01-HLA-DQB1*02:06',
         'HLA-DQA1*03:01-HLA-DQB1*03:01', 'HLA-DQA1*03:01-HLA-DQB1*03:02', 'HLA-DQA1*03:01-HLA-DQB1*03:03', 'HLA-DQA1*03:01-HLA-DQB1*03:04',
         'HLA-DQA1*03:01-HLA-DQB1*03:05', 'HLA-DQA1*03:01-HLA-DQB1*03:06',
         'HLA-DQA1*03:01-HLA-DQB1*03:07', 'HLA-DQA1*03:01-HLA-DQB1*03:08', 'HLA-DQA1*03:01-HLA-DQB1*03:09', 'HLA-DQA1*03:01-HLA-DQB1*03:10',
         'HLA-DQA1*03:01-HLA-DQB1*03:11', 'HLA-DQA1*03:01-HLA-DQB1*03:12',
         'HLA-DQA1*03:01-HLA-DQB1*03:13', 'HLA-DQA1*03:01-HLA-DQB1*03:14', 'HLA-DQA1*03:01-HLA-DQB1*03:15', 'HLA-DQA1*03:01-HLA-DQB1*03:16',
         'HLA-DQA1*03:01-HLA-DQB1*03:17', 'HLA-DQA1*03:01-HLA-DQB1*03:18',
         'HLA-DQA1*03:01-HLA-DQB1*03:19', 'HLA-DQA1*03:01-HLA-DQB1*03:20', 'HLA-DQA1*03:01-HLA-DQB1*03:21', 'HLA-DQA1*03:01-HLA-DQB1*03:22',
         'HLA-DQA1*03:01-HLA-DQB1*03:23', 'HLA-DQA1*03:01-HLA-DQB1*03:24',
         'HLA-DQA1*03:01-HLA-DQB1*03:25', 'HLA-DQA1*03:01-HLA-DQB1*03:26', 'HLA-DQA1*03:01-HLA-DQB1*03:27', 'HLA-DQA1*03:01-HLA-DQB1*03:28',
         'HLA-DQA1*03:01-HLA-DQB1*03:29', 'HLA-DQA1*03:01-HLA-DQB1*03:30',
         'HLA-DQA1*03:01-HLA-DQB1*03:31', 'HLA-DQA1*03:01-HLA-DQB1*03:32', 'HLA-DQA1*03:01-HLA-DQB1*03:33', 'HLA-DQA1*03:01-HLA-DQB1*03:34',
         'HLA-DQA1*03:01-HLA-DQB1*03:35', 'HLA-DQA1*03:01-HLA-DQB1*03:36',
         'HLA-DQA1*03:01-HLA-DQB1*03:37', 'HLA-DQA1*03:01-HLA-DQB1*03:38', 'HLA-DQA1*03:01-HLA-DQB1*04:01', 'HLA-DQA1*03:01-HLA-DQB1*04:02',
         'HLA-DQA1*03:01-HLA-DQB1*04:03', 'HLA-DQA1*03:01-HLA-DQB1*04:04',
         'HLA-DQA1*03:01-HLA-DQB1*04:05', 'HLA-DQA1*03:01-HLA-DQB1*04:06', 'HLA-DQA1*03:01-HLA-DQB1*04:07', 'HLA-DQA1*03:01-HLA-DQB1*04:08',
         'HLA-DQA1*03:01-HLA-DQB1*05:01', 'HLA-DQA1*03:01-HLA-DQB1*05:02',
         'HLA-DQA1*03:01-HLA-DQB1*05:03', 'HLA-DQA1*03:01-HLA-DQB1*05:05', 'HLA-DQA1*03:01-HLA-DQB1*05:06', 'HLA-DQA1*03:01-HLA-DQB1*05:07',
         'HLA-DQA1*03:01-HLA-DQB1*05:08', 'HLA-DQA1*03:01-HLA-DQB1*05:09',
         'HLA-DQA1*03:01-HLA-DQB1*05:10', 'HLA-DQA1*03:01-HLA-DQB1*05:11', 'HLA-DQA1*03:01-HLA-DQB1*05:12', 'HLA-DQA1*03:01-HLA-DQB1*05:13',
         'HLA-DQA1*03:01-HLA-DQB1*05:14', 'HLA-DQA1*03:01-HLA-DQB1*06:01',
         'HLA-DQA1*03:01-HLA-DQB1*06:02', 'HLA-DQA1*03:01-HLA-DQB1*06:03', 'HLA-DQA1*03:01-HLA-DQB1*06:04', 'HLA-DQA1*03:01-HLA-DQB1*06:07',
         'HLA-DQA1*03:01-HLA-DQB1*06:08', 'HLA-DQA1*03:01-HLA-DQB1*06:09',
         'HLA-DQA1*03:01-HLA-DQB1*06:10', 'HLA-DQA1*03:01-HLA-DQB1*06:11', 'HLA-DQA1*03:01-HLA-DQB1*06:12', 'HLA-DQA1*03:01-HLA-DQB1*06:14',
         'HLA-DQA1*03:01-HLA-DQB1*06:15', 'HLA-DQA1*03:01-HLA-DQB1*06:16',
         'HLA-DQA1*03:01-HLA-DQB1*06:17', 'HLA-DQA1*03:01-HLA-DQB1*06:18', 'HLA-DQA1*03:01-HLA-DQB1*06:19', 'HLA-DQA1*03:01-HLA-DQB1*06:21',
         'HLA-DQA1*03:01-HLA-DQB1*06:22', 'HLA-DQA1*03:01-HLA-DQB1*06:23',
         'HLA-DQA1*03:01-HLA-DQB1*06:24', 'HLA-DQA1*03:01-HLA-DQB1*06:25', 'HLA-DQA1*03:01-HLA-DQB1*06:27', 'HLA-DQA1*03:01-HLA-DQB1*06:28',
         'HLA-DQA1*03:01-HLA-DQB1*06:29', 'HLA-DQA1*03:01-HLA-DQB1*06:30',
         'HLA-DQA1*03:01-HLA-DQB1*06:31', 'HLA-DQA1*03:01-HLA-DQB1*06:32', 'HLA-DQA1*03:01-HLA-DQB1*06:33', 'HLA-DQA1*03:01-HLA-DQB1*06:34',
         'HLA-DQA1*03:01-HLA-DQB1*06:35', 'HLA-DQA1*03:01-HLA-DQB1*06:36',
         'HLA-DQA1*03:01-HLA-DQB1*06:37', 'HLA-DQA1*03:01-HLA-DQB1*06:38', 'HLA-DQA1*03:01-HLA-DQB1*06:39', 'HLA-DQA1*03:01-HLA-DQB1*06:40',
         'HLA-DQA1*03:01-HLA-DQB1*06:41', 'HLA-DQA1*03:01-HLA-DQB1*06:42',
         'HLA-DQA1*03:01-HLA-DQB1*06:43', 'HLA-DQA1*03:01-HLA-DQB1*06:44', 'HLA-DQA1*03:02-HLA-DQB1*02:01', 'HLA-DQA1*03:02-HLA-DQB1*02:02',
         'HLA-DQA1*03:02-HLA-DQB1*02:03', 'HLA-DQA1*03:02-HLA-DQB1*02:04',
         'HLA-DQA1*03:02-HLA-DQB1*02:05', 'HLA-DQA1*03:02-HLA-DQB1*02:06', 'HLA-DQA1*03:02-HLA-DQB1*03:01', 'HLA-DQA1*03:02-HLA-DQB1*03:02',
         'HLA-DQA1*03:02-HLA-DQB1*03:03', 'HLA-DQA1*03:02-HLA-DQB1*03:04',
         'HLA-DQA1*03:02-HLA-DQB1*03:05', 'HLA-DQA1*03:02-HLA-DQB1*03:06', 'HLA-DQA1*03:02-HLA-DQB1*03:07', 'HLA-DQA1*03:02-HLA-DQB1*03:08',
         'HLA-DQA1*03:02-HLA-DQB1*03:09', 'HLA-DQA1*03:02-HLA-DQB1*03:10',
         'HLA-DQA1*03:02-HLA-DQB1*03:11', 'HLA-DQA1*03:02-HLA-DQB1*03:12', 'HLA-DQA1*03:02-HLA-DQB1*03:13', 'HLA-DQA1*03:02-HLA-DQB1*03:14',
         'HLA-DQA1*03:02-HLA-DQB1*03:15', 'HLA-DQA1*03:02-HLA-DQB1*03:16',
         'HLA-DQA1*03:02-HLA-DQB1*03:17', 'HLA-DQA1*03:02-HLA-DQB1*03:18', 'HLA-DQA1*03:02-HLA-DQB1*03:19', 'HLA-DQA1*03:02-HLA-DQB1*03:20',
         'HLA-DQA1*03:02-HLA-DQB1*03:21', 'HLA-DQA1*03:02-HLA-DQB1*03:22',
         'HLA-DQA1*03:02-HLA-DQB1*03:23', 'HLA-DQA1*03:02-HLA-DQB1*03:24', 'HLA-DQA1*03:02-HLA-DQB1*03:25', 'HLA-DQA1*03:02-HLA-DQB1*03:26',
         'HLA-DQA1*03:02-HLA-DQB1*03:27', 'HLA-DQA1*03:02-HLA-DQB1*03:28',
         'HLA-DQA1*03:02-HLA-DQB1*03:29', 'HLA-DQA1*03:02-HLA-DQB1*03:30', 'HLA-DQA1*03:02-HLA-DQB1*03:31', 'HLA-DQA1*03:02-HLA-DQB1*03:32',
         'HLA-DQA1*03:02-HLA-DQB1*03:33', 'HLA-DQA1*03:02-HLA-DQB1*03:34',
         'HLA-DQA1*03:02-HLA-DQB1*03:35', 'HLA-DQA1*03:02-HLA-DQB1*03:36', 'HLA-DQA1*03:02-HLA-DQB1*03:37', 'HLA-DQA1*03:02-HLA-DQB1*03:38',
         'HLA-DQA1*03:02-HLA-DQB1*04:01', 'HLA-DQA1*03:02-HLA-DQB1*04:02',
         'HLA-DQA1*03:02-HLA-DQB1*04:03', 'HLA-DQA1*03:02-HLA-DQB1*04:04', 'HLA-DQA1*03:02-HLA-DQB1*04:05', 'HLA-DQA1*03:02-HLA-DQB1*04:06',
         'HLA-DQA1*03:02-HLA-DQB1*04:07', 'HLA-DQA1*03:02-HLA-DQB1*04:08',
         'HLA-DQA1*03:02-HLA-DQB1*05:01', 'HLA-DQA1*03:02-HLA-DQB1*05:02', 'HLA-DQA1*03:02-HLA-DQB1*05:03', 'HLA-DQA1*03:02-HLA-DQB1*05:05',
         'HLA-DQA1*03:02-HLA-DQB1*05:06', 'HLA-DQA1*03:02-HLA-DQB1*05:07',
         'HLA-DQA1*03:02-HLA-DQB1*05:08', 'HLA-DQA1*03:02-HLA-DQB1*05:09', 'HLA-DQA1*03:02-HLA-DQB1*05:10', 'HLA-DQA1*03:02-HLA-DQB1*05:11',
         'HLA-DQA1*03:02-HLA-DQB1*05:12', 'HLA-DQA1*03:02-HLA-DQB1*05:13',
         'HLA-DQA1*03:02-HLA-DQB1*05:14', 'HLA-DQA1*03:02-HLA-DQB1*06:01', 'HLA-DQA1*03:02-HLA-DQB1*06:02', 'HLA-DQA1*03:02-HLA-DQB1*06:03',
         'HLA-DQA1*03:02-HLA-DQB1*06:04', 'HLA-DQA1*03:02-HLA-DQB1*06:07',
         'HLA-DQA1*03:02-HLA-DQB1*06:08', 'HLA-DQA1*03:02-HLA-DQB1*06:09', 'HLA-DQA1*03:02-HLA-DQB1*06:10', 'HLA-DQA1*03:02-HLA-DQB1*06:11',
         'HLA-DQA1*03:02-HLA-DQB1*06:12', 'HLA-DQA1*03:02-HLA-DQB1*06:14',
         'HLA-DQA1*03:02-HLA-DQB1*06:15', 'HLA-DQA1*03:02-HLA-DQB1*06:16', 'HLA-DQA1*03:02-HLA-DQB1*06:17', 'HLA-DQA1*03:02-HLA-DQB1*06:18',
         'HLA-DQA1*03:02-HLA-DQB1*06:19', 'HLA-DQA1*03:02-HLA-DQB1*06:21',
         'HLA-DQA1*03:02-HLA-DQB1*06:22', 'HLA-DQA1*03:02-HLA-DQB1*06:23', 'HLA-DQA1*03:02-HLA-DQB1*06:24', 'HLA-DQA1*03:02-HLA-DQB1*06:25',
         'HLA-DQA1*03:02-HLA-DQB1*06:27', 'HLA-DQA1*03:02-HLA-DQB1*06:28',
         'HLA-DQA1*03:02-HLA-DQB1*06:29', 'HLA-DQA1*03:02-HLA-DQB1*06:30', 'HLA-DQA1*03:02-HLA-DQB1*06:31', 'HLA-DQA1*03:02-HLA-DQB1*06:32',
         'HLA-DQA1*03:02-HLA-DQB1*06:33', 'HLA-DQA1*03:02-HLA-DQB1*06:34',
         'HLA-DQA1*03:02-HLA-DQB1*06:35', 'HLA-DQA1*03:02-HLA-DQB1*06:36', 'HLA-DQA1*03:02-HLA-DQB1*06:37', 'HLA-DQA1*03:02-HLA-DQB1*06:38',
         'HLA-DQA1*03:02-HLA-DQB1*06:39', 'HLA-DQA1*03:02-HLA-DQB1*06:40',
         'HLA-DQA1*03:02-HLA-DQB1*06:41', 'HLA-DQA1*03:02-HLA-DQB1*06:42', 'HLA-DQA1*03:02-HLA-DQB1*06:43', 'HLA-DQA1*03:02-HLA-DQB1*06:44',
         'HLA-DQA1*03:03-HLA-DQB1*02:01', 'HLA-DQA1*03:03-HLA-DQB1*02:02',
         'HLA-DQA1*03:03-HLA-DQB1*02:03', 'HLA-DQA1*03:03-HLA-DQB1*02:04', 'HLA-DQA1*03:03-HLA-DQB1*02:05', 'HLA-DQA1*03:03-HLA-DQB1*02:06',
         'HLA-DQA1*03:03-HLA-DQB1*03:01', 'HLA-DQA1*03:03-HLA-DQB1*03:02',
         'HLA-DQA1*03:03-HLA-DQB1*03:03', 'HLA-DQA1*03:03-HLA-DQB1*03:04', 'HLA-DQA1*03:03-HLA-DQB1*03:05', 'HLA-DQA1*03:03-HLA-DQB1*03:06',
         'HLA-DQA1*03:03-HLA-DQB1*03:07', 'HLA-DQA1*03:03-HLA-DQB1*03:08',
         'HLA-DQA1*03:03-HLA-DQB1*03:09', 'HLA-DQA1*03:03-HLA-DQB1*03:10', 'HLA-DQA1*03:03-HLA-DQB1*03:11', 'HLA-DQA1*03:03-HLA-DQB1*03:12',
         'HLA-DQA1*03:03-HLA-DQB1*03:13', 'HLA-DQA1*03:03-HLA-DQB1*03:14',
         'HLA-DQA1*03:03-HLA-DQB1*03:15', 'HLA-DQA1*03:03-HLA-DQB1*03:16', 'HLA-DQA1*03:03-HLA-DQB1*03:17', 'HLA-DQA1*03:03-HLA-DQB1*03:18',
         'HLA-DQA1*03:03-HLA-DQB1*03:19', 'HLA-DQA1*03:03-HLA-DQB1*03:20',
         'HLA-DQA1*03:03-HLA-DQB1*03:21', 'HLA-DQA1*03:03-HLA-DQB1*03:22', 'HLA-DQA1*03:03-HLA-DQB1*03:23', 'HLA-DQA1*03:03-HLA-DQB1*03:24',
         'HLA-DQA1*03:03-HLA-DQB1*03:25', 'HLA-DQA1*03:03-HLA-DQB1*03:26',
         'HLA-DQA1*03:03-HLA-DQB1*03:27', 'HLA-DQA1*03:03-HLA-DQB1*03:28', 'HLA-DQA1*03:03-HLA-DQB1*03:29', 'HLA-DQA1*03:03-HLA-DQB1*03:30',
         'HLA-DQA1*03:03-HLA-DQB1*03:31', 'HLA-DQA1*03:03-HLA-DQB1*03:32',
         'HLA-DQA1*03:03-HLA-DQB1*03:33', 'HLA-DQA1*03:03-HLA-DQB1*03:34', 'HLA-DQA1*03:03-HLA-DQB1*03:35', 'HLA-DQA1*03:03-HLA-DQB1*03:36',
         'HLA-DQA1*03:03-HLA-DQB1*03:37', 'HLA-DQA1*03:03-HLA-DQB1*03:38',
         'HLA-DQA1*03:03-HLA-DQB1*04:01', 'HLA-DQA1*03:03-HLA-DQB1*04:02', 'HLA-DQA1*03:03-HLA-DQB1*04:03', 'HLA-DQA1*03:03-HLA-DQB1*04:04',
         'HLA-DQA1*03:03-HLA-DQB1*04:05', 'HLA-DQA1*03:03-HLA-DQB1*04:06',
         'HLA-DQA1*03:03-HLA-DQB1*04:07', 'HLA-DQA1*03:03-HLA-DQB1*04:08', 'HLA-DQA1*03:03-HLA-DQB1*05:01', 'HLA-DQA1*03:03-HLA-DQB1*05:02',
         'HLA-DQA1*03:03-HLA-DQB1*05:03', 'HLA-DQA1*03:03-HLA-DQB1*05:05',
         'HLA-DQA1*03:03-HLA-DQB1*05:06', 'HLA-DQA1*03:03-HLA-DQB1*05:07', 'HLA-DQA1*03:03-HLA-DQB1*05:08', 'HLA-DQA1*03:03-HLA-DQB1*05:09',
         'HLA-DQA1*03:03-HLA-DQB1*05:10', 'HLA-DQA1*03:03-HLA-DQB1*05:11',
         'HLA-DQA1*03:03-HLA-DQB1*05:12', 'HLA-DQA1*03:03-HLA-DQB1*05:13', 'HLA-DQA1*03:03-HLA-DQB1*05:14', 'HLA-DQA1*03:03-HLA-DQB1*06:01',
         'HLA-DQA1*03:03-HLA-DQB1*06:02', 'HLA-DQA1*03:03-HLA-DQB1*06:03',
         'HLA-DQA1*03:03-HLA-DQB1*06:04', 'HLA-DQA1*03:03-HLA-DQB1*06:07', 'HLA-DQA1*03:03-HLA-DQB1*06:08', 'HLA-DQA1*03:03-HLA-DQB1*06:09',
         'HLA-DQA1*03:03-HLA-DQB1*06:10', 'HLA-DQA1*03:03-HLA-DQB1*06:11',
         'HLA-DQA1*03:03-HLA-DQB1*06:12', 'HLA-DQA1*03:03-HLA-DQB1*06:14', 'HLA-DQA1*03:03-HLA-DQB1*06:15', 'HLA-DQA1*03:03-HLA-DQB1*06:16',
         'HLA-DQA1*03:03-HLA-DQB1*06:17', 'HLA-DQA1*03:03-HLA-DQB1*06:18',
         'HLA-DQA1*03:03-HLA-DQB1*06:19', 'HLA-DQA1*03:03-HLA-DQB1*06:21', 'HLA-DQA1*03:03-HLA-DQB1*06:22', 'HLA-DQA1*03:03-HLA-DQB1*06:23',
         'HLA-DQA1*03:03-HLA-DQB1*06:24', 'HLA-DQA1*03:03-HLA-DQB1*06:25',
         'HLA-DQA1*03:03-HLA-DQB1*06:27', 'HLA-DQA1*03:03-HLA-DQB1*06:28', 'HLA-DQA1*03:03-HLA-DQB1*06:29', 'HLA-DQA1*03:03-HLA-DQB1*06:30',
         'HLA-DQA1*03:03-HLA-DQB1*06:31', 'HLA-DQA1*03:03-HLA-DQB1*06:32',
         'HLA-DQA1*03:03-HLA-DQB1*06:33', 'HLA-DQA1*03:03-HLA-DQB1*06:34', 'HLA-DQA1*03:03-HLA-DQB1*06:35', 'HLA-DQA1*03:03-HLA-DQB1*06:36',
         'HLA-DQA1*03:03-HLA-DQB1*06:37', 'HLA-DQA1*03:03-HLA-DQB1*06:38',
         'HLA-DQA1*03:03-HLA-DQB1*06:39', 'HLA-DQA1*03:03-HLA-DQB1*06:40', 'HLA-DQA1*03:03-HLA-DQB1*06:41', 'HLA-DQA1*03:03-HLA-DQB1*06:42',
         'HLA-DQA1*03:03-HLA-DQB1*06:43', 'HLA-DQA1*03:03-HLA-DQB1*06:44',
         'HLA-DQA1*04:01-HLA-DQB1*02:01', 'HLA-DQA1*04:01-HLA-DQB1*02:02', 'HLA-DQA1*04:01-HLA-DQB1*02:03', 'HLA-DQA1*04:01-HLA-DQB1*02:04',
         'HLA-DQA1*04:01-HLA-DQB1*02:05', 'HLA-DQA1*04:01-HLA-DQB1*02:06',
         'HLA-DQA1*04:01-HLA-DQB1*03:01', 'HLA-DQA1*04:01-HLA-DQB1*03:02', 'HLA-DQA1*04:01-HLA-DQB1*03:03', 'HLA-DQA1*04:01-HLA-DQB1*03:04',
         'HLA-DQA1*04:01-HLA-DQB1*03:05', 'HLA-DQA1*04:01-HLA-DQB1*03:06',
         'HLA-DQA1*04:01-HLA-DQB1*03:07', 'HLA-DQA1*04:01-HLA-DQB1*03:08', 'HLA-DQA1*04:01-HLA-DQB1*03:09', 'HLA-DQA1*04:01-HLA-DQB1*03:10',
         'HLA-DQA1*04:01-HLA-DQB1*03:11', 'HLA-DQA1*04:01-HLA-DQB1*03:12',
         'HLA-DQA1*04:01-HLA-DQB1*03:13', 'HLA-DQA1*04:01-HLA-DQB1*03:14', 'HLA-DQA1*04:01-HLA-DQB1*03:15', 'HLA-DQA1*04:01-HLA-DQB1*03:16',
         'HLA-DQA1*04:01-HLA-DQB1*03:17', 'HLA-DQA1*04:01-HLA-DQB1*03:18',
         'HLA-DQA1*04:01-HLA-DQB1*03:19', 'HLA-DQA1*04:01-HLA-DQB1*03:20', 'HLA-DQA1*04:01-HLA-DQB1*03:21', 'HLA-DQA1*04:01-HLA-DQB1*03:22',
         'HLA-DQA1*04:01-HLA-DQB1*03:23', 'HLA-DQA1*04:01-HLA-DQB1*03:24',
         'HLA-DQA1*04:01-HLA-DQB1*03:25', 'HLA-DQA1*04:01-HLA-DQB1*03:26', 'HLA-DQA1*04:01-HLA-DQB1*03:27', 'HLA-DQA1*04:01-HLA-DQB1*03:28',
         'HLA-DQA1*04:01-HLA-DQB1*03:29', 'HLA-DQA1*04:01-HLA-DQB1*03:30',
         'HLA-DQA1*04:01-HLA-DQB1*03:31', 'HLA-DQA1*04:01-HLA-DQB1*03:32', 'HLA-DQA1*04:01-HLA-DQB1*03:33', 'HLA-DQA1*04:01-HLA-DQB1*03:34',
         'HLA-DQA1*04:01-HLA-DQB1*03:35', 'HLA-DQA1*04:01-HLA-DQB1*03:36',
         'HLA-DQA1*04:01-HLA-DQB1*03:37', 'HLA-DQA1*04:01-HLA-DQB1*03:38', 'HLA-DQA1*04:01-HLA-DQB1*04:01', 'HLA-DQA1*04:01-HLA-DQB1*04:02',
         'HLA-DQA1*04:01-HLA-DQB1*04:03', 'HLA-DQA1*04:01-HLA-DQB1*04:04',
         'HLA-DQA1*04:01-HLA-DQB1*04:05', 'HLA-DQA1*04:01-HLA-DQB1*04:06', 'HLA-DQA1*04:01-HLA-DQB1*04:07', 'HLA-DQA1*04:01-HLA-DQB1*04:08',
         'HLA-DQA1*04:01-HLA-DQB1*05:01', 'HLA-DQA1*04:01-HLA-DQB1*05:02',
         'HLA-DQA1*04:01-HLA-DQB1*05:03', 'HLA-DQA1*04:01-HLA-DQB1*05:05', 'HLA-DQA1*04:01-HLA-DQB1*05:06', 'HLA-DQA1*04:01-HLA-DQB1*05:07',
         'HLA-DQA1*04:01-HLA-DQB1*05:08', 'HLA-DQA1*04:01-HLA-DQB1*05:09',
         'HLA-DQA1*04:01-HLA-DQB1*05:10', 'HLA-DQA1*04:01-HLA-DQB1*05:11', 'HLA-DQA1*04:01-HLA-DQB1*05:12', 'HLA-DQA1*04:01-HLA-DQB1*05:13',
         'HLA-DQA1*04:01-HLA-DQB1*05:14', 'HLA-DQA1*04:01-HLA-DQB1*06:01',
         'HLA-DQA1*04:01-HLA-DQB1*06:02', 'HLA-DQA1*04:01-HLA-DQB1*06:03', 'HLA-DQA1*04:01-HLA-DQB1*06:04', 'HLA-DQA1*04:01-HLA-DQB1*06:07',
         'HLA-DQA1*04:01-HLA-DQB1*06:08', 'HLA-DQA1*04:01-HLA-DQB1*06:09',
         'HLA-DQA1*04:01-HLA-DQB1*06:10', 'HLA-DQA1*04:01-HLA-DQB1*06:11', 'HLA-DQA1*04:01-HLA-DQB1*06:12', 'HLA-DQA1*04:01-HLA-DQB1*06:14',
         'HLA-DQA1*04:01-HLA-DQB1*06:15', 'HLA-DQA1*04:01-HLA-DQB1*06:16',
         'HLA-DQA1*04:01-HLA-DQB1*06:17', 'HLA-DQA1*04:01-HLA-DQB1*06:18', 'HLA-DQA1*04:01-HLA-DQB1*06:19', 'HLA-DQA1*04:01-HLA-DQB1*06:21',
         'HLA-DQA1*04:01-HLA-DQB1*06:22', 'HLA-DQA1*04:01-HLA-DQB1*06:23',
         'HLA-DQA1*04:01-HLA-DQB1*06:24', 'HLA-DQA1*04:01-HLA-DQB1*06:25', 'HLA-DQA1*04:01-HLA-DQB1*06:27', 'HLA-DQA1*04:01-HLA-DQB1*06:28',
         'HLA-DQA1*04:01-HLA-DQB1*06:29', 'HLA-DQA1*04:01-HLA-DQB1*06:30',
         'HLA-DQA1*04:01-HLA-DQB1*06:31', 'HLA-DQA1*04:01-HLA-DQB1*06:32', 'HLA-DQA1*04:01-HLA-DQB1*06:33', 'HLA-DQA1*04:01-HLA-DQB1*06:34',
         'HLA-DQA1*04:01-HLA-DQB1*06:35', 'HLA-DQA1*04:01-HLA-DQB1*06:36',
         'HLA-DQA1*04:01-HLA-DQB1*06:37', 'HLA-DQA1*04:01-HLA-DQB1*06:38', 'HLA-DQA1*04:01-HLA-DQB1*06:39', 'HLA-DQA1*04:01-HLA-DQB1*06:40',
         'HLA-DQA1*04:01-HLA-DQB1*06:41', 'HLA-DQA1*04:01-HLA-DQB1*06:42',
         'HLA-DQA1*04:01-HLA-DQB1*06:43', 'HLA-DQA1*04:01-HLA-DQB1*06:44', 'HLA-DQA1*04:02-HLA-DQB1*02:01', 'HLA-DQA1*04:02-HLA-DQB1*02:02',
         'HLA-DQA1*04:02-HLA-DQB1*02:03', 'HLA-DQA1*04:02-HLA-DQB1*02:04',
         'HLA-DQA1*04:02-HLA-DQB1*02:05', 'HLA-DQA1*04:02-HLA-DQB1*02:06', 'HLA-DQA1*04:02-HLA-DQB1*03:01', 'HLA-DQA1*04:02-HLA-DQB1*03:02',
         'HLA-DQA1*04:02-HLA-DQB1*03:03', 'HLA-DQA1*04:02-HLA-DQB1*03:04',
         'HLA-DQA1*04:02-HLA-DQB1*03:05', 'HLA-DQA1*04:02-HLA-DQB1*03:06', 'HLA-DQA1*04:02-HLA-DQB1*03:07', 'HLA-DQA1*04:02-HLA-DQB1*03:08',
         'HLA-DQA1*04:02-HLA-DQB1*03:09', 'HLA-DQA1*04:02-HLA-DQB1*03:10',
         'HLA-DQA1*04:02-HLA-DQB1*03:11', 'HLA-DQA1*04:02-HLA-DQB1*03:12', 'HLA-DQA1*04:02-HLA-DQB1*03:13', 'HLA-DQA1*04:02-HLA-DQB1*03:14',
         'HLA-DQA1*04:02-HLA-DQB1*03:15', 'HLA-DQA1*04:02-HLA-DQB1*03:16',
         'HLA-DQA1*04:02-HLA-DQB1*03:17', 'HLA-DQA1*04:02-HLA-DQB1*03:18', 'HLA-DQA1*04:02-HLA-DQB1*03:19', 'HLA-DQA1*04:02-HLA-DQB1*03:20',
         'HLA-DQA1*04:02-HLA-DQB1*03:21', 'HLA-DQA1*04:02-HLA-DQB1*03:22',
         'HLA-DQA1*04:02-HLA-DQB1*03:23', 'HLA-DQA1*04:02-HLA-DQB1*03:24', 'HLA-DQA1*04:02-HLA-DQB1*03:25', 'HLA-DQA1*04:02-HLA-DQB1*03:26',
         'HLA-DQA1*04:02-HLA-DQB1*03:27', 'HLA-DQA1*04:02-HLA-DQB1*03:28',
         'HLA-DQA1*04:02-HLA-DQB1*03:29', 'HLA-DQA1*04:02-HLA-DQB1*03:30', 'HLA-DQA1*04:02-HLA-DQB1*03:31', 'HLA-DQA1*04:02-HLA-DQB1*03:32',
         'HLA-DQA1*04:02-HLA-DQB1*03:33', 'HLA-DQA1*04:02-HLA-DQB1*03:34',
         'HLA-DQA1*04:02-HLA-DQB1*03:35', 'HLA-DQA1*04:02-HLA-DQB1*03:36', 'HLA-DQA1*04:02-HLA-DQB1*03:37', 'HLA-DQA1*04:02-HLA-DQB1*03:38',
         'HLA-DQA1*04:02-HLA-DQB1*04:01', 'HLA-DQA1*04:02-HLA-DQB1*04:02',
         'HLA-DQA1*04:02-HLA-DQB1*04:03', 'HLA-DQA1*04:02-HLA-DQB1*04:04', 'HLA-DQA1*04:02-HLA-DQB1*04:05', 'HLA-DQA1*04:02-HLA-DQB1*04:06',
         'HLA-DQA1*04:02-HLA-DQB1*04:07', 'HLA-DQA1*04:02-HLA-DQB1*04:08',
         'HLA-DQA1*04:02-HLA-DQB1*05:01', 'HLA-DQA1*04:02-HLA-DQB1*05:02', 'HLA-DQA1*04:02-HLA-DQB1*05:03', 'HLA-DQA1*04:02-HLA-DQB1*05:05',
         'HLA-DQA1*04:02-HLA-DQB1*05:06', 'HLA-DQA1*04:02-HLA-DQB1*05:07',
         'HLA-DQA1*04:02-HLA-DQB1*05:08', 'HLA-DQA1*04:02-HLA-DQB1*05:09', 'HLA-DQA1*04:02-HLA-DQB1*05:10', 'HLA-DQA1*04:02-HLA-DQB1*05:11',
         'HLA-DQA1*04:02-HLA-DQB1*05:12', 'HLA-DQA1*04:02-HLA-DQB1*05:13',
         'HLA-DQA1*04:02-HLA-DQB1*05:14', 'HLA-DQA1*04:02-HLA-DQB1*06:01', 'HLA-DQA1*04:02-HLA-DQB1*06:02', 'HLA-DQA1*04:02-HLA-DQB1*06:03',
         'HLA-DQA1*04:02-HLA-DQB1*06:04', 'HLA-DQA1*04:02-HLA-DQB1*06:07',
         'HLA-DQA1*04:02-HLA-DQB1*06:08', 'HLA-DQA1*04:02-HLA-DQB1*06:09', 'HLA-DQA1*04:02-HLA-DQB1*06:10', 'HLA-DQA1*04:02-HLA-DQB1*06:11',
         'HLA-DQA1*04:02-HLA-DQB1*06:12', 'HLA-DQA1*04:02-HLA-DQB1*06:14',
         'HLA-DQA1*04:02-HLA-DQB1*06:15', 'HLA-DQA1*04:02-HLA-DQB1*06:16', 'HLA-DQA1*04:02-HLA-DQB1*06:17', 'HLA-DQA1*04:02-HLA-DQB1*06:18',
         'HLA-DQA1*04:02-HLA-DQB1*06:19', 'HLA-DQA1*04:02-HLA-DQB1*06:21',
         'HLA-DQA1*04:02-HLA-DQB1*06:22', 'HLA-DQA1*04:02-HLA-DQB1*06:23', 'HLA-DQA1*04:02-HLA-DQB1*06:24', 'HLA-DQA1*04:02-HLA-DQB1*06:25',
         'HLA-DQA1*04:02-HLA-DQB1*06:27', 'HLA-DQA1*04:02-HLA-DQB1*06:28',
         'HLA-DQA1*04:02-HLA-DQB1*06:29', 'HLA-DQA1*04:02-HLA-DQB1*06:30', 'HLA-DQA1*04:02-HLA-DQB1*06:31', 'HLA-DQA1*04:02-HLA-DQB1*06:32',
         'HLA-DQA1*04:02-HLA-DQB1*06:33', 'HLA-DQA1*04:02-HLA-DQB1*06:34',
         'HLA-DQA1*04:02-HLA-DQB1*06:35', 'HLA-DQA1*04:02-HLA-DQB1*06:36', 'HLA-DQA1*04:02-HLA-DQB1*06:37', 'HLA-DQA1*04:02-HLA-DQB1*06:38',
         'HLA-DQA1*04:02-HLA-DQB1*06:39', 'HLA-DQA1*04:02-HLA-DQB1*06:40',
         'HLA-DQA1*04:02-HLA-DQB1*06:41', 'HLA-DQA1*04:02-HLA-DQB1*06:42', 'HLA-DQA1*04:02-HLA-DQB1*06:43', 'HLA-DQA1*04:02-HLA-DQB1*06:44',
         'HLA-DQA1*04:04-HLA-DQB1*02:01', 'HLA-DQA1*04:04-HLA-DQB1*02:02',
         'HLA-DQA1*04:04-HLA-DQB1*02:03', 'HLA-DQA1*04:04-HLA-DQB1*02:04', 'HLA-DQA1*04:04-HLA-DQB1*02:05', 'HLA-DQA1*04:04-HLA-DQB1*02:06',
         'HLA-DQA1*04:04-HLA-DQB1*03:01', 'HLA-DQA1*04:04-HLA-DQB1*03:02',
         'HLA-DQA1*04:04-HLA-DQB1*03:03', 'HLA-DQA1*04:04-HLA-DQB1*03:04', 'HLA-DQA1*04:04-HLA-DQB1*03:05', 'HLA-DQA1*04:04-HLA-DQB1*03:06',
         'HLA-DQA1*04:04-HLA-DQB1*03:07', 'HLA-DQA1*04:04-HLA-DQB1*03:08',
         'HLA-DQA1*04:04-HLA-DQB1*03:09', 'HLA-DQA1*04:04-HLA-DQB1*03:10', 'HLA-DQA1*04:04-HLA-DQB1*03:11', 'HLA-DQA1*04:04-HLA-DQB1*03:12',
         'HLA-DQA1*04:04-HLA-DQB1*03:13', 'HLA-DQA1*04:04-HLA-DQB1*03:14',
         'HLA-DQA1*04:04-HLA-DQB1*03:15', 'HLA-DQA1*04:04-HLA-DQB1*03:16', 'HLA-DQA1*04:04-HLA-DQB1*03:17', 'HLA-DQA1*04:04-HLA-DQB1*03:18',
         'HLA-DQA1*04:04-HLA-DQB1*03:19', 'HLA-DQA1*04:04-HLA-DQB1*03:20',
         'HLA-DQA1*04:04-HLA-DQB1*03:21', 'HLA-DQA1*04:04-HLA-DQB1*03:22', 'HLA-DQA1*04:04-HLA-DQB1*03:23', 'HLA-DQA1*04:04-HLA-DQB1*03:24',
         'HLA-DQA1*04:04-HLA-DQB1*03:25', 'HLA-DQA1*04:04-HLA-DQB1*03:26',
         'HLA-DQA1*04:04-HLA-DQB1*03:27', 'HLA-DQA1*04:04-HLA-DQB1*03:28', 'HLA-DQA1*04:04-HLA-DQB1*03:29', 'HLA-DQA1*04:04-HLA-DQB1*03:30',
         'HLA-DQA1*04:04-HLA-DQB1*03:31', 'HLA-DQA1*04:04-HLA-DQB1*03:32',
         'HLA-DQA1*04:04-HLA-DQB1*03:33', 'HLA-DQA1*04:04-HLA-DQB1*03:34', 'HLA-DQA1*04:04-HLA-DQB1*03:35', 'HLA-DQA1*04:04-HLA-DQB1*03:36',
         'HLA-DQA1*04:04-HLA-DQB1*03:37', 'HLA-DQA1*04:04-HLA-DQB1*03:38',
         'HLA-DQA1*04:04-HLA-DQB1*04:01', 'HLA-DQA1*04:04-HLA-DQB1*04:02', 'HLA-DQA1*04:04-HLA-DQB1*04:03', 'HLA-DQA1*04:04-HLA-DQB1*04:04',
         'HLA-DQA1*04:04-HLA-DQB1*04:05', 'HLA-DQA1*04:04-HLA-DQB1*04:06',
         'HLA-DQA1*04:04-HLA-DQB1*04:07', 'HLA-DQA1*04:04-HLA-DQB1*04:08', 'HLA-DQA1*04:04-HLA-DQB1*05:01', 'HLA-DQA1*04:04-HLA-DQB1*05:02',
         'HLA-DQA1*04:04-HLA-DQB1*05:03', 'HLA-DQA1*04:04-HLA-DQB1*05:05',
         'HLA-DQA1*04:04-HLA-DQB1*05:06', 'HLA-DQA1*04:04-HLA-DQB1*05:07', 'HLA-DQA1*04:04-HLA-DQB1*05:08', 'HLA-DQA1*04:04-HLA-DQB1*05:09',
         'HLA-DQA1*04:04-HLA-DQB1*05:10', 'HLA-DQA1*04:04-HLA-DQB1*05:11',
         'HLA-DQA1*04:04-HLA-DQB1*05:12', 'HLA-DQA1*04:04-HLA-DQB1*05:13', 'HLA-DQA1*04:04-HLA-DQB1*05:14', 'HLA-DQA1*04:04-HLA-DQB1*06:01',
         'HLA-DQA1*04:04-HLA-DQB1*06:02', 'HLA-DQA1*04:04-HLA-DQB1*06:03',
         'HLA-DQA1*04:04-HLA-DQB1*06:04', 'HLA-DQA1*04:04-HLA-DQB1*06:07', 'HLA-DQA1*04:04-HLA-DQB1*06:08', 'HLA-DQA1*04:04-HLA-DQB1*06:09',
         'HLA-DQA1*04:04-HLA-DQB1*06:10', 'HLA-DQA1*04:04-HLA-DQB1*06:11',
         'HLA-DQA1*04:04-HLA-DQB1*06:12', 'HLA-DQA1*04:04-HLA-DQB1*06:14', 'HLA-DQA1*04:04-HLA-DQB1*06:15', 'HLA-DQA1*04:04-HLA-DQB1*06:16',
         'HLA-DQA1*04:04-HLA-DQB1*06:17', 'HLA-DQA1*04:04-HLA-DQB1*06:18',
         'HLA-DQA1*04:04-HLA-DQB1*06:19', 'HLA-DQA1*04:04-HLA-DQB1*06:21', 'HLA-DQA1*04:04-HLA-DQB1*06:22', 'HLA-DQA1*04:04-HLA-DQB1*06:23',
         'HLA-DQA1*04:04-HLA-DQB1*06:24', 'HLA-DQA1*04:04-HLA-DQB1*06:25',
         'HLA-DQA1*04:04-HLA-DQB1*06:27', 'HLA-DQA1*04:04-HLA-DQB1*06:28', 'HLA-DQA1*04:04-HLA-DQB1*06:29', 'HLA-DQA1*04:04-HLA-DQB1*06:30',
         'HLA-DQA1*04:04-HLA-DQB1*06:31', 'HLA-DQA1*04:04-HLA-DQB1*06:32',
         'HLA-DQA1*04:04-HLA-DQB1*06:33', 'HLA-DQA1*04:04-HLA-DQB1*06:34', 'HLA-DQA1*04:04-HLA-DQB1*06:35', 'HLA-DQA1*04:04-HLA-DQB1*06:36',
         'HLA-DQA1*04:04-HLA-DQB1*06:37', 'HLA-DQA1*04:04-HLA-DQB1*06:38',
         'HLA-DQA1*04:04-HLA-DQB1*06:39', 'HLA-DQA1*04:04-HLA-DQB1*06:40', 'HLA-DQA1*04:04-HLA-DQB1*06:41', 'HLA-DQA1*04:04-HLA-DQB1*06:42',
         'HLA-DQA1*04:04-HLA-DQB1*06:43', 'HLA-DQA1*04:04-HLA-DQB1*06:44',
         'HLA-DQA1*05:01-HLA-DQB1*02:01', 'HLA-DQA1*05:01-HLA-DQB1*02:02', 'HLA-DQA1*05:01-HLA-DQB1*02:03', 'HLA-DQA1*05:01-HLA-DQB1*02:04',
         'HLA-DQA1*05:01-HLA-DQB1*02:05', 'HLA-DQA1*05:01-HLA-DQB1*02:06',
         'HLA-DQA1*05:01-HLA-DQB1*03:01', 'HLA-DQA1*05:01-HLA-DQB1*03:02', 'HLA-DQA1*05:01-HLA-DQB1*03:03', 'HLA-DQA1*05:01-HLA-DQB1*03:04',
         'HLA-DQA1*05:01-HLA-DQB1*03:05', 'HLA-DQA1*05:01-HLA-DQB1*03:06',
         'HLA-DQA1*05:01-HLA-DQB1*03:07', 'HLA-DQA1*05:01-HLA-DQB1*03:08', 'HLA-DQA1*05:01-HLA-DQB1*03:09', 'HLA-DQA1*05:01-HLA-DQB1*03:10',
         'HLA-DQA1*05:01-HLA-DQB1*03:11', 'HLA-DQA1*05:01-HLA-DQB1*03:12',
         'HLA-DQA1*05:01-HLA-DQB1*03:13', 'HLA-DQA1*05:01-HLA-DQB1*03:14', 'HLA-DQA1*05:01-HLA-DQB1*03:15', 'HLA-DQA1*05:01-HLA-DQB1*03:16',
         'HLA-DQA1*05:01-HLA-DQB1*03:17', 'HLA-DQA1*05:01-HLA-DQB1*03:18',
         'HLA-DQA1*05:01-HLA-DQB1*03:19', 'HLA-DQA1*05:01-HLA-DQB1*03:20', 'HLA-DQA1*05:01-HLA-DQB1*03:21', 'HLA-DQA1*05:01-HLA-DQB1*03:22',
         'HLA-DQA1*05:01-HLA-DQB1*03:23', 'HLA-DQA1*05:01-HLA-DQB1*03:24',
         'HLA-DQA1*05:01-HLA-DQB1*03:25', 'HLA-DQA1*05:01-HLA-DQB1*03:26', 'HLA-DQA1*05:01-HLA-DQB1*03:27', 'HLA-DQA1*05:01-HLA-DQB1*03:28',
         'HLA-DQA1*05:01-HLA-DQB1*03:29', 'HLA-DQA1*05:01-HLA-DQB1*03:30',
         'HLA-DQA1*05:01-HLA-DQB1*03:31', 'HLA-DQA1*05:01-HLA-DQB1*03:32', 'HLA-DQA1*05:01-HLA-DQB1*03:33', 'HLA-DQA1*05:01-HLA-DQB1*03:34',
         'HLA-DQA1*05:01-HLA-DQB1*03:35', 'HLA-DQA1*05:01-HLA-DQB1*03:36',
         'HLA-DQA1*05:01-HLA-DQB1*03:37', 'HLA-DQA1*05:01-HLA-DQB1*03:38', 'HLA-DQA1*05:01-HLA-DQB1*04:01', 'HLA-DQA1*05:01-HLA-DQB1*04:02',
         'HLA-DQA1*05:01-HLA-DQB1*04:03', 'HLA-DQA1*05:01-HLA-DQB1*04:04',
         'HLA-DQA1*05:01-HLA-DQB1*04:05', 'HLA-DQA1*05:01-HLA-DQB1*04:06', 'HLA-DQA1*05:01-HLA-DQB1*04:07', 'HLA-DQA1*05:01-HLA-DQB1*04:08',
         'HLA-DQA1*05:01-HLA-DQB1*05:01', 'HLA-DQA1*05:01-HLA-DQB1*05:02',
         'HLA-DQA1*05:01-HLA-DQB1*05:03', 'HLA-DQA1*05:01-HLA-DQB1*05:05', 'HLA-DQA1*05:01-HLA-DQB1*05:06', 'HLA-DQA1*05:01-HLA-DQB1*05:07',
         'HLA-DQA1*05:01-HLA-DQB1*05:08', 'HLA-DQA1*05:01-HLA-DQB1*05:09',
         'HLA-DQA1*05:01-HLA-DQB1*05:10', 'HLA-DQA1*05:01-HLA-DQB1*05:11', 'HLA-DQA1*05:01-HLA-DQB1*05:12', 'HLA-DQA1*05:01-HLA-DQB1*05:13',
         'HLA-DQA1*05:01-HLA-DQB1*05:14', 'HLA-DQA1*05:01-HLA-DQB1*06:01',
         'HLA-DQA1*05:01-HLA-DQB1*06:02', 'HLA-DQA1*05:01-HLA-DQB1*06:03', 'HLA-DQA1*05:01-HLA-DQB1*06:04', 'HLA-DQA1*05:01-HLA-DQB1*06:07',
         'HLA-DQA1*05:01-HLA-DQB1*06:08', 'HLA-DQA1*05:01-HLA-DQB1*06:09',
         'HLA-DQA1*05:01-HLA-DQB1*06:10', 'HLA-DQA1*05:01-HLA-DQB1*06:11', 'HLA-DQA1*05:01-HLA-DQB1*06:12', 'HLA-DQA1*05:01-HLA-DQB1*06:14',
         'HLA-DQA1*05:01-HLA-DQB1*06:15', 'HLA-DQA1*05:01-HLA-DQB1*06:16',
         'HLA-DQA1*05:01-HLA-DQB1*06:17', 'HLA-DQA1*05:01-HLA-DQB1*06:18', 'HLA-DQA1*05:01-HLA-DQB1*06:19', 'HLA-DQA1*05:01-HLA-DQB1*06:21',
         'HLA-DQA1*05:01-HLA-DQB1*06:22', 'HLA-DQA1*05:01-HLA-DQB1*06:23',
         'HLA-DQA1*05:01-HLA-DQB1*06:24', 'HLA-DQA1*05:01-HLA-DQB1*06:25', 'HLA-DQA1*05:01-HLA-DQB1*06:27', 'HLA-DQA1*05:01-HLA-DQB1*06:28',
         'HLA-DQA1*05:01-HLA-DQB1*06:29', 'HLA-DQA1*05:01-HLA-DQB1*06:30',
         'HLA-DQA1*05:01-HLA-DQB1*06:31', 'HLA-DQA1*05:01-HLA-DQB1*06:32', 'HLA-DQA1*05:01-HLA-DQB1*06:33', 'HLA-DQA1*05:01-HLA-DQB1*06:34',
         'HLA-DQA1*05:01-HLA-DQB1*06:35', 'HLA-DQA1*05:01-HLA-DQB1*06:36',
         'HLA-DQA1*05:01-HLA-DQB1*06:37', 'HLA-DQA1*05:01-HLA-DQB1*06:38', 'HLA-DQA1*05:01-HLA-DQB1*06:39', 'HLA-DQA1*05:01-HLA-DQB1*06:40',
         'HLA-DQA1*05:01-HLA-DQB1*06:41', 'HLA-DQA1*05:01-HLA-DQB1*06:42',
         'HLA-DQA1*05:01-HLA-DQB1*06:43', 'HLA-DQA1*05:01-HLA-DQB1*06:44', 'HLA-DQA1*05:03-HLA-DQB1*02:01', 'HLA-DQA1*05:03-HLA-DQB1*02:02',
         'HLA-DQA1*05:03-HLA-DQB1*02:03', 'HLA-DQA1*05:03-HLA-DQB1*02:04',
         'HLA-DQA1*05:03-HLA-DQB1*02:05', 'HLA-DQA1*05:03-HLA-DQB1*02:06', 'HLA-DQA1*05:03-HLA-DQB1*03:01', 'HLA-DQA1*05:03-HLA-DQB1*03:02',
         'HLA-DQA1*05:03-HLA-DQB1*03:03', 'HLA-DQA1*05:03-HLA-DQB1*03:04',
         'HLA-DQA1*05:03-HLA-DQB1*03:05', 'HLA-DQA1*05:03-HLA-DQB1*03:06', 'HLA-DQA1*05:03-HLA-DQB1*03:07', 'HLA-DQA1*05:03-HLA-DQB1*03:08',
         'HLA-DQA1*05:03-HLA-DQB1*03:09', 'HLA-DQA1*05:03-HLA-DQB1*03:10',
         'HLA-DQA1*05:03-HLA-DQB1*03:11', 'HLA-DQA1*05:03-HLA-DQB1*03:12', 'HLA-DQA1*05:03-HLA-DQB1*03:13', 'HLA-DQA1*05:03-HLA-DQB1*03:14',
         'HLA-DQA1*05:03-HLA-DQB1*03:15', 'HLA-DQA1*05:03-HLA-DQB1*03:16',
         'HLA-DQA1*05:03-HLA-DQB1*03:17', 'HLA-DQA1*05:03-HLA-DQB1*03:18', 'HLA-DQA1*05:03-HLA-DQB1*03:19', 'HLA-DQA1*05:03-HLA-DQB1*03:20',
         'HLA-DQA1*05:03-HLA-DQB1*03:21', 'HLA-DQA1*05:03-HLA-DQB1*03:22',
         'HLA-DQA1*05:03-HLA-DQB1*03:23', 'HLA-DQA1*05:03-HLA-DQB1*03:24', 'HLA-DQA1*05:03-HLA-DQB1*03:25', 'HLA-DQA1*05:03-HLA-DQB1*03:26',
         'HLA-DQA1*05:03-HLA-DQB1*03:27', 'HLA-DQA1*05:03-HLA-DQB1*03:28',
         'HLA-DQA1*05:03-HLA-DQB1*03:29', 'HLA-DQA1*05:03-HLA-DQB1*03:30', 'HLA-DQA1*05:03-HLA-DQB1*03:31', 'HLA-DQA1*05:03-HLA-DQB1*03:32',
         'HLA-DQA1*05:03-HLA-DQB1*03:33', 'HLA-DQA1*05:03-HLA-DQB1*03:34',
         'HLA-DQA1*05:03-HLA-DQB1*03:35', 'HLA-DQA1*05:03-HLA-DQB1*03:36', 'HLA-DQA1*05:03-HLA-DQB1*03:37', 'HLA-DQA1*05:03-HLA-DQB1*03:38',
         'HLA-DQA1*05:03-HLA-DQB1*04:01', 'HLA-DQA1*05:03-HLA-DQB1*04:02',
         'HLA-DQA1*05:03-HLA-DQB1*04:03', 'HLA-DQA1*05:03-HLA-DQB1*04:04', 'HLA-DQA1*05:03-HLA-DQB1*04:05', 'HLA-DQA1*05:03-HLA-DQB1*04:06',
         'HLA-DQA1*05:03-HLA-DQB1*04:07', 'HLA-DQA1*05:03-HLA-DQB1*04:08',
         'HLA-DQA1*05:03-HLA-DQB1*05:01', 'HLA-DQA1*05:03-HLA-DQB1*05:02', 'HLA-DQA1*05:03-HLA-DQB1*05:03', 'HLA-DQA1*05:03-HLA-DQB1*05:05',
         'HLA-DQA1*05:03-HLA-DQB1*05:06', 'HLA-DQA1*05:03-HLA-DQB1*05:07',
         'HLA-DQA1*05:03-HLA-DQB1*05:08', 'HLA-DQA1*05:03-HLA-DQB1*05:09', 'HLA-DQA1*05:03-HLA-DQB1*05:10', 'HLA-DQA1*05:03-HLA-DQB1*05:11',
         'HLA-DQA1*05:03-HLA-DQB1*05:12', 'HLA-DQA1*05:03-HLA-DQB1*05:13',
         'HLA-DQA1*05:03-HLA-DQB1*05:14', 'HLA-DQA1*05:03-HLA-DQB1*06:01', 'HLA-DQA1*05:03-HLA-DQB1*06:02', 'HLA-DQA1*05:03-HLA-DQB1*06:03',
         'HLA-DQA1*05:03-HLA-DQB1*06:04', 'HLA-DQA1*05:03-HLA-DQB1*06:07',
         'HLA-DQA1*05:03-HLA-DQB1*06:08', 'HLA-DQA1*05:03-HLA-DQB1*06:09', 'HLA-DQA1*05:03-HLA-DQB1*06:10', 'HLA-DQA1*05:03-HLA-DQB1*06:11',
         'HLA-DQA1*05:03-HLA-DQB1*06:12', 'HLA-DQA1*05:03-HLA-DQB1*06:14',
         'HLA-DQA1*05:03-HLA-DQB1*06:15', 'HLA-DQA1*05:03-HLA-DQB1*06:16', 'HLA-DQA1*05:03-HLA-DQB1*06:17', 'HLA-DQA1*05:03-HLA-DQB1*06:18',
         'HLA-DQA1*05:03-HLA-DQB1*06:19', 'HLA-DQA1*05:03-HLA-DQB1*06:21',
         'HLA-DQA1*05:03-HLA-DQB1*06:22', 'HLA-DQA1*05:03-HLA-DQB1*06:23', 'HLA-DQA1*05:03-HLA-DQB1*06:24', 'HLA-DQA1*05:03-HLA-DQB1*06:25',
         'HLA-DQA1*05:03-HLA-DQB1*06:27', 'HLA-DQA1*05:03-HLA-DQB1*06:28',
         'HLA-DQA1*05:03-HLA-DQB1*06:29', 'HLA-DQA1*05:03-HLA-DQB1*06:30', 'HLA-DQA1*05:03-HLA-DQB1*06:31', 'HLA-DQA1*05:03-HLA-DQB1*06:32',
         'HLA-DQA1*05:03-HLA-DQB1*06:33', 'HLA-DQA1*05:03-HLA-DQB1*06:34',
         'HLA-DQA1*05:03-HLA-DQB1*06:35', 'HLA-DQA1*05:03-HLA-DQB1*06:36', 'HLA-DQA1*05:03-HLA-DQB1*06:37', 'HLA-DQA1*05:03-HLA-DQB1*06:38',
         'HLA-DQA1*05:03-HLA-DQB1*06:39', 'HLA-DQA1*05:03-HLA-DQB1*06:40',
         'HLA-DQA1*05:03-HLA-DQB1*06:41', 'HLA-DQA1*05:03-HLA-DQB1*06:42', 'HLA-DQA1*05:03-HLA-DQB1*06:43', 'HLA-DQA1*05:03-HLA-DQB1*06:44',
         'HLA-DQA1*05:04-HLA-DQB1*02:01', 'HLA-DQA1*05:04-HLA-DQB1*02:02',
         'HLA-DQA1*05:04-HLA-DQB1*02:03', 'HLA-DQA1*05:04-HLA-DQB1*02:04', 'HLA-DQA1*05:04-HLA-DQB1*02:05', 'HLA-DQA1*05:04-HLA-DQB1*02:06',
         'HLA-DQA1*05:04-HLA-DQB1*03:01', 'HLA-DQA1*05:04-HLA-DQB1*03:02',
         'HLA-DQA1*05:04-HLA-DQB1*03:03', 'HLA-DQA1*05:04-HLA-DQB1*03:04', 'HLA-DQA1*05:04-HLA-DQB1*03:05', 'HLA-DQA1*05:04-HLA-DQB1*03:06',
         'HLA-DQA1*05:04-HLA-DQB1*03:07', 'HLA-DQA1*05:04-HLA-DQB1*03:08',
         'HLA-DQA1*05:04-HLA-DQB1*03:09', 'HLA-DQA1*05:04-HLA-DQB1*03:10', 'HLA-DQA1*05:04-HLA-DQB1*03:11', 'HLA-DQA1*05:04-HLA-DQB1*03:12',
         'HLA-DQA1*05:04-HLA-DQB1*03:13', 'HLA-DQA1*05:04-HLA-DQB1*03:14',
         'HLA-DQA1*05:04-HLA-DQB1*03:15', 'HLA-DQA1*05:04-HLA-DQB1*03:16', 'HLA-DQA1*05:04-HLA-DQB1*03:17', 'HLA-DQA1*05:04-HLA-DQB1*03:18',
         'HLA-DQA1*05:04-HLA-DQB1*03:19', 'HLA-DQA1*05:04-HLA-DQB1*03:20',
         'HLA-DQA1*05:04-HLA-DQB1*03:21', 'HLA-DQA1*05:04-HLA-DQB1*03:22', 'HLA-DQA1*05:04-HLA-DQB1*03:23', 'HLA-DQA1*05:04-HLA-DQB1*03:24',
         'HLA-DQA1*05:04-HLA-DQB1*03:25', 'HLA-DQA1*05:04-HLA-DQB1*03:26',
         'HLA-DQA1*05:04-HLA-DQB1*03:27', 'HLA-DQA1*05:04-HLA-DQB1*03:28', 'HLA-DQA1*05:04-HLA-DQB1*03:29', 'HLA-DQA1*05:04-HLA-DQB1*03:30',
         'HLA-DQA1*05:04-HLA-DQB1*03:31', 'HLA-DQA1*05:04-HLA-DQB1*03:32',
         'HLA-DQA1*05:04-HLA-DQB1*03:33', 'HLA-DQA1*05:04-HLA-DQB1*03:34', 'HLA-DQA1*05:04-HLA-DQB1*03:35', 'HLA-DQA1*05:04-HLA-DQB1*03:36',
         'HLA-DQA1*05:04-HLA-DQB1*03:37', 'HLA-DQA1*05:04-HLA-DQB1*03:38',
         'HLA-DQA1*05:04-HLA-DQB1*04:01', 'HLA-DQA1*05:04-HLA-DQB1*04:02', 'HLA-DQA1*05:04-HLA-DQB1*04:03', 'HLA-DQA1*05:04-HLA-DQB1*04:04',
         'HLA-DQA1*05:04-HLA-DQB1*04:05', 'HLA-DQA1*05:04-HLA-DQB1*04:06',
         'HLA-DQA1*05:04-HLA-DQB1*04:07', 'HLA-DQA1*05:04-HLA-DQB1*04:08', 'HLA-DQA1*05:04-HLA-DQB1*05:01', 'HLA-DQA1*05:04-HLA-DQB1*05:02',
         'HLA-DQA1*05:04-HLA-DQB1*05:03', 'HLA-DQA1*05:04-HLA-DQB1*05:05',
         'HLA-DQA1*05:04-HLA-DQB1*05:06', 'HLA-DQA1*05:04-HLA-DQB1*05:07', 'HLA-DQA1*05:04-HLA-DQB1*05:08', 'HLA-DQA1*05:04-HLA-DQB1*05:09',
         'HLA-DQA1*05:04-HLA-DQB1*05:10', 'HLA-DQA1*05:04-HLA-DQB1*05:11',
         'HLA-DQA1*05:04-HLA-DQB1*05:12', 'HLA-DQA1*05:04-HLA-DQB1*05:13', 'HLA-DQA1*05:04-HLA-DQB1*05:14', 'HLA-DQA1*05:04-HLA-DQB1*06:01',
         'HLA-DQA1*05:04-HLA-DQB1*06:02', 'HLA-DQA1*05:04-HLA-DQB1*06:03',
         'HLA-DQA1*05:04-HLA-DQB1*06:04', 'HLA-DQA1*05:04-HLA-DQB1*06:07', 'HLA-DQA1*05:04-HLA-DQB1*06:08', 'HLA-DQA1*05:04-HLA-DQB1*06:09',
         'HLA-DQA1*05:04-HLA-DQB1*06:10', 'HLA-DQA1*05:04-HLA-DQB1*06:11',
         'HLA-DQA1*05:04-HLA-DQB1*06:12', 'HLA-DQA1*05:04-HLA-DQB1*06:14', 'HLA-DQA1*05:04-HLA-DQB1*06:15', 'HLA-DQA1*05:04-HLA-DQB1*06:16',
         'HLA-DQA1*05:04-HLA-DQB1*06:17', 'HLA-DQA1*05:04-HLA-DQB1*06:18',
         'HLA-DQA1*05:04-HLA-DQB1*06:19', 'HLA-DQA1*05:04-HLA-DQB1*06:21', 'HLA-DQA1*05:04-HLA-DQB1*06:22', 'HLA-DQA1*05:04-HLA-DQB1*06:23',
         'HLA-DQA1*05:04-HLA-DQB1*06:24', 'HLA-DQA1*05:04-HLA-DQB1*06:25',
         'HLA-DQA1*05:04-HLA-DQB1*06:27', 'HLA-DQA1*05:04-HLA-DQB1*06:28', 'HLA-DQA1*05:04-HLA-DQB1*06:29', 'HLA-DQA1*05:04-HLA-DQB1*06:30',
         'HLA-DQA1*05:04-HLA-DQB1*06:31', 'HLA-DQA1*05:04-HLA-DQB1*06:32',
         'HLA-DQA1*05:04-HLA-DQB1*06:33', 'HLA-DQA1*05:04-HLA-DQB1*06:34', 'HLA-DQA1*05:04-HLA-DQB1*06:35', 'HLA-DQA1*05:04-HLA-DQB1*06:36',
         'HLA-DQA1*05:04-HLA-DQB1*06:37', 'HLA-DQA1*05:04-HLA-DQB1*06:38',
         'HLA-DQA1*05:04-HLA-DQB1*06:39', 'HLA-DQA1*05:04-HLA-DQB1*06:40', 'HLA-DQA1*05:04-HLA-DQB1*06:41', 'HLA-DQA1*05:04-HLA-DQB1*06:42',
         'HLA-DQA1*05:04-HLA-DQB1*06:43', 'HLA-DQA1*05:04-HLA-DQB1*06:44',
         'HLA-DQA1*05:05-HLA-DQB1*02:01', 'HLA-DQA1*05:05-HLA-DQB1*02:02', 'HLA-DQA1*05:05-HLA-DQB1*02:03', 'HLA-DQA1*05:05-HLA-DQB1*02:04',
         'HLA-DQA1*05:05-HLA-DQB1*02:05', 'HLA-DQA1*05:05-HLA-DQB1*02:06',
         'HLA-DQA1*05:05-HLA-DQB1*03:01', 'HLA-DQA1*05:05-HLA-DQB1*03:02', 'HLA-DQA1*05:05-HLA-DQB1*03:03', 'HLA-DQA1*05:05-HLA-DQB1*03:04',
         'HLA-DQA1*05:05-HLA-DQB1*03:05', 'HLA-DQA1*05:05-HLA-DQB1*03:06',
         'HLA-DQA1*05:05-HLA-DQB1*03:07', 'HLA-DQA1*05:05-HLA-DQB1*03:08', 'HLA-DQA1*05:05-HLA-DQB1*03:09', 'HLA-DQA1*05:05-HLA-DQB1*03:10',
         'HLA-DQA1*05:05-HLA-DQB1*03:11', 'HLA-DQA1*05:05-HLA-DQB1*03:12',
         'HLA-DQA1*05:05-HLA-DQB1*03:13', 'HLA-DQA1*05:05-HLA-DQB1*03:14', 'HLA-DQA1*05:05-HLA-DQB1*03:15', 'HLA-DQA1*05:05-HLA-DQB1*03:16',
         'HLA-DQA1*05:05-HLA-DQB1*03:17', 'HLA-DQA1*05:05-HLA-DQB1*03:18',
         'HLA-DQA1*05:05-HLA-DQB1*03:19', 'HLA-DQA1*05:05-HLA-DQB1*03:20', 'HLA-DQA1*05:05-HLA-DQB1*03:21', 'HLA-DQA1*05:05-HLA-DQB1*03:22',
         'HLA-DQA1*05:05-HLA-DQB1*03:23', 'HLA-DQA1*05:05-HLA-DQB1*03:24',
         'HLA-DQA1*05:05-HLA-DQB1*03:25', 'HLA-DQA1*05:05-HLA-DQB1*03:26', 'HLA-DQA1*05:05-HLA-DQB1*03:27', 'HLA-DQA1*05:05-HLA-DQB1*03:28',
         'HLA-DQA1*05:05-HLA-DQB1*03:29', 'HLA-DQA1*05:05-HLA-DQB1*03:30',
         'HLA-DQA1*05:05-HLA-DQB1*03:31', 'HLA-DQA1*05:05-HLA-DQB1*03:32', 'HLA-DQA1*05:05-HLA-DQB1*03:33', 'HLA-DQA1*05:05-HLA-DQB1*03:34',
         'HLA-DQA1*05:05-HLA-DQB1*03:35', 'HLA-DQA1*05:05-HLA-DQB1*03:36',
         'HLA-DQA1*05:05-HLA-DQB1*03:37', 'HLA-DQA1*05:05-HLA-DQB1*03:38', 'HLA-DQA1*05:05-HLA-DQB1*04:01', 'HLA-DQA1*05:05-HLA-DQB1*04:02',
         'HLA-DQA1*05:05-HLA-DQB1*04:03', 'HLA-DQA1*05:05-HLA-DQB1*04:04',
         'HLA-DQA1*05:05-HLA-DQB1*04:05', 'HLA-DQA1*05:05-HLA-DQB1*04:06', 'HLA-DQA1*05:05-HLA-DQB1*04:07', 'HLA-DQA1*05:05-HLA-DQB1*04:08',
         'HLA-DQA1*05:05-HLA-DQB1*05:01', 'HLA-DQA1*05:05-HLA-DQB1*05:02',
         'HLA-DQA1*05:05-HLA-DQB1*05:03', 'HLA-DQA1*05:05-HLA-DQB1*05:05', 'HLA-DQA1*05:05-HLA-DQB1*05:06', 'HLA-DQA1*05:05-HLA-DQB1*05:07',
         'HLA-DQA1*05:05-HLA-DQB1*05:08', 'HLA-DQA1*05:05-HLA-DQB1*05:09',
         'HLA-DQA1*05:05-HLA-DQB1*05:10', 'HLA-DQA1*05:05-HLA-DQB1*05:11', 'HLA-DQA1*05:05-HLA-DQB1*05:12', 'HLA-DQA1*05:05-HLA-DQB1*05:13',
         'HLA-DQA1*05:05-HLA-DQB1*05:14', 'HLA-DQA1*05:05-HLA-DQB1*06:01',
         'HLA-DQA1*05:05-HLA-DQB1*06:02', 'HLA-DQA1*05:05-HLA-DQB1*06:03', 'HLA-DQA1*05:05-HLA-DQB1*06:04', 'HLA-DQA1*05:05-HLA-DQB1*06:07',
         'HLA-DQA1*05:05-HLA-DQB1*06:08', 'HLA-DQA1*05:05-HLA-DQB1*06:09',
         'HLA-DQA1*05:05-HLA-DQB1*06:10', 'HLA-DQA1*05:05-HLA-DQB1*06:11', 'HLA-DQA1*05:05-HLA-DQB1*06:12', 'HLA-DQA1*05:05-HLA-DQB1*06:14',
         'HLA-DQA1*05:05-HLA-DQB1*06:15', 'HLA-DQA1*05:05-HLA-DQB1*06:16',
         'HLA-DQA1*05:05-HLA-DQB1*06:17', 'HLA-DQA1*05:05-HLA-DQB1*06:18', 'HLA-DQA1*05:05-HLA-DQB1*06:19', 'HLA-DQA1*05:05-HLA-DQB1*06:21',
         'HLA-DQA1*05:05-HLA-DQB1*06:22', 'HLA-DQA1*05:05-HLA-DQB1*06:23',
         'HLA-DQA1*05:05-HLA-DQB1*06:24', 'HLA-DQA1*05:05-HLA-DQB1*06:25', 'HLA-DQA1*05:05-HLA-DQB1*06:27', 'HLA-DQA1*05:05-HLA-DQB1*06:28',
         'HLA-DQA1*05:05-HLA-DQB1*06:29', 'HLA-DQA1*05:05-HLA-DQB1*06:30',
         'HLA-DQA1*05:05-HLA-DQB1*06:31', 'HLA-DQA1*05:05-HLA-DQB1*06:32', 'HLA-DQA1*05:05-HLA-DQB1*06:33', 'HLA-DQA1*05:05-HLA-DQB1*06:34',
         'HLA-DQA1*05:05-HLA-DQB1*06:35', 'HLA-DQA1*05:05-HLA-DQB1*06:36',
         'HLA-DQA1*05:05-HLA-DQB1*06:37', 'HLA-DQA1*05:05-HLA-DQB1*06:38', 'HLA-DQA1*05:05-HLA-DQB1*06:39', 'HLA-DQA1*05:05-HLA-DQB1*06:40',
         'HLA-DQA1*05:05-HLA-DQB1*06:41', 'HLA-DQA1*05:05-HLA-DQB1*06:42',
         'HLA-DQA1*05:05-HLA-DQB1*06:43', 'HLA-DQA1*05:05-HLA-DQB1*06:44', 'HLA-DQA1*05:06-HLA-DQB1*02:01', 'HLA-DQA1*05:06-HLA-DQB1*02:02',
         'HLA-DQA1*05:06-HLA-DQB1*02:03', 'HLA-DQA1*05:06-HLA-DQB1*02:04',
         'HLA-DQA1*05:06-HLA-DQB1*02:05', 'HLA-DQA1*05:06-HLA-DQB1*02:06', 'HLA-DQA1*05:06-HLA-DQB1*03:01', 'HLA-DQA1*05:06-HLA-DQB1*03:02',
         'HLA-DQA1*05:06-HLA-DQB1*03:03', 'HLA-DQA1*05:06-HLA-DQB1*03:04',
         'HLA-DQA1*05:06-HLA-DQB1*03:05', 'HLA-DQA1*05:06-HLA-DQB1*03:06', 'HLA-DQA1*05:06-HLA-DQB1*03:07', 'HLA-DQA1*05:06-HLA-DQB1*03:08',
         'HLA-DQA1*05:06-HLA-DQB1*03:09', 'HLA-DQA1*05:06-HLA-DQB1*03:10',
         'HLA-DQA1*05:06-HLA-DQB1*03:11', 'HLA-DQA1*05:06-HLA-DQB1*03:12', 'HLA-DQA1*05:06-HLA-DQB1*03:13', 'HLA-DQA1*05:06-HLA-DQB1*03:14',
         'HLA-DQA1*05:06-HLA-DQB1*03:15', 'HLA-DQA1*05:06-HLA-DQB1*03:16',
         'HLA-DQA1*05:06-HLA-DQB1*03:17', 'HLA-DQA1*05:06-HLA-DQB1*03:18', 'HLA-DQA1*05:06-HLA-DQB1*03:19', 'HLA-DQA1*05:06-HLA-DQB1*03:20',
         'HLA-DQA1*05:06-HLA-DQB1*03:21', 'HLA-DQA1*05:06-HLA-DQB1*03:22',
         'HLA-DQA1*05:06-HLA-DQB1*03:23', 'HLA-DQA1*05:06-HLA-DQB1*03:24', 'HLA-DQA1*05:06-HLA-DQB1*03:25', 'HLA-DQA1*05:06-HLA-DQB1*03:26',
         'HLA-DQA1*05:06-HLA-DQB1*03:27', 'HLA-DQA1*05:06-HLA-DQB1*03:28',
         'HLA-DQA1*05:06-HLA-DQB1*03:29', 'HLA-DQA1*05:06-HLA-DQB1*03:30', 'HLA-DQA1*05:06-HLA-DQB1*03:31', 'HLA-DQA1*05:06-HLA-DQB1*03:32',
         'HLA-DQA1*05:06-HLA-DQB1*03:33', 'HLA-DQA1*05:06-HLA-DQB1*03:34',
         'HLA-DQA1*05:06-HLA-DQB1*03:35', 'HLA-DQA1*05:06-HLA-DQB1*03:36', 'HLA-DQA1*05:06-HLA-DQB1*03:37', 'HLA-DQA1*05:06-HLA-DQB1*03:38',
         'HLA-DQA1*05:06-HLA-DQB1*04:01', 'HLA-DQA1*05:06-HLA-DQB1*04:02',
         'HLA-DQA1*05:06-HLA-DQB1*04:03', 'HLA-DQA1*05:06-HLA-DQB1*04:04', 'HLA-DQA1*05:06-HLA-DQB1*04:05', 'HLA-DQA1*05:06-HLA-DQB1*04:06',
         'HLA-DQA1*05:06-HLA-DQB1*04:07', 'HLA-DQA1*05:06-HLA-DQB1*04:08',
         'HLA-DQA1*05:06-HLA-DQB1*05:01', 'HLA-DQA1*05:06-HLA-DQB1*05:02', 'HLA-DQA1*05:06-HLA-DQB1*05:03', 'HLA-DQA1*05:06-HLA-DQB1*05:05',
         'HLA-DQA1*05:06-HLA-DQB1*05:06', 'HLA-DQA1*05:06-HLA-DQB1*05:07',
         'HLA-DQA1*05:06-HLA-DQB1*05:08', 'HLA-DQA1*05:06-HLA-DQB1*05:09', 'HLA-DQA1*05:06-HLA-DQB1*05:10', 'HLA-DQA1*05:06-HLA-DQB1*05:11',
         'HLA-DQA1*05:06-HLA-DQB1*05:12', 'HLA-DQA1*05:06-HLA-DQB1*05:13',
         'HLA-DQA1*05:06-HLA-DQB1*05:14', 'HLA-DQA1*05:06-HLA-DQB1*06:01', 'HLA-DQA1*05:06-HLA-DQB1*06:02', 'HLA-DQA1*05:06-HLA-DQB1*06:03',
         'HLA-DQA1*05:06-HLA-DQB1*06:04', 'HLA-DQA1*05:06-HLA-DQB1*06:07',
         'HLA-DQA1*05:06-HLA-DQB1*06:08', 'HLA-DQA1*05:06-HLA-DQB1*06:09', 'HLA-DQA1*05:06-HLA-DQB1*06:10', 'HLA-DQA1*05:06-HLA-DQB1*06:11',
         'HLA-DQA1*05:06-HLA-DQB1*06:12', 'HLA-DQA1*05:06-HLA-DQB1*06:14',
         'HLA-DQA1*05:06-HLA-DQB1*06:15', 'HLA-DQA1*05:06-HLA-DQB1*06:16', 'HLA-DQA1*05:06-HLA-DQB1*06:17', 'HLA-DQA1*05:06-HLA-DQB1*06:18',
         'HLA-DQA1*05:06-HLA-DQB1*06:19', 'HLA-DQA1*05:06-HLA-DQB1*06:21',
         'HLA-DQA1*05:06-HLA-DQB1*06:22', 'HLA-DQA1*05:06-HLA-DQB1*06:23', 'HLA-DQA1*05:06-HLA-DQB1*06:24', 'HLA-DQA1*05:06-HLA-DQB1*06:25',
         'HLA-DQA1*05:06-HLA-DQB1*06:27', 'HLA-DQA1*05:06-HLA-DQB1*06:28',
         'HLA-DQA1*05:06-HLA-DQB1*06:29', 'HLA-DQA1*05:06-HLA-DQB1*06:30', 'HLA-DQA1*05:06-HLA-DQB1*06:31', 'HLA-DQA1*05:06-HLA-DQB1*06:32',
         'HLA-DQA1*05:06-HLA-DQB1*06:33', 'HLA-DQA1*05:06-HLA-DQB1*06:34',
         'HLA-DQA1*05:06-HLA-DQB1*06:35', 'HLA-DQA1*05:06-HLA-DQB1*06:36', 'HLA-DQA1*05:06-HLA-DQB1*06:37', 'HLA-DQA1*05:06-HLA-DQB1*06:38',
         'HLA-DQA1*05:06-HLA-DQB1*06:39', 'HLA-DQA1*05:06-HLA-DQB1*06:40',
         'HLA-DQA1*05:06-HLA-DQB1*06:41', 'HLA-DQA1*05:06-HLA-DQB1*06:42', 'HLA-DQA1*05:06-HLA-DQB1*06:43', 'HLA-DQA1*05:06-HLA-DQB1*06:44',
         'HLA-DQA1*05:07-HLA-DQB1*02:01', 'HLA-DQA1*05:07-HLA-DQB1*02:02',
         'HLA-DQA1*05:07-HLA-DQB1*02:03', 'HLA-DQA1*05:07-HLA-DQB1*02:04', 'HLA-DQA1*05:07-HLA-DQB1*02:05', 'HLA-DQA1*05:07-HLA-DQB1*02:06',
         'HLA-DQA1*05:07-HLA-DQB1*03:01', 'HLA-DQA1*05:07-HLA-DQB1*03:02',
         'HLA-DQA1*05:07-HLA-DQB1*03:03', 'HLA-DQA1*05:07-HLA-DQB1*03:04', 'HLA-DQA1*05:07-HLA-DQB1*03:05', 'HLA-DQA1*05:07-HLA-DQB1*03:06',
         'HLA-DQA1*05:07-HLA-DQB1*03:07', 'HLA-DQA1*05:07-HLA-DQB1*03:08',
         'HLA-DQA1*05:07-HLA-DQB1*03:09', 'HLA-DQA1*05:07-HLA-DQB1*03:10', 'HLA-DQA1*05:07-HLA-DQB1*03:11', 'HLA-DQA1*05:07-HLA-DQB1*03:12',
         'HLA-DQA1*05:07-HLA-DQB1*03:13', 'HLA-DQA1*05:07-HLA-DQB1*03:14',
         'HLA-DQA1*05:07-HLA-DQB1*03:15', 'HLA-DQA1*05:07-HLA-DQB1*03:16', 'HLA-DQA1*05:07-HLA-DQB1*03:17', 'HLA-DQA1*05:07-HLA-DQB1*03:18',
         'HLA-DQA1*05:07-HLA-DQB1*03:19', 'HLA-DQA1*05:07-HLA-DQB1*03:20',
         'HLA-DQA1*05:07-HLA-DQB1*03:21', 'HLA-DQA1*05:07-HLA-DQB1*03:22', 'HLA-DQA1*05:07-HLA-DQB1*03:23', 'HLA-DQA1*05:07-HLA-DQB1*03:24',
         'HLA-DQA1*05:07-HLA-DQB1*03:25', 'HLA-DQA1*05:07-HLA-DQB1*03:26',
         'HLA-DQA1*05:07-HLA-DQB1*03:27', 'HLA-DQA1*05:07-HLA-DQB1*03:28', 'HLA-DQA1*05:07-HLA-DQB1*03:29', 'HLA-DQA1*05:07-HLA-DQB1*03:30',
         'HLA-DQA1*05:07-HLA-DQB1*03:31', 'HLA-DQA1*05:07-HLA-DQB1*03:32',
         'HLA-DQA1*05:07-HLA-DQB1*03:33', 'HLA-DQA1*05:07-HLA-DQB1*03:34', 'HLA-DQA1*05:07-HLA-DQB1*03:35', 'HLA-DQA1*05:07-HLA-DQB1*03:36',
         'HLA-DQA1*05:07-HLA-DQB1*03:37', 'HLA-DQA1*05:07-HLA-DQB1*03:38',
         'HLA-DQA1*05:07-HLA-DQB1*04:01', 'HLA-DQA1*05:07-HLA-DQB1*04:02', 'HLA-DQA1*05:07-HLA-DQB1*04:03', 'HLA-DQA1*05:07-HLA-DQB1*04:04',
         'HLA-DQA1*05:07-HLA-DQB1*04:05', 'HLA-DQA1*05:07-HLA-DQB1*04:06',
         'HLA-DQA1*05:07-HLA-DQB1*04:07', 'HLA-DQA1*05:07-HLA-DQB1*04:08', 'HLA-DQA1*05:07-HLA-DQB1*05:01', 'HLA-DQA1*05:07-HLA-DQB1*05:02',
         'HLA-DQA1*05:07-HLA-DQB1*05:03', 'HLA-DQA1*05:07-HLA-DQB1*05:05',
         'HLA-DQA1*05:07-HLA-DQB1*05:06', 'HLA-DQA1*05:07-HLA-DQB1*05:07', 'HLA-DQA1*05:07-HLA-DQB1*05:08', 'HLA-DQA1*05:07-HLA-DQB1*05:09',
         'HLA-DQA1*05:07-HLA-DQB1*05:10', 'HLA-DQA1*05:07-HLA-DQB1*05:11',
         'HLA-DQA1*05:07-HLA-DQB1*05:12', 'HLA-DQA1*05:07-HLA-DQB1*05:13', 'HLA-DQA1*05:07-HLA-DQB1*05:14', 'HLA-DQA1*05:07-HLA-DQB1*06:01',
         'HLA-DQA1*05:07-HLA-DQB1*06:02', 'HLA-DQA1*05:07-HLA-DQB1*06:03',
         'HLA-DQA1*05:07-HLA-DQB1*06:04', 'HLA-DQA1*05:07-HLA-DQB1*06:07', 'HLA-DQA1*05:07-HLA-DQB1*06:08', 'HLA-DQA1*05:07-HLA-DQB1*06:09',
         'HLA-DQA1*05:07-HLA-DQB1*06:10', 'HLA-DQA1*05:07-HLA-DQB1*06:11',
         'HLA-DQA1*05:07-HLA-DQB1*06:12', 'HLA-DQA1*05:07-HLA-DQB1*06:14', 'HLA-DQA1*05:07-HLA-DQB1*06:15', 'HLA-DQA1*05:07-HLA-DQB1*06:16',
         'HLA-DQA1*05:07-HLA-DQB1*06:17', 'HLA-DQA1*05:07-HLA-DQB1*06:18',
         'HLA-DQA1*05:07-HLA-DQB1*06:19', 'HLA-DQA1*05:07-HLA-DQB1*06:21', 'HLA-DQA1*05:07-HLA-DQB1*06:22', 'HLA-DQA1*05:07-HLA-DQB1*06:23',
         'HLA-DQA1*05:07-HLA-DQB1*06:24', 'HLA-DQA1*05:07-HLA-DQB1*06:25',
         'HLA-DQA1*05:07-HLA-DQB1*06:27', 'HLA-DQA1*05:07-HLA-DQB1*06:28', 'HLA-DQA1*05:07-HLA-DQB1*06:29', 'HLA-DQA1*05:07-HLA-DQB1*06:30',
         'HLA-DQA1*05:07-HLA-DQB1*06:31', 'HLA-DQA1*05:07-HLA-DQB1*06:32',
         'HLA-DQA1*05:07-HLA-DQB1*06:33', 'HLA-DQA1*05:07-HLA-DQB1*06:34', 'HLA-DQA1*05:07-HLA-DQB1*06:35', 'HLA-DQA1*05:07-HLA-DQB1*06:36',
         'HLA-DQA1*05:07-HLA-DQB1*06:37', 'HLA-DQA1*05:07-HLA-DQB1*06:38',
         'HLA-DQA1*05:07-HLA-DQB1*06:39', 'HLA-DQA1*05:07-HLA-DQB1*06:40', 'HLA-DQA1*05:07-HLA-DQB1*06:41', 'HLA-DQA1*05:07-HLA-DQB1*06:42',
         'HLA-DQA1*05:07-HLA-DQB1*06:43', 'HLA-DQA1*05:07-HLA-DQB1*06:44',
         'HLA-DQA1*05:08-HLA-DQB1*02:01', 'HLA-DQA1*05:08-HLA-DQB1*02:02', 'HLA-DQA1*05:08-HLA-DQB1*02:03', 'HLA-DQA1*05:08-HLA-DQB1*02:04',
         'HLA-DQA1*05:08-HLA-DQB1*02:05', 'HLA-DQA1*05:08-HLA-DQB1*02:06',
         'HLA-DQA1*05:08-HLA-DQB1*03:01', 'HLA-DQA1*05:08-HLA-DQB1*03:02', 'HLA-DQA1*05:08-HLA-DQB1*03:03', 'HLA-DQA1*05:08-HLA-DQB1*03:04',
         'HLA-DQA1*05:08-HLA-DQB1*03:05', 'HLA-DQA1*05:08-HLA-DQB1*03:06',
         'HLA-DQA1*05:08-HLA-DQB1*03:07', 'HLA-DQA1*05:08-HLA-DQB1*03:08', 'HLA-DQA1*05:08-HLA-DQB1*03:09', 'HLA-DQA1*05:08-HLA-DQB1*03:10',
         'HLA-DQA1*05:08-HLA-DQB1*03:11', 'HLA-DQA1*05:08-HLA-DQB1*03:12',
         'HLA-DQA1*05:08-HLA-DQB1*03:13', 'HLA-DQA1*05:08-HLA-DQB1*03:14', 'HLA-DQA1*05:08-HLA-DQB1*03:15', 'HLA-DQA1*05:08-HLA-DQB1*03:16',
         'HLA-DQA1*05:08-HLA-DQB1*03:17', 'HLA-DQA1*05:08-HLA-DQB1*03:18',
         'HLA-DQA1*05:08-HLA-DQB1*03:19', 'HLA-DQA1*05:08-HLA-DQB1*03:20', 'HLA-DQA1*05:08-HLA-DQB1*03:21', 'HLA-DQA1*05:08-HLA-DQB1*03:22',
         'HLA-DQA1*05:08-HLA-DQB1*03:23', 'HLA-DQA1*05:08-HLA-DQB1*03:24',
         'HLA-DQA1*05:08-HLA-DQB1*03:25', 'HLA-DQA1*05:08-HLA-DQB1*03:26', 'HLA-DQA1*05:08-HLA-DQB1*03:27', 'HLA-DQA1*05:08-HLA-DQB1*03:28',
         'HLA-DQA1*05:08-HLA-DQB1*03:29', 'HLA-DQA1*05:08-HLA-DQB1*03:30',
         'HLA-DQA1*05:08-HLA-DQB1*03:31', 'HLA-DQA1*05:08-HLA-DQB1*03:32', 'HLA-DQA1*05:08-HLA-DQB1*03:33', 'HLA-DQA1*05:08-HLA-DQB1*03:34',
         'HLA-DQA1*05:08-HLA-DQB1*03:35', 'HLA-DQA1*05:08-HLA-DQB1*03:36',
         'HLA-DQA1*05:08-HLA-DQB1*03:37', 'HLA-DQA1*05:08-HLA-DQB1*03:38', 'HLA-DQA1*05:08-HLA-DQB1*04:01', 'HLA-DQA1*05:08-HLA-DQB1*04:02',
         'HLA-DQA1*05:08-HLA-DQB1*04:03', 'HLA-DQA1*05:08-HLA-DQB1*04:04',
         'HLA-DQA1*05:08-HLA-DQB1*04:05', 'HLA-DQA1*05:08-HLA-DQB1*04:06', 'HLA-DQA1*05:08-HLA-DQB1*04:07', 'HLA-DQA1*05:08-HLA-DQB1*04:08',
         'HLA-DQA1*05:08-HLA-DQB1*05:01', 'HLA-DQA1*05:08-HLA-DQB1*05:02',
         'HLA-DQA1*05:08-HLA-DQB1*05:03', 'HLA-DQA1*05:08-HLA-DQB1*05:05', 'HLA-DQA1*05:08-HLA-DQB1*05:06', 'HLA-DQA1*05:08-HLA-DQB1*05:07',
         'HLA-DQA1*05:08-HLA-DQB1*05:08', 'HLA-DQA1*05:08-HLA-DQB1*05:09',
         'HLA-DQA1*05:08-HLA-DQB1*05:10', 'HLA-DQA1*05:08-HLA-DQB1*05:11', 'HLA-DQA1*05:08-HLA-DQB1*05:12', 'HLA-DQA1*05:08-HLA-DQB1*05:13',
         'HLA-DQA1*05:08-HLA-DQB1*05:14', 'HLA-DQA1*05:08-HLA-DQB1*06:01',
         'HLA-DQA1*05:08-HLA-DQB1*06:02', 'HLA-DQA1*05:08-HLA-DQB1*06:03', 'HLA-DQA1*05:08-HLA-DQB1*06:04', 'HLA-DQA1*05:08-HLA-DQB1*06:07',
         'HLA-DQA1*05:08-HLA-DQB1*06:08', 'HLA-DQA1*05:08-HLA-DQB1*06:09',
         'HLA-DQA1*05:08-HLA-DQB1*06:10', 'HLA-DQA1*05:08-HLA-DQB1*06:11', 'HLA-DQA1*05:08-HLA-DQB1*06:12', 'HLA-DQA1*05:08-HLA-DQB1*06:14',
         'HLA-DQA1*05:08-HLA-DQB1*06:15', 'HLA-DQA1*05:08-HLA-DQB1*06:16',
         'HLA-DQA1*05:08-HLA-DQB1*06:17', 'HLA-DQA1*05:08-HLA-DQB1*06:18', 'HLA-DQA1*05:08-HLA-DQB1*06:19', 'HLA-DQA1*05:08-HLA-DQB1*06:21',
         'HLA-DQA1*05:08-HLA-DQB1*06:22', 'HLA-DQA1*05:08-HLA-DQB1*06:23',
         'HLA-DQA1*05:08-HLA-DQB1*06:24', 'HLA-DQA1*05:08-HLA-DQB1*06:25', 'HLA-DQA1*05:08-HLA-DQB1*06:27', 'HLA-DQA1*05:08-HLA-DQB1*06:28',
         'HLA-DQA1*05:08-HLA-DQB1*06:29', 'HLA-DQA1*05:08-HLA-DQB1*06:30',
         'HLA-DQA1*05:08-HLA-DQB1*06:31', 'HLA-DQA1*05:08-HLA-DQB1*06:32', 'HLA-DQA1*05:08-HLA-DQB1*06:33', 'HLA-DQA1*05:08-HLA-DQB1*06:34',
         'HLA-DQA1*05:08-HLA-DQB1*06:35', 'HLA-DQA1*05:08-HLA-DQB1*06:36',
         'HLA-DQA1*05:08-HLA-DQB1*06:37', 'HLA-DQA1*05:08-HLA-DQB1*06:38', 'HLA-DQA1*05:08-HLA-DQB1*06:39', 'HLA-DQA1*05:08-HLA-DQB1*06:40',
         'HLA-DQA1*05:08-HLA-DQB1*06:41', 'HLA-DQA1*05:08-HLA-DQB1*06:42',
         'HLA-DQA1*05:08-HLA-DQB1*06:43', 'HLA-DQA1*05:08-HLA-DQB1*06:44', 'HLA-DQA1*05:09-HLA-DQB1*02:01', 'HLA-DQA1*05:09-HLA-DQB1*02:02',
         'HLA-DQA1*05:09-HLA-DQB1*02:03', 'HLA-DQA1*05:09-HLA-DQB1*02:04',
         'HLA-DQA1*05:09-HLA-DQB1*02:05', 'HLA-DQA1*05:09-HLA-DQB1*02:06', 'HLA-DQA1*05:09-HLA-DQB1*03:01', 'HLA-DQA1*05:09-HLA-DQB1*03:02',
         'HLA-DQA1*05:09-HLA-DQB1*03:03', 'HLA-DQA1*05:09-HLA-DQB1*03:04',
         'HLA-DQA1*05:09-HLA-DQB1*03:05', 'HLA-DQA1*05:09-HLA-DQB1*03:06', 'HLA-DQA1*05:09-HLA-DQB1*03:07', 'HLA-DQA1*05:09-HLA-DQB1*03:08',
         'HLA-DQA1*05:09-HLA-DQB1*03:09', 'HLA-DQA1*05:09-HLA-DQB1*03:10',
         'HLA-DQA1*05:09-HLA-DQB1*03:11', 'HLA-DQA1*05:09-HLA-DQB1*03:12', 'HLA-DQA1*05:09-HLA-DQB1*03:13', 'HLA-DQA1*05:09-HLA-DQB1*03:14',
         'HLA-DQA1*05:09-HLA-DQB1*03:15', 'HLA-DQA1*05:09-HLA-DQB1*03:16',
         'HLA-DQA1*05:09-HLA-DQB1*03:17', 'HLA-DQA1*05:09-HLA-DQB1*03:18', 'HLA-DQA1*05:09-HLA-DQB1*03:19', 'HLA-DQA1*05:09-HLA-DQB1*03:20',
         'HLA-DQA1*05:09-HLA-DQB1*03:21', 'HLA-DQA1*05:09-HLA-DQB1*03:22',
         'HLA-DQA1*05:09-HLA-DQB1*03:23', 'HLA-DQA1*05:09-HLA-DQB1*03:24', 'HLA-DQA1*05:09-HLA-DQB1*03:25', 'HLA-DQA1*05:09-HLA-DQB1*03:26',
         'HLA-DQA1*05:09-HLA-DQB1*03:27', 'HLA-DQA1*05:09-HLA-DQB1*03:28',
         'HLA-DQA1*05:09-HLA-DQB1*03:29', 'HLA-DQA1*05:09-HLA-DQB1*03:30', 'HLA-DQA1*05:09-HLA-DQB1*03:31', 'HLA-DQA1*05:09-HLA-DQB1*03:32',
         'HLA-DQA1*05:09-HLA-DQB1*03:33', 'HLA-DQA1*05:09-HLA-DQB1*03:34',
         'HLA-DQA1*05:09-HLA-DQB1*03:35', 'HLA-DQA1*05:09-HLA-DQB1*03:36', 'HLA-DQA1*05:09-HLA-DQB1*03:37', 'HLA-DQA1*05:09-HLA-DQB1*03:38',
         'HLA-DQA1*05:09-HLA-DQB1*04:01', 'HLA-DQA1*05:09-HLA-DQB1*04:02',
         'HLA-DQA1*05:09-HLA-DQB1*04:03', 'HLA-DQA1*05:09-HLA-DQB1*04:04', 'HLA-DQA1*05:09-HLA-DQB1*04:05', 'HLA-DQA1*05:09-HLA-DQB1*04:06',
         'HLA-DQA1*05:09-HLA-DQB1*04:07', 'HLA-DQA1*05:09-HLA-DQB1*04:08',
         'HLA-DQA1*05:09-HLA-DQB1*05:01', 'HLA-DQA1*05:09-HLA-DQB1*05:02', 'HLA-DQA1*05:09-HLA-DQB1*05:03', 'HLA-DQA1*05:09-HLA-DQB1*05:05',
         'HLA-DQA1*05:09-HLA-DQB1*05:06', 'HLA-DQA1*05:09-HLA-DQB1*05:07',
         'HLA-DQA1*05:09-HLA-DQB1*05:08', 'HLA-DQA1*05:09-HLA-DQB1*05:09', 'HLA-DQA1*05:09-HLA-DQB1*05:10', 'HLA-DQA1*05:09-HLA-DQB1*05:11',
         'HLA-DQA1*05:09-HLA-DQB1*05:12', 'HLA-DQA1*05:09-HLA-DQB1*05:13',
         'HLA-DQA1*05:09-HLA-DQB1*05:14', 'HLA-DQA1*05:09-HLA-DQB1*06:01', 'HLA-DQA1*05:09-HLA-DQB1*06:02', 'HLA-DQA1*05:09-HLA-DQB1*06:03',
         'HLA-DQA1*05:09-HLA-DQB1*06:04', 'HLA-DQA1*05:09-HLA-DQB1*06:07',
         'HLA-DQA1*05:09-HLA-DQB1*06:08', 'HLA-DQA1*05:09-HLA-DQB1*06:09', 'HLA-DQA1*05:09-HLA-DQB1*06:10', 'HLA-DQA1*05:09-HLA-DQB1*06:11',
         'HLA-DQA1*05:09-HLA-DQB1*06:12', 'HLA-DQA1*05:09-HLA-DQB1*06:14',
         'HLA-DQA1*05:09-HLA-DQB1*06:15', 'HLA-DQA1*05:09-HLA-DQB1*06:16', 'HLA-DQA1*05:09-HLA-DQB1*06:17', 'HLA-DQA1*05:09-HLA-DQB1*06:18',
         'HLA-DQA1*05:09-HLA-DQB1*06:19', 'HLA-DQA1*05:09-HLA-DQB1*06:21',
         'HLA-DQA1*05:09-HLA-DQB1*06:22', 'HLA-DQA1*05:09-HLA-DQB1*06:23', 'HLA-DQA1*05:09-HLA-DQB1*06:24', 'HLA-DQA1*05:09-HLA-DQB1*06:25',
         'HLA-DQA1*05:09-HLA-DQB1*06:27', 'HLA-DQA1*05:09-HLA-DQB1*06:28',
         'HLA-DQA1*05:09-HLA-DQB1*06:29', 'HLA-DQA1*05:09-HLA-DQB1*06:30', 'HLA-DQA1*05:09-HLA-DQB1*06:31', 'HLA-DQA1*05:09-HLA-DQB1*06:32',
         'HLA-DQA1*05:09-HLA-DQB1*06:33', 'HLA-DQA1*05:09-HLA-DQB1*06:34',
         'HLA-DQA1*05:09-HLA-DQB1*06:35', 'HLA-DQA1*05:09-HLA-DQB1*06:36', 'HLA-DQA1*05:09-HLA-DQB1*06:37', 'HLA-DQA1*05:09-HLA-DQB1*06:38',
         'HLA-DQA1*05:09-HLA-DQB1*06:39', 'HLA-DQA1*05:09-HLA-DQB1*06:40',
         'HLA-DQA1*05:09-HLA-DQB1*06:41', 'HLA-DQA1*05:09-HLA-DQB1*06:42', 'HLA-DQA1*05:09-HLA-DQB1*06:43', 'HLA-DQA1*05:09-HLA-DQB1*06:44',
         'HLA-DQA1*05:10-HLA-DQB1*02:01', 'HLA-DQA1*05:10-HLA-DQB1*02:02',
         'HLA-DQA1*05:10-HLA-DQB1*02:03', 'HLA-DQA1*05:10-HLA-DQB1*02:04', 'HLA-DQA1*05:10-HLA-DQB1*02:05', 'HLA-DQA1*05:10-HLA-DQB1*02:06',
         'HLA-DQA1*05:10-HLA-DQB1*03:01', 'HLA-DQA1*05:10-HLA-DQB1*03:02',
         'HLA-DQA1*05:10-HLA-DQB1*03:03', 'HLA-DQA1*05:10-HLA-DQB1*03:04', 'HLA-DQA1*05:10-HLA-DQB1*03:05', 'HLA-DQA1*05:10-HLA-DQB1*03:06',
         'HLA-DQA1*05:10-HLA-DQB1*03:07', 'HLA-DQA1*05:10-HLA-DQB1*03:08',
         'HLA-DQA1*05:10-HLA-DQB1*03:09', 'HLA-DQA1*05:10-HLA-DQB1*03:10', 'HLA-DQA1*05:10-HLA-DQB1*03:11', 'HLA-DQA1*05:10-HLA-DQB1*03:12',
         'HLA-DQA1*05:10-HLA-DQB1*03:13', 'HLA-DQA1*05:10-HLA-DQB1*03:14',
         'HLA-DQA1*05:10-HLA-DQB1*03:15', 'HLA-DQA1*05:10-HLA-DQB1*03:16', 'HLA-DQA1*05:10-HLA-DQB1*03:17', 'HLA-DQA1*05:10-HLA-DQB1*03:18',
         'HLA-DQA1*05:10-HLA-DQB1*03:19', 'HLA-DQA1*05:10-HLA-DQB1*03:20',
         'HLA-DQA1*05:10-HLA-DQB1*03:21', 'HLA-DQA1*05:10-HLA-DQB1*03:22', 'HLA-DQA1*05:10-HLA-DQB1*03:23', 'HLA-DQA1*05:10-HLA-DQB1*03:24',
         'HLA-DQA1*05:10-HLA-DQB1*03:25', 'HLA-DQA1*05:10-HLA-DQB1*03:26',
         'HLA-DQA1*05:10-HLA-DQB1*03:27', 'HLA-DQA1*05:10-HLA-DQB1*03:28', 'HLA-DQA1*05:10-HLA-DQB1*03:29', 'HLA-DQA1*05:10-HLA-DQB1*03:30',
         'HLA-DQA1*05:10-HLA-DQB1*03:31', 'HLA-DQA1*05:10-HLA-DQB1*03:32',
         'HLA-DQA1*05:10-HLA-DQB1*03:33', 'HLA-DQA1*05:10-HLA-DQB1*03:34', 'HLA-DQA1*05:10-HLA-DQB1*03:35', 'HLA-DQA1*05:10-HLA-DQB1*03:36',
         'HLA-DQA1*05:10-HLA-DQB1*03:37', 'HLA-DQA1*05:10-HLA-DQB1*03:38',
         'HLA-DQA1*05:10-HLA-DQB1*04:01', 'HLA-DQA1*05:10-HLA-DQB1*04:02', 'HLA-DQA1*05:10-HLA-DQB1*04:03', 'HLA-DQA1*05:10-HLA-DQB1*04:04',
         'HLA-DQA1*05:10-HLA-DQB1*04:05', 'HLA-DQA1*05:10-HLA-DQB1*04:06',
         'HLA-DQA1*05:10-HLA-DQB1*04:07', 'HLA-DQA1*05:10-HLA-DQB1*04:08', 'HLA-DQA1*05:10-HLA-DQB1*05:01', 'HLA-DQA1*05:10-HLA-DQB1*05:02',
         'HLA-DQA1*05:10-HLA-DQB1*05:03', 'HLA-DQA1*05:10-HLA-DQB1*05:05',
         'HLA-DQA1*05:10-HLA-DQB1*05:06', 'HLA-DQA1*05:10-HLA-DQB1*05:07', 'HLA-DQA1*05:10-HLA-DQB1*05:08', 'HLA-DQA1*05:10-HLA-DQB1*05:09',
         'HLA-DQA1*05:10-HLA-DQB1*05:10', 'HLA-DQA1*05:10-HLA-DQB1*05:11',
         'HLA-DQA1*05:10-HLA-DQB1*05:12', 'HLA-DQA1*05:10-HLA-DQB1*05:13', 'HLA-DQA1*05:10-HLA-DQB1*05:14', 'HLA-DQA1*05:10-HLA-DQB1*06:01',
         'HLA-DQA1*05:10-HLA-DQB1*06:02', 'HLA-DQA1*05:10-HLA-DQB1*06:03',
         'HLA-DQA1*05:10-HLA-DQB1*06:04', 'HLA-DQA1*05:10-HLA-DQB1*06:07', 'HLA-DQA1*05:10-HLA-DQB1*06:08', 'HLA-DQA1*05:10-HLA-DQB1*06:09',
         'HLA-DQA1*05:10-HLA-DQB1*06:10', 'HLA-DQA1*05:10-HLA-DQB1*06:11',
         'HLA-DQA1*05:10-HLA-DQB1*06:12', 'HLA-DQA1*05:10-HLA-DQB1*06:14', 'HLA-DQA1*05:10-HLA-DQB1*06:15', 'HLA-DQA1*05:10-HLA-DQB1*06:16',
         'HLA-DQA1*05:10-HLA-DQB1*06:17', 'HLA-DQA1*05:10-HLA-DQB1*06:18',
         'HLA-DQA1*05:10-HLA-DQB1*06:19', 'HLA-DQA1*05:10-HLA-DQB1*06:21', 'HLA-DQA1*05:10-HLA-DQB1*06:22', 'HLA-DQA1*05:10-HLA-DQB1*06:23',
         'HLA-DQA1*05:10-HLA-DQB1*06:24', 'HLA-DQA1*05:10-HLA-DQB1*06:25',
         'HLA-DQA1*05:10-HLA-DQB1*06:27', 'HLA-DQA1*05:10-HLA-DQB1*06:28', 'HLA-DQA1*05:10-HLA-DQB1*06:29', 'HLA-DQA1*05:10-HLA-DQB1*06:30',
         'HLA-DQA1*05:10-HLA-DQB1*06:31', 'HLA-DQA1*05:10-HLA-DQB1*06:32',
         'HLA-DQA1*05:10-HLA-DQB1*06:33', 'HLA-DQA1*05:10-HLA-DQB1*06:34', 'HLA-DQA1*05:10-HLA-DQB1*06:35', 'HLA-DQA1*05:10-HLA-DQB1*06:36',
         'HLA-DQA1*05:10-HLA-DQB1*06:37', 'HLA-DQA1*05:10-HLA-DQB1*06:38',
         'HLA-DQA1*05:10-HLA-DQB1*06:39', 'HLA-DQA1*05:10-HLA-DQB1*06:40', 'HLA-DQA1*05:10-HLA-DQB1*06:41', 'HLA-DQA1*05:10-HLA-DQB1*06:42',
         'HLA-DQA1*05:10-HLA-DQB1*06:43', 'HLA-DQA1*05:10-HLA-DQB1*06:44',
         'HLA-DQA1*05:11-HLA-DQB1*02:01', 'HLA-DQA1*05:11-HLA-DQB1*02:02', 'HLA-DQA1*05:11-HLA-DQB1*02:03', 'HLA-DQA1*05:11-HLA-DQB1*02:04',
         'HLA-DQA1*05:11-HLA-DQB1*02:05', 'HLA-DQA1*05:11-HLA-DQB1*02:06',
         'HLA-DQA1*05:11-HLA-DQB1*03:01', 'HLA-DQA1*05:11-HLA-DQB1*03:02', 'HLA-DQA1*05:11-HLA-DQB1*03:03', 'HLA-DQA1*05:11-HLA-DQB1*03:04',
         'HLA-DQA1*05:11-HLA-DQB1*03:05', 'HLA-DQA1*05:11-HLA-DQB1*03:06',
         'HLA-DQA1*05:11-HLA-DQB1*03:07', 'HLA-DQA1*05:11-HLA-DQB1*03:08', 'HLA-DQA1*05:11-HLA-DQB1*03:09', 'HLA-DQA1*05:11-HLA-DQB1*03:10',
         'HLA-DQA1*05:11-HLA-DQB1*03:11', 'HLA-DQA1*05:11-HLA-DQB1*03:12',
         'HLA-DQA1*05:11-HLA-DQB1*03:13', 'HLA-DQA1*05:11-HLA-DQB1*03:14', 'HLA-DQA1*05:11-HLA-DQB1*03:15', 'HLA-DQA1*05:11-HLA-DQB1*03:16',
         'HLA-DQA1*05:11-HLA-DQB1*03:17', 'HLA-DQA1*05:11-HLA-DQB1*03:18',
         'HLA-DQA1*05:11-HLA-DQB1*03:19', 'HLA-DQA1*05:11-HLA-DQB1*03:20', 'HLA-DQA1*05:11-HLA-DQB1*03:21', 'HLA-DQA1*05:11-HLA-DQB1*03:22',
         'HLA-DQA1*05:11-HLA-DQB1*03:23', 'HLA-DQA1*05:11-HLA-DQB1*03:24',
         'HLA-DQA1*05:11-HLA-DQB1*03:25', 'HLA-DQA1*05:11-HLA-DQB1*03:26', 'HLA-DQA1*05:11-HLA-DQB1*03:27', 'HLA-DQA1*05:11-HLA-DQB1*03:28',
         'HLA-DQA1*05:11-HLA-DQB1*03:29', 'HLA-DQA1*05:11-HLA-DQB1*03:30',
         'HLA-DQA1*05:11-HLA-DQB1*03:31', 'HLA-DQA1*05:11-HLA-DQB1*03:32', 'HLA-DQA1*05:11-HLA-DQB1*03:33', 'HLA-DQA1*05:11-HLA-DQB1*03:34',
         'HLA-DQA1*05:11-HLA-DQB1*03:35', 'HLA-DQA1*05:11-HLA-DQB1*03:36',
         'HLA-DQA1*05:11-HLA-DQB1*03:37', 'HLA-DQA1*05:11-HLA-DQB1*03:38', 'HLA-DQA1*05:11-HLA-DQB1*04:01', 'HLA-DQA1*05:11-HLA-DQB1*04:02',
         'HLA-DQA1*05:11-HLA-DQB1*04:03', 'HLA-DQA1*05:11-HLA-DQB1*04:04',
         'HLA-DQA1*05:11-HLA-DQB1*04:05', 'HLA-DQA1*05:11-HLA-DQB1*04:06', 'HLA-DQA1*05:11-HLA-DQB1*04:07', 'HLA-DQA1*05:11-HLA-DQB1*04:08',
         'HLA-DQA1*05:11-HLA-DQB1*05:01', 'HLA-DQA1*05:11-HLA-DQB1*05:02',
         'HLA-DQA1*05:11-HLA-DQB1*05:03', 'HLA-DQA1*05:11-HLA-DQB1*05:05', 'HLA-DQA1*05:11-HLA-DQB1*05:06', 'HLA-DQA1*05:11-HLA-DQB1*05:07',
         'HLA-DQA1*05:11-HLA-DQB1*05:08', 'HLA-DQA1*05:11-HLA-DQB1*05:09',
         'HLA-DQA1*05:11-HLA-DQB1*05:10', 'HLA-DQA1*05:11-HLA-DQB1*05:11', 'HLA-DQA1*05:11-HLA-DQB1*05:12', 'HLA-DQA1*05:11-HLA-DQB1*05:13',
         'HLA-DQA1*05:11-HLA-DQB1*05:14', 'HLA-DQA1*05:11-HLA-DQB1*06:01',
         'HLA-DQA1*05:11-HLA-DQB1*06:02', 'HLA-DQA1*05:11-HLA-DQB1*06:03', 'HLA-DQA1*05:11-HLA-DQB1*06:04', 'HLA-DQA1*05:11-HLA-DQB1*06:07',
         'HLA-DQA1*05:11-HLA-DQB1*06:08', 'HLA-DQA1*05:11-HLA-DQB1*06:09',
         'HLA-DQA1*05:11-HLA-DQB1*06:10', 'HLA-DQA1*05:11-HLA-DQB1*06:11', 'HLA-DQA1*05:11-HLA-DQB1*06:12', 'HLA-DQA1*05:11-HLA-DQB1*06:14',
         'HLA-DQA1*05:11-HLA-DQB1*06:15', 'HLA-DQA1*05:11-HLA-DQB1*06:16',
         'HLA-DQA1*05:11-HLA-DQB1*06:17', 'HLA-DQA1*05:11-HLA-DQB1*06:18', 'HLA-DQA1*05:11-HLA-DQB1*06:19', 'HLA-DQA1*05:11-HLA-DQB1*06:21',
         'HLA-DQA1*05:11-HLA-DQB1*06:22', 'HLA-DQA1*05:11-HLA-DQB1*06:23',
         'HLA-DQA1*05:11-HLA-DQB1*06:24', 'HLA-DQA1*05:11-HLA-DQB1*06:25', 'HLA-DQA1*05:11-HLA-DQB1*06:27', 'HLA-DQA1*05:11-HLA-DQB1*06:28',
         'HLA-DQA1*05:11-HLA-DQB1*06:29', 'HLA-DQA1*05:11-HLA-DQB1*06:30',
         'HLA-DQA1*05:11-HLA-DQB1*06:31', 'HLA-DQA1*05:11-HLA-DQB1*06:32', 'HLA-DQA1*05:11-HLA-DQB1*06:33', 'HLA-DQA1*05:11-HLA-DQB1*06:34',
         'HLA-DQA1*05:11-HLA-DQB1*06:35', 'HLA-DQA1*05:11-HLA-DQB1*06:36',
         'HLA-DQA1*05:11-HLA-DQB1*06:37', 'HLA-DQA1*05:11-HLA-DQB1*06:38', 'HLA-DQA1*05:11-HLA-DQB1*06:39', 'HLA-DQA1*05:11-HLA-DQB1*06:40',
         'HLA-DQA1*05:11-HLA-DQB1*06:41', 'HLA-DQA1*05:11-HLA-DQB1*06:42',
         'HLA-DQA1*05:11-HLA-DQB1*06:43', 'HLA-DQA1*05:11-HLA-DQB1*06:44', 'HLA-DQA1*06:01-HLA-DQB1*02:01', 'HLA-DQA1*06:01-HLA-DQB1*02:02',
         'HLA-DQA1*06:01-HLA-DQB1*02:03', 'HLA-DQA1*06:01-HLA-DQB1*02:04',
         'HLA-DQA1*06:01-HLA-DQB1*02:05', 'HLA-DQA1*06:01-HLA-DQB1*02:06', 'HLA-DQA1*06:01-HLA-DQB1*03:01', 'HLA-DQA1*06:01-HLA-DQB1*03:02',
         'HLA-DQA1*06:01-HLA-DQB1*03:03', 'HLA-DQA1*06:01-HLA-DQB1*03:04',
         'HLA-DQA1*06:01-HLA-DQB1*03:05', 'HLA-DQA1*06:01-HLA-DQB1*03:06', 'HLA-DQA1*06:01-HLA-DQB1*03:07', 'HLA-DQA1*06:01-HLA-DQB1*03:08',
         'HLA-DQA1*06:01-HLA-DQB1*03:09', 'HLA-DQA1*06:01-HLA-DQB1*03:10',
         'HLA-DQA1*06:01-HLA-DQB1*03:11', 'HLA-DQA1*06:01-HLA-DQB1*03:12', 'HLA-DQA1*06:01-HLA-DQB1*03:13', 'HLA-DQA1*06:01-HLA-DQB1*03:14',
         'HLA-DQA1*06:01-HLA-DQB1*03:15', 'HLA-DQA1*06:01-HLA-DQB1*03:16',
         'HLA-DQA1*06:01-HLA-DQB1*03:17', 'HLA-DQA1*06:01-HLA-DQB1*03:18', 'HLA-DQA1*06:01-HLA-DQB1*03:19', 'HLA-DQA1*06:01-HLA-DQB1*03:20',
         'HLA-DQA1*06:01-HLA-DQB1*03:21', 'HLA-DQA1*06:01-HLA-DQB1*03:22',
         'HLA-DQA1*06:01-HLA-DQB1*03:23', 'HLA-DQA1*06:01-HLA-DQB1*03:24', 'HLA-DQA1*06:01-HLA-DQB1*03:25', 'HLA-DQA1*06:01-HLA-DQB1*03:26',
         'HLA-DQA1*06:01-HLA-DQB1*03:27', 'HLA-DQA1*06:01-HLA-DQB1*03:28',
         'HLA-DQA1*06:01-HLA-DQB1*03:29', 'HLA-DQA1*06:01-HLA-DQB1*03:30', 'HLA-DQA1*06:01-HLA-DQB1*03:31', 'HLA-DQA1*06:01-HLA-DQB1*03:32',
         'HLA-DQA1*06:01-HLA-DQB1*03:33', 'HLA-DQA1*06:01-HLA-DQB1*03:34',
         'HLA-DQA1*06:01-HLA-DQB1*03:35', 'HLA-DQA1*06:01-HLA-DQB1*03:36', 'HLA-DQA1*06:01-HLA-DQB1*03:37', 'HLA-DQA1*06:01-HLA-DQB1*03:38',
         'HLA-DQA1*06:01-HLA-DQB1*04:01', 'HLA-DQA1*06:01-HLA-DQB1*04:02',
         'HLA-DQA1*06:01-HLA-DQB1*04:03', 'HLA-DQA1*06:01-HLA-DQB1*04:04', 'HLA-DQA1*06:01-HLA-DQB1*04:05', 'HLA-DQA1*06:01-HLA-DQB1*04:06',
         'HLA-DQA1*06:01-HLA-DQB1*04:07', 'HLA-DQA1*06:01-HLA-DQB1*04:08',
         'HLA-DQA1*06:01-HLA-DQB1*05:01', 'HLA-DQA1*06:01-HLA-DQB1*05:02', 'HLA-DQA1*06:01-HLA-DQB1*05:03', 'HLA-DQA1*06:01-HLA-DQB1*05:05',
         'HLA-DQA1*06:01-HLA-DQB1*05:06', 'HLA-DQA1*06:01-HLA-DQB1*05:07',
         'HLA-DQA1*06:01-HLA-DQB1*05:08', 'HLA-DQA1*06:01-HLA-DQB1*05:09', 'HLA-DQA1*06:01-HLA-DQB1*05:10', 'HLA-DQA1*06:01-HLA-DQB1*05:11',
         'HLA-DQA1*06:01-HLA-DQB1*05:12', 'HLA-DQA1*06:01-HLA-DQB1*05:13',
         'HLA-DQA1*06:01-HLA-DQB1*05:14', 'HLA-DQA1*06:01-HLA-DQB1*06:01', 'HLA-DQA1*06:01-HLA-DQB1*06:02', 'HLA-DQA1*06:01-HLA-DQB1*06:03',
         'HLA-DQA1*06:01-HLA-DQB1*06:04', 'HLA-DQA1*06:01-HLA-DQB1*06:07',
         'HLA-DQA1*06:01-HLA-DQB1*06:08', 'HLA-DQA1*06:01-HLA-DQB1*06:09', 'HLA-DQA1*06:01-HLA-DQB1*06:10', 'HLA-DQA1*06:01-HLA-DQB1*06:11',
         'HLA-DQA1*06:01-HLA-DQB1*06:12', 'HLA-DQA1*06:01-HLA-DQB1*06:14',
         'HLA-DQA1*06:01-HLA-DQB1*06:15', 'HLA-DQA1*06:01-HLA-DQB1*06:16', 'HLA-DQA1*06:01-HLA-DQB1*06:17', 'HLA-DQA1*06:01-HLA-DQB1*06:18',
         'HLA-DQA1*06:01-HLA-DQB1*06:19', 'HLA-DQA1*06:01-HLA-DQB1*06:21',
         'HLA-DQA1*06:01-HLA-DQB1*06:22', 'HLA-DQA1*06:01-HLA-DQB1*06:23', 'HLA-DQA1*06:01-HLA-DQB1*06:24', 'HLA-DQA1*06:01-HLA-DQB1*06:25',
         'HLA-DQA1*06:01-HLA-DQB1*06:27', 'HLA-DQA1*06:01-HLA-DQB1*06:28',
         'HLA-DQA1*06:01-HLA-DQB1*06:29', 'HLA-DQA1*06:01-HLA-DQB1*06:30', 'HLA-DQA1*06:01-HLA-DQB1*06:31', 'HLA-DQA1*06:01-HLA-DQB1*06:32',
         'HLA-DQA1*06:01-HLA-DQB1*06:33', 'HLA-DQA1*06:01-HLA-DQB1*06:34',
         'HLA-DQA1*06:01-HLA-DQB1*06:35', 'HLA-DQA1*06:01-HLA-DQB1*06:36', 'HLA-DQA1*06:01-HLA-DQB1*06:37', 'HLA-DQA1*06:01-HLA-DQB1*06:38',
         'HLA-DQA1*06:01-HLA-DQB1*06:39', 'HLA-DQA1*06:01-HLA-DQB1*06:40',
         'HLA-DQA1*06:01-HLA-DQB1*06:41', 'HLA-DQA1*06:01-HLA-DQB1*06:42', 'HLA-DQA1*06:01-HLA-DQB1*06:43', 'HLA-DQA1*06:01-HLA-DQB1*06:44',
         'HLA-DQA1*06:02-HLA-DQB1*02:01', 'HLA-DQA1*06:02-HLA-DQB1*02:02',
         'HLA-DQA1*06:02-HLA-DQB1*02:03', 'HLA-DQA1*06:02-HLA-DQB1*02:04', 'HLA-DQA1*06:02-HLA-DQB1*02:05', 'HLA-DQA1*06:02-HLA-DQB1*02:06',
         'HLA-DQA1*06:02-HLA-DQB1*03:01', 'HLA-DQA1*06:02-HLA-DQB1*03:02',
         'HLA-DQA1*06:02-HLA-DQB1*03:03', 'HLA-DQA1*06:02-HLA-DQB1*03:04', 'HLA-DQA1*06:02-HLA-DQB1*03:05', 'HLA-DQA1*06:02-HLA-DQB1*03:06',
         'HLA-DQA1*06:02-HLA-DQB1*03:07', 'HLA-DQA1*06:02-HLA-DQB1*03:08',
         'HLA-DQA1*06:02-HLA-DQB1*03:09', 'HLA-DQA1*06:02-HLA-DQB1*03:10', 'HLA-DQA1*06:02-HLA-DQB1*03:11', 'HLA-DQA1*06:02-HLA-DQB1*03:12',
         'HLA-DQA1*06:02-HLA-DQB1*03:13', 'HLA-DQA1*06:02-HLA-DQB1*03:14',
         'HLA-DQA1*06:02-HLA-DQB1*03:15', 'HLA-DQA1*06:02-HLA-DQB1*03:16', 'HLA-DQA1*06:02-HLA-DQB1*03:17', 'HLA-DQA1*06:02-HLA-DQB1*03:18',
         'HLA-DQA1*06:02-HLA-DQB1*03:19', 'HLA-DQA1*06:02-HLA-DQB1*03:20',
         'HLA-DQA1*06:02-HLA-DQB1*03:21', 'HLA-DQA1*06:02-HLA-DQB1*03:22', 'HLA-DQA1*06:02-HLA-DQB1*03:23', 'HLA-DQA1*06:02-HLA-DQB1*03:24',
         'HLA-DQA1*06:02-HLA-DQB1*03:25', 'HLA-DQA1*06:02-HLA-DQB1*03:26',
         'HLA-DQA1*06:02-HLA-DQB1*03:27', 'HLA-DQA1*06:02-HLA-DQB1*03:28', 'HLA-DQA1*06:02-HLA-DQB1*03:29', 'HLA-DQA1*06:02-HLA-DQB1*03:30',
         'HLA-DQA1*06:02-HLA-DQB1*03:31', 'HLA-DQA1*06:02-HLA-DQB1*03:32',
         'HLA-DQA1*06:02-HLA-DQB1*03:33', 'HLA-DQA1*06:02-HLA-DQB1*03:34', 'HLA-DQA1*06:02-HLA-DQB1*03:35', 'HLA-DQA1*06:02-HLA-DQB1*03:36',
         'HLA-DQA1*06:02-HLA-DQB1*03:37', 'HLA-DQA1*06:02-HLA-DQB1*03:38',
         'HLA-DQA1*06:02-HLA-DQB1*04:01', 'HLA-DQA1*06:02-HLA-DQB1*04:02', 'HLA-DQA1*06:02-HLA-DQB1*04:03', 'HLA-DQA1*06:02-HLA-DQB1*04:04',
         'HLA-DQA1*06:02-HLA-DQB1*04:05', 'HLA-DQA1*06:02-HLA-DQB1*04:06',
         'HLA-DQA1*06:02-HLA-DQB1*04:07', 'HLA-DQA1*06:02-HLA-DQB1*04:08', 'HLA-DQA1*06:02-HLA-DQB1*05:01', 'HLA-DQA1*06:02-HLA-DQB1*05:02',
         'HLA-DQA1*06:02-HLA-DQB1*05:03', 'HLA-DQA1*06:02-HLA-DQB1*05:05',
         'HLA-DQA1*06:02-HLA-DQB1*05:06', 'HLA-DQA1*06:02-HLA-DQB1*05:07', 'HLA-DQA1*06:02-HLA-DQB1*05:08', 'HLA-DQA1*06:02-HLA-DQB1*05:09',
         'HLA-DQA1*06:02-HLA-DQB1*05:10', 'HLA-DQA1*06:02-HLA-DQB1*05:11',
         'HLA-DQA1*06:02-HLA-DQB1*05:12', 'HLA-DQA1*06:02-HLA-DQB1*05:13', 'HLA-DQA1*06:02-HLA-DQB1*05:14', 'HLA-DQA1*06:02-HLA-DQB1*06:01',
         'HLA-DQA1*06:02-HLA-DQB1*06:02', 'HLA-DQA1*06:02-HLA-DQB1*06:03',
         'HLA-DQA1*06:02-HLA-DQB1*06:04', 'HLA-DQA1*06:02-HLA-DQB1*06:07', 'HLA-DQA1*06:02-HLA-DQB1*06:08', 'HLA-DQA1*06:02-HLA-DQB1*06:09',
         'HLA-DQA1*06:02-HLA-DQB1*06:10', 'HLA-DQA1*06:02-HLA-DQB1*06:11',
         'HLA-DQA1*06:02-HLA-DQB1*06:12', 'HLA-DQA1*06:02-HLA-DQB1*06:14', 'HLA-DQA1*06:02-HLA-DQB1*06:15', 'HLA-DQA1*06:02-HLA-DQB1*06:16',
         'HLA-DQA1*06:02-HLA-DQB1*06:17', 'HLA-DQA1*06:02-HLA-DQB1*06:18',
         'HLA-DQA1*06:02-HLA-DQB1*06:19', 'HLA-DQA1*06:02-HLA-DQB1*06:21', 'HLA-DQA1*06:02-HLA-DQB1*06:22', 'HLA-DQA1*06:02-HLA-DQB1*06:23',
         'HLA-DQA1*06:02-HLA-DQB1*06:24', 'HLA-DQA1*06:02-HLA-DQB1*06:25',
         'HLA-DQA1*06:02-HLA-DQB1*06:27', 'HLA-DQA1*06:02-HLA-DQB1*06:28', 'HLA-DQA1*06:02-HLA-DQB1*06:29', 'HLA-DQA1*06:02-HLA-DQB1*06:30',
         'HLA-DQA1*06:02-HLA-DQB1*06:31', 'HLA-DQA1*06:02-HLA-DQB1*06:32',
         'HLA-DQA1*06:02-HLA-DQB1*06:33', 'HLA-DQA1*06:02-HLA-DQB1*06:34', 'HLA-DQA1*06:02-HLA-DQB1*06:35', 'HLA-DQA1*06:02-HLA-DQB1*06:36',
         'HLA-DQA1*06:02-HLA-DQB1*06:37', 'HLA-DQA1*06:02-HLA-DQB1*06:38',
         'HLA-DQA1*06:02-HLA-DQB1*06:39', 'HLA-DQA1*06:02-HLA-DQB1*06:40', 'HLA-DQA1*06:02-HLA-DQB1*06:41', 'HLA-DQA1*06:02-HLA-DQB1*06:42',
         'HLA-DQA1*06:02-HLA-DQB1*06:43', 'HLA-DQA1*06:02-HLA-DQB1*06:44',
         'H2-IAd', 'H2-IAb'])

    __version = "3.1"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def command(self):
        """
        Defines the commandline call for external tool
        """
        return self.__command

    @property
    def supportedLength(self):
        """
        A list of supported :class:`~Fred2.Core.Peptide.Peptide` lengths
        """
        return self.__supported_length

    @property
    def supportedAlleles(self):
        """A list of valid :class:`~Fred2.Core.Allele.Allele` models"""
        return self.__alleles

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name


class PickPocket_1_1(AExternalEpitopePrediction):
    """
    Implementation of PickPocket adapter.

    .. note::

        Zhang, H., Lund, O., & Nielsen, M. (2009). The PickPocket method for predicting binding specificities
        for receptors based on receptor pocket similarities: application to MHC-peptide binding.
        Bioinformatics, 25(10), 1293-1299.

    """
    __name = "pickpocket"
    __supported_length = frozenset([8, 9, 10, 11])
    __command = 'PickPocket -p {peptides} -a {alleles} {options} | grep -v "#" > {out}'
    __supported_alleles = frozenset(['HLA-A*01:01', 'HLA-A*01:02', 'HLA-A*01:03', 'HLA-A*01:06', 'HLA-A*01:07', 'HLA-A*01:08', 'HLA-A*01:09',
                                     'HLA-A*01:10', 'HLA-A*01:12', 'HLA-A*01:13', 'HLA-A*01:14', 'HLA-A*01:17', 'HLA-A*01:19', 'HLA-A*01:20',
                                     'HLA-A*01:21', 'HLA-A*01:23', 'HLA-A*01:24',
                                     'HLA-A*01:25', 'HLA-A*01:26', 'HLA-A*01:28', 'HLA-A*01:29', 'HLA-A*01:30', 'HLA-A*01:32', 'HLA-A*01:33',
                                     'HLA-A*01:35', 'HLA-A*01:36', 'HLA-A*01:37',
                                     'HLA-A*01:38', 'HLA-A*01:39', 'HLA-A*01:40', 'HLA-A*01:41', 'HLA-A*01:42', 'HLA-A*01:43', 'HLA-A*01:44',
                                     'HLA-A*01:45', 'HLA-A*01:46', 'HLA-A*01:47',
                                     'HLA-A*01:48', 'HLA-A*01:49', 'HLA-A*01:50', 'HLA-A*01:51', 'HLA-A*01:54', 'HLA-A*01:55', 'HLA-A*01:58',
                                     'HLA-A*01:59', 'HLA-A*01:60', 'HLA-A*01:61',
                                     'HLA-A*01:62', 'HLA-A*01:63', 'HLA-A*01:64', 'HLA-A*01:65', 'HLA-A*01:66', 'HLA-A*02:01', 'HLA-A*02:02',
                                     'HLA-A*02:03', 'HLA-A*02:04', 'HLA-A*02:05',
                                     'HLA-A*02:06', 'HLA-A*02:07', 'HLA-A*02:08', 'HLA-A*02:09', 'HLA-A*02:10', 'HLA-A*02:11', 'HLA-A*02:12',
                                     'HLA-A*02:13', 'HLA-A*02:14', 'HLA-A*02:16',
                                     'HLA-A*02:17', 'HLA-A*02:18', 'HLA-A*02:19', 'HLA-A*02:20', 'HLA-A*02:21', 'HLA-A*02:22', 'HLA-A*02:24',
                                     'HLA-A*02:25', 'HLA-A*02:26', 'HLA-A*02:27',
                                     'HLA-A*02:28', 'HLA-A*02:29', 'HLA-A*02:30', 'HLA-A*02:31', 'HLA-A*02:33', 'HLA-A*02:34', 'HLA-A*02:35',
                                     'HLA-A*02:36', 'HLA-A*02:37', 'HLA-A*02:38',
                                     'HLA-A*02:39', 'HLA-A*02:40', 'HLA-A*02:41', 'HLA-A*02:42', 'HLA-A*02:44', 'HLA-A*02:45', 'HLA-A*02:46',
                                     'HLA-A*02:47', 'HLA-A*02:48', 'HLA-A*02:49',
                                     'HLA-A*02:50', 'HLA-A*02:51', 'HLA-A*02:52', 'HLA-A*02:54', 'HLA-A*02:55', 'HLA-A*02:56', 'HLA-A*02:57',
                                     'HLA-A*02:58', 'HLA-A*02:59', 'HLA-A*02:60',
                                     'HLA-A*02:61', 'HLA-A*02:62', 'HLA-A*02:63', 'HLA-A*02:64', 'HLA-A*02:65', 'HLA-A*02:66', 'HLA-A*02:67',
                                     'HLA-A*02:68', 'HLA-A*02:69', 'HLA-A*02:70',
                                     'HLA-A*02:71', 'HLA-A*02:72', 'HLA-A*02:73', 'HLA-A*02:74', 'HLA-A*02:75', 'HLA-A*02:76', 'HLA-A*02:77',
                                     'HLA-A*02:78', 'HLA-A*02:79', 'HLA-A*02:80',
                                     'HLA-A*02:81', 'HLA-A*02:84', 'HLA-A*02:85', 'HLA-A*02:86', 'HLA-A*02:87', 'HLA-A*02:89', 'HLA-A*02:90',
                                     'HLA-A*02:91', 'HLA-A*02:92', 'HLA-A*02:93',
                                     'HLA-A*02:95', 'HLA-A*02:96', 'HLA-A*02:97', 'HLA-A*02:99', 'HLA-A*02:101', 'HLA-A*02:102', 'HLA-A*02:103',
                                     'HLA-A*02:104', 'HLA-A*02:105',
                                     'HLA-A*02:106', 'HLA-A*02:107', 'HLA-A*02:108', 'HLA-A*02:109', 'HLA-A*02:110', 'HLA-A*02:111', 'HLA-A*02:112',
                                     'HLA-A*02:114', 'HLA-A*02:115',
                                     'HLA-A*02:116', 'HLA-A*02:117', 'HLA-A*02:118', 'HLA-A*02:119', 'HLA-A*02:120', 'HLA-A*02:121', 'HLA-A*02:122',
                                     'HLA-A*02:123', 'HLA-A*02:124',
                                     'HLA-A*02:126', 'HLA-A*02:127', 'HLA-A*02:128', 'HLA-A*02:129', 'HLA-A*02:130', 'HLA-A*02:131', 'HLA-A*02:132',
                                     'HLA-A*02:133', 'HLA-A*02:134',
                                     'HLA-A*02:135', 'HLA-A*02:136', 'HLA-A*02:137', 'HLA-A*02:138', 'HLA-A*02:139', 'HLA-A*02:140', 'HLA-A*02:141',
                                     'HLA-A*02:142', 'HLA-A*02:143',
                                     'HLA-A*02:144', 'HLA-A*02:145', 'HLA-A*02:146', 'HLA-A*02:147', 'HLA-A*02:148', 'HLA-A*02:149', 'HLA-A*02:150',
                                     'HLA-A*02:151', 'HLA-A*02:152',
                                     'HLA-A*02:153', 'HLA-A*02:154', 'HLA-A*02:155', 'HLA-A*02:156', 'HLA-A*02:157', 'HLA-A*02:158', 'HLA-A*02:159',
                                     'HLA-A*02:160', 'HLA-A*02:161',
                                     'HLA-A*02:162', 'HLA-A*02:163', 'HLA-A*02:164', 'HLA-A*02:165', 'HLA-A*02:166', 'HLA-A*02:167', 'HLA-A*02:168',
                                     'HLA-A*02:169', 'HLA-A*02:170',
                                     'HLA-A*02:171', 'HLA-A*02:172', 'HLA-A*02:173', 'HLA-A*02:174', 'HLA-A*02:175', 'HLA-A*02:176', 'HLA-A*02:177',
                                     'HLA-A*02:178', 'HLA-A*02:179',
                                     'HLA-A*02:180', 'HLA-A*02:181', 'HLA-A*02:182', 'HLA-A*02:183', 'HLA-A*02:184', 'HLA-A*02:185', 'HLA-A*02:186',
                                     'HLA-A*02:187', 'HLA-A*02:188',
                                     'HLA-A*02:189', 'HLA-A*02:190', 'HLA-A*02:191', 'HLA-A*02:192', 'HLA-A*02:193', 'HLA-A*02:194', 'HLA-A*02:195',
                                     'HLA-A*02:196', 'HLA-A*02:197',
                                     'HLA-A*02:198', 'HLA-A*02:199', 'HLA-A*02:200', 'HLA-A*02:201', 'HLA-A*02:202', 'HLA-A*02:203', 'HLA-A*02:204',
                                     'HLA-A*02:205', 'HLA-A*02:206',
                                     'HLA-A*02:207', 'HLA-A*02:208', 'HLA-A*02:209', 'HLA-A*02:210', 'HLA-A*02:211', 'HLA-A*02:212', 'HLA-A*02:213',
                                     'HLA-A*02:214', 'HLA-A*02:215',
                                     'HLA-A*02:216', 'HLA-A*02:217', 'HLA-A*02:218', 'HLA-A*02:219', 'HLA-A*02:220', 'HLA-A*02:221', 'HLA-A*02:224',
                                     'HLA-A*02:228', 'HLA-A*02:229',
                                     'HLA-A*02:230', 'HLA-A*02:231', 'HLA-A*02:232', 'HLA-A*02:233', 'HLA-A*02:234', 'HLA-A*02:235', 'HLA-A*02:236',
                                     'HLA-A*02:237', 'HLA-A*02:238',
                                     'HLA-A*02:239', 'HLA-A*02:240', 'HLA-A*02:241', 'HLA-A*02:242', 'HLA-A*02:243', 'HLA-A*02:244', 'HLA-A*02:245',
                                     'HLA-A*02:246', 'HLA-A*02:247',
                                     'HLA-A*02:248', 'HLA-A*02:249', 'HLA-A*02:251', 'HLA-A*02:252', 'HLA-A*02:253', 'HLA-A*02:254', 'HLA-A*02:255',
                                     'HLA-A*02:256', 'HLA-A*02:257',
                                     'HLA-A*02:258', 'HLA-A*02:259', 'HLA-A*02:260', 'HLA-A*02:261', 'HLA-A*02:262', 'HLA-A*02:263', 'HLA-A*02:264',
                                     'HLA-A*02:265', 'HLA-A*02:266',
                                     'HLA-A*03:01', 'HLA-A*03:02', 'HLA-A*03:04', 'HLA-A*03:05', 'HLA-A*03:06', 'HLA-A*03:07', 'HLA-A*03:08',
                                     'HLA-A*03:09', 'HLA-A*03:10', 'HLA-A*03:12',
                                     'HLA-A*03:13', 'HLA-A*03:14', 'HLA-A*03:15', 'HLA-A*03:16', 'HLA-A*03:17', 'HLA-A*03:18', 'HLA-A*03:19',
                                     'HLA-A*03:20', 'HLA-A*03:22', 'HLA-A*03:23',
                                     'HLA-A*03:24', 'HLA-A*03:25', 'HLA-A*03:26', 'HLA-A*03:27', 'HLA-A*03:28', 'HLA-A*03:29', 'HLA-A*03:30',
                                     'HLA-A*03:31', 'HLA-A*03:32', 'HLA-A*03:33',
                                     'HLA-A*03:34', 'HLA-A*03:35', 'HLA-A*03:37', 'HLA-A*03:38', 'HLA-A*03:39', 'HLA-A*03:40', 'HLA-A*03:41',
                                     'HLA-A*03:42', 'HLA-A*03:43', 'HLA-A*03:44',
                                     'HLA-A*03:45', 'HLA-A*03:46', 'HLA-A*03:47', 'HLA-A*03:48', 'HLA-A*03:49', 'HLA-A*03:50', 'HLA-A*03:51',
                                     'HLA-A*03:52', 'HLA-A*03:53', 'HLA-A*03:54',
                                     'HLA-A*03:55', 'HLA-A*03:56', 'HLA-A*03:57', 'HLA-A*03:58', 'HLA-A*03:59', 'HLA-A*03:60', 'HLA-A*03:61',
                                     'HLA-A*03:62', 'HLA-A*03:63', 'HLA-A*03:64',
                                     'HLA-A*03:65', 'HLA-A*03:66', 'HLA-A*03:67', 'HLA-A*03:70', 'HLA-A*03:71', 'HLA-A*03:72', 'HLA-A*03:73',
                                     'HLA-A*03:74', 'HLA-A*03:75', 'HLA-A*03:76',
                                     'HLA-A*03:77', 'HLA-A*03:78', 'HLA-A*03:79', 'HLA-A*03:80', 'HLA-A*03:81', 'HLA-A*03:82', 'HLA-A*11:01',
                                     'HLA-A*11:02', 'HLA-A*11:03', 'HLA-A*11:04',
                                     'HLA-A*11:05', 'HLA-A*11:06', 'HLA-A*11:07', 'HLA-A*11:08', 'HLA-A*11:09', 'HLA-A*11:10', 'HLA-A*11:11',
                                     'HLA-A*11:12', 'HLA-A*11:13', 'HLA-A*11:14',
                                     'HLA-A*11:15', 'HLA-A*11:16', 'HLA-A*11:17', 'HLA-A*11:18', 'HLA-A*11:19', 'HLA-A*11:20', 'HLA-A*11:22',
                                     'HLA-A*11:23', 'HLA-A*11:24', 'HLA-A*11:25',
                                     'HLA-A*11:26', 'HLA-A*11:27', 'HLA-A*11:29', 'HLA-A*11:30', 'HLA-A*11:31', 'HLA-A*11:32', 'HLA-A*11:33',
                                     'HLA-A*11:34', 'HLA-A*11:35', 'HLA-A*11:36',
                                     'HLA-A*11:37', 'HLA-A*11:38', 'HLA-A*11:39', 'HLA-A*11:40', 'HLA-A*11:41', 'HLA-A*11:42', 'HLA-A*11:43',
                                     'HLA-A*11:44', 'HLA-A*11:45', 'HLA-A*11:46',
                                     'HLA-A*11:47', 'HLA-A*11:48', 'HLA-A*11:49', 'HLA-A*11:51', 'HLA-A*11:53', 'HLA-A*11:54', 'HLA-A*11:55',
                                     'HLA-A*11:56', 'HLA-A*11:57', 'HLA-A*11:58',
                                     'HLA-A*11:59', 'HLA-A*11:60', 'HLA-A*11:61', 'HLA-A*11:62', 'HLA-A*11:63', 'HLA-A*11:64', 'HLA-A*23:01',
                                     'HLA-A*23:02', 'HLA-A*23:03', 'HLA-A*23:04',
                                     'HLA-A*23:05', 'HLA-A*23:06', 'HLA-A*23:09', 'HLA-A*23:10', 'HLA-A*23:12', 'HLA-A*23:13', 'HLA-A*23:14',
                                     'HLA-A*23:15', 'HLA-A*23:16', 'HLA-A*23:17',
                                     'HLA-A*23:18', 'HLA-A*23:20', 'HLA-A*23:21', 'HLA-A*23:22', 'HLA-A*23:23', 'HLA-A*23:24', 'HLA-A*23:25',
                                     'HLA-A*23:26', 'HLA-A*24:02', 'HLA-A*24:03',
                                     'HLA-A*24:04', 'HLA-A*24:05', 'HLA-A*24:06', 'HLA-A*24:07', 'HLA-A*24:08', 'HLA-A*24:10', 'HLA-A*24:13',
                                     'HLA-A*24:14', 'HLA-A*24:15', 'HLA-A*24:17',
                                     'HLA-A*24:18', 'HLA-A*24:19', 'HLA-A*24:20', 'HLA-A*24:21', 'HLA-A*24:22', 'HLA-A*24:23', 'HLA-A*24:24',
                                     'HLA-A*24:25', 'HLA-A*24:26', 'HLA-A*24:27',
                                     'HLA-A*24:28', 'HLA-A*24:29', 'HLA-A*24:30', 'HLA-A*24:31', 'HLA-A*24:32', 'HLA-A*24:33', 'HLA-A*24:34',
                                     'HLA-A*24:35', 'HLA-A*24:37', 'HLA-A*24:38',
                                     'HLA-A*24:39', 'HLA-A*24:41', 'HLA-A*24:42', 'HLA-A*24:43', 'HLA-A*24:44', 'HLA-A*24:46', 'HLA-A*24:47',
                                     'HLA-A*24:49', 'HLA-A*24:50', 'HLA-A*24:51',
                                     'HLA-A*24:52', 'HLA-A*24:53', 'HLA-A*24:54', 'HLA-A*24:55', 'HLA-A*24:56', 'HLA-A*24:57', 'HLA-A*24:58',
                                     'HLA-A*24:59', 'HLA-A*24:61', 'HLA-A*24:62',
                                     'HLA-A*24:63', 'HLA-A*24:64', 'HLA-A*24:66', 'HLA-A*24:67', 'HLA-A*24:68', 'HLA-A*24:69', 'HLA-A*24:70',
                                     'HLA-A*24:71', 'HLA-A*24:72', 'HLA-A*24:73',
                                     'HLA-A*24:74', 'HLA-A*24:75', 'HLA-A*24:76', 'HLA-A*24:77', 'HLA-A*24:78', 'HLA-A*24:79', 'HLA-A*24:80',
                                     'HLA-A*24:81', 'HLA-A*24:82', 'HLA-A*24:85',
                                     'HLA-A*24:87', 'HLA-A*24:88', 'HLA-A*24:89', 'HLA-A*24:91', 'HLA-A*24:92', 'HLA-A*24:93', 'HLA-A*24:94',
                                     'HLA-A*24:95', 'HLA-A*24:96', 'HLA-A*24:97',
                                     'HLA-A*24:98', 'HLA-A*24:99', 'HLA-A*24:100', 'HLA-A*24:101', 'HLA-A*24:102', 'HLA-A*24:103', 'HLA-A*24:104',
                                     'HLA-A*24:105', 'HLA-A*24:106',
                                     'HLA-A*24:107', 'HLA-A*24:108', 'HLA-A*24:109', 'HLA-A*24:110', 'HLA-A*24:111', 'HLA-A*24:112', 'HLA-A*24:113',
                                     'HLA-A*24:114', 'HLA-A*24:115',
                                     'HLA-A*24:116', 'HLA-A*24:117', 'HLA-A*24:118', 'HLA-A*24:119', 'HLA-A*24:120', 'HLA-A*24:121', 'HLA-A*24:122',
                                     'HLA-A*24:123', 'HLA-A*24:124',
                                     'HLA-A*24:125', 'HLA-A*24:126', 'HLA-A*24:127', 'HLA-A*24:128', 'HLA-A*24:129', 'HLA-A*24:130', 'HLA-A*24:131',
                                     'HLA-A*24:133', 'HLA-A*24:134',
                                     'HLA-A*24:135', 'HLA-A*24:136', 'HLA-A*24:137', 'HLA-A*24:138', 'HLA-A*24:139', 'HLA-A*24:140', 'HLA-A*24:141',
                                     'HLA-A*24:142', 'HLA-A*24:143',
                                     'HLA-A*24:144', 'HLA-A*25:01', 'HLA-A*25:02', 'HLA-A*25:03', 'HLA-A*25:04', 'HLA-A*25:05', 'HLA-A*25:06',
                                     'HLA-A*25:07', 'HLA-A*25:08', 'HLA-A*25:09',
                                     'HLA-A*25:10', 'HLA-A*25:11', 'HLA-A*25:13', 'HLA-A*26:01', 'HLA-A*26:02', 'HLA-A*26:03', 'HLA-A*26:04',
                                     'HLA-A*26:05', 'HLA-A*26:06', 'HLA-A*26:07',
                                     'HLA-A*26:08', 'HLA-A*26:09', 'HLA-A*26:10', 'HLA-A*26:12', 'HLA-A*26:13', 'HLA-A*26:14', 'HLA-A*26:15',
                                     'HLA-A*26:16', 'HLA-A*26:17', 'HLA-A*26:18',
                                     'HLA-A*26:19', 'HLA-A*26:20', 'HLA-A*26:21', 'HLA-A*26:22', 'HLA-A*26:23', 'HLA-A*26:24', 'HLA-A*26:26',
                                     'HLA-A*26:27', 'HLA-A*26:28', 'HLA-A*26:29',
                                     'HLA-A*26:30', 'HLA-A*26:31', 'HLA-A*26:32', 'HLA-A*26:33', 'HLA-A*26:34', 'HLA-A*26:35', 'HLA-A*26:36',
                                     'HLA-A*26:37', 'HLA-A*26:38', 'HLA-A*26:39',
                                     'HLA-A*26:40', 'HLA-A*26:41', 'HLA-A*26:42', 'HLA-A*26:43', 'HLA-A*26:45', 'HLA-A*26:46', 'HLA-A*26:47',
                                     'HLA-A*26:48', 'HLA-A*26:49', 'HLA-A*26:50',
                                     'HLA-A*29:01', 'HLA-A*29:02', 'HLA-A*29:03', 'HLA-A*29:04', 'HLA-A*29:05', 'HLA-A*29:06', 'HLA-A*29:07',
                                     'HLA-A*29:09', 'HLA-A*29:10', 'HLA-A*29:11',
                                     'HLA-A*29:12', 'HLA-A*29:13', 'HLA-A*29:14', 'HLA-A*29:15', 'HLA-A*29:16', 'HLA-A*29:17', 'HLA-A*29:18',
                                     'HLA-A*29:19', 'HLA-A*29:20', 'HLA-A*29:21',
                                     'HLA-A*29:22', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*30:03', 'HLA-A*30:04', 'HLA-A*30:06', 'HLA-A*30:07',
                                     'HLA-A*30:08', 'HLA-A*30:09', 'HLA-A*30:10',
                                     'HLA-A*30:11', 'HLA-A*30:12', 'HLA-A*30:13', 'HLA-A*30:15', 'HLA-A*30:16', 'HLA-A*30:17', 'HLA-A*30:18',
                                     'HLA-A*30:19', 'HLA-A*30:20', 'HLA-A*30:22',
                                     'HLA-A*30:23', 'HLA-A*30:24', 'HLA-A*30:25', 'HLA-A*30:26', 'HLA-A*30:28', 'HLA-A*30:29', 'HLA-A*30:30',
                                     'HLA-A*30:31', 'HLA-A*30:32', 'HLA-A*30:33',
                                     'HLA-A*30:34', 'HLA-A*30:35', 'HLA-A*30:36', 'HLA-A*30:37', 'HLA-A*30:38', 'HLA-A*30:39', 'HLA-A*30:40',
                                     'HLA-A*30:41', 'HLA-A*31:01', 'HLA-A*31:02',
                                     'HLA-A*31:03', 'HLA-A*31:04', 'HLA-A*31:05', 'HLA-A*31:06', 'HLA-A*31:07', 'HLA-A*31:08', 'HLA-A*31:09',
                                     'HLA-A*31:10', 'HLA-A*31:11', 'HLA-A*31:12',
                                     'HLA-A*31:13', 'HLA-A*31:15', 'HLA-A*31:16', 'HLA-A*31:17', 'HLA-A*31:18', 'HLA-A*31:19', 'HLA-A*31:20',
                                     'HLA-A*31:21', 'HLA-A*31:22', 'HLA-A*31:23',
                                     'HLA-A*31:24', 'HLA-A*31:25', 'HLA-A*31:26', 'HLA-A*31:27', 'HLA-A*31:28', 'HLA-A*31:29', 'HLA-A*31:30',
                                     'HLA-A*31:31', 'HLA-A*31:32', 'HLA-A*31:33',
                                     'HLA-A*31:34', 'HLA-A*31:35', 'HLA-A*31:36', 'HLA-A*31:37', 'HLA-A*32:01', 'HLA-A*32:02', 'HLA-A*32:03',
                                     'HLA-A*32:04', 'HLA-A*32:05', 'HLA-A*32:06',
                                     'HLA-A*32:07', 'HLA-A*32:08', 'HLA-A*32:09', 'HLA-A*32:10', 'HLA-A*32:12', 'HLA-A*32:13', 'HLA-A*32:14',
                                     'HLA-A*32:15', 'HLA-A*32:16', 'HLA-A*32:17',
                                     'HLA-A*32:18', 'HLA-A*32:20', 'HLA-A*32:21', 'HLA-A*32:22', 'HLA-A*32:23', 'HLA-A*32:24', 'HLA-A*32:25',
                                     'HLA-A*33:01', 'HLA-A*33:03', 'HLA-A*33:04',
                                     'HLA-A*33:05', 'HLA-A*33:06', 'HLA-A*33:07', 'HLA-A*33:08', 'HLA-A*33:09', 'HLA-A*33:10', 'HLA-A*33:11',
                                     'HLA-A*33:12', 'HLA-A*33:13', 'HLA-A*33:14',
                                     'HLA-A*33:15', 'HLA-A*33:16', 'HLA-A*33:17', 'HLA-A*33:18', 'HLA-A*33:19', 'HLA-A*33:20', 'HLA-A*33:21',
                                     'HLA-A*33:22', 'HLA-A*33:23', 'HLA-A*33:24',
                                     'HLA-A*33:25', 'HLA-A*33:26', 'HLA-A*33:27', 'HLA-A*33:28', 'HLA-A*33:29', 'HLA-A*33:30', 'HLA-A*33:31',
                                     'HLA-A*34:01', 'HLA-A*34:02', 'HLA-A*34:03',
                                     'HLA-A*34:04', 'HLA-A*34:05', 'HLA-A*34:06', 'HLA-A*34:07', 'HLA-A*34:08', 'HLA-A*36:01', 'HLA-A*36:02',
                                     'HLA-A*36:03', 'HLA-A*36:04', 'HLA-A*36:05',
                                     'HLA-A*43:01', 'HLA-A*66:01', 'HLA-A*66:02', 'HLA-A*66:03', 'HLA-A*66:04', 'HLA-A*66:05', 'HLA-A*66:06',
                                     'HLA-A*66:07', 'HLA-A*66:08', 'HLA-A*66:09',
                                     'HLA-A*66:10', 'HLA-A*66:11', 'HLA-A*66:12', 'HLA-A*66:13', 'HLA-A*66:14', 'HLA-A*66:15', 'HLA-A*68:01',
                                     'HLA-A*68:02', 'HLA-A*68:03', 'HLA-A*68:04',
                                     'HLA-A*68:05', 'HLA-A*68:06', 'HLA-A*68:07', 'HLA-A*68:08', 'HLA-A*68:09', 'HLA-A*68:10', 'HLA-A*68:12',
                                     'HLA-A*68:13', 'HLA-A*68:14', 'HLA-A*68:15',
                                     'HLA-A*68:16', 'HLA-A*68:17', 'HLA-A*68:19', 'HLA-A*68:20', 'HLA-A*68:21', 'HLA-A*68:22', 'HLA-A*68:23',
                                     'HLA-A*68:24', 'HLA-A*68:25', 'HLA-A*68:26',
                                     'HLA-A*68:27', 'HLA-A*68:28', 'HLA-A*68:29', 'HLA-A*68:30', 'HLA-A*68:31', 'HLA-A*68:32', 'HLA-A*68:33',
                                     'HLA-A*68:34', 'HLA-A*68:35', 'HLA-A*68:36',
                                     'HLA-A*68:37', 'HLA-A*68:38', 'HLA-A*68:39', 'HLA-A*68:40', 'HLA-A*68:41', 'HLA-A*68:42', 'HLA-A*68:43',
                                     'HLA-A*68:44', 'HLA-A*68:45', 'HLA-A*68:46',
                                     'HLA-A*68:47', 'HLA-A*68:48', 'HLA-A*68:50', 'HLA-A*68:51', 'HLA-A*68:52', 'HLA-A*68:53', 'HLA-A*68:54',
                                     'HLA-A*69:01', 'HLA-A*74:01', 'HLA-A*74:02',
                                     'HLA-A*74:03', 'HLA-A*74:04', 'HLA-A*74:05', 'HLA-A*74:06', 'HLA-A*74:07', 'HLA-A*74:08', 'HLA-A*74:09',
                                     'HLA-A*74:10', 'HLA-A*74:11', 'HLA-A*74:13',
                                     'HLA-A*80:01', 'HLA-A*80:02', 'HLA-B*07:02', 'HLA-B*07:03', 'HLA-B*07:04', 'HLA-B*07:05', 'HLA-B*07:06',
                                     'HLA-B*07:07', 'HLA-B*07:08', 'HLA-B*07:09',
                                     'HLA-B*07:10', 'HLA-B*07:11', 'HLA-B*07:12', 'HLA-B*07:13', 'HLA-B*07:14', 'HLA-B*07:15', 'HLA-B*07:16',
                                     'HLA-B*07:17', 'HLA-B*07:18', 'HLA-B*07:19',
                                     'HLA-B*07:20', 'HLA-B*07:21', 'HLA-B*07:22', 'HLA-B*07:23', 'HLA-B*07:24', 'HLA-B*07:25', 'HLA-B*07:26',
                                     'HLA-B*07:27', 'HLA-B*07:28', 'HLA-B*07:29',
                                     'HLA-B*07:30', 'HLA-B*07:31', 'HLA-B*07:32', 'HLA-B*07:33', 'HLA-B*07:34', 'HLA-B*07:35', 'HLA-B*07:36',
                                     'HLA-B*07:37', 'HLA-B*07:38', 'HLA-B*07:39',
                                     'HLA-B*07:40', 'HLA-B*07:41', 'HLA-B*07:42', 'HLA-B*07:43', 'HLA-B*07:44', 'HLA-B*07:45', 'HLA-B*07:46',
                                     'HLA-B*07:47', 'HLA-B*07:48', 'HLA-B*07:50',
                                     'HLA-B*07:51', 'HLA-B*07:52', 'HLA-B*07:53', 'HLA-B*07:54', 'HLA-B*07:55', 'HLA-B*07:56', 'HLA-B*07:57',
                                     'HLA-B*07:58', 'HLA-B*07:59', 'HLA-B*07:60',
                                     'HLA-B*07:61', 'HLA-B*07:62', 'HLA-B*07:63', 'HLA-B*07:64', 'HLA-B*07:65', 'HLA-B*07:66', 'HLA-B*07:68',
                                     'HLA-B*07:69', 'HLA-B*07:70', 'HLA-B*07:71',
                                     'HLA-B*07:72', 'HLA-B*07:73', 'HLA-B*07:74', 'HLA-B*07:75', 'HLA-B*07:76', 'HLA-B*07:77', 'HLA-B*07:78',
                                     'HLA-B*07:79', 'HLA-B*07:80', 'HLA-B*07:81',
                                     'HLA-B*07:82', 'HLA-B*07:83', 'HLA-B*07:84', 'HLA-B*07:85', 'HLA-B*07:86', 'HLA-B*07:87', 'HLA-B*07:88',
                                     'HLA-B*07:89', 'HLA-B*07:90', 'HLA-B*07:91',
                                     'HLA-B*07:92', 'HLA-B*07:93', 'HLA-B*07:94', 'HLA-B*07:95', 'HLA-B*07:96', 'HLA-B*07:97', 'HLA-B*07:98',
                                     'HLA-B*07:99', 'HLA-B*07:100', 'HLA-B*07:101',
                                     'HLA-B*07:102', 'HLA-B*07:103', 'HLA-B*07:104', 'HLA-B*07:105', 'HLA-B*07:106', 'HLA-B*07:107', 'HLA-B*07:108',
                                     'HLA-B*07:109', 'HLA-B*07:110',
                                     'HLA-B*07:112', 'HLA-B*07:113', 'HLA-B*07:114', 'HLA-B*07:115', 'HLA-B*08:01', 'HLA-B*08:02', 'HLA-B*08:03',
                                     'HLA-B*08:04', 'HLA-B*08:05',
                                     'HLA-B*08:07', 'HLA-B*08:09', 'HLA-B*08:10', 'HLA-B*08:11', 'HLA-B*08:12', 'HLA-B*08:13', 'HLA-B*08:14',
                                     'HLA-B*08:15', 'HLA-B*08:16', 'HLA-B*08:17',
                                     'HLA-B*08:18', 'HLA-B*08:20', 'HLA-B*08:21', 'HLA-B*08:22', 'HLA-B*08:23', 'HLA-B*08:24', 'HLA-B*08:25',
                                     'HLA-B*08:26', 'HLA-B*08:27', 'HLA-B*08:28',
                                     'HLA-B*08:29', 'HLA-B*08:31', 'HLA-B*08:32', 'HLA-B*08:33', 'HLA-B*08:34', 'HLA-B*08:35', 'HLA-B*08:36',
                                     'HLA-B*08:37', 'HLA-B*08:38', 'HLA-B*08:39',
                                     'HLA-B*08:40', 'HLA-B*08:41', 'HLA-B*08:42', 'HLA-B*08:43', 'HLA-B*08:44', 'HLA-B*08:45', 'HLA-B*08:46',
                                     'HLA-B*08:47', 'HLA-B*08:48', 'HLA-B*08:49',
                                     'HLA-B*08:50', 'HLA-B*08:51', 'HLA-B*08:52', 'HLA-B*08:53', 'HLA-B*08:54', 'HLA-B*08:55', 'HLA-B*08:56',
                                     'HLA-B*08:57', 'HLA-B*08:58', 'HLA-B*08:59',
                                     'HLA-B*08:60', 'HLA-B*08:61', 'HLA-B*08:62', 'HLA-B*13:01', 'HLA-B*13:02', 'HLA-B*13:03', 'HLA-B*13:04',
                                     'HLA-B*13:06', 'HLA-B*13:09', 'HLA-B*13:10',
                                     'HLA-B*13:11', 'HLA-B*13:12', 'HLA-B*13:13', 'HLA-B*13:14', 'HLA-B*13:15', 'HLA-B*13:16', 'HLA-B*13:17',
                                     'HLA-B*13:18', 'HLA-B*13:19', 'HLA-B*13:20',
                                     'HLA-B*13:21', 'HLA-B*13:22', 'HLA-B*13:23', 'HLA-B*13:25', 'HLA-B*13:26', 'HLA-B*13:27', 'HLA-B*13:28',
                                     'HLA-B*13:29', 'HLA-B*13:30', 'HLA-B*13:31',
                                     'HLA-B*13:32', 'HLA-B*13:33', 'HLA-B*13:34', 'HLA-B*13:35', 'HLA-B*13:36', 'HLA-B*13:37', 'HLA-B*13:38',
                                     'HLA-B*13:39', 'HLA-B*14:01', 'HLA-B*14:02',
                                     'HLA-B*14:03', 'HLA-B*14:04', 'HLA-B*14:05', 'HLA-B*14:06', 'HLA-B*14:08', 'HLA-B*14:09', 'HLA-B*14:10',
                                     'HLA-B*14:11', 'HLA-B*14:12', 'HLA-B*14:13',
                                     'HLA-B*14:14', 'HLA-B*14:15', 'HLA-B*14:16', 'HLA-B*14:17', 'HLA-B*14:18', 'HLA-B*15:01', 'HLA-B*15:02',
                                     'HLA-B*15:03', 'HLA-B*15:04', 'HLA-B*15:05',
                                     'HLA-B*15:06', 'HLA-B*15:07', 'HLA-B*15:08', 'HLA-B*15:09', 'HLA-B*15:10', 'HLA-B*15:11', 'HLA-B*15:12',
                                     'HLA-B*15:13', 'HLA-B*15:14', 'HLA-B*15:15',
                                     'HLA-B*15:16', 'HLA-B*15:17', 'HLA-B*15:18', 'HLA-B*15:19', 'HLA-B*15:20', 'HLA-B*15:21', 'HLA-B*15:23',
                                     'HLA-B*15:24', 'HLA-B*15:25', 'HLA-B*15:27',
                                     'HLA-B*15:28', 'HLA-B*15:29', 'HLA-B*15:30', 'HLA-B*15:31', 'HLA-B*15:32', 'HLA-B*15:33', 'HLA-B*15:34',
                                     'HLA-B*15:35', 'HLA-B*15:36', 'HLA-B*15:37',
                                     'HLA-B*15:38', 'HLA-B*15:39', 'HLA-B*15:40', 'HLA-B*15:42', 'HLA-B*15:43', 'HLA-B*15:44', 'HLA-B*15:45',
                                     'HLA-B*15:46', 'HLA-B*15:47', 'HLA-B*15:48',
                                     'HLA-B*15:49', 'HLA-B*15:50', 'HLA-B*15:51', 'HLA-B*15:52', 'HLA-B*15:53', 'HLA-B*15:54', 'HLA-B*15:55',
                                     'HLA-B*15:56', 'HLA-B*15:57', 'HLA-B*15:58',
                                     'HLA-B*15:60', 'HLA-B*15:61', 'HLA-B*15:62', 'HLA-B*15:63', 'HLA-B*15:64', 'HLA-B*15:65', 'HLA-B*15:66',
                                     'HLA-B*15:67', 'HLA-B*15:68', 'HLA-B*15:69',
                                     'HLA-B*15:70', 'HLA-B*15:71', 'HLA-B*15:72', 'HLA-B*15:73', 'HLA-B*15:74', 'HLA-B*15:75', 'HLA-B*15:76',
                                     'HLA-B*15:77', 'HLA-B*15:78', 'HLA-B*15:80',
                                     'HLA-B*15:81', 'HLA-B*15:82', 'HLA-B*15:83', 'HLA-B*15:84', 'HLA-B*15:85', 'HLA-B*15:86', 'HLA-B*15:87',
                                     'HLA-B*15:88', 'HLA-B*15:89', 'HLA-B*15:90',
                                     'HLA-B*15:91', 'HLA-B*15:92', 'HLA-B*15:93', 'HLA-B*15:95', 'HLA-B*15:96', 'HLA-B*15:97', 'HLA-B*15:98',
                                     'HLA-B*15:99', 'HLA-B*15:101', 'HLA-B*15:102',
                                     'HLA-B*15:103', 'HLA-B*15:104', 'HLA-B*15:105', 'HLA-B*15:106', 'HLA-B*15:107', 'HLA-B*15:108', 'HLA-B*15:109',
                                     'HLA-B*15:110', 'HLA-B*15:112',
                                     'HLA-B*15:113', 'HLA-B*15:114', 'HLA-B*15:115', 'HLA-B*15:116', 'HLA-B*15:117', 'HLA-B*15:118', 'HLA-B*15:119',
                                     'HLA-B*15:120', 'HLA-B*15:121',
                                     'HLA-B*15:122', 'HLA-B*15:123', 'HLA-B*15:124', 'HLA-B*15:125', 'HLA-B*15:126', 'HLA-B*15:127', 'HLA-B*15:128',
                                     'HLA-B*15:129', 'HLA-B*15:131',
                                     'HLA-B*15:132', 'HLA-B*15:133', 'HLA-B*15:134', 'HLA-B*15:135', 'HLA-B*15:136', 'HLA-B*15:137', 'HLA-B*15:138',
                                     'HLA-B*15:139', 'HLA-B*15:140',
                                     'HLA-B*15:141', 'HLA-B*15:142', 'HLA-B*15:143', 'HLA-B*15:144', 'HLA-B*15:145', 'HLA-B*15:146', 'HLA-B*15:147',
                                     'HLA-B*15:148', 'HLA-B*15:150',
                                     'HLA-B*15:151', 'HLA-B*15:152', 'HLA-B*15:153', 'HLA-B*15:154', 'HLA-B*15:155', 'HLA-B*15:156', 'HLA-B*15:157',
                                     'HLA-B*15:158', 'HLA-B*15:159',
                                     'HLA-B*15:160', 'HLA-B*15:161', 'HLA-B*15:162', 'HLA-B*15:163', 'HLA-B*15:164', 'HLA-B*15:165', 'HLA-B*15:166',
                                     'HLA-B*15:167', 'HLA-B*15:168',
                                     'HLA-B*15:169', 'HLA-B*15:170', 'HLA-B*15:171', 'HLA-B*15:172', 'HLA-B*15:173', 'HLA-B*15:174', 'HLA-B*15:175',
                                     'HLA-B*15:176', 'HLA-B*15:177',
                                     'HLA-B*15:178', 'HLA-B*15:179', 'HLA-B*15:180', 'HLA-B*15:183', 'HLA-B*15:184', 'HLA-B*15:185', 'HLA-B*15:186',
                                     'HLA-B*15:187', 'HLA-B*15:188',
                                     'HLA-B*15:189', 'HLA-B*15:191', 'HLA-B*15:192', 'HLA-B*15:193', 'HLA-B*15:194', 'HLA-B*15:195', 'HLA-B*15:196',
                                     'HLA-B*15:197', 'HLA-B*15:198',
                                     'HLA-B*15:199', 'HLA-B*15:200', 'HLA-B*15:201', 'HLA-B*15:202', 'HLA-B*18:01', 'HLA-B*18:02', 'HLA-B*18:03',
                                     'HLA-B*18:04', 'HLA-B*18:05',
                                     'HLA-B*18:06', 'HLA-B*18:07', 'HLA-B*18:08', 'HLA-B*18:09', 'HLA-B*18:10', 'HLA-B*18:11', 'HLA-B*18:12',
                                     'HLA-B*18:13', 'HLA-B*18:14', 'HLA-B*18:15',
                                     'HLA-B*18:18', 'HLA-B*18:19', 'HLA-B*18:20', 'HLA-B*18:21', 'HLA-B*18:22', 'HLA-B*18:24', 'HLA-B*18:25',
                                     'HLA-B*18:26', 'HLA-B*18:27', 'HLA-B*18:28',
                                     'HLA-B*18:29', 'HLA-B*18:30', 'HLA-B*18:31', 'HLA-B*18:32', 'HLA-B*18:33', 'HLA-B*18:34', 'HLA-B*18:35',
                                     'HLA-B*18:36', 'HLA-B*18:37', 'HLA-B*18:38',
                                     'HLA-B*18:39', 'HLA-B*18:40', 'HLA-B*18:41', 'HLA-B*18:42', 'HLA-B*18:43', 'HLA-B*18:44', 'HLA-B*18:45',
                                     'HLA-B*18:46', 'HLA-B*18:47', 'HLA-B*18:48',
                                     'HLA-B*18:49', 'HLA-B*18:50', 'HLA-B*27:01', 'HLA-B*27:02', 'HLA-B*27:03', 'HLA-B*27:04', 'HLA-B*27:05',
                                     'HLA-B*27:06', 'HLA-B*27:07', 'HLA-B*27:08',
                                     'HLA-B*27:09', 'HLA-B*27:10', 'HLA-B*27:11', 'HLA-B*27:12', 'HLA-B*27:13', 'HLA-B*27:14', 'HLA-B*27:15',
                                     'HLA-B*27:16', 'HLA-B*27:17', 'HLA-B*27:18',
                                     'HLA-B*27:19', 'HLA-B*27:20', 'HLA-B*27:21', 'HLA-B*27:23', 'HLA-B*27:24', 'HLA-B*27:25', 'HLA-B*27:26',
                                     'HLA-B*27:27', 'HLA-B*27:28', 'HLA-B*27:29',
                                     'HLA-B*27:30', 'HLA-B*27:31', 'HLA-B*27:32', 'HLA-B*27:33', 'HLA-B*27:34', 'HLA-B*27:35', 'HLA-B*27:36',
                                     'HLA-B*27:37', 'HLA-B*27:38', 'HLA-B*27:39',
                                     'HLA-B*27:40', 'HLA-B*27:41', 'HLA-B*27:42', 'HLA-B*27:43', 'HLA-B*27:44', 'HLA-B*27:45', 'HLA-B*27:46',
                                     'HLA-B*27:47', 'HLA-B*27:48', 'HLA-B*27:49',
                                     'HLA-B*27:50', 'HLA-B*27:51', 'HLA-B*27:52', 'HLA-B*27:53', 'HLA-B*27:54', 'HLA-B*27:55', 'HLA-B*27:56',
                                     'HLA-B*27:57', 'HLA-B*27:58', 'HLA-B*27:60',
                                     'HLA-B*27:61', 'HLA-B*27:62', 'HLA-B*27:63', 'HLA-B*27:67', 'HLA-B*27:68', 'HLA-B*27:69', 'HLA-B*35:01',
                                     'HLA-B*35:02', 'HLA-B*35:03', 'HLA-B*35:04',
                                     'HLA-B*35:05', 'HLA-B*35:06', 'HLA-B*35:07', 'HLA-B*35:08', 'HLA-B*35:09', 'HLA-B*35:10', 'HLA-B*35:11',
                                     'HLA-B*35:12', 'HLA-B*35:13', 'HLA-B*35:14',
                                     'HLA-B*35:15', 'HLA-B*35:16', 'HLA-B*35:17', 'HLA-B*35:18', 'HLA-B*35:19', 'HLA-B*35:20', 'HLA-B*35:21',
                                     'HLA-B*35:22', 'HLA-B*35:23', 'HLA-B*35:24',
                                     'HLA-B*35:25', 'HLA-B*35:26', 'HLA-B*35:27', 'HLA-B*35:28', 'HLA-B*35:29', 'HLA-B*35:30', 'HLA-B*35:31',
                                     'HLA-B*35:32', 'HLA-B*35:33', 'HLA-B*35:34',
                                     'HLA-B*35:35', 'HLA-B*35:36', 'HLA-B*35:37', 'HLA-B*35:38', 'HLA-B*35:39', 'HLA-B*35:41', 'HLA-B*35:42',
                                     'HLA-B*35:43', 'HLA-B*35:44', 'HLA-B*35:45',
                                     'HLA-B*35:46', 'HLA-B*35:47', 'HLA-B*35:48', 'HLA-B*35:49', 'HLA-B*35:50', 'HLA-B*35:51', 'HLA-B*35:52',
                                     'HLA-B*35:54', 'HLA-B*35:55', 'HLA-B*35:56',
                                     'HLA-B*35:57', 'HLA-B*35:58', 'HLA-B*35:59', 'HLA-B*35:60', 'HLA-B*35:61', 'HLA-B*35:62', 'HLA-B*35:63',
                                     'HLA-B*35:64', 'HLA-B*35:66', 'HLA-B*35:67',
                                     'HLA-B*35:68', 'HLA-B*35:69', 'HLA-B*35:70', 'HLA-B*35:71', 'HLA-B*35:72', 'HLA-B*35:74', 'HLA-B*35:75',
                                     'HLA-B*35:76', 'HLA-B*35:77', 'HLA-B*35:78',
                                     'HLA-B*35:79', 'HLA-B*35:80', 'HLA-B*35:81', 'HLA-B*35:82', 'HLA-B*35:83', 'HLA-B*35:84', 'HLA-B*35:85',
                                     'HLA-B*35:86', 'HLA-B*35:87', 'HLA-B*35:88',
                                     'HLA-B*35:89', 'HLA-B*35:90', 'HLA-B*35:91', 'HLA-B*35:92', 'HLA-B*35:93', 'HLA-B*35:94', 'HLA-B*35:95',
                                     'HLA-B*35:96', 'HLA-B*35:97', 'HLA-B*35:98',
                                     'HLA-B*35:99', 'HLA-B*35:100', 'HLA-B*35:101', 'HLA-B*35:102', 'HLA-B*35:103', 'HLA-B*35:104', 'HLA-B*35:105',
                                     'HLA-B*35:106', 'HLA-B*35:107',
                                     'HLA-B*35:108', 'HLA-B*35:109', 'HLA-B*35:110', 'HLA-B*35:111', 'HLA-B*35:112', 'HLA-B*35:113', 'HLA-B*35:114',
                                     'HLA-B*35:115', 'HLA-B*35:116',
                                     'HLA-B*35:117', 'HLA-B*35:118', 'HLA-B*35:119', 'HLA-B*35:120', 'HLA-B*35:121', 'HLA-B*35:122', 'HLA-B*35:123',
                                     'HLA-B*35:124', 'HLA-B*35:125',
                                     'HLA-B*35:126', 'HLA-B*35:127', 'HLA-B*35:128', 'HLA-B*35:131', 'HLA-B*35:132', 'HLA-B*35:133', 'HLA-B*35:135',
                                     'HLA-B*35:136', 'HLA-B*35:137',
                                     'HLA-B*35:138', 'HLA-B*35:139', 'HLA-B*35:140', 'HLA-B*35:141', 'HLA-B*35:142', 'HLA-B*35:143', 'HLA-B*35:144',
                                     'HLA-B*37:01', 'HLA-B*37:02',
                                     'HLA-B*37:04', 'HLA-B*37:05', 'HLA-B*37:06', 'HLA-B*37:07', 'HLA-B*37:08', 'HLA-B*37:09', 'HLA-B*37:10',
                                     'HLA-B*37:11', 'HLA-B*37:12', 'HLA-B*37:13',
                                     'HLA-B*37:14', 'HLA-B*37:15', 'HLA-B*37:17', 'HLA-B*37:18', 'HLA-B*37:19', 'HLA-B*37:20', 'HLA-B*37:21',
                                     'HLA-B*37:22', 'HLA-B*37:23', 'HLA-B*38:01',
                                     'HLA-B*38:02', 'HLA-B*38:03', 'HLA-B*38:04', 'HLA-B*38:05', 'HLA-B*38:06', 'HLA-B*38:07', 'HLA-B*38:08',
                                     'HLA-B*38:09', 'HLA-B*38:10', 'HLA-B*38:11',
                                     'HLA-B*38:12', 'HLA-B*38:13', 'HLA-B*38:14', 'HLA-B*38:15', 'HLA-B*38:16', 'HLA-B*38:17', 'HLA-B*38:18',
                                     'HLA-B*38:19', 'HLA-B*38:20', 'HLA-B*38:21',
                                     'HLA-B*38:22', 'HLA-B*38:23', 'HLA-B*39:01', 'HLA-B*39:02', 'HLA-B*39:03', 'HLA-B*39:04', 'HLA-B*39:05',
                                     'HLA-B*39:06', 'HLA-B*39:07', 'HLA-B*39:08',
                                     'HLA-B*39:09', 'HLA-B*39:10', 'HLA-B*39:11', 'HLA-B*39:12', 'HLA-B*39:13', 'HLA-B*39:14', 'HLA-B*39:15',
                                     'HLA-B*39:16', 'HLA-B*39:17', 'HLA-B*39:18',
                                     'HLA-B*39:19', 'HLA-B*39:20', 'HLA-B*39:22', 'HLA-B*39:23', 'HLA-B*39:24', 'HLA-B*39:26', 'HLA-B*39:27',
                                     'HLA-B*39:28', 'HLA-B*39:29', 'HLA-B*39:30',
                                     'HLA-B*39:31', 'HLA-B*39:32', 'HLA-B*39:33', 'HLA-B*39:34', 'HLA-B*39:35', 'HLA-B*39:36', 'HLA-B*39:37',
                                     'HLA-B*39:39', 'HLA-B*39:41', 'HLA-B*39:42',
                                     'HLA-B*39:43', 'HLA-B*39:44', 'HLA-B*39:45', 'HLA-B*39:46', 'HLA-B*39:47', 'HLA-B*39:48', 'HLA-B*39:49',
                                     'HLA-B*39:50', 'HLA-B*39:51', 'HLA-B*39:52',
                                     'HLA-B*39:53', 'HLA-B*39:54', 'HLA-B*39:55', 'HLA-B*39:56', 'HLA-B*39:57', 'HLA-B*39:58', 'HLA-B*39:59',
                                     'HLA-B*39:60', 'HLA-B*40:01', 'HLA-B*40:02',
                                     'HLA-B*40:03', 'HLA-B*40:04', 'HLA-B*40:05', 'HLA-B*40:06', 'HLA-B*40:07', 'HLA-B*40:08', 'HLA-B*40:09',
                                     'HLA-B*40:10', 'HLA-B*40:11', 'HLA-B*40:12',
                                     'HLA-B*40:13', 'HLA-B*40:14', 'HLA-B*40:15', 'HLA-B*40:16', 'HLA-B*40:18', 'HLA-B*40:19', 'HLA-B*40:20',
                                     'HLA-B*40:21', 'HLA-B*40:23', 'HLA-B*40:24',
                                     'HLA-B*40:25', 'HLA-B*40:26', 'HLA-B*40:27', 'HLA-B*40:28', 'HLA-B*40:29', 'HLA-B*40:30', 'HLA-B*40:31',
                                     'HLA-B*40:32', 'HLA-B*40:33', 'HLA-B*40:34',
                                     'HLA-B*40:35', 'HLA-B*40:36', 'HLA-B*40:37', 'HLA-B*40:38', 'HLA-B*40:39', 'HLA-B*40:40', 'HLA-B*40:42',
                                     'HLA-B*40:43', 'HLA-B*40:44', 'HLA-B*40:45',
                                     'HLA-B*40:46', 'HLA-B*40:47', 'HLA-B*40:48', 'HLA-B*40:49', 'HLA-B*40:50', 'HLA-B*40:51', 'HLA-B*40:52',
                                     'HLA-B*40:53', 'HLA-B*40:54', 'HLA-B*40:55',
                                     'HLA-B*40:56', 'HLA-B*40:57', 'HLA-B*40:58', 'HLA-B*40:59', 'HLA-B*40:60', 'HLA-B*40:61', 'HLA-B*40:62',
                                     'HLA-B*40:63', 'HLA-B*40:64', 'HLA-B*40:65',
                                     'HLA-B*40:66', 'HLA-B*40:67', 'HLA-B*40:68', 'HLA-B*40:69', 'HLA-B*40:70', 'HLA-B*40:71', 'HLA-B*40:72',
                                     'HLA-B*40:73', 'HLA-B*40:74', 'HLA-B*40:75',
                                     'HLA-B*40:76', 'HLA-B*40:77', 'HLA-B*40:78', 'HLA-B*40:79', 'HLA-B*40:80', 'HLA-B*40:81', 'HLA-B*40:82',
                                     'HLA-B*40:83', 'HLA-B*40:84', 'HLA-B*40:85',
                                     'HLA-B*40:86', 'HLA-B*40:87', 'HLA-B*40:88', 'HLA-B*40:89', 'HLA-B*40:90', 'HLA-B*40:91', 'HLA-B*40:92',
                                     'HLA-B*40:93', 'HLA-B*40:94', 'HLA-B*40:95',
                                     'HLA-B*40:96', 'HLA-B*40:97', 'HLA-B*40:98', 'HLA-B*40:99', 'HLA-B*40:100', 'HLA-B*40:101', 'HLA-B*40:102',
                                     'HLA-B*40:103', 'HLA-B*40:104',
                                     'HLA-B*40:105', 'HLA-B*40:106', 'HLA-B*40:107', 'HLA-B*40:108', 'HLA-B*40:109', 'HLA-B*40:110', 'HLA-B*40:111',
                                     'HLA-B*40:112', 'HLA-B*40:113',
                                     'HLA-B*40:114', 'HLA-B*40:115', 'HLA-B*40:116', 'HLA-B*40:117', 'HLA-B*40:119', 'HLA-B*40:120', 'HLA-B*40:121',
                                     'HLA-B*40:122', 'HLA-B*40:123',
                                     'HLA-B*40:124', 'HLA-B*40:125', 'HLA-B*40:126', 'HLA-B*40:127', 'HLA-B*40:128', 'HLA-B*40:129', 'HLA-B*40:130',
                                     'HLA-B*40:131', 'HLA-B*40:132',
                                     'HLA-B*40:134', 'HLA-B*40:135', 'HLA-B*40:136', 'HLA-B*40:137', 'HLA-B*40:138', 'HLA-B*40:139', 'HLA-B*40:140',
                                     'HLA-B*40:141', 'HLA-B*40:143',
                                     'HLA-B*40:145', 'HLA-B*40:146', 'HLA-B*40:147', 'HLA-B*41:01', 'HLA-B*41:02', 'HLA-B*41:03', 'HLA-B*41:04',
                                     'HLA-B*41:05', 'HLA-B*41:06', 'HLA-B*41:07',
                                     'HLA-B*41:08', 'HLA-B*41:09', 'HLA-B*41:10', 'HLA-B*41:11', 'HLA-B*41:12', 'HLA-B*42:01', 'HLA-B*42:02',
                                     'HLA-B*42:04', 'HLA-B*42:05', 'HLA-B*42:06',
                                     'HLA-B*42:07', 'HLA-B*42:08', 'HLA-B*42:09', 'HLA-B*42:10', 'HLA-B*42:11', 'HLA-B*42:12', 'HLA-B*42:13',
                                     'HLA-B*42:14', 'HLA-B*44:02', 'HLA-B*44:03',
                                     'HLA-B*44:04', 'HLA-B*44:05', 'HLA-B*44:06', 'HLA-B*44:07', 'HLA-B*44:08', 'HLA-B*44:09', 'HLA-B*44:10',
                                     'HLA-B*44:11', 'HLA-B*44:12', 'HLA-B*44:13',
                                     'HLA-B*44:14', 'HLA-B*44:15', 'HLA-B*44:16', 'HLA-B*44:17', 'HLA-B*44:18', 'HLA-B*44:20', 'HLA-B*44:21',
                                     'HLA-B*44:22', 'HLA-B*44:24', 'HLA-B*44:25',
                                     'HLA-B*44:26', 'HLA-B*44:27', 'HLA-B*44:28', 'HLA-B*44:29', 'HLA-B*44:30', 'HLA-B*44:31', 'HLA-B*44:32',
                                     'HLA-B*44:33', 'HLA-B*44:34', 'HLA-B*44:35',
                                     'HLA-B*44:36', 'HLA-B*44:37', 'HLA-B*44:38', 'HLA-B*44:39', 'HLA-B*44:40', 'HLA-B*44:41', 'HLA-B*44:42',
                                     'HLA-B*44:43', 'HLA-B*44:44', 'HLA-B*44:45',
                                     'HLA-B*44:46', 'HLA-B*44:47', 'HLA-B*44:48', 'HLA-B*44:49', 'HLA-B*44:50', 'HLA-B*44:51', 'HLA-B*44:53',
                                     'HLA-B*44:54', 'HLA-B*44:55', 'HLA-B*44:57',
                                     'HLA-B*44:59', 'HLA-B*44:60', 'HLA-B*44:62', 'HLA-B*44:63', 'HLA-B*44:64', 'HLA-B*44:65', 'HLA-B*44:66',
                                     'HLA-B*44:67', 'HLA-B*44:68', 'HLA-B*44:69',
                                     'HLA-B*44:70', 'HLA-B*44:71', 'HLA-B*44:72', 'HLA-B*44:73', 'HLA-B*44:74', 'HLA-B*44:75', 'HLA-B*44:76',
                                     'HLA-B*44:77', 'HLA-B*44:78', 'HLA-B*44:79',
                                     'HLA-B*44:80', 'HLA-B*44:81', 'HLA-B*44:82', 'HLA-B*44:83', 'HLA-B*44:84', 'HLA-B*44:85', 'HLA-B*44:86',
                                     'HLA-B*44:87', 'HLA-B*44:88', 'HLA-B*44:89',
                                     'HLA-B*44:90', 'HLA-B*44:91', 'HLA-B*44:92', 'HLA-B*44:93', 'HLA-B*44:94', 'HLA-B*44:95', 'HLA-B*44:96',
                                     'HLA-B*44:97', 'HLA-B*44:98', 'HLA-B*44:99',
                                     'HLA-B*44:100', 'HLA-B*44:101', 'HLA-B*44:102', 'HLA-B*44:103', 'HLA-B*44:104', 'HLA-B*44:105', 'HLA-B*44:106',
                                     'HLA-B*44:107', 'HLA-B*44:109',
                                     'HLA-B*44:110', 'HLA-B*45:01', 'HLA-B*45:02', 'HLA-B*45:03', 'HLA-B*45:04', 'HLA-B*45:05', 'HLA-B*45:06',
                                     'HLA-B*45:07', 'HLA-B*45:08', 'HLA-B*45:09',
                                     'HLA-B*45:10', 'HLA-B*45:11', 'HLA-B*45:12', 'HLA-B*46:01', 'HLA-B*46:02', 'HLA-B*46:03', 'HLA-B*46:04',
                                     'HLA-B*46:05', 'HLA-B*46:06', 'HLA-B*46:08',
                                     'HLA-B*46:09', 'HLA-B*46:10', 'HLA-B*46:11', 'HLA-B*46:12', 'HLA-B*46:13', 'HLA-B*46:14', 'HLA-B*46:16',
                                     'HLA-B*46:17', 'HLA-B*46:18', 'HLA-B*46:19',
                                     'HLA-B*46:20', 'HLA-B*46:21', 'HLA-B*46:22', 'HLA-B*46:23', 'HLA-B*46:24', 'HLA-B*47:01', 'HLA-B*47:02',
                                     'HLA-B*47:03', 'HLA-B*47:04', 'HLA-B*47:05',
                                     'HLA-B*47:06', 'HLA-B*47:07', 'HLA-B*48:01', 'HLA-B*48:02', 'HLA-B*48:03', 'HLA-B*48:04', 'HLA-B*48:05',
                                     'HLA-B*48:06', 'HLA-B*48:07', 'HLA-B*48:08',
                                     'HLA-B*48:09', 'HLA-B*48:10', 'HLA-B*48:11', 'HLA-B*48:12', 'HLA-B*48:13', 'HLA-B*48:14', 'HLA-B*48:15',
                                     'HLA-B*48:16', 'HLA-B*48:17', 'HLA-B*48:18',
                                     'HLA-B*48:19', 'HLA-B*48:20', 'HLA-B*48:21', 'HLA-B*48:22', 'HLA-B*48:23', 'HLA-B*49:01', 'HLA-B*49:02',
                                     'HLA-B*49:03', 'HLA-B*49:04', 'HLA-B*49:05',
                                     'HLA-B*49:06', 'HLA-B*49:07', 'HLA-B*49:08', 'HLA-B*49:09', 'HLA-B*49:10', 'HLA-B*50:01', 'HLA-B*50:02',
                                     'HLA-B*50:04', 'HLA-B*50:05', 'HLA-B*50:06',
                                     'HLA-B*50:07', 'HLA-B*50:08', 'HLA-B*50:09', 'HLA-B*51:01', 'HLA-B*51:02', 'HLA-B*51:03', 'HLA-B*51:04',
                                     'HLA-B*51:05', 'HLA-B*51:06', 'HLA-B*51:07',
                                     'HLA-B*51:08', 'HLA-B*51:09', 'HLA-B*51:12', 'HLA-B*51:13', 'HLA-B*51:14', 'HLA-B*51:15', 'HLA-B*51:16',
                                     'HLA-B*51:17', 'HLA-B*51:18', 'HLA-B*51:19',
                                     'HLA-B*51:20', 'HLA-B*51:21', 'HLA-B*51:22', 'HLA-B*51:23', 'HLA-B*51:24', 'HLA-B*51:26', 'HLA-B*51:28',
                                     'HLA-B*51:29', 'HLA-B*51:30', 'HLA-B*51:31',
                                     'HLA-B*51:32', 'HLA-B*51:33', 'HLA-B*51:34', 'HLA-B*51:35', 'HLA-B*51:36', 'HLA-B*51:37', 'HLA-B*51:38',
                                     'HLA-B*51:39', 'HLA-B*51:40', 'HLA-B*51:42',
                                     'HLA-B*51:43', 'HLA-B*51:45', 'HLA-B*51:46', 'HLA-B*51:48', 'HLA-B*51:49', 'HLA-B*51:50', 'HLA-B*51:51',
                                     'HLA-B*51:52', 'HLA-B*51:53', 'HLA-B*51:54',
                                     'HLA-B*51:55', 'HLA-B*51:56', 'HLA-B*51:57', 'HLA-B*51:58', 'HLA-B*51:59', 'HLA-B*51:60', 'HLA-B*51:61',
                                     'HLA-B*51:62', 'HLA-B*51:63', 'HLA-B*51:64',
                                     'HLA-B*51:65', 'HLA-B*51:66', 'HLA-B*51:67', 'HLA-B*51:68', 'HLA-B*51:69', 'HLA-B*51:70', 'HLA-B*51:71',
                                     'HLA-B*51:72', 'HLA-B*51:73', 'HLA-B*51:74',
                                     'HLA-B*51:75', 'HLA-B*51:76', 'HLA-B*51:77', 'HLA-B*51:78', 'HLA-B*51:79', 'HLA-B*51:80', 'HLA-B*51:81',
                                     'HLA-B*51:82', 'HLA-B*51:83', 'HLA-B*51:84',
                                     'HLA-B*51:85', 'HLA-B*51:86', 'HLA-B*51:87', 'HLA-B*51:88', 'HLA-B*51:89', 'HLA-B*51:90', 'HLA-B*51:91',
                                     'HLA-B*51:92', 'HLA-B*51:93', 'HLA-B*51:94',
                                     'HLA-B*51:95', 'HLA-B*51:96', 'HLA-B*52:01', 'HLA-B*52:02', 'HLA-B*52:03', 'HLA-B*52:04', 'HLA-B*52:05',
                                     'HLA-B*52:06', 'HLA-B*52:07', 'HLA-B*52:08',
                                     'HLA-B*52:09', 'HLA-B*52:10', 'HLA-B*52:11', 'HLA-B*52:12', 'HLA-B*52:13', 'HLA-B*52:14', 'HLA-B*52:15',
                                     'HLA-B*52:16', 'HLA-B*52:17', 'HLA-B*52:18',
                                     'HLA-B*52:19', 'HLA-B*52:20', 'HLA-B*52:21', 'HLA-B*53:01', 'HLA-B*53:02', 'HLA-B*53:03', 'HLA-B*53:04',
                                     'HLA-B*53:05', 'HLA-B*53:06', 'HLA-B*53:07',
                                     'HLA-B*53:08', 'HLA-B*53:09', 'HLA-B*53:10', 'HLA-B*53:11', 'HLA-B*53:12', 'HLA-B*53:13', 'HLA-B*53:14',
                                     'HLA-B*53:15', 'HLA-B*53:16', 'HLA-B*53:17',
                                     'HLA-B*53:18', 'HLA-B*53:19', 'HLA-B*53:20', 'HLA-B*53:21', 'HLA-B*53:22', 'HLA-B*53:23', 'HLA-B*54:01',
                                     'HLA-B*54:02', 'HLA-B*54:03', 'HLA-B*54:04',
                                     'HLA-B*54:06', 'HLA-B*54:07', 'HLA-B*54:09', 'HLA-B*54:10', 'HLA-B*54:11', 'HLA-B*54:12', 'HLA-B*54:13',
                                     'HLA-B*54:14', 'HLA-B*54:15', 'HLA-B*54:16',
                                     'HLA-B*54:17', 'HLA-B*54:18', 'HLA-B*54:19', 'HLA-B*54:20', 'HLA-B*54:21', 'HLA-B*54:22', 'HLA-B*54:23',
                                     'HLA-B*55:01', 'HLA-B*55:02', 'HLA-B*55:03',
                                     'HLA-B*55:04', 'HLA-B*55:05', 'HLA-B*55:07', 'HLA-B*55:08', 'HLA-B*55:09', 'HLA-B*55:10', 'HLA-B*55:11',
                                     'HLA-B*55:12', 'HLA-B*55:13', 'HLA-B*55:14',
                                     'HLA-B*55:15', 'HLA-B*55:16', 'HLA-B*55:17', 'HLA-B*55:18', 'HLA-B*55:19', 'HLA-B*55:20', 'HLA-B*55:21',
                                     'HLA-B*55:22', 'HLA-B*55:23', 'HLA-B*55:24',
                                     'HLA-B*55:25', 'HLA-B*55:26', 'HLA-B*55:27', 'HLA-B*55:28', 'HLA-B*55:29', 'HLA-B*55:30', 'HLA-B*55:31',
                                     'HLA-B*55:32', 'HLA-B*55:33', 'HLA-B*55:34',
                                     'HLA-B*55:35', 'HLA-B*55:36', 'HLA-B*55:37', 'HLA-B*55:38', 'HLA-B*55:39', 'HLA-B*55:40', 'HLA-B*55:41',
                                     'HLA-B*55:42', 'HLA-B*55:43', 'HLA-B*56:01',
                                     'HLA-B*56:02', 'HLA-B*56:03', 'HLA-B*56:04', 'HLA-B*56:05', 'HLA-B*56:06', 'HLA-B*56:07', 'HLA-B*56:08',
                                     'HLA-B*56:09', 'HLA-B*56:10', 'HLA-B*56:11',
                                     'HLA-B*56:12', 'HLA-B*56:13', 'HLA-B*56:14', 'HLA-B*56:15', 'HLA-B*56:16', 'HLA-B*56:17', 'HLA-B*56:18',
                                     'HLA-B*56:20', 'HLA-B*56:21', 'HLA-B*56:22',
                                     'HLA-B*56:23', 'HLA-B*56:24', 'HLA-B*56:25', 'HLA-B*56:26', 'HLA-B*56:27', 'HLA-B*56:29', 'HLA-B*57:01',
                                     'HLA-B*57:02', 'HLA-B*57:03', 'HLA-B*57:04',
                                     'HLA-B*57:05', 'HLA-B*57:06', 'HLA-B*57:07', 'HLA-B*57:08', 'HLA-B*57:09', 'HLA-B*57:10', 'HLA-B*57:11',
                                     'HLA-B*57:12', 'HLA-B*57:13', 'HLA-B*57:14',
                                     'HLA-B*57:15', 'HLA-B*57:16', 'HLA-B*57:17', 'HLA-B*57:18', 'HLA-B*57:19', 'HLA-B*57:20', 'HLA-B*57:21',
                                     'HLA-B*57:22', 'HLA-B*57:23', 'HLA-B*57:24',
                                     'HLA-B*57:25', 'HLA-B*57:26', 'HLA-B*57:27', 'HLA-B*57:29', 'HLA-B*57:30', 'HLA-B*57:31', 'HLA-B*57:32',
                                     'HLA-B*58:01', 'HLA-B*58:02', 'HLA-B*58:04',
                                     'HLA-B*58:05', 'HLA-B*58:06', 'HLA-B*58:07', 'HLA-B*58:08', 'HLA-B*58:09', 'HLA-B*58:11', 'HLA-B*58:12',
                                     'HLA-B*58:13', 'HLA-B*58:14', 'HLA-B*58:15',
                                     'HLA-B*58:16', 'HLA-B*58:18', 'HLA-B*58:19', 'HLA-B*58:20', 'HLA-B*58:21', 'HLA-B*58:22', 'HLA-B*58:23',
                                     'HLA-B*58:24', 'HLA-B*58:25', 'HLA-B*58:26',
                                     'HLA-B*58:27', 'HLA-B*58:28', 'HLA-B*58:29', 'HLA-B*58:30', 'HLA-B*59:01', 'HLA-B*59:02', 'HLA-B*59:03',
                                     'HLA-B*59:04', 'HLA-B*59:05', 'HLA-B*67:01',
                                     'HLA-B*67:02', 'HLA-B*73:01', 'HLA-B*73:02', 'HLA-B*78:01', 'HLA-B*78:02', 'HLA-B*78:03', 'HLA-B*78:04',
                                     'HLA-B*78:05', 'HLA-B*78:06', 'HLA-B*78:07',
                                     'HLA-B*81:01', 'HLA-B*81:02', 'HLA-B*81:03', 'HLA-B*81:05', 'HLA-B*82:01', 'HLA-B*82:02', 'HLA-B*82:03',
                                     'HLA-B*83:01', 'HLA-C*01:02', 'HLA-C*01:03',
                                     'HLA-C*01:04', 'HLA-C*01:05', 'HLA-C*01:06', 'HLA-C*01:07', 'HLA-C*01:08', 'HLA-C*01:09', 'HLA-C*01:10',
                                     'HLA-C*01:11', 'HLA-C*01:12', 'HLA-C*01:13',
                                     'HLA-C*01:14', 'HLA-C*01:15', 'HLA-C*01:16', 'HLA-C*01:17', 'HLA-C*01:18', 'HLA-C*01:19', 'HLA-C*01:20',
                                     'HLA-C*01:21', 'HLA-C*01:22', 'HLA-C*01:23',
                                     'HLA-C*01:24', 'HLA-C*01:25', 'HLA-C*01:26', 'HLA-C*01:27', 'HLA-C*01:28', 'HLA-C*01:29', 'HLA-C*01:30',
                                     'HLA-C*01:31', 'HLA-C*01:32', 'HLA-C*01:33',
                                     'HLA-C*01:34', 'HLA-C*01:35', 'HLA-C*01:36', 'HLA-C*01:38', 'HLA-C*01:39', 'HLA-C*01:40', 'HLA-C*02:02',
                                     'HLA-C*02:03', 'HLA-C*02:04', 'HLA-C*02:05',
                                     'HLA-C*02:06', 'HLA-C*02:07', 'HLA-C*02:08', 'HLA-C*02:09', 'HLA-C*02:10', 'HLA-C*02:11', 'HLA-C*02:12',
                                     'HLA-C*02:13', 'HLA-C*02:14', 'HLA-C*02:15',
                                     'HLA-C*02:16', 'HLA-C*02:17', 'HLA-C*02:18', 'HLA-C*02:19', 'HLA-C*02:20', 'HLA-C*02:21', 'HLA-C*02:22',
                                     'HLA-C*02:23', 'HLA-C*02:24', 'HLA-C*02:26',
                                     'HLA-C*02:27', 'HLA-C*02:28', 'HLA-C*02:29', 'HLA-C*02:30', 'HLA-C*02:31', 'HLA-C*02:32', 'HLA-C*02:33',
                                     'HLA-C*02:34', 'HLA-C*02:35', 'HLA-C*02:36',
                                     'HLA-C*02:37', 'HLA-C*02:39', 'HLA-C*02:40', 'HLA-C*03:01', 'HLA-C*03:02', 'HLA-C*03:03', 'HLA-C*03:04',
                                     'HLA-C*03:05', 'HLA-C*03:06', 'HLA-C*03:07',
                                     'HLA-C*03:08', 'HLA-C*03:09', 'HLA-C*03:10', 'HLA-C*03:11', 'HLA-C*03:12', 'HLA-C*03:13', 'HLA-C*03:14',
                                     'HLA-C*03:15', 'HLA-C*03:16', 'HLA-C*03:17',
                                     'HLA-C*03:18', 'HLA-C*03:19', 'HLA-C*03:21', 'HLA-C*03:23', 'HLA-C*03:24', 'HLA-C*03:25', 'HLA-C*03:26',
                                     'HLA-C*03:27', 'HLA-C*03:28', 'HLA-C*03:29',
                                     'HLA-C*03:30', 'HLA-C*03:31', 'HLA-C*03:32', 'HLA-C*03:33', 'HLA-C*03:34', 'HLA-C*03:35', 'HLA-C*03:36',
                                     'HLA-C*03:37', 'HLA-C*03:38', 'HLA-C*03:39',
                                     'HLA-C*03:40', 'HLA-C*03:41', 'HLA-C*03:42', 'HLA-C*03:43', 'HLA-C*03:44', 'HLA-C*03:45', 'HLA-C*03:46',
                                     'HLA-C*03:47', 'HLA-C*03:48', 'HLA-C*03:49',
                                     'HLA-C*03:50', 'HLA-C*03:51', 'HLA-C*03:52', 'HLA-C*03:53', 'HLA-C*03:54', 'HLA-C*03:55', 'HLA-C*03:56',
                                     'HLA-C*03:57', 'HLA-C*03:58', 'HLA-C*03:59',
                                     'HLA-C*03:60', 'HLA-C*03:61', 'HLA-C*03:62', 'HLA-C*03:63', 'HLA-C*03:64', 'HLA-C*03:65', 'HLA-C*03:66',
                                     'HLA-C*03:67', 'HLA-C*03:68', 'HLA-C*03:69',
                                     'HLA-C*03:70', 'HLA-C*03:71', 'HLA-C*03:72', 'HLA-C*03:73', 'HLA-C*03:74', 'HLA-C*03:75', 'HLA-C*03:76',
                                     'HLA-C*03:77', 'HLA-C*03:78', 'HLA-C*03:79',
                                     'HLA-C*03:80', 'HLA-C*03:81', 'HLA-C*03:82', 'HLA-C*03:83', 'HLA-C*03:84', 'HLA-C*03:85', 'HLA-C*03:86',
                                     'HLA-C*03:87', 'HLA-C*03:88', 'HLA-C*03:89',
                                     'HLA-C*03:90', 'HLA-C*03:91', 'HLA-C*03:92', 'HLA-C*03:93', 'HLA-C*03:94', 'HLA-C*04:01', 'HLA-C*04:03',
                                     'HLA-C*04:04', 'HLA-C*04:05', 'HLA-C*04:06',
                                     'HLA-C*04:07', 'HLA-C*04:08', 'HLA-C*04:10', 'HLA-C*04:11', 'HLA-C*04:12', 'HLA-C*04:13', 'HLA-C*04:14',
                                     'HLA-C*04:15', 'HLA-C*04:16', 'HLA-C*04:17',
                                     'HLA-C*04:18', 'HLA-C*04:19', 'HLA-C*04:20', 'HLA-C*04:23', 'HLA-C*04:24', 'HLA-C*04:25', 'HLA-C*04:26',
                                     'HLA-C*04:27', 'HLA-C*04:28', 'HLA-C*04:29',
                                     'HLA-C*04:30', 'HLA-C*04:31', 'HLA-C*04:32', 'HLA-C*04:33', 'HLA-C*04:34', 'HLA-C*04:35', 'HLA-C*04:36',
                                     'HLA-C*04:37', 'HLA-C*04:38', 'HLA-C*04:39',
                                     'HLA-C*04:40', 'HLA-C*04:41', 'HLA-C*04:42', 'HLA-C*04:43', 'HLA-C*04:44', 'HLA-C*04:45', 'HLA-C*04:46',
                                     'HLA-C*04:47', 'HLA-C*04:48', 'HLA-C*04:49',
                                     'HLA-C*04:50', 'HLA-C*04:51', 'HLA-C*04:52', 'HLA-C*04:53', 'HLA-C*04:54', 'HLA-C*04:55', 'HLA-C*04:56',
                                     'HLA-C*04:57', 'HLA-C*04:58', 'HLA-C*04:60',
                                     'HLA-C*04:61', 'HLA-C*04:62', 'HLA-C*04:63', 'HLA-C*04:64', 'HLA-C*04:65', 'HLA-C*04:66', 'HLA-C*04:67',
                                     'HLA-C*04:68', 'HLA-C*04:69', 'HLA-C*04:70',
                                     'HLA-C*05:01', 'HLA-C*05:03', 'HLA-C*05:04', 'HLA-C*05:05', 'HLA-C*05:06', 'HLA-C*05:08', 'HLA-C*05:09',
                                     'HLA-C*05:10', 'HLA-C*05:11', 'HLA-C*05:12',
                                     'HLA-C*05:13', 'HLA-C*05:14', 'HLA-C*05:15', 'HLA-C*05:16', 'HLA-C*05:17', 'HLA-C*05:18', 'HLA-C*05:19',
                                     'HLA-C*05:20', 'HLA-C*05:21', 'HLA-C*05:22',
                                     'HLA-C*05:23', 'HLA-C*05:24', 'HLA-C*05:25', 'HLA-C*05:26', 'HLA-C*05:27', 'HLA-C*05:28', 'HLA-C*05:29',
                                     'HLA-C*05:30', 'HLA-C*05:31', 'HLA-C*05:32',
                                     'HLA-C*05:33', 'HLA-C*05:34', 'HLA-C*05:35', 'HLA-C*05:36', 'HLA-C*05:37', 'HLA-C*05:38', 'HLA-C*05:39',
                                     'HLA-C*05:40', 'HLA-C*05:41', 'HLA-C*05:42',
                                     'HLA-C*05:43', 'HLA-C*05:44', 'HLA-C*05:45', 'HLA-C*06:02', 'HLA-C*06:03', 'HLA-C*06:04', 'HLA-C*06:05',
                                     'HLA-C*06:06', 'HLA-C*06:07', 'HLA-C*06:08',
                                     'HLA-C*06:09', 'HLA-C*06:10', 'HLA-C*06:11', 'HLA-C*06:12', 'HLA-C*06:13', 'HLA-C*06:14', 'HLA-C*06:15',
                                     'HLA-C*06:17', 'HLA-C*06:18', 'HLA-C*06:19',
                                     'HLA-C*06:20', 'HLA-C*06:21', 'HLA-C*06:22', 'HLA-C*06:23', 'HLA-C*06:24', 'HLA-C*06:25', 'HLA-C*06:26',
                                     'HLA-C*06:27', 'HLA-C*06:28', 'HLA-C*06:29',
                                     'HLA-C*06:30', 'HLA-C*06:31', 'HLA-C*06:32', 'HLA-C*06:33', 'HLA-C*06:34', 'HLA-C*06:35', 'HLA-C*06:36',
                                     'HLA-C*06:37', 'HLA-C*06:38', 'HLA-C*06:39',
                                     'HLA-C*06:40', 'HLA-C*06:41', 'HLA-C*06:42', 'HLA-C*06:43', 'HLA-C*06:44', 'HLA-C*06:45', 'HLA-C*07:01',
                                     'HLA-C*07:02', 'HLA-C*07:03', 'HLA-C*07:04',
                                     'HLA-C*07:05', 'HLA-C*07:06', 'HLA-C*07:07', 'HLA-C*07:08', 'HLA-C*07:09', 'HLA-C*07:10', 'HLA-C*07:11',
                                     'HLA-C*07:12', 'HLA-C*07:13', 'HLA-C*07:14',
                                     'HLA-C*07:15', 'HLA-C*07:16', 'HLA-C*07:17', 'HLA-C*07:18', 'HLA-C*07:19', 'HLA-C*07:20', 'HLA-C*07:21',
                                     'HLA-C*07:22', 'HLA-C*07:23', 'HLA-C*07:24',
                                     'HLA-C*07:25', 'HLA-C*07:26', 'HLA-C*07:27', 'HLA-C*07:28', 'HLA-C*07:29', 'HLA-C*07:30', 'HLA-C*07:31',
                                     'HLA-C*07:35', 'HLA-C*07:36', 'HLA-C*07:37',
                                     'HLA-C*07:38', 'HLA-C*07:39', 'HLA-C*07:40', 'HLA-C*07:41', 'HLA-C*07:42', 'HLA-C*07:43', 'HLA-C*07:44',
                                     'HLA-C*07:45', 'HLA-C*07:46', 'HLA-C*07:47',
                                     'HLA-C*07:48', 'HLA-C*07:49', 'HLA-C*07:50', 'HLA-C*07:51', 'HLA-C*07:52', 'HLA-C*07:53', 'HLA-C*07:54',
                                     'HLA-C*07:56', 'HLA-C*07:57', 'HLA-C*07:58',
                                     'HLA-C*07:59', 'HLA-C*07:60', 'HLA-C*07:62', 'HLA-C*07:63', 'HLA-C*07:64', 'HLA-C*07:65', 'HLA-C*07:66',
                                     'HLA-C*07:67', 'HLA-C*07:68', 'HLA-C*07:69',
                                     'HLA-C*07:70', 'HLA-C*07:71', 'HLA-C*07:72', 'HLA-C*07:73', 'HLA-C*07:74', 'HLA-C*07:75', 'HLA-C*07:76',
                                     'HLA-C*07:77', 'HLA-C*07:78', 'HLA-C*07:79',
                                     'HLA-C*07:80', 'HLA-C*07:81', 'HLA-C*07:82', 'HLA-C*07:83', 'HLA-C*07:84', 'HLA-C*07:85', 'HLA-C*07:86',
                                     'HLA-C*07:87', 'HLA-C*07:88', 'HLA-C*07:89',
                                     'HLA-C*07:90', 'HLA-C*07:91', 'HLA-C*07:92', 'HLA-C*07:93', 'HLA-C*07:94', 'HLA-C*07:95', 'HLA-C*07:96',
                                     'HLA-C*07:97', 'HLA-C*07:99', 'HLA-C*07:100',
                                     'HLA-C*07:101', 'HLA-C*07:102', 'HLA-C*07:103', 'HLA-C*07:105', 'HLA-C*07:106', 'HLA-C*07:107', 'HLA-C*07:108',
                                     'HLA-C*07:109', 'HLA-C*07:110',
                                     'HLA-C*07:111', 'HLA-C*07:112', 'HLA-C*07:113', 'HLA-C*07:114', 'HLA-C*07:115', 'HLA-C*07:116', 'HLA-C*07:117',
                                     'HLA-C*07:118', 'HLA-C*07:119',
                                     'HLA-C*07:120', 'HLA-C*07:122', 'HLA-C*07:123', 'HLA-C*07:124', 'HLA-C*07:125', 'HLA-C*07:126', 'HLA-C*07:127',
                                     'HLA-C*07:128', 'HLA-C*07:129',
                                     'HLA-C*07:130', 'HLA-C*07:131', 'HLA-C*07:132', 'HLA-C*07:133', 'HLA-C*07:134', 'HLA-C*07:135', 'HLA-C*07:136',
                                     'HLA-C*07:137', 'HLA-C*07:138',
                                     'HLA-C*07:139', 'HLA-C*07:140', 'HLA-C*07:141', 'HLA-C*07:142', 'HLA-C*07:143', 'HLA-C*07:144', 'HLA-C*07:145',
                                     'HLA-C*07:146', 'HLA-C*07:147',
                                     'HLA-C*07:148', 'HLA-C*07:149', 'HLA-C*08:01', 'HLA-C*08:02', 'HLA-C*08:03', 'HLA-C*08:04', 'HLA-C*08:05',
                                     'HLA-C*08:06', 'HLA-C*08:07', 'HLA-C*08:08',
                                     'HLA-C*08:09', 'HLA-C*08:10', 'HLA-C*08:11', 'HLA-C*08:12', 'HLA-C*08:13', 'HLA-C*08:14', 'HLA-C*08:15',
                                     'HLA-C*08:16', 'HLA-C*08:17', 'HLA-C*08:18',
                                     'HLA-C*08:19', 'HLA-C*08:20', 'HLA-C*08:21', 'HLA-C*08:22', 'HLA-C*08:23', 'HLA-C*08:24', 'HLA-C*08:25',
                                     'HLA-C*08:27', 'HLA-C*08:28', 'HLA-C*08:29',
                                     'HLA-C*08:30', 'HLA-C*08:31', 'HLA-C*08:32', 'HLA-C*08:33', 'HLA-C*08:34', 'HLA-C*08:35', 'HLA-C*12:02',
                                     'HLA-C*12:03', 'HLA-C*12:04', 'HLA-C*12:05',
                                     'HLA-C*12:06', 'HLA-C*12:07', 'HLA-C*12:08', 'HLA-C*12:09', 'HLA-C*12:10', 'HLA-C*12:11', 'HLA-C*12:12',
                                     'HLA-C*12:13', 'HLA-C*12:14', 'HLA-C*12:15',
                                     'HLA-C*12:16', 'HLA-C*12:17', 'HLA-C*12:18', 'HLA-C*12:19', 'HLA-C*12:20', 'HLA-C*12:21', 'HLA-C*12:22',
                                     'HLA-C*12:23', 'HLA-C*12:24', 'HLA-C*12:25',
                                     'HLA-C*12:26', 'HLA-C*12:27', 'HLA-C*12:28', 'HLA-C*12:29', 'HLA-C*12:30', 'HLA-C*12:31', 'HLA-C*12:32',
                                     'HLA-C*12:33', 'HLA-C*12:34', 'HLA-C*12:35',
                                     'HLA-C*12:36', 'HLA-C*12:37', 'HLA-C*12:38', 'HLA-C*12:40', 'HLA-C*12:41', 'HLA-C*12:43', 'HLA-C*12:44',
                                     'HLA-C*14:02', 'HLA-C*14:03', 'HLA-C*14:04',
                                     'HLA-C*14:05', 'HLA-C*14:06', 'HLA-C*14:08', 'HLA-C*14:09', 'HLA-C*14:10', 'HLA-C*14:11', 'HLA-C*14:12',
                                     'HLA-C*14:13', 'HLA-C*14:14', 'HLA-C*14:15',
                                     'HLA-C*14:16', 'HLA-C*14:17', 'HLA-C*14:18', 'HLA-C*14:19', 'HLA-C*14:20', 'HLA-C*15:02', 'HLA-C*15:03',
                                     'HLA-C*15:04', 'HLA-C*15:05', 'HLA-C*15:06',
                                     'HLA-C*15:07', 'HLA-C*15:08', 'HLA-C*15:09', 'HLA-C*15:10', 'HLA-C*15:11', 'HLA-C*15:12', 'HLA-C*15:13',
                                     'HLA-C*15:15', 'HLA-C*15:16', 'HLA-C*15:17',
                                     'HLA-C*15:18', 'HLA-C*15:19', 'HLA-C*15:20', 'HLA-C*15:21', 'HLA-C*15:22', 'HLA-C*15:23', 'HLA-C*15:24',
                                     'HLA-C*15:25', 'HLA-C*15:26', 'HLA-C*15:27',
                                     'HLA-C*15:28', 'HLA-C*15:29', 'HLA-C*15:30', 'HLA-C*15:31', 'HLA-C*15:33', 'HLA-C*15:34', 'HLA-C*15:35',
                                     'HLA-C*16:01', 'HLA-C*16:02', 'HLA-C*16:04',
                                     'HLA-C*16:06', 'HLA-C*16:07', 'HLA-C*16:08', 'HLA-C*16:09', 'HLA-C*16:10', 'HLA-C*16:11', 'HLA-C*16:12',
                                     'HLA-C*16:13', 'HLA-C*16:14', 'HLA-C*16:15',
                                     'HLA-C*16:17', 'HLA-C*16:18', 'HLA-C*16:19', 'HLA-C*16:20', 'HLA-C*16:21', 'HLA-C*16:22', 'HLA-C*16:23',
                                     'HLA-C*16:24', 'HLA-C*16:25', 'HLA-C*16:26',
                                     'HLA-C*17:01', 'HLA-C*17:02', 'HLA-C*17:03', 'HLA-C*17:04', 'HLA-C*17:05', 'HLA-C*17:06', 'HLA-C*17:07',
                                     'HLA-C*18:01', 'HLA-C*18:02', 'HLA-C*18:03',
                                     'HLA-G*01:01', 'HLA-G*01:02', 'HLA-G*01:03', 'HLA-G*01:04', 'HLA-G*01:06', 'HLA-G*01:07', 'HLA-G*01:08',
                                     'HLA-G*01:09', 'HLA-E*01:01',
                                     'H2-Db', 'H2-Dd', 'H2-Kb', 'H2-Kd', 'H2-Kk', 'H2-Ld'])
    __version = "1.1"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    @property
    def command(self):
        """
        Defines the commandline call for external tool
        """
        return self.__command

    @property
    def supportedLength(self):
        """
        A list of supported :class:`~Fred2.Core.Peptide.Peptide` lengths
        """
        return self.__supported_length

    @property
    def supportedAlleles(self):
        """
        A list of supported :class:`~Fred2.Core.Allele.Allele`
        """
        return self.__supported_alleles

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    def _represent(self, allele):
        """
        Internal function transforming an allele object into its representative string
        :param allele: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: str
        """
        if isinstance(allele, MouseAllele):
            return "H-2-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)
        else:
            return "HLA-%s%s:%s" % (allele.locus, allele.supertype, allele.subtype)

    def convert_alleles(self, alleles):
        """
        Converts :class:`~Fred2.Core.Allele.Allele` into the internal :class:`~Fred2.Core.Allele.Allele` representation
        of the predictor and returns a string representation

        :param alleles: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`
        :rtype: list(str)
        """
        return [self._represent(a) for a in alleles]

    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        result = defaultdict(defaultdict)
        with open(file, "r") as f:
            for row in f:
                if row[0] in ["#", "-"] or row.strip() == "" or "pos" in row:
                    continue
                else:
                    s = row.split()
                    result[s[1].replace("*", "")][s[2]] = float(s[4])
        return result

    def get_external_version(self, path=None):
        """
        Returns the external version of the tool by executing
        >{command} --version

        might be dependent on the method and has to be overwritten
        therefore it is declared abstract to enforce the user to
        overwrite the method. The function in the base class can be called
        with super()

        :param str path: Optional specification of executable path if deviant from :attr:`elf.__command`
        :return: The external version of the tool or None if tool does not support versioning
        :rtype: str
        """
        return None

    def prepare_input(self, input, file):
        """
        Prepares input for external tools and writes them to file in the specific format

        No return value!

        :param: list(str) input: The :class:`~Fred2.Core.Peptide.Peptide` sequences to write into _file
        :param File file: File-handler to input file for external tool
        """
        file.write("\n".join(input))


class NetCTLpan_1_1(AExternalEpitopePrediction):
    """
    Interface for NetCTLpan 1.1.

    .. note::

        NetCTLpan - Pan-specific MHC class I epitope predictions Stranzl T., Larsen M. V., Lundegaard C., Nielsen M.
        Immunogenetics. 2010 Apr 9. [Epub ahead of print]
    """
    __name = "netctlpan"
    __command = "netCTLpan -f {peptides} -a {alleles} {options} > {out}"
    __supported_length = frozenset([8, 9, 10, 11])
    __alleles = frozenset(
        ['HLA-A*01:01', 'HLA-A*01:02', 'HLA-A*01:03', 'HLA-A*01:06', 'HLA-A*01:07', 'HLA-A*01:08', 'HLA-A*01:09', 'HLA-A*01:10', 'HLA-A*01:12',
         'HLA-A*01:13', 'HLA-A*01:14', 'HLA-A*01:17', 'HLA-A*01:19', 'HLA-A*01:20', 'HLA-A*01:21', 'HLA-A*01:23', 'HLA-A*01:24', 'HLA-A*01:25',
         'HLA-A*01:26', 'HLA-A*01:28', 'HLA-A*01:29', 'HLA-A*01:30', 'HLA-A*01:32', 'HLA-A*01:33', 'HLA-A*01:35', 'HLA-A*01:36', 'HLA-A*01:37',
         'HLA-A*01:38', 'HLA-A*01:39', 'HLA-A*01:40', 'HLA-A*01:41', 'HLA-A*01:42', 'HLA-A*01:43', 'HLA-A*01:44', 'HLA-A*01:45', 'HLA-A*01:46',
         'HLA-A*01:47', 'HLA-A*01:48', 'HLA-A*01:49', 'HLA-A*01:50', 'HLA-A*01:51', 'HLA-A*01:54', 'HLA-A*01:55', 'HLA-A*01:58', 'HLA-A*01:59',
         'HLA-A*01:60', 'HLA-A*01:61', 'HLA-A*01:62', 'HLA-A*01:63', 'HLA-A*01:64', 'HLA-A*01:65', 'HLA-A*01:66', 'HLA-A*02:01', 'HLA-A*02:02',
         'HLA-A*02:03', 'HLA-A*02:04', 'HLA-A*02:05', 'HLA-A*02:06', 'HLA-A*02:07', 'HLA-A*02:08', 'HLA-A*02:09', 'HLA-A*02:10', 'HLA-A*02:101',
         'HLA-A*02:102', 'HLA-A*02:103', 'HLA-A*02:104', 'HLA-A*02:105', 'HLA-A*02:106', 'HLA-A*02:107', 'HLA-A*02:108', 'HLA-A*02:109',
         'HLA-A*02:11', 'HLA-A*02:110', 'HLA-A*02:111', 'HLA-A*02:112', 'HLA-A*02:114', 'HLA-A*02:115', 'HLA-A*02:116', 'HLA-A*02:117',
         'HLA-A*02:118', 'HLA-A*02:119', 'HLA-A*02:12', 'HLA-A*02:120', 'HLA-A*02:121', 'HLA-A*02:122', 'HLA-A*02:123', 'HLA-A*02:124',
         'HLA-A*02:126', 'HLA-A*02:127', 'HLA-A*02:128', 'HLA-A*02:129', 'HLA-A*02:13', 'HLA-A*02:130', 'HLA-A*02:131', 'HLA-A*02:132',
         'HLA-A*02:133', 'HLA-A*02:134', 'HLA-A*02:135', 'HLA-A*02:136', 'HLA-A*02:137', 'HLA-A*02:138', 'HLA-A*02:139', 'HLA-A*02:14',
         'HLA-A*02:140', 'HLA-A*02:141', 'HLA-A*02:142', 'HLA-A*02:143', 'HLA-A*02:144', 'HLA-A*02:145', 'HLA-A*02:146', 'HLA-A*02:147',
         'HLA-A*02:148', 'HLA-A*02:149', 'HLA-A*02:150', 'HLA-A*02:151', 'HLA-A*02:152', 'HLA-A*02:153', 'HLA-A*02:154', 'HLA-A*02:155',
         'HLA-A*02:156', 'HLA-A*02:157', 'HLA-A*02:158', 'HLA-A*02:159', 'HLA-A*02:16', 'HLA-A*02:160', 'HLA-A*02:161', 'HLA-A*02:162',
         'HLA-A*02:163', 'HLA-A*02:164', 'HLA-A*02:165', 'HLA-A*02:166', 'HLA-A*02:167', 'HLA-A*02:168', 'HLA-A*02:169', 'HLA-A*02:17',
         'HLA-A*02:170', 'HLA-A*02:171', 'HLA-A*02:172', 'HLA-A*02:173', 'HLA-A*02:174', 'HLA-A*02:175', 'HLA-A*02:176', 'HLA-A*02:177',
         'HLA-A*02:178', 'HLA-A*02:179', 'HLA-A*02:18', 'HLA-A*02:180', 'HLA-A*02:181', 'HLA-A*02:182', 'HLA-A*02:183', 'HLA-A*02:184',
         'HLA-A*02:185', 'HLA-A*02:186', 'HLA-A*02:187', 'HLA-A*02:188', 'HLA-A*02:189', 'HLA-A*02:19', 'HLA-A*02:190', 'HLA-A*02:191',
         'HLA-A*02:192', 'HLA-A*02:193', 'HLA-A*02:194', 'HLA-A*02:195', 'HLA-A*02:196', 'HLA-A*02:197', 'HLA-A*02:198', 'HLA-A*02:199',
         'HLA-A*02:20', 'HLA-A*02:200', 'HLA-A*02:201', 'HLA-A*02:202', 'HLA-A*02:203', 'HLA-A*02:204', 'HLA-A*02:205', 'HLA-A*02:206',
         'HLA-A*02:207', 'HLA-A*02:208', 'HLA-A*02:209', 'HLA-A*02:21', 'HLA-A*02:210', 'HLA-A*02:211', 'HLA-A*02:212', 'HLA-A*02:213',
         'HLA-A*02:214', 'HLA-A*02:215', 'HLA-A*02:216', 'HLA-A*02:217', 'HLA-A*02:218', 'HLA-A*02:219', 'HLA-A*02:22', 'HLA-A*02:220',
         'HLA-A*02:221', 'HLA-A*02:224', 'HLA-A*02:228', 'HLA-A*02:229', 'HLA-A*02:230', 'HLA-A*02:231', 'HLA-A*02:232', 'HLA-A*02:233',
         'HLA-A*02:234', 'HLA-A*02:235', 'HLA-A*02:236', 'HLA-A*02:237', 'HLA-A*02:238', 'HLA-A*02:239', 'HLA-A*02:24', 'HLA-A*02:240',
         'HLA-A*02:241', 'HLA-A*02:242', 'HLA-A*02:243', 'HLA-A*02:244', 'HLA-A*02:245', 'HLA-A*02:246', 'HLA-A*02:247', 'HLA-A*02:248',
         'HLA-A*02:249', 'HLA-A*02:25', 'HLA-A*02:251', 'HLA-A*02:252', 'HLA-A*02:253', 'HLA-A*02:254', 'HLA-A*02:255', 'HLA-A*02:256',
         'HLA-A*02:257', 'HLA-A*02:258', 'HLA-A*02:259', 'HLA-A*02:26', 'HLA-A*02:260', 'HLA-A*02:261', 'HLA-A*02:262', 'HLA-A*02:263',
         'HLA-A*02:264', 'HLA-A*02:265', 'HLA-A*02:266', 'HLA-A*02:27', 'HLA-A*02:28', 'HLA-A*02:29', 'HLA-A*02:30', 'HLA-A*02:31', 'HLA-A*02:33',
         'HLA-A*02:34', 'HLA-A*02:35', 'HLA-A*02:36', 'HLA-A*02:37', 'HLA-A*02:38', 'HLA-A*02:39', 'HLA-A*02:40', 'HLA-A*02:41', 'HLA-A*02:42',
         'HLA-A*02:44', 'HLA-A*02:45', 'HLA-A*02:46', 'HLA-A*02:47', 'HLA-A*02:48', 'HLA-A*02:49', 'HLA-A*02:50', 'HLA-A*02:51', 'HLA-A*02:52',
         'HLA-A*02:54', 'HLA-A*02:55', 'HLA-A*02:56', 'HLA-A*02:57', 'HLA-A*02:58', 'HLA-A*02:59', 'HLA-A*02:60', 'HLA-A*02:61', 'HLA-A*02:62',
         'HLA-A*02:63', 'HLA-A*02:64', 'HLA-A*02:65', 'HLA-A*02:66', 'HLA-A*02:67', 'HLA-A*02:68', 'HLA-A*02:69', 'HLA-A*02:70', 'HLA-A*02:71',
         'HLA-A*02:72', 'HLA-A*02:73', 'HLA-A*02:74', 'HLA-A*02:75', 'HLA-A*02:76', 'HLA-A*02:77', 'HLA-A*02:78', 'HLA-A*02:79', 'HLA-A*02:80',
         'HLA-A*02:81', 'HLA-A*02:84', 'HLA-A*02:85', 'HLA-A*02:86', 'HLA-A*02:87', 'HLA-A*02:89', 'HLA-A*02:90', 'HLA-A*02:91', 'HLA-A*02:92',
         'HLA-A*02:93', 'HLA-A*02:95', 'HLA-A*02:96', 'HLA-A*02:97', 'HLA-A*02:99', 'HLA-A*03:01', 'HLA-A*03:02', 'HLA-A*03:04', 'HLA-A*03:05',
         'HLA-A*03:06', 'HLA-A*03:07', 'HLA-A*03:08', 'HLA-A*03:09', 'HLA-A*03:10', 'HLA-A*03:12', 'HLA-A*03:13', 'HLA-A*03:14', 'HLA-A*03:15',
         'HLA-A*03:16', 'HLA-A*03:17', 'HLA-A*03:18', 'HLA-A*03:19', 'HLA-A*03:20', 'HLA-A*03:22', 'HLA-A*03:23', 'HLA-A*03:24', 'HLA-A*03:25',
         'HLA-A*03:26', 'HLA-A*03:27', 'HLA-A*03:28', 'HLA-A*03:29', 'HLA-A*03:30', 'HLA-A*03:31', 'HLA-A*03:32', 'HLA-A*03:33', 'HLA-A*03:34',
         'HLA-A*03:35', 'HLA-A*03:37', 'HLA-A*03:38', 'HLA-A*03:39', 'HLA-A*03:40', 'HLA-A*03:41', 'HLA-A*03:42', 'HLA-A*03:43', 'HLA-A*03:44',
         'HLA-A*03:45', 'HLA-A*03:46', 'HLA-A*03:47', 'HLA-A*03:48', 'HLA-A*03:49', 'HLA-A*03:50', 'HLA-A*03:51', 'HLA-A*03:52', 'HLA-A*03:53',
         'HLA-A*03:54', 'HLA-A*03:55', 'HLA-A*03:56', 'HLA-A*03:57', 'HLA-A*03:58', 'HLA-A*03:59', 'HLA-A*03:60', 'HLA-A*03:61', 'HLA-A*03:62',
         'HLA-A*03:63', 'HLA-A*03:64', 'HLA-A*03:65', 'HLA-A*03:66', 'HLA-A*03:67', 'HLA-A*03:70', 'HLA-A*03:71', 'HLA-A*03:72', 'HLA-A*03:73',
         'HLA-A*03:74', 'HLA-A*03:75', 'HLA-A*03:76', 'HLA-A*03:77', 'HLA-A*03:78', 'HLA-A*03:79', 'HLA-A*03:80', 'HLA-A*03:81', 'HLA-A*03:82',
         'HLA-A*11:01', 'HLA-A*11:02', 'HLA-A*11:03', 'HLA-A*11:04', 'HLA-A*11:05', 'HLA-A*11:06', 'HLA-A*11:07', 'HLA-A*11:08', 'HLA-A*11:09',
         'HLA-A*11:10', 'HLA-A*11:11', 'HLA-A*11:12', 'HLA-A*11:13', 'HLA-A*11:14', 'HLA-A*11:15', 'HLA-A*11:16', 'HLA-A*11:17', 'HLA-A*11:18',
         'HLA-A*11:19', 'HLA-A*11:20', 'HLA-A*11:22', 'HLA-A*11:23', 'HLA-A*11:24', 'HLA-A*11:25', 'HLA-A*11:26', 'HLA-A*11:27', 'HLA-A*11:29',
         'HLA-A*11:30', 'HLA-A*11:31', 'HLA-A*11:32', 'HLA-A*11:33', 'HLA-A*11:34', 'HLA-A*11:35', 'HLA-A*11:36', 'HLA-A*11:37', 'HLA-A*11:38',
         'HLA-A*11:39', 'HLA-A*11:40', 'HLA-A*11:41', 'HLA-A*11:42', 'HLA-A*11:43', 'HLA-A*11:44', 'HLA-A*11:45', 'HLA-A*11:46', 'HLA-A*11:47',
         'HLA-A*11:48', 'HLA-A*11:49', 'HLA-A*11:51', 'HLA-A*11:53', 'HLA-A*11:54', 'HLA-A*11:55', 'HLA-A*11:56', 'HLA-A*11:57', 'HLA-A*11:58',
         'HLA-A*11:59', 'HLA-A*11:60', 'HLA-A*11:61', 'HLA-A*11:62', 'HLA-A*11:63', 'HLA-A*11:64', 'HLA-A*23:01', 'HLA-A*23:02', 'HLA-A*23:03',
         'HLA-A*23:04', 'HLA-A*23:05', 'HLA-A*23:06', 'HLA-A*23:09', 'HLA-A*23:10', 'HLA-A*23:12', 'HLA-A*23:13', 'HLA-A*23:14', 'HLA-A*23:15',
         'HLA-A*23:16', 'HLA-A*23:17', 'HLA-A*23:18', 'HLA-A*23:20', 'HLA-A*23:21', 'HLA-A*23:22', 'HLA-A*23:23', 'HLA-A*23:24', 'HLA-A*23:25',
         'HLA-A*23:26', 'HLA-A*24:02', 'HLA-A*24:03', 'HLA-A*24:04', 'HLA-A*24:05', 'HLA-A*24:06', 'HLA-A*24:07', 'HLA-A*24:08', 'HLA-A*24:10',
         'HLA-A*24:100', 'HLA-A*24:101', 'HLA-A*24:102', 'HLA-A*24:103', 'HLA-A*24:104', 'HLA-A*24:105', 'HLA-A*24:106', 'HLA-A*24:107',
         'HLA-A*24:108', 'HLA-A*24:109', 'HLA-A*24:110', 'HLA-A*24:111', 'HLA-A*24:112', 'HLA-A*24:113', 'HLA-A*24:114', 'HLA-A*24:115',
         'HLA-A*24:116', 'HLA-A*24:117', 'HLA-A*24:118', 'HLA-A*24:119', 'HLA-A*24:120', 'HLA-A*24:121', 'HLA-A*24:122', 'HLA-A*24:123',
         'HLA-A*24:124', 'HLA-A*24:125', 'HLA-A*24:126', 'HLA-A*24:127', 'HLA-A*24:128', 'HLA-A*24:129', 'HLA-A*24:13', 'HLA-A*24:130',
         'HLA-A*24:131', 'HLA-A*24:133', 'HLA-A*24:134', 'HLA-A*24:135', 'HLA-A*24:136', 'HLA-A*24:137', 'HLA-A*24:138', 'HLA-A*24:139',
         'HLA-A*24:14', 'HLA-A*24:140', 'HLA-A*24:141', 'HLA-A*24:142', 'HLA-A*24:143', 'HLA-A*24:144', 'HLA-A*24:15', 'HLA-A*24:17', 'HLA-A*24:18',
         'HLA-A*24:19', 'HLA-A*24:20', 'HLA-A*24:21', 'HLA-A*24:22', 'HLA-A*24:23', 'HLA-A*24:24', 'HLA-A*24:25', 'HLA-A*24:26', 'HLA-A*24:27',
         'HLA-A*24:28', 'HLA-A*24:29', 'HLA-A*24:30', 'HLA-A*24:31', 'HLA-A*24:32', 'HLA-A*24:33', 'HLA-A*24:34', 'HLA-A*24:35', 'HLA-A*24:37',
         'HLA-A*24:38', 'HLA-A*24:39', 'HLA-A*24:41', 'HLA-A*24:42', 'HLA-A*24:43', 'HLA-A*24:44', 'HLA-A*24:46', 'HLA-A*24:47', 'HLA-A*24:49',
         'HLA-A*24:50', 'HLA-A*24:51', 'HLA-A*24:52', 'HLA-A*24:53', 'HLA-A*24:54', 'HLA-A*24:55', 'HLA-A*24:56', 'HLA-A*24:57', 'HLA-A*24:58',
         'HLA-A*24:59', 'HLA-A*24:61', 'HLA-A*24:62', 'HLA-A*24:63', 'HLA-A*24:64', 'HLA-A*24:66', 'HLA-A*24:67', 'HLA-A*24:68', 'HLA-A*24:69',
         'HLA-A*24:70', 'HLA-A*24:71', 'HLA-A*24:72', 'HLA-A*24:73', 'HLA-A*24:74', 'HLA-A*24:75', 'HLA-A*24:76', 'HLA-A*24:77', 'HLA-A*24:78',
         'HLA-A*24:79', 'HLA-A*24:80', 'HLA-A*24:81', 'HLA-A*24:82', 'HLA-A*24:85', 'HLA-A*24:87', 'HLA-A*24:88', 'HLA-A*24:89', 'HLA-A*24:91',
         'HLA-A*24:92', 'HLA-A*24:93', 'HLA-A*24:94', 'HLA-A*24:95', 'HLA-A*24:96', 'HLA-A*24:97', 'HLA-A*24:98', 'HLA-A*24:99', 'HLA-A*25:01',
         'HLA-A*25:02', 'HLA-A*25:03', 'HLA-A*25:04', 'HLA-A*25:05', 'HLA-A*25:06', 'HLA-A*25:07', 'HLA-A*25:08', 'HLA-A*25:09', 'HLA-A*25:10',
         'HLA-A*25:11', 'HLA-A*25:13', 'HLA-A*26:01', 'HLA-A*26:02', 'HLA-A*26:03', 'HLA-A*26:04', 'HLA-A*26:05', 'HLA-A*26:06', 'HLA-A*26:07',
         'HLA-A*26:08', 'HLA-A*26:09', 'HLA-A*26:10', 'HLA-A*26:12', 'HLA-A*26:13', 'HLA-A*26:14', 'HLA-A*26:15', 'HLA-A*26:16', 'HLA-A*26:17',
         'HLA-A*26:18', 'HLA-A*26:19', 'HLA-A*26:20', 'HLA-A*26:21', 'HLA-A*26:22', 'HLA-A*26:23', 'HLA-A*26:24', 'HLA-A*26:26', 'HLA-A*26:27',
         'HLA-A*26:28', 'HLA-A*26:29', 'HLA-A*26:30', 'HLA-A*26:31', 'HLA-A*26:32', 'HLA-A*26:33', 'HLA-A*26:34', 'HLA-A*26:35', 'HLA-A*26:36',
         'HLA-A*26:37', 'HLA-A*26:38', 'HLA-A*26:39', 'HLA-A*26:40', 'HLA-A*26:41', 'HLA-A*26:42', 'HLA-A*26:43', 'HLA-A*26:45', 'HLA-A*26:46',
         'HLA-A*26:47', 'HLA-A*26:48', 'HLA-A*26:49', 'HLA-A*26:50', 'HLA-A*29:01', 'HLA-A*29:02', 'HLA-A*29:03', 'HLA-A*29:04', 'HLA-A*29:05',
         'HLA-A*29:06', 'HLA-A*29:07', 'HLA-A*29:09', 'HLA-A*29:10', 'HLA-A*29:11', 'HLA-A*29:12', 'HLA-A*29:13', 'HLA-A*29:14', 'HLA-A*29:15',
         'HLA-A*29:16', 'HLA-A*29:17', 'HLA-A*29:18', 'HLA-A*29:19', 'HLA-A*29:20', 'HLA-A*29:21', 'HLA-A*29:22', 'HLA-A*30:01', 'HLA-A*30:02',
         'HLA-A*30:03', 'HLA-A*30:04', 'HLA-A*30:06', 'HLA-A*30:07', 'HLA-A*30:08', 'HLA-A*30:09', 'HLA-A*30:10', 'HLA-A*30:11', 'HLA-A*30:12',
         'HLA-A*30:13', 'HLA-A*30:15', 'HLA-A*30:16', 'HLA-A*30:17', 'HLA-A*30:18', 'HLA-A*30:19', 'HLA-A*30:20', 'HLA-A*30:22', 'HLA-A*30:23',
         'HLA-A*30:24', 'HLA-A*30:25', 'HLA-A*30:26', 'HLA-A*30:28', 'HLA-A*30:29', 'HLA-A*30:30', 'HLA-A*30:31', 'HLA-A*30:32', 'HLA-A*30:33',
         'HLA-A*30:34', 'HLA-A*30:35', 'HLA-A*30:36', 'HLA-A*30:37', 'HLA-A*30:38', 'HLA-A*30:39', 'HLA-A*30:40', 'HLA-A*30:41', 'HLA-A*31:01',
         'HLA-A*31:02', 'HLA-A*31:03', 'HLA-A*31:04', 'HLA-A*31:05', 'HLA-A*31:06', 'HLA-A*31:07', 'HLA-A*31:08', 'HLA-A*31:09', 'HLA-A*31:10',
         'HLA-A*31:11', 'HLA-A*31:12', 'HLA-A*31:13', 'HLA-A*31:15', 'HLA-A*31:16', 'HLA-A*31:17', 'HLA-A*31:18', 'HLA-A*31:19', 'HLA-A*31:20',
         'HLA-A*31:21', 'HLA-A*31:22', 'HLA-A*31:23', 'HLA-A*31:24', 'HLA-A*31:25', 'HLA-A*31:26', 'HLA-A*31:27', 'HLA-A*31:28', 'HLA-A*31:29',
         'HLA-A*31:30', 'HLA-A*31:31', 'HLA-A*31:32', 'HLA-A*31:33', 'HLA-A*31:34', 'HLA-A*31:35', 'HLA-A*31:36', 'HLA-A*31:37', 'HLA-A*32:01',
         'HLA-A*32:02', 'HLA-A*32:03', 'HLA-A*32:04', 'HLA-A*32:05', 'HLA-A*32:06', 'HLA-A*32:07', 'HLA-A*32:08', 'HLA-A*32:09', 'HLA-A*32:10',
         'HLA-A*32:12', 'HLA-A*32:13', 'HLA-A*32:14', 'HLA-A*32:15', 'HLA-A*32:16', 'HLA-A*32:17', 'HLA-A*32:18', 'HLA-A*32:20', 'HLA-A*32:21',
         'HLA-A*32:22', 'HLA-A*32:23', 'HLA-A*32:24', 'HLA-A*32:25', 'HLA-A*33:01', 'HLA-A*33:03', 'HLA-A*33:04', 'HLA-A*33:05', 'HLA-A*33:06',
         'HLA-A*33:07', 'HLA-A*33:08', 'HLA-A*33:09', 'HLA-A*33:10', 'HLA-A*33:11', 'HLA-A*33:12', 'HLA-A*33:13', 'HLA-A*33:14', 'HLA-A*33:15',
         'HLA-A*33:16', 'HLA-A*33:17', 'HLA-A*33:18', 'HLA-A*33:19', 'HLA-A*33:20', 'HLA-A*33:21', 'HLA-A*33:22', 'HLA-A*33:23', 'HLA-A*33:24',
         'HLA-A*33:25', 'HLA-A*33:26', 'HLA-A*33:27', 'HLA-A*33:28', 'HLA-A*33:29', 'HLA-A*33:30', 'HLA-A*33:31', 'HLA-A*34:01', 'HLA-A*34:02',
         'HLA-A*34:03', 'HLA-A*34:04', 'HLA-A*34:05', 'HLA-A*34:06', 'HLA-A*34:07', 'HLA-A*34:08', 'HLA-A*36:01', 'HLA-A*36:02', 'HLA-A*36:03',
         'HLA-A*36:04', 'HLA-A*36:05', 'HLA-A*43:01', 'HLA-A*66:01', 'HLA-A*66:02', 'HLA-A*66:03', 'HLA-A*66:04', 'HLA-A*66:05', 'HLA-A*66:06',
         'HLA-A*66:07', 'HLA-A*66:08', 'HLA-A*66:09', 'HLA-A*66:10', 'HLA-A*66:11', 'HLA-A*66:12', 'HLA-A*66:13', 'HLA-A*66:14', 'HLA-A*66:15',
         'HLA-A*68:01', 'HLA-A*68:02', 'HLA-A*68:03', 'HLA-A*68:04', 'HLA-A*68:05', 'HLA-A*68:06', 'HLA-A*68:07', 'HLA-A*68:08', 'HLA-A*68:09',
         'HLA-A*68:10', 'HLA-A*68:12', 'HLA-A*68:13', 'HLA-A*68:14', 'HLA-A*68:15', 'HLA-A*68:16', 'HLA-A*68:17', 'HLA-A*68:19', 'HLA-A*68:20',
         'HLA-A*68:21', 'HLA-A*68:22', 'HLA-A*68:23', 'HLA-A*68:24', 'HLA-A*68:25', 'HLA-A*68:26', 'HLA-A*68:27', 'HLA-A*68:28', 'HLA-A*68:29',
         'HLA-A*68:30', 'HLA-A*68:31', 'HLA-A*68:32', 'HLA-A*68:33', 'HLA-A*68:34', 'HLA-A*68:35', 'HLA-A*68:36', 'HLA-A*68:37', 'HLA-A*68:38',
         'HLA-A*68:39', 'HLA-A*68:40', 'HLA-A*68:41', 'HLA-A*68:42', 'HLA-A*68:43', 'HLA-A*68:44', 'HLA-A*68:45', 'HLA-A*68:46', 'HLA-A*68:47',
         'HLA-A*68:48', 'HLA-A*68:50', 'HLA-A*68:51', 'HLA-A*68:52', 'HLA-A*68:53', 'HLA-A*68:54', 'HLA-A*69:01', 'HLA-A*74:01', 'HLA-A*74:02',
         'HLA-A*74:03', 'HLA-A*74:04', 'HLA-A*74:05', 'HLA-A*74:06', 'HLA-A*74:07', 'HLA-A*74:08', 'HLA-A*74:09', 'HLA-A*74:10', 'HLA-A*74:11',
         'HLA-A*74:13', 'HLA-A*80:01', 'HLA-A*80:02', 'HLA-B*07:02', 'HLA-B*07:03', 'HLA-B*07:04', 'HLA-B*07:05', 'HLA-B*07:06', 'HLA-B*07:07',
         'HLA-B*07:08', 'HLA-B*07:09', 'HLA-B*07:10', 'HLA-B*07:100', 'HLA-B*07:101', 'HLA-B*07:102', 'HLA-B*07:103', 'HLA-B*07:104',
         'HLA-B*07:105', 'HLA-B*07:106', 'HLA-B*07:107', 'HLA-B*07:108', 'HLA-B*07:109', 'HLA-B*07:11', 'HLA-B*07:110', 'HLA-B*07:112',
         'HLA-B*07:113', 'HLA-B*07:114', 'HLA-B*07:115', 'HLA-B*07:12', 'HLA-B*07:13', 'HLA-B*07:14', 'HLA-B*07:15', 'HLA-B*07:16', 'HLA-B*07:17',
         'HLA-B*07:18', 'HLA-B*07:19', 'HLA-B*07:20', 'HLA-B*07:21', 'HLA-B*07:22', 'HLA-B*07:23', 'HLA-B*07:24', 'HLA-B*07:25', 'HLA-B*07:26',
         'HLA-B*07:27', 'HLA-B*07:28', 'HLA-B*07:29', 'HLA-B*07:30', 'HLA-B*07:31', 'HLA-B*07:32', 'HLA-B*07:33', 'HLA-B*07:34', 'HLA-B*07:35',
         'HLA-B*07:36', 'HLA-B*07:37', 'HLA-B*07:38', 'HLA-B*07:39', 'HLA-B*07:40', 'HLA-B*07:41', 'HLA-B*07:42', 'HLA-B*07:43', 'HLA-B*07:44',
         'HLA-B*07:45', 'HLA-B*07:46', 'HLA-B*07:47', 'HLA-B*07:48', 'HLA-B*07:50', 'HLA-B*07:51', 'HLA-B*07:52', 'HLA-B*07:53', 'HLA-B*07:54',
         'HLA-B*07:55', 'HLA-B*07:56', 'HLA-B*07:57', 'HLA-B*07:58', 'HLA-B*07:59', 'HLA-B*07:60', 'HLA-B*07:61', 'HLA-B*07:62', 'HLA-B*07:63',
         'HLA-B*07:64', 'HLA-B*07:65', 'HLA-B*07:66', 'HLA-B*07:68', 'HLA-B*07:69', 'HLA-B*07:70', 'HLA-B*07:71', 'HLA-B*07:72', 'HLA-B*07:73',
         'HLA-B*07:74', 'HLA-B*07:75', 'HLA-B*07:76', 'HLA-B*07:77', 'HLA-B*07:78', 'HLA-B*07:79', 'HLA-B*07:80', 'HLA-B*07:81', 'HLA-B*07:82',
         'HLA-B*07:83', 'HLA-B*07:84', 'HLA-B*07:85', 'HLA-B*07:86', 'HLA-B*07:87', 'HLA-B*07:88', 'HLA-B*07:89', 'HLA-B*07:90', 'HLA-B*07:91',
         'HLA-B*07:92', 'HLA-B*07:93', 'HLA-B*07:94', 'HLA-B*07:95', 'HLA-B*07:96', 'HLA-B*07:97', 'HLA-B*07:98', 'HLA-B*07:99', 'HLA-B*08:01',
         'HLA-B*08:02', 'HLA-B*08:03', 'HLA-B*08:04', 'HLA-B*08:05', 'HLA-B*08:07', 'HLA-B*08:09', 'HLA-B*08:10', 'HLA-B*08:11', 'HLA-B*08:12',
         'HLA-B*08:13', 'HLA-B*08:14', 'HLA-B*08:15', 'HLA-B*08:16', 'HLA-B*08:17', 'HLA-B*08:18', 'HLA-B*08:20', 'HLA-B*08:21', 'HLA-B*08:22',
         'HLA-B*08:23', 'HLA-B*08:24', 'HLA-B*08:25', 'HLA-B*08:26', 'HLA-B*08:27', 'HLA-B*08:28', 'HLA-B*08:29', 'HLA-B*08:31', 'HLA-B*08:32',
         'HLA-B*08:33', 'HLA-B*08:34', 'HLA-B*08:35', 'HLA-B*08:36', 'HLA-B*08:37', 'HLA-B*08:38', 'HLA-B*08:39', 'HLA-B*08:40', 'HLA-B*08:41',
         'HLA-B*08:42', 'HLA-B*08:43', 'HLA-B*08:44', 'HLA-B*08:45', 'HLA-B*08:46', 'HLA-B*08:47', 'HLA-B*08:48', 'HLA-B*08:49', 'HLA-B*08:50',
         'HLA-B*08:51', 'HLA-B*08:52', 'HLA-B*08:53', 'HLA-B*08:54', 'HLA-B*08:55', 'HLA-B*08:56', 'HLA-B*08:57', 'HLA-B*08:58', 'HLA-B*08:59',
         'HLA-B*08:60', 'HLA-B*08:61', 'HLA-B*08:62', 'HLA-B*13:01', 'HLA-B*13:02', 'HLA-B*13:03', 'HLA-B*13:04', 'HLA-B*13:06', 'HLA-B*13:09',
         'HLA-B*13:10', 'HLA-B*13:11', 'HLA-B*13:12', 'HLA-B*13:13', 'HLA-B*13:14', 'HLA-B*13:15', 'HLA-B*13:16', 'HLA-B*13:17', 'HLA-B*13:18',
         'HLA-B*13:19', 'HLA-B*13:20', 'HLA-B*13:21', 'HLA-B*13:22', 'HLA-B*13:23', 'HLA-B*13:25', 'HLA-B*13:26', 'HLA-B*13:27', 'HLA-B*13:28',
         'HLA-B*13:29', 'HLA-B*13:30', 'HLA-B*13:31', 'HLA-B*13:32', 'HLA-B*13:33', 'HLA-B*13:34', 'HLA-B*13:35', 'HLA-B*13:36', 'HLA-B*13:37',
         'HLA-B*13:38', 'HLA-B*13:39', 'HLA-B*14:01', 'HLA-B*14:02', 'HLA-B*14:03', 'HLA-B*14:04', 'HLA-B*14:05', 'HLA-B*14:06', 'HLA-B*14:08',
         'HLA-B*14:09', 'HLA-B*14:10', 'HLA-B*14:11', 'HLA-B*14:12', 'HLA-B*14:13', 'HLA-B*14:14', 'HLA-B*14:15', 'HLA-B*14:16', 'HLA-B*14:17',
         'HLA-B*14:18', 'HLA-B*15:01', 'HLA-B*15:02', 'HLA-B*15:03', 'HLA-B*15:04', 'HLA-B*15:05', 'HLA-B*15:06', 'HLA-B*15:07', 'HLA-B*15:08',
         'HLA-B*15:09', 'HLA-B*15:10', 'HLA-B*15:101', 'HLA-B*15:102', 'HLA-B*15:103', 'HLA-B*15:104', 'HLA-B*15:105', 'HLA-B*15:106',
         'HLA-B*15:107', 'HLA-B*15:108', 'HLA-B*15:109', 'HLA-B*15:11', 'HLA-B*15:110', 'HLA-B*15:112', 'HLA-B*15:113', 'HLA-B*15:114',
         'HLA-B*15:115', 'HLA-B*15:116', 'HLA-B*15:117', 'HLA-B*15:118', 'HLA-B*15:119', 'HLA-B*15:12', 'HLA-B*15:120', 'HLA-B*15:121',
         'HLA-B*15:122', 'HLA-B*15:123', 'HLA-B*15:124', 'HLA-B*15:125', 'HLA-B*15:126', 'HLA-B*15:127', 'HLA-B*15:128', 'HLA-B*15:129',
         'HLA-B*15:13', 'HLA-B*15:131', 'HLA-B*15:132', 'HLA-B*15:133', 'HLA-B*15:134', 'HLA-B*15:135', 'HLA-B*15:136', 'HLA-B*15:137',
         'HLA-B*15:138', 'HLA-B*15:139', 'HLA-B*15:14', 'HLA-B*15:140', 'HLA-B*15:141', 'HLA-B*15:142', 'HLA-B*15:143', 'HLA-B*15:144',
         'HLA-B*15:145', 'HLA-B*15:146', 'HLA-B*15:147', 'HLA-B*15:148', 'HLA-B*15:15', 'HLA-B*15:150', 'HLA-B*15:151', 'HLA-B*15:152',
         'HLA-B*15:153', 'HLA-B*15:154', 'HLA-B*15:155', 'HLA-B*15:156', 'HLA-B*15:157', 'HLA-B*15:158', 'HLA-B*15:159', 'HLA-B*15:16',
         'HLA-B*15:160', 'HLA-B*15:161', 'HLA-B*15:162', 'HLA-B*15:163', 'HLA-B*15:164', 'HLA-B*15:165', 'HLA-B*15:166', 'HLA-B*15:167',
         'HLA-B*15:168', 'HLA-B*15:169', 'HLA-B*15:17', 'HLA-B*15:170', 'HLA-B*15:171', 'HLA-B*15:172', 'HLA-B*15:173', 'HLA-B*15:174',
         'HLA-B*15:175', 'HLA-B*15:176', 'HLA-B*15:177', 'HLA-B*15:178', 'HLA-B*15:179', 'HLA-B*15:18', 'HLA-B*15:180', 'HLA-B*15:183',
         'HLA-B*15:184', 'HLA-B*15:185', 'HLA-B*15:186', 'HLA-B*15:187', 'HLA-B*15:188', 'HLA-B*15:189', 'HLA-B*15:19', 'HLA-B*15:191',
         'HLA-B*15:192', 'HLA-B*15:193', 'HLA-B*15:194', 'HLA-B*15:195', 'HLA-B*15:196', 'HLA-B*15:197', 'HLA-B*15:198', 'HLA-B*15:199',
         'HLA-B*15:20', 'HLA-B*15:200', 'HLA-B*15:201', 'HLA-B*15:202', 'HLA-B*15:21', 'HLA-B*15:23', 'HLA-B*15:24', 'HLA-B*15:25', 'HLA-B*15:27',
         'HLA-B*15:28', 'HLA-B*15:29', 'HLA-B*15:30', 'HLA-B*15:31', 'HLA-B*15:32', 'HLA-B*15:33', 'HLA-B*15:34', 'HLA-B*15:35', 'HLA-B*15:36',
         'HLA-B*15:37', 'HLA-B*15:38', 'HLA-B*15:39', 'HLA-B*15:40', 'HLA-B*15:42', 'HLA-B*15:43', 'HLA-B*15:44', 'HLA-B*15:45', 'HLA-B*15:46',
         'HLA-B*15:47', 'HLA-B*15:48', 'HLA-B*15:49', 'HLA-B*15:50', 'HLA-B*15:51', 'HLA-B*15:52', 'HLA-B*15:53', 'HLA-B*15:54', 'HLA-B*15:55',
         'HLA-B*15:56', 'HLA-B*15:57', 'HLA-B*15:58', 'HLA-B*15:60', 'HLA-B*15:61', 'HLA-B*15:62', 'HLA-B*15:63', 'HLA-B*15:64', 'HLA-B*15:65',
         'HLA-B*15:66', 'HLA-B*15:67', 'HLA-B*15:68', 'HLA-B*15:69', 'HLA-B*15:70', 'HLA-B*15:71', 'HLA-B*15:72', 'HLA-B*15:73', 'HLA-B*15:74',
         'HLA-B*15:75', 'HLA-B*15:76', 'HLA-B*15:77', 'HLA-B*15:78', 'HLA-B*15:80', 'HLA-B*15:81', 'HLA-B*15:82', 'HLA-B*15:83', 'HLA-B*15:84',
         'HLA-B*15:85', 'HLA-B*15:86', 'HLA-B*15:87', 'HLA-B*15:88', 'HLA-B*15:89', 'HLA-B*15:90', 'HLA-B*15:91', 'HLA-B*15:92', 'HLA-B*15:93',
         'HLA-B*15:95', 'HLA-B*15:96', 'HLA-B*15:97', 'HLA-B*15:98', 'HLA-B*15:99', 'HLA-B*18:01', 'HLA-B*18:02', 'HLA-B*18:03', 'HLA-B*18:04',
         'HLA-B*18:05', 'HLA-B*18:06', 'HLA-B*18:07', 'HLA-B*18:08', 'HLA-B*18:09', 'HLA-B*18:10', 'HLA-B*18:11', 'HLA-B*18:12', 'HLA-B*18:13',
         'HLA-B*18:14', 'HLA-B*18:15', 'HLA-B*18:18', 'HLA-B*18:19', 'HLA-B*18:20', 'HLA-B*18:21', 'HLA-B*18:22', 'HLA-B*18:24', 'HLA-B*18:25',
         'HLA-B*18:26', 'HLA-B*18:27', 'HLA-B*18:28', 'HLA-B*18:29', 'HLA-B*18:30', 'HLA-B*18:31', 'HLA-B*18:32', 'HLA-B*18:33', 'HLA-B*18:34',
         'HLA-B*18:35', 'HLA-B*18:36', 'HLA-B*18:37', 'HLA-B*18:38', 'HLA-B*18:39', 'HLA-B*18:40', 'HLA-B*18:41', 'HLA-B*18:42', 'HLA-B*18:43',
         'HLA-B*18:44', 'HLA-B*18:45', 'HLA-B*18:46', 'HLA-B*18:47', 'HLA-B*18:48', 'HLA-B*18:49', 'HLA-B*18:50', 'HLA-B*27:01', 'HLA-B*27:02',
         'HLA-B*27:03', 'HLA-B*27:04', 'HLA-B*27:05', 'HLA-B*27:06', 'HLA-B*27:07', 'HLA-B*27:08', 'HLA-B*27:09', 'HLA-B*27:10', 'HLA-B*27:11',
         'HLA-B*27:12', 'HLA-B*27:13', 'HLA-B*27:14', 'HLA-B*27:15', 'HLA-B*27:16', 'HLA-B*27:17', 'HLA-B*27:18', 'HLA-B*27:19', 'HLA-B*27:20',
         'HLA-B*27:21', 'HLA-B*27:23', 'HLA-B*27:24', 'HLA-B*27:25', 'HLA-B*27:26', 'HLA-B*27:27', 'HLA-B*27:28', 'HLA-B*27:29', 'HLA-B*27:30',
         'HLA-B*27:31', 'HLA-B*27:32', 'HLA-B*27:33', 'HLA-B*27:34', 'HLA-B*27:35', 'HLA-B*27:36', 'HLA-B*27:37', 'HLA-B*27:38', 'HLA-B*27:39',
         'HLA-B*27:40', 'HLA-B*27:41', 'HLA-B*27:42', 'HLA-B*27:43', 'HLA-B*27:44', 'HLA-B*27:45', 'HLA-B*27:46', 'HLA-B*27:47', 'HLA-B*27:48',
         'HLA-B*27:49', 'HLA-B*27:50', 'HLA-B*27:51', 'HLA-B*27:52', 'HLA-B*27:53', 'HLA-B*27:54', 'HLA-B*27:55', 'HLA-B*27:56', 'HLA-B*27:57',
         'HLA-B*27:58', 'HLA-B*27:60', 'HLA-B*27:61', 'HLA-B*27:62', 'HLA-B*27:63', 'HLA-B*27:67', 'HLA-B*27:68', 'HLA-B*27:69', 'HLA-B*35:01',
         'HLA-B*35:02', 'HLA-B*35:03', 'HLA-B*35:04', 'HLA-B*35:05', 'HLA-B*35:06', 'HLA-B*35:07', 'HLA-B*35:08', 'HLA-B*35:09', 'HLA-B*35:10',
         'HLA-B*35:100', 'HLA-B*35:101', 'HLA-B*35:102', 'HLA-B*35:103', 'HLA-B*35:104', 'HLA-B*35:105', 'HLA-B*35:106', 'HLA-B*35:107',
         'HLA-B*35:108', 'HLA-B*35:109', 'HLA-B*35:11', 'HLA-B*35:110', 'HLA-B*35:111', 'HLA-B*35:112', 'HLA-B*35:113', 'HLA-B*35:114',
         'HLA-B*35:115', 'HLA-B*35:116', 'HLA-B*35:117', 'HLA-B*35:118', 'HLA-B*35:119', 'HLA-B*35:12', 'HLA-B*35:120', 'HLA-B*35:121',
         'HLA-B*35:122', 'HLA-B*35:123', 'HLA-B*35:124', 'HLA-B*35:125', 'HLA-B*35:126', 'HLA-B*35:127', 'HLA-B*35:128', 'HLA-B*35:13',
         'HLA-B*35:131', 'HLA-B*35:132', 'HLA-B*35:133', 'HLA-B*35:135', 'HLA-B*35:136', 'HLA-B*35:137', 'HLA-B*35:138', 'HLA-B*35:139',
         'HLA-B*35:14', 'HLA-B*35:140', 'HLA-B*35:141', 'HLA-B*35:142', 'HLA-B*35:143', 'HLA-B*35:144', 'HLA-B*35:15', 'HLA-B*35:16', 'HLA-B*35:17',
         'HLA-B*35:18', 'HLA-B*35:19', 'HLA-B*35:20', 'HLA-B*35:21', 'HLA-B*35:22', 'HLA-B*35:23', 'HLA-B*35:24', 'HLA-B*35:25', 'HLA-B*35:26',
         'HLA-B*35:27', 'HLA-B*35:28', 'HLA-B*35:29', 'HLA-B*35:30', 'HLA-B*35:31', 'HLA-B*35:32', 'HLA-B*35:33', 'HLA-B*35:34', 'HLA-B*35:35',
         'HLA-B*35:36', 'HLA-B*35:37', 'HLA-B*35:38', 'HLA-B*35:39', 'HLA-B*35:41', 'HLA-B*35:42', 'HLA-B*35:43', 'HLA-B*35:44', 'HLA-B*35:45',
         'HLA-B*35:46', 'HLA-B*35:47', 'HLA-B*35:48', 'HLA-B*35:49', 'HLA-B*35:50', 'HLA-B*35:51', 'HLA-B*35:52', 'HLA-B*35:54', 'HLA-B*35:55',
         'HLA-B*35:56', 'HLA-B*35:57', 'HLA-B*35:58', 'HLA-B*35:59', 'HLA-B*35:60', 'HLA-B*35:61', 'HLA-B*35:62', 'HLA-B*35:63', 'HLA-B*35:64',
         'HLA-B*35:66', 'HLA-B*35:67', 'HLA-B*35:68', 'HLA-B*35:69', 'HLA-B*35:70', 'HLA-B*35:71', 'HLA-B*35:72', 'HLA-B*35:74', 'HLA-B*35:75',
         'HLA-B*35:76', 'HLA-B*35:77', 'HLA-B*35:78', 'HLA-B*35:79', 'HLA-B*35:80', 'HLA-B*35:81', 'HLA-B*35:82', 'HLA-B*35:83', 'HLA-B*35:84',
         'HLA-B*35:85', 'HLA-B*35:86', 'HLA-B*35:87', 'HLA-B*35:88', 'HLA-B*35:89', 'HLA-B*35:90', 'HLA-B*35:91', 'HLA-B*35:92', 'HLA-B*35:93',
         'HLA-B*35:94', 'HLA-B*35:95', 'HLA-B*35:96', 'HLA-B*35:97', 'HLA-B*35:98', 'HLA-B*35:99', 'HLA-B*37:01', 'HLA-B*37:02', 'HLA-B*37:04',
         'HLA-B*37:05', 'HLA-B*37:06', 'HLA-B*37:07', 'HLA-B*37:08', 'HLA-B*37:09', 'HLA-B*37:10', 'HLA-B*37:11', 'HLA-B*37:12', 'HLA-B*37:13',
         'HLA-B*37:14', 'HLA-B*37:15', 'HLA-B*37:17', 'HLA-B*37:18', 'HLA-B*37:19', 'HLA-B*37:20', 'HLA-B*37:21', 'HLA-B*37:22', 'HLA-B*37:23',
         'HLA-B*38:01', 'HLA-B*38:02', 'HLA-B*38:03', 'HLA-B*38:04', 'HLA-B*38:05', 'HLA-B*38:06', 'HLA-B*38:07', 'HLA-B*38:08', 'HLA-B*38:09',
         'HLA-B*38:10', 'HLA-B*38:11', 'HLA-B*38:12', 'HLA-B*38:13', 'HLA-B*38:14', 'HLA-B*38:15', 'HLA-B*38:16', 'HLA-B*38:17', 'HLA-B*38:18',
         'HLA-B*38:19', 'HLA-B*38:20', 'HLA-B*38:21', 'HLA-B*38:22', 'HLA-B*38:23', 'HLA-B*39:01', 'HLA-B*39:02', 'HLA-B*39:03', 'HLA-B*39:04',
         'HLA-B*39:05', 'HLA-B*39:06', 'HLA-B*39:07', 'HLA-B*39:08', 'HLA-B*39:09', 'HLA-B*39:10', 'HLA-B*39:11', 'HLA-B*39:12', 'HLA-B*39:13',
         'HLA-B*39:14', 'HLA-B*39:15', 'HLA-B*39:16', 'HLA-B*39:17', 'HLA-B*39:18', 'HLA-B*39:19', 'HLA-B*39:20', 'HLA-B*39:22', 'HLA-B*39:23',
         'HLA-B*39:24', 'HLA-B*39:26', 'HLA-B*39:27', 'HLA-B*39:28', 'HLA-B*39:29', 'HLA-B*39:30', 'HLA-B*39:31', 'HLA-B*39:32', 'HLA-B*39:33',
         'HLA-B*39:34', 'HLA-B*39:35', 'HLA-B*39:36', 'HLA-B*39:37', 'HLA-B*39:39', 'HLA-B*39:41', 'HLA-B*39:42', 'HLA-B*39:43', 'HLA-B*39:44',
         'HLA-B*39:45', 'HLA-B*39:46', 'HLA-B*39:47', 'HLA-B*39:48', 'HLA-B*39:49', 'HLA-B*39:50', 'HLA-B*39:51', 'HLA-B*39:52', 'HLA-B*39:53',
         'HLA-B*39:54', 'HLA-B*39:55', 'HLA-B*39:56', 'HLA-B*39:57', 'HLA-B*39:58', 'HLA-B*39:59', 'HLA-B*39:60', 'HLA-B*40:01', 'HLA-B*40:02',
         'HLA-B*40:03', 'HLA-B*40:04', 'HLA-B*40:05', 'HLA-B*40:06', 'HLA-B*40:07', 'HLA-B*40:08', 'HLA-B*40:09', 'HLA-B*40:10', 'HLA-B*40:100',
         'HLA-B*40:101', 'HLA-B*40:102', 'HLA-B*40:103', 'HLA-B*40:104', 'HLA-B*40:105', 'HLA-B*40:106', 'HLA-B*40:107', 'HLA-B*40:108',
         'HLA-B*40:109', 'HLA-B*40:11', 'HLA-B*40:110', 'HLA-B*40:111', 'HLA-B*40:112', 'HLA-B*40:113', 'HLA-B*40:114', 'HLA-B*40:115',
         'HLA-B*40:116', 'HLA-B*40:117', 'HLA-B*40:119', 'HLA-B*40:12', 'HLA-B*40:120', 'HLA-B*40:121', 'HLA-B*40:122', 'HLA-B*40:123',
         'HLA-B*40:124', 'HLA-B*40:125', 'HLA-B*40:126', 'HLA-B*40:127', 'HLA-B*40:128', 'HLA-B*40:129', 'HLA-B*40:13', 'HLA-B*40:130',
         'HLA-B*40:131', 'HLA-B*40:132', 'HLA-B*40:134', 'HLA-B*40:135', 'HLA-B*40:136', 'HLA-B*40:137', 'HLA-B*40:138', 'HLA-B*40:139',
         'HLA-B*40:14', 'HLA-B*40:140', 'HLA-B*40:141', 'HLA-B*40:143', 'HLA-B*40:145', 'HLA-B*40:146', 'HLA-B*40:147', 'HLA-B*40:15',
         'HLA-B*40:16', 'HLA-B*40:18', 'HLA-B*40:19', 'HLA-B*40:20', 'HLA-B*40:21', 'HLA-B*40:23', 'HLA-B*40:24', 'HLA-B*40:25', 'HLA-B*40:26',
         'HLA-B*40:27', 'HLA-B*40:28', 'HLA-B*40:29', 'HLA-B*40:30', 'HLA-B*40:31', 'HLA-B*40:32', 'HLA-B*40:33', 'HLA-B*40:34', 'HLA-B*40:35',
         'HLA-B*40:36', 'HLA-B*40:37', 'HLA-B*40:38', 'HLA-B*40:39', 'HLA-B*40:40', 'HLA-B*40:42', 'HLA-B*40:43', 'HLA-B*40:44', 'HLA-B*40:45',
         'HLA-B*40:46', 'HLA-B*40:47', 'HLA-B*40:48', 'HLA-B*40:49', 'HLA-B*40:50', 'HLA-B*40:51', 'HLA-B*40:52', 'HLA-B*40:53', 'HLA-B*40:54',
         'HLA-B*40:55', 'HLA-B*40:56', 'HLA-B*40:57', 'HLA-B*40:58', 'HLA-B*40:59', 'HLA-B*40:60', 'HLA-B*40:61', 'HLA-B*40:62', 'HLA-B*40:63',
         'HLA-B*40:64', 'HLA-B*40:65', 'HLA-B*40:66', 'HLA-B*40:67', 'HLA-B*40:68', 'HLA-B*40:69', 'HLA-B*40:70', 'HLA-B*40:71', 'HLA-B*40:72',
         'HLA-B*40:73', 'HLA-B*40:74', 'HLA-B*40:75', 'HLA-B*40:76', 'HLA-B*40:77', 'HLA-B*40:78', 'HLA-B*40:79', 'HLA-B*40:80', 'HLA-B*40:81',
         'HLA-B*40:82', 'HLA-B*40:83', 'HLA-B*40:84', 'HLA-B*40:85', 'HLA-B*40:86', 'HLA-B*40:87', 'HLA-B*40:88', 'HLA-B*40:89', 'HLA-B*40:90',
         'HLA-B*40:91', 'HLA-B*40:92', 'HLA-B*40:93', 'HLA-B*40:94', 'HLA-B*40:95', 'HLA-B*40:96', 'HLA-B*40:97', 'HLA-B*40:98', 'HLA-B*40:99',
         'HLA-B*41:01', 'HLA-B*41:02', 'HLA-B*41:03', 'HLA-B*41:04', 'HLA-B*41:05', 'HLA-B*41:06', 'HLA-B*41:07', 'HLA-B*41:08', 'HLA-B*41:09',
         'HLA-B*41:10', 'HLA-B*41:11', 'HLA-B*41:12', 'HLA-B*42:01', 'HLA-B*42:02', 'HLA-B*42:04', 'HLA-B*42:05', 'HLA-B*42:06', 'HLA-B*42:07',
         'HLA-B*42:08', 'HLA-B*42:09', 'HLA-B*42:10', 'HLA-B*42:11', 'HLA-B*42:12', 'HLA-B*42:13', 'HLA-B*42:14', 'HLA-B*44:02', 'HLA-B*44:03',
         'HLA-B*44:04', 'HLA-B*44:05', 'HLA-B*44:06', 'HLA-B*44:07', 'HLA-B*44:08', 'HLA-B*44:09', 'HLA-B*44:10', 'HLA-B*44:100', 'HLA-B*44:101',
         'HLA-B*44:102', 'HLA-B*44:103', 'HLA-B*44:104', 'HLA-B*44:105', 'HLA-B*44:106', 'HLA-B*44:107', 'HLA-B*44:109', 'HLA-B*44:11',
         'HLA-B*44:110', 'HLA-B*44:12', 'HLA-B*44:13', 'HLA-B*44:14', 'HLA-B*44:15', 'HLA-B*44:16', 'HLA-B*44:17', 'HLA-B*44:18', 'HLA-B*44:20',
         'HLA-B*44:21', 'HLA-B*44:22', 'HLA-B*44:24', 'HLA-B*44:25', 'HLA-B*44:26', 'HLA-B*44:27', 'HLA-B*44:28', 'HLA-B*44:29', 'HLA-B*44:30',
         'HLA-B*44:31', 'HLA-B*44:32', 'HLA-B*44:33', 'HLA-B*44:34', 'HLA-B*44:35', 'HLA-B*44:36', 'HLA-B*44:37', 'HLA-B*44:38', 'HLA-B*44:39',
         'HLA-B*44:40', 'HLA-B*44:41', 'HLA-B*44:42', 'HLA-B*44:43', 'HLA-B*44:44', 'HLA-B*44:45', 'HLA-B*44:46', 'HLA-B*44:47', 'HLA-B*44:48',
         'HLA-B*44:49', 'HLA-B*44:50', 'HLA-B*44:51', 'HLA-B*44:53', 'HLA-B*44:54', 'HLA-B*44:55', 'HLA-B*44:57', 'HLA-B*44:59', 'HLA-B*44:60',
         'HLA-B*44:62', 'HLA-B*44:63', 'HLA-B*44:64', 'HLA-B*44:65', 'HLA-B*44:66', 'HLA-B*44:67', 'HLA-B*44:68', 'HLA-B*44:69', 'HLA-B*44:70',
         'HLA-B*44:71', 'HLA-B*44:72', 'HLA-B*44:73', 'HLA-B*44:74', 'HLA-B*44:75', 'HLA-B*44:76', 'HLA-B*44:77', 'HLA-B*44:78', 'HLA-B*44:79',
         'HLA-B*44:80', 'HLA-B*44:81', 'HLA-B*44:82', 'HLA-B*44:83', 'HLA-B*44:84', 'HLA-B*44:85', 'HLA-B*44:86', 'HLA-B*44:87', 'HLA-B*44:88',
         'HLA-B*44:89', 'HLA-B*44:90', 'HLA-B*44:91', 'HLA-B*44:92', 'HLA-B*44:93', 'HLA-B*44:94', 'HLA-B*44:95', 'HLA-B*44:96', 'HLA-B*44:97',
         'HLA-B*44:98', 'HLA-B*44:99', 'HLA-B*45:01', 'HLA-B*45:02', 'HLA-B*45:03', 'HLA-B*45:04', 'HLA-B*45:05', 'HLA-B*45:06', 'HLA-B*45:07',
         'HLA-B*45:08', 'HLA-B*45:09', 'HLA-B*45:10', 'HLA-B*45:11', 'HLA-B*45:12', 'HLA-B*46:01', 'HLA-B*46:02', 'HLA-B*46:03', 'HLA-B*46:04',
         'HLA-B*46:05', 'HLA-B*46:06', 'HLA-B*46:08', 'HLA-B*46:09', 'HLA-B*46:10', 'HLA-B*46:11', 'HLA-B*46:12', 'HLA-B*46:13', 'HLA-B*46:14',
         'HLA-B*46:16', 'HLA-B*46:17', 'HLA-B*46:18', 'HLA-B*46:19', 'HLA-B*46:20', 'HLA-B*46:21', 'HLA-B*46:22', 'HLA-B*46:23', 'HLA-B*46:24',
         'HLA-B*47:01', 'HLA-B*47:02', 'HLA-B*47:03', 'HLA-B*47:04', 'HLA-B*47:05', 'HLA-B*47:06', 'HLA-B*47:07', 'HLA-B*48:01', 'HLA-B*48:02',
         'HLA-B*48:03', 'HLA-B*48:04', 'HLA-B*48:05', 'HLA-B*48:06', 'HLA-B*48:07', 'HLA-B*48:08', 'HLA-B*48:09', 'HLA-B*48:10', 'HLA-B*48:11',
         'HLA-B*48:12', 'HLA-B*48:13', 'HLA-B*48:14', 'HLA-B*48:15', 'HLA-B*48:16', 'HLA-B*48:17', 'HLA-B*48:18', 'HLA-B*48:19', 'HLA-B*48:20',
         'HLA-B*48:21', 'HLA-B*48:22', 'HLA-B*48:23', 'HLA-B*49:01', 'HLA-B*49:02', 'HLA-B*49:03', 'HLA-B*49:04', 'HLA-B*49:05', 'HLA-B*49:06',
         'HLA-B*49:07', 'HLA-B*49:08', 'HLA-B*49:09', 'HLA-B*49:10', 'HLA-B*50:01', 'HLA-B*50:02', 'HLA-B*50:04', 'HLA-B*50:05', 'HLA-B*50:06',
         'HLA-B*50:07', 'HLA-B*50:08', 'HLA-B*50:09', 'HLA-B*51:01', 'HLA-B*51:02', 'HLA-B*51:03', 'HLA-B*51:04', 'HLA-B*51:05', 'HLA-B*51:06',
         'HLA-B*51:07', 'HLA-B*51:08', 'HLA-B*51:09', 'HLA-B*51:12', 'HLA-B*51:13', 'HLA-B*51:14', 'HLA-B*51:15', 'HLA-B*51:16', 'HLA-B*51:17',
         'HLA-B*51:18', 'HLA-B*51:19', 'HLA-B*51:20', 'HLA-B*51:21', 'HLA-B*51:22', 'HLA-B*51:23', 'HLA-B*51:24', 'HLA-B*51:26', 'HLA-B*51:28',
         'HLA-B*51:29', 'HLA-B*51:30', 'HLA-B*51:31', 'HLA-B*51:32', 'HLA-B*51:33', 'HLA-B*51:34', 'HLA-B*51:35', 'HLA-B*51:36', 'HLA-B*51:37',
         'HLA-B*51:38', 'HLA-B*51:39', 'HLA-B*51:40', 'HLA-B*51:42', 'HLA-B*51:43', 'HLA-B*51:45', 'HLA-B*51:46', 'HLA-B*51:48', 'HLA-B*51:49',
         'HLA-B*51:50', 'HLA-B*51:51', 'HLA-B*51:52', 'HLA-B*51:53', 'HLA-B*51:54', 'HLA-B*51:55', 'HLA-B*51:56', 'HLA-B*51:57', 'HLA-B*51:58',
         'HLA-B*51:59', 'HLA-B*51:60', 'HLA-B*51:61', 'HLA-B*51:62', 'HLA-B*51:63', 'HLA-B*51:64', 'HLA-B*51:65', 'HLA-B*51:66', 'HLA-B*51:67',
         'HLA-B*51:68', 'HLA-B*51:69', 'HLA-B*51:70', 'HLA-B*51:71', 'HLA-B*51:72', 'HLA-B*51:73', 'HLA-B*51:74', 'HLA-B*51:75', 'HLA-B*51:76',
         'HLA-B*51:77', 'HLA-B*51:78', 'HLA-B*51:79', 'HLA-B*51:80', 'HLA-B*51:81', 'HLA-B*51:82', 'HLA-B*51:83', 'HLA-B*51:84', 'HLA-B*51:85',
         'HLA-B*51:86', 'HLA-B*51:87', 'HLA-B*51:88', 'HLA-B*51:89', 'HLA-B*51:90', 'HLA-B*51:91', 'HLA-B*51:92', 'HLA-B*51:93', 'HLA-B*51:94',
         'HLA-B*51:95', 'HLA-B*51:96', 'HLA-B*52:01', 'HLA-B*52:02', 'HLA-B*52:03', 'HLA-B*52:04', 'HLA-B*52:05', 'HLA-B*52:06', 'HLA-B*52:07',
         'HLA-B*52:08', 'HLA-B*52:09', 'HLA-B*52:10', 'HLA-B*52:11', 'HLA-B*52:12', 'HLA-B*52:13', 'HLA-B*52:14', 'HLA-B*52:15', 'HLA-B*52:16',
         'HLA-B*52:17', 'HLA-B*52:18', 'HLA-B*52:19', 'HLA-B*52:20', 'HLA-B*52:21', 'HLA-B*53:01', 'HLA-B*53:02', 'HLA-B*53:03', 'HLA-B*53:04',
         'HLA-B*53:05', 'HLA-B*53:06', 'HLA-B*53:07', 'HLA-B*53:08', 'HLA-B*53:09', 'HLA-B*53:10', 'HLA-B*53:11', 'HLA-B*53:12', 'HLA-B*53:13',
         'HLA-B*53:14', 'HLA-B*53:15', 'HLA-B*53:16', 'HLA-B*53:17', 'HLA-B*53:18', 'HLA-B*53:19', 'HLA-B*53:20', 'HLA-B*53:21', 'HLA-B*53:22',
         'HLA-B*53:23', 'HLA-B*54:01', 'HLA-B*54:02', 'HLA-B*54:03', 'HLA-B*54:04', 'HLA-B*54:06', 'HLA-B*54:07', 'HLA-B*54:09', 'HLA-B*54:10',
         'HLA-B*54:11', 'HLA-B*54:12', 'HLA-B*54:13', 'HLA-B*54:14', 'HLA-B*54:15', 'HLA-B*54:16', 'HLA-B*54:17', 'HLA-B*54:18', 'HLA-B*54:19',
         'HLA-B*54:20', 'HLA-B*54:21', 'HLA-B*54:22', 'HLA-B*54:23', 'HLA-B*55:01', 'HLA-B*55:02', 'HLA-B*55:03', 'HLA-B*55:04', 'HLA-B*55:05',
         'HLA-B*55:07', 'HLA-B*55:08', 'HLA-B*55:09', 'HLA-B*55:10', 'HLA-B*55:11', 'HLA-B*55:12', 'HLA-B*55:13', 'HLA-B*55:14', 'HLA-B*55:15',
         'HLA-B*55:16', 'HLA-B*55:17', 'HLA-B*55:18', 'HLA-B*55:19', 'HLA-B*55:20', 'HLA-B*55:21', 'HLA-B*55:22', 'HLA-B*55:23', 'HLA-B*55:24',
         'HLA-B*55:25', 'HLA-B*55:26', 'HLA-B*55:27', 'HLA-B*55:28', 'HLA-B*55:29', 'HLA-B*55:30', 'HLA-B*55:31', 'HLA-B*55:32', 'HLA-B*55:33',
         'HLA-B*55:34', 'HLA-B*55:35', 'HLA-B*55:36', 'HLA-B*55:37', 'HLA-B*55:38', 'HLA-B*55:39', 'HLA-B*55:40', 'HLA-B*55:41', 'HLA-B*55:42',
         'HLA-B*55:43', 'HLA-B*56:01', 'HLA-B*56:02', 'HLA-B*56:03', 'HLA-B*56:04', 'HLA-B*56:05', 'HLA-B*56:06', 'HLA-B*56:07', 'HLA-B*56:08',
         'HLA-B*56:09', 'HLA-B*56:10', 'HLA-B*56:11', 'HLA-B*56:12', 'HLA-B*56:13', 'HLA-B*56:14', 'HLA-B*56:15', 'HLA-B*56:16', 'HLA-B*56:17',
         'HLA-B*56:18', 'HLA-B*56:20', 'HLA-B*56:21', 'HLA-B*56:22', 'HLA-B*56:23', 'HLA-B*56:24', 'HLA-B*56:25', 'HLA-B*56:26', 'HLA-B*56:27',
         'HLA-B*56:29', 'HLA-B*57:01', 'HLA-B*57:02', 'HLA-B*57:03', 'HLA-B*57:04', 'HLA-B*57:05', 'HLA-B*57:06', 'HLA-B*57:07', 'HLA-B*57:08',
         'HLA-B*57:09', 'HLA-B*57:10', 'HLA-B*57:11', 'HLA-B*57:12', 'HLA-B*57:13', 'HLA-B*57:14', 'HLA-B*57:15', 'HLA-B*57:16', 'HLA-B*57:17',
         'HLA-B*57:18', 'HLA-B*57:19', 'HLA-B*57:20', 'HLA-B*57:21', 'HLA-B*57:22', 'HLA-B*57:23', 'HLA-B*57:24', 'HLA-B*57:25', 'HLA-B*57:26',
         'HLA-B*57:27', 'HLA-B*57:29', 'HLA-B*57:30', 'HLA-B*57:31', 'HLA-B*57:32', 'HLA-B*58:01', 'HLA-B*58:02', 'HLA-B*58:04', 'HLA-B*58:05',
         'HLA-B*58:06', 'HLA-B*58:07', 'HLA-B*58:08', 'HLA-B*58:09', 'HLA-B*58:11', 'HLA-B*58:12', 'HLA-B*58:13', 'HLA-B*58:14', 'HLA-B*58:15',
         'HLA-B*58:16', 'HLA-B*58:18', 'HLA-B*58:19', 'HLA-B*58:20', 'HLA-B*58:21', 'HLA-B*58:22', 'HLA-B*58:23', 'HLA-B*58:24', 'HLA-B*58:25',
         'HLA-B*58:26', 'HLA-B*58:27', 'HLA-B*58:28', 'HLA-B*58:29', 'HLA-B*58:30', 'HLA-B*59:01', 'HLA-B*59:02', 'HLA-B*59:03', 'HLA-B*59:04',
         'HLA-B*59:05', 'HLA-B*67:01', 'HLA-B*67:02', 'HLA-B*73:01', 'HLA-B*73:02', 'HLA-B*78:01', 'HLA-B*78:02', 'HLA-B*78:03', 'HLA-B*78:04',
         'HLA-B*78:05', 'HLA-B*78:06', 'HLA-B*78:07', 'HLA-B*81:01', 'HLA-B*81:02', 'HLA-B*81:03', 'HLA-B*81:05', 'HLA-B*82:01', 'HLA-B*82:02',
         'HLA-B*82:03', 'HLA-B*83:01', 'HLA-C*01:02', 'HLA-C*01:03', 'HLA-C*01:04', 'HLA-C*01:05', 'HLA-C*01:06', 'HLA-C*01:07', 'HLA-C*01:08',
         'HLA-C*01:09', 'HLA-C*01:10', 'HLA-C*01:11', 'HLA-C*01:12', 'HLA-C*01:13', 'HLA-C*01:14', 'HLA-C*01:15', 'HLA-C*01:16', 'HLA-C*01:17',
         'HLA-C*01:18', 'HLA-C*01:19', 'HLA-C*01:20', 'HLA-C*01:21', 'HLA-C*01:22', 'HLA-C*01:23', 'HLA-C*01:24', 'HLA-C*01:25', 'HLA-C*01:26',
         'HLA-C*01:27', 'HLA-C*01:28', 'HLA-C*01:29', 'HLA-C*01:30', 'HLA-C*01:31', 'HLA-C*01:32', 'HLA-C*01:33', 'HLA-C*01:34', 'HLA-C*01:35',
         'HLA-C*01:36', 'HLA-C*01:38', 'HLA-C*01:39', 'HLA-C*01:40', 'HLA-C*02:02', 'HLA-C*02:03', 'HLA-C*02:04', 'HLA-C*02:05', 'HLA-C*02:06',
         'HLA-C*02:07', 'HLA-C*02:08', 'HLA-C*02:09', 'HLA-C*02:10', 'HLA-C*02:11', 'HLA-C*02:12', 'HLA-C*02:13', 'HLA-C*02:14', 'HLA-C*02:15',
         'HLA-C*02:16', 'HLA-C*02:17', 'HLA-C*02:18', 'HLA-C*02:19', 'HLA-C*02:20', 'HLA-C*02:21', 'HLA-C*02:22', 'HLA-C*02:23', 'HLA-C*02:24',
         'HLA-C*02:26', 'HLA-C*02:27', 'HLA-C*02:28', 'HLA-C*02:29', 'HLA-C*02:30', 'HLA-C*02:31', 'HLA-C*02:32', 'HLA-C*02:33', 'HLA-C*02:34',
         'HLA-C*02:35', 'HLA-C*02:36', 'HLA-C*02:37', 'HLA-C*02:39', 'HLA-C*02:40', 'HLA-C*03:01', 'HLA-C*03:02', 'HLA-C*03:03', 'HLA-C*03:04',
         'HLA-C*03:05', 'HLA-C*03:06', 'HLA-C*03:07', 'HLA-C*03:08', 'HLA-C*03:09', 'HLA-C*03:10', 'HLA-C*03:11', 'HLA-C*03:12', 'HLA-C*03:13',
         'HLA-C*03:14', 'HLA-C*03:15', 'HLA-C*03:16', 'HLA-C*03:17', 'HLA-C*03:18', 'HLA-C*03:19', 'HLA-C*03:21', 'HLA-C*03:23', 'HLA-C*03:24',
         'HLA-C*03:25', 'HLA-C*03:26', 'HLA-C*03:27', 'HLA-C*03:28', 'HLA-C*03:29', 'HLA-C*03:30', 'HLA-C*03:31', 'HLA-C*03:32', 'HLA-C*03:33',
         'HLA-C*03:34', 'HLA-C*03:35', 'HLA-C*03:36', 'HLA-C*03:37', 'HLA-C*03:38', 'HLA-C*03:39', 'HLA-C*03:40', 'HLA-C*03:41', 'HLA-C*03:42',
         'HLA-C*03:43', 'HLA-C*03:44', 'HLA-C*03:45', 'HLA-C*03:46', 'HLA-C*03:47', 'HLA-C*03:48', 'HLA-C*03:49', 'HLA-C*03:50', 'HLA-C*03:51',
         'HLA-C*03:52', 'HLA-C*03:53', 'HLA-C*03:54', 'HLA-C*03:55', 'HLA-C*03:56', 'HLA-C*03:57', 'HLA-C*03:58', 'HLA-C*03:59', 'HLA-C*03:60',
         'HLA-C*03:61', 'HLA-C*03:62', 'HLA-C*03:63', 'HLA-C*03:64', 'HLA-C*03:65', 'HLA-C*03:66', 'HLA-C*03:67', 'HLA-C*03:68', 'HLA-C*03:69',
         'HLA-C*03:70', 'HLA-C*03:71', 'HLA-C*03:72', 'HLA-C*03:73', 'HLA-C*03:74', 'HLA-C*03:75', 'HLA-C*03:76', 'HLA-C*03:77', 'HLA-C*03:78',
         'HLA-C*03:79', 'HLA-C*03:80', 'HLA-C*03:81', 'HLA-C*03:82', 'HLA-C*03:83', 'HLA-C*03:84', 'HLA-C*03:85', 'HLA-C*03:86', 'HLA-C*03:87',
         'HLA-C*03:88', 'HLA-C*03:89', 'HLA-C*03:90', 'HLA-C*03:91', 'HLA-C*03:92', 'HLA-C*03:93', 'HLA-C*03:94', 'HLA-C*04:01', 'HLA-C*04:03',
         'HLA-C*04:04', 'HLA-C*04:05', 'HLA-C*04:06', 'HLA-C*04:07', 'HLA-C*04:08', 'HLA-C*04:10', 'HLA-C*04:11', 'HLA-C*04:12', 'HLA-C*04:13',
         'HLA-C*04:14', 'HLA-C*04:15', 'HLA-C*04:16', 'HLA-C*04:17', 'HLA-C*04:18', 'HLA-C*04:19', 'HLA-C*04:20', 'HLA-C*04:23', 'HLA-C*04:24',
         'HLA-C*04:25', 'HLA-C*04:26', 'HLA-C*04:27', 'HLA-C*04:28', 'HLA-C*04:29', 'HLA-C*04:30', 'HLA-C*04:31', 'HLA-C*04:32', 'HLA-C*04:33',
         'HLA-C*04:34', 'HLA-C*04:35', 'HLA-C*04:36', 'HLA-C*04:37', 'HLA-C*04:38', 'HLA-C*04:39', 'HLA-C*04:40', 'HLA-C*04:41', 'HLA-C*04:42',
         'HLA-C*04:43', 'HLA-C*04:44', 'HLA-C*04:45', 'HLA-C*04:46', 'HLA-C*04:47', 'HLA-C*04:48', 'HLA-C*04:49', 'HLA-C*04:50', 'HLA-C*04:51',
         'HLA-C*04:52', 'HLA-C*04:53', 'HLA-C*04:54', 'HLA-C*04:55', 'HLA-C*04:56', 'HLA-C*04:57', 'HLA-C*04:58', 'HLA-C*04:60', 'HLA-C*04:61',
         'HLA-C*04:62', 'HLA-C*04:63', 'HLA-C*04:64', 'HLA-C*04:65', 'HLA-C*04:66', 'HLA-C*04:67', 'HLA-C*04:68', 'HLA-C*04:69', 'HLA-C*04:70',
         'HLA-C*05:01', 'HLA-C*05:03', 'HLA-C*05:04', 'HLA-C*05:05', 'HLA-C*05:06', 'HLA-C*05:08', 'HLA-C*05:09', 'HLA-C*05:10', 'HLA-C*05:11',
         'HLA-C*05:12', 'HLA-C*05:13', 'HLA-C*05:14', 'HLA-C*05:15', 'HLA-C*05:16', 'HLA-C*05:17', 'HLA-C*05:18', 'HLA-C*05:19', 'HLA-C*05:20',
         'HLA-C*05:21', 'HLA-C*05:22', 'HLA-C*05:23', 'HLA-C*05:24', 'HLA-C*05:25', 'HLA-C*05:26', 'HLA-C*05:27', 'HLA-C*05:28', 'HLA-C*05:29',
         'HLA-C*05:30', 'HLA-C*05:31', 'HLA-C*05:32', 'HLA-C*05:33', 'HLA-C*05:34', 'HLA-C*05:35', 'HLA-C*05:36', 'HLA-C*05:37', 'HLA-C*05:38',
         'HLA-C*05:39', 'HLA-C*05:40', 'HLA-C*05:41', 'HLA-C*05:42', 'HLA-C*05:43', 'HLA-C*05:44', 'HLA-C*05:45', 'HLA-C*06:02', 'HLA-C*06:03',
         'HLA-C*06:04', 'HLA-C*06:05', 'HLA-C*06:06', 'HLA-C*06:07', 'HLA-C*06:08', 'HLA-C*06:09', 'HLA-C*06:10', 'HLA-C*06:11', 'HLA-C*06:12',
         'HLA-C*06:13', 'HLA-C*06:14', 'HLA-C*06:15', 'HLA-C*06:17', 'HLA-C*06:18', 'HLA-C*06:19', 'HLA-C*06:20', 'HLA-C*06:21', 'HLA-C*06:22',
         'HLA-C*06:23', 'HLA-C*06:24', 'HLA-C*06:25', 'HLA-C*06:26', 'HLA-C*06:27', 'HLA-C*06:28', 'HLA-C*06:29', 'HLA-C*06:30', 'HLA-C*06:31',
         'HLA-C*06:32', 'HLA-C*06:33', 'HLA-C*06:34', 'HLA-C*06:35', 'HLA-C*06:36', 'HLA-C*06:37', 'HLA-C*06:38', 'HLA-C*06:39', 'HLA-C*06:40',
         'HLA-C*06:41', 'HLA-C*06:42', 'HLA-C*06:43', 'HLA-C*06:44', 'HLA-C*06:45', 'HLA-C*07:01', 'HLA-C*07:02', 'HLA-C*07:03', 'HLA-C*07:04',
         'HLA-C*07:05', 'HLA-C*07:06', 'HLA-C*07:07', 'HLA-C*07:08', 'HLA-C*07:09', 'HLA-C*07:10', 'HLA-C*07:100', 'HLA-C*07:101', 'HLA-C*07:102',
         'HLA-C*07:103', 'HLA-C*07:105', 'HLA-C*07:106', 'HLA-C*07:107', 'HLA-C*07:108', 'HLA-C*07:109', 'HLA-C*07:11', 'HLA-C*07:110',
         'HLA-C*07:111', 'HLA-C*07:112', 'HLA-C*07:113', 'HLA-C*07:114', 'HLA-C*07:115', 'HLA-C*07:116', 'HLA-C*07:117', 'HLA-C*07:118',
         'HLA-C*07:119', 'HLA-C*07:12', 'HLA-C*07:120', 'HLA-C*07:122', 'HLA-C*07:123', 'HLA-C*07:124', 'HLA-C*07:125', 'HLA-C*07:126',
         'HLA-C*07:127', 'HLA-C*07:128', 'HLA-C*07:129', 'HLA-C*07:13', 'HLA-C*07:130', 'HLA-C*07:131', 'HLA-C*07:132', 'HLA-C*07:133',
         'HLA-C*07:134', 'HLA-C*07:135', 'HLA-C*07:136', 'HLA-C*07:137', 'HLA-C*07:138', 'HLA-C*07:139', 'HLA-C*07:14', 'HLA-C*07:140',
         'HLA-C*07:141', 'HLA-C*07:142', 'HLA-C*07:143', 'HLA-C*07:144', 'HLA-C*07:145', 'HLA-C*07:146', 'HLA-C*07:147', 'HLA-C*07:148',
         'HLA-C*07:149', 'HLA-C*07:15', 'HLA-C*07:16', 'HLA-C*07:17', 'HLA-C*07:18', 'HLA-C*07:19', 'HLA-C*07:20', 'HLA-C*07:21', 'HLA-C*07:22',
         'HLA-C*07:23', 'HLA-C*07:24', 'HLA-C*07:25', 'HLA-C*07:26', 'HLA-C*07:27', 'HLA-C*07:28', 'HLA-C*07:29', 'HLA-C*07:30', 'HLA-C*07:31',
         'HLA-C*07:35', 'HLA-C*07:36', 'HLA-C*07:37', 'HLA-C*07:38', 'HLA-C*07:39', 'HLA-C*07:40', 'HLA-C*07:41', 'HLA-C*07:42', 'HLA-C*07:43',
         'HLA-C*07:44', 'HLA-C*07:45', 'HLA-C*07:46', 'HLA-C*07:47', 'HLA-C*07:48', 'HLA-C*07:49', 'HLA-C*07:50', 'HLA-C*07:51', 'HLA-C*07:52',
         'HLA-C*07:53', 'HLA-C*07:54', 'HLA-C*07:56', 'HLA-C*07:57', 'HLA-C*07:58', 'HLA-C*07:59', 'HLA-C*07:60', 'HLA-C*07:62', 'HLA-C*07:63',
         'HLA-C*07:64', 'HLA-C*07:65', 'HLA-C*07:66', 'HLA-C*07:67', 'HLA-C*07:68', 'HLA-C*07:69', 'HLA-C*07:70', 'HLA-C*07:71', 'HLA-C*07:72',
         'HLA-C*07:73', 'HLA-C*07:74', 'HLA-C*07:75', 'HLA-C*07:76', 'HLA-C*07:77', 'HLA-C*07:78', 'HLA-C*07:79', 'HLA-C*07:80', 'HLA-C*07:81',
         'HLA-C*07:82', 'HLA-C*07:83', 'HLA-C*07:84', 'HLA-C*07:85', 'HLA-C*07:86', 'HLA-C*07:87', 'HLA-C*07:88', 'HLA-C*07:89', 'HLA-C*07:90',
         'HLA-C*07:91', 'HLA-C*07:92', 'HLA-C*07:93', 'HLA-C*07:94', 'HLA-C*07:95', 'HLA-C*07:96', 'HLA-C*07:97', 'HLA-C*07:99', 'HLA-C*08:01',
         'HLA-C*08:02', 'HLA-C*08:03', 'HLA-C*08:04', 'HLA-C*08:05', 'HLA-C*08:06', 'HLA-C*08:07', 'HLA-C*08:08', 'HLA-C*08:09', 'HLA-C*08:10',
         'HLA-C*08:11', 'HLA-C*08:12', 'HLA-C*08:13', 'HLA-C*08:14', 'HLA-C*08:15', 'HLA-C*08:16', 'HLA-C*08:17', 'HLA-C*08:18', 'HLA-C*08:19',
         'HLA-C*08:20', 'HLA-C*08:21', 'HLA-C*08:22', 'HLA-C*08:23', 'HLA-C*08:24', 'HLA-C*08:25', 'HLA-C*08:27', 'HLA-C*08:28', 'HLA-C*08:29',
         'HLA-C*08:30', 'HLA-C*08:31', 'HLA-C*08:32', 'HLA-C*08:33', 'HLA-C*08:34', 'HLA-C*08:35', 'HLA-C*12:02', 'HLA-C*12:03', 'HLA-C*12:04',
         'HLA-C*12:05', 'HLA-C*12:06', 'HLA-C*12:07', 'HLA-C*12:08', 'HLA-C*12:09', 'HLA-C*12:10', 'HLA-C*12:11', 'HLA-C*12:12', 'HLA-C*12:13',
         'HLA-C*12:14', 'HLA-C*12:15', 'HLA-C*12:16', 'HLA-C*12:17', 'HLA-C*12:18', 'HLA-C*12:19', 'HLA-C*12:20', 'HLA-C*12:21', 'HLA-C*12:22',
         'HLA-C*12:23', 'HLA-C*12:24', 'HLA-C*12:25', 'HLA-C*12:26', 'HLA-C*12:27', 'HLA-C*12:28', 'HLA-C*12:29', 'HLA-C*12:30', 'HLA-C*12:31',
         'HLA-C*12:32', 'HLA-C*12:33', 'HLA-C*12:34', 'HLA-C*12:35', 'HLA-C*12:36', 'HLA-C*12:37', 'HLA-C*12:38', 'HLA-C*12:40', 'HLA-C*12:41',
         'HLA-C*12:43', 'HLA-C*12:44', 'HLA-C*14:02', 'HLA-C*14:03', 'HLA-C*14:04', 'HLA-C*14:05', 'HLA-C*14:06', 'HLA-C*14:08', 'HLA-C*14:09',
         'HLA-C*14:10', 'HLA-C*14:11', 'HLA-C*14:12', 'HLA-C*14:13', 'HLA-C*14:14', 'HLA-C*14:15', 'HLA-C*14:16', 'HLA-C*14:17', 'HLA-C*14:18',
         'HLA-C*14:19', 'HLA-C*14:20', 'HLA-C*15:02', 'HLA-C*15:03', 'HLA-C*15:04', 'HLA-C*15:05', 'HLA-C*15:06', 'HLA-C*15:07', 'HLA-C*15:08',
         'HLA-C*15:09', 'HLA-C*15:10', 'HLA-C*15:11', 'HLA-C*15:12', 'HLA-C*15:13', 'HLA-C*15:15', 'HLA-C*15:16', 'HLA-C*15:17', 'HLA-C*15:18',
         'HLA-C*15:19', 'HLA-C*15:20', 'HLA-C*15:21', 'HLA-C*15:22', 'HLA-C*15:23', 'HLA-C*15:24', 'HLA-C*15:25', 'HLA-C*15:26', 'HLA-C*15:27',
         'HLA-C*15:28', 'HLA-C*15:29', 'HLA-C*15:30', 'HLA-C*15:31', 'HLA-C*15:33', 'HLA-C*15:34', 'HLA-C*15:35', 'HLA-C*16:01', 'HLA-C*16:02',
         'HLA-C*16:04', 'HLA-C*16:06', 'HLA-C*16:07', 'HLA-C*16:08', 'HLA-C*16:09', 'HLA-C*16:10', 'HLA-C*16:11', 'HLA-C*16:12', 'HLA-C*16:13',
         'HLA-C*16:14', 'HLA-C*16:15', 'HLA-C*16:17', 'HLA-C*16:18', 'HLA-C*16:19', 'HLA-C*16:20', 'HLA-C*16:21', 'HLA-C*16:22', 'HLA-C*16:23',
         'HLA-C*16:24', 'HLA-C*16:25', 'HLA-C*16:26', 'HLA-C*17:01', 'HLA-C*17:02', 'HLA-C*17:03', 'HLA-C*17:04', 'HLA-C*17:05', 'HLA-C*17:06',
         'HLA-C*17:07', 'HLA-C*18:01', 'HLA-C*18:02', 'HLA-C*18:03', 'HLA-E*01:01', 'HLA-G*01:01', 'HLA-G*01:02', 'HLA-G*01:03', 'HLA-G*01:04',
         'HLA-G*01:06', 'HLA-G*01:07', 'HLA-G*01:08', 'HLA-G*01:09',
         'H2-Db', 'H2-Dd', 'H2-Kb', 'H2-Kd', 'H2-Kk', 'H2-Ld'])
    __version = "1.1"

    @property
    def version(self):
        """The version of the predictor"""
        return self.__version

    def _represent(self, allele):
        """
        Internal function transforming an allele object into its representative string
        :param allele: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: str
        """
        if isinstance(allele, MouseAllele):
            return "H-2-%s%s%s" % (allele.locus, allele.supertype, allele.subtype)
        else:
            return "HLA-%s%s:%s" % (allele.locus, allele.supertype, allele.subtype)

    def convert_alleles(self, alleles):
        """
        Converts :class:`~Fred2.Core.Allele.Allele` into the internal :class:`~Fred2.Core.Allele.Allele` representation
        of the predictor and returns a string representation

        :param alleles: The :class:`~Fred2.Core.Allele.Allele` for which the internal predictor representation is
                        needed
        :type alleles: :class:`~Fred2.Core.Allele.Allele`
        :return: Returns a string representation of the input :class:`~Fred2.Core.Allele.Allele`
        :rtype: list(str)
        """
        return [self._represent(a) for a in alleles]

    @property
    def supportedAlleles(self):
        """
        A list of supported :class:`~Fred2.Core.Allele.Allele`
        """
        return self.__alleles

    @property
    def name(self):
        """The name of the predictor"""
        return self.__name

    @property
    def command(self):
        """
        Defines the commandline call for external tool
        """
        return self.__command

    @property
    def supportedLength(self):
        """
        A list of supported :class:`~Fred2.Core.Peptide.Peptide` lengths
        """
        return self.__supported_length

    def parse_external_result(self, file):
        """
        Parses external results and returns the result

        :param str file: The file path or the external prediction results
        :return: A dictionary containing the prediction results
        :rtype: dict
        """
        result = defaultdict(defaultdict)
        with open(file, "r") as f:
            for l in f:
                if l.startswith("#") or l.startswith("-") or l.strip() == "":
                    continue
                row = l.strip().split()
                if not row[0].isdigit():
                    continue

                epitope, allele, comb_score = row[3], row[2], row[7]
                result[allele.replace("*", "")][epitope] = float(comb_score)
        return result

    def get_external_version(self, path=None):
        """
        Returns the external version of the tool by executing
        >{command} --version

        might be dependent on the method and has to be overwritten
        therefore it is declared abstract to enforce the user to
        overwrite the method. The function in the base class can be called
        with super()

        :param str path: Optional specification of executable path if deviant from :attr:`self.__command`
        :return: The external version of the tool or None if tool does not support versioning
        :rtype: str
        """
        return None

    def prepare_input(self, input, file):
        """
        Prepares input for external tools and writes them to file in the specific format

        No return value!

        :param: list(str) input: The :class:`~Fred2.Core.Peptide.Peptide` sequences to write into file
        :param File file: File-handler to input file for external tool
        """
        file.write("\n".join(">pepe_%i\n%s" % (i, p) for i, p in enumerate(input)))