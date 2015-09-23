# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: HLAtyping.External
   :synopsis: This module contains all classes for HLA typing with external methods.
.. moduleauthor:: schubert

"""

import os.path
import subprocess
import csv
import abc
import shutil
import warnings

from tempfile import NamedTemporaryFile

from Fred2.Core import AHLATyping, AExternal
from Fred2.Core import Allele


class AExternalHLATyping(AHLATyping, AExternal):

    def predict(self, ngsFile, output, command=None, options=None, delete=True, **kwargs):
        """
        Implementation of prediction

        :param str ngsFile: The path to the NGS file of interest
        :param str output: The path to the output file or directory
        :param str command: The path to a alternative binary (if binary is not globally executable)
        :param str options: A string with additional options that is directly past to the tool
        :param bool delete: Boolean indicator whether generated files should be deleted afterwards
        :return: list(Allele) - A list of Allele objects representing the most likely HLA genotype
        """

        if not self.is_in_path() and "path" not in kwargs:
            raise RuntimeError("{name} {version} could not be found in PATH".format(name=self.name,
                                                                                    version=self.version))
        external_version = self.get_external_version(command)
        if self.version != external_version and external_version is not None:
            raise RuntimeError("Internal version {internal_version} does "
                               "not match external version {external_version}".format(internal_version=self.version,
                                                                                      external_version=external_version))

        if not os.path.exists(ngsFile):
            raise ValueError("Specified file or directory {fil} does not exist".format(fil=ngsFile))

        #allowe customary executable specification
        if command is not None:
            exe = self.command.split()[0]
            _command = self.command.replace(exe, command)
        else:
            _command = self.command

        if output is None:
            tmp_output = NamedTemporaryFile(delete=False)
            output = tmp_output.name
            tmp_output.close()

        try:
            stdo = None
            stde = None
            cmd = _command.format(file=ngsFile, options="" if options is None else options, out=output)
            print cmd
            p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait() #block the rest
            stdo, stde = p.communicate()
            stdr = p.returncode
            print stdo
            print stde
            if stdr > 0:
                raise RuntimeError("Unsuccessful execution of " + cmd + " (EXIT!=0) with error: " + stde)
        except Exception as e:
            raise RuntimeError(e)

        genotype = self.parse_external_result(output)

        if delete:
            self.clean_up(output)

        return genotype

    @abc.abstractmethod
    def clean_up(self, _output):
        """
        Cleans the generated files after prediction

        :param str output: The path to the output file or directory
        """
        raise NotImplementedError


class OptiType_1_0(AExternalHLATyping):
    """
    Wrapper of OptiType v1.0

    Szolek, A., Schubert, B., Mohr, C., Sturm, M., Feldhahn, M., & Kohlbacher, O. (2014).
    OptiType: precision HLA typing from next-generation sequencing data. Bioinformatics, 30(23), 3310-3316.
    """
    __name = "OptiType"
    __version = "1.0"
    __command = "OptiTypePipeline.py -i {file} {options} -o {out}"

    @property
    def version(self):
        return self.__version

    @property
    def name(self):
        return self.__name

    @property
    def command(self):
        return self.__command

    def get_external_version(self, path=None):
        return None

    def parse_external_result(self, _output):
        """
        Searches within the defined dir _file for the newest dir and reads
        the prediction file from there

        :param str _output: the path to the output dir
        :return: list(Allele) - The predicted HLA genotype
        """
        all_subdirs = [os.path.join(_output,d) for d in os.listdir(_output) if os.path.isdir(os.path.join(_output,d))]
        latest_subdir = max(all_subdirs, key=os.path.getmtime)
        result_file = latest_subdir+"/"+os.path.basename(os.path.normpath(latest_subdir))+"_result.tsv"
        with open(result_file, "r") as f:
            row = csv.DictReader(f, delimiter="\t").next()
            return map(lambda x: Allele("HLA-"+x), [ row[k] for k in ["A1","A2","B1","B2","C1","C2"]])

    def clean_up(self, _output):
        """
        Searches within the defined dir _file for the newest dir and deletes it.
        This should be the one OptiType had created

        This could cause some terrible site effects if someone or something also writes in that directory!!
        OptiType should change the way it writes its output!

        :param str _output: the path to the output file or directory of the programme
        """
        all_subdirs = [os.path.join(_output, d) for d in os.listdir(_output) if os.path.isdir(os.path.join(_output,d))]
        latest_subdir = max(all_subdirs, key=os.path.getmtime)
        shutil.rmtree(latest_subdir)


class Seq2HLA_2_2(AExternalHLATyping):
    """
    Wrapper of seq2HLA v2.2


    Boegel, S., Scholtalbers, J., Loewer, M., Sahin, U., & Castle, J. C. (2015).
    In Silico HLA Typing Using Standard RNA-Seq Sequence Reads. Molecular Typing of Blood Cell Antigens, 247.
    """
    __name = "seq2HLA"
    __version = "2.2"
    __command = "seq2HLA.py -1 {file} {options} -r {out}"

    @property
    def version(self):
        return self.__version

    @property
    def name(self):
        return self.__name

    @property
    def command(self):
        return self.__command

    def get_external_version(self, path=None):
        try:
            stdo = None
            stde = None
            cmd = self.command.split()[0]+" --version"
            p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait() #block the rest
            stdo, stde = p.communicate()
            stdr = p.returncode
            if stdr > 0:
                raise RuntimeError("Unsuccessful execution of " + cmd + " (EXIT!=0) with error: " + stde)
            _, version = stdo.split()
            return version
        except Exception as e:
            raise RuntimeError(e)

    def predict(self, ngsFile, output, command=None, options=None, delete=True, **kwargs):
        if "-2" not in options:
            ValueError("Seq2HLA only supports paired-end inputs. Please use the options "
                       "parameter to specify the second ngs file (e.g. options='-2 /path/to/ngs2.fq'")
        return super(Seq2HLA_2_2, self).predict(ngsFile, output, command=command,
                                                options=options, delete=delete, **kwargs)

    def parse_external_result(self, _file):
        alleles = []
        try:
            with open(_file+"-ClassI.HLAgenotype4digits") as c1:
                for row in csv.DictReader(c1, delimiter="\t"):
                    alleles.extend([Allele("HLA-"+row["Allele 1"]), Allele("HLA-"+row["Allele 2"])])
        except IOError as e:
            warnings.warn("Output file {c1} for HLA-I could not be found. {error}".format(
                c1=_file + "-ClassI.HLAgenotype4digits"), error=e)

        try:
            with open(_file+"-ClassII.HLAgenotype4digits") as c2:
                for row in csv.DictReader(c2, delimiter="\t"):
                    alleles.extend([Allele("HLA-"+row["Allele 1"]), Allele("HLA-"+row["Allele 2"])])
        except IOError as e:
            warnings.warn("Output file {c2} for HLA-I could not be found. {error}".format(
                c2=_file + "-ClassII.HLAgenotype4digits"), error=e)

        return alleles

    def clean_up(self, _output):
        if os.path.isdir(_output):
            #if _output was mistakenly set to a directory all seq2HLA files will start with -ClassI or -ClassII
            for f in os.listdir(_output):
                if f.startswith("-ClassI") or f.startswith("-ClassII"):
                    os.remove(os.path.join(_output, f))
        else:
            basedir = os.path.dirname(_output)
            prefix = os.path.basename(_output)
            for f in os.listdir(basedir):
                if f.startswith(prefix):
                    os.remove(os.path.join(basedir, f))


class ATHLATES_1_0(AExternalHLATyping):
    """
    Wrapper for ATHLATES

    C. Liu, X. Yang, B. Duffy, T. Mohanakumar, R.D. Mitra, M.C. Zody, J.D. Pfeifer (2012) ATHLATES:
    accurate typing of human leukocyte antigen through exome sequencing, Nucl. Acids Res.
    (2013)
    """

    __name = "athlates"
    __version = "1.0"
    __command = "typing -bam {file} {options} -o {out}"

    @property
    def version(self):
        return self.__version

    @property
    def name(self):
        return self.__name

    @property
    def command(self):
        return self.__command

    def get_external_version(self, path=None):
        return None

    def parse_external_result(self, _output):
        """
        Searches within the defined dir _file for the newest dir and reads
        the prediction file from there

        :param str _output: the path to the output dir
        :return: list(Allele) - The predicted HLA genotype
        """
        alleles = []
        if os.path.isdir(_output):
            _file = os.path.join(_output, ".typing.txt")
        else:
            _file = _output+".typing.txt"

        typing = False
        with open(_file, "r") as f:

            for l in f:
                if typing and l.strip() != "":
                    a1, a2, _ = l.split()
                    alleles.append(Allele("HLA-"+":".join(a1.split(":")[:2])))
                    alleles.append(Allele("HLA-"+":".join(a2.split(":")[:2])))
                if "------------ Inferred Allelic Pairs -------------" in l:
                    typing = True
        return alleles

    def clean_up(self, _output):
        """
        deletes files created by ATHLATES within _output

        :param str _output: the path to the output file or directory of the programme
        """
        if os.path.isdir(_output):
            for f in os.listdir(_output):
                if f.startswith("."):
                    os.remove(os.path.join(_output, f))
        else:
            basedir = os.path.dirname(_output)
            prefix = os.path.basename(_output)
            for f in os.listdir(basedir):
                if f.startswith(prefix):
                    os.remove(os.path.join(basedir, f))


class Polysolver(AExternalHLATyping):
    """
    Wrapper for Polysolver


    Shukla, Sachet A., Rooney, Michael S., Rajasagi, Mohini, Tiao, Grace, et al. (2015).
    Comprehensive analysis of cancer-associated somatic mutations in class I HLA genes.
    Nat Biotech, advance online publication. doi: 10.1038/nbt.3344
    """

    __name = "polysolver"
    __version = "1.0"
    __command = "shell_call_hla_type {file} {options} {out}"

    @property
    def version(self):
        return self.__version

    @property
    def name(self):
        return self.__name

    @property
    def command(self):
        return self.__command

    def get_external_version(self, path=None):
        return None

    def parse_external_result(self, _output):
        """
        Searches within the defined dir _file for the newest dir and reads
        the prediction file from there

        :param str _output: the path to the output dir
        :return: list(Allele) - The predicted HLA genotype
        """
        alleles = []
        try:
            with open(os.path.join(_output, "winner.hla.txt"), "r") as f:
                for l in f:
                    try:
                        _, a1, a2 = l.replace("-n", "").replace("-e", "").strip().split()
                        a1 = a1.split("_")
                        a2 = a2.split("_")
                        alleles.extend([Allele("HLA-"+a1[1].upper()+"*"+a1[2]+":"+a1[3]),
                                        Allele("HLA-"+a2[1].upper()+"*"+a2[2]+":"+a2[3])])
                    except ValueError:
                        IOError(
                            "Output format seems incorrect:\n{line}\n. Please check if Polysolver ran correctly.".format(
                                lines=l))
                return alleles
        except IOError:
            raise IOError("File {out} could not be found. Please check your specified output folder".format(
                out=os.path.join(_output, "winner.hla.txt")))

    def clean_up(self, _output):
        """
        deletes files created by Polysolver within _output

        :param str _output: the path to the output file or directory of the programme
        """
        os.remove(os.path.join(_output, "winner.hla.txt"))