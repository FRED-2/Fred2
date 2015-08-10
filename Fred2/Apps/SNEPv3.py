__author__ = 'schubert'
"""
Generates Polymorphic epitopes for human proteins

input is either a list of gene identifier or a vcf file

"""
import argparse, sys, os
from Fred2.IO import FileReader
from Fred2.Core.Generator import generate_peptides_from_protein
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import AEpitopePrediction
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.IO.MartsAdapter import MartsAdapter


def main():
    parser = argparse.ArgumentParser(description="Reads protein or peptide sequences and predicts peptides "+
                                                 "for a specified prediction method and HLA alleles.")
    parser.add_argument("-i", "--input",
                        nargs="+",
                        required=True,
                        help="Input data can be RefSeq ID, UniProt ID, fasta file, peptide file (one peptide per line),"
                             +" or peptide sequences as sequences (max 50)"
                        )
    input_types = parser.add_mutually_exclusive_group(required=True)
    input_types.add_argument("-p","--protein",
                             action="store_true",
                             help= "Specifies if IDs are protein IDs")
    input_types.add_argument("-r","--rna",
                             action="store_true",
                             help= "Specifies the input as rna IDs")
    input_types.add_argument("-f","--fasta",
                             action="store_true",
                             help= "Specifies the input as protein (multi-)Fasta file")
    parser.add_argument("-a", "--alleles",
                        nargs="+",
                        required=True,
                        help="Specifies for which alleles prediction should be made. " +
                             "Input either can be alleles as string (new nomenclature), or a file with one allele per line.")
    allele_types = parser.add_mutually_exclusive_group(required=True)
    allele_types.add_argument("-af", "--allelefile",
                               action="store_true",
                               help="Specifies the allele input as allele file.")
    allele_types.add_argument("-as", "--allelestring",
                               action="store_true",
                               help="Specifies the allele input as allele string.")
    parser.add_argument("-m", "--method",
                       required=True,
                       nargs="+",
                       help="Specifies the method used for prediction.")
    parser.add_argument("-l", "--length",
                        required=False,
                        type=int,
                        default=9,
                        help="Specifies the length of the peptides (default=9).")
    parser.add_argument("-o", "--output",
                        required=True,
                        help="Specifies the output path. Results will be written to CSV")
    parser.add_argument("-am", "--available",
                        required=False,
                        action="store_true",
                        help="Returns all available methods and their allele models.")

    #COMMENT: These options are hidden and only used for ETK2
    parser.add_argument("-html", "--html",
                        required=False,
                        action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("-od", "--outdir",
                        required=False,
                        default="",
                        help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.available:
        for pred, obj in AEpitopePrediction.registry.iteritems():
            if pred not in ["AEpitopePrediction", "APSSMEpitopePredictor", "ANetMHC", "ASVMEpitopePrediction"]:
                print "Method: ",pred
                print "Supported Alleles: ", " ".join(getattr(obj, "_"+pred+"__alleles" ))
                print "Supported Length: ", " ".join(map(str, getattr(obj,  "_"+pred+"__supported_length")))
                print
        sys.exit(0)


    mart = MartsAdapter()
    #RefSeq
    transcripts = []
    if args.protein:
        for pid in args.input:
            ids = mart.get_transcript_information_from_protein_id(pid)

    #UniProt
    elif args.uniprot:
        pass