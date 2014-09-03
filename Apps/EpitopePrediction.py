__author__ = 'schubert'

import argparse, sys
from Fred2.IO import FileReader
from Fred2.Core.Generator import generate_peptides_from_protein
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import AEpitopePrediction
from Fred2.EpitopePrediction import EpitopePredictorFactory

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reads protein or peptide sequences and predicts peptides "+
                                                 "for a specified prediction method and HLA alleles.")
    parser.add_argument("-i", "--input",
                        nargs="+",
                        requried=True,
                        help="Input data can be RefSeq ID, UniProt ID, fasta file, peptide file (one peptide per line),"
                             +" or peptide sequences as sequences (max 50)"
                        )
    input_types = argparse.add_mutually_exclusive_group(required=True)
    input_types.add_argument("-r","--refseq",
                             action="store_true",
                             help= "Specifies the input as RefSeq IDs")
    input_types.add_argument("-u","--uniprot",
                             action="store_true",
                             help= "Specifies the input as UniProt IDs")
    input_types.add_argument("-f","--fasta",
                             action="store_true",
                             help= "Specifies the input as protein (multi-)Fasta file")
    input_types.add_argument("-pf","--pepfile",
                             action="store_true",
                             help= "Specifies the input as peptide file")
    input_types.add_argument("-p","--peptide",
                             action="store_true",
                             help= "Specifies the input as peptide sequences")
    parser.add_argument("-a", "--alleles",
                        nargs="+",
                        required=True,
                        help="Specifies for which alleles prediction should be made. " +
                             "Input either can be alleles as string (new nomenclature), or a file with one allele per line.")
    allele_types = argparse.add_mutually_exclusive_group(required=True)
    allele_types.add_arguments("-af", "--allelefile",
                               action="store_true",
                               help="Specifies the allele input as allele file.")
    allele_types.add_arguments("-as", "--allelestring",
                               action="store_true",
                               help="Specifies the allele input as allele string.")
    parser.add_argument("-m", "--method",
                       requried=True,
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
    args = parser.parse_args()

    if args.available:
        for pred, obj in AEpitopePrediction.registry.iteritems():
            if pred != "AEpitopePrediction":
                print "Method: ",pred
                print "Supported Alleles: ", " ".join(pred.supportedAlleles)
                print "Supported Length: ", " ".join(map(str, pred.supportedLength))
        sys.exit(0)


    '''
    Parser Input
    '''
    #RefSeq
    if args.refseq:
        pass

    #UniProt
    elif args.uniprot:
        pass

    #fasta protein
    elif args.fasta:
        proteins = FileReader.get_sequence(args.input, type="Protein")
        peptides = generate_peptides_from_protein(proteins, args.length)

    elif args.pepfile:
        peptides = FileReader.read_line(args.input, type="Peptide")

    elif args.peptide:
        peptides = [Peptide(s) for s in args.input]

    #read in alleles
    if args.allelefile:
        alleles = FileReader.read_line(args.input, type="Allele")
    else:
        alleles = [Allele(a) for a in args.allelestring]

    result = EpitopePredictorFactory(args.method).predict(peptides, alleles)

    with open(args.output) as out:
        result.to_csv(out)