__author__ = 'schubert'

import argparse, sys, os
from Fred2.IO import FileReader
from Fred2.Core.Generator import generate_peptides_from_protein
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import AEpitopePrediction
from Fred2.EpitopePrediction import EpitopePredictorFactory


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
        proteins = FileReader.read_fasta(args.input, type="Protein")
        peptides = generate_peptides_from_protein(proteins, args.length)

    elif args.pepfile:
        peptides = FileReader.read_lines(args.input, type="Peptide")

    elif args.peptide:
        peptides = [Peptide(s) for s in args.input]

    #read in alleles
    if args.allelefile:
        alleles = FileReader.read_lines(args.alleles, type="Allele")
    else:
        alleles = [Allele(a.upper()) for a in args.alleles]

    result = [EpitopePredictorFactory(m).predict(peptides, alleles) for m in args.method]
    r_df = result.pop()
    for r in result:
        r_df_a, r_a = r_df.align(r, fill_value=0)
        r_df = r_df_a + r_a

    output = args.output if args.outdir == "" else args.outdir + os.path.basename(args.output)
    with open(output, "w") as out:
        r_df.to_csv(out)



    #generate Galaxy HTML output
    if args.html:
        begin_html = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <link rel="stylesheet" href="/static/style/blue/etk.css" type="text/css" />
    <script type="text/javascript" src="/static/scripts/packed/libs/jquery/jquery.js"></script>
    <script type="text/javascript" src="/static/scripts/packed/libs/jquery/jquery.tablesorter.js"></script>
    <script type="text/javascript" src="/static/scripts/libs/etk.js"></script>
</head>

<body>
    <div class="document">"""

        setting = """  <h2 class="etk-heading">Epitope Prediction Results</h2>

        <table class="etk-parameterT">
            <tr> <th class ="etk-innerHeading" colspan="2"> Parameters </th></tr>
            <tr>
                <th>Prediction Method:</th>
                <td>%s</td>
            </tr>
        </table>"""%args.method



        table="""

        <input id="etk-search" placeholder="  filter">
        <table class="etk-sortT etk-resultsT etk-filterT">

            <thead>
                <tr>
                    <th>Peptide</th>"""+"".join("<th>%s</th>"%str(a) for a in result.columns) \
            +"""
                </tr>
            </thead>"""+"".join("<tr><td>%s<td>%s</tr>"%(r[0] ,"".join("<td align='right'>%s</td>"%str(result.loc[r, c])))
                                for r in result.index for c in result.columns)+"</table>"

        end_html = "</div></body></html>"

        html_out = ".".join(output.split(".")[:-1])+".html"
        with open(html_out, "w") as html_o:
            html_o.write(begin_html+setting+table+end_html)

if __name__ == "__main__":
    main()
