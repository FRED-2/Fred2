__author__ = 'schubert'


import argparse, sys, os

from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import ACleavageSitePrediction
from Fred2.IO.FileReader import read_lines
from Fred2.EpitopeAssembly.EpitopeAssembly import EpitopeAssembly
from Fred2.CleavagePrediction import CleavageSitePredictorFactory

def main():
    parser = argparse.ArgumentParser(description="Reads peptide sequences and predicts best arrangement "+
                                                 "for a string-of-beats peptide vaccine based on proteasomal cleavage prediction.")
    parser.add_argument("-i", "--input",
                        nargs="+",
                        required=True,
                        help="peptide file (one peptide per line),"
                             +" or peptide sequences as sequences (max 50)"
                        )
    input_types = parser.add_mutually_exclusive_group(required=True)
    input_types.add_argument("-pf","--pepfile",
                             action="store_true",
                             help= "Specifies the input as peptide file")
    input_types.add_argument("-p","--peptide",
                             action="store_true",
                             help= "Specifies the input as peptide sequences")
    parser.add_argument("-m", "--method",
                        default="PCM",
                        help="Specifies the Cleavage Site prediction tool to use - default PCM."
                      )
    parser.add_argument("-s", "--solver",
                        default="glpk",
                        help="Specifies the ILP solver to use (must be installed) - default glpk"
                      )
    parser.add_argument("-o", "--output",
                        required=True,
                        help="Specifies the output path. Results will be written to CSV")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Specifies verbose run."
                      )
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
        for pred, obj in ACleavageSitePrediction.registry.iteritems():
            if pred not in ["ACleavageSitePrediction", "APSSMCleavageSitePredictor"]:
                print "Method: ",pred
                print "Supported Length: ", " ".join(map(str, getattr(obj,  "_"+pred+"__supported_length")))
                print
        sys.exit(0)


    if args.pepfile:
        peptides = read_lines(args.input, type="Peptide")

    else:
        peptides = [Peptide(s) for s in args.input]

    cleav_pred = CleavageSitePredictorFactory(args.method)
    assembler = EpitopeAssembly(peptides, cleav_pred, solver=args.solver, verbosity=1 if args.verbose else 0)
    result = assembler.solve()
    output = args.output if args.outdir == "" else args.outdir + os.path.basename(args.output)
    with open(output, "w") as out:
        out.write("Order,Peptide\n")
        out.write("\n".join("%i,%s"%(i+1,str(p)) for i, p in enumerate(result)))


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

        setting = """  <h2 class="etk-heading">Peptide String-of-Beats Design</h2>

        <table class="etk-parameterT">
            <tr> <th class ="etk-innerHeading" colspan="2"> Parameters </th></tr>
            <tr>
                <th>Cleave Site Prediction Method:</th>
                <td>%s</td>
            </tr>
        </table>"""%args.method


        table="""

        <input id="etk-search" placeholder="  filter">
        <table class="etk-sortT etk-resultsT etk-filterT">

            <thead>
                <tr>
                    <th>Order</th><th>Peptide</th>
                </tr>
            </thead>"""+"".join("<tr><td>%i</td><td>%s</td></tr>"%(i+1,p) for i, p in enumerate(result))+"</table>"

        end_html = "</div></body></html>"

        html_out = ".".join(output.split(".")[:-1])+".html"
        with open(html_out, "w") as html_o:
            html_o.write(begin_html+setting+table+end_html)

if __name__ == "__main__":
    main()
