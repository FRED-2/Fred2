"""Defines globally some variants and matation syntaxes 
for testing purposes.
"""

from Fred2.Core.Variant import Variant, VariationType, MutationSyntax

# DATA:     M U T A T I O N   S Y N T A X
# ======================================================================
mut_syn1_1 = MutationSyntax("tsc_1", # transcript_id
                            2,       # pos in transc
                            1,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

mut_syn1_2 = MutationSyntax("tsc_1", # transcript_id
                            5,       # pos in transc
                            2,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

mut_syn1_3 = MutationSyntax("tsc_1", # transcript_id
                            7,       # pos in transc
                            3,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

mut_syn1_4 = MutationSyntax("tsc_1", # transcript_id
                            9,      # pos in transc
                            3,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

mut_syn1_5 = MutationSyntax("tsc_1", # transcript_id
                            10,      # pos in transc
                            4,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

mut_syn1_6 = MutationSyntax("tsc_1", # transcript_id
                            10,      # pos in transc
                            7,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

# DATA:     V A R I A N T S
# ======================================================================
# assuming:
#   genome size:            30 bp, 10 aa
#   transcripts:            1 (30 bp, 10 aa)
# (as genomic position I actually used the protein-position)

# SNP, HOMOZYGOUS
var_1 = Variant(
    "var_1", VariationType.SNP, "chr1",
    20,                    # Genomic Position 
    "C", "T",             # reference , observed
    {"tsc_1":mut_syn1_1}, # dict of (transcrip_id : mutSnytaxes)
    True, False)          # isHomozygous?, isSynonymous?

# SNP, HOMOZYGOUS
var_2 = Variant(
    "var_2", VariationType.SNP, "chr1",
    5,                    # Genomic Position 
    "C", "T",             # reference , observed
    {"tsc_1":mut_syn1_2}, # dict of (transcrip_id : mutSnytaxes)
    True, False)         # isHomozygous?, isSynonymous?

# +F-SHIFT, HOMOZYGOUS
var_3 = Variant(
    "var_3", VariationType.FSINS, "chr1",
    7,                    # Genomic Position 
    "", "TT",             # reference , observed
    {"tsc_1":mut_syn1_3}, # dict of (transcrip_id : mutSnytaxes)
    True, False)          # isHomozygous?, isSynonymous?

# -F-SHIFT, HOMOZYGOUS
var_4 = Variant(
    "var_4", VariationType.FSDEL, "chr1",
    9,                    # Genomic Position 
    "CCCCC", "",          # reference , observed
    {"tsc_1":mut_syn1_4}, # dict of (transcrip_id : mutSnytaxes)
    True, False)          # isHomozygous?, isSynonymous?

# INSERTION, HOMOZYGOUS
var_5 = Variant(
    "var_5", VariationType.INS, "chr1",
    10,                   # Genomic Position 
    "", "TTT",            # reference , observed
    {"tsc_1":mut_syn1_5}, # dict of (transcrip_id : mutSnytaxes)
    True, False)          # isHomozygous?, isSynonymous?

# -F-SHIFT, HOMOZYGOUS, INVALID
var_6 = Variant(
    "var_6", VariationType.FSDEL, "chr1",
    10,                   # Genomic Position 
    "CC", "",             # reference , observed
    {"tsc_1":mut_syn1_6}, # dict of (transcrip_id : mutSnytaxes)
    True, False)          # isHomozygous?, isSynonymous?


# SNP, HOMOZYGOUS
mut_syn1_7 = MutationSyntax("tsc_1", # transcript_id
                            9,      # pos in transc
                            3,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

var_7 = Variant(
    "var_7", VariationType.SNP, "chr1",
    9,                    # Genomic Position 
    "C", "T",             # reference , observed
    {"tsc_1":mut_syn1_7}, # dict of (transcrip_id : mutSnytaxes)
    True, False)         # isHomozygous?, isSynonymous?


# INSERTION, HOMOZYGOUS
mut_syn1_8 = MutationSyntax("tsc_1", # transcript_id
                            14,      # pos in transc
                            4,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

var_8 = Variant(
    "var_8", VariationType.INS, "chr1",
    14,                   # Genomic Position 
    "", "TTT",            # reference , observed
    {"tsc_1":mut_syn1_8}, # dict of (transcrip_id : mutSnytaxes)
    True, False)          # isHomozygous?, isSynonymous?
