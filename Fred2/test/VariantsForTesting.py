"""Defines globally some variants and matation syntaxes 
for testing purposes.
"""

from Fred2.Core.Variant import Variant, VariationType, MutationSyntax

# DATA:     M U T A T I O N   S Y N T A X
# ======================================================================
mut_syn1_1 = MutationSyntax("tsc_1", # transcript_id
                            1,       # pos in transc
                            0,       # pos in protein
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
                            5,       # pos in transc
                            3,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

mut_syn1_5 = MutationSyntax("tsc_1", # transcript_id
                            8,      # pos in transc
                            2,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

mut_syn1_6 = MutationSyntax("tsc_1", # transcript_id
                            8,      # pos in transc
                            2,       # pos in protein
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
    1,                    # Genomic Position
    "A", "T",             # reference , observed
    {"tsc_1":mut_syn1_1}, # dict of (transcrip_id : mutSnytaxes)
    True, False, experimentalDesign="tumor")          # isHomozygous?, isSynonymous?

# SNP, HOMOZYGOUS
var_2 = Variant(
    "var_2", VariationType.SNP, "chr1",
    5,                    # Genomic Position 
    "C", "T",             # reference , observed
    {"tsc_1":mut_syn1_2}, # dict of (transcrip_id : mutSnytaxes)
    True, False, experimentalDesign="normal")         # isHomozygous?, isSynonymous?

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
    5,                    # Genomic Position
    "CCCCC", "",          # reference , observed
    {"tsc_1":mut_syn1_4}, # dict of (transcrip_id : mutSnytaxes)
    True, False)          # isHomozygous?, isSynonymous?

# INSERTION, HOMOZYGOUS
var_5 = Variant(
    "var_5", VariationType.INS, "chr1",
    8,                   # Genomic Position
    "", "TTT",            # reference , observed
    {"tsc_1":mut_syn1_5}, # dict of (transcrip_id : mutSnytaxes)
    True, False)          # isHomozygous?, isSynonymous?

# -F-SHIFT, HOMOZYGOUS, INVALID
var_6 = Variant(
    "var_6", VariationType.FSDEL, "chr1",
    8,                   # Genomic Position
    "CC", "",             # reference , observed
    {"tsc_1":mut_syn1_6}, # dict of (transcrip_id : mutSnytaxes)
    True, False)          # isHomozygous?, isSynonymous?


# SNP, HOMOZYGOUS
mut_syn1_7 = MutationSyntax("tsc_1", # transcript_id
                            10,      # pos in transc
                            3,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

var_7 = Variant(
    "var_7", VariationType.SNP, "chr1",
    10,                    # Genomic Position
    "G", "T",             # reference , observed
    {"tsc_1":mut_syn1_7}, # dict of (transcrip_id : mutSnytaxes)
    True, False)         # isHomozygous?, isSynonymous?


# INSERTION+3, HOMOZYGOUS, pos 14
mut_syn1_8 = MutationSyntax("tsc_1", # transcript_id
                            15,      # pos in transc
                            5,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

var_8 = Variant(
    "var_8", VariationType.INS, "chr1",
    15,                   # Genomic Position
    "", "TTT",            # reference , observed
    {"tsc_1":mut_syn1_8}, # dict of (transcrip_id : mutSnytaxes)
    True, False)          # isHomozygous?, isSynonymous?


# INSERTION+2, HOMOZYGOUS, pos 2
mut_syn1_9 = MutationSyntax("tsc_1", # transcript_id
                            2,      # pos in transc
                            1,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

var_9 = Variant(
    "var_9", VariationType.INS, "chr1",
    2,                   # Genomic Position 
    "", "TT",            # reference , observed
    {"tsc_1":mut_syn1_9}, # dict of (transcrip_id : mutSnytaxes)
    True, False)          # isHomozygous?, isSynonymous?


################################################################################
#      H E T E R O Z Y G O U S     V A R I A N T S
################################################################################

# 3-DEL(-2), HETEROZYGOUS, pos 3
mut_syn1_10 = MutationSyntax("tsc_1", # transcript_id
                             0,       # pos in transc
                             0,       # pos in protein
                             "",      # cdsMutationSyntax - irrelevant
                             "")      # aaMutationSyntax - irrelevant
mut_syn2_10 = MutationSyntax("tsc_2", # transcript_id
                             0,       # pos in transc
                             0,       # pos in protein
                             "",      # cdsMutationSyntax - irrelevant
                             "")      # aaMutationSyntax - irrelevant

var_10 = Variant(
    "var_10", VariationType.DEL, "chr1",
    0,                     # Genomic Position
    "AA", "",              # reference , observed
    {"tsc_1":mut_syn1_10,
     "tsc_2":mut_syn2_10}, # dict of (transcrip_id : mutSnytaxes)
    False, False)          # isHomozygous?, isSynonymous?


# 5-INS(+3), HOMOZYGOUS, pos 5
mut_syn1_11 = MutationSyntax("tsc_1", # transcript_id
                             5,       # pos in transc
                             1,       # pos in protein
                             "",      # cdsMutationSyntax - irrelevant
                             "")      # aaMutationSyntax - irrelevant
mut_syn2_11 = MutationSyntax("tsc_2", # transcript_id
                             5,       # pos in transc
                             1,       # pos in protein
                             "",      # cdsMutationSyntax - irrelevant
                             "")      # aaMutationSyntax - irrelevant
var_11 = Variant(
    "var_11", VariationType.INS, "chr1",
    5,                     # Genomic Position
    "", "TTT",             # reference , observed
    {"tsc_1":mut_syn1_11,
     "tsc_2":mut_syn2_11}, # dict of (transcrip_id : mutSnytaxes)
    True, False)           # isHomozygous?, isSynonymous?


# 7-DEL(-4), HETEROZYGOUS, pos 7
mut_syn1_12 = MutationSyntax("tsc_1", # transcript_id
                             5,       # pos in transc
                             1,       # pos in protein
                             "",      # cdsMutationSyntax - irrelevant
                             "")      # aaMutationSyntax - irrelevant
mut_syn2_12 = MutationSyntax("tsc_2", # transcript_id
                             5,       # pos in transc
                             1,       # pos in protein
                             "",      # cdsMutationSyntax - irrelevant
                             "")      # aaMutationSyntax - irrelevant
var_12 = Variant(
    "var_12", VariationType.DEL, "chr1",
    5,                     # Genomic Position
    "CCCCC", "",            # reference , observed
    {"tsc_1":mut_syn1_12,
     "tsc_2":mut_syn2_12}, # dict of (transcrip_id : mutSnytaxes)
    False, False, experimentalDesign="normal")          # isHomozygous?, isSynonymous?


# -INS(+2), HETEROZYGOUS, pos 15
mut_syn2_13 = MutationSyntax("tsc_1", # transcript_id
                             15,      # pos in transc
                             5,       # pos in protein
                             "",      # cdsMutationSyntax - irrelevant
                             "")      # aaMutationSyntax - irrelevant
var_13 = Variant(
    "var_13", VariationType.INS, "chr1",
    15,                    # Genomic Position 
    "", "CC",              # reference , observed
    {"tsc_1":mut_syn2_13}, # dict of (transcrip_id : mutSnytaxes)
    False, True)          # isHomozygous?, isSynonymous?

# -INS(+1), HETEROZYGOUS, pos 1
mut_syn2_14 = MutationSyntax("tsc_1", # transcript_id
                             1,       # pos in transc
                             1,       # pos in protein
                             "",      # cdsMutationSyntax - irrelevant
                             "")      # aaMutationSyntax - irrelevant
var_14 = Variant(
    "var_14", VariationType.INS, "chr1",
    1,                     # Genomic Position 
    "", "C",              # reference , observed
    {"tsc_1":mut_syn2_14}, # dict of (transcrip_id : mutSnytaxes)
    False, False)          # isHomozygous?, isSynonymous?


mut_syn1_15 = MutationSyntax("tsc_1", # transcript_id
                            2,      # pos in transc
                            1,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant


##### Tumor tests
mut_syn1_n = MutationSyntax("tsc_1", # transcript_id
                            5,      # pos in transc
                            1,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

mut_syn1_t = MutationSyntax("tsc_1", # transcript_id
                            11,      # pos in transc
                            5,       # pos in protein
                            "",      # cdsMutationSyntax - irrelevant
                            "")      # aaMutationSyntax - irrelevant

var_n1 = Variant(
    "var_n1", VariationType.SNP, "chr1",
    5,                   # Genomic Position
    "A", "T",            # reference , observed
    {"tsc_1":mut_syn1_n}, # dict of (transcrip_id : mutSnytaxes)
    False, False, experimentalDesign="normal")          # isHomozygous?, isSynonymous?

var_n2 = Variant(
    "var_n2", VariationType.SNP, "chr1",
    5,                   # Genomic Position
    "A", "T",            # reference , observed
    {"tsc_1":mut_syn1_n}, # dict of (transcrip_id : mutSnytaxes)
    True, False, experimentalDesign="normal")

var_t1 = Variant(
    "var_t1", VariationType.SNP, "chr1",
    11,                   # Genomic Position
    "G", "T",            # reference , observed
    {"tsc_1":mut_syn1_t}, # dict of (transcrip_id : mutSnytaxes)
    False, False, experimentalDesign="tumor")

var_t2 = Variant(
    "var_t2", VariationType.SNP, "chr1",
    11,                   # Genomic Position
    "G", "T",            # reference , observed
    {"tsc_1":mut_syn1_t}, # dict of (transcrip_id : mutSnytaxes)
    True, False, experimentalDesign="tumor")