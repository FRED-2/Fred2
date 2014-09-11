__author__ = 'schubert'

"""
pcm_matrices: Matrix file for proteasomal cleavage prediction.
Reference to original method: Doennes, P. and Kohlbacher, O. (2005) Integrated modeling of the major events in the
MHC class I antigen processing pathway. Protein Sci
"""

pcm_6 = {0: {'A': 0.4731237565819792, 'C': -0.5025268209512956, 'D': -0.5025268209512956, 'E': -0.21691300156357377,
           'F': -0.3495574761698686, 'G': -2.2537949288246137, 'H': 1.1647120903726331, 'I': 0.004987541511038968,
           'K': -0.3495574761698686, 'L': 0.09984533496971612, 'M': 0.004987541511038968, 'N': 0.18647956694261839,
           'P': -0.3495574761698686, 'Q': -0.6831968497067774, 'R': 0.2662030407746567, 'S': 0.2662030407746567,
           'T': -0.21691300156357377, 'V': 0.004987541511038968, 'W': -0.903868211875598, 'Y': 0.2662030407746567},
       1: {'A': 0.4731237565819792, 'C': -0.5025268209512956, 'D': -0.0998203352822109, 'E': -0.0998203352822109,
           'F': -0.5025268209512956, 'G': -5.298317366548036, 'H': 0.7907275088988094, 'I': -0.0998203352822109,
           'K': -0.21691300156357377, 'L': 0.2662030407746567, 'M': -0.21691300156357377, 'N': 0.7907275088988094,
           'P': -0.5025268209512956, 'Q': -0.903868211875598, 'R': 0.2662030407746567, 'S': 0.2662030407746567,
           'T': -0.3495574761698686, 'V': -0.3495574761698686, 'W': -0.903868211875598, 'Y': 0.5335651107354802},
       2: {'A': 0.5335651107354802, 'C': 0.004987541511038968, 'D': -0.6831968497067774, 'E': 0.18647956694261839,
           'F': 0.40879289820083897, 'G': -2.2537949288246137, 'H': 1.1002775679871708, 'I': 0.09984533496971612,
           'K': -0.6831968497067774, 'L': -0.21691300156357377, 'M': -1.1874435023747254, 'N': 0.3400373027857091,
           'P': -0.0998203352822109, 'Q': -0.6831968497067774, 'R': 0.4731237565819792, 'S': -1.1874435023747254,
           'T': 0.004987541511038968, 'V': -0.0998203352822109, 'W': -2.2537949288246137, 'Y': 0.40879289820083897},
       3: {'A': 0.7443154671343447, 'C': 0.09984533496971612, 'D': 0.3400373027857091, 'E': -0.903868211875598,
           'F': -0.21691300156357377, 'G': -5.298317366548036, 'H': 0.5335651107354802, 'I': 0.2662030407746567,
           'K': -0.5025268209512956, 'L': 0.6444820085786643, 'M': -0.3495574761698686, 'N': -0.3495574761698686,
           'P': 0.004987541511038968, 'Q': -0.6831968497067774, 'R': -0.3495574761698686, 'S': -0.903868211875598,
           'T': -0.6831968497067774, 'V': 0.3400373027857091, 'W': -0.0998203352822109, 'Y': 0.5905605917848442},
       4: {'A': 0.8775499035577246, 'C': 0.004987541511038968, 'D': -0.21691300156357377, 'E': 0.18647956694261839,
           'F': -0.21691300156357377, 'G': -5.298317366548036, 'H': 1.0314035389746596, 'I': -0.0998203352822109,
           'K': -0.6831968497067774, 'L': 0.004987541511038968, 'M': -1.584745299843729, 'N': 0.40879289820083897,
           'P': -2.2537949288246137, 'Q': -0.5025268209512956, 'R': 0.4731237565819792, 'S': 0.2662030407746567,
           'T': 0.004987541511038968, 'V': 0.18647956694261839, 'W': -0.6831968497067774, 'Y': -0.6831968497067774},
       5: {'A': 0.3400373027857091, 'C': -0.21691300156357377, 'D': 0.3400373027857091, 'E': 0.5905605917848442,
           'F': 0.18647956694261839, 'G': -5.298317366548036, 'H': 0.9182887345368281, 'I': 0.004987541511038968,
           'K': -0.21691300156357377, 'L': 0.40879289820083897, 'M': 0.004987541511038968, 'N': -0.903868211875598,
           'P': -0.21691300156357377, 'Q': -1.584745299843729, 'R': -0.21691300156357377, 'S': -0.3495574761698686,
           'T': 0.18647956694261839, 'V': -0.21691300156357377, 'W': -0.5025268209512956, 'Y': 0.004987541511038968}}

ginodi_11 = {
0: {'A': -1.511391, 'C': -0.360375, 'E': -0.940226, 'D': -1.359827, 'G': -3.095267, 'F': -1.355572, 'I': -0.499191,
    'H': -2.819591, 'K': -3.05951, 'M': -0.423525, 'L': -0.733422, 'N': -0.266152, 'Q': -1.241387, 'P': -3.237045,
    'S': -1.239708, 'R': -0.164552, 'T': -3.817465, 'W': -0.370045, 'V': -1.016888, 'Y': -1.278611},
1: {'A': 0.194364, 'C': -3.96114, 'E': 0.434893, 'D': 0.939321, 'G': 0.740418, 'F': -0.109587, 'I': 0.662333,
    'H': 0.080279, 'K': 0.71949, 'M': -0.564399, 'L': 0.028457, 'N': 1.702132, 'Q': -0.630433, 'P': 0.804683,
    'S': 0.435486, 'R': -0.159664, 'T': 0.325815, 'W': -1.631427, 'V': 0.497904, 'Y': -0.067161},
2: {'A': -0.561323, 'C': -0.05371, 'E': -1.533472, 'D': -1.963718, 'G': -1.471855, 'F': -0.985792, 'I': -3.80482,
    'H': -1.581678, 'K': 0.289743, 'M': 1.08178, 'L': -4.341714, 'N': -0.368058, 'Q': -0.540673, 'P': -2.052335,
    'S': -0.815012, 'R': -1.955067, 'T': -1.949179, 'W': 0.731348, 'V': -3.413529, 'Y': -2.227598},
3: {'A': -0.028153, 'C': -2.76344, 'E': -1.236879, 'D': -1.510436, 'G': 0.472955, 'F': -1.151859, 'I': -1.788194,
    'H': -0.673313, 'K': -1.527173, 'M': -1.456656, 'L': -2.048402, 'N': 1.173897, 'Q': -2.466156, 'P': -0.761433,
    'S': 0.454867, 'R': 0.440716, 'T': -1.036751, 'W': -1.112729, 'V': -1.721837, 'Y': -0.387162},
4: {'A': 0.208343, 'C': -2.144223, 'E': 0.05252, 'D': -0.590682, 'G': -3.14597, 'F': 1.531049, 'I': -0.312761,
    'H': -1.420491, 'K': -0.792268, 'M': 2.490757, 'L': 2.348379, 'N': -2.899709, 'Q': -0.547701, 'P': -1.549403,
    'S': -0.208625, 'R': 0.016924, 'T': -6.807949, 'W': 0.124181, 'V': 1.287458, 'Y': 2.106059}}