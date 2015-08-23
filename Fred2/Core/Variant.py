# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Variant
   :synopsis: Contains relevant classes describing variants.
.. moduleauthor:: schubert, szolek, walzer

"""
import math

from Fred2.Core.Base import MetadataLogger


VariationType = (lambda **enums: type('Enum', (), enums))(SNP=0,
                                                          DEL=1,
                                                          INS=2,
                                                          FSDEL=3,
                                                          FSINS=4,
                                                          UNKNOWN=5)
"""
Enum for variation types:
type.SNP, type.DEL, type.INS, type.FSDEL, type.FSINS, type.UNKNOWN
"""


class MutationSyntax():
    """
    This class represents the mutation syntax of a variant and stores its 
    transcript and protein position

    :param str transID: the transcript id
    :param int transPos: the position of the variant within the transcript
    :param int protPos: the protein position of the variant within the 
                        transcript
    :param str cds: the complete cds_mutation_syntax string
    :param str aas: the complete protein_mutation_syntax string
    """
    def __init__(self, transID, transPos, protPos, cds, aas):
        #TODO: is protPos always given? what about synonymous variants?
        self.transID = transID
        self.tranPos = transPos
        self.protPos = protPos
        self.cdsMutationSyntax = cds  #c. ...
        self.aaMutationSyntax = aas  #p. ...


class Variant(MetadataLogger):
    """
    A Variant contains information about a single genetic modification of
    the reference genome.

    :param str id: Variant id
    :param Enum.VariationType type: An Enum type of the variant either SNP, 
                                    DEL, or INS
    :param str chrom: The chromosome on which the variant lies
    :param int genomePos: The genomic position of the variant
    :param str ref: The reference seq at the genomic position
    :param str obs: The observed variation at the genomic position
    :param dict(str,MutationSyntax) coding: A dictionary of associated 
                                            transcripts. Key=transcript_id, 
                                            value=MutationSyntax
    :param bool isHomozygous: defines if variant is homozygous or not
    :param bool isSynonymous: defines if variant is a synonymous mutation 
                              or not
    :param dict(str,int) offsets: the position offset this variant has in a 
                                  specific transcript-variant
                                  key=transcript-variant-id (xxx:FRED2_nn)
                                  value=offset
    :param defaultdict(list) metadata: meta information (not relevant for core
                                       functionality of Fred2)

    """
    def __init__(self, id, type, chrom, genomePos, ref, obs, coding,
                 isHomozygous, isSynonymous, experimentalDesign=None, metadata=None):
        """
        Constructor for a variant, see init-types in class parameters
        """
        MetadataLogger.__init__(self)
        self.id = id
        self.type = type
        self.chrom = chrom
        self.genomePos = genomePos
        self.ref = ref
        self.obs = obs
        self.gene = None
        self.isHomozygous = isHomozygous
        self.isSynonymous = isSynonymous
        self.offsets = {}
        self.coding = coding  # dict transcript_id:MutationSyntax
        self.experimentalDesign = "" if experimentalDesign is None else experimentalDesign

        if metadata is not None:
            for meta in metadata:
                self.log_metadata(meta, metadata[meta])

    def __repr__(self):
        return "Variant(g.%i%s>%s):%s" % (self.genomePos, self.ref, self.obs, self.experimentalDesign) \
            if self.experimentalDesign else "Variant(g.%i%s>%s)" % (self.genomePos, self.ref, self.obs)

    def get_transcript_offset(self):
        return len(self.obs) - len(self.ref)

    def get_shift(self):
        return self.get_transcript_offset() % 3

    def get_transcript_position(self, trans_variant_id):
        """
        .. note:: May only be used for transcript variants that were created from this variant via :func:
        `Generator.generate_transcripts_from_variants`.
        ..TODO:: what? that should not be encoded in transcript id but in separate member

        returns the specific transcript position of a given transcript_id. 
        If variant is not associated with the given transcript id the function 
        throws a KeyError

        :param str transcriptId: A transcript_id
        :return: (int) -- transcript position
        :raises: KeyError
        """
        trans_id = trans_variant_id.split(":FRED2_")[0]
        try:
            return self.coding[trans_id].tranPos + \
                   self.offsets.get(trans_variant_id, 0)
        except KeyError:
            raise KeyError("Transcript ID %s not associated with variant %s"%
                           (trans_variant_id, str(self)))

    def get_protein_position(self, trans_variant_id):
        """
        .. note:: May only be used for transcript variants that were created from this variant via :func:
        `Generator.generate_transcripts_from_variants`.
        ..TODO:: what? that should not be encoded in transcript id but in separate member

        returns the specific protein position of a given transcript_id. If 
        variant is not associated with the given transcript id the function 
        throws a KeyError

        :param str transcriptId: A transcript_id
        :return: (int) -- the protein position of the variant
        :raises: KeyError
        """
        trans_id = trans_variant_id.split(":FRED2_")[0]
        try: 
            # get actual transcript position
            tpos = self.coding[trans_id].tranPos + \
                   self.offsets.get(trans_variant_id, 0)
        except KeyError:
            raise KeyError("Transcript ID %s not associated with variant %s"%
                           (str(trans_variant_id), self.__str__))
        return int(math.ceil(tpos/3.0)) # generate protein pos from transcript pos
