from Fred2.EpitopePrediction import EpitopePredictorFactory, Peptide

peptides = [Peptide("SYFPEITHI"), Peptide("FIASNGVKL"), Peptide("LLGATCMFV")]
print EpitopePredictorFactory("mhcflurry").predict(peptides)
print EpitopePredictorFactory("mhcnuggets-class-1").predict(peptides)
print EpitopePredictorFactory("mhcnuggets-class-2").predict(peptides)
