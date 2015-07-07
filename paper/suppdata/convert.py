from cno.io import *
s = sif.SIF("PKN-YEAST.sif")
sbml = s.to_sbmlqual("PKN-YEAST.xml")
c = cnograph.CNOGraph(s, "MD-YEAST.csv")
c2 = cnograph.CNOGraph("PKN-YEAST.sif", "MD-YEAST.csv")
print c2 == c

