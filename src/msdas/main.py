


def getfile(filename):
    from easydev import get_share_file
    return get_share_file("msdas", "data/upstream/", filename)

from cellnopt.wrapper.wrapper_cnor import *
from cellnopt.wrapper import cnor

cnolist0 = CNOlist(getfile("a0Up.csv"))
cnolist1 = CNOlist(getfile("a1Up.csv"))
cnolist2 = CNOlist(getfile("a5Up.csv"))
cnolist3 = CNOlist(getfile("a10Up.csv"))
cnolist4 = CNOlist(getfile("a20Up.csv"))
cnolist5 = CNOlist(getfile("a40Up.csv"))
salt = CNOlist(getfile("salt0Up.csv"))
alldata = [cnolist0, cnolist1, cnolist2, cnolist3, cnolist4, cnolist5, salt]

alldata = [normaliseCNOlist(x, detection=-1000, mode="raw") for x in alldata]


pknmodel = getfile("TopoUpstreamX.sif")


pknmodel = getfile("TopoUpstreamXNeg.sif")
o = cnor.CNORxy(pknmodel, alldata)


o.run(maxteme=10, ndiverse=10, dim_refset=10)
o.plotFit()
o.plotFitness(1)
o.boxplot()


o.plot_k(fontsize=8)
o.plot_n(fontsize=8)




cnolist0 = CNOlist("a0UNSCALEDmidas.csv")
cnolist1 = CNOlist("a1UNSCALEDmidas.csv")
cnolist2 = CNOlist("a5UNSCALEDmidas.csv")
cnolist3 = CNOlist("a10UNSCALEDmidas.csv")
cnolist4 = CNOlist("a20UNSCALEDmidas.csv")
cnolist5 = CNOlist("a40UNSCALEDmidas.csv")
salt = CNOlist("salt0UNSCALEDmidas.csv")
alldata = [cnolist0, cnolist1, cnolist2, cnolist3, cnolist4, cnolist5, salt]
pknmodel = getfile("../Topo2.sif")

o = cnor.CNORxy(pknmodel, alldata)
o.scaling(method="max")

