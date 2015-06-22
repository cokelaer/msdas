import pandas as pd
from msdas import *

from easydev import gsf

def test_alignmentbase():
    m = alignment.MassSpecAlignmentBase()
    m.df = pd.DataFrame({'Psite':[1,2]})
    m.df = pd.DataFrame({'Protein':["ZAP70"]})
    print(m) 

def test_peptides_tcell():
    filename = get_tcell_filenames()
    #gsf("msdas", "data", "tcell_sample.csv")
    m = alignment.MassSpecAlignmentTCell(filename)
    m.check_format()
    m.df.Psite.ix[0] = "A3"
    try:
        m.check_format()
        assert False
    except:
        assert True




def test_alignment_yeast():

    filenames = get_yeast_filenames()
    m = alignment.MassSpecAlignmentYeast(yeast.get_yeast_filenames(), verbose=False)
    m.check_format()
    m.mode
