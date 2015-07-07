import pandas as pd
from msdas import *

from easydev import gsf

def test_alignmentbase():
    m = alignment.MassSpecAlignmentBase()
    m.df = pd.DataFrame({'Psite':[1,2]})
    m.df = pd.DataFrame({'Protein':["ZAP70"]})
    print(m) 




def test_alignment_yeast():

    filenames = get_yeast_filenames()
    m = alignment.MassSpecAlignmentYeast(yeast.get_yeast_filenames(), verbose=False)
    m.check_format()
    m.mode
