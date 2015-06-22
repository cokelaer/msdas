
from easydev import gsf

__all__ = ['get_tcell_filenames']

def get_tcell_filenames():
    x1 = gsf("msdas", "data", "donor_1_processed.csv")
    x2 = gsf("msdas", "data", "donor_2_processed.csv")
    x3 = gsf("msdas", "data", "donor_3_processed.csv")
    return [x1,x2,x3]

