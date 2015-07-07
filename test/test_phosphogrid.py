from msdas import *
from easydev import TempFile


def test_phosphogrid():
    m = MassSpecReader(get_yeast_small_data(), verbose=False)
    gene_names = set(list(m.df.Protein))
    p = phosphogrid.PhosphoGRID(directory = "../share/data")
    p.run(gene_names=gene_names)
    fh = TempFile(suffix='.sif')
    p.export2sif(filename=fh.name)
    p.plot()
    #p.run()
    
    
    fh.delete()
