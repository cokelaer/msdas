from msdas import *



def test_phosphogrid():
    m = MassSpecReader(get_yeast_small_data(), verbose=False)
    gene_names = set(list(m.df.Protein))
    p = phosphogrid.PhosphoGRID(directory = "../share/data")
    p.run(gene_names=gene_names)
    p.export2sif()
    p.plot()

    p.run()
