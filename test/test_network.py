from msdas import *
from easydev import gsf

def test_network_uniprot():
    annotations = gsf("msdas", "data", "YEAST_annotations_small.pkl")
    r = MassSpecReader(get_yeast_small_data(), mode="YEAST", verbose=False)
    r.read_annotations(annotations)
    # test constructor 1
    n = network.NetworkFromUniProt(r.annotations)
    # test constructor 2
    n = network.NetworkFromUniProt(r)
    #m = MassSpecReader(get_yeast_small_data)
    n = network.NetworkFromUniProt(annotations)
    c = n.get_cnograph_intact(label="entry_id")
    c = n.get_cnograph_intact()
    c.plot()
    c.to_sif("PKN-uniprot.sif")

    names = list(set(r.df.Protein))
    n = network.CombineNetworks( 
        {"Curated": gsf("msdas", "data", "PKN-yeast.sif"),  
         "UniProt": "PKN-uniprot.sif",  
         #"PhosPho": "PKN-phospho.sif"},
         },
         signals=names[:], stimuli=["a", "NaCl"]) 


    #n.plot_multiedge_graph()
    #n.get_multiedge_graph()
    n.get_digraph()

    import os
    os.remove("PKN-uniprot.sif")




