from msdas import *
from easydev import TempFile
import os
import pandas as pd

verbose=False
df = pd.DataFrame(
     {'Protein':['DIG1', 'LEU1', 'ASC1'],
      'Sequence_Phospho':['S(Phospho)APAQVTQHSK', 'VEVTS(Phospho)EDEK', 'DS(Phospho)VTIISAGNDK'],
      'Psite':['S142','S495', 'S166']})

def test_constructor():
    a = annotations.Annotations(df, "YEAST", verbose=verbose)
    try:
        a = annotations.Annotations(df)
        assert False
    except:
        assert True

def test_annotations_update():
    a = annotations.Annotations(df, "YEAST", verbose=verbose)

    try:
        a.plot_goid_histogram()
        assert False
    except:
        assert True

    a.get_uniprot_entries()
    a.set_annotations()
    a._mapping['LEU1_YEAST'] = ['P07264']
    a.get_uniprot_entries() # to update the main df with new entries
    a.set_annotations() # to retrieve the sequence of LEUC_YEAST
    a.set_annotations() # nothing to do here
    a.check_entries_versus_sequence()

def test_annotations_simple():
    r = MassSpecReader(get_yeast_small_data(), verbose=verbose)
    a = annotations.Annotations(r, "YEAST", verbose=verbose)

    a.get_uniprot_entries()
    a.set_annotations()

    try:
        a.to_pickle("test")
    except:
        if os.path.exists("YEAST_annotations_test.pkl"):
            os.remove("YEAST_annotations_test.pkl")
        a.to_pickle("test")
    finally:
            os.remove("YEAST_annotations_test.pkl")

    a.plot_goid_histogram()
    a.hist_most_relevant_goids()
    a.check_entries_versus_sequence()
    #a.find_sequence_blast(a.df.Sequence[0], "test@ebi.ac.uk")


def test_yeast_annotations():
    from easydev import gsf
    filename = gsf('msdas', "data", "YEAST_raw_sample.csv")
    r = MassSpecReader(filename, verbose=verbose)
    a = AnnotationsYeast(r, verbose=verbose)
    a.df = a.df.ix[0:200] # 200 is enough to get gene name cases and ambiguous gene names cases
    # e.g., ALD3_YEAST ['P54114', 'P40047']
    a.get_uniprot_entries()
    a.update_mapping()
    a.set_annotations()
    a.annotations.Sequence

    t = TempFile()
    a.to_csv(t.name)
    t.delete()

    t = TempFile()
    a.to_pickle("test", overwrite=True)
    try:
        a.to_pickle("test", overwrite=False)
        assert False
    except IOError:
        assert True
    a.read_pickle("YEAST_annotations_test.pkl")

    # create constructor given the annotations
    a = AnnotationsYeast(r, verbose=verbose, annotations="YEAST_annotations_test.pkl")
    a.get_uniprot_entries() # populate entry and entry_names in the df
    a.plot_goid_histogram()


    #cleanup
    os.remove("YEAST_annotations_test.pkl")




