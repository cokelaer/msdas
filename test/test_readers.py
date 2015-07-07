from msdas import MassSpecReader, Cleaner, yeast, get_yeast_small_data
import pandas as pd
from easydev import  gsf

verbose=False

def test_Cleaner():
    c = Cleaner()
    c.df = pd.DataFrame({'A': [0, 1] , 'B':[4,5]})
    assert pd.isnull(c.df).sum().sum() == 0
    c.set_zero_to_na()
    assert pd.isnull(c.df).sum().sum() == 1

def test_MSReader():

    # we can just create an instance
    r = MassSpecReader(verbose=verbose)

    # fails if wrong file
    try:
        r = MassSpecReader("dummy.csv", verbose=verbose)
        assert False
    except:
        assert True


    filename = yeast.get_yeast_filenames()[0]
    r = MassSpecReader(filename, verbose=verbose)
    print(r)
    r.mode
    r.N
    r.df
    r.measurements
    r.metadata


    try:
        r.mode = None
        assert False
    except:
        assert True

    r.sort_psites_ors_only()
    r['DIG1']
    r['DIG1',"S142"]
    r['DIG1_S142']
    try:
        r['DIG1', 'S142', 'dummy']
        assert False
    except:
        assert True
    r.sequences
    r.psites

    from easydev import TempFile
    f = TempFile()
    r.to_csv(f.name)
    f.delete()


def test_plot():
    filename = yeast.get_yeast_filenames()[0]
    r = MassSpecReader(filename, verbose=verbose)
    r.pcolor("DIG1")
    r.status()
    r.plot_timeseries("DIG1_S142")

    r = MassSpecReader(get_yeast_small_data(), verbose=verbose)
    r.plot_experiments("DIG1_S142")
    r.read_annotations(gsf("msdas", "data", "YEAST_annotations_small.pkl"))
    r.show_sequence("DIG1")

def test_reader_yeast_small():
    filename = gsf("msdas", "data", "alpha0.csv")
    r = MassSpecReader(filename, verbose=verbose)
    assert len(r.df) == 57
    r.plot_phospho_stats()
    r.merge_peptides()
    assert len(r.df) == 57

    r.pcolor("DIG1", "t")
    r.hist_peptide_sequence_length()

def _test_reader_yeast_large():
    filename = gsf("msdas", "data", "Yeast_all_raw.csv")

    r = MassSpecReader(filename, cleanup=False,verbose=verbose)
    assert len(r.df) == 8895

    r = MassSpecReader(filename, cleanup=True,verbose=verbose)
    assert len(r.df) == 8570

    r.drop_oxidation()

    r.get_na_count()
    r.drop_na_count(50, inplace=False) # low values means that merge_peptides will be faster
    r.drop_na_count(50, inplace=True)
    r.merge_peptides()


def test_constructor():
    filename = gsf("msdas", "data", "alpha0.csv")
    r = MassSpecReader(filename,verbose=verbose)
    r2 = MassSpecReader(r,verbose=verbose)
    assert r == r2
    r3 = MassSpecReader(verbose=verbose)
    r3.read_csv(filename)
    r3.cleanup()
    assert r3 == r

    r4 = MassSpecReader(r.df, verbose=verbose)
    assert r4 == r




