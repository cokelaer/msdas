from msdas import *
from easydev import gsf

data = yeast.YEAST2MIDAS(get_yeast_small_data(), get_yeast_raw_data())


def _test_files():
    assert len(get_yeast_filenames()) == 6
    #assert len(get_yeast_filenames("all")) == 6
    try:
        get_yeast_filenames("dummy")
        assert False
    except:
        assert True

def _test_yeast_pkn():

    filename = gsf("msdas", "data", "PKN-yeast.sif")
    df = data.get_df_exemplars(-1)
    cnograph = data.get_expanded_cnograph(filename, df)
    c,m,e = data.export_pkn_and_midas(filename)



def test_yeast_june():
    #y = yeast.YEAST2MIDAS(get_yeast_small_data(), get_yeast_raw_data(),  verbose=False)
    #y.cleanup_june()
    #y.cleanup_june()
    #len(y.df)<100
    filename = gsf("msdas", "data", "PKN-yeastScaffold.sif")
    data.cleanup_june()
    c,m,e = data.export_pkn_and_midas_june(filename)

    from easydev import TempFile
    f = TempFile()
    data.to_midas(f.name)
    f.delete()

    cv = data.get_cv()
    m = data.get_midas()
    data.pcolor_na()
    data.plot_timeseries("DIG1_S126+S127")

def _test_yeast_cluster():
    data.get_group_psite("DIG1")
    data.get_group_psite_transposed("DIG1")
    psites = data.get_psites_exemplars(preference=-30)
    psites = data.get_psites_mapping(preference=-30)
    g = data.groups
    data.cluster.run_affinity_propagation_clustering()
    g = data.groups_cluster()

def _test_yest_Errors():

    data.get_errors()
    data.pcolor_errors()



