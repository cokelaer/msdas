from msdas import  *
from easydev import gsf

filename = gsf("msdas", "data", "YEAST_raw_sample.csv")


def test_replicates_plotting():
    r = replicates.Replicates(filename, verbose=False)

    r.boxplot_replicates_by_experiment("a0_t0")
    r.boxplot_replicates_by_experiment("a0_t0", indices=r.df.index[0:20])

    r.hist_coefficient_variation()
    r.plot_mu_sigma()

    r.hist2d_mu_versus_cv()
    r.hist2d_errors_versus_cv()

    r.hist_na_per_experiments()
    r.pcolor_na()
    assert r == r.copy()


def test_yeast():
    r = replicates.ReplicatesYeast(get_yeast_raw_data(), verbose=False)
    assert r.get_mu_df().mean().mean()>1
    r.normalise()
    # should be small value below 0
    assert r.get_mu_df().mean().mean()<1
    r.plot_timeseries_midas("ABF1_S720")
    r.get_na_per_experiment_after_averaging()
    r.get_na_per_salt_experiment_after_averaging()
    r.clean_na_paper()
    r.clean_na_paper2()


def test_replicates_others():
    r = replicates.Replicates(filename, verbose=False)
    assert r.df.shape == (198,23)
    r.average()
    assert r.df.shape == (198,11)
    r.reset()
    assert r.df.shape == (198,23)

    assert r == r.copy()

    r.set_irrelevant_replicates_to_na()
    r.reset()

    r.get_average_sigma()
    r.get_average_mu()
    r.get_standarised_residuals()


