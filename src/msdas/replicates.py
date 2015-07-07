# -*- python -*-
#
#  This file is part of MS-DAS software
#
#  Copyright (c) 2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#
##############################################################################
import pandas as pd
import numpy as np
import pylab
from matplotlib import colors


from msdas import readers

__all__ = ["ReplicatesYeast", "Replicates"]



class Replicates(readers.MassSpecReader):
    """Investigate the replicates in the raw dat sets

    The purpose of this class is to analyse a full data set with replicates,
    to remove rows that are meaningless, to compute mean/sigma on
    replicates and identify rows that have large errors.

    Replicates are identified as follows. If the input is a CSV file,
    columns with the same name are considered replicates to each other. Internally,
    in the dataframe, an index is appended. The first replicate does not have
    an index appended though. For instance, 3 columns called **X1** will be
    renamed internally as::

        X1
        X1.1
        X1.2

    If you save the dataframe into a CSV file, these convention is kept
    and names will be kept as it is.

    ::

        from msdas import *
        r = ReplicatesYeast(get_yeast_raw_data(), verbose=True)

    Since Replicates inherits from :class:`msdas.readers.MassSpecReader`, some
    functionalities can be re-used such as plotting, selection of data, annotations...

    Here are some notations used in the documentation.

    * :math:`i` represent a row that corresponds to a combination
      of protein name and psites.
    * There are possibly more than 1 experiment, which number is :math:`E` and indexed
      thanks to the :math:`e` letter.
    * There are :math:`R` replicates for each measurement. A replicate is
      indexed with :math:`r`.

    Using the notations above, one datum in the dataframe can be denoted as
    :math:`X_{i,e,r}`. We denote :math:`\mu` the mean and :math:`\sigma` the
    standard deviation.

    * :math:`\mu_{i,e} = E(X_{i,e,r})_r` is the mean over replicates and is a
      NxE matrix.
    * If you then want to average the resulting mean over experiment (to get the
      **grand mean**), we denote it as :math:`\mu_i = E(\mu_{i,e})_e`. This
      is equivalent to :math:`\mu_i = E(X_{i,e,r})_{e,r}`.

    Note however, that for standard deviation, the distinction and order of
    the indices matters.

    * :math:`\sigma_i = \sqrt{V(X_{i,e,r})_{r,e}}` is a N-length vector of standard
      deviations for each combi of protein and psite.
    * :math:`\sigma_{i,e} = \sqrt{V(X_{i,e,r})_r}` is a NxE matrix of standard deviations.

    There are quite a few methods to extract these information either as full
    dataframe or set of dataframe inside a dictionary. For instance, :meth:`get_mus`
    returns the mean over each


    """
    def __init__(self, data, verbose=False, cleanup=True):
        """.. rubric:: constructor

        :param data: a valid input to :class:`msdas.readers.MassSpecReader`
        :param cleanup: see :class:`msdas.readers.MassSpecReader` for details.


        """
        super(Replicates, self).__init__(data, verbose=verbose, cleanup=cleanup,
            merge_peptides=False)

        # hack. should be in readers but could not find where.
        self.df.fillna(np.nan, inplace=True)

        if "count_phospho" in self.df.columns:
            del self.df['count_phospho']

        # saves annotations
        self._rawdf = self.df.copy()

    def average(self):
        """Average the data over replicates and replace the dataframe inplace"""
        mu = pd.DataFrame(self.get_mus())
        self.df[mu.columns] = mu
        self._drop_replicates()

    def boxplot_replicates_by_experiment(self, tag, indices=None, fontsize=8):
        """BoxPlot of the replicates for a given experiment

        :param str tag: a valid measurement name to look for replicates.
        :param indices: list of indices to look at. If data set is large,
            you will want to use this parameter. Note that indices is not a range
            i,e., range(0,20) may not work because index 10 is not present.
        :param fontsize: fontsize used for the x-axis

        .. plot::
            :include-source:
            :width: 80%

            from msdas import *
            r = replicates.ReplicatesYeast(get_yeast_raw_data())
            r.boxplot_replicates_by_experiment("a0_t0", indices=r.df.index[0:20])

        """
        if tag not in self.get_unique_measurement_name():
            raise ValueError("invalid tag. call get_unique_measurement_name)")

        if indices is None:
            indices = self.df.index
        if len(indices) > 200:
            self.logging.info("boxplot with more than 200 variables will be slow. You can stop the plotting pressing CTRL+C")

        df = self.get_replicates_from_one_unique_measurement(tag)

        # tranpose to make protein the variable (column)
        df = df.ix[indices].transpose()
        # df.columns is equivalent to indices is principle
        # does not work; use xticks instead : df.columns = self.df.ix[df.columns].Identifier.values
        colnames = self.df.ix[df.columns].Identifier.values
        pylab.clf()
        self.dfdf = df
        df.boxplot(rot=90, return_type='axes')
        pylab.semilogy()
        pylab.xticks(range(0,len(indices)), colnames, rotation=90,fontsize=fontsize)
        pylab.tight_layout()
        return df

    def copy(self):
        """Return a copy of Replicates instance. """
        r = Replicates(self)
        return r

    def get_unique_measurement_name(self):
        """Returns unique column names without replicates.

        If there are replicates, they have the same column name in the CSV file
        but in the dataframe, there are encoded differently with an appended suffix.
        For instance, three replicates of the colum *t0* will be renamed *t0*, *t0.1*, *t0.2*.

        This method returns all unique names before the dot charater.
        For instance, with a dataframe that contains the following columns:
        *t0*, *t0.1*, *t0.2*, *t5*, *t5.1*, *t5.2*, this methods returns a list that
        contains *t0* and *t5* strings.

        .. note:: names are sorted alphabetically

        .. seealso:: :meth:`get_replicates_from_one_unique_measurement`
        """
        return sorted(list(set([x.split(".")[0] for x in self.measurements.columns])))

    def get_replicates_from_one_unique_measurement(self, name):
        """Given a unique name, returns a dataframe with only its replicates

        .. seealso:: :meth:`get_unique_measurement_name`
        """
        if name not in self.get_unique_measurement_name():
            raise ValueError("invalid tag")
        else:
            names = [c for c in self.df.columns if c.split(".")[0]==name]
            return self.df[names]


    def get_mu_df(self):
        """Version of :meth:`get_mus` returning a dataframe """

        mu = pd.DataFrame(self.get_mus())
        mu = mu[[x for x in self.measurements.columns if "." not in x]]
        return mu


    def get_mus(self):
        r"""Return mean of replicates as a dictionary

        Mathematical, it returns a NxE matrix computed as follows:

        .. math::

            \mu_{i,e} = E(X_{i,e,r})_r =  \frac{1}{R} \sum_{r=1}^{R} X_{i,e,r}


        :return: a dictionary. Keys are the E unique measurement names
            Values are each vector :math:`\mu_i` of length N protein/psites

        """
        tags = self.get_unique_measurement_name()
        mus = {}
        for tag in tags:
           mu = self.get_replicates_from_one_unique_measurement(tag).mean(axis=1)
           mus[tag] = mu.copy()
        return mus

    def get_sigmas(self):
        r"""Return standard deviation of replicates as a dictionary

        Mathematical, it returns a NxE matrix computed as follows :

        .. math::

            \sigma_{i,e} = V(X_{i,e,r})_r

        :return: a dictionary. Keys are the E unique measurement names
            Values are each vector :math:`\mu_i` of length N protein/psites

        """
        tags = self.get_unique_measurement_name()
        sigmas = {}
        for tag in tags:
           sigma = self.get_replicates_from_one_unique_measurement(tag).std(axis=1)
           sigmas[tag] = sigma.copy()
        return sigmas

    def get_mean_errors_df(self):
        """Returns dataframe with errors averaged (dataframe version)

        This is taking the results of :meth:`get_errors`, compute the average
        for each key/value and builds up a dataframe.

        :return: dataframe with averaged errors.

        .. seealso:: :meth:`get_errors`


        """
        errors = self.get_errors()
        errors = dict([(tag,errors[tag].mean(axis=1)) for tag in errors.keys()])
        errors = pd.DataFrame(errors)
        return errors

    def get_errors(self):
        r"""Returns dataframe with errors (dictionary version)

        Errors are computed as follows:

        .. math:: \epsilon_{i,e,r} =  \frac{\left| X_{i,e,r}- \mu_{i,e} \right|}{\mu_{i,e}}

        with

        .. math:: \mu_{i,e} = \frac{1}{R} \sum_{r=1}^{R} X_{i,e,r}


        In the dictionary, keys correspond to the name of a measurement as returned
        by :meth:`get_unique_measurement_name`. Values are dataframe with the
        errors for each replicate.

        """
        errors = {}
        mus = self.get_mus()
        tags = self.get_unique_measurement_name()
        for time in tags:
            data = self.get_replicates_from_one_unique_measurement(time)
            errors[time] = np.abs(data.sub(mus[time], axis="index")).divide(mus[time],
                axis="index")*100
        return errors

    def set_irrelevant_replicates_to_na(self):
        """Set unique replicate to NAs

        If an experiment has no replicates (0 or 1), you may want to set the
        experiment to NA. This is not relevant if there are zero replicates
        since the value may already be an NA but may make sense when there
        is only one replicate, for which no errors can be obtained.


        .. plot::
            :include-source:
            :width: 70%

            from msdas import *
            r = ReplicatesYeast(get_yeast_raw_data(), verbose=False)
            r.set_irrelevant_replicates_to_na()
            r.hist_na_per_experiments(color="r", alpha=0.5)
            r.reset()
            r.hist_na_per_experiments(color="g", alpha=0.5)

        """
        tags = self.get_unique_measurement_name()
        for tag in tags:
            df = self.get_replicates_from_one_unique_measurement(tag)
            indices = pd.notnull(df).sum(axis=1)<=1
            indices = [k for k,v in indices.iteritems() if v]
            colnames = [c for c in self.df.columns if c.split(".")[0]==tag]
            self.df.ix[indices, colnames] = np.nan




    def _drop_replicates(self):
        """Remove replicates from the dataframe"""
        todrop = [x for x in self.measurements.columns if "." in x]
        self.df.drop(todrop, inplace=True, axis=1)


    def get_average_mu(self):
        r"""Return grand mean

        :return: a vector with the grand mean over protein/psites

        The grand mean is the mean of the mean over replicates and experiments
        that is:

        .. math:: \mu_i = E(X_{i,e,r})_{e,r}

        This is equivalent to taking the mean over the experiment of the output
        of :meth:`get_mu_df`

        .. warning:: may contain NAs

        """
        mu = self.get_mu_df().mean(axis=1)
        return mu

    def get_average_sigma(self):
        r"""average of individual sigma on data measurements replicates

        :return: a vector with average sigma computed for each experiment.

        For instance, if you have 3 replicates (R=3) on 5 psites (N=5) for
        E=10 experiments, you get 5 values (average of 10 sigmas computed on 3 replicates)

        .. math:: \bar{\sigma}_{i} = E( \sigma_{e,i})_e

        .. note:: this is not the standard deviation over experiment and
            replicates :math:`\sigma_i = V( X_{i,e,r})_{e,r}` but the average
            of the standard deviation over replicates. See figure below.

        .. plot::
            :include-source:
            :width: 70%

            >>> import pylab
            >>> from msdas import *
            >>> r = replicates.ReplicatesYeast(get_yeast_raw_data())
            >>> x = pylab.linspace(1e8, 1e9, 100)
            >>> r.get_average_sigma().hist(bins=x, alpha=0.5, normed=True, label=r"$\bar{\sigma}_i$")
            >>> r.df.std(axis=1).hist(alpha=0.5, bins=x, normed=True, label=r"$\sigma_{i}$")
            >>> pylab.legend()

        .. warning:: may contain NAs
        """
        tags = self.get_unique_measurement_name()
        sigma = pd.concat([self.get_sigmas()[tag] for tag in tags], axis=1).mean(axis=1)
        return sigma

    def hist_na_per_experiments(self, bins=None, **kargs):
        """Histogram of the number of NA per rows

        :param kargs: optional parameter accepted by matplotlib.hist function.
        :param bins: binning of the histogram

        .. plot::
            :include-source:
            :width: 80%

            from msdas import *
            from easydev import gsf
            filename= gsf("msdas", "data", "YEAST_raw_sample.csv")
            r = Replicates(filename, verbose=False)
            r.hist_na_per_experiments()


        """
        if bins == None:
            bins = len(self.measurements.columns)
        data = self.get_na_count()
        data.hist(bins=bins, **kargs)
        pylab.xticks([x+0.5 for x in range(0,bins)], self.measurements.columns, rotation=90)

    def hist_coefficient_variation(self, merge=True, tags=None, **kargs):
        r"""Plots histogram of the coefficient variation

        That is, the quantity:

        .. math::

            \frac{\sigma_{i,e}}{\mu_{i,e}}

        .. plot::
            :include-source:
            :width: 70%

            import pylab
            import scipy.stats
            from msdas import *
            r = replicates.ReplicatesYeast(get_yeast_raw_data())
            r.hist_coefficient_variation(normed=True)
            x = pylab.linspace(0, 1, 200)
            pylab.plot(x, scipy.stats.gamma.pdf(x,2, scale=0.068), lw=2, color='r')
            pylab.xlim([0,1])

        """
        kargs['bins'] = kargs.get('bins', 100)
        if tags == None:
            tags = self.get_unique_measurement_name()

        cv = self.get_coefficient_variation()

        coeffs = {}
        if merge == False:
            for i, key in enumerate(tags):
                pylab.hist(cv[key].dropna(), label=key, **kargs)
            pylab.legend()
        else:
            coeffs = []
            for i, key in enumerate(tags):
                coeffs.extend(list(cv[key].dropna()))
            pylab.hist(coeffs, **kargs)
        pylab.grid()
        pylab.xlabel("Coefficient of Variation")
        pylab.ylabel("#")
        return coeffs


    def hist2d_mu_versus_cv(self, bins=100, cmap="hot_r", fontsize=10, Nlevels=4,
                            contour=True,**kargs):
        """plots histogram of mean across replicates versus coefficient variation

        :param int bins: binning for the 2D histogram
        :param fontsize: fontsize for the labels
        :param contour: show some contours
        :param int Nlevels: must be more than 2

        .. plot::
            :include-source:
            :width: 50%

            >>> from msdas import *
            >>> r = replicates.ReplicatesYeast(get_yeast_raw_data())
            >>> r.drop_na_count(54) # to speed up the plot creation
            >>> r.hist2d_mu_versus_cv()

        """
        cv = self.get_coefficient_variation_df().unstack()
        mu = self.get_mu_df().unstack()

        tokeep = (pd.isnull(cv)==False) & (pd.isnull(mu)==False)

        if len(cv)>10000:
            self.logging.info("Computing 2D histogram. Please wait")

        pylab.clf()
        res = pylab.hist2d( cv[tokeep]*100, pylab.log10(mu[tokeep]), bins=bins,
                cmap=cmap, norm=colors.LogNorm())
        pylab.colorbar()

        X, Y = pylab.meshgrid(res[1][0:bins], res[2][0:bins])
        if contour:
            levels = [round(x) for x in pylab.logspace(0, pylab.log10(res[0].max().max()),Nlevels)]
            pylab.contour(X,Y,res[0].transpose(), levels[2:], color="g")
            #pylab.clabel(C, fontsize=fontsize, inline=1)
        pylab.xlabel(r"mean of errors, $E\left[ \epsilon_{i,e,r} \right]_{e,r}$", fontsize=fontsize)
        pylab.ylabel(r"Coefficient of Variation,  $\frac{\sigma_{i,e}}{\mu_{i,e}}$", fontsize=fontsize)

        pylab.grid(True)
        pylab.xlabel("Coefficient Variation (percent)")
        pylab.ylabel("Average (log10)")
        return res

    def hist2d_errors_versus_cv(self, bins=100, cmap="hot_r", fontsize=10,
                                Nlevels=4, contour=True,**kargs):

        r"""Show hist 2D of mean of errors across replicates and experiment versus
        Coefficient variation.

        :param int bins: binning for the 2D histogram
        :param fontsize: fontsize for the labels
        :param contour: show some contours
        :param int Nlevels: must be more than 2



        .. plot::
            :include-source:
            :width: 50%

            >>> from msdas import *
            >>> r = replicates.ReplicatesYeast(get_yeast_raw_data())
            >>> r.drop_na_count(54) # to speed up the plot creation
            >>> r.hist2d_errors_versus_cv()

        """

        cv = self.get_coefficient_variation_df().unstack()
        grand_mu_errors = self.get_mean_errors_df().unstack()
        tokeep = (pd.isnull(cv)==False) & (pd.isnull(grand_mu_errors)==False)



        if len(cv)>10000:
            self.logging.info("Computing 2D histogram. Please wait")
        pylab.clf()
        res = pylab.hist2d(grand_mu_errors[tokeep], cv[tokeep]*100, bins=bins,
                cmap=cmap, norm=colors.LogNorm())
        pylab.colorbar()

        X, Y = pylab.meshgrid(res[1][0:bins], res[2][0:bins])
        if contour:
            levels = [round(x) for x in pylab.logspace(0, pylab.log10(res[0].max().max()),Nlevels)]
            pylab.contour(X,Y,res[0].transpose(), levels[2:], color="g")
            #pylab.clabel(C, fontsize=fontsize, inline=1)
        pylab.xlabel(r"mean of errors, $E\left[ \epsilon_{i,e,r} \right]_{e,r}$", fontsize=fontsize)
        pylab.ylabel(r"Coefficient of Variation,  $\frac{\sigma_{i,e}}{\mu_{i,e}}$", fontsize=fontsize)

        pylab.grid(True)
        return res

    def get_coefficient_variation_df(self):
        r"""Return dataframe containing the coefficient of variation

        .. math:: CV_{i,e} = \frac{\sigma_{i,e} }{\mu_{i,e}}

        """
        cv = pd.DataFrame(self.get_coefficient_variation())
        cv = cv[[x for x in self.measurements.columns if "." not in x]]
        return cv

    def get_coefficient_variation(self):
        r"""Return dictionary containing the coefficient of variation for each measurement

        .. math:: CV_{i,e} = \frac{\sigma_{i,e} }{\mu_{i,e}}

        In the dictionary, keys correspond to the name of a measurement as returned
        by :meth:`get_unique_measurement_name`. Values are time series with the
        coefficient of variation.

        """
        tags = self.get_unique_measurement_name()
        coeffs = {}
        mus = self.get_mus()
        sigmas = self.get_sigmas()
        coeffs = {}
        for i, time in enumerate(tags):
            coeffs[time] = (sigmas[time]/mus[time])
        return coeffs

    def get_standarised_residuals(self):
        r"""Returns dictionary (times being the keys)  with residual errors

        .. math:: \epsilon_{i,e} =  \frac{\left|X_{i,e,r}-\mu_{i,e}\right|}{\sigma_{i,e}}

        """
        tags = self.get_unique_measurement_name()
        errors = {}
        mus = self.get_mus()
        sigmas = self.get_sigmas()

        for i, tag in enumerate(tags):
            df = self.get_replicates_from_one_unique_measurement(tag)
            errors[tag] = np.abs(df.sub(mus[tag], axis="index")).divide(sigmas[tag],
                axis="index")
        return errors

    def _get_na_count_per_experiment(self):
        nas = {}
        for tag in self.get_unique_measurement_name():
           df = self.get_replicates_from_one_unique_measurement(tag)
           R = len(df.columns)
           na = R  - pd.notnull(df).sum(axis=1)
           nas[tag] = na.copy()
        return nas

    def pcolor_na(self, sort_index=False, noxticks=True, cmap="jet"):
        """plots the number of NA versus experiment and protein

        :param sort_index:
        :param noxticks: can be slow if you have thousands of protein names
        :return: dataframe used with NAs per row and experiment


        .. plot::
            :include-source:
            :width: 70%

            from msdas import *
            from easydev import gsf
            filename= gsf("msdas", "data", "YEAST_raw_sample.csv")
            r = Replicates(filename, verbose=False)
            r.pcolor_na()

        """
        df = pd.DataFrame(self._get_na_count_per_experiment()).transpose()
        df.columns = self.df.Identifier

        if sort_index:
            dfna = self.get_na_count()
            dfna.sort(axis=1)
            df = df[dfna.index]
            xtickslabel = self.metadata.ix[dfna.index].Protein + "_" + self.metadata.Psite
        else:
            xtickslabel = self.metadata.Protein + "_" + self.metadata.Psite

        pylab.clf()
        #pcolor in at least 15 times faster than pcolor but is not in pandas, so convert data to matrix before
        pylab.pcolormesh(df.as_matrix(), cmap=cmap)
        pylab.yticks([0.5+x for x in range(0,len(df.index))], df.index)
        if noxticks==False:
            pylab.xticks(range(0, len(self.metadata)), xtickslabel, rotation=90)
        pylab.colorbar()
        pylab.title("Number of NA")
        pylab.ylim([0,len(df.index)])
        pylab.xlim([0,len(self.metadata)])
        return df

    def plot_na_per_experiment(self, percent=False):
        """plot curve showing number of NAs per experiment

        .. plot::
            :include-source:
            :width: 70%

            from msdas import *
            from easydev import gsf
            filename= gsf("msdas", "data", "YEAST_raw_sample.csv")
            r = Replicates(filename, verbose=False)
            r.plot_na_per_experiment()



        """
        N = len(self.df)
        pylab.clf()
        if percent == False:
            count = pd.isnull(self.measurements).sum(axis=0)
            count.plot()
            pylab.ylabel("Number of NAs per experiment")
            pylab.ylim([0, N])
        else:
            count = pd.isnull(self.measurements).sum(axis=0)
            count/=float(N)
            count *= 100
            count.plot()
            pylab.ylabel("Number of NAs (percentage) per experiment")
            pylab.ylim([0,101])

        pylab.xticks(range(0, self.N), self.measurements.columns, rotation=90)



        pylab.tight_layout()


    def plot_mu_sigma(self, tag_list=None, histograms=True, loglog=True, fontsize_legend=10):
        """Plots coefficient variation versus sigma as a scatter plot with histograms

        :param list tag_list: list of measurement to look at (Defaults to all)

        .. plot::
            :include-source:
            :width: 90%

            from msdas import *
            r = replicates.ReplicatesYeast(get_yeast_raw_data())
            r.plot_mu_sigma(tag_list=["a0_t0", "a0_t45"], loglog=True, histograms=True)


        """
        mus = self.get_mus()
        sigmas = self.get_coefficient_variation()

        if tag_list == None:
            tag_list = self.get_unique_measurement_name()

        f = pylab.figure()
        f.clf()
        if histograms is True:
            a1 = f.add_axes([.12,.12,.5,.5])
            a_top = f.add_axes([.12,.7,.5,.2])
            a_right = f.add_axes([.7,.12,.2,.5])
        else:
            a1 = pylab.axes()

        for time in tag_list:
            pylab.sca(a1)
            if loglog:
                a1.loglog(mus[time], sigmas[time], 'o', alpha=0.5, label=time)
            else:
                a1.plot(mus[time], sigmas[time], 'o', alpha=0.5, label=time)
            pylab.grid(True)
            if histograms:
                pylab.sca(a_top)
                if loglog:
                    pylab.log10(mus[time]).hist(alpha=0.5, bins=20)
                    a = pylab.gca()
                    a.set_xlim(pylab.log10(a1.get_xlim()))
                else:
                    mus[time].hist(alpha=0.5, bins=20)
                    a = pylab.gca()
                    a.set_xlim(a1.get_xlim())
                a.grid(True)
                pylab.sca(a_right)
                if loglog:
                    pylab.log10(sigmas[time]).hist(alpha=0.5, bins=20, orientation="horizontal")
                    a = pylab.gca()
                    a.set_ylim(pylab.log10(a1.get_ylim()))
                else:
                    sigmas[time].hist(alpha=0.5, orientation="horizontal", bins=20, xrot=True)
                    a = pylab.gca()
                    a.set_ylim(a1.get_ylim())
                a.grid(True)
        #pylab.legend(loc="upper left", fontsize=fontsize_legend)
        #pylab.tight_layout()
        pylab.xlabel("mu")
        pylab.ylabel("sigma")
        pylab.legend(tag_list, fontsize=fontsize_legend, loc="lower right")

    def reset(self):
        """Reset the dataframe to its original contenst.

        Since loading large data set takes some time to initialise, this
        could be useful for debugging if you cut the dataframe using for
        instance :meth:`drop_na_count`.
        """
        self.df = self._rawdf.copy()
        #self.metadata = self._metadata.copy()



class ReplicatesYeast(Replicates):
    """This class is a specialisation of :class:`Replicates` dedicated to the YEAST data


    The data set comprises 36 experiment that are combination of
    alpha hormone at different times (0,1,5,10,20,45)
    """

    #: times
    times = [0,1,5,10,20,45]
    #: input tags to be found
    yeast_tags = ["Ma0s04", "Ma1s04", "Ma5s04", "Ma10s04", "Ma20s04","Ma45s04"]

    def __init__(self, data, verbose=False, cleanup=True):
        """.. rubric:: constructor

        :param data: a valid input to :class:`msdas.readesr.MassSpecReader`
        :param cleanup: see :class:`msdas.readesr.MassSpecReader` for details.


        """
        super(ReplicatesYeast, self).__init__(data, verbose=verbose, cleanup=cleanup)

        # renaming columns
        for i, tag in enumerate(self.yeast_tags):
            self.df.columns = [c if c.startswith(tag)==False
                else "a" + str(self.times[i]) + "_t" + c.split("_")[1] for c in self.df.columns]

        self._rawdf = self.df.copy()

    def copy(self):
        """Return a copy of Replicates instance. """
        r = ReplicatesYeast(self)
        return r

    def normalise(self):

        from easydev import gsf
        self.tic = pd.read_csv(gsf("msdas", "data", "YEAST_TIC.csv"))
        # FIXME we assume that replicates are sorted and rename them to differentiate their name
        # and to agree with data. Should therefore be ordered in the same way.
        self.tic.columns = ['name', 'TIC', 'dummy']
        def func_rename(x,i):
            if i%3==0:
                return x
            elif i%3==1:
                return x+".1"
            elif i%3==2:
                return x+".2"
                # rename to add the replicate index
        self.tic['name'] = [func_rename(x,i) for i,x in enumerate(self.tic.name)]


        # rename to replcae Ma10s04_20 into convention a10_t20
        def func_rename2(x):
            res = x.split(".")
            name = res[0]
            if len(res) == 1:
                replicate = None
            else:
                replicate = res[1]
            hormone, salt = name.split("_")
            hormone = hormone.replace("Ma", "")
            hormone = hormone.replace("s04", "")
            if replicate:
                return "a" + hormone + "_t" + salt + "." + replicate
            else:
                return "a" + hormone + "_t" + salt
        self.tic['name'] = [func_rename2(x) for x in self.tic['name']]
        # we do not need these informatoion anymore
        self.tic.drop(['dummy'], axis=1, inplace=True)

        # let us have the name as column names to be consistent with the data
        self.tic = self.tic.set_index("name").transpose()

        # let us reorder the columns similarly to the data
        self.tic = self.tic[self.measurements.columns]

        # FIXME there is probably a way to divide enture matrix by the tic row
        # but could not figure it out so let us do it column by column
        for col in self.measurements.columns:
            self.df[col] /= self.tic[col].values[0]


    # CLAUDIA this is an attemp to build scaling function using as template scale function from clustering.py
    def scale(self):

        # the scale function uses unbiased variance
        for col in self.measurements.columns:
            X = self.df[col]
            self.df[col]=(X-X.mean())/X.std(ddof=1)


     # CLAUDIA this is an attemp to build log2 function
    def log_b2(self):

        X=self.df[self.measurements.columns]
        self.df[self.measurements.columns]=np.log2(X)


    def plot_timeseries_midas(self, identifier):

         measurements = [u'a0_t0', u'a0_t1', u'a0_t5', u'a0_t10', u'a0_t20', u'a0_t45',
         u'a1_t1', u'a5_t5', u'a10_t10', u'a20_t20', u'a45_t45',
         u'a1_t0', u'a5_t0', u'a10_t0', u'a20_t0', u'a45_t0']

         l = [0,1,5,10,20,45,46,50,55,65,90,91,95,100,120,145]
         data= self.get_mu_df()[measurements]

         index = self.df[self.df.Identifier == identifier].index[0]
         data = data.ix[index]
         print data
         pylab.plot(l, data, 'ro-')

    def get_na_per_experiment_after_averaging(self):
        """average data and figure out number of NA per alpha experiment"""
        mu = self.get_mu_df()

        na = pd.DataFrame()
        for exp in ['a0', 'a1', 'a5', 'a10', 'a20', 'a45']:
            sum_na = pd.isnull(mu[[x for x in mu.columns if x.startswith(exp+"_")]]).sum(axis=1)
            na[exp] = sum_na
        return na

    def get_na_per_salt_experiment_after_averaging(self):
        """average data and figure out number of NA per alpha experiment"""
        mu = self.get_mu_df()

        na = pd.DataFrame()
        for exp in ['t0', 't1', 't5', 't10', 't20', 't45']:
            sum_na = pd.isnull(mu[[x for x in mu.columns if x.endswith(exp)]]).sum(axis=1)
            na[exp] = sum_na
        return na

    #def is_na_count_exceeds_minnonzero_in_one_experiment(self, min_non_zero=4):
    #    na = self.get_na_per_experiment_after_averaging()
    #    return na.max(axis=1) >= min_non_zero

    def is_na_count_exceeds_minnonzero_in_all_experiment(self, min_non_zero=4):
        na = self.get_na_per_experiment_after_averaging()

        N = len(self.get_unique_measurement_name())

        # <= like in Stephania's code so that a row with only 4 non-na values are removed.
        ts =  N - na.sum(axis=1) <=  min_non_zero
        return ts

    def drop_na_exceeds_minnonzero(self, min_non_zero=4):
        """drop rows where Number of number of non-null values is below 4

        .. note:: Not used anymore. Maybe dropped.

        """
        ts = self.is_na_count_exceeds_minnonzero_in_all_experiment(min_non_zero=min_non_zero)
        self.logging.info("Dropping %s rows that have less than %s non-na values" % (sum(ts),min_non_zero))
        self.df.drop(ts[ts].index, inplace=True)

    def clean_na_paper(self, inplace=True,N=2):
        """To be used with YEAST data only

        :param N: max number of NAs per experiment

        The YEAST data set is made of 6 experiments of 6 data points times 3 replicates.
        We first average the data over replicates so we have 6 experiments of 6 data points.
        In each experiments, if there have less than 6 - N values, then that experiment is reset:
        all data (including replicates) are set to NAs.

        """
        na = self.get_na_per_experiment_after_averaging()
        filterna = na>N # means that there less than 4 measurements

        labels = list(na.columns)
        # Let us exoande the na dataframe by copying each column (e.g. a0)
        # into expected replicates e.f. (a0_t0, a0_t0.1, ...a0_t45)
        for column in filterna.columns:
            expanded_columns = [x for x in self.measurements.columns if x.startswith(column+"_")]
            #col1 = [this+".1" for this in expanded_columns]
            #col2 = [this+".2" for this in expanded_columns]
            #expanded_columns.extend(col1)
            #expanded_columns.extend(col2)
            for x in expanded_columns:
                filterna[x] = filterna[column]

        # let us remove the labels a0,a1,a5,a10,a20,a45
        filterna.drop(labels, axis=1, inplace=True)

        for this in self.metadata.columns:
            filterna[this] = [False] * len(filterna)
        self._filterna = filterna

        if inplace:
            self.df = self.df.mask(filterna)
        else:
            df = self.df.mask(filterna)
            return df

    def clean_na_paper2(self, inplace=True, N=2, M=2):
        """To be used with YEAST data only

        Identify rows with at least 2 experiments (M=2) with more than 6-N values
        and crosses the 6 x replicates measurements in the experiment.
        """
        na = self.get_na_per_salt_experiment_after_averaging()
        filterna = na>N # means that there less than 4 measurements

        labels = list(na.columns)
        # Let us exoande the na dataframe by copying each column (e.g. a0)
        # into expected replicates e.f. (a0_t0, a0_t0.1, ...a0_t45)
        for column in filterna.columns:
            expanded_columns = [x for x in self.measurements.columns if x.split(".")[0].endswith(column)]
            for x in expanded_columns:
                filterna[x] = filterna[column]

        # let us remove the labels a0,a1,a5,a10,a20,a45
        filterna.drop(labels, axis=1, inplace=True)

        for this in self.metadata.columns:
            filterna[this] = [False] * len(filterna)
        self._filterna = filterna

        if inplace:
            self.df = self.df.mask(filterna)
        else:
            df = self.df.mask(filterna)
            return df

