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
"""Clustering utilities"""
import pylab
import numpy as np
import scipy.cluster
import scipy

from msdas import readers

__all__ = ["MSClustering", "Affinity"]



class MSClustering(readers.MassSpecReader):
    """Clustering utilities

    This class is used to perform clustering analysis on measurements to
    be found in a :class:`msdas.readers.MassSpecReader` instance.

    The data extraction consists in selecting a subset of the initial dataframe
    provided by the Reader classes. The selection is made on (1) the protein
    name that can be given by the user or not, and (2) the time measurements are
    kept only. This is achieved via the :meth:`get_group`.

    Currently, only the Affinity propagation algorithm is provided. It can
    be used via the :meth:`run_affinity_propagation_clustering` method or
    directly with the :class:`Affinity` class.

    Here is an example that reads data, creates a MSClustering
    instance and finally performs a clustering with the affinity propagation
    algorithm.

    ::

        from msdas import *
        y = MassSpecReader(get_yeast_small_data(), mode="YEAST")
        c = MSClustering(y)
        c.run_affinity_propagation_clustering(preference=-120)


    """
    def __init__(self, data, mode="default", fillna=True, verbose=True, cleanup=False):
        """.. rubric:: constructor


        :param data: can be either a dataframe obtained from
            :class:`~msdas.readers.MassSpecReader` or an instance of
            :class:`~msdas.readers.MassSpecReader`
            itself.
        :param str mode: must be **yeast**. 
        """
        super(MSClustering,self).__init__(data, verbose=verbose, mode=mode, cleanup=cleanup)
        self._metadata_names.append("cluster")
        if fillna:
            self.logging.info("Filling Nas with zeros if any")
            self.df.fillna(0, inplace=True)

    def _get_groups(self):
        return self.df.groupby("Protein").groups
    groups = property(_get_groups,
            doc="Returns all group for each protein. Alias to dataframe.groupby('Protein')")

    def get_group(self, name=None, tag="onlyCD3"):
        """Returns a subset of the original dataframe (stored in :attr:`df`)

        If name is provided, only the protein corresponding to that name is
        selected, otherwise all proteins are kept. Then, only columns related to
        time measurements are kept (ie. protein columns and non numeric columns
        are ignored).

        In the **yeast** case, the criteria to select time measurements is to
        keep columns that starts with the letter **a**.

        In the tcell case, "Nor_Unstimulated" column is selected as as as
        columns starting with "Nor" and containin the tag provided(e.g., onlyCD3)
        indices on time are renamed as 00, 05, 10, 30.
        """
        if name in self.groups.keys():
            indices = self.groups[name]
            subdf = self.df.ix[indices]
            subdf = subdf[[x for x in self.df.columns if x in self.measurements.columns]]
            subdf = subdf.transpose()
            subdf.drop_duplicates(inplace=True)
            subdf.columns = self.df.Identifier[indices]
            return subdf
        elif name==None:

            subdf = self.df[[x for x in self.df.columns if x in self.measurements.columns]]
            subdf = subdf.transpose()
            subdf.drop_duplicates(inplace=True)
            subdf.columns = self.df.Identifier
            return subdf
        else:
            raise KeyError("name not found")


    def scale(self, df):
        """scale a df and returns the scaled version.

        This can be used on the numeric dataframe returned by :meth:`get_group`.

        ::

            from msdas import *
            df = MassSpecReader(get_yeast_small_data())
            c = MSClustering(df)
            subdf = c.get_group()
            scaled_subdf = c.scale(subdf)


        .. math:: \hat{X} = {X-\mu}{\sigma (X)}

        where the standard deviation is the biased one computed with
        numpy as X.std(ddof=1).

        """
        # the scale function uses biased variance
        #from sklearn.preprocessing import scale
        #newdf = pd.DataFrame(scale(df, axis=axis, with_mean=with_mean,
        #    with_std=with_std), index=df.index, columns=df.columns)
        print("rescaling over %s columns" % len(df.columns))
        newdf = df.copy()
        for i in newdf.columns:
            X = newdf.ix[:,i]
            newdf.ix[:,i] = (X-X.mean())/X.std(ddof=1)

        return newdf

    def _normalise(self, subdf):
        """in place normalisation"""
        for this in subdf.columns:
            subdf.ix[:,this] /= pylab.sqrt(subdf.ix[:,this].dot(subdf.ix[:,this]))
        return subdf

    def _get_euclidean_distance_matrix(self, name=None):
        """

        .. plot::
            :include-source:
            :width: 50%

            clf()
            pcolor(c.get_dot_matrix("DIG1"))
            colorbar()


        .. seealso:: used by :meth:`dendogram`

        .. note:: The data is normalised before calling scipy.spatial.distance.eucllidean.
        may not be required actually.

        """
        import scipy.spatial.distance
        subdf = self.get_group(name)
        self._normalise(subdf)

        N = len(subdf.columns)

        m = np.zeros((N, N))
        for i, ix in enumerate(subdf.columns):
            for j, jx in enumerate(subdf.columns):
                m[i,j] = scipy.spatial.distance.euclidean(subdf[ix], subdf[jx])
        return m

    def _get_dot_matrix(self, name=None):
        """

        ::
            clf()
            pcolor(c.get_dot_matrix("DIG1"))
            colorbar()

        .. seealso:: used by :meth:`dendogram`

        .. note:: The data is normalised before calling scipy.spatial.distance.eucllidean.
            may not be required actually.
        """
        subdf = self.get_group(name)
        self._normalise(subdf)
        N = len(subdf.columns)
        m = np.zeros((N,N))
        for i,ix in enumerate(subdf.columns):
            for j,jx in enumerate(subdf.columns):
                m[i,j] = subdf.ix[:,ix].dot(subdf.ix[:,jx])
        return m

    def plot_vectors(self, name, normalise=False, scale=True, legend=True,
                     subplots=False, fontsize=8, **kargs):
        """Plot time series or data of the time measurements for a group of
        protein (i.e. individual peptides).

        :param name: must be a valid protein name (see :attr:`df.Protein`) or
            None to select all proteins.
        :param normalise, scale: mutually exclusive pre processing option to
            normalise or scale the data before plotting.
        :param bool legend: add the legend or not
        :param fontsize: used for the legend
        :param kargs: any option accepted y the plot method used by the
            pandas dataframe. location and fontsize are used for the legend

        .. plot::
            :include-source:
            :width: 50%

            >>> from msdas import *
            >>> y = MassSpecReader(get_yeast_small_data())
            >>> c = MSClustering(y)
            >>> # c.run_affinity_propagation_clustering(-90)
            >>> c.plot_vectors("DIG1", normalise=False, scale=True, legend=False)


        """
        subdf = self.get_group(name)
        assert (normalise == True and scale == False) or (normalise == False and
                scale == True) or (normalise == False and scale == False)

        if normalise:
            self._normalise(subdf)
        if scale:
            subdf = self.scale(subdf)


        if self.mode == "YEAST":
            kargs['xticks'] = kargs.get("xticks", range(0, self.N))
            kargs['rot'] = kargs.get("rot", 90)
        for this in subdf.columns:
            subdf.ix[:, this].plot(legend=False, subplots=subplots, **kargs)

        if legend:
            ax = pylab.gca()
            patches, labels = ax.get_legend_handles_labels()
            ax.legend(patches,labels, fontsize=kargs.get("fontsize",fontsize),
                      loc=kargs.get('location',"best"))

    def dendogram(self, name=None, method="euclidean", metric="euclidean",p=2,
                  pdist=False, rotation=90, **kargs):
        """Dendogram. In progress, please do not use or with care.

        :param name: name of the protein to perform dendogram on. IF none, all proteins are selected
        :param method: euclidean or dot
        :param pdist: pdist computes the Euclidean distance between pairs of
            objects in m-by-n data matrix X. Rows of X correspond to observations,
            and columns correspond to variables.
        :param rotation: rotation of xlabels (protein name + psite)
        :param kargs: any valid parameter of the scipy dendogram function

        .. plot::
            :include-source:
            :width: 80%

            import pylab
            from msdas import *
            y = MassSpecReader(get_yeast_small_data())
            c = MSClustering(y, mode="YEAST")
            c.dendogram("DIG2", pdist="euclidean")
            # equivalent to c.dendogram("DIG2", method="euclidean")
            pylab.figure(2)
            c.plot_vectors("DIG2")

        .. note:: pdist has many metrics and we recommand to use the pdist argument
            in place of method parameter, which maybe deprecated in the future.

        """
        pylab.clf()

        if method == "euclidean":
            X = self._get_euclidean_distance_matrix(name)
            names = list(self.get_group(name).columns)
            linkage = scipy.cluster.hierarchy.linkage(X)
        elif method == "dot":
            X = self._get_dot_matrix(name)
            names = list(self.get_group(name).columns)
            linkage = scipy.cluster.hierarchy.linkage(X)
        elif pdist == True:

            X = scipy.spatial.distance.pdist(self.get_group(name).transpose(), metric=metric,p=p).transpose()
            subdf = list(self.get_group(name).columns)
            names = list(subdf.columns)
            linkage = scipy.cluster.hierarchy.linkage(pdist)
        else:
            raise ValueError("method must be euclidean or dot or pdist must be provided")

        # note that
        # linkage(pdist(c.get_group("DIG2").transpose()).transpose())
        # is equivalent to
        # linkage(X)

        scipy.cluster.hierarchy.dendrogram(linkage, labels=names,
                                               leaf_rotation=rotation, **kargs)
        pylab.title(name)
        pylab.tight_layout()

        self._linkage = linkage

    def run_affinity_propagation_clustering(self, protein=None, preference=-30):
        """Run the affinity propagation algorithm

        :param str name: a protein name. If None, consider all proteins.
            (Defaults to None)
        :param float preference: the affinity propagation preference parameter.
            should be negative. Large negative values gives less clusters.
        :return: an affinity instance that can be manipulated to extract labels,
            cluster, or plot the clusters. See :class:`Affinity` for more
            details

        The original data stored in :attr:`df` will contain anew column
        called `cluster` will the label of each cluster.

        .. warning:: NA are filled with zero for now.
        """
        data = self.get_group(protein)
        data.fillna(0, inplace=True)
        affinity = Affinity(data, method="euclidean",
                preference=preference, transpose=True)

        # clustering
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # this raises annoying warnings
            res = self.df.query("Identifier in @data.columns", engine="python")

        indices = res.index

        if "cluster" not in self.df.columns:
            self.df['cluster'] = self.df.shape[0] * [None]


        self.affinity_results = affinity
        self.df['cluster'].ix[indices] = affinity.labels
        return affinity

    def plotClusters_and_protein(self, name, **kargs):
        """Alias to :meth:`Affinity.plotClusters_and_protein`"""
        self.affinity_results.plotClusters_psites(name, **kargs)

    def plotClusters(self, **kargs):
        """Alias to :meth:`Affinity.plotClusters`"""
        self.affinity_results.plotClusters(**kargs)


    def plotClusteredProtein(self, name, scale=True, fontsize=8, **kargs):
        """For a given protein, plot time series by group/label


        If a group does not have the protein in it, then nothing is plotted.

        .. note:: the affinity propagation algorithm must have been run before.
        """



        groups = self.df.groupby(["Protein", "cluster"]).groups
        #groups = [g for g in groups if g[0]==name]
        #columns = [c for c in self.df.columns if c.startswith("a")]
        for group in groups.iteritems():
            if group[0][0] == name:
                index = self.metadata.ix[group[1]].Identifier
                data = self.measurements.ix[group[1]].transpose()
                data.columns = list(index)
                if scale:
                    data = self.scale(data)
                data.plot()
                pylab.title("Group %s" % group[0][1])

                #self.df.ix[group[1], columns].transpose().plot()
        pylab.legend(fontsize=fontsize)





class Affinity(object):
    """Apply affinity propagation clustering on MS data set.

    Get some data first. It must be a dataframe that contains the
    measurements as columns and the Psites as indices. It can be the
    output of :meth:`MSClustering.get_group` method.

    .. plot::
        :include-source:
        :width: 30%

        >>> from msdas import *
        >>> y = MassSpecReader(get_yeast_small_data(), mode="YEAST")
        >>> c = clustering.MSClustering(y)
        >>> af = c.run_affinity_propagation_clustering(preference=-90)

        >>> # OR
        >>> data = c.get_group()
        >>> af = Affinity(data, preference=-90, method="euclidean")

        >>> af.plotClusters()

        >>> af.nclusters
        7
        >>> labels = af.labels
        >>> center_indices = af.cluster_centers_indices

    The data will be scaled using the
    :meth:`~msdas.clustering.MSClustering.scale` method.

    Similarly for the cell except that the method should be correlation instead
    of euclidean. In which case, the correlation is performed internally.

    """
    def __init__(self, data, method, affinity_matrix=None, preference=-30,
            damping=0.9, transpose=False, verbose=True):
        """.. rubric:: contructor

        :param DataFrame data:
        :param method: either euclidean or correlation.
        :param affinity_matrix:
        :param preference: the affinity propagation parameter that affects the number
            of clusters. (default is -30)
        :param float damping: values between 0.5 and 1 (stricly). (default is
            0.9 to be similar to what is used in R affinity propagation routine)
        :param bool transpose: transpose the affinity matrix before calling the
            affinity propagation algorith. Should be used if euclidean distance
            method is selected.


        .. note:: the dataframe is scaled internally using
            :math:`\hat{X} = (X-\mu)/\sigma`

        """
        self.verbose = verbose
        self.df = data.copy()

        # TODO what about another scaling methods ?
        self.df = self._scale(self.df)
        #self.df = self._normalise(self.df)
        self.method = method
        if affinity_matrix!=None:
            self.method = "other"
            self.affinity_matrix = affinity_matrix
        elif self.method == "correlation":
            self.affinity_matrix = self.df.corr()
        elif self.method=="euclidean":
            self.affinity_matrix = self.df.copy()
        else:
            raise ValueError("method incorrect or affinity_matrix not provided")

        self.affinity_results = None
        self._preference = preference
        self._damping = damping

        ## TODO this has to be done in yeast case
        if transpose==True:
            self.affinity_matrix = self.affinity_matrix.transpose()
        self.AffinityPropagation()

    def _scale(self, df):
        if self.verbose:
            print("rescaling over %s columns" % len(df.columns))
        newdf = df.copy()
        for i in newdf.columns:
            X = newdf.ix[:,i]
            newdf.ix[:,i] = (X-X.mean())/X.std(ddof=1)
        return newdf

    def _normalise(self, df):
        """in place normalisation"""
        newdf = df.copy()
        for this in df.columns:
            newdf.ix[:,this] /= np.sqrt(newdf.ix[:,this].dot(newdf.ix[:,this]))
        return newdf

    def AffinityPropagation(self, preference=None, damping=None, max_iter=200):
        """The core of the Affinity Propagation algorithm.

        This method uses scikit-learn

        :param preference: the affinity propagation parameter that affects the number
            of clusters. (default is -30)
        :param float damping: values between 0.5 and 1 (stricly). (default is
            0.9 to be similar to what is used in R affinity propagation routine)
        :param int max_iter: 200 by default


        Populates :attr:`affinity_results` with the results of the clustering.

        """
        if preference == None:
            preference = self._preference
        if damping == None:
            damping = self._damping

        if self.verbose:
            print("Running affinity propagation using preference=%s" % preference)
        from sklearn.cluster import AffinityPropagation
        if self.method != "euclidean":
            self.affinity_results = AffinityPropagation(preference=preference, max_iter=max_iter,
                    damping=damping, affinity="precomputed").fit(self.affinity_matrix)
        else:
            self.affinity_results = AffinityPropagation(preference=preference, max_iter=max_iter,
                    damping=damping).fit(self.affinity_matrix)
        if self.verbose:
            print("Found %s clusters" % self.nclusters)
        return self.affinity_results

    def plotClusters(self, legend=True, fontsize=8, **kargs):
        """Plots N figures related to the N clusters found with all time series

        See class documentation for an example


        .. seealso:: :class:`~msdas.clustering.MSClustering`


        All timeseries for a given cluster are plotted in blue.
        The exampler in red.

        """
        if self.affinity_results == None:
            raise ValueError("Call AffinityPropagation first.")
        else:
            for icluster, cluster in enumerate(self.affinity_results.cluster_centers_indices_):
                # plot all time series from this cluster
                self.df.ix[:,self.affinity_results.labels_==icluster].plot(legend=legend,
                        color="b")
                # and the exampler:
                self.df.ix[:,cluster].plot(legend=legend, color="r", lw=4)
                if legend:
                    pylab.legend(fontsize=fontsize)

    def plotClusters_and_protein(self, name, fontsize=8, legend=True):
        """Same as `plotClusters` but superimposed a specific protein if present

        All timeseries for a given cluster are plotted in blue.
        The exampler is in red and the time series matching the name provided are
        in green.


        .. plot::
            :include-source:
            :width: 30%

            >>> from msdas import *
            >>> y = MassSpecReader(get_yeast_small_data(), mode="yeast")
            >>> c = clustering.MSClustering(y)
            >>> af = c.run_affinity_propagation_clustering(preference=-90)
            >>> af.plotClusters_and_protein("DIG1", legend=True)

        """
        if self.affinity_results == None:
            raise ValueError("Call AffinityPropagation first.")

        for label in set(self.labels):
            psites = self.df.columns[np.logical_and(self.labels==label,
                [x.startswith(name) for x in self.df.columns])]
            cluster_index = self.cluster_centers_indices[label]

            # all timeseries
            self.df.ix[:,self.labels==label].plot(legend=legend, color="b")

            # the exampler
            self.df.ix[:,cluster_index].plot(legend=legend, color="r", lw=4)

            # the psites
            for psite in psites:
                self.df.ix[:,psite].plot(legend=legend, color="g", lw=2)
            if legend:
                pylab.legend(fontsize=fontsize)


    def _get_nclusters(self):
        return len(self.affinity_results.cluster_centers_indices_)
    nclusters = property(_get_nclusters)

    def _get_preference(self):
        return self.affinity_results.preference
    preference = property(_get_preference)

    def _get_cluster_center_indices(self):
        return self.affinity_results.cluster_centers_indices_
    cluster_centers_indices = property(_get_cluster_center_indices,
            doc="Returns indices of each cluster center")

    def _get_labels(self):
        return self.affinity_results.labels_
    labels = property(_get_labels, doc="Returns labels of the clustering")

    def __str__(self):
        str_ = "Number of clusters found %s" % self.nclusters
        str_+="\n"
        for i in range(0,self.nclusters):
            str_+="Cluster %s is made of %s items" % (i+1,
                    sum(self.affinity_results.labels_==i))
            str_+="\n"
        return str_

