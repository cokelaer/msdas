from __future__  import division

import clustering
import numpy as np
import pylab
import pandas as pd
from msdas import replicates
from readers import MassSpecReader
from easydev import get_share_file as gsf
from cno import CNOGraph

__all__ = ["YEAST", "get_yeast_filenames", "YEAST2MIDAS",
           "get_yeast_raw_data", "get_yeast_small_data"]

"""

    sim = plotLBodeFitness(cnolist, pknmodel,res, plotParams=list(F=.3, cmap_scale=2, cex=0.5))

    Ntimes = length(sim)
    N = nrow(sim[[1]]) * ncol(sim[[1]]) * Ntimes
    signals = cnolist@signal
    sim2 = list()
    for(i in seq_along(sim)){
         sim2[[i]] = matrix(colMeans((signals[[i]]-sim[[i]])^2, na.rm=T), nrow=1)
    }
    sim2 = matrix(unlist(sim2), nrow=Ntimes, byrow=T)
    sim2 = colMeans(sim2, na.rm=T)
    d = data.frame(sim2, row.names=colnames(cnolist@signals[[1]]))
    write.table(d, col.names=F, quote=F, sep=",", file="test.csv")


"""


class YEAST(object):
    times = [0,1,5,10,20,45]
    def __init__(self):
        pass


def get_yeast_small_data():
    """Return filename of data related to YEAST

    This is the CSV file resulting from the alignment of 6 small data files
    provided in the package. Here is the procedure followed to create that file.

    ::

        from msdas import *
        m = MassSpecAlignmentYeast(get_yeast_small_files())
        # data is now aligned, we can save it
        # We also wanted the annotations:
        a = AnnotationsYeast(m) # creates the identifiers
        a.get_uniprot_entries() # Creates the entry
        a.set_annotations()     # Creates Entry names
        a.check_entries_versus_sequence()
        # columns are rename with prefix "a"
        a.to_csv("YEAST_small_all_test.csv")


    """
    return gsf("msdas", "data", "YEAST_small_all.csv")

def get_yeast_raw_data():
    """Return filename to the YEAST raw data

    See :func:`get_yeast_small_data` for procedure on how to generate this
    data file. The input filenames used are :meth:`get_yeast_filenames`
    """
    return gsf("msdas", "data", "Yeast_all_raw.csv")

def get_yeast_filenames(mode="subset"):
    """Returns the list of filenames related to the YEAST data

    :param str mode: valid values are *subset* or *all*. If subset,
        filenames looked for are files called alphaX.csv . If mode
        is set to *all*,  filenames looked for all files
        called MaX_annotated_CLEAN_COMPACT_Normedian_FINAL.csv
    :return: list of 6 full path names related to the YEAST data set (6 files)
        The corresponding files are sorted by times (ie, 0, 1, 5, 10, 20, 45).

    The 6 filenames must be provided as input to the :class:`~msdas.alignment.MassSpecAlignmentYeast`.

    .. seealso:: :class:`~msdas.alignment.MassSpecAlignmentYeast`.
    """
    times = [0, 1, 5, 10, 20, 45]
    if mode == "subset":
        filenames = [gsf("msdas", "data","alpha%s.csv" % str(x)) for x in times]
        return filenames
    elif mode == "all":
        prefix = "Ma"
        suffix = "s04_annotated_CLEAN_COMPACT_NORMmedian_FINAL.csv"
        filenames = [gsf("msdas", "data","{}{}{}".format(prefix,str(x),suffix)) for x in times]
        return filenames
    else:
        raise ValueError("mode must be in [subset, all]. Defaults to *subset*")


class YEAST2MIDAS(MassSpecReader, YEAST):
    """Class pipeline to read/cluster/write yeast data into MIDAS/PKN.

    The constructor of this class performs the following tasks:

        #. Reads small data set (with the interpolated data)
        #. Reads the raw data set (with replicates)
        #. Cleanup and merging similar to what was done in the small dat set (see below for details)
        #. Creates an instance of :class:`msdas.clustering.MSClustering`

    Tools are provided to

        #. Reads a PKN and build a MIDAS from the data. A new PKN is generated
           with new names containing the psites.
        #. plot time series
        #. perform clustering

    YEAST2MIDAS is also a :class:`msdas.readers.MassSpecReader` so you can get
    lots of plotting tools in addition to those provided within this class.

    Here is an example::

        >>> from msdas import yeast
        >>> y = yeast.YEAST2MIDAS(get_yeast_small_data(), get_yeast_raw_data())
        >>> # get a new data frame containing only exemplars
        >>> df = y.get_df_exemplars(preference=-30)
        >>> c = y.get_expanded_cnograph("PKN.sif", df)

        >>> # create the new MIDAS and a CNOGraph for the new PKN
        >>> c, m, e = y.export_pkn_and_midas("../share/data/PKN-yeast.sif")

    :Details: The dataframe :attr:`df` contains all measurements and extra
        information such as the protein/psites/sequence/uniprot_entry. The
        measurements itself is made of a subset of the 36 measurements, which
        are combination of 6 times and 6 alpha experiments.

        You can imagine the data as a matrix 6 by 6. The measurements that are
        kept for MIDAS are made of the first row, first column and main diagonal.



    The raw data is preprocessed using the :class:`msdas.replicates.ReplicatesYeast`
    and follows those steps:

    #. drop rows where peptide contains an oxidation :meth:`msdas.readers.Cleaner.drop_oxidation`
    #. Merge peptide with identical peptide sequence (but different psites)
    #. clean rows with too much NAs using :meth:`~msdas.replicates.ReplicatesYeast.clean_na_paper`
       and :meth:`~msdas.replicates.ReplicatesYeast.clean_na_paper`
    #. Drop rows with only NAs
    #. normalise the data using the TIC values
    #. Call :meth:`~msdas.readers.Cleaner.fix_duplicated_identifier` (rename
       duplicated identifiers)


    """
    def __init__(self, data=None, rawdata=None, verbose=False) :
        """.. rubric:: Constructor

       :param data: optional data set (default to :func:`get_yeast_small_data`)
       :param rawdata: optional data set (default to :func:`get_yeast_raw_data`)

        """
        if data is None:
            data = get_yeast_small_data()
        if rawdata is None:
            rawdata = get_yeast_raw_data()
        super(YEAST2MIDAS, self).__init__(data, cleanup=True,
            merge_peptides=True, verbose=verbose)

        self.logging.info("Reading raw data set")
        self.replicates = replicates.ReplicatesYeast(rawdata)
        # because of the next drop_na_exceeds function.

        # one way to drop NA: 4 values required over the 36 cells
        #self.replicates.drop_na_exceeds_minnonzero()

        # another way, which looks like what was done in the paper is
        # to set to drop rows where one experiment does not have  values
        self.replicates.drop_oxidation()  # to match small data set
        self.replicates.merge_identical_peptide_sequence()
        #self.replicates._rebuild_identifier()


        self.replicates.clean_na_paper() # cross (replaces with NA) alpha experiment
        # if not enough values (at least 4 per set of 6 experiment in alpha direction)
        self.replicates.clean_na_paper2() # cross (replaces with NA) alpha experiment
        # if not enough values (at least 4 per set of 6 experiment in NaCl direction
        self.replicates.drop_na_count(108) # to get closer to small data set

        self.replicates.normalise() # can be before or after
        #self.replicates.merge_peptides()
        #self.replicates.drop_oxidation()  # to match small data set
        self.replicates.sort_psites_ors_only()  # ???
        self.replicates.fix_duplicated_identifier()

        # called once for all
        self._cv_buf = pd.DataFrame(self.replicates.get_coefficient_variation())


        # here, we set na to zero, just for the clustering
        self.cluster = clustering.MSClustering(self, fillna=True, verbose=False)

        # FIXME: yeast.df is a reference so if changed, it affects the user parameter

        # select only data with the 3 conditions we are interested in

        self._measures = list(self.measurements.columns)
        self._measures_salt = ["a"+str(a)+"_t"+str(salt) for salt in self.times for a in self.times]
        # probably not needed since done in the instanciation above but good check

        #self._drop_psites_with_bad_quality()
        #self._drop_visual_clustering()

        #self.run_clustering()

        #print("Scaled dataframe as in the clustering is available in  self.df_scaled")
        #self.df_scaled = self.cluster.scale(self.cluster.get_group()).transpose()
        self.logging.info("data is in `df` attribute and full raw data with replicates in `replicates` attribute")
        self._cleanup = False

        self.mapping = self.mapping_small_to_raw_june()

    def _get_groups(self):
        return self.cluster.groups
    groups = property(_get_groups,
        doc="get the dataframe indices grouped by proteins (alias to cluster.groups)")

    def cleanup_june(self):
        """Cleanup the raw and small data set.

        The small data set contains labels that are incorrect and do not represent
        what pre-processing was performed.

        Besides, interpolation (and extrapolation) were applied.

        In order to replace the interpolated values by NA, we need to know
        the correct label. See :meth:`mapping_small_to_raw_june`. This method
        replaces the interpolated values by NAs in the :attr:`df` dataframe, relabel
        corerctly the identifiers and removes rows in the big dataframe :attr:`replicates`
        that are now useless

        FUS3_S177^T180+S177^Y182+T180^Y182', 'HOG1_T174^Y176' time zero are missing
        for now, we keep the data with interpolated values.


        #. FAR1_S114 has no t0 values. Downstream of PKN. Removed for now
        #. 3 RCK2 have no t0 values. Downstream of PKN. Removed for now

        """

        if self._cleanup == True:
            print("already cleanup. nothing to do")

        mapping = self.mapping_small_to_raw_june()

        values = [x for x in set(mapping.values()) if x]

        # should replace HOG1 t0 a 0.0001 because there is 1 missing value, which
        # is the only missing value.
        #index = self.replicates['HOG1', "T174^Y176"].index[0]
        #self.replicates.df.ix[index, 'a0_t0'] = 0.0001

        # Selection of the data in the replicates that we wawnt to look at
        selection = self.replicates.df.Identifier.apply(lambda x: x in values)
        print(sum(selection))
        self.replicates.df = self.replicates.df[selection]

        # replace the data with the average of the replicates so that now we have the NAs
        # first, let us get the average, keeping track of the identifier as the index of the df
        mu = self.replicates.get_mu_df()
        mu.index = self.replicates.df.Identifier

        # now let us drop identifier that are not in the mapping (keys)
        keys = [k for k,v in mapping.iteritems() if v!=None]
        # the tokeep should be be touched:
        self.tokeep = self.df.copy()

        self.df = self.df[self.df.Identifier.apply(lambda x : x in keys)]

        # replace identifiers by the correct ones
        self.df.Identifier = [mapping[x] for x in self.df.Identifier]
        self.df.Psite = [x.split("_")[1] for x in self.df.Identifier]
        # not need to redo protein column, which should already be correct.

        #replace the data
        for identifier in mu.index:
            index = self.df[self.df.Identifier == identifier].index[0]
            self.df.ix[index, mu.columns] = mu.ix[identifier].values


        index = self.tokeep.query("Identifier == 'HOG1_T174^Y176'").index[0]
        self.df.ix[index] = self.tokeep.ix[index]

        index = self.tokeep.query("Identifier == 'FUS3_S177+T180^Y182'").index[0]
        self.df.ix[index] = self.tokeep.ix[index]
        self.df.ix[index, 'Identifier'] = 'FUS3_S177^T180+S177^Y182+T180^Y182'
        self.df.ix[index, 'Psite'] = 'S177^T180+S177^Y182+T180^Y182'
        self._cleanup == True

        self.cluster = clustering.MSClustering(self, fillna=False, cleanup=False)

    def cleanup_june2(self):
        # 'GPD1'  'S23+S24+S27', 'S23+S24+S25+S27', 'S24^S27'
        #can be clustered. We keep S14^S27
        print("WARNING: should be fixed to not use indices but names")
        self.df.drop([15, 16], axis=0, inplace=True) # keep GPD1_S24^S27 only (paper)
        self.df.drop([8, 9], axis=0, inplace=True) #DIG2 clustered into S225 (paper)
        
        # DIG1 has 4 clusters: 1      we pick up randomly one of them when there is
        # more than 1. no NAs exccept 2 in S272^S275+T277^S279
        #S126+S127
        #S142 and S330 and S395+S397. kee the latest with ;ean CV 13% against 33 and 16%
        self.df.drop([5,1], axis=0, inplace=True)
        #S272  AND   S272^S275. first has CV 10% and 0 NA; second has 15%.later has 1 NA
        # drop S272^S275
        self.df.drop(3, axis=0, inplace=True)
        #S272^S275+T277^S279

        # NO a0_t0 point so removed
        self.df.drop(self['FAR1_S114'].index[0], inplace=True)

        # 4 NAs
        self.df.drop(self['SSK1_S110'].index[0], inplace=True)
        self.df.drop(self['SSK1_S673'].index[0], inplace=True)


        # STE20 clustering is difficult despite 9 peptides
        #We remove T413+S418 because it contains lots of 16 NAs anyway ot of 36
        #and 6 out of 16 (MIDAS)
        index = self.df[self.df.Identifier == "STE20_T413+S418"].index[0]
        self.df.drop(index, inplace=True)

        # 35% errors
        index = self.df[self.df.Identifier == "STE20_T546+S547"].index[0]
        self.df.drop(index, inplace=True)

        # we can drop RCK2. It is going to be removed in the export of MIDAS
        # isnce not in the PKN, but we can do it here as well
        # true for FPS1, SIC, TEC
        nodes = ['FPS1', 'RCK2', 'SIC1', 'TEC1']
        self.df = self.df[self.df.Protein.apply(lambda x: x not in nodes)]

    def corr2d(self, name1, name2, show=True):
        d1 = self.get_data_matrix(name1).fillna(0).as_matrix()
        d2 = self.get_data_matrix(name2).fillna(0).as_matrix()
        d1 /= pylab.sqrt(pylab.sum(d1**2))
        d2 /= pylab.sqrt(pylab.sum(d2**2))
        from  scipy.signal import correlate2d
        im = correlate2d(d1, d2)
        if show:
            pylab.clf();
            pylab.subplot(2,2,1)
            pylab.imshow(d1, interpolation="None")
            pylab.subplot(2,2,2)
            pylab.imshow(d2, interpolation="None")
            pylab.subplot(2,2,3)
            pylab.imshow(im, interpolation="None")
            pylab.colorbar()
        return im

    def corr2d_all(self, mode="max"):
        ident = self.df.Identifier
        N = len(ident)
        df = pd.DataFrame(np.zeros((N,N)), columns=ident.values)
        df.index = df.columns
        for i1 in self.df.Identifier.index:
            print(i1)
            for i2 in self.df.Identifier.index:
                corr = self.corr2d(ident.ix[i1], ident.ix[i2], show=False)
                if mode=="max":
                    M = corr.max()
                elif mode=="integral":
                    M = corr.sum()/36.
                elif mode=="center":
                    M = corr[5,5]
                else:
                    raise ValueError("mode must be max, integral or center")
                df.ix[ident.ix[i1], ident.ix[i2]] = M
        return df


    def plot_group(self, protein):
        df = self.measurements.ix[self.groups[protein]].transpose()
        df.columns = self.df.Identifier.ix[self.groups[protein]]
        df.plot()

    def remove_protein(self, identifier):
        self.df.drop(self[identifier].index[0], inplace=True)

    def remove_species_not_in_pkn(self, filename):
        c = CNOGraph(filename)
        # let us get rid of the data that will not be available
        nodes = [y for y in set(self.df.Protein) if y not in sorted(set([x.split("_")[0] for x in c.nodes()]))]
        # ['FPS1', 'RCK2', 'SIC1', 'TEC1']
        self.df = self.df[self.df.Protein.apply(lambda x: x not in nodes)]
        return nodes

    def get_group_psite_transposed(self, name):
        """Get a dataframe containing a group of protein

        :param str name: a valid protein name. See :attr:`df.Protein` to get a
            list
        :return: transposed dataframe (indices are psites, columns are
            measurements for ech combination of alpha/NaCl

        .. seealso:: :meth:`get_group_psite`

        """
        return self.cluster.get_group(name).transpose()

    def get_group_psite(self, name):
        """Get a dataframe containing a groupd of psites for one protein

        :param str name: a valid protein name. See :attr:`df.Protein` to get a
            list
        :return: dataframe (indices are psites, columns are
            measurements for ech combination of alpha/time (see class
            documentation)

        .. plot::
            :include-source:
            :width: 80%

            >>> from msdas import yeast
            >>> y = yeast.YEAST2MIDAS()
            >>> y.get_group_psite("DIG1").plot()

        .. seealso:: :meth:`get_group_psite_transposed`
        """
        return self.cluster.get_group(name)

    def groups_cluster(self):
        """Returns dataframe grouped by protein and clustered.

        You must run :meth:`run_affinity_propagation_clustering` before hand
        to populate the **cluster** column.

        ::

            >>> y.cluster = y['DIG1']
            >>> y.cluster.run_affinity_propagation_clustering(preference=-30)
            >>> y.groups_cluster()

        """
        if "cluster" not in self.cluster.df.columns:
            raise ValueError("cluster column not found in the dataframe. See documentation example.")


        return self.cluster.df.groupby(["Protein", "cluster"]).groups

    def get_psites_exemplars(self, preference=-30):
        """Returns list of psites corresponding to the exemplars found in the clustering

        :param float preference: the parameter of the Affinity Propagation
            algorithm. See :class:`msdas.cluster.Affinity`

        Affinity propagation algorithm is run on each protein. For each of them,
        clusters and their exemplars are found. This method returns the list of
        exemplars

            >>> psites = y.get_psites_exemplars(preference=-30)

        .. seealso:: :meth:`get_psites_mapping`.

        """
        psites_tokeep = []
        proteins = list(set(self.df.Protein))
        self.af_results = {}
        for protein in proteins:
            self.logging.debug(protein),
            if len(self.groups[protein]) == 1:
                self.logging.debug(" : no clustering required.")
                psite = list(self.df.ix[self.groups[protein]]['Identifier'])
                psites_tokeep.append(psite[0])
            else:
                
                af = clustering.Affinity(self.get_group_psite(protein),
                    method="euclidean", transpose=True, preference=preference,
                    verbose=False)
                # Get the indices in the entire dataframe of the exemplars:
                indices = np.array(self.groups[protein])[af.cluster_centers_indices]
                self.logging.debug("Found %s clusters " % len(indices))
                psites_tokeep.extend(self.df.Identifier[indices])
                self.af_results[protein] = af
        return psites_tokeep

    def get_psites_mapping(self, preference=-30):
        """Returns exemplars and the list of proteins in the same cluster

        Affinity propagation algorithm is run on each protein. For each of them,
        clusters and their exemplars are found. This method returns a dictionary
        where each key is an exemplar and its value contains the list of
        protein/psites that belongs to this cluster.

        :param float preference: the parameter of the Affinity Propagation
            algorithm. See :class:`msdas.clustering.Affinity`
        :return: dictionary as explained in the documentation here above.

        .. seealso:: :meth:`get_psites_exemplars`.

        """
        # this call populates the af_results that can be used here below
        self.logging.debug("Entering get_psites_mapping----------------------------")
        psites = self.get_psites_exemplars(preference=preference)

        all_proteins = list(set(self.df.Protein))
        res = {}
        for protein in all_proteins:
            if protein in self.af_results.keys():
                self.logging.debug("in clustered protein")
                af = self.af_results[protein]
                names = af.df.columns
                d = dict([(names[af.cluster_centers_indices[i]] ,
                    list(names[af.labels==i])) for i in  set(af.labels)])
                res.update(d)
            else:
                name = [x for x in psites if x.split("_")[0] == protein]
                res.update({protein:name})
        for name in sorted(res.keys()):
            self.logging.debug("%20s" %  name + "\t" + ", ".join(res[name]))
        return res


    def get_df_exemplars(self, preference=-30, normalise=None):
        """Returns a normalised dataframe that corresponds to the exemplars

        The psites are obtained from :meth:`get_psites_exemplars`.

        :param float preference: the parameter of the Affinity Propagation
            algorithm. See :class:`msdas.cluster.Affinity`
        :param normalise: minmax / maxonly / timezero. Set to None to ignore
            normalisation


        """
        psites = self.get_psites_exemplars(preference=preference)
        newdf = self.cluster.get_group()[psites].transpose()

        if normalise=="minmax":
            for index in newdf.index:
                M = newdf.ix[index].max()
                m = newdf.ix[index].min()
                newdf.ix[index] = (newdf.ix[index] - m)/(M-m)
        elif normalise=="maxonly":
            for index in newdf.index:
                M = newdf.ix[index].max()
                newdf.ix[index] = newdf.ix[index] / M
        elif normalise == "timezero":
            for index in newdf.index:
                m = newdf.ix[index][0]
                newdf.ix[index] = newdf.ix[index] - m
                M = newdf.ix[index].max()
                newdf.ix[index] /= M
        else:
            pass
        return newdf

    def to_midas(self, filename, preference=-1):
        """Given a dataframe, export measurements into MIDAS format

        :param str filename:
        :param int preference: a negative parameter used in the clustering. See
             :class:`msdas.clustering.Affinity`
        :return: the MIDAS object.

        Internally, a new dataframe is created by selecting the
        examplars from the clustering (Default preference of -1 is in principle
        equivalen of not having any clustering). Then, the MIDAS file is
        created and saved into a file.


        ::

            y = yeast.YEAST2MIDAS(get_yeast_small_data(), get_yeast_raw_data() )
            m = y.to_midas(df, "MD-test.csv")

            from cno import XMIDAS
            m = XMIDAS("MD-test.csv")

        """
        # create a midas builder
        df = self.get_df_exemplars(preference=preference)
        mb = self._get_midas_builder_from_df(df)
        # and save to
        m = mb.xmidas
        m.to_midas(filename)
        return m

    def _get_midas_builder_from_df(self, df):
        # exports alpha=0, as one condition, NaCl=0 as 1 conditions and
        # alpha=NaCl= on (where t_alpha = t_NaCl) as a thrid condition
        from cno import MIDASBuilder
        m = MIDASBuilder()
        from cno.io.measurements import Measurement as Experiment

        _measures = ['a0_t0', 'a0_t1', 'a0_t5', 'a0_t10', 'a0_t20', 'a0_t45',
                'a1_t1', 'a5_t5', 'a10_t10', 'a20_t20', 'a45_t45',
                'a1_t0', 'a5_t0', 'a10_t0', 'a20_t0', 'a45_t0']
        df = df[_measures]

        inhibitors = {}
        for psite in df.index:
            for col in df.columns:
                value = df.ix[psite][col]

                alpha, nacl = col.split("_")
                if col == "a0_t0":
                    # special case that need to be added 3 times i.e. a0_t0
                    time = 0
                    stimuli = {"a":1, "NaCl":1}
                    e = Experiment(psite, time, stimuli, inhibitors, value)
                    m.add_measurements([e])

                    stimuli = {"a":1, "NaCl":0}
                    e = Experiment(psite, time, stimuli, inhibitors, value)
                    m.add_measurements([e])

                    stimuli = {"a":0, "NaCl":1}
                    e = Experiment(psite, time, stimuli, inhibitors, value)
                    m.add_measurements([e])
                elif col.startswith("a0"):
                    # case alpha=0, NaCl=1 e.g., a0_t5
                    time = int(nacl[1:])
                    stimuli = {"a":0, "NaCl":1}
                    e = Experiment(psite, time, stimuli, inhibitors, value)
                    m.add_measurements([e])
                elif col.endswith("t0"):
                    # case alpha=1, NaCl=0 e.g., a5_t0
                    time = int(alpha[1:])
                    stimuli = {"a":1, "NaCl":0}
                    e = Experiment(psite, time, stimuli, inhibitors, value)
                    m.add_measurements([e])
                else:
                    # case alpha=1, NaCl=1 e.g. a5_t5
                    time = int(alpha[1:])
                    # could be time = int(nacl[1:])
                    stimuli = {"a":1, "NaCl":1}
                    e = Experiment(psite, time, stimuli, inhibitors, value)
                    m.add_measurements([e])
        return m


    def get_expanded_cnograph(self, filename_pkn_no_sites, df):
        """Expands node in a SIF into their psites

        :param str filename_pkn_no_sites: a valid filename to a PKN in SIF
            format.
        :param dataframe df: a dataframe with protein/psites names as indices
            like the output of :meth:`get_df_exemplars`.


        df = y.get_df_exemplars(-1)
        y.get_expanded_cnograph("../../share/data/PKN-yeastScaffold.sif", df)

        :return: a CNOGraph instance

        """
        from cno import CNOGraph
        c = CNOGraph(filename_pkn_no_sites)

        proteins = [this.split("_")[0] for this in df.index]
        proteins = list(set(proteins))
        c._signals = [x.replace("^", "~") for x in proteins[:]]
        for protein in proteins:
            psites = [this for this in df.index if this.split("_")[0] ==  protein]
            psites = [psite.replace("^", "~") for psite in psites]
            psites = [psite.replace("+", "-") for psite in psites]
            self.logging.debug(protein ,psites)
            if psites is not None and len(psite) == 1:
                # just rename the node in the PKN if found
                if protein in c.nodes():
                    #print("Warning. renaming %s into %s" % (protein, psites[0]))
                    c = c.relabel_nodes({protein: psites[0]})

                else:
                    self.debug("Warning. %s not found in PKN" % protein)
            else:
                if protein in c.nodes():
                    self.debug("split {} in {} nodes".format(protein, psites))
                    c.split_node(protein, psites)
                else:
                    self.debug("Warnig. %s not found in PKN" % protein)

        c._stimuli = ["a", "NaCl"]
        #c._signals = [df.index

        return c

    def export_pkn_and_midas_june(self, pkn, tag="undefined"):
        df = self.measurements
        df.index = self.df.Identifier
        c = self.get_expanded_cnograph(pkn, df)

        # midas
        m = self.get_midas()
        m.df.columns = [x.replace("^", "~") for x in m.df.columns]
        m.df.columns = [x.replace("+", "-") for x in m.df.columns]
        m.sim.columns = m.df.columns[:]
        m.errors = m.sim.copy()
        #FIXME: bug in cellnopt.core.xmidas need to set cellline manually if created from mifdasbuilder
        m._cellLine = "undefined"
        # FIXME: make sure the order is correct
        #m.experiments.columns = ["NaCl", "a"]

        df = self.get_cv()
        errors = self._get_midas_builder_from_df(df).xmidas
        errors.df.columns = [x.replace("^", "~") for x in errors.df.columns]
        errors.df.columns = [x.replace("+", "-") for x in errors.df.columns]

        # FIXME: should use the mapping here
        errors.df = errors.df[m.df.columns]
        errors.create_empty_simulation() # FIXME required to update sim and have corret plotting

        # rename
        return c,m, errors

    def export_pkn_and_midas(self, pkn_filename,
            preference=-30, tag="undefined"):
        """Creates the expanded PKN and MIDAS file given PKN and normalise method

        Saves new PKN and MIDAS into PKN-Yeast_psites.sif and MD-Yeast_psites.cv

        :param str pkn_filename:

        :param float preference: the parameter of the Affinity Propagation
            algorithm. See :class:`msdas.cluster.Affinity`

        :return: a tuple with cnograph and midas instances.

        """

        output_pkn_filename = "PKN-Yeast_psites_%s.sif" % tag
        output_midas_filename = "MD-Yeast_%s.csv" % tag

        normalise = "None"

        print("Generating new dataframe given preference of {} for the clustering".format(preference))
        newdf = self.get_df_exemplars(normalise=normalise, preference=preference)
        #newdf = newdf[self.df.columns]
        #print newdf.columns

        print("Expanding the original PKN and saving into {}".format(output_pkn_filename))
        c = self.get_expanded_cnograph(pkn_filename, newdf)

        print("Creating  the MIDAS file into {} ".format(output_midas_filename))
        # raw data contains all measurements
        m = self._get_midas_builder_from_df(newdf)
        m = m.xmidas

        m.df.columns = [x.replace("^", "~") for x in m.df.columns]
        m.df.columns = [x.replace("+", "-") for x in m.df.columns]
        m.sim.columns = m.df.columns[:]
        m.errors = m.sim.copy()

        # need to remove species in MIDAS that are not in the PKN

        #pkn = [x for x in c.nodes() if x not in c._find_and_nodes()]
        not_found = []
        #print m.names_species
        for name in m.names_species:
            if name not in [x for x in c.nodes() if x not in c._find_and_nodes()]:
                print("{} not found in PKN. Removing from MIDAS file".format(name))
                # FUS3_S177+T180~Y182 not found in PKN
                not_found.append(name)
        print("to be removed")
        print(not_found)
        #print m.df.index
        #print m.df.columns
        print(m.df.columns)
        m.remove_species(not_found)

        #FIXME: bug in cellnopt.core.xmidas need to set cellline manually if created from mifdasbuilder
        m._cellLine = "undefined"
        # FIXME: make sure the order is correct
        #m.save(output_midas_filename)

        # get the errors, filter with respect to the species in m
        e = self.xmidas_errors()

        for this in e.df.columns:
            if this not in m.df.columns:
                e.df.drop(this,axis=1, inplace=True)
        # need to reset the simulation. This is usefule for the layout
        e.create_empty_simulation()
        e._cellLine = "undefined"
        # FIXME: make sure the order is correct
        e.experiments.columns = ["NaCl", "a"]
        #m.save(output_midas_filename)

        return (c, m, e)

    def get_errors(self,  default_cv=0.5):
        """Return errors (coefficient variation) for each identifier

        To get the actual errors scaled, you need to multiply by the mean, that is
        the data itself.

        .. todo:: could be move to replicates module
        """
        errors = []
        self.logging.info("filling NA with {}".format(default_cv))

        for this in self.df.Identifier:
            error = self.get_coefficient_variation(this,  default_cv=default_cv)
            errors.append(error)
        df = pd.DataFrame(errors)
        df.index = self.df.Identifier  # set psites as indices
        df = df[self._measures] # rearrange columns
        return df

    def xmidas_errors(self, default_cv=0.5):
        """Return coefficient of variation into a XMIDAS object"""
        df = self.get_errors(default_cv=default_cv)
        errors = self._get_midas_builder_from_df(df).xmidas
        errors.df.columns = [x.replace("^","~").replace("+", "-") for x in errors.df.columns]
        try:
            errors.cellLine = "undefined"
        except:
            pass
        return errors

    def mapping_small_to_raw_june(self):

        mapping = {}
        for k in self.df.Identifier:
            mapping[k] = k

        mapping['DIG1_S272^T277^S279'] = 'DIG1_S272^S275+T277^S279'
        mapping['GPD1_S24+S25+S27'] = 'GPD1_S23+S24+S25+S27'
        mapping['GPD1_S23+S24'] = 'GPD1_S23+S24+S27'
        mapping["RCK2_S33+T35+T44+S46"] = "RCK2_S32+S33+T35+T44+S46"
        mapping["RCK2_S45"] = "RCK2_S45+S46"
        mapping['RCK2_S32+S33'] = "RCK2_S32+S33+T35"
        mapping['DIG1_S395'] = 'DIG1_S395+S397'
        mapping['GPA1_S199'] = 'GPA1_T189+S199+S200'
        mapping['PTP2_S258'] = 'PTP2_Y257+S258'
        mapping['PBS2_S68'] = 'PBS2_S68+S71+S83'
        mapping['FUS3_T180'] = 'FUS3_T173+S177+T180+Y182'
        mapping["FUS3_S177+T180^Y182"] = 'FUS3_S177^T180+S177^Y182+T180^Y182'
        mapping["SKO1_S94^S108^T113"] = "SKO1_S94+S96^S108^T113"
        mapping['SSK1_S351'] = "SSK1_S350+S351"
        mapping['STE11_S326'] = "STE11_S323+S326+S326"  # TODO fix this psite name
        mapping['STE20_S192'] = "STE20_S192+S195"
        mapping['STE20_S195'] = "STE20_S192+S195+S196"
        mapping['STE20_S196^T197'] =  "STE20_S192^S195+S196^T197"
        mapping['STE20_S418'] = "STE20_T413+S418"
        mapping["STE20_T170^T172"] = "STE20_S169+T170^T172"
        mapping["STE20_T203^T207"] = "STE20_T203^T207+T217+T218"
        mapping["STE20_T573"] = "STE20_T573+T575"
        #mapping["STE20_T511"] = "STE20_T511_1"
        mapping['DIG2_T83'] = 'DIG2_T83+S84'
        mapping['DIG2_S84'] = 'DIG2_S84+T83'
        mapping['SIC1_S201'] = 'SIC1_S198+S201'
        mapping['SIC1_T173'] = 'SIC1_T173+S175'
        mapping['SIC1_S191'] = 'SIC1_S191_1'
        mapping['DIG2_S84'] = 'DIG2_T83+S84_2'
        mapping['DIG2_T83'] = 'DIG2_T83+S84_1'

        """# for now, let us ignore DIG1, DIG2, and STE12, downstream of FUS3
        for k,v in mapping.iteritems():
            for this in ['DIG1', 'DIG2', 'STE12']:
                if k.startswith(this):
                    mapping[k] = None
        """
        return mapping



    def get_midas(self):
        """Return subset of the small data as a MIDAS object"""
        data = self.measurements
        data.index = self.df.Identifier
        m = self._get_midas_builder_from_df(data).xmidas
        m.df.columns = [x.replace("^", "~") for x in m.df.columns]
        m.df.columns = [x.replace("+", "-") for x in m.df.columns]
        return m

    def get_cv(self):
        """Return coefficient of variation with proper indices.

        .. todo:: could be in replicates
        """
        cv = pd.DataFrame(self.replicates.get_coefficient_variation())
        cv.index = self.replicates.df.Identifier
        return cv

    def get_coefficient_variation(self, identifier, default_cv=0.5):
        """Return CV for a particular identifier

        .. todo:: could be in replicates
        """
        # no more duplicatesd rows to average
        #mapping = self.mapping_small_to_raw_june()

        cv = self._cv_buf

        indices = list(self.replicates.df[self.replicates.df.Identifier == identifier].index)
        if len(indices)>1:
            self.logging.warning("get_coefficient_variation on several rows...")
        errors = cv.ix[indices]
        errors.index = self.replicates.metadata.ix[indices]['Identifier']
        errors = errors.mean()  # takes the mean of the errors.

        errors.fillna(default_cv, inplace=True)
        return errors

    def pcolor_errors(self, vmax=None, cmap="hot_r", default_cv=np.nan, fontsize=8, vmin=0,
                      salt=False):
        """plot coefficient of variation for the small daa set using replicates
        found in the raw data set.


        :param salt:  if set to True, re-arrange the columns based on salt
            rather than alpha
        :return: errors in a dataframe

        .. plot::
            :include-source:
            :width: 80%

            from msdas import *
            y  = YEAST2MIDAS(get_yeast_small_data(), get_yeast_raw_data())
            y.cleanup_june() # replaces data with regenerated data including NAs
            errors = y.pcolor_errors(vmax=1,vmin=.2,fontsize=7)


        """
        errors = self.get_errors(default_cv=default_cv)
        errors = errors.ix[sorted(errors.index)]
        if salt:
            errors = errors[self._measures_salt]
        mask = np.ma.array(errors.as_matrix(), mask=np.isnan(errors.as_matrix()))

        pylab.clf();
        if vmax == None:
            vmax = errors.max().max()

        cmap = pylab.get_cmap(cmap)
        cmap.set_bad("grey", 1)
        pylab.pcolormesh(mask, vmin=vmin, vmax=vmax, cmap=cmap)
        pylab.colorbar()
        N, M = errors.shape
        pylab.ylim([0, N])
        pylab.xlim(0, M)
        pylab.xticks([0.5 +x for x in range(0, M)], errors.columns,
                      rotation=90, fontsize=fontsize)
        pylab.yticks([0.5+x for x in range(0, N)], errors.index, fontsize=fontsize)
        # adding the errors to be conservative
        pylab.tight_layout()
        return errors


    def pcolor_na(self, raw=False, fontsize=8):
        """Plot number of NA for protein that are in the small data set

        Final number of rows is therefore larger thatn in the small data set.

        overlap between small and raw data is 32 rows. the remaining 25 are combination from the
        raw data set.

        .. plot::
            :include-source:
            :width: 80%

            from msdas import *
            import pylab
            y  = YEAST2MIDAS(get_yeast_small_data(), get_yeast_raw_data())
            y.cleanup_june()
            errors = y.pcolor_na()

        """
        proteins = list(set(self.df.Protein))
        r = self.replicates
        psites = r.metadata.ix[r.metadata.query("Protein in proteins", engine="python").index].Identifier
        tags = self._measures

        NAs = {}
        for tag in tags:
            df = r.get_replicates_from_one_unique_measurement(tag)
            nas = 3 - pd.notnull(df.ix[psites.index]).sum(axis=1)
            NAs[tag] = nas.values

        NAs = pd.DataFrame(NAs, index=psites)
        NAs = NAs[self._measures]
        NAs = NAs.ix[sorted(NAs.index)]
        pylab.clf()
        pylab.pcolor(NAs)
        N, M = NAs.shape

        pylab.ylim([0,N])
        pylab.xlim(0, M)
        pylab.xticks([0.5 +x for x in range(0,M)], NAs.columns, rotation=90, fontsize=fontsize)
        pylab.yticks([0.5+x for x in range(0,N)], NAs.index, fontsize=fontsize)
        pylab.colorbar()
        pylab.tight_layout()
        return NAs

    def plot_psites_mapping(self, preference=-30):
        """For each protein, plot all data with the exemplar highlighted.

        The clustering is made with :mod:`msdas.clustering` and its affinity
        propagation algorithm. In the yeast case, clustering is performed on
        euclidean distance.


        ::

            from msdas import *
            y = yeast.YEAST2MIDAS(get_yeast_small_data(), get_yeast_raw_data(), drop_non_midas=False
            y.plot_psites_mapping()


        """
        mapping = self.get_psites_mapping(preference=preference)
        for k,v in mapping.iteritems():
            if len(v)>1:
                protein = k.split("_")[0]
                df = self.df.query("Protein==protein", engine="python")

                df = df.set_index("Identifier").drop(self.metadata.drop("Identifier", axis=1), axis=1)
                df = df.transpose()

                v.remove(k)
                ax = df.apply(lambda x:x/np.sqrt((x.abs()**2).sum()))[v].plot(marker='o', lw=1)
                df.apply(lambda x:x/np.sqrt((x.abs()**2).sum()))[k].plot(lw=2, ax=ax, marker='o')

    def plot_timeseries_midas(self, psite="DIG1_S126+S127",  sigma=2, hold=False):
        """Plot data related to a given identifier restricted to the MIDAS-compatible set

        Data that can be used in the ODE package must be time-series. This
        correspond to salt=0, or NaCl=0 or salt==NaCl; IN the 6x6 matrix,
        this correspond to first row, first colum ad diagonal.

        .. plot::
            :include-source:
            :width: 80%

            from msdas import *
            y = YEAST2MIDAS(get_yeast_small_data(), get_yeast_raw_data())
            y.plot_timeseries("DIG1_S126+S127")

        Errors (2 sigmas) are also shown
        """

        l = [0,1,5,10,20,45,46,50,55,65,90,91,95,100,120,145]

        _measures = ['a0_t0', 'a0_t1', 'a0_t5', 'a0_t10', 'a0_t20', 'a0_t45',
                'a1_t1', 'a5_t5', 'a10_t10', 'a20_t20', 'a45_t45',
                'a1_t0', 'a5_t0', 'a10_t0', 'a20_t0', 'a45_t0']
        df = self.df[_measures].transpose()

        df.columns = self.df.Identifier
        df['time'] = l
        if hold==False:
            pylab.clf()
        df.plot(x="time", y=psite, marker="o")
        pylab.legend([psite])
        errors = self.get_errors()

        normerr = df[psite]
        pylab.errorbar(df['time'],df[psite], yerr=errors.ix[psite].ix[_measures]*normerr*sigma)
        pylab.axvline(45.5, alpha=0.5)
        pylab.axvline(90.5, alpha=0.5)
        Y0 = df[psite]['a0_t0']
        pylab.axhline(Y0, alpha=0.5)
        pylab.text(20,Y0, "a0_tX")
        pylab.text(70,Y0, "aX_tX")
        pylab.text(120,Y0, "aX_t0")
        pylab.xlim([-.5,146])

    def boxplot_errors(self):
      """ Show errors (coefficient variation) as a boxplot for each identifier"""
      pylab.clf()
      self.get_errors(default_cv=np.nan).transpose().boxplot(vert=False)
      pylab.tight_layout()


    def _find_identifiers_with_missing_time0(self):
        """Return identifiers for which a0_t0 is missing"""
        identifiers = self.df[pd.isnull(self.df['a0_t0'])].Identifier
        return identifiers.values

    def interpolate_time0(self, identifier, cmap="hot_r", method="cubic"):
        """


        :param method: 'linear', 'nearest', 'cubic'

        """
        # the data may contains NA so interp2d does not work. we use griddata
        from scipy.interpolate import griddata
        grid_x = np.array([[0]*6,[1]*6,[5]*6,[10]*6,[20]*6,[45]*6]);
        grid_y = grid_x.T

        data = self.get_data_matrix(identifier).as_matrix()

        points = []
        values = []
        for i in range(0,6):
            for j in range(0,6):
                if np.isnan(data[i,j]) == False:
                    points.append([self.times[i], self.times[j]])
                    values.append(data[i,j])

        grid_z0 = griddata(np.array(points), np.array(values), (grid_x, grid_y),
                           method=method)

        pylab.figure(1)
        pylab.clf()
        pylab.subplot(1,2,1)

        mask = np.ma.array(data, mask=np.isnan(data))
        cmap = pylab.get_cmap(cmap)
        cmap.set_bad("grey", 1)
        #pylab.pcolormesh(pylab.flipud(mask), vmin=vmin, vmax=vmax, cmap=cmap)


        pylab.imshow(data, interpolation="None", cmap=cmap)



        pylab.subplot(1,2,2)
        pylab.imshow(pylab.flipud(grid_z0), origin="lower", cmap=cmap, interpolation='None')
        pylab.suptitle(identifier)
        return data, grid_z0

    def plot_figure6(self, identifier,vmin=0,vmax=2):
        from easydev import colors
        c = colors.ColorMapTools()
        d = {'blue': [0,1,1],
             'green':[0,1,0],
             'red':  [1,1,0]}
        cmap = c.get_cmap(d, reverse=False)


        pylab.clf(); 
        m = self.get_data_matrix(identifier)

        im1 = (1./(m.transpose().divide(m[0], axis=1))).transpose()
        pylab.subplot(1,2,1)
        pylab.imshow(pylab.flipud(im1), origin="lower", interpolation="None",
                cmap=cmap,vmin=vmin,vmax=vmax)
        pylab.colorbar(); 
        pylab.xlim([.5,5.5]); 
        pylab.yticks(range(0,6),  ['45','20','10', '5','1', 'Phe0'])
        for i,label in enumerate([1,2,3,4,5]):
            for j,x in enumerate(im1[label].values):
                pylab.text(1+i, 6-j-1, int(100*x)/100.)


        pylab.subplot(1,2,2)
        im2 = (m.ix[0]/m).transpose()
        pylab.imshow(im2, origin="lower", interpolation="None",  cmap=cmap,
                vmin=vmin,vmax=vmax); 
        pylab.colorbar(); 
        pylab.xlim([.5,5.5]); 
        pylab.yticks(range(0,6),  ['45','20','10', '5','1', 'NaCl0'])
        for i,label in enumerate(im2.columns[1:]):
            for j,x in im2[label].iterkv():
                pylab.text(1+i, 6-j-1, int(100*x)/100.)



        pylab.suptitle("%s" % identifier)


        return im1, im2



def load_sim_data(filename, times=[0,1,5,10,20,45]):
    """
    signals = colnames(cnolist@signals$`0`)
    sim = t(sapply(sim, unlist))
    colnames(sim) = signals
    write.csv(sim, "sim.csv")

     sim = yeast.load_sim_data("sim.csv")
     m = midas.XMIDAS("MD-Yeast_test_alpha0oneexp_maxonly.csv")
     m.sim = sim.copy()
     m.plot(mode="trend")
     m.plot(mode="mse")


    c = CNOGraph("PKN-", "MD-")
    for name in c.nodes():
    if name in diffs.columns:
        c.node[name]['mse'] = diffs[name][0]
    else:
        c.node[name]['mse'] = None


    """
    sim = pd.read_csv(filename, index_col=0)
    N = len(sim.columns)


    sim['time'] = times
    sim['experiment'] = ['experiment_0'] * 6
    sim['cellLine'] = ['undefined'] * 6

    sim = sim.set_index(['cellLine', 'experiment', 'time'])
    return sim





# function for the linear fit to automate the process
def plotfit(x, y, yerr=None, order=1):
    w = None if (yerr is None or np.sum(yerr)==0) else 1/yerr
    p, cov = np.polyfit(x, y, order, w=w, cov=True)  # coefficients and covariance matrix
    yfit = np.polyval(p, x)          # evaluate the polynomial at x

    perr = np.sqrt(np.diag(cov))     # standard-deviation estimates for each coefficient
    R2 = np.corrcoef(x, y)[0, 1]**2  # coefficient of determination between x and y
    resid = y - yfit
    chi2red = np.sum((resid/yerr)**2)/(y.size - 2) if w is not None else np.nan

    #return yfit, p, R2, chi2red, perr, resid
    pylab.errorbar(x, y, yerr=yerr, fmt = 'bo', ecolor='b', capsize=0, elinewidth=2)
    pylab.plot(x,y)
    pylab.xlim([-.5,46])
    pylab.plot(x, yfit, 'r', linewidth=3, color=[1, 0, 0, .5])


