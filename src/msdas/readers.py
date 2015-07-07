# -*- python -*-
#
#  This file is part of MS-DAS software
#
#  Copyright (c) 2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>, Marti Bernardo Faura
#  bernardo@ebi.ac.uk
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#
##############################################################################
"""Mass Spectrometry Readers (raw data)"""
import os
import pandas as pd
from easydev import Logging
import numpy as np
from msdas.psites import PSites
from msdas.tools import SequenceTools, Modification

import pylab

__all__ = ["MassSpecReader", "Cleaner"]


class Cleaner(object):
    """Tools inherited by :class:`msdas.readers.MassSpecReader` to perform some
    tasks on grouping the peptides. Allows also some manipulation on the NAs.

    .. note:: In principle you do not need to use this class, which is
        available via :class:`msdas.readers.MassSpecReader`. However, it can be
        useful as a standalone class.

    Requires the dataframe stored  in :attr:`df` to be populated with a dataframe
    similar to those used in :class:`msdas.readers.MassSpecReader`.


    .. rubric:: handling NA

    Here, we create a simple dataframe with two columns. The last row contains
    1 NA and 1 zero.::

        from msdas import Cleaner
        import pandas as pd

        c = Cleaner()
        c.df = pd.DataFrame({'A':[1,2,3,np.nan], 'B':[4,5,6,0]})

    In Mass Spectrometry, zeros have no meaning. There is always some
    measurements of abundance being made. Zeros are sometimes used to represent a
    NA. The :class:`Cleaner` class has a method :meth:`set_zero_to_na` to replace all
    zeros by NAs::

        >>> c.set_zero_to_na()
        >>> c.df
            A   B
        0   1   4
        1   2   5
        2   3   6
        3 NaN NaN

    Rows that contains all NAs are useless and one may want to remove them::

        c.drop_na_count(2)

    This statement removes all rows that contain 2 NAs. You can get the count
    of all NAs through the entire dataframe with the method :meth:`get_na_count`.

    """
    def __init__(self):
        self.df = None

    def drop_oxidation(self):
        """Drop rows where the Sequence_Phospho contains the Oxidation string


        That the the :meth:`merge_identical_peptide_sequence` renames
        some sequences and add the "MERGED" prefix. Those rows are not to be removed
        since they contains measurements that correspond to a peptide with an oxidation
        but also possible another peptide without oxidation. There is no way currently
        to know that. So the MERGED peptide are ignored. Ideally, :meth:`drop_oxidation`
        should be called before any merging.

        """
        indices = self.df.Sequence_Phospho.apply(lambda x: "Oxidation" in x and "MERGED" not in x)
        try:
            self.logging.warning("Removing {} rows that contain Oxidation string".format(sum(indices)))
        except:
            print("Removing {} rows that contain Oxidation string".format(sum(indices)))

        indices = self.df.Sequence_Phospho.apply(lambda x: "Oxidation" not in x and "MERGED" not in x)
        self.df = self.df[indices]

    def set_zero_to_na(self):
        """Replace all zeros with NA

        If you are manipulating raw mass spec data, a zero most probably means
        NA. Indeed, each measurement is the ion current measured, that is greater
        than zero. So, if you want to replace zeros by NA, just call this method.
        """
        self.df = self.df.where(self.df!=0, None)

    def get_na_count(self):
        """Return vector with number of NA per row (psite/protein)"""
        return len(self.df.columns) - pd.notnull(self.df).sum(axis=1)

    def drop_na_count(self, N=100000, inplace=True):
        """removes rows with at least N NA values

        A row may have replicates. In the YEAST data case, we have 3 replicates
        and 36 experiments; setting N=54 means removes
        all rows for which half of the measurements are NAs.

        Here below, we use :class:`msdas.replicates.ReplicatesYeast`, which
        inherit from :class:`MassSpecReader` and therefore from :class:`Cleaner`.
        Consequently, we can simplify the data by calling :meth:`drop_oxidation`
        and :meth:`merge_identical_peptides_sequence`. Then, we drop
        rows that have 108 NAs and plot the number of remaining rows. Then, we drop
        rows that have 2 NAs and so on until there is no more data. The plot shows
        the size of the dataframe as a function of the maximal number of NAs that
        is acceptable.

        .. plot::
            :include-source:
            :width: 80%

            >>> import pylab
            >>> from msdas import ReplicatesYeast, get_yeast_raw_data
            >>> r = ReplicatesYeast(get_yeast_raw_data())
            >>> r.drop_oxidation()
            >>> r.merge_identical_peptide_sequence()
            >>> size = []
            >>> for n in range(108,-1,-1):
            ...     r.drop_na_count(n)
            ...     size.append(r.df.shape[0])
            >>> pylab.plot(size, 'o-')
            >>> pylab.grid()
            >>> pylab.xlabel("Max number of NA per row")
            >>> pylab.ylabel("Number of rows")


        """
        data = self.get_na_count()
        try:
            self.logging.info("Removing {} rows where number of NA is larger than {}".format(len(data[data>=N]), N))
        except:
            print("Removing {} rows where number of NA is larger than {}".format(len(data[data>=N]), N))

        if inplace:
            self.df.drop(data[data>=N].index, inplace=inplace)
        else:
            return self.df.drop(data[data>N].index, inplace=False)

    # TODO rename into merge_same_psites
    def merge_peptides_with_same_psites(self):
        """Merge rows that are different peptides but same psites

        in other words same psites are presumably peptides including in each other.


        Rows may have same protein and psite but different peptides. This
        may happen when an enzyme cuts the peptide at different places. Consider
        for instance this case::

            KVASMNSASLQDEAEPYDS(Phospho)DEAISK,S49,TOP1,
             VASMNSASLQDEAEPYDS(Phospho)DEAISK,S49,TOP1

        In such cases, we add the data measurements and replace
        the 2 rows by a single one. The Peptide stored is the first one (arbitrary choice).


        Actually, you may also have the Oxidation present, which may be 2 different rows but
        there are added as expected::

                DS(Phospho)SLLFSPAAVAM(Oxidation)R       S515
                DS(Phospho)SLLFSPAAVAMR       S515

        In this case, the 2 rows are added.

        Identifier are unchanged

        .. note:: if two peptides are identical, there are also taken into account here and being added.

        .. seealso:: :class:`msdas.psites.PSites`
        """
        try:
            self.logging.info("Merging ambiguous peptides (different peptides but same psites")
        except:
            print("Merging ambiguous peptides (different peptides but same psites")
        tobeadded = []

        g = self.df.groupby(["Protein", "Psite"])
        groups = g.groups
        count = 0
        count_group = 0

        S = g.aggregate(sum)

        try:
            self.logging.info("Merging {} rows into {} groups.".format(len(self.df), len(groups)))
        except:
            pass

        todrop = []
        for k,v in groups.iteritems():
            if len(v)>1:
                # add the measurements in the group and replace the entire
                # group by this new data
                group = self.df.ix[groups[k]]
                # this does not work :  newrow = group.sum() because some
                # columns contains strings. Consequently, the entire row is
                # casted into strings. So, we need to use apply
                newrow = S.ix[k]
                metadata = group.ix[v[0]].copy()
                newrow['Psite'] = metadata['Psite']
                newrow['Sequence'] = metadata['Sequence']
                newrow['Protein'] = metadata['Protein']
                newrow['Sequence_Phospho'] = metadata['Sequence_Phospho']
                newrow['Identifier'] = metadata['Identifier']
                newrow['Entry'] = metadata['Entry']
                newrow['Entry_name'] = metadata['Entry_name']

                tobeadded.append(newrow)
                todrop.extend(v)

                count += len(v)
                count_group += 1

        self.df.drop(todrop, inplace=True)
        if len(tobeadded):
            self.df = self.df.append(tobeadded, ignore_index=True)

        try:
            self.logging.info("Merged {} groups (with more than 1 peptide) into {} new rows.".format(count, len(tobeadded)))
            self.logging.info("New data frame has {} rows".format(len(self.df)))
        except:
            print("Merged {} groups (with more than 1 peptide) into {} new rows.".format(count, len(tobeadded)))
            print("New data frame has {} rows".format(len(self.df)))


        self._rebuild_identifier()

    def merge_identical_peptide_sequence(self):
        """Merge (by adding) rows that have same protein and peptide sequence

        They may have different psite locations but their number must agree.
        A new name for the psite has to   be created to reflect what addition
        was performed. Instead of concatenating  the psites, we factorise if
        possible. For instance, consider the simple case with no common location
        but a single location::

            DVSQIT(Phospho)SSPK        T44
            DVSQITSS(Phospho)PK        S46

        This is a simple OR case  and the new name is **T44+S46**. Let us now consider a
        case with 2 psites, one being common::

            TRS(Phospho)S(Phospho)GTNNKDVSQITSSPK    S32^S33
            TRSS(Phospho)GT(Phospho)NNKDVSQITSSPK    S33^T35

        S33 is common to two sequence, so the new psite name is ::

            S32+T35^S33

        and so on.

        .. seealso:: :class:`msdas.psites.PSites`
        """
        try:
            self.logging.info("Merging identical peptides (but with psites at different locations")
        except:
            print("Merging identical peptides (but with psites at different locations")
        tobeadded = []
        # parentheses in (Phospho) are required to not clash with prefix added to the
        # sequence_original when there ambiguities
        self.df['count_phospho'] = self.df.Sequence_Phospho.apply(lambda x: x.count("(Phospho)"))

        # ignore sequence without peptide (if already result of a merge)
        df = self.df[self.df['count_phospho']>0]
        g = df.groupby(["Protein", "Sequence", "count_phospho"])
        groups = g.groups

        try:
            self.logging.info("grouping in progress")
        except:
            print("grouping in progress")
        count_group = 0
        count = 0

        pfactory = PSites(verbose=self.level)

        todrop = []

        S = g.aggregate(sum) # compute some before the loop. This is much 2,3 times faster than a apply function
        self.logging.info("Merging {} rows into {} groups.".format(len(self.df), len(g.groups)))

        #print("##############")
        #print(len(groups))
        NN = len(groups)
        count = 0
        for k, indices in groups.iteritems():
            #print(k, count)
            count += 1


            # if there is only one row for a given protein/sequence/number of phosphos,
            # then there is nothing to do otherwise, there are several cases as explained
            # below
            if len(indices)<=1:
                continue

            # S142^S150
            # S142^S151
            #
            # in this case, S142 is not ambiguous and final psite could be named S142^150+S151
            # note that there is no paraenthese. Yet, it means S142^S150 OR S142^S151
            #
            # k[2] conains the counting of phosphos
            #

            # TODO those 3 lines takes 3 seconds ...
            group = self.df.ix[groups[k]]
            Npsite = pfactory.remove_duplicated("^".join(list(self.df.ix[indices].Psite)))
            Npsite = len(Npsite.split("^"))

            if k[2] == 1:
                psite_name = "+".join(self.df.ix[indices].Psite)
                #psite_name = pfactory.sorted(psite_name.split("^"))
                # the apply/sum function is slow...
                newrow = S.ix[k].copy()
                #newrow = group.apply(lambda x: x.sum())
                metadata = group.ix[indices[0]]

                newrow['Psite'] = psite_name
                newrow['Sequence'] = metadata['Sequence']
                newrow['Protein'] = metadata['Protein']

                try:
                    newrow['Identifier'] = metadata['Identifier']
                    newrow['Entry'] = metadata['Entry']
                    newrow['Entry_name'] = metadata['Entry_name']
                except Exception:
                    raise Exception

                newrow['Sequence_Phospho'] = 'MERGED_{}Phosphos_{}locations_'.format(k[2], Npsite) + newrow['Sequence']
                newrow['count_phospho'] = k[2]
                todrop.extend(indices)
                tobeadded.append(newrow)
                count += len(indices)
                count_group += 1
            else:
                psites = list(self.df.ix[indices].Psite)
                common_psites = pfactory.get_common_psites(psites)
                if len(common_psites) == k[2]-1:
                    psite_name = pfactory.get_factorised_psites(psites)
                else:
                    # concatenate with +
                    psite_name  = "+".join(sorted(psites))

                newrow = S.ix[k].copy()
                #newrow = group.apply(lambda x: x.sum())
                newrow['Psite'] = psite_name
                metadata = group.ix[indices[0]]
                newrow['Sequence'] = metadata['Sequence']
                newrow['Protein'] = metadata['Protein']
                try:
                    newrow['Identifier'] = metadata['Identifier']
                    newrow['Entry'] = metadata['Entry']
                    newrow['Entry_name'] = metadata['Entry_name']
                except Exception:
                    raise Exception

                newrow['Sequence_Phospho'] = 'MERGED_{}Phosphos_{}locations_'.format(k[2], Npsite) + newrow['Sequence']
                newrow['count_phospho'] = k[2]
                todrop.extend(indices)
                #self.df.drop(indices, inplace=True)
                tobeadded.append(newrow)
                count += len(indices)
                count_group += 1

        self.df.drop(todrop, inplace=True)
        if len(tobeadded):
            self.df = self.df.append(tobeadded, ignore_index=True)
        try:
            self.logging.info("Merged {} groups (with more than 1 peptide) into {} new rows.".format(count, len(tobeadded)))
        except:
            print("Merged {} groups (with more than 1 peptide) into {} new rows.".format(count, len(tobeadded)))
        del self.df["count_phospho"]
        self._rebuild_identifier()

    def fix_duplicated_identifier(self):
        """Make sure identifierse are unique

        Some identifiers may be duplicates (same protein, same psites
        but different peptides) e.g., ::

             Protein                Sequence    Psite     Identifier
                DIG2  LSQKEEDHSGKPPTITTSPAEK  T83+S84   DIG2_T83+S84
                DIG2      EEDHSGKPPTITTSPAEK  T83+S84   DIG2_T83+S84

        In this case, we append an additional id, which is an integer
        starting from 1. First row found has the id 1.::

             Protein                Sequence    Psite     Identifier
                DIG2  LSQKEEDHSGKPPTITTSPAEK  T83+S84   DIG2_T83+S84_1
                DIG2      EEDHSGKPPTITTSPAEK  T83+S84   DIG2_T83+S84_2

        """
        groups = self.df.groupby("Identifier").groups
        for k, indices in groups.iteritems():
            if len(indices) > 1:
                identifiers = list(self.df.ix[indices,'Identifier'])
                identifiers = [identifier + "_" + str(i+1) for i,identifier in enumerate(identifiers)]
                self.df.ix[indices,'Identifier'] = identifiers





class MassSpecReader(Logging, SequenceTools, Cleaner):
    """Read a simple Mass Spec data CSV file

    First, the MassSpecReader reads a CSV file and stores it in the
    :attr:`df` attribute as a dataframe. At this stage, there is no
    manipulation except from the fact that columns found in the header are
    cleanup (e.g., spaces are removed).

    Second, the header is interpreted. The User Guide provides a draft
    convention of how the header should look like. It should contains the
    following columns:

    #. **Protein**: a valid uniprot protein name (e.g., DIG1)
    #. **Sequence**: a peptide sequence without annotations. If not provided,
       Sequence_Phospho must be provided
    #. **Sequence_Phospho**: a sequence with the phsphorylation position tagged
       as "(Phospho)" e.g., AVS(Phospho)KK means there was a Phophorylation
       on S at position 3 and possibly oxidation as (Oxidation)
    #. **Psite**: The absolute position of the phosphorylations. For instance S132^S140

    and possibly:

    #. **Entry_name**: a valid uniprot entry name (e.g., DIG1_YEAST)
    #. **Entry**: a valid uniprot entry (e.g., P23321)
    #. **Identifier**: concatenation of Protein and Psite. Maybe overwritten if
       :meth:`_rebuild_identifier` is called

    If you have another name (e.g., PeptideSequence instead of Sequence_Phospho),
    you can populate the :attr:`header_mapping` attribute.


    If "Accession" and "Modifications" are found in the header, the accession
    column will be expected to contain these kind of strings (FASTA header)::

        sp|Q8IYB3-2|SRRM1_HUMAN

    from which the Protein column can be extract. Besides, the Modifications
    will contain information that could be used to
    extract the Psite positions. Modifications entries look like::

        [5] S+7777| [1] T+5555

    Once the data is intrepreted. Protein, Sequence and Psite should be available.

    Finally, a third step performed by this class is to performs some cleanup.
    See :meth:`cleanup` for details.

    .. warning:: Psites are separated by the ^ character (AND). If there is
        an ambiguity, the + sign can be used as a logical OR.


    The creation of a MassSpecReader is quite versatile. You can read a CSV
    file as described above, or provide an instance of MassSpecReader, or an object
    that contains an attribute :attr:`df`, which is compatible.


    ::

        # Let us create a dummy dataframe
        import pandas as pd
        df = pd.DataFrame({
            'Psite': ['S1+S2', 'T3'],
            'Sequence_Phospho': ["S(Phospho)S(Phospho)TTT", "SST(Phospho)TT"],
            'Protein': ['DIG1', 'DIG2'],
            'X1':[1,2],
            'X2':[5,2],
            'X3':[4,5],
            'dummy'[1,2]
        }

        from msdas import *
        r = readers.MassSpecReader(df)

    There are some sanity checks performed. In the example above, the case Psite, S1+S2
    is acutally ambiguous: it means that there a psite on S1 OR S2 (i.e., there is one psite).
    Yet the sequence contains 2 phophoosites. So, there is a typo. If you check the content of
    **r.df** you would see only 1 row. Let us fix it

    ::

        >>> df['Psite'] = ['S1^S2', 'T3']
        >>> r = readers.MassSpecReader(df)
        >>> len(r.df)
        2

    .. seealso:: the User Guide and notebooks provided in the MSDAS documentation


    """
    #: symbol used to encode psites (AND case)
    AND = "^"
    #: symbol used to encode psites (OR case)
    OR = "+"

    header_mapping = {
        'Protein': ['Standard name', 'ProteinName'],
        'Psite': ['Phosphosites'],
        'Sequence_Phospho': ['Peptide sequence', 'PeptideSequence']
    }

    valid_header = ["Identifier", "Protein", "Sequence", "Sequence_Phospho",
                    "Psite", "Entry", "Entry_name"]
    def __init__(self, data=None, mode="default", verbose=True, quotechar="'", merge_peptides=False,
                 index_col=None, sep=",",cleanup=True):
        """.. rubric:: constructor

        :param str data: Can be (i) None, (ii) a valid filename with the appropiate format (CSV), (iii)
            an instance of MassSpecReader itself, (iv) a dataframe with the correct column names (v)
            any object that contains an attribute **df** that is a dataframe with the correct column names.
        :param bool verbose: Set to True/False or valid string (DEBUG,INFO,WARNING,ERROR)
        :param quotechar: Ignore the quotechar argument if found in the CSV file.
            Defaults to ' character. For Tcell data, set to " since
            the CSV file contains " characters.
        :param bool merge_peptides: calls :meth:`merge_peptides` is True (default to False)
        :param int index_col: may be required to read some CSV files
        :param str sep: the seprator to be used when reading the CSV file
        :param cleanup: calls :meth:`cleanup` method (default to True)

        """
        super(MassSpecReader, self).__init__(level=verbose)
        self._quotechar = quotechar
        self._index_col = index_col
        self._sep = sep
        self._info = {}

        self._merge_peptides = merge_peptides

        self.removed = {}
        # must be set to zero now
        self._mode = None
        self._df = None

        if isinstance(data, MassSpecReader):
            self.df = data.df.copy()
            self.mode = data.mode

        elif isinstance(data, pd.DataFrame):
            self.df = data.copy()
            data = True # to force the cleaup and avoid co,paring dataframe with None
            self.mode = mode
        elif hasattr(data,"df"):
            self.df = data.df.copy()
            self.mode = mode
        # comparison with None should be after dataframe case
        elif data == None:
            self.mode = mode
        else:
            # mode must be set now, before calling read_csv
            self.mode = mode
            if data != None and os.path.exists(data) == False:
                raise IOError("The file {} does not exist".format(data))
            if data!= None:
                self.read_csv(data, quotechar=self._quotechar,
                              index_col=self._index_col, sep=self._sep)

        # if this is failing, we want to raise an error. no try/except
        if cleanup == True and data!=None:
            self.cleanup()

        self._metadata_names = self.valid_header[:]

    def _set_df(self, df):
        self._df = df
        self._interpret_header()
    def _get_df(self):
        return self._df
    df = property(_get_df, _set_df, doc="getter/setter for the main dataframe")

    def _get_N(self):
        return len(self.measurements.columns)
    N = property(_get_N, doc="Return number of measurements.")

    def _get_metadata(self):
        return self.df[[x for x in self._metadata_names if x in self.df.columns]]
    metadata = property(_get_metadata,
                        doc="""read-only attribute. Returns subset of the :attr:`df` attribute, which contains
                        metadata (not measurements). The subset is selected based on the :attr:`_metadata_names`
                        attribute, which by default contains Sequence, Protein, Identifier, Entry, Entry_name,...
                        look at  `_metadata_names` for up-to-date list of column names""")

    def _get_measurements(self):
        col = [x for x in self.df.columns if x not in self._metadata_names]
        return self.df[col]
    measurements = property(_get_measurements,
        doc="getter to subset of the dataframe :attr:`df` that is not in :attr:`metadata`")

    def _get_mode(self):
        if self._mode == None:
            # try to figure out
            self.logging.warning("mode not provided. Trying to figure out")
            if "Entry_name" in self.df.columns:
                mode = self.df.Entry_name.ix[0].split("_")[1]
                self.mode = mode
            else:
                self.logging.warning("Could not figure out. please provide YEAST")
        return self._mode
    def _set_mode(self, mode):
        if mode == None:
            raise ValueError("mode cannot be None")
        assert mode.upper() in ["DEFAULT", "YEAST"], "mode can be yeast "
        if mode.upper() == "DEFAULT":
            mode = "YEAST"
        self._mode = mode.upper()
        if self._mode == "YEAST":
            self._quotechar = "\'"
    mode = property(_get_mode, _set_mode,
                    doc="mode of the Merger (e.g., yeast)")

    def _set_sequence(self):
        self.debug("set_sequence starting")
        if "Sequence_Phospho" not in self.df.columns:
            return

        tags = ["(Phospho)", "(Oxidation)"]

        if "Sequence" in self.df.columns:
            self.df.loc[:,'Sequence'] = list(self.df['Sequence_Phospho'].values)
        else:
            self.df['Sequence'] = list(self.df['Sequence_Phospho'].values)

        for index in self.df.index:
            seq = self.df.ix[index, 'Sequence']
            for tag in tags:
                seq = seq.replace(tag, "")
            if "(" in seq:
                raise ValueError("Found invalid pattern in the sequence: %s" % seq)
            self.df.ix[index, 'Sequence'] = seq
        self.debug("set_sequence done")

    def _get_sequences(self):
        return self.df.Sequence
    sequences = property(_get_sequences, doc="returns list of Sequences found in the DataFrame")

    def _get_psites(self):
        return self.df.Psite
    psites = property(_get_psites, doc="returns list of Psites found in the DataFrame")

    def _rebuild_identifier(self):
        self.logging.warning("Rebuilding identifier in the dataframe. MERGED prefixes will be lost")

        self.merged = self.df.Sequence_Phospho.apply(lambda x: "MERGED" in x)

        if "Identifier" in self.df.columns:
            del self.df['Identifier']

        if "Identifier" not in self.df.columns:
            self.df['Identifier'] = self.df.Protein + "_" + self.df.Psite

        #for k,v in self.merged.iterkv():
        #    if v == True:
        #        self.df.loc[k,'Identifier'] = "MERGED_" + self.df.loc[k,'Identifier']

        if len(set(self.df.Identifier)) != len(self.df.Identifier):
            self.logging.warning("Identifiers are not unique. Have you called merge_peptides() ?")

    def _check_sequence_name_header(self):

        if "Sequence_Phospho" not in self.df.columns and "Sequence" not in self.df.columns:
            raise ValueError("Expecting either Sequence of Sequence_Phospho in the header")

        # maybe sequence columns is not named properly:
        if "Sequence_Phospho" not in self.df.columns:
            if sum(self.df.Sequence.apply(lambda x:"Phospho" in x)) == 0:
                self._build_sequence_phospho()
            else:
                self.logging.warning("Some Phospho strings found in Sequence column. No Sequence_Phospho column found.Renaming Sequence into Sequence_Phospho")
                self.df.columns = [x if x!="Sequence" else "Sequence_Phospho" for x in self.df.columns]

        # if sequence_phospho privded, all sequence must have at least one phospho
        if "Sequence_Phospho" in self.df.columns:
            if sum(self.df.Sequence_Phospho.apply(lambda x:"Phospho" in x)) != len(self.df):
                self.logging.warning("Some Sequence in Sequence_Phospho column do not have nay phosphos.")

    def _build_sequence_phospho(self):
        """

        .. warning:: not implemented, for now just a copy of Sequence into Sequence_Phospho,
            but should be rebuilt from the protein Sequence and psites location.
        """
        self.df['Sequence_Phospho'] = self.df.Sequence

    def _interpret_header(self):
        """Interpret the header of the CSV file.

        #. Check that the Sequence_Phospho column and/or Sequence column are available
        #. Check that the expected column from :attr:`header_mapping` are present
        #. Set Sequence_Phospho and Sequence properly if this is not the case (see :meth:`_set_sequence`)
        #. Finally check that the header is valid

        """
        self._check_sequence_name_header()


        columns = self.df.columns
        # rename columns if needed
        for valid, invalids in self.header_mapping.iteritems():
            for invalid in invalids:
                if valid not in self.df.columns and invalid in columns:
                    txt = "Column *{}* not found. Please use the correct header. Trying to replace possible match called **{}**"
                    self.logging.warning(txt.format(valid, invalid))
                    columns = [valid if c==invalid else c for c in columns]
        self.df.columns = columns


        self._set_sequence()


        self._check_header()


    def _check_header(self):
        assert "Protein" in self.df.columns, "Protein must be found in header"
        assert "Psite" in self.df.columns, "Psite must be found in header"
        assert "Sequence_Phospho" in self.df.columns or "Sequence" in self.df.columns, \
            "either Sequence or Sequence_Phospho must be provided"

    def boxplot(self, index=None, logy=True):
        """A simple boxplot of the data using logy
     
        :param index: you can rearrange the ordering of the axis column or select a subset
            by providing a list of indices.
        :param logy: set y-axis with log scale
     
        .. plot::
    
            from msdas import *
            m = MassSpecReader(get_yeast_small_data())
            m.boxplot()
   
  
        """
        if index == None:
            self.df.boxplot(rot=90, return_type='axes')
        else:
            self.df[index].boxplot(rot=90, return_type='axes')

        if logy:
            pylab.semilogy()
 


    def cleanup(self):
        """Call functions that clean up the dataframe

        Here are the steps followed:

        #. set zeros to NA :meth:`msdas.readers.Cleaner.set_zero_to_na`
        #. Rename psites using MSDAS convention (see :meth:`rename_psites`).
           Here just get rid of spaces
        #. calls :meth:`remove_ambiguities`
        #. calls :meth:`merge_peptides` if :attr:`_merge_peptides` is True
        #. rebuild the identifier based on Protein and Psite columns.

        """
        #self.debugLevel = "ERROR"
        if self.df is None:
            self.logging.warning("No data loaded yet. Please populate the attribute *df*, or call read_csv")

        # psites seprated by spaces are renamed so that there are now seprated
        # by the ^ (and ) character
        # oxidation psite are also interpreted.
        self.logging.info("Renaming psites with ^ character")
        self.rename_psites()

        # get rid of zeros in the data
        # should be used for raw data, not normalised so should be called by user
        self.logging.info("Replacing zeros with NAs")
        self.set_zero_to_na()

        # see documentatoion of the method itself
        self.remove_ambiguities()

        if self._merge_peptides:
            self.merge_peptides()

        self._rebuild_identifier()

    def merge_peptides(self):
        """Merge data using :class:`Cleaner` class

        #. Merge peptide with same psites. See :meth:`~msdas.readers.Cleaner.merge_peptide_with_same_psites`
        #. Merge identical peptide sequences. See :meth:`~msdas.readers.Cleaner.merge_identical_peptide_sequence`

        Note that in the yeast case, we use only the merge_identical_peptide_sequence

        """
        self.merge_peptides_with_same_psites()
        self.merge_identical_peptide_sequence()

    def read_csv(self, filename=None, quotechar=None, index_col=None, sep=None):
        """Read a CSV file file Mass Spec data

        :param str filename: if not provided, uses the class attribute
            :attr:`filename`


        This function is used by the constructor of MassSpecReader if a filename
        is provided. Yet, you can either create an empty instance and read the file
        later on. If you read the file later, the difference with the constructor is
        that by default the :meth:`cleanup` function is called. So these two codes
        are equivalent only if you call :meth:`cleanup`::

            from msdas import *
            m1 = MassSpecReader(get_yeast_small_data())

        ::

            from msdas import *
            m2 = MassSpecReader()
            m2.read_csv(get_yeast_small_data())
            m2.cleanup()

        You can then check that the 2 files are equivalent::

            >>> m1 == m2
            True


        """
        self.level = self.level
        if self.mode is None:
            raise ValueError("mode is not set. mode must be provided as YEAST")
        self.logging.info("Reading %s" % filename)

        # Read the data. Must work
        if quotechar == None:
            quotechar = self._quotechar

        #print index_col, sep, quotechar, self.mode
        #sep=None
        df = pd.read_csv(filename, quotechar=quotechar, index_col=index_col, sep=sep, engine='python')
        df.dropna(how="all", inplace=True) # ignore empty lines
        
        # replaces None by NA
        import numpy as np
        df.fillna(np.NaN, inplace=True)

        # spaces need to be removed in the header
        df.columns = [col.strip() for col in df.columns]

        self.df = df
        self.rawdf = self.df.copy() # for bookkeeping

        # and columns that are made of strings should also be stripped
        def clean_func(x):
            try:
                return x.strip()
            except:
                return x
        for col in df.columns:
            self.df[col] = self.df[col].map(clean_func)

        self._interpret_header()

    def remove_ambiguities(self):
        """

        #. Remove rows where psites are ambiguous or could not be interpreted.
           See :meth:`remove_ambiguous_psites`
        #. Remove rows where number of proteins is > 2 :meth:`remove_ambiguous_proteins`

        """
        # remove psites where number of psites does not match phospho+oxidation
        # rename_psites (append_oxidation) must be called before otherwise all
        # sequences with oxidation are removed. converns 203 rows in yeast
        self.remove_ambiguous_psites()

        # not redundant with drop above. remove rows where protein names are > 1
        self.remove_ambiguous_proteins()


    def rename_psites(self):
        """Rename Psites using AND/OR convention.


        To store the Psite column in the header, we use the following conventions:

        * individual psite are separated by AND using the **^** character
        * ambiguous psite when position is not certain are separated by OR that is **+** character

        See :class:`msdas.psites.PSites` class documentation for details.

        Here, the CSV file may contain Psites stored with spaces, which are considered
        as ANDs. So, we first remove spaces, replace them by ANDs. There are also
        extra ; characters that are ignored. If there are duplicated psites, which occurs
        in the raw data set (yeast case), we remove them. Finally, in the yeast case,
        some sequence contain the (Oxidation) string but the corresponding Psite is not
        provided. This function extract this information and stores it in the Psite string.

        Psites are also sorted by position e.g., S5^S3 is renamed as S3^S5.
        """
        ptools = PSites(verbose=self.level) # important to provide level

        # removes spaces
        psites = [ptools.remove_spaces(p) for p in self.df.Psite]

        # add AND character
        psites = [self.AND.join(psite.split()) for psite in psites]

        # what to do with ; character? assume an AND for now
        psites = [x.replace(";", self.AND) for x in psites]

        psites = [ptools.remove_duplicated(x) for x in psites]
        self.df['Psite'] = psites

    def remove_ambiguous_proteins(self):
        """Remove rows where proteins are not uniquely defined

        Some rows may have the protein field incorrectly defined with several
        protein names, in which case the row is discarded by this function


        """
        # Removes all rows where we cannot understand what is going on.
        # e.g. protein or sites names are not unique
        toremove = self.df.index[self.df.Protein.apply(lambda x: " " in x)]
        self.logging.info("-- Removing {} rows with ambigous protein names:".format(len(toremove)))
        self._removed_ambiguous_proteins = toremove[:]
        #FIXME: here we use self.df.ix whereas in remove_ambigous_psites, we use self.df directly why?
        self.removed['ambiguous_proteins'] = self.df.ix[self._removed_ambiguous_proteins].copy()

        self._ambiguous_proteins_df = self.df.ix[toremove==True][['Protein', 'Sequence', 'Psite', 'Sequence_Phospho']].copy()

        for index in toremove:
            self.logging.debug(" {}:{}".format(index+1,self.df.Protein[index]))
        self.logging.info("--------------------------------------------------")
        self._info["N_unknown"] = len(toremove)
        self.df.drop(toremove, inplace=True)

    def remove_ambiguous_psites(self):
        """Removes Psites where number of Psite does not match content of Sequence

        Some sequences are different number of phospho as compared to the Psites
        columns. In such case, we remove the rows. The field used for the sequence
        is the Sequence_Phospho where psites are encoded with the (Phsospho) string

        ::

          ALY1,    S118^S216,   TPLPSSS(Phospho)R

        .. note:: Oxidation are ignored

        """
        if "Sequence_Phospho" not in self.df.columns:
            return
        # count number of phospho
        # IF MERGED, should be skipped
        l3 = self.df.Sequence_Phospho.apply(lambda x: "MERGED" not in x)

        count_phosphos = self.df.Sequence_Phospho.apply(lambda x: x.count("Phospho"))
        #count_oxis = self.df.Sequence_original.apply(lambda x: x.count("Oxidation"))
        count_psites = self.df.Psite.apply(lambda x: len(x.split(self.AND)))

        #toremove = count_phosphos+count_oxis != count_psites
        toremove = (count_phosphos != count_psites) & l3

        self._removed_ambiguous_psites = toremove[:]

        self.removed['ambiguous_psites'] = self.df.ix[self._removed_ambiguous_psites].copy()

        if sum(toremove) > 0:
            self._ambiguous_psites_df = self.df.ix[toremove==True][['Protein', 'Sequence', 'Psite', 'Sequence_Phospho']].copy()
            self.logging.info("-- %s rows have ambiguous psites and are removed" % sum(toremove))
            self.logging.info("save data in attribute _ambiguous_psites_df")
            self.df.drop(toremove.index[toremove==True], inplace=True)
            self.logging.info("--------------------------------------------------")

        #self._info["N_ambiguous_psites"] = sum(toremove==True)

    def _rename_measurements(self, tag):
        """Append prefix to all measurement columns


        Measurements columns are those that starts with the letter "t" for now.

        """
        oldcols = self.df.columns
        columns = [tag+"_" + c if c.startswith("t") else c for c in oldcols]
        self.df.columns = columns

    def to_csv(self, filename):
        """Save dataframe into a CSV file

        """
        #FIXME: why a copy ??
        df = self.df.copy()

        #df['Sequence'] = df['Sequence_Phospho']
        #del df['Sequence_original']
        df.to_csv(filename, index=False, sep=",")

    def sort_psites_ors_only(self):
        """Sort the psites location in psites string that contains OR, only

        e.g., S3+S1 is renamed as S1+S3 but S3+S1^S7 is unchanged.

        """
        func = PSites(verbose=False).sort_psites_ors_only
        self.df['Psite'] = self.df.Psite.apply(lambda x: func(x))
        self._rebuild_identifier()

    def plot_phospho_stats(self, **kargs):
        """Plots descriptive statistics about the Psites

        The data used comes from the DataFrame.

        There are 3 plots created:

            #. pie charts of number of sites of each category (e.g., S, T, ..)
            #. Histogram of number of phospho sites per protein
            #. Histogram of number of phospho sites per peptides

        .. plot::
            :include-source:
            :width: 45%

            from msdas import *
            filenames = get_yeast_filenames()
            p1 = readers.MassSpecReader(filenames[0], mode="YEAST")
            p1.plot_phospho_stats()

        """
        self._global_count = {}
        self._number_sites_per_peptide = []
        self._number_sites_per_protein = {}
        for sites in self.df.Psite:

            # then split on "-" in case there are ambiguities and select the
            # first letter in any case.
            sites = [this[0] for this in sites.split(self.AND)]
            for s in sites:
                if s not in self._global_count.keys():
                    self._global_count[s] = 0
                else:
                    self._global_count[s] += 1
            # how many phospho sites per peptide ?
            self._number_sites_per_peptide.append(len(sites))

        pylab.figure(1)
        pylab.clf()
        pylab.pie(self._global_count.values(), labels=[k+" ("+str(v)+")" for k,v in
            self._global_count.iteritems()])
        pylab.title("Type of of phospho sites (over all peptides)")

        pylab.figure(2)
        pylab.clf()
        pylab.hist(self._number_sites_per_peptide, [x+0.5 for x in
            range(0,1+max(self._number_sites_per_peptide))])
        pylab.title("Number of phospho sites per peptide")
        pylab.grid()

        g = self.df[['Protein', 'Psite']].groupby("Protein")
        groups = g.Psite.groups

        pfact = PSites()

        for group in groups:
            sites = [x for x in self.df.Psite[g.Psite.groups[group]]]
            #sites = list(pylab.flatten([[y[0] for y in x.split("_")[1:]] for x in sites]))
            # FIXME: if merge or not, the results will be different !
            sites = pfact.get_unique_psites(sites)
            self._number_sites_per_protein[group] = len(sites)

        pylab.figure(3)
        pylab.clf()
        bins = [x+0.5 for x in
            range(0,1+max(self._number_sites_per_protein.values()))]
        pylab.hist(self._number_sites_per_protein.values(), bins)
        #pylab.xticks(bins, )
        pylab.title("Number of phospho sites per protein")
        pylab.grid()


    def hist_number_peptides(self):
        """Create a bar plot with number of peptides per protein


        """
        groups = [(len(v),k) for k,v in self.df.groupby("Protein").groups.iteritems()]
        l = sorted(groups)
        l.reverse()
        y = [this[0] for this in l]
        xlabel = [this[1] for this in l]
        pylab.clf()
        pylab.bar(range(0,len(y)), y)
        pylab.xticks([this+.5 for this in range(0,len(y))], xlabel, rotation=90)
        for i in range(0,len(y)):
            pylab.text(i+0.5, y[i]+.5, y[i])
        pylab.ylabel("Number of peptide per protein")
        pylab.ylim([0,max(y)+1])


    def hist_peptide_sequence_length(self, drop_duplicates=True, **kargs):
        """Plot histogram of the peptide sequences length.

        :param bool drop_duplicates: duplicated peptide sequence are dropped by default
            Indeed, you may have two identical peptides with different psites. Since we are
            interested in the peptide sequence, we drop duplicated ones.
        :param kargs: optional keywords are passed to the dataframe.hist method.

        .. plot::

            >>> from msdas import *
            >>> m = MassSpecReader(get_yeast_filenames()[0], mode="YEAST")
            >>> m.hist_peptide_sequence_length()

        """
        if drop_duplicates:
            data = self.df.Sequence.drop_duplicates().apply(lambda x: len(x))
        else:
            data = self.df.Sequence.apply(lambda x: len(x))
        data.hist(**kargs)
        pylab.title("")
        pylab.xlabel("Sequence length")
        pylab.ylabel("#")

    def pcolor(self, protein=None, tag="t", fillna=None, vmax=None, vmin=None,
               cmap="hot_r", fontsize_yticks=8):
        """

        :param str protein: a protein name. If not provided, plot all the data.
        :param str tag: only column that starts with this tag are used
        :param fillna: if na are found, there are plot as grey pixels. You can set
            all NAs to a scalar if required.
        :param float vmin: vmin parameters used by the matplotlib.pcolor function
        :param float vmax: vmax parameters used by the matplotlib.pcolor function
        :param cmap: a valid color map from matplotlib

        .. plot::

            from msdas import *
            m = MassSpecReader(get_yeast_small_data(), verbose=False)
            m.pcolor("DIG1", tag="a") # pick up measurements with name starting with "a" letter


        """
        measures = [c for c in self.df.columns if c.startswith(tag)]
        if len(measures)==0:
            print("no column matches the provided tag (%s)"%tag)
            return
        pylab.clf()

        if protein:
            data = self[protein][measures]
            psites = list(self.df.Psite[self[protein].index])
        else:
            data = self.df[measures]
            psites = list(self.df.Psite)
            proteins = list(self.df.Protein)
            psites = [x+"_"+y for x,y in zip(proteins, psites)]

        if vmax == None:
            vmax = data.max().max()
        if vmin == None:
            vmin = data.min().min()
        mask = np.ma.array(data.as_matrix(), mask=data.isnull())

        cmap = pylab.get_cmap(cmap)
        cmap.set_bad("grey", 1)
        pylab.pcolormesh(mask, vmin=vmin, vmax=vmax, cmap=cmap)
        pylab.colorbar()
        pylab.xlabel("measurements")
        pylab.ylabel("peptides")
        pylab.xlim([0, len(measures)])
        pylab.ylim([0, data.shape[0]])
        pylab.yticks([x+0.5 for x in range(0,data.shape[0])], psites,
                      fontsize=fontsize_yticks)
        if protein:
            pylab.title(protein)
        pylab.tight_layout()
        return data

    def read_annotations(self, filename):
        """Read annotations saved in a dataframe.

        :param str filename: filename where annotations are saved (pickle file)

        See :class:`msdas.annotations.Annotations`

        This function is required to populate the annotatoins attribute used for example
        by :meth:`show_sequence`

        """
        self.annotations = pd.read_pickle(filename)

    def status(self):
        """Some basic sanity checks on the dataframe that contains the data"""
        # check header
        for this in self.valid_header:
            if this not in self.df.columns:
                print("Could not find %s in the dataframe columns" % this)

        # check Psites
        import psites
        pfactory = psites.PSites()
        if "Psite" in self.df.columns:
            for k,v in self.df.Psite.iteritems():
                if pfactory.isvalid(v) == False:
                    print("Invalid psite (%s) at index %s" % (k,v))
                    print("There are maybe others. break for now the psite checking")
                    break

        return True

    def plot_timeseries(self, identifier, replicate=None, **kargs):
        """This is for the yeast case with replicates

        :param str identifier:
        :param replicates: None means use the first replicate (the one without prefix).
            otherwise, provide the suffix corresponding to the replicate (see example below)
        :param kargs: any paraemter accepted by the matplotlib.plot function

        .. plot::
            :include-source:
            :width: 80%


            from msdas import *
            m = readers.MassSpecReader(get_yeast_small_data(), verbose=False)
            m.plot_timeseries("DIG1_S142")

        If you have replicates, only the first replicate is taken into account and
        you therefore need to call plot_timeseries several time providing the
        tag of the replicate. The first replicate has no tag. The second replicate
        has a suffix (**.1**), the third replicates has the tag **.2** and so on

        .. plot::
            :include-source:
            :width: 80%

            from msdas import *
            M = readers.MassSpecReader(get_yeast_raw_data(), verbose=False)
            M.plot_timeseries("DIG1_S142", replicate=None, color="r")
            M.plot_timeseries("DIG1_S142", replicate=1, color="g")
            M.plot_timeseries("DIG1_S142", replicate=2, color="b")

        """
        if replicate == None:
            measures = [x for x in self.measurements if "." not in x]
        else:
            measures = [x for x in self.measurements if x.endswith("."+str(replicate))]
        if len(measures)==0:
            raise ValueError("No measurement data found that correspond to the replicate provided" )

        try:
            # for the replicates; get_mu_df is required
            data = self.get_mu_df()[measures]
        except:
            data= self.df[measures]

        indices = self.df[self.df.Identifier == identifier].index
        data = data.ix[indices]
        kargs['color'] = kargs.get("color", "r")
        kargs['marker'] = kargs.get("marker", "o")
        kargs['markersize'] = kargs.get("markersize", 10)
        kargs['lw'] = kargs.get("lw", "2")
        if len(indices)>1:
            self.logging.warning("More than 1 row found. Consider calling merging_ambiguous_peptides method")
        for index in indices:
            pylab.plot(data.ix[index], '--', **kargs)

        N = len(measures)
        pylab.title(identifier)
        pylab.xticks(range(0,N), measures, rotation=90)
        pylab.grid()
        return data

    def get_data_matrix(self, identifier):
        """YEAST case with replicates

        .. todo:: to be moved to yeast module
        """
        alpha = ['a0', 'a1', 'a5', 'a10', 'a20', 'a45']
        salt = ['t0', 't1', 't5', 't10', 't20', 't45']

        try:
            data = self.get_mu_df()
            data = data[self.df.Identifier == identifier]
            if len(data)>1:
                self.logging.error("Found more than 1 row matching the psite")
            print('a')
        except:
            data = self.df[self.df.Identifier == identifier]
        exps = np.zeros((6,6))

        for i,a in enumerate(alpha):
            for j,s in enumerate(salt):
                exps[i,j] = data[a+"_"+s]

        exps = pd.DataFrame(exps)
        #TODO : is this the correct order
        exps.index = alpha
        exps.column = salt
        return exps

    def plot_experiments(self, identifier, cmap="hot_r", vmin=None,vmax=None):
        """

        should be for yeast only. Put inside a YEAST class ?

        if not merged, you could end up with several rows.
        may fail in such case...


        """
        #        fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
        #ax.zaxis.set_major_locator(LinearLocator(10))
        #ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        #fig.colorbar(surf, shrink=0.5, aspect=5)

        self.logging.warning("Works with yeast data set only")

        exps = self.get_data_matrix(identifier)
        mask = np.ma.array(exps, mask=np.isnan(exps))
        pylab.clf()
        if vmax == None:
            vmax = np.nanmax(exps)
        if vmin == None:
            vmin = np.nanmin(exps)
        cmap = pylab.get_cmap(cmap)
        cmap.set_bad("grey", 1)

        pylab.pcolormesh(pylab.flipud(mask), vmin=vmin, vmax=vmax, cmap=cmap)
        alpha = ['a0', 'a1', 'a5', 'a10', 'a20', 'a45']
        salt = ['t0', 't1', 't5', 't10', 't20', 't45']
        pylab.yticks([x+0.5 for x in range(0, 6)], alpha[::-1])
        pylab.xticks([x+0.5 for x in range(0, 6)], salt)
        pylab.colorbar()
        pylab.title(identifier)
        pylab.ylabel("Time after alpha-factor stimulation")
        pylab.xlabel("Time after NaCl stimulation")

        return exps

    def show_sequence(self, protein, N=40, vmax=None,cmap="hot_r", fontsize=10):
        """Show peptide and phosphorylation distribution on a sequence

        .. plot::
            :include-source:

            import os
            from msdas import *
            r = readers.MassSpecReader(get_yeast_raw_data())
            if os.path.exists("YEAST_annotations.pkl")==False:
                a = annotations.Annotations(r, "YEAST")
                a.get_uniprot_entries()
                a.set_annotations()
                a.to_pickle()
            r.read_annotations("YEAST_annotations.pkl")
            r.show_sequence("DIG1",fontsize=8)

        .. warning:: once merged, some phospho are tagged ith merged name so we use the colum Sequence (not Sequence_Phospho)
        """
        import textwrap
        self.warning("Requires annotations")
        indices = self.metadata.groupby("Protein").groups[protein]
        psites = self.metadata.Psite.ix[indices]
        peptide_sequences = self.metadata.Sequence.ix[indices]
        entry = self.metadata.Entry.ix[indices[0]]
        sequence = self.annotations.ix[entry].Sequence

        psites_location = set([w[1:] for x in psites for z in x.split("^") for w in z.split("+")])
        psites_location = [int(x) for x in psites_location]
        # build the matrix to plot. The matrix will have the following length and text
        text = textwrap.wrap(sequence, N)
        ncol = N
        nrow = len(text)
        # but for now, let us work with the data as an array. We will want a matrix,
        # so let us complete the sequence with zeros such that length is the same as a
        # matrix
        values = [0] * ncol * nrow

        # let us look at peptide coverage
        for peptide in peptide_sequences:
            i = sequence.index(peptide)
            j = len(peptide)
            for pos in range(i, i+j):
                values[pos] += 1
        if vmax == None:
            M = max(values)
        else:
            M = vmax
        # now we can build the matrix
        m = np.array(values).reshape(nrow, ncol)
        m = pylab.flipud(m)

        pylab.clf();
        pylab.pcolor(m, edgecolors="k", vmin=0, vmax=M, cmap=cmap, alpha=0.8);
        pylab.colorbar()
        for i,x in enumerate(sequence):
            if i+1 in psites_location:
                pylab.text(i%N+.2, -0.5+nrow-i/N, x, color="green",
                           bbox=dict(facecolor='white', alpha=0.5), fontsize=fontsize)
            else:
                pylab.text(i%N+.2, -0.5+nrow-i/N, x, fontsize=fontsize)

            pylab.text(i%N+.6, -.9+nrow-i/N, i+1, fontsize=fontsize/1.5)
        pylab.title(protein)
        pylab.xticks([],[])
        pylab.yticks([],[])
        return text

    def __eq__(self, this):
        metadata = self.metadata.fillna(0).eq(this.metadata.fillna(0)).all().all()
        if metadata == False:
            return False
        measures1 = self.measurements.fillna(0).apply(np.round, args=[10]) # 10 decimals
        measures2 = this.measurements.fillna(0).apply(np.round, args=[10])
        return measures1.eq(measures2).all().all()

    def __getitem__(self, name):
        if isinstance(name, str):
            name = [name]
        if len(name)==1:
            name = name[0]
            df = self.df[self.df.Protein == name].copy()
            if len(df)==0:
                # try someting else based on identifier instead of protein name
                df =  self.df[self.df.Identifier == name].copy()
        elif len(name) == 2:
            protein = name[0]
            psite= name[1]
            df = self.df[(self.df.Protein == protein) & (self.df.Psite == psite)]
        else:
            raise ValueError("Expects only 1 or 2 parameters.")

        return df

    def __str__(self):
        import textwrap
        N = len([x for x in self.df.columns if x not in ['Sequence', 'Psite', 'Protein']])
        msg = "This dataframe contains {} columns in addition to the ".format(N) + \
            "standard columns Protein, Sequence, Psite"
        txt = "\n".join(textwrap.wrap(msg, 80))
        txt += "\nYour data contains {} unique proteins\n".format(len(set(self.df.Protein)))
        txt += "Your data contains {} combination of psites/proteins\n".format(self.df.shape[0])
        return txt
