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
"""Module to align several files together.

The contents of the files in term of protein name and peptide sequence can
overlap or not. Each file is taken as a different experiment and therefore column
names that are not standard (Protein, Sequence, Entry, ...) are renamed with a
prefix (the filenames or prefixes provided by the user.)


"""
import os

import pylab
from easydev import Logging
from msdas.readers import MassSpecReader
from msdas.readers import PSites
from msdas.tools import SequenceTools


__all__ = ["MassSpecAlignmentYeast", "MassSpecAlignmentBase"]



# ALIGNEMENT not merging
class MassSpecAlignmentBase(Logging):
    """Base Class related to Reading/Writing Mass Spectrometry data sets

    This class gathers common functionalities for more specialised
    classes :class:`MassSpecAlignmentYeast`.

    Inputs can be a filename of a instance of MassSpecReader::

        from msdas import *
        m = MassSpecAlignmentBase()
        r1 = MassSpecReader(get_yeast_filenames()[0])
        r2 = MassSpecReader(get_yeast_filenames()[1])
        df = m.merge([r1.df,r2.df])

    The contents of r1 and r2 is similar (same column names). So columns must be
    renamed (except for the Proiten, Sequence and so on i.e. metadata).


    .. seealso:: :class:`MassSpecAlignment`, :class:`MassSpecAlignmentYeast`,
       :class:`MassSpecAlignmentTCell`.

    """
    AND = "^"
    def __init__(self, filenames=None, verbose=True,
                 ):
        """.. rubric:: Constructor

        :param str filename: Not required but used by sub classes.
        :param bool verbose: verbosity set on by default


        """
        super(MassSpecAlignmentBase, self).__init__(level=verbose)

        # create this class to help us reading single CSV later on
        self.ms_reader = MassSpecReader(verbose=verbose)

        if isinstance(filenames, list):
            self.filenames = filenames[:]
        if isinstance(filenames, str):
            self.filenames = [filenames]
        else:
            self.filenames = filenames

        #: the dataframe, set to None by default
        self.df = None


        self._error_messages = {
                'df': "df attribute not found. Please call " +
                      "read(filename) method or use proper MassSpecAlignment class",
                }

        #: list of column names to performed the merge on 'Protein', "Sequence", 'Psite', 'Sequence_Phospho'
        self.merge_on = ['Protein', "Sequence",'Psite', 'Sequence_Phospho']

    def _get_mode(self):
        return self.ms_reader.mode
    def _set_mode(self, mode):
        self.ms_reader.mode = mode
    mode = property(_get_mode, _set_mode, doc="get/set mode of the MassSpecReader")

    def check_format(self):
        ptools = PSites(verbose=self.level)
        for psite in self.df.Psite:
            if ptools.isvalid(psite) == False:
                raise ValueError("found invalid psite %s" % psite)

    def merge(self, dfs, on=None):
        """Merge/align several dataframes

        :param list dfs: list of dataframes to align
        :param list on: list of column to perform the merge on. If not
            provided, use the :attr:`merge_on` attribute (Protein, Sequence,
            Psite, Sequence_Phospho).

        .. note:: This is a merge using the merge function from pandas library.
            the merge is performed using the **how** parameter set to **outer**.

        """
        if on == None:
            on = self.merge_on[:]

        if len(dfs)==0:
            raise ValueError("No data provided. Please provide a list of dataframe")

        if len(dfs) == 1:
            # nothing to do
            return dfs[0].copy()

        df = dfs[0].copy()

        # length may be different
        for index, thisdf in enumerate(dfs[1:]):
            df = df.merge(thisdf, how="outer", on=on)

        # rename columns so that all Xt are at the end
        df = df[on + [c for c in df.columns if c not in on]]
        return df

    def __str__(self):
        txt = "This is MassSpecAlignment instance\n"
        txt += "It contains %s combination of protein names and psites\n" % len(self.df.Protein)
        txt += "Data is contained in the attribute df (dataframe)"
        return txt


class MassSpecAlignmentTCell(MassSpecAlignmentBase, SequenceTools):
    """Align several Mass Spectrometry data files.

    See :class:`MassSpecAlignmentBase` for the list of columns to be found to
    perform the alignment/merging and :class:`MassSpecAlignmentYeast` for more
    details

    MS-DAS provides a set of 3 files to play with::

        from msdas import *
        m = MassSpecAlignmentTCell(get_tcell_filenames(),

    .. todo:: check specific cases such as Q8IYB3-2 isoform 2 why the isoform 1 picked up

    The original files contains extra information such as Retention_time__min_,
    Charge, measured mass that are kept in rawdf attribute.

    .. warning:: not fully tested

    """
    def __init__(self, filenames=None, verbose=True, prefixes=None):
        """.. rubric:: contructor

        :param str filename: See documentation of the class for details about
            the format
        :param bool verbose: verbosity set on by default
        """
        super(MassSpecAlignmentTCell, self).__init__(filenames, verbose=verbose)
        self.prefixes = prefixes
        if self.filenames:
            if prefixes == None:
                self.logging.warning("No prefixes provided. We will use the filename")
                prefixes = [os.path.splitext(os.path.split(f)[-1])[0] for f in self.filenames]
                self.prefixes =prefixes
            if len(prefixes)!=len(self.filenames):
                raise ValueError("Number of prefixes must match number of filenames")

            self._init_tcell(prefixes)

        self.logging.error("MassSpecAlignmentTCell has not been tested thoroughly. Use with care")

    def _update_df_tcell(self):

        # we select only some columns
        selected_columns = [x for x in self.rawdf.columns if x.startswith("Nor") or
            x.startswith("Raw") or x in ["Max_fold_change","Psite",
            "Protein", "ProteinID", "Unstimulated", "Sequence"]]
        self.df = self.rawdf.ix[:, selected_columns].copy()


    def _read_tcell_csv(self, filename, tag):
        # Use an external reader that will gather all kind of conventions
        self.ms_reader.read_csv(filename,  sep=",")
        # if you merge the single file, if may not work. to be tested
        self.ms_reader.rename_psites()
        self.ms_reader.set_zero_to_na()
        # renaming columns
        self.ms_reader._rename_measurements(tag)

        return self.ms_reader.df

    def _init_tcell(self, prefixes):
        # read first one
        self._dfs = []

        df = self._read_tcell_csv(self.filenames[0], tag=prefixes[0])
        self._dfs.append(df)

        # read other filename and merge into original data frame
        if self.filenames>1:
            for i, filename in enumerate(self.filenames[1:]):
                thisdf = self._read_tcell_csv(filename, tag=prefixes[i+1])
                # if we set names in read_csv, floats are converted to strings.
                # so we set the column names afterwards.
                self._dfs.append(thisdf)

            self.df = self.merge(self._dfs).copy()





class MassSpecAlignmentYeast(MassSpecAlignmentBase):
    """Align several Mass Spectrometry data files.

    MS-DAS provides a set of 6 files to play with. Their paths can be obtained
    using function :func:`get_yeast_filenames`::

        from msdas import *
        filenames = alignment.get_yeast_filenames()
        yeast = alignment.MassSpecAlignmentYeast(filenames)

    See :class:`MassSpecAlignmentBase` for the list of columns to be found to
    perform the alignment/merging.

    Since measurements may have the same names from one file to another, then need
    to be renamed.

    Internallty, a prefix is added. It is populated with prefixes provided by
    the use; otherwise with the filenames themselves.

    ::

        from msdas import *
        m = MassSpecAlignmentYeast(get_yeast_filenames(),
                            prefixes = ["a0", "a1", "a5", "a10" "a20", "a45])

    You can then create a MassSpecReader instance from the dataframe just
    created and save it to a file                              ::

        r = MassSpecReader(m)
        r.to_csv("test_align.csv")

    """
    def __init__(self, filenames=None, verbose=True,  prefixes=None):
        """.. rubric:: constructor

        :param list filename: list of filenames (readable by MassSpecReader).
        :param bool verbose: verbosity set on by default

        """
        super(MassSpecAlignmentYeast, self).__init__(filenames, verbose=verbose)
        self.mode = "YEAST"

        self.prefixes = prefixes
        if self.filenames:
            if prefixes == None:
                self.logging.warning("No prefixes provided. We will use the filename")
                prefixes = [os.path.splitext(os.path.split(f)[-1])[0] for f in self.filenames]
            if len(prefixes)!=len(self.filenames):
                raise ValueError("Number of prefixes must match number of filenames")

            self._init_yeast(prefixes)

    def _read_yeast_csv(self, filename, tag):
        # Use an external reader that will gather all kind of conventions
        self.ms_reader.read_csv(filename,  sep=",")
        # if you merge the single file, if may not work. to be tested
        self.ms_reader.rename_psites()
        self.ms_reader.set_zero_to_na()
        # renaming columns
        self.ms_reader._rename_measurements(tag)

        return self.ms_reader.df

    def _init_yeast(self, prefixes):
        # read first one
        self._dfs = []

        df = self._read_yeast_csv(self.filenames[0], tag=prefixes[0])
        self._dfs.append(df)

        # read other filename and merge into original data frame
        if self.filenames>1:
            for i, filename in enumerate(self.filenames[1:]):
                thisdf = self._read_yeast_csv(filename, tag=prefixes[i+1])
                # if we set names in read_csv, floats are converted to strings.
                # so we set the column names afterwards.
                self._dfs.append(thisdf)

            self.df = self.merge(self._dfs).copy()




