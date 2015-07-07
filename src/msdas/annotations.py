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
import os
import textwrap
import types

import pandas as pd
import pylab
from bioservices import UniProt

from msdas.readers import MassSpecReader



__all__ = ["Annotations", "AnnotationsYeast"]


class Annotations(MassSpecReader):
    """Create/store/read annotations from uniprot and figure out entry names

    The Annotations classes allows one to populate the dataframe attribute
    :attr:`df` with the **Entry** and **Entry_name**  columns (UniProt entries).
    This is not strictly speaking required columns but provide more tools if
    available. The annotations also creates a new dataframe called :attr:`annotations`
    that stores in particular the protein sequence and the GO terms. The former
    being used to check the peptide sequence and the latter to plot relevant
    histogram about GO terms.

    This class inherits from :class:`msdas.readers.MassSpecReader`. Consequently,
    input can be a MassSpecReader instance, or a filename or even nothing (data
    can be read at a later stage). The dataframe must contain the Protein column.

    One reason to fetch the entries from UniProt is that the protein column name
    may contain typos or non-uniprot entries, therefore it is quite useful to
    fetch all entries from uniprot based on the protein name provided. This can
    be done thanks to the :meth:`get_uniprot_entries`. This method fills a dictionary called
    :attr:`_mapping` (note the underscore), which is used to populate a new column
    in the dataframe called **Entry**.

    If your initial dataframe contains the columns "Entry" with all valid UniProt
    entries (e.g., P23300) then the :attr:`_mapping` attribute is populated during the
    initialisation and the call to :meth:`get_uniprot_entries` can be skipped.
    If called, it will also be faster but will overwrite the content of the Entry column.
    You can also fill/correct/complete the :attr:`_mapping` attribute before calling
    :meth:`get_uniprot_entries`

    .. doctest::

        >>> from msdas import annotations
        >>> import pandas as pd
        >>> df = pd.DataFrame({
            'Protein':['DIG1'],
            'Sequence_Phospho':['SS(Phospho)T'],
            'Psite':['S2']})
        >>> a = annotations.Annotations(df, "YEAST")
        >>> a._mapping
        {}
        >>> a.get_uniprot_entries()
        {'DIG1_YEAST': ['Q03063']}
        >>> a.df.Entry
        0    Q03063
        Name: Entry, dtype: object

    Then, call :meth:`set_annotations`, which will fetch all annotations from
    uniprot and store them in a new dataframe in the :attr:`annotations` attribute

    ::

        a.set_annotations()
        a.annotations

    A new field called **Entry_name** is also added to the dataframe itself::

        a.df.Entry_name

    On a big data set, it may take a few minutes to fetch all information from uniprot.
    So, we also provide tools to save and read back the relevant information (
    :meth:`read_annotations`, :meth:`to_pickle`, :meth:`to_csv`  ) ::

        from msdas import *
        r = readers.MassSpecReader(get_yeast_raw_data())
        # this takes about 10 minutes depending on the connection for 1600 unique protein names
        r.get_uniprot_entries()
        r.set_annotations()
        r.to_pickle(tag="test") # creates a file called YEAST_annotations_test.pkl
        r.to_csv("data.csv")

    Next time, just type::

        from msdas import *
        a = annotations.Annotations("data.csv", "YEAST")
        a.read_annotations("YEAST_annotations_test.pkl")

    To check that the entries are correct, one thing that can be done is to
    look for the peptide sequence into the FASTA sequence found in the annotations::

        a.check_entries_versus_sequence()

    This is a very good sanity check to verify that the entry names found
    correspond to the peptide provided. If not, the protein name was probably wrong
    or was a gene name that could not be mapped correctly to the correct protein.

    If some entries are not found or mapping was not found, you need to
    manually check the issues and update the :attr:`_mapping` attribute,
    update the uniprot entries and annotations::

        a._mapping[entry] = ['entry name']
        a.get_uniprot_entries()
        a.set_annotations()
        a.check_entries_versus_sequence()

    if you cannot find a mapping, we would recommend to delete the item from
    the dataframe :attr:`df`.
    """
    def __init__(self, data, organism=None, verbose=True, annotations=None, **kargs):
        """.. rubric:: Constructor

        :param data: a MassSpecReader compatible input (e.g., CSV file, None, a
            MassSpecReader instance). See :class:`msdas.readers.MassSpecReader`
            documentation for details
        :param organism: valid uniprot identifier for the organism e.g., HUMAN
            YEAST.
        :param annotations: a pickled file containing the annotations saved
            using :meth:`to_pickle`.
        :param kargs: valid parameter recognised by :class:`msdas.readers.MassSpecReader`

        """
        super(Annotations, self).__init__(data=data, verbose=verbose, **kargs)
        if organism is None:
            raise ValueError("organism must be provided e.g. YEAST, HUMAN")

        self.organism = organism

        #: the dataframe where annotations from uniprot will be stored.
        self.annotations = None

        self._mapping = {}
        self.build_mapping_from_df()# if Entry is provided

        self._init_uniprot()

        if annotations:
            self.read_pickle(annotations)

    def _init_uniprot(self):
        if hasattr(self, "_uniprot") == False:
            self._uniprot = UniProt(verbose=self.debugLevel)

    def _update_species_to_find(self):
        entry_names = [x + "_" + self.organism for x in self.df.Protein]
        #unique_entry_names = list(set(entry_names))

        species_to_find = [k for k in entry_names if k not in self._mapping.keys()]
        species_to_find = list(set(species_to_find))
        self._species_to_find = list(set(species_to_find))

    def build_mapping_from_df(self):
        """Populate the _mapping dictionary using the Uniprot Entry column"""
        if "Entry" in self.df.columns:
            for index in self.df.index:
                k = self.df.Protein.ix[index]
                if k.endswith("_"+self.organism) == False:
                    k += "_" + self.organism
                v = self.df.Entry.ix[index]
                self._mapping[k] = [v]
        else:
            self.warning("Entry column not found in the dataframe. call get_uniprot_entries")

    def get_uniprot_entries(self, Nmax=50):
        """Search for the uniprot entries and entry names given protein column.

        Protein names from the dataframe are first used to feed uniprot mapping tool.
        Some protein names won't be found as a uniprot entry because there are
        not uniprot entry name but gene names.  We therefore also scan
        missing entries by looking for gene names. Once found, the proposed
        items that contain the gene names and organism are candidates for the
        entry names. There may be several solutions though, which explain why
        the values in the :attr:`_mapping` dictionary  are made of lists. If several
        candidates are found, warning and raised.

        Results are stored in :attr:`_mapping` and in the dataframe itself.

        Let us show one example with 3 protein names that cover all cases:

            * DIG1, is a valid uniprot entry
            * ASC1 is not a uniprot entry. It is a gene name from which the entry
              may be retrieved automatically.
            * LEU1 is a gene name AND a uniprot entry. This is an ambiguous
              case. The default is to use the uniprot entry but if you call
              :meth:`check_entries_versus_sequence` (after meth:`set_annotations`)
              you will see that there is a mismatch meaning that LEU1_YEAST provided
              in the protein column is catually not the protein name but the gene name


        ::

            >>> import pandas as pd
            >>> from msdas import *
            >>> df = pd.DataFrame({'Protein':['DIG1', 'LEU1', 'ASC1'],
                'Sequence_Phospho':['S(Phospho)APAQVTQHSK', 'VEVTS(Phospho)EDEK',
                                    'DS(Phospho)VTIISAGNDK'],
                'Psite':['S142','S495', 'S166']})
            >>> a = Annotations(df, "YEAST)
            >>> a.get_uniprot_entries()
            >>> a._mapping
            {'ASC1_YEAST': ['P38011', 'P01120'],
             'DIG1_YEAST': ['Q03063'],
             'LEU1_YEAST': ['P06208-1', 'P06208']}

        Here, DIG1 has one unique entry. This is expected because DIG1 is in fact
        an entry name (unique by definition). ASC1 is a gene name. This method
        figures out that it correspond to either P38011 or P01120. There are
        several entries because mapping from gene to protein is not unique.
        By default, the entry with highest score appears first. There is no 100% guarantee
        that this mapping is correct and :meth:`check_entries_versus_sequence`
        should be called to check that the peptide sequence is contained in
        this entry sequence. The last case (LEU1) is even more problematic because
        it is a valid entry name even though the protein name provided is actually
        a gene name... again call :meth:`check_entries_versus_sequence`.


            >>> a.set_annotations()
            >>> a.check_entries_versus_sequence()
            P06208-1 not found in the annotations index

        So, here we are told that amongst the 3 entries, P06208-1 is not found.
        This the LEU1 case. If you were to use batch tool, you would figure out
        given the peptide sequence that this is actually LEUC_YEAST entry with
        uniprot entry LEUC_YEAST/P07264.

        So, you need to manually update the mapping:

            >>> a._mapping['LEU1_YEAST'] = ['P07264']
            >>> a.get_uniprot_entries() # to update the main df with new entries
            >>> a.set_annotations() # to retrieve the sequence of LEUC_YEAST
            >>> a.check_entries_versus_sequence()

        .. seealso:: :meth:`set_annotations`

        """
        # get the mapping using bioservices.uniprot

        # apply function is 3 times slower than list...
        # entry_names = self.df.Protein.apply(lambda x: x + "_" + self.organism)
        self._update_species_to_find()
        if len(self._species_to_find)>0:
            self.logging.info("Fetching uniprot accession numbers for %s entries" % len(self.df.Protein))
            self.logging.info("Fetching uniprot accession numbers for %s unique entries" % len(self.df.Protein.unique()))
            mapping = self._uniprot.multi_mapping(fr="ID", to="ACC",Nmax=Nmax,
                                              query=self._species_to_find)
            for k,v in mapping.iteritems():
                if k not in self._mapping.keys():
                    self._mapping[k] = v

        # some species may not be found (secondary accession number) if _human
        # appended in tcell case. so we may need to call again the mapping but
        # without the appended organism string.
        self._update_species_to_find()
        if len(self._species_to_find):
            self.logging.info("Some species were not found ({}). Using secondary accession:".format(len(self._species_to_find)))
            self.logging.info("Fetching uniprot without trailing species")
            self.logging.info("Fetching  %s new ones " % len(self._species_to_find))
            mapping = self._uniprot.multi_mapping(fr="ID", to="ACC",
                    query=[x.split("_")[0] for x in self._species_to_find], Nmax=Nmax)
            for k,v in mapping.iteritems():
                self._mapping[k+ "_" + self.organism] = v

        # Some are not yet found. this could be because the provided protein name is actually a
        # gene name...
        def func(x, tag):
            if len(x)==0:
                return False
            else:
                return tag in x[0].split()

        self._genes = {}
        self._update_species_to_find()
        if len(self._species_to_find):

            self.logging.info("Some species are still not found {}. Trying to use gene names".format(len(self._species_to_find)))
            self.logging.info("Fetching uniprot accession numbers for those without _species appended")
            self.logging.info("Fetching  %s new ones " % len(self._species_to_find))
            for i,this in enumerate(self._species_to_find):
                if " " in this:
                    continue
                self.logging.info("Searching for entry {}/{} for gene names".format(i+1,len(self._species_to_find)))
                df = self._uniprot.get_df(this.split("_")[0], organism=self.organism)
                l1 = df['Gene names'].apply(lambda x : func(x,this.split("_")[0] ))
                l2 = df['Entry name'].apply(lambda x: x.endswith(self.organism))

                if sum(l1&l2) >= 1:
                    k = list(df.ix[l2&l1]['Entry name'])
                    v = list(df.ix[l2&l1]['Entry'])
                    self.logging.debug(k, v, this)
                    if k in self._mapping.keys():
                        raise ValueError("!!!!!!!%s Already in the dictionary " % k)
                    #self._mapping[k] = v
                    self._mapping[this] = v
                else:
                    print("skipping %s... sum=%s" % (this, sum(l1&l2)))

        self._update_species_to_find()
        if len(self._species_to_find):
            self.logging.info("Some species were not found. Using gene names")
            self.logging.info("Fetching uniprot accession numbers for those without _species appended")
            self.logging.info("Fetching  %s new ones " % len(self._species_to_find))

        self._append_uniprot_entries_to_df()

    def _append_uniprot_entries_to_df(self):

        if "Entry" in self.df.columns:
            self.logging.warning("Overwritting column called Entry in the dataframe")

        # get list of unique entry names.
        entry_names = self.df.Protein.apply(lambda x: x + "_" + self.organism)
        #remapping = [(k,v[0]) for k,va.g in a._mapping.iteritems()]
        # add into dataframe the uniprot entries but must have same order as
        # in the dataframe (entry_names)
        uniprot_entries = []
        for name in entry_names:
            # if not found, let us use unknown as a label. could use NA?
            uniprot_entry = self._mapping.get(name, "")
            if uniprot_entry == "":
                uniprot_entries.append("")
                print("!! ", name, " not found")
            else:
                if len(uniprot_entry)>1:
                    self.logging.info("Found entry with several matches: %s %s . Only first one is selected (highest uniprot score)" % (name,
                                      uniprot_entry))
                uniprot_entries.append(uniprot_entry[0])

        # index=df.index is important to use the join afterwards
        thisdf = pd.DataFrame({'Entry': uniprot_entries}, index=self.df.index)
        if "Entry" in self.df.columns:
            del self.df['Entry']
        self.df = self.df.join(thisdf)

    def _append_uniprot_entry_names_to_df(self):

        if isinstance(self.annotations, types.NoneType) == True:
            self.error("must call set_annotations first")
            return
        # let us add the Entry_name column as well
        entry_names = [self.annotations.ix[e]['Entry name'] if e in self.annotations.index else "" for e in self.df.Entry]
        self.df['Entry_name'] = entry_names

    def plot_goid_histogram(self, drop_duplicates=True):
        """Histrogram of the number of GO terms per petide or protein

        :param drop_duplicates: ignore duplicates entries

        .. plot::
            :width: 80%
            :include-source:

            from msdas import *
            m = Annotations(get_yeast_small_data(), "YEAST", verbose=False)
            m.set_annotations()
            m.plot_goid_histogram()

        .. todo:: is this functional process or not

        """
        if self.annotations is False:
            raise AttributeError(self._error_messages['annotations'])
        if drop_duplicates:
            entries = self.df.Entry.drop_duplicates()
        counter = self.annotations.ix[entries]['Gene ontology IDs'].apply(lambda x: len(x))

        M = counter.max()
        # if we want the GO per peptides, then we need to look at the original
        # dataframe that contains several psites per peptide. UniProt_entry is
        # not a set so values from counter may be duplicated, which is what we
        # want for this first figure
        duplicated_counter = [counter[x] for x in self.df.Entry]
        pylab.figure(1)
        pylab.clf()
        pylab.hist(duplicated_counter, bins=[x+.5 for x in range(0,M+1)])
        pylab.title("Distribution of number of GO id terms per peptide")
        pylab.grid()

        # annotations contains the unique protein entry, so here we get the number of GO terms per protein
        counter = self.annotations['Gene ontology IDs'].apply(lambda x: len(x))
        M = counter.max()
        pylab.figure(2)
        pylab.clf()
        pylab.hist(counter, bins = [x+.5 for x in range(0,M+1)])
        pylab.title("Distribution of number of GO id terms per protein")
        pylab.grid()

    def set_annotations(self, Nmax=100):
        """Fetched all information from uniprot and set :attr:`annotations`
        as a pandas dataframe.


        Look into the dataframe Entry column and update the annotations dataframe
        to populate missing entries. The Entry column in the :attr:`df` should have been
        populated by :meth:`get_uniprot_entries` with valid entries from Uniprot.

        If you have thousand of entries, this is taking a few minutes. You can
        save the annotations and read them back using :meth:`msdas.MassSpecReader.read_annotations`
        and :meth:`to_pickle`.


        """
        self.logging.info("Fectching information from uniprot. Takes some time")

        #could split if too long
        entries = [this for this in list(set(self.df.Entry)) if this]

        # not need to search again if already present in the attribute
        if self.annotations is not None:
            entries = [x for x in entries if x not in list(self.annotations.index)]

        if len(entries)==0:
            self.warning("No new entries found. Your annotations dataframe is already up-to-date")
            self.annotations.drop_duplicates(subset="Entry name", inplace=True)
            self._append_uniprot_entry_names_to_df()
            return

        annotations = self._uniprot.get_df(entries, nChunk=Nmax)
        annotations = annotations[annotations.Entry.apply(lambda x: x in entries)]

        if len(annotations) == 0:
            raise ValueError("your list of protein is empty")
        self.logging.info("Fectching {}".format(len(annotations)))
        annotations.set_index(["Entry"], inplace=True)
        if self.annotations is None:
            self.annotations = annotations
        else:
            self.annotations = self.annotations.append(annotations)
        #self.annotations.set_index(["Entry"], inplace=True)
        self.logging.info("Annotations have been loaded. You can save the annotations" +
            " dataframe attribute using x.to_pickle('annotations.pkl') " +
            " Next time, you could just load if using \n\n" +
            "     >>> m = readers.MassSpecReader(filename, mode='yeast')\n" +
            "     >>>  m.read_annotations('annotations.pkl')")

        #indices are the uniprot entry. Some may be identical with slightly different columns
        # but the entry name should be unique. Here, we keep the first instance of each entry
        self.annotations.drop_duplicates(subset="Entry name", inplace=True)
        self._append_uniprot_entry_names_to_df()

    def to_pickle(self, tag=None, overwrite=False):
        """Save annotations dataframe as a pickle

        :param tag: a tag to append to the name of the annotations file.
        :param overwrite: overwrite file if it exists

        filename is going to be organism_annotations_tag.pkl
        """
        filename = self.organism + "_annotations"
        if tag != None and isinstance(tag, str):
            filename += "_" + tag
        filename += ".pkl"
        if overwrite == False:
            if os.path.exists(filename):
                raise IOError("file %s already exists" % filename)
        self.annotations.to_pickle(filename)

    def read_pickle(self, filename):
        """Read annotations in pickled format as saved by :meth:`to_pickle`

        :param str filename: filename to read
        """
        try:
            self.annotations = pd.read_pickle(filename)

            # update the mapping dictionary
            for k,v in self.annotations['Entry name'].iteritems():
                if k not in self._mapping.keys():
                    self._mapping[v] = [k]
        except:
            self.logging.error("Could not read your file. Expected a pkl \
            containing a dataframe with Entry name and index being uniprot \
            indices. ")

    def hist_most_relevant_goids(self, N=10, tight_layout=True, wrap_length=40,
                                 drop_duplicates=True, **kargs):
        """Plot histogram of the GO identifiers found in all proteins.

        :param int N: restrict histogram to terms that appear at least N times
        :param int wrap_length:  wrap text on the y-axis by wrap_length (defaults to 40)
        :param drop_duplicates: drop the duplicated entries
        :param kargs: pandas.plot arguments accepted.

        .. plot::
            :include-source:
            :width: 80%

            from msdas import *
            m = Annotations(get_yeast_small_data(), "YEAST", verbose=False)
            m.set_annotations()
            m.hist_most_relevant_goids(N=5)


        .. todo:: this is made on the annotations dataframe. Should be
            done based on the entry names in the dataframe

        """
        if self.annotations is False:
            raise AttributeError(self._error_messages['annotations'])
        kargs['legend'] = kargs.get("legend", False)


        if drop_duplicates:
            entries = self.df.Entry.drop_duplicates()

        goids = [y for x in self.annotations.ix[entries]['Gene ontology (GO)'] for y in x]
        uniq_goids = set(goids)
        names = [x for x in uniq_goids]

        # let us wrap the string by 40 character max to avoid long labels in the figure
        names = ["\n".join(textwrap.wrap(name, width=wrap_length)) for name in names]

        count = [goids.count(x) for x in uniq_goids]
        df = pd.DataFrame({'name':names, 'size':count}, index=range(0, len(uniq_goids)))
        if N:
            subdf = df[df['size']>N].set_index("name")
        subdf.sort("size").plot(kind="barh", **kargs)
        if tight_layout:
            pylab.tight_layout()

    def check_entries_versus_sequence(self):
        """Check that peptide sequence are contained in uniprot sequence

        This is a very good sanity check on the validity of the uniprot entry names
        found by :meth:`get_uniprot_entries` method

        If a peptide sequence is not found, it means that the protein name is
        not correct.

        See AnnotationsYeast class where the :meth:`AnnotationsYeast.update_mapping` is used to
        update the incorrect mapping.

        .. seealso:: :meth:`find_sequence_blast`

        """

        self.logging.info("Comparing peptide sequence in the attribute df with sequences in the annotations")
        self.logging.info("row index, protein name, uniprot entry")
        if isinstance(self.annotations, types.NoneType):
            raise Exception("annotations not set. call set_annotations")
        found = False
        for i in self.df.index:
            entry = self.df.ix[i].Entry
            if entry not in self.annotations.index:
                print("{} not found in the annotations index".format(entry))
                continue
            if self.df.ix[i].Sequence not in self.annotations.ix[entry].Sequence:
                if found == False:
                    print("Found unknown entries\nindex, protein name, uniprot entry ")
                    found = True
                print(i, self.df.ix[i].Protein, self.df.ix[i].Entry)

    def find_sequence_blast(self, seq, email):
        """Utility to search for a sequence using BLAST via bioservices

        :param str seq: the sequence
        :param email: a valid email address

        .. note:: This is using NCIBlast web service via
            `BioServices <https://pypi.python.org/pypi/bioservices>`_.
        """
        from bioservices import NCBIblast
        s = NCBIblast(verbose=self.level)
        jobid = s.run(program="blastp", sequence=seq, stype="protein",
                      database="uniprotkb", email=email)
        return s.getResult(jobid, "out")

    def to_csv(self, filename):
        """Export the dataframe with data and annotations into a CSV file

        :meth:`set_annotations` and :meth:`get_uniprot_entries` must have been called.


        """
        if "Entry" not in self.df.columns or "Entry_name" not in self.df.columns:
            raise ValueError("Entry or Entry_name missing in dataframe. You must call get_entries_uniprot and set_annotations methods")
        self.df.Identifier = self.df.Protein + "_" + self.df.Psite
        self.df.to_csv(filename, index=False, sep=",")




class AnnotationsYeast(Annotations):
    """Class dedicated to the YEAST data analysis

    This class is almost identical to :class:`Annotations`. It contains
    extra code to cleanup the mapping based on further manual investigations
    of the gene mapping to protein, which may be ambiguous (see
    :meth:`Annotations.get_uniprot_entries` for details).

    ::

        from msdas import *
        from easydev import gsf
        filename = gsf("msdas", "data", "YEAST_raw_sample.csv")
        r = MassSpecReader(r)
        a = AnnotationsYEAST(r)
        a.get_uniprot()
        a.set_annotations()


    Only 80% of the protein names are found directly using UniProt.
    The 20% remaining are actually gene names on which a mapping to protein has to be done.
    Yet, sometimes, there is an ambiguity that remains either because the gene name
    is also a valid entry name or because the gene maps to several entry names.
    This list gives some of these ambiguities. The first one is used by default
    (based on highest score) but may not be correct.
    See for instance :meth:`Annotations.check_entries_versus_sequence` to help
    you figuring out which one is the correct one.

        * ALD3_YEAST ['P54114', 'P40047']
        * ALD4_YEAST ['P46367', 'P54114']
        * CKI1_YEAST ['P20485', 'P23292']
        * CPR1_YEAST ['P14832', 'P16603']
        * PRS5_YEAST ['Q12265', 'P23638']
        * RPL16B_YEAST ['P26785', 'Q3E757', 'P05739']
        * RPL32_YEAST ['P38061', 'P14120']
        * RPL6A_YEAST ['Q02326', 'P05737']
        * RPS7A_YEAST ['P26786', 'P0CX36']
        * RPS7B_YEAST ['P48164', 'P0CX35']
        * CPR1_YEAST ['P14832', 'P16603']
        * ECM17_YEAST ['P40458', 'P47169']
        * RPL16B_YEAST ['P26785', 'Q3E757', 'P05739']
        * RPS7B_YEAST ['P48164', 'P0CX35']
        * ASC1_YEAST ['P38011', 'P01120']
        * ECM17_YEAST ['P40458', 'P47169']


     :Notes on the data:

        * NPL3_356^S349 has a wrong psite name. Given the sequence, it should be NPL3_S349
        * Same issue with TIF3 (IF4B_YEAST) where trailing number without phospho was removed by hand.
        * FLO9 1004^554^464^374^T329^T779, T(Phospho)GTFTSTSTEM(Oxidation)TTVTGTNGQPTDETVIVI should be T779
        * One entry is AD5;7 which is wrongly named to not clash with CSV format. The proper
          name is indeed AD5,7. We renamed it in the file as PUR2_YEAST. We checked the sequence
        * To find the mapping, we used blast from bioservices to figure out the sequence of
          the protein and checked on uniprot. See update_mapping function
        * Typo in the original code for ABP1 peptide: small e was found.
          KEPVKT eSP APAAK should be KEPVKT PSP APAAK
        * Possible typo is STE11_S323 location  should be S326
        * Also possible typoe IMP2 has 2 rows called IMP2' (note the quote)


    Here are proteins names provided in the Yeast_data_set that are actually gene
    names. Using bioservices, we figure out possible uniprot entries but similarlly
    to the YEF3 there is maybe an ambiguity on the name::

        77 P54114 ALD4   could be ALDH4_YEAST ([') whereas P54114 is ALD3
        2077 P23638 PRS5    KPR5_YEAST  '
        6664 P01120 ASC1  GBLP_YEAST  P38011

    there are 6 peptide labelled CTR9. In fact, it is 2 different peptides. first 4 are CTR9
    2 last are EAMAISEHNVKDDSDLSDKDNEYDEEQPR . This is defintely CTR9 but is not exactly in
    the sequence. Missing K at the end

    SIR1_YEAST exists but the peptide sequence  cannot be found either using blast or manual searh on uniprot
    possibly sirp1 but not yeast organism; it is yeasb. Actually may still be SIR1
    but Uniprot changed the sequence. See persona communicatoin with uniprot.


    """
    def __init__(self, data,  verbose=True, annotations=None, **kargs):
        """.. rubric:: Constructor

        Same as :class:`Annotations` except that organism must not be provided.
        """


        super(AnnotationsYeast, self).__init__(data=data, organism="YEAST",
            verbose=verbose, annotations=annotations, **kargs)


        index = list(self.df[self.df.Protein == "ADE5;7"].index)
        if len(index)==1:
            self.warning("renaming ADE5;7 (should be ADE5,7 anyway) into PUR2")
            self.df  = self.df.set_value(index[0], "Protein", "PUR2")
            self._rebuild_identifier()
        if len(index)>1:
            raise NotImplementedError

    def update_mapping(self):
        """Update the mapping with known keys

        There are issues in the naming because of a mixing of protein and gene
        names. Methods in :class:`Annotations` found most of the mapping.
        However, there are some ambiguities and the mappind dictionary is corrected
        as follows. Checked with uniprot and blast.

        we remove "ADE5,7_YEAST" (replaced by PUR2 if found)

        ================ ======== ============================================
        ================ ======== ============================================
        ADE5,7_YEAST     P07244   PUR2_YEAST
        YEF3_YEAST       P16521
        FEN1_YEAST       P25358   FEN1 is actually ELO2_YEAST
        HIS7_YEAST       P33734   HIS7 is actually HIS5_YEAST
        LEU1_YEAST       P07264   this is LEUC_YEAST P07264
        NTH1_YEAST       P32356   TREA_YEAST
        NTH2_YEAST       P35172   TREB_YEAST
        ECM17_YEAST      P47169   MET5_YEAST
        RPL6A_YEAST      Q02326   RL6A_YEAST
        RPS7A_YEAST      P26786   RS7A_YEAST
        YJU2_YEAST       P28320   CWC16_YEAST no blast results but checked
        ASC1_YEAST       P38011   GBLP_YEAST
        PSA1_YEAST       P41940   no blast result; gene name MPG1_YEAST
        CTR9_YEAST       P89105   CTR9_YEAST correct  but see note below
        IMP2_YEAST       P32351   IMPX_YEAST no blast result;IMP2 is gene name
        ================ ======== ============================================


        .. note:: there is also a protein called IMP2' (note the quote), which
            presumably is also IMP2. Kept as it is fow now


         """
        if len(self._mapping) == 0:
            raise ValueError("You should call set_annotations and get_uniprot first")
        #self._mapping['ADE5,7_YEAST'] = ['P07244'] # PUR2_YEAST
        self._mapping['YEF3_YEAST'] = ["P16521"]  #
        self._mapping['FEN1_YEAST'] = ["P25358"]  # FEN1 is actually ELO2_YEAST
        self._mapping['HIS7_YEAST'] = ["P33734"]   # HIS5 is actually HIS5_YEAST
        self._mapping['LEU1_YEAST'] = ["P07264"]  # this is LEUC_YEAST P07264
        self._mapping['NTH1_YEAST'] = ['P32356']  # TREA_YEAST
        self._mapping['NTH2_YEAST'] = ['P35172']  # TREB_YEAST
        self._mapping['ECM17_YEAST'] = ['P47169'] # MET5_YEAST
        self._mapping['RPL6A_YEAST'] = ['Q02326'] # RL6A_YEAST
        self._mapping['RPS7A_YEAST'] = ['P26786'] # RS7A_YEAST
        self._mapping['YJU2_YEAST'] = ['P28320']  # CWC16_YEAST no blast results but checked
        self._mapping['ASC1_YEAST'] =  ['P38011'] # GBLP_YEAST
        self._mapping['PSA1_YEAST'] = ['P41940']  # no blast result; a gene name corresponding to MPG1_YEAST
        self._mapping['CTR9_YEAST'] = ['P89105']  # CTR9_YEAST  P89105 this is correct  but see note documetnation
        self._mapping['IMP2_YEAST'] = ['P32351'] # IMPX_YEAST no blast result; IMP" is gene name

