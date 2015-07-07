# -*- coding: utf-8 -*-
# -*- python -*-
#
#  This file is part of MS-DAS software
#
#  Copyright (c) 2011-2012 - EBI-EMBL
#
#  File author(s): Claudia Hernand, Marti Bernardo,
#       Thomas Cokelaer <cokelaer@ebi.ac.uk>
##
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website:
#
##############################################################################
"""
Created on Fri Feb  7 16:58:10 2014

@author: chernand
"""
import os
import pandas as pd
from tools import Requires


__all__ = ["PhosphoGRID"]


class PhosphoGRID(Requires):
    """Build PKN from phosphogrid

    You need first a list of protein names that you are interested in.
    This list can be obtained from the merger instance as follows::

        from msdas import *
        m = MassSpecReader(get_yeast_small_data, "yeast")
        gene_names = set(list(m.df.Protein))

    Then, you can use this class as follows to visualise and save into a SIF file the
    network built from phosphogrid database.

    .. plot::
        :include-source:
        :width: 70%

        >>> from msdas import *
        >>> r = MassSpecReader(get_yeast_small_data())
        >>> gene_names = set(list(r.df.Protein))
        >>> p = phosphogrid.PhosphoGRID()
        >>> p.run(gene_names=gene_names)
        >>> p.export2sif()
        >>> p.plot()

    2 databases are used. Their filenames are hard coded as:

        #. BIOGRID-PTM-RELATIONSHIPS-3.1.93.ptmrel.txt
        #. BIOGRID-PTM-15-3.1.93.ptmtab.txt

    You can change the name by changing :attr:`dbR_filename` and
    :attr:`dbP_filename`.

    """

    def __init__(self, directory="../../share/data/"):
        """.. rubric:: constructor

        :param str directory: directory where to find the databases



        They can be found at Claudia's Hernand home directory.

        """
        super(PhosphoGRID, self).__init__()
        self.directory = directory

        self.dbR_filename = "BIOGRID-PTM-RELATIONSHIPS-3.1.93.ptmrel.txt"
        # original file was huged. We removed sequence, author,pumed that are not used.
        #self.dfP_filename = "BIOGRID-PTM-15-3.1.93.ptmtab.txt"
        self.dfP_filename = "BIOGRID-PTM-15-3.1.93_SMALL.ptmtab.txt"

    def run(self, gene_names=None):
        """Build dataframe from the relations found in the PhosphoGRID databases


        relations are saved in :attr:`dfSIF`

        """
        #===PhosphoGRID Database
        #Open file with relationships PhosphoGRID tab delimited files
        import gzip
        dbr_arch = gzip.GzipFile(self.directory+os.sep+self.dbR_filename+'.gz')
        dbp_arch = gzip.GzipFile(self.directory+os.sep+self.dfP_filename+'.gz')
        import StringIO

        dfR = pd.read_csv(StringIO.StringIO(dbr_arch.read(self.dbR_filename)),
                          sep='\t', header=0)
        dfP = pd.read_csv(StringIO.StringIO(dbp_arch.read(self.dfP_filename)),
                          sep='\t', header=0)

        #===Merge PhosphoGRID data frame using PTMID as matching column
        dfRP = dfR.merge(dfP, on="PTMID")
        dfRP.columns = [this.replace("_x","_reg") for this in dfRP.columns.map(str)]
        dfRP.columns = [this.replace("_y","_tgt") for this in dfRP.columns.map(str)]
        dfRP["RelationshipCode"] = ""

        #===Asign interaction type according to  enzymatic activity
        dfRP.RelationshipCode[dfRP.Relationship == "kinase"] = 1
        dfRP.RelationshipCode[dfRP.Relationship == "phosphatase"] = -1
        dfRP["Site"] = dfRP.Residue+dfRP.Position.map(str)

        #Yeast dataframe
        #Select unique gene names
        if gene_names == None:
            filename = self.directory + os.sep + "RawData.csv"
            df = pd.read_csv(filename)
            gene_names = df["Protein"].drop_duplicates()
        else:
            pass

        #===Merge Yeast experimental dataset with PhosphoGRID
        # Intersection between both datasets querying  by Standard 
        # name == OfficialSymbol both in regulators and substrates
        columns=["OfficialSymbol_reg","RelationshipCode","OfficialSymbol_tgt",
                 "Relationship","Identity","Site"]
        out = pd.DataFrame(columns=columns)
        for gene in gene_names:
            this = dfRP[dfRP.OfficialSymbol_reg == gene][columns].drop_duplicates()
            out = out.merge(this, "outer")
            this = dfRP[dfRP.OfficialSymbol_tgt == gene][columns].drop_duplicates()
            out = out.merge(this, "outer")

        #Create SIF file
        self.dfSIF = out[["OfficialSymbol_reg","RelationshipCode",
                        "OfficialSymbol_tgt"]].drop_duplicates()

    def export2sif(self, filename="matchYeastPhosGRID_protein.sif"):
        """Export the found relations into a SIF file

        :param str filename: Defaults to matchYeastPhosGRID_protein.sif

        """
        # could be a decorator but argument are hidden, so just call it for now.
        self._requires("dfSIF", "Call run() method first")
        # TODO : check that dfSIF is present otherwise call run()
        self.dfSIF.to_csv(filename, index=False, header=False, sep="\t")
        self.sif_filename = filename

    def plot(self):
        """Plot the relations using CNOGraph

        .. warning:: consider using a faster visualisation tools if the number
            of edges/nodes exceeds ~ 300

        .. todo:: build the SIF without saving it and import directly in the
            cnograph structure.
        """
        self._requires("sif_filename", "sif_filename not set yet. Please, Call export2sif() method first.")
        from cno import CNOGraph
        c = CNOGraph(self.sif_filename)
        c.plot()










