# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 16:02:50 2014

@author: cokelaer
"""
from readers import MassSpecReader
from cno import CNOGraph
from cno.io.multigraph import  CNOGraphMultiEdges
import pylab
import pandas as pd

__all__ = ["NetworkFromUniProt", "CombineNetworks"]


class NetworkFromUniProt(object):
    """Build a PKN based on the uniprot **interact with** **annotations**

    The annotations dataframe can be obtained from the MassSpecMerger instance

    ::

        from msdas import *
        a = annotations.Annotations(get_yeast_small_data(), "YEAST")
        a.get_uniprot_entries()
        a.set_annotations()

    Then, you can create this object::

        n = network.NetworkFromUniProt(a.annotations)

    And finally get a graph structure that extract all relations found in
    the annotations dataframe based on the  uniprot field called "Interacts with".

    """
    def __init__(self, annotations, verbose=True):
        """.. rubric:: constructor

        :param df annotations: a dataframe similar to the :attr:`annotations`
            found in :meth:`MassSpecReader`
        """
        if isinstance(annotations, str):
            self.annotations = pd.read_pickle(annotations)
        elif isinstance(annotations, MassSpecReader):
            self.annotations = annotations.df
        else: # trying as a dataframe
            self.annotations = annotations.copy()
        self.verbose=verbose

    def get_cnograph_intact(self, label="entry_name"):
        """Return cnograph made of the protein names found in the interactions
        of the annotations.


        .. plot::
            :include-source:
            :width: 50%

            from msdas import *
            a = annotations.Annotations(get_yeast_small_data(), "YEAST")
            a.get_uniprot_entries()
            a.set_annotations()
            n = network.NetworkFromUniProt(a.annotations)
            c = n.get_cnograph_intact()
            c.plotdot()


        """
        assert label in ["entry_id", "entry_name"]
        c = CNOGraph()
        interactions = self.annotations["Interacts with"]
        
        # add all nodes
        c.add_nodes_from(interactions.index)

        # some have no interactions in which case, it is filled with NaN. let us drop those
        # entries. 
        interactions = interactions.dropna()
        indices = interactions.index
        for i, index in enumerate(indices):
            print("{}/{}".format(i+1, len(indices)))
            these_interactions = interactions.ix[index].split(';')
            these_interactions = [x.strip() for x in these_interactions]
            for interaction in these_interactions:
                if interaction == "Itself":
                    interaction = index
                c.add_reaction("{}={}".format(index, interaction))

        if label == "entry_id":
            c._signals = list(self.annotations.index)
        else:
            # bioservices required because interacting species may not be part
            # of the list of measurements,
            from bioservices import UniProt
            u = UniProt(verbose=self.verbose)
            mapping = u.multi_mapping(fr="ACC", to="ID", query=c.nodes())
            for k, v in mapping.iteritems():
                if len(mapping[k])>1:
                    print("ambigous case {} with more than 1 mapping. will take only first".format(k))
                mapping[str(k)] = str(v[0].split("_")[0])
            c.relabel_nodes(mapping)

            measured = [x.split("_")[0] for x in self.annotations['Entry name']]
            c._signals = measured

        return c





class CombineNetworks(object):
    """Combine several PKN from different methods


    THis class serves as an example on how to combine several PKNs into
    a common one. The input PKN used may come from:

    #. In this example, you need to build a uniprot PKN using
    :class:`NetworkFromUniProt`, a PKN using PhosphoGrid
    :class:`msdas.phospho.PhosphoGrid` and a list of names to indicates nodes
    where you have measurements.

    .. plot::
        :include-source:
        :width: 50%

        # Get list of names
        from msdas import *
        a = annotations.Annotations(get_yeast_small_data(), "YEAST")
        a.get_uniprot_entries()
        a.set_annotations()
        n = network.NetworkFromUniProt(a.annotations)

        names = list(set(a.df.Protein))

        from easydev import get_share_file as gsf
        n = network.CombineNetworks(
            {"Curated": gsf("msdas", "data", "PKN-yeast.sif"),
            "UniProt": "PKN-uniprot.sif",
            "PhosPho": "PKN-phospho.sif"},
            signals=names[:], stimuli=["a", "NaCl"])

        c = n.get_digraph()

        c.plot()
        #c.export2sif("PKN-combined.sif")

    """


    def __init__(self, dict_network, stimuli=[], signals=[], cutnonc=True,
                 remove_self_loops=True):
        """

        :param dict dict_network: a dictionary of network. keys are used for labelling
            values must be a SIF filename
        :param list stimuli: list of stimuli
        :param list signals: list of signals
        :param bool cutnonc: remove useless nodes, not measured or without influence on signals
        :param bool remove_self_loops: remove self loop from the network.

        """
        self.filenames = []
        self.labels = []
        for k,v in dict_network.iteritems():
            self.filenames.append(v)
            self.labels.append(k)
        self.stimuli = stimuli[:]
        self.signals = signals[:]
        self.cutnonc = cutnonc
        self.remove_self_loops = remove_self_loops

    def plot_multiedge_graph(self, cmap="jet"):
        """Creates a multiedge graph and plots it

        :param cmap: a valid color map from matplotlib. jet, spring, hot, ...
        :return: CNOGraphMultiEdges object


        .. plot::
            :include-source:
            :width: 50%

            # Get list of names
            from msdas import *
            from easydev import gsf

            m = MassSpecReader()
            m.read_annotations(gsf("msdas", "data", "YEAST_annotations_small.pkl"))
            n = network.NetworkFromUniProt(a.annotations)

            names = list(set(m.df.Protein))


            n = network.CombineNetworks(
                {"Curated": gsf("msdas", "data", "PKN-yeastScaffold.sif"),
                 "UniProt": "PKN-uniprot.sif",
                 "PhosPho": "PKN-phospho.sif"},
                 signals=names[:], stimuli=["a", "NaCl"])

            c = n.plot_multiedge_graph()
            c.plot()


        """
        N = len(self.labels)
        values = pylab.linspace(.1,.9, N)

        # build network
        c = self.get_multiedge_graph()
        c.plot(edge_attribute="edgecolor", edge_attribute_labels=False, cmap=cmap)

        # #build legend
        for i, label in enumerate(self.labels):
            print label, c._get_hex_color_from_value(values[i], cmap)
            pylab.barh(0,0,1,color=c._get_hex_color_from_value(values[i], cmap), label=label)
        pylab.legend(title="edge legend", fontsize="small", loc="lower right")

        return c

    def get_multiedge_graph(self):
        """Creates a multiedge graph from the input networks

        :return: CNOGraphMultiEdges object


        """
        # build network
        N = len(self.labels)
        values = pylab.linspace(.1,.9, N)

        c = CNOGraphMultiEdges()
        for i,filename in enumerate(self.filenames):
            print("Reading {}".format(filename))
            graph = CNOGraph(filename)
            for e in graph.edges(data=True):
                c.add_edge(e[0], e[1], source=self.labels[i], edgecolor=values[i], **e[2])

        c._signals = self.signals[:]
        c._stimuli = self.stimuli[:]

        if self.signals and self.stimuli and self.cutnonc:
            c.cutnonc()
        if self.remove_self_loops:
            c.remove_self_loops()

        return c


    def get_digraph(self, sources_priority=["Curated", "PhosPho", "UniProt"]):
        """Creates a directed graph from the input networks"""
        multic = self.get_multiedge_graph()

        # We could cast the multiedge into a normal digraph but when edges have different
        # meaning, (e.g., activation/inhibition), then there is an ambiguity and
        # the remaining edge may not be the one we want. There is no clear answer
        # on what is the best one to keep but we could prioritise the
        # edge based on the source (e.g., curated is most trustful that a automatic
        # database).


        # let us merge edges with same input/output nodes if they have different
        # link base on the priority of the source.
        if sources_priority==None:
            labels = self.labels[:]
        else:
            labels = sources_priority[:]

        for node1 in multic.nodes():
            for node2 in multic.edge[node1]:
                # figure out the one to keep and remove others
                keys = multic.edge[node1][node2].keys()
                sources = [multic[node1][node2][this]['source'] for this in keys]

                # FIXME: could have a + and - from same source ?
                #links = [multic[node1][node2][this]['link'] for this in multic.edge[node1][node2].keys()]

                orders = sorted(zip(sources, keys),
                               cmp=lambda x,y: cmp(labels.index(x[0]), labels.index(y[0])))
                # keep first one, get rid of others

                if len(orders)>1:
                    if orders[0][0] == orders[1][0]:
                        print("looks like {}-{ edge is ambiguous}".format(node1, node2) )
                    # remove all except first one
                    for order in orders[1:]                         :
                        key = order[1]  # order is a 2-length tuple n e.g., ("Curated", 0)
                        print("removing {}-{}, key={} (source={})".format(node1, node2, key, order[0]))
                        multic.remove_edge(node1, node2,key)
                    print("Keeping {}-{}, source={}, key={}\n".format(node1, node2, orders[0][0], orders[0][1]))
                multic.edge[node1][node2][orders[0][1]]['source'] = sources

        c = CNOGraph(multic)
        c._signals = self.signals[:]
        c._stimuli = self.stimuli[:]
        return c

















