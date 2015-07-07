from msdas import *

verbose = False




def test_clustering_yeast():
    yeast = MassSpecReader(get_yeast_small_data(), mode="yeast")
    c = clustering.MSClustering(yeast)
    af = c.run_affinity_propagation_clustering(preference=-90)
    af.plotClusters_and_protein("DIG1")

    yeast = MassSpecReader(get_yeast_filenames()[0], mode="yeast")
    c = clustering.MSClustering(yeast)
    af = c.run_affinity_propagation_clustering(preference=-90)
    af.plotClusters_and_protein("DIG1")


    c.dendogram(method="euclidean")
    c.dendogram(method="dot")
    c.dendogram(pdist="euclidean")

    c.plotClusteredProtein("DIG1", scale=True)

def test_clustering_affinity():
    yeast = MassSpecReader(get_yeast_small_data(), mode="yeast" ,
            verbose=verbose)
    c = clustering.MSClustering(yeast)
    data = c.get_group()
    af = Affinity(data, preference=-90, method="euclidean")
    print(af)
    #af.preferences
    af.plotClusters()
    c.plot_vectors("DIG1")
    c.plot_vectors("DIG1", normalise=False,  scale=True, legend=False)





def test_affinity():

    yeast = MassSpecReader(get_yeast_small_data(), mode="yeast" , verbose=verbose)
    c = clustering.MSClustering(yeast)
    a = Affinity(c.get_group(), "euclidean")
    a.preference == -30

    try:
        a.plotClusters_and_protein
        assert False
    except:
        assert True

