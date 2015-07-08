

#################################
MS-DAS documentation
#################################

Motivation 
###########

**MS-DAS** stands for **M**\ass **S**\pectrometry **D**\ata **A**\nalysi\ **S**\. MS-DAS is a Python package that provides tools to read, manipulate and visualise Mass Spectrometry measurements. MS-DAS is not about pre-processing (at the spectra level) but rather focuses on providing tools to get more insights about relation between the peptides (e.g., for logical modeling).

MSDAS was designed to analysis a specific data set (YEAST) with a specific format. In this package, we try to enforce a common data structure so as to write all the pertinent tools only once. 

MSDAS allows one to 

#. read some Mass Spectrometry data stored in excel files
#. Align several files in a single data structure (dataframe).
#. Add annotations using UniProt service (sequence, protein, go terms)
#. Cluster the data if required
#. Creation of prior knowledge networks using phosphogrid and uniprot services
#. Export of the data into MIDAS format
#. Understanding the Quality of the data (if replicates are available)
#. Descriptive analysis is also possible at different level.

Installation
###################

From the source file, type::

    python setup.py install

If installing globally, you may need root permission (e.g., sudo access)

MS-DAS relies on standard libraries (e.g., numpy, matplotlib, pandas, 
networkx) and bioservices, which is also available on pypi website.

If the command above does not succeed, it may be related to one of the
dependencies that failed to be installed, in which case you will need to install
them manually. We would recommend to install NumPy first. Others should then work out of the box.
You can install all those packages using **pip**, for instance::

    pip install numpy


User guide
##################


.. toctree::
    :maxdepth: 2
    :numbered:

    quickstart.rst
    pkn.rst
    notebooks.rst

References
##################


.. toctree::
    :maxdepth: 2
    :numbered:

    references




.. .. toctree::

..    biblio.rst


..    ChangeLog.rst
    authors.rst
    contribute

Others
##########

.. toctree::

    report/index.rst

