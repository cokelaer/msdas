

#################################
MS-DAS documentation
#################################

Motivation 
###########

**MS-DAS** stands for **M**\ass **S**\pectrometry **D**\ata **A**\nalysi\ **S**\. MS-DAS is a Python package that provides tools to read, manipulate and visualise Mass Spectrometry measurements. MS-DAS is not about pre-processing (at the spectra level) but rather focuses on providing tools to get more insights about relation between the peptides (e.g., for logical modeling).

MSDAS was designed based on two data sets with different formats and conventions. Usually, one converts its input data set (usually excel documents) into its own data structure, and then writes tools (e.g., merging of different data sets, clustering of time courses, averaging, normalisation, plotting...) to extract relevant information. In this package, we try to enforce a common data structure so as to write all the pertinent tools only once. This format may change so as to use more standard format (e.g. mzML); see references and quickstart section for the format definition. For mzML format, one can use existing libraries such as **pymzml** or **pyteomics**.

For now, MSDAS allows one to 

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

MS-DAS relies on standard libraries (e.g., numpy, matplotlib, pandas, scikit-learn,
networkx) and bioservices, which is also available on pypi website.

If the command above does not succeed, it may be related to one of the
dependencies that failed to be installed, in which case you will need to install
it manually. The trickiest to install is numpy that may require (in some rare cases) a manual installation; 
we would recommend to install numpy first. Others should then work out of the box.
You can install all those packages using **pip**, for instance::

    pip install pandas


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

