Temporary directory to provide PKN and MIDAS files related to the YEAST
analysis (and fitness plots).

--------------

Please download the supplementary data archive in `tar gz file
format <suppdata.tar.gz>`_ or `zip format <suppdata.zip>`_, which
contains 8 directories. Each directory contains

-  a logic model in SIF and SBML-qual formats called PKN-YEAST.sig and
   PKN-YEAST.xml respectively
-  a MIDAS file containing the data called MD-YEASY.csv
-  a R script that loads the PKN, data and use a set of ODE parameters
   to simulate the model
-  a PDF called Fitness.pdf, which is the result of the simulation
   fitted to the data

--------------

Instructions under Linux:

::

    - tar xvfz suppdata.tar.gz
    - Browse the directory
    - In a given directory (e.g., STE20_All), type "sh create_fitness_plot.sh", which create the PDF Fitness.pdf (already provided but can be reproduced). data.R contains the decision vector for the median model.

You can also view the fitness plot individually here below:

-  STE20\_All: `fitness
   plot <http://www.cellnopt.org/data/yeast/STE20_All/Fitness.pdf>`_
-  STE20\_All\_HOG1\_GPD1: `fitness
   plot <http://www.cellnopt.org/data/yeast/STE20_All_HOG1_GPD1/Fitness.pdf>`_
-  STE20\_All\_pruned\_edges: `fitness
   plot <http://www.cellnopt.org/data/yeast/STE20_All_pruned_edges/Fitness.pdf>`_
-  STE20\_All\_HOG1\_GPD1\_pruned\_edges: `fitness
   plot <http://www.cellnopt.org/data/yeast/STE20_All_HOG1_GPD1_pruned_edges/Fitness.pdf>`_
-  STE20\_T511: `fitness
   plot <http://www.cellnopt.org/data/yeast/STE20_T511/Fitness.pdf>`_
-  STE20\_T511\_HOG1\_GPD1: `fitness
   plot <http://www.cellnopt.org/data/yeast/STE20_T511_HOG1_GPD1/Fitness.pdf>`_
-  STE20\_T511\_pruned\_edges: `fitness
   plot <http://www.cellnopt.org/data/yeast/STE20_T511_pruned_edges/Fitness.pdf>`_
-  STE20\_T511\_HOG1\_GPD1\_pruned\_edges: `fitness
   plot <http://www.cellnopt.org/data/yeast/STE20_T511_HOG1_GPD1_pruned_edges/Fitness.pdf>`_

You will need CNORode version 1.5.2 and CellNOptR version 1.11.2
available on
`www.cellnopt.org <http://www.cellnopt.org/downloads.html>`_

--------------

author: TC 2014
