Getting Started
================

First, clone the repository

.. code-block:: shell

   git clone https://gitlab.com/LMSAL_HUB/aia_hub/aiapy.git

Next, install the needed dependencies using Anaconda and activate 
the environment,

.. code-block:: shell

   cd aiapy
   conda env create -f environment.yml
   source activate aiapy-dev

Alternatively, you can download all of the packages listed in
`environment.yml` using `pip` if you're not using the Anaconda
installation. Finally, install the aiapy package,

.. code-block:: shell

   python setup.py install

Once development is a bit more stable, releases will be made on PyPI and conda-forge