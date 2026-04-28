.. _aiapy-installing:

******************
aiapy Installation
******************

The easiest way to install aiapy is to follow the instructions on
:ref:`sunpy-tutorial-installing`. The only difference is that you will
need to install aiapy instead of sunpy.

If you are not using Miniforge:

.. code-block:: bash

    $ conda config --add channels conda-forge
    $ conda config --set channel_priority strict

And then you can do:

.. code-block:: bash

    $ conda install aiapy

Alternatively, aiapy can be installed with pip:

.. code-block:: bash

    $ pip install aiapy