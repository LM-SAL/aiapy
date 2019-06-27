aiapy
-------------------------------

.. image:: http://img.shields.io/badge/powered%20by-SunPy-orange.svg?style=flat 
    :target: http://www.sunpy.org                                               
    :alt: Powered by SunPy Badge
.. image:: https://readthedocs.org/projects/aiapy/badge/?version=latest
    :target: https://aiapy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

aiapy is a Python package for analyzing data from the Atmospheric Imaging Assembly instrument
onboard the Solar Dynamics Observatory spacecraft.

Installation
------------

First, clone the repository

.. code-block:: shell

   git clone https://gitlab.com/LMSAL_HUB/aia_hub/aiapy.git

Next, install the needed dependencies using Anaconda and activate 
the environment,

.. code-block:: shell

   cd aiapy
   conda env create -f environment.yml
   source activate aiapy-dev

Finally, install the aiapy package,

.. code-block:: shell

   python setup.py install

Once development is a bit more stable, releases will be made on PyPI and conda-forge

License
-------

This project is Copyright (c) AIA Instrument Team and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.

Contributing
------------

We love contributions! aiapy is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

Note: This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
aiapy based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.