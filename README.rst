aiapy
=====

aiapy is a Python package for analyzing data from the Atmospheric Imaging Assembly (AIA) instrument onboard NASA's Solar Dynamics Observatory spacecraft.
For more information, see the `aiapy documentation <https://aiapy.readthedocs.io/en/latest/>`__.
For some examples of using aiapy, see `our gallery <https://aiapy.readthedocs.io/en/latest/generated/gallery/index.html>`__.

Installation
------------

The current stable version of aiapy is available through the Python Package Index and can be installed via ``pip``

.. code-block:: shell

   pip install aiapy

or through the Anaconda distribution via ``conda-forge``,

.. code-block:: shell

   conda install -c conda-forge aiapy

These are the recommended ways to obtain and install aiapy.
Alternatively, you can install the current development version directly from GitLab,

.. code-block:: shell

   git clone https://gitlab.com/LMSAL_HUB/aia_hub/aiapy.git
   cd aiapy
   pip install -e .[all]

If you'll be developing aiapy, see the `development setup guide <https://aiapy.readthedocs.io/en/latest/develop.html>`__.

Testing
-------

If you want to run the test suite, first install the dev requirements,

.. code-block:: shell

   pip install -e .[dev]

and then run

.. code-block:: shell

   pytest --remote-data=any

If an internet connection is not available, exclude the ``--remote-data`` flag.

A valid install of IDL and SSW are required to run the tests that compare results from aiapy and SSW.
If one is not available, these tests are automatically skipped.

The entire test suite can also be run using tox.
For additional instructions, please see the `SunPy development guide on testing <https://docs.sunpy.org/en/latest/dev_guide/tests.html>`__.

Citing
------

If you use aiapy in your scientific work, we would appreciate you citing it in your publications.
The latest citation information can be found in the `documentation <https://aiapy.readthedocs.io/en/latest/about.html>`__, or obtained with ``aiapy.__citation__``.

Contributing
------------

We love contributions! aiapy is open source, built on open source, and we'd love to have you hang out in our community.

If you can write code at all, you can contribute code to open source.
Contributing to open source projects is a fantastic way to advance one's coding skills ; it's trying to create something, making mistakes, and learning from those mistakes.
That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either.
You can help out by writing documentation, tests, or even giving feedback about the project (and yes - that includes giving feedback about the contribution process).
Some of these contributions may be the most valuable to the project as a whole, because you're coming to the project with fresh eyes, so you can see the errors and assumptions that seasoned contributors have glossed over.
