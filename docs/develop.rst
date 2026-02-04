.. _aiapy_dev-guide:

============
Contributing
============

Contributing to open source projects is a good way to advance one's coding skills; it's trying to create something, making mistakes, and learning from those mistakes.
That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either.
You can help out by writing documentation, tests, or even giving feedback about the project (and yes - that includes giving feedback about the contribution process).
Some of these contributions may be the most valuable to the project as a whole, because you're coming to the project with fresh eyes, you can see the errors and assumptions that seasoned contributors have glossed over.

Issue Tracking
--------------

All bugs, feature requests, and other issues related to ``aiapy`` should be recorded using the GitHub issue tracker.
You can find instructions for how to create an issue `here <https://github.com/LM-SAL/aiapy/issues>`__.

All conversation regarding these issues should take place on the issue tracker.
When a merge request resolves an issue, the issue will be closed and the appropriate merge request will be referenced.
Issues will not be closed without a reason given.

Code
----

If you would like to contribute to ``aiapy``, you will first need to setup your development environment.

We suggest reading through the `SunPy developer's guide`_ for a more detailed description of the development process, as we follow there process for ``aiapy``.

This will hopefully make it easier for you to contribute to ``aiapy`` and other SunPy affiliated packages in the future.
If you encounter any problems, please don't hesitate to ask for help on any of our communication channels (found on the landing page of the documentation).

Testing
-------

There are several tests that require a working installation of `sswidl <http://www.lmsal.com/solarsoft/>`__ in order to compare results from IDL and Python.
This is managed via the `hissw <https://github.com/wtbarnes/hissw/>`__ package.
If you'd like to run these tests, you must first tell ``hissw`` where to find your IDL and SSW installations by placing the following lines in the file: ``$HOME/.hissw/hisswrc``,

.. code:: yaml

    [hissw]
    ssw_home=/path/to/ssw
    idl_home=/another/path/to/idl

where ``ssw_home`` is the path to the top of the sswidl tree and ``idl_home`` is the path to a working installation of IDL.
For more details, see the `hissw documentation <https://wtbarnes.github.io/hissw/>`__.
If a working installation is not available, these tests are automatically skipped.

Standards
---------

All contributors to the ``aiapy`` codebase should follow the `SunPy developer's guide`_.
This guide lays out a set of best practices for contributing, reviewing, testing, and documenting code.
All contributions to ``aiapy`` must adhere to the `Python in Heliophysics Community Standards <https://doi.org/10.5281/zenodo.2529130>`__.

.. _`SunPy developer's guide`: https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html
