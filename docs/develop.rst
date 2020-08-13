.. _dev-guide:

Developing aiapy
================

This page outlines how to setup your Python environment to develop
aiapy as well as some basic best practices for collaborative development.
These instructions are adapted from the
`SunPy developer's guide`_.
This guide is by no means complete or definitive. If you feel like something
is missing or incorrect, please
`create an issue <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/issues>`_

Issue Tracking
--------------
All bugs, feature requests, and other issues related to aiapy should be
recorded using the GitLab issue tracker. You can find instructions for
how to create an issue
`here <https://docs.gitlab.com/ee/user/project/issues/create_new_issue.html>`_.
All conversation regarding these issues should take place on the issue tracker.

Labels for issues should be created and used as appropriate. When a merge
request resolves an issue, the issue should be closed and the appropriate
merge request should be referenced. Issues should not be closed without a
reason.

Creating a Fork
---------------

If you would like to contribute to aiapy, you will first need to setup your
development environment. First create a fork of the main aiapy repository under
your GitLab username. You can find the instructions for how to do this
`here <https://docs.gitlab.com/ee/gitlab-basics/fork-project.html>`__.
If you don't already have an account on GitLab, you'll need to create one. You
can also sign into GitLab using your GitHub username.

Next, clone your fork of aiapy to your local machine,

.. code:: shell

    git clone https://gitlab.com/<your_username>/aiapy.git

Now add the main aiapy repository as an upstream repository,

.. code:: shell

    git remote add upstream https://gitlab.com/LMSAL_HUB/aia_hub/aiapy.git

You can now keep your fork up to date with main repository by running,

.. code:: shell

    git pull upstream master

Installation
-------------

If you're using the `Anaconda Python distribution <https://www.anaconda.com/distribution/#download-section>`_ (recommended),
create a new environment for aiapy development,

.. code-block:: shell

   conda create -n aiapy-dev python
   conda activate aiapy-dev

If you're using an alternate python installation, you can also use
`virtual environments <https://docs.python.org/3/tutorial/venv.html>`_.
Next, install the needed dependencies,

.. code-block:: shell

   cd aiapy
   pip install -e .[test,docs]

This includes all of the dependencies for the package as well as `pytest` for
running tests and `sphinx` for building the docs.

To make sure everything is working alright, you can run the tests,

.. code-block:: shell

    pytest --remote-data=any

See testing_ for more details regarding running the tests.

Making a Contribution
---------------------

If you want to add a feature or bugfix to aiapy, start by first making sure the
master branch of your fork is up to date with the master branch of the main
repository (see above, this will help to prevent potential file conflicts).
Next, create a new branch and switch to it,

.. code:: shell

    git checkout -b my-new-feature

After you've made your changes, commit and push them up to GitLab,

.. code:: shell

    git add changed_file_1.py changed_file_2.py
    git commit -m "short description of my change"
    git push origin my-new-feature

Once you see the changes in GitLab, create a merge request against the main
aiapy repository. You can find instructions for how to do this
`here <https://docs.gitlab.com/ee/gitlab-basics/add-merge-request.html>`__.
Others will likely have comments and suggestions regarding your proposed
changes. You can make these changes using the instructions listed above.

At least one other aiapy developer must approve your changes before the code
can be merged. Additionally, all automated tests should pass and all
conversations should be resolved. Once these steps are complete, the code can
be merged and you can delete  your branch `my-new-feature`.

.. _testing:

Testing
-------

Before committing any changes, you should ensure that the all of the tests pass
locally. To run the tests,

.. code:: shell

    pytest --remote-data=any

This will generate report showing which tests passed and which failed (if any).
Dropping the `--remote-data` flag will skip tests that require a network
connection. aiapy uses the `pytest <https://pytest.org/en/latest/>`_ framework
for discovering and running all of the tests.

Additions to the codebase should be accompanied by appropriate
tests such that the test coverage of the entire package does not decrease.
You can check the test coverage by running,

.. code:: shell

    pytest --remote-data=any --cov aiapy

Additionally, the test suite, including the documentation build and code
style checks can be run with `tox <https://tox.readthedocs.io/en/latest/>`_. See the
`SunPy developer's guide`_ for more information on running the test suite with `tox`.

Tests should be added to the directory in the appropriate subpackage, e.g.
for  `calibrate`, the tests should be placed in `calibrate/tests`. Your
tests can be added to an existing file or placed in a new file following
the naming convention `test_*.py`. This organization allows the tests to
be automatically discovered by pytest.

There are several tests that require a working installation of `sswidl <http://www.lmsal.com/solarsoft/>`_
in order to compare results from IDL and Python. This is managed
via the `hissw <https://github.com/wtbarnes/hissw/>`_ package.
If you'd like to run these tests, you must first tell `hissw`
where to find your IDL and SSW installations by placing the
following lines in the file `$HOME/.hissw/hisswrc`,

.. code:: yaml

    [hissw]
    ssw_home=/path/to/ssw
    idl_home=/another/path/to/idl

where `ssw_home` is the path to the top of the sswidl tree and `idl_home` is
the path to a working installation of IDL. For more details, see
the `hissw documentation <https://wtbarnes.github.io/hissw/>`_. If a working
installation is not available, these tests are automatically skipped.

Documentation
--------------

All documentation is written in `reStructuredText <https://docutils.sourceforge.io/rst.html>`_
and rendered using `Sphinx <https://www.sphinx-doc.org/en/master/>`_.
Docstrings are automatically pulled from all modules, functions and classes to
create the API documentation.
You can build and test the documentation locally by running,

.. code:: shell

    cd docs
    make html

This will run Sphinx on the restructured text files in order to create the
HTML version of the documentation.
The built documentation, in HTML format, is in `docs/_build/html`.

Best Practices
--------------

All contributors to the aiapy codebase should follow the `SunPy developer's guide`_.
This guide lays out a set of best practices for contributing, reviewing,
testing, and documenting code. All contributions to aiapy must adhere to the
`Python in Heliophysics Community Standards <https://doi.org/10.5281/zenodo.2529130>`_.

.. _`SunPy developer's guide`: https://docs.sunpy.org/en/latest/dev_guide/index.html
