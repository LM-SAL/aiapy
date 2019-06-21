Developing aiapy
================

Issue Tracking
--------------
All bugs, feature requests, and other issues related to aiapy should be
recorded using the GitLab issue tracker. You can find instructions for
how to create an issue `here <https://docs.gitlab.com/ee/user/project/issues/create_new_issue.html>`_.
All conversation regarding these issues should take place on the issue tracker.

Labels for issues should be created and used as appropriate. When a merge request
resolves an issue, the issue should be closed and the appropriate merge request
should be referenced. Issues should not be closed without a reason.

Creating a Fork
---------------

If you would like to contribute to aiapy, you will first need to setup your
development environment. First create a fork of the main aiapy repository under
your GitLab username. You can find the instructions for how to do this
`here <https://docs.gitlab.com/ee/gitlab-basics/fork-project.html>`_.
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

Making a Contribution
---------------------

If you want to add a feature or bugfix to aiapy, start by first making sure the
master branch of your fork is up to date with the master branch of the main
repository (see above, this will help to prevent potential file conflicts). Next,
create a new branch and switch to it,

.. code:: shell

    git checkout -b my-new-feature

After you've made your changes, commit and push them up to GitLab,

.. code:: shell

    git add changed_file_1.py changed_file_2.py
    git commit -m "short description of my change"
    git push origin my-new-feature

Once you see the changes in GitLab, create a merge request against the main
aiapy repository. You can find instructions for how to do this `here <https://docs.gitlab.com/ee/gitlab-basics/add-merge-request.html>`_.
Others will likely have comments and suggestions regarding your proposed changes.
You can make these changes using the above "add-commit-push" instructions listed
above. 

At least one other aiapy developer must approve your changes before the code can be
merged. Additionally, all automated tests should pass and all conversations should be
resolved. Once these steps are complete, the code can be merged and you can delete 
your branch `my-new-feature`.

Best Practices
--------------

All contributors to the aiapy codebase should follow
`SunPy developer's guide <https://docs.sunpy.org/en/latest/dev_guide/index.html>`_.
This guide lays out a set of best practices for contributing, reviewing, testing,
and documenting code. All contributors to aiapy must also follow the
`Python in Heliophysics Community Standards <https://doi.org/10.5281/zenodo.2529130>`_.
