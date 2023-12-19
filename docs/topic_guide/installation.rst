.. _aiapy-topic-guide-installing:

************
Installation
************

This topic guide details how to get a working installation of Python and ``aiapy``.

Installing Python
=================

There are many ways to install Python, but even if you have Python installed somewhere on your computer we recommend following these instructions anyway.
That's because we will create a new Python environment.
As well as containing a Python installation, this environment provides an isolated place to install Python packages (like ``aiapy``) without affecting any other current Python installation.
If you already have Python and either ``conda`` or ``pip`` working you can skip the next section.

Installing miniforge
--------------------

If you don't already have a Python installation then we recommend installing Python with `miniforge <https://github.com/conda-forge/miniforge/#miniforge>`__.
This will install ``conda`` and automatically configure the default channel (a channel is a remote software repository) to be ``conda-forge``, which is where ``aiapy`` is available.

First, download the installer for your system and architecture from the links below:

.. grid:: 3

    .. grid-item-card:: Linux

        `x86-64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh>`__

        `aarch64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh>`__

        `ppc64le <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-ppc64le.sh>`__

    .. grid-item-card:: Windows
        :link: https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe

        `x86-64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe>`__

    .. grid-item-card:: Mac

        `x86-64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh>`__

        `arm64 (Apple
        Silicon) <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh>`__

Then select your platform to install miniforge:

.. tab-set::

    .. tab-item:: Linux & Mac
        :sync: platform

        Linux & Mac Run the script downloaded above, with
        ``bash <filename>``. The following should work:

        .. code-block:: console

            bash Miniforge3-$(uname)-$(uname -m).sh

        Once the installer has completed, close and reopen your terminal.

    .. tab-item:: Windows
        :sync: platform

        Double click the executable file downloaded from
        the links above.

        Once the installer has completed you should have a new "miniforge
        Prompt" entry in your start menu.

In a new terminal (miniforge Prompt on Windows) run ``conda list`` to test that the install has worked.

Installing aiapy
----------------

To install ``aiapy``, start by launching a terminal (under a UNIX-like system) or the miniforge Prompt (under Windows).
Now we will create and activate new virtual environment to install ``aiapy`` into:

.. code-block:: bash

    $ conda create --name aiapy
    # Only run the following two lines
    # if you have NOT installed miniforge or added conda-forge to your channels
    $ conda config --add channels conda-forge
    $ conda config --set channel_priority strict
    $ conda activate aiapy

In this case the environment is named 'aiapy'.
Feel free to change this to a different environment name.

The benefit of using a virtual environment is that it allows you to install packages without affecting any other Python installation on your system.
This also means you can work on multiple projects (research or coding) with different package requirements without them interfering with each other.

Now we have a fresh environment we can install ``aiapy``:

.. code-block:: bash

    $ conda install aiapy

This will install ``aiapy`` and all of its dependencies.
If you want to install another package later, you can run ``conda install <package_name>``.
