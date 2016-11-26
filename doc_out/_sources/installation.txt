Installing the Spatial Transcriptomics pipeline
-----------------------------------------------

These are the general instructions for installing the st_pipeline from scratch
on you compute environment. All the commands can be performed as a user with no
elevated permissions.

Check that you have Anaconda installed, and that it is using python 2.7

	``python -c "import sys; print sys.version"``

Example output:

	``2.7.12 |Anaconda custom (64-bit)| (default, Jul  2 2016,
	17:42:40) [GCC 4.4.7 20120313 (Red Hat 4.4.7-1)]``

Here one can see that the python version is 2.7.12, and that Anaconda is
installed.  This shows an example of output when Anaconda is not installed:

	``2.7.10 (default, Oct 23 2015, 19:19:21) [GCC 4.2.1 Compatible Apple LLVM 7.0.0
	(clang-700.0.59.5)]``

If Anaconda is not installed on your system, it can be downloaded from the
following link: Important! Download the python 2.7 version

https://www.continuum.io/downloads#linux

Direct link for wget or curl (at the time of writing)
https://repo.continuum.io/archive/Anaconda2-4.1.1-Linux-x86_64.sh

	``wget https://repo.continuum.io/archive/Anaconda2-4.1.1-Linux-x86_64.sh``

To run the installer, type the following command (the numbers might change
depending on the current version at the time of download):

	``bash Anaconda2-4.1.1-Linux-x86_64.sh``

During the installation you will be asked if you want to pre-append the Anaconda
path to your path. Choose yes.  If you do not choose yes, run the following
command:

	``echo export PATH="~/anaconda2/bin:$PATH" > ~/.bash_profile``

**OR**

Edit your ~/.bash_profile directly by entering the following text using your
favourite text editor:

	``export PATH="~/anaconda2/bin:$PATH"``

To make sure that the changes have been regestired by the shell you can either:

Log out from the server and then log back in

**OR**

Type
	``source ~/.bash_profile``

After doing this, we need to install a couple of packages using anaconda to use
the pipeline. Type the following commands:

	``conda install nomkl``

	``conda install pysam``

We then create a virtual environment from which we will run the pipeline in.
Type the following command:

	``conda create -n pipeline python=2.7 anaconda``

The name for the virtual environment that we have just created is specified by
the -n flag. Here is is called pipeline, but this can be anything that you want
to name it. To run the pipeline, this virtual environment must be activated. To
activate the virtual environment, enter the following command:

	``source activate pipeline``

Where pipeline is the name of your virtual environment (here the virtual
environment is called pipeline). To deactivate the virtual environment, type the
following command:

	``source deactivate``

You need to obtain the pipeline from github to use it. The following steps will
tell you how to perform this.

Change to your home directory

	``cd``

Clone the repository from github

	``git clone git://github.com/SpatialTranscriptomicsResearch/st_pipeline.git``

Change into the st_pipeline directory

	``cd st_pipeline``

Activate the virtual environment (if not already active)

	``source activate pipeline``

Install the pipeline

	``python setup.py build``

	``python setup.py install``

Now the pipeline is installed and ready to be run.
