Installing the Spatial Transcriptomics pipeline
-----------------------------------------------

These are the general instructions for installing the st_pipeline from scratch
on you compute environment. All the commands can be performed as a user with no
elevated permissions.

We recommend to download and install Anaconda (https://www.anaconda.com/products/individual)

We then create a virtual environment from which we will run the pipeline in.
Type the following command:

	``conda create -n stpipeline python=3.9``

The name for the virtual environment that we have just created is specified by
the -n flag. Here is is called stpipeline, but this can be anything that you want
to name it. To run the pipeline, this virtual environment must be activated. To
activate the virtual environment, enter the following command:

	``source activate stpipeline``

Where stpipeline is the name of your virtual environment (here the virtual
environment is called stpipeline). To deactivate the virtual environment, type the
following command:

	``source deactivate``

You need to obtain the pipeline from github to use it. The following steps will
tell you how to perform this.

Change to your home directory

	``cd``

Clone the repository from github

	``git clone git://github.com/jfnavarro/st_pipeline.git``

Change into the st_pipeline directory

	``cd st_pipeline``

Activate the virtual environment (if not already active)

	``source activate pipeline``

Install the pipeline

	``python setup.py build``

	``python setup.py install``
	
Alternatively, you can simply install the pipeline using PyPy:

	``pip install stpipeline``

Now the pipeline is installed and ready to be run.
