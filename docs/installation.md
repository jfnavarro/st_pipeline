# Installation

Python 3.10, 3.11 or 3.12 is required.

## Requirements

The ST Pipeline requires [STAR][] installed in the system (minimum version 2.5.4 if you use a ST Pipeline version >= 1.6.0).

If you use anaconda you can install STAR with:

```bash
conda install -c bioconda star
```

The ST Pipeline needs a computer with at least 32GB of RAM (depending on the size of the genome) and 8 cpu cores.

## Dependencies

The ST Pipeline depends on some Python packages that will
be automatically installed during the installation process.
You can see them in the file `requirements.txt`

## From source

The source for `ST Pipeline` can be downloaded from the [Github repo][].

You can either clone the public repository:

``` console
git clone https://github.com/jfnavarro/stpipeline
```

Or download the [tarball][]:

``` console
curl -OJL https://github.com/jfnavarro/stpipeline/tarball/master
```

Once you have a copy of the source, you can install it with:

### Using Poetry

If you don't have [Poetry](https://python-poetry.org/docs/) installed
you can use the following command:

``` console
curl -sSL https://install.python-poetry.org | python -
```

Install the package:

``` console
poetry install
```

Now you can run the ST Pipeline:

``` console
poetry run st_pipeline_run --help
```

### Using Pip

If you don't have [pip][] installed, this [Python installation guide][]
can guide you through the process.

Install the package:

``` console
pip install .
```

You can also use the official PyPy repositories:

``` console
pip install stpipeline
```

Now you can run ST Pipeline:

``` console
st_pipeline_run --help
```

### Using Docker

Before installing, ensure that [Docker](https://www.docker.com/) is installed
in your environment.

First, build a Docker image:

``` console
docker buildx build --platform linux/amd64 -t stpipeline .
```

Then, you can run ST Pipeline using Docker:

To run `ST Pipeline` commands:

``` console
docker run --rm stpipeline st_pipeline_run --help
```

### Using Anaconda

Before installing, ensure you have either [Anaconda](https://www.anaconda.com/)
or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed in your environment.

First, create the environment:

``` console
conda env create -n stpipeline python=3.10
```

Then, activate the environment:

``` console
conda activate stpipeline
```

Install the package:

``` console
pip install .
```

Now you can run ST Pipeline:

``` console
st_pipeline_run --help
```

  [STAR]: https://github.com/alexdobin/STAR
  [pip]: https://pip.pypa.io
  [Python installation guide]: http://docs.python-guide.org/en/latest/starting/installation/
  [Github repo]: https://github.com/jfnavarro/st_pipeline/
  [tarball]: https://github.com/jfnavarro/st_pipeline/releases
