#!/bin/bash 
set -e
echo "Installing ST EMR Pipeline and dependencies"
#sudo apt-get -y install python-dev python-numpy python-matplotlib python-setuptools build-essential
sudo apt-get install -y libboost-dev cmake
sudo easy_install -U pip
sudo easy_install -U cython
sudo easy_install -U distribute
#sudo pip install --upgrade pip #do not think we need this 
sudo pip install HTSeq
sudo pip install pysam
sudo pip install argparse
if [ ! -d ~/bin ]; then
        mkdir ~/bin
fi
#hadoop fs -copyToLocal s3n://stpipelinedev/bin/* ~/bin
hadoop fs -copyToLocal s3n://stpipelinedev/src/bowtie2-2.1.0-source.zip .
hadoop fs -copyToLocal s3n://stpipelinedev/src/findIndexes.tar.gz .
hadoop fs -copyToLocal s3n://stpipelinedev/src/pipeline.tar.gz .
unzip bowtie2-2.1.0-source.zip
cd bowtie2-2.1.0
make
cp bowtie2* ~/bin
cd ..
tar -zxvf findIndexes.tar.gz
cd findIndexes
mkdir build
cd build
cmake ../
make
cp findIndexes ~/bin
cd ..
cd ..
tar -zxvf pipeline.tar.gz
cd pipeline
sudo python setup.py build
sudo python setup.py install
export PATH=$PATH:~/bin
export PYTHONPATH=$PATH:~/pipeline
echo "Done"
