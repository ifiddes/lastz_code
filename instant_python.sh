#!/usr/bin/env bash
# instant_python.sh: instantly install Python 2.7 with pip

# Die on error
set -e

# Where do we want Python?
PYTHON_PREFIX=$HOME/python

# Make a directory to build in
BUILD_DIR=`mktemp -d`

cd $BUILD_DIR

# Install Python
wget "http://python.org/ftp/python/2.7.5/Python-2.7.5.tgz" --no-check-certificate
tar -xvzf Python-2.7.5.tgz
cd Python-2.7.5
./configure --prefix=$PYTHON_PREFIX
make
make install

# Add it to the PATH
echo >> ~/.bashrc
echo "# Instant Python Installation" >> ~/.bashrc
echo "export PATH=$PYTHON_PREFIX/bin:\$PATH" >> ~/.bashrc

# Add it to the current script's path
export PATH=$PYTHON_PREFIX/bin:$PATH

cd ..

# Install setuptools
wget "https://pypi.python.org/packages/source/s/setuptools/setuptools-0.9.8.tar.gz" --no-check-certificate
tar -xvzf setuptools-0.9.8.tar.gz
cd setuptools-0.9.8
python setup.py install

cd ..

# Install pip
wget "https://raw.github.com/pypa/pip/master/contrib/get-pip.py" --no-check-certificate
python get-pip.py

# Add .local/bin to PATH for pip things
echo "export PATH=$HOME/.local/bin:\$PATH" >> ~/.bashrc


# Get rid of the build directory
cd
rm -Rf "$BUILD_DIR"

echo "Installed Python"

