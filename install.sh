# This works on Ubuntu 18.04

# Try to install python/pip/cython only if they are not already installed
python --version || sudo apt-get -y install python
pip --version || sudo apt -y install python-pip
cython --version || sudo apt -y install cython

# Get fastPHASE executable
curl http://scheet.org/code/Linuxfp.tar.gz --output Linuxfp.tar.gz
tar -xvzf Linuxfp.tar.gz
chmod a+x fastPHASE
rm Linuxfp.tar.gz

# Pip install everything else, but run pip as user, not superuser
su -c "pip install -r snpko/requirements.txt --user --upgrade" $SUDO_USER
