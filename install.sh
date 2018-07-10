sudo apt-get -y install python
sudo apt -y install python-pip
sudo apt -y install cython
curl http://scheet.org/code/Linuxfp.tar.gz --output Linuxfp.tar.gz
tar -xvzf Linuxfp.tar.gz
chmod a+x fastPHASE
mv fastPHASE snpko/
rm Linuxfp.tar.gz
