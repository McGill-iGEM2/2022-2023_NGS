# CONFIGURING PCGR USING CVMFS
echo -e "\n## GenPipes/MUGQIC genomes and modules" >> ~/.bashrc
echo "export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6" >> ~/.bashrc
echo "module use $MUGQIC_INSTALL_HOME/modulefiles" >> ~/.bashrc

# INSTALLING DATA BUNDLE
cd pcgr

GENOME="grch38" # or "grch37"
BUNDLE_VERSION="20220203"
BUNDLE="pcgr.databundle.${GENOME}.${BUNDLE_VERSION}.tgz"

wget http://insilico.hpc.uio.no/pcgr/${BUNDLE}
gzip -dc ${BUNDLE} | tar xvf -