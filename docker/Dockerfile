# This is the official rocker R Docker image.
FROM rocker/r-base

LABEL container.base.image="rocker/r-base"
LABEL software.name="UCSC Cell Browser"
LABEL software.version="2018-04-23.0"
LABEL software.description="UCSC Cell Browser"
LABEL software.website="https://github.com/maximilianh/cellBrowser"
LABEL software.documentation="https://github.com/maximilianh/cellBrowser"
LABEL software.license="Open Source"
LABEL tags="Genomics"

WORKDIR /home/

run apt-get --yes --force-yes update
run apt-get --yes --force-yes install libssh2-1-dev
run apt-get --yes --force-yes install libssl-dev
run apt-get --yes --force-yes install libcurl4-openssl-dev
run apt-get --yes --force-yes install bash
run apt-get --yes --force-yes install git
run apt-get --yes --force-yes install python
run apt-get --yes --force-yes install time

run echo 'install.packages("devtools")' >> /home/installRpackages.R
run echo 'library(devtools)' >> /home/installRpackages.R
run echo 'install_github("satijalab/seurat", ref = "3bd092a")' >> /home/installRpackages.R

run R CMD BATCH /home/installRpackages.R

# Pull the source code for cellBrowser 
run git clone https://github.com/maximilianh/cellBrowser.git

ENTRYPOINT ["/bin/bash"]
