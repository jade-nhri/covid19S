# covid19S
A bioinformatic protocol for calling spike protein variants of SARS-CoV-2

**To run with Docker**

``git clone https://github.com/jade-nhri/covid19S.git``

``cd covid19S``

``docker build -t "covid19s:v1" ./``

``docker run -h covidS --name covidS -i -t -v /:/MyData covid19s:v1 /bin/bash``

Please note that you need to run “vdb-config --interactive” before using sra toolkit.

Installation
------------
**Installation from source**

``cd /opt``

``git clone https://github.com/jade-nhri/covid19S.git``

``cd covid19S/covid19S``

``chmod +x *.py``

``export PATH="$PATH:/opt/covid19S/covid19S"``


## Dependencies

- [pyspoa-0.0.3](https://github.com/nanoporetech/pyspoa)
- [medaka-1.2.6](https://github.com/nanoporetech/medaka)
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](http://github.com/samtools/)
- [bcftools](https://github.com/samtools/bcftools)
- [seqkit=0.15](https://github.com/shenwei356/seqkit)
- [sratoolkit](https://github.com/ncbi/sra-tools)


## Usage
- [readme](https://www.dropbox.com/s/vz1xb8ywsotgcyw/Manual%20of%20covidS.pdf?dl=0)


