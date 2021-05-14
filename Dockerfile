FROM nvidia/cuda:11.1-devel-ubuntu20.04
ENV http_proxy=http://proxy.nhri.org.tw:3128
ENV https_proxy=http://proxy.nhri.org.tw:3128
ENV TZ=Asia/Taipei
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y \
wget dos2unix \
python3 python3-pip \
python3-dev \
libbz2-dev \
liblzma-dev \
libncurses5-dev \
libncursesw5-dev \
zlib1g-dev \
libcurl4-gnutls-dev libssl-dev \
cmake unzip git wget libz-dev vim autoconf python

#Download covid19S
WORKDIR /opt
RUN git clone https://github.com/jade-nhri/covid19S.git
WORKDIR /opt/covid19S/covid19S
RUN chmod +x *.py


#Install pyspoa & medaka1.2.6
WORKDIR /opt
RUN wget https://github.com/nanoporetech/pyspoa/releases/download/v0.0.3/pyspoa-0.0.3.tar.gz
RUN tar -xzvf pyspoa-0.0.3.tar.gz
WORKDIR /opt/pyspoa
RUN make
RUN python3 -m pip install --upgrade pip
RUN pip install tensorflow==2.2.0
RUN pip3 install biopython pandas lxml six ete3 medaka==1.2.6

#Install minimap2 
WORKDIR /opt 
RUN git clone https://github.com/lh3/minimap2 
WORKDIR /opt/minimap2 
RUN make 
#Install bgzip and tabix 
WORKDIR /opt 
RUN git clone https://github.com/samtools/htslib.git 
WORKDIR /opt/htslib 
RUN git submodule update --init --recursive 
RUN autoreconf -i 
RUN ./configure 
RUN make 
RUN make install 

WORKDIR /opt 
RUN git clone http://github.com/samtools/samtools.git 
WORKDIR samtools 
RUN make
RUN make install

#Install bcftools 
WORKDIR /opt 
RUN git clone https://github.com/samtools/bcftools.git 
WORKDIR /opt/bcftools 
RUN make 
RUN make install

#Install seqkit 0.15
WORKDIR /opt 
RUN wget https://github.com/shenwei356/seqkit/releases/download/v0.15.0/seqkit_linux_amd64.tar.gz 
RUN tar -zxvf seqkit_linux_amd64.tar.gz 
RUN cp seqkit /usr/local/bin

#Install sratoolkit
WORKDIR /opt
RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
RUN tar zxvf sratoolkit.current-ubuntu64.tar.gz

WORKDIR /
ENV PATH $PATH:/opt/covid19S/covid19S/:/opt:/opt/minimap2:/opt/htslib:/opt/bcftools:/opt/samtools-1.12:/opt/sratoolkit.2.11.0-ubuntu64/bin

