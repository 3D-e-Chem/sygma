FROM continuumio/miniconda

MAINTAINER Lars Ridder <l.ridder@esciencecenter.nl>

RUN /opt/conda/bin/conda install -y -q -c https://conda.anaconda.org/rdkit rdkit && \
/opt/conda/bin/conda install -y nose && \
/opt/conda/bin/conda clean -y -s -p -t -l -i

ENV PATH /opt/conda/bin:$PATH

ADD . /sygma

WORKDIR /sygma
RUN /opt/conda/bin/python setup.py install

WORKDIR /

ENTRYPOINT ["sygma"]
