FROM ubuntu:20.04

# Install basic dependencies
RUN apt-get update && apt-get install -y \
    openjdk-11-jre-headless \
    curl \
    unzip \
    wget \
    python3 \
    python3-pip \
    git

# Install SnpEff
RUN mkdir -p /opt/snpEff && \
    wget -O /opt/snpEff/snpEff.zip https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    cd /opt/snpEff && unzip snpEff.zip && \
    chmod +x snpEff/snpEff.jar

# Install VEP
RUN mkdir -p /opt/vep && \
    cd /opt/vep && \
    wget -O vep.tar.gz https://ftp.ensembl.org/pub/release-104/variation/scripts/variant_effect_predictor/variant_effect_predictor.tar.gz && \
    tar xzf vep.tar.gz && \
    cd variant_effect_predictor && \
    perl INSTALL.pl -a a -s homo_sapiens -y GRCh38

# Install Python dependencies
RUN pip3 install pandas seaborn matplotlib

# Set environment variables
ENV PATH="/opt/snpEff:$PATH"
ENV PATH="/opt/vep:$PATH"

# Copy analysis scripts
COPY scripts/ /opt/analysis_scripts/

# Set working directory
WORKDIR /opt/analysis_scripts

# Default command
CMD ["/bin/bash"]
