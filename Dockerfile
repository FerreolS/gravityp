FROM fedora:latest
LABEL maintainer "Ferreol Soulez <ferreol.soulez@univ-lyon1.fr>"
ENV GRAVITYKIT=gravity-kit-1.11.0
ENV GRAVITYVERSION=1.11.0
RUN dnf install dnf-plugins-core  libffi-devel java-latest-openjdk-devel wget  perl bzip2 gnuplot awk zlib-ng-compat-devel  libcurl-devel  git -y  && \
    cd $HOME  && \
    wget -nv  https://ftp.eso.org/pub/dfs/pipelines/instruments/gravity/$GRAVITYKIT.tar.gz && \
    tar xzf $GRAVITYKIT.tar.gz && rm -f $GRAVITYKIT.tar.gz  && \
    cd $GRAVITYKIT  && \
    sed -i '911d' install_pipeline && sed -i 's/\&confirm(/0 and \&confirm(/g' install_pipeline && \
    sed -i '0,/   -t STDIN ||/s/    -t STDIN ||/  -t STDIN;/g' install_pipeline  && \
    echo "n"  |  ./install_pipeline   && \
    rm -rf gravity-calib-$GRAVITYVERSION && rm -rf cfitsio-* esorex-* fftw-* cpl-* erfa-* gsl-* wcslib-* &&\
    cd gravity-$GRAVITYVERSION/ && make clean && \
    mkdir -p /work/data && ln -s /usr/local/calib/gravity-$GRAVITYVERSION /work/common_calibration && \
    cd $HOME  &&  wget  -nv https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  && \
    bash ./Miniconda3-latest-Linux-x86_64.sh -b  && rm -f Miniconda3-latest-Linux-x86_64.sh && \
    export PATH=$PATH:$HOME/miniconda3/bin  && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda install -y  -c conda-forge numpy matplotlib scipy astropy astroquery scikit-learn reportlab pdfrw svglib pypdf   && \
    cd $HOME && git clone https://gitlab.obspm.fr/gravity-devs/gravi_tools.git   && \
    export PATH=$PATH:$HOME/gravi_tools:$HOME/gravi_tools/gravi_shell:$HOME/gravi_tools/gravi_quicklook && \
    export PYTHONPATH=$HOME/gravi_tools:$PYTHONPATH && \
    dnf remove  git -y
WORKDIR /work/data
ENV PATH=/root/.local/bin:/root/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/root/gravi_tools:/root/gravi_tools/gravi_shell:/root/gravi_tools/gravi_quicklook:/root/miniconda3/bin
ENV PYTHONPATH=/root/gravi_tools
ENTRYPOINT bash
