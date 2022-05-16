FROM ferreol/gravitypipeline:release
LABEL maintainer "Ferreol Soulez <ferreol.soulez@univ-lyon1.fr>"
ARG BRANCH
RUN cd $HOME && wget -nv https://github.com/FerreolS/gravityp/archive/refs/heads/$BRANCH.zip  && \
    unzip $BRANCH.zip && rm -f $BRANCH.zip && \
    rm -f gravity-kit-1.5.0-6/gravity-1.5.0/gravi/*.{h,c} && \
    rm -f gravity-kit-1.5.0-6/gravity-1.5.0/recipes/*.{h,c} && \
    mv -f gravityp-$BRANCH/gravi/*.{h,c} gravity-kit-1.5.0-6/gravity-1.5.0/gravi/. &&\
    mv -f gravityp-$BRANCH/recipes/*.c gravity-kit-1.5.0-6/gravity-1.5.0/recipes/. &&\
    cd gravity-kit-1.5.0-6/gravity-1.5.0/ && make && make install && make clean
