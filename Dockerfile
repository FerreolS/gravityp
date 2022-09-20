FROM ferreol/gravitypipeline:release
LABEL maintainer "Ferreol Soulez <ferreol.soulez@univ-lyon1.fr>"
ARG BRANCH
RUN cd $HOME && wget -nv https://github.com/FerreolS/gravityp/archive/refs/heads/$BRANCH.zip  && \
    unzip $BRANCH.zip && rm -f $BRANCH.zip && \
    rm -f $GRAVITYKIT/gravity$GRAVITYVERSION/gravi/*.{h,c} && \
    rm -f $GRAVITYKIT/gravity-$GRAVITYVERSION/recipes/*.{h,c} && \
    mv -f gravityp-$BRANCH/gravi/*.{h,c} $GRAVITYKIT/gravity-$GRAVITYVERSION/gravi/. &&\
    mv -f gravityp-$BRANCH/recipes/*.c $GRAVITYKIT/gravity-$GRAVITYVERSION/recipes/. &&\
    cd $GRAVITYKIT/gravity-$GRAVITYVERSION/ && make && make install && make clean
