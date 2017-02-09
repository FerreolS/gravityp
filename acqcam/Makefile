CPL_INCLUDES = -I/home/aamorim/gvpipeline-9.8/include 

all: libacqcam.a


libacqcam.a:
	gcc -c -O2 -DCCS_LIGHT -DCCS_LIGHT -fPIC -Wno-format-security -Wno-deprecated-declarations -I. ${CPL_INCLUDES} -Iinclude -Iclipm \
                         gvacqPupilTracker.c \
                         gvacqAberrationTracker.c \
                         gvacqFieldTracker.c \
                         gvoProcessImageAC.c \
                         gvoacqWrite.c \
                         gvoacqPreprocess.c \
                         gvacqBicubicSpline.c \
                         clipm/clipm_image_signal.c \
                         clipm/clipm_priv_optimize.c \
                         clipm/clipm_priv_matrix.c\
                         clipm/clipm_priv_array.c \
                         clipm/clipm_math.c \
                         clipm/clipm_priv_image.c\
                         clipm/clipm_priv_image_signal.c\
                         clipm/clipm_priv_error.c \
                         clipm/clipm_priv_vector.c \
                         clipm/clipm_centroiding.c \
                         clipm/clipm_priv_checks.c
	ar -r libacqcam.a  *.o
	ranlib libacqcam.a

clean:
	rm -f *.o
	rm libacqcam.a
