# GravityPipeline
Dockerfile to build the latest VLTI/Gravity pipeline with the python tool

## Running pipeline
```
docker run -it --rm --pull always  -v $FOLDER:/work/data ferreol/gravitypipeline:release
 ```
  where:
  - `$FOLDER` is the folder were all the raw files are (`$PWD` for the current directory).
