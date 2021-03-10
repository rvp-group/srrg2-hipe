# HiPE Docker

This directory contains the docker file used to generate the *hipe* image. We also provide a *shell* script that you can use to build the image (the first time that you run it) and run the container. You can provide to the script a path to the folder that contains your datasets:
```
$ . docker.sh <path/to/datasets>
```

*Warning* The image contain extra packages which are not required by *srrg2_hipe*. We add those to provide some additional functionality within the container (editing, uncompressing ecc)

If you are wondering what a *docker* is, you can check it out [here](https://docs.docker.com/) 
