This docker image is designed to process data from various quality control tools, such as multi-qc, and generate simplified reports and alerts. 

It comes pre-installed with popular Python libraries for data analysis and visualization, including matplotlib, numpy, and pandas.


## Local Testing Commands

To build it locally, you can use the following command from the project root directory:

```
docker build containers/python -t qc-pipeline-python:dev
```


To open a shell on a Docker container created from this image, run the following command:

```
docker run -it --rm qc-pipeline-python:dev /bin/bash
```

Since the `--rm` option is specified, the container will be automatically deleted once you exit the shell.

## Update python dependencies

To update the Python dependencies, you can follow the instructions provided in the [how_to_update_python_dependencies.md](how_to_update_python_dependencies.md) file.