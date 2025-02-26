To update versions in the `requirements_full.txt` file from a `requirements_minimal.txt` file containing only package names without versions, you can follow these steps:

1. Create a and activate a new virtual environment using a tool like `virtualenv` or `venv`. This will ensure a clean and isolated environment for the update process.  Here is an example with venv:

```
python3 -m venv ./venv

source ./venv/bin/activate
```

Make sure to use the same python version as in the dockerfile (here 3.11)

2. Install the packages listed in the `requirements_minimal.txt` file using the command `pip install -r requirements_minimal.txt`. This will install the latest versions of the packages.

3. Use the `pip freeze` command to generate a list of installed packages and their versions in the virtual environment:
 `pip freeze > requirements_full.txt`.

4. Test that scripts relying on this image still works.

5. Deactivate the virtual environment. For venv you can simply use the command `deactivate`. 

6. Make sure to commit and push your changes.

