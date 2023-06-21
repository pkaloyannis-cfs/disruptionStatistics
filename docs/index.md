# Welcome to Disruption Statistics

Disruption Statistics is a software package by Panagiotis Kaloyannis for exploring the operational spaces of tokamaks and producing disruptivity maps.
This documentation page will be filled out more in the future. 😃

<!--- TODO: FIX THE CODE DUPLICATION IN THIS MD FILE AND THE DOCS INDEX SITE. CAN I DUPE THE README AS THE INDEX? SYMBOLIC LINK? -->

## Getting Started

The first step to getting started is ensuring that an adequate version of python (>3.8) is installed.
```
sudo apt install python-is-python3
```
Next, ensure you have installed poetry, as is seen [here](https://python-poetry.org/docs/). Once that is good and ready to go, install all the dependencies in one fell swoop using the following command.
```
poetry install
```
Finally, to ensure the code will work properly, feel free to run the unit tests.
```
poetry run pytest
```
If all those pass, you are ready to take on the world! The first plan of action would be to launch the jupyter notebook tutorials. For this run
```
poetry run jupyter notebook
```
and using the web GUI, open the `getting_started.ipynb` file. One could optionally also interact with this file by using another jupyter compatible editor (e.g. vscode), but just ensure that the python environment is the poetry environment generated during the poetry install step.

## Documentation

Generating code documentation is easy after this step, run the line below and click on the link to view docs!

```
poetry run mkdocs serve
```
Alternatively, you could build the website and browse it freely.
```
poetry run mkdocs build
```
To get the data used for data manipulation, consider looking at the `.mat` files found in [this google drive](https://drive.google.com/drive/folders/1ZeCJhrHs7nmj6mKEgLdV-hCJNjFZb2kk).



## For Developers

To contribute to this project, you will need access to `pre-commit` to handle linting of the python code. The above installation steps are the same, but you must initialize pre-commit before making any changes to the code.
```
pre-commit install
```
Now, on each commit the linter will run. You can also manually run the linter using.
```
pre-commit run --all-files
```


## Project layout

    .pre-commit-config.yaml         # The pre-commit hooks.
    .gitignore                      # Git Ignore
    .github/                        # Git Actions
    mkdocs.yml                      # The documentation configuration file.
    conftest.py                     # Global Unit test configuation.
    poetry.lock                     # Poetry Lock File
    pyproject.toml                  # Poetry TOML File
    README.md                       # The Repository README.
    docs/                           # Documentation Directory
    src/                            # Source Code Directory
        vis/                        # Visualization Package
        tokamaks/                   # Tokamak Configuration Directory
    tests/                          # Unit Test Directory
    notebooks/                      # Jupyter Notebook Directory

