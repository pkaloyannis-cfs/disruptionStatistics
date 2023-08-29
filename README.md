# Disruption Statistics

Disruption Statistics is a codebase that aims to provide a tokamak agnostic computing framework for computing disruptibity values from time series data. Along with these computational routines, a plethora of plotting tools are available to conveniently visualize, save and share high quality figures. The repository is written in python with the intention of using it in jupyter notebooks or other similar apps.

This work is the result of an EPFL Master's Project by Peter Kaloyannis in colaboration with MIT PSFC. The project was supervised by Prof. Ambrogio Fasoli (EPFL), Dr. Cristina Rea (PSFC), and Dr. Alessandro Pau (EPFL) and is built off the work of Nadezda Yalynskaya of IAEA.

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
If all those pass, you are ready to take on the world! Generating documentation is easy after this step, run the line below and click on the link to view docs!
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


