# Building the LFRic Apps Documentation

This README describes how to build the lfric_apps documentation.

For first time users, please create a virtual environment to build the docs.

From the `lfric_apps/documentation` folder of the repository run:

```bash
# Create a virtual environment and install dependencies (in ./documentation):
cd documentation
python3.12 -m venv .venv
.venv/bin/pip install .

# Activate the environment:
source .venv/bin/activate

# Build HTML documentation (in ./build/html/):
make clean html
```

At the Met Office you can also run the following to deploy the html documents
directly into `~/public_html/lfric_apps/<branch>/`:

```bash
make clean deploy
```
