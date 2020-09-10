[![Generic badge](https://img.shields.io/badge/docs-passing-<green>.svg)](https://prcurran.github.io/pdb_superimposer/)
[![prcurran](https://circleci.com/gh/prcurran/pdb_superimposer.svg?style=shield)](https://app.circleci.com/pipelines/github/prcurran)



# Non Standard Dependencies


[BioPython](https://biopython.org/)

[NumPy](https://numpy.org/)

[SciPy](https://www.scipy.org/)

[tqdm](https://pypi.org/project/tqdm/)



# Test the Code


1. Clone this repository

        git clone https://github.com/prcurran/pdb_superimposer.git

2. Create a new environment

        conda env create -f environment.yml

3. Install this package

        conda activate super
        pip install .

4. Run the CDK2 example

        python example.py

5. Remove environment and clean up

        conda env remove --name super --all
        rm -rf pdb_superimposer
