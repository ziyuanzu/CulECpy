# CulECpy

The process for enzyme-constrained model construction.

## About

The pipeline was written and tested on Windows Visual Studio Code. The core libraries essential for the pipeline including: cobra, Pyomo and related packages. 

## Installation

1. create CulECpy environment using conda:

```shell
$ conda create -n ECMpy python=3.6.5
```

2. install related packages using pip:

```shell
$ conda activate CulECpy
$ pip install cobra==0.13.3
$ pip install -r requirements.txt
$ python -m ipykernel install --user --name CulECpy --display-name "CulECpy"
```

## Verifying multiple instances of synthetic consortia construction

**1.Construction of synthetic consortia for sakuranetin with different pathway assignment strategies.**

**code:**CulECpy_TestFlow_1.ipynb

**data:**new_idea_SKM_pre_meta

**2.Construction of synthetic consortia for curcumin with different pathway assignment strategies.**

**code:**CulECpy_TestFlow_2.ipynb

**data:**cur

**3.Construction of synthetic consortium of rosmarinic acid with different two/three-strain pathway design ideas.**

**code:**CulECpy_TestFlow_3.ipynb.ipynb

**data:**3-mass

**4.Construction of synthetic consortia for salidroside.**

**code:**CulECpy_TestFlow_4.ipynb

**data:**new_idea_AGGD

**5.Construction of synthetic consortia for cadaverine.**

**code:**CulECpy_TestFlow_5.ipynb

**data:**iCW773+iml1515

