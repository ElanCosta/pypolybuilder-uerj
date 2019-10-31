# pyPolyBuilder

PyPolyBuilder is a python tool that generates MTFs for linear and branched polymers, macromolecules and supramolecular structures

based on predefined BBs. The program considers two main cases: (1) the preparation of MTBs for radial polymers such as dendrimers; and

(2) a general procedure for connecting BBs based on a BB connectivity file. Currently, the output formats are compatible with GROMOS

and GROMACS simulation packages.

### 1- Installation:
The recommended way to install pyPolyBuilder is using pip:
```sh
pip install pypolybuilder
```

Another option is to clone this repository by simply:

```sh
git clone https://vitor_horta@bitbucket.org/vitor_horta/pypolybuilder.git
```

### 2- Basic Usage

pyPolyBuilder currently supports two building modes: Dendrimer and Polymer.

The required inputs for each mode are detailed in the following subsections.

Also the help command can be used to list all available options:

```sh
pypolybuilder --help
```

#### 2.1 Dendrimer

The dendrimer building process is based on the following inputs:

--core (Core filename)

--ter (Terminal filename)

--inter (Intermediate filename)

--core (Core filename)

--params (Force Field parameters filename)

--ngen (Number of generations to be built)

#### 2.2 Polymer

### TODO

### How to cite pyPolyBuilder

Authors of scientific papers including results generated using pyPolyBuilder are encouraged to cite the following paper.
```xml

@article{pyPolyBuilder2017,
    author    = " Vitor A. C. Horta, Mayk Caldas Ramos, Bruno A. C. Horta ",
    title     = { {Pypolybuilder}: A Python Program for the Creation of Molecular Topologies for Branched Macromolecules, Polymers and Supramolecular Structures },
    pages    = { PAGENUMBER },
    volume    = { VOLUMENUMBER },
    month     = { MONTH },
    year      = { 2017 },
    journal   = { JOURNAL }
}
```