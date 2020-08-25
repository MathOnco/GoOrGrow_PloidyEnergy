# Go or grow: Analysis of interaction of ploidy in nutrient-limiting environment

We investigate the impact of a nutrient-limiting environment on the co-evoling polyploidy cells in a petri dish using coupled partial differential equations (PDEs) that govern the spatio-temporal dynamics of the nutrient (energy), low ploidy (energy efficient) and high ploidy (energy inefficient).

## Getting Started

Before being able to execute the code you will need to install the armadillo package (http://arma.sourceforge.net/download.html) and have -std=c++11 or higher.

### Prerequisites

To run the code you will need to install:

```
Armadillo (version 9.86)
```

### Installing
On Mac using homebrew:

```
1. brew install armadillo
```

## Running tests
To make sure everything is running properly, go to src directory and type ```make```. If everything compiles, here is a working example that uses default parameters:

```
./go_or_grow -v 1
```

To give user defined values:

```
./go_or_grow -v 1 -i <inputfile>
```

where ```<inputfile>``` is a tab separated file that contains the parameter names and values:
```
# example input.in contents:
outFile	out
statisticsFile	statsfile
tFinal	100.0
R	10.0
R0	5.0
dt	0.01
dr	0.1
u0	0.1
v0	0.1
eta	100.0
a	3.4
xi_u	0.7
phi_u	0.3
xi_v	0.5
phi_v	0.75
dim	2
growOutsideR0	true
saveFiles	false
randomCellMotion	1e-2
```

## Running the R Markdown script "ModelCalibration_MEP-LINCS.Rmd" requires:

  1) downloading the contents of Git folder RMarkdown
  2) downloading SupplementaryData1,2 and 3 from the corresponding article published in Cancer Research (Kimmel et al., 2020)
  3) installing package "rmarkdown" in R: install.packages("rmarkdown")


```Calling `rmarkdown::render("ModelCalibration_MEP-LINCS.Rmd")` should then reproduce RMarkdown/ModelCalibration_MEP-LINCS.html locally.```

## Authors

* **Gregory Kimmel**
* **Philipp Altrock**
* **Noemi Andor**

## License

## Acknowledgements