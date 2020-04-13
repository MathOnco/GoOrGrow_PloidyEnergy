# Go or grow: Analysis of interaction of ploidy in nutrient-limiting environment

We investigate the impact of a nutrient-limiting environment on the co-evoling polyploidy cells in a petri dish using coupled partial differential equations (PDEs) that govern the spatio-temporal dynamics of the nutrient (energy), low ploidy (energy efficient and high ploidy (energy inefficient).

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
To make sure everything is running properly, go to src directory and type ```make```. If everything compiles, here is a working example

```
./go_or_grow 100.0 10.0 1.0 0.005 0.1 0.2 0.4 10.0 2.3 0.5 0.3 0.1 0.8 t.txt r.txt E.txt u.txt v.txt statistics.txt 1 1
```

## Authors

* **Gregory Kimmel**
* **Philipp Altrock**
* **Noemi Andor**

## License

## Acknowledgements