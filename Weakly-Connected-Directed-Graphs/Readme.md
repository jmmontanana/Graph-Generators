# Weakly-Connected-Directed-Graphs

This project contains a Graph Generator for Benchmarking Algorithms of Detection of Communities

## Motivation 

The existing community detection algorithms have difficulties to detect the communities in graphs like the SED-graph [[1]](#ref1).
Thus, we consider that a benchmarking graph generator with the properties of the SED-graph 
is needed in order to develop a detection algorithm. 

## Prerequisites

  Before you can proceed, please clone the repository:

  ```bash
  git clone git://github.com/jmmontanana/Graph-Generators.git
  ```

## Compilation Instructions

Then compile with gcc as appear in the "compilation-script.sh":

  ```bash
  gcc -o duo duo.c write_r.c grafos.c tamano_com.c -lm -Wall;
  gcc -o dbo dbo.c write_r.c grafos.c tamano_com.c -lm -Wall;
  ```
 
### Dependencies

    This project has been tested with Linux Ubuntu 16.04 and Ubuntu 17.10  and the following installed packages:

    | Component        | Homepage                 | Version   |
    |----------------- |------------------------  |---------  |
    | gcc              | https://gcc.gnu.org      | 5.4.0     | 
 
# Example of generation of Synthetic Graphs:

Examples of generation of graphs from the command line:

### with the DUO generator

```bash 
	for i in 100 101 102 103 104 105; do #different seeds
		./duo -c 30 -n 0 -s $i > /dev/null; 
	done;
```
### with the DBO generator

```bash
	for i in 100 101 102 103 104 105; do #different seeds
		./dbo -c 30 -n 0 -s $i > /dev/null; 
	done;
``` 

## Acknowledgment 

    This work has been financed by the Project "Complex Networks" from the IUMM of the U.P.V.

## Contributing
    Find a bug? Have a feature request?
    Please [create](https://github.com/jmmontanana/Graph-Generators/issues) an issue.

## Main Contributors

**Montañana, José Miguel HLRS**
+ [github/jmmontanana](https://github.com/jmmontanana)

**Hervas, Antonio, UPV**

**Soriano, Pedro Pablo, UPV**

## License

The description of the algorithms is in the process of revision before being published.
The authors kindly request that publication be referenced by those who use this generator.
We will give details of the publication as soon as it is accepted

Copyright (C) 2018 J.M. Montañana

[Apache License v2](LICENSE).

## Release History

| Date        | Version | Comment          |
| ----------- | ------- | ---------------- |
| 2018-02-02  | 1.0     | First public functional version | 

## References

  <a name="ref1"></a>1. A. Hervas, P. P. Soriano, A. Jimenez, J. Peinado, R. Capilla, J. M. Montañana, 
  Modeling Human Behavior: Individuals and Organizations, Nova Science Publishers, Inc., New York, 2017, Ch. 
  Applying a graph model for the Spanish public university system, pp. 9–24.

 
