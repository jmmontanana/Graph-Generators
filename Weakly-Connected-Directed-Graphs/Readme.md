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

### Usage

```
  Usage: dbo [options] 
  
  Options:
 
 -n number_of_vertices(nodes) 
 -c number_of_communities 
 -s value_of_rand_seed
 -h Shows this usage information

   number_of_vertices = 0, implies that size of the communities based on the ratio of the SED-graph
   number_of_vertices > 0, implies a fixed size for all the communities.

 The use of anyone of the parameters is optional, it is used the default value when they are not provided

 The default values are:
   number_of_vertices = 0, number_of_communities = 8, seed = 100
```

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
**Montañana, José Miguel, HLRS**
+ [github/jmmontanana](https://github.com/jmmontanana)

**Hervas, Antonio, UPV**

**Soriano, Pedro Pablo, UPV**

## License
The description of the algorithms is in the process of revision before being published.
The authors kindly request that publication be referenced by those who use this generator.
We will give details of the publication as soon as it is accepted.

Copyright (C) 2018-19 J.M. Montañana

[Apache License v2](LICENSE).

## Release History
| Date        | Version | Comment          |
| ----------- | ------- | ---------------- |
| 2018-02-02  | 1.0     | First public functional version | 

## References
  <a name="ref1"></a>1. A. Hervas, P. P. Soriano, A. Jimenez, J. Peinado, R. Capilla, and J. M. Montañana. (2017) 
  Modeling Human Behavior: Individuals and Organizations, Nova Science Publishers, Inc., New York,  Ch. 
  Applying a graph model for the Spanish public university system, pp. 9–24.
  
  
  <a name="ref2"></a>2. J.M. Montanana, A. Hervas, and P.P. Soriano-Jiménez. (2020) A Proposal for a Benchmark Generator of Weakly Connected Directed Graphs. Open Journal of Modelling and Simulation. Vol.8 No.1, January 2020, Article ID:  97276 
https://doi.org/10.4236/ojmsi.2020.81002

