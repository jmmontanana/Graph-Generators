//Author: J.M. Montañana, October 2019
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include "grafos.h"
#include "write_r.h"
#include "tamano_com.h"
char method[] = "duo";
char folder_base_name[] = "data_duo";

//#include <config.h>
#include <math.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

//if windows #include <iostream>

unsigned int rand_popular_node(const unsigned int num_nodes, const unsigned long long int max_rand, const unsigned int last, unsigned int* vector_nodes_available) {
	if (max_rand == 0) {
		printf(" %lld \n", max_rand); fflush(stdout);
		exit(0);
	}
	//default definition of max_rand is:
	// unsigned int max_rand = 0;//(2<<(num_nodes+1))-1;// es el sumatorio desde i=0 hasta num_nodes-1 de (2<<i)
	// for(node=0;node<num_nodes;node++){
	// max_rand+= (2<<(num_nodes-node));
	// vector_nodes_available[node]=node;
	// }
	unsigned long long int opt = rand() % max_rand;
	//buscamos a que nodo corresponde
	unsigned int node = 0;
	unsigned long long int accum = 0;
	unsigned int encontrado = false;
	while ((node < last) && (encontrado == false)) {
		accum += ((unsigned long long int)2 << (num_nodes - vector_nodes_available[node]));
		if (opt < accum) {
			encontrado = true;
			//printf(" acum %i node %i max_rand%i\n", accum, vector_nodes_available[node], max_rand);
		}
		else {
			node++;
		}
	}
	if (node == last) {
		unsigned long long int mio = 0;
		for (unsigned int i = 0; i < last; i++)
			mio += ((unsigned long long int)2 << (num_nodes - vector_nodes_available[i]));
		printf("valor correcto de max_rand seria %llu != max_rand %llu\n", mio, max_rand);
		printf("Error no encotramos nodo node %i=last %i\n", node, last);
		printf(" opt %llu max_rand %llu\n", opt, max_rand);
		exit(1);
	}
	return node;
}

/* The probability density function of a Weibull random variable is:
   f(x,L,k) = (k/L) (x/L)^(k-1) exp(-(x/L)^k) with x>=0, 0 if x<0
   https://en.wikipedia.org/wiki/Weibull_distribution
 */
double weibull(const double x, const double L, const double k) {
	if (x < 0) {
		return 0;
	}
	else if (x == 0) {
		if (k == 1)
			return 1 / L;
		else
			return 0;
	}
	else if (k == 1) {
		return exp(-x / L) / L;
	}
	else {
		//      return (k/L) * exp (-pow (x/L, k) + (k - 1) * log (x/L));
		return (k / L) * pow(x / L, k - 1) * exp(-pow(x / L, k));
	}
}

unsigned long long int popularidad_de(const unsigned int node, const unsigned int num_nodes) {
	return ((unsigned long long int)2 << (num_nodes - node));
}


int main(int argc, char* argv[]) // argv[i] from i = 0; to i<argcv;
{
	unsigned int  enlace;
	unsigned int num_nodes = 0;
	unsigned int num_communities = 8;
	//float weight_factor_intra_comm=1.0;
	unsigned int myseed = 100;
	unsigned int unnormalized = false;
	for (int i = 1; i < argc; i++) {//i=0 es el nombre del programa
		if (argv[i][0] == '-') {
			if (strcmp(argv[i], "-n") == 0) {
				if (argc > i) {
					printf("parametro nodos es '%s'\n", argv[i + 1]);
					num_nodes = atoi(argv[i + 1]);
				}
				else {
					printf("error parametro de numero de nodos sin valor asignado\n"); exit(1);
				}
				//} else if(strcmp(argv[i],"-l")==0){
				//	if(argc>i){
				//		printf("parametro enlaces es '%s'\n",argv[i+1]);
				//	}else{
				//		printf("error parametro de numero de enlaces sin valor asignado\n");exit(1);
				//	}
			}
			else if (strcmp(argv[i], "-c") == 0) {
				if (argc > i) {
					printf("parametro comunidades es '%s'\n", argv[i + 1]);
					num_communities = atoi(argv[i + 1]);
				}
				else {
					printf("error parametro de comunidades sin valor asignado\n"); exit(1);
				}
			}
			else if (strcmp(argv[i], "-s") == 0) {
				if (argc > i) {
					printf("parametro semilla es '%s'\n", argv[i + 1]);
					myseed = atoi(argv[i + 1]);
				}
				else {
					printf("error parametro semilla sin valor asignado\n"); exit(1);
				}
			// } else if(strcmp(argv[i],"-w")==0){
			// 	if(argc>i){
			// 		printf("parametro Weight_factor on the vertices between communities is: '%s'\n",argv[i+1]);
			// 		weight_factor_intra_comm=atof(argv[i+1]);
			// 	}else{
			// 		printf("error: missing value for the Weight_factor.\n");exit(1);
			// 	}	
			}else if (strcmp(argv[i], "-u") == 0) {
				printf("set mode as unnormalized.\n");
				unnormalized = true;
			}
			else if (strcmp(argv[i], "-h") == 0) {
				printf(" Usage: duo [options]\n\n Options:\n");
				printf("   -n number_of_vertices(nodes)\n   -c number_of_communities\n   -s value_of_rand_seed\n");
				printf("   -u Unnormalized output weights\n\n");
				printf("   -h Shows this usage information\n\n");
				printf("   number_of_vertices = 0, implies that size of the communities based on the ratio of the SED-graph\n");
				printf("   number_of_vertices > 0, implies that all the communia fixed size for all the communities.\n\n");

				printf(" The use of anyone of the parameters is optional, it is used the default value when they are not provided\n\n");
				printf(" The default values are:\n");
				printf("   number_of_vertices = 0, number_of_communities = 8, seed = 100\n");
				exit(1);
			}
		}
	}

	char folder_name[255];
	char mode[2] = { unnormalized ? 'u' : 'n', '\0' };

	sprintf_s(folder_name, "%s_%d_%d_%s", folder_base_name, num_communities, num_nodes, mode);

	srand(myseed);//((unsigned) time(&t));
	//definimos la estructura de comunidades:
	//una comunidad tiene un total de enlaces de: num_links*num_nodes
	_link*** comunidades; //es un array de comunidades, cada comunidad es un array de enlaces.
	comunidades = (_link***)malloc(num_communities * sizeof(_link**));
	if (comunidades == NULL) {
		fprintf(stderr, "Error: No se pudo asignar memoria para comunidades.\n");
		exit(EXIT_FAILURE);
	}
	unsigned int* comunity_num_nodes = (unsigned int*)malloc(num_communities * sizeof(unsigned int));
	if (comunity_num_nodes == NULL) {
		fprintf(stderr, "Error: No se pudo asignar memoria para comunity_num_nodes.\n");
		exit(EXIT_FAILURE);
	}
	unsigned int** comunity_size_node_links = (unsigned int**)malloc(num_communities * sizeof(unsigned int*));
	if (comunity_size_node_links == NULL) {
		fprintf(stderr, "Error: No se pudo asignar memoria para comunity_size_node_links.\n");
		exit(EXIT_FAILURE);
	}
	unsigned int** total_used_links = (unsigned int**)malloc(num_communities * sizeof(unsigned int*));
	if (total_used_links == NULL) {
		fprintf(stderr, "Error: No se pudo asignar memoria para total_used_links.\n");
		exit(EXIT_FAILURE);
	}
	unsigned int comm, node;
	unsigned int total_nodes = 0;
	for (comm = 0; comm < num_communities; comm++) {
		if (num_nodes == 0) {
			comunity_num_nodes[comm] = get_new_community_size(num_communities);
			if (comunity_num_nodes[comm] == 0) {
				printf("error comunidad comm %i de tamano 0\n", comm);
				exit(1);
			}
		}
		else {
			comunity_num_nodes[comm] = num_nodes;
		}
		total_nodes += comunity_num_nodes[comm];
		total_used_links[comm] = (unsigned int*)malloc(comunity_num_nodes[comm] * sizeof(unsigned int));
		if (total_used_links[comm] == NULL) {
			fprintf(stderr, "Error: No se pudo asignar memoria para total_used_links[comm].\n");
			exit(EXIT_FAILURE);
		}
		comunity_size_node_links[comm] = (unsigned int*)malloc(comunity_num_nodes[comm] * sizeof(unsigned int));
		if (comunity_size_node_links[comm] == NULL) {
			fprintf(stderr, "Error: No se pudo asignar memoria para comunity_size_node_links[comm].\n");
			exit(EXIT_FAILURE);
		}
		comunidades[comm] = (_link**)malloc(comunity_num_nodes[comm] * sizeof(_link*));
		if (comunidades[comm] == NULL) {
			fprintf(stderr, "Error: No se pudo asignar memoria para comunidades[comm].\n");
			exit(EXIT_FAILURE);
		}
		for (node = 0; node != comunity_num_nodes[comm]; node++) {
			comunidades[comm][node] = NULL;
			total_used_links[comm][node] = 0;
			printf("asigned comm %i node %i\n", comm, node);
			comunity_size_node_links[comm][node] = 0;
		}
		//comunidades[comm][node].link_node = UINT_MAX;// UINT_MAX stands for undefined destination node
	}

	unsigned int num_estudiantes = 150;
	unsigned int long_estud_preferen = 5;
	unsigned int step, stud;
	unsigned int last;
	bool encontrado;
	unsigned int** estudiante_preferencias;
	estudiante_preferencias = (unsigned int**)malloc(num_estudiantes * sizeof(unsigned int*));
	if (estudiante_preferencias == NULL) { 
		fprintf(stderr, "Error: No se pudo asignar memoria para estudiante_preferencias.\n");
		exit(EXIT_FAILURE);
	}
	for (stud = 0; stud < num_estudiantes; stud++) {
		estudiante_preferencias[stud] = (unsigned int*)malloc(long_estud_preferen * sizeof(unsigned int));
		if(estudiante_preferencias[stud] == NULL) {
			fprintf(stderr, "Error: No se pudo asignar memoria para estudiante_preferencias[stud].\n");
			exit(EXIT_FAILURE);
		}
	}
	// la probabilidad del nodo x es p(x)
	// donde x es un valor entero entre 0 y numero-de-nodos-menos-1
	// Queremos definir una funcion que nos de aleatoriamente un valor para seleccionar un nodo
	// por ejemplo cuando el valor aleatorio este en el intervalo [0,p(0)[ corresponderá al nodo 0
	// cuando el valor aleatorio este en el intervalo [p(0),p(0)+p(1)[ corresponderá al nodo 1
	// cuando el valor aleatorio este en el intervalo [p(0)+p(1),p(0)+p(1)+p(2)[ corresponderá al nodo 2
	// definimos la funcion de probabilidad como: p(x)=cte*(1/2)^x
	// donde x es un valor entre 0 y numnodes-1
	// cte= 2^(num_nodes)
	// por ello la funcion de probabilidad queda como: p(x)=2^(num_nodes-x)
	// Con esta definicion, todos los valores de p(x) son valores enteros mayores de 1.
	// el valor aleatorio maximo a calcular es: SUM(desde i=0, hasta i=num_nodes-1) 2^(num_nodes-i) = 2^(num_nodes+1)-1
	// el valor aleatorio minimo a calcular es 0
	// vamos a usar una funcion aleatoria de numeros enteros
	// cuando seleccionemos un nodo N: quitaremos el nodo N de la lista-vector, y restaremos p(N) al valor aleatorio maximo antes de buscar el siguiente nodo.
	// ATENCION: el valor maximo 2^(num_nodes)-1 con 64 bits, implica que el valor maximo de nodos por comunidad es 63
	for (comm = 0; comm < num_communities; comm++) {
		unsigned int* vector_nodes_available = (unsigned int*)malloc(comunity_num_nodes[comm] * sizeof(unsigned int));
		if (vector_nodes_available == NULL) {
			fprintf(stderr, "Error: No se pudo asignar memoria para vector_nodes_available.\n");
			exit(EXIT_FAILURE);
		}
		for (stud = 0; stud < num_estudiantes; stud++)
			for (step = 0; step < long_estud_preferen; step++)
				estudiante_preferencias[stud][step] = UINT_MAX;// UINT_MAX stands for undefined destination node
		for (stud = 0; stud < num_estudiantes; stud++) {
			unsigned long long int max_rand = 0; //es el sumatorio desde i=0 hasta num_nodes-1 de (2<<i) no es (2<<(num_nodes+1))-1
			for (node = 0; node < comunity_num_nodes[comm]; node++) {
				max_rand += popularidad_de(node, comunity_num_nodes[comm]);//(2<<(num_nodes-node));
				vector_nodes_available[node] = node;
			}
			last = comunity_num_nodes[comm];

			for (step = 0; step != long_estud_preferen; step++) {
				if (max_rand > 0 && last > 0) {//elementos disponibles para conectar
					node = rand_popular_node(comunity_num_nodes[comm], max_rand, last, vector_nodes_available);
					if (node >= comunity_num_nodes[comm]){
						printf("Error node >=comunity_num_nodes[comm]\n");
						exit(1);
					}
					max_rand -= popularidad_de(vector_nodes_available[node], comunity_num_nodes[comm]);//(2<<(num_nodes-vector_nodes_available[node]));
					estudiante_preferencias[stud][step] = vector_nodes_available[node];// UINT_MAX stands for undefined destination node
					if (vector_nodes_available[node] >= comunity_num_nodes[comm]) {
						printf("Error node >=numnodes\n");
						exit(1);
					}
					if (node != last - 1 && node < comunity_num_nodes[comm])
						vector_nodes_available[node] = vector_nodes_available[last - 1];
					last--;
				}
				// printf("Estudiante %i preferencias:",stud);
				// for (step = 0; step != long_estud_preferen; step++)
				// 	printf(" %i",estudiante_preferencias[stud][step]);
				// printf("\n");
			}
		}
		//ahora ponemos las listas de los estudiantes en la estructura de comunidades
		for (stud = 0; stud < num_estudiantes; stud++) {
			if (long_estud_preferen > 1) {
				for (step = 0; step < long_estud_preferen - 1; step++) { //para cada vector de 5 nodos, tenemos en total 4 enlaces, contado cada par de nodos
					if (estudiante_preferencias[stud][step] != UINT_MAX) { // el vector es mas corto que long_estud_preferen
						if (estudiante_preferencias[stud][step] == estudiante_preferencias[stud][step + 1]) {
							printf("error, nodo repetido %i stu %i step %i+1 k1\n", estudiante_preferencias[stud][step + 1], stud, step);
							exit(1);
						}
						enlace = 0;
						encontrado = false;
						unsigned int nodo = estudiante_preferencias[stud][step];
						if (comm < num_communities && nodo < comunity_num_nodes[comm]) {
							while ((encontrado == false) && (enlace < comunity_size_node_links[comm][nodo])
								) {
								//Warning, comunidades[comm][nodo][enlace] is not defined if enlace >= comunity_size_node_links[comm][nodo]
								if ((comunidades[comm][nodo][enlace].link_node == estudiante_preferencias[stud][step + 1]) &&
									(comunidades[comm][nodo][enlace].link_comunity == comm)) {
									encontrado = true;
									comunidades[comm][nodo][enlace].link_weight_out += ((unsigned long long int)2 << (long_estud_preferen - step)) / 16.0;
								}
								else if (comunidades[comm][nodo][enlace].link_node == UINT_MAX) {
									encontrado = true;
									if (estudiante_preferencias[stud][step + 1] != UINT_MAX) {
										total_used_links[comm][nodo] = total_used_links[comm][nodo] + 1;
										comunidades[comm][nodo][enlace].link_node = estudiante_preferencias[stud][step + 1];
										comunidades[comm][nodo][enlace].link_comunity = comm;//esto es la comunidad del nodo destino
										comunidades[comm][nodo][enlace].link_weight_out = ((unsigned long long int)2 << (long_estud_preferen - step)) / 16.0;//el peso va por el numero de salto
									}
								}
								else {
									enlace++;
								}
							}
						}
						if (encontrado == false) {
							unsigned int total_links = comunity_size_node_links[comm][nodo];
							_link* previo = comunidades[comm][nodo];
							comunidades[comm][nodo] = (_link*)malloc(total_links + 50 * sizeof(_link));
							if (comunidades[comm][nodo]==NULL) {
								fprintf(stderr, "Error: No se pudo asignar memoria para comunidades[comm][nodo].\n");
								exit(EXIT_FAILURE);
							}
							for (unsigned int temp_link = 0; temp_link < total_links; temp_link++) {								
								comunidades[comm][nodo][temp_link].link_node = previo[temp_link].link_node;
								comunidades[comm][nodo][temp_link].link_comunity = previo[temp_link].link_comunity;
								comunidades[comm][nodo][temp_link].link_weight_out = previo[temp_link].link_weight_out;
							}
							for (unsigned int temp_link = total_links; temp_link < total_links + 50; temp_link++) {
								comunidades[comm][nodo][temp_link].link_node = UINT_MAX;
								comunidades[comm][nodo][temp_link].link_comunity = UINT_MAX;
							}
							free(previo);
							comunity_size_node_links[comm][nodo] = total_links + 50;
							enlace = total_links;
							if (estudiante_preferencias[stud][step + 1] != UINT_MAX)
								total_used_links[comm][nodo] = total_used_links[comm][nodo] + 1;
							comunidades[comm][nodo][enlace].link_node = estudiante_preferencias[stud][step + 1];
							comunidades[comm][nodo][enlace].link_comunity = comm;//esto es la comunidad del nodo destino
							comunidades[comm][nodo][enlace].link_weight_out = ((unsigned long long int)2 << (long_estud_preferen - step)) / 16.0;//el peso va por el numero de salto
						}
					}//
				}//for step
			}//long_estud_preferen
		}//for stud
		free(vector_nodes_available);
	}//for comm

	for (stud = 0; stud < num_estudiantes; stud++)
		free(estudiante_preferencias[stud]);
	free(estudiante_preferencias);


	//buscamos el nodo mas popular
	for (comm = 0; comm < num_communities; comm++) {
		double maxpeso = 0.0;
		unsigned int nodo_maxpeso_link = UINT_MAX;
		unsigned int nodo_maxpeso = 0;//si la comunidad con solo un nodo, el de mas peso sera el 0, y no se podra realizar el bucle for
		for (node = 0; node != comunity_num_nodes[comm]; node++) {
			double totalpeso = 0.0;
			for (unsigned int enlace = 0; enlace < comunity_size_node_links[comm][node]; enlace++) {
				if (comunidades[comm][node][enlace].link_node != UINT_MAX) {
					totalpeso += comunidades[comm][node][enlace].link_weight_out;
					if (totalpeso >= maxpeso) {
						maxpeso = totalpeso;
						nodo_maxpeso = node;
						nodo_maxpeso_link = enlace;
					}
				}
			}
		}
		printf("nodo con mas peso de salida es nodo %i enlace %i peso %f comm %i\n", nodo_maxpeso, nodo_maxpeso_link, maxpeso, comm);
		//buscamos nodos que han quedado completamente sin conectar, y los conectamos
		for (node = 0; node != comunity_num_nodes[comm]; node++) {
			if (comunity_size_node_links[comm][node] == 0) {
				printf(" error nodo %i sin conexion. comm %i\n", node, comm);
				unsigned int total_links = comunity_size_node_links[comm][node];
				_link* previo = comunidades[comm][node];
				comunidades[comm][node] = (_link*)malloc(total_links + 50 * sizeof(_link));
				for (unsigned int temp_link = total_links; temp_link < total_links + 50; temp_link++)
					comunidades[comm][node][temp_link].link_node = UINT_MAX;
				free(previo);
				comunity_size_node_links[comm][node] = total_links + 50;
				enlace = total_links;
				total_used_links[comm][node] = total_used_links[comm][node] + 1;
				comunidades[comm][node][enlace].link_node = nodo_maxpeso;
				comunidades[comm][node][enlace].link_comunity = comm;//esto es la comunidad del nodo destino
				comunidades[comm][node][enlace].link_weight_out = ((unsigned long long int)2 << (long_estud_preferen - 0)) / 16.0;//el peso va por el numero de salto
			}else {
				unsigned int conectado = 0;
				for (unsigned int enlace = 0; enlace < comunity_size_node_links[comm][node]; enlace++)
					if (comunidades[comm][node][enlace].link_node != UINT_MAX)
						conectado = 1;
				if (conectado == 0) {
					printf(" node %i[%i] no conected, create connection to node %i\n", node, comm, nodo_maxpeso);
					comunidades[comm][node][0].link_node = nodo_maxpeso;
					float newpeso = 2.0;
					//newpeso=newpeso/64.0;
					comunidades[comm][node][0].link_weight_out = newpeso;
				}
			}
		}//for node
	}//for comm

	unsigned int max_used_links = 0;
	unsigned int** comunity_used_links = (unsigned int**)malloc(num_communities * sizeof(unsigned int*));
	for (comm = 0; comm < num_communities; comm++)
		for (node = 0; node != comunity_num_nodes[comm]; node++)
			for (unsigned int enlace = 0; enlace < comunity_size_node_links[comm][node]; enlace++)
				if (comunidades[comm][node][enlace].link_node != UINT_MAX)
					if (max_used_links < enlace + 1) //if used the link number 0 then means we using at least 1 link
						max_used_links = enlace + 1;
	printf("max_used_links %i\n", max_used_links);
	for (comm = 0; comm < num_communities; comm++) {
		comunity_used_links[comm] = (unsigned int*)malloc(comunity_num_nodes[comm] * sizeof(unsigned int));
		if (comunity_used_links[comm] == NULL) {
			printf("Error Out of Memory comunity_used_links[comm]\n");
			exit(1);
		}
		for (node = 0; node != comunity_num_nodes[comm]; node++) {
			comunity_used_links[comm][node] = (comunity_size_node_links[comm][node] == 0)? 0: max_used_links;
		}
	}

	if (unnormalized == false)
		normaliza_pesos_salida(num_communities, comunidades, NULL, comunity_num_nodes, comunity_size_node_links, NULL);

	print_graph_double(num_communities, comunidades, NULL, comunity_num_nodes, comunity_used_links, NULL);
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// ahora conectamos enlaces entre comunidades

	char comm2[50];
	char net2[50];
	char nodes_by_comm_file[50];
	sprintf_s(comm2, "%s%s%s_%d.dat", folder_name, "/community_", method, myseed);
	sprintf_s(net2, "%s%s%s_%d_0.dat", folder_name, "/network_", method, myseed);
	sprintf_s(nodes_by_comm_file, "%s%s%s_%d.dat", folder_name, "/nodes_by_comm_", method, myseed);
	makedir(folder_name);
	write_r_file(num_communities, comunity_num_nodes, comm2, net2, nodes_by_comm_file, comunity_size_node_links, comunidades, NULL, NULL);

	_link*** comunidades_otra;
	comunidades_otra = (_link***)malloc(num_communities * sizeof(_link**));
	unsigned int** comunity_size_node_links_otra = (unsigned int**)malloc(num_communities * sizeof(unsigned int*));
	unsigned int** total_used_links_at_node = (unsigned int**)malloc(num_communities * sizeof(unsigned int*));
	double* comm_external_links = (double*)malloc(num_communities * sizeof(double));
	for (unsigned int comm = 0; comm < num_communities; comm++) {
		comunidades_otra[comm] = (_link**)malloc(comunity_num_nodes[comm] * sizeof(_link*));
		comunity_size_node_links_otra[comm] = (unsigned int*)malloc(comunity_num_nodes[comm] * sizeof(unsigned int));
		total_used_links_at_node[comm] = (unsigned int*)malloc(comunity_num_nodes[comm] * sizeof(unsigned int));
	}

	for (unsigned int comm = 0; comm < num_communities; comm++) {
		for (unsigned int nodo = 0; nodo < comunity_num_nodes[comm]; nodo++) {
			comunity_size_node_links_otra[comm][nodo] = 0; //tamaño de memoria reservada para comunidades_otra[comm][nodo]
			comunidades_otra[comm][nodo] = NULL;
		}
	}

	//se crean enlaces diferentes en cada caso de num enlaces entre comunidades
	unsigned int max_mu_extern_links = 9;
	unsigned int mu_extern_links;
	float extern_links;
	for (mu_extern_links = 1; mu_extern_links < max_mu_extern_links; mu_extern_links++) { //vamos de 10% 20% 30% ... 60% de enlaces extra
		printf(" mu_extern_links %.1f%%\n", mu_extern_links * 10.0);
		for (unsigned int comm = 0; comm < num_communities; comm++) {
			extern_links = 0.0;
			for (unsigned int nodo = 0; nodo < comunity_num_nodes[comm]; nodo++) {
				//int locallinks=0;
				//for(enlace=0;enlace <comunity_size_node_links[comm][nodo];enlace++)
				//	if(comunidades[comm][nodo][enlace].link_node !=UINT_MAX)
				//		locallinks++;
				extern_links += total_used_links[comm][nodo];//locallinks;
				//printf ("links from %u[%u] are: %u=%u\n", nodo, comm, locallinks,total_used_links[comm][nodo]);
				total_used_links_at_node[comm][nodo] = 0;
				//comunidades_otra[comm][nodo]=( _link *) malloc( (extern_links)*sizeof( _link ));
				//for(int enlace = 0; enlace < extern_links; enlace++){
				//	comunidades_otra[comm][nodo][enlace].link_node = UINT_MAX;// UINT_MAX stands for undefined destination node
				//	comunidades_otra[comm][nodo][enlace].link_comunity = UINT_MAX;// UINT_MAX stands for undefined destination node
				//}
			}
			printf(" community %i total linkxs %f, ext links %f\n", comm, extern_links, (float)(extern_links * (float)mu_extern_links) / 10.0);
			comm_external_links[comm] = (extern_links * (float)mu_extern_links) / 10.0;//total de enlaces nuevos para la comunidad, no sabemos cuantos tiene cada comunidad
		}
		// mu_extern_links stands for mu from 0.1 ..0.2 .. 1.0.   0.1 means add a 10% of links of a community with the other comunities
		interconnecta_comunidades_new(num_communities,
			comunidades, comunity_size_node_links,
			comunidades_otra, comunity_size_node_links_otra, comunity_num_nodes, total_nodes, comm_external_links,
			total_used_links_at_node, mu_extern_links);
		//interconnecta_comunidades( num_communities, comunidades_otra, comunity_size_node_links_otra, comunity_num_nodes, total_nodes);
		if (unnormalized == false)
			normaliza_pesos_salida(num_communities, comunidades, comunidades_otra, comunity_num_nodes, comunity_size_node_links, comunity_size_node_links_otra);
		print_graph_double(num_communities, comunidades, comunidades_otra, comunity_num_nodes, comunity_used_links, comunity_size_node_links_otra);
		sprintf_s(comm2, "%s%s%s_%d.dat", folder_name, "/community_", method, myseed);
		sprintf_s(net2, "%s%s%s_%d_%d.dat", folder_name, "/network_", method, myseed, mu_extern_links);
		sprintf_s(nodes_by_comm_file, "%s%s%s_%d.dat", folder_name, "/nodes_by_comm_", method, myseed);
		write_r_file(num_communities, comunity_num_nodes, comm2, net2, nodes_by_comm_file,
			comunity_size_node_links, comunidades, comunity_size_node_links_otra, comunidades_otra);
	}

	free(comm_external_links);
	for (unsigned int comm = 0; comm < num_communities; comm++) {
		//comunidades_otra[comm_src][node_src]=( _link *) malloc(total_links+50*sizeof( _link ));
		//for (int nodo = 0; nodo != comunity_num_nodes[comm]; nodo++){
		//comunity_size_node_links_otra[comm_src][node_src]=total_links+50;


		for (unsigned int nodo = 0; nodo < comunity_num_nodes[comm]; nodo++) {
			free(comunidades_otra[comm][nodo]);
			free(comunidades[comm][nodo]);
		}
		free(total_used_links[comm]);
		free(comunity_used_links[comm]);
		free(comunity_size_node_links[comm]);
		free(comunidades[comm]);
		free(comunidades_otra[comm]);
		free(comunity_size_node_links_otra[comm]);
		free(total_used_links_at_node[comm]);
	}
	free(total_used_links_at_node);
	free(total_used_links);
	free(comunidades_otra);
	free(comunity_size_node_links_otra);
	free(comunity_size_node_links);
	free(comunidades);

	free(comunity_used_links);
	free(comunity_num_nodes);
	printf("fin de programa\n");
	// Do not forget to free the memory when you are done:
	return 0;
}//end main
