//Author: J.M. Monta�ana, April 2017
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include "grafos.h"



//ordenados de mayor a menor
void bubble_sort_linkout_accumulated(long list[], double** linkout_accumulated, const unsigned int n, const unsigned int comm) {
	long c, d, t;
	for (c = 0; c < (n - 1); c++) {
		for (d = 0; d < n - c - 1; d++) {
			if (linkout_accumulated[comm][list[d]] < linkout_accumulated[comm][list[d + 1]]) {
				/* Swapping */
				t = list[d];
				list[d] = list[d + 1];
				list[d + 1] = t;
			}
		}
	}
}


////ordenados de mayor a menor
//void bubble_sort(long list[], long n){
//  long c, d, t;
//  for (c = 0 ; c < ( n - 1 ); c++)  {
//    for (d = 0 ; d < n - c - 1; d++)    {
//      if (list[d] < list[d+1])      {
//        /* Swapping */
//        t         = list[d];
//        list[d]   = list[d+1];
//        list[d+1] = t;
//      }
//    }
//  }
//}

void normaliza_pesos_salida(const unsigned int num_communities, _link*** comunidades,
	_link*** comunidades_otra, unsigned int* comunity_num_nodes,
	unsigned int** comunity_size_node_links, unsigned int** comunity_size_node_links_otra) {
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//ahora ordenamos los nodos en cada comunidad por peso de los enlaces de salida, de mayor a menor.
	double** linkout_accumulated;
	long** list_pos;
	linkout_accumulated = (double**)malloc(num_communities * sizeof(double*));
	list_pos = (long**)malloc(num_communities * sizeof(long*));
	for (int comm = 0; comm < num_communities; comm++) {
		linkout_accumulated[comm] = (double*)malloc(comunity_num_nodes[comm] * sizeof(double));
		list_pos[comm] = (long*)malloc(comunity_num_nodes[comm] * sizeof(long));
	}
	for (int comm = 0; comm < num_communities; comm++) {
		for (int nodo = 0; nodo < comunity_num_nodes[comm]; nodo++) {
			linkout_accumulated[comm][nodo] = 0;
			for (int enlace = 0; enlace < comunity_size_node_links[comm][nodo]; enlace++) {
				if (comunidades[comm][nodo][enlace].link_node != UINT_MAX) {//	conectado
					linkout_accumulated[comm][nodo] += comunidades[comm][nodo][enlace].link_weight_out;
				}
			}
			list_pos[comm][nodo] = nodo;
		}
		if (comunidades_otra != NULL) {
			for (int nodo = 0; nodo < comunity_num_nodes[comm]; nodo++) {
				for (int enlace = 0; enlace < comunity_size_node_links_otra[comm][nodo]; enlace++) {
					if (comunidades_otra[comm][nodo][enlace].link_node != UINT_MAX) {//	conectado
						linkout_accumulated[comm][nodo] += comunidades_otra[comm][nodo][enlace].link_weight_out;
					}
				}
			}
		}
		bubble_sort_linkout_accumulated(list_pos[comm], linkout_accumulated, comunity_num_nodes[comm], comm);
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//escalamos los pesos para que el que tiene el mayor peso de enlaces de salida sea 75
	for (int comm = 0; comm < num_communities; comm++) {
		double factor_escala = 75.0 / linkout_accumulated[comm][list_pos[comm][0]];
		//printf("\n Communidad comm %i, total nodos: %i   factor_escala %f\n",comm,comunity_num_nodes[comm], factor_escala);
		for (int nodo = 0; nodo < comunity_num_nodes[comm]; nodo++) {
			linkout_accumulated[comm][list_pos[comm][nodo]] = linkout_accumulated[comm][list_pos[comm][nodo]] * factor_escala;
			//printf("   Nodo %li: link_out weight %f\n",list_pos[comm][nodo],linkout_accumulated[comm][list_pos[comm][nodo]]);
			for (int enlace = 0; enlace < comunity_size_node_links[comm][nodo]; enlace++) {
				if (comunidades[comm][nodo][enlace].link_node != UINT_MAX) {
					comunidades[comm][nodo][enlace].link_weight_out = comunidades[comm][nodo][enlace].link_weight_out * factor_escala;
				}
			}
			if (comunidades_otra != NULL) {
				for (int enlace = 0; enlace < comunity_size_node_links_otra[comm][nodo]; enlace++) {
					if (comunidades_otra[comm][nodo][enlace].link_node != UINT_MAX) {
						comunidades_otra[comm][nodo][enlace].link_weight_out = comunidades_otra[comm][nodo][enlace].link_weight_out * factor_escala;
					}
				}
			}
		}
	}
}



float peso_link(_link* nodo, const unsigned int enlace, const unsigned int dest, const int max_link) {
	if ((nodo == NULL) || (enlace >= max_link) || (dest == UINT_MAX)) {
		return 0.0;
	}
	else {
		return nodo[enlace].link_weight_out;
	}
}



float print_link(_link* nodo, const unsigned int enlace, const unsigned int dest, const int max_link, const int tabla_abreviada) {
	if (tabla_abreviada == false) {
		if ((nodo == NULL) || (enlace >= max_link)) {
			printf("     |         |");
			return 0.0;
		}
		else {
			if (dest == UINT_MAX) {
				printf(" -[-]|    (-)  |");
				return 0.0;
			}
			else if (nodo[enlace].link_weight_out == 0.0) {
				printf(" %u[%u]| ----.-- |",
					dest,
					nodo[enlace].link_comunity);
				return 0.0;
			}
			else {
				printf(" %u[%u]| %7.2f |",
					dest,
					nodo[enlace].link_comunity,
					nodo[enlace].link_weight_out);
				return nodo[enlace].link_weight_out;
			}
		}
	}
	else {
		if ((nodo == NULL) || (enlace >= max_link)) {
			printf("          |");
			return 0.0;
		}
		else {
			if (dest == UINT_MAX) {
				printf("      (-)  |");
				return 0.0;
			}
			else if (nodo[enlace].link_weight_out == 0.0) {
				printf("  ----.-- |");
				return 0.0;
			}
			else {
				printf("  %7.2f |",
					nodo[enlace].link_weight_out);
				return nodo[enlace].link_weight_out;
			}
		}
	}
}

//this function prints in screen the nodes in each community
void print_graph_double(const unsigned int num_communities, _link*** comunidades, _link*** comunidades_otra, unsigned int* comunity_num_nodes,
	unsigned int** comunity_size_node_links, unsigned int** comunity_size_node_links_otra) {
	unsigned int dest, comm, enlace, nodo, tabla_abreviada;

	tabla_abreviada = (comunidades_otra == NULL) ? true : false;

	//reservamos memoria para las metricas en total_peso_entradas y total_peso_salidas
	float** total_peso_entradas = (float**)malloc(num_communities * sizeof(float*));
	float** total_peso_salidas = (float**)malloc(num_communities * sizeof(float*));
	unsigned int dest_comm;
	for (comm = 0; comm < num_communities; comm++) {
		total_peso_entradas[comm] = (float*)malloc(comunity_num_nodes[comm] * sizeof(float));
		total_peso_salidas[comm] = (float*)malloc(comunity_num_nodes[comm] * sizeof(float));
		for (nodo = 0; nodo < comunity_num_nodes[comm]; nodo++) {
			total_peso_entradas[comm][nodo] = 0.0;//para los destinos
			total_peso_salidas[comm][nodo] = 0.0;
		}
	}
	//generamos las metricas en total_peso_entradas y total_peso_salidas
	for (comm = 0; comm < num_communities; comm++) {
		for (nodo = 0; nodo != comunity_num_nodes[comm]; nodo++) {
			if (comunidades != NULL) {
				for (dest = 0; dest < comunity_num_nodes[comm]; dest++) {
					enlace = 0;
					unsigned int found = false;
					while ((found == false) && (enlace < comunity_size_node_links[comm][nodo])) {
						if (dest == comunidades[comm][nodo][enlace].link_node) {
							found = true;
							float peso = peso_link(comunidades[comm][nodo], enlace, dest, comunity_size_node_links[comm][nodo]);
							total_peso_salidas[comm][nodo] += peso;
							dest_comm = comunidades[comm][nodo][enlace].link_comunity;
							if (dest >= comunity_num_nodes[dest_comm]) {
								printf("error dest %i[%i] > max %i\n", dest, dest_comm, comunity_num_nodes[dest_comm]); fflush(stdout);  exit(1);
							}
							total_peso_entradas[dest_comm][dest] += peso;
						}
						else {
							enlace++;
						}
					}
					if (found == false) {
						peso_link(NULL, 0, UINT_MAX, comunity_size_node_links[comm][nodo]);
					}
				}
				// for (enlace=0;enlace < comunity_size_node_links[comm][nodo];enlace++){
				// float peso=peso_link( comunidades[comm][nodo], enlace, dest, comunity_size_node_links[comm][nodo] );
				// total_peso_salidas[comm][nodo]+= peso;
				// total_peso_entradas[dest_comm][dest]+= peso;
				// }
			}
			if (comunidades_otra != NULL) {
				for (enlace = 0; enlace < comunity_size_node_links_otra[comm][nodo]; enlace++) {
					if (comunity_size_node_links[comm][nodo] != 0) { //not have defined links, and not reserved memory for it
						dest = comunidades_otra[comm][nodo][enlace].link_node;
						if ((dest != UINT_MAX) && (comm == comunidades_otra[comm][nodo][enlace].link_comunity)) {
							float peso = peso_link(comunidades_otra[comm][nodo], enlace, dest, comunity_size_node_links_otra[comm][nodo]);
							total_peso_salidas[comm][nodo] += peso;
							dest_comm = comunidades[comm][nodo][enlace].link_comunity;
							if (dest >= comunity_num_nodes[comm]) {
								printf("ERROR dest %i > max %i\n", dest, comunity_num_nodes[comm]); fflush(stdout); exit(1);
							}
							total_peso_entradas[dest_comm][dest] += peso;
						}
					}
				}
			}
		}//for nodo
	}//for comm

	// imprimimos las metricas por pantalla
	for (comm = 0; comm < num_communities; comm++) {
		printf(" COMUNIDAD :%i dest[comm](peso)\n", comm);
		if (comunidades_otra == NULL) {
			printf("                     |");
			for (dest = 0; dest != comunity_num_nodes[comm]; dest++) {
				printf("     %u[%u] |", dest, comm);
			}printf(" peso de salida\n");
		}

		printf("total dest =%i\n", comunity_num_nodes[comm]);
		for (nodo = 0; nodo != comunity_num_nodes[comm]; nodo++) {
			if (comunidades != NULL) {
				printf("   NODO %i (links %i): |", nodo, comunity_size_node_links[comm][nodo]); fflush(stdout);
				for (dest = 0; dest < comunity_num_nodes[comm]; dest++) {
					enlace = 0;
					unsigned int found = false;
					while ((found == false) && (enlace < comunity_size_node_links[comm][nodo])) {
						if (dest == comunidades[comm][nodo][enlace].link_node) {
							found = true;
							print_link(comunidades[comm][nodo], enlace, dest, comunity_size_node_links[comm][nodo], tabla_abreviada);
						}
						else {
							enlace++;
						}
					}
					if (found == false) {
						print_link(NULL, 0, UINT_MAX, comunity_size_node_links[comm][nodo], tabla_abreviada);
					}fflush(stdout);
				}
			}
			if (comunidades_otra != NULL) {
				//printf("Nodo %i %i\n", nodo, comunity_size_node_links_otra[comm][nodo]);fflush(stdout);
				for (enlace = 0; enlace < comunity_size_node_links_otra[comm][nodo]; enlace++) {
					if (comunity_size_node_links[comm][nodo] == 0) { //not have defined links, and not reserved memory for it
						print_link(NULL, 0, UINT_MAX, comunity_size_node_links_otra[comm][nodo], false);
					}
					else {
						dest = comunidades_otra[comm][nodo][enlace].link_node;
						if ((dest != UINT_MAX) && (comm == comunidades_otra[comm][nodo][enlace].link_comunity)) {
							print_link(comunidades_otra[comm][nodo], enlace, dest, comunity_size_node_links_otra[comm][nodo], false);
							if (dest >= comunity_num_nodes[comm]) {
								printf("ERROR dest %i > max %i\n", dest, comunity_num_nodes[comm]); fflush(stdout); exit(1);
							}
						}
					}
				}
			}
			if (tabla_abreviada == false) printf(" peso de salida|");
			if (total_peso_salidas[comm][nodo] == 0.0) {
				printf("  ---.--\n");
			}
			else
				printf(" %7.2f\n", total_peso_salidas[comm][nodo]);
		}//for nodo
		printf("peso de entrada     :|");
		for (dest = 0; dest != comunity_num_nodes[comm]; dest++)
			if (total_peso_entradas[comm][dest] == 0.0) {
				printf("  ---.-- |"); fflush(stdout);
			}
			else
				printf("  %7.2f |", total_peso_entradas[comm][dest]);
		printf("\n--------------------------------------\n"); fflush(stdout);
	}//for comm

	for (comm = 0; comm < num_communities; comm++) {
		free(total_peso_entradas[comm]);
		free(total_peso_salidas[comm]);
	}
	free(total_peso_entradas);
	free(total_peso_salidas);
}

//this function prints in screen the nodes in each community
// void print_graph (const unsigned int num_communities, _link*** comunidades, unsigned int *comunity_num_nodes,
// unsigned int** comunity_size_node_links)
// {
// print_graph_double(num_communities, comunidades, NULL, comunity_num_nodes, comunity_size_node_links, NULL);
// }
void verify_comunities(const unsigned int num_communities, _link*** comunidades, unsigned int* comunity_num_nodes, unsigned int** comunity_size_node_links)
{
	//vamos a contar cuantos enlaces usa cada nodo, alguno puede haber quedado sin conectar. (cuenta_local)
	//vamos a contar cuantos nodos se conectan con �l. (cuenta_conectados)
	//ambos cuentas tienen que ser iguales.
	unsigned int total_nodes = 0;
	for (unsigned int comm = 0; comm < num_communities; comm++) {
		total_nodes += comunity_num_nodes[comm];
	}
	unsigned int* cuenta_local = (unsigned int*)malloc(total_nodes * sizeof(unsigned int));
	if (cuenta_local == NULL) {
		fprintf(stderr, "Error: No se pudo asignar memoria para cuenta_local.\n");
		exit(EXIT_FAILURE);
	}
	unsigned int* cuenta_conectados = (unsigned int*)malloc(total_nodes * sizeof(unsigned int));
	if (cuenta_conectados == NULL) {
		fprintf(stderr, "Error: No se pudo asignar memoria para cuenta_conectados.\n");
		exit(EXIT_FAILURE);
	}


	for (unsigned int comm = 0; comm < num_communities; comm++) {
		for (unsigned int nodo = 0; nodo < comunity_num_nodes[comm]; nodo++) {
			cuenta_local[nodo] = 0;//iniciamos el contador, es el siguiente nodo
			cuenta_conectados[nodo] = 0;
		}
		for (unsigned int nodo = 0; nodo != comunity_num_nodes[comm]; nodo++) {
			for (unsigned int enlace = 0; enlace < comunity_size_node_links[comm][nodo]; enlace++) {
				if (comunidades[comm][nodo][enlace].link_node != UINT_MAX) {
					cuenta_local[nodo] = cuenta_local[nodo] + 1;
					cuenta_conectados[comunidades[comm][nodo][enlace].link_node] = cuenta_conectados[comunidades[comm][nodo][enlace].link_node] + 1;
				}
			}
		}
		//int valido=true;
		for (int nodo = 0; nodo < comunity_num_nodes[comm]; nodo++)
			if (cuenta_local[nodo] != cuenta_conectados[nodo]) {
				//	valido=false;
				printf("comunidad %u erronea nodo %i tot1 %i tot2 %i\n", comm, nodo, cuenta_local[nodo], cuenta_conectados[nodo]);
				exit(ERR_COMUNIDAD);
			}
		//if (valido==true)
		printf("comunidad %u validada\n", comm);
	}
	free(cuenta_local);
	free(cuenta_conectados);
}



void remove_connection_comm(const unsigned int nodo, const unsigned int link, const unsigned int comm,
	_link*** comunidades, unsigned int** comunity_size_node_links) {
	unsigned int enlace = link;
	if (enlace >= comunity_size_node_links[comm][nodo]) {
		printf("Error remove_connection_comm_ini %u\n", enlace);
		exit(ERR_REM_CONNECTION);
	}
	do {// Eliminamos el enlace "enlace", con sobre-escribiendo con otro de indice menor, al final querar� libre la posicion 0
		if (enlace != 0) {
			comunidades[comm][nodo][enlace].link_node = comunidades[comm][nodo][enlace - 1].link_node;
			comunidades[comm][nodo][enlace].link_comunity = comunidades[comm][nodo][enlace - 1].link_comunity;
			comunidades[comm][nodo][enlace - 1].link_node = UINT_MAX;
			enlace--;
		}
		else
			comunidades[comm][nodo][enlace].link_node = UINT_MAX;
	} while (enlace != 0);
}


void print_vector_com(unsigned int** vector_com_nodes, unsigned int** vector_com_links, const unsigned int* last, const unsigned int num_communities)
{
	for (unsigned int community = 0; community < num_communities; community++) {
		printf("comunidad %u vector last %i:\n ", community, last[community]);
		for (int i = 0; i < last[community]; i++) {
			printf(" %u  ", vector_com_nodes[community][i]);
			printf(" links(%u);", vector_com_links[community][i]);
		}
		printf(";last=%u\n", last[community]);
	}
}
void remove_from_vector_comm(unsigned int** vector_com_nodes, unsigned int** vector_com_links, unsigned int* last,
	const unsigned int node_a, const unsigned int comm_a,
	const unsigned int node_b, const unsigned int comm_b) {
	//quitamos los valores en la posicion "node_a" y posicion "node_b" del vector temporal
	// atencion, si borramos un elemento de la cadena, los elementos siguientes avanzan una posicion.
	//por eso (en el caso que sean de la misma comunidad) borramos primero el elemento mas al final de la cadena,
	// y despues el mas cercano al comienzo
	unsigned int mylast;
	unsigned int node_earlier, node_later, comm_earlier, comm_later;
	if ((comm_a != comm_b) || (node_a < node_b)) {
		node_earlier = node_a;
		comm_earlier = comm_a;
		node_later = node_b;
		comm_later = comm_b;
	}
	else if (node_a > node_b) {
		node_earlier = node_b;
		comm_earlier = comm_b;
		node_later = node_a;
		comm_later = comm_a;
	}
	else {
		printf("Error, no pueden ser iguales, node_a and node_b\n");
		exit(ERR_REM_VEC_COM);
	}
	if (vector_com_links[comm_later][node_later] == 1) {
		mylast = last[comm_later];
		vector_com_nodes[comm_later][node_later] = vector_com_nodes[comm_later][mylast - 1];
		vector_com_links[comm_later][node_later] = vector_com_links[comm_later][mylast - 1];
		last[comm_later] = mylast - 1;
	}
	else
		vector_com_links[comm_later][node_later] = vector_com_links[comm_later][node_later] - 1;
	if (vector_com_links[comm_earlier][node_earlier] == 1) {
		mylast = last[comm_earlier];
		vector_com_nodes[comm_earlier][node_earlier] = vector_com_nodes[comm_earlier][mylast - 1];
		vector_com_links[comm_earlier][node_earlier] = vector_com_links[comm_earlier][mylast - 1];
		last[comm_earlier] = mylast - 1;
	}
	else
		vector_com_links[comm_earlier][node_earlier] = vector_com_links[comm_earlier][node_earlier] - 1;
}


//buscamos nodo disponible en otra comunidad para conectarlo
//los enlaces libres en las comunidades son los mas altos: total_enlaces-externos hasta total_enlaces
// para facilitar posible depuracion, (por inspecion visual saber que tipo de enlace esta libre) vamos a llenar desde total_enlaces-externos
//el ultimo enlace en llenar es total_enlaces
//Es posible que se hayan llenado todas las otras comunidades y todavia tengamos mas de un nodo que conectar
//en ese caso romperemos algun enlace entre las otras comunidades para conectar estos nodos.
//vector_nodes contiene los nodos de esta comunidad a conectar, vector_links indica cuantos enlaces tienen libres
//devuelve: la posicion del destino a conectar en dest_comm, queremos conectarlo con source-source_comm
unsigned int busca_dest_otra_comunidad(
	unsigned int* dest_comm, //extra for otra_comunidad
	const unsigned int source, //always equal to 0, because it is in the first position of vector_com_nodes
	const unsigned int misma_comunidad,
	unsigned int* last,
	unsigned int** vector_com_nodes, //los vectores van a ir haciendose mas cortos, a medida que sean asignados los enlaces
	unsigned int** vector_com_links, //los vectores van a ir haciendose mas cortos, a medida que sean asignados los enlaces
	_link*** comunidades,
	const unsigned int source_comm, //identifica esta comunidad
	unsigned int* comunity_num_nodes,
	const unsigned int total_nodes,
	//const unsigned int num_links, --> ahora es comunity_size_node_links[source_comm][vector_com_nodes[source_comm][source]]
	//const unsigned int extern_links,
	const unsigned int num_communities, //extra for otra_comunidad
	unsigned int** comunity_size_node_links) {
	const int debug = false;
	unsigned int i, enlace, dest_node;
	unsigned int community;
	if (debug) print_vector_com(vector_com_nodes, vector_com_links, last, num_communities);
	if (debug) printf("BUSCA: estamos en node %i commm %u total enlaces entre comunicades %i\n", vector_com_nodes[source_comm][source], source_comm,
		comunity_size_node_links[source_comm][vector_com_nodes[source_comm][source]]);
	unsigned int valido = false;
	//comprobamos que queda algun enlace libre en el nodo
	for (enlace = 0; enlace < comunity_size_node_links[source_comm][vector_com_nodes[source_comm][source]]; enlace++) {
		if (comunidades[source_comm][vector_com_nodes[source_comm][source]][enlace].link_node == UINT_MAX)
			valido = true;
	}
	if (valido == false) {
		printf("Error, todos los enlaces del nodo %u[%u] estan ocupados\n", vector_com_nodes[source_comm][source], source_comm);
		exit(ERR_BUSC_COMM);
	}
	unsigned int** possible_nodes = (unsigned int**)malloc(num_communities * sizeof(unsigned int*));
	for (community = 0; community < num_communities; community++) {
		possible_nodes[community] = (unsigned int*)malloc(comunity_num_nodes[community] * sizeof(unsigned int));
		for (i = 0; i < comunity_num_nodes[community]; i++)
			possible_nodes[community][i] = false;
	}

	unsigned int* local_vector = (unsigned int*)malloc(total_nodes * num_communities * sizeof(unsigned int));
	if (local_vector == NULL) {
		printf("error reservando mem\n"); fflush(stdout);
		exit(1);
	}
	unsigned int* local_vector_comm = (unsigned int*)malloc(total_nodes * num_communities * sizeof(unsigned int));
	if (local_vector_comm == NULL) {
		printf("error reservando mem\n"); fflush(stdout);
		exit(1);
	}
	//buscamos cuales de esos destinos tienen todavia un enlace libre !!

	if (debug)
		for (community = 0; community < num_communities; community++)
			printf(" last[%i]=%i ", community, last[community]);
	if (debug) printf("\n");
	if (misma_comunidad == false) {
		for (community = 0; community < num_communities; community++) {
			if (community != source_comm) {//en otra comunidad distinta del source
				for (i = 0; i < last[community]; i++) {//consideramos todos los nodos, porque estan en otra comunidad
					possible_nodes[community][vector_com_nodes[community][i]] = true;
					if (debug) {
						printf("BUSCA: source_comm %u !!possible node %u[%u] ->", source_comm, vector_com_nodes[community][i], community);
						for (int kkjj = 0; kkjj < comunity_size_node_links[community][i]; kkjj++) {
							printf("  %u[%u]", comunidades[community][i][kkjj].link_node, comunidades[community][i][kkjj].link_comunity);
						}
						printf("\n");
					}
				}
			}
		}
		possible_nodes[source_comm][vector_com_nodes[source_comm][source]] = false;// no queremos que se conecte consigo mismo
	}
	else {
		for (i = 1; i < last[source_comm]; i++)//i==0 es el source, no queremos que se conecte consigo mismo
			possible_nodes[source_comm][vector_com_nodes[source_comm][i]] = true;
		if (possible_nodes[source_comm][vector_com_nodes[source_comm][source]] == true) {
			printf("BUSCA: Error the source should be a possible destination node\n");
			exit(ERR_BUSC_SAME);
		}
	}

	//buscamos los posibles destinos, no conectados con source
	for (enlace = 0; enlace < comunity_size_node_links[source_comm][vector_com_nodes[source_comm][source]]; enlace++) {
		unsigned int n = comunidades[source_comm][vector_com_nodes[source_comm][source]][enlace].link_node;
		unsigned int c = comunidades[source_comm][vector_com_nodes[source_comm][source]][enlace].link_comunity;
		if (n != UINT_MAX) {
			possible_nodes[c][n] = false;
			//printf("no vale, ya conectado %u con source\n", k);
		}
	}

	//hacemos un vector local de posibles destinos
	unsigned int local_last = 0;
	if (misma_comunidad == false) {
		for (community = 0; community < num_communities; community++) {
			for (i = 0; i < comunity_num_nodes[community]; i++) {
				if (possible_nodes[community][i] == true) {
					local_vector[local_last] = i;
					local_vector_comm[local_last] = community;
					local_last++;
				}
			}
		}
	}
	else {
		for (i = 0; i < comunity_num_nodes[source_comm]; i++) {
			if (possible_nodes[source_comm][i] == true) {
				local_vector[local_last] = i;
				local_vector_comm[local_last] = source_comm;
				local_last++;
			}
		}
	}
	if (debug) printf("BUSCA: num_communities %u local_last %u\n", num_communities, local_last);
	if (debug)
		for (int ii = 0; ii < local_last; ii++)
			printf("local_vector[]=%u[%u]\n", local_vector[ii], local_vector_comm[ii]);
	if (local_last > 0) {//existe algun destino posible
		//elegimos uno random
		unsigned int lll = rand() % local_last;
		dest_node = local_vector[lll];
		*dest_comm = local_vector_comm[lll];
		if (debug) printf("BUSCA: dest_node%u dest_comm%u\n", dest_node, *dest_comm);
		i = 0;//entonces tenemos que devolver su posicion en el vector_nodes
		while ((vector_com_nodes[*dest_comm][i] != dest_node) && (i < last[*dest_comm]))
			i++;
		if (i == last[*dest_comm]) {
			printf("BUSCA: Error, no encontrado posible destino nodo %u, en vector_com_nodes, comunity %u\n", dest_node, *dest_comm);
			exit(ERR_BUSC_COM_NODE);
		}
		if (debug) printf("BUSCA: posicion del destino %u es %u, en comunidad %u\n", dest_node, i, *dest_comm);

		for (unsigned int enlace = 0; enlace < comunity_size_node_links[source_comm][vector_com_nodes[source_comm][source]]; enlace++) {
			//printf(" %i:-- %u\n", enlace, comunidades[source_comm][vector_nodes[source]*num_links+enlace].link_node);
			if ((comunidades[source_comm][vector_com_nodes[source_comm][source]][enlace].link_node == vector_com_nodes[*dest_comm][i])
				&& (comunidades[source_comm][vector_com_nodes[source_comm][source]][enlace].link_comunity == *dest_comm)) {
				printf("xyxy emparejamos %u[%u] con %u[%u]\n", vector_com_nodes[source_comm][source], source_comm, vector_com_nodes[*dest_comm][i], *dest_comm);
				printf("xyxy Error source ya conectado con dest, enlace %i\n", enlace);
				exit(ERR_MAIN_CONF_B);
			}
		}
		for (community = 0; community < num_communities; community++)
			free(possible_nodes[community]);
		free(possible_nodes);
		free(local_vector);
		free(local_vector_comm);
		return i;
	}
	else {
		if (debug)
			printf("BUSCA: source %u no se puede conectar con alguien!!\n", vector_com_nodes[source_comm][source]);
		// source no se puede conectar con alguien, entonces buscamos un nodo (modify_node) que no este conectado con source,
		// evidentemente, todos los posibles destinos tienen todos los enlaces ocupados, porque no encontrabamos ninguno candidato valido con enlaces libres
		// rompemos aleatoriamente uno de los enlaces, los enlaces rotos los quitamos del vector comunidades[][], y los a�adiremos al vector_nodes[] (modify_node conectado con modify_dest)
		//	pero es posible que el nodo modify_dest puede que ya estuviese conectado con source, por eso no lo habriamos considerado aunque tuviese enlaces disponibles, en ese caso todavia estara el nodo modify_dest en el vector vector_nodes[]
		if (misma_comunidad == false) {
			for (community = 0; community < num_communities; community++)
				if (community != source_comm) //en otra comunidad distinta del source
					for (i = 0; i < comunity_num_nodes[community]; i++)
						possible_nodes[community][i] = true;
		}
		else {
			for (i = 0; i < comunity_num_nodes[source_comm]; i++)
				possible_nodes[source_comm][i] = true;
		}
		possible_nodes[source_comm][vector_com_nodes[source_comm][source]] = false;//no queremos que se conecte consigo mismo

		//buscamos los posibles destinos, no conectados con source
		for (enlace = 0; enlace < comunity_size_node_links[source_comm][vector_com_nodes[source_comm][source]]; enlace++) {
			int n = comunidades[source_comm][vector_com_nodes[source_comm][source]][enlace].link_node;
			int c = comunidades[source_comm][vector_com_nodes[source_comm][source]][enlace].link_comunity;
			if (n != UINT_MAX) {
				possible_nodes[c][n] = false;
				//printf("no vale, ya conectado %u con source\n", k);
			}
		}
		//hacemos un vector local de posibles destinos QUE NO ESTE YA CONECTADO CON SOURCE !!***
		int local_last = 0;
		if (misma_comunidad == false) {
			for (community = 0; community < num_communities; community++) {
				for (i = 0; i < comunity_num_nodes[community]; i++) {
					if (possible_nodes[community][i] == true) {
						if (debug)
							printf(" posible node %i[%i] =%i[%i]", i, community, comunidades[community][i][0].link_node, comunidades[community][i][0].link_comunity);
						if (debug)
							for (int jjkk = 1; jjkk < comunity_size_node_links[community][i]; jjkk++)
								printf("                   =%i[%i]", comunidades[community][i][jjkk].link_node, comunidades[community][i][jjkk].link_comunity);
						local_vector[local_last] = i;
						local_vector_comm[local_last] = community;
						local_last++;
						if (debug) printf("\n");
					}
				}
			}
		}
		else {
			for (i = 0; i < comunity_num_nodes[source_comm]; i++) {
				if (possible_nodes[source_comm][i] == true) {
					local_vector[local_last] = i;
					local_vector_comm[local_last] = source_comm;
					local_last++;
				}
			}
		}
		if (debug) printf("z1\n");
		if (local_last > 0) {// queda algun nodo sin conectarse con source-source_comm
			if (debug) printf("z2  %i\n", local_last);
			//elegimos uno random, tiene todos los enlaces ocupados, PERO OJO QUE NO ESTE YA CONECTADO CON SOURCE !!***
			unsigned int modify_poss = rand() % local_last;
			unsigned int modify_node = local_vector[modify_poss];
			unsigned int modify_comm = local_vector_comm[modify_poss];

			if (debug) printf("z7 %u[%u] \n", modify_node, modify_comm);
			if (debug) printf("z7 %u \n", comunity_size_node_links[modify_comm][modify_node]);


			//elegimos un enlace random
			unsigned int modify_node_link = rand() % comunity_size_node_links[modify_comm][modify_node]; // num_links;
			//unsigned int modify_node_link= rand() %
			if (debug) printf("z3 33\n");
			unsigned int modify_dest = comunidades[modify_comm][modify_node][modify_node_link].link_node;
			if (debug) printf("z3 44\n");
			unsigned int modify_dest_comm = comunidades[modify_comm][modify_node][modify_node_link].link_comunity;
			if (debug) printf("z3  55\n");
			unsigned int modify_dest_link = UINT_MAX;
			//printf("buscamos en la comunidad %u el enlace de dest=%u que lleva a %u\n", source_comm, modify_dest, modify_node);
			//if(debug==true) print_graph ( num_communities, comunidades, comunity_num_nodes, omunity_size_node_links);
			if (modify_dest == UINT_MAX) {
				printf("ERROR modify_dest=UINT_MAX\n");
				exit(ERR_BUSC_COM_MOD_A);
			}
			if (debug) printf("z3\n");
			unsigned int enlace = 0;
			do {
				if ((modify_node == comunidades[modify_dest_comm][modify_dest][enlace].link_node)
					&& (modify_comm == comunidades[modify_dest_comm][modify_dest][enlace].link_comunity)) {
					modify_dest_link = enlace;
				}
				enlace++;
			} while ((enlace < comunity_size_node_links[modify_dest_comm][vector_com_nodes[modify_dest_comm][modify_dest]]) && (modify_dest_link == UINT_MAX));
			if (modify_dest_link == UINT_MAX) {
				printf("BUSCA: ERROR, no encontrado enlace de entrada en nodo %u desde nodo %u\n", modify_dest, modify_node);
				exit(ERR_BUSC_COM_MOD_B);
			}
			if (debug) printf("z4\n");
			// queremos romper la conexion modify_node-modify_node_link con modify_dest-modify_dest_link
			if (debug) printf("queremos romper la conexion node:%u[%u]-link:%u con node:%u[%u]-link:%u\n",
				modify_node, modify_comm, modify_node_link,
				modify_dest, modify_dest_comm, modify_dest_link);
			// para borrar la conexion source-dest a la vector de comunidades:
			// se desplazan las conexiones existentes, dejando sin definir la ULTIMA de la conexiones.
			// primero para el nodo origen
			remove_connection_comm(modify_node, modify_node_link, modify_comm, comunidades, comunity_size_node_links);
			// segundo para el nodo destino
			remove_connection_comm(modify_dest, modify_dest_link, modify_dest_comm, comunidades, comunity_size_node_links);
			//ahora añadimos al vector estos dos nuevos enlaces modify_node y modify_dest
			// el nodo modify_dest puede que si este en vector_com_nodes[]
			unsigned int j = 0;
			unsigned int found = false;
			while ((j < last[modify_dest_comm]) && (found == false)) {
				if (vector_com_nodes[modify_dest_comm][j] == modify_dest) {
					found = true;
				}
				else
					j++;
			}
			if (found == false) {
				vector_com_nodes[modify_dest_comm][last[modify_dest_comm]] = modify_dest;
				vector_com_links[modify_dest_comm][last[modify_dest_comm]] = 1;
				last[modify_dest_comm] = last[modify_dest_comm] + 1;
			}
			else
				vector_com_links[modify_dest_comm][j] = vector_com_links[modify_dest_comm][j] + 1;
			//el nodo modify_node no esta en vector_com_nodes[]
			vector_com_nodes[modify_comm][last[modify_comm]] = modify_node;
			vector_com_links[modify_comm][last[modify_comm]] = 1;
			last[modify_comm] = last[modify_comm] + 1;
			*dest_comm = modify_comm;
			i = last[modify_comm] - 1;
			for (unsigned int enlace = 0; enlace < comunity_size_node_links[source_comm][vector_com_nodes[source_comm][source]]; enlace++) {
				//printf(" %i:-- %u\n", enlace, comunidades[source_comm][vector_nodes[source]*num_links+enlace].link_node);
				if ((comunidades[source_comm][vector_com_nodes[source_comm][source]][enlace].link_node == vector_com_nodes[*dest_comm][i])
					&& (comunidades[source_comm][vector_com_nodes[source_comm][source]][enlace].link_comunity == *dest_comm)) {
					printf("xyxy emparejamos %u[%u] con %u[%u]\n", vector_com_nodes[source_comm][source], source_comm, vector_com_nodes[*dest_comm][i], *dest_comm);
					printf("xyxy ERROR source ya conectado con dest, enlace %i\n", enlace);
					exit(ERR_MAIN_CONF_B);
				}
			}
			for (community = 0; community < num_communities; community++)
				free(possible_nodes[community]);
			free(possible_nodes);
			free(local_vector);
			free(local_vector_comm);
			return (last[modify_comm] - 1);
		}
		else {
			for (community = 0; community < num_communities; community++)
				free(possible_nodes[community]);
			free(possible_nodes);
			free(local_vector);
			free(local_vector_comm);
			return (UINT_MAX);
		}
	}
}//busca_dest_otra_comunidad


 // comunity_size_node_links_otra contains mu_extern_links (stands for mu from 0.1 ..0.2 .. 1.0.   0.1 means add a 10% of links) links of a community with the other comunities
int interconnecta_comunidades(const unsigned int num_communities,
	_link*** comunidades_otra, unsigned int** comunity_size_node_links_otra,
	unsigned int* comunity_num_nodes, const unsigned int total_nodes)
{
	unsigned int source_comm, dest_comm;
	unsigned int dest_pos = UINT_MAX;
	unsigned int total_last = 0;
	unsigned int total_links = 0;
	unsigned int misma_comunidad = false;
	const unsigned int debug = false;
	printf("\nConectando enlaces entre comunidades .......\n");
	//los vectores van a ir haciendose mas cortos, a medida que sean asignados los enlaces
	unsigned int** vector_com_nodes = (unsigned int**)malloc(num_communities * sizeof(unsigned int*));
	unsigned int** vector_com_links = (unsigned int**)malloc(num_communities * sizeof(unsigned int*));
	unsigned int* vector_com_last = (unsigned int*)malloc(num_communities * sizeof(unsigned int));
	for (int comm = 0; comm < num_communities; comm++) {
		vector_com_nodes[comm] = (unsigned int*)malloc(comunity_num_nodes[comm] * sizeof(unsigned int));
		vector_com_links[comm] = (unsigned int*)malloc(comunity_num_nodes[comm] * sizeof(unsigned int));
	}

	//ahora conectamos las comunidades entre ellas.
	//vamos a conectar el enlace libre de cada nodo con otro nodo de otra comunidad
	// if(misma_comunidad==false){
	unsigned int contador;
	for (unsigned int comm = 0; comm < num_communities; comm++) {	//this is the community
		contador = 0;
		for (unsigned int nodo = 0; nodo != comunity_num_nodes[comm]; nodo++) {
			vector_com_nodes[comm][contador] = nodo;
			vector_com_links[comm][contador] = comunity_size_node_links_otra[comm][nodo];//links que nos guardamos para conectar cada nodo con otra comunidad
			if (comunity_size_node_links_otra[comm][nodo] > 0)
				contador++;
		}
		vector_com_last[comm] = contador;//comunity_num_nodes[comm];
	}
	// }else{
	// for(int comm=0;comm<num_communities;comm++)
	// vector_com_last[comm]=0;
	// for (int nodo = 0;nodo != comunity_num_nodes[source_comm];nodo++) {
	// vector_com_nodes[source_comm][nodo] = nodo;
	// vector_com_links[source_comm][nodo] = num_links;//nos guardamos un link para conectar cada nodo con otra comunidad
	// }
	// vector_com_last[source_comm]=comunity_num_nodes[source_comm];
	// }

	//for(int comm=0;comm<num_communities;comm++)
	//	for (int nodo = 0; nodo < comunity_num_nodes[comm]; nodo++)
	//		printf(" DEFF: %u[%u]-> links %u\n", nodo, comm, comunity_size_node_links_otra[comm][nodo]);

	unsigned int source_pos = 0;
	source_comm = UINT_MAX;
	for (unsigned int comm = 0; comm < num_communities; comm++) { //this is the community
		total_last += vector_com_last[comm];
		for (unsigned int pos = 0; pos < vector_com_last[comm]; pos++) {
			total_links += vector_com_links[comm][pos];
			printf(" links of node %u[%u] are %u\n", vector_com_nodes[comm][pos], comm, vector_com_links[comm][pos]);
		}
		if ((vector_com_last[comm] != 0) && (source_comm == UINT_MAX))
			source_comm = comm;
	}

	//	while(comunity_size_node_links[source_comm][vector_com_nodes[source_comm][source_pos]]==0){
	//		source_pos++;
	//		if(source_pos>=comunity_num_nodes[source_comm]){
	//			source_comm++;
	//			source_pos=0;
	//		}
	//	}
	if (source_comm >= num_communities) {
		printf("Error,  inicio interconnecta_comunidades\n");
	}
	misma_comunidad = false;
	while ((total_last > 1) || (total_links > 1)) {//at least "last" is 
		if (source_comm >= num_communities || source_pos >= vector_com_last[source_comm]) {
			printf("Error \n");
			exit(1);
		}
		if (debug) printf("Empezamos vector_com_last[comm] %i total_last %i total_links %i\n", vector_com_last[source_comm], total_last, total_links);
		//const source_pos=0;procesamos el primer nodo del vector
		dest_pos = busca_dest_otra_comunidad(&dest_comm, source_pos, misma_comunidad, vector_com_last, vector_com_nodes, vector_com_links,
			comunidades_otra, source_comm, comunity_num_nodes, total_nodes, num_communities, comunity_size_node_links_otra);
		if (debug) printf("source_comm %i num_nodes[comm]%i dest_comm%i   source_pos %i dest_pos %i  extern_links%i\n",
			source_comm, comunity_num_nodes[source_comm], dest_comm, source_pos, dest_pos,
			comunity_size_node_links_otra[source_comm][0]);
		if (dest_pos == UINT_MAX) {//not possible to find dest_pos
			//tenemos que quitar el source_pos source_comm de la lista vector_com_nodes[source_comm][source]
			unsigned int lastpos = vector_com_last[source_comm];
			vector_com_nodes[source_comm][source_pos] = vector_com_nodes[source_comm][lastpos - 1];
			vector_com_links[source_comm][source_pos] = vector_com_links[source_comm][lastpos - 1];
			vector_com_last[source_comm] = vector_com_last[source_comm] - 1;
		}
		else {
			if (vector_com_links[dest_comm][dest_pos] > comunity_size_node_links_otra[dest_comm][dest_pos]) {
				printf("error link demasiado grande\n");
				printf(" link %i > max=%i \n", vector_com_links[dest_comm][dest_pos], comunity_size_node_links_otra[dest_comm][dest_pos]);
				exit(1);
			}
			if (debug) printf("source %u[%u] dest %u[%u]\n", vector_com_nodes[source_comm][source_pos],
				source_comm, vector_com_nodes[dest_comm][dest_pos], dest_comm);
			for (unsigned int enlace = 0; enlace < comunity_size_node_links_otra[source_comm][vector_com_nodes[source_comm][source_pos]]; enlace++) {
				//printf(" %i:-- %u\n", enlace, comunidades[source_comm][vector_nodes[source_pos]*num_links+enlace].link_node);
				if ((comunidades_otra[source_comm][vector_com_nodes[source_comm][source_pos]][enlace].link_node == vector_com_nodes[dest_comm][dest_pos])
					&& (comunidades_otra[source_comm][vector_com_nodes[source_comm][source_pos]][enlace].link_comunity == dest_comm)) {
					printf("emparejamos %u[%u] con %u[%u]\n", vector_com_nodes[source_comm][source_pos], source_comm, vector_com_nodes[dest_comm][dest_pos], dest_comm);
					printf("Error source ya conectado con dest, enlace %i\n", enlace);
					exit(ERR_MAIN_CONF_B);
				}
			}
			if ((vector_com_nodes[source_comm][source_pos] != vector_com_nodes[dest_comm][dest_pos]) || (source_comm != dest_comm)) {
				if (debug)
					printf("emparejamos %u con %u\n", vector_com_nodes[source_comm][source_pos], vector_com_nodes[dest_comm][dest_pos]);
				//añadimos la conexion source-dest a la vector de comunidades
				comunidades_otra[source_comm][vector_com_nodes[source_comm][source_pos]][vector_com_links[source_comm][source_pos] - 1].link_node
					= vector_com_nodes[dest_comm][dest_pos];
				comunidades_otra[source_comm][vector_com_nodes[source_comm][source_pos]][vector_com_links[source_comm][source_pos] - 1].link_comunity
					= dest_comm;
				comunidades_otra[source_comm][vector_com_nodes[source_comm][source_pos]][vector_com_links[source_comm][source_pos] - 1].link_weight_out = 1;

				comunidades_otra[dest_comm][vector_com_nodes[dest_comm][dest_pos]][vector_com_links[dest_comm][dest_pos] - 1].link_node
					= vector_com_nodes[source_comm][source_pos];
				comunidades_otra[dest_comm][vector_com_nodes[dest_comm][dest_pos]][vector_com_links[dest_comm][dest_pos] - 1].link_comunity
					= source_comm;
				comunidades_otra[dest_comm][vector_com_nodes[dest_comm][dest_pos]][vector_com_links[dest_comm][dest_pos] - 1].link_weight_out = 1;
				//quitamos los valores source y dest del vector temporal
				if (debug)
					printf("quitamos de vector la posicion %u[%u] y posicion %u[%u]\n", source_pos, source_comm, dest_pos, dest_comm);
				if (debug)  printf("quedan vector_com_last[comm] %i vector_com_last[comm] %i\n", vector_com_last[source_comm], vector_com_last[dest_comm]);
				//if(debug) print_vector_com( vector_com_nodes, vector_com_links, vector_com_last, num_communities);
				remove_from_vector_comm(vector_com_nodes, vector_com_links, vector_com_last, source_pos, source_comm, dest_pos, dest_comm);
				if (debug)  printf("quedan vector_com_last[comm] %i vector_com_last[comm] %i\n", vector_com_last[source_comm], vector_com_last[dest_comm]);

				// if(debug) print_graph_double ( num_communities, NULL, comunidades_otra, comunity_num_nodes, NULL, comunity_size_node_links_otra );
				// print_graph ( num_communities, comunidades_otra, comunity_num_nodes, comunity_size_node_links_otra);
				//if(debug) print_vector( vector_nodes, vector_links, last);
			}
		}
		//*****************************************************************
		//Antencion	fix !!!
		//*****************************************************************
		//si queremos conectar extern_links de cada nodo dentro de una comunidad, nos hacen falta una posible cantidad de destinos: extern_links*num_nodos
		// max_destinos= sum(comunity_nodes[different-comm]
		//pero si la comunidad es muy grande puede no existir tantos destinos por lo que el bucle nunca terminaria !!
		//lo limitamos a max_destinos/2 en el caso que necesitmos romper conexiones entre otros enlaces externos
		//contamos max_posibles_dest:
		unsigned int max_posibles_dest = 0;
		for (int temp_comm = 0; temp_comm < num_communities; temp_comm++) { //this is the community
			if (temp_comm != source_comm)
				max_posibles_dest += comunity_num_nodes[temp_comm];
		}
		//contamos cuantos quedan sin definir:
		unsigned int total_links_sin_definir = 0;
		for (unsigned int j = 0; j < num_communities; j++) {	//this is the community
			if (j != source_comm) {
				for (unsigned int pos = 0; pos < vector_com_last[j]; pos++)
					total_links_sin_definir += vector_com_links[j][pos];
			}
		}
		//total de enlaces en la comunidad source_comm que queriamos usar:
		unsigned int total_links_queriamos_usar = 0;
		for (int ni = 0; ni < comunity_num_nodes[source_comm]; ni++)
			total_links_queriamos_usar += comunity_size_node_links_otra[source_comm][ni];
		//contamos cuantos faltan por usar en source_comm:
		unsigned int links_por_usar_source_comm = 0;
		for (unsigned int j = 0; j < num_communities; j++) {	//this is the community
			if (j == source_comm) {
				for (unsigned int pos = 0; pos < vector_com_last[j]; pos++)
					links_por_usar_source_comm += vector_com_links[j][pos];
			}
		}
		unsigned int total_usados = 0;
		if (total_links_queriamos_usar > links_por_usar_source_comm)
			total_usados = total_links_queriamos_usar - links_por_usar_source_comm;
		//si no quedan sin definir, y hemos usado mas de max_posibles_dest/2 entonces hemos terminado con esta communidad source_comm
		//quitamos los que quedan de esta comunidad source_comm en vector_com_links
		//con la tranquilidad de que si rompemos algun enlace de source_comm, se a�adiria nuevamente a la lista
		if ((total_usados > max_posibles_dest / 2) && (total_links_sin_definir == 0)) {
			vector_com_last[source_comm] = 0;
		}
		//*****************************************************************
		//end fix
		//*****************************************************************
		total_last = 0;
		total_links = 0;
		source_comm = UINT_MAX;
		for (unsigned int j = 0; j < num_communities; j++) {	//this is the community
			total_last += vector_com_last[j];
			for (unsigned int pos = 0; pos < vector_com_last[j]; pos++)
				total_links += vector_com_links[j][pos];
			if (vector_com_last[j] != 0) {
				if (source_comm == UINT_MAX)
					source_comm = j;
			}
		}
	}
	for (int comm = 0; comm < num_communities; comm++) {
		free(vector_com_nodes[comm]);
		free(vector_com_links[comm]);
	}
	free(vector_com_nodes);
	free(vector_com_links);
	free(vector_com_last);
	printf(" fin conectando\n");
	if (dest_pos != UINT_MAX) {
		return 0;
	}
	else {
		return 1;// stands for not possible to find dest_pos
	}
}


int interconnecta_comunidades_new(const unsigned int num_communities, //constante total de comunidades
	_link*** comunidades_orig, unsigned int** comunity_size_node_links_orig,
	_link*** comunidades_otra, unsigned int** comunity_size_node_links_otra, //arrays para los resultados, inicialmente vacios
	unsigned int* comunity_num_nodes, const unsigned int total_nodes, //constantes de numero de nodos
	double* comm_external_links, unsigned int** total_used_links_at_node, unsigned int mu_extern_links) { //comm_external_links es el numero de enlaces de salida de cada comunidad
	//const int debug=true;
	//posibles destinos para cada node-comm
	unsigned int*** posible_dst_nodes = (unsigned int***)malloc(num_communities * sizeof(unsigned int**));
	if (posible_dst_nodes == NULL) {
		printf("error reservando mem posible_dst_nodes\n"); fflush(stdout);
		exit(1);
	}
	unsigned int*** posible_dst_comm = (unsigned int***)malloc(num_communities * sizeof(unsigned int**));
	if (posible_dst_comm == NULL) {
		printf("error reservando mem\n"); fflush(stdout);
		exit(1);
	}
	unsigned int** posible_dst_total = (unsigned int**)malloc(num_communities * sizeof(unsigned int*));
	if (posible_dst_total == NULL) {
		printf("error reservando mem posible_dst_total\n"); fflush(stdout);
		exit(1);
	}
	for (int comm_src = 0; comm_src < num_communities; comm_src++) {
		posible_dst_nodes[comm_src] = (unsigned int**)malloc(comunity_num_nodes[comm_src] * sizeof(unsigned int*));
		if (posible_dst_nodes[comm_src] == NULL) {
			printf("error reservando mem\n"); fflush(stdout);
			exit(1);
		}
		posible_dst_comm[comm_src] = (unsigned int**)malloc(comunity_num_nodes[comm_src] * sizeof(unsigned int*));
		if (posible_dst_comm[comm_src] == NULL) {
			printf("error reservando mem\n"); fflush(stdout);
			exit(1);
		}
		posible_dst_total[comm_src] = (unsigned int*)malloc(comunity_num_nodes[comm_src] * sizeof(unsigned int));
		for (int node = 0; node < comunity_num_nodes[comm_src]; node++) {
			posible_dst_nodes[comm_src][node] = (unsigned int*)malloc((total_nodes - comunity_num_nodes[comm_src]) * sizeof(unsigned int));
			if (posible_dst_nodes[comm_src][node] == NULL) {
				printf("error reservando mem\n"); fflush(stdout);
				exit(1);
			}
			posible_dst_comm[comm_src][node] = (unsigned int*)malloc((total_nodes - comunity_num_nodes[comm_src]) * sizeof(unsigned int));
			if (posible_dst_comm[comm_src][node] == NULL) {
				printf("error reservando mem\n"); fflush(stdout);
				exit(1);
			}
		}
		for (int node_src = 0; node_src < comunity_num_nodes[comm_src]; node_src++) {
			//for(int i=0; i<total_nodes-comunity_num_nodes[comm_src]; i++){
			//	posible_dst_nodes[comm_src][node_src][i]=0;
			//	posible_dst_comm[comm_src][node_src][i]=0;
			//	posible_src_nodes[comm_src][node_src][i]=0;
			//	posible_src_comm[comm_src][node_src][i]=0;
			//}
			posible_dst_total[comm_src][node_src] = total_nodes - comunity_num_nodes[comm_src];
			unsigned int contador = 0;
			for (int comm_dst = 0; comm_dst < num_communities; comm_dst++) {
				if (comm_src != comm_dst) {
					for (int node_dest = 0; node_dest < comunity_num_nodes[comm_dst]; node_dest++) {
						posible_dst_nodes[comm_src][node_src][contador] = node_dest;
						posible_dst_comm[comm_src][node_src][contador] = comm_dst;
						contador++;
					}
				}
			}//contador must be equal to posible_dst_total[comm_src][node_src]
		}
	}

	for (unsigned int comm_src = 0; comm_src < num_communities; comm_src++) {
		unsigned int* lista_posibles_src = (unsigned int*)malloc(comunity_num_nodes[comm_src] * sizeof(unsigned int));
		for (unsigned int i = 0; i < comunity_num_nodes[comm_src]; i++)
			lista_posibles_src[i] = i;
		int total_lista_posibles_src = comunity_num_nodes[comm_src];
		//if(debug)
		printf("peticion de enlaces de %3.1f enlaces de salida para la comunidad %2i\n", comm_external_links[comm_src], comm_src);
		int comm_external_links_this = comm_external_links[comm_src];

		//while(( comm_external_links_this>0) &&(node_src< comunity_num_nodes[comm_src])){
		while ((comm_external_links_this > 0) && (total_lista_posibles_src > 0)) {
			int pos_node_src = rand() % total_lista_posibles_src;
			int node_src = lista_posibles_src[pos_node_src];
			unsigned int total = posible_dst_total[comm_src][node_src];
			if (total == 0) {
				printf("Error node %i comm %i total nodes %i\n", node_src, comm_src, comunity_num_nodes[comm_src]);
				exit(1);
			}
			int cont_rand = rand() % posible_dst_total[comm_src][node_src];
			int node_dst = posible_dst_nodes[comm_src][node_src][cont_rand];
			int comm_dst = posible_dst_comm[comm_src][node_src][cont_rand];
			if (node_dst == UINT_MAX) {
				printf("node_dst UINT_MAX, EOP\n");
				fflush(stdout);
				exit(1);
			}
			//comunity_size_node_links_otra[comm_src][node_src]; total de memoria reservada para los enlaces.
			int enlace = total_used_links_at_node[comm_src][node_src];
			if (enlace >= comunity_size_node_links_otra[comm_src][node_src]) { //tenemos que reservar mas memoria
				unsigned int total_links = comunity_size_node_links_otra[comm_src][node_src];
				_link* previo = comunidades_otra[comm_src][node_src];
				comunidades_otra[comm_src][node_src] = (_link*)malloc(total_links + 50 * sizeof(_link));
				for (int temp_link = 0; temp_link < total_links; temp_link++) {
					comunidades_otra[comm_src][node_src][temp_link].link_node = previo[temp_link].link_node;
					comunidades_otra[comm_src][node_src][temp_link].link_comunity = previo[temp_link].link_comunity;
					comunidades_otra[comm_src][node_src][temp_link].link_weight_out = previo[temp_link].link_weight_out;
				}
				for (int temp_link = total_links; temp_link < total_links + 50; temp_link++) {
					comunidades_otra[comm_src][node_src][temp_link].link_node = UINT_MAX;
					comunidades_otra[comm_src][node_src][temp_link].link_comunity = UINT_MAX;
				}
				free(previo);
				comunity_size_node_links_otra[comm_src][node_src] = total_links + 50;
			}
			comunidades_otra[comm_src][node_src][enlace].link_node = node_dst;
			comunidades_otra[comm_src][node_src][enlace].link_comunity = comm_dst;//esto es la comunidad del nodo destino
			comunidades_otra[comm_src][node_src][enlace].link_weight_out = 1;
			// en lugar de peso 1, vamos a poner el 90% del enlace mas debil
			// de entre los enlaces que llegan o salen del node_src o node_dst
			unsigned int total_linksxx = comunity_size_node_links_orig[comm_src][node_src];
			//if(comm_src ==0 )
			//printf("***** new link node_src %d node_dst %d com_dst %i\n", node_src, node_dst, node_dst);
			float www = -1;
			// unsigned int node_dst_www =0;
			for (int temp_link = 0; temp_link < total_linksxx; temp_link++) {
				if (comunidades_orig[comm_src][node_src][temp_link].link_node != UINT_MAX) {
					float wgwg = comunidades_orig[comm_src][node_src][temp_link].link_weight_out;
					// 	if(comm_src==0 || comm_src==7)
					// 	printf(" www com_src %d node_src %d node_dst %d   peso %f \n",comm_src,node_src,
					// comunidades_orig[comm_src][node_src][temp_link].link_node, wgwg);
					if ((wgwg > 0 && wgwg < www) || www == -1) {
						www = wgwg;
						// node_dst_www= comunidades_orig[comm_src][node_src][temp_link].link_node;
						// printf(" mayor peso comm_src%d node_src %d comm_dst %d node_dst %d temp_link %d peso %f\n", comm_src, node_src, comm_dst, comunidades_orig[comm_src][node_src][temp_link].link_node, temp_link, www);
					}
				}
			}

			if (www > 0) {
				comunidades_otra[comm_src][node_src][enlace].link_weight_out = www * .9;
				//if (comm_src == 0 )
					//printf("***** new link comm_src%d node_src %d node_dst %d com_dst %i peso %f\n", comm_src,node_src,node_dst, node_dst, comunidades_otra[comm_src][node_src][enlace].link_weight_out);
				// unsigned int node_dst_www= comunidades_orig[comm_src][node_src][enlace].link_node;
				// printf("---xxx comm_src%d node_src %d node_dst %d linknum %d peso %f\n", comm_src, node_src, node_dst_www, enlace, www);
			}

			total_used_links_at_node[comm_src][node_src] = enlace + 1;
			//printf(" new ext link %u[%u] -> %u[%u]\n", node_src, comm_src, node_dst, comm_dst);
			total--;
			posible_dst_nodes[comm_src][node_src][cont_rand] = posible_dst_nodes[comm_src][node_src][total];
			posible_dst_comm[comm_src][node_src][cont_rand] = posible_dst_comm[comm_src][node_src][total];
			posible_dst_total[comm_src][node_src] = total;
			comm_external_links_this--;
			//printf(" %u[%u] elegimos node_dest %u comm_%u quedan en esa cmm %u pendientes %u\n", node_src,comm_src, node_dst, comm_dst ,total, comm_external_links_this);
			if (total == 0) {
				lista_posibles_src[pos_node_src] = lista_posibles_src[total_lista_posibles_src - 1];
				total_lista_posibles_src--;
				//int dst_total=posible_dst_total[comm_src][node_src];
				//while((dst_total==0)&&(node_src<comunity_num_nodes[comm_src])){
				//	node_src++;
				//	if( node_src < comunity_num_nodes[comm_src] )
				//		dst_total=posible_dst_total[comm_src][node_src];
				//}
			}
		}
		//if((node_src>= comunity_num_nodes[comm_src]) &&(comm_external_links_this>0))
		if ((total_lista_posibles_src == 0) && (comm_external_links_this > 0))
			printf("Warning, se acabaron los destinos de la comm %u y quedaban por conectar %u de un total de %.1f:  mu %i%%\n",
				comm_src, comm_external_links_this, comm_external_links[comm_src], mu_extern_links * 10);
		//ENLACES DE SALIDA
		//para cada comunidad comm asignamos 1/2 de los enlaces extra
		//tomamos un nodo aleatorio y lo conectamos con alguno de su lista node_dst=posible_dst_nodes[comm][node][rand]
		//then posible_dst_total[comm][node]--
		//comunidades_otra[comm][nodo][enlace].link_node       = node_dst
		//comunidades_otra[comm][node][enlace].link_comunity   = comm_dst;//esto es la comunidad del nodo destino
		//comunidades_otra[comm][node][enlace].link_weight_out = (2<<(long_estud_preferen-0))/16.0;
		free(lista_posibles_src);
	}
	//ENLACES DE ENTRADA
	//para cada comunidad comm asignamos 1/2 de los enlaces extra
	//tomamos un nodo aleatorio y lo conectamos con alguno de su lista nodo_src=posible_src_nodes[comm][node][rand]
	//then posible_src_total[comm][node]--
	//comunidades_otra[comm_src][nodo_src][enlace].link_node=node
	//comunidades_otra[comm_src][nodo_src][enlace].link_comunity =comm;//esto es la comunidad del nodo destino
	//comunidades_otra[comm_src][nodo_src][enlace].link_weight_out = (2<<(long_estud_preferen-0))/16.0;

	for (unsigned int comm = 0; comm < num_communities; comm++) {
		for (unsigned int node = 0; node < comunity_num_nodes[comm]; node++) {
			free(posible_dst_nodes[comm][node]);
			free(posible_dst_comm[comm][node]);
		}
		free(posible_dst_nodes[comm]);
		free(posible_dst_comm[comm]);
		free(posible_dst_total[comm]);
	}
	free(posible_dst_nodes);
	free(posible_dst_comm);
	free(posible_dst_total);
	return 0;
}