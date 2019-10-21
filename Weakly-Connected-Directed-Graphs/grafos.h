#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
//typedef enum {FALSE = 0, TRUE} boolean;
#define ERR_COMUNIDAD 		100
#define ERR_REM_VECTOR 		101
#define ERR_REM_VEC_COM 	102
#define ERR_REM_CONNECTION 	103
#define ERR_REM_CONN_COM 	104
#define ERR_BUSC_COMM 		205
#define ERR_BUSC_COM_NODE 	206
#define ERR_BUSC_COM_MOD_A 	207
#define ERR_BUSC_COM_MOD_B 	208
#define ERR_BUSC_SAME 		305
#define ERR_BUSC_SAM_NODE 	306
#define ERR_BUSC_SAM_MOD_A 	307
#define ERR_BUSC_SAM_MOD_B 	308
#define ERR_MAIN_CONF_A  	401
#define ERR_MAIN_CONF_B  	402




typedef struct{
		unsigned int link_node;		//id of the destination node
		unsigned int link_comunity;	// community of the destination node 
		double link_weight_out;		// community of the destination node 
} _link;


void normaliza_pesos_salida(const unsigned int num_communities, _link*** comunidades, _link*** comunidades_otra, unsigned int *comunity_num_nodes, 
unsigned int** comunity_size_node_links, unsigned int** comunity_size_node_links_otra);



//this function prints in screen the nodes in each community 	 
void print_graph (const unsigned int num_comunities, _link*** comunidades, unsigned int *comunity_num_nodes, 
		unsigned int** comunity_size_node_links);
	
void print_graph_double (const unsigned int num_comunities, _link*** comunidades, _link*** comunidades_otra, unsigned int *comunity_num_nodes, 
		unsigned int** comunity_size_node_links,
		unsigned int** comunity_size_node_links_otra);

void verify_comunities(const unsigned int num_comunities, _link*** comunidades, unsigned int *comunity_num_nodes, unsigned int** comunity_size_node_links) ;
	
void remove_connection_comm(const unsigned int nodo, const unsigned int link, const unsigned int comm, 
		_link*** comunidades,
		unsigned int** comunity_size_node_links);

void remove_from_vector_comm(unsigned int** vector_com_nodes,unsigned int** vector_com_links, unsigned int *last, 
		const unsigned int node_a, const unsigned int comm_a, 
		const unsigned int node_b, const unsigned int comm_b);

unsigned int busca_dest_otra_comunidad(
			unsigned int *dest_comm, //extra for otra_comunidad
			const unsigned int source, //always equal to 0, because it is in the first position of vector_com_nodes
			const unsigned int misma_comunidad, 
			unsigned int *last, 
			unsigned int** vector_com_nodes, //los vectores van a ir haciendose mas cortos, a medida que sean asignados los enlaces
			unsigned int** vector_com_links, //los vectores van a ir haciendose mas cortos, a medida que sean asignados los enlaces
			_link*** comunidades, 
			const unsigned int source_comm, //identifica esta comunidad
			unsigned int *comunity_num_nodes, 
			const unsigned int total_nodes,
			//const unsigned int num_links, 
			//const unsigned int extern_links, 
			const unsigned int num_comunities, //extra for otra_comunidad
			unsigned int **comunity_size_node_links );
			 
int interconnecta_comunidades(  const unsigned int num_comunities, 
_link*** comunidades_otra, unsigned int **comunity_size_node_links_otra, unsigned int *comunity_num_nodes, const unsigned int total_nodes ) ;

int interconnecta_comunidades_new( const unsigned int num_communities, //constante total de comunidades
	_link*** comunidades_otra, unsigned int **comunity_size_node_links_otra, //arrays para los resultados, inicialmente vacios
	unsigned int *comunity_num_nodes, const unsigned int total_nodes, //constantes de numero de nodos
	float *comm_external_links, unsigned int **total_used_links_at_node, unsigned int mu_extern_links);

