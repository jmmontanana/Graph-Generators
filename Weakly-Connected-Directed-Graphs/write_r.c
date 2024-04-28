#include <stdio.h>
#include <errno.h>
#include "grafos.h"

#if defined _MSC_VER
#include <string>
#include <direct.h>
#include <stdio.h>
#include <io.h> // Para _access
#include <stdlib.h> // Para _access y _CRT_ERROR
#include <iostream>
#include <fstream>
#include <vector>
#define access _access
/* Values for the second argument to access.
These may be OR'd together.  */
#define R_OK    4       /* Test for read permission.  */
#define W_OK    2       /* Test for write permission.  */
//#define   X_OK    1       /* execute permission - unsupported in windows*/
#define F_OK    0       /* Test for existence.  */
#elif defined __GNUC__
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

void makedir(const char path[]) { 
	struct stat st = { 0 };
	if (stat(path, &st) == -1) {
#if defined _MSC_VER
		if (_mkdir(path) != 0) {
			perror("_mkdir failed");
			exit(EXIT_FAILURE);
		}
#elif defined __GNUC__
		if (mkdir(path, 0700) != 0) {
			perror("mkdir failed");
			exit(EXIT_FAILURE);
		}
#endif
	}
}

char* mitoa(int value, char* result, int base) {
	// check that the base if valid
	if (base < 2 || base > 36) { *result = '\0'; return result; }
	char* ptr = result, * ptr1 = result, tmp_char;
	int tmp_value;
	do {
		tmp_value = value;
		value /= base;
		*ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz"[35 + (tmp_value - value * base)];
	} while (value);
	// Apply negative sign
	if (tmp_value < 0) *ptr++ = '-';
	*ptr-- = '\0';
	while (ptr1 < ptr) {
		tmp_char = *ptr;
		*ptr-- = *ptr1;
		*ptr1++ = tmp_char;
	}
	return result;
}

   //return 1 if error
int write_r_file(const unsigned int num_comunities,
	unsigned int* comunity_num_nodes,
	char* comunity_file,//for instance "community.dat"
	char* network_file, //for instance "network.dat"
	char* nodes_by_comm_file, //for instance "nodes_by_comm.dat"
	unsigned int** comunity_size_node_links, _link*** comunidades,
	unsigned int** comunity_size_node_links_otra, _link*** comunidades_otra)
{
	unsigned int* acumulado_id = (unsigned int*)malloc(num_comunities * sizeof(unsigned int));
	unsigned int comm, nodo, dest, enlace, dest_comm;
	FILE* fpcom, * fpnet, * fpnod;
	// Verifica si alguno de los archivos existe
#if defined _MSC_VER
 	if ((_access(nodes_by_comm_file, F_OK) != -1) || (_access(comunity_file, F_OK) != -1)) {
		printf(" Not generating community file because it already exists\n");
	}
#elif defined __GNUC__
	if ((access(nodes_by_comm_file, F_OK) != -1) || (access(comunity_file, F_OK) != -1)) { //other values that can be tested are R_OK W_OK X_OK, from unistd.h
		printf(" Not generating community file because it already exists\n");
	}
#endif	
	else {//file not exists
#if defined _MSC_VER
		fopen_s(&fpnod,nodes_by_comm_file, "w");
		fopen_s(&fpcom,comunity_file, "w");
#elif defined __GNUC__
		fpnod = fopen(nodes_by_comm_file, "w");
		fpcom = fopen(comunity_file, "w");
#endif
		if (fpnod != NULL && fpcom != NULL) {
			//los nodos dentro de cada comunidad los hemos numerado desde 0.
			//ahora queremos que el id de un nodo de la comunidad "comm+1" sea el id+1 del ultimo nodo de la comunidad "comm"
			for (comm = 0; comm < num_comunities; comm++)
				acumulado_id[comm] = (comm == 0) ? 0 : acumulado_id[comm - 1] + comunity_num_nodes[comm - 1];
			for (comm = 0; comm < num_comunities; comm++) {
				for (nodo = 0; nodo < comunity_num_nodes[comm]; nodo++) {
					fprintf(fpcom, "%i\t%i\n", 1 + nodo + acumulado_id[comm], comm + 1);
					if (nodo == 0) {
						fprintf(fpnod, "%i", 1 + nodo + acumulado_id[comm]);
					}
					else {
						fprintf(fpnod, "\t%i", 1 + nodo + acumulado_id[comm]);
					}
				}
				fprintf(fpnod, "\n"); 
			}
			fclose(fpnod);
			fclose(fpcom);
		}
		else {
			if (fpnod == NULL) printf("Error opening nodes_comm_file  %s\n", network_file);
			if (fpcom == NULL) printf("Error opening comunity_file %s\n", comunity_file);
			return 1;
		}
	}
#if defined _MSC_VER
	fopen_s(&fpnet,network_file, "w");
#elif defined __GNUC__
	fpnet = fopen(network_file, "w");
#endif
//debugging
	 //check if links between community 0 and others have biggest weight than the links internal !!
	 for (comm = 0; comm < num_comunities; comm++) {
		 acumulado_id[comm] = (comm == 0) ? 0 : acumulado_id[comm - 1] + comunity_num_nodes[comm - 1];
	 }
	 for (unsigned int comm_src = 0; comm_src < num_comunities; comm_src++) {
		 for (unsigned int node_src = 0; node_src < comunity_num_nodes[comm_src]; node_src++) {
			 if (comunidades != NULL) 
			 for (unsigned int temp_link = 0; temp_link < comunity_size_node_links[comm_src][node_src]; temp_link++) {
				 if (comunidades[comm_src][node_src][temp_link].link_node != UINT_MAX
					 ||
					 comunidades[comm_src][node_src][temp_link].link_comunity != UINT_MAX) {
					if(comunidades[comm_src][node_src][temp_link].link_comunity>= num_comunities){
						printf("Error comm_src %d node_src %d temp_link %d\n", comm_src, node_src, temp_link);fflush(stdout);
						printf("comunidades[comm_src][node_src][temp_link].link_comunity>= num_comunities\n");
						printf("link_comunity %d >= %d\n",comunidades[comm_src][node_src][temp_link].link_comunity, num_comunities);
						printf("link_node %d\n",comunidades[comm_src][node_src][temp_link].link_node);
						fflush(stdout);exit(1);
					}
					printf("ORIG %d[%d]->%d[%d] %f (%d)", node_src, comm_src,
					 	comunidades[comm_src][node_src][temp_link].link_node,
					 	comunidades[comm_src][node_src][temp_link].link_comunity,
					 	comunidades[comm_src][node_src][temp_link].link_weight_out,
					 	temp_link);
					printf(" ==> %d=%d\n", 1 + node_src + acumulado_id[comm_src],
						 1 + comunidades[comm_src][node_src][temp_link].link_node + acumulado_id[comunidades[comm_src][node_src][temp_link].link_comunity]);
				}
			 }
			 if (comunidades_otra != NULL) 
			 //int comm_external_links_this = comm_external_links[comm_src];				
			 for (unsigned int temp_link = 0; temp_link < comunity_size_node_links_otra[comm_src][node_src]; temp_link++) {
				 if (comunidades_otra[comm_src][node_src][temp_link].link_node != UINT_MAX
					 ||
					 comunidades_otra[comm_src][node_src][temp_link].link_comunity != UINT_MAX) {
					 printf("OTRA %d[%d]->%d[%d] %f (%d) ==> %d=%d\n", node_src, comm_src,
					 	comunidades_otra[comm_src][node_src][temp_link].link_node,
					 	comunidades_otra[comm_src][node_src][temp_link].link_comunity,
					 	comunidades_otra[comm_src][node_src][temp_link].link_weight_out,
					 	temp_link,
						 1 + node_src + acumulado_id[comm_src],
						 1 + comunidades_otra[comm_src][node_src][temp_link].link_node + acumulado_id[comunidades_otra[comm_src][node_src][temp_link].link_comunity]);
				 }
			 }
			 printf("\n");
		 }
	 }

	if (fpnet != NULL) {
		//los nodos dentro de cada comunidad los hemos numerado desde 0.
		//ahora queremos que el id de un nodo de la comunidad "comm+1" sea el id+1 del ultimo nodo de la comunidad "comm"
		for (comm = 0; comm < num_comunities; comm++) {
			acumulado_id[comm] = (comm == 0) ? 0 : acumulado_id[comm - 1] + comunity_num_nodes[comm - 1];
		}
		for (comm = 0; comm < num_comunities; comm++) {
			for (nodo = 0; nodo < comunity_num_nodes[comm]; nodo++) {
				if (comunidades != NULL) {
					for (enlace = 0; enlace != comunity_size_node_links[comm][nodo]; enlace++) {
						if (comunidades[comm][nodo][enlace].link_node != UINT_MAX) {
							dest = comunidades[comm][nodo][enlace].link_node;
							dest_comm = comunidades[comm][nodo][enlace].link_comunity;
							fprintf(fpnet, "%i\t%i\t%7.2f\n", 1 + nodo + acumulado_id[comm],
							1 + dest + acumulado_id[dest_comm],
							comunidades[comm][nodo][enlace].link_weight_out);
							//el nodo id minimo en R es 1
						}
					}
				}
				if (comunidades_otra != NULL) {
					for (enlace = 0; enlace != comunity_size_node_links_otra[comm][nodo]; enlace++) {
						if (comunidades_otra[comm][nodo][enlace].link_node != UINT_MAX) {
							dest = comunidades_otra[comm][nodo][enlace].link_node;
							dest_comm = comunidades_otra[comm][nodo][enlace].link_comunity;
							fprintf(fpnet, "%i \t%i\t%7.2f\n", 1 + nodo + acumulado_id[comm], 
							1 + dest + acumulado_id[dest_comm],
							comunidades_otra[comm][nodo][enlace].link_weight_out);
							//el nodo id minimo en R es 1
						}
					}
				}
			}
		}
		fclose(fpnet);
	}
	else {
		if (fpnet == NULL) printf("Error opening network_file  %s\n", network_file);
		return 1;
	}
	free(acumulado_id);
	return 0;
}
