#include <stdio.h>
#include <errno.h>
#include "grafos.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
void makedir(const char path[]){ //"/some/directory"
	struct stat st = {0};
	if (stat(path, &st) == -1) 
		mkdir(path, 0700); 
}

char* itoa (int value, char *result, int base){
	// check that the base if valid
	if (base < 2 || base > 36) { *result = '\0'; return result; }
	char* ptr = result, *ptr1 = result, tmp_char;
	int tmp_value; 
	do {
		tmp_value = value;
		value /= base;
		*ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
	} while ( value ); 
	// Apply negative sign
	if (tmp_value < 0) *ptr++ = '-';
	*ptr-- = '\0';
	while (ptr1 < ptr) {
		tmp_char = *ptr;
		*ptr--= *ptr1;
		*ptr1++ = tmp_char;
	}
	return result;
}

//return 1 if error
int write_r_file( const unsigned int num_comunities,
	unsigned int *comunity_num_nodes, 
	char *comunity_file,//for instance "community.dat"
	char *network_file, //for instance "network.dat"
	char *nodes_by_comm_file, //for instance "nodes_by_comm.dat"
	unsigned int ** comunity_size_node_links, _link*** comunidades,
	unsigned int ** comunity_size_node_links_otra, _link*** comunidades_otra)
{
	unsigned int *acumulado_id= (unsigned int *) malloc(num_comunities*sizeof(unsigned int )); 
	unsigned int comm,nodo,dest,enlace,dest_comm;
	FILE *fpcom,*fpnet, *fpnod;
if( (access(nodes_by_comm_file, F_OK)!=-1) || (access(comunity_file, F_OK)!=-1) ){ //other values that can be tested are R_OK W_OK X_OK, from unistd.h
	printf(" Not generating community file because it already exists\n");
}else{//file not exists
	fpnod = fopen(nodes_by_comm_file,"w");
	fpcom = fopen(comunity_file,"w");
	if((fpnod!=NULL)&&(fpcom!=NULL)) { 
		//los nodos dentro de cada comunidad los hemos numerado desde 0.
		//ahora queremos que el id de un nodo de la comunidad "comm+1" sea el id+1 del ultimo nodo de la comunidad "comm"
		for(comm=0;comm<num_comunities;comm++)
			acumulado_id[comm]= (comm==0)? 0 : acumulado_id[comm-1]+comunity_num_nodes[comm-1];
		for(comm=0;comm<num_comunities;comm++){ 
			for(nodo = 0;nodo < comunity_num_nodes[comm];nodo++){
				fprintf(fpcom,"%i\t%i\n",1+nodo+acumulado_id[comm],comm+1);
				if(nodo==0){
					fprintf(fpnod,"%i",1+nodo+acumulado_id[comm]); 
				}else{
					fprintf(fpnod,"\t%i",1+nodo+acumulado_id[comm]); 				
				} 
			}
			fprintf(fpnod,"\n"); 
		} 
		fclose(fpnod);
		fclose(fpcom); 
	}else {
		if(fpnod==NULL) printf("Error opening nodes_comm_file  %s\n",network_file);
		if(fpcom==NULL) printf("Error opening comunity_file %s\n",comunity_file);
		return 1;
	}
}
	fpnet = fopen(network_file,"w");
	if(fpnet!=NULL) { 
		//los nodos dentro de cada comunidad los hemos numerado desde 0.
		//ahora queremos que el id de un nodo de la comunidad "comm+1" sea el id+1 del ultimo nodo de la comunidad "comm"
		for(comm=0;comm<num_comunities;comm++)
			acumulado_id[comm]= (comm==0)? 0 : acumulado_id[comm-1]+comunity_num_nodes[comm-1];
		for(comm=0;comm<num_comunities;comm++){ 
			for(nodo = 0;nodo < comunity_num_nodes[comm];nodo++){ 
				if(comunidades!=NULL){
					for(enlace = 0;enlace != comunity_size_node_links[comm][nodo];enlace++){ 
						if(comunidades[comm][nodo][enlace].link_node !=UINT_MAX){ 
							dest=comunidades[comm][nodo][enlace].link_node; 
							dest_comm=comunidades[comm][nodo][enlace].link_comunity;
							fprintf(fpnet,"%i\t%i\t%7.2f\n",1+nodo+acumulado_id[comm],
								1+dest+acumulado_id[dest_comm], 
								comunidades[comm][nodo][enlace].link_weight_out); 
							//el nodo id minimo en R es 1 
						}
					}
				} 
				if(comunidades_otra!=NULL){
					for(enlace = 0;enlace != comunity_size_node_links_otra[comm][nodo];enlace++){
						if(comunidades_otra[comm][nodo][enlace].link_node !=UINT_MAX){
							dest=comunidades_otra[comm][nodo][enlace].link_node;
							dest_comm=comunidades_otra[comm][nodo][enlace].link_comunity;
							fprintf(fpnet,"%i\t%i\t%7.2f\n",1+nodo+acumulado_id[comm], 
								1+dest+acumulado_id[dest_comm], 
								comunidades_otra[comm][nodo][enlace].link_weight_out);
							//el nodo id minimo en R es 1 
						}
					}
				} 
			} 
		}   
		fclose(fpnet);
	}else {
		if(fpnet==NULL) printf("Error opening network_file  %s\n",network_file); 
		return 1;
	}
	free(acumulado_id);
	return 0;
}
