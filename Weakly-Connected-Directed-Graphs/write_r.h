#if defined _MSC_VER
#pragma once
#endif
#include <stdio.h>
#include <errno.h>
#if defined _MSC_VER
#include "grafos.h"
#endif
void makedir(const char path[]);

char* mitoa(int value, char* result, int base);

int write_r_file(const unsigned int num_comunities,
	unsigned int* comunity_num_nodes,
	char* comunity_file,//for instance "community.dat"
	char* network_file, //for instance "network.dat"
	char* nodes_by_comm_file, //for instance "nodes_by_comm.dat"
	unsigned int** comunity_size_node_links, _link*** comunidades,
	unsigned int** comunity_size_node_links_otra, _link*** comunidades_otra);
