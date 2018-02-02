//Author: J.M. Montañana, April 2017
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <math.h>
#include "grafos.h"
#include "write_r.h"
#include "tamano_com.h"

char method[]="dbol";
char folder_name[]="data_dbol";

//node es un valor entre 0 y num_nodes-1
const unsigned int max_num_links = 20;

//function para calcular el maximo de total_enlaces_salida de un nodo
unsigned int popularidad(const float node, const float num_nodes){
	const unsigned int num_links = 1.2*sqrt( num_nodes);
	//numero de links (aristas) es f(node)= ((1-k)/num_nodes) *node + k
	// el maximo es num_links
	const float k=num_links;
	float i = ((1.0-k)/num_nodes) *node + k;
	return i;	
}

//funcion para elegir nodos visitados
unsigned long long int popularidad_de (const unsigned int node, const unsigned num_nodes){
//	return (2<<(num_nodes-node)); //funcion parabolica, problema si (num_nodes-nodes)>64
//	f(x)=d* [ a/b * (x/b)^(a−1) * exp(− (x/b)^a) ] +c Funcion Weibull
	double a = 1.22737455473989;
	double b =176.985725394169;
	double c =1.21672489858613;
	double d =1536.90835967452;
	//calculada para un total puntos es: 57 para un intervalo de x=[0:840], por lo que toca aproximadamente un punto por cada x=14
	double x= node *14.0;
	//#include <math.h> printf("9 elevado a 15 es %lf\n", pow (9, 15));
	double termino1 = (a/b) * pow( x/b, a-1.0);
	double termino2 = -1.0 *pow( x/b,a);
	unsigned long long int result = (unsigned int) (d* (termino1 * exp(termino2))+c);
	return result;
}


 
	
unsigned int rand_popular_node(const unsigned int num_nodes, const unsigned long long int max_rand, const unsigned int last, unsigned int* vector_nodes_available){
	//default definition of max_rand is:
	// unsigned int max_rand = 0;//(2<<(num_nodes+1))-1;// es el sumatorio desde i=0 hasta num_nodes-1 de (2<<i) 
	// for(node=0;node<num_nodes;node++){
	// 		max_rand+= (2<<(num_nodes-node));
	// 		vector_nodes_available[node]=node;
	// } 
	unsigned long long int opt= rand() % max_rand;
	//buscamos a que nodo corresponde
	unsigned int node=0;
	unsigned long long int accum=0;
	unsigned int encontrado=false;
	while((node<last)&&(encontrado==false)){
		accum+=popularidad_de(vector_nodes_available[node], num_nodes);// (unsigned long long int ) (2<<(num_nodes-vector_nodes_available[node]));
		if(opt<accum){
			encontrado=true;
		}else{
			node++;
		}		
	} 
	if(node==last){
		unsigned long long int mio=0;
		for(unsigned int i=0;i<last;i++) 
		mio+= popularidad_de(vector_nodes_available[i], num_nodes);// (unsigned long long int) (2<<(num_nodes-vector_nodes_available[i]));
		printf("valor correcto de max_rand seria %llu != max_rand %llu\n", mio, max_rand);
		printf("Error no encotramos nodo node %i=last %i\n", node, last);
		printf(" opt %llu max_rand %llu\n", opt, max_rand);
		exit(1);
	} 		
	return node;
}




int main(int argc, char *argv[]) {// argv[i] from i = 0; to i<argcv;
	unsigned int i;
	unsigned int num_nodes=0;
	unsigned int num_communities=8;
	// unsigned int extern_links=1;
	const int debug=false; 
	unsigned int myseed = 100;
	for(i=1;i<argc;i++){//i=0 es el nombre del programa
		if(argv[i][0]=='-'){
			if(strcmp(argv[i],"-n")==0){ 
				if(argc>i){ 
					printf("parametro nodos es '%s'\n",argv[i+1]);
					num_nodes=atoi(argv[i+1]);
				}else{
					printf("error parametro de numero de nodos sin valor asignado\n");exit(1);
				}
			//} else if(strcmp(argv[i],"-l")==0){
			//	if(argc>i){					
			//		printf("parametro enlaces es '%s'\n",argv[i+1]);
			//	}else{
			//		printf("error parametro de numero de enlaces sin valor asignado\n");exit(1);
			//	}
			} else if(strcmp(argv[i],"-c")==0){
				if(argc>i){
					printf("parametro comunidades es '%s'\n",argv[i+1]);
					num_communities=atoi(argv[i+1]);
				}else{
					printf("error parametro de comunidades sin valor asignado\n");exit(1);
				}
			} else if(strcmp(argv[i],"-s")==0){
				if(argc>i){
					printf("parametro semilla es '%s'\n",argv[i+1]);
					myseed=atoi(argv[i+1]);
				}else{
					printf("error parametro semilla sin valor asignado\n");exit(1);
				}
			} else if(strcmp(argv[i],"-h")==0){ 
				printf(" uso de la funcion:\n %s -n number_of_vertices(nodes) -c number_of_communities -s value_of_rand_seed\n",argv[0]);
				printf(" el uso de cualquiera de los parametros es opcional, los valores por defecto son:\n");
				printf(" number_of_vertices = 0 , number_of_communities = 8, seed = 100\n");
				printf("    number of vertices = 0 implica una seleccion aleatoria del tamaño de cada comunidad");
				printf("    number of vertices > 0 asigna ese valor como tamaño para todas las comunidades ");
				exit(1); 
			}
		}
	}
	srand(myseed);//((unsigned) time(&t));
	//definimos la estructura de comunidades:
	//una comunidad tiene un total de enlaces de: num_links*num_nodes
	// el total de enlaces de salida de un nodo, lo define la funcion popularidad
	_link*** comunidades;//es un array de comunidades, cada comunidad es un array de enlaces.
	comunidades = ( _link ***) malloc(num_communities*sizeof( _link **));
	unsigned int **comunity_size_node_links= (unsigned int **) malloc(num_communities*sizeof(unsigned int *));	
	unsigned int enlace,comm,node, total_enlaces_salida;
	unsigned int *comunity_num_nodes= (unsigned int *) malloc(num_communities*sizeof(unsigned int ));
	unsigned int total_nodes=0;	
	for(comm=0;comm<num_communities;comm++){ 
		if (num_nodes==0){
			comunity_num_nodes[comm]= get_new_community_size( num_communities);
		}else{
			comunity_num_nodes[comm]=num_nodes;
		}	
		total_nodes+=comunity_num_nodes[comm];
		comunity_size_node_links[comm] = (unsigned int *) malloc(comunity_num_nodes[comm]*sizeof(unsigned int ));
		if(comunity_size_node_links[comm] == NULL){ 
			printf("Out of memory\n");
			exit(-1);
		}
		comunidades[comm]= ( _link **) malloc(comunity_num_nodes[comm]*sizeof( _link *));
		if(comunidades[comm] == NULL){ 
			printf("Out of memory\n");
			exit(-1);
		}		
		for(node = 0;node != comunity_num_nodes[comm]; node++){ 
			total_enlaces_salida= popularidad(node,comunity_num_nodes[comm]);	 
			comunity_size_node_links[comm][node]=total_enlaces_salida;		
			comunidades[comm][node]=( _link *) malloc(total_enlaces_salida*sizeof( _link ));
			if(comunidades[comm][node] == NULL){ 
				printf("Out of memory\n");
				exit(-1);
			} 
		} 
		//ahora asignamos valores a las posiciones del array comunidades
		for(node = 0; node != comunity_num_nodes[comm]; node++)
			for(enlace=0;enlace<comunity_size_node_links[comm][node];enlace++) {
				comunidades[comm][node][enlace].link_node=UINT_MAX;
				comunidades[comm][node][enlace].link_weight_out=0;
			}
		//comunidades[comm][node].link_node = UINT_MAX;// UINT_MAX stands for undefined destination node
	} 		
	//ahora hay que generar rutas sobre los nodos
	//las primeras rutas definen donde van los enlaces de salida
	//las ultimas rutas tienen que usar los enlaces ya definidos
	unsigned int ruta,step,nodo;
	//unsigned int num_rutas=500;
	unsigned int long_estud_preferen=7;
	unsigned int* path_nodes;
	unsigned int* path_links;	
	unsigned int next_node;
	//hacemos un vector de los posibles destinos
	unsigned int* nodos_no_visitados;
	unsigned int* enlaces_a_nodos_no_visitados;
	enlaces_a_nodos_no_visitados = (unsigned int *) malloc(max_num_links*sizeof(unsigned int ));
	unsigned int total_enlaces_a_nodos_no_visitados;
	//float newpeso;
        // El peso pasa a ser la probabilidad de que un alumno sea admitido en la opción i, i de 1 a 10
        //float peso[]={0.570464,0.164352,0.090634,0.069855,0.03752,0.026522,0.015767,0.011521,0.007588,0.005778};
          float peso[]={1.05 ,1, 1, 1.21, 1.21, 1.21, 1.21, 1.21};
  	//float peso[]={1.4,0.01,0.01,1.1,1.21,1.31,1.41,1.51};
	for(comm=0;comm<num_communities;comm++){  
		nodos_no_visitados = (unsigned int *) malloc(comunity_num_nodes[comm]*sizeof(unsigned int )); 
		unsigned int* vector_nodes_available=( unsigned int *) malloc(comunity_num_nodes[comm]*sizeof(unsigned int ));
		path_nodes= (unsigned int *) malloc(comunity_num_nodes[comm]*sizeof(unsigned int ));
		path_links= (unsigned int *) malloc(comunity_num_nodes[comm]*sizeof(unsigned int )); 
		unsigned long long int max_rand,new_rand;
printf(" total numrutas %i para la comm %i\n", comunity_num_nodes[comm]*20,comm);fflush(stdout);
		for(ruta=0;ruta<200+comunity_num_nodes[comm]*10;ruta++){
			if(debug) printf(" ruta %i max_nodes=%i\n",ruta, comunity_num_nodes[comm]);fflush(stdout);
			max_rand = 0; //(2<<(num_nodes+1))-1;// es el sumatorio desde i=0 hasta num_nodes-1 de (2<<i)  
			//sumamos la probabilidad de cada nodo segun su popularidad
			for(node=0;node<comunity_num_nodes[comm];node++){
				new_rand = popularidad_de(node, comunity_num_nodes[comm]);// 
				if( max_rand > ULLONG_MAX-new_rand){
					printf("Error, network too big\n");fflush(stdout);
					exit(1);
				}
				max_rand+=new_rand;
				vector_nodes_available[node]=node; 
			}  
			unsigned int last =comunity_num_nodes[comm];		
			//next_node= rand() % comunity_num_nodes[comm];//el primer nodo de la ruta lo elegimos aleatoriamente
			next_node= rand_popular_node( comunity_num_nodes[comm], max_rand, last, vector_nodes_available);
			max_rand-= popularidad_de(vector_nodes_available[next_node], comunity_num_nodes[comm]);//   
			unsigned int sin_salida_posible=false; 
			//newpeso=2.0;
			unsigned int total_nodos_no_visitados;
			unsigned int path_length=0;
			path_nodes[path_length]=next_node;
			if(long_estud_preferen >= (sizeof(peso)/sizeof(peso[0]))){
				printf("Error long_estud_preferen %u > than weight list %lu\n",long_estud_preferen, (sizeof(peso)/sizeof(peso[0])));
				fflush(stdout);exit(1);
			}
			for(step=0;step<long_estud_preferen; step++){
				node=next_node; 
				//quitamos de la lista completa, aquellos nodos que estan visitados
                		unsigned int visited;
                		total_nodos_no_visitados=0;
                		for(nodo=0;nodo<comunity_num_nodes[comm];nodo++){
							visited=false;
							i=0;
							while((i<=path_length)&&(visited==false)){
								if(nodo==path_nodes[i]){
									visited=true;
								}else{
									i++;
								}
							}
							if(visited==false){
								nodos_no_visitados[total_nodos_no_visitados]=nodo; //lista completa de todos los nodos
								total_nodos_no_visitados++;
							}
                		}    
				//la ruta ya contiene un nodo inicial, que es el next_node
				if(debug) printf(" step %i,total_nodos_no_visitados =%i total nodes %i\n",step,total_nodos_no_visitados, comunity_num_nodes[comm]);
				if((sin_salida_posible==false)&&(total_nodos_no_visitados>0)){
					//si existe algun enlace de salida no conectado
					//entonces podemos elegir cualquier destino aletoriamente no visitado en esta ruta
					enlace=0;
					unsigned int encontrado_enlace_libre=false;
					// printf("buscando enlaces libres node %i total enlaces %i\n",node, comunity_size_node_links[comm][node]);
					while((encontrado_enlace_libre==false)&&(enlace<comunity_size_node_links[comm][node])) { 
						if(comunidades[comm][node][enlace].link_node==UINT_MAX){
							encontrado_enlace_libre=true;
							// printf("enlace encontrado libre node=%i enlace=%i\n",node,enlace);
						}else{
							// printf("nodo %i enlace %i ocupado dest=%i\n",node,enlace,comunidades[comm][node][enlace].link_node);
							enlace++;
						} 
					}
					// el primer enlace no conectado, si existe alguno es "enlace" 
					if(encontrado_enlace_libre==true){
						if(debug) printf(" nodo encontrado\n");
						//buscamos un nodo aleatoriamente, de entre los no visitados 
						max_rand = 0; //(2<<(comunity_num_nodes[comm]+1))-1;// es el sumatorio desde i=0 hasta comunity_num_nodes[comm]-1 de (2<<i) 
						if (debug) printf("total nodos no visitados %i\n",total_nodos_no_visitados);
						for(i=0;i<total_nodos_no_visitados;i++){
						//if (debug) printf("comunity_num_nodes[comm] %i nodos_no_visitados[i]%i\n",comunity_num_nodes[comm],nodos_no_visitados[i]);
							//printf("-- max_rand %u\n",max_rand);							
							new_rand= popularidad_de(nodos_no_visitados[i], comunity_num_nodes[comm]);//2<<(comunity_num_nodes[comm]-nodos_no_visitados[i])
							if( max_rand > ULLONG_MAX-new_rand){
								printf("Error, network too Big\n");fflush(stdout);
								exit(1);
							}
							max_rand+=new_rand;
						}						
						if (debug) printf(" max_rand %llu\n",max_rand);
						if(max_rand==0) {
							printf("error dim too large, max_rand=0\n");exit(1);
						}
						i= rand_popular_node( comunity_num_nodes[comm], max_rand, total_nodos_no_visitados,nodos_no_visitados); 
						if(debug) printf(" rand_popular_node %i\n",i);
						next_node=nodos_no_visitados[i]; 
                        if(next_node==path_nodes[path_length]){
                            printf("Error, repitiendo nodo %i\n",next_node);
                            exit(1);
                        } 
						//next_node=nodos_no_visitados[rand() % total_nodos_no_visitados]; 
						//buscamos si tiene conexion previa, en caso contrario usamos el "enlace" no definido
						if(next_node==node){
							printf("error next node=%i ,total_nodos_no_visitados %i\n",next_node,total_nodos_no_visitados);
							for(unsigned int i=0;i< total_nodos_no_visitados;i++)
								printf(" %i",nodos_no_visitados[i]);
							printf("\n");
							exit(1);
						}
						unsigned int temp_link=0;
						unsigned int encontrada_conexion_previa=false;
						while((temp_link<comunity_size_node_links[comm][node])&&(encontrada_conexion_previa==false)){
							if(comunidades[comm][node][temp_link].link_node==next_node){
								encontrada_conexion_previa=true;
							}else{
								temp_link++;
							}
						}
						if(encontrada_conexion_previa==true){
							enlace=temp_link;
						}else if(comunidades[comm][node][enlace].link_node!=UINT_MAX){
							printf("error enlace already defined\n");exit(1);
						} 
						if(debug) printf(" definiendo com %i nodo %i enlace %i\n",comm,node,enlace);
						comunidades[comm][node][enlace].link_node=next_node;
						comunidades[comm][node][enlace].link_comunity=comm;					
						path_links[path_length]=enlace;	
						path_length+=1;
						path_nodes[path_length]=next_node; 
						if(comunidades[comm][node][enlace].link_weight_out==0) 
							comunidades[comm][node][enlace].link_weight_out=1;
						//long long int jajaa=popularidad_de(next_node, comunity_num_nodes[comm]);
						comunidades[comm][node][enlace].link_weight_out*= peso[step];//newpeso;							
						//newpeso=newpeso/2.0;
					}else{ 
						// printf("node %i no quedan enlaces indefinidos, total enlaces %i\n",node,comunity_size_node_links[comm][node]);
						//no quedan enlaces libres, tenemos que escoger uno de los enlaces del nodo, que vaya a un nodo no visitado
						total_enlaces_a_nodos_no_visitados=0;
						//vamos a hacer un vector con los enlaces definidos que no llevan a ningun nodo visitado
						for(enlace=0;enlace<comunity_size_node_links[comm][node];enlace++){
							unsigned int next_node=comunidades[comm][node][enlace].link_node;
							if(next_node==UINT_MAX){
								printf("error enlace=%i\n",enlace);exit(1);
							}else{ 
								unsigned int i=0;
								while((i<total_nodos_no_visitados)&&(nodos_no_visitados[i]!=next_node)){
									i++;
								}
								if(i<total_nodos_no_visitados){
									enlaces_a_nodos_no_visitados[total_enlaces_a_nodos_no_visitados]=enlace;
									total_enlaces_a_nodos_no_visitados++;
								} 
							}
						}
						// printf("total_enlaces_a_nodos_no_visitados %i\n",total_enlaces_a_nodos_no_visitados);
						// printf("enlace 0 lleva a %i\n",comunidades[comm][node][enlace].link_node);
						//si existe algun enlace posible, tomamos uno aleatoriamente
						if(total_enlaces_a_nodos_no_visitados!=0){ 
							enlace=enlaces_a_nodos_no_visitados[rand()% total_enlaces_a_nodos_no_visitados]; 
							next_node=comunidades[comm][node][enlace].link_node;
							path_links[path_length]=enlace; 
							if(next_node==path_nodes[path_length]){
								printf("Error, repited node %i ...\n",next_node);
								exit(1);
							} 
							path_length+=1;
							path_nodes[path_length]=next_node; 
							path_length=step+1;				
							if (comunidades[comm][node][enlace].link_weight_out==0) comunidades[comm][node][enlace].link_weight_out=1;			
							//long long int jajab= popularidad_de(next_node, comunity_num_nodes[comm]);
							comunidades[comm][node][enlace].link_weight_out*= peso[step];//newpeso;							
							//newpeso=newpeso/2.0;
						}else{
							sin_salida_posible=true;
						}					
					} 
				}//if sin_salida_posible==false
			}//for steps
			// printf(" Community %i Ruta %i:",comm, ruta);
			for(unsigned int sstep=0;sstep<=path_length; sstep++){
				// printf(" %i",path_nodes[sstep]);
				if(sstep>0){
					if (path_nodes[sstep]==path_nodes[sstep-1]){
						printf("error enlace repetido step %i path lenght %i\n",sstep, path_length);
						printf("path_nodes[0] =%i, path_nodes[1]=%i\n", path_nodes[0],path_nodes[1]);
						exit(1);
					}
				}
			}
			// printf("\n");	 
		} //for rutas
		free(vector_nodes_available);
		free(path_nodes);
		free(path_links);	
		free(nodos_no_visitados);
	}//for comm 
     free(enlaces_a_nodos_no_visitados);
	
	//buscamos el nodo mas popular
	for(comm=0;comm<num_communities;comm++){
		float maxpeso=0.0;
		unsigned int nodo_maxpeso=0;//si la comunidad solo tiene un nodo, el de mas peso sera el 0, y no se podra realizar el bucle for
		for(node = 0; node != comunity_num_nodes[comm]; node++){
			float totalpeso=0.0;
			for(int enlace=0;enlace <comunity_size_node_links[comm][node];enlace++){
				if(comunidades[comm][node][enlace].link_node != UINT_MAX){
					totalpeso+=comunidades[comm][node][enlace].link_weight_out;
					if(totalpeso>=maxpeso){
						maxpeso=totalpeso;
						nodo_maxpeso=node;
					}
				}
			}
		}
		//printf("nodo con mas peso de salida es nodo %i peso %f comm %i\n", nodo_maxpeso, maxpeso, comm);
        //buscamos nodos que han quedado completamente sin conectar, y los conectamos
        for(node = 0; node != comunity_num_nodes[comm]; node++){
			if(comunity_size_node_links[comm][node]==0){
				printf(" Error nodo %i sin conexion. comm %i\n",node,comm); 
				printf("this code only for duol");exit(1);
				//unsigned int total_links=comunity_size_node_links[comm][node];
				//_link *previo=comunidades[comm][node];
				//comunidades[comm][node]=( _link *) malloc(total_links+50*sizeof( _link ));
				//for(int temp_link=total_links;temp_link<total_links+50;temp_link++)
				//	comunidades[comm][node][temp_link].link_node=UINT_MAX;
				//free(previo);
				//comunity_size_node_links[comm][node]=total_links+50;
				//enlace=total_links;
				//comunidades[comm][node][enlace].link_node =nodo_maxpeso;
				//comunidades[comm][node][enlace].link_comunity =comm;//esto es la comunidad del nodo destino
				//comunidades[comm][node][enlace].link_weight_out = (2<<(long_estud_preferen-0))/16.0;//el peso va por el numero de salto
			}else{   
				unsigned int conectado=0;
				for(int enlace=0;enlace <comunity_size_node_links[comm][node];enlace++)
					if(comunidades[comm][node][enlace].link_node != UINT_MAX)
						conectado=1;
				if(conectado==0){   
					//if (debug)
						printf(" node %i[%i] no conected, create connection to node %i\n",node,comm, nodo_maxpeso);
					comunidades[comm][node][0].link_node = nodo_maxpeso;
					comunidades[comm][node][0].link_comunity =comm;//esto es la comunidad del nodo destino
					//float newpeso=2.0;
					//newpeso=newpeso/64.0;
				        //long long int jajac= popularidad_de(nodo_maxpeso, comunity_num_nodes[comm]);	
					comunidades[comm][node][0].link_weight_out=peso[1];//newpeso;
				}
			}
        }//for node
    }//for comm 
	

 

normaliza_pesos_salida(num_communities, comunidades, NULL, comunity_num_nodes, comunity_size_node_links, NULL);

	print_graph_double(num_communities, comunidades, NULL, comunity_num_nodes, comunity_size_node_links, NULL);
	/////////////////////////////////////////////////////////////////////////////////////////////////////// 
	// ahora conectamos enlaces entre comunidades 
	char comm2[50];
	char net2[50];
	char nodes_by_comm_file[50];
	sprintf(comm2, "%s%s%s_%d.dat", folder_name,"/community_",method,myseed);
	sprintf(net2 , "%s%s%s_%d_0.dat", folder_name,"/network_",method,myseed);
	sprintf(nodes_by_comm_file, "%s%s%s_%d.dat", folder_name,"/nodes_by_comm_",method,myseed);
 	makedir(folder_name);	 
	write_r_file(num_communities, comunity_num_nodes, comm2, net2, nodes_by_comm_file,comunity_size_node_links, comunidades, NULL, NULL);
	
	_link*** comunidades_otra; 
	comunidades_otra = ( _link ***) malloc(num_communities*sizeof( _link **));	 
	unsigned int **comunity_size_node_links_otra= (unsigned int **) malloc(num_communities*sizeof(unsigned int *)); 
	unsigned int **total_used_links_at_node= (unsigned int **) malloc(num_communities*sizeof(unsigned int *));	
	for(int comm=0;comm<num_communities;comm++){	 
		comunidades_otra[comm] = ( _link **) malloc(comunity_num_nodes[comm]*sizeof( _link *)); 
		comunity_size_node_links_otra[comm] = (unsigned int *) malloc(comunity_num_nodes[comm]*sizeof(unsigned int )); 
		total_used_links_at_node[comm] = (unsigned int *) malloc(comunity_num_nodes[comm]*sizeof(unsigned int )); 
	}

	for(int comm=0;comm<num_communities;comm++){  
		for(int nodo = 0; nodo < comunity_num_nodes[comm]; nodo++){						
			comunity_size_node_links_otra[comm][nodo]=0; //tamaño de memoria reservada para comunidades_otra[comm][nodo] 
			comunidades_otra[comm][nodo]=NULL; 
		}
	}
	float *comm_external_links= (float *) malloc(num_communities*sizeof(float ));

	//se crean enlaces diferentes en cada caso de n enlaces entre comunidades
	int max_mu_extern_links=9; 
	unsigned int mu_extern_links;
	float extern_links; 
	for(mu_extern_links =1;mu_extern_links< max_mu_extern_links;mu_extern_links++){ //vamos de 10% 20% 30% ... 60% de enlaces extra
		printf(" mu_extern_links %.1f%%\n",mu_extern_links*10.0);
		for(int comm=0;comm<num_communities;comm++){	 
			extern_links=0.0;
			for(int nodo = 0; nodo < comunity_num_nodes[comm]; nodo++){ 
				//printf("links node %u[%u] are %u\n", nodo, comm, comunity_size_node_links[comm][nodo]);
				extern_links+= (comunity_size_node_links[comm][nodo] * mu_extern_links);   
				total_used_links_at_node[comm][nodo]=0; 
  				//comunidades_otra[comm][nodo]=( _link *) malloc( (comunity_size_node_links_otra[comm][nodo])*sizeof( _link )); 
				//for(int enlace = 0; enlace != comunity_size_node_links_otra[comm][nodo]; enlace++) {
				//	comunidades_otra[comm][nodo][enlace].link_node = UINT_MAX;// UINT_MAX stands for undefined destination node	
				//	comunidades_otra[comm][nodo][enlace].link_comunity = UINT_MAX;// UINT_MAX stands for undefined destination node
				//} 
			} 
			comm_external_links[comm] =extern_links / 10.0;//total de enlaces nuevos para la comunidad, no sabemos cuantos tiene cada comunidad
			//printf(" *** Enlaces extra para la comunidad [%u] is %.1f\n", comm, comm_external_links[comm]);
		}  	 	
		// mu_extern_links stands for mu from 0.1 ..0.2 .. 1.0.   0.1 means add a 10% of links of a community with the other comunities
 		interconnecta_comunidades_new( num_communities, comunidades_otra, comunity_size_node_links_otra, comunity_num_nodes, total_nodes, comm_external_links, total_used_links_at_node, mu_extern_links); 
		//interconnecta_comunidades( num_communities, comunidades_otra, comunity_size_node_links_otra, comunity_num_nodes, total_nodes); 
		
		normaliza_pesos_salida(num_communities, comunidades, comunidades_otra, comunity_num_nodes, comunity_size_node_links, comunity_size_node_links_otra);
		print_graph_double(num_communities, comunidades, comunidades_otra, comunity_num_nodes, comunity_size_node_links, comunity_size_node_links_otra);
		sprintf(comm2, "%s%s%s_%d.dat", folder_name,"/community_",method,myseed);
		sprintf(net2 , "%s%s%s_%d_%d.dat", folder_name,"/network_",method, myseed, mu_extern_links);
		sprintf(nodes_by_comm_file, "%s%s%s_%d.dat", folder_name,"/nodes_by_comm_",method,myseed);
		write_r_file(num_communities, comunity_num_nodes, comm2,net2,nodes_by_comm_file,
			comunity_size_node_links, comunidades, comunity_size_node_links_otra, comunidades_otra);
 	}

	free(comm_external_links);
	for(int comm=0;comm<num_communities;comm++){	 
		for (int nodo = 0; nodo != comunity_num_nodes[comm]; nodo++){   
			free(comunidades_otra[comm][nodo]);
			free(comunidades[comm][nodo]);
		}
		free(total_used_links_at_node[comm]); 
		free(comunity_size_node_links[comm]); 
		free(comunidades[comm]);  
		free(comunidades_otra[comm]);
		free(comunity_size_node_links_otra[comm]);
	}
	free(total_used_links_at_node);
	free(comunidades_otra );
	free(comunity_size_node_links_otra);
	free(comunity_size_node_links); 
	free(comunidades); 
	
	//free(comunity_used_links);
	free(comunity_num_nodes);
	printf("fin de programa\n");
	// Do not forget to free the memory when you are done: 
	return 0;
}//end main
 
