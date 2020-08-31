#include <omp.h>
# include<stdio.h>
# include<stdlib.h>
# include<string.h>
# include<time.h>
# include<math.h>
int *InDegree, *OutDegree, **Edges, **NeighborsIn, **NeighborsOut ,M, N;
float **Pij_weight, *Pij , *Influence_star , *Influence, *Influence1 ,*Susceptibility_star , *Susceptibility , *Susceptibility1 , *Fstar, *Gstar ; 
float deleta = 0.0 ;
int **Neighbor_total_try , **Neighbor_sucess_Infect ; 
void Get_size_of_vertex_and_edges(char **argv) ;
void Init()  ;
void Calculate_degree_neighbors(char **argv) ;
void Input_Fstar_and_Gstar(char **argv) ;
void Input_Pij(char **argv) ;
void Init_Influence_and_Susceptibility(char **argv) ; 
void Iterate_Influence_label(int Iteration,char **argv) ;
void Parallel_Influence_and_Susceptibility(int i) ;
void Exchange_Influence_and_Susceptibility() ;
void Conserve_the_lasttime_I_and_S() ;
float Generate_pij_weight() ;
void Output(char **argv) ;
void FreeMemory() ;
void Output_detaildata(char **argv) ;
float I_S_error(char **argv,int iteration) ; 
int main(int agrs, char **argv){
	srand(time(NULL)) ; 
	printf("read dataset is starting\n") ;
	Get_size_of_vertex_and_edges(argv) ; 	// input the dataset, here the dataset consist of two-tuple.
	printf("read dataset end\n") ;
	Init()  ; 								//initialize the parameters, allocate the memory.
	Calculate_degree_neighbors(argv) ; 		// calculate the degrees(in/out) and set of neighbors(in/out) of all individuals
	Input_Fstar_and_Gstar(argv) ;      		// input the f and g 
	Init_Influence_and_Susceptibility(argv);// initialize the initial condition of I and S.
	Iterate_Influence_label(30000,argv) ;   // start to iterate the I and S, here the first parameter of 
											// the function is the maxmium iteration time if the error of I
											// and S don't converge to given condition.							
	Output(argv) ; // output results
	Output_detaildata(argv) ; // output results
	FreeMemory() ; // free memory
	return 0 ;
	//Note 1: the program can run in parellel when complie it with the command: gcc -fopenmp -o main main.c
	//Note 2: if you don't want to run it in parallel, use the command: gcc -o main main.c. 
	//Note 3: For incompatible reason, you may also need to specify the version of  C, with the command:
	// gcc -std=c99 -fopenmp -o main main.c 
	
}
void Get_size_of_vertex_and_edges(char **argv)
{
	
	char *file_name ;
	char f1[10] = "datasets/";
	char f2[5] = ".txt" ;
	file_name = (char *) malloc((strlen(argv[1])+strlen(argv[5])+strlen(f1)+strlen(f2)+2)*sizeof(char)) ;
	char one_line[100] ;
	sprintf(file_name,"%s%s_%s%s",f1,argv[1],argv[5],f2) ;
	FILE *fout ;
	Edges = (int**)calloc(1,sizeof(int*)) ; 

	int ii = 0 , i , j , maxID = 0 , minID = 100  ,flag = 1 ,count_lines = 0;
	if ((fout = fopen(file_name, "r")) == NULL) { printf("file %sdoes not exist\n",argv[1]); exit(1); }   
	while (fgets(one_line, sizeof(one_line), fout) != NULL)
	{
		if (one_line[0] != '#') {
			Edges[ii] = (int*)calloc(2,sizeof(int)) ;
			flag = sscanf(one_line, "%d\t%d\n", &i, &j);
		//	printf("%d %d\n",i,j)  ;
			Edges[ii][0] = i ;
			Edges[ii][1] = j ;
			ii++ ;
			Edges = (int**)realloc(Edges,(ii+1)*sizeof(int*)) ;
			if (flag != 2) { printf("failed to parse line %s\n", one_line); exit(1); }
			if (i < minID) minID = i;
			if (i > maxID) maxID = i;
			if (j < minID) minID = j;
			if (j > maxID) maxID = j;
			count_lines++;
		}
	}
	fclose(fout);
	Edges = (int**)realloc(Edges,ii*sizeof(int*)) ;
	N  = maxID  ; //the size of vertex
	M = count_lines ;  //the size of edges
	
}

void Init() {
	InDegree =(int*)calloc(N,sizeof(int)) ;	
	OutDegree =(int*)calloc(N,sizeof(int)) ;	
	Influence = (float*)calloc(N,sizeof(float)) ;
	Susceptibility = (float*)calloc(N,sizeof(float)) ;
	Influence1 = (float*)calloc(N,sizeof(float)) ;
	Susceptibility1 = (float*)calloc(N,sizeof(float)) ;
	Fstar = (float*)calloc(N,sizeof(float)) ;
	Gstar = (float*)calloc(N,sizeof(float)) ;
	Influence_star = (float*)calloc(N,sizeof(float)) ;
	Susceptibility_star = (float*)calloc(N,sizeof(float)) ;
	Pij = (float*)calloc(M,sizeof(float)) ;
}
void Calculate_degree_neighbors(char **argv)
{
	int *countOut, *countIn  ;
	for(int i = 0 ; i < M ;i++)
	{	Edges[i][0]-- ;
	    Edges[i][1]-- ;
		OutDegree[Edges[i][0]]++ ;
		InDegree[Edges[i][1]]++ ;
	} 
	NeighborsOut = (int**)calloc(N,sizeof(int*)) ;
	NeighborsIn = (int**)calloc(N,sizeof(int*)) ;
	countOut = (int*)calloc(N,sizeof(int)) ;
	countIn = (int*)calloc(N,sizeof(int)) ;
	for(int i = 0 ; i < N ; i++)
	{
		NeighborsOut[i] = (int*)calloc(1,sizeof(int)) ;
		NeighborsIn[i] = (int*)calloc(1,sizeof(int)) ;
		NeighborsOut[i][0] = -100 ;
		NeighborsIn[i][0] = -100 ;
	}
	for(int i = 0 ; i < M ; i++)
	{//	printf("%d\n",i) ;
		NeighborsOut[Edges[i][0]][countOut[Edges[i][0]]] = Edges[i][1]  ;
		NeighborsIn[Edges[i][1]][countIn[Edges[i][1]]] = Edges[i][0] ;
		countOut[Edges[i][0]]++ ; 
		countIn[Edges[i][1]]++ ;
		//printf("%d %d %d\n",i , Edges[i][0],Edges[i][1]) ;
		NeighborsOut[Edges[i][0]] = (int*)realloc(NeighborsOut[Edges[i][0]],(countOut[Edges[i][0]]+1)*sizeof(int)) ;
		NeighborsIn[Edges[i][1]] = (int*)realloc(NeighborsIn[Edges[i][1]],(countIn[Edges[i][1]]+1)*sizeof(int)) ;
	} 
}
void Input_Fstar_and_Gstar(char **argv) 
{
	char *file_name1 , *file_name2 ;
	char f1[16] = "datasets/Fstar/" ;
	char f3[16] = "datasets/Gstar/" ;
	char f2[5] = ".txt" ;
	file_name1 = (char *)malloc((strlen(argv[1])+strlen(f1)+strlen(argv[5])+strlen(f2)+2)*sizeof(char)) ;
	file_name2 = (char *)malloc((strlen(argv[1])+strlen(f3)+strlen(argv[5])+strlen(f2)+2)*sizeof(char)) ;
	sprintf(file_name1,"%s%s_%s%s",f1,argv[1],argv[5],f2) ;
	sprintf(file_name2,"%s%s_%s%s",f3,argv[1],argv[5],f2) ;
	printf("filename1 is:%s\n",file_name1) ;
	FILE *F1, *F2 ;
	F1 = fopen(file_name1,"r") ;
	F2 = fopen(file_name2,"r") ;
	for(int i = 0 ; i < N ; i++)
	{
		fscanf(F1,"%f",&Fstar[i]) ;
		fscanf(F2,"%f",&Gstar[i]) ;
	}
	printf("filename1 is:%s\n",file_name1) ;

	deleta = atof(argv[2])  ;
	fclose(F1) ;
	fclose(F2) ;
	Pij_weight = (float**)calloc(N,sizeof(float*)) ;
	for(int i = 0 ; i < N ; i++)
	{	
		Pij_weight[i] = (float*)calloc(OutDegree[i],sizeof(float)) ;
	} 
}
void Input_Pij(char **argv){
	char *file_name1  ;
	char f1[10] = "datasets/" ;
	char f2[11] = "_edges.txt" ;
	file_name1 = (char*)malloc((strlen(argv[1])+strlen(f1)+strlen(f2)+1)*sizeof(char)) ;
	sprintf(file_name1,"%s%s%s",f1,argv[1],f2) ;
	FILE *F1 ;
	F1 = fopen(file_name1,"r") ;
	int index1 , index2 ;
	for(int i = 0 ; i < M ; i++)
	{
		fscanf(F1,"%d %d %f\n",&index1, &index2 ,&Pij[i]) ;
	}
	fclose(F1) ;
	
	
}
void Init_Influence_and_Susceptibility(char **argv)
{
	FILE *fin ;
	printf("argv3 is %s\n",argv[3]) ;
	if(strcmp(argv[3],"NO")==0){
		for(int i = 0 ; i < N ; i++)
		{
			Influence[i] = atof(argv[4]);
			Susceptibility[i] = atof(argv[4]) ;
		}
	}
	else
	{
		fin = fopen("datasets/I_S.txt","r") ;
		for(int i = 0 ; i < N ; i++){
			fscanf(fin,"%f %f\n",&Influence[i],&Susceptibility[i]) ;
			//printf("I and S are:%f %f\n",Influence[i],Susceptibility[i]) ;
		}
	}
}
void Iterate_Influence_label(int Iteration,char **argv)
{
	int i = 0  ,jj = 0;
	float value1 , value2 ;
	FILE *fout ;
	char *file_name ; 
	char f1[9] = "results/" ;
	char f2[4] = "I_S" ;
	char f3[5] = ".txt" ;
	file_name = (char*)calloc((strlen(f1)+strlen(f2)+strlen(f3)+strlen(argv[2])+strlen(argv[4])+1),sizeof(char)) ;
	sprintf(file_name,"%s%s_%s_%s%s",f1,f2,argv[2],argv[4],f3) ;
	while(jj < Iteration)
	{  
		#pragma omp parallel for 
		for(i = 0 ; i < N ; i++)
		{	
			Parallel_Influence_and_Susceptibility(i) ;
		}
		value2 = I_S_error(argv,jj) ;
		if(value2<0.01)
		{	
			break ;
		}
		else
		{
			Exchange_Influence_and_Susceptibility(Influence,Influence1) ;
			Exchange_Influence_and_Susceptibility(Susceptibility,Susceptibility1) ;
			jj++ ;
		}		
	}
	
}
void Parallel_Influence_and_Susceptibility(int i)
{
			float part1, part2 ,part3, part4 ;
			part2 = 0.0 ;	
			for(int j = 0 ; j < OutDegree[i] ; j++)
			{	part1 = 0.0 ;
			    if(NeighborsOut[i][j]>=0){
			    	for(int k = 0 ; k < InDegree[NeighborsOut[i][j]] ; k++)
					{
						part1 =  part1 + Influence[NeighborsIn[NeighborsOut[i][j]][k]] ;	
							
					}
					part1 = (float)Gstar[NeighborsOut[i][j]]/(part1+deleta) ;
					part2 = part2 + part1 ; 
				}					
			}
			Influence1[i] = (float)Fstar[i]/(part2+deleta) ;		
			part4 = 0.0 ;
			
			for(int j = 0 ; j < InDegree[i] ; j++)
			{	
			    part3 = 0.0 ;
			    if(NeighborsIn[i][j]>=0)
			    {
				    for(int k = 0 ; k < OutDegree[NeighborsIn[i][j]] ; k++)
					{
						part3 =  part3 + Susceptibility[NeighborsOut[NeighborsIn[i][j]][k]] ;
							
					}
					part3 = (float)Fstar[NeighborsIn[i][j]]/(part3+deleta) ;
					part4 = part4 + part3 ; 
				}	
			}
			Susceptibility1[i] = (float)Gstar[i]/(part4+deleta) ;
	
}
void Exchange_Influence_and_Susceptibility(int *arr1, int *arr2)
{
	for(int i = 0 ; i < N ; i++)
	{
		arr1[i] = arr2[i]  ;
	}	
}
float Generate_pij_weight()
{	
	float sum = 0 ;
	for(int i = 0 ; i < N ; i++)
	{	
		for(int j = 0 ; j < OutDegree[i] ; j++)
		{	
			Pij_weight[i][j] = Influence1[i]*Susceptibility1[NeighborsOut[i][j]] ;
		}
	}
	return sum ;
}
void Output(char **argv)
{	
	char f1[9] = "results/";
	char f2[5] = ".txt" ;
	char *file_name = (char *) malloc((strlen(argv[1])+strlen(argv[2])+strlen(f1)+strlen(argv[5])+strlen(argv[4])+strlen(f2)+4)*sizeof(char)) ;
	sprintf(file_name,"%s%s_%s_%s_%s%s",f1,argv[1],argv[2],argv[4],argv[5],f2) ;
	int i ; //output filename results/filename_deleta_initial_solution_randomremovededge.txt
	FILE *fout ;
	fout =fopen(file_name,"w") ; 
	for(int i = 0 ; i < N ; i++)
	{	
		for(int j = 0 ; j < OutDegree[i] ; j++)
		{
			Pij_weight[i][j] = Influence1[i]*Susceptibility1[NeighborsOut[i][j]] ;
			fprintf(fout,"%d %d %f\n",i,NeighborsOut[i][j],Pij_weight[i][j]) ;
		}
	} 
	fclose(fout);	
}
void Output_detaildata(char **argv)
{
	char f1[20] = "results/detaildata/";
	char f2[5] = ".txt" ;
	char *file_name = (char *) malloc((strlen(argv[1])+strlen(argv[2])+strlen(argv[4])+strlen(argv[5])+strlen(f1)+strlen(f2)+4)*sizeof(char)) ;
	sprintf(file_name,"%s%s_%s_%s_%s%s",f1,argv[1],argv[2],argv[4],argv[5],f2) ;
	int i  ;
	FILE *fout ;
	fout = fopen(file_name,"w") ;
	fprintf(fout,"%s %s %s %s\n","Inf_REC","Sus_REC","F_IC","G_IC") ;
	for(i = 0 ; i < N ; i++)
		fprintf(fout,"%f %f %f %f\n",Influence1[i],Susceptibility1[i],Fstar[i],Gstar[i]) ;
	fclose(fout)   ;	
}
void FreeMemory() 
{	
	free(InDegree) ;
	free(OutDegree) ;
	for(int i = 0 ; i < M ; i++)
	{
		free(Edges[i]) ;
	}
	free(Edges) ;
}
void Conserve_the_lasttime_I_and_S(){
	
	FILE *fout ;
	fout = fopen("results/I_S_N_N_1.txt","w") ;
	fprintf(fout,"%s %s %s %s\n","I_n","I_n+1","S_n","S_n+1") ;
	for(int i = 0 ; i < N ;i++)
	{
		
		fprintf(fout,"%f %f %f %f\n",Influence[i],Influence1[i],Susceptibility[i],Susceptibility1[i]) ;
	}
	fclose(fout) ;
} 
float I_S_error(char **argv,int iteration){
	FILE *fout ;
	char *file_name ; 
	char f1[11] = "Iteration/" ;
	char f2[4] = "I_S" ;
	char f3[5] = ".txt" ;
	file_name = (char*)calloc((strlen(f1)+strlen(f2)+strlen(f3)+strlen(argv[1])+strlen(argv[5])+3),sizeof(char)) ;
	sprintf(file_name,"%s%s_%s_%s%s",f1,f2,argv[1],argv[5],f3) ;
	float sum1 = 0.0 , sum2 = 0.0 ;
	for(int i = 0 ; i < N ; i++)
	{	
		sum1 = sum1 + fabs(Influence1[i]-Influence[i]) ;
		sum2 = sum2 + fabs(Susceptibility1[i]-Susceptibility[i]) ;
	}
	fout = fopen(file_name,"a+")  ;
	fprintf(fout,"%d %f %f\n",iteration,sum1,sum2) ;
	printf("error(I) is:%f ,error(S) is:%f\n",sum1, sum2) ;
	fclose(fout)  ;	
	return sum1>sum2?sum1:sum2 ;		
}

