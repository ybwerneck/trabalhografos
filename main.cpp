#include <limits.h>
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <queue>

int numA, inf = INT_MAX, f = 0, r = -1;

struct Aresta {//Node para a lista encadeadas de Arestas. Representa uma Aresta no grafo 
	int destino;
	int origem;
	int peso;

	struct Aresta* next;
	struct Aresta* prev;
};
struct Vertice { //Node para a lista encadeadas de Vertice. Representa um Vertice no grafo 
	int numero;// identificador
	int grau;
	struct Vertice* next;
	struct Vertice* prev;


	struct Aresta* arestas; //Lista encadeada contendo as arestas que partem desse vértice
};

struct Grafo { // Estrutura que representa um grafo
	struct Vertice* vertices; // lista encadeada de vertices nos grafos
	int nArestas = 0;  // para facilitar a obtenção dessa informação e nao precisar rodar o grafo inteiro para contar toda vez que for usa-la
	int nVertices=0;
	int* a;
};

typedef struct Grafo Grafo;
typedef struct Vertice Vertice;
typedef struct Aresta Aresta;

Vertice* getVertice(Grafo** a, int num) {
	Vertice* v = (*a)->vertices;
	while (v != NULL && v->numero != num)
		v = v->next;
	return v;
}

bool temVertice(Grafo** graf, int num) {
	return(getVertice(graf, num) != NULL);

}
bool temTodosVertices(Grafo** graf, Grafo** graf2) {
	Vertice* vertices = (*graf)->vertices;
	while (vertices != NULL)
	{
		if (temVertice(graf2, vertices->numero))
			return false;
		vertices = vertices->next;
	}
	return true;
}


Aresta* getAresta(Grafo** a, int numV, int numD) {
	Vertice* v = getVertice(a, numV);
	Aresta* ar = (v)->arestas;
	while (ar != NULL && ar->destino != numD)
		ar = ar->next;
	return ar;
}
bool temAresta(Grafo** a, int numV, int numD) {
	return !((getAresta(a, numV, numD) == NULL) && (getAresta(a, numV, numD) == NULL)) ;
}

Vertice* initVertice(int num, Vertice* next, Vertice* prev) {// Inicia um vértice 
	Vertice* v = (Vertice*)malloc(sizeof(Vertice));
	v->numero = num;
	return v;
};



Aresta* initAresta(int destino, int origem, int peso, Aresta* next, Aresta* prev) { //inicia a uma Aresta
	Aresta* a = (Aresta*)malloc(sizeof(Aresta));
	a->origem = origem;
	a->destino = destino;
	a->peso = peso;
	a->next = next;
	a->prev = prev;

	return a;
};
void inserirVertice(Grafo** a, int num) // insere Vertice no Grafo
{
	Vertice* v;
	if ((*a)->vertices == NULL) {
		v = initVertice(num, NULL, NULL);
		v->next = NULL;
		v->prev = NULL;
		(*a)->vertices = v;

	}
	else {
		Vertice* oldhead = (*a)->vertices;
		Vertice* newhead = initVertice(num, oldhead, NULL);
		v = newhead;
		oldhead->prev = newhead;
		newhead->next = oldhead;
		newhead->prev = NULL;
		(*a)->vertices = newhead;
	}
	v->arestas = NULL;
}


int inserirAresta(Grafo** a, int origem, int destino, int peso) { //função que insere as arestas no grafo
	Vertice* v = getVertice(a, origem);
	if (v == NULL)
		return 1;
	if (getAresta(a, origem, destino) != NULL)
		return 1;
	Vertice* v2 = getVertice(a, destino);
	if (v2 == NULL)
		return 1;
	if ((v)->arestas == NULL) {
		(v)->arestas = initAresta(destino, v->numero, peso, NULL, NULL);
	}
	else {
		Aresta* oldhead = (v)->arestas;
		Aresta* newhead = initAresta(destino, v->numero, peso, oldhead, NULL);
		initAresta(destino, v->numero, peso, oldhead, NULL);
		oldhead->prev = newhead;
		(v)->arestas = newhead;
	}

	if ((v2)->arestas == NULL)
		(v2)->arestas = initAresta(origem, v2->numero, peso, NULL, NULL);
	else {
		Aresta* oldhead = (v2)->arestas;
		Aresta* newhead = initAresta(origem, v2->numero, peso, oldhead, NULL);
		initAresta(origem, v2->numero, peso, oldhead, NULL);
		oldhead->prev = newhead;
		(v2)->arestas = newhead;
	}
	(*a)->nArestas++;
	(*a)->nArestas++;

	return 0;
}
int inserirArestaHard(Grafo** g, int origem, int destino, int peso) { //Inseri a aresta, se alguma das pontas não estiver no grafo, insere elas também
	Vertice* v = getVertice(g, origem);
	if (v == NULL)
		inserirVertice(g, origem);
	Vertice* v2 = getVertice(g, destino);
	if (v2 == NULL)
		inserirVertice(g, destino);
	return inserirAresta(g, origem, destino, peso);

}

Grafo* initGrafo(int n) {
	Grafo* g = (Grafo*)malloc(sizeof(Grafo));
	g->vertices = NULL;
	g->nArestas = 0;
	for (int i = 0; i < n; i++) {
		inserirVertice(&g, i);
	}
	g->nVertices = n;
	return g;

}
void printArestas(Aresta* aresta) { // função que percorre as arestas e as imprime
	Aresta* a = aresta;
	while (a != NULL) {
			printf("(No:%d Peso: ", a->destino);
			printf("%d", a->peso);
			printf(") ");
		

		a = a->next;
	}
}
void printVertices(Grafo** g, bool printGrau) { //função que percorre os vértices e os imprime
	Vertice* vertice = (*g)->vertices;
	if (vertice != NULL) {
		Vertice* v = vertice;
		while (v != NULL) {
			printf("No  (%d)",v->numero);		
			if(printGrau)printf("grau (%) ", v->numero, v->grau);
			printArestas(v->arestas);
			printf("\n");
			v = v->next;
		}
	}
}
void printGrafo(Grafo** graf, bool printfGrau) { //função que imprime o grafo
	printf("\n\n");
	printVertices(graf, printfGrau);
	printf("\n\n");
}
void calculaGrau(Grafo** graf) //essa função calcula o grau de cada vértice do grafo, tem como parametro o grafo
{
	int nV = (*graf)->nVertices;
	int** matriz = (int**)malloc(sizeof(int*) * nV); //é criada uma matriz através de alocação dinamica
	for (int i = 0; i < nV; i++) {
		matriz[i] = (int*)calloc(sizeof(int), nV);
	}

	//são percorridos todos os vértices e todas as arestas, os pesos são colocados em uma matriz 
	Vertice* vertice = (*graf)->vertices;
	if (vertice != NULL) {
		Vertice* v = vertice;
		while (v != NULL) {
			Aresta* a = v->arestas;
			while (a != NULL) {
					matriz[v->numero][a->destino] = a->peso;
				
				a = a->next; //vai para a próxima aresta
			}
			printf("\n");
			v = v->next; //vai para o próximo vértice
		}
	}
	//-------------
	Vertice* v = vertice;
	int cont = 0;
	for (int i = 0; i < nV; i++) //percorre a matriz que auxiliara na contagem do grau de cada vertice
	{
		cont = 0;
		for (int j = 0; j < nV; j++)
		{
			if (matriz[i][j] != 0) //se houver um peso, ou seja, uma aresta entre os vértices conta o grau
			{
				cont++; // utilizado para calcular o grau de cada vertice
			}

		}
		v->grau = cont; //atualiza o grau
		v = v->next; //vai pro próximo vértice
	}
}
float grauMedioGrafo(Grafo** graf) //essa função é utilizada para calcular o grau médio do grafo, ou seja, realiza a soma
//do grau de todos os vértices, e divide pela quantidade de vértices. A função recebe como parâmetro o grafo no qual se deseja
//saber o grau médio.
{
	int nV = (*graf)->nVertices;
	int somaGrau = 0;
	Vertice* vertice = (*graf)->vertices;
	if (vertice != NULL) {
		Vertice* v = vertice;
		while (v != NULL)
		{
			somaGrau = somaGrau + v->grau; //soma-se o grau de todos os vértices, até que todos os vértices tenham sido percorridos.
			v = v->next;
		}

	}
	float grauMedio = 0;
	grauMedio = float(somaGrau) / nV; //divide-se a soma total dos graus dos vértices pela quantidade de vértices para obter o grau médio.
	return grauMedio; //ao final, a função retorna esse valor (grau médio do grafo) para quem a chamou.
}
Grafo* carrega( char* _arq1, char* _arq2) 
//essa função é utilizada para realizar o carregamento dos arquivos de
//entrada e saída. Recebe como parâmetros o grafo, um ponteiro para o arquivo de entrada e outro para o arquivo de saída
{

	int nV;
	FILE* f, * f2;
	int origem, destino, peso;
	f = fopen(_arq1, "r"); //arquivo de entrada
	f2 = fopen(_arq2, "w"); //arquivo de saída

	if (f == NULL) {
		printf("Erro na leitura do arquivo");
		return 0;
	};

	fscanf(f, "%d", &nV); //utilizado para ler do arquivo de entrada o número de vértices do grafo
	Grafo* graf = initGrafo(nV);
	fprintf(f2, "%d\n", nV); //utilizado para escrever no arquivo de saída o número de vértices do grafo
    
	while (!feof(f)) {
		fscanf(f, "%d %d %d", &origem, &destino, &peso); //utilizado para ler do arquivo de entrada a origem, destino e peso de cada aresta
		Aresta* b = (Aresta*)malloc(sizeof(Aresta));
		b->peso = peso;
		b->origem = origem;
		b->destino = destino;
		b->next = nullptr;
		inserirAresta(&graf, origem, b->destino, b->peso);
		Aresta* a = getAresta(&graf, origem, destino);
	};
	calculaGrau(&graf);
	fprintf(f2, "%d\n", graf->nArestas); //utilizado para escrever no arquivo de saída o número de arestas do grafo
	fprintf(f2, "%f\n", grauMedioGrafo(&graf)); //utilizado para escrever no arquivo de saída o grau médio do grafo
	for (int i = 0; i < nV; i++) //essa estrutura de repetição é utilizada para calcular a frequência relativa de cada grau do grafo
	{
		Vertice* vertice = graf->vertices;
		Vertice* v = vertice;
		float frequencia = 0;
		int cont = 0;
		while (v != NULL) { //utilizado para percorrer todos os vértices
			if (v->grau == i) // verifica se o grau do vértice atual é o mesmo do qual está fazendo a contagem
			{
				cont++; //utilizado para fazer o somatório de quantos vértices possuem grau "d"
			}
			v = v->next;
		}
		frequencia = float(cont) / nV; //o cálculo é realizado pela divisão da quantidade de vértices com grau "d" pela quantidade de vértices do grafo
		fprintf(f2, "%f\n", frequencia); //utilizado para escrever no arquivo de saída a frequência relativa de cada grau "d" do grafo
	}

	fclose(f);
	fclose(f2);
	return graf;
}
void caminhamentoLargura(int** matriz, int visitado[], int q[], int verticeInicial,int nV)
{
	for (int i = 0; i < nV; i++)
		if (matriz[verticeInicial][i] && !visitado[i])
			q[++r] = i;
	if (f <= r)
	{
		visitado[q[f]] = 1;
		caminhamentoLargura(matriz, visitado, q, q[f++],nV);
	}


}
void buscaLargura_Grafo(int ini, int* visitado, int** matrizAdj,int nV) {
	int i, vert, cont = 1;
	int* fila, IF = 0, FF = 0;
	for (i = 0; i < nV; i++)
		visitado[i] = 0;
	fila = (int*)malloc(nV * sizeof(int));
	FF++;
	fila[FF] = ini;
	visitado[ini] = cont;
	while (IF != FF) {
		IF = (IF + 1) % nV;
		vert = fila[IF];
		cont++;
		for (i = 0; i < nV; i++) {
			if (!visitado[matrizAdj[vert][i]]) {
				FF = (FF + 1) % nV;
				fila[FF] = matrizAdj[vert][i];
				visitado[matrizAdj[vert][i]] = cont;
			}
		}
	}
	free(fila);
	for (i = 0; i < nV; i++)
		printf("%d -> %d\n", i, visitado[i]);
}
void caminhamentoProfundidade(int** matrizAdj, int visitado[], int verticeInicial,int nV)
{
	visitado[verticeInicial] = 1;
	for (int i = 0; i < nV; i++) {
		if (matrizAdj[verticeInicial][i] && !visitado[i])
		{
			printf("%d->%d\n", verticeInicial, i);
			caminhamentoProfundidade(matrizAdj, visitado, i,nV);

		}
	}
}
int distanciaMinima(int dist[], bool vetor_b[],int nV) {
	int min = INT_MAX, indice_minimo;

	for (int v = 0; v < nV; v++)
		if (vetor_b[v] == false && dist[v] <= min)
			min = dist[v], indice_minimo = v;

	return indice_minimo;
}
void imprimeSolucao(int vetorDist[],int nV) {
	// printf("Vertice Distancia da origem\n");
	for (int i = 0; i < nV; i++)
		printf(
			"\nCusto do vértice origem para o vértice [%d] = %d\n",
			i,
			vetorDist[i]);
}
void dijkstra(int** matrizAdj, int origem,int nV) {
	int* dist = (int*)calloc(sizeof(int), nV);
	bool* vetor_b = (bool*)calloc(sizeof(bool), nV);
	for (int i = 0; i < nV; i++) {
		dist[i] = INT_MAX, vetor_b[i] = false;
	}
	dist[origem] = 0;

	for (int count = 0; count < nV - 1; count++) {
		int u = distanciaMinima(dist, vetor_b,nV);
		vetor_b[u] = true;
		for (int v = 0; v < nV; v++)
			if (!vetor_b[v] && matrizAdj[u][v] && dist[u] != INT_MAX &&
				dist[u] + matrizAdj[u][v] < dist[v]) {
				dist[v] = dist[u] + matrizAdj[u][v];
			}
	}

	imprimeSolucao(dist,nV);
}
void floyd(int** matrizAdj, int nV) {
	
	int** matrizDist = (int**)malloc(sizeof(int*) * nV);
	for (int i = 0; i < nV; i++) {
		matrizDist[i] = (int*)calloc(sizeof(int), nV);
	}
	for (int i = 0; i < nV; i++)
		for (int j = 0; j < nV; j++)
			matrizDist[i][j] = matrizAdj[i][j];
	for (int k = 0; k < nV; k++) {
		for (int i = 0; i < nV; i++) {
			for (int j = 0; j < nV; j++) {
				if (matrizDist[i][k] + matrizDist[k][j] < matrizDist[i][j]) {
					matrizDist[i][j] = matrizDist[i][k] + matrizDist[k][j];
				}
			}
		}
	}
	for (int i = 0; i < nV; i++) {
		for (int j = 0; j < nV; j++) {
			printf(
				"\nCusto do vértice [%d] para o vértice [%d] = %d\n",
				i,
				j,
				matrizDist[i][j]);
		}
	}
}
void listAresta(Grafo** graf, Aresta** arestas) {

	int i=0;
	Vertice* vertice = (*graf)->vertices;
	if (vertice != NULL) {
		Vertice* v = vertice;
		while (v != NULL) {
			Aresta* a = v->arestas;
			while (a != NULL) {
				arestas[i++] = a;
				a = a->next;
			}
			v = v->next;
		}
	}

}
void swap(Aresta* a1, Aresta* a2)
{
	Aresta temp = *a1;
	*a1 = *a2;
	*a2 = temp;
}
void sortAresta(Aresta** arestas, int n)
{
	int i, j;

	for (i = 0; i < n - 1; i++)
		for (j = 0; j < n - i - 1; j++)
			if (arestas[j]->peso > arestas[j + 1]->peso)
				swap(arestas[j], arestas[j + 1]);
}
void printAresta(Aresta** a, int n) {
	int i;
	printf("\n{");
	for (i = 0; i < n ; i++)
		printf(" %d->%d (%d)  ", a[i]->origem, a[i]->destino,a[i]->peso);

	printf("}\n");
}
void listVertice(Grafo* graf, Vertice** vertices) {

	int i = 0;
	Vertice* vertice = (graf)->vertices;
	if (vertice != NULL) {
		Vertice* v = vertice;
		while(v!=NULL){
		vertices[i++] = v;
		v = v->next;
		}
	}


}
void printVertice(Vertice** v, int n) {
	int i;
	printf("\n{");
	for (i = 0; i < n; i++)
		printf(" %d (%d) ", v[i]->numero, v[i]->grau);

	printf("}\n");
}
void sortVertice(Vertice** vertices, int n)
{
	int i, j;
	for (i = 0; i < n - 1; i++)
		for (j = 0; j < n - i - 1; j++)
			if (vertices[j]->grau < vertices[j + 1]->grau)
			{
				Vertice* temp = vertices[j];
				vertices[j] = vertices[j + 1];
				vertices[j + 1] = temp;
			}
}
void MSTprim(Grafo** graf) {

		Grafo* grafoAux = initGrafo(0);
		int nV = (*graf)->nVertices;
		/*int nA = (*graf)->nArestas;
		Vertice** vertices = NULL;
		(vertices) = (Vertice**)malloc(sizeof(Vertice) * nV);
		listVertice(graf, vertices);
		sortVertice(vertices, nV);
		Aresta** arestas = NULL;
		(arestas) = (Aresta**)malloc(sizeof(Aresta) * nA);
		listAresta(graf, arestas);
		sortAresta(arestas, nA);*/
		int valorCaminho, totalCaminho = 0;


		for (int i = 0; i < nV; i++) {
			valorCaminho = 0;
			for (int j = 0; j < nV; j++) {

				if (i == j)
					j++;
				Aresta* ar = getAresta(graf, i, j);
				if (valorCaminho == 0) {
					valorCaminho = ar->peso;
				}
				else {
					if (valorCaminho > ar->peso)
						valorCaminho = ar->peso;
				}
			}
			totalCaminho = totalCaminho + valorCaminho;
		}

		printf("%d", totalCaminho);
		//printGrafo(&grafoAux, false);

	
}
bool visitado(Vertice** visitadas,int k,Vertice* vertice){
	for (int i = 0; i < k; i++) {
		if (visitadas[i]->numero == vertice->numero)
			return true;
	}
	return false;
}
bool hasPath(Grafo** graf,int origem,int destino) {
	int k=0;
	Vertice* v=getVertice(graf,origem);
	Vertice** visitadas = (Vertice**)malloc((sizeof(Vertice) *(*graf)->nVertices));
	std::queue <Vertice*> vertices;
	vertices.push(v);
	visitadas[k++] = v;
	while (!vertices.empty())
	{
		if (v->numero == destino)
			return true;
		v=vertices.front();
		vertices.pop();
		Aresta* a = v->arestas;
		while (a != NULL)
		{
			Vertice* aux = getVertice(graf, a->destino);
			if (!visitado(visitadas, k, aux)) {
			vertices.push(getVertice(graf, a->destino));
			visitadas[k++] = aux;
			}
			a = a->next;
	}

	}
	return false;
}
void MSTKruskal(Grafo** graf) { 
	Grafo* grafoAux = initGrafo(5);			//inicializa a solução sem arestas
	int nA = (*graf)->nArestas;
	Aresta** arestas = NULL;
	(arestas) = (Aresta**)malloc(sizeof(Aresta*) * nA);
	listAresta(graf, arestas);  //cria um conjunto com as arestas do grafo
	sortAresta(arestas,nA);		//ordena as aretas do grafo
	
	printAresta(arestas, nA);
	for (int i = 0; i < nA; i++) { //para cada elemento do conjunto começando pelo de menor peso
		if ((!hasPath(&grafoAux, arestas[i]->origem, arestas[i]->destino))&&
			!hasPath(&grafoAux, arestas[i]->destino, arestas[i]->origem)) //se a aresta liga duas arvores que ainda nao estão ligadas
		{	
			inserirAresta(&grafoAux, arestas[i]->origem, arestas[i]->destino, arestas[i]->peso); //liga as duas
			printf("\nnao tem caminho de %d %d \n", arestas[i]->origem, arestas[i]->destino);

		}
		else{
			printf("\nja tem caminho de %d %d \n", arestas[i]->origem, arestas[i]->destino);
		}
		}
	printGrafo(&grafoAux,false);
}

Grafo* grafoVizinho(Grafo** graf,Grafo** GrafoReferencia) {
	Grafo* grafoaux = initGrafo(0);
	Aresta** arestas;
	int nA = (*GrafoReferencia)->nArestas;
	(arestas) = (Aresta**)malloc(sizeof(Aresta*) * nA);
	listAresta(GrafoReferencia, arestas);  //cria um conjunto com as arestas do grafo

	for (int i = 0; i < (*GrafoReferencia)->nArestas; i++)
		if(temVertice(graf,arestas[i]->origem))
		if (!temVertice(graf, arestas[i]->destino))
			inserirVertice(&grafoaux, arestas[i]->destino);
		
	return grafoaux;

}
bool eDominante(Grafo** graf,Grafo** sol) {
	Grafo* grafoV;
	grafoV = grafoVizinho(sol,graf);
	for (Vertice* vertices = (*graf)->vertices;vertices != NULL;vertices=vertices->next)
	{	if (temVertice(sol, vertices->numero))
			continue;
		if (temVertice(&grafoV, vertices->numero))
			continue;
		return false;
	}
	return true;
}
void Guloso(Grafo** graf) {
	Grafo* grafoV;
	int nV = (*graf)->nVertices;
	Vertice** vertices =(Vertice**) malloc(sizeof(Vertice) * nV);
	int k = 0;

	listVertice(*graf, vertices);
	sortVertice(vertices, nV);
	grafoV=initGrafo(0);
	
	while (!eDominante(graf, &grafoV))
	{
		Vertice* v = vertices[k++];
		Grafo* gV = grafoVizinho(&grafoV,graf);
		if(!temVertice(&grafoV,v->numero))
		if(!temVertice(&gV,v->numero)){
			inserirVertice(&grafoV,v->numero);
			
			
			
		}
	}
	printf("O conjunto gerador minimo é ");
	printGrafo(&grafoV,false);

}

void chama_caminhamentoLargura(int** matriz, int verticeOrigem,int nV)
{
	printf("\n\nVertice origem: %d \n", verticeOrigem);
	int* visitado = (int*)calloc(sizeof(int), nV);
	int* q = (int*)calloc(sizeof(int), nV);
	for (int j = 0; j < nV; j++)
	{
		visitado[j] = 0;
		q[j] = 0;
	}
	caminhamentoLargura(matriz, visitado, q, verticeOrigem,nV);
	for (int k = 0; k < nV; k++) {
		if (visitado[k]) {
			printf("%d ", k);
		}
	}
}
void chama_caminhamentoProfundidade(int** matriz, int verticeOrigem,int nV)
{
	printf("\n\nVertice origem: %d \n", verticeOrigem);
	int* visitado = (int*)calloc(sizeof(int), nV);

	for (int j = 0; j < nV; j++)
	{
		visitado[j] = 0;
	}
	caminhamentoProfundidade(matriz, visitado, verticeOrigem,nV);
	int cont = 0;
	for (int k = 0; k < nV; k++)
	{
		if (visitado[k])
			cont++;
	}
	if (cont == nV)
		printf("O grafo é conexo.");
	else
		printf("O grafo não é conexo.");
}
//int main(int argc, char **argv)
int main()
{
	//as duas linhas abaixo foram utilizadas para garantir o padrão para a execução a ser utilizado pelo professor
	//char * arquivo1 = argv[1];
	//char * arquivo2 = argv[2];
	char arquivo1[] = "dados.txt";
	char arquivo2[] = "saida.txt";

	Grafo* graf = carrega(arquivo1, arquivo2); //carrega o grafo
	int nV = graf->nVertices;


	int** matrizAdj = (int**)malloc(sizeof(int*) * nV);
	for (int i = 0; i < nV; i++) {
		matrizAdj[i] = (int*)calloc(sizeof(int), nV);
	}

	int** matriz = (int**)malloc(sizeof(int*) * nV);
	for (int i = 0; i < nV; i++) {
		matriz[i] = (int*)calloc(sizeof(int), nV);
	}
	// Gera a matrix de adjacencia usada em algumas operações 
	Vertice* vertice = graf->vertices;
	if (vertice != NULL) {
		Vertice* v = vertice;
		while (v != NULL) {
			Aresta* a = v->arestas;
			while (a != NULL) {
					matrizAdj[v->numero][a->destino] = a->peso;
					matriz[v->numero][a->destino] = 1;
					matriz[v->numero][v->numero] = 1;
				
				a = a->next;
			}
			printf("\n");
			v = v->next;
		}
	}
	printf("-----------TRABALHO DE GRAFOS----------------- ");
	printf("\n-----------GRUPO 3---------------------------- ");
	printf("\nAna Beatriz Simoes Goncalves");
	printf("\nIgor Coimbra Vargas Lorenzeto");
	printf("\nLuisa Silva Ribeiro");
	printf("\nMatheus Fajardo Galvao");
	printf("\nYan Barbosa Werneck");
	printf("\n\n----------GRAFO---------------------------- ");
	printGrafo(&graf, true);
	printf("\n ----------MENU-----------------------------\n");
	int op = 0, sair = 0;
	do {
		printf("\n1 - Caminhamento em largura\n");
		printf("2 - Caminhamento em profundidade \n");
		printf("3 - Dijkstra: caminho minimo \n");
		printf("4 - Floyd: caminho minimo \n");
		printf("5 - Prim: Arvore Geradora\n");
		printf("6 - Kruskal: Arvore Geradora Minima \n");
		printf("7 - Sair \n");
		printf("Escolha sua opcao: ");
		scanf("%d", &op);
		switch (op) {
		case 1:
			chama_caminhamentoLargura(matriz, 0, graf->nVertices);
			break;
		case 2:
			chama_caminhamentoProfundidade(matriz, 0, graf->nVertices);
			break;
		case 3:
			printf("\nOrigem: %d \n", 0);
			dijkstra(matrizAdj, 0, graf->nVertices);
			break;
		case 4:
			floyd(matrizAdj,graf->nVertices);
			break;
		case 7:
			Guloso(&graf);
			break;
		case 6:
			MSTKruskal(&graf);
			break;
		default:
			printf("Valor invalido!\n");
		}
	} while (!sair);
	return 0;
}