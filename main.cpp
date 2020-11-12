#include <limits.h>
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int nV, numA, inf = INT_MAX, f = 0, r = -1;

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

Aresta* getAresta(Grafo** a, int numV, int numD) {
	Vertice* v = getVertice(a, numV);
	Aresta* ar = (v)->arestas;
	while (ar != NULL && ar->destino != numD)
		ar = ar->next;
	return ar;
}

Vertice* initVertice(int num, Vertice* next, Vertice* prev) {// Inicia um vértice 
	Vertice* v = (Vertice*)malloc(sizeof(Vertice));
	v->next = next;
	v->prev = prev;
	v->numero = num;
	return v;
};

void inserirVertice(Grafo** a, int num) // insere Vertice no Grafo
{
	Vertice* v;
	if ((*a)->vertices == NULL) {
		v = initVertice(num, NULL, NULL);

		(*a)->vertices = v;
	}
	else {
		Vertice* oldhead = (*a)->vertices;
		Vertice* newhead = initVertice(num, oldhead, NULL);
		v = newhead;
		oldhead->prev = newhead;
		(*a)->vertices = newhead;
	}
	v->arestas = NULL;
}

Aresta* initAresta(int destino, int origem, int peso, Aresta* next, Aresta* prev) { //inicia a uma Aresta
	Aresta* a = (Aresta*)malloc(sizeof(Aresta));
	a->origem = origem;
	a->destino = destino;
	a->peso = peso;
	a->next = next;
	a->prev = prev;

	return a;
};

int inserirAresta(Grafo** a, int origem, int destino, int peso) { //função que insere as arestas no grafo
	Vertice* v = getVertice(a, origem);
	if (v == NULL)
		return 1;
	if (getAresta(a, origem, destino) != NULL)
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
	Vertice* ax = getVertice(a, origem);

	if (v == ax)
		printf("asd");
	Vertice* v2 = getVertice(a, destino);
	if (v2 == NULL)
		return 1;

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
	return 0;
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
void printVertices(Grafo** g, bool printSupers) { //função que percorre os vértices e os imprime
	Vertice* vertice = (*g)->vertices;
	if (vertice != NULL) {
		Vertice* v = vertice;
		while (v != NULL) {
			printf("No");
			printf(" (%d) grau (%d) -> ", v->numero, v->grau);
			printArestas(v->arestas);
			printf("\n");
			v = v->next;
		}
	}
}

void printGrafo(Grafo** graf, bool printfSupers) { //função que imprime o grafo
	printf("\n\n");
	printVertices(graf, printfSupers);
	printf("\n\n");
}

void calculaGrau(Grafo** graf) //essa função calcula o grau de cada vértice do grafo, tem como parametro o grafo
{
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

int carrega(Grafo** graf, char* _arq1, char* _arq2) 
//essa função é utilizada para realizar o carregamento dos arquivos de
//entrada e saída. Recebe como parâmetros o grafo, um ponteiro para o arquivo de entrada e outro para o arquivo de saída
{
	FILE* f, * f2;
	int origem, destino, peso;
	f = fopen(_arq1, "r"); //arquivo de entrada
	f2 = fopen(_arq2, "w"); //arquivo de saída

	if (f == NULL) {
		printf("Erro na leitura do arquivo");
		return 0;
	};

	fscanf(f, "%d", &nV); //utilizado para ler do arquivo de entrada o número de vértices do grafo
	fprintf(f2, "%d\n", nV); //utilizado para escrever no arquivo de saída o número de vértices do grafo
	for (int i = 0; i < nV; i++) {
		inserirVertice(graf, i);
	}
	numA = 0;
	while (!feof(f)) {
		fscanf(f, "%d %d %d", &origem, &destino, &peso); //utilizado para ler do arquivo de entrada a origem, destino e peso de cada aresta
		Aresta* b = (Aresta*)malloc(sizeof(Aresta));
		b->peso = peso;
		b->origem = origem;
		b->destino = destino;
		b->next = nullptr;
		inserirAresta(graf, origem, b->destino, b->peso);
		Aresta* a = getAresta(graf, origem, destino);
	};
	calculaGrau(graf);
	fprintf(f2, "%d\n", (*graf)->nArestas); //utilizado para escrever no arquivo de saída o número de arestas do grafo
	fprintf(f2, "%f\n", grauMedioGrafo(graf)); //utilizado para escrever no arquivo de saída o grau médio do grafo
	for (int i = 0; i < nV; i++) //essa estrutura de repetição é utilizada para calcular a frequência relativa de cada grau do grafo
	{
		Vertice* vertice = (*graf)->vertices;
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
	return 1;
}
void caminhamentoLargura(int** matriz, int visitado[], int q[], int verticeInicial)
{
	for (int i = 0; i < nV; i++)
		if (matriz[verticeInicial][i] && !visitado[i])
			q[++r] = i;
	if (f <= r)
	{
		visitado[q[f]] = 1;
		caminhamentoLargura(matriz, visitado, q, q[f++]);
	}


}
void buscaLargura_Grafo(int ini, int* visitado, int** matrizAdj) {
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
void caminhamentoProfundidade(int** matrizAdj, int visitado[], int verticeInicial)
{
	visitado[verticeInicial] = 1;
	for (int i = 0; i < nV; i++) {
		if (matrizAdj[verticeInicial][i] && !visitado[i])
		{
			printf("%d->%d\n", verticeInicial, i);
			caminhamentoProfundidade(matrizAdj, visitado, i);

		}
	}
}

int distanciaMinima(int dist[], bool vetor_b[]) {
	int min = INT_MAX, indice_minimo;

	for (int v = 0; v < nV; v++)
		if (vetor_b[v] == false && dist[v] <= min)
			min = dist[v], indice_minimo = v;

	return indice_minimo;
}

void imprimeSolucao(int vetorDist[]) {
	// printf("Vertice Distancia da origem\n");
	for (int i = 0; i < nV; i++)
		printf(
			"\nCusto do vértice origem para o vértice [%d] = %d\n",
			i,
			vetorDist[i]);
}

void dijkstra(int** matrizAdj, int origem) {
	int* dist = (int*)calloc(sizeof(int), nV);
	bool* vetor_b = (bool*)calloc(sizeof(bool), nV);
	for (int i = 0; i < nV; i++) {
		dist[i] = INT_MAX, vetor_b[i] = false;
	}
	dist[origem] = 0;

	for (int count = 0; count < nV - 1; count++) {
		int u = distanciaMinima(dist, vetor_b);
		vetor_b[u] = true;
		for (int v = 0; v < nV; v++)
			if (!vetor_b[v] && matrizAdj[u][v] && dist[u] != INT_MAX &&
				dist[u] + matrizAdj[u][v] < dist[v]) {
				dist[v] = dist[u] + matrizAdj[u][v];
			}
	}

	imprimeSolucao(dist);
}

void floyd(int** matrizAdj) {
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
		printf(" %d ", a[i]->peso);

	printf("}\n");
}

void listVertice(Grafo** graf, Vertice** vertices) {

	int i = 0;
	Vertice* vertice = (*graf)->vertices;
	if (vertice != NULL) {
		Vertice* v = vertice;
		while(v!=NULL){
		vertices[i++] = v;
		v = v->next;
		}
	}


}
void swap(Vertice* v1, Vertice* v2)
{
	Vertice temp = *v1;
	*v1 = *v2;
	*v2 = temp;
}
void sortVertice(Vertice** vertices, int n)
{
	int i, j;

	for (i = 0; i < n - 1; i++)
		for (j = 0; j < n - i - 1; j++)
			if (vertices[j]->grau > vertices[j + 1]->grau)
				swap(vertices[j], vertices[j + 1]);
}
void printVertice(Vertice** v, int n) {
	int i;
	printf("\n{");
	for (i = 0; i < n ; i++)
		printf(" %d ", v[i]->grau);

	printf("}\n");
}

void teste(Grafo** graf) {
	Vertice** vertices = NULL;
	(vertices) = (Vertice**)malloc(sizeof(Vertice) * nV);
	listVertice(graf, vertices);
	printVertice(vertices, nV);
	sortVertice(vertices, nV);
	printVertice(vertices, nV);
}

void MSTprim(Grafo** graf) {

}
void MSTKruskal(Grafo** graf) {


}
void chama_caminhamentoLargura(int** matriz, int verticeOrigem)
{
	printf("\n\nVertice origem: %d \n", verticeOrigem);
	int* visitado = (int*)calloc(sizeof(int), nV);
	int* q = (int*)calloc(sizeof(int), nV);
	for (int j = 0; j < nV; j++)
	{
		visitado[j] = 0;
		q[j] = 0;
	}
	caminhamentoLargura(matriz, visitado, q, verticeOrigem);
	for (int k = 0; k < nV; k++) {
		if (visitado[k]) {
			printf("%d ", k);
		}
	}
}
void chama_caminhamentoProfundidade(int** matriz, int verticeOrigem)
{
	printf("\n\nVertice origem: %d \n", verticeOrigem);
	int* visitado = (int*)calloc(sizeof(int), nV);

	for (int j = 0; j < nV; j++)
	{
		visitado[j] = 0;
	}
	caminhamentoProfundidade(matriz, visitado, verticeOrigem);
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
	char arquivo1[] = "grafo_125.txt";
	char arquivo2[] = "saida.txt";

	Grafo* graf = (Grafo*)malloc(sizeof(Grafo));
	graf->vertices = NULL;
	graf->nArestas=0;

	carrega(&graf, arquivo1, arquivo2); //carrega o grafo

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
		printf("5 - Prim: Arvore Geradora \n");
		printf("6 - Kruskal: Arvore Geradora Minima (POR ENQUANTO USE ESSE PARA TESTAR A FUNÇÃO QUE IMPLEMENTA A ORDENACAO DE VERTICES PELO GRAU -TESTE-) \n");
		printf("7 - Sair \n");
		printf("Escolha sua opcao: ");
		scanf("%d", &op);
		switch (op) {
		case 1:
			chama_caminhamentoLargura(matriz, 0);
			break;
		case 2:
			chama_caminhamentoProfundidade(matriz, 0);
			break;
		case 3:
			printf("\nOrigem: %d \n", 0);
			dijkstra(matrizAdj, 0);
			break;
		case 4:
			floyd(matrizAdj);
			break;
		case 7:
			sair = 1;
			break;
		case 6:
			teste(&graf);
			break;
		default:
			printf("Valor invalido!\n");
		}
	} while (!sair);
	return 0;
}