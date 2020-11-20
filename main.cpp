#include <limits.h>
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <queue>

int numA, inf = INT_MAX;
int r = -1, f = 0;
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
	int nVertices = 0;
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
	return !((getAresta(a, numV, numD) == NULL) && (getAresta(a, numV, numD) == NULL));
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
	printf("{");
	while (a != NULL) {
		printf(" %d(P=%d) ", a->destino, a->peso);


		a = a->next;
	}printf("}");
}
void printVertices(Grafo** g, bool printGrau) { //função que percorre os vértices e os imprime
	Vertice* vertice = (*g)->vertices;
	if (vertice != NULL) {
		Vertice* v = vertice;
		while (v != NULL) {
			printf("\n");
			printf("No  %d", v->numero);
			if (printGrau)printf(" G=%d  ", v->grau);
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

void calculaGrau(Grafo** grafo) { //calcula grau do grafo
	Vertice* v = (*grafo)->vertices;
	while (v != NULL) {
		int k = 0;
		Aresta* a = v->arestas;
		while (a != NULL) {
			if (a->destino == v->numero)
				k++; //self loops são contados duas vezes
			k++;
			a = a->next;
		}
		v->grau = k;
		v = v->next;
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
Grafo* carrega(char* _arq1, char* _arq2)
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
		inserirAresta(&graf, origem, destino, peso);
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
}int distanciaMinima(int dist[], bool vetor_b[], int nV) {//essa função encontra o vértice do conjunto de vértices que ainda não está na solução do caminho mais curto e possui o menor valor de distância 
	int min = INT_MAX, indice_minimo;

	for (int v = 0; v < nV; v++)
		if (vetor_b[v] == false && dist[v] <= min) //verifica se o vértice não pertence ao caminho mais curto(vetor_b[v] == false) e se sua distância é menor que o mínimo
			min = dist[v], indice_minimo = v; //caso positivo, min recebe dist[v] e o índice minimo passa a ser o vértice v

	return indice_minimo;
}
void imprimeSolucao(int vetorDist[], int nV) {//essa função imprime a solução do algoritmo de dijkstra
	for (int i = 0; i < nV; i++)
		printf(
			"\nCusto do vértice origem para o vértice [%d] = %d\n",
			i,
			vetorDist[i]);
}
void dijkstra(int** matriz, int origem, int nV) {//o algoritmo de dijkstra recebe como parâmetro a matriz de adjacência do grafo, o vértice de origem e o número total de vértices

	int* dist = (int*)calloc(sizeof(int), nV); //alocando memória para dist 
	bool* vetor_b = (bool*)calloc(sizeof(bool), nV); //alocando memória para vetor_b
	for (int i = 0; i < nV; i++) {
		dist[i] = INT_MAX, vetor_b[i] = false;//inicializando o vetor de distância como infinito em todas as posições e vetor_b com false
	}
	dist[origem] = 0; //a distância do vértice de origem é sempre 0

	for (int count = 0; count < nV - 1; count++) {
		int u = distanciaMinima(dist, vetor_b, nV);//escolhendo o vértice de menor distância dentre os que ainda não foram percorridos
		vetor_b[u] = true;//atualizando o estado do vértice escolhido como visitado

		for (int v = 0; v < nV; v++)
			if (!vetor_b[v] && matriz[u][v] && dist[u] != INT_MAX && dist[u] + matriz[u][v] < dist[v]) {//se v não estiver em vetor_b, se existir uma aresta de u para v e se o peso do caminho do vértice de origem até v através de u for menor que o valor atual de dist[v], atualiza dist[v]
				dist[v] = dist[u] + matriz[u][v];
			}
	}

	imprimeSolucao(dist, nV);
}
void floyd(int** matriz, int nV) {//o algoritmo de floyd recebe como parâmetro a matriz de adjacência do grafo e o número de vértices

	int** matrizDist = (int**)malloc(sizeof(int*) * nV); //matrizDist é a matriz que possuirá a menor distância entre cada par de vértices

	for (int i = 0; i < nV; i++) {
		matrizDist[i] = (int*)calloc(sizeof(int), nV);//fazendo alocação de memória
	}
	for (int i = 0; i < nV; i++)
		for (int j = 0; j < nV; j++)
			matrizDist[i][j] = matriz[i][j];  //inicializando a matriz de distâncias 


	for (int k = 0; k < nV; k++) {//adicionando ao conjunto de vértices intermediários cada um dos vértices
		for (int i = 0; i < nV; i++) {//utilizando todos os vértices como origem
			for (int j = 0; j < nV; j++) {//para cada vértice escolhido como origem acima, escolhe-se todos os outros como destino

				if (matrizDist[i][k] + matrizDist[k][j] < matrizDist[i][j]) {//se o vértice k pertencer ao caminho mais curto de i para j, o valor de matrizDist[i][j] é atualizado
					matrizDist[i][j] = matrizDist[i][k] + matrizDist[k][j];
				}
			}
		}
	}
	//imprimindo a solução
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

	//Gera um array com todas arestas do grafo
  //Faciita na hora de realizar operacoes, como ordenacao
  int i = 0;
	Vertice* vertice = (*graf)->vertices;
	if (vertice != NULL) {
		Vertice* v = vertice;
		while (v != NULL) {
			Aresta* a = v->arestas;
			while (a != NULL) {
				Aresta* aux = (Aresta*)malloc(sizeof(Aresta));
				memcpy(aux, a, sizeof(Aresta));
				arestas[i++] = aux;
        a = a->next;
        //O array gerado contém ponteiros para copias das arestas originais e nao para os objetos Aresta originais em si
        //Isso foi feito para que o objeto possa ser manipulado sem modificar o grafo original

			}
			v = v->next;
		}
	}

}
void swap(Aresta* a1, Aresta* a2) // realiza a troca
{
	Aresta temp = *a1;
	*a1 = *a2;
	*a2 = temp;
}
void sortAresta(Aresta** arestas, int n) //ordena um array de  arestas do grafo
{
	int i, j;

	for (i = 0; i < n - 1; i++)
		for (j = 0; j < n - i - 1; j++)
			if (arestas[j]->peso > arestas[j + 1]->peso)
				swap(arestas[j], arestas[j + 1]);
}
void printAresta(Aresta** a, int n) { //imprime as arestas
	int i;
	printf("\n{");
	for (i = 0; i < n; i++)
		printf(" %d->%d (%d)  ", a[i]->origem, a[i]->destino, a[i]->peso);

	printf("}\n");
}
void printVertice(Vertice** v, int n) { //imprime os vértices com seus respectivos graus
	int i;
	printf("\n{");
	for (i = 0; i < n; i++)
		printf(" %d (%d) ", v[i]->numero, v[i]->grau);

	printf("}\n");
}
void sortVertice(Vertice** vertices, int n)//ordena os vértices do grafo
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
	//função que calcula a árvore geradora de custo mínimo utilizando o algoritimo de Prim.

	int nV = (*graf)->nVertices;//numero de vertices do grafo.

	int valorCaminho, totalCaminho = 0, cont = 0, i = 0, k = 0;
	bool* visitado = (bool*)calloc(sizeof(bool), nV);//vetor pra controle dos vertices que ja foram visitados.

	for (int vis = 1; vis < nV; vis++) {// preenche o vetor visitado com falso.
		visitado[vis] = false;
	}
	visitado[0] = true;// define o indice 0 do visitado como verdadeiro ja que o caminho começa nele.

	while (cont < nV) {// vai gerar um loop enquanto nao percorrer todos os vertices.
		valorCaminho = 0;
		for (int j = 0; j < nV; j++) {//percorre todos os vertices do grafo.
			if (visitado[j]) {// verifica se o vertice ja foi visitado, se ja foi vai para o proximo vertice.
				j++;
			}
			Aresta* ar = getAresta(graf, i, j);// cria uma aresta e define outra aresta com o i e o j como vertices.
			if (ar != NULL) {//verifica se a aresta é valida.
				if (valorCaminho == 0) {//verifica se o valor do caminho ja foi anteriormente instanciado.
					valorCaminho = ar->peso;// se nao foi instanciado define o valor do caminho com o peso da aresta
				}
				else {
					if (valorCaminho > ar->peso) {// se o valor ja foi instanciado e o peso da aresta for menor que o valor anterior define o novo valor com o peso da aresta
						valorCaminho = ar->peso;
						k = j;// k terá o valor da aresta que possui o menor peso;
					}
				}
			}
		}
		i = k;
		visitado[k] = true;// o vertice k para para o valor visitado.
		totalCaminho = totalCaminho + valorCaminho;// o total do caminho é adicionado do valor atual encontrado.
		cont++;// adiciona 1 ao contador
	}
	free(visitado);

	printf("Custo minimo para percorrer o grafo: %d", totalCaminho);//imprime o  custo mínimo para percorrer o grafo.

}
bool visitado(Vertice** visitadas, int k, Vertice* vertice) {
	for (int i = 0; i < k; i++) {
		if (visitadas[i]->numero == vertice->numero)
			return true;
	}
	return false;
}
bool hasPath(Grafo** graf, int origem, int destino) { // verifica se a aresta liga duas arvores que ainda nao estão ligadas
	int k = 0;
	Vertice* v = getVertice(graf, origem);

  //Faz isso fazendo uma busca de largura a partir do vertice de origem tentando achar o vertice de destino

	Vertice** visitadas = (Vertice**)malloc((sizeof(Vertice) * (*graf)->nVertices));//vertices visitados
	std::queue <Vertice*> vertices;//fila de vertices, ordem da busca por profundidade
	vertices.push(v);
	visitadas[k++] = v; 
	while (!vertices.empty()) //busca em profundiadade
	{
		if (v->numero == destino)
			return true; //tem caminho
		v = vertices.front();//pega o proximo vertice na primeira posicao fila
		vertices.pop();//tira ele da fila
		Aresta* a = v->arestas;
		while (a != NULL)
		{
			Vertice* aux = getVertice(graf, a->destino);
			if (!visitado(visitadas, k, aux)) {
				vertices.push(getVertice(graf, a->destino)); //adiciona todo os seus vizinho no final da fila
				visitadas[k++] = aux;
			}
			a = a->next;
		}

	}
  //percorreu todos os caminhos a partir da origem e nao achou o destino
	free(visitadas);
	return false; //nao tem caminho
}
void MSTKruskal(Grafo** graf) {
	Grafo* grafoAux = initGrafo((*graf)->nVertices);			//inicializa a solução sem arestas
	int nA = (*graf)->nArestas;
	Aresta** arestas = NULL;
	(arestas) = (Aresta**)malloc(sizeof(Aresta*) * nA);
	listAresta(graf, arestas);  //cria um conjunto com as arestas do grafo
	sortAresta(arestas, nA);		//ordena as arestas do grafo

	for (int i = 0; i < nA; i++) { //para cada elemento do conjunto começando pelo de menor peso
		if ((!hasPath(&grafoAux, arestas[i]->origem, arestas[i]->destino))) //se a aresta liga duas arvores que ainda nao estão ligadas
		{
			inserirAresta(&grafoAux, arestas[i]->origem, arestas[i]->destino, arestas[i]->peso); //liga as duas

		}
	}
	calculaGrau(&grafoAux);
	printGrafo(&grafoAux, true);
	free(arestas);
}

Grafo* grafoVizinho(Grafo** graf, Grafo** GrafoReferencia) {
	Grafo* grafoaux = initGrafo(0);
	Aresta** arestas;
	int nA = (*GrafoReferencia)->nArestas;
	(arestas) = (Aresta**)malloc(sizeof(Aresta*) * nA);
	listAresta(GrafoReferencia, arestas);  //cria um conjunto com as arestas do grafo

	for (int i = 0; i < (*GrafoReferencia)->nArestas; i++)
		if (temVertice(graf, arestas[i]->origem))
			if (!temVertice(graf, arestas[i]->destino))
				inserirVertice(&grafoaux, arestas[i]->destino);

	return grafoaux;

}

void caminhamentoLargura(int** matriz, int visitado[], int q[], int verticeInicial, int nV)
{
	for (int i = 0; i < nV; i++)
		if (matriz[verticeInicial][i] && !visitado[i])//vertifica se existe a aresta através da matriz de adjacencia e se o nó não foi visitado
			q[++r] = i; //coloca no vetor q o indice do vértice (r variável global)
	if (f <= r) //verifica  de f é menor ou igual a r (f e r são variáveis globais) para realizar a visitação do vértice
	{
		visitado[q[f]] = 1; //marca o vértice como visitado
		caminhamentoLargura(matriz, visitado, q, q[f++], nV); //recursão
	}


}

void caminhamentoProfundidade(int** matrizAdj, int visitado[], int verticeInicial, int nV)
{
	visitado[verticeInicial] = 1; //marca o vertice inicial como visitado
	for (int i = 0; i < nV; i++) { //percorre todos os vértices
		if (matrizAdj[verticeInicial][i] && !visitado[i]) //vertifica se existe a aresta através da matriz de adjacencia e se o nó não foi visitado, se sim visita o nó
		{
			printf("%d->%d\n", verticeInicial, i);
			caminhamentoProfundidade(matrizAdj, visitado, i, nV); //recursividade

		}
	}
}

void chama_caminhamentoLargura(int** matriz, int verticeOrigem, int nV) //função utilizada para chamar a função caminhamentoLargura
{
	printf("\n\nVertice origem: %d \n", verticeOrigem);
	int* visitado = (int*)malloc(sizeof(int)* nV);
	int* q = (int*)malloc(sizeof(int)* nV*nV);
	for (int j = 0; j < nV; j++)
	{
		visitado[j] = 0; //marca todos os nós como não visitado
		q[j] = 0; //vetor auxiliar
	}
	caminhamentoLargura(matriz, visitado, q, verticeOrigem, nV);
	for (int k = 0; k < nV; k++) {
		if (visitado[k]) {
			printf("%d ", k); //imprime nós visitados 
		}
	}
}
bool chama_caminhamentoProfundidade(int** matriz, int verticeOrigem, int nV) //função utilizada para chamar a função caminhamentoProfundidade
{
	printf("\n\nVertice origem: %d \n", verticeOrigem);
	int* visitado = (int*)calloc(sizeof(int), nV);

	for (int j = 0; j < nV; j++)
	{
		visitado[j] = 0;
	}
	caminhamentoProfundidade(matriz, visitado, verticeOrigem, nV);
	int cont = 0;
	for (int k = 0; k < nV; k++)
	{
		if (visitado[k])
			cont++;
	}
	if (cont == nV) //vertifica se o grafo é conexo
		{
    printf("O grafo é conexo.");
    return true;
    }
	else
		printf("O grafo não é conexo.");
    return false;
}
//int main(int argc, char **argv)
int main()
{
	//as duas linhas abaixo foram utilizadas para garantir o padrão para a execução a ser utilizado pelo professor
	//char * arquivo1 = argv[1];
	//char * arquivo2 = argv[2];
	char arquivo1[] = "grafo_125.txt";
	char arquivo2[] = "saida.txt";

	Grafo* graf = carrega(arquivo1, arquivo2); //carrega o grafo
	int nV = graf->nVertices;


	int** matriz = (int**)malloc(sizeof(int*) * nV);
	for (int i = 0; i < nV; i++) {
		matriz[i] = (int*)calloc(sizeof(int), nV);
	}

	int** matrizAdj = (int**)malloc(sizeof(int*) * nV);
	for (int i = 0; i < nV; i++) {
		matrizAdj[i] = (int*)calloc(sizeof(int), nV);
	}
	// Gera a matrizes usadas em algumas operações 
	Vertice* vertice = graf->vertices;
	if (vertice != NULL) {
		Vertice* v = vertice;
		while (v != NULL) {
			Aresta* a = v->arestas;
			while (a != NULL) {
				matriz[v->numero][a->destino] = a->peso;
				matrizAdj[v->numero][a->destino] = 1;
				matrizAdj[v->numero][v->numero] = 1;

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
			chama_caminhamentoLargura(matrizAdj, 0, graf->nVertices);
			break;
		case 2:
			chama_caminhamentoProfundidade(matrizAdj, 0, graf->nVertices);
			break;
		case 3:
			printf("\nOrigem: %d \n", 0);
			dijkstra(matriz, 0, graf->nVertices);
			break;
		case 4:
			floyd(matriz, graf->nVertices);
			break;
		case 5:
			MSTprim(&graf);
			break;
		case 6:
			MSTKruskal(&graf);
			break;
		case 7:
			exit(0);
		default:
			printf("Valor invalido!\n");
		}
	} while (!sair);
	return 0;
}