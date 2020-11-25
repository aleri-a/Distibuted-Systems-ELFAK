//#include "stdafx.h"
#include<iostream>

#include<mpi.h>
#include<string.h>
#include <climits> 
#include<cfloat>
#pragma warning(disable : 4996)
void FunkcijaA_2zad(int argc, char**argv);
void FunkcijaB_2zad(int argc, char** argv);
void Zad6A(int argc, char**argv);
void ZadSaLabVezbe(int argc, char **argv);
void zad6A_ja_septembar(int argc, char **argv);
void zadSaIspita_slika(int argc, char **argv);
void zadSaIspita_slika_grupneKom(int argc, char ** argv);
void zad1_labvezba(int argc, char** argv);
//void main(int argc, char**argv)
//{
	//FunkcijaA_2zad(argc, argv);
	//Zad6A(argc, argv);
	//zadSaIspita_slika(argc, argv);
	//zadSaIspita_slika_grupneKom(argc, argv);
	//zad6A_ja_septembar(argc, argv);
	//FunkcijaB_2zad(argc, argv);
	//ZadSaLabVezbe(argc, argv);
	//zad1_labvezba(argc, argv);
//}

#pragma region MyRegion



void zad1_labvezba(int argc, char** argv) {
	int rank, size;
	int i;
	int IdProcesa;
	char prvaPoruka[20], drugaPoruka[20];
	MPI_Status st;

	MPI_Init(&argc, &argv);//inicajlizacija MPI kruzenja
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		strcpy_s(prvaPoruka, 20, "Prva poruka");
		strcpy_s(drugaPoruka, 20, "Druga Poruka");
		for (i = 1; i < size; i++)
		{
			MPI_Send(prvaPoruka, 20, MPI_CHAR, i, 1, MPI_COMM_WORLD);
			MPI_Send(drugaPoruka, 20, MPI_CHAR, i, 2, MPI_COMM_WORLD);

			MPI_Recv(&IdProcesa, 1, MPI_INT, i, 55, MPI_COMM_WORLD, &st);
			printf("Proces %d je primio poruke\n", IdProcesa);
		}
	}
	else
	{
		MPI_Recv(prvaPoruka, 20, MPI_CHAR, 0, 2, MPI_COMM_WORLD, &st);
		MPI_Recv(drugaPoruka, 20, MPI_CHAR, 0, 1, MPI_COMM_WORLD, &st);
		printf("Proces %d je primio poruke: %s    %s  \n", rank, prvaPoruka, drugaPoruka);

		MPI_Send(&rank, 1, MPI_INT, 0, 55, MPI_COMM_WORLD);

	}






	MPI_Finalize();
}


void FunkcijaB_2zad(int argc, char** argv) //zad2 pod B
{
	int rank, size;
	int i;
	int IdProcesa;
	int n; int np;
	float niz[100]; //0ti ga salje
	float niz2[100];
	MPI_Status st;
	float min;
	float minMinimuma = FLT_MAX;

	MPI_Init(&argc, &argv);//inicajlizacija MPI kruzenja
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		printf("Broj elementaka niza  ");
		fflush(stdout);//
		scanf("%d", &n);
		for (i = 0; i < n; i++)
			scanf("%f", &niz[i]);
		np = n / (size - 1);
	}
	MPI_Bcast(&np, 1, MPI_INT, 0, MPI_COMM_WORLD);//svi primau iy niz u niz2 2 np karaktera
	MPI_Scatter(niz, np, MPI_FLOAT, niz2, np, MPI_FLOAT, 0, MPI_COMM_WORLD); // 


	for (i = 0; i < np; i++)
		if (niz2[i] < minMinimuma)
			minMinimuma = niz2[i];

	MPI_Reduce(&minMinimuma, &min, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		printf("min je %f", min);
	}
	
	MPI_Finalize();
}

void FunkcijaA_2zad(int argc, char**argv)
{
	int rank, size;
	int i;
	int IdProcesa;
	int n; int np;
	float niz[100];
	MPI_Status st;
	float min;
	float minMinimuma = FLT_MAX;

	MPI_Init(&argc, &argv);//inicajlizacija MPI kruzenja
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		printf("Broj elementaka niza  ");
		fflush(stdout);//
		scanf("%d", &n);
		for (i = 0; i < n; i++)
			scanf("%f", &niz[i]);
		np = n / (size - 1); //zato sto je jedan proces master i njega ne racunas 
		for (i = 1; i < size; i++)
		{

			MPI_Send(&np, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(&niz[(i - 1)*np], np, MPI_FLOAT, i, 2, MPI_COMM_WORLD);

			MPI_Recv(&min, 1, MPI_FLOAT, i, 55, MPI_COMM_WORLD, &st);
			if (min < minMinimuma)
				minMinimuma = min;

		}
		printf("Proces %f minimumMinimuma\n", minMinimuma);
	}
	else
	{
		MPI_Recv(&np, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &st);
		MPI_Recv(niz, np, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &st);
		for (i = 0; i < np; i++)
			if (niz[i] < minMinimuma)
				minMinimuma = niz[i];

		MPI_Send(&minMinimuma, 1, MPI_FLOAT, 0, 55, MPI_COMM_WORLD);

	}






	MPI_Finalize();

}

void Zad1_A(int argc, char **argv)
{
	char prvaPoruka[100], drugaPoruka[100];
	int rank, size, i, povratnaPoruka;
	MPI_Status st;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		//MASTER PROCES
		strcpy_s(prvaPoruka, 100 * sizeof(char), "PRVA PORUKA");
		strcpy_s(drugaPoruka, 100 * sizeof(char), "DRUGA PORUKA");
		for (i = 1; i < size; i++) {
			MPI_Send(prvaPoruka, 100, MPI_CHAR, i, 1, MPI_COMM_WORLD);
		}
		for (i = 1; i < size; i++) {
			MPI_Send(drugaPoruka, 100, MPI_CHAR, i, 2, MPI_COMM_WORLD);
		}

		for (i = 1; i < size; i++) {
			MPI_Recv(&povratnaPoruka, 1, MPI_INT, i, 3, MPI_COMM_WORLD, &st);
			printf("Povratna poruka: %d\n", povratnaPoruka);
			fflush(stdout);
		}
	}
	else {
		//RADNICI PROCES
		MPI_Recv(drugaPoruka, 100, MPI_CHAR, 0, 2, MPI_COMM_WORLD, &st);
		printf("Proces %d je primio poruku %s\n", rank, drugaPoruka);
		fflush(stdout);
		MPI_Recv(prvaPoruka, 100, MPI_CHAR, 0, 1, MPI_COMM_WORLD, &st);
		printf("Proces %d je primio poruku %s\n", rank, prvaPoruka);
		fflush(stdout);
		MPI_Send(&rank, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);// procesu 0, id poruke je 3
	}

	MPI_Finalize();

}

void Zad2_A(int argc, char **argv)
{
	int n, rank, size, n_p, i;
	float niz[100], loklaniMin, Min;
	MPI_Status st;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		printf("Unesite velicinu niza: ");
		fflush(stdout);
		scanf("%d", &n);
		//niz = (float*)(malloc(n * sizeof(float)));
		printf("Unesite niz: ");
		fflush(stdout);
		for (i = 0; i < n; i++) {
			scanf("%f", &niz[i]);
		}

		n_p = n / (size - 1);
		for (i = 1; i < size; i++) {
			MPI_Send(&n_p, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
			MPI_Send(&niz[(i - 1) * n_p], n_p, MPI_FLOAT, i, 1, MPI_COMM_WORLD);
		}
		MPI_Recv(&Min, 1, MPI_FLOAT, 1, 2, MPI_COMM_WORLD, &st);

		for (i = 2; i < size; i++) {
			MPI_Recv(&loklaniMin, 1, MPI_FLOAT, i, 2, MPI_COMM_WORLD, &st);
			//printf("lokalni min: %f", loklaniMin);
			if (loklaniMin < Min) {
				Min = loklaniMin;
			}
		}
		printf("Minimum je %f", Min);
		fflush(stdout);
		//free(niz);
	}
	else {
		MPI_Recv(&n_p, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &st);
		//niz = (float*)(malloc(n_p * sizeof(float)));
		MPI_Recv(&niz, n_p, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &st);
		loklaniMin = niz[0];
		for (i = 1; i < n_p; i++) {

			if (loklaniMin > niz[i]) {
				loklaniMin = niz[i];
			}
		}
		MPI_Send(&loklaniMin, 1, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
		//free(niz);
	}
	MPI_Finalize();
}

/*
void Zad3(int argc, char **argv)
{
	
	int rank, size;
	int i;
	int IdProcesa;
	int pom2;//od kog da salje;
	int n; int pom;
	int sum; int niz3[100];
	MPI_Status st;
	float min;
	float minMinimuma = FLT_MAX;

	MPI_Init(&argc, &argv);//inicajlizacija MPI kruzenja
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		printf("Broj elementaka niza  ");
		fflush(stdout);//
		scanf("%d", &n);
		for (i = 0; i < n; i++)
			scanf("%d", &niz3[i]);
		pom = 2; //zato sto 1ti prima 2 broja prvi 
		pom2 = 0;//salje se od 0tog el niza
		for (i = 1; i < size; i++)
		{

			MPI_Send(&niz3[pom2] , pom, MPI_INT, i, 1, MPI_COMM_WORLD);

			//ja sam ovako napisala  -> MPI_Send(&niz[(i-1) * 2], i * 2 , MPI_INT, i, 22, MPI_COMM_WORLD);
			pom2 += pom;
			pom += 2;
			
		}
	}
	else
	{
		
		MPI_Recv(niz3, rank*2, MPI_INT, 0, 1, MPI_COMM_WORLD, &st);
		
		for (i = 0; i < rank * 2; i++)			
				sum += niz3[i];

		printf("Proces %d, suma je %d \n", rank, sum);

	}

	MPI_Finalize();
	
}

*/



////Napisati MPI program koji realizuje množenje matrice A reda n i vektora x, čime se dobija rezultujući vektor c. 
//Množenje se obavlja tako što master proces (sa rankom 0) inicijalizuje vrednosti za matricu A i vektor x i nakon toga šalje 
//svakom procesu po jednu
//vrstu prve matrice (direktno iz matrice A) i ceo vektor x. Svi procesi učestvuju u izračunavanju. Konačni rezultat se
//generiše i prikazuje u master procesu. Zadatak rešiti korišćenjem isključivo Point-to-Point operacija.
void ZadSaLabVezbe(int argc, char **argv)
{
	int rank, size;
	int i, j;
	int n;	
	int  X[100];
	int  C[100];
	int  A[100][100];
	int pomNiz[100];
	int rez = 0;
	
	MPI_Status st;


	MPI_Init(&argc, &argv);//inicajlizacija MPI kruzenja
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		printf("Broj elementaka vektora X  ");
		
		scanf("%d", &n);
		printf("Elementi niza su : \n");
		for (i = 0; i < n; i++)
			scanf("%d", &X[i]);
		printf("Unesite elemente matrice");
		for (i = 0; i < n; i++)
			for(j=0;j<n;j++)
			scanf("%d", &A[i][j]);
		fflush(stdout);
		for (i = 1; i < size; i++)
		{

			MPI_Send(&n, 1, MPI_INT, i, 1, MPI_COMM_WORLD);//koliko podataka ima niz koji primaju 
			MPI_Send(&A[n*i], n, MPI_INT, i, 2, MPI_COMM_WORLD);
			MPI_Send(&X, n, MPI_INT, i, 3, MPI_COMM_WORLD);

			MPI_Recv(&C[i], 1, MPI_INT, i, 55, MPI_COMM_WORLD, &st);	
			
		}
		printf("Proces %s Konacni niz je \n", C);
	}
	else
	{
		MPI_Recv(&n, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &st);
		MPI_Recv(pomNiz, n, MPI_INT, 0, 2, MPI_COMM_WORLD, &st);
		MPI_Recv(X, n, MPI_INT, 0, 3, MPI_COMM_WORLD, &st);
		for (i = 0; i < n; i++)
		{
			rez += pomNiz[i] * X[rank+1];
		}

		MPI_Send(&rez, 1, MPI_INT, 0, 55, MPI_COMM_WORLD);

	}


	
}

void Zad33(int argc, char** argv)
{
	int rank, size, *niz, sum, n, i, pom;
	MPI_Status st;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		//master
		n = size * (size - 1);
		niz = (int*)malloc(n * sizeof(int));
		printf("Unesite niz: ");
		fflush(stdout);

		for (i = 0; i < n; i++) {
			//std::cin >> niz[i];
			scanf("%d", &niz[i]);
		}
		printf("Ucitano1\n");
		fflush(stdout);
		pom = 2;
		sum = 0;
		for (i = 1; i < size; i++) {
			MPI_Send(&niz[sum], pom, MPI_INT, i, 1, MPI_COMM_WORLD);
			sum += pom;
			pom += 2;
		}

		free(niz);
	}
	else {
		//radnik
		pom = rank << 1;
		niz = (int*)malloc(pom * sizeof(int));
		MPI_Recv(niz, pom, MPI_INT, 0, 1, MPI_COMM_WORLD, &st);
		sum = 0;
		for (i = 0; i < pom; i++) {
			sum += niz[i];
		}
		printf("Suma procesa %d je : %d", rank, sum);
		fflush(stdout);
		free(niz);
	}

	MPI_Finalize();

}
void domaci(int argc, char**argv)
{
	int rank, size, sum = 0, recv;
	int send = 1;
	int shoudlSend = 2;
	int pair = 1;
	MPI_Status st;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank % shoudlSend != 0) {
		MPI_Send(&rank, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD);
	}
	else {

		MPI_Recv(&recv, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD, &st);
		sum = rank + recv;
		size = size >> 2;
		while (size > 0) {
			shoudlSend = shoudlSend << 1;
			pair = pair << 1;
			if (rank % shoudlSend != 0) {
				MPI_Send(&sum, 1, MPI_INT, rank - pair, 1, MPI_COMM_WORLD);
				break;
			}
			MPI_Recv(&recv, 1, MPI_INT, rank + pair, 1, MPI_COMM_WORLD, &st);
			sum += recv;
			size = size >> 1;
		}
		if (rank == 0) {
			printf("Suma brojeva je %d %d", sum, size);
		}
	}

	MPI_Finalize();
}

void zad4(int argc, char**argv)
{
	int rank, size, niz[20], n, n_p, lokalniNiz[20], sum, status;
	MPI_Status st;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		//Master
		//srand(time(nullptr));
		n = 20;
		printf("Array: ");
		for (int i = 0; i < n; i++) {
			niz[i] = rand() % 100;
			printf("%d ", niz[i]);
		}
		fflush(stdout);
		n_p = n / (size);
	}//sve ovo sluzi da napravi niz i svima posalje elemente neke da imaju u svom nizu 
	MPI_Bcast(&n_p, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(&niz, n_p, MPI_INT, &lokalniNiz, n_p, MPI_INT, 0, MPI_COMM_WORLD);

	sum = 0;

	for (int i = 0; i < n_p; i++) {
		sum += lokalniNiz[i];
	}
	if (rank == size - 1) 
	{
		if (sum % 2 == 0) {
			status = 1;
		}
		else {
			status = 0;
		}
		MPI_Send(&status, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD);//poslace ovom nizem od sebeg
	}
	else if (rank != 0) 
	{
		MPI_Recv(&status, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD, &st);//primi je od ovog viseg od njega 
		status <<= 1;
		if (sum % 2 == 0) {
			status |= 1;
		}
		MPI_Send(&status, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD);//nastavlja lanac da salje ovom ispod sebe 
	}
	else {
		MPI_Recv(&status, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &st);//master 
		status <<= 1;
		if (sum % 2 == 0) {
			status |= 1;
		}
		printf("\nProcesi sa parnim brojem: ");
		for (int i = 0; i < size; i++) {
			if (status & 1) {//ispitace samo poslednji bit 
				printf("%d ", i);
				fflush(stdout);
			}
			status >>= 1;
		}
	}
	MPI_Finalize();
}


void zad5(int argc, char**argv)
{
	int rank, size, sum, niz[10], uslov, i;
	MPI_Status st;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//srand(time(nullptr));


	if (rank == 0) {
		for (i = 1; i < size; i++) {
			MPI_Recv(&uslov, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &st);
			if (!uslov) { //u uslovu ce biti 0 ako je deljivo 
				printf("%d ", i);
				fflush(stdout);
			}
		}
	}
	else {
		sum = 0;
		for (i = 0; i < 10; i++) {
			niz[i] = rand() % 100;
			sum += niz[i];
		}
		uslov = (sum % rank);
		MPI_Send(&uslov, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}


	MPI_Finalize();
}

#pragma endregion


void ScanFja_prekoPTP(int argc, char **argv) {//pokusaj 6og zad ali nije to to 
	int size, rank, sum, parcSum, j, i;

	MPI_Status st;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	sum = rank;

	j = rank + 1;
	//printf("%d\n", j);
	while (j < size)
	{
		//printf("ja sam %d, i ja cu poslati %d \n", rank, j);
		MPI_Send(&rank, 1, MPI_INT, j, 33, MPI_COMM_WORLD);
		//printf("%d\n", rank);
		j++;
	}

	i = rank - 1;
	sum = 0;
	while (i >= 0) {
		MPI_Recv(&parcSum, 1, MPI_INT, i, 33, MPI_COMM_WORLD, &st);
		sum += parcSum + 1;
		i--;

	}
	sum += rank + 1;
	printf("Ovo je proces %d, i suma prethodnika je %d\n ", rank, sum);
	MPI_Finalize();
}
void Zad6A(int argc, char**argv)  //????
{
	int rank, size, sum = 0, numOfSteps = 0, i, lokalniBafer, baferZaPrimanje, dest, stepen2 = 1;
	MPI_Status st;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	for (i = 1; i < size; i <<= 1) { //ako je size =6 
		numOfSteps++;	//i ovde stoji 	printf("%d ",i); onda ce izbaciti 2 6 kao rezultat 
	}
	lokalniBafer = rank;
	sum = rank;
	for (i = 0; i < numOfSteps; i++) {
		dest = rank ^ (1 << i); //ovo je ^ XOR tj  dve iste 0 dve razlicite 1  ==> tj 00001, ako je i=2 to je 00010 tj 2
		//i=3 je 000100 tj 4 itd 
		printf("ja sam %d i saljem %d", rank, dest);
		MPI_Send(&lokalniBafer, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
		MPI_Recv(&baferZaPrimanje, 1, MPI_INT, dest, 1, MPI_COMM_WORLD, &st);
		lokalniBafer += baferZaPrimanje;
		if (dest < rank) {
			sum += baferZaPrimanje;
		}
	}

	printf("Proces %d ima sumu: %d", rank, sum);
	fflush(stdout);
	MPI_Finalize();
}


void zad6B(int argc, char**argv)
{
	int rank, size, sum;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Scan(&rank, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	printf("Proces %d sa sumom: %d", rank, sum);

	MPI_Finalize();
}
void zad6A_by_Masa(int argc, char **argv) {

	int rank, size, b, niz[8], log, s, pom1, pom2;
	double niz2[8], suma, max;
	struct {
		double val;
		int rank;
	} in, out;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	
	log = log2(size); //3
	if (rank == 0)
	{
		for (int i = 0; i < size; i++)
			niz[i] = i;
	}
	s = 1;
	int pom3 = size;
	b = size / 2;
	for (int i = 0; i < log; i++)
	{
	
		pom1 = 0;
		for (int j = 0; j < s; j++) {//bilo je j<s

			if (rank == pom1) {
				MPI_Send(&niz[b], b, MPI_INT, pom1 + b, 15, MPI_COMM_WORLD);
			}
			if (rank == pom1 + b)
			{
				MPI_Recv(niz, b, MPI_INT, pom1, 15, MPI_COMM_WORLD, &status);
			}
			pom1 = pom1 + pom3;
		}
		s = s * 2;
		b = b / 2;
		pom3 = pom3 / 2;

	}

	printf("ja sam %d, i moj br je %d", rank, niz[0]);
	MPI_Finalize();

}


void zad7(int argc, char **argv)
{
	//Napisati MPI program koji izračunava sumu N celih brojeva (N je stepen 2)
	//korišćenjem Point-to-point komunikacije. U prvom korako procesi se grupišu u parove (P0, P1),
	//(P2, P3),…, (Pp-2, Pp-1). Zatim se izačunavaju parcijalne sume u svim parovima korišćenjem
	//P-to-P komunikacije i akumuliraju u procesima (P0, P2,…, Pp-2). Npr. process Pi (i-parno) 
	//izračunava parcijalne sume za par procesa (Pi , Pi+1). U sledećem koraku razmatraju se parovi
	//procesa (P0, P2), (P4, P6),…, (Pp-4, Pp-2) pronalaze parcijalne sume i akumuliraju u
	//(P0, P4,…, Pp-4). Postupak se ponavlja dok ne ostanu 2 procesa 
	//i rezultat se akumulira u P0. Zadatak rešiti i korišćenjem jedne grupne operacije.
	int rank, size, sum;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Reduce(&rank, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		printf("Suma prvih %d brojeva je: %d", size, sum);
	

	}

	MPI_Finalize();
}




//void zadSaStablon(int argc, char **argv)
//{
//	// Masin kod 
//	int rank, size, b, i, pom = 1;
//	MPI_Status status;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//	b = 1;
//
//	if (rank == 0)
//		pom = 5;
//	/*
//
//	ja bih ovde stavila
//	int numObilzaka=0;
//	for(i=1;i<size;i<<=1)
//		numObilazka++;
//	*/
//	for (int j = 0; j < numObilazka; j++)
//	{
//
//	}
//
//
//	while (b < size)
//	{
//		for (i = 0; i < b; i++)
//		{
//			if (rank == i)
//			{
//				MPI_Send(&pom, 1, MPI_INT, rank + b, 15, MPI_COMM_WORLD);
//			}
//			if (rank == i + b)
//			{
//				MPI_Recv(&pom, 1, MPI_INT, i, 15, MPI_COMM_WORLD, &status);
//				printf("ja sam proces %d i moj pom je %d", rank, pom);
//			}
//
//		}
//		b *= 2;
//	}
//	MPI_Finalize();
//}

//Napisati MPI program koji izračunava i prikazuje skalarni proizvod dva vektora a i b dužine m, 
//korišćenjem podele vektora na blokove, sa cikličnom uniformnom distribucijom podataka vektora,
//kao na slici.Vektori a i b se inicijalizuju u master procesu.Na slici je prikazan raspodela podataka 
//vektora a po procesima  na primeru za m = 16 i p = 4 
//(p je broj procesa, m deljivo sa p).Ista raspodela podataka se podrazumeva i za vektor b.
//Nakon raspodele elemenata vektora a i b po procesima, obaviti odgovarajuća izračunavanja,
//generisati i prikazati skalarni proizvod.Zadatak rešiti korišćenjem :
//a) Isključivo grupnih operacija
//b) Point - to - Point operacija za deo koda kojim se vrši slanje / primanje podataka
//ja radila ova 2 zadatka i rade 
void zadSaIspita_slika(int argc, char **argv) {
	int size, rank, np, n, parcSum, sum = 0;
	int niz1[100], niz2[100];
	MPI_Status st;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
		printf("Unesite br n");
		fflush(stdout);
		scanf("%d", &n);
		for (int i = 0; i < n; i++)
		{
			niz1[i] = rand()%10 +1;
			printf("%d ", niz1[i]);
			niz2[i] = rand() % 10+1;
	    }
		printf("\n");
		for (int i = 0; i < n; i++)
			printf("%d ", niz2[i]);
		printf("\n");
		fflush(stdout);
		
		np = n / (size-1);
		for (int i = 1; i < size; i++)
			MPI_Send(&np, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
		
		int j;
			for (int i = 1; i < size; i++)
			{
				
				for (j=i-1 ; j < n; j=(size-1) +j)
				{
					MPI_Send(&niz1[j], 1, MPI_INT, i, 11, MPI_COMM_WORLD);
					MPI_Send(&niz2[j], 1, MPI_INT, i, 22, MPI_COMM_WORLD);
				}
				
				MPI_Recv(&parcSum, 1, MPI_INT, i, 55, MPI_COMM_WORLD, &st);
				sum += parcSum;
			}
			printf("rez je %d", sum);

	}
	else
	{
		parcSum = 0;
		MPI_Recv(&np, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &st);
		for (int i = 0; i < np; i++)
		{
			MPI_Recv(&niz1[i], 1, MPI_INT, 0, 11, MPI_COMM_WORLD, &st);
			MPI_Recv(&niz2[i], 1, MPI_INT, 0, 22, MPI_COMM_WORLD, &st);
			printf("Ja sam %d i primio sam niz1: %d  niz2: %d \n", rank, niz1[i], niz2[i]);
			parcSum += niz1[i] * niz2[i];
		}
		MPI_Send(&parcSum, 1, MPI_INT, 0, 55, MPI_COMM_WORLD);


	}
	MPI_Finalize();


}

void zadSaIspita_slika_grupneKom(int argc, char ** argv) {
	int size, rank, np, n, parcSum = 0 , sum = 0;
	int niz1[100], niz2[100], niz13[100], niz23[100];
	MPI_Status st;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
		printf("Unesite br n");
		fflush(stdout);
		scanf("%d", &n);
		for (int i = 0; i < n; i++)
		{
			niz1[i] = rand() % 10 + 1;
			printf("%d ", niz1[i]);
			niz2[i] = rand() % 10 + 1;
		}
		printf("\n");
		for (int i = 0; i < n; i++)
			printf("%d ", niz2[i]);
		printf("\n");
		fflush(stdout);

		int j = 0;
		int pomJ = 0;
		int i = 0;
		while (j < size) {		
			printf("%d ", j);
			for ( j ; j < n; j = size + j)
			{
				niz13[i] = niz1[j];
				niz23[i] = niz2[j];
				i++;
			}
			pomJ++;
			j = pomJ;
		}
		printf("elementi nizova 13 i 23 \n");
		for (int i = 0; i < n; i++)
			printf("%d ", niz13[i]);
		printf("\n");
		for (int i = 0; i < n; i++)
			printf("%d ", niz23[i]);
		fflush(stdout);



		np = n / size;
	}
	MPI_Bcast(&np, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(niz13, np, MPI_INT, niz1, np, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(niz23, np, MPI_INT, niz2, np, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < np; i++) {
		printf("ja sam proces %d i moji niz1: %d  niz2: %d \n",rank, niz1[i],niz2[i]);
	}

	for (int i = 0; i < np; i++)
		parcSum += niz1[i] * niz2[i];
	MPI_Reduce(&parcSum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0)
		printf("rez je %d \n ", sum);
		

	MPI_Finalize();
}

/*
void stabloByMasa(int argc, char **argv)
{
	int rank, size, b, i, pom = 1;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	b = 1;

	if (rank == 0)
		pom = 5;
	
	//ja bih ovde stavila
	//int numObilzaka=0;
	//for(i=1;i<size;i<<=1)
	//	numObilazka++;
	
	for (int j = 0; j < numObilazka; j++)
	{

	}


	while (b < size)
	{
		for (i = 0; i < b; i++)
		{
			if (rank == i)
			{
				MPI_Send(&pom, 1, MPI_INT, rank + b, 15, MPI_COMM_WORLD);
			}
			if (rank == i + b)
			{
				MPI_Recv(&pom, 1, MPI_INT, i, 15, MPI_COMM_WORLD, &status);
				printf("ja sam proces %d i moj pom je %d", rank, pom);
			}

		}
		b *= 2;
	}
	MPI_Finalize();
}*/


void stabloByMe(int argc, char **argv) 

{

}



//by Danica 7 zad lab 
void fja(int argc, char* argv[])
{
	int rank, size, sum = 0, sumrcv = 0, j, x, a, b, n = 16, odnos, brojiter, l, m;
	int niz[] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
	int pom[5];
	MPI_Status st1;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);//br trenutnog procesa
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	odnos = n / size;
	brojiter = log2(size);
	b = 2;
	l = 1;
	for (j = 0; j < brojiter; j++)
	{
		if (j == 0)
		{
			if (rank % 2 == 1){ // tj neparno salje 
				a = 0;
				for (x = rank * odnos; x < rank * odnos + odnos; x++)
					a += niz[x];
				MPI_Send(&a, 1, MPI_INT, rank - 1, rank - 1, MPI_COMM_WORLD);
				memset(pom, 0, 5);}
			else{
				sum = 0;
				MPI_Recv(&m, 1, MPI_INT, rank + 1, rank, MPI_COMM_WORLD, &st1);
				sum += m;
				for (x = rank * odnos; x < rank * odnos + odnos; x++)
					sum += niz[x];
			}
		}
		else
		{
			if (rank % l == 0)
			{
				if (rank % b == 0)
				{
					MPI_Recv(&sumrcv, 1, MPI_INT, rank + l, rank + l, MPI_COMM_WORLD, &st1);
					sum += sumrcv;
				}
				else if (rank % b == l)
				{
					MPI_Send(&sum, 1, MPI_INT, rank - l, rank, MPI_COMM_WORLD);
				}
			}
		}
		b *= 2; l *= 2;
	}
	if (rank == 0)
		std::cout << sum;

	MPI_Finalize();
}
