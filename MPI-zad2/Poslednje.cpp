
#include<iostream>

#include<mpi.h>
#include<string.h>
#include <climits> 
#include<cfloat>
#include <math.h>
#pragma warning(disable : 4996)

void lab1(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int rank, size;
	char prvaPoruka[20], drugaPor[20];
	MPI_Comm comm= MPI_COMM_WORLD;
	MPI_Status st;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	if (rank == 0)
	{
		strcpy(prvaPoruka, "Prva poruka");
		strcpy(drugaPor, "Druga poruka");

		for (int i = 1; i < size; i++)
		{
			MPI_Send(prvaPoruka, 20, MPI_CHAR, i, 1, comm);
			MPI_Send(drugaPor, 20, MPI_CHAR, i, 2, comm);
			int pom;
			MPI_Recv(&pom, 1, MPI_INT, i, 3, comm, &st);
			printf("mater primio: %d", pom);
		}


	}
	else
	{

		MPI_Recv(drugaPor, 20, MPI_CHAR, 0, 2, comm, &st);
		MPI_Recv(prvaPoruka, 20, MPI_CHAR, 0, 1, comm, &st);
		printf("ja sam %d, primio: %s, %s", rank, drugaPor, prvaPoruka);
		MPI_Send(&rank, 1, MPI_INT, 0, 3, comm);
	}
	MPI_Finalize();

}

void lab2a(int argc, char** argv)
{
	int rank, size;
	int n; int min = INT_MAX;
	int* niz = new int[20];
	MPI_Init(&argc, &argv);
	MPI_Status st;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);


	if (rank == 0)
	{
		n = 8;
		for (int i = 0; i < n; i++)
			niz[i] = i;
		int np = n / size;
		for (int i = 1; i < size; i++)
		{
			MPI_Send(&np, 1, MPI_INT, i, 1, comm);
			MPI_Send(&niz[(i-1) * np], np, MPI_INT, i, 2, comm);
			int pom;
			MPI_Recv(&pom, 1, MPI_INT, i, 3, comm, &st);
			if (min > pom)
				min = pom;
		}
		printf("Konacni min %d", min);
	}
	else
	{
		MPI_Recv(&n, 1, MPI_INT, 0, 1, comm, &st);		
		MPI_Recv(niz , n, MPI_INT, 0, 2, comm, &st);
		for (int i = 0; i < n; i++)
		{
			if (min > niz[i])			
				min = niz[i];
			
		}
		MPI_Send(&min, 1, MPI_INT, 0, 3, comm);
		printf("ja sam %d, moj min %d", rank, min);

	}
	MPI_Finalize();
}

void lab2b(int argc, char** argv) // bila greska za min sam stavila int_min i uvek je to manje 
{
	MPI_Init(&argc, &argv);
	MPI_Status st;
	MPI_Comm comm = MPI_COMM_WORLD;
	int rank, size;
	int n, min=INT_MAX, minKonacno;
	int *niz=new int[10];
	int *niz2=new int[10];
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	if (rank == 0)
	{
		n = 8;
		for (int i = 0; i < n; i++)
			niz[i] = i;
		n= n / size;
		//printf("Ja sma 0 nizz je %d %d %d %d", niz[0], niz[1], niz[2], niz[3]);
	}


	MPI_Bcast(&n, 1, MPI_INT, 0, comm);
	MPI_Scatter(niz, n, MPI_INT, niz2, n, MPI_INT, 0, comm);
	printf("ja sam %d , niz moj je: %d %d ", rank, niz2[0], niz2[1]);

	for (int i = 0; i < n; i++)
	{
		if (min > niz2[i])
			min = niz2[i];
	}
	printf("(%d) min je %d", rank, min);
		//MPI_Scan(&min, &minKonacno, 1, MPI_INT, MPI_MIN, comm);
	MPI_Reduce(&min, &minKonacno, 1, MPI_INT, MPI_MIN, 0, comm);
	if (rank == 0)
		printf("Konacno min je %d", minKonacno);

	MPI_Finalize();
}

void lab3(int argc, char** argv)
{
	int rank, size;
	int n;
	MPI_Init(&argc, &argv);
	MPI_Status st;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	int niz[100];
	if (rank == 0)
	{
		n = size * (size + 1);
		for (int i = 0; i < n; i++)
			niz[i] = i;
		for (int i = 1; i < size; i++)
			MPI_Send(&niz[(i - 1) * i], 2 * i, MPI_INT, i, 1, comm);
	}
	else
	{
		MPI_Recv(&niz, rank * 2, MPI_INT, 0, 1, comm, &st);
		int sum = 0;
		for (int i = 0; i < 2 * rank; i++)
		{
			sum += niz[i];
		}
		printf("ja sam %d, sum= %d", rank, sum);
	}
	MPI_Finalize();
}

void lab4(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int size, rank;
	MPI_Status st;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	int niz[10], sum = 0;
	int n = 5, rez = 0, rezOld;
	
		for (int i = 0; i < n; i++)
		{
			niz[i] = rank  + i;
			sum += niz[i];
	
		}
		
		if (sum % 2 == 0)
		{
			rez = rez | 1;
		}
		if (rank < size - 1)		
		{
			MPI_Recv(&rezOld, 1, MPI_INT, rank + 1, 1, comm, &st);
			rezOld= (rezOld << 1);
			rez = rez | rezOld;
		}
		if (rank == 0)
		{
			int i = 0;
			while (i < size)
			{
				printf("proces %d je paran = %d \n ", i, rez & 1);
				rez >>= 1;
				i++;
			}
		}
		else
		{
			MPI_Send(&rez, 1, MPI_INT, rank - 1, 1, comm);
		}
		MPI_Finalize();


	
}


void lab5(int argc, char** argv)
{
	struct { int sum; int id; } radnik;
	MPI_Init(&argc, &argv);
	int size, rank;
	int niz[100];
	int n, sum = 0;
	MPI_Status st;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	int niz2[100];
	if (rank == 0)	
	{
		n = 20;
		for (int i = 0; i < n; i++)
		{
			niz[i] = rand()%100;
		}
		n = n / (size );
		printf("n=%d  ", n);
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, comm);
	MPI_Scatter(&niz, n, MPI_INT, &niz2, n, MPI_INT, 0, comm);
	
	for (int i = 0; i < n; i++)
		sum += niz2[i];
	if (rank == 0)
	{
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&radnik, 1, MPI_2INT, i, 1, comm, &st);
			if (radnik.sum % radnik.id == 0)
				printf("radnik %d ima deljivu sumu %d", radnik.id, radnik.sum);
		}
	}
	else
	{
		radnik.id = rank;
		radnik.sum = sum;
		MPI_Send(&radnik, 1, MPI_2INT, 0, 1, comm);
	}

	MPI_Finalize();

}




void main(int argc, char** argv)
{
	//lab1(argc, argv);
	//lab2a(argc, argv);
	//lab2b(argc, argv);
	//lab3(argc, argv);
	//lab4(argc, argv);
	lab5(argc, argv);
}