#include<iostream>

#include<mpi.h>
#include<string.h>
#include <climits> 
#include<cfloat>
#include <math.h>
#pragma warning(disable : 4996)

//0 salje 2 poruke, ostali primaju prvo 2. pa prvu i stampaju ih
//i vracaju svoj id 0li
void zad1_lab(int argc, char** argv)
{
	int rank, size;
	int i, idProcesa;
	char prvaPoruka[20], drugaPoruka[20];
	MPI_Status st;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0)
	{
		strcpy_s(prvaPoruka, 20, "Prva poruka");
		strcpy_s(drugaPoruka, 20, "Druga Poruka");
		for (int i = 1; i < size; i++)
		{
			MPI_Send(prvaPoruka, 20, MPI_CHAR, i, 1, MPI_COMM_WORLD);
			MPI_Send(drugaPoruka, 20, MPI_CHAR, i, 2, MPI_COMM_WORLD);
			MPI_Recv(&idProcesa, 1, MPI_INT, i, 3, MPI_COMM_WORLD, &st);
			printf("0ti primio potvrdu od %d \n", idProcesa);
		
		}
	}
	else
	{
		MPI_Recv(drugaPoruka, 20, MPI_CHAR, 0, 2, MPI_COMM_WORLD, &st);		
		MPI_Recv(prvaPoruka, 20, MPI_CHAR, 0, 1, MPI_COMM_WORLD, &st);
		printf("Proces %d je primio poruke: %s   %s  \n", rank, drugaPoruka, prvaPoruka);
		MPI_Send(&rank, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
	}

}


//0 podeli svakom procexu N/p el iz niza , svako najde min u sebi i vrati masteru
// master odredi ptravi min i odstampa ga 
void zad2a_lab(int argc, char** argv) //Ptp
{
	int rank, size;
	MPI_Status st;
	int niz[100];
	int min = INT_MAX;
	int n;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0)
	{
		// n = (size-1)*2;
		n = 9;
		for (int i = 0; i <= n; i++)
			niz[i] = i+1;
		for (int i = 1; i < size; i++)
		{
			int pomMin;
			int startIndex = n / (size - 1) * (i-1);
			int brElzaslanje = n / (size - 1);
			
			MPI_Send(&brElzaslanje, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(&niz[startIndex], brElzaslanje, MPI_INT, i, 2, MPI_COMM_WORLD);
			MPI_Recv(&pomMin, 1, MPI_INT, i, 3, MPI_COMM_WORLD, &st);
			if (min > pomMin)
				min = pomMin;
		}

		printf("Cao od 0, konacan min je %d", min);
	}
	else
	{
		MPI_Recv(&n, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &st);
		MPI_Recv(niz, n, MPI_INT, 0, 2, MPI_COMM_WORLD, &st);
		for (int i = 0; i < n; i++)
			if (min > niz[i])
				min = niz[i];
		MPI_Send(&min, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
		printf("Cao od %d, moj min je %d", rank, min);

	}

	MPI_Finalize();

}

void zab2b_lab(int argc, char** argv)
{
	int rank, size, min=INT_MAX, minMinumuma,n ;
	MPI_Status st;
	int niz[100];
	int niz2[100];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0)
	{
		n = size*3;
		for (int i = 0; i <= n; i++)
			niz[i] = size-i + 1;		
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int np = n / (size );
	
	MPI_Scatter(&niz, np, MPI_INT, &niz2, np,MPI_INT, 0, MPI_COMM_WORLD);
	
		//printf("ja sam %d, niz je %d, %d, %d\n", rank, niz[rank*3],niz[rank*3+1],niz[rank*3+2]);
		//printf("ja sam %d, niz2 je %d, %d, %d\n", rank, niz2[0],niz2[1],niz2[2]);
		fflush(stdout);

		for (int i = 0; i < np; i++)
			if (min > niz2[i])
				min = niz2[i];

		printf("ja sam %d, moj min %d\n", rank, min);

	

	MPI_Reduce(&min, &minMinumuma, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	if (rank == 0)
		printf("ja sam %d i krajnji min je %d", rank, minMinumuma);


	MPI_Finalize();


}


// niz n=(p+1)*p, svaki radnik dobije niz duzine 2*i tj
// size=4, n=12, p=3 jer je toj br procesa radnika ,   1. dobije 2 el, 2. -> 4, 3. ->6 
//odredi sumu svakog radnika 
void zad3_lab(int argc, char** argv)
{
	int rank, size, n;
	MPI_Status st;
	int niz[100];

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	if (rank == 0)
	{
		n = (size - 1) * (size - 1 + 1);
		for (int i = 0; i < n; i++)
		{
			niz[i] = i + 1;
		}

		for (int i = 1; i < size; i++)
		{
			int sendNum = i * 2;
			MPI_Send(&sendNum, 1, MPI_INT, i,1, MPI_COMM_WORLD);
			MPI_Send(&niz[(i - 1) * i], sendNum, MPI_INT, i, 2, MPI_COMM_WORLD);
		
		}

	}
	else
	{
		MPI_Recv(&n, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &st);
		MPI_Recv(&niz, n, MPI_INT, 0, 2, MPI_COMM_WORLD, &st);
		int sum = 0;
		for (int i = 0; i < n; i++)
		{
			sum += niz[i];
		}
		printf("ja sam %d, suma je %d", rank, sum);
	}

	MPI_Finalize();

}


//master nadje i ispise ID procesa cija suma el u nizu parna, trazi od poslednjeg ka masteru
// svaki proc prosledjuje info od ovog prethodnog da li je prethodni imao parnu osobinu i onda master odtsmpa taj niz ili kako vec 
// prosledjuje se podatak velicine int gde je jedan bit za jedan proces znaci 8 bit za 8. proces, ako je 1 onda je suma u tom procesu parna
//RESENJE: u check pamtim te bitove, kad primim zavrsni cek:
//(i-1) zato sto sam na bitu 0 upisivala informaccije za proces 1, i onda sa >>(i-1) dovedem na poslednji bit 
// bit koji oznacava za taj proces i da li je paran, sa &1 proveravam taj zadnji bit, ako je jedan odtsmpace 1 ako ne onda je 0 i taj proces nije imao parnu sumu
//MILAN radio na drugi nacin 
void zad4_lab(int argc, char** argv)
{

	int rank, size, n;
	MPI_Status st;
	int niz[100];
	int check;


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0)
	{
		MPI_Recv(&check, 1, MPI_INT, 1, 2, MPI_COMM_WORLD, &st);
		
		for (int i = 1; i < size; i++)
		{
			int checkCopy = check;
			printf("Proces %d je %d\n", i, (checkCopy >>= (i-1)) & 1);
		}
		printf("konacan check %d", check);
	}
	else
	{
		n = 6;
		for (int i = 0; i < n; i++)
			niz[i] = i * rank;

		int sum = 0;
		for (int i = 0; i < n; i++)
			sum += niz[i];
		
		printf("Ja sam %d i sum= %d", rank, sum);
		fflush(stdout);

		if (rank != size - 1)
			MPI_Recv(&check, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD, &st);
		else
			check = 0; //all bits 0

		int mask = 1;
		mask <<=(rank- 1);
		if (sum % 2 == 0)
		{
			check |= mask;
			printf("Ja sam %d i postaljam maska je %d", rank, mask);
		}
		if (rank != 1)
			MPI_Send(&check, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD);
		else
			MPI_Send(&check, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);

	}
	MPI_Finalize();

}

//master nadje id onog cija suma el u nizu deljiva sa idijem tog radika, kom izmedju mastera i radnika pojedinacnom 
void zad5_lab(int argc, char** argv)
{

	int size, rank, n,np, sum=0;
	int niz[100], niz2[100];
	struct { int sum;  int id; } radnik;
	MPI_Status st;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// rank 0 kreira niz i salje svim korisnicima, ni ne mora tako
	if (rank == 0)
	{
		n = 20;
		for (int i = 0; i < n; i++)
		{
			niz[i] = rand() % 100;
			printf("%d ", niz[i]);
		}
		 np = n / (size ); // nije size-1 jer ce bcast i 0tom dodeliti deo niza
	}
	MPI_Bcast(&np, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(&niz, np, MPI_INT, &niz2, np, MPI_INT, 0, MPI_COMM_WORLD);


	for (int i = 0; i < np; i++)
	{
		sum += niz2[i];
	}


	if (rank == 0)
	{
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&radnik, 1, MPI_2INT, i, 1, MPI_COMM_WORLD, &st);			
			if (radnik.sum % radnik.id == 0)
				printf("\n%d / %d", radnik.sum, radnik.id);
			fflush(stdout);
		}
	}
	else
	{
		radnik.id = rank;
		radnik.sum = sum;
		MPI_Send(&radnik, 1, MPI_2INT, 0, 1, MPI_COMM_WORLD);
		
	}

	MPI_Finalize();
}

//suma N brojeva, N je stepen 2, 
// rekla sma n=size 
// 1 proces 1 broj tj radim zbir brojeva od 1 do n 
/*void zad7a_lab(int argc, char** argv)
{
	int size, rank, n = 16;
	int sum = 0; int br, odDrugog;
	//int niz[] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
	MPI_Status st;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	br = rank + 1;
	int brIteracija = log2(n); // 4

	for (int iter = 0; iter < brIteracija; iter++)
	{
		int rcvFrom = rank + iter+1; 
		if (rank % (iter + 1) == 0)
		{
			if (iter == 0 && rank % 2 != 0) // (p0, P1) (P2,P3)-> neparno salje
			{
				MPI_Send(&br, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD);
			}
			else
			{
				MPI_Recv(&br, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD, &st);
			}
		

		}

		MPI_Recv(&odDrugog, 1, MPI_INT, rcvFrom, 1, MPI_COMM_WORLD);

	}
	

}  

void zad7b_lab(int argc, char** argv)
{
	int size, rank, n, br[100];
	MPI_Status st;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);



}*/

////Napisati MPI program koji realizuje množenje matrice A reda n i vektora x, čime se dobija rezultujući vektor c. 
//Množenje se obavlja tako što master proces (sa rankom 0) inicijalizuje vrednosti za matricu A i vektor x i nakon toga šalje 
//svakom procesu po jednu
//vrstu prve matrice (direktno iz matrice A) i ceo vektor x. Svi procesi učestvuju u izračunavanju. Konačni rezultat se
//generiše i prikazuje u master procesu. Zadatak rešiti korišćenjem isključivo Point-to-Point operacija.
void mojSaLab(int argc, char** argv) // matrica salje vrste, treba 4 procesa RADIIII
{
	int size, rank, n;
	int mat[10][10], vector[5]  ;
	int deoMatrice[10];
	MPI_Status st;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
		n = 3;
		int mat2[3][3] = { {3,2,4},{1,8,9},{1,1,2} };
		int vec2[3] = { 3,4,5 };
		for (int i = 1; i < size; i++)
		{
			int pom;
			MPI_Send(&n, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(mat2[i-1], n, MPI_INT, i, 2, MPI_COMM_WORLD);  
			MPI_Send(&vec2, n, MPI_INT, i, 3, MPI_COMM_WORLD);
			MPI_Recv(&pom, 1, MPI_INT, i, 4, MPI_COMM_WORLD, &st);
			vector[i - 1] = pom;
			printf("Primio sam od procesa %d da je suma %d\n", i,pom);
		}

	}
	else
	{
		MPI_Recv(&n, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &st);		
		MPI_Recv(deoMatrice, n, MPI_INT, 0, 2, MPI_COMM_WORLD, &st);
		MPI_Recv(&vector, n, MPI_INT, 0, 3, MPI_COMM_WORLD, &st);		
		int sum = 0;
		for (int i = 0; i < n; i++)
			sum += vector[i] * deoMatrice[i];		
		MPI_Send(&sum, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);


	}
	
	MPI_Finalize();


}

void matSaljeKolone(int argc, char** argv) // matrica salje kolone, treba 4 procesa RADIIII
{
	int size, rank, n;
	int mat[10][10], vector[5];
	int deoMatrice[10];
	MPI_Status st;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//&mat2[(i-1)*n] salju se vrste 
	if (rank == 0)
	{
		//printf("Primio sam od procesa %d da je suma %d\n", i, pom);
		n = 3;
		int mat2[9] = { 1,2,3,4,5,6,7,8,9 };
		int vec2[3] = { 3,4,5 };
		for (int i = 1; i < size; i++)
		{
			int pom;
			MPI_Send(&n, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			int kolona[10];
			for (int jj = 0; jj < n; jj++)
			{
				kolona[jj] = mat2[(i - 1) +jj* n];
			}
			//kolona{mat2[0], mst2[3] mat2[6] }
			MPI_Comm comm = MPI_COMM_WORLD; //da skratimo vreme pisanja
			MPI_Send(kolona, n, MPI_INT, i, 44, comm);
			//MPI_Send(&mat2[(i-1)*n], n, MPI_INT, i, 2, MPI_COMM_WORLD); // Ovde salje vrste
			MPI_Send(&vec2, n, MPI_INT, i, 3, MPI_COMM_WORLD);
			MPI_Recv(&pom, 1, MPI_INT, i, 4, MPI_COMM_WORLD, &st);
			
			vector[i - 1] = pom;
			//printf("Primio sam od procesa %d da je suma %d\n", i, pom);
		}
		

	}
	else
	{
		MPI_Recv(&n, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &st);
		MPI_Recv(deoMatrice, n, MPI_INT, 0, 44, MPI_COMM_WORLD, &st);
		//MPI_Recv(deoMatrice, n, MPI_INT, 0, 2, MPI_COMM_WORLD, &st);
		MPI_Recv(&vector, n, MPI_INT, 0, 3, MPI_COMM_WORLD, &st);
		printf("ja sam %d: kolonaMat< %d %d %d >", rank, deoMatrice[0], deoMatrice[1], deoMatrice[2]);
		int sum = 0;
		for (int i = 0; i < n; i++)
			sum += vector[i] * deoMatrice[i];
		MPI_Send(&sum, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);


	}

	MPI_Finalize();


}


void ispitKocka(int argc, char** argv)//28.4.2017 MNOGO DOBRO RADI
{
	int size, rank;

	MPI_Status st;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	MPI_Request request;
	
	int brIter = log2(size);
	int niz[8];

	// 1. nacin da odredis kome saljes je : poww=pow(2,brIteracije-iter-1) i petlja od 0 do brIteracije
	// 2. je samo unazad 
	if (rank == 0)
	{
		for (int i = 0; i < 8; i++)
			niz[i] = i ;
	}

	for (int iter = brIter-1; iter >= 0; iter--)
	{
		int fleg;
		int poww = pow(2, iter);
		int komeSaljem = rank ^ poww; //0iter: 1-2, 3-4  sad : 
		MPI_Isend(&niz[poww], poww, MPI_INT, komeSaljem, 0, comm, &request);	
		MPI_Test(&request, &fleg, &st);
		/*if (rank > komeSaljem)
			MPI_Recv();*/
			//MPI_Irecv(niz, poww, MPI_INT, komeSaljem, 0, comm, &request);
		if (fleg)
		{
			if (rank > komeSaljem)
				MPI_Irecv(niz, poww, MPI_INT, komeSaljem, 0, comm, &request);
			if (rank < komeSaljem) {
				MPI_Irecv(&niz[poww], poww, MPI_INT, komeSaljem, 0, comm, &request);
			}

		}

		MPI_Wait(&request, &st);

		

		//ovo je da izbrisemo one elemente koje smo poslali 
		
		if (iter == 2) {
			printf(" Ja sam %d, niz koji imam je : ", rank);
			for (int i = 0; i < 8; i++)
				printf("%d ", niz[i]);
			printf("\n");
		}

	}
	//printf("KONACNO ja sam %d i ja imam %d \n", rank, niz[0]);
	/*printf(" Ja sam %d, niz koji imam je : ", rank);
	for (int i = 0; i < 8; i++)
		printf("%d ", niz[i]);
	printf("\n")*/

	/*
	
		int poww = pow(2, iter);
		int komeSaljem = rank ^ poww; //0iter: 1-2, 3-4  sad : 
		MPI_Isend(&niz[poww], poww, MPI_INT, komeSaljem, 0, comm, &request);
		if(rank>komeSaljem) 
			MPI_Irecv(niz, poww, MPI_INT, komeSaljem, 0, comm, &request);

		MPI_Wait(&request, &st);*/




	MPI_Finalize();
}



void lab7Milosija(int argc, char** argv) // MNOGO DOBRO RADI
{

	int size, rank;
	int sum = 0;
	MPI_Status st;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	MPI_Request request;
	sum = rank;

	int brIteracije = log2(size);
	for (int it = 0; it < brIteracije; it++)
	{

		
			int pomSum = 0;
			// i=0 , razlika je 1, odnosno 2^0=1, tj 2^it
			// i=1, razlika je 2, odnosno 2^1
			
			// i=0 rank%2==1 oni salju\ =1
			// i=1 rank%2==0  // 3->1 i 7->5 =2
			// i=2	rank%2==0			// 4->0 =4
			int poww = pow(2, it);
			if (rank % poww == 0) // ovo omogucava da aktiviramo samo one procese koji treba da rader
				// tj u 0toj iteraciji se aktiviraju svi, u 1oj 0,2,4,6 i 2goj 0,4 
			{

				//odredimo ko nam je par tj u 1oj iteraciji (0,1),(2,3),(4,5) u 2goj je (0,2) i (4,6)
				//trecoj (0,6)
				int par = rank ^ poww; //(0,2) i (4,6)
				// i:0 par= rank ^1 -> 0
				// i:1 par= rank^2
				if (rank < par)
					MPI_Recv(&pomSum, 1, MPI_INT, par, 0, comm, &st);
				else
					MPI_Send(&sum, 1, MPI_INT, par, 0, comm);
				sum += pomSum;
				printf("%d Ja sam proces %d, moja suma je %d \n", it, rank, sum);
			}
		
			
	}


	MPI_Finalize();
}


/*

*/
void lab6a(int argc, char** argv)
{

	int size, rank ;
	
	MPI_Status st;
	MPI_Comm comm=MPI_COMM_WORLD;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	MPI_Request request;
	 int p = size, k = rank;
	int brIter = log2(p);
	
	//int* nizZaSumu = (int*)malloc((10) * sizeof(int));
	int nizZaSumu[100];
	int nizZaSlanje[8]; //max elemnta = 4??
	nizZaSlanje[0] = rank; //okrugle
	nizZaSumu[0] = rank; //kockaste
	int j = 0;
	for (int iter = 0; iter < brIter; iter++)
	{
		int poww = pow(2, iter);
		int l = k^poww  ; //xor

		MPI_Isend(&nizZaSlanje[0], poww, MPI_INT, l, 0, comm, &request);

		MPI_Irecv(&nizZaSlanje[poww], poww, MPI_INT, l, 0, comm, &request);

		

		MPI_Wait(&request, &st);
		
		
	}
	for (int index = 0; index < 8; index++)
	{
		if (rank >= nizZaSlanje[index])
		{
			nizZaSumu[j] = nizZaSlanje[index];
			j++;
		}
	}
	
	
	printf(" ja sam %d Niz za slanje: ", rank);
	for (int i = 0; i < 8; i++)
		printf("%d ", nizZaSlanje[i]);

	printf(" ja sam %d Niz za sumu: ", rank);
	for (int i = 0; i < rank; i++)
		printf("%d ", nizZaSumu[i]);
	




	MPI_Finalize();
}



/*realizacija sumiranja 
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			s+=i+j
ravnomernom ciklicnom raspodelom posla izmedju p procesa. k je indeks procesa,
tj Pk(0<=k<=p-1) treba da izvrsi k,k+p, k+2p ,... , k+N^2-p sumiranje po redu u sekvencijalnom izvrsenju programa.
Pretpostaviti da je N vece od broja procesa p i da je N deljivo sa p.
Rezultat prikazati u procesu koji sadrzi najmanji broj sabiraka koji su prosti brojevi.
Nije dozvoljeno koriscenje indeksiranim promenljivih
*/
void ispit_jun2020(int argc, char** argv) //?? nzm sta trazi uopst 
{

	int rank, size;
	MPI_Status st;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);




}



/*
p=2
mat:		vector = 
1,2,3,4		9
5,6,7,8		6
			4
			3
*/
void ispit_oktobar2019(int argc, char** argv)
{
	//size=4
	//mat je nxk, vector je kx1
	int rank, size, k = 8, n = 16;
	MPI_Status st;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int mat[16][8];
	int vec[8], pomNiz[2], kolonaZaRoot[16];
	int s = k / size;

	if (rank == 0)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < k; j++)
				mat[i][j] = j ;
		for (int j = 0; j < k; j++)
			vec[j] = j;
	}
	MPI_Scatter(vec, s, MPI_INT, pomNiz, s, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		for (int kprocesu = 1; kprocesu < size; kprocesu++) // p1, p2
		{
			// svakom saljes po 2 kol -> p1: 0,1 p2:2,3
			for (int index = 0; index < n; index++)
			{
				kolonaZaRoot[index] = mat[index][(kprocesu - 1) * s];

			}
			MPI_Send(kolonaZaRoot, n, MPI_INT, kprocesu, 1, MPI_COMM_WORLD);
		}
	}
	else {
		MPI_Recv(kolonaZaRoot, n, MPI_INT, 0, 1, MPI_COMM_WORLD,&st);
	}

	printf(" \n Ja sam %d \n", rank);
	for (int z = 0; z < n; z++)
		printf(" %d ", kolonaZaRoot[z]);

	MPI_Finalize();
}


/*
Mnozenje matrice A(nxk) i vektoraB(k) dobija se vektor C
Pronalazi max vrednost u A
Proizvod elemenata svake vrste mat A
Master salje svakom procesu po s kolona matrice A (s je zadata konst) i po s elemenata vektora B.
Svi elementi kolone matrice A se salju odjednom. 
Svi procesi ucestvuju u izracunavanjima. Rezultati se prikazuju u procesu koji sadrzi max vrednost 
elemenata u matrici A, nakon raspodele kolona po procesima
Resiti pomocu grupnih operacija, osim deo za slanje kolona matrice A za koje se koristi point-to-point

*/
void ispit_avgust2019(int argc, char** argv)  //by Milos ali ne radi lepo treba da izbaci 2016
{
	//sabrca brojeve od 1 do 64 

	//n=8, 64 sabiranja 
	//p=size= 4, to znaci proces P3 tj k=3 treba da odradi sabiranja: 3ce,3+4=7, 3+8=11 ,..., 3+64-4=63
	//0 do 64 sabiranja,poslednji je 63-njega treba da odradi P3(jer su procesi od P0 do P3)
	//P1 radi sve do 1+64-4=61
	//i+j..s+=i+j..s += 3+4=7...s+=3+8=7+11=18
	int rank, size,N=8,k,broj,s=0,prostiBrojevi=0, result;
	MPI_Status st;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	bool prost = true;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{

			//i*N je prvih 8, drugih 8 itd 
			k = i * N + j; //1,2,3,4,5,6,7,8,9,0, -> 0: k=0 svakom 4tom tj svakom size 
			if ((k - rank) % size == 0) // da li to k treba ja da uzmem
			{
				broj = i + j;
				s += broj;
				/*for (int z = 2; z <= broj / 2; ++z)
				{
					if (broj % z == 0)
						prost = false;
				}
				if (prost)
					prostiBrojevi++;*/
			}
		}
	}
	MPI_Reduce(&s, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0)
		printf("%d", result);
	MPI_Finalize();

}




void isput_avgust2020(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Status st;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);



}




/*
void main(int argc, char** argv)
{

	//svi koji su ovde rade super 

	//zad1_lab(argc, argv);
	//zad2a_lab(argc, argv);
	//zab2b_lab(argc, argv);
	//zad3_lab(argc, argv);
	//zad4_lab(argc, argv);
	//zad5_lab(argc, argv);
	//mojSaLab(argc, argv);
	//matSaljeKolone(argc, argv); // ne salje kolone NE RADI LEPO
	//lab6a(argc, argv);
	//ispitKocka(argc, argv);
	//lab7Milosija(argc, argv);
	//ispit_avgust2019(argc, argv);
	/*int sum = 0;
	for (int i = 0; i < 64; i++)
		sum += i;
	printf("  %d", sum);
	*/

//	ispit_oktobar2019(argc, argv);

//}