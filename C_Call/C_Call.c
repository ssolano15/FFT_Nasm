#include <stdio.h>
#include <stdint.h>	
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define BILLION 1000000000L

extern double *fft(const double *in_data, const int size);//Declaración externa de fft_asm

//Documentacion para caso a 4 muestras
double *fftc(const double *in_data, const int size) //Declara el puntero de fft
{
	int i;
	int k = 0;

	while ((1 << k) < size) //Bucle de desplazamiento de bits - Representacion muestras maximo
		k++;
	
	//size = 4 = 100
  	//1er -> 1 para k=0 -> k=1
  	//2do -> 10 para k=1 -> k=2
  	//3er -> 100 = 100 -> k=2

	int *rev = (int*) calloc(size, sizeof(int)); //Crea un espacio en memoria y retorna un puntero rev

	rev[0] = 0; //Iguala la matriz rev en posicion cero a cero
	
	//rev = [0,x,x,x]

	int high1 = -1;

	// Metodo de Bit reversal

	for (i = 1; i < size; ++i) // i<4
	{
		//Cuando i es una potencia de dos highl aumenta en 1
		if ((i & (i - 1)) == 0)		
			high1++; //highl = 1
		
		rev[i] = rev[i ^ (1 << high1)];
		rev[i] |= (1 << (k - high1 - 1));
	}

	//Para i=1 -> highl = 0 -> rev[1]=rev[0] -> rev[1] |= (1<<(2-0-1)) rev[1] |= 2
    	//Para i=2 -> highl = 1 -> rev[2]=rev[0] -> rev[2] |= (1<<(2-1-1)) rev[2] |= 1
   	//Para i=3 -> highl = 1 -> rev[3]=rev[1] -> rev[3] |= (1<<(2-1-1)) rev[3] |= 3

	double *roots = (double*) calloc(2 * size, sizeof(double)); //Crea un espacio en memoria y retorna el puntero roots
	
	// Twiddle Factor
	double alpha = 2 * M_PI / size;	
	for (i = 0; i < size; ++i)
	{
		roots[2 * i] = cos(i * alpha);
		roots[2 * i + 1] = sin(i * alpha);
	}

	double *cur = (double*) calloc(2 * size, sizeof(double)); //Crea un espacio en memoria y retorna el puntero cur

	for (i = 0; i < size; ++i) //Carga y acomodo de los datos en el orden correspondiente. Ejemplo: 0-2 1-3
	{
		int ni = rev[i];
				
		cur[2 * i] = in_data[2 * ni];
		cur[2 * i + 1] = in_data[2 * ni + 1];
	}

	free(rev);

	int len;

	for (len = 1; len < size; len <<= 1) //for de desplazamiento de bits en potencias de dos -> len:1,10
	{
		double *ncur = (double*) calloc(2 * size, sizeof(double)); //Crea un espacio en memoria y retorna el puntero ncur
		
		int rstep = size / (2 * len); //rstep= 2,1
		int p1;
		for (p1 = 0; p1 < size; p1 += len) //Para len=1 p1=0,2 - Para len=2 p1=0
		{
			for (i = 0; i < len; ++i) //i<1,2
			{
				double val_r = roots[2 * (i * rstep)] * cur[2 * (p1 + len)] - roots[2 * (i * rstep) + 1] * cur[2 * (p1 + len) + 1];
				double val_i = roots[2 * (i * rstep)] * cur[2 * (p1 + len) + 1] + roots[2 * (i * rstep) + 1] * cur[2 * (p1 + len)];
				
				ncur[2 * p1] = cur[2 * p1] + val_r;
				ncur[2 * p1 + 1] = cur[2 * p1 + 1] + val_i;
				ncur[2 * (p1 + len)] = cur[2 * p1] - val_r;
				ncur[2 * (p1 + len) + 1] = cur[2 * p1 + 1] - val_i;
				
				p1++;
			}
		}
		
		// Se apuntan los resultados con cur para liberar ncur utilizando un puntero temporal
		double *tmp = ncur;
		ncur = cur;
		cur = tmp;

		free(ncur);
	}	
	free(roots);
	
	return cur;
}

int main(int argc, char* argv[]){
	
	uint64_t diff;
	struct timespec start, end;

	uint64_t diff_asm;
	struct timespec start_asm, end_asm;
	//Generar Matriz ********************************************************************************************
	//***********************************************************************************************************
	
	
	int m_max;
	printf("\nNumero de n para un total de 2^n muestras:");
	scanf ("%d", &m_max);
	printf("\n*****************************************************************************************************************************");
	printf("\nMuestras	|	Tiempo CPU en C (ns)		|	Tiempo CPU en asm (ns)		|	Redimiento del asm (%%)");

	for (int m=1; m <= m_max; ++m){
		
		FILE *in = fopen("Matrices.in", "w");

		fprintf(in,"%d\n", (int) pow (2,m)); 

		for (int cont = 0; cont < (int) pow (2,m); ++cont)
			fprintf(in,"%lF + %lFi\n",(double)rand()/(double)(RAND_MAX/10), (double)rand()/(double)(RAND_MAX/10));
		
		fclose(in);
		
	//***********************************************************************************************************
	//***********************************************************************************************************

		//Guardar Datos en Memoria **********************************************************************************
		//***********************************************************************************************************
		
		in = fopen("Matrices.in", "r"); //Crea un puntero *in hacia las matrices de entrada

		if (in == 0) //Si no contiene datos
		{
			fprintf(stderr, "El archivo 'Matrices.in' esta vacio.\n");
			return 1;
		}
		
		int size, i; //Se definen dos variables
		
	  	fscanf(in, "%d", &size); //Lee la cantidad de muestras del archivo de entrada
		double *in_data = (double*) calloc(2 * size, sizeof(double)); //Crea un espacio en memoria para la entrada y retorna el puntero

		for (i = 0; i < 2 * size; ++i) //Almacena los datos de entrada en la memoria in_data
			fscanf(in, "%lF + %lFi", &in_data[2 * i], &in_data[2 * i + 1]); //Leer los num imaginarios del archivo de entrada

		fclose(in);

		//***********************************************************************************************************	
		//***********************************************************************************************************

		//Imprimir Entradas *****************************************************************************************
		//***********************************************************************************************************		
		
		/*
		printf("\nCantidad de Muestras: %d\n", size);
		printf("Entrada:\n");
		for (int s = 0; s < size; ++s) //Almacena los datos de entrada en la memoria in_data
			printf("%lF + %lFi\n", *(in_data+2*s), *(in_data+2*s+1));
		*/
		
		//***********************************************************************************************************	
		//***********************************************************************************************************

		//Calcular FFT en C *****************************************************************************************
		//***********************************************************************************************************
		
		if (size == 2)
			fftc(in_data, size);

		double *out_datac = (double*) calloc(2 * size, sizeof(double)); //Crea un espacio en memoria para la salida y retorna un puntero
  		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start); //Punto de partida del clock
	  	out_datac = fftc(in_data, size); // Almacena en out_data_c los resultados de la FFT
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end); //Punto final de clock
		
		//Calcular diferencias entre tiempo de inicio y final
		diff = BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
		
		FILE *outc = fopen("fft_c.out", "w"); //Crea un puntero *out hacia el archivo con los resultados

		fprintf(outc, "Resultados de la FFT:\n"); 
		//Imprime en el archivo .out_c 
		for (i = 0; i < size; ++i) //crea un bucle para la impresión de archivos
			fprintf(outc, "%lF + %lFi\n", out_datac[2 * i], out_datac[2 * i + 1]);//Escribe los numeros imaginarios de los resultados en el archivo de salida
			
		fclose(outc);

		
		//***********************************************************************************************************	
		//***********************************************************************************************************

		//Calcular FFT en ASM ***************************************************************************************
		//***********************************************************************************************************

		
		double *out_data = (double*) calloc(2 * size, sizeof(double)); //Crea un espacio en memoria para la salida y retorna un puntero
  	
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_asm); //Punto de partida del clock
	  	out_data = fft(in_data, size); // Almacena en out_data los resultados de la FFT
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_asm); //Punto final de clock
		
		//Calcular diferencias entre tiempo de inicio y final
		diff_asm = BILLION * (end_asm.tv_sec - start_asm.tv_sec) + end_asm.tv_nsec - start_asm.tv_nsec;
		
		double rend_asm = (((float) diff - (float) diff_asm)/ (float) diff_asm)*100;

		FILE *out = fopen("fft_asm.out", "w"); //Crea un puntero *out hacia el archivo con los resultados

		fprintf(out, "Resultados de la FFT:\n"); 
		//Imprime en el archivo .out 
		for (i = 0; i < size; ++i) //crea un bucle para la impresión de archivos
			fprintf(out, "%lF + %lFi\n", out_data[2 * i], out_data[2 * i + 1]);//Escribe los numeros imaginarios de los resultados en el archivo de salida
			
		fclose(out);
		
		
		//***********************************************************************************************************	
		//***********************************************************************************************************


		//Imprimir Resultados C *************************************************************************************
		//***********************************************************************************************************

		/*
		printf("\nResultados en C:\n");
		for (i = 0; i < size; ++i)
			printf("%lF + %lFi\n", *(out_datac+2*i), *(out_datac+2*i+1));
		*/

		//***********************************************************************************************************
		//***********************************************************************************************************

		//Imprimir Resultados asm ***********************************************************************************
		//***********************************************************************************************************

		/*
		printf("\nResultados en asm:\n");
		for (i = 0; i < size; ++i)
			printf("%lF + %lFi\n", *(out_data+2*i), *(out_data+2*i+1));
		*/

		//***********************************************************************************************************
		//***********************************************************************************************************


		//Imprimir Tiempos ***********************************************************************************
		//***********************************************************************************************************
		
		printf("\n%d		|	%llu				|	%llu				|	+%f %%", size, (long long unsigned int) diff, (long long unsigned int) diff_asm, rend_asm);

		//***********************************************************************************************************
		//***********************************************************************************************************
	  	

		//Libera los punteros
		free(in_data);
		free(out_data);
		free(out_datac);

	}

	printf("\n*****************************************************************************************************************************\n");

	return 0;
}