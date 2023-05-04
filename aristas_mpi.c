#include <stdlib.h>
#include <mpi.h>

// Estructura que contiene lo necesario para almacenar una imagen en formato PGM
typedef struct
{
	int row;	  // número de filas en la imagen
	int col;	  // número de columnas en la imagen
	int max_gray; // maximo valor gray
	int **matrix; // matriz de pixeles entre 0 y 255
} PGMData;

int **CrearArray2D_int(int, int);
void LiberarArray2D_int(int, double **);
void readPGM(char *, PGMData *);
void writePGM(char *, PGMData *);
void Filtro_Laplace(int **, int **, int, int);
void crea_pgm(int, int, int, int **, PGMData *);

int main(int argc, char **argv)
{
	char archivo_imagen_ori[100] = "logo.pgm";
	char archivo_imagen_aristas[100] = "logo_edge_paralelo.pgm";

	PGMData img_data, img_edge;
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0)
	{
		switch (argc)
		{
		case 2:
			strcpy(archivo_imagen_ori, argv[1]);
			strcat(archivo_imagen_ori, ".pgm");
			strcpy(archivo_imagen_aristas, argv[1]);
			strcat(archivo_imagen_aristas, "_edge.pgm");
			break;
		case 1:
			break;
		default:
			printf("Demasiados parametros\n");
			return 0;
		}

		printf("\n  *************** DATOS DE LA EJECUCION ***************************\n");
		printf("  * Archivo imagen original   : %25s         *\n", archivo_imagen_ori);
		printf("  * Archivo imagen con aristas: %25s         *\n", archivo_imagen_aristas);
		printf("  *****************************************************************\n\n");
		printf("  Leyendo imagen \"%s\" ... \n", archivo_imagen_ori);

		readPGM(archivo_imagen_ori, &img_data);
		printf("  Dimension de la imagen: %d x %d\n", img_data.row, img_data.col);
	}

	// Enviamos los datos de la imagen a todos los procesos
	MPI_Bcast(&img_data.row, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&img_data.col, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&img_data.max_gray, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// A cada proceso enviamos una parte de la imagen usando MPI_Send.
	// El proceso 0 tendrá las filas que le corresponden más una fila extra por abajo
	// y el último proceso tendrá las filas que le corresponden más una fila extra por arriba
	// Los demás procesos tendrán las filas que les corresponden más una fila extra por arriba y otra por abajo
	int filas_por_proceso = img_data.row / size;
	int indice = filas_por_proceso + img_data.row % size;
	int rowlocal = 0;


	if (rank == 0)
	{
		rowlocal = indice + 1;
		for (int i = 1; i < size - 1; i++)
		{
			
			for (int j = (indice - 1); j < (indice + filas_por_proceso + 1); j++)
			{
				MPI_Send(img_data.matrix[j], img_data.col, MPI_INT, i, 0, MPI_COMM_WORLD);
			}
			indice += filas_por_proceso;
		}
		for (int i = (indice - 1); i < img_data.row; i++)
		{
			MPI_Send(img_data.matrix[i], img_data.col, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
		}
		img_edge.matrix = CrearArray2D_int(img_data.row, img_data.col);
	}
	else
	{
		if (rank == size - 1) {
			rowlocal = filas_por_proceso + 1;
		}
		else {
			rowlocal = filas_por_proceso + 2;
		}
		img_data.matrix = CrearArray2D_int(rowlocal, img_data.col);
		for (int i = 0; i < rowlocal; i++)
		{
			MPI_Recv(img_data.matrix[i], img_data.col, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		img_edge.matrix = CrearArray2D_int(rowlocal, img_data.col);
	}

	printf("  Aplicando el filtro de Laplace proceso %d ...\n", rank);
	img_edge.row = img_data.row;
	img_edge.col = img_data.col;
	img_edge.max_gray = img_data.max_gray;

	Filtro_Laplace(img_data.matrix, img_edge.matrix, rowlocal, img_data.col);

	if (rank != 0 && rank != size - 1) {
		for (int i = 1; i < rowlocal - 1; i++) {
			MPI_Send(img_edge.matrix[i], img_data.col, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		LiberarArray2D_int(rowlocal, img_data.matrix);
		LiberarArray2D_int(rowlocal, img_edge.matrix);
	}
	else if (rank == size - 1) {
		for (int i = 1; i < rowlocal; i++) {
			MPI_Send(img_edge.matrix[i], img_data.col, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		LiberarArray2D_int(rowlocal, img_data.matrix);
		LiberarArray2D_int(rowlocal, img_edge.matrix);
	}
	else {
		// Recibimos los datos de los demás procesos
		indice = filas_por_proceso + img_data.row % size;
		for (int i = 1; i < size; i++) {
			if (i != size - 1) {
				for (int j = indice; j < indice + filas_por_proceso; j++) {
					MPI_Recv(img_edge.matrix[j], img_data.col, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			}
			else {
				for (int j = indice; j < img_data.row; j++) {
					MPI_Recv(img_edge.matrix[j], img_data.col, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			}
			indice += filas_por_proceso;
		}
		printf("  Guardando la imagen con la detección de aristas en \"%s\"\n\n", archivo_imagen_aristas);
		writePGM(archivo_imagen_aristas, &img_edge);
		LiberarArray2D_int(img_data.row, img_data.matrix);
		LiberarArray2D_int(img_edge.row, img_edge.matrix);
	}
	MPI_Finalize();
	return 0;
}