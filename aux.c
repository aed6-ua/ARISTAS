#include <stdlib.h>
#include <malloc.h>

// Estructira para definir una imagen en tonos de gris en formato PGM
typedef struct {
    int row;
    int col;
    int max_gray;
    int **matrix;
} PGMData;


int **CrearArray2D_int(int, int);

// Funciones para leer y escribir imagenes en formato PGM
#define HI(num) (((num) & 0x0000FF00) >> 8)
#define LO(num) ((num) & 0x000000FF)

void SkipComments(FILE *fp)
{
    int ch;
    char line[100];
 
    while ((ch = fgetc(fp)) != EOF && isspace(ch))
        ;
    if (ch == '#') {
        fgets(line, sizeof(line), fp);
        SkipComments(fp);
    } else
        fseek(fp, -1, SEEK_CUR);
}

void readPGM(char *file_name, PGMData *data)
{
    FILE *pgmFile;
    char version[3];
    int i, j;
    int lo, hi;

    pgmFile = fopen(file_name, "rb");
    if (pgmFile == NULL) {
        perror("No se puede abrir el archivo para leer");
        exit(EXIT_FAILURE);
    }

    fgets(version, sizeof(version), pgmFile);
    if (strcmp(version, "P5")) {
        fprintf(stderr, "Tipo de archivo incorrecto!\n");
        exit(EXIT_FAILURE);
    }

    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->col);
    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->row);
    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->max_gray);
    fgetc(pgmFile);
 
    data->matrix = CrearArray2D_int(data->row, data->col);
    if (data->max_gray > 255)
        for (i = 0; i < data->row; ++i)
            for (j = 0; j < data->col; ++j) {
                hi = fgetc(pgmFile);
                lo = fgetc(pgmFile);
                data->matrix[i][j] = (hi << 8) + lo;
            }
    else
        for (i = 0; i < data->row; ++i)
            for (j = 0; j < data->col; ++j) {
                lo = fgetc(pgmFile);
                data->matrix[i][j] = lo;
            }
 
    fclose(pgmFile);
    return;
 
}
 
void writePGM(char *filename, PGMData *data)
{
    FILE *pgmFile;
    int i, j;
    int hi, lo;
 
    pgmFile = fopen(filename, "wb");
    if (pgmFile == NULL) {
        perror("No se puede abrir el archivo para escribir");
        exit(EXIT_FAILURE);
    }
 
    fprintf(pgmFile, "P5 ");
    fprintf(pgmFile, "%d %d ", data->col, data->row);
    fprintf(pgmFile, "%d ", data->max_gray);
 
    if (data->max_gray > 255) {
        for (i = 0; i < data->row; ++i) {
            for (j = 0; j < data->col; ++j) {
                hi = HI(data->matrix[i][j]);
                lo = LO(data->matrix[i][j]);
                fputc(hi, pgmFile);
                fputc(lo, pgmFile);
            }
 
        }
    } else {
        for (i = 0; i < data->row; ++i)
            for (j = 0; j < data->col; ++j) {
                lo = LO(data->matrix[i][j]);
                fputc(lo, pgmFile);
            }
    }
 
    fclose(pgmFile);
}

// Aplica el filtro de Laplace a una matriz de enteros (imagen_in) y la devuelve en otra matriz de enteros (imagen_out)
// row y col son, respectivamente el numero de filas y columnas de la imagen
void Filtro_Laplace(int **imagen_in, int **imagen_out, int row, int col)
{
	int i,j,suma;

	for (i=0; i<row;i++) {
  		for (j=0; j<col; j++){
			if (i==0 || i==row-1)
			        suma = 255;
			else if (j==0 || j==col-1)
				suma = 255;
			else
				suma = 8*imagen_in[i][j] -
					 imagen_in[i-1][j-1] - imagen_in[i-1][j] - imagen_in[i-1][j+1] -
					 imagen_in[i][j-1] - imagen_in[i][j+1] -
					 imagen_in[i+1][j-1] - imagen_in[i+1][j] - imagen_in[i+1][j+1];
			if (suma < 0) suma = 0;
			if (suma > 255) suma = 255;
			imagen_out[i][j] = suma;
		}
       }
	return; 
}

// Dimensiona un array 2D de enteros de dimension Filas x Columnas
int **CrearArray2D_int(int Filas, int Columnas)
{
	// crea un array de 2 dimensiones
 	int **Array2D;
	Array2D = (int **) malloc( Filas * sizeof(int *) );
	if (Array2D == NULL)
	{
		perror("Problemas al asignar memoria");
		exit(EXIT_FAILURE);
	}
	for (int i=0;i<Filas;i++)
	{
		Array2D[i] = (int *) malloc( Columnas * sizeof(int) );
		if (Array2D[i] == NULL)
	        {
		        perror("Problemas al asignar memoria");
		        exit(EXIT_FAILURE);
	        }
	}
 	return Array2D;
}

// Libera memoria de un array 2D de dimension Filas x Columnas
void LiberarArray2D_int(int Filas, double **Array2D)
{
	// Liberamos la memoria asignada a cada una de las filas 
        for ( int i = 0; i < Filas; i++ )
		free(Array2D[i]);
	//i Liberamos la memoria asignada al array de punteros 
	free(Array2D);
	return;
}

