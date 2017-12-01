#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

#define GENE_NUMBER 4549
#define DISEASED_MAX 8
#define NORMAL_MAX 52
#define DISTRIBUTION_SIZE 1000
#define FILE_NAME "NCI-60.csv"
#define MASTER 0
#define EOL "\n"
#define SEPARATOR ","


gsl_rng *r;

struct Gene {
    char name[10];
    double d_value;
};

typedef struct Gene Gene;
int compare_genes(const void * value1, const void * value2) {
    Gene *gene1 = (Gene*) value1;
    Gene *gene2 = (Gene*) value2;
    if (gene1->d_value > gene2->d_value) {
        return -1;
    } else if (gene2->d_value > gene1->d_value) {
        return 1;
    } else {
        return 0;
    }
}

double t_stat(double *dataset1, int n1, double *dataset2, int n2) {
    double mean1 = gsl_stats_mean(dataset1, 1, n1);
    double variance1 = gsl_stats_variance_m(dataset1, 1, n1, mean1);

    double mean2 = gsl_stats_mean(dataset2, 1, n2);
    double variance2 = gsl_stats_variance_m(dataset2, 1, n2, mean2);

    return (mean1 - mean2) / sqrt((variance1 / n1) + (variance2 / n2));
}


double* duplicate_array(double *array, int n) {
    int i;
    double *duplicate = malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        duplicate[i] = array[i];
    }

    return duplicate;
}

double* concat_array(double *array1, int n1, double *array2, int n2) {
    int concat_size = n1 + n2;
    double *concatenated = malloc(concat_size * sizeof(double));

    int i;
    for (i = 0; i < n1; i++) {
        concatenated[i] = array1[i];
    }

    for (i = 0; i < n2; i++ ) {
        concatenated[i + n1] = array2[i];
    }

    return concatenated;
}

// Fisher-Yates shuffle
double* shuffle_array(double *array, int n) {
    int i;
    double *shuffled = duplicate_array(array, n);


    for (i = 0; i < n; i++) {
        int j = gsl_rng_get(r) % n;
        double temp = shuffled[i];
        shuffled[i] = shuffled[j];
        shuffled[j] = temp;
    }

    return shuffled;
}

void random_permutation(double *dataset1, int n1, double *dataset2, int n2, double *permutation1, double *permutation2) {
    int i;
    int concat_size = n1 + n2;
    double *concatenated = concat_array(dataset1, n1, dataset2, n2);
    double *shuffled = shuffle_array(concatenated, concat_size);

    for (i = 0; i < n1; i++) {
        permutation1[i] = shuffled[i];
    }

    for (i = 0; i < n2; i++) {
        permutation2[i] = shuffled[i + n1];
    }
}

double random_permutation_t(double *dataset1, int n1, double *dataset2, int n2) {
    double permutation1[n1];
    double permutation2[n2];

    random_permutation(dataset1, n1, dataset2, n2, permutation1, permutation2);

    return t_stat(permutation1, n1, permutation2, n2);
}

int main(int argc, char* argv[]) {
    r = gsl_rng_alloc(gsl_rng_taus);
    if (r == NULL) {
        printf("Could not allocate rng");
        exit(1);
    }

    // Seed the rng with time
    gsl_rng_set(r, time(NULL));

    int my_rank, source, num_nodes;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);

    int genes_per_process = ceil((float)GENE_NUMBER/num_nodes);
    int values_per_gene = (NORMAL_MAX + DISEASED_MAX);
    int total_values = values_per_gene * GENE_NUMBER;
    double *genes_values = malloc(total_values * sizeof(double));
    if (my_rank == MASTER) {
        char * buffer = 0;
        long length;
        FILE *f = fopen(FILE_NAME, "rb");
        if (!f) {
            printf("Could not open file");
            return 1;
        }

        fseek(f, 0, SEEK_END);
        length = ftell(f);
        fseek(f, 0, SEEK_SET);
        buffer = malloc(length + 1);

        if (!buffer) {
            printf("Could not allocate memory for the file");
            return 2;
        }

        fread(buffer, 1, length, f);
        buffer[length] = '\0';
        fclose(f);


        char *line, *value;
        int line_index = 0;
        int value_index = 0;

        //Headers line
        line = strsep(&buffer, EOL);
        while ((line = strsep(&buffer, EOL)) != NULL) {
            char *name = strsep(&line, SEPARATOR);
            while ((value = strsep(&line, SEPARATOR)) != NULL) {
                genes_values[value_index] = strtod(value, NULL);
                value_index++;
            }
            line_index++;
        }

        free(buffer);
    }

    int values_per_process = values_per_gene * genes_per_process;
    int receive_buffer_size = values_per_process * sizeof(double);
    double *receive_buffer = malloc(receive_buffer_size);
    if (!receive_buffer) {
        printf("Could not allocate memory for receive_buffer");
        return 3;
    }

    MPI_Scatter(
        genes_values, values_per_process, MPI_DOUBLE,
        receive_buffer, values_per_process, MPI_DOUBLE,
        MASTER, MPI_COMM_WORLD);

    int i, j;
    for(i = 0; i < genes_per_process; i++) {
        double diseased_genes[DISEASED_MAX];
        int diseased_n = 0;
        for (j = 0; j < DISEASED_MAX; j++) {
            if (receive_buffer[i * values_per_gene + j] != 0.0 ) {
                diseased_genes[diseased_n] = receive_buffer[i * values_per_gene + j];
                diseased_n++;
            }
        }

        double normal_genes[NORMAL_MAX];
        int normal_n = 0;
        for (j = DISEASED_MAX; j < values_per_gene; j++) {
            if (receive_buffer[i * values_per_gene + j] != 0.0) {
                normal_genes[normal_n] = receive_buffer[i * values_per_gene + j];
                normal_n++;
            }
        }

        double distribution[DISTRIBUTION_SIZE];
        double reference_t = t_stat(normal_genes, normal_n, diseased_genes, diseased_n);

        for (i = 0; i< DISTRIBUTION_SIZE; i++) {
            distribution[i] = random_permutation_t(normal_genes, normal_n, diseased_genes, diseased_n);
        }

        double distribution_mean = gsl_stats_mean(distribution, 1, DISTRIBUTION_SIZE);
        double distribution_sd = gsl_stats_sd_m(distribution, 1, DISTRIBUTION_SIZE, distribution_mean);
        double d_value = fabs(reference_t - distribution_mean) / distribution_sd;
    }


    free(genes_values);
    free(receive_buffer);
    MPI_Finalize();

    return 0;
}
