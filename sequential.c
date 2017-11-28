#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define GENE_NUMBER 4549
#define DISEASED_MAX 8
#define NORMAL_MAX 52
#define DISTRIBUTION_SIZE 1000
#define FILE_NAME "NCI-60.csv"

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

int main() {
    r = gsl_rng_alloc(gsl_rng_taus);
    if (r == NULL) {
        printf("Could not allocate rng");
        exit(1);
    }

    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    size_t read;
    size_t line_index = 0;
    // Seed the rng with time
    gsl_rng_set(r, time(NULL));

    fp = fopen(FILE_NAME, "r");
    if (fp == NULL) {
        printf("Error opening file");
        exit(2);
    }


    Gene genes[GENE_NUMBER];

    while ((read = getline(&line, &len, fp)) != -1) {
        int normal_n = 0;
        int diseased_n = 0;
        char *normal_tokens[NORMAL_MAX];
        char *diseased_tokens[DISEASED_MAX];

        if (line_index == 0) {
            line_index++;
            continue;
        }

        char *token;
        int token_index = 0;
        Gene gene;

        while ((token = strsep(&line, ",")) != NULL) {
            if (token_index == 0) {
                strncpy(gene.name, token, 10);
                token_index++;
                continue;
            }

            if (strncmp(token, "", 1) == 0) {
                token_index++;
                continue;
            }

            if (token_index <= DISEASED_MAX) {
                diseased_tokens[diseased_n] = token;
                diseased_n++;
            } else {
                normal_tokens[normal_n] = token;
                normal_n++;
            }

            token_index++;
        }

        double diseased_genes[diseased_n];
        double normal_genes[normal_n];
        int i;

        for (i = 0; i < diseased_n; i++) {
            diseased_genes[i] = strtod(diseased_tokens[i], NULL);
        }

        for (i = 0; i < normal_n; i++) {
            normal_genes[i] = strtod(normal_tokens[i], NULL);
        }

        double distribution[DISTRIBUTION_SIZE];
        double reference_t = t_stat(normal_genes, normal_n, diseased_genes, diseased_n);

        for (i = 0; i< DISTRIBUTION_SIZE; i++) {
            distribution[i] = random_permutation_t(normal_genes, normal_n, diseased_genes, diseased_n);
        }

        double distribution_mean = gsl_stats_mean(distribution, 1, DISTRIBUTION_SIZE);
        double distribution_sd = gsl_stats_sd_m(distribution, 1, DISTRIBUTION_SIZE, distribution_mean);
        double d_value = fabs(reference_t - distribution_mean) / distribution_sd;
        gene.d_value = d_value;
        genes[line_index - 1] = gene;

        line_index++;
    }

    fclose(fp);
    if (line) {
        free(line);
    }

    qsort(genes, GENE_NUMBER, sizeof(Gene), compare_genes);
    int i;
    for (i = 0; i < 10; i++) {
        Gene gene = genes[i];
        printf("%s, %.2f\n", gene.name, gene.d_value);
    }

    gsl_rng_free(r);
    return 0;
}
