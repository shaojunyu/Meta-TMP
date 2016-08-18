#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hash.h"
#include <errno.h>
#include <unistd.h>
#include <math.h>

typedef struct col_data {
    char col_name[30];
    int col_index;
    float abd_array[300];
    int col_size;
    struct col_data *next;
} col_data;

void col_insert(col_data **head, col_data *node) {
    if (*head == NULL) {
        *head = node;
        return;
    }
    // printf("%s\n", (*head)->col_name);
    col_data *tmp;
    tmp = *head;
    while (tmp->next) {
        tmp = tmp->next;
    }
    tmp->next = node;
}

float compute_variance(float *array, int size) {
    float dist_sum = 0;
    float mean = 0;
    for (int i = 0; i < size; ++i) {
        mean = mean + array[i];
    }
    mean = mean / size;

    for (int i = 0; i < size; ++i) {
        dist_sum += (array[i] - mean) * (array[i] - mean);
    }
    dist_sum = dist_sum / (size - 1);
    return dist_sum;
}

float compute_covariance(float *array1, int size1, float *array2, int size2) {
    float cov = 0;

    if (size1 != size2) {
        printf("%s\n", "size must be the same!!");
        return 0;
    }

    float m1 = 0, m2 = 0;
    for (int i = 0; i < size1; ++i) {
        m1 += array1[i];
        m2 += array2[i];
    }
    m1 = m1 / size1;
    m2 = m2 / size2;
    // printf("%f\n", m1);
    // printf("%f\n", m2);
    for (int i = 0; i < size1; ++i) {
        cov += ((array1[i] - m1) * (array2[i] - m2));
    }
    cov = cov / (size1 - 1);
    return cov;
}

float determinant(float *martrix, int row, int col) {
    if (row != col) {
        printf("col and row must be the same!!\n");
        return 0;
    }

    for (int i = 0; i < row - 1; ++i) {
        for (int j = 0; j < col - 1; ++j) {
        }
    }
    return 0;
}

double **array_to_covariance_matrix(float *array1, float *array2, float *array3,
                                    int size) {
    if (array3 == NULL) {
        double **m = (double **)malloc(sizeof(double *) * 2);
        for (int i = 0; i < 2; ++i) {
            m[i] = (double *)malloc(sizeof(double) * 2);
        }

        m[0][0] = compute_variance(array1, size);
        m[0][1] = compute_covariance(array1, size, array2, size);
        m[1][0] = compute_covariance(array2, size, array1, size);
        m[1][1] = compute_variance(array2, size);

        return m;
    } else {
        double **m = (double **)malloc(sizeof(double *) * 3);
        for (int i = 0; i < 3; ++i) {
            m[i] = (double *)malloc(sizeof(double) * 3);
        }

        m[0][0] = compute_variance(array1, size);
        m[0][1] = compute_covariance(array1, size, array2, size);
        m[0][2] = compute_covariance(array1, size, array3, size);
        m[1][0] = compute_covariance(array2, size, array1, size);
        m[1][1] = compute_variance(array2, size);
        m[1][2] = compute_covariance(array2, size, array3, size);
        m[2][0] = compute_covariance(array3, size, array1, size);
        m[2][1] = compute_covariance(array3, size, array2, size);
        m[2][2] = compute_variance(array3, size);
        return m;
    }
}

void free_matrix(double **m) {}

double Determinant(double **a, int n) {
    int i, j, j1, j2;  // general loop and matrix subscripts
    double det = 0;    // init determinant
    double **m = NULL; // pointer to pointers to implement 2d
                       // square array

    if (n < 1) {
    } // error condition, should never get here

    else if (n == 1) { // should not get here
        det = a[0][0];
    }

    else if (n == 2) { // basic 2X2 sub-matrix determinate
                       // definition. When n==2, this ends the
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1]; // the recursion series
    }

    // recursion continues, solve next sub-matrix
    else {       // solve the next minor by building a
                 // sub matrix
        det = 0; // initialize determinant of sub-matrix

        // for each column in sub-matrix
        for (j1 = 0; j1 < n; j1++) {
            // get space for the pointer list
            m = (double **)malloc((n - 1) * sizeof(double *));

            for (i = 0; i < n - 1; i++)
                m[i] = (double *)malloc((n - 1) * sizeof(double));
            //     i[0][1][2][3]  first malloc
            //  m -> +  +  +  +   space for 4 pointers
            //       |  |  |  |          j  second malloc
            //       |  |  |  +-> _ _ _ [0] pointers to
            //       |  |  +----> _ _ _ [1] and memory for
            //       |  +-------> _ a _ [2] 4 doubles
            //       +----------> _ _ _ [3]
            //
            //                   a[1][2]
            // build sub-matrix with minor elements excluded
            for (i = 1; i < n; i++) {
                j2 = 0; // start at first sum-matrix column position
                        // loop to copy source matrix less one column
                for (j = 0; j < n; j++) {
                    if (j == j1)
                        continue; // don't copy the minor column element

                    m[i - 1][j2] =
                        a[i][j]; // copy source element into new sub-matrix
                                 // i-1 because new sub-matrix is one row
                                 // (and column) smaller with excluded minors
                    j2++;        // move to next sub-matrix column position
                }
            }

            det += pow(-1.0, 1.0 + j1 + 1.0) * a[0][j1] * Determinant(m, n - 1);
            // sum x raised to y power
            // recursively get determinant of next
            // sub-matrix which is now one
            // row & column smaller

            for (i = 0; i < n - 1; i++)
                free(m[i]); // free the storage allocated to
                            // to this minor's set of pointers
            free(m);        // free the storage for the original
                            // pointer to pointer
        }
    }
    return (det);
}

double compute_MI(float *array1, float *array2, int size) {
    double MI = compute_variance(array1, size) * compute_variance(array2, size);
    MI = MI /
         Determinant(array_to_covariance_matrix(array1, array2, NULL, size), 2);
    MI = 0.5 * log(MI);
    return MI;
}

double compute_CMI(float *array1, float *array2, float *array3,
                   int size) { // CMI(array1, array2 | array3)
    double CMI =
        Determinant(array_to_covariance_matrix(array1, array3, NULL, size), 2) *
        Determinant(array_to_covariance_matrix(array2, array3, NULL, size), 2);
    CMI = CMI / compute_variance(array3, size) /
          Determinant(array_to_covariance_matrix(array1, array2, array3, size),
                      3);
    CMI = 0.5 * log(CMI);
    return CMI;
}

int main(int argc, char const *argv[]) {
    FILE *fp;
    hash_t *Sample_hash = hash_new();

    fp = fopen(argv[1], "r");
    if (fp == NULL) {
        printf("%s \n", strerror(errno));
        return 0;
    }

    char line[1000000];
    col_data *head = NULL;
    while (fgets(line, 1000000, fp)) {
        line[strlen(line) - 1] = '\0';
        // printf("%s \n", strerror(errno));
        if (strstr(line, "Sample")) {
            char *c;
            c = (char *)line;
            c = strtok(c, "\t");
            c = strtok(NULL, "\t");
            int col_index = 1;

            while (c != NULL) {
                col_data *node = (col_data *)malloc(sizeof(col_data));
                node->next = NULL;
                strcpy(node->col_name, c);
                node->col_index = col_index;
                node->col_size = 0;
                col_insert(&head, node);

                col_index++;
                c = strtok(NULL, "\t");
            }
        } else {
            char *c;
            c = (char *)line;
            c = strtok(c, "\t");
            c = strtok(NULL, "\t");
            int col_index = 1;
            while (c != NULL) {
                col_data *tmp = head;
                while (tmp && tmp->col_index != col_index) {
                    tmp = tmp->next;
                }
                tmp->abd_array[tmp->col_size] = atof(c);
                tmp->col_size++;
                col_index++;
                c = strtok(NULL, "\t");
            }
        }
    }
    // return 0;
    // printf("%s\n", head->col_name);
    // PCA reduce
    hash_t *MI_hash = hash_new();
    col_data *tmp = head;
    while (tmp) {
        col_data *next;
        if (tmp->next) {
            next = tmp->next;
            while (next) {
                if (compute_MI(tmp->abd_array, next->abd_array, tmp->col_size) >
                    0.0000009) {
                    char *key = (char *)malloc(100 * sizeof(char));
                    sprintf(key, "%s\t%s", tmp->col_name, next->col_name);
                    double *MI = (double *)malloc(sizeof(double));
                    *MI = compute_MI(tmp->abd_array, next->abd_array,
                                     tmp->col_size);
                    hash_set(MI_hash, key, MI);
                }

                next = next->next;
            }
        }
        tmp = tmp->next;
    }

    fclose(fp);
    char filename[100];
    sprintf(filename, "%s-0-order.out", argv[1]);
    fp = fopen(filename,"w");
    hash_each(MI_hash,{
    	double *MI = (double *)val;
    	fprintf(fp, "%s\t%f\n", key, *MI);
    });
    fclose(fp);

    int order = 1;
    int hash_diff_size = hash_size(MI_hash);
    while (order <= 10 && hash_diff_size > 0) {
        printf("%d-order PCA reducing...\n", order);
        //pca reduce
        hash_each(MI_hash, {
            double *MI = (double *)val;
            col_data *col_1 = NULL;
            col_data *col_2 = NULL;
            col_data *tmp = head;
            while (tmp) {
                if (strstr(key, tmp->col_name) == key) {
                    col_1 = tmp;
                }
                if (strstr(key, tmp->col_name) - key > 0) {
                    col_2 = tmp;
                }
                tmp = tmp->next;
            }
            tmp = head;
            while (tmp) {
                if (strcmp(tmp->col_name, col_1->col_name) != 0 &&
                    strcmp(tmp->col_name, col_2->col_name)) {
                    double CMI;
                    CMI = compute_CMI(col_1->abd_array, col_2->abd_array,
                                      tmp->abd_array, col_1->col_size);
                    if (CMI == 0) {
                        hash_del(MI_hash, (char*)key);
                        printf("%s\t%s\t%f\n", key, tmp->col_name, CMI);
                    }
                }
                tmp = tmp->next;
            }
        });
        hash_diff_size = hash_diff_size - hash_size(MI_hash);
        printf("%d edges reduced\n", hash_diff_size);

        if (hash_diff_size > 0) {
        	
        	sprintf(filename, "%s-%d-order-reduced.out", argv[1], order);
            fp = fopen(filename, "w");
            hash_each(MI_hash, {
                double *MI = (double *)val;
                fprintf(fp, "%s\t%f\n", key, *MI);
            });
            fclose(fp);
        }

        order++;
    }

    // float a[] = {1, 2, 3, 4, 5, 6, 7};
    // float b[] = {1, 2, 3, 4, 5, 6, 8};
    // float c[] = {4, 1, 6, 7, 8, 9, 10};
    // double** cov;
    // cov = array_to_covariance_matrix(a,b,NULL,7);
    // printf("%f\n", compute_covariance(a,7,b,7));
    // printf("%f\n", cov[0][1]);
    // double MI;
    // MI = compute_MI(a,b,7);
    // printf("MI:%f\n", MI);

    // double CMI = compute_CMI(a, b, c, 7);
    // printf("CMI:%f\n", CMI);
    // return 0;

    // double m[2][2] = {2, 3, 1, 4};
    // double **mm = (double **)malloc(sizeof(double *) * 2);
    // mm[0] = (double *)malloc(sizeof(double) * 2);
    // mm[1] = (double *)malloc(sizeof(double) * 2);
    // mm[0][0] = 2;
    // mm[0][1] = 3;
    // mm[1][0] = 1;
    // mm[1][1] = 4;
    // mm = &m[0];
    // printf("%f\n", m[0][1]);
    // printf("%f\n", mm[0][0]);
    // printf("%f\n", Determinant(mm, 2));
    return 0;
}