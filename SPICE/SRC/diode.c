#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spice.h"
#include <time.h>

typedef struct{
	struct circuit* c;
	struct list* head;
	struct HashTable *t;
	int size;
}my_params;

int compute_F(const gsl_vector* x, void* params, gsl_vector* f) {
    my_params *p = (my_params *)params;

	for (int i = 0; i <p->size; i++) {
		gsl_vector_set(f,i,0.0);
	}

	for (int i = 0; i < p->size; i++) {
        double linear_contrib = 0.0;
        for (int j = 0; j < p->size; j++) {
            linear_contrib += p->c->A[i][j] * gsl_vector_get(x,j);
        }
        gsl_vector_set(f, i, linear_contrib);
    }

    struct list* curr= p->head;
	while(curr!=NULL){
		double v1;
		double v2;
		double vd;
		double current = 0.0;
		int pos = hash_find(p->t, curr->node_pos) - 1;
		int neg = hash_find(p->t, curr->node_neg) - 1;

		if(pos <0){
			v1 = 0;
			v2 = gsl_vector_get(x,neg);
		}
		else if(neg<0){
			v1 = gsl_vector_get(x, pos);
			v2 = 0;
		}
		else{
			v1 = gsl_vector_get(x, pos);
			v2 = gsl_vector_get(x,neg);
		}

		vd = v1-v2;

        if (curr->type == 'd') {
            current = curr->model.is * (exp(vd / (curr->model.n * curr->model.vt)) - 1);
        }

        if(pos <0){
			gsl_vector_set(f, neg, gsl_vector_get(f,neg) - current);
		}
		else if(neg<0){
			gsl_vector_set(f, pos,gsl_vector_get(f,pos) + current);
		}
		else{
			gsl_vector_set(f, pos, gsl_vector_get(f,pos) + current);
			gsl_vector_set(f, neg, gsl_vector_get(f,neg) - current);
		}

		curr=curr->next;
	}

	for (int i = 0; i <p->size; i++) {
        gsl_vector_set(f, i, gsl_vector_get(f, i) - p->c->b[i]);
    }

	return GSL_SUCCESS;

}

void newtons_method(struct circuit *c, struct HashTable *t, struct list **head) {
    int size = c->m2 + t->unique_ids;
	my_params params;

	gsl_vector *x = gsl_vector_alloc(size);

	 for (int i = 0; i < size; i++) {
        gsl_vector_set(x,i,0.0);
    }

    params.t = t;
	params.c = c;
	params.head = *head;
	params.size = size;

    gsl_multiroot_function f = {&compute_F, size, &params};

	gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_dnewton, size);
	if (solver == NULL) {
		printf("Failed to allocate memory for the GSL solver.\n");
		return;  // Handle the error appropriately
	}

	gsl_multiroot_fsolver_set(solver, &f, x);

	int status;
	size_t iter = 0;

    do{
        iter++;

        status = gsl_multiroot_fsolver_iterate(solver);

		if(status){
			printf("Iteration %zu failed with status: %d\n", iter, status);
			break;
		}

		status = gsl_multiroot_test_residual (solver->f, 1e-7);

    }while(status == GSL_CONTINUE && iter < 1000);

	if (status == GSL_SUCCESS) {
		for (int i = 0; i < size; i++) {
			c->x[i] = 0;
			c->x[i] = gsl_vector_get(solver->x, i);
		}
    } else {
        printf("Solver failed to converge.\n");
    }

    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(x);

}


