#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "spice.h"
#include <time.h>

typedef struct {
    struct circuit *c;
    struct list *head;
    struct HashTable *hash;
    int size;
    double gmin;
    double x;
} memristor_params;

int memristor_compute_F(const gsl_vector *x, void *params, gsl_vector *f) {
    memristor_params *p = (memristor_params *)params;

    for (int i = 0; i < p->size; i++) {
        gsl_vector_set(f, i, 0.0);
    }

    for (int i = 0; i < p->size; i++) {
        double linear_contrib = 0.0;
        for (int j = 0; j < p->size; j++) {
            linear_contrib += p->c->A[i][j] * gsl_vector_get(x, j);
        }
        gsl_vector_set(f, i, linear_contrib);
    }

    struct list* curr = p->head;
    while (curr != NULL) {
        double v1, v2, current=0.0, vd;
        int pos = hash_find(p->hash, curr->node_pos) - 1;
        int neg = hash_find(p->hash, curr->node_neg) - 1;

        if (pos < 0) {
            v1 = 0;
            v2 = gsl_vector_get(x, neg);
        } else if (neg < 0) {
            v1 = gsl_vector_get(x, pos);
            v2 = 0;
        } else {
            v1 = gsl_vector_get(x, pos);
            v2 = gsl_vector_get(x, neg);
        }

        vd = v1 - v2;
        if (curr->type == 'x') {

            if (vd >= 0) {
                current = curr->params.a1 * p->x * sinh(curr->params.b * vd);
            } else {
                current = curr->params.a2 * p->x * sinh(curr->params.b * vd);
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
        }
        curr = curr->next;
    }

    for (int i = 0; i < p->size; i++) {
        gsl_vector_set(f, i, gsl_vector_get(f, i) - p->c->b[i]);
    }

    for (int i = 0; i < p->size; i++) {
        double xi = gsl_vector_get(x, i);
        gsl_vector_set(f, i, gsl_vector_get(f, i) + p->gmin * xi);
    }

    return GSL_SUCCESS;
}

void memristor_newtons_method(struct circuit *c, struct HashTable *t, struct list **head,double x) {
	int size = c->m2 + t->unique_ids;
	memristor_params params;
    gsl_set_error_handler_off();

    const double GMIN_START = 1e-3;
    const double GMIN_END   = 1e-9;
    const double GMIN_FACTOR = 10.0;

    double gmin = GMIN_START;
	double best_x[size];
    int converged = 0;

    while (gmin >= GMIN_END) {
		params.c = c;
        params.hash = t;
        params.head = *head;
		params.x = x;
        params.gmin = gmin;
		params.size = size;

        gsl_vector *x = gsl_vector_alloc(size);
        for (int i = 0; i < size; i++) {
            gsl_vector_set(x, i, 0.0);
        }

        gsl_multiroot_function f = {&memristor_compute_F, size, &params};
        gsl_multiroot_fsolver *solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_dnewton, size);
        gsl_multiroot_fsolver_set(solver, &f, x);

        int status = GSL_CONTINUE;
        size_t iter = 0;
        const int max_iter = 100;

        while (status == GSL_CONTINUE && iter < max_iter) {
            iter++;
            status = gsl_multiroot_fsolver_iterate(solver);
            if (status) break;

            status = gsl_multiroot_test_residual(solver->f, 1e-10);
        }

        if (status == GSL_SUCCESS) {
			for (int i = 0; i < size; i++) {
                best_x[i] = gsl_vector_get(solver->x, i);
            }
            converged = 1;
            gsl_vector_free(x);
            gsl_multiroot_fsolver_free(solver);
            gmin /= GMIN_FACTOR;
        } else {
            printf("Newton failed at Gmin = %.2e\n", gmin);
            gsl_vector_free(x);
            gsl_multiroot_fsolver_free(solver);
            break;
        }
    }

    if (converged) {
		for (int i = 0; i < size; i++) {
            c->x[i] = best_x[i];
        }
    } else {
        printf("Newton failed even with Gmin stepping.\n");
    }
}

int memristor_ode_op(double t, const double y[], double dydt[], void *params) {
	memristor_params *p = (memristor_params *)params;
    double V=p->c->b[p->hash->unique_ids];

	struct list* curr = p->head;
	while (curr != NULL) {
		if (curr->type == 'x') {
			break;
		}
		curr = curr->next;
	}

    dydt[0] = curr->params.eta*F(V, y[0],curr) * G(V,curr);
    return GSL_SUCCESS;
}

void memristor_op(struct circuit *c, struct HashTable *t, struct list **head) {
    int size = c->m2 + t->unique_ids;
    memristor_params params;
    gsl_set_error_handler_off();

    double x0;
    struct list* n = *head;
    while (n != NULL) {
        if (n->type == 'x') {
            x0 = n->params.x0;
            break;
        }
        n = n->next;
    }

    params.c = c;
    params.head = *head;
    params.hash = t;
    params.size = size;

    gsl_odeiv2_system sys = { memristor_ode_op, NULL, 1, &params };
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-10, 1e-10);

    double max_time = 0.01;
    double time = 0.0;
    double x = x0;
    double dx_thresh = 452.970661 + -456.663947*x + 3.255518*pow(x,2) + 0.500980*pow(x,3) + -0.068750*pow(x,4);
    double dt = 1e-10;
    while (time < max_time) {
        double prev_x = x;
        gsl_odeiv2_driver_apply(d, &time, time + dt, &x);
        double dxdt = (x - prev_x) / dt;
        if (fabs(dxdt) < dx_thresh) break;
    }
    printf("V(xsv)   = %.6f V\n", x);

    memristor_newtons_method(c,t,head,x);

    gsl_odeiv2_driver_free(d);
}

