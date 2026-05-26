#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "spice.h"
#include <time.h>

void print_sparse_matrix(cs* A) {
    printf("Sparse Matrix A:\n");

    for (int j = 0; j < A->n; ++j) {
        printf("col %d: locations %d to %d\n", j, A->p[j], A->p[j + 1] - 1);

        for (int p = A->p[j]; p < A->p[j + 1]; ++p) {
            int i = A->i[p];
            double x = A->x[p];
            printf("  A[%d][%d] = %f\n", i, j, x);
        }
    }
}

int nz_exists(cs* A,int i,int j,int n){
	for(int z=0;z<n;z++){
		if(A->i[z] == i && A->p[z] == j){
			return z;
		}
	}
	return -1;
}




void create_sparse_t(struct list **head,struct circuit *c, struct HashTable *t){
   long int size = c->m2+t->unique_ids;
	struct list * curr = *head;
    cs *A_s = cs_spalloc(size,size,c->nz,1,1);
	
    long int p = 0;
	int k =0;
	for(curr=*head;curr!=NULL;curr=curr->next){
		if(curr->type=='r'){

			if(strcmp(curr->node_pos,"0")==0){
				int ind = hash_find(t,curr->node_neg)-1;

				A_s->i[p] = ind;
				A_s->p[p] = ind;
				A_s->x[p] = 1/curr->value;
				p++;

			}
			else if(strcmp(curr->node_neg,"0")==0){
				int ind = hash_find(t,curr->node_pos)-1;
				
				A_s->i[p] = ind;
				A_s->p[p] = ind;
				A_s->x[p] = 1/curr->value;
				p++;

			}
			else{
				int ind1 = hash_find(t,curr->node_pos)-1;
				int ind2 = hash_find(t,curr->node_neg)-1;
				
				A_s->i[p] = ind1;
				A_s->p[p] = ind1;
				A_s->x[p] = 1/curr->value;
				p++;

				A_s->i[p] = ind2;
				A_s->p[p] = ind2;
				A_s->x[p] = 1/curr->value;
				p++;

				A_s->i[p] = ind1;
				A_s->p[p] = ind2;
				A_s->x[p] = -1/curr->value;
				p++;

				A_s->i[p] = ind2;
				A_s->p[p] = ind1;
				A_s->x[p] = -1/curr->value;
				p++;

			}
		}
		else if(curr->type=='v' || curr->type == 'l'){
			if(strcmp(curr->node_pos,"0")==0){
				int ind = hash_find(t,curr->node_neg)-1;
	
				A_s->i[p] = t->unique_ids+k;
				A_s->p[p] = ind;
				A_s->x[p] = -1;
				p++;

				A_s->i[p] = ind;
				A_s->p[p] = t->unique_ids+k;
				A_s->x[p] = -1;
				p++;

			}
			else if(strcmp(curr->node_neg,"0")==0){
				int ind = hash_find(t,curr->node_pos)-1;
				
				A_s->i[p] = t->unique_ids+k;
				A_s->p[p] = ind;
				A_s->x[p] = 1;
				p++;

				A_s->i[p] = ind;
				A_s->p[p] = t->unique_ids+k;
				A_s->x[p] = 1;
				p++;

			}
			else{
				int ind1 = hash_find(t,curr->node_pos)-1;
				int ind2 = hash_find(t,curr->node_neg)-1;
				
				A_s->i[p] = t->unique_ids+k;
				A_s->p[p] = ind2;
				A_s->x[p] = -1;
				p++;

				A_s->i[p] = ind2;
				A_s->p[p] = t->unique_ids+k;
				A_s->x[p] = -1;
				p++;

				A_s->i[p] = t->unique_ids+k;
				A_s->p[p] = ind1;
				A_s->x[p] = 1;
				p++;

				A_s->i[p] = ind1;
				A_s->p[p] = t->unique_ids+k;
				A_s->x[p] = 1;
				p++;
			}
			
			k++;
		}
	}
	
	A_s->nz = p;
    c->A_sparse = cs_compress(A_s);
    cs_spfree(A_s);
	cs_dupl(c->A_sparse);

}

void sparse_lu_solve(struct circuit *c, struct HashTable *t){
    long int size = c->m2+t->unique_ids;
    css* S;
    csn* N;
    double *y;

    y = (double*)calloc(size,sizeof(double));
	
    S = cs_sqr(2,c->A_sparse,0);
    N = cs_lu(c->A_sparse,S,1);
    cs_spfree(c->A_sparse);

    cs_ipvec(N->pinv,c->b,y,size);
    cs_lsolve(N->L,y);
    cs_usolve(N->U,y);
    cs_ipvec(S->q,y,c->x,size);

    cs_sfree(S);
    cs_nfree(N);
    free(y);
}

int sparse_chol_solve(struct circuit *c, struct HashTable *t){
    int size = c->m2+t->unique_ids;
    css* S;
    csn* N;
    double *y;

    y = (double*)calloc(size,sizeof(double));

    S = cs_schol(1,c->A_sparse);
    N = cs_chol(c->A_sparse,S);
    cs_spfree(c->A_sparse);
		
	if (N == NULL) {
		free(y);
		cs_sfree(S);
		cs_nfree(N);
		return -1;
	}
	
    cs_ipvec(S->pinv,c->b,y,size);
    cs_lsolve(N->L,y);
    cs_ltsolve(N->L,y);
    cs_pvec(S->pinv,y,c->x,size);

    cs_sfree(S);
    cs_nfree(N);
    free(y);
	return 0;
}

void sparse_preconditioner(cs* A, int size, double *r, double *z) {
    double *m = (double *)calloc(size, sizeof(double));

    for (int j = 0; j < size; j++) {
        for (int p = A->p[j]; p < A->p[j + 1]; p++) {
            if (A->i[p] == j) {
                m[j] = A->x[p];
                break;
            }
        }
    }

    for (int i = 0; i < size; i++) {
        if (m[i] != 0) {
            z[i] = r[i] / m[i];
        } else {
            z[i] = r[i];
        }
    }

    free(m);
}

void sparse_conjugate_gradients(struct circuit* c, struct HashTable* t) {
	int size = c->m2 + t->unique_ids;
    double *r = (double *)calloc(size,sizeof(double*));
    double *p = (double *)calloc(size,sizeof(double*));
    double *z = (double *)calloc(size,sizeof(double*));
    double *Ap = (double *)calloc(size,sizeof(double*));
	double *Ax = (double *)calloc(size,sizeof(double*));
	double *q = (double *)calloc(size,sizeof(double*));
	double rho=0, rho1=0, beta=0, alpha=0, pq=0;

    for (int i = 0; i < size; i++) {
        c->x[i] = 0;
    }

    cs_gaxpy(c->A_sparse, c->x, Ax);

    for (int i = 0; i < size; i++){
        r[i] = c->b[i] - Ax[i];
    }
	double norm_r = 0;
	double normSquared = 0;

    for(int i = 0; i < size; i++){
        normSquared += pow(r[i], 2);
    }
    norm_r = sqrt(normSquared);

	double norm_b = 0;
	normSquared = 0;

    for(int i = 0; i < size;++i){
        normSquared += pow(c->b[i], 2);
    }
    norm_b = sqrt(normSquared);
	if(norm_b == 0){
		norm_b = 1;
	}
    int iteration = 0;
    while ((norm_r/norm_b) > itol && iteration<(size*2)) {
		pq = 0;
		rho = 0;
		sparse_preconditioner(c->A_sparse,size,r,z);
		for (int i = 0; i < size; i++) {
			rho += r[i]*z[i];
		}
		if(iteration == 0){
			for (int i = 0; i < size; i++) {
				p[i] = z[i];
			}
		}
		else{
			beta = rho/rho1;
			for (int i = 0; i < size; i++) {
				p[i] = z[i] + beta*p[i];

			}
		}
		rho1 = rho;
		for (int i = 0; i < size; i++) {
            Ap[i] = 0;
		}
        cs_gaxpy(c->A_sparse, p, Ap);
		for (int i = 0; i < size; i++) {
			q[i] = Ap[i];
		}

		for (int i = 0; i < size; i++) {
			pq += p[i]*q[i];
		}
		alpha = rho/pq;
		for (int i = 0; i < size; i++) {
			c->x[i] += alpha*p[i];
		}
		for (int i = 0; i < size; i++) {
			r[i] -= alpha*q[i];
		}
		normSquared = 0;
		for(int i = 0; i < size; i++){
			normSquared += pow(r[i], 2);
		}
		norm_r = sqrt(normSquared);
		iteration++;

	}
	cs_spfree(c->A_sparse);
    free(r);
    free(p);
    free(z);
    free(Ap);
	free(Ax);
	free(q);

}

void sparse_bi_conjugate_gradients(struct circuit* c, struct HashTable* t) {
    int size = c->m2 + t->unique_ids;
    double *r = (double *)calloc(size,sizeof(double*));
    double *p = (double *)calloc(size,sizeof(double*));
    double *z = (double *)calloc(size,sizeof(double*));
    double *Ap = (double *)calloc(size,sizeof(double*));
	double *ATp = (double *)calloc(size,sizeof(double*));
	double *Ax = (double *)calloc(size,sizeof(double*));
	double *q = (double *)calloc(size,sizeof(double*));
	double *bi_r = (double *)calloc(size,sizeof(double*));
	double *bi_p = (double *)calloc(size,sizeof(double*));
	double *bi_z = (double *)calloc(size,sizeof(double*));
	double *bi_q = (double *)calloc(size,sizeof(double*));
	double rho=0, rho1=0, beta=0, alpha=0, bi_pq=0, omega=0;

    for (int i = 0; i < size; i++) {
        c->x[i] = 0;
        Ax[i] = 0;
    }
    cs_gaxpy(c->A_sparse,c->x,Ax);

    for (int i = 0; i < size; i++){
        r[i] = c->b[i] - Ax[i];
		bi_r[i] = r[i];
    }
	double norm_r = 0;
	double normSquared = 0;

    for(int i = 0; i < size; i++){
        normSquared += pow(r[i], 2);
    }
    norm_r = sqrt(normSquared);

	double norm_b = 0;
	normSquared = 0;

    for(int i = 0; i < size;++i){
        normSquared += pow(c->b[i], 2);
    }
    norm_b = sqrt(normSquared);
	if(norm_b == 0){
		norm_b = 1;
	}
    int iteration = 0;
    while ((norm_r/norm_b) > itol && iteration<(size*2)) {
		bi_pq = 0;
		rho = 0;
		sparse_preconditioner(c->A_sparse,size,r,z);
		sparse_preconditioner(c->A_sparse,size,bi_r,bi_z);
		for (int i = 0; i < size; i++) {
			rho += bi_r[i]*z[i];
		}

		if(fabs(rho) < 1e-14){
			printf("Algorithm Failed\n");
			exit(1);
		}
		if(iteration == 0){
			for (int i = 0; i < size; i++) {
				p[i] = z[i];
				bi_p[i] = bi_z[i];
			}
		}
		else{
			beta = rho/rho1;
			for (int i = 0; i < size; i++) {
				p[i] = z[i] + beta*p[i];
				bi_p[i] = bi_z[i] + beta*bi_p[i];
			}
		}
		rho1 = rho;
		for (int i = 0; i < size; i++) {
			Ap[i] = 0;
		}
		cs_gaxpy(c->A_sparse,p,Ap);

        for (int i = 0; i < size; i++) {
            ATp[i] = 0;
            for (int p = c->A_sparse->p[i]; p < c->A_sparse->p[i + 1]; p++) {
                ATp[i] += c->A_sparse->x[p] * bi_p[c->A_sparse->i[p]];
            }
        }

		for (int i = 0; i < size; i++) {
			q[i] = Ap[i];
			bi_q[i] = ATp[i];
		}

		for (int i = 0; i < size; i++) {
			bi_pq += bi_p[i]*q[i];
		}
		omega = bi_pq;
		if(fabs(omega) < 1e-14){
			printf("Algorithm Failed\n");
			exit(1);
		}
		alpha = rho/omega;
		for (int i = 0; i < size; i++) {
			c->x[i] += alpha*p[i];
		}
		for (int i = 0; i < size; i++) {
			r[i] -= alpha*q[i];
		}
		for (int i = 0; i < size; i++) {
			bi_r[i] -= alpha*bi_q[i];
		}
		normSquared = 0;
		for(int i = 0; i < size; i++){
			normSquared += pow(r[i], 2);
		}
		norm_r = sqrt(normSquared);
		iteration++;
	}
    free(r);
    free(p);
    free(z);
    free(Ap);
	free(ATp);
	free(Ax);
	free(q);
	free(bi_r);
	free(bi_p);
	free(bi_z);
	free(bi_q);
    cs_spfree(c->A_sparse);
}


