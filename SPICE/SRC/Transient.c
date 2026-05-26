#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include"spice.h"

typedef struct {
	struct circuit *c;
    struct tran_circuit *tran;
    struct list *head;
    struct HashTable *hash;
	double x0;
	double gmin;
	double V;
	int size;
} ode_params;

struct transient* createNode_transient(double time_step, double fin_time, char* node_p) {
    struct transient* newnode = (struct transient*)malloc(sizeof(struct transient));

    newnode->node_p = (char*)malloc((strlen(node_p)+1)*sizeof(char));

	if (newnode == NULL) {
        printf("Memory allocation failed\n");
        exit(EXIT_FAILURE);
    };

    strcpy(newnode->node_p,node_p);
	newnode->time_step = time_step;
    newnode->fin_time = fin_time;
	newnode->next = NULL;
    return newnode;
}

void create_transient(struct transient** head, double time_step, double fin_time, char* node_p) {

	struct transient* newnode = createNode_transient(time_step,fin_time,node_p);

	if (*head == NULL) {
        *head = newnode;
    } else {
        struct transient* current = *head;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = newnode;
    }
}

void clear_transient(struct transient** head) {
    struct transient* current = *head;
    while (current != NULL) {
        struct transient* next = current->next;
        free(current->node_p);
        free(current);
        current = next;
    }
    *head = NULL;
}

double pwl_function(double time, double* pwl, int num_points) {
    int segment = 0;
    while (segment < num_points - 1 && time >= pwl[2 * (segment + 1)]) {
        segment++;
    }

    if (segment < num_points - 1) {
        double t1 = pwl[2 * segment];
        double t2 = pwl[2 * (segment + 1)];
        double i1 = pwl[2 * segment + 1];
        double i2 = pwl[2 * (segment + 1) + 1];

        return i1 + (i2 - i1) * (time - t1) / (t2 - t1);
    } else {
        return pwl[2 * (num_points - 1) + 1];
    }
}


void tran_spec(struct list* node,double time){
	double epsilon = 1e-10;
	if(node->tran!=NULL){
		if(strcmp(node->tran,"exp")==0){
			if(time>=0 && time<=node->exp[2]){
				node->value = node->exp[0];
			}
			else if(time>=node->exp[2] && time<=node->exp[4]){
				node->value = node->exp[0] + ((node->exp[1] - node->exp[0])*(1-exp(-(time-node->exp[2])/node->exp[3])));
			}
			else if(time>=node->exp[4] && time<=fin_time+epsilon){
				node->value = node->exp[0] + ((node->exp[1] - node->exp[0])*(exp(-(time-node->exp[4])/node->exp[5])-exp(-(time-node->exp[2])/node->exp[3])));
			}
			
		}
		else if(strcmp(node->tran,"sin")==0){
			if(time>=0 && time<=node->sin[3]){
				node->value = node->sin[0] + (node->sin[1]*sin((2*M_PI)*node->sin[5]/360));
			}
			else if(time>node->sin[3] && time<=fin_time+epsilon){
				node->value = node->sin[0] + (node->sin[1]*sin((2*M_PI*node->sin[2]*(time-node->sin[3]))+(2*M_PI*node->sin[5]/360))*exp(-(time-node->sin[3])*node->sin[4]));
			}
		}
		else if (strcmp(node->tran, "pulse") == 0) {
            double i1 = node->pulse[0];
            double i2 = node->pulse[1];
            double td = node->pulse[2];
            double tr = node->pulse[3];
            double tf = node->pulse[4];
            double pw = node->pulse[5];
            double per = node->pulse[6];

            double k = floor((double)time / per);
            double time_in_period = (double)time - k * per; 

            if (time_in_period <= td) {
                node->value = i1; 
            } else if (time_in_period <= td + tr) {
                node->value = i1 + (i2 - i1) / tr * (time_in_period - td);
            } else if (time_in_period <= td + tr + pw) {
                node->value = i2;  
            } else if (time_in_period <= td + tr + pw + tf) {
                node->value = i2 + (i1 - i2) / tf * (time_in_period - (td + tr + pw));
            } else if (time_in_period <= td + per) {
                node->value = i1;
            } 
        }
        if (strcmp(node->tran, "pwl") == 0) {
			int num_points = node->pos/2;
			node->value = pwl_function(time,node->pwl, num_points);
		}
	}
}

void update_values(struct list** head,struct tran_circuit* tran,struct HashTable* t,double time){
	struct list* curr = *head;

	int k=0;
	for(curr=*head;curr!=NULL;curr=curr->next){

		if(curr->type == 'v' || curr->type == 'l'){
			tran_spec(curr,time);
			tran->b[t->unique_ids+k] = curr->value;
			k++;
		}
		else if(curr->type == 'i'){
			tran_spec(curr,time);
			if(strcmp(curr->node_pos,"0")==0){
				int ind = hash_find(t,curr->node_neg)-1;
				tran->b[ind] = curr->value;
			}
			else if(strcmp(curr->node_neg,"0")==0){
				int ind = hash_find(t,curr->node_pos)-1;
				tran->b[ind] = -curr->value;
			}
			else{
				int ind1 = hash_find(t,curr->node_pos)-1;
				int ind2 = hash_find(t,curr->node_neg)-1;
				tran->b[ind1] = -curr->value;
				tran->b[ind2] = curr->value;
			}
		}

	}
}

void tran_sparse_tables(struct list **head,struct circuit *c,struct tran_circuit* tran,struct HashTable *t){

	long int size = c->m2+t->unique_ids;
	struct list * curr = *head;
    cs *A_s = cs_spalloc(size,size,c->nz,1,1);
	cs *C_s = cs_spalloc(size,size,c->nz,1,1);

	tran->b = (double*)calloc(c->m2+t->unique_ids,sizeof(double));
	tran->x = (double*)calloc(c->m2+t->unique_ids,sizeof(double));

	for(int i=0; i<(c->m2+t->unique_ids); i++){
		tran->b[i] = c->b[i];
		tran->x[i] = c->x[i];
	}
    long int p = 0,tr = 0;
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

			if(curr->type == 'l' && tran_flg ==1){
				C_s->i[tr] = t->unique_ids+k;
				C_s->p[tr] = t->unique_ids+k;
				C_s->p[tr] = -curr->value;
				tr++;
			}

			k++;
		}
		else if(curr->type=='c'){
			if(strcmp(curr->node_pos,"0")==0){
				int ind = hash_find(t,curr->node_neg)-1;
				C_s->i[tr] = ind;
				C_s->p[tr] = ind;
				C_s->x[tr] = curr->value;
				tr++;
			}
			else if(strcmp(curr->node_neg,"0")==0){
				int ind = hash_find(t,curr->node_pos)-1;
				C_s->i[tr] = ind;
				C_s->p[tr] = ind;
				C_s->x[tr] = curr->value;
				tr++;
			}
			else{
				int ind1 = hash_find(t,curr->node_pos)-1;
				int ind2 = hash_find(t,curr->node_neg)-1;
				C_s->i[tr] = ind1;
				C_s->p[tr] = ind1;
				C_s->x[tr] = curr->value;
				tr++;
				C_s->i[tr] = ind2;
				C_s->p[tr] = ind2;
				C_s->x[tr] = curr->value;
				tr++;
				C_s->i[tr] = ind1;
				C_s->p[tr] = ind2;
				C_s->x[tr] = -curr->value;
				tr++;
				C_s->i[tr] = ind2;
				C_s->p[tr] = ind1;
				C_s->x[tr] = -curr->value;
				tr++;
			}
		}
	}
	C_s->nz = tr;
	A_s->nz = p;
    tran->A_sparse = cs_compress(A_s);
	tran->C_sparse = cs_compress(C_s);
    cs_spfree(A_s);
	cs_spfree(C_s);
	cs_dupl(tran->A_sparse);
	cs_dupl(tran->C_sparse);

	if(trapezoidal==0){
		for (int j = 0; j < tran->C_sparse->n; j++) {
			for (int i = tran->C_sparse->p[j]; i < tran->C_sparse->p[j+1]; i++) {
				tran->C_sparse->x[i] /= time_step;
			}
		}
		tran->A_sparse = cs_add(tran->A_sparse,tran->C_sparse, 1, 1);
	}else{
		for (int j = 0; j < tran->C_sparse->n; j++) {
			for (int i = tran->C_sparse->p[j]; i < tran->C_sparse->p[j+1]; i++) {
				tran->C_sparse->x[i] = (2*tran->C_sparse->x[i])/time_step;
			}
		}
		tran->A_sparse = cs_add(tran->A_sparse,tran->C_sparse, 1, 1);
	}

}

void tran_create_tables(struct list** head,struct HashTable* t,struct tran_circuit* tran,struct circuit* c){
	struct list* curr2 = *head;
	tran->A = (double**)calloc(c->m2+t->unique_ids,sizeof(double*));
	tran->old_A = (double**)calloc(c->m2+t->unique_ids,sizeof(double*));
	tran->C = (double**)calloc(c->m2+t->unique_ids,sizeof(double*));

	for(int i=0;i<(c->m2+t->unique_ids);i++){
		tran->A[i]  = (double*)calloc(c->m2 + t->unique_ids, sizeof(double));
		tran->old_A[i]  = (double*)calloc(c->m2 + t->unique_ids, sizeof(double));
		tran->C[i]  = (double*)calloc(c->m2 + t->unique_ids, sizeof(double));
	}
	tran->b = (double*)calloc(c->m2+t->unique_ids,sizeof(double));
	tran->x = (double*)calloc(c->m2+t->unique_ids,sizeof(double));

	for(int i=0; i<(c->m2+t->unique_ids); i++){
		tran->b[i] = c->b[i];
		for(int j=0; j<(c->m2+t->unique_ids); j++){
			tran->old_A[i][j] = c->A[i][j];
			tran->A[i][j] = c->A[i][j];
		}
	}

	int k=0;
	for(curr2=*head;curr2!=NULL;curr2=curr2->next){
		if(curr2->type=='v' || curr2->type=='l'){
			if(curr2->type == 'l' && tran_flg ==1){
				tran->C[t->unique_ids+k][t->unique_ids+k] -= curr2->value;
			}
			k++;
		}
		else if(curr2->type=='c' && tran_flg ==1){
			if(strcmp(curr2->node_pos,"0")==0){
				int ind = hash_find(t,curr2->node_neg)-1;
				tran->C[ind][ind] += curr2->value;
			}
			else if(strcmp(curr2->node_neg,"0")==0){
				int ind = hash_find(t,curr2->node_pos)-1;

				tran->C[ind][ind] += curr2->value;
			}
			else{
				int ind1 = hash_find(t,curr2->node_pos)-1;
				int ind2 = hash_find(t,curr2->node_neg)-1;
				tran->C[ind1][ind1] += curr2->value;
				tran->C[ind2][ind2] += curr2->value;
				tran->C[ind1][ind2] -= curr2->value;
				tran->C[ind2][ind1] -= curr2->value;
			}
		}
	}

	if(trapezoidal==0){
		for(int i=0; i<(c->m2+t->unique_ids); i++){
			tran->x[i] = c->x[i];
			for(int j=0; j<(c->m2+t->unique_ids); j++){
				tran->A[i][j] += (tran->C[i][j]/time_step);
			}
		}
	}
	else{
		for(int i=0; i<(c->m2+t->unique_ids); i++){
			tran->x[i] = c->x[i];
			for(int j=0; j<(c->m2+t->unique_ids); j++){
				tran->A[i][j] += (2*tran->C[i][j])/time_step;
			}
		}
	}
}

void tran_sparse_lu_solve(struct circuit *c,struct tran_circuit* tran, struct HashTable *t){
    int size = c->m2+t->unique_ids;
    css* S;
    csn* N;
    double *y;

    y = (double*)calloc(size,sizeof(double));

    S = cs_sqr(2,tran->A_sparse,0);
    N = cs_lu(tran->A_sparse,S,1);

    cs_ipvec(N->pinv,tran->b,y,size);
    cs_lsolve(N->L,y);
    cs_usolve(N->U,y);
    cs_ipvec(S->q,y,tran->x,size);

    cs_sfree(S);
    cs_nfree(N);
    free(y);
}


void tran_LU_func(struct circuit* c, struct tran_circuit* tran, struct HashTable* t) {
	int size = c->m2 + t->unique_ids,m=0;
	double x=0;

    tran->l = (double**)calloc(size,sizeof(double*));
    tran->u = (double**)calloc(size,sizeof(double*));
    tran->p = (int*)malloc(size * sizeof(int));

    for (int i = 0; i < size; i++) {
        tran->l[i] = (double*)calloc(size,sizeof(double));
        tran->u[i] = (double*)calloc(size,sizeof(double));
        tran->p[i] = i;
    }
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
				tran->u[i][j] = tran->A[i][j];
			}
	}
	for (int i = 0; i < size; i++) {
	}

    for (int k = 0; k < size; k++) {
		x = fabs(tran->u[k][k]);
		m = k;
        for (int i = k+1; i < size; i++) {
            if(fabs(tran->u[i][k]) > x){
				x = fabs(tran->u[i][k]);
				m = i;
			}
        }
		int temp_p = tran->p[k];
		tran->p[k] = tran->p[m];
		tran->p[m] = temp_p;

        double* temp = tran->u[k];
		tran->u[k] = tran->u[m];
		tran->u[m] = temp;
		temp = tran->l[k];
		tran->l[k] = tran->l[m];
		tran->l[m] = temp;

		for (int i = k+1; i < size; i++) {
			tran->l[i][k] = (tran->u[i][k] / tran->u[k][k]);
			for (int j = k+1; j < size; j++) {
				tran->u[i][j] = tran->u[i][j] - (tran->l[i][k]*tran->u[k][j]);
			}
		}
    }
    tran->l[0][0] = 1;

    for(int i=1;i<size;i++){
		for (int j = 0; j <i; j++) {
			tran->u[i][j] = 0;
		}
		tran->l[i][i] = 1;
	}
		//printf("l:\n");
		//printMatrix(tran->l,c->m2+t->unique_ids);
		//printf("u:\n");
		//printMatrix(tran->u,c->m2+t->unique_ids);
}

void tran_back_for_sub(struct circuit* c, struct tran_circuit* tran, struct HashTable* t) {
    int size = c->m2 + t->unique_ids;

    double* y = (double*)malloc(size * sizeof(double));

    if (choleski == 0) {
        for (int i = 0; i < size; i++) {
            y[i] = tran->b[tran->p[i]];
        }
    } else {
        for (int i = 0; i < size; i++) {
            y[i] = tran->b[i];
        }
    }

    for (int k = 0; k < size; k++) {
        for (int j = 0; j < k; j++) {
            y[k] -= tran->l[k][j] * y[j];
        }
        y[k] /= tran->l[k][k];
    }

    for (int k = size - 1; k >= 0; k--) {
        for (int j = k + 1; j < size; j++) {
            y[k] -= tran->u[k][j] * y[j];
        }
        y[k] /= tran->u[k][k];
    }

    for (int i = 0; i < size; i++) {
        tran->x[i] = y[i];
    }

    free(y);
}

void solve_system(struct circuit* c,struct tran_circuit * tran, struct HashTable* t, int size,double b_tk_1[c->m2+t->unique_ids]) {
	
	double temp[c->m2+t->unique_ids];
		
	for(int i=0; i<(c->m2+t->unique_ids);i++){
		temp[i]=0.0;
		temp[i] = tran->x[i];
		//printf("TMP %f\n",temp[i]);
		tran->x[i] = 0;
		//printf("X %f\n",tran->x[i]);
	}
	
	
	
	if(trapezoidal == 0){

		double temp1[c->m2+t->unique_ids];

		if(sparse==1){
			for (int j = 0; j <c->m2+t->unique_ids; j++) {
				temp1[j] = 0;
				for (int p = tran->C_sparse->p[j]; p < tran->C_sparse->p[j + 1]; p++) {
					temp1[tran->C_sparse->i[p]] += tran->C_sparse->x[p] * temp[j];
				}
			}
			for (int i = 0; i < c->m2+t->unique_ids; i++) {
				tran->b[i] += temp1[i];
			}
		}
		else{

			for (int i = 0; i < c->m2+t->unique_ids; i++) {
				temp1[i] = 0;
				for (int j = 0; j < c->m2+t->unique_ids; j++) {
					temp1[i] += tran->C[i][j] * temp[j];
				}
			}

			for (int i = 0; i < c->m2+t->unique_ids; i++) {
				tran->b[i] += temp1[i]/time_step;
			}
		}
	}
	else{

		if(sparse==1){
			cs* G_2C_h = cs_spalloc(c->m2+t->unique_ids, c->m2+t->unique_ids,tran->A_sparse->nz, 1, 1);

			for (int i = 0; i < c->m2+t->unique_ids; i++) {
				for (int j = 0; j < c->m2+t->unique_ids; j++) {
					double value_A = 0.0;
					double value_C = 0.0;
					for (int p = tran->A_sparse->p[j]; p < tran->A_sparse->p[j + 1]; p++) {
						if (tran->A_sparse->i[p] == i) {
							value_A = tran->A_sparse->x[p];
							break;
						}
					}
					for (int p = tran->C_sparse->p[j]; p < tran->C_sparse->p[j + 1]; p++) {
						if (tran->C_sparse->i[p] == i) {
							value_C = tran->C_sparse->x[p] / time_step;
							break;
						}
					}
					double value_G_2C_h = value_A - (2 * value_C);

					if (value_G_2C_h != 0.0) {
						cs_entry(G_2C_h, i, j, value_G_2C_h);
					}
				}
			}
			double temp1[c->m2+t->unique_ids];
			for (int j = 0; j < G_2C_h->n; j++) {
				for (int p = G_2C_h->p[j]; p < G_2C_h->p[j + 1]; p++) {
					temp1[G_2C_h->i[p]] += G_2C_h->x[p] * temp[j];
				}
			}
			for (int i = 0; i < c->m2+t->unique_ids; i++) {
				tran->b[i] +=  b_tk_1[i]-temp1[i];
			}
			cs_free(G_2C_h);
		}
		else{
			double G_2C_h[c->m2+t->unique_ids][c->m2+t->unique_ids];

			for(int i = 0; i < size; i++) {
				for(int j = 0; j < size; j++) {
					G_2C_h[i][j] = tran->old_A[i][j] - (2.0 * tran->C[i][j] / time_step);
				}
			}

			// Update the right-hand side vector
			double temp1[size];
			for (int i = 0; i < size; i++) {
				temp1[i] = 0;
				for (int j = 0; j < size; j++) {
					temp1[i] += G_2C_h[i][j] * temp[j];
				}
			}

			for (int i = 0; i < size; i++) {
				tran->b[i] += b_tk_1[i] - temp1[i];
			}

		}
	}
	if (sparse==1){
		tran_sparse_lu_solve(c,tran,t);
		
		printf("---New X----\n");
		
		for (int i = 0; i < c->m2+t->unique_ids; i++) {
			printf("%f\t",tran->x[i]);
		}
		printf("\n");
	}
	else{
		tran_LU_func(c,tran,t);
		tran_back_for_sub(c,tran,t);
	}
}

void free_tran_circ(struct HashTable* t,struct circuit* c, struct tran_circuit* tran){
	if(sparse==0){
		
		for(int i=0;i<(c->m2+t->unique_ids);i++){
			free(tran->A[i]);
			free(tran->old_A[i]);
			free(tran->C[i]);
		}
		free(tran->C);
		free(tran->A);
		free(tran->old_A);
	}
	else{
		cs_spfree(tran->A_sparse);
		cs_spfree(tran->C_sparse);
	}
	free(tran->x);
	free(tran->b);
}

int memristor_compute_F_ode(const gsl_vector *x, void *params, gsl_vector *f) {
    ode_params *p = (ode_params *)params;

    for (int i = 0; i < p->size; i++) {
        gsl_vector_set(f, i, 0.0);
    }

    for (int i = 0; i < p->size; i++) {
        double linear_contrib = 0.0;
        for (int j = 0; j < p->size; j++) {
            linear_contrib += p->tran->A[i][j] * gsl_vector_get(x, j);
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
			p->V = vd;

            if (vd >= 0) {
                current = curr->params.a1 * p->x0 * sinh(curr->params.b * vd);
            } else {
                current = curr->params.a2 * p->x0 * sinh(curr->params.b * vd);
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
        gsl_vector_set(f, i, gsl_vector_get(f, i) - p->tran->b[i]);
    }

    for (int i = 0; i < p->size; i++) {
        double xi = gsl_vector_get(x, i);
        gsl_vector_set(f, i, gsl_vector_get(f, i) + p->gmin * xi);
    }

    /*printf("F\n");
    for (int i = 0; i < p->size; i++) {
        printf("%g\t",gsl_vector_get(f, i));
    }
    printf("\n");*/

    return GSL_SUCCESS;
}

void memristor_newtons_ode(struct circuit *c, struct tran_circuit *tran, struct HashTable *t, struct list **head,double x0, double *V) {
	int size = c->m2 + t->unique_ids;
	ode_params params;
    gsl_set_error_handler_off();

    const double GMIN_START = 1e-3;
    const double GMIN_END   = 1e-9;
    const double GMIN_FACTOR = 10.0;

    double gmin = GMIN_START;
	double best_x[size];
    int converged = 0;

    while (gmin >= GMIN_END) {
		params.c = c;
        params.tran = tran;
        params.hash = t;
        params.head = *head;
		params.x0 = x0;
        params.gmin = gmin;
		params.size = size;

        gsl_vector *x = gsl_vector_alloc(size);
        for (int i = 0; i < size; i++) {
            gsl_vector_set(x, i, 0.0);
        }

        gsl_multiroot_function f = {&memristor_compute_F_ode, size, &params};
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
            tran->x[i] = best_x[i];
        }
		*V = params.V;
    } else {
        printf("Newton failed even with Gmin stepping.\n");
    }
}

double G(double V, struct list* curr) {
    if (V > curr->params.Vp)
        return curr->params.Ap * (exp(V) - exp(curr->params.Vp));
    else if (V < -curr->params.Vn)
        return -curr->params.An * (exp(-V) - exp(curr->params.Vn));
    else
        return 0.0;
}

double fp(double x, struct list* curr) {
    if (x >= curr->params.xp)
        return exp(-curr->params.alphap * (x - curr->params.xp)) * ((curr->params.xp - x) / (1 - curr->params.xp) + 1);
    else
        return 1.0;
}

double fn(double x, struct list* curr) {
    if (x <= (1 - curr->params.xn))
        return exp(curr->params.alphan * (x + curr->params.xn - 1)) * (x / (1 - curr->params.xn));
    else
        return 1.0;
}

double F(double V, double x,struct list* curr) {
    return (curr->params.eta * V >= 0) ? fp(x,curr) : fn(x,curr);
}

// dx/dt = F(V, x) * G(V)
int memristor_ode(double t, const double y[], double dydt[], void *params) {
	ode_params *p = (ode_params *)params;
	double V=0;

	update_values(&p->head,p->tran,p->hash,t);

	memristor_newtons_ode(p->c,p->tran,p->hash,&p->head,p->x0,&V);

	struct list* curr = p->head;
	while (curr != NULL) {
		if (curr->type == 'x') {
			break;
		}
		curr = curr->next;
	}

    dydt[0] = curr->params.eta*F(V, y[0],curr) * G(V,curr);
	p->x0 = y[0];
    return GSL_SUCCESS;
}


void tran_analysis(struct transient** head1,struct list** head,struct circuit* c,struct tran_circuit* tran,struct HashTable* t){
	struct transient* curr = *head1;
	double epsilon = 1e-10;
	FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
	int first=0,plot_size=0;
	double *x = NULL;
	double **y = NULL;
	struct list* this;
	struct list* that;

	for(curr=*head1;curr!=NULL;curr=curr->next){
		
		char* str = (char*) malloc(sizeof(char*)*(strlen(curr->node_p)+8));
		int i;
		strcpy(str,"Results/");
		for(i=0;i<(strlen(curr->node_p));i++){
			str[i+8] = curr->node_p[i];
		}
		str[i+8] = '\0';
		str = (char*) realloc(str, (strlen(curr->node_p)+10)*sizeof(char*));

		str[i+8] = '_';
		str[i+9] ='t';
		str[i+10] = 'r';
		str[i+11] = 'a';
		str[i+12] = 'n';
		str[i+13] = '.';
		str[i+14] = 't';
		str[i+15] = 'x';
		str[i+16] = 't';
		str[i+17] = '\0';
		
		
		curr->file = fopen(str,"w");

		plot_size++;
		free(str);
	}

	if(memristor_c!=0){
		y = (double **) calloc(plot_size+4,sizeof(double*));
		x = (double *) calloc((fin_time/time_step)+1,sizeof(double));
		for(int e=0; e<plot_size+4;e++){
			y[e] = (double *) calloc((fin_time/time_step)+1,sizeof(double));
		}
		int position;

		ode_params myparams;
		double x0;
		struct list* n = *head;
		while (n != NULL) {
			if (n->type == 'x') {
				x0 = n->params.x0;
				that = n;
			}
			if (n->type == 'v') {
				position = hash_find(t, n->node_pos) - 1;
				this = n;
			}
			n = n->next;
		}

		myparams.c = c;
		myparams.tran = tran;
		myparams.head = *head;
		myparams.hash = t;
		myparams.x0 = x0;

		gsl_odeiv2_system sys = {memristor_ode, NULL, 1, &myparams};
		gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, time_step, 1e-8, 1e-8);

		double time = 0.0;
		double xsv = x0;

		while (time < fin_time + 1e-10) {
			plot_size=0;
			x[first] = time;
			for(curr=*head1;curr!=NULL;curr=curr->next){
				int index = hash_find(t,curr->node_p)-1;
				fprintf(curr->file,"%f  %.6e\n",time,tran->x[index]);
				y[plot_size][first] = tran->x[index];
				plot_size++;
			}
			y[plot_size][first] = tran->x[position];
			y[plot_size+1][first] = tran->x[t->unique_ids];
			y[plot_size+2][first] = xsv;
			first++;

			gsl_odeiv2_driver_apply(d, &time, time + time_step, &xsv);
		}

	}
	else{
		y = (double **) calloc(plot_size,sizeof(double*));
		x = (double *) calloc((fin_time/time_step)+1,sizeof(double));
		for(int e=0; e<plot_size;e++){
			y[e] = (double *) calloc((fin_time/time_step)+1,sizeof(double));
		}

		for(double i=0;i<fin_time+epsilon;i+=time_step){

			if(trapezoidal==1){
				double b_tk_1[c->m2+t->unique_ids];

				for(int z=0;z<c->m2+t->unique_ids;z++){
					b_tk_1[z] = tran->b[z];
					tran->b[z] = 0;
				}
				update_values(head,tran,t,i);

				solve_system(c,tran,t,c->m2+t->unique_ids,b_tk_1);

				plot_size=0;
				x[first] = i;
				for(curr=*head1;curr!=NULL;curr=curr->next){
					int index = hash_find(t,curr->node_p)-1;
					fprintf(curr->file,"%f  %.6e\n",i,tran->x[index]);
					y[plot_size][first] = tran->x[index];
					plot_size++;
				}
				first++;

				if (sparse==0){
					for(int i=0;i<(c->m2+t->unique_ids);i++){
						free(tran->u[i]);
						free(tran->l[i]);
					}
					free(tran->l);
					free(tran->u);
					if (choleski ==0){
						free(tran->p);
					}
				}
			}
			else{
				double b_tk_1[c->m2+t->unique_ids];

				update_values(head,tran,t,i);

				solve_system(c,tran,t,c->m2+t->unique_ids,b_tk_1);

				for(curr=*head1;curr!=NULL;curr=curr->next){
					int index = hash_find(t,curr->node_p)-1;
					fprintf(curr->file,"%f  %.6e\n",i,tran->x[index]);
				}
				for(int z=0;z<c->m2+t->unique_ids;z++){
					tran->b[z] = 0;
				}
				if (sparse==0){
					for(int i=0;i<(c->m2+t->unique_ids);i++){
						free(tran->u[i]);
						free(tran->l[i]);
					}
					free(tran->l);
					free(tran->u);
					if (choleski ==0){
						free(tran->p);
					}
				}
			}
		}
	}
	if(memristor_c!=0){

		for(int j=0;j<first;j++){
			if(y[plot_size+1][j]==0){
				y[plot_size+3][j]=0;
			}
			else{
				y[plot_size+3][j] = -y[plot_size+1][j];
			}
		}
	}
	char* time_val = (char*) malloc(3*sizeof(char));
	if (fin_time >= 1) {
		time_val[0] = 's';
		time_val[1] = '\0';
	} else if (fin_time >= 1e-3 && fin_time < 1) {
		for (int i = 0; i < first; i++) {
			x[i] = 1e3*x[i];
		}
		time_val[0] = 'm';
		time_val[1] = 's';
		time_val[2] = '\0';
	} else if (fin_time >= 1e-6 && fin_time < 1e-3) {
		for (int i = 0; i < first; i++) {
			x[i] = 1e6*x[i];
		}
		time_val[0] = 'u';
		time_val[1] = 's';
		time_val[2] = '\0';
	} else if (fin_time >= 1e-9 && fin_time < 1e-6) {
		for (int i = 0; i < first; i++) {
			x[i] = 1e9*x[i];
		};
		time_val[0] = 'n';
		time_val[1] = 's';
		time_val[2] = '\0';
	} else if(fin_time >= 1e-12 && fin_time < 1e-9){
		for (int i = 0; i < first; i++) {
			x[i] = 1e12*x[i];
		}
		time_val[0] = 'p';
		time_val[1] = 's';
		time_val[2] = '\0';
	}

	char* ampere_val = (char*) malloc(3*sizeof(char));
	char* vol_val = (char*) malloc(3*sizeof(char));
	char* xsv_val = (char*) malloc(3*sizeof(char));
	double min_vol = 0;
	double max_vol = 0;
	double min_i = 0;
	double max_i = 0;
	double min_xsv = that->params.x0;
	double max_xsv = 0;
	if(memristor_c!=0){
		for (int i = 0; i < first; i++) {
			double original_val = y[plot_size][i];

			if (original_val < min_vol) min_vol = original_val;
			if (original_val > max_vol) max_vol = original_val;

			original_val = y[plot_size+1][i];

			if (original_val < min_i) min_i = original_val;
			if (original_val > max_i) max_i = original_val;

			original_val = y[plot_size+2][i];

			if (original_val < min_xsv) min_xsv = original_val;
			if (original_val > max_xsv) max_xsv = original_val;

		}
		double max_abs_vol = fabs(min_vol) > fabs(max_vol) ? fabs(min_vol) : fabs(max_vol);
		double max_abs_i = fabs(min_i) > fabs(max_i) ? fabs(min_i) : fabs(max_i);
		double max_abs_xsv = fabs(min_xsv) > fabs(max_xsv) ? fabs(min_xsv) : fabs(max_xsv);
		double min_vol_temp;
		double max_vol_temp;
		min_vol_temp = fabs(min_vol);
		max_vol_temp = fabs(max_vol);
		if(max_abs_vol >= 1){
			vol_val[0] = 'V';
			vol_val[1] = '\0';
			min_vol_temp = min_vol_temp + 0.5;
			max_vol_temp = max_vol_temp + 0.5;
		}
		else if(max_abs_vol >= 1e-3 && max_abs_vol < 1){
			for (int i = 0; i < first; i++) {
				y[plot_size][i] = 1e3*y[plot_size][i];
			}
			vol_val[0] = 'm';
			vol_val[1] = 'V';
			vol_val[2] = '\0';
			min_vol_temp = 1e3*min_vol_temp;
			max_vol_temp = 1e3*max_vol_temp;
			if(min_vol_temp<10){
				min_vol_temp = min_vol_temp + 0.5;
			}
			else if(min_vol_temp>=10 && min_vol_temp<100){
				min_vol_temp = min_vol_temp + 5;
			}
			else if(min_vol_temp>=100 && min_vol_temp<1000){
				min_vol_temp = min_vol_temp + 50;
			}
			if(max_vol_temp<10){
				max_vol_temp = max_vol_temp + 0.5;
			}
			else if(max_vol_temp>=10 && max_vol_temp<100){
				max_vol_temp = max_vol_temp + 5;
			}
			else if(max_vol_temp>=100 && max_vol_temp<1000){
				max_vol_temp = max_vol_temp + 50;
			}
		}
		else if(max_abs_vol >= 1e-6 && max_abs_vol < 1e-3){
			for (int i = 0; i < first; i++) {
				y[plot_size][i] = 1e6*y[plot_size][i];
			}
			vol_val[0] = 'u';
			vol_val[1] = 'V';
			vol_val[2] = '\0';
			min_vol_temp = 1e6*min_vol_temp;
			max_vol_temp = 1e6*max_vol_temp;
			if(min_vol_temp<10){
				min_vol_temp = min_vol_temp + 0.5;
			}
			else if(min_vol_temp>=10 && min_vol_temp<100){
				min_vol_temp = min_vol_temp + 5;
			}
			else if(min_vol_temp>=100 && min_vol_temp<1000){
				min_vol_temp = min_vol_temp + 50;
			}
			if(max_vol_temp<10){
				max_vol_temp = max_vol_temp + 0.5;
			}
			else if(max_vol_temp>=10 && max_vol_temp<100){
				max_vol_temp = max_vol_temp + 5;
			}
			else if(max_vol_temp>=100 && max_vol_temp<1000){
				max_vol_temp = max_vol_temp + 50;
			}
		}
		else if(max_abs_vol >= 1e-9 && max_abs_vol < 1e-6){
			for (int i = 0; i < first; i++) {
				y[plot_size][i] = 1e9*y[plot_size][i];
			}
			vol_val[0] = 'n';
			vol_val[1] = 'V';
			vol_val[2] = '\0';
			min_vol_temp = 1e9*min_vol_temp;
			max_vol_temp = 1e9*max_vol_temp;
			if(min_vol_temp<10){
				min_vol_temp = min_vol_temp + 0.5;
			}
			else if(min_vol_temp>=10 && min_vol_temp<100){
				min_vol_temp = min_vol_temp + 5;
			}
			else if(min_vol_temp>=100 && min_vol_temp<1000){
				min_vol_temp = min_vol_temp + 50;
			}
			if(max_vol_temp<10){
				max_vol_temp = max_vol_temp + 0.5;
			}
			else if(max_vol_temp>=10 && max_vol_temp<100){
				max_vol_temp = max_vol_temp + 5;
			}
			else if(max_vol_temp>=100 && max_vol_temp<1000){
				max_vol_temp = max_vol_temp + 50;
			}
		}
		else if(max_abs_vol >= 1e-12 && max_abs_vol < 1e-9){
			for (int i = 0; i < first; i++) {
				y[plot_size][i] = 1e12*y[plot_size][i];
			}
			vol_val[0] = 'p';
			vol_val[1] = 'V';
			vol_val[2] = '\0';
			min_vol_temp = 1e12*min_vol_temp;
			max_vol_temp = 1e12*max_vol_temp;
			if(min_vol_temp<10){
				min_vol_temp = min_vol_temp + 0.5;
			}
			else if(min_vol_temp>=10 && min_vol_temp<100){
				min_vol_temp = min_vol_temp + 5;
			}
			else if(min_vol_temp>=100 && min_vol_temp<1000){
				min_vol_temp = min_vol_temp + 50;
			}
			if(max_vol_temp<10){
				max_vol_temp = max_vol_temp + 0.5;
			}
			else if(max_vol_temp>=10 && max_vol_temp<100){
				max_vol_temp = max_vol_temp + 5;
			}
			else if(max_vol_temp>=100 && max_vol_temp<1000){
				max_vol_temp = max_vol_temp + 50;
			}
		}
		min_vol = (min_vol < 0) ? -min_vol_temp : min_vol_temp;
		max_vol = (max_vol < 0) ? -max_vol_temp : max_vol_temp;

		if(max_abs_i >= 1){
			ampere_val[0] = 'A';
			ampere_val[1] = '\0';
		}
		else if(max_abs_i >= 1e-3 && max_abs_i < 1){
			for (int i = 0; i < first; i++) {
				y[plot_size+1][i] = 1e3*y[plot_size+1][i];
				y[plot_size+3][i] = 1e3*y[plot_size+3][i];
			}
			ampere_val[0] = 'm';
			ampere_val[1] = 'A';
			ampere_val[2] = '\0';
		}
		else if(max_abs_i >= 1e-6 && max_abs_i < 1e-3){
			for (int i = 0; i < first; i++) {
				y[plot_size+1][i] = 1e6*y[plot_size+1][i];
				y[plot_size+3][i] = 1e6*y[plot_size+3][i];
			}
			ampere_val[0] = 'u';
			ampere_val[1] = 'A';
			ampere_val[2] = '\0';
		}
		else if(max_abs_i >= 1e-9 && max_abs_i < 1e-6){
			for (int i = 0; i < first; i++) {
				y[plot_size+1][i] = 1e9*y[plot_size+1][i];
				y[plot_size+3][i] = 1e9*y[plot_size+3][i];
			}
			ampere_val[0] = 'n';
			ampere_val[1] = 'A';
			ampere_val[2] = '\0';
		}
		else if(max_abs_i >= 1e-12 && max_abs_i < 1e-9){
			for (int i = 0; i < first; i++) {
				y[plot_size+1][i] = 1e12*y[plot_size+1][i];
				y[plot_size+3][i] = 1e12*y[plot_size+3][i];
			}
			ampere_val[0] = 'p';
			ampere_val[1] = 'A';
			ampere_val[2] = '\0';
		}

		if(max_abs_xsv >= 1){
			xsv_val[0] = 'V';
			xsv_val[1] = '\0';
			max_xsv = max_xsv + 0.2;
			if(min_xsv>0.1){
				min_xsv = min_xsv - 0.03;
			}
		}
		else if(max_abs_xsv >= 1e-3 && max_abs_xsv < 1){
			for (int i = 0; i < first; i++) {
				y[plot_size+2][i] = 1e3*y[plot_size+2][i];
			}
			xsv_val[0] = 'm';
			xsv_val[1] = 'V';
			xsv_val[2] = '\0';
			max_xsv = 1e3*max_xsv;
			min_xsv = 1e3*min_xsv;
			if(min_xsv>1 && min_xsv<10){
				min_xsv = min_xsv - 0.3;
			}
			else if(min_xsv>=10 && min_xsv<100){
				min_xsv = min_xsv - 3;
			}
			else if(min_xsv>=100 && min_xsv<1000){
				min_xsv = min_xsv - 30;
			}
			if(max_xsv<10){
				max_xsv = max_xsv + 0.2;
			}
			else if(max_xsv>=10 && max_xsv<100){
				max_xsv = max_xsv + 2;
			}
			else if(max_xsv>=100 && max_xsv<1000){
				max_xsv = max_xsv + 20;
			}
		}
		else if(max_abs_xsv >= 1e-6 && max_abs_xsv < 1e-3){
			for (int i = 0; i < first; i++) {
				y[plot_size+2][i] = 1e6*y[plot_size+2][i];
			}
			xsv_val[0] = 'u';
			xsv_val[1] = 'V';
			xsv_val[2] = '\0';
			max_xsv = 1e6*max_xsv;
			min_xsv = 1e6*min_xsv;
			if(min_xsv>1 && min_xsv<10){
				min_xsv = min_xsv - 0.3;
			}
			else if(min_xsv>=10 && min_xsv<100){
				min_xsv = min_xsv - 3;
			}
			else if(min_xsv>=100 && min_xsv<1000){
				min_xsv = min_xsv - 30;
			}
			if(max_xsv<10){
				max_xsv = max_xsv + 0.2;
			}
			else if(max_xsv>=10 && max_xsv<100){
				max_xsv = max_xsv + 2;
			}
			else if(max_xsv>=100 && max_xsv<1000){
				max_xsv = max_xsv + 20;
			}
		}
		else if(max_abs_xsv >= 1e-9 && max_abs_xsv < 1e-6){
			for (int i = 0; i < first; i++) {
				y[plot_size+2][i] = 1e9*y[plot_size+2][i];
			}
			xsv_val[0] = 'n';
			xsv_val[1] = 'V';
			xsv_val[2] = '\0';
			max_xsv = 1e9*max_xsv;
			min_xsv = 1e9*min_xsv;
			if(min_xsv>1 && min_xsv<10){
				min_xsv = min_xsv - 0.3;
			}
			else if(min_xsv>=10 && min_xsv<100){
				min_xsv = min_xsv - 3;
			}
			else if(min_xsv>=100 && min_xsv<1000){
				min_xsv = min_xsv - 30;
			}
			if(max_xsv<10){
				max_xsv = max_xsv + 0.2;
			}
			else if(max_xsv>=10 && max_xsv<100){
				max_xsv = max_xsv + 2;
			}
			else if(max_xsv>=100 && max_xsv<1000){
				max_xsv = max_xsv + 20;
			}
		}
		else if(max_abs_xsv >= 1e-12 && max_abs_xsv < 1e-9){
			for (int i = 0; i < first; i++) {
				y[plot_size+2][i] = 1e12*y[plot_size+2][i];
			}
			xsv_val[0] = 'p';
			xsv_val[1] = 'V';
			xsv_val[2] = '\0';
			max_xsv = 1e12*max_xsv;
			min_xsv = 1e12*min_xsv;
			if(min_xsv>1 && min_xsv<10){
				min_xsv = min_xsv - 0.3;
			}
			else if(min_xsv>=10 && min_xsv<100){
				min_xsv = min_xsv - 3;
			}
			else if(min_xsv>=100 && min_xsv<1000){
				min_xsv = min_xsv - 30;
			}
			if(max_xsv<10){
				max_xsv = max_xsv + 0.2;
			}
			else if(max_xsv>=10 && max_xsv<100){
				max_xsv = max_xsv + 2;
			}
			else if(max_xsv>=100 && max_xsv<1000){
				max_xsv = max_xsv + 20;
			}
		}
	}

	plot_size=0;
	for(curr=*head1;curr!=NULL;curr=curr->next){
		fprintf(gnuplotPipe, "set title 'V(%s)'\n", curr->node_p);
		fprintf(gnuplotPipe, "set xlabel 'Time (%s)'\n",time_val);
		fprintf(gnuplotPipe, "set ylabel 'Voltage (V)'\n");
		fprintf(gnuplotPipe, "set xrange [0:%lf]\n", x[first - 1]);
		fprintf(gnuplotPipe, "plot '-' with lines notitle linewidth 2 linecolor rgb 'blue'\n");

		for (int i = 0; i < first; i++) {
			fprintf(gnuplotPipe, "%lf %lf\n", x[i], y[plot_size][i]);
		}
		fprintf(gnuplotPipe, "e\n");
		plot_size++;

		fprintf(gnuplotPipe, "set term png\n");
		fprintf(gnuplotPipe, "set output 'Results/%s_for_tran.png'\n", curr->node_p);
		fprintf(gnuplotPipe, "replot\n");
		fprintf(gnuplotPipe, "unset xrange\n");
	}

	if(memristor_c!=0){
		//V IN
		fprintf(gnuplotPipe, "set title 'Graph of V(%s) vs time'\n",this->node_pos);
		fprintf(gnuplotPipe, "set xlabel 'Time (%s)'\n",time_val);
		fprintf(gnuplotPipe, "set ylabel 'Voltage (%s)'\n",vol_val);
		fprintf(gnuplotPipe, "set xrange [0:%lf]\n", x[first - 1]);


		fprintf(gnuplotPipe, "plot '-' with lines notitle linewidth 2 linecolor rgb 'blue'\n");

		for (int i = 0; i < first; i++) {
			fprintf(gnuplotPipe, "%lf %lf\n", x[i], y[plot_size][i]);
		}
		fprintf(gnuplotPipe, "e\n");

		fprintf(gnuplotPipe, "set term png\n");
		fprintf(gnuplotPipe, "set output 'Results/V(%s).png'\n",this->node_pos);
		fprintf(gnuplotPipe, "replot\n");
		fprintf(gnuplotPipe, "unset xrange\n");

		//CURRENT
		fprintf(gnuplotPipe, "set title 'Graph of I(V%s) vs time'\n",this->name);
		fprintf(gnuplotPipe, "set xlabel 'Time (%s)'\n",time_val);
		fprintf(gnuplotPipe, "set ylabel 'Current (%s)'\n",ampere_val);
		fprintf(gnuplotPipe, "set xrange [0:%lf]\n", x[first - 1]);


		fprintf(gnuplotPipe, "plot '-' with lines notitle linewidth 2 linecolor rgb 'red'\n");

		for (int i = 0; i < first; i++) {
			fprintf(gnuplotPipe, "%lf %lf\n", x[i], y[plot_size+3][i]);
		}
		fprintf(gnuplotPipe, "e\n");

		fprintf(gnuplotPipe, "set term png\n");
		fprintf(gnuplotPipe, "set output 'Results/I(V%s).png'\n",this->name);
		fprintf(gnuplotPipe, "replot\n");
		fprintf(gnuplotPipe, "unset xrange\n");


		//XSV
		fprintf(gnuplotPipe, "set title 'Graph of V(%s) vs time'\n",that->params.xsv);
		fprintf(gnuplotPipe, "set xlabel 'Time (%s)'\n",time_val);
		fprintf(gnuplotPipe, "set ylabel 'Voltage (%s)'\n",xsv_val);
		fprintf(gnuplotPipe, "set xrange [0:%lf]\n", x[first - 1]);
		fprintf(gnuplotPipe, "set yrange [%lf:%lf]\n",min_xsv,max_xsv);

		fprintf(gnuplotPipe, "plot '-' with lines notitle linewidth 2 linecolor rgb 'green'\n");

		for (int i = 0; i < first; i++) {
			fprintf(gnuplotPipe, "%lf %lf\n", x[i], y[plot_size+2][i]);
		}
		fprintf(gnuplotPipe, "e\n");

		fprintf(gnuplotPipe, "set term png\n");
		fprintf(gnuplotPipe, "set output 'Results/V(%s).png'\n",that->params.xsv);
		fprintf(gnuplotPipe, "replot\n");
		fprintf(gnuplotPipe, "unset xrange\n");
		fprintf(gnuplotPipe, "unset yrange\n");

		//I-V
		fprintf(gnuplotPipe, "set title 'Graph of I(V%s)-V(%s)'\n",this->name, this->node_pos);
		fprintf(gnuplotPipe, "set xlabel 'Voltage (%s)'\n",vol_val);
		fprintf(gnuplotPipe, "set ylabel 'Current (%s)'\n",ampere_val);
		fprintf(gnuplotPipe, "set xrange [%lf:%lf]\n", min_vol ,max_vol);

		fprintf(gnuplotPipe, "plot '-' with lines notitle linewidth 2 linecolor rgb 'black'\n");

		for (int i = 0; i < first; i++) {
			fprintf(gnuplotPipe, "%lf %lf\n", y[plot_size][i], y[plot_size+3][i]);
		}
		fprintf(gnuplotPipe, "e\n");

		fprintf(gnuplotPipe, "set term png\n");
		fprintf(gnuplotPipe, "set output 'Results/I(V%s)-V(%s).png'\n",this->name, this->node_pos);
		fprintf(gnuplotPipe, "replot\n");
		fprintf(gnuplotPipe, "unset xrange\n");
	}

	free(time_val);
	free(ampere_val);
	free(xsv_val);
	free(vol_val);
	fflush(gnuplotPipe);
	if(memristor_c!=0){
		for(int e=0;e<plot_size+4;e++){
			free(y[e]);
		}
	}
	else{
		for(int e=0;e<plot_size;e++){
			free(y[e]);
		}
	}
	free(x);
	free(y);
	
	for(curr=*head1;curr!=NULL;curr=curr->next){
		fclose(curr->file);
	}
	
	free_tran_circ(t,c,tran);
	
}
