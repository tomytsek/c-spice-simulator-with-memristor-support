#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include"spice.h"
#include <time.h>

void create_tables(struct list** head,struct HashTable* t,struct circuit* c){
	int k=0;
	struct list* curr2 = *head;
	c->nz=0;
	if(sparse==0){
	
		c->A = (double**)calloc(c->m2+t->unique_ids,sizeof(double*));

		for(int i=0;i<(c->m2+t->unique_ids);i++){
			c->A[i]  = (double*)calloc(c->m2 + t->unique_ids, sizeof(double));
		}
	}
	c->b = (double*)calloc(c->m2+t->unique_ids,sizeof(double));
	c->x = (double*)calloc(c->m2+t->unique_ids,sizeof(double));
	
	for(curr2=*head;curr2!=NULL;curr2=curr2->next){
		if(curr2->type=='r'){

			if(strcmp(curr2->node_pos,"0")==0){
				
				if(sparse==0){
					int ind = hash_find(t,curr2->node_neg)-1;
					c->A[ind][ind] = (1/curr2->value) + c->A[ind][ind];
				}
				else{
					c->nz++;
				}
			}
			else if(strcmp(curr2->node_neg,"0")==0){
				
				if(sparse==0){
					int ind = hash_find(t,curr2->node_pos)-1;
					c->A[ind][ind] = (1/curr2->value) + c->A[ind][ind];
				}
				else{
					c->nz++;
				}
				
			}
			else{
				if(sparse==0){
					int ind1 = hash_find(t,curr2->node_pos)-1;
					int ind2 = hash_find(t,curr2->node_neg)-1;
					c->A[ind1][ind1] = (1/curr2->value) + c->A[ind1][ind1];
					c->A[ind2][ind2] = (1/curr2->value) + c->A[ind2][ind2];
					c->A[ind1][ind2] = c->A[ind1][ind2] - (1/curr2->value);
					c->A[ind2][ind1] = c->A[ind2][ind1] - (1/curr2->value);
				}
				else{
					c->nz = c->nz+4;
				}
			}
		}
		else if(curr2->type=='i'){

			if(strcmp(curr2->node_pos,"0")==0){
				int ind = hash_find(t,curr2->node_neg)-1;
				c->b[ind] = curr2->value + c->b[ind];
			}
			else if(strcmp(curr2->node_neg,"0")==0){
				int ind = hash_find(t,curr2->node_pos)-1;
				c->b[ind] = c->b[ind] - curr2->value;
			}
			else{
				int ind1 = hash_find(t,curr2->node_pos)-1;
				int ind2 = hash_find(t,curr2->node_neg)-1;
				c->b[ind1] = c->b[ind1] - curr2->value;
				c->b[ind2] = curr2->value + c->b[ind2];
			}
		}
		else if(curr2->type=='v' || curr2->type == 'l'){
			if(strcmp(curr2->node_pos,"0")==0){
				if(sparse==0){
					int ind = hash_find(t,curr2->node_neg)-1;
					c->A[t->unique_ids+k][ind] = -1 ;
					c->A[ind][t->unique_ids+k] = -1 ;
				}
				else{
					c->nz = c->nz+2;
				}
			}
			else if(strcmp(curr2->node_neg,"0")==0){
				if(sparse==0){
					int ind = hash_find(t,curr2->node_pos)-1;
					c->A[t->unique_ids+k][ind] = 1 ;
					c->A[ind][t->unique_ids+k] = 1 ;
				}
				else{
					c->nz = c->nz+2;
				}
			}
			else{
				if(sparse==0){
					int ind1 = hash_find(t,curr2->node_pos)-1;
					int ind2 = hash_find(t,curr2->node_neg)-1;
					c->A[t->unique_ids+k][ind2] = -1 ;
					c->A[ind2][t->unique_ids+k] = -1 ;
					c->A[t->unique_ids+k][ind1] = 1 ;
					c->A[ind1][t->unique_ids+k] = 1 ;
				}
				else{
					c->nz = c->nz+4;
				}
			}
			if(curr2->type =='v'){
				c->b[t->unique_ids+k] = curr2->value;
			}
			
			k++;
		}
	}

	printf("---- A ----\n");
	for(int i=0; i<(c->m2+t->unique_ids); i++){
		for(int j=0; j<(c->m2+t->unique_ids); j++){
			printf("%lf\t",c->A[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("---- b ----\n");
	for(int i=0; i<(c->m2+t->unique_ids); i++){
		printf("%lf\t",c->b[i]);
	}
	printf("\n");


}

void LU_func(struct circuit* c, struct HashTable* t) {
   int size = c->m2 + t->unique_ids,m=0;
	double x=0;
	
    c->l = (double**)calloc(size,sizeof(double*));
    c->u = (double**)calloc(size,sizeof(double*));
    c->p = (int*)malloc(size * sizeof(int));

    for (int i = 0; i < size; i++) {
        c->l[i] = (double*)calloc(size,sizeof(double));
        c->u[i] = (double*)calloc(size,sizeof(double));
        c->p[i] = i;
    }
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
				c->u[i][j] = c->A[i][j];
			}
	}
	for (int i = 0; i < size; i++) {
	}
	
    for (int k = 0; k < size; k++) {
		x = fabs(c->u[k][k]);
		m = k;
        for (int i = k+1; i < size; i++) {
            if(fabs(c->u[i][k]) > x){
				x = fabs(c->u[i][k]);
				m = i;
			}
        }
		int temp_p = c->p[k];
		c->p[k] = c->p[m];
		c->p[m] = temp_p;
		
        double* temp = c->u[k];
		c->u[k] = c->u[m];
		c->u[m] = temp;
		temp = c->l[k];
		c->l[k] = c->l[m];
		c->l[m] = temp;
		
		for (int i = k+1; i < size; i++) {
			c->l[i][k] = (c->u[i][k] / c->u[k][k]);
			for (int j = k+1; j < size; j++) {
				c->u[i][j] = c->u[i][j] - (c->l[i][k]*c->u[k][j]);
			}
		}
    }
    c->l[0][0] = 1;
    for(int i=1;i<size;i++){
		for (int j = 0; j <i; j++) {
			c->u[i][j] = 0;
		}
		c->l[i][i] = 1;
	}
}


int Choleski(struct circuit* c, struct HashTable* t){
	int size = c->m2 + t->unique_ids;
	double sum1,sum2;
	
	c->l = (double**)calloc(size,sizeof(double*));
    c->u = (double**)calloc(size,sizeof(double*));

    for (int i = 0; i < size; i++) {
        c->l[i] = (double*)calloc(size,sizeof(double));
        c->u[i] = (double*)calloc(size,sizeof(double));
    }
	
	 for (int k = 0; k < size; k++) {
		 sum1=0;
		 for(int j=0;j<k;j++){
			 sum1 = (c->l[k][j]*c->l[k][j]) + sum1;
		}
		
		if(sum1>c->A[k][k]){
			return -1;
		}
		
		c->l[k][k] = sqrt(c->A[k][k] - sum1);
		 
        for (int i = k+1; i < size; i++) {
			sum2=0;
			for(int j=0;j<k;j++){
				sum2 = (c->l[i][j]*c->l[k][j]) + sum2;
			}
			c->l[i][k] = (1/c->l[k][k])*(c->A[k][i] - sum2);
			
		}
	}
	
	for (int i=0; i<size; i++){
		for (int j=0; j<size; j++){
			c->u[i][j] = c->l[j][i];
		}	
	}
	return 1;
	
}

void back_for_sub(struct circuit* c, struct HashTable* t) {
    int size = c->m2 + t->unique_ids;

    double* y = (double*)malloc(size * sizeof(double));

    if (choleski == 0) {
        for (int i = 0; i < size; i++) {
            y[i] = c->b[c->p[i]];
        }
    } else {
        for (int i = 0; i < size; i++) {
            y[i] = c->b[i];
        }
    }

    for (int k = 0; k < size; k++) {
        for (int j = 0; j < k; j++) {
            y[k] -= c->l[k][j] * y[j];
        }
        y[k] /= c->l[k][k];
    }

    for (int k = size - 1; k >= 0; k--) {
        for (int j = k + 1; j < size; j++) {
            y[k] -= c->u[k][j] * y[j];
        }
        y[k] /= c->u[k][k];
    }

    for (int i = 0; i < size; i++) {
        c->x[i] = y[i];
    }

    free(y);
}

void preconditioner(double** A, int size, double *r, double *z){
	for (int i = 0; i < size; i++){
		if(A[i][i] == 0){
			z[i] = r[i];
		}
		else{
			z[i] = r[i]*(1/A[i][i]);
		}
	}
}

void conjugate_gradients(struct circuit* c, struct HashTable* t) {
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
    for (int i = 0; i < size; i++) {
        Ax[i] = 0;
        for (int j = 0; j < size; j++) {
            Ax[i] += c->A[i][j] * c->x[j];
        }
    }
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
		preconditioner(c->A,size,r,z);
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
			for (int j = 0; j < size; j++) {
				Ap[i] += c->A[i][j] * p[j];
			}
		}

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
    free(r);
    free(p);
    free(z);
    free(Ap);
	free(Ax);
	free(q);
}

void bi_conjugate_gradients(struct circuit* c, struct HashTable* t) {
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
	double** AT = (double**)calloc(size,sizeof(double*));

    for (int i = 0; i < size; i++) {
        AT[i] = (double*)calloc(size,sizeof(double));
    }

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			AT[j][i] = c->A[i][j];
		}
	}

    for (int i = 0; i < size; i++) {
        c->x[i] = 0;
    }
    for (int i = 0; i < size; i++) {
        Ax[i] = 0;
        for (int j = 0; j < size; j++) {
            Ax[i] += c->A[i][j] * c->x[j];
        }
    }
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
		preconditioner(c->A,size,r,z);
		preconditioner(AT,size,bi_r,bi_z);
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
			for (int j = 0; j < size; j++) {
				Ap[i] += c->A[i][j] * p[j];
			}
		}
		for (int i = 0; i < size; i++) {
			ATp[i] = 0;
			for (int j = 0; j < size; j++) {
				ATp[i] += AT[i][j] * bi_p[j];
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
	for(int i = 0; i < size; i++){
		free(AT[i]);
	}
	free(AT);
	
}

void print_dc_op(struct list ** head,struct HashTable* t,struct circuit* c,char* filename){
	FILE* file;
	char* str;
	char *filePath;
	int i;

	const char *lastSlash = strrchr(filename, '/');

    if (lastSlash != NULL) {
        size_t length = strlen(lastSlash + 1);
        char benchmark[length + 1];
        strncpy(benchmark, lastSlash + 1, length);
        benchmark[length] = '\0';
		str = (char*) malloc(sizeof(char*)*(strlen(benchmark)-4));

		for(i=0;i<(strlen(benchmark)-4);i++){
			str[i] = benchmark[i];
		}

		str = (char*) realloc(str, ((strlen(benchmark)-4)+8)*sizeof(char*));

		str[i] = '_';
		str[i+1] ='o';
		str[i+2] = 'p';
		str[i+3] = '.';
		str[i+4] = 't';
		str[i+5] = 'x';
		str[i+6] = 't';
		str[i+7] = '\0';

    } else {
		str = (char*) malloc(sizeof(char*)*(strlen(filename)-4));

		for(i=0;i<(strlen(filename)-4);i++){
			str[i] = filename[i];
		}

		str = (char*) realloc(str, ((strlen(filename)-4)+8)*sizeof(char*));

		str[i] = '_';
		str[i+1] ='o';
		str[i+2] = 'p';
		str[i+3] = '.';
		str[i+4] = 't';
		str[i+5] = 'x';
		str[i+6] = 't';
		str[i+7] = '\0';
    }

	filePath = (char*) malloc(sizeof(char*)*(strlen(str)+8));
	sprintf(filePath, "%s/%s", "Results", str);

	file = fopen(filePath,"w");
	fprintf(file,"g  0.00000e+00\n");

	struct list* curr = *head;
	char ch = '\0';
	for(int i=0;i<t->unique_ids;i++){
		while (curr != NULL) {
            if (hash_find(t,curr->node_neg) == i+1) {
				ch = 'n';
               break;
            }
            else if(hash_find(t,curr->node_pos) == i+1){
				ch = 'p';
				break;
			}
            curr = curr->next;
        }
        
        if(ch == 'p'){
			fprintf(file,"%s  %.6e\n",curr->node_pos,c->x[i]);
		}
        else if(ch == 'n'){
			fprintf(file,"%s  %.6e\n",curr->node_neg,c->x[i]);
		}
		
	}
	if(memristor_c!=0){
		struct list* curre = *head;
		while (curre != NULL) {
			if(curre->type == 'x'){
				ch = 'x';
				break;
			}
			curre = curre->next;
		}
		if(ch == 'x'){
			fprintf(file,"%s  %lf\n",curre->params.xsv,curre->params.x0);
		}
	}

	fclose(file);
	free(str);
	free(filePath);
}

void free_circ(struct HashTable* t,struct circuit* c){
	if(choleski == 1 && iter == 1 && sparse==0){
		for(int i=0;i<(c->m2+t->unique_ids);i++){
			free(c->A[i]);
		}
		free(c->A);
		free(c->x);
		free(c->b);

	}
	else if(iter == 1 && sparse==0){
		for(int i=0;i<(c->m2+t->unique_ids);i++){
			free(c->A[i]);
		}
		free(c->A);
		free(c->x);
		free(c->b);
	}
	else if(sparse == 1){
		free(c->x);
		free(c->b);
	}
	else{
		for(int i=0;i<(c->m2+t->unique_ids);i++){
			free(c->A[i]);
			if(diode==0 && memristor_c==0){
				free(c->u[i]);
				free(c->l[i]);
			}
		}
		if(diode==0 && memristor_c==0){
			free(c->l);
			free(c->u);
			
			if (choleski == 0){
				free(c->p);
			}
			
		}
		free(c->A);
		free(c->x);
		free(c->b);
	}
	
}

void printMatrix(double** matrix,int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%f\t", matrix[i][j]);
        }
        printf("\n");
    }
}
