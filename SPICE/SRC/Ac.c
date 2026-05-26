#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include"spice.h"

struct ac* createNode_ac(char* sweep, double points, double start_freq, double end_freq, char* node_p) {
    struct ac* newnode = (struct ac*)malloc(sizeof(struct ac));

    newnode->sweep = (char*)malloc((strlen(sweep)+1)*sizeof(char));
    newnode->node_p = (char*)malloc((strlen(node_p)+1)*sizeof(char));

	if (newnode == NULL) {
        printf("Memory allocation failed\n");
        exit(EXIT_FAILURE);
    };

    strcpy(newnode->sweep,sweep);
    strcpy(newnode->node_p,node_p);
	newnode->points = points;
    newnode->start_freq = start_freq;
    newnode->end_freq = end_freq;
	newnode->next = NULL;
    return newnode;
}

void create_ac(struct ac** head, char* sweep, double points, double start_freq, double end_freq, char* node_p) {

	struct ac* newnode = createNode_ac(sweep,points,start_freq,end_freq,node_p);

	if (*head == NULL) {
        *head = newnode;
    } else {
        struct ac* current = *head;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = newnode;
    }
}

void clear_ac(struct ac** head) {
    struct ac* current = *head;
    while (current != NULL) {
        struct ac* next = current->next;
        free(current->node_p);
        free(current->sweep);
        free(current);
        current = next;
    }
    *head = NULL;
}

void create_ac_t(struct list** head,struct HashTable* t,struct circuit* c, double f){
	int k=0;
	double omega = 2*M_PI*f;
	struct list* curr2 = *head;
	c->nz=0;
	if(sparse==0){
	
		c->C_A = (double complex**)calloc(c->m2+t->unique_ids,sizeof(double complex*));

		for(int i=0;i<(c->m2+t->unique_ids);i++){
			c->C_A[i]  = (double complex*)calloc(c->m2 + t->unique_ids, sizeof(double complex));
		}
	}
	c->C_b = (double complex*)calloc(c->m2+t->unique_ids,sizeof(double complex));
	c->C_x = (double complex*)calloc(c->m2+t->unique_ids,sizeof(double complex));
	
	for(curr2=*head;curr2!=NULL;curr2=curr2->next){
		if(curr2->type=='r'){

			if(strcmp(curr2->node_pos,"0")==0){
				
				if(sparse==0){
					int ind = hash_find(t,curr2->node_neg)-1;
					c->C_A[ind][ind] = (1/CMPLX(curr2->value,0)) + c->C_A[ind][ind];
				}
				else{
					c->nz++;
				}
			}
			else if(strcmp(curr2->node_neg,"0")==0){
				
				if(sparse==0){
					int ind = hash_find(t,curr2->node_pos)-1;
					c->C_A[ind][ind] = (1/CMPLX(curr2->value,0)) + c->C_A[ind][ind];
				}
				else{
					c->nz++;
				}
				
			}
			else{
				if(sparse==0){
					int ind1 = hash_find(t,curr2->node_pos)-1;
					int ind2 = hash_find(t,curr2->node_neg)-1;
					c->C_A[ind1][ind1] = (1/CMPLX(curr2->value,0)) + c->C_A[ind1][ind1];
					c->C_A[ind2][ind2] = (1/CMPLX(curr2->value,0)) + c->C_A[ind2][ind2];
					c->C_A[ind1][ind2] = c->C_A[ind1][ind2] - (1/CMPLX(curr2->value,0));
					c->C_A[ind2][ind1] = c->C_A[ind2][ind1] - (1/CMPLX(curr2->value,0));
				}
				else{
					c->nz = c->nz+4;
				}
			}
		}
		else if(curr2->type=='i' && curr2->ac!=NULL){

			complex double z = curr2->ac_array[0] * cexp(I*curr2->ac_array[1]);
			if(strcmp(curr2->node_pos,"0")==0){
				int ind = hash_find(t,curr2->node_neg)-1;
				c->C_b[ind] = z + c->C_b[ind];
			}
			else if(strcmp(curr2->node_neg,"0")==0){
				int ind = hash_find(t,curr2->node_pos)-1;
				c->C_b[ind] = z - curr2->value;
			}
			else{
				int ind1 = hash_find(t,curr2->node_pos)-1;
				int ind2 = hash_find(t,curr2->node_neg)-1;
				c->C_b[ind1] = c->C_b[ind1] - z;
				c->C_b[ind2] = z + c->C_b[ind2];
			}
		}
		else if(curr2->type=='v' || curr2->type == 'l'){
			if(strcmp(curr2->node_pos,"0")==0){
				if(sparse==0){
					int ind = hash_find(t,curr2->node_neg)-1;
					c->C_A[t->unique_ids+k][ind] = -1 ;
					c->C_A[ind][t->unique_ids+k] = -1 ;
				}
				else{
					c->nz = c->nz+2;
				}
			}
			else if(strcmp(curr2->node_neg,"0")==0){
				if(sparse==0){
					int ind = hash_find(t,curr2->node_pos)-1;
					c->C_A[t->unique_ids+k][ind] = 1 ;
					c->C_A[ind][t->unique_ids+k] = 1 ;
				}
				else{
					c->nz = c->nz+2;
				}
			}
			else{
				if(sparse==0){
					int ind1 = hash_find(t,curr2->node_pos)-1;
					int ind2 = hash_find(t,curr2->node_neg)-1;
					c->C_A[t->unique_ids+k][ind2] = -1 ;
					c->C_A[ind2][t->unique_ids+k] = -1 ;
					c->C_A[t->unique_ids+k][ind1] = 1 ;
					c->C_A[ind1][t->unique_ids+k] = 1 ;
				}
				else{
					c->nz = c->nz+4;
				}
			}
			if(curr2->type =='v' && curr2->ac!=NULL){
				complex double z = curr2->ac_array[0] * cexp(I*curr2->ac_array[1]);
				c->C_b[t->unique_ids+k] = z;
			}
			
			if(curr2->type =='l'){
				c->C_A[t->unique_ids][k] = -I*omega*curr2->value;
			}
			
			k++;
		}
		else if(curr2->type=='c'){
			if(strcmp(curr2->node_pos,"0")==0){
				int ind = hash_find(t,curr2->node_neg)-1;
				c->C_A[ind][ind] += I*omega*curr2->value;
			}
			else if(strcmp(curr2->node_neg,"0")==0){
				int ind = hash_find(t,curr2->node_pos)-1;
				
				c->C_A[ind][ind] += I*omega*curr2->value;
			}
			else{
				int ind1 = hash_find(t,curr2->node_pos)-1;
				int ind2 = hash_find(t,curr2->node_neg)-1;
				c->C_A[ind1][ind1] += I*omega*curr2->value;
				c->C_A[ind2][ind2] += I*omega*curr2->value;
				c->C_A[ind1][ind2] -= I*omega*curr2->value;
				c->C_A[ind2][ind1] -= I*omega*curr2->value;
			}
		}
	}
	 printf("--A MATRIX--\n");
	for (int i = 0; i < c->m2+t->unique_ids; i++) {
        for (int j = 0; j < c->m2+t->unique_ids; j++) {
            printf("(%f + %fi)\t", creal(c->C_A[i][j]), cimag(c->C_A[i][j]));
        }
        printf("\n");
    }
    printf("\n");
    printf("--B MATRIX--\n");
    for (int i = 0; i < c->m2+t->unique_ids; i++) {
		printf("(%f + %fi)\t", creal(c->C_b[i]), cimag(c->C_b[i]));
	}
	printf("\n");
}

void AC_lu(struct circuit* c,struct HashTable* t){
	int size = c->m2 + t->unique_ids,m=0;
	double x=0;
	
    c->C_l = (double complex**)calloc(size,sizeof(double complex*));
    c->C_u = (double complex**)calloc(size,sizeof(double complex*));
    c->p = (int*)malloc(size * sizeof(int));

    for (int i = 0; i < size; i++) {
        c->C_l[i] = (double complex*)calloc(size,sizeof(double complex));
        c->C_u[i] = (double complex*)calloc(size,sizeof(double complex));
        c->p[i] = i;
    }
    
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
				c->C_u[i][j] = c->C_A[i][j];
			}
	}
	for (int i = 0; i < size; i++) {
	}
	
    for (int k = 0; k < size; k++) {
		x = cabs(c->C_u[k][k]);
		m = k;
        for (int i = k+1; i < size; i++) {
            if(cabs(c->C_u[i][k]) > x){
				x = cabs(c->C_u[i][k]);
				m = i;
			}
        }
		int temp_p = c->p[k];
		c->p[k] = c->p[m];
		c->p[m] = temp_p;
		
        double complex* temp = c->C_u[k];
		c->C_u[k] = c->C_u[m];
		c->C_u[m] = temp;
		temp = c->C_l[k];
		c->C_l[k] = c->C_l[m];
		c->C_l[m] = temp;
		
		for (int i = k+1; i < size; i++) {
			c->C_l[i][k] = (c->C_u[i][k] / c->C_u[k][k]);
			for (int j = k+1; j < size; j++) {
				c->C_u[i][j] = c->C_u[i][j] - (c->C_l[i][k]*c->C_u[k][j]);
			}
		}
    }
    c->C_l[0][0] = 1;
    for(int i=1;i<size;i++){
		for (int j = 0; j <i; j++) {
			c->C_u[i][j] = 0;
		}
		c->C_l[i][i] = 1;
	}
}

void back_AC(struct circuit* c, struct HashTable* t) {
    int size = c->m2 + t->unique_ids;

    double complex* y = (double complex*)malloc(size * sizeof(double complex));

    if (choleski == 0) {
        for (int i = 0; i < size; i++) {
            y[i] = c->C_b[c->p[i]];
        }
    } else {
        for (int i = 0; i < size; i++) {
            y[i] = c->C_b[i];
        }
    }

    for (int k = 0; k < size; k++) {
        for (int j = 0; j < k; j++) {
            y[k] -= c->C_l[k][j] * y[j];
        }
        y[k] /= c->C_l[k][k];
    }

    for (int k = size - 1; k >= 0; k--) {
        for (int j = k + 1; j < size; j++) {
            y[k] -= c->C_u[k][j] * y[j];
        }
        y[k] /= c->C_u[k][k];
    }

    for (int i = 0; i < size; i++) {
        c->C_x[i] = y[i];
    }

    free(y);
}

void free_AC(struct HashTable* t,struct circuit* c){
	
	for(int i=0;i<c->m2+t->unique_ids;i++){
		if(iter == 0){
			free(c->C_l[i]);
			free(c->C_u[i]);
		}
		free(c->C_A[i]);

	}
	if(iter == 0){
		free(c->C_l);
		free(c->C_u);
		free(c->p);
	}
	free(c->C_A);
	free(c->C_b);
	free(c->C_x);
}


void AC_analysis(struct list** head,struct HashTable* t,struct circuit* c,struct ac *head_ac){
	struct ac* curr = head_ac;
	
	
	while (curr != NULL){
		
		char* str = (char*) malloc(sizeof(char*)*(strlen(curr->node_p)+8));
		int i;
		strcpy(str,"Results/");
		for(i=0;i<(strlen(curr->node_p));i++){
			str[i+8] = curr->node_p[i];
		}
		str[i+8] = '\0';
		str = (char*) realloc(str, (strlen(curr->node_p)+8)*sizeof(char*));
		str[i+8] = '_';
		str[i+9] ='a';
		str[i+10] = 'c';
		str[i+11] = '.';
		str[i+12] = 't';
		str[i+13] = 'x';
		str[i+14] = 't';
		str[i+15] = '\0';
		
		
		
		FILE* file;
		file = fopen(str,"w");
		free(str);
		if(strcmp(curr->sweep,"lin")==0){
			
			for(int i=0;i<curr->points;i++){
				double f = curr->start_freq + i * (curr->end_freq - curr->start_freq)/(curr->points - 1);
				
				create_ac_t(head,t,c,f);
				AC_lu(c,t);
				back_AC(c,t);
				int index = hash_find(t,curr->node_p)-1;
				fprintf(file,"%.6e\t%.6e\t%.6e\n",f,cabs(c->C_x[index]),carg(c->C_x[index])*180.0/M_PI);
				
				free_AC(t,c);
				
			}
		}
		else{
			
			for(int i=0;i<curr->points;i++){
				double f = curr->start_freq * exp(i * log(curr->end_freq / curr->start_freq) / (curr->points - 1));
				
				create_ac_t(head,t,c,f);
				AC_lu(c,t);
				back_AC(c,t);
				int index = hash_find(t,curr->node_p)-1;
				fprintf(file,"%.6e\t%.6e\t%.6e\n",f,cabs(c->C_x[index]),carg(c->C_x[index]));
				
				free_AC(t,c);
			}
		}
		curr = curr->next;
		fclose(file);
	}
}

