#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include"spice.h"

struct dc_sweep* createNode_dc_sweep(char type_dc, char* var, double start, double end, double incr, char* node_p) {
    struct dc_sweep* newnode = (struct dc_sweep*)malloc(sizeof(struct dc_sweep));

	newnode->var = (char*)malloc((strlen(var)+1)*sizeof(char));
	newnode->node_p = (char*)malloc((strlen(node_p)+1)*sizeof(char));

	if (newnode == NULL) {
        printf("Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    strcpy(newnode->node_p,node_p);
	newnode->type_dc = type_dc;
	strcpy(newnode->var, var);
	newnode->start = start;
    newnode->end = end;
	newnode->incr = incr;
	newnode->next = NULL;
    return newnode;
}

void create_dc_sweep(struct dc_sweep** head, char type_dc,char* var, double start, double end, double incr, char* node_p) {

	struct dc_sweep* newnode = createNode_dc_sweep(type_dc,var,start,end,incr,node_p);

	if (*head == NULL) {
        *head = newnode;
    } else {
        struct dc_sweep* current = *head;
        while (current->next != NULL) {
            current = current->next;
        }
        current->next = newnode;
    }
}

int printDc_sweep(struct dc_sweep* dc,struct list** head,struct HashTable *t,struct circuit* c) {
	struct dc_sweep* current2 = dc;
	double epsilon = 1e-10;
	FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
	int first=0;
	double *x = NULL;
	double *y = NULL;

	
	free_circ(t,c);

    while (current2 != NULL){
		first = 0;
		char* str = (char*) malloc(sizeof(char*)*(strlen(current2->node_p)+8));
		int i,j;
		strcpy(str,"Results/");
		for(i=0;i<(strlen(current2->node_p));i++){
			str[i+8] = current2->node_p[i];
		}
		str[i+8] = '\0';
		str = (char*) realloc(str, (strlen(current2->node_p)+strlen(current2->var)+18)*sizeof(char*));

		str[i+8] = '_';
		str[i+9] ='d';
		str[i+10] = 'c';
		str[i+11] = '_';
		str[i+12] = current2->type_dc;
		for(j=0;j<(strlen(current2->var));j++){
			str[j+i+13] = current2->var[j];
		}
		str[j+i+13] = '\0';
		str[j+i+13] = '.';
		str[j+i+14] = 't';
		str[j+i+15] = 'x';
		str[j+i+16] = 't';
		str[j+i+17] = '\0';

		double start = current2->start;
		int w=0;
		FILE* file = NULL;
		while(start< current2->end + epsilon){
			double temp=0;
			struct list* current1 = *head;
			while (current1 != NULL){
				if(current1->type=='v' || current1->type=='i'){
					if(current1->type==current2->type_dc && strcmp(current1->name,current2->var)==0){
						temp = current1->value;
						current1->value = start;
						break;
					}
				}
				current1 = current1->next;
			}
			
			create_tables(head,t,c);
			
			if(choleski == 1 && iter == 1){
				if(w==0){
					file = fopen(str,"w");
				}
				if(sparse==1){
					create_sparse_t(head,c,t);
					sparse_conjugate_gradients(c,t);
				}
				else{
					conjugate_gradients(c,t);
				}

			}
			else if(iter == 1){

				if(w==0){
					file = fopen(str,"w");
				}

				if(sparse==1){
					create_sparse_t(head,c,t);
					sparse_bi_conjugate_gradients(c,t);
				}
				else{
					bi_conjugate_gradients(c,t);
				}
			}
			else if(choleski==1){
				if(w==0){
					file = fopen(str,"w");
				}
				
				if(sparse==1){
					create_sparse_t(head,c,t);
					sparse_chol_solve(c,t);
				}
				else{
					Choleski(c,t);
					back_for_sub(c,t);
				}
			}
			else{
				if(w==0){
					file = fopen(str,"w");
				}

				if(sparse==1){
					create_sparse_t(head,c,t);
					sparse_lu_solve(c,t);
				}
				else if(diode != 0){
					newtons_method(c,t,head);
				}
				else{
					LU_func(c,t);
					back_for_sub(c,t);
				}
			}
			int index;

			index = hash_find(t,current2->node_p)-1;

			fprintf(file,"%f V(%s): %.6e V\n",start,current2->node_p,c->x[index]);

			if(first==0){
				x = (double *) malloc(1*sizeof(double));
				y = (double *) malloc(1*sizeof(double));
				x[first] = start;
				y[first] = c->x[index];
				first++;
			}
			else{
				x = (double *) realloc(x, (first+1)*sizeof(double));
				y = (double *) realloc(y, (first+1)*sizeof(double));
				x[first] = start;
				y[first] = c->x[index];
				first++;
			}

			current1->value = temp;

			free_circ(t,c);

			start = start + current2->incr;
			w++;
		}

        fprintf(gnuplotPipe, "set title 'V(%s)'\n", current2->node_p);
		fprintf(gnuplotPipe, "set xlabel 'X-Axis A'\n");
		fprintf(gnuplotPipe, "set ylabel 'Y-Axis V'\n");

		fprintf(gnuplotPipe, "plot '-' with lines linewidth 2 linecolor rgb 'blue', \
			'-' with points pt 7 ps 1 linecolor rgb 'red' notitle\n");

		for (int i = 0; i < first; i++) {
			fprintf(gnuplotPipe, "%lf %lf\n", x[i], y[i]);
		}
		fprintf(gnuplotPipe, "e\n");
		for (int i = 0; i < first; i++) {
			fprintf(gnuplotPipe, "%lf %lf\n", x[i], y[i]);
		}

		// End of marks data
		fprintf(gnuplotPipe, "e\n");

		fprintf(gnuplotPipe, "set term png\n");
		fprintf(gnuplotPipe, "set output 'Results/%s_for_dc_%c%s.png'\n", current2->node_p, current2->type_dc ,current2->var);
		fprintf(gnuplotPipe, "replot\n");

		fflush(gnuplotPipe);
        free(x);
        free(y);

        current2 = current2->next;
        fclose(file);
        free(str);
    }
	pclose(gnuplotPipe);
    return 1;
}

void clear_dc_sweep(struct dc_sweep** head) {
    struct dc_sweep* current = *head;
    while (current != NULL) {
        struct dc_sweep* next = current->next;
		free(current->var);
		free(current->node_p);
        free(current);
        current = next;
    }
    *head = NULL;
}
