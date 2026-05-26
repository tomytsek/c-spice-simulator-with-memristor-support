#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include"spice.h"
#include <time.h>

int choleski = 0;
int iter = 0;
double itol = 1e-3;
int sparse = 0;
int trapezoidal = 1;
int backward_euler = 0;
int tran_flg = 0;
int ac_flg = 0;
double time_step=0;
double fin_time=0;
int diode = 0;
int memristor_c = 0;

struct list* createnode(char type,char* name,char* node_pos,char* node_neg,double value,char* tran,double* tran_array, int pos, char* ac,double* ac_array,char* model, char* xsv,char* memristor) {

    struct list* newnode = (struct list*)malloc(sizeof(struct list));

	if (newnode == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

	newnode->name = (char*)malloc((strlen(name)+1)*sizeof(char));
	newnode->node_pos = (char*)malloc((strlen(node_pos)+1)*sizeof(char));
	newnode->node_neg = (char*)malloc((strlen(node_neg)+1)*sizeof(char));


	if (newnode->name == NULL || newnode->node_pos == NULL || newnode->node_neg == NULL) {
        free(newnode->name);
        free(newnode->node_pos);
        free(newnode->node_neg);
        free(newnode);
        printf("Memory allocation failed.\n");
        exit(1);
    }


    newnode->type = type;
	strcpy(newnode->name, name);
	strcpy(newnode->node_pos, node_pos);
	strcpy(newnode->node_neg, node_neg);
	if(type == 'd'){
		newnode->model.modelname = (char*)malloc((strlen(model)+1)*sizeof(char));
		strcpy(newnode->model.modelname,model);
		newnode->model.vt =0;
		newnode->model.is =0;
		newnode->model.n =0;
	}else if(type == 'x'){
		newnode->params.memristorname = (char*)malloc((strlen(memristor)+1)*sizeof(char));
		strcpy(newnode->params.memristorname,memristor);
		newnode->params.xsv = (char*)malloc((strlen(xsv)+1)*sizeof(char));
		strcpy(newnode->params.xsv,xsv);
		newnode->params.a1 =0;
		newnode->params.a2 =0;
		newnode->params.b =0;
		newnode->params.Vp =0;
		newnode->params.Vn =0;
		newnode->params.Ap =0;
		newnode->params.xp =0;
		newnode->params.xn =0;
		newnode->params.alphap =0;
		newnode->params.alphan =0;
		newnode->params.x0 =0;
		newnode->params.eta =0;
	}else{
		newnode->value = value;
	}
    newnode->next = NULL;

	if(tran != NULL){
		if(strcmp(tran, "ac")==0){
			newnode->ac = (char*)malloc((strlen(tran)+1)*sizeof(char));
			strcpy(newnode->ac,tran);
			newnode->ac_array = (double*)malloc(pos*sizeof(double));
			for(int m=0;m<pos;m++){
				newnode->ac_array[m] = tran_array[m];
			}
			newnode->tran = NULL;
		}
		else{
			newnode->tran = (char*)malloc((strlen(tran)+1)*sizeof(char));
			strcpy(newnode->tran,tran);
			if(strcmp(tran, "exp")==0){
				newnode->exp = (double*)malloc(pos*sizeof(double));
				for(int m=0;m<pos;m++){
					newnode->exp[m] = tran_array[m];
				}
				newnode->pos = pos;
			}
			else if(strcmp(tran, "sin")==0){
				newnode->sin = (double*)malloc(pos*sizeof(double));
				for(int m=0;m<pos;m++){
					newnode->sin[m] = tran_array[m];
				}
				newnode->pos = pos;
			}
			else if(strcmp(tran, "pulse")==0){
				newnode->pulse = (double*)malloc(pos*sizeof(double));
				for(int m=0;m<pos;m++){
					newnode->pulse[m] = tran_array[m];
				}
				newnode->pos = pos;
			}
			else if(strcmp(tran, "pwl")==0){
				newnode->pwl = (double*)malloc(pos*sizeof(double));
				for(int m=0;m<pos;m++){
					newnode->pwl[m] = tran_array[m];
				}
				newnode->pos = pos;
			}
		}
		if(ac != NULL){
			if(strcmp(ac, "ac")==0){
				newnode->ac = (char*)malloc((strlen(ac)+1)*sizeof(char));
				strcpy(newnode->ac,ac);
				newnode->ac_array = (double*)malloc(2*sizeof(double));
				for(int m=0;m<2;m++){
					newnode->ac_array[m] = ac_array[m];
				}
			}
		}
		else if(strcmp(tran, "ac")!=0){
			newnode->ac = ac;
		}
	}
	else{
		newnode->ac = ac;
		newnode->tran = tran;
	}

    return newnode;
}

void createlist(struct list** head,char type,char* name,char* node_pos,char* node_neg,double value,char* tran,double* tran_array,int pos, char* ac,double* ac_array,char* model, char* xsv,char* memristor) {


    struct list* newnode = createnode(type,name,node_pos,node_neg,value,tran,tran_array,pos,ac,ac_array,model,xsv,memristor);

	if (*head == NULL) {
        *head = newnode;
		(*head)->last = *head;
    } else {
        struct list* current = (*head)->last;
        current->next = newnode;
		(*head)->last = newnode;
    }
}

void printLinkedList(struct list* head) {
    struct list* current = head;
    while (current != NULL) {
		if (current->type == 'd'){
			printf("Type: %c, Name: %s, Node1: %s, Node2: %s, Model: %s", current->type, current->name, current->node_pos, current->node_neg, current->model.modelname);
			printf("\nModel params: Is=%.5e N=%lf Vt=%lf", current->model.is, current->model.n, current->model.vt);
		}
		else if (current->type == 'x'){
			printf("Type: %c, Name: %s, Node1: %s, Node2: %s, XSV Name: %s, Memristor: %s", current->type, current->name, current->node_pos, current->node_neg, current->params.xsv, current->params.memristorname);
			printf("\nMemristor params: a1=%lf a2=%lf b=%lf Vp=%lf Vn=%lf Ap=%lf An=%lf xp=%lf xn=%lf alphap=%lf alphan=%lf x0=%lf eta=%lf", current->params.a1, current->params.a2, current->params.b, current->params.Vp, current->params.Vn, current->params.Ap, current->params.An, current->params.xp, current->params.xn, current->params.alphap, current->params.alphan, current->params.x0, current->params.eta);
		}else{
			printf("Type: %c, Name: %s, Node1: %s, Node2: %s, Value: %lf", current->type, current->name, current->node_pos, current->node_neg, current->value);
		}
		if(current->tran!=NULL){
			if(strcmp(current->tran,"exp")==0){
				printf(", Tran: %s Values: ",current->tran);
				for(int i=0; i<current->pos; i++){
					printf("%lf ",current->exp[i]);
				}
			}
			else if(strcmp(current->tran,"sin")==0){
				printf(", Tran: %s Values: ",current->tran);
				for(int i=0; i<current->pos; i++){
					printf("%lf ",current->sin[i]);
				}
			}
			else if(strcmp(current->tran,"pulse")==0){
				printf(", Tran: %s Values: ",current->tran);
				for(int i=0; i<current->pos; i++){
					printf("%lf ",current->pulse[i]);
				}
			}
			else if(strcmp(current->tran,"pwl")==0){
				printf(", Tran: %s Values: ",current->tran);
				for(int i=0; i<current->pos; i++){
					printf("%lf ",current->pwl[i]);
				}
			}
		}
		if(current->ac != NULL){
			if(strcmp(current->ac,"ac")==0){
				printf(", Tran: %s Values: ",current->ac);
				for(int i=0; i<2; i++){
					printf("%lf ",current->ac_array[i]);
				}
			}
		}
		printf("\n");
		current = current->next;
    }
}

void clear_list(struct list** head) {
    struct list* current = *head;
    while (current != NULL) {
        struct list* next = current->next;
        free(current->name);
        free(current->node_pos);
        free(current->node_neg);
		if(current->type =='d'){
			free(current->model.modelname);
		}
		if(current->type =='x'){
			free(current->params.xsv);
			free(current->params.memristorname);
		}
		if(current->tran != NULL){
			if(strcmp(current->tran,"exp")==0){
				free(current->exp);
				free(current->tran);
			}
			else if(strcmp(current->tran,"sin")==0){
				free(current->sin);
				free(current->tran);
			}
			else if(strcmp(current->tran,"pulse")==0){
				free(current->pulse);
				free(current->tran);
			}
			else if(strcmp(current->tran,"pwl")==0){
				free(current->pwl);
				free(current->tran);
			}
		}
		if(current->ac != NULL){
			if(strcmp(current->ac,"ac")==0){
				free(current->ac_array);
				free(current->ac);
			}
		}
		free(current);
		current = next;
    }
    *head = NULL;
}

void parser_input(struct list **head,char* filename,struct dc_sweep **head_dc, struct transient **head_tran, struct ac **head_ac, struct circuit *c){

	char curr_ch = '\0';
	char type ='\0';
	char type_dc = '\0';
	char *name = NULL;
	char *node_pos = NULL;
	char *node_neg = NULL;
	char *val = NULL;
	char *dot_buff=NULL;
	char *var = NULL;
	char *tran = NULL;
	char *ac = NULL;
	char *sweep = NULL;
	char* model = NULL;
	char* xsv = NULL, *memristor = NULL;
	double *tran_array = NULL;
	double *ac_array = NULL;
	double value=0, start=0, end=0, incr=0, points=0, start_freq=0, end_freq=0;
	int i,h=0,pos=0, hh=0;
	int dc_en=0, tran_en=0, ac_en=0;
	size_t buffer_size = 2;
	FILE* file;
	c->m2 = 0;

	file = fopen(filename,"r");

	if (file == NULL) {
        perror("Error opening file");
        return ;
    }

	while((curr_ch = tolower(getc(file))) != EOF){

		if(curr_ch == '*'){                           //elegxos gia sxolia
			while((curr_ch = tolower(getc(file))) !='\n'){}
		}
		else if(curr_ch == '\n'){						//elexgos gia newline
		}
		else if(curr_ch == '\r'){
		}
		else if(curr_ch == ' ' || curr_ch == '\t'){  //elegxos gia tabs h spaces sthn arxh
		}
		else if(curr_ch =='.'){						//elegxos gia teleia

			curr_ch = tolower(getc(file));

			buffer_size = 2;
			dot_buff = (char*) malloc(buffer_size);

			for(i=0;curr_ch !=' ' && curr_ch !='\t' && curr_ch !='\n' && curr_ch !=EOF;curr_ch = tolower(getc(file)),i++){
				dot_buff[i] = curr_ch;
				if(i>=buffer_size-1){
					buffer_size*=2;
					dot_buff = (char *)realloc(dot_buff,buffer_size);
				}
			}
			dot_buff[i] = '\0';

			if(strcmp(strtok(dot_buff, " \t\n\r"),"options")==0 && curr_ch!='\n' && curr_ch!=EOF){
				free(dot_buff);
				while(curr_ch != '\n' && curr_ch != EOF && curr_ch != '\r'){
					while(curr_ch == ' ' || curr_ch == '\t'){
						curr_ch = tolower(getc(file));
					}
					buffer_size = 2;
					dot_buff = (char*) malloc(buffer_size);

					for(i=0;curr_ch !='\n' && curr_ch != EOF && curr_ch != ' ' &&  curr_ch !='\t';curr_ch = tolower(getc(file)),i++){
						dot_buff[i] = curr_ch;
						if(i>=buffer_size-1){
							buffer_size*=2;
							dot_buff = (char *)realloc(dot_buff,buffer_size);
						}
					}

					dot_buff[i] = '\0';
					if (strcmp(strtok(dot_buff, " \t\n\r"), "spd") == 0) {
							choleski = 1;
					}
					else if(strcmp(strtok(dot_buff, " \t\n\r"), "iter") == 0) {
							iter = 1;
					}
					else if(strstr(strtok(dot_buff, " \t\n\r"),"itol") != NULL){
							if(strcmp(strtok(dot_buff, " \t\n\r"),"itol")==0){
								itol = 1e-3;
							}
							else{
								buffer_size = 2;
								char* dot_buff1 = (char*) malloc(buffer_size);
								for(int n=5;n<strlen(dot_buff);n++){
									dot_buff1[n-5] = dot_buff[n];
									if(n>=buffer_size-1){
										buffer_size*=2;
										dot_buff1 = (char *)realloc(dot_buff1,buffer_size);
									}
								}
								itol = strtod(dot_buff1, NULL);
								free(dot_buff1);
							}
					}
					else if (strcmp(strtok(dot_buff, " \t\n\r"), "sparse") == 0){
						sparse = 1;
					}
					else if(strcmp(strtok(dot_buff, " \t\n\r"), "method=be") == 0){
						backward_euler = 1;
						trapezoidal = 0;
					}
					while(curr_ch == ' ' || curr_ch == '\t'){
						curr_ch = tolower(getc(file));
					}
					free(dot_buff);
				}
			}
			else if(strcmp(strtok(dot_buff, " \t\n\r"),"params")==0 && curr_ch!='\n' && curr_ch!=EOF){
				free(dot_buff);
				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}
				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\n' && curr_ch != EOF && curr_ch != ' ' &&  curr_ch !='\t';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}

				dot_buff[i] = '\0';
				struct list* curr = *head;

				char* memname;
				memname = (char*)malloc(strlen(dot_buff)+1);
				strcpy(memname,dot_buff);

				long int position = ftell(file);
				if (position == -1L) {
					perror("Error getting file position");
					fclose(file);
					return;
				}

				int memristor_exists=0;

				while(curr != NULL){
					if(curr->type == 'x' && (strcmp(strtok(memname, " \t\n\r"),curr->params.memristorname)==0)){
						fseek(file, position, SEEK_SET);
						curr_ch = tolower(getc(file));
						while(curr_ch == '+'){
							free(dot_buff);
							curr_ch = tolower(getc(file));
							while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '\n'){
								curr_ch = tolower(getc(file));
							}
							buffer_size = 2;
							dot_buff = (char*) malloc(buffer_size);

							for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch !='=';curr_ch = tolower(getc(file)),i++){
								dot_buff[i] = curr_ch;
								if(i>=buffer_size-1){
									buffer_size*=2;
									dot_buff = (char *)realloc(dot_buff,buffer_size);
								}
							}
							dot_buff[i] = '\0';

							if(strcmp(strtok(dot_buff, " \t\n\r"),"a1") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.a1 = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"a2") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.a2 = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"b") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.b = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"vp") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.Vp = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"vn") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.Vn = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"ap") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.Ap = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"an") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.An = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"xp") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.xp = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"xn") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.xn = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"alphap") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.alphap = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"alphan") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.alphan = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"x0") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.x0 = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"eta") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->params.eta = strtod(dot_buff, NULL);
							}
							else{
								free(dot_buff);
							}
							while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '\n' || curr_ch == '\r'){
								curr_ch = tolower(getc(file));
							}
							if(curr_ch=='.'){
								fseek(file,-1,SEEK_CUR);
							}
						}
						memristor_exists++;
					}
					curr = curr->next;
				}
				free(dot_buff);
				free(memname);
				if(memristor_exists==0){
					while(curr_ch != EOF){
						long int pos = ftell(file);
						if(position == -1L){
							perror("Error getting file position");
							fclose(file);
							return;
						}
						int next_ch = tolower(getc(file));
						if(curr_ch == '.' && isalpha(next_ch)){
							fseek(file,pos,SEEK_SET);
							break;
						}
						fseek(file,pos,SEEK_SET);
						curr_ch = tolower(getc(file));
					}

					if(curr_ch == '.'){
						fseek(file,-1,SEEK_CUR);
					}
				}
			}
			else if(strcmp(strtok(dot_buff, " \t\n\r"),"model")==0 && curr_ch!='\n' && curr_ch!=EOF){
				free(dot_buff);
				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}
				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\n' && curr_ch != EOF && curr_ch != ' ' &&  curr_ch !='\t';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}

				dot_buff[i] = '\0';
				struct list* curr = *head;

				char* modname;
				modname = (char*)malloc(strlen(dot_buff)+1);
				strcpy(modname,dot_buff);

				long int position = ftell(file);
				if (position == -1L) {
					perror("Error getting file position");
					fclose(file);
					return;
				}

				int diode_exists=0;

				while(curr != NULL){
					if(curr->type == 'd' && (strcmp(strtok(modname, " \t\n\r"),curr->model.modelname)==0)){
						fseek(file, position, SEEK_SET);
						curr_ch = tolower(getc(file));
						while(curr_ch == '+'){
							free(dot_buff);
							curr_ch = tolower(getc(file));
							while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '\n'){
								curr_ch = tolower(getc(file));
							}
							buffer_size = 2;
							dot_buff = (char*) malloc(buffer_size);

							for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch !='=';curr_ch = tolower(getc(file)),i++){
								dot_buff[i] = curr_ch;
								if(i>=buffer_size-1){
									buffer_size*=2;
									dot_buff = (char *)realloc(dot_buff,buffer_size);
								}
							}
							dot_buff[i] = '\0';
							if(strcmp(strtok(dot_buff, " \t\n\r"),"is") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->model.is = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"n") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->model.n = strtod(dot_buff, NULL);
							}
							else if(strcmp(dot_buff,"vt") == 0){
								free(dot_buff);
								while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '='){
									curr_ch = tolower(getc(file));
								}
								buffer_size = 2;
								dot_buff = (char*) malloc(buffer_size);

								for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch!='\n';curr_ch = tolower(getc(file)),i++){
									dot_buff[i] = curr_ch;
									if(i>=buffer_size-1){
										buffer_size*=2;
										dot_buff = (char *)realloc(dot_buff,buffer_size);
									}
								}
								dot_buff[i] = '\0';
								curr->model.vt = strtod(dot_buff, NULL);
							}
							else{
								free(dot_buff);
							}
							while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '\n' || curr_ch == '\r'){
								curr_ch = tolower(getc(file));
							}
							if(curr_ch=='.'){
								fseek(file,-1,SEEK_CUR);
							}
						}
						diode_exists++;
					}
					curr = curr->next;
				}
				free(dot_buff);
				free(modname);
				if(diode_exists==0){
					while(curr_ch != EOF){
						long int pos = ftell(file);
						if(position == -1L){
							perror("Error getting file position");
							fclose(file);
							return;
						}
						int next_ch = tolower(getc(file));
						if(curr_ch == '.' && isalpha(next_ch)){
							fseek(file,pos,SEEK_SET);
							break;
						}
						fseek(file,pos,SEEK_SET);
						curr_ch = tolower(getc(file));
					}

					if(curr_ch == '.'){
						fseek(file,-1,SEEK_CUR);
					}
				}
			}
			else if(strcmp(strtok(dot_buff, " \t\n\r"),"dc")==0){
				if(var != NULL){
					free(var);
				}
				free(dot_buff);

				while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '\n'){
					curr_ch = tolower(getc(file));
				}
				type_dc = curr_ch;
				curr_ch = tolower(getc(file));
				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\t' && curr_ch !=' ' && curr_ch !='\n';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}
				dot_buff[i] = '\0';

				var = (char*)malloc((strlen(dot_buff)+1)*sizeof(char));
				strcpy(var,dot_buff);

				free(dot_buff);

				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}

				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\n' && curr_ch !=EOF && curr_ch !=' ';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}
				dot_buff[i] = '\0';

				start = strtod(dot_buff, NULL);

				free(dot_buff);

				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}

				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\n' && curr_ch !=EOF && curr_ch !=' ';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}
				dot_buff[i] = '\0';

				end = strtod(dot_buff, NULL);

				free(dot_buff);

				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}

				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\n' && curr_ch !=EOF && curr_ch !=' ';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}
				dot_buff[i] = '\0';

				incr = strtod(dot_buff, NULL);

				free(dot_buff);

				dc_en = 1;

			}
			else if(strcmp(strtok(dot_buff, " \t\n\r"),"tran")==0){
				free(dot_buff);
				tran_flg = 1;
				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}

				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\n' && curr_ch !=EOF && curr_ch !=' ';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}
				dot_buff[i] = '\0';

				char *endptr;
				time_step = strtod(dot_buff, &endptr);
				if (strcmp(endptr, "s") == 0 || *endptr == '\0') {
					time_step = time_step;  // Seconds
				} else if (strcmp(endptr, "ms") == 0) {
					time_step = time_step * 1e-3;  // Convert milliseconds to seconds
				} else if (strcmp(endptr, "us") == 0) {
					time_step = time_step * 1e-6;  // Convert microseconds to seconds
				} else if (strcmp(endptr, "ns") == 0) {
					time_step = time_step * 1e-9;  // Convert nanoseconds to seconds
				} else if (strcmp(endptr, "ps") == 0) {
					time_step = time_step * 1e-12; // Convert picoseconds to seconds
				}

				free(dot_buff);

				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}

				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\n' && curr_ch !=EOF && curr_ch !=' ';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}
				dot_buff[i] = '\0';

				char *endptr1;
				fin_time = strtod(dot_buff, &endptr1);
				if (strcmp(endptr1, "s") == 0 || *endptr1 == '\0') {
					fin_time = fin_time;  // Seconds
				} else if (strcmp(endptr1, "ms") == 0) {
					fin_time = fin_time * 1e-3;  // Convert milliseconds to seconds
				} else if (strcmp(endptr1, "us") == 0) {
					fin_time = fin_time * 1e-6;  // Convert microseconds to seconds
				} else if (strcmp(endptr1, "ns") == 0) {
					fin_time = fin_time * 1e-9;  // Convert nanoseconds to seconds
				} else if (strcmp(endptr1, "ps") == 0) {
					fin_time = fin_time * 1e-12; // Convert picoseconds to seconds
				}

				free(dot_buff);

				tran_en = 1;

			}
			else if(strcmp(strtok(dot_buff, " \t\n\r"),"ac")==0){
				if(sweep != NULL){
					free(var);
				}
				free(dot_buff);
				ac_flg = 1;
				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}

				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\n' && curr_ch !=EOF && curr_ch !=' ';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}
				dot_buff[i] = '\0';

				sweep = (char*)malloc((strlen(dot_buff)+1)*sizeof(char));
				strcpy(sweep,dot_buff);

				free(dot_buff);

				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}

				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\n' && curr_ch !=EOF && curr_ch !=' ';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}
				dot_buff[i] = '\0';

				points = strtod(dot_buff, NULL);

				free(dot_buff);

				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}

				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\n' && curr_ch !=EOF && curr_ch !=' ';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}
				dot_buff[i] = '\0';

				start_freq = strtod(dot_buff, NULL);

				free(dot_buff);

				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}

				buffer_size = 2;
				dot_buff = (char*) malloc(buffer_size);

				for(i=0;curr_ch !='\n' && curr_ch !=EOF && curr_ch !=' ';curr_ch = tolower(getc(file)),i++){
					dot_buff[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size*=2;
						dot_buff = (char *)realloc(dot_buff,buffer_size);
					}
				}
				dot_buff[i] = '\0';

				end_freq = strtod(dot_buff, NULL);

				free(dot_buff);

				ac_en = 1;

			}else if(strcmp(strtok(dot_buff, " \t\n\r"),"plot")==0 || strcmp(strtok(dot_buff, " \t\n\r"),"print")==0){
				free(dot_buff);
				while(curr_ch != '\n' && curr_ch != EOF){
					while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '\n'){
						curr_ch = tolower(getc(file));
					}
					curr_ch = tolower(getc(file));
					curr_ch = tolower(getc(file));
					buffer_size = 2;
					dot_buff = (char*) malloc(buffer_size);

					for(i=0;curr_ch !=')';curr_ch = tolower(getc(file)),i++){
						dot_buff[i] = curr_ch;
						if(i>=buffer_size-1){
							buffer_size*=2;
							dot_buff = (char *)realloc(dot_buff,buffer_size);
						}
					}
					dot_buff[i] = '\0';
					if(dc_en == 1){
						create_dc_sweep(head_dc,type_dc,var,start,end,incr,dot_buff);
					}
					if(tran_en == 1){
						create_transient(head_tran,time_step,fin_time,dot_buff);
					}
					if(ac_en == 1){
						create_ac(head_ac,sweep,points,start_freq,end_freq,dot_buff);
					}
					free(dot_buff);
					while(curr_ch == ')' || curr_ch == ' ' || curr_ch == '\t' || curr_ch == '\r'){
						curr_ch = tolower(getc(file));
					}
				}
			}
			else{
				free(dot_buff);
			}
		}
		else{											//eisagwgh twn stoixeiwn sthn lista
			h=0;
			hh=0;
			while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == '\n'){
				curr_ch = tolower(getc(file));
			}
			type = curr_ch;
			if(type =='l' || type =='v'){
				c->m2++;
			}
			curr_ch = tolower(getc(file));

			buffer_size = 2;
			name = (char*)malloc(buffer_size);
			name[0] = curr_ch;
			curr_ch = tolower(getc(file));

			for(i=1; curr_ch != ' ' && curr_ch != '\t'; curr_ch = tolower(getc(file)), i++){
				name[i] = curr_ch;
				if(i>=buffer_size-1){
					buffer_size *= 2;
					name = (char *)realloc(name,buffer_size);
				}
			}
			name[i] = '\0';

			while(curr_ch == ' ' || curr_ch == '\t'){
				curr_ch = tolower(getc(file));
			}

			buffer_size = 2;
			node_pos = (char*)malloc(buffer_size);
			node_pos[0] = curr_ch;
			curr_ch = tolower(getc(file));

			for(i=1; curr_ch != ' ' && curr_ch != '\t'; curr_ch = tolower(getc(file)), i++){
				node_pos[i] = curr_ch;
				if(i>=buffer_size-1){
					buffer_size *= 2;
					node_pos = (char *)realloc(node_pos,buffer_size);
				}
			}
			node_pos[i] = '\0';

			while(curr_ch == ' ' || curr_ch == '\t'){
				curr_ch = tolower(getc(file));
			}

			buffer_size = 2;
			node_neg = (char*)malloc(buffer_size);
			node_neg[0] = curr_ch;
			curr_ch = tolower(getc(file));

			for(i=1; curr_ch != ' ' && curr_ch != '\t'; curr_ch = tolower(getc(file)), i++){
				node_neg[i] = curr_ch;
				if(i>=buffer_size-1){
					buffer_size *= 2;
					node_neg = (char *)realloc(node_neg,buffer_size);
				}
			}
			node_neg[i] = '\0';

			while(curr_ch == ' ' || curr_ch == '\t'){
				curr_ch = tolower(getc(file));
			}

			if(type == 'd'){
				diode++;
				buffer_size = 2;
				model = (char*)malloc(buffer_size);
				model[0] = curr_ch;
				curr_ch = tolower(getc(file));

				for(i=1; curr_ch != ' ' && curr_ch != '\t' && curr_ch != '\n' && curr_ch != EOF; curr_ch = tolower(getc(file)), i++){
					model[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size *= 2;
						model = (char *)realloc(model,buffer_size);
					}
				}
				model[i] = '\0';


			}else if(type == 'x'){
				memristor_c++;
				buffer_size = 2;
				xsv = (char*)malloc(buffer_size);
				xsv[0] = curr_ch;
				curr_ch = tolower(getc(file));

				for(i=1; curr_ch != ' ' && curr_ch != '\t' && curr_ch != '\n' && curr_ch != EOF; curr_ch = tolower(getc(file)), i++){
					xsv[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size *= 2;
						xsv = (char *)realloc(xsv,buffer_size);
					}
				}
				xsv[i] = '\0';

				while(curr_ch == ' ' || curr_ch == '\t'){
					curr_ch = tolower(getc(file));
				}

				buffer_size = 2;
				memristor = (char*)malloc(buffer_size);
				memristor[0] = curr_ch;
				curr_ch = tolower(getc(file));

				for(i=1; curr_ch != ' ' && curr_ch != '\t' && curr_ch != '\n' && curr_ch != EOF; curr_ch = tolower(getc(file)), i++){
					memristor[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size *= 2;
						memristor = (char *)realloc(memristor,buffer_size);
					}
				}
				memristor[i] = '\0';
			}
			else if(isdigit(curr_ch)){
				buffer_size = 2;
				val = (char*)malloc(buffer_size);
				val[0] = curr_ch;
				curr_ch = tolower(getc(file));

				for(i=1; curr_ch != ' ' && curr_ch != '\t' && curr_ch != '\n' && curr_ch != EOF; curr_ch = tolower(getc(file)), i++){
					val[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size *= 2;
						val = (char *)realloc(val,buffer_size);
					}
				}
				val[i] = '\0';
				value = strtod(val, NULL);
			}

			while(curr_ch != '\n' && (curr_ch == ' ' || curr_ch == '\t')){
				curr_ch = tolower(getc(file));
			}

			if(curr_ch != '\n'){
				h=1;
				buffer_size = 2;
				tran = (char*)malloc(buffer_size);
				tran[0] = curr_ch;
				curr_ch = tolower(getc(file));

				for(i=1; curr_ch != ' ' && curr_ch != '\t' && curr_ch != '('; curr_ch = tolower(getc(file)), i++){
					tran[i] = curr_ch;
					if(i>=buffer_size-1){
						buffer_size *= 2;
						tran = (char *)realloc(tran,buffer_size);
					}
				}
				tran[i] = '\0';

				if(strcmp(strtok(tran, " \t\n\r"),"exp")==0 || strcmp(strtok(tran, " \t\n\r"),"sin")==0){
					pos=0;
					char* temp = NULL;
					tran_array = (double *) malloc(6*sizeof(double*));
					while(curr_ch == ' ' || curr_ch == '\t'){
						curr_ch = tolower(getc(file));
					}
					curr_ch = tolower(getc(file));
					while(curr_ch != ')'){
						while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == ','){
							curr_ch = tolower(getc(file));
						}
						buffer_size = 2;
						temp = (char*)malloc(buffer_size);
						temp[0] = curr_ch;
						curr_ch = tolower(getc(file));

						for(i=1; curr_ch != ' ' && curr_ch != '\t' && curr_ch != ')' && curr_ch != ','; curr_ch = tolower(getc(file)), i++){
							temp[i] = curr_ch;
							if(i>=buffer_size-1){
								buffer_size *= 2;
								temp = (char *)realloc(temp,buffer_size);
							}
						}
						temp[i] = '\0';
						tran_array[pos] = strtod(temp, NULL);
						free(temp);
						pos++;
					}
				}
				else if(strcmp(strtok(tran, " \t\n\r"),"pulse")==0){
					pos=0;
					char* temp = NULL;
					tran_array = (double *) malloc(7*sizeof(double));
					while(curr_ch == ' ' || curr_ch == '\t'){
						curr_ch = tolower(getc(file));
					}
					curr_ch = tolower(getc(file));
					while(curr_ch != ')'){
						while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == ','){
							curr_ch = tolower(getc(file));
						}
						buffer_size = 2;
						temp = (char*)malloc(buffer_size);
						temp[0] = curr_ch;
						curr_ch = tolower(getc(file));

						for(i=1; curr_ch != ' ' && curr_ch != '\t' && curr_ch != ')' && curr_ch != ','; curr_ch = tolower(getc(file)), i++){
							temp[i] = curr_ch;
							if(i>=buffer_size-1){
								buffer_size *= 2;
								temp = (char *)realloc(temp,buffer_size);
							}
						}
						temp[i] = '\0';
						tran_array[pos] = strtod(temp, NULL);
						free(temp);
						pos++;
					}
				}
				else if(strcmp(strtok(tran, " \t\n\r"),"pwl")==0){
					pos=0;
					int ke=0;
					char* temp = NULL;
					tran_array = (double *) malloc(2*sizeof(double));
					while(curr_ch != '\n'){
						while(curr_ch == ' ' || curr_ch == '\t'){
							curr_ch = tolower(getc(file));
						}
						curr_ch = tolower(getc(file));
						while(curr_ch != ')'){
							while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == ','){
								curr_ch = tolower(getc(file));
							}
							buffer_size = 2;
							temp = (char*)malloc(buffer_size);
							temp[0] = curr_ch;
							curr_ch = tolower(getc(file));

							for(i=1; curr_ch != ' ' && curr_ch != '\t' && curr_ch != ')' && curr_ch != ','; curr_ch = tolower(getc(file)), i++){
								temp[i] = curr_ch;
								if(i>=buffer_size-1){
									buffer_size *= 2;
									temp = (char *)realloc(temp,buffer_size);
								}
							}
							temp[i] = '\0';

							char *endptr;
							tran_array[pos] = strtod(temp, &endptr);
							if (strcmp(endptr, "s") == 0 || *endptr == '\0') {
								tran_array[pos] = tran_array[pos];  // Seconds
							} else if (strcmp(endptr, "ms") == 0) {
								tran_array[pos] = tran_array[pos] * 1e-3;  // Convert milliseconds to seconds
							} else if (strcmp(endptr, "us") == 0) {
								tran_array[pos] = tran_array[pos] * 1e-6;  // Convert microseconds to seconds
							} else if (strcmp(endptr, "ns") == 0) {
								tran_array[pos] = tran_array[pos] * 1e-9;  // Convert nanoseconds to seconds
							} else if (strcmp(endptr, "ps") == 0) {
								tran_array[pos] = tran_array[pos] * 1e-12; // Convert picoseconds to seconds
							}

							free(temp);
							pos++;
							ke=1;
						}
						curr_ch = tolower(getc(file));
						if(ke==1){
							tran_array = (double *) realloc(tran_array,(pos+2)*sizeof(double));
						}
					}
				}
				else if(strcmp(strtok(tran, " \t\n\r"),"ac")==0){
					pos=0;
					char* temp = NULL;
					tran_array = (double *) malloc(2*sizeof(double*));
					while(curr_ch == ' ' || curr_ch == '\t'){
						curr_ch = tolower(getc(file));
					}
					while(pos != 2){
						while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == ','){
							curr_ch = tolower(getc(file));
						}
						buffer_size = 2;
						temp = (char*)malloc(buffer_size);
						temp[0] = curr_ch;
						curr_ch = tolower(getc(file));

						for(i=1; curr_ch != ' ' && curr_ch != '\t' && curr_ch != ',' && curr_ch != '\n'; curr_ch = tolower(getc(file)), i++){
							temp[i] = curr_ch;
							if(i>=buffer_size-1){
								buffer_size *= 2;
								temp = (char *)realloc(temp,buffer_size);
							}
						}
						temp[i] = '\0';
						tran_array[pos] = strtod(temp, NULL);
						free(temp);
						pos++;
					}
				}


				while(curr_ch != '\n' && (curr_ch == ' ' || curr_ch == '\t' || curr_ch == ')')){
					curr_ch = tolower(getc(file));
				}
				if(curr_ch == 'a'){
					hh=1;
					buffer_size = 2;
					ac = (char*)malloc(buffer_size);
					ac[0] = curr_ch;
					curr_ch = tolower(getc(file));

					for(i=1; curr_ch != ' ' && curr_ch != '\t' && curr_ch != '('; curr_ch = tolower(getc(file)), i++){
						ac[i] = curr_ch;
						if(i>=buffer_size-1){
							buffer_size *= 2;
							ac = (char *)realloc(ac,buffer_size);
						}
					}
					ac[i] = '\0';

					if(strcmp(strtok(ac, " \t\n\r"),"ac")==0){
						int pospos=0;
						char* temp = NULL;
						ac_array = (double *) malloc(2*sizeof(double*));
						while(curr_ch == ' ' || curr_ch == '\t'){
							curr_ch = tolower(getc(file));
						}
						while(pospos != 2){
							while(curr_ch == ' ' || curr_ch == '\t' || curr_ch == ','){
								curr_ch = tolower(getc(file));
							}
							buffer_size = 2;
							temp = (char*)malloc(buffer_size);
							temp[0] = curr_ch;
							curr_ch = tolower(getc(file));

							for(i=1; curr_ch != ' ' && curr_ch != '\t' && curr_ch != ',' && curr_ch != '\n'; curr_ch = tolower(getc(file)), i++){
								temp[i] = curr_ch;
								if(i>=buffer_size-1){
									buffer_size *= 2;
									temp = (char *)realloc(temp,buffer_size);
								}
							}
							temp[i] = '\0';
							ac_array[pospos] = strtod(temp, NULL);
							free(temp);
							pospos++;
						}
					}
				}
				else{
					ac = NULL;
					ac_array = NULL;
				}

				while(curr_ch != '\n' && (curr_ch == ' ' || curr_ch == '\t')){
					curr_ch = tolower(getc(file));
				}
			}
			else{
				tran = NULL;
				tran_array = NULL;
				ac = NULL;
				ac_array = NULL;
			}

			if(curr_ch =='*'){
				while((curr_ch = tolower(getc(file))) !='\n'){}
			}

			createlist(head,type,name,node_pos,node_neg,value,tran,tran_array,pos,ac,ac_array,model,xsv,memristor);
			free(name);
			free(node_pos);
			free(node_neg);
			if(type !='d' && type != 'x'){
				free(val);
			}
			else if(type =='d'){
				free(model);
			}
			else if(type =='x'){
				free(memristor);
				free(xsv);
			}
			if(h==1){
				free(tran);
				free(tran_array);
			}
			if(hh==1){
				free(ac);
				free(ac_array);
			}
		}
	}
	if(var != NULL){
		free(var);
	}
	if(sweep != NULL){
		free(sweep);
	}
	fclose(file);
	printLinkedList(*head);
	return;
}
