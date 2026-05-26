#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spice.h"
#include <sys/stat.h>
#include <time.h>

int in_tran=0;

int main(int argc, char* argv[]){
	char directoryName[] = "Results";
	int f=0;
	
	if(argc !=2){
		printf("Wrong number of arguements!\n");
		return 1;
	}

	struct list *head = NULL;
	struct dc_sweep *head_dc = NULL;
	struct transient *head_tran = NULL;
	struct ac *head_ac = NULL;
	struct tran_circuit head_tran_circuit;
    struct HashTable hash_table;
	struct circuit c;

	parser_input(&head,argv[1],&head_dc,&head_tran,&head_ac,&c);
	create_hashtable(&head,&hash_table);
	
	struct stat st;
    if (stat(directoryName, &st) != 0) {
        if (mkdir(directoryName, 0777) != 0) {
            perror("Error creating directory");
            exit(EXIT_FAILURE);
        }
    }

	create_tables(&head,&hash_table,&c);
		
	if(sparse == 1){
		create_sparse_t(&head,&c,&hash_table);
		printf("____SPARSE____\n");
	}

	if(choleski == 1 && iter == 1){

		printf("SPD ITER\n");
		if(sparse==1){
			
			if(sparse_chol_solve(&c,&hash_table)==-1){
				printf("Matrix A is not spd!\n");
				free_circ(&hash_table,&c);
				free_hashtable(&hash_table);
				clear_dc_sweep(&head_dc);
				clear_list(&head);
				return -1;
			}
			create_sparse_t(&head,&c,&hash_table);
			sparse_conjugate_gradients(&c,&hash_table);
			print_dc_op(&head,&hash_table,&c,argv[1]);
		}
		else{
			if(Choleski(&c,&hash_table)==-1){
				printf("Matrix A is not spd!\n");
				free_circ(&hash_table,&c);
				free_hashtable(&hash_table);
				clear_dc_sweep(&head_dc);
				clear_list(&head);
				cs_spfree(c.A_sparse);
				return -1;
			}
			
			conjugate_gradients(&c,&hash_table);
			print_dc_op(&head,&hash_table,&c,argv[1]);
		}

		if(head_dc !=NULL){
			f=1;
			printDc_sweep(head_dc,&head,&hash_table,&c);
		}
	}
	else if(iter == 1){
		printf("ITER\n");

		if(sparse==1){
			sparse_bi_conjugate_gradients(&c,&hash_table);
			print_dc_op(&head,&hash_table,&c,argv[1]);
		}
		else{
			bi_conjugate_gradients(&c,&hash_table);
			print_dc_op(&head,&hash_table,&c,argv[1]);
		}

		if(head_dc !=NULL){
			f=1;
			printDc_sweep(head_dc,&head,&hash_table,&c);
		}
	}
	else if(choleski == 1){

		printf("SPD\n");
		if(sparse==1){
			if(sparse_chol_solve(&c,&hash_table)==-1){
				printf("Matrix A is not spd!\n");
				free_circ(&hash_table,&c);
				free_hashtable(&hash_table);
				clear_dc_sweep(&head_dc);
				clear_list(&head);
				return -1;
			}
			print_dc_op(&head,&hash_table,&c,argv[1]);
		}else{
			
			if(Choleski(&c,&hash_table)==-1){
				printf("Matrix A is not spd!\n");
				free_circ(&hash_table,&c);
				free_hashtable(&hash_table);
				clear_dc_sweep(&head_dc);
				clear_list(&head);
				return -1;
			}
			back_for_sub(&c,&hash_table);
			print_dc_op(&head,&hash_table,&c,argv[1]);
		}


		if(head_dc !=NULL){
			f=1;
			printDc_sweep(head_dc,&head,&hash_table,&c);
		}
	}
	else{
		printf("NUN\n");
		if(sparse==1){
			sparse_lu_solve(&c,&hash_table);
			print_dc_op(&head,&hash_table,&c,argv[1]);
	
		}
		else{
			if(diode!=0){
				newtons_method(&c ,&hash_table,&head);
			}
			/*else if(memristor_c!=0){
				memristor_op(&c,&hash_table,&head);
			}*/
			else if(memristor_c!=1){
				LU_func(&c,&hash_table);
				back_for_sub(&c,&hash_table);
			}
			print_dc_op(&head,&hash_table,&c,argv[1]);
		}

		if(head_dc != NULL){
			f=1;
			printDc_sweep(head_dc,&head,&hash_table,&c);
		}
	}
	
	if(tran_flg){
		in_tran =1;
		if(sparse==1){
			tran_sparse_tables(&head,&c,&head_tran_circuit,&hash_table);
		}
		else{
			tran_create_tables(&head,&hash_table,&head_tran_circuit,&c);
		}
		if(f!=1){
			free_circ(&hash_table,&c);
		}
		tran_analysis(&head_tran,&head,&c,&head_tran_circuit,&hash_table);
	}
	else{
		if(f!=1){
			free_circ(&hash_table,&c);
		}
	}
	
	if(head_ac!=NULL){
		AC_analysis(&head,&hash_table,&c,head_ac);
	}
	
	free_hashtable(&hash_table);
	clear_dc_sweep(&head_dc);
	clear_transient(&head_tran);
	clear_ac(&head_ac);
	clear_list(&head);

	return 0;
}


			/*printf("U:\n");
			printMatrix(c.u,c.m2+hash_table.unique_ids);
			printf("\n");
			printf("L:\n");
			printMatrix(c.l,c.m2+hash_table.unique_ids);
			printf("\n");
			printf("A:\n");
			printMatrix(c.A,c.m2+hash_table.unique_ids);
			printf("\n[");
			for (int i = 0; i < c.m2+hash_table.unique_ids; i++) {
				printf("%f\t",c.x[i]);
			}
			printf("]\n");
			printf("\n[");
			for (int i = 0; i < c.m2+hash_table.unique_ids; i++) {
				printf("%d\t",c.p[i]);
			}
			printf("]\n");*/

