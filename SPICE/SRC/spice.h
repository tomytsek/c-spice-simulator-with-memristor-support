#ifndef SPICE_H
#define SPICE_H

#include "/usr/include/suitesparse/cs.h"
#include <math.h>
#include <complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

extern int LU;
extern int choleski;
extern int iter;
extern double itol;
extern int sparse;
extern int trapezoidal;
extern int backward_euler;
extern int tran_flg;
extern double time_step;
extern double fin_time;
extern int in_tran;
extern int diode;
extern int memristor_c;

struct dc_sweep{
	char type_dc;
	char* var;
	double start;
	double end;
	double incr;
	char *node_p;
	struct dc_sweep* next;
};

struct transient{
	double time_step;
	double fin_time;
	char *node_p;
	struct transient* next;
	FILE * file;
};

struct ac{
	char *sweep;
	double points;
	double start_freq;
	double end_freq;
	char *node_p;
	struct ac* next;
};

struct model{
	char* modelname;
	double is;
	double n;
	double vt;
};

struct params{
	char* memristorname;
	char* xsv;
	double a1;
	double a2;
	double b;
	double Vp;
	double Vn;
	double Ap;
	double An;
	double xp;
	double xn;
	double alphap;
	double alphan;
	double x0;
	double eta;
};

struct list{
	char type;
	char* name;
	char* node_pos;
	char* node_neg;
	struct model model;
	struct params params;
	double value;
	char* tran;
	char* ac;
	double* exp;
	double* sin;
	double* pulse;
	double* pwl;
	double* ac_array;
	int pos;
	struct list* next;
	struct list* last;

};

struct circuit{
	double** A;
	double* b;
	double* x;
	int m2;
	double** l;
	double** u;
	int* p;
	long int nz;
	cs* A_sparse;
	double complex** C_A;
	double complex* C_b;
	double complex* C_x;
	double complex** C_l;
	double complex** C_u;
	double** h;
	double* gx;
	double** J;
};

struct tran_circuit{
	double** A;
	double** old_A;
	double* b;
	double* x;
	double** C;
	double** l;
	double** u;
	int* p;
	cs* A_sparse;
	cs* C_sparse;
};

typedef struct Node {
    char* key;
    long int id;
    struct Node* next;
} Node;

struct HashTable {
    Node **table;
    long int hash_size;
    long int unique_ids;
};

//PARSER

struct list* createnode(char type,char* name,char* node_pos,char* node_neg,double value,char* tran,double* tran_array, int pos, char* ac,double* ac_array,char* model,char* xsv,char* memristor);

void createlist(struct list** head,char type,char* name,char* node_pos,char* node_neg,double value,char* tran,double* tran_array, int pos, char* ac,double* ac_array,char* model,char* xsv,char* memristor);

void printLinkedList(struct list* head);

void clear_list(struct list** head);

void parser_input(struct list **head,char* filename,struct dc_sweep **head_dc, struct transient **head_tran, struct ac **head_ac, struct circuit *c);


// DC_SWEEP

struct dc_sweep* createNode_dc_sweep(char type_dc, char* var, double start, double end, double incr, char* node_p);

void create_dc_sweep(struct dc_sweep** head, char type_dc,char* var, double start, double end, double incr, char* node_p);

void clear_dc_sweep(struct dc_sweep** head);

int printDc_sweep(struct dc_sweep* dc,struct list** head,struct HashTable *t,struct circuit* c);


// HASH_TABLE

unsigned int hash_function(const char* key,int size);

void create_hashtable(struct list **head,struct HashTable *t);

int hash_find(struct HashTable *t,const char* key);

char* find_string(struct HashTable* t,int id);

void free_hashtable(struct HashTable* t);


// CIRCUIT

void create_tables(struct list** head,struct HashTable* t,struct circuit* c);

void LU_func(struct circuit* c, struct HashTable* t);

int Choleski(struct circuit* c, struct HashTable* t);

void back_for_sub(struct circuit* c, struct HashTable* t);

void conjugate_gradients(struct circuit* c, struct HashTable* t);

void bi_conjugate_gradients(struct circuit* c, struct HashTable* t);

void print_dc_op(struct list ** head,struct HashTable* t,struct circuit* c,char* filename);

void free_circ(struct HashTable* t,struct circuit* c);

void printMatrix(double** matrix,int size);

// SPARSE

void print_sparse_matrix(cs* A);

void create_sparse_t(struct list **head,struct circuit *c, struct HashTable *t);

void sparse_lu_solve(struct circuit *c, struct HashTable *t);

int sparse_chol_solve(struct circuit *c, struct HashTable *t);

void sparse_conjugate_gradients(struct circuit* c, struct HashTable* t);

void sparse_bi_conjugate_gradients(struct circuit* c, struct HashTable* t);

// TRANSIENT

struct transient* createNode_transient(double time_step, double fin_time, char* node_p);

void create_transient(struct transient** head, double time_step, double fin_time, char* node_p);

void tran_LU_func(struct circuit* c, struct tran_circuit* tran, struct HashTable* t);

void tran_back_for_sub(struct circuit* c, struct tran_circuit* tran, struct HashTable* t);

void clear_transient(struct transient** head);

void tran_create_tables(struct list** head,struct HashTable* t,struct tran_circuit* tran,struct circuit* c);

void tran_analysis(struct transient** head1,struct list** head,struct circuit* c,struct tran_circuit* tran,struct HashTable* t);

void tran_sparse_tables(struct list **head,struct circuit *c,struct tran_circuit* tran,struct HashTable *t);

#endif

// AC

struct ac* createNode_ac(char* sweep, double points, double start_freq, double end_freq, char* node_p);

void create_ac(struct ac** head, char* sweep, double points, double start_freq, double end_freq, char* node_p);

void AC_analysis(struct list** head,struct HashTable* t,struct circuit* c,struct ac *head_ac);

void clear_ac(struct ac** head);


//DIODE
void newtons_method(struct circuit *c, struct HashTable *t, struct list **head);

//MEMRISTOR
void memristor_op(struct circuit *c, struct HashTable *t, struct list **head);

double F(double V, double x,struct list* curr);

double G(double V, struct list* curr);


