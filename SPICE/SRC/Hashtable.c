#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include"spice.h"

unsigned int hash_function(const char* key,int size) {
    unsigned int hash = 0;
    while (*key) {
        hash = (hash * 31) + (*key++);
    }
    return hash % size;
}

void create_hashtable(struct list **head,struct HashTable *t) {
    int found = 0;
    long int size, id = 1;
    struct list *current = *head;


    size = 0;
    while(current != NULL){
        size++;
        current = current->next;
    }
    t->hash_size = size/2;
    t->table = (Node**)malloc(t->hash_size * sizeof(Node*));
    for (size_t i = 0; i < t->hash_size; i++) {
        t->table[i] = NULL;
    }

    struct list *cur = *head;

    while(cur != NULL){
        if(strcmp(cur->node_pos,"0")!=0){
            found = 0;
            unsigned int index = hash_function(cur->node_pos,t->hash_size);

            Node* current1 = t->table[index];
            while (current1 != NULL) {
                if (strcmp(current1->key, cur->node_pos) == 0) {
                    found = 1;
                }
                current1 = current1->next;
            }

            if(found != 1){
                Node* newNode = malloc(sizeof(Node));
                newNode->key = (char*)malloc((strlen(cur->node_pos)+1)*sizeof(char));
                strcpy(newNode->key, cur->node_pos);
                newNode->id = id;
                id++;
                newNode->next = t->table[index];
                t->table[index] = newNode;
            }
        }
        if(strcmp(cur->node_neg,"0")!=0){
            found = 0;
            unsigned int index = hash_function(cur->node_neg,t->hash_size);

            Node* current2 = t->table[index];
            while (current2 != NULL) {
                if (strcmp(current2->key, cur->node_neg) == 0) {
                    found = 1;
                }
                current2 = current2->next;
            }

            if(found != 1){
                Node* newNode = malloc(sizeof(Node));
                newNode->key = (char*)malloc((strlen(cur->node_neg)+1)*sizeof(char));
				strcpy(newNode->key, cur->node_neg);
                newNode->id = id;
                id++;
                newNode->next = t->table[index];
                t->table[index] = newNode;
            }
        }
        /*if(cur->type == 'x'){
            found = 0;
            unsigned int index = hash_function(cur->params.xsv,t->hash_size);

            Node* current1 = t->table[index];
            while (current1 != NULL) {
                if (strcmp(current1->key, cur->params.xsv) == 0) {
                    found = 1;
                }
                current1 = current1->next;
            }

            if(found != 1){
                Node* newNode = malloc(sizeof(Node));
                newNode->key = (char*)malloc((strlen(cur->params.xsv)+1)*sizeof(char));
                strcpy(newNode->key, cur->params.xsv);
                newNode->id = id;
                id++;
                newNode->next = t->table[index];
                t->table[index] = newNode;
            }
        }*/
        cur = cur->next;
    }

    t->unique_ids = id-1;

	/*printf("Hash Table Contents:\n");
    for (int i = 0; i < t->hash_size; i++) {
        printf("[%d]: ",i);
        Node* currentt = t->table[i];
        while (currentt != NULL) {
            printf("(Key: %s, Value: %ld) -> ", currentt->key, currentt->id);
            currentt = currentt->next;
        }
        printf("\n");
    }*/
}

int hash_find(struct HashTable *t,const char* key){
	int index = hash_function(key,t->hash_size);

	Node* current = t->table[index];
	while (current != NULL) {
        if (strcmp(current->key, key) == 0) {
            return current->id;
        }
        current = current->next;
    }

	return -1;
}

char* find_string(struct HashTable* t,int id) {
    for (int i = 0; i < t->hash_size; i++) {
        Node* current = t->table[i];
        while (current != NULL) {
            if (current->id == id) {
                return current->key;
            }
            current = current->next;
        }
    }
    return NULL;
}

void free_hashtable(struct HashTable* t) {

	for (int i = 0; i < t->hash_size; i++) {
		Node* currentt = t->table[i];
		while (currentt != NULL) {
			Node* temp = currentt;
			currentt = currentt->next;
			free(temp->key);
			free(temp);
		}
	}
	free(t->table);
} 
