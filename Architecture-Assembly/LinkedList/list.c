#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "list.h"

void panic(char *s) {
    fprintf(stderr, "%s\n", s);
    exit(1);
}

struct list *new() {
    struct list *newlist = (struct list *)malloc(sizeof(struct list));
    if(newlist == NULL){
      panic("Out of memory! Cannot create new list.");
    }
    newlist->size = 0;
    newlist->head = NULL;
    newlist->tail = NULL;
    //deal with assertion
    return newlist;
}

void destroy(struct list *l) {
    assert(l != NULL);
    
    while(l->size != 0){
      remove_at(l, 0);
    }
    free(l);
}

int size(struct list *l) {
    assert(l != NULL);
    return l->size;
}

void add_tail(struct list *l, long val) {
    assert(l != NULL);
    //deal with assertion
    struct node *newNode = (struct node *)malloc(sizeof(struct node));
    if(newNode == NULL){
      panic("Out of memory! Cannot add tail.");
    }
    //deal with assertion
    newNode->val = val;
    newNode->next = NULL;
    if(l->tail != NULL){
      l->tail->next = newNode;
      newNode->prev = l->tail;
      l->tail = newNode;
    }
    else{
      newNode->prev = NULL;
      l->tail = newNode;
      l->head = newNode;
    }
    l->size++;
}

long remove_at(struct list *l, int index) {
    assert(l != NULL && index >= 0 && index < l->size);

    struct node *currentNode = l->head;
    for(int i = 1; i <= index; i++){
      currentNode = currentNode->next;
    }

    if(currentNode->prev == NULL){
      l->head = currentNode->next;
      if(l->head != NULL)
	l->head->prev = NULL;
    }
    else if(currentNode->next == NULL){
      l->tail = currentNode->prev;
      if(l->tail != NULL)
	l->tail->next = NULL;
    }
    if(currentNode->prev != NULL && currentNode->next != NULL){
      currentNode->prev->next = currentNode->next;
      currentNode->next->prev = currentNode->prev;
    }
    
    long value = currentNode->val;
    free(currentNode);
    l->size--;
    return value;
}

long get(struct list *l, int index) {
    assert(l != NULL && index >= 0 && index < l->size);
    struct node *currentNode = l->head;
    for(int i = 1; i <= index; i++){
      currentNode = currentNode->next;
    }
    return currentNode->val;
}
