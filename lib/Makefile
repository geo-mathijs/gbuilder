SHELL  = /bin/sh
CC     = gcc
FLAGS  = -std=gnu99 -Iinclude
CFLAGS = -fPIC -shared -Wall

TARGET  = libpredicates.so
SOURCES = predicates.c
HEADERS = predicates.h
OBJECTS = $(SOURCES:.c=.o)

PREFIX = /usr/local

all: $(TARGET)

install:
	cp ./$(TARGET) $(PREFIX)/lib
	cp ./$(HEADERS) $(PREFIX)/include

$(TARGET): $(OBJECTS)
	$(CC) $(FLAGS) $(CFLAGS) -o $(TARGET) $(OBJECTS)
