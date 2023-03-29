# this makefile is intended for g++ on Linux

CC = g++
CFLAGS = -c -Wall -O3
# CFLAGS = -c
LDFLAGS =
INCLUDE = -I../include/
SOURCES =  filter.cpp iofunc.cpp PLL.cpp project.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = project

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	-rm $(OBJECTS)
