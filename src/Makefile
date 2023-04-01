# this makefile is intended for g++ on Linux

CC = g++
CFLAGS = -c -Wall -O3 -mcpu=cortex-a72+crypto -pthread
# CFLAGS = -c
LDFLAGS = -pthread
INCLUDE = -I../include/
SOURCES =  filter.cpp iofunc.cpp PLL.cpp Final_file.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = project

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	-rm $(OBJECTS)
