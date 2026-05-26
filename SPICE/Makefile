CC = gcc
CFLAGS = -Wall -g -O0
LDFLAGS = -lgsl -lgslcblas -lm -lcxsparse
SRCDIR = SRC
OBJDIR = OBJ

SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJECTS = $(patsubst $(SRCDIR)/%.c, $(SRCDIR)/%.o, $(SOURCES))
EXECUTABLE = spice

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(SRCDIR)/%.o: $(SRCDIR)/%.c | (SRCDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf $(EXECUTABLE) $(OBJECTS)
