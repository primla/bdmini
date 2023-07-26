cc = g++
CFLAGS = -Wall -Wextra -O3 -std=c++20 -pedantic 

SRCDIR = ./src/
BINDIR = ./bin/

TASKS = pre_process spac fkanaly dspac
all: $(TASKS)

pre_process:
	$(cc) $(CFLAGS) -o $(BINDIR)$@ $(SRCDIR)util.cpp $(SRCDIR)obs_site.cpp $(SRCDIR)io.cpp $(SRCDIR)obs_site_group.cpp $(SRCDIR)pre_process.cpp

spac:
	$(cc) $(CFLAGS) -o $(BINDIR)$@ $(SRCDIR)$@.cpp

fkanaly:
	$(cc) $(CFLAGS) -o $(BINDIR)$@ $(SRCDIR)$@.cpp

dspac:
	$(cc) $(CFLAGS) -o $(BINDIR)$@ $(SRCDIR)pso.cpp $(SRCDIR)$@.cpp

clean:
	rm -f *.o $(BINDIR)*
