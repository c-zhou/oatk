CC=			gcc
CFLAGS=		-g -Wall -O3 -Wextra -Wno-unused-result -Wunused-parameter -fno-strict-aliasing
CPPFLAGS=
INCLUDES=	
OBJS=
PROG=		run_syncasm hmm_annotation path_finder path_to_fasta
PROG_EXTRA=
LIBS=		-lm -lz -lpthread

.PHONY:all extra clean depend
.SUFFIXES:.c .o

ifneq ($(asan),)
		CFLAGS+=-fsanitize=address
		LIBS+=-fsanitize=address
endif

all: $(PROG)

extra: all $(PROG_EXTRA)

debug: $(PROG)
debug: CFLAGS += -DDEBUG

run_syncasm: run_syncasm.c syncasm.c syncmer.c graph.c alignment.c sstream.c cov.c misc.c MurmurHash3.c kalloc.c kopen.c
		$(CC) $(CFLAGS) run_syncasm.c syncasm.c syncmer.c graph.c alignment.c sstream.c cov.c misc.c MurmurHash3.c kalloc.c kopen.c -o $@ -L. $(LIBS) $(INCLUDES)

hmm_annotation: annotation.c hmmannot.c misc.c kalloc.c kthread.c
		$(CC) $(CFLAGS) annotation.c hmmannot.c misc.c kalloc.c kthread.c -o $@ -L. $(LIBS) $(INCLUDES)

path_finder: path_finder.c path.c graph.c hmmannot.c misc.c kalloc.c kopen.c
		$(CC) $(CFLAGS) path_finder.c path.c graph.c hmmannot.c misc.c kalloc.c kopen.c -o $@ -L. $(LIBS) $(INCLUDES)

path_to_fasta: path_to_fasta.c path.c graph.c misc.c kalloc.c kopen.c
		$(CC) $(CFLAGS) path_to_fasta.c path.c graph.c misc.c kalloc.c kopen.c -o $@ -L. $(LIBS) $(INCLUDES)

clean:
		rm -fr *.o a.out $(PROG) $(PROG_EXTRA)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

