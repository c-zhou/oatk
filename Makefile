CC=			gcc
#CFLAGS=		-g -Wall -O3 -Wextra -Wno-unused-result -Wunused-parameter -fno-strict-aliasing
CFLAGS=		-g -Wall -O3 -Wno-unused-function
CPPFLAGS=
INCLUDES=	
OBJS=
PROG=		syncasm hmm_annotation pathfinder path_to_fasta oatk
PROG_EXTRA=
LIBS=		-lm -lz -lpthread
DESTDIR=	~/bin

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

syncasm: run_syncasm.c syncasm.c syncmer.c syncerr.c levdist.c graph.c alignment.c sstream.c misc.c kthread.c kalloc.c kopen.c
		$(CC) $(CFLAGS) -DSYNCASM_MAIN run_syncasm.c syncasm.c syncmer.c syncerr.c levdist.c graph.c alignment.c sstream.c misc.c kthread.c kalloc.c kopen.c -o $@ -L. $(LIBS) $(INCLUDES)

hmm_annotation: hmm_annotation.c hmmannot.c misc.c kalloc.c kthread.c
		$(CC) $(CFLAGS) -DANNOTATION_MAIN hmm_annotation.c hmmannot.c misc.c kalloc.c kthread.c -o $@ -L. $(LIBS) $(INCLUDES)

pathfinder: path_finder.c syncasm.c syncmer.c syncerr.c levdist.c path.c graph.c hmmannot.c alignment.c sstream.c misc.c kalloc.c kopen.c kthread.c
		$(CC) $(CFLAGS) -DPATHFINDER_MAIN path_finder.c syncasm.c syncmer.c syncerr.c levdist.c path.c graph.c hmmannot.c alignment.c sstream.c misc.c kalloc.c kopen.c kthread.c -o $@ -L. $(LIBS) $(INCLUDES)

path_to_fasta: path_to_fasta.c path.c graph.c hmmannot.c misc.c kalloc.c kopen.c
		$(CC) $(CFLAGS) path_to_fasta.c path.c graph.c hmmannot.c misc.c kalloc.c kopen.c -o $@ -L. $(LIBS) $(INCLUDES)

oatk: oatk.c run_syncasm.c hmm_annotation.c path_finder.c hmmannot.c syncasm.c syncmer.c syncerr.c levdist.c path.c graph.c alignment.c sstream.c misc.c kalloc.c kopen.c kthread.c
		$(CC) $(CFLAGS) oatk.c run_syncasm.c hmm_annotation.c path_finder.c hmmannot.c syncasm.c syncmer.c syncerr.c levdist.c path.c graph.c alignment.c sstream.c misc.c kalloc.c kopen.c kthread.c -o $@ -L. $(LIBS) $(INCLUDES)

clean:
		rm -fr *.o a.out $(PROG) $(PROG_EXTRA)

install:
		cp $(PROG) $(DESTDIR)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

