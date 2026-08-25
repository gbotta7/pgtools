CC=			gcc
CFLAGS=		-g -Wall -O2
CPPFLAGS=
INCLUDES=	
OBJS=		bed.o count.o shtab.o kthread.o parser.o sys.o utils.o
PROG=		pgkmc
LIBS=		-lm -lz -lpthread

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

$(PROG):$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

main.o:    main.c    shtab.h utils.h ketopt.h sys.h
bed.o:     bed.c     bed.h khashl.h utils.h
count.o:   count.c   bed.h shtab.h utils.h kseq.h kthread.h parser.h
shtab.o:   shtab.c   shtab.h utils.h khashl.h parser.h kseq.h
kthread.o: kthread.c kthread.h
parser.o:  parser.c  parser.h kseq.h utils.h
sys.o:     sys.c     sys.h
utils.o:   utils.c   shtab.h utils.h