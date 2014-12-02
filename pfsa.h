/*
 * pfsa.h
 * anand, 22 mar 95.
 * templates, declarations and common includes for pfsa programs
 *
 * About the structure of the PFSA:
 *
 * The pfsa is a linked list of nodes of type NODE.  The first node is a dummy
 * header node.  In each node other than the first node, the state field stores
 * the state number, nsymbols stores the number of distinct symbols from that
 * state, ntrans stores the total number of transitions from that state, and
 * nvisits stores the total number of times this state is targeted by other
 * states.  Translist is a linked list of transitions again with a dummy
 * header.  Srclist (not implemented) is a list of states for which this
 * state is a target.  The mark field is used by various algorithms that use
 * this PFSA, see for instance beams.c.
 * state_list stores an ordered list of states from which this node was constructed.
 *   used only in beam search
 *   does not include the lowest state, which is stored in the state field.
 * In the root node (dummy header) of the pfsa, the state field should be
 * set to -1, ntrans stores the total number of states in the pfsa and
 * nsymbols stores the number of the maximum numbered state.
 *
 * 24 Sep 96: I fixed this so that the transition count is not affected
 * by transitions on delimiter symbols.  The trancnt is mainly used in
 * the Wallace formula to compute the cost of specifying one of N target
 * states, which is not necessary for the delimiter, since the target is
 * implicitly understood to be state 0.
 */
#ifndef PFSA_H
#define PFSA_H
#include <sys/types.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

/*extern char yytext[];*/

typedef struct {
   int type;
   union {
      int num;
      char *sym;
   } elem;
} *ELEM;

/* PONDY list of states, attached to each node, for comparing pfsa */
typedef struct stateListNode {
   int state;
   struct stateListNode *next;
} STATE;

/*
 * When creating transitions for lookahead automata, ensure that 
 * transitions on the same symbol from a state must be consecutive
 * links in the transition list. This is guaranteed if the list of
 * transitions is sorted in order of sym.
 */
typedef struct trans {
   struct node *target;
   int sym;			/* Really an index into Symtab where sym is */
   int freq;
   struct trans *next_tran;
} TRANS;

typedef struct source {
   struct node *source;
   int sym;
   int freq;
   struct source *next_src;
} SOURCE;

/*
 * ntrans below should be equal to nvisits by kirchoff's law.  The
 * field is preserved for historical reasons.
 * Note: transitions on delim are not counted for ntrans (MML reasons).
 */
typedef struct node {
   int state;
   int nsymbols;		/* # of distinct symbols from this state */
   int ntrans;			/* # of transitions from this state */
   int nvisits;			/* # of times this state is visited */
   u_char mark;			/* To mark the node as visited in traversals */
   TRANS *translist;            /* llist of targetnode, symbol (index to), and freq */
   SOURCE *srclist;		/* To get at all source nodes */
   STATE *state_list;     /* PONDY Ordered list of states merged into this node */
   struct node *nextnode;
} NODE;

/*
 * For efficiency reasons, we use the ntrans field of the root node
 * of the pfsa to store its number of states and the nsymbols field
 * to store the maximum state number.  We also store the total number
 * of unique arcs (NOT ON DELIMITER SYMBOLS) in the whole pfsa
 * in the nvisits field of the root node.
 *
 * 23/9/96: For the purposes of MML, the transition count only includes
 * the total number of unique transitions not on the delimiter symbol.
 */
#define incr_nodecnt(p) ((p)->ntrans++)
#define decr_nodecnt(p) ((p)->ntrans--)
#define nodecnt(p) ((p)->ntrans)
#define nstates(p) ((p)->ntrans)
#define incr_trancnt(p) ((p)->nvisits++)
#define decr_trancnt(p) ((p)->nvisits--)
#define trancnt(p) ((p)->nvisits)
#define setmaxstatenum(p,n) ((p)->nsymbols=n)
#define getmaxstatenum(p) ((p)->nsymbols)

#ifndef MAXNODES
# define MAXNODES 4096		/* Max nodes our dfa program can handle */
#endif

#define MAXSYMS 256
#define MAXSYMSIZE 64
typedef struct symbol {
   char label[MAXSYMSIZE];
   int freq;
} SYMBOL;
#define Sym(x)  (Symtab[x].label[0]=='\n'? "\\n" : Symtab[x].label)
#define Freq(x) (Symtab[x].freq)
#define DELIMITER 1		/* index of Delimiter symbol */
/*
 * Allow MAXSYMS possible symbols. Set of symbols is initially empty.
 */
#ifdef MAIN
   char *Author = "Author: Anand Raman, Department of Computer Science,\n"
                  "Massey University, Palmerston North, New Zealand.\n"
                  "Email: A.Raman@massey.ac.nz.\n"
                  "No warranties of any kind provided.\n"
                  "Studies using this program must cite it.\n";

   SYMBOL Symtab[MAXSYMS] = {
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, 
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0},
     {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0}, {"",0} 
  };
  int Lineno = 1;
  NODE *Pfsa = (NODE *) 0;
  int Debug = 0;
  int Graphplace = 0;
  char Delim = '\n';
  int Verbose = 0;
#else
  extern SYMBOL Symtab[];
  extern int Lineno, Debug, Graphplace, Verbose;
  extern char Delim;
  extern NODE *Pfsa;
#endif

extern int optind;
extern char option, *optarg;

extern int isatty(int);
/*
 * Return types of misc.c functions
 */
void buildpfsa(char []);
void resetmaxstatenum(NODE *pfsa);
void memerr(void);
void Perror(char *); 
void printpfsa(NODE *);
void writepfsa(FILE *, NODE *);
void output_pfsa(NODE *, char []);
void delpfsa(NODE *);
NODE *sortpfsa(NODE *);
NODE *createnode(void);
void statelimiterror(void);
NODE *addnode(NODE *, int);
void addtrans(NODE *, NODE *, int, int);
NODE *findnode(NODE *, int);
int addsym(char []);
int findsym(char []);
void printsyms(int *);
NODE *renumber(NODE *);
NODE *bf_renumber(NODE *pfsa);
NODE *copypfsa(NODE *);
void merge(NODE *, NODE *, NODE *);
int mealymerge(NODE *p1, NODE *p2);
NODE *mergecopy(NODE *, NODE *, NODE *);
NODE *newnode(NODE *);
TRANS *findtrans(TRANS *, int); 
int matchlen(NODE *, int *);
TRANS *lfindtrans(TRANS *, int *); 
int acceptable(NODE *, int *);
char *mkfname(char *, char *);
int *toks2syms(char *);
char *syms2toks(int *);
void onusr1(int);
void clearmarks(NODE *);
int isequiv(NODE *, NODE *);
   /* PONDY added this */
int isequiv_unrealised(NODE *proot, NODE *p1, NODE *p2, NODE *qroot, NODE *q1, NODE *q2);
NODE *trim(NODE *);
   /* PONDY added this */
STATE *consStateNode(int state, STATE *next);

/* int * functions equiv to the str * functions: Historical reasons */

int *intcpy(int *, int *);
int *intcat(int *, int);
int *intdup(int *);
int intcmp(int *, int *);
int intlen(int *);

/*
 * mml.c
 */
#define mml(p,q) mml_nfa(p,q)
#define mergedmml(p,q,r,s) mergedmml_nfa(p,q,r,s)

double mml_dfa(NODE *, double *);
double mml_nfa(NODE *, double *);
double mergedmml_dfa(NODE *pfsa, NODE *p, NODE *q, double oldmml);
double mergedmml_nfa(NODE *pfsa, NODE *p, NODE *q, double oldmml);

/*
 * opt.c
 */
void setfilenames(char *arg);

/*
 * dfa.c
 */
NODE *determinise(NODE *);
NODE *minimise(NODE *);

/*
 * Various optimisation functions
 */
NODE *ktail(int, char **);		/* ktail.c */
NODE *skstr(int, char **);		/* skstr.c */
NODE *beams(int, char **);		/* beams.c */
NODE *simba(int, char **);		/* simba.c */
#endif /*#ifndef PFSA_H*/