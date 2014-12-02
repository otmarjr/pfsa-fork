/*
 * Miscellaneous useful functions for pfsa programs.
 * anand, 22 March 95
 * 1 May 1996, added getannos() function to prettify graphplace output.
 * 8 Oct 1996, added bf_renumber() for breadth first renumbering of nodes.
 */
#include "pfsa.h"
#include <string.h>
#include <signal.h>
#ifndef MISC_C
#define MISC_C

/*
 * The following 2 structures and four queue manipulation functions are
 * for the breadth-first renumbering function bf_renumber().  They are
 * not used anywhere else at the moment.
 */
typedef struct cell {
    NODE *nodep;
    struct cell *next;
} CELL;

typedef struct {
    CELL *front, *rear;
} NODEQ;

static NODEQ *createnodeq(void);
static void add2q(NODEQ *nq, NODE *p);
static NODE *dequeue(NODEQ *nq);
static void deleteq(NODEQ *nq);
NODE *isort(NODE *p);
STATE *consStateNode(int state, STATE *next);

void buildpfsa(char specsfile[]) {
    extern char Infile[];
    char temp[BUFSIZ];
    extern FILE *yyin;
    int parse_err, yyparse();

    strcpy(temp, Infile);
    strcpy(Infile, specsfile);

    if (strcmp(Infile, "-")) {
        yyin = fopen(Infile, "r");
        if (!yyin)
            Perror(Infile);
    }

    /*
     * Make delimiter the first symbol in the symbol table, Note Symtab[0]
     * is a sentinel, can't use it, so symtab[1] is the first
     */
    sprintf(Symtab[DELIMITER].label, "%c", Delim);

    /*
     * Build the pfsa from the specs in the file.
     */
    Pfsa = (NODE *) calloc(1, sizeof (NODE)); /* set up header */
    if (!Pfsa)
        memerr();
    Pfsa->state = -1;
    if ((parse_err = yyparse())) {
        fprintf(stderr, "yyparse() returned %d\n", parse_err);
        exit(parse_err);
    }
    if (yyin != stdin)
        (void) fclose(yyin);
    strcpy(Infile, temp);
}

/* This function - superceded by macro, 16/5/96
   
int nstates(pfsa)
NODE *pfsa;
{
   int i = 0;

   while ((pfsa = pfsa->nextnode))
      i++;
   return i;
}

 */

void resetmaxstatenum(NODE *pfsa) /* 16/6/96 */ {
    NODE *p;

    setmaxstatenum(pfsa, 0);
    p = pfsa;
    while ((p = p->nextnode))
        if (p->state > getmaxstatenum(pfsa))
            setmaxstatenum(pfsa, p->state);
}

NODE *createnode() {
    NODE *instance;

    instance = (NODE *) calloc(1, sizeof (NODE));
    if (!instance)
        memerr();
    instance->translist = (TRANS *) calloc(1, sizeof (TRANS));
    instance->srclist = (SOURCE *) calloc(1, sizeof (SOURCE));
    if (!instance->translist || !instance->srclist)
        memerr();
    instance->translist->sym = instance->srclist->sym = -1;
    return instance;
}

void statelimiterror() {
    fprintf(stderr, "More than %d nodes in this PFSA!\n", MAXNODES);
    fprintf(stderr, "Change set size in the DFA conversion program\n");
    fprintf(stderr, "and the statelist array size in copypfsa()\n");
    exit(1);
}

NODE *newnode(NODE *pfsa) /* Create a new node at the end of list pfsa */ {
    NODE *p;

    p = pfsa;
    while (p->nextnode)
        p = p->nextnode;

    p->nextnode = createnode();
    p->nextnode->state = p->state + 1;

    if (p->nextnode->state == MAXNODES)
        statelimiterror();
    incr_nodecnt(pfsa);
    return p->nextnode;
}

/*
 * add the node for the state numbered "state". If it already exists
 * then return a pointer to it. Nodes are kept sorted in order of state 
 * number. Remember, the first node is a dummy header guaranteed to exist.
 */
NODE *addnode(NODE *pfsa,
        int state) {
    NODE *p, *newInstance;

    p = pfsa;
    while (p) {
        if (!p->nextnode || p->nextnode->state > state) {
            if (!p->nextnode)
                setmaxstatenum(pfsa, state);
            if (incr_nodecnt(pfsa) >= MAXNODES)
                statelimiterror();
            newInstance = createnode();
            newInstance->state = state;
            newInstance->nextnode = p->nextnode;
            p->nextnode = newInstance;
            return newInstance;
        } else if (p->nextnode->state == state)
            return p->nextnode;
        else
            p = p->nextnode;
    }
    return (NODE *) 0; /* Not reached */
}

/* I realise now that I should have ordered both translist and srclist
 * using the state number as a secondary key in addition to using sym
 * as the primary key.  Hopefully, I should re-implement this soon
 * when I get some time. 15/01/97
 */
void addtrans(NODE *src, NODE *dst,
        int sym,
        int freq) {
    TRANS *tp, *newtp;
    SOURCE *sp, *newsp;
    int newsym = 1;

    tp = src->translist;
    while (tp) {
        if (!tp->next_tran || tp->next_tran->sym > sym) {
            if (tp->sym == sym)
                newsym = 0;
            newtp = (TRANS *) calloc(1, sizeof (*tp->next_tran));
            if (!newtp)
                memerr();
            newtp->sym = sym;
            newtp->freq = freq;
            newtp->target = dst;
            newtp->next_tran = tp->next_tran;
            tp->next_tran = newtp;
            if (Symtab[sym].label[0] != Delim)
                incr_trancnt(Pfsa);
            break;
        } else if (tp->next_tran->sym == sym && tp->next_tran->target == dst) {
            newsym = 0;
            tp->next_tran->freq += freq;
            break;
        } else
            tp = tp->next_tran;
    }
    /* Add a back pointer to the source node from the dest node on sym
     */
    sp = dst->srclist;
    while (sp) {
        if (!sp->next_src || sp->next_src->sym > sym) {
            newsp = (SOURCE *) calloc(1, sizeof (*sp->next_src));
            if (!newsp)
                memerr();
            newsp->sym = sym;
            newsp->freq = freq;
            newsp->source = src;
            newsp->next_src = sp->next_src;
            sp->next_src = newsp;
            break;
        } else if (sp->next_src->sym == sym && sp->next_src->source == src) {
            sp->next_src->freq += freq;
            break;
        } else
            sp = sp->next_src;
    }

    if (newsym)
        src->nsymbols++;
    src->ntrans += freq;
    dst->nvisits += freq;
    Symtab[sym].freq += freq;
}

NODE *findnode(NODE *p,
        int state) /* find a node "state" in a nodelist p */ {
    while (p && p->state != state)
        p = p->nextnode;
    return p;
}

TRANS *findtrans(TRANS *tp,
        int sym) /* find 1st transition on "sym" in  list t */ {
    while (tp && tp->sym != sym)
        tp = tp->next_tran;
    return tp;
}

/*
 * If the label is already known return its index in the symbol table.
 * Otherwise, create a new label by that name and return its index.
 * Henceforth, that label will be known by its index.
 * Note: Symbol zero is reserved as sentinel.
 */
int addsym(char label[]) {
    int i;

    for (i = 1; i < MAXSYMS; i++) {
        if (!Symtab[i].label[0]) {
            strncpy(Symtab[i].label, label, MAXSYMSIZE - 1);
            return i;
        }
        if (!strncmp(Symtab[i].label, label, MAXSYMSIZE))
            return i;
    }
    fprintf(stderr, "Too many symbols: Consider increasing Sym table size\n");
    exit(1);
}

void printsyms(int *syms) {
    while (*syms) {
        printf("%s", Symtab[*syms].label);
        if (Symtab[*syms].label[0] == Delim)
            break;
        printf(":");
        syms++;
    }
    if (!*syms)
        printf("*UNDELIMITED STRING*\n");
}

/*
 * This is the same as addsym, but it is an error to not find sym
 * already in the table
 * Note: Symbol zero is reserved as sentinel.
 */
int findsym(char label[]) {
    int i;

    for (i = 1; i < MAXSYMS && Symtab[i].label[0]; i++)
        if (!strncmp(Symtab[i].label, label, MAXSYMSIZE))
            return i;
    fprintf(stderr, "Couldn't find symbol %s in table\n",
            label[0] == '\n' ? "\\n" : label);
    exit(1);
}

/*
 * How far into node p in the pfsa can you insert the string of symbols s?
 */
int matchlen(NODE *p, int *s) {
    TRANS *t;
    int len, bestlen;

    /*
     * It seems lfindtrans should be used here, but findtrans is actually
     * the one - the code is correct if you check it through.
     */
    t = findtrans(p->translist->next_tran, *s);
    if (!t || !*s)
        return 0;
    if (Symtab[*s].label[0] == Delim)
        return 1;
    bestlen = 1 + matchlen(t->target, s + 1);
    while (t->next_tran && t->next_tran->sym == *s) {
        t = t->next_tran;
        len = 1 + matchlen(t->target, s + 1);
        if (bestlen < len)
            bestlen = len;
    }
    return bestlen;
}

/*
 * lfindtrans() is the lookahead version of findtrans. There might
 * be more than one transition on the same symbol from a state, so the
 * one which promises the most consumption of symbols in s is chosen
 */
TRANS *lfindtrans(TRANS *t, int *s) {
    TRANS *best;

    while (t && t->sym != *s) /* find first transition on sym *s */
        t = t->next_tran;
    if (!t || !*s)
        return (TRANS *) 0;

    /*
     * If the current symbol is the delimiter, we have found the only
     * transition - return it.
     */
    if (Symtab[*s].label[0] == Delim)
        return t;

    best = t;
    while (t->next_tran && t->next_tran->sym == *s) {
        t = t->next_tran;
        if (matchlen(best->target, s + 1) < matchlen(t->target, s + 1))
            best = t;
    }
    return best;
}

int acceptable(NODE *q, int *s) {
    TRANS *tp;

    while (*s) {
        tp = lfindtrans(q->translist, s);
        if (!tp)
            return 0;
        else
            q = tp->target;
        ++s;
    }
    return 1;
}

/*
 * Renumber the nodes of the pfsa such that the states are sequentially
 * numbered. Call this when there a lot of holes in the pfsa after a bout
 * of optimisation. This function became a lot smaller as of 3 Aug 1995
 * when t->targets were changed to be node pointers. Now we don't have to
 * renumber the targets, since they all point to nodes.
 */
NODE *renumber(NODE *pfsa) {
    NODE *p;
    int prev = -1;

    if (Debug)
        fprintf(stderr, "Renumbering...\n");
    for (p = pfsa; p; p = p->nextnode)
        p->state = prev++;
    setmaxstatenum(pfsa, prev - 1);
    return pfsa;
}

static NODEQ *createnodeq() {
    NODEQ *nq;

    nq = (NODEQ *) calloc(1, sizeof (NODEQ));
    if (!nq)
        memerr();
    nq->front = nq->rear = (CELL *) NULL;
    return nq;
}

static void add2q(NODEQ *nq, NODE *p) {
    CELL *nc;

    nc = (CELL *) calloc(1, sizeof (CELL));
    if (!nc)
        memerr();
    nc->nodep = p;
    if (!nq->rear)
        nq->rear = nq->front = nc;
    else {
        nq->rear->next = nc;
        nq->rear = nc;
    }
}

static NODE *dequeue(NODEQ *nq) {
    CELL *nc;
    NODE *p;

    if (!nq->front)
        return (NODE *) NULL;
    nc = nq->front;
    nq->front = nq->front->next;
    if (!nq->front)
        nq->rear = nq->front;
    p = nc->nodep;
    free((void *) nc);
    return p;
}

static void deleteq(NODEQ *nq) {
    CELL *nc;

    while (nq->front) {
        nc = nq->front;
        nq->front = nq->front->next;
        free((void *) nc);
    }
    free((void *) nq);
}

NODE *isort(NODE *p) /* insertion sort to improve output */ {
    NODE *q, *r, *s;

    if (!p->nextnode)
        return p;
    p->nextnode = isort(p->nextnode);
    if (p->nextnode->state >= p->state)
        return p;
    r = p;
    q = p->nextnode;
    while (r->nextnode && r->nextnode->state < p->state)
        r = r->nextnode;
    s = r->nextnode;
    r->nextnode = p;
    p->nextnode = s;
    return q;
}

void printq(NODEQ *nodeq) {
    CELL *nc;

    nc = nodeq->front;

    while (nc) {
        fprintf(stderr, "%d->", nc->nodep->state);
        nc = nc->next;
    }
    fprintf(stderr, "[END]\n");
}

NODE *bf_renumber(NODE *pfsa) {
    int prev = -1;
    NODEQ *nodeq;
    NODE *p;
    TRANS *tp;

    clearmarks(pfsa);
    nodeq = createnodeq();
    pfsa->nextnode->mark = 1;
    add2q(nodeq, pfsa->nextnode); /* state 0 */

    while ((p = dequeue(nodeq)) != NULL) {
        p->state = ++prev;
        for (tp = p->translist->next_tran; tp; tp = tp->next_tran)
            if (!tp->target->mark) {
                tp->target->mark = 1;
                add2q(nodeq, tp->target);
            }
        /* printq(nodeq); */
    }
    deleteq(nodeq);
    pfsa->nextnode = isort(pfsa->nextnode);
    return pfsa;
}

void delpfsa(NODE *pfsa) {
    NODE *p;
    TRANS *tp;
    SOURCE *sp;
    STATE *stp;

    while (pfsa) {
        while (pfsa->translist) {
            tp = pfsa->translist;
            pfsa->translist = tp->next_tran;
            free((void *) tp);
        }
        while (pfsa->srclist) {
            sp = pfsa->srclist;
            pfsa->srclist = sp->next_src;
            free((void *) sp);
        }
        /* PONDY remove the state list */
        while (pfsa->state_list) {
            stp = pfsa->state_list;
            pfsa->state_list = stp->next;
            free((void *) stp);
        }
        p = pfsa;
        pfsa = pfsa->nextnode;
        free((void *) p);
    }
}

void ps(NODE *p) {
    char *label;
    SOURCE *sp;

    while (p && p->nextnode) {
        sp = p->nextnode->srclist->next_src;
        fprintf(stderr, "%d:", p->nextnode->state);
        while (sp) {
            label = Symtab[sp->sym].label;
            if (label[0] == '\n')
                label = (char*)"\\n";
            fprintf(stderr, "->(%d, %s,%d)", sp->source->state, label, sp->freq);
            sp = sp->next_src;
        }
        fprintf(stderr, ".\n");
        sp = p->nextnode->srclist->next_src;
        fprintf(stderr, "%d:", p->nextnode->state);
        while (sp) {
            label = Symtab[sp->sym].label;
            if (label[0] == '\n')
                label = (char*)"\\n";
            fprintf(stderr, "->(%p)", sp);
            sp = sp->next_src;
        }
        fprintf(stderr, ".\n");
        p = p->nextnode;
    }
}

/*
 * After a change to the structure of the pfsa, this function turned out
 * not to be as trivial as it was since the target pointers in the trans
 * lists needed to be updated to point into nodes in the new pfsa.  Thus
 * I have now changed this into a two-pass function.
 */
NODE *copypfsa(NODE *oldpfsa) {
    NODE *newpfsa, *newp, *oldp, *statelist[MAXNODES];
    TRANS *newtp, *oldtp;
    SOURCE *newsp, *oldsp;
    STATE *rest, *tail; /*PONDY for copying the statelist */

    /*
     * First create all the state nodes.  Then fill in the transitions
     */
    newpfsa = (NODE *) calloc(1, sizeof (NODE));
    if (!newpfsa)
        memerr();
    memcpy(newpfsa, oldpfsa, sizeof (NODE));

    oldp = oldpfsa;
    newp = newpfsa;

    /* PASS 1: Create the list of nodes and remember their addresses
     */
    while (oldp->nextnode) {
        newp->nextnode = (NODE *) calloc(1, sizeof (NODE));
        if (!newp->nextnode)
            memerr();
        memcpy(newp->nextnode, oldp->nextnode, sizeof (NODE));
        statelist[newp->nextnode->state] = newp->nextnode;
        /*
         * PONDY Create newp->state_list list of state, if there is one there.
         */
        rest = oldp->nextnode->state_list;
        if (rest) {
            newp->nextnode->state_list = consStateNode(rest->state, NULL);
            tail = newp->nextnode->state_list;
            rest = rest->next;
            while (rest) {
                tail->next = consStateNode(rest->state, NULL);
                tail = tail->next;
                rest = rest->next;
            }
        }
        oldp = oldp->nextnode;
        newp = newp->nextnode;
    }


    /* PASS 2: Now go through and create the transition and source lists
     *         using node addresses from the NEW pfsa node list.
     */
    oldp = oldpfsa;
    newp = newpfsa;
    while (oldp->nextnode) {
        /*
         * Create newp->nextnode->translist
         */
        newp->nextnode->translist = (TRANS *) calloc(1, sizeof (TRANS));
        if (!newp->nextnode->translist)
            memerr();
        oldtp = oldp->nextnode->translist;
        newtp = newp->nextnode->translist;
        while (oldtp->next_tran) {
            newtp->next_tran = (TRANS *) calloc(1, sizeof (TRANS));
            if (!newtp->next_tran)
                memerr();
            memcpy(newtp->next_tran, oldtp->next_tran, sizeof (TRANS));
            newtp->next_tran->target = statelist[oldtp->next_tran->target->state];
            oldtp = oldtp->next_tran;
            newtp = newtp->next_tran;
        }

        /*
         * Create newp->nextnode->srclist
         */
        newp->nextnode->srclist = (SOURCE *) calloc(1, sizeof (SOURCE));
        if (!newp->nextnode->srclist)
            memerr();
        oldsp = oldp->nextnode->srclist;
        newsp = newp->nextnode->srclist;
        while (oldsp->next_src) {
            newsp->next_src = (SOURCE *) calloc(1, sizeof (SOURCE));
            if (!newsp->next_src)
                memerr();
            memcpy(newsp->next_src, oldsp->next_src, sizeof (SOURCE));
            newsp->next_src->source = statelist[oldsp->next_src->source->state];
            oldsp = oldsp->next_src;
            newsp = newsp->next_src;
        }

        oldp = oldp->nextnode;
        newp = newp->nextnode;
    }
    return newpfsa;
}

/*
 * This returns 1 if the merge of nodes p1 and p2 in the pfsa will result
 * in a mealy type machine - i.e, non-lookahead automaton.  This is
 * critical for the Wallace & Georgeff 1984 formula to continue to apply.
 */
int mealymerge(NODE *p1, NODE *p2) {
    TRANS *tp1, *tp2;

    if (p1 == p2)
        return 1;
    /*
     * A merge will result in a non-mealy machine if there there are two
     * transitions t1 and t2 from nodes p1 and p2 on the same symbol which
     * go to different targets *after the merge*
     */
    tp1 = p1->translist->next_tran;
    tp2 = p2->translist->next_tran;
    while (tp1 && tp2) {
        if (tp1->sym < tp2->sym)
            tp1 = tp1->next_tran;
        else if (tp2->sym < tp1->sym)
            tp2 = tp2->next_tran;
        else { /* same symbols */
            if (tp1->target == tp2->target ||
                    ((tp1->target == p1 && tp2->target == p2) ||
                    (tp1->target == p2 && tp2->target == p1))) {
                tp1 = tp1->next_tran;
                tp2 = tp2->next_tran;
            } else
                return 0;
        }
    }
    return 1;
}

/* Heuristic to determine if two PFSA are equal in O(n) time.
 * 1. They have the same number of states
 * 2. The states have the same numbers in the same order
 * 3. States with the same numbers have the same transitions
 */
int isequiv(NODE *pfsa1, NODE *pfsa2) {
    NODE *p, *q;
    TRANS *tp, *tq, *t;
    int foundtrans;

    if (nstates(pfsa1) != nstates(pfsa2))
        return 0;
    p = pfsa1->nextnode;
    q = pfsa2->nextnode;
    while (p && q) {
        if (p->state != q->state)
            return 0;
        p = p->nextnode;
        q = q->nextnode;
    }
    p = pfsa1->nextnode;
    q = pfsa2->nextnode;
    while (p && q) {
        tp = p->translist->next_tran;
        tq = q->translist->next_tran;

        foundtrans = 0;
        for (t = tp; t && t->sym == tq->sym; t = t->next_tran) {
            if (tq->target->state == t->target->state && tq->freq == t->freq) {
                foundtrans = 1;
                break;
            }
        }
        if (!foundtrans)
            return 0;

        foundtrans = 0;
        for (t = tq; t && t->sym == tp->sym; t = t->next_tran) {
            if (tp->target->state == t->target->state && tp->freq == t->freq) {
                foundtrans = 1;
                break;
            }
        }
        if (!foundtrans)
            return 0;
        p = p->nextnode;
        q = q->nextnode;
    }
    return 1;
}

/* PONDY
 * Heuristic to determine if two ***unrealised*** PFSA are equal in O(n) time.
 * p1, p2 are states that should be merged in the pfsa starting at proot
 * q1, q2 are states that should be merged in the pfsa starting at qroot
 * 1. They have the same number of states
 * 2. The states have the same numbers in the same order
 * 3. States with the same numbers have the same transitions
 * Assumes they were constructed by merging from the same original PFSA
 * code is long, with repeated bits for each different case - more efficient
 */
int isequiv_unrealised(NODE *proot, NODE *p1, NODE *p2, NODE *qroot, NODE *q1, NODE *q2) {
    NODE *p, *q, *temp;
    STATE *restp, *restq, *restp2, *restq2;
    int pstate, qstate, p2state, q2state;

    if (p1->state > p2->state) { /*shouldn't be in reverse order */
        temp = p1;
        p1 = p2;
        p2 = temp;
    }

    if (q1->state > q2->state) { /*shouldn't be in reverse order */
        temp = q1;
        q1 = q2;
        q2 = temp;
    }
    if (nstates(proot) != nstates(qroot))
        return 0;
    if (proot == qroot) { /* both pfsa derived from the same root*/
        if (p1 != p2 || q1 != q2) /* merge nodes must be the same pair of nodes */
            return 0;
    }
    /* Assumes that nodes are stored in a canonical order of state number */
    p = proot->nextnode;
    q = qroot->nextnode;
    while (p && q) {
        if (p->state != q->state)
            return 0;
        restp = p->state_list;
        restq = q->state_list;
        if (p != p1 && q != q1) { /* easy case - no merge nodes involved */
            while (restp && restq) { /* just step along all the states */
                if (restp->state != restq->state)
                    return 0;
                restp = restp->next;
                restq = restq->next;
            }
        } else if (p == p1 && q != q1) { /* only p has a merge node */
            restp2 = p2->state_list;
            p2state = p2->state;
            /* While still some states left to compare */
            while ((restp || p2state || restp2) && restq) {
                /* Find the next state in the p node */
                if (p2state && (!restp || p2state < restp->state)) {
                    if (p2state != restq->state)
                        return 0;
                    p2state = 0;
                } else if (p2state || !restp2 || (restp && (restp->state < restp2->state))) {
                    /* p2state > rest1, or comparing the two lists, and restp1 is first */
                    if (restp->state != restq->state)
                        return 0;
                    restp = restp->next;
                } else { /* comparing lists: step along p2->state_list */
                    if (restp2->state != restq->state)
                        return 0;
                    restp2 = restp2->next;
                }
                restq = restq->next;
            }
            /* If there are any states left on either side, then lose */
            if (restp || p2state || restp2 || restq)
                return 0;
        } else if (p != p1 && q == q1) { /*  only q has a merge node */
            restq2 = q2->state_list;
            q2state = q2->state;
            /* While still some states left to compare */
            while (restp && (restq || q2state || restq2)) {
                pstate = restp->state;
                /* Find the next state in the q node */
                if (q2state && (!restq || q2state < restq->state)) {
                    if (restp->state != q2state)
                        return 0;
                    q2state = 0; /* gobble q2->state */
                } else if (q2state || !restq2 || (restq && (restq->state < restq2->state))) {
                    /* q2state > rest1, or comparing the two lists, and restq1 is first */
                    if (restp->state != restq->state)
                        return 0;
                    restq = restq->next;
                } else { /* comparing lists: step along q2->state_list */
                    if (restp->state != restq2->state)
                        return 0;
                    restq2 = restq2->next;
                }
            }
            /* If there are any states left on either side, then lose */
            if (restp || restq || q2state || restq2)
                return 0;
        } else { /*  p and q have merge nodes */
            restp2 = p2->state_list;
            restq2 = q2->state_list;
            p2state = p2->state;
            q2state = q2->state;
            /* While still some states left to compare */
            while ((restp || p2state || restp2) && (restq || q2state || restq2)) {
                /* Find the next state in the p node */
                if (p2state && (!restp || p2state < restp->state)) {
                    pstate = p2state; /* gobble p2->state */
                    p2state = 0;
                } else if (p2state || !restp2 || (restp && (restp->state < restp2->state))) {
                    /* p2state > rest1, or comparing the two lists, and restp1 is first */
                    pstate = restp->state;
                    restp = restp->next;
                } else { /* comparing lists: step along p2->state_list */
                    pstate = restp2->state;
                    restp2 = restp2->next;
                }
                /* Find the next state in the q node */
                if (q2state && (!restq || q2state < restq->state)) {
                    qstate = q2state; /* gobble q2->state */
                    q2state = 0;
                } else if (q2state || !restq2 || (restq && (restq->state < restq2->state))) {
                    /* q2state > rest1, or comparing the two lists, and restq1 is first */
                    qstate = restq->state;
                    restq = restq->next;
                } else { /* comparing lists: step along q2->state_list */
                    qstate = restq2->state;
                    restq2 = restq2->next;
                }
                /* Check the two states */
                if (qstate != pstate)
                    return 0;
            }
            /* If there are any states left on either side, then lose */
            if (restp || p2state || restp2 || restq || q2state || restq2)
                return 0;
        }
        p = p->nextnode;
        q = q->nextnode;

        /* skip over p2 and q2 - they were checked before */
        if (p == p2)
            p = p->nextnode;
        if (q == q2)
            q = q->nextnode;

    }

    return 1;
}

/* PONDY Constructs a node in a linked list of states */
STATE *consStateNode(int state, STATE *next) {
    STATE *ans;
    ans = (STATE *) calloc(1, sizeof (STATE));
    if (!ans)
        memerr();
    ans->state = state;
    ans->next = next;
    return ans;
}

/*
 * Merge nodes p and q in the pfsa.
 * This is currently an O(n) algorithm. As soon as I get some time I should
 * convert this to an O(1) algorithm.
 *
 * This is currently the only place where an unrealised pfsa is realised.
 *
 * NOTE: p2 is merged into p1 not the other way around.  After the merge
 *       you can't refer to p2.  The process is symmetric, but the 
 *       nomenclature is not.
 */
void merge(NODE *pfsa,NODE *p1,NODE *p2) {
    NODE *p;
    TRANS *tp, *tp1, *tp2, *temptp;
    SOURCE *sp, *sp1, *sp2, *tempsp;
    if (p1 == p2)
        return;

    /* Any node that has p2 as a target needs to be changed to have p1
     * instead.  These nodes are in p2->srclist.  Likewise, any node
     * that has p2 as a source needs to be changed to have p1 instead.
     * These nodes are in p2->translist
     *
     * (Actually works out faster to traverse the node list rather than
     * sourcelist and translist since they tend to have duplicates.)
     */
    for (p = pfsa->nextnode; p; p = p->nextnode) {
        for (tp = p->translist->next_tran; tp; tp = tp->next_tran)
            if (tp->target == p2)
                tp->target = p1;
        for (sp = p->srclist->next_src; sp; sp = sp->next_src)
            if (sp->source == p2)
                sp->source = p1;
    }

    /* PONDY  Merge the state list of p2 into p1 must be new copy
       (can't share because of need to dispose)  */

    if (!p1->state_list) { /* p1 was a one state node */
        p1->state_list = consStateNode(p2->state, p2->state_list);
    } else {
        STATE *list2, *rest1, *temp; /* PONDY for merging the state lists */

        list2 = consStateNode(p2->state, p2->state_list);

        /* p1 has a statelist of its own, but we know that
         * p1->state < list2->state.
         * Unlink each node in list2 and link it into p1->state_list
         * in the right place sequentially.
         */
        rest1 = p1->state_list;
        while (list2) {
            while (rest1->next && rest1->next->state < list2->state)
                rest1 = rest1->next; /* found either gap or End */
            if (!rest1->next) {
                rest1->next = list2;
                break;
            }
            temp = list2->next;
            list2->next = rest1->next;
            rest1->next = list2;
            list2 = temp;
        }
    }

    /* Merge the transition list of p2 into p1.  This may result in
     * duplicate transitions - eg (a->4)..(a->3)..(a,4), but this will
     * be collapsed in the next step where we merge duplicate
     * transitions.
     */
    tp1 = p1->translist;
    tp2 = p2->translist;
    while (tp2->next_tran) {
        temptp = tp2->next_tran;
        tp2->next_tran = temptp->next_tran; /* unlink */
        while (tp1->next_tran && tp1->next_tran->sym < temptp->sym)
            tp1 = tp1->next_tran;
        if (tp1->next_tran && tp1->next_tran->sym == temptp->sym &&
                tp1->next_tran->target == temptp->target) {
            tp1->next_tran->freq += temptp->freq;
            free((void *) temptp);
            if (Symtab[tp1->next_tran->sym].label[0] != Delim)
                decr_trancnt(pfsa);
        } else {
            if (!tp1->next_tran ||
                    (tp1->next_tran && tp1->next_tran->sym > temptp->sym))
                p1->nsymbols++;
            temptp->next_tran = tp1->next_tran;
            tp1->next_tran = temptp;
        }
    }
    free((void *) p2->translist);

    /* Likewise for p2's source list.
     */
    sp1 = p1->srclist;
    sp2 = p2->srclist;
    while (sp2->next_src) {
        tempsp = sp2->next_src;
        sp2->next_src = tempsp->next_src; /* unlink */
        while (sp1->next_src && sp1->next_src->sym < tempsp->sym)
            sp1 = sp1->next_src;
        if (sp1->next_src && sp1->next_src->sym == tempsp->sym &&
                sp1->next_src->source == tempsp->source) {
            sp1->next_src->freq += tempsp->freq;
            free((void *) tempsp);
        } else {
            tempsp->next_src = sp1->next_src;
            sp1->next_src = tempsp;
        }
    }
    free((void *) p2->srclist);

    /* The previous steps would have caused duplicate transitions and
     * sources in the lists, eg, (1,a)->(2,a)->... may become
     * (1,a)->(1,a)->..  on merging 1 and 2.  Merge all such duplicate
     * items.
     */
    for (p = pfsa->nextnode; p; p = p->nextnode) {
        for (tp = p->translist; tp; tp = tp->next_tran) {
            tp1 = tp;
            while (tp1->next_tran && tp->sym == tp1->next_tran->sym) {
                if (tp->target == tp1->next_tran->target) {
                    tp->freq += tp1->next_tran->freq;
                    if (Symtab[tp1->next_tran->sym].label[0] != Delim)
                        decr_trancnt(pfsa);
                    temptp = tp1->next_tran;
                    tp1->next_tran = temptp->next_tran;
                    free((void *) temptp);
                    continue;
                }
                tp1 = tp1->next_tran;
            }
        }
        for (sp = p->srclist; sp; sp = sp->next_src) {
            sp1 = sp;
            while (sp1->next_src && sp->sym == sp1->next_src->sym) {
                if (sp->source == sp1->next_src->source) {
                    sp->freq += sp1->next_src->freq;
                    tempsp = sp1->next_src;
                    sp1->next_src = tempsp->next_src;
                    free((void *) tempsp);
                    continue;
                }
                sp1 = sp1->next_src;
            }
        }
    }

    p1->nvisits += p2->nvisits;
    p1->ntrans += p2->ntrans;

    /*
     * Step 3: Unlink p2
     */
    for (p = pfsa; p->nextnode; p = p->nextnode) {
        if (p->nextnode == p2) {
            p->nextnode = p2->nextnode;
            free((void *) p2);
            break;
        }
    }
    if (p2->state == getmaxstatenum(pfsa))
        resetmaxstatenum(pfsa);
    decr_nodecnt(pfsa);
}

/*
 * mergecopy is the same as merge - in fact it calls merge() to do the job,
 * but it returns a copy of the merged pfsa, leaving the original pfsa
 * unchanged.
 * Effect: Merges nodes q and r in a copy of the PFSA p and returns it.
 */
NODE *mergecopy(NODE *p, NODE *q, NODE *r) {
    int state1, state2;
    NODE *newpfsa;

    state1 = q->state;
    state2 = r->state;
    newpfsa = copypfsa(p);
    q = findnode(newpfsa, state1);
    r = findnode(newpfsa, state2);
    merge(newpfsa, q, r);
    return newpfsa;
}

/*
 * Printpfsa prints out the pfsa in its internal representation:
 *
 *   Node -> list of transitions out of that node
 *    |
 *    V
 * nextnode
 *   ...
 * Note that the first element of every list is a dummy for the header.
 */
void printlist(TRANS *tp) {
    char *label;

    while (tp) {
        label = Symtab[tp->sym].label;
        if (label[0] == '\n')
            label = (char*)"\\n";
        fprintf(stderr, "->(%d, %s,%d)", tp->target->state, label, tp->freq);
        tp = tp->next_tran;
    }
    fprintf(stderr, ".\n");
}

void printpfsa(NODE *pfsa) {
    NODE *p;

    if (!pfsa)
        return;
    fprintf(stderr, "nstates = %d, narcs = %d, maxstate = %d\n",
            nstates(pfsa), trancnt(pfsa), getmaxstatenum(pfsa));
    for (p = pfsa->nextnode; p; p = p->nextnode) {
        fprintf(stderr, "(%d,%ds,%dt)",
                p->state, p->nsymbols, p->ntrans);
        printlist(p->translist->next_tran);
    }
}

/*
 * getanno() is used to make a character string to annotate an edge
 * for the graphplace program.  This is used because if there is more
 * than one transition to the same target from the same source, graphplace
 * will print them one over the other, thus making all but one illegible.
 * getanno() will order all such transitions as a single continuous string.
 * The way it works is that each time getanno() is called with a transition
 * pointer, it will return a string with annotations of all transitions to
 * that target.  If the current transition is to a target that has already
 * been handled, it will return a NULL string;
 * Because it is prettier to have the arc frequency in superscript, which is
 * only possible if there is a single transition to a state, we build both
 * possibilities (Singleanno requires "(freq) (label)", multianno requires
 * () (label^freq,label^freq,label^freq,...)", but discard the former if
 * there is more than one transition to the same state.
 * Also, to make things tidier, frequency counts of 1 are not displayed.
 */
char *getanno(TRANS *trans, TRANS *tp) {
    extern char *Prog;
    TRANS *t;
    static char annos[512], singleanno[64];
    char buf[256];
    int singletrans = 1;

    /*
     * First go through all the targets in translist before the current
     * one (tp->target) and see if we have already handled this target on
     * non-delimiter symbols.
     */
    for (t = trans; t != tp; t = t->next_tran)
        if (t->target == tp->target && Symtab[t->sym].label[0] != Delim)
            return (char *) 0;

    if (tp->freq == 1)
        sprintf(singleanno, "() (%s)", Symtab[tp->sym].label);
    else
        sprintf(singleanno, "(%d) (%s)", tp->freq, Symtab[tp->sym].label);

    if (tp->freq == 1)
        sprintf(annos, "() (%s", Symtab[tp->sym].label);
    else
        sprintf(annos, "() (%s^%d", Symtab[tp->sym].label, tp->freq);

    while ((t = t->next_tran)) {
        if (t->target == tp->target) {
            if (strlen(annos) > 500) {
                fprintf(stderr, "%s: Anno too long on arc to %d\n", Prog,
                        tp->target->state);
                exit(1);
            }
            singletrans = 0;
            if (t->freq > 1)
                sprintf(buf, ",%s^%d", Symtab[t->sym].label, t->freq);
            else
                sprintf(buf, ",%s", Symtab[t->sym].label);
            strcat(annos, buf);
        }
    }
    if (singletrans)
        return singleanno;
    strcat(annos, ")");
    return annos;
}

void writepfsa(FILE *fp,
        NODE *pfsa) {
    NODE *p;
    TRANS *tp;
    char *label, *s;
    int ndelims; /* no. of delims from a state (for graphplace) */

    if (Graphplace) {
        for (p = pfsa->nextnode; p; p = p->nextnode) {
            ndelims = 0;
            tp = p->translist->next_tran;
            while (tp) {
                label = Symtab[tp->sym].label;
                if (label[0] == '\n')
                    ndelims += tp->freq;
                else if ((s = getanno(p->translist->next_tran, tp)))
                    fprintf(fp, "%s %d %d edge\n", s, p->state, tp->target->state);
                tp = tp->next_tran;
            }
            if (ndelims)
                fprintf(fp, "(!^%d) (%d) () %d node\n", ndelims, p->state, p->state);
            else
                fprintf(fp, "(%d) () %d node\n", p->state, p->state);
        }
    } else {
        for (p = pfsa->nextnode; p; p = p->nextnode) {
            tp = p->translist->next_tran;
            while (tp) {
                label = Symtab[tp->sym].label;
                fprintf(fp, "%d\t%d\t%20s\t%d\n", p->state,
                        tp->target->state,
                        label[0] == '\n' ? "\\n" : label,
                        tp->freq);
                tp = tp->next_tran;
            }
        }
    }
}

/*
 * The following are non-pfsa related general functions
 */

void memerr() {
    Perror((char*)"calloc");
}

void Perror(char *s) /* Fatal perror() */
{
    extern int errno;

    perror(s);
    exit(errno);
}

/*
 * Normalise filename. If fname already has the extension ext, it is left
 * alone. If there is no extension, then the extension ext is tacked on.
 * If there is an extension, then it must be .ext or else mkfname will
 * complain and die.
 */
char *mkfname(char *fname, char *ext)
{
    static char buf[BUFSIZ], *p;

    strcpy(buf, fname);
    if (strcmp(fname, "-")) {
        p = strrchr(fname, '.');
        if (!p || *(p + 1) == '/')
            strcat(buf, ext);
        else if (strcmp(p, ext)) {
            fprintf(stderr, "Unexpected file type %s when expecting %s\n",
                    fname, ext + 1);
            exit(1);
        }
    }
    return buf;
}

/*
 * This converts a buffer of delimited tokens into a string of
 * integers which index into the Symbol table
 */
int *toks2syms(char *buf)
{
    static int symbols[BUFSIZ];
    char *label;
    int sym, i = 0;

    label = strtok(buf, ":");
    while (label && i < BUFSIZ - 1) {
        sym = findsym(label);
        symbols[i++] = sym;
        label = strtok((char *) NULL, ":");
    }
    symbols[i] = 0;
    return symbols;
}

/*
 * This converts an array of indexes into Symtab into a : delimited
 * character tokens
 */
char *syms2toks(int *p)
{
    static char buf[BUFSIZ];
    char *label = (char *) 0;

    buf[0] = '\0';
    while (*p) {
        label = Symtab[*p].label;
        if (label[0] == Delim) {
            strcat(buf, "\\n\n");
            break;
        } else
            strcat(buf, label);
        strcat(buf, ":");
        ++p;
    }
    if (label && label[0] != Delim)
        strcat(buf, "*UNDELIMITED*\n");
    return buf;
}

/*
 * Turns debug on and off during a long execution
 */
void onusr1(int par)
// int par; /* UNUSED */
{
    Debug = Debug ? 0 : 1;
    signal(SIGUSR1, onusr1);
}

void clearmarks(NODE *pfsa)
{
    while (pfsa) {
        pfsa->mark = 0;
        pfsa = pfsa->nextnode;
    }
}

/*
 * Historical reasons:  The program was first implemented using single
 * character symbols, so an input string for induction can be just a
 * character string.  Later versions used tokens for arc labels.  But
 * to retain the rest of the code, the label names were entered in an
 * integer array and the labels were still indexed by a single number.
 * However, it became necessary to create the following integer
 * versions of the str() functions.
 */
int *intcpy(int *p, int *q)
{
    int *r = p;

    while (*q)
        *r++ = *q++;
    *r = 0;
    return p;
}

int *intcat(int p[], int n)
{
    int *r = p;

    while (*r)
        r++;
    *r++ = n;
    *r = 0;
    return p;
}

int intcmp(int *p, int *q) {
    while (*p && *q) {
        if (*p != *q)
            break;
        ++p;
        ++q;
    }
    return *p - *q;
    /*
       return strcmp(Symtab[*p].label, Symtab[*q].label);
     */
}

int intlen(int *p)
{
    int len = 0;

    while (*p++)
        len++;
    return len;
}

int *intdup(int *p)
{
    int len, *s;

    len = intlen(p);
    s = (int*)calloc(len + 1, sizeof (int));
    if (!s)
        return (int *) NULL;
    return intcpy(s, p);
}

/* Remove 0 frequency arcs from the pfsa
 */
NODE * trim(NODE *pfsa) {
    NODE *p, *temp_p;
    TRANS *tp, *temp_tp;

    p = pfsa;
    while (p->nextnode) {
        tp = p->nextnode->translist;
        while (tp->next_tran) {
            if (!tp->next_tran->freq) {
                if (tp->next_tran->sym != DELIMITER)
                    decr_trancnt(pfsa);
                temp_tp = tp->next_tran;
                tp->next_tran = tp->next_tran->next_tran;
                free((void *) temp_tp);
            } else
                tp = tp->next_tran;
        }

        /* Remove node if all trans have gone 
         */
        if (!p->nextnode->translist->next_tran) {
            temp_p = p->nextnode;
            p->nextnode = p->nextnode->nextnode;
            free((void *) temp_p);
            decr_nodecnt(pfsa);
        } else
            p = p->nextnode;
    }
    return pfsa;
}

#endif /* MISC_C */