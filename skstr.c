/**
 * This is the stochastic k-tails algorithm. This is the same as the
 * k-tails algorithm (see k-tails.c) except that only the top s% of the
 * distribution of k-strings (not k-tails) is considered. 
 *
 * This is an ineffecient implementation, but it was done in a hurry and
 * I never got around to redoing it.
 *
 * 16/8/95: &
 * NOTE: low values of k give best results mostly.
 *
 * 29/9/95: Reverted to this version after having tested the version where
 * the top n% of the strings of each node p and q must agree.
 *
 * 01/07/96: This is the long awaited rethink of the implementation of this
 * algorithm.  Includes heuristics to bound the number of strings generated.
 *
 * 22/08/96: Now includes the other three options I discuss in the Raman &
 * Patrick 1995 paper:  skstr_or, skstr_lax and skstr_strict
 *
 * The number of strings from a state that need to be considered is
 * exponential in tail size.  Thus for tail sizes of 2 and above, the
 * space complexity increases rapidly and makes the exact
 * implementation unfeasible.  Ideally, we should build strings in
 * decreasing order of probability and stop when the sum of
 * probabilities of strings generated thus far exceeds Agreepct.  But
 * this is not possible unless all strings have been examined.  Thus
 * we need to ignore some strings at least.  The following heuristics
 * were looked at:
 *
 * Heuristic 1: At each state, choose transitions in decreasing order of
 * probability and stop when the maximum number of strings have been
 * built up from a state. This has the following disadvantage.  Suppose
 * that a state 0 has a transition of prob 0.9 to state 1 and 0.1 to
 * state 2.  And that state 1 has final transition with a prob of 0.99
 * and 10 other final transitions of probability 0.001 each.  Suppose
 * also that state 2 has one transition with prob 0.9 and one with prob
 * 0.1.  Then the string list should include one string with a
 * probability of 0.9 * 0.99, one string with probability 0.1 * 0.9 and
 * one with 0.1 * 0.1.  However this will not be the case since all the
 * 0.9 * 0.001 strings will be added to the string list before the 0.1 *
 * 0.1 strings are.  If string list should fill up before considering
 * state 2, then some undeserving candidates have been included while
 * some deserving candidates have been left out.
 *
 * Heuristic 2: Reject strings less than a certain threshold of
 * probability. Now consider a state with one or more low probability (<
 * 1%) transitions. Let's say that the threshold is 1%.  The heuristic is
 * to ignore these transitions in building strings to be considered from
 * this state.  We must see if this is safe.  Now either this state has
 * other higher probability transitions, or not.  If it has other high
 * probability transitions, then there will either be many high
 * probability transitions or one very high probability transition. In
 * either of these cases, these high-prob transition(s) determine the
 * character of this state and not the low prob ones.  Thus it is safe to
 * ignore the low probability transitions.  However, if there are no high
 * probability transitions from this state, then there are many (at least
 * 100 in this case) low prob transitions.  These strings will all have
 * to be considered in merging this with another similar state.  However,
 * this is an unlikely situation, thus this heuristic can be said to
 * perform reasonably.
 *
 * Heuristic 3: Reject strings that are less than n (say n=10) times
 * as probable as the most probable k-string from that state.  This
 * is a refinement of the previous heuristic in which case low prob
 * strings are rejected, but the rejection is based on the relative
 * probability of a string as compared with its neighbours as opposed
 * to its absolute probability.  The complexity of finding the most
 * probable string from a state is exponential again in the length
 * of the string, since all possible strings need to be considered.
 * It won't do to just follow the maximum freq trans at each state
 * because of the possibility of bad fragmentation down the line.  Eg.
 * on (a,40) and (b,60) we may choose b, but a may proceed unfragmented
 * to a delimiter and thus derive a string with prob .4, but b may fragment
 * into 10 paths of prob .1 each, in which case, our initial guess would
 * only yield a string of prob .06 at best!  This was tried and found
 * to be empirically unsuitable too.
 */
#ifndef SKSTR_C
#define SKSTR_C
#include "pfsa.h"
#include <math.h>

#define TAILSIZE 1
#define AGREEPCT 50      /* What % of strings must agree before merging */
#define PREC 1000        /* How many decimal places (/10) in percentages */
#define MINPROB (1*PREC) /* reject strings less than 1% probable */
#define MINENTROPY 0.5

/*
 * The MAXSTR constant determines how many generated strings can be held
 * in the comparison buffer.  The default value is 1K. If MINPROB is set
 * to (100/MAXSTR)%,it ensures that the maximum number of strings generated
 * from a state is MAXSTR.
 */
#define MAXSTR 1000

typedef struct {
    int *kstr;
    u_long prob;
} kstring;

struct kstrList {
    kstring *ks;
    int nstr;
};


/*
 * Externals
 */
extern char *Prog, Outfile[], Infile[], Callstring[];
extern int Tailsize; /* Defined in ktail.c */

/*
 * Globals:
 *
 * Agreepct is the percentage of the distribution of strings between
 * two states that must agree if they are to be merged.  The default is
 * 100% or all strings generated from state 1 must agree with state 2
 * and vice versa.  For the actual meaning of the word "agree" in this
 * context, see the Raman&Patrick paper on the sk-strings algorithm.
 * In this implementation, all conditions are used and are selectable
 * using the -H switch (25/9/96), although the default is AND (state 2
 * must accept all of the top s% state 1's strings AND vice versa)
 *
 * Strings lower than Minprob from a state will be ignored (not added
 * to the string list for consideration)  This helps because the no.
 * of strings is exponential with tailsize.
 *
 * Ksv_cache is the cache of strings generated from a given state.
 * We need this since we may repeatedly require to access strings
 * from a state.  When a merge is done anywhere in the pfsa, the
 * cache needs to be flushed.
 */
int Agreepct = AGREEPCT;
u_long Minprob = MINPROB;
double MinEntropy = -1;
struct kstrList **Ksv_cache = (struct kstrList **) NULL;
int Cache_size = 0;
char Heuristic[128] = "AND";

/*
 * local function prototypes
 */
static int (*Sk_compare)(const void *, const void *);
static int sk_compare_byProb(const void *, const void *);
static int sk_compare_byStr(const void *, const void *);

static NODE *do_skstrings(NODE *);
int sk_distinguishable(NODE *p1, NODE *p2);
int acceptlist(struct kstrList *, NODE *);
void addstring(int [], u_long, struct kstrList *);
void get_kstrList(int, NODE *, int [], u_long, struct kstrList *);
void dispose_strs(struct kstrList *);
void flush_cache(void);
struct kstrList *get_sorted_kstrList(int k, NODE *p, int syms[], u_long prob);
static void usage_skstr(char *);

/* sk-string algorithms */
static int (*Sk_mergeable)(int, NODE *, NODE *);
static int skstr_or(int k, NODE *p, NODE *q);
static int skstr_and(int k, NODE *p, NODE *q);
static int skstr_lax(int k, NODE *p, NODE *q);
static int skstr_strict(int k, NODE *p, NODE *q);
static int skstr_xentropic(int k, NODE *p, NODE *q);
static int skstr_vardist(int k, NODE *p, NODE *q);
static void onusr2(int);

NODE *skstr(int argc,
        char **argv) {
    int c;

    setbuf(stderr, (char *) NULL);
    Tailsize = TAILSIZE;
    while ((c = getopt(argc, argv, "dvgD:o:m:t:p:e:hH:")) != EOF) {
        switch (c) {
            case 'H':
                strcpy(Heuristic, optarg);
                break;
            case 'D':
                Delim = optarg[0];
                break;
            case 'd':
                ++Debug;
                break;
            case 'v':
                ++Verbose;
                break;
            case 'g':
                ++Graphplace;
                break;
            case 'o':
                strcpy(Outfile, optarg);
                break;
            case 'm':
                Minprob = (u_long) (atof(optarg) * PREC);
                if (Minprob <= 0 || Minprob > 100 * PREC) {
                    fprintf(stderr, "Illegal -m optarg reset to %d%%\n", MINPROB / PREC);
                    Minprob = MINPROB;
                }
                break;
            case 'e':
                MinEntropy = (double) atof(optarg);
                if (MinEntropy < 0 || MinEntropy > 1) {
                    fprintf(stderr, "Illegal -e optarg reset to %.2f\n", MINENTROPY);
                    MinEntropy = MINENTROPY;
                }
                break;
            case 't':
                Tailsize = atoi(optarg);
                if (Tailsize < 0) {
                    fprintf(stderr, "Illegal -t optarg reset to %d\n", TAILSIZE);
                    Tailsize = TAILSIZE;
                }
                break;
            case 'p':
                Agreepct = atoi(optarg);
                if (Agreepct < 0 || Agreepct > 100) {
                    fprintf(stderr, "Invalid Sk-tails probability value (%d)"
                            "reset to %d%%.\n", Agreepct, AGREEPCT);
                    Agreepct = AGREEPCT;
                }
                break;
            case 'h':
            default:
                usage_skstr(Prog);
                exit(1);
                break;
        }
    }
    if (argc > optind)
        setfilenames(argv[optind]);

    Sk_compare = sk_compare_byProb;
    if (!strcasecmp(Heuristic, "or")) {
        sprintf(Callstring, "%s -H %s %s%s-t %d -p %d -m %.2f -o %s %s", Prog,
                Heuristic, Verbose ? "-v " : "", Debug ? "-d " : "", Tailsize,
                Agreepct, ((double) Minprob) / PREC, Outfile, Infile);
        Sk_mergeable = skstr_or;
    } else if (!strcasecmp(Heuristic, "and")) {
        Sk_mergeable = skstr_and;
        sprintf(Callstring, "%s -H %s %s%s-t %d -p %d -m %.2f -o %s %s", Prog,
                Heuristic, Verbose ? "-v " : "", Debug ? "-d " : "", Tailsize,
                Agreepct, ((double) Minprob) / PREC, Outfile, Infile);
    } else if (!strcasecmp(Heuristic, "lax")) {
        Sk_mergeable = skstr_lax;
        sprintf(Callstring, "%s -H %s %s%s-t %d -p %d -m %.2f -o %s %s", Prog,
                Heuristic, Verbose ? "-v " : "", Debug ? "-d " : "", Tailsize,
                Agreepct, ((double) Minprob) / PREC, Outfile, Infile);
    } else if (!strcasecmp(Heuristic, "strict")) {
        Sk_mergeable = skstr_strict;
        sprintf(Callstring, "%s -H %s %s%s-t %d -p %d -m %.2f -o %s %s", Prog,
                Heuristic, Verbose ? "-v " : "", Debug ? "-d " : "", Tailsize,
                Agreepct, ((double) Minprob) / PREC, Outfile, Infile);
    } else if (!strcasecmp(Heuristic, "xentropic")) {
        Agreepct = 100;
        Sk_mergeable = skstr_xentropic;
        Sk_compare = sk_compare_byStr;
        if (MinEntropy < 0)
            MinEntropy = MINENTROPY;
        sprintf(Callstring, "%s -H %s %s%s-t %d -e %.2f -m %.2f -o %s %s", Prog,
                Heuristic, Verbose ? "-v " : "", Debug ? "-d " : "", Tailsize,
                (double) MinEntropy, ((double) Minprob) / PREC, Outfile, Infile);
    } else if (!strcasecmp(Heuristic, "vardist")) {
        Agreepct = 100;
        Sk_mergeable = skstr_vardist;
        Sk_compare = sk_compare_byStr;
        if (MinEntropy < 0)
            MinEntropy = MINENTROPY;
        sprintf(Callstring, "%s -H %s %s%s-t %d -e %.2f -m %.2f -o %s %s", Prog,
                Heuristic, Verbose ? "-v " : "", Debug ? "-d " : "", Tailsize,
                (double) MinEntropy, ((double) Minprob) / PREC, Outfile, Infile);
    } else {
        usage_skstr(Prog);
        exit(1);
    }

    signal(SIGUSR2, onusr2);
    buildpfsa(Infile);
    Ksv_cache = (struct kstrList **) calloc(nstates(Pfsa),
            sizeof (struct kstrList *));
    if (!Ksv_cache)
        memerr();
    Cache_size = nstates(Pfsa);
    return do_skstrings(Pfsa);
    signal(SIGUSR2, SIG_DFL);
}

/*
 * Two states are indistinguishable if they produce exactly the same strings
 * with exactly the same probabilities.  We need this to make the skstrings
 * function more efficient.  If two indistinguishable states are merged
 * this will not affect the output strings of up-stream states.  Thus there
 * is no need to flush the cache of strings or to restart comparisons from
 * state 0.
 * Note that the get_sorted_kstrList() call is of O(1) complexity in this
 * function as it should just be a cache access.
 */
int sk_distinguishable(NODE *p1, NODE *p2) {
    int i, syms[128];
    struct kstrList *ksv1, *ksv2;

    syms[0] = 0;
    ksv1 = get_sorted_kstrList(Tailsize, p1, syms, 100 * PREC);
    syms[0] = 0;
    ksv2 = get_sorted_kstrList(Tailsize, p2, syms, 100 * PREC);

    if (ksv1->nstr != ksv2->nstr)
        return 1;
    for (i = 0; i < ksv1->nstr; i++) {
        if (intcmp(ksv1->ks[i].kstr, ksv2->ks[i].kstr) ||
                ksv1->ks[i].prob != ksv2->ks[i].prob)
            return 1;
    }
    return 0;
}

static NODE *do_skstrings(NODE *pfsa) {
    NODE *p1, *p2;
    int restart = 0;

    for (p1 = pfsa->nextnode; p1; p1 = restart ? pfsa->nextnode : p1->nextnode) {
        restart = 0;
        for (p2 = p1->nextnode; p2; p2 = p2->nextnode) {
            if (Debug)
                fprintf(stderr, "%d-equiv(%04d,%04d)?%s",
                    Tailsize, p1->state, p2->state, isatty(2) ? "\r" : "\n");
            if ((*Sk_mergeable)(Tailsize, p1, p2)) {
                if (Debug)
                    fprintf(stderr, "\nMerging %d & %d\n\n", p1->state, p2->state);
                if (sk_distinguishable(p1, p2)) {
                    merge(pfsa, p1, p2);
                    flush_cache();
                    restart = 1;
                    break;
                } else {
                    dispose_strs(Ksv_cache[p2->state]);
                    free((void *) Ksv_cache[p2->state]);
                    Ksv_cache[p2->state] = (struct kstrList *) NULL;
                    merge(pfsa, p1, p2);
                }
            }
        }
    }
    pfsa = renumber(pfsa);
    return pfsa;
}

static int sk_compare_byProb(const void *p, const void *q) /* Str is the secondary key */ {
    kstring *p1, *q1;

    p1 = (kstring *) p;
    q1 = (kstring *) q;

    if (p1->prob < q1->prob) /* Reverse sort by probability */
        return 1;
    else if (p1->prob > q1->prob)
        return -1;

    return intcmp(p1->kstr, q1->kstr);
}

static int sk_compare_byStr(const void *p, const void *q) /* Prob is the secondary key */ {
    kstring *p1, *q1;
    int diff;

    p1 = (kstring *) p;
    q1 = (kstring *) q;

    diff = intcmp(p1->kstr, q1->kstr);
    if (!diff)
        return (p1->prob == q1->prob ? 0 : (p1->prob < q1->prob ? 1 : -1));
    else
        return diff;
}

int acceptlist(struct kstrList *ksv, NODE *p) /* Is the top Agreepct% of ksv acceptable at p? */ {
    u_long cutoff;
    int i;

    cutoff = 0;
    for (i = 0; i < ksv->nstr; i++) {
        cutoff += ksv->ks[i].prob;
        if (!acceptable(p, ksv->ks[i].kstr))
            return 0;
        if (cutoff > (Agreepct * PREC))
            break;
    }
    return 1;
}

void addstring(int s[],
        u_long prob,
        struct kstrList *ksv) {
    int *r;

    if (ksv->nstr > 0 && !intcmp(ksv->ks[ksv->nstr - 1].kstr, s)) {
        ksv->ks[ksv->nstr - 1].prob += prob;
        return;
    }
    r = intdup(s);
    if (!r)
        memerr();
    ksv->ks[ksv->nstr].kstr = r;
    ksv->ks[ksv->nstr].prob = prob;
    ksv->nstr++;
    if (ksv->nstr == MAXSTR) {
        fprintf(stderr, "Fatal Error:\n");
        fprintf(stderr, "Too many strings are being considered (%d)\n", MAXSTR);

        fprintf(stderr,
                "Consider limiting this by increasing the value of the Minprob parameter\n"
                "which determines the least probability a string should have, to be\n"
                "considered.  A setting of -m %.2f ensures that a maximum of %d strings\n"
                "will be considered\n\n", 100.0 / MAXSTR, MAXSTR);

        fprintf(stderr,
                "Alternatively, you may try recompiling this program with a larger value\n"
                "of MAXSTR (currently %d) to allow more strings to be considered, but be\n"
                "aware that the number of strings is exponential in tailsize and must be\n"
                "limited with the MINPROB parameter (currently %.2f)\n",
                MAXSTR, ((double) Minprob) / PREC);

        exit(1);
    }
    return;
}

void get_kstrList(int k,
        NODE *p,
        int s[],
        u_long prob,
        struct kstrList *ksv) {
    u_long newprob;
    int syms[128];
    TRANS *tp;

    if (k == 0) {
        if (s[0])
            addstring(s, prob, ksv);
        return;
    }

    /*
     * Refer to heuristic 2 above.  We reject paths that are less than
     * 1% probable.  This gives us a maximum of 100 paths to handle
     * from any one state.  Tweak this number using the Minprob parameter
     * above.
     */
    for (tp = p->translist->next_tran; tp; tp = tp->next_tran) {
        intcpy(syms, s);
        intcat(syms, tp->sym);
        newprob = prob * tp->freq / p->ntrans;
        if (newprob < Minprob)
            return;
        if (Symtab[tp->sym].label[0] == Delim)
            addstring(syms, newprob, ksv);
        else
            get_kstrList(k - 1, tp->target, syms, newprob, ksv);
    }
    return;
}

void printstrings(struct kstrList *ksv) {
    int i = 0;
    u_long cutoff = 0;

    for (i = 0; i < ksv->nstr; i++) {
        fprintf(stderr, "%3ld.%-3ld:%s %s", ksv->ks[i].prob / PREC,
                ksv->ks[i].prob % PREC, cutoff > Agreepct * PREC ? "-" : " ",
                syms2toks(ksv->ks[i].kstr));
        cutoff += ksv->ks[i].prob;
    }
}

struct kstrList *get_sorted_kstrList(int k, NODE *p, int syms[], u_long prob) {
    struct kstrList *ksv;

    if (!Ksv_cache) {
        Ksv_cache = (struct kstrList **) calloc(MAXNODES, sizeof (struct kstrList *));
        if (!Ksv_cache)
            memerr();
        Cache_size = MAXNODES;
    }

    if (Ksv_cache && Ksv_cache[p->state])
        return Ksv_cache[p->state];

    ksv = (struct kstrList *) calloc(1, sizeof (struct kstrList));
    if (!ksv)
        memerr();
    ksv->ks = (kstring *) calloc((int) (100L * PREC / Minprob),
            sizeof (kstring));
    if (!ksv->ks)
        memerr();
    ksv->nstr = 0;
    get_kstrList(k, p, syms, prob, ksv);
    qsort((void *) ksv->ks, ksv->nstr, sizeof (kstring), Sk_compare);
    Ksv_cache[p->state] = ksv;
    if (Debug > 1) {
        fprintf(stderr, "Strings from state %d\n", p->state);
        printstrings(ksv);
    }
    return ksv;
}

void dispose_strs(struct kstrList *ksv) {
    while (ksv->nstr > 0) {
        --ksv->nstr;
        free((void *) ksv->ks[ksv->nstr].kstr);
    }
    free((void *) ksv->ks);
    ksv->ks = (kstring *) NULL;
}

void flush_cache(void) {
    int i;

    for (i = 0; i < Cache_size; i++) {
        if (!Ksv_cache[i])
            continue;
        dispose_strs(Ksv_cache[i]);
        free(Ksv_cache[i]);
        Ksv_cache[i] = (struct kstrList *) NULL;
    }
}

/*
 * This handles the case when p's strings must be acceptable at q AND
 * q's strings must be acceptable at p.
 */
static int skstr_and(int k, NODE *p, NODE *q) {
    int syms[128];
    struct kstrList *ksv_p, *ksv_q;

    /*
     * Get all the k-strings that can be generated from p with their probs.
     * Sort the strings in descending order of probability.
     * Check if the top Agreepct% of the list is acceptable at q.
     */
    syms[0] = 0;
    ksv_p = get_sorted_kstrList(k, p, syms, 100 * PREC);
    if (!acceptlist(ksv_p, q))
        return 0;

    /*
     * Now the other way
     */
    syms[0] = 0;
    ksv_q = get_sorted_kstrList(k, q, syms, 100 * PREC);

    if (!acceptlist(ksv_q, p))
        return 0;

    /*
     * At this point we know the top s% of p's k-strings are acceptable
     * by q and vice versa.  So the two states must be sk-equiv.
     */
    return 1;
}

/*
 * This handles the case when p's strings must be acceptable at q OR
 * q's strings must be acceptable at p.
 */
static int skstr_or(int k, NODE *p, NODE *q) {
    int syms[128];
    struct kstrList *ksv_p, *ksv_q;

    /*
     * Get all the k-strings that can be generated from p with their probs.
     * Sort the strings in descending order of probability.
     * Check if the top Agreepct% of the list is acceptable at q.
     */
    syms[0] = 0;
    ksv_p = get_sorted_kstrList(k, p, syms, 100 * PREC);
    if (acceptlist(ksv_p, q))
        return 1;

    /*
     * Now the other way
     */
    syms[0] = 0;
    ksv_q = get_sorted_kstrList(k, q, syms, 100 * PREC);

    if (acceptlist(ksv_q, p))
        return 1;

    /*
     * At this point we know that the top s% of neither state's
     * strings are acceptable at the other state. So the 2 states are
     * not sk-mergeable under this criterion.
     */
    return 0;
}

/*
 * This handles the case where the two states must agree in the
 * first Agreepct% of their strings, but not necessarily with the
 * same probabilities.
 */
static int skstr_lax(int k, NODE *p, NODE *q) {
    int i, syms[128];
    struct kstrList *ksv_p, *ksv_q;
    u_long cutoffp = 0, cutoffq = 0;

    syms[0] = 0;
    ksv_p = get_sorted_kstrList(k, p, syms, 100 * PREC);
    syms[0] = 0;
    ksv_q = get_sorted_kstrList(k, q, syms, 100 * PREC);
    for (i = 0; i < ksv_p->nstr && i < ksv_q->nstr; i++) {
        if (intcmp(ksv_p->ks[i].kstr, ksv_q->ks[i].kstr))
            return 0;
        cutoffp += ksv_p->ks[i].prob;
        cutoffq += ksv_q->ks[i].prob;
        if (cutoffp >= Agreepct && cutoffq >= Agreepct)
            return 1;
    }
    return 0;
}

/*
 * This handles the case where the two states must agree in the
 * first Agreepct% of their strings, AND with the same probabilities.
 */
static int skstr_strict(int k, NODE *p, NODE *q) {
    int i, syms[128];
    struct kstrList *ksv_p, *ksv_q;
    u_long cutoffp = 0, cutoffq = 0;

    syms[0] = 0;
    ksv_p = get_sorted_kstrList(k, p, syms, 100 * PREC);
    syms[0] = 0;
    ksv_q = get_sorted_kstrList(k, q, syms, 100 * PREC);
    for (i = 0; i < ksv_p->nstr && i < ksv_q->nstr; i++) {
        if (ksv_p->ks[i].prob != ksv_q->ks[i].prob ||
                intcmp(ksv_p->ks[i].kstr, ksv_q->ks[i].kstr))
            return 0;
        cutoffp += ksv_p->ks[i].prob;
        cutoffq += ksv_q->ks[i].prob;
        if (cutoffp >= Agreepct && cutoffq >= Agreepct)
            return 1;
    }
    return 0;
}

/*
 * This handles the case when the cross-entropic coding of strings
 * (Kullback)in each state must be less than a certain threshold.
 * i.e.  Sigma { (pi*log(qi)+qi*log(pi) - (pi*log(pi)+qi*log(qi))) }
 * where pi and qi are probabilities of the i'th string at p and q resp.
 * For computational purposes, where pi or qi is zero, we assume
 * the value MINPROB which is the least probability a string can have.
 * If the two states are identical, this measure will be zero.
 * It gives the expected number of bits to encode the difference
 * between them.
 */
static int skstr_xentropic(int k, NODE *p, NODE *q) {
    int i, j, diff, syms[128];
    struct kstrList *ksv_p, *ksv_q;
    double xentropy = 0, epsilon, pi, qi;

    syms[0] = 0;
    ksv_p = get_sorted_kstrList(k, p, syms, 100 * PREC);
    syms[0] = 0;
    ksv_q = get_sorted_kstrList(k, q, syms, 100 * PREC);
    i = j = 0;
    epsilon = (double) Minprob / 100.0 / PREC;
    while (i < ksv_p->nstr && j < ksv_q->nstr) {
        diff = intcmp(ksv_p->ks[i].kstr, ksv_q->ks[j].kstr);
        if (!diff) {
            pi = (double) ksv_p->ks[i++].prob / 100 / PREC;
            qi = (double) ksv_q->ks[j++].prob / 100 / PREC;
        } else if (diff < 0) {
            pi = (double) ksv_p->ks[i++].prob / 100 / PREC;
            qi = epsilon;
        } else {
            pi = epsilon;
            qi = (double) ksv_q->ks[j++].prob / 100 / PREC;
        }
        xentropy += (pi - qi) * log(pi / qi);
    }
    while (i < ksv_p->nstr) {
        pi = (double) ksv_p->ks[i++].prob / 100 / PREC;
        qi = epsilon;
        xentropy += (pi - qi) * log(pi / qi);
    }
    while (j < ksv_q->nstr) {
        pi = epsilon;
        qi = (double) ksv_q->ks[j++].prob / 100 / PREC;
        xentropy += (pi - qi) * log(pi / qi);
    }
    xentropy /= -2.0 * (1.0 - epsilon) * log(epsilon);

    if (Debug > 1)
        fprintf(stderr, "Nodes %d and %d: xentropy = %.3f\n",
            p->state, q->state, xentropy);
    return (xentropy <= MinEntropy);
}

/*
 * This handles the case when the sum of the absolute differences of
 * string probabilities at each state must be less than a certain threshold.
 * i.e.  Sigma { Abs(pi-qi) } < e
 * where pi and qi are probabilities of the i'th string at p and q resp.
 */
static int skstr_vardist(int k, NODE *p, NODE *q) {
    int i, j, diff, syms[128];
    struct kstrList *ksv_p, *ksv_q;
    double vardist = 0, pi, qi;

    syms[0] = 0;
    ksv_p = get_sorted_kstrList(k, p, syms, 100 * PREC);
    syms[0] = 0;
    ksv_q = get_sorted_kstrList(k, q, syms, 100 * PREC);
    i = j = 0;
    while (i < ksv_p->nstr && j < ksv_q->nstr) {
        diff = intcmp(ksv_p->ks[i].kstr, ksv_q->ks[j].kstr);
        if (!diff) {
            pi = (double) ksv_p->ks[i++].prob / 100 / PREC;
            qi = (double) ksv_q->ks[j++].prob / 100 / PREC;
        } else if (diff < 0) {
            pi = (double) ksv_p->ks[i++].prob / 100 / PREC;
            qi = 0;
        } else {
            pi = 0;
            qi = (double) ksv_q->ks[j++].prob / 100 / PREC;
        }
        vardist += fabs(pi - qi);
    }
    while (i < ksv_p->nstr)
        vardist += (double) ksv_p->ks[i++].prob / 100 / PREC;
    while (j < ksv_q->nstr)
        vardist += (double) ksv_q->ks[j++].prob / 100 / PREC;
    vardist /= 2.0;

    if (Debug > 1)
        fprintf(stderr, "Nodes %d and %d: vardist = %.3f\n",
            p->state, q->state, vardist);
    return (vardist <= MinEntropy);
}

static void usage_skstr(char *prog) {
    char *usagestring = (char*)
            "This program optimises the given minimal canonical pfsa by successively\n"
            "merging pairs of states it deems equivalent.  Six types of equivalence\n"
            "relations can be invoked by calling this program with the following\n"
            "options:\n"
            "\n"
            "Progname      Heuristic\n"
            "--------      ---------\n"
            "-H and        (Default) P and Q are mergeable if both accept the top\n"
            "              Agreepct% of the distribution of strings of the other.\n"
            "\n"
            "-H or         P and Q are mergeable if either accepts the top Agreepct%\n"
            "              of the distribution of strings of the other.\n"
            "\n"
            "-H lax        P and Q are mergeable if the top Agreepct% of their\n"
            "              distribution of strings is the same, in the same order.\n"
            "\n"
            "-H strict     P and Q are mergeable if the top Agreepct% of their\n"
            "              distribution of strings is the same, in the same order,\n"
            "              with the same probabilities.\n"
            "\n"
            "-H xentropic  P and Q are mergeable if the cross-entropic coding of\n"
            "              their strings is less than 'e'. The x-entropic coding is:\n"
            "              Sigma (pi*log(qi/pi) + qi*log(pi/qi)), where pi and qi are\n"
            "              probabilities of the i'th string at p and q respectively.\n"
            "\n"
            "-H vardist    P and Q are mergeable if the variation distance between\n"
            "              their string probability distributions is less than 'e'\n"
            "              That is Sigma (Abs(pi - qi)) where pi and qi are as above.\n"
            "              This is basically a hack version of the xentropic scheme.\n"
            "\n"
            "              xentropy and vardist are normalised to be in the range\n"
            "              0..1 before comparison with e.  Thus 0 <= e <= 1.\n"
            "\n"
            "If the strings are in the file f1.pfsa, the output is written to the\n"
            "file f1.opfsa.\n"
            "\n"
            "PFSA are specified as a list of source target symbol frequency' values\n"
            "\n"
            "Options: (Defaults shown in square brackets)\n"
            "Options unnecessary for the selected algorithm will be ignored.\n\n"
            "\n"
            "-p num    Set Agreepct% to num [50%]\n"
            "-t num    Consider output strings of size <= num for every state [1]\n"
            "-m num    String must be at least num% probable to be considered [1%]\n"
            "-e num    Set minimum entropy [0.5] or minimum vardist [0.5], see above\n"
            "-d        Debug mode: prints miscellaneous info while executing [0]\n"
            "-v        Verbose mode: prints extra information in result [0]\n"
            "-D char   Set delimiter to 'char' [\\n]\n"
            "-g        Output PFSA in Graphplace format [0]\n"
            "-o file   Output PFSA in `file' [`infile.opfsa' or stdout]\n"
            "\n"
            "Minprob (set using -m) determines the least probability a string must\n"
            "have in order to be considered in comparison with strings from another\n"
            "state.  Use it to place an upper bound on the number of strings that\n"
            "are generated from a state.  Otherwise, the number of strings is\n"
            "exponential in tailsize.  The default value of Minprob is 1%.  For\n"
            "practical reasons, Minprob canot be 0. Since the precision of operation\n"
            "is 3 decimal places, setting Minprob to less than .001% is the same as\n"
            "setting it to zero, causing it to be reset to 1%.\n";
    fprintf(stderr, "usage: skstr [options] [input file]\n");
    fprintf(stderr, "%s", usagestring);
    return;
}

static void onusr2(int par) {
    fprintf(stderr, "No progress report available.  Try USR1\n");
}
#endif /*#ifndef SKSTR_C*/