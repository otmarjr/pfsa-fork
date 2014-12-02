/* 
 * File:   main.cpp
 * Author: otmar
 *
 * Created on 1 de Dezembro de 2014, 22:25
 */

#include <iostream>
#include "pfsa.h"
#include "misc.c"
#include "skstr.c"


/*
 * 
 */
/* Globals */
char *Prog;
char Outfile[BUFSIZ]="-", Infile[BUFSIZ]="-";
char Callstring[128] = "opt";

int main(int argc, char** argv) {

    
    
    return 0;
}

void setfilenames(char *arg)
{
   strcpy(Infile, mkfname(arg, (char *)".pfsa"));
   if (strchr(Infile, '.') && !strcmp(Outfile, "-")) {
      strcpy(Outfile, Infile);
      strcpy(strrchr(Outfile, '.'), ".opfsa");
   }
   return;
}   

/*
 * output_pfsa() is used to output the PFSA in the same format we expect
 * to read them in, except if the Graphplace flag is given.
 */
void output_pfsa(NODE *pfsa, char outfile[])
{
   FILE *fp;
   char buf[128];

   if (!strcmp(outfile, "-"))
      fp = stdout;
   else
      fp = fopen(outfile, "w");
   if (!fp)
      Perror(outfile);
   if (Graphplace) {
      fprintf (fp, "%%!PS\n/makearrows true def\n");
      if (Verbose)
	 fprintf (fp, "/Times-Roman findfont 18 scalefont setfont\n"
		  "300 810 (%s) centertext\n", Callstring);
      writepfsa(fp, pfsa);
      if (Verbose) {
	 sprintf(buf, "MML = %.2f bits", mml(Pfsa, (double *) 0));
	 fprintf (fp, "/Times-Roman findfont 12 scalefont setfont\n"
		  "300 36 (%s) centertext\n", buf);
      }
   } else {
      if (Verbose)
	 fprintf(fp, "# %s\n", Callstring);
      writepfsa(fp, pfsa);
      if (Verbose) {
	 fprintf(fp, "# nstates = %d, narcs = %d, maxstate = %d\n",
		 nodecnt(pfsa), trancnt(pfsa), getmaxstatenum(pfsa));
	 fprintf(fp, "# MML = %.2f bits\n", mml(pfsa, (double *) 0));
      }
   }
   if (fp != stdout)
      (void) fclose(fp);
}

static void usage(char *prog)
{
   char *usagestring = (char *)
      "To perform an algorithm on the canonical machine, program name should\n"
      "be one of the following:\n\n" 
      "\t beams:  Do a breadth first beam search\n"
      "\t simba:  Do a breadth first simba search\n"
      "\t ktail:  Do Biermann & Feldman's (1979) k-tails algorithm\n"
      "\t skstr:  Do Raman & Patrick's (1995) sk-strings algorithm\n\n"
      "For further information on each of the algorithm's options, invoke\n"
      "the appropriate program with the -h option\n\n";
   fprintf(stderr, "This program was called with the name: %s\n", prog);
   fprintf(stderr, "%s", usagestring);
}
