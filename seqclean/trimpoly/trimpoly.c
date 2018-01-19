/* takes EST sequences from stdin only
*/
  
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "trimpoly.h"
#include "color_c.h"
#include "ctype.h"

#define USAGE  "Usage: \n\
   trimpoly [-h] [-n <percN>] [-y <Nperc>] [-m <percdist>] [-PN] [-e <end>] [-r <minrun>] \n\
    [-s <minscore>] [-u <maxgap>] [-l <percdist>|-L <range>] [-X][-D] \n\
    \n\
   Trims polyA/polyT and Ns from the end of the FASTA DNA sequences read\n\
   from standard input.\n\
   Outputs the following columns, tab-delimited:\n\
          seq_name, perc_N, end5, end3, initial_length \n\
   \n\
   -n try the N-trimming if percent of N is above <percN>\n\
      this option MUST be present if you want to enable N-trimming!\n\
      Unless -A or -N switches are present, N trimming is done before \n\
      polyA/T trimming\n\
   -y try to extract from a sequence the maximum range having \n\
      percent of N lower than <Nperc> (default: 3.00)\n\
   -m <percdist> - incremental distance from ends, in percentage of\n\
       sequence length, where N-trimming is done (default:12) \n\
      The resulting distance is auto-limited to 30nt. max\n\
   -N if -n flag was specified, do only N-trimming (skip the \n\
      polyA/T trimming)\n\
   -P do the polyA/T trimming first and then the N trimming;\n\
      only makes sense if -n option is present.\n\
 Poly A/T trimming:\n\
   -e <end> restrict the poly trimming to one end only or to the \n\
            maximum scoring poly match\n\
        (default is: trim both ends)\n\
        if <end> is 3 - trim only polyA from end3\n\
        if <end> is 5 - trim only polyT from end5\n\
        if <end> is x - only the maximum scoring of end3/end5 will be trimmed\n\
       (if -X option is used, polyT/polyA will be swapped)\n\
   -s <minscore> minimum score for poly A/T trimming (default 7)\n\
      it is dynamically increased by 1 for every 10nt distance from the ends\n\
      (scoring: +1 for each A/T, -2 for each intervening nucleotide) \n\
   -r <minrun> minimum count of consecutive A/T (default 5)\n\
      to start extending and computing score from\n\
   -u <extstop> the count of non-matching nucleotides\n\
      that will determine a stop in extending a match (default 4)\n\
   -l incremental distance from ends, in percentage of length,\n\
      where polyA/T trimming is allowed (default:15)\n\
   -L fixed distance (in nucleotides) from both ends where\n\
      polyA/T trimming is allowed (overrides -l)\n\
   -o <file> outputs a multi FASTA file containing only the trimmed\n\
      sequences \n\
   -R relative shifting for end3, end5 will be displayed \n\
      in the tabbed output instead of absolute coordinates\n\
   -H parse starting end5, end3 coordinates from the FASTA file \n\
      (expect to find them after sequence name, tab delimited)\n\
   -X inverted trimming: trim polyA from end5 and polyT from end3\n\
   -Y try inverted trimming if no poly(A/T) was trimmed in the first pass\n\
   -D debug display mode: sequences are shown trimmed, colored\n\
      using console codes (trimmed parts have a red background)\n"      

#define MIN_PERC_LEN 30
/*
 trimpoly -n2 -m5 
  seems to be a good way to trim ESTs with possible 'dirty' ends
*/

typedef struct feature_t {
 char* seq_name;
 char* defline;
 char* seq; /* sequence, upper case, spaces, tabs, newlines are removed */
 int seqlen;
 double perc_N;
 int end5;
 int end3;
 int* NPos; /* position of Ns in the sequence  - dynamically alocated array */
 int NCount;
} feature_t;

struct feature_t feat; /* current sequence data being analysed */

/* global options/flags: */
int skip_poly=0;  /* -N */
int poly_first=0; /* -P */
int perc_lenN=12; /* -m */
double perc_N=0; /* -n */
double nsave_limit=3.00; /* -y */
int trim_end=0;   /* -e  */
int min_score=7;  /* -s */
int min_run=5;    /* -r  */
int max_gap=4;    /* -u */
int perc_len=15;  /* -l */
int swap_poly=0;  /* -X, -Y */
int c_show=0;     /* -D */
char outfile[256]; /* -o */
int do_write=0;   /* -o */
int relative_shift=0;  /* -R */
int init_ends=0; /* -H */
int fixed_len=0;
char* seed_5=SEED_T;
char* seed_3=SEED_A;
char  chr_5='T';
char  chr_3='A';
int first_pass=1;

void optErr(char opt);
void showUsage();
char *seq_upper(char *dest, char *src, int *seqsize);
char* lstrstr(char* lstart, char *rend, char* substr);
char* rstrstr(char* rstart, char *lend, char* substr);

/* these funcs are acting directly on feat structure */
int analyse3();
int analyse5();
void N_trim();
void N_analyse(int l5, int l3, int p5, int p3);
void N_save(int l5, int l3);
double N_calc(int e5, int e3, int* percN);

int main(int argc, char * const argv[]) {
  int c,wpos, init_end5, init_end3, score1, score2, save5, save3;
  int alloc_size;
  int seq_size=0, readlen;
  char *p, *p_end, *bufpos;
  char seq_name[256];
  char def_line[512];
  char line[121];
  FILE * fout;
  char inbuf[BUF_LEN]; /* incoming buffer for sequence lines. */
  char * seq; /* working buffer */  
  alloc_size=BUF_LEN;
  seq=malloc(BUF_LEN); /* initial allocation for sequence buffer */
  while ((c=getopt(argc, argv, "hNPXYDRHm:n:e:r:s:l:L:o:u:y:")) > 0) {
     switch (c)
      { case 'N': skip_poly=1;break;
        case 'P': poly_first=1;break;      
        case 'm': perc_lenN=strtol(optarg, (char **)NULL, 10);
                  if (perc_lenN==0) optErr('m');
                  break;      
        case 'n': perc_N=strtod(optarg, (char **)NULL);
                  if (perc_N==0) optErr('n');
                  break;      
        case 'y': nsave_limit=strtod(optarg, (char **)NULL);
                  if (nsave_limit==0) optErr('y');
                  break;      
        case 'e': if (optarg[0]=='x') {
                         trim_end=1;
                         break;
                         }
                  trim_end=strtol(optarg, (char **)NULL, 10);       
                  if ((trim_end!=3 && trim_end!=5) || trim_end==0) 
                         optErr('e');
                  break;
        case 'r': min_run=strtol(optarg, (char **)NULL, 10);
                  if (min_run<=1) optErr('r');
                  break;
        case 'u': max_gap=strtol(optarg, (char **)NULL, 10);
                  if (max_gap<=1) optErr('u');
                  break;
        case 's': min_score=strtol(optarg, (char **)NULL, 10);
                  if (min_score==0) optErr('s');
                  break;
        case 'l': perc_len=strtol(optarg, (char **)NULL, 10);
                  if (perc_len==0) optErr('l');
                  break;
        case 'L': fixed_len=strtol(optarg, (char **)NULL, 10);
                  if (fixed_len==0) optErr('L');
                  break;
        case 'o': strcpy(outfile, optarg);
                  if (strlen(outfile)==0) 
                          optErr('o');
                     else do_write=1;
                  break;
        case 'X': swap_poly=1;break;
        case 'Y': swap_poly=2;break;
        case 'H': init_ends=1;break;
        case 'D': c_show=1;break;
        case 'h': showUsage();exit(1);
        case 'R': relative_shift=1;
       }
  } /* while options */

  if (swap_poly==1) {
     seed_5=SEED_A;
     seed_3=SEED_T;
     chr_5='A';
     chr_3='T';     
     }
  bufpos=seq;
  if (do_write) {
    if ((fout=fopen(outfile,"w"))==NULL) optErr('o');
    }
  for (;;) {
   /* is this a sequence name? */
   p_end=fgets(inbuf, BUF_LEN-1,stdin);
   /* basic sanity check: */
   if (!p_end) {   
     p=strchr(inbuf, '\n');
     if (!p || *(p+1)!='\0') {
        fprintf(stderr, "fgets error! Line too long?!\n");
        exit(1);
        }
     }
   p=strchr(inbuf, '>');
   if (p || !p_end) { /* starts new sequence, accumulated old one */
      if (bufpos!=seq) { /* close and analyse existing previous sequence */
        /* printf("Final: seqsize=%d\n",seq_size); */
        feat.seq=seq;
        feat.seqlen=seq_size;
        feat.seq_name=seq_name;
        if (init_ends) {/* parsed initial coordinates */
           feat.end5=init_end5;
           feat.end3=init_end3;                     
          }
          else {
           feat.end5=1;
           feat.end3=feat.seqlen;
           }
        if (skip_poly)
               N_trim(); /* only test for N_trim  */
          else {  
            if (poly_first==0) N_trim();
            save5=feat.end5;
            save3=feat.end3;
            if (trim_end!=5) score1=analyse3();
            if (trim_end!=3) score2=analyse5();
            if (swap_poly==2 && save5==feat.end5 &&
                save3==feat.end3) {
                   /* printf("here: Score1=%d, Score2=%d\n", score1, score2);       */
                   seed_5=SEED_A;
                   seed_3=SEED_T;
                   chr_5='A';
                   chr_3='T';
                   if (trim_end!=5) score1=analyse3();
                   if (trim_end!=3) score2=analyse5();
                   seed_5=SEED_T;
                   seed_3=SEED_A;
                   chr_5='T';
                   chr_3='A';
                   }
            /* printf("Score1=%d, Score2=%d\n", score1, score2);       */
            if (trim_end==1) {
              /* printf("score5=%d, score3=%d\n", score2,score1); */
                  if (score1<score2) feat.end3=save3;
                                else feat.end5=save5;
                  }
            if (poly_first!=0) N_trim();
            }
        if (relative_shift) {
           if (init_ends)
               printf("%s\t\%4.2f\t%d\t%d\t%d\n", 
                 seq_name, feat.perc_N, feat.end5-init_end5,init_end3-feat.end3, 
                   feat.seqlen);
              else     
                printf("%s\t\%4.2f\t%d\t%d\t%d\n", 
                 seq_name, feat.perc_N, feat.end5-1,feat.seqlen-feat.end3, 
                   feat.seqlen);
              }
           else
             printf("%s\t\%4.2f\t%d\t%d\t%d\n", 
                seq_name, feat.perc_N, feat.end5,feat.end3, 
                  feat.seqlen);                   
        if (c_show) 
          c_print_trim(seq, feat.end5, feat.end3);
          
        if (do_write && feat.perc_N<3.00 && feat.end3-feat.end5>=99) { 
           /* write only those having length>=100 and perc_N<3*/
           if (init_ends) fprintf(fout,">%s\n",seq_name);
                     else fprintf(fout,"%s",def_line);
           wpos=feat.end5-1;
           while (feat.end3-wpos>=60) {
             strncpy(line,seq+wpos,60);line[60]='\0';
             fprintf(fout,"%s\n",line);
             wpos+=60;
             }
           if (feat.end3-wpos>0) {
             strncpy(line,seq+wpos,feat.end3-wpos);
             line[feat.end3-wpos]='\0';
             fprintf(fout,"%s\n",line);
             }
           }  
        }
      if (!p_end) break; /* exit on end-of-file */  
      /* parse sequence name and, eventually, init coordinates */
      strcpy(def_line, p); /* keep defline there for -o output */
      for (p_end=p+1;(*p_end)!='\0' && !isspace(*p_end); p_end++);
      strncpy(seq_name, p+1,p_end-p-1); 
      seq_name[p_end-p-1]='\0';
      if (init_ends) { /* tab delimited fields in def_line*/
       /* after p_end (\t) there are init values */
       p_end=strchr(p,'\t');
       if (p_end==NULL) optErr('H');
       p_end++;
       p=strchr(p_end,'\t');
       if (p==NULL) optErr('H');
       *p='\0';
       init_end5=strtol(p_end, (char **)NULL, 10);
       if (init_end5==0) optErr('H');       
       p++;
       p_end=strchr(p,'\n');
       if (!p_end) p_end=strchr(p,'\t');
       if (p_end==NULL) optErr('H');
       *p_end='\0';
       init_end3=strtol(p, (char **)NULL, 10);
       if (init_end3==0 || init_end3<init_end5) optErr('H');
       }
      seq_size=0;
      bufpos=seq; /* no need to reallocate the buffer for each sequence..
                     let it stay at maximum all along */
    }
   else { /* extending current sequence : */
      readlen=strlen(inbuf);
      if (seq_size+readlen>alloc_size) { /* must realloc, it's gonna blow */
        alloc_size+=BUF_LEN;
        seq=realloc(seq, alloc_size);
        /* printf("Realloc! seqsize=%d, readlen=%d\n",seq_size, readlen); */
        bufpos=seq+seq_size;
        }
      bufpos=seq_upper(bufpos,inbuf, &seq_size); 
         /* append, removing newline chars and incrementing seq_size accordingly */
      }
   }/* for input */
free(seq);
if (do_write) fclose(fout);
return 0;
}

char* seq_upper(char *dest, char *src, int *seq_size) {
/* assume destination already has room for all the chars in src! 
  skip any blanks or newlines on copying 
  RETURN: the position immediately after the last char copied!
  */
int i,j;
j=0;
for (i=0;(src[i]!='\0');i++) {
 if (src[i]==' ' || src[i]=='\n'|| src[i]=='\t') continue;
 dest[j]=toupper(src[i]);j++;
 }
dest[j]='\0'; 
(*seq_size) += j;
return &dest[j]; 
}



/*********************************************************************/
/*                        Helper functions                           */
/*********************************************************************/



void optErr(char opt) {
 showUsage();
 fprintf(stderr, "\nIncorrect value specified for -%c option!\n", opt);
 exit(1);
}

void showUsage() {
 fprintf(stderr, USAGE);
}


char* lstrstr(char* lstart, char *rend, char* substr) {  /*like strstr, but starts searching
 from lstart, going up to rend and returns a pointer to the first (left) 
 matching character in str */
 char *p;
 int l,i;
 l=strlen(substr);
 p=lstart;
 while (p+l<rend) {
     for (i=0;i<l;i++) if (*(p+i) != *(substr+i)) break;
     if (i==l) return p;
     p++;
     }
 return NULL;
 }


char* rstrstr(char* rstart, char *lend, char* substr) {  /*like strstr, but starts searching
 from right end, going up to lend and returns a pointer to the last (right) 
 matching character in str */
 char *p;
 int l,i;
 l=strlen(substr);
 p=rstart-l+1;
 while (p>=lend) {
    for (i=0;i<l;i++) if (*(p+i) != *(substr+i)) break;
    if (i==l) return p+l-1;
    p--;
    }
 return NULL;
 }

/*=================================================*/

int analyse5() { /*==  end5 trim (polyT if no -X) ======*/
 int gaplen, locmaxpos, 
     locmaxscore, locminscore, run, locmaxrun, v;
 char *p;
 int score;
 int first_pass=0;
 int l5=0; /* offset from current end5 */
 int norange=1;
 int maxhit_score=0;
 int maxhit_pos=0;
 int startpos=0;
 int pos5=l5;  /* distance to end5 */
 char *seqstart=feat.seq+feat.end5-1;
 char *seqend=feat.seq+feat.end3-1;
 while (pos5<feat.end3-feat.end5 && 
         (p=lstrstr(seqstart+pos5, seqend, seed_5))!=NULL &&
           (p-seqstart)<(feat.end3-feat.end5)/2 /* while hit is in 5' half! */
          ) { /* big loop, going through all possible hits */
   startpos=p-seqstart; /*offset from end5 */
   locmaxscore=0;
   locmaxpos=l5;
   run=0;
   locmaxrun=0;
   p+=SEED_LEN;
   /*
   printf("Seed found: pos5=%d, (p-feat.seq)=%d, *p=%c%c%c%c%c\n",
          pos5, p-feat.seq, *p,*(p+1),*(p+2),*(p+3),*(p+4));   
   */            
   score=SEED_LEN;
   run=SEED_LEN;
   gaplen=0;
   while (gaplen<max_gap && score>0 && p+gaplen<=seqend) {
    if (p[gaplen]==chr_5) {
          if (gaplen) { /* new poly run after a gap  */
             p+=gaplen;
             score-=gaplen*2; /* penalty */
             gaplen=0;
             }
            p++;
            run++;
            score++;
            }
         else { /* non-match */
          if (gaplen==0) { /* new gap starting or end of sequence*/
             if (run>locmaxrun) locmaxrun=run;
             run=0;          
             if (score>=locmaxscore) { /* check locmax! */
              locmaxscore=score;
              locmaxpos=(p-seqstart); /* offset from end5 */
              }
             }
          gaplen++;
          }   
    /*                  
      printf("T loop: l5=%d, gaplen=%d, score=%d, locmaxpos=%d, locmaxscore=%d locmaxrun=%d\n   (p-feat.seq)=%d, *p=%c%c%c%c%c\n",
       l5, gaplen, score, locmaxpos, locmaxscore, locmaxrun, p-feat.seq, *p,*(p+1),*(p+2),*(p+3),*(p+4));         
     */         
    }/* while */
    if (run>locmaxrun) locmaxrun=run;
    if (score>=locmaxscore) { /* check locmax! */
        locmaxscore=score;
        locmaxpos=(p-seqstart); /* offset from end5 */
        }
   pos5=locmaxpos;
   locminscore=min_score+startpos/10; /* min_score is dynamically adjusted
         based on the hit distance to the ends */
   /*
   printf("feat.end3=%d, feat.end5=%d, l5=%d, startpos=%d, locmaxpos=%d, locmaxscore=%d, locmaxrun=%d\n",
           feat.end3,    feat.end5,    l5,   startpos, locmaxpos, locmaxscore, locmaxrun);
   
   */     
   if (locmaxscore>=locminscore && locmaxrun>=min_run) { /* potential hit */
         norange=1;
         if (fixed_len) {
                 if (startpos<=fixed_len) {
                                  l5=locmaxpos;
                                  maxhit_score=locmaxscore;
                                  maxhit_pos=startpos;
                                  norange=0;
                                  }
                                  else break; /* strict fixed range criteria */
                 }
               else { /* percentual distance test: assess distance to previous hit */                   
                  v=(feat.end3-feat.end5-l5)*perc_len/100;
                  if ((first_pass && v<MIN_PERC_LEN)) {
                        v=MIN_PERC_LEN;
                        first_pass=0;
                        }                       
                  if (startpos-l5<v) { /* hit in allowed range */
                                     l5=locmaxpos;
                                     maxhit_score=locmaxscore;
                                     maxhit_pos=startpos;
                                     norange=0;
                                     }
                  }
         /*  test for any huge hits in the first half */
         if (norange && l5!=locmaxpos && 
                   ((locmaxscore>30 && startpos+(locmaxpos-startpos)/2 < (feat.end3-feat.end5)/2)
                     || (feat.end3-feat.end5)-(locmaxpos-startpos)<100) ) 
                            {
                            l5=locmaxpos;
                            maxhit_score=locmaxscore;
                            maxhit_pos=startpos;
                            }
       }
   } /* big while loop */
 feat.end5+=l5;
 return maxhit_score*10-maxhit_pos;
}    



int analyse3() { /*===  end3 search: polyA (if no -X) =============== */
 int gaplen, score, locmaxpos, 
     locmaxscore, locminscore, run, locmaxrun, v;
 char *p;
 int first_pass=0;
 int norange;
 char *seqstart=feat.seq+feat.end5-1;
 char *seqend=feat.seq+feat.end3-1;
 int l3=0; /*offset from end3 */
 int pos3=l3; /*distance to end3 */
 int startpos=0;
 int maxhit_score=0;
 int maxhit_pos=0;
 while (pos3<feat.end3-feat.end5 &&
         (p=rstrstr(seqend-pos3, seqstart, seed_3))!=NULL &&
           (seqend-p)<(feat.end3-feat.end5)/2 /* while hit is in 3' half! */
          ) { /* big loop, going through all possible hits */
   /* printf("feat.end3=%d, pos3=%d, startpos=%d\n",feat.end3,pos3, seqend-p); */
   startpos=seqend-p; /*distance to seqend */
   locmaxscore=0;
   locmaxpos=l3;
   run=0;
   locmaxrun=0;
   /*
   printf("Seed found *p=%c%c%c%c%c%c%c\n",
               *p, *(p+1), *(p+2), *(p+3), *(p+4), *(p+5),*(p+6));
   */
   p-=SEED_LEN; 
   /*
   printf("After seed *p=%c%c%c%c%c%c%c\n",
               *p, *(p+1), *(p+2), *(p+3), *(p+4), *(p+5),*(p+6));
   */
   /*
   printf("Seed found: l5=%d, (p-feat.seq)=%d, *p=%c%c%c%c%c\n",
          l5, p-feat.seq, *p,*(p+1),*(p+2),*(p+3),*(p+4) );   
   */           
   score=SEED_LEN;
   run=SEED_LEN;
   gaplen=0;
   while (gaplen<max_gap && score>0 && p-gaplen>=seqstart) {   
     /*printf("A loop: l3=%d, gaplen=%d, score=%d, locmaxpos=%d, locmaxscore=%d, locmaxrun=%d\n \
       (p-feat.seq)=%d, *p=%c%c%c%c%c%c%c\n",
       l3, gaplen, score, locmaxpos, locmaxscore, locmaxrun, 
        p-feat.seq, *(p-gaplen),*(p-gaplen+1),*(p-gaplen+2),*(p-gaplen+3),*(p-gaplen+4),*(p-gaplen+5),
         *(p-gaplen+6));         
     */
     if (*(p-gaplen)==chr_3) { 
          if (gaplen) { /* new poly run after a gap  */
             p-=gaplen;
             score-=(gaplen<<1); /* penalty */
             gaplen=0;
             }
            p--;
            run++;
            score++;
            /*
            printf("match: score=%d, *p=%c%c%c%c%c%c%c\n",score,
               *p, *(p+1), *(p+2), *(p+3), *(p+4), *(p+5),*(p+6));
               */
            }
         else { /* non-match*/
          if (gaplen==0) { /* new gap */
             if (run>locmaxrun) locmaxrun=run;
             run=0;          
             if (score>=locmaxscore) { /* check locmax! */
              locmaxscore=score;
              locmaxpos=(seqend-p); /* distance to end3 */
              }
             }
          gaplen++; /*extend gap */
          }
                      
              
     }/* while */
     if (run>locmaxrun) locmaxrun=run;
     if (score>=locmaxscore) { /* check locmax! */
              locmaxscore=score;
              locmaxpos=(seqend-p); /* distance to end3 */
              }
   pos3=locmaxpos;
   locminscore=min_score+startpos/10; /* min_score is dynamically adjusted
         based on the hit distance to the ends */
   if (locmaxscore>=locminscore && locmaxrun>=min_run) { /* potential hit */
         norange=1;
         if (fixed_len) {
                 if (startpos<=fixed_len) {
                                 l3=locmaxpos;
                                 maxhit_score=locmaxscore;
                                 maxhit_pos=startpos;
                                 norange=0;
                                 }
                                else break; /* strict fixed range criteria */
                 }
               else {
                   v=((feat.end3-l3-feat.end5)*perc_len/100);
                   if (first_pass && v<MIN_PERC_LEN) {
                        v=MIN_PERC_LEN;
                        first_pass=0;
                        }
                   if (startpos-l3<v) {
                                l3=locmaxpos;
                                maxhit_score=locmaxscore;
                                maxhit_pos=startpos;
                                norange=0;
                                }
                  }
                 /* instead, just aim for big hits in the 3' half */
                 if (norange && l3!=locmaxpos && 
                      ( (locmaxscore>30 && startpos+(locmaxpos-startpos)/2 < (feat.end3-feat.end5)/2)
                         || (feat.end3-feat.end5)-(locmaxpos-startpos)<100) ) {
                      l3=locmaxpos;
                      maxhit_score=locmaxscore;
                      maxhit_pos=startpos;
                      }
         } /* potential hit analysis */
   } /* big while loop */
 feat.end3-=l3;
 return maxhit_score*10-maxhit_pos;
}


double N_calc(int e5, int e3, int* Ncount) { 
/* computes the percent and count of Ns; updates feat structure */
int i;
double result;
*Ncount=0;
for (i=e5;i<=e3;i++)
  if (feat.seq[i]!='A' && feat.seq[i]!='C' && feat.seq[i]!='G' && feat.seq[i]!='T')
      (*Ncount)++;
result = ((double)(100*(*Ncount)))/(e3-e5+1);
return result;
}

void N_trim() {
/* computes initial perc_N and fill N_Pos
   call N_analyse as apropriate (according to switches)
  */
int i,j; 
feat.perc_N=N_calc(feat.end5-1, feat.end3-1, &feat.NCount);
/* printf("entering N trim: %f compare to %f\n",feat.perc_N, perc_N); */
if (perc_N!=0 && feat.perc_N>perc_N) {
  j=0;
  feat.NPos=malloc(feat.NCount*sizeof(int));
  for (i=feat.end5-1;i<feat.end3;i++) 
     if (feat.seq[i]!='A' && feat.seq[i]!='C' 
        && feat.seq[i]!='G' && feat.seq[i]!='T') {
           feat.NPos[j]=i;  
           j++;
           }
  N_analyse(feat.end5-1, feat.end3-1, 0, feat.NCount-1);
  /* compute perc_N again after trimming! */
  feat.perc_N=N_calc(feat.end5-1, feat.end3-1, &feat.NCount);
  if (feat.perc_N >= nsave_limit && feat.end3-feat.end5>=99) { 
     /* try to save this at any cost: */
     N_save(feat.end5-1, feat.end3-1); 
     /* this will also update feat.perc_N and feat.NCount */
     }
  free(feat.NPos);
  }
}

void N_analyse(int l5, int l3, int p5, int p3) {
/* assumes feat was filled properly */
 int old_dif, t5,t3,v;
 if (l3-l5<100 || p5>p3 ) {
   feat.end5=l5+1;
   feat.end3=l3+1;
   return; 
   }

 t5=feat.NPos[p5]-l5;
 t3=l3-feat.NPos[p3];
 old_dif=p3-p5;
 v=(int)((((double)(l3-l5))*perc_lenN)/100);
 if (v>30) v=30; /* enforce limit for long ESTs */ 
 if (t5 < v ) {
   l5=feat.NPos[p5]+1;
   p5++;
   }
 if (t3 < v) {
   l3=feat.NPos[p3]-1;
   p3--;
   }
 /* restNs=p3-p5; number of Ns in the new CLR */
 if (p3-p5==old_dif) { /* no change, return */
           /* don't trim if this may shorten < 100 */
           if (l5-l3<100) return;
           feat.end5=l5+1;
           feat.end3=l3+1;
           }
    else 
      N_analyse(l5,l3, p5,p3);           
}

void N_save(int l5, int l3) {
/* find the maximum range having percN<3 */
/* needs feat.NPos[], feat.NCount to be set before*/
 int countN,i,p5,p3, max_extent, start3, max_5, max_3, e3, e5, max_countN;
 double percN,max_percN;
 for (i=0;feat.NPos[i]<=feat.end5-1;i++);
 p5=i-1; /* found first N's position in the selected range */
 for (i=feat.NCount-1;feat.NPos[i]>=feat.end3-1;i--);
 p3=i+1; /* found last N's position in the selected range */ 
 if (l3-l5<100 || p5>p3 ) return;
 e5=p5;e3=p3;
 start3=l3;
 max_extent=99;max_3=0;max_5=0;
 while (l3-l5>max_extent /* && e3>=e5 */) { /* outer end5 loop */
    e3=p3;    
    percN=N_calc(l5,l3, &countN);
    while (percN>=nsave_limit && l3-l5>max_extent) { /* e3 loop -- decreases until percN<3 */
       /* printf("l5=%d, l3=%d; e5=%d (%d), e3=%d(%d), percN=%f, %d ? %d\n",l5,l3,e5,
          feat.NPos[e5], e3, feat.NPos[e3], percN, l3-l5, max_extent); 
          */
      e3--;
      while (e3>e5 && feat.NPos[e3-1]==feat.NPos[e3]-1) e3--; 
               /* skip consecutive Ns */
      l3=feat.NPos[e3]-1;
      percN=N_calc(l5,l3, &countN);
      }
    if (l3-l5>max_extent) {
      max_extent=l3-l5;
      max_5=l5;max_3=l3;
      max_percN=percN;
      max_countN=countN;
      }
    while (e3>e5 && feat.NPos[e5+1]==feat.NPos[e5]+1) e5++; 
              /* skip consecutive Ns */
    l5=feat.NPos[e5]+1;
    /* restart end3: */
    l3=start3;e3=p3;
    e5++;
  } /* outer loop */
 if (max_3!=0) {
    feat.end5=max_5+1;
    feat.end3=max_3+1;
    feat.perc_N=max_percN;
    feat.NCount=max_countN;
  }
}


