#include "color_c.h"
#include <stdio.h>
#include <string.h>

char buf[6]={0x1B,'[', 'n','m','m','\0'};

void c_setfg(int c) {
 sprintf((char *)(&buf[2]),"%dm",c+30);
 fwrite(buf,1,strlen(buf), stdout);
}

void c_setbg(int c) {
 sprintf((char *)(&buf[2]),"%dm",c+40);
 fwrite(buf,1,strlen(buf), stdout);
};

void c_resetfg() {
 sprintf((char *)(&buf[2]),"39m");
 fwrite(buf,1,strlen(buf), stdout);
};
void c_resetbg() {
 sprintf((char *)(&buf[2]),"49m");
 fwrite(buf,1,strlen(buf), stdout);
}
void c_reset() {
 sprintf((char *)(&buf[2]),"0m");
 fwrite(buf,1,strlen(buf), stdout);
};

void c_setbold() {};

void c_setbright() {};

void c_setblink() {};

void c_setreverse() {};

void c_setunderline() {};

void c_setnormal() {
 sprintf((char *)(&buf[2]),"22m");
 fwrite(buf,1,strlen(buf), stdout);
};

/*
sub txc_resetfg {
print chr(27),"[39m";
}

sub txc_resetbg {
print chr(27),"[49m";
}

sub txc_reset {
print chr(27),"[0m";
}

sub txc_setbold {
print chr(27),"[1m";
}

sub txc_setbright {
print chr(27),"[2m";
}

sub txc_setblink {
print chr(27),"[5m";
}

sub txc_setreverse {
print chr(27),"[7m";
}

sub txc_setunderline {
print chr(27),"[4m";
}

sub txc_setnormal {
print chr(27),"[22m";
}
*/

void c_print_trim (char * seq, int end5, int end3) 
{
 int pos;
 int len, ex;
 char buf[81];
 if (end5>0) {
   c_reset(); 
   c_setbg(c_red); 
   }
 /*-----------first colored part: */
 len=end5-1;pos=0;
 while (len-pos>=60) {
     strncpy(buf, seq+pos, 60);buf[60]=0;
     printf("%s\n",buf);
     pos+=60;
    }
 if (len-pos>0) {
     ex=len-pos;
     strncpy(buf, seq+pos,ex);buf[ex]=0;
     printf("%s",buf);
     pos+=len-pos; 
    }  
 /* -----------middle part:  */
 c_reset(); /* restore normal color */
 len=end3;
 if (len/60 > pos/60) {
    ex=(pos/60+1)*60-pos;
    strncpy(buf, seq+pos, ex);buf[ex]=0;
    printf("%s\n",buf);
    pos=(pos/60+1)*60; 
    while (len-pos>=60) {
       strncpy(buf, seq+pos, 60);buf[60]=0;
       printf("%s\n",buf);
       pos+=60;
       }
    if (len-pos>0) {
       ex=len-pos;
       strncpy(buf, seq+pos, ex);buf[ex]=0;
       printf("%s",buf);
       pos+=len-pos; 
       }   
    }
 else {
    ex=len-pos;
    strncpy(buf, seq+pos, ex);buf[ex]=0;
    printf("%s\n",buf);
    pos+=len-pos;
    }
/*----------end colored part */
 len=strlen(seq);
 if (end3<len) {
    c_setbg(c_red); /* red again */
    if (len/60 > pos/60) {
      ex=(pos/60+1)*60-pos;
      strncpy(buf, seq+pos, ex);buf[ex]=0;
      printf("%s\n",buf);
      pos=(pos/60+1)*60; 
      while (len-pos>=60) {
         strncpy(buf, seq+pos, 60);buf[60]=0;
         printf("%s\n",buf);
         pos+=60;
         }
      if (len-pos>0)  {
         ex=len-pos; 
         strncpy(buf, seq+pos, ex);buf[ex]=0;
         printf("%s",buf);
         pos+=len-pos;
         }  
      }
     else {
      ex=len-pos;
      strncpy(buf, seq+pos, ex);buf[ex]=0;
      printf("%s",buf);
      }
    c_reset();
    }
printf("\n"); 
}
