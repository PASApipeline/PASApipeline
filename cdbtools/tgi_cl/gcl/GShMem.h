#ifndef _GShMem_h_
#define _GShMem_h_
#define _XOPEN_SOURCE 1
#include "gcl/GBase.h"
#include <string.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>


class GShMem {
   int shmsize;
   int shmid;
   char* shname;
   char* mem;
   shmid_ds shminfo;
   bool isNew;
   int logpos;
 protected:
   //
   key_t getKey(const char* name) {
      /*
      char* homepath=getenv("HOME");
      if (homepath==NULL)
        GError("Error: cannot get HOME environment variable.");
      char* path;
      GMALLOC(path, (strlen(homepath)+strlen(name)+8));
      strcpy(path, homepath);
      strcat(path, "/.SHM_");
      strcat(path,name);
      key_t key=ftok(path,'M');      
      GMessage("path sent to ftok is: %s; key=%p\n",path, key);
      
      if (key<0) {
        perror("ftok");
        exit(1);      
        }
        
      GFREE(path);
      */
      int key=strhash(name);
      if (key<0) key=abs(key);
      //right now just don't care about collisions!
      //TO DO: a static member could keep a real hashlist for these
      return (key_t) key;
      }
 public:
  GShMem(const char* name, int size=0, bool destroy=false) {
     key_t key=getKey(name);
     logpos=0;
     shmid=-1;
     mem=NULL;
     shname=NULL;
     if (size>0) { //new shared segment creation
      shmid=shmget(key, size, IPC_CREAT | 0666);
      if (shmid<0) {
         GMessage("Error creating the shared memory block '%s'!\n", 
              name);
         perror("shmget");
         exit(1);
         }
      //shmsize=size;
      isNew=true;
      }
    else { //block retrieval, for reading 
      //GMessage("Size given=%d\n",size);
      shmid=shmget(key, size, 0666);      
      if (shmid<0) {
        GMessage("Error opening shared memory block '%s'!\n",
              name);
        perror("shmget");
        exit(1);
        }
      isNew=false;
      }
     shname=Gstrdup(name);
     //attach to it 
     if ((mem = (char*)shmat(shmid, 0, 0))<(char*)0) {
        perror("shmat");
        exit(1);
        }
          
     int r= shmctl(shmid, IPC_STAT, &shminfo); //get info
     if (r<0) {
       GMessage("Error getting shared memory info for block '%s'!\n", 
              name);
       perror("shmctl");
       exit(1);
       }
     //GMessage("shmize=%d\n", shminfo.shm_segsz);
     shmsize=shminfo.shm_segsz;
     if (destroy) toDestroy();
     }
     
  ~GShMem() {
     GFREE(shname);        
     if (shmid>0 && mem!=NULL) { //detach
        if ( shmdt(mem) < 0) {
          perror("shmdt");
          }
        }
     }
 int getSize() { return shmsize;}
 const char* getPtr() { return mem; }
 const char* getName() { return shname; }

 void toDestroy() {
     //mark the shared block to be destroyed after the last detach to it
     if (shmctl(shmid, IPC_RMID, &shminfo)<0) {
         GMessage("Error setting shared memory as destroyable (block '%s')!\n", 
                shname);
         perror("shmctl");
         }
     }
 //use the shared memory as a communication buffer
 //loggin operations: must be used exclusively!
 void log(const char* msg) {
    //appends to the memory buffer, if enough room
    int msglen=strlen(msg);
    if (msglen>=shmsize) {msglen=shmsize-1;logpos=0;}
    if (msglen>=shmsize-logpos) {//"scroll" the buffer up at least one line    
       int need=msglen-(shmsize-logpos);
       int to=0;
       while (*(mem+need+to)!='\n' && need+to<logpos) to++;
       if (need+to==logpos) to=+1; //couldn't find a line boundary, scroll chars
                       else to++;
       //GMessage("memmove: from %d \n", need+to);
       memmove((void*)mem, (const void*) (mem+need+to), logpos-need-to);
       logpos-=need+to;
       }
     //logpos is set safely; append the given string  
     //if (strcmp(shname, "sumcclust")==0)
     //  GMessage("shmsize=%d: at logpos=%d (s=\"%s\", msglen=%d)\n",shmsize, logpos, msg, msglen);
     memmove((void*)&mem[logpos], (const void*)msg, msglen);
     logpos+=msglen;
     mem[logpos]='\0';
     }
};



#endif
