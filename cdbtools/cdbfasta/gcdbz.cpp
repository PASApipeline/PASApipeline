#include "gcdbz.h"
#include "gcl/gcdb.h"

GCdbz::GCdbz(FILE* azf, bool uc, int zrsize) {
 uncompress=uc;
 zrecsize=-1;
 zpos=0;
 defline_cap=1024;
 begin_defline();
 GMALLOC(defline, defline_cap);
 zf=azf;
 // FULL_FLUSH method instead of finish:
 if (uncompress)
     decomp_start(zrsize);
   else
     compress_start();
}

GCdbz::~GCdbz() {
 //if (zf!=NULL && zf!=stdout && zf!=stdin) fclose(zf);
 // FULL_FLUSH method instead of finish
  if (uncompress) decomp_end();
   else
     if (!zclosed) compress_end();
 GFREE(defline);
}



void GCdbz::extend_defline(int ch) {
 if (defline_len+1 >= defline_cap) {
   defline_cap+=(defline_cap>>2);
   GREALLOC(defline, defline_cap);
   }
 defline[defline_len]= ch;
 defline_len++;
 }


#define DUMMY_ZREC ">AA1234567890 DNA protein\n\
ACGTTGCTAGCT\n\
NRMTPYYHEIEP\n\
RTASNTSPTPNS\n\
IKSAHPAEPPKR\n"

void GCdbz::compress_start() {
 //initialize zstream compression
 zstream.zalloc = (alloc_func)0; //no alloc function to use
 zstream.zfree = (free_func)0;   //no free function to use
 zstream.opaque = (voidpf)0;     //no private object to pass to zalloc/zfree

 int err=deflateInit(&zstream, Z_DEFAULT_COMPRESSION);
 if (err!=Z_OK)
     GError("GCdbz error: deflateInit failed!(err=%d)\n",err);
 zclosed=false;
 //write a dummy record as the first record, 
 //so we can use random access (FULL_FLUSH style) later
 char ztag[5];strcpy(ztag, "CDBZ");
 uint32 zsize=0;
 zstream.next_in = (Bytef*)sbuf;
 strcpy(sbuf, DUMMY_ZREC);
 zstream.avail_in=strlen(sbuf);
 zstream.next_out = (Bytef*)lbuf;
 zstream.avail_out = GCDBZ_LBUF_LEN;
 uLong t_out=zstream.total_out;
 err = deflate(&zstream, Z_FULL_FLUSH);
 zsize=zstream.total_out-t_out;
 if ((err !=Z_OK && err!=Z_STREAM_END) || zsize<=0)
       GError("GCdbz error: deflate 1st record failed! (err=%d)\n", err);
 //now write the header and the dummy record
     //in case this was not done before:
 gcvt_uint=(endian_test())? &uint32_sun : &uint32_x86;
 uint32 zfv = gcvt_uint(&zsize);
 if (fwrite(ztag, 1, 4, zf)<4 || 
       fwrite(&zfv,1,sizeof(uint32), zf) < sizeof(uint32) ||
         fwrite(lbuf, 1, zsize, zf) < zsize)
       GError("Error writing 1st deflated record!\n");
 zpos+=4+sizeof(uint32)+zsize;
 }

void GCdbz::compress_end() {
 zstream.next_out = (Bytef*)lbuf;
 zstream.avail_out = GCDBZ_LBUF_LEN;
 zstream.avail_in = 0;
 uLong t_out=zstream.total_out;
 int err = deflate(&zstream, Z_FINISH);
 if (err != Z_STREAM_END) {
   GError("GCdbz error: deflate/Z_FINISH() failed! (err=%d) \n", err);
   }
 uLong toWrite=zstream.total_out-t_out;
 if (toWrite>0) {
   if (fwrite(lbuf, 1, toWrite, zf)<toWrite)
        GError("Error writing FINISH deflate chunk!\n");
   //GError("GCdbz error: out data after Z_FINISH (%d bytes)\n",
   //    zstream.total_out-t_out);
   }
 err=deflateEnd(&zstream);
 if (err!=Z_OK)
   GError("GCdbz error: deflateEnd() failed! (err=%d) \n", err);   
 zclosed=true;   
}

char* GCdbz::compress(GReadBuf *readbuf, char* delim) {
  //compress everything coming from the input stream inf
  //until \n is encountered followed by delim
  //returns this->defline or NULL if error encountered

  //-- WARNING: this subrutine assumes that inf file position
  // is at the beginning of the record, right AFTER the delim
  // (exactly as left after a previous call)
 if (zf==NULL || uncompress)
    GError("GCdbz Error: cannot use compress() method !\n");
 unsigned int total_out=0;
 int c=0;
 bool in_rec=true;
 int delimlen=strlen(delim);
 zrecsize=0;
 if ((c=readbuf->peekCmp(delim, delimlen))!=0) {
    if (c<-1) return NULL; //end of file reached
    GError("GCdbZ::compress error: delimiter '%s' expected at record start!\n",
          delim);
    }
 bool bol=false; //beginning of line flag
 int deflate_flag=0;
 begin_defline();
 int rec_pos=0;
 int err=0;
 while (in_rec) { // main read loop
     int bytes_read=0;
     while ((c=readbuf->getch())>=0) {
       sbuf[bytes_read++]=c;
       if (c=='\n' || c=='\r') { //beginning of line
          bol = true;
          if (in_defline) end_defline();
          //look_ahead for record delimiter:
          if (readbuf->peekCmp(delim, delimlen)==0) {
              in_rec=false;
              break;
              }
          }
        else bol = false;
       if (rec_pos>delimlen-1 && in_defline)
            extend_defline(c);
       rec_pos++;
       if (bytes_read == GCDBZ_SBUF_LEN) break;     
       }//while not EOF or space in buffer
     /*if (bytes_read==0)
           return NULL;*/
     if (c==EOF) {
        in_rec=false;
        if (in_defline) end_defline();
        }
     zstream.next_in = (Bytef*)sbuf;
     zstream.avail_in = bytes_read;
     //deflate_flag = in_rec ? 0 : Z_FINISH;
     deflate_flag = in_rec ? 0 : Z_FULL_FLUSH;
     do { //compression loop
        zstream.next_out = (Bytef*)lbuf;
        zstream.avail_out = GCDBZ_LBUF_LEN;
        uLong t_out=zstream.total_out;
        err = deflate(&zstream, deflate_flag);
        if (err !=Z_OK && err!=Z_STREAM_END) 
             GError("GCdbz error: deflate failed! (err=%d)\n", err);
        uLong toWrite=zstream.total_out-t_out;
        if (toWrite>0) {
             if (fwrite(lbuf, 1, toWrite, zf)<toWrite)
                GError("Error writing deflate chunk!\n");
             total_out+=toWrite;
             zrecsize+=toWrite;
             zpos+=toWrite;
             }
       } while (err!=Z_STREAM_END && zstream.avail_out==0);//compression loop
   } //read loop
  //if (deflate_flag!=Z_FINISH)
  if (deflate_flag!=Z_FULL_FLUSH)
     GError("Deflate flag not set to FINISH!\n");
  return defline;
}


void GCdbz::decomp_start(int zrsize) {
 zstream.zalloc = (alloc_func)0;
 zstream.zfree = (free_func)0;
 zstream.opaque = (voidpf)0;
 zstream.next_in  = (Bytef*)sbuf;
 zstream.avail_in = 0;
 zstream.next_out = (Bytef*)lbuf;
 int err = inflateInit(&zstream);
 if (err!=Z_OK)
     GMessage("Error at inflateInit()\n");
 //-- now read and discard the first record, so we can use random access later
 // (needed by zlib)
 int bytes_read=fread(sbuf, 1, zrsize, zf);
 if (bytes_read<zrsize)
     GError("Error reading 1st record from zrec file\n");
 zstream.next_in = (Bytef*)sbuf;
 zstream.avail_in = bytes_read;
//decompress first chunk
 zstream.next_out = (Bytef*)lbuf;
 zstream.avail_out = GCDBZ_LBUF_LEN;
 err = inflate(&zstream, Z_SYNC_FLUSH);
 if (err !=Z_OK && err!=Z_STREAM_END)
     GError("GCdbz error: 1st record inflate failed! (err=%d)\n",err);     
}

void GCdbz::decomp_end() {
  int err = inflateEnd(&zstream);
  if (err!=Z_OK)
     GError("Error at inflateEnd() (err=%d)\n", err);

}


//record decompress
//returns: the number of bytes decompressed
int GCdbz::decompress(FILE* outf, int csize, int zfofs) {
 if (zfofs>=0) {
    if (fseek(zf, zfofs, 0))
      GError("GCdbz::decompress: error fseek() to %d\n", zfofs);
    }
  else
     if (feof(zf)) return 0;
 bool in_rec=true;
 int err=0;
 int total_read=0;
 int total_written=0;
 while (in_rec) { // main read loop
     int to_read=0;
     int bytes_read=0;
     if (csize<=0) { //read one byte at a time
        to_read=1;
        int c;
        if ((c =fgetc(zf))!=EOF) {
           bytes_read = 1;
           sbuf[0]=c;
           }
          else {
            //bytes_read=0;
            return 0; //eof
            }
        total_read+=bytes_read;
        }
      else {
        to_read = csize-total_read>GCDBZ_SBUF_LEN ?
                                 GCDBZ_SBUF_LEN : csize-total_read;
       // check for csize vs bytes_read match:
        if (to_read==0) return 0;
        bytes_read=fread(sbuf, 1, to_read, zf);
        if (bytes_read!=to_read)
            GError("Error reading from zrec file\n");
        total_read+=bytes_read;
        in_rec=(total_read<csize);
        }
     if (bytes_read==0) {
        //GMessage("bytes_read = 0\n");
        return 0;
        }
     if (in_rec && bytes_read<to_read) in_rec=false;
     zstream.next_in = (Bytef*)sbuf;
     zstream.avail_in = bytes_read;

     do { //decompression loop
        zstream.next_out = (Bytef*)lbuf;
        zstream.avail_out = GCDBZ_LBUF_LEN;
        uLong t_out=zstream.total_out;
        err = inflate(&zstream, Z_SYNC_FLUSH);
        uLong toWrite=zstream.total_out-t_out;
        if (toWrite>0) {
             if (fwrite(lbuf, 1, toWrite, outf)<toWrite) {
               GError("Error writing inflated chunk!\n");
               }
             total_written+=toWrite;
             }
        if (err==Z_STREAM_END) {
              in_rec=false;
              if (total_written==0) {
                GMessage("Z_STREAM_END found but total_written=0!\n");
                }
              break;
              }
         else if (err !=Z_OK)
                GError("GCdbz error: inflate failed! (err=%d)\n",err);
        } while (zstream.avail_in!=0); //decompression loop
   } //read loop
 /*if (err!=Z_STREAM_END) {
   GError("decompress: Z_STREAM_END not found!\n");
   }*/
  return total_written;
}

