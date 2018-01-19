#ifndef _color_c_h_
#define _color_c_h_
           
enum {c_black=0,
      c_red, c_green,c_brown,c_blue,c_magenta,c_cyan,c_white
      };
      
void c_print_trim(char *seq, int end5, int end3);

void c_setfg(int c);
void c_setbg(int c);
void c_resetfg();
void c_resetbg();
void c_reset();
void c_setbold();

void c_setbright();

void c_setblink();

void c_setreverse();

void c_setunderline();

void c_setnormal();

#endif
