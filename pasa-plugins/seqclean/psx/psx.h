#define True 1
#define False 0

#define MAX_PROCS 16
#define IN_LINE_SIZE 4096


/*
#define WCONTFLG                0177777
#define WWORD(stat)             ((int)((stat))&0177777)
#define WIFCONTINUED(stat)      (WWORD(stat) == WCONTFLG)

#define  _W_INT(i) (*(int *)&(i))
*/



/* evaluates to the low-order 8 bits of the child exit status   
#define WEXITSTATUS(x)  (WIFEXITED(x) ? ((_W_INT(x) >> 8) & 0xff00) : -1)

#define WEXITSTATUS(x) -1
*/
