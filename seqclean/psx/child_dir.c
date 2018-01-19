#include "psx.h"
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
/* #include <sys/exec.h> */
#include <signal.h>

/* Child process ID array */
static pid_t childID[MAX_PROCS];

/* Child process directories */
static char *childdir[MAX_PROCS];

/* Number of procs */
static int stat_numprocs;

extern void wait_then_exit();

void
init_cur_dir(execdir, numprocs)

char	*execdir;
int	numprocs;
{
int	i, dirlen;
FILE	*testfp;

	stat_numprocs = numprocs;
	dirlen = strlen(execdir) + 20;
	for (i = 0; i < numprocs; i++) {
		childID[i] = (pid_t)0;
		if ((childdir[i] = (char *)malloc((unsigned)(2*dirlen)))== NULL) {
			fprintf(stderr, "ERROR: Could not allocate memory for childdir! (errno = %d)\n", errno);
			wait_then_exit(4);
		}
		sprintf(childdir[i], "%s_%d", execdir, (i + 1));
		if (mkdir(childdir[i], (mode_t)0750) < 0) {
			if (errno == EEXIST) {
				strcat(childdir[i], "/test_psx");
				if ((testfp = fopen(childdir[i], "w")) != NULL) {
					fclose(testfp);
					unlink(childdir[i]);
					sprintf(childdir[i], "%s_%d", execdir, (i + 1));
					fprintf(stderr, "WARNING: directory %s already exists!\n", childdir[i]);
					continue;
				}
			}
			fprintf(stderr, "ERROR: Could not make directory %s! (errno = %d)\n", childdir[i], errno);
			wait_then_exit(4);
		}
	}
	return;
}

char *
get_cur_dir(child)

pid_t	child;
{
int	i;

	for (i = 0; i < stat_numprocs; i++) {
		if (childID[i] == child) {
			return (childdir[i]);
		}
	}
	return ((char *)NULL);
}

void
change_cur_dir(child, curdir)

pid_t	child;
char	*curdir;
{
int	i;

	for (i = 0; i < stat_numprocs; i++) {
		if (strcmp(childdir[i], curdir) == 0) {
			childID[i] = child;
			return;
		}
	}
	fprintf(stderr, "ERROR: Nonexistent directory %s !\n", curdir);
	wait_then_exit(4);
}
