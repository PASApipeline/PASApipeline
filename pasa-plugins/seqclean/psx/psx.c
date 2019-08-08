/*
$Header: /usr/local/devel/psx-linux/src/RCS/psx.c,v 1.4 1998/09/25 14:10:57 fliang Exp fliang $

$Log: psx.c,v $
Revision 1.4 1998/09/25 14:10:57 fliang
modified to run on Dec alpha

Revision 1.3 1998/09/25 14:09:15 fliang
modified to run in Linux box

Revision 1.2 1998/09/25 14:07:19 fliang
sx program was modified by Granger Sutton to run in mutiple processors

Revision 1.1 1998/09/25 14:04:28 fliang
Initial revision


*/

/*
 * ===========================================================================
 *
 * PUBLIC DOMAIN NOTICE
 * National Institutes of Health
 *
 * This software/database is a "United States Government Work" under the
 * terms of the United States Copyright Act. It was written as part of
 * the author's official duties as a United States Government employee and
 * thus cannot be copyrighted. This software/database is freely available
 * to the public for use. The National Institutes of Health and the U.S.
 * Government have not placed any restriction on its use or reproduction.
 *
 * Although all reasonable efforts have been taken to ensure the accuracy
 * and reliability of the software and data, the NIH and the U.S.
 * Government do not and cannot warrant the performance or results that
 * may be obtained by using this software or data. The NIH and the U.S.
 * Government disclaim all warranties, express or implied, including
 * warranties of performance, merchantability or fitness for any particular
 * purpose.
 *
 * Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * File Name: sx.c
 *
 * Author: Joshua M. Yulish
 *
 * Version Creation Date: June 14, 1992
 *
 * $Revision: 1.4 $
 *
 * File Description:
 *
 * This program is a program which will do bulk batching. Based on an input
 * file this program will execute a specified program or script continuously
 * until all sequences specified have been run through.
 *
 * This program takes as input six flags. -n = number of sequences to group
 * together to run with (def 1). -i = name of input file otherwise default to
 * stdin. -a = the switch to append the environment variables onto the command
 * line to be executed. -s = number of sequences in input file to skip before
 * seqs are sent to executable (def 0). -t = total num of seqs to do in the
 * input file (def ALL). -c = in " " is the command to execute, the name of a 
 * script to run with the sequences. The input file must be in standard FASTA 
 * format, each sequence seperated by a comment line starting with a '>', '#'.
 *
 * Modifications:
 * --------------------------------------------------------------------------
 * Date Name Description of modification
 * ------- ---------- -----------------------------------------------------
 *
 *
 * ==========================================================================
 */

#include <fcntl.h>
#include "psx.h"
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>

/*added by ray to copy file. link function does not work on the cluster bc hard link is not permitted*/

int copy_file(const char * from,
  const char * to) {
  int fd_to, fd_from;
  char buf[4096];
  ssize_t nread;
  int saved_errno;

  fd_from = open(from, O_RDONLY);
  if (fd_from < 0)
    return -1;

  fd_to = open(to, O_WRONLY | O_CREAT | O_EXCL, 0666);
  if (fd_to < 0)
    goto out_error;

  while (nread = read(fd_from, buf, sizeof buf), nread > 0) {
    char * out_ptr = buf;
    ssize_t nwritten;

    do {
      nwritten = write(fd_to, out_ptr, nread);

      if (nwritten >= 0) {
        nread -= nwritten;
        out_ptr += nwritten;
      } else if (errno != EINTR) {
        goto out_error;
      }
    } while (nread > 0);
  }

  if (nread == 0) {
    if (close(fd_to) < 0) {
      fd_to = -1;
      goto out_error;
    }
    close(fd_from);

    /* Success! */
    return 0;
  }

  out_error:
    saved_errno = errno;

  close(fd_from);
  if (fd_to >= 0)
    close(fd_to);

  errno = saved_errno;
  return -1;
}

/* #include <sys/exec.h> */
/* for sun4 and Dec alpha */

void wait_then_exit();

/* A function used to print a standard error message, along with specifying
where they made there mistake on the command line. */

usagerr(err)

int err;

{
  if (err == 1)
    fprintf(stderr, "\nCheck your t switch\n\n");
  else if (err == 2)
    fprintf(stderr, "\nCheck your n switch\n\n");
  else if (err == 3)
    fprintf(stderr, "\nCheck your s switch\n\n");
  else if (err == 4)
    fprintf(stderr, "\nCheck your c switch\n\n");
  else if (err == 6)
    fprintf(stderr, "\nCheck your p switch\n\n");
  else if (err == 5)
    fprintf(stderr, "\nIllegal switch\n\n");

  fprintf(stderr, "Usage: psx [-i infile] [-n numseq] [-s skipnum] [-t total] [-p procs]\n");
  fprintf(stderr, " [-d dir] [-C 'args'] -c executable\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "-i infile Name of file with seqs, default = (stdin)\n");
  fprintf(stderr, "-n numseq Number of seqs per iteration of executable, num_seq > 0 (1)\n");
  fprintf(stderr, "-a Appends environment variables to the call to executable\n");
  fprintf(stderr, "-s skipnum Skips skipnum num of seqs before sending to executable >= 0 (0)\n");
  fprintf(stderr, "-t total Total number of sequences to use from infile, total > 0 (ALL)\n");
  fprintf(stderr, "-p procs number of processes to spawn, 0 < procs <= %d (1)\n", MAX_PROCS);
  fprintf(stderr, "-d dir root name of directories to spawn processes from, (tmp)\n");
  fprintf(stderr, "-C 'args' command line arguments to pass to executable\n");
  fprintf(stderr, "-c executable Name of the program or script which sx will call\n");
  exit(2);
}

pid_t
handle_wait()

{
  pid_t child_pid;
  int child_stat, child_exit;

  do {
    if ((child_pid = wait( & child_stat)) < 0) {
      if (errno == ECHILD) {
        return (child_pid);
      }
      if (errno == EINTR) {
        fprintf(stderr, "WARNING: psx: wait failed due to a signal interrupt\n");
      } else {
        fprintf(stderr, "ERROR: psx: wait failed for unknown reason - exiting!!\n");
        exit(4);
      }
    } else if (WIFEXITED(child_stat)) {
      if (child_exit = WEXITSTATUS(child_stat)) {
        fprintf(stderr, " ERROR: psx: child terminated with nonzero exit status (%d) - exiting!!\n", child_exit);
        wait_then_exit(child_exit);
      }
      return (child_pid);
    } else if (WIFSIGNALED(child_stat)) {
      fprintf(stderr, "ERROR: psx: child terminated by a signal (%d) - exiting!!\n", WTERMSIG(child_stat));
      wait_then_exit(4);
    } else if (WIFSTOPPED(child_stat)) {
      fprintf(stderr, "ERROR: psx: child stopped by a signal (%d) - exiting!!\n", WSTOPSIG(child_stat));
      wait_then_exit(4);
    }
    /*else if (WIFCONTINUED(child_stat)) {
fprintf(stderr, "ERROR: psx: child has been continued - exiting!!\n");
wait_then_exit(4);
} */
    else {
      fprintf(stderr, "ERROR: psx: wait returned unknown status - exiting!!\n");
      exit(4);
    }
  } while (True);
}

void
wait_then_exit(exit_status)

int exit_status;

{
  int status;

  fprintf(stderr, "WAITING for all children to finish before exiting!\n");
  while (handle_wait( & status) >= (pid_t) 0);
  exit(exit_status);
}

main(argc, argv)

int argc;
char ** argv;

{
  extern void init_cur_dir();
  extern char * get_cur_dir();
  extern void change_cur_dir();

  /* Used to hold the -c switch for later use */
  char * executable = NULL;

  /* Used to hold the -d switch for later use */
  char * execdir = "tmp";

  /* Used to hold the -C switch for later use */
  char * cmdargs = "";

  /* Arrays used for holding the expressions to be put in the enviornment */
  char newfilename[1024], filename[512], fileroot[256], filestring[256], numstring[128];
  char skipstring[128], totalstring[128], cntstring[128], endstr[128];
  /*char ls_string[512], cat_string[512];*/

  /* Used for holding the command line sent to exec */
  char exec_line[5120];

  /* Name of the file with the -i switch */
  char * infile_name;

  /* Used with the wait() function after an exec */
  int status = True;

  /* Counter */
  int i;

  /* Number after the -n switch */
  int num = 1;

  /* Counter to do number requested with -t switch */
  int total = 1;

  /* Number after -s switch */
  int skipnum = 0, skipnum_env;

  /* Number after -t switch */
  int tlimit = 0;

  /* Number after -p switch and current # of procs */
  int numprocs = 1, curprocs = 0;

  /* Counter for iterations */
  int cnt = 1;

  /* The temporary file which hold the groups of sequences, SXFILE */
  FILE * fd = NULL;

  /* Holds stdin if not used by -i switch */
  FILE * in = NULL;

  /* Holds each line of read-in sequence from file */
  char s[IN_LINE_SIZE];

  /* Signifies if last sequence or not */
  int finish = False;

  /* Child process ID number */
  pid_t child;

  /* current working directory */
  char * curdir;

  /* Goes through all of the commands in argv[] */
  while (argc > 1 && argv[1][0] == '-') {
    argv++;
    argc--;
    switch (argv[0][1]) {
    case 'p':
      i = sscanf(argv[1], "%d", & numprocs);
      argv++;
      argc--;
      if ((numprocs <= 0) || (numprocs > MAX_PROCS)) {
        usagerr(6);
      }
      break;
    case 't':
      i = sscanf(argv[1], "%d", & tlimit);
      argv++;
      argc--;
      if (tlimit <= 0) {
        usagerr(1);
      }
      break;
    case 's':
      i = sscanf(argv[1], "%d", & skipnum);
      argv++;
      argc--;
      if (skipnum < 0) {
        usagerr(3);
      }
      break;
    case 'n':
      i = sscanf(argv[1], "%d", & num);
      argv++;
      argc--;
      if (num <= 0) {
        usagerr(2);
      }
      break;
    case 'i':
      infile_name = argv[1]; in = fopen(argv[1], "r");
      if ( in == NULL) {
        fprintf(stderr, "Error opening file %s for input \n", infile_name);
        exit(2);
      }
      argv++;
      argc--;
      break;
    case 'd':
      execdir = argv[1];
      argv++;
      argc--;
      break;
    case 'C':
      cmdargs = argv[1];
      argv++;
      argc--;
      break;
    case 'c':
      executable = argv[1];
      argv++;
      argc--;
      break;

      /* If a switch other than the ones above were present default is run */
    default:
      usagerr(5);
      break;
    }
  }

  /* Gives error and exits if no -c was on command line */
  if (executable == NULL) {
    usagerr(4);
  }

  if ( in == NULL)
    in = stdin;
  else
    freopen("/dev/null", "r", stdin);

  /* Checks to see if input file is empty */
  if (fgets(s, (IN_LINE_SIZE - 1), in ) == NULL) {
    fprintf(stderr, "The input file appears to be empty.\n");
    exit(2);
  }

  /* Checks the first character of the file, has to be either > or # */
  if ((strncmp(s, ">", 1) != 0) && (strncmp(s, "#", 1) != 0)) {
    fprintf(stderr, "The file is not in the correct FASTA format \n");
    exit(2);
  }

  /* Set up the format of SXFILE */
  sprintf(fileroot, "%s.%d.sx_file_", getenv("USER"), getpid());
  sprintf(filename, "%s%d", fileroot, cnt);
  /*sprintf(ls_string,"ls -lg %s",filename);*/
  /*sprintf(cat_string,"cat %s",filename);*/

  /* Sets up a temp file also as an environment var for the repeated batch */
  /* Sets up SXNUM and SXOFF with their respective numbers */
  sprintf(numstring, "SXNUM=%d", num);
  sprintf(skipstring, "SXSKIP=%d", skipnum);
  sprintf(totalstring, "SXTOTAL=%d", tlimit);
  sprintf(endstr, "SXOFF=%d", finish);
  sprintf(cntstring, "SXCNT=%d", cnt);
  sprintf(filestring, "SXFILE=%s", filename);

  /* Puts them into the environment */
  if ((putenv(filestring) != 0) || (putenv(numstring) != 0) ||
    (putenv(skipstring) != 0) || (putenv(totalstring) != 0) ||
    (putenv(endstr) != 0) || (putenv(cntstring) != 0)) {
    fprintf(stderr, "ERROR: putenv failed!\n");
    exit(4);
  }

  /* Opens SXFILE for writing */
  fd = fopen(filename, "w");
  if (fd == NULL) {
    fprintf(stderr, "The file - %s - is not writable.\n", filename);
    exit(3);
  }

  /* initialize working directories for procs */
  init_cur_dir(execdir, numprocs);

  skipnum_env = skipnum;
  while (skipnum > 0) {
    if (fgets(s, (IN_LINE_SIZE - 1), in ) == NULL)
      exit(0);
    if ((strncmp(s, ">", 1) == 0) || (strncmp(s, "#", 1) == 0))
      skipnum--;
  }
  i = 0;
  curprocs = 0;
  do {
    if (fputs(s, fd) == EOF) {
      fprintf(stderr, "ERROR: fputs failed for %s!\n", filename);
      wait_then_exit(4);
    }
    /*system(ls_string);*/
    /*system(cat_string);*/
    if (fgets(s, (IN_LINE_SIZE - 1), in ) == NULL)
      finish = True;
    if (finish || (strncmp(s, ">", 1) == 0) || (strncmp(s, "#", 1) == 0)) {
      /* Continues until tlimit has been reached */
      if (total == tlimit) {
        finish = True;
      }
      i++;
      total++;
      if (finish || (i >= num)) {
        if (fclose(fd) == EOF) {
          fprintf(stderr, "WARNING: fclose failed for %s! (errno = %d)\n", filename, errno);
        }
        /*system(ls_string);*/
        /*system(cat_string);*/
        if (curprocs < numprocs) {
          curdir = get_cur_dir(0);
          curprocs++;
        } else {
          if ((child = handle_wait( & status)) < (pid_t) 0) {
            fprintf(stderr, "ERROR: wait failed! (errno = %d)\n", errno);
            exit(4);
          }
          curdir = get_cur_dir(child);
        }
        if (finish) {
          fprintf(stderr, "WAITING for all children to finish before starting last child!\n");
          while (handle_wait( & status) >= (pid_t) 0);
        }
        if (curdir == (char * ) NULL) {
          fprintf(stderr, "ERROR: get_cur_dir failed!\n");
          wait_then_exit(4);
        }
        /*system(ls_string);*/
        /*system(cat_string);*/
        sprintf(newfilename, "%s/%s", curdir, filename);
        if ((unlink(newfilename) < 0) && (errno != ENOENT)) {
          fprintf(stderr, "ERROR: unlink failed for %s ! (errno = %d)\n", newfilename, errno);
          wait_then_exit(4);
        }
        if (copy_file(filename, newfilename) < 0) //changed by Ray to copy_file
        {
          fprintf(stderr, "ERROR: link failed for %s to %s! (errno = %d)\n", filename, newfilename, errno);
          wait_then_exit(4);
        }
        if (unlink(filename) < 0) {
          fprintf(stderr, "ERROR: unlink failed for %s ! (errno = %d)\n", filename, errno);
          wait_then_exit(4);
        }
        /* Sets up SXNUM and SXOFF with their respective numbers */
        sprintf(numstring, "SXNUM=%d", i);
        sprintf(endstr, "SXOFF=%d", finish);
        if ((child = fork()) == (pid_t) 0) {
          /* Spawns off shell with command after -c flag with
           * the environment variables appended to the command line
           */
          if (chdir(curdir) < 0) {
            fprintf(stderr, "ERROR: chdir failed! (errno = %d)\n", errno);
            wait_then_exit(4);
          }
          freopen("/dev/null", "r", stdin);
          sprintf(exec_line, "%s %s %d %d %d %d %d %s", executable, filename, i, cnt, finish, skipnum_env, tlimit, cmdargs);
          execl("/bin/sh", "sh", "-c", exec_line, NULL);
        } else if (child < (pid_t) 0) {
          fprintf(stderr, "ERROR: fork failed! (errno = %d)\n", errno);
          wait_then_exit(4);
        }
        change_cur_dir(child, curdir);
        cnt++;
        sprintf(cntstring, "SXCNT=%d", cnt);
        if (!finish) {
          sprintf(filename, "%s%d", fileroot, cnt);
          fd = fopen(filename, "w");
          if (fd == NULL) {
            fprintf(stderr, "The file - %s - is not rewritable. (errno = %d)\n", filename, errno);
            wait_then_exit(3);
          }
          /*system(ls_string);*/
          /*system(cat_string);*/
        } else {
          fprintf(stderr, "WAITING for all children to finish!\n");
          while (handle_wait( & status) >= (pid_t) 0);
        }
        i = 0;
      }
    }
  }
  while (!finish);

  exit(0);
}
