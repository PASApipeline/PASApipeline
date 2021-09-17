#!/usr/bin/env perl

package main;
our ($DB_SEE, $DEBUG);


package DB_connect;

require 5.6.0;
require Exporter;
use Carp;
use strict;
use Data::Dumper;

use threads;
use threads::shared;
my $LOCKVAR :shared;


our @ISA = qw(Exporter);

## export the following for general use
our @EXPORT = qw (do_sql_2D connect_to_db RunMod first_result_sql very_first_result_sql get_last_insert_id );

my $QUERYFAIL;

my $CONNECTED_FLAG = 0; # set to 1 on first db connection. Mostly for logging info.

############### DATABASE CONNECTIVITY ################################
####

sub configure_db_driver {
    my ($db) = @_;
    
    if ($db =~ /^\.*\//) {
        # have fully qualified path:
        $ENV{DBI_DRIVER} = 'SQLite';
        print STDERR "-connecting to SQLite db: $db\n" unless $CONNECTED_FLAG;
    }
    else {
        $ENV{DBI_DRIVER} = 'mysql';
        print STDERR "-connecting to MySQL db: $db\n" unless $CONNECTED_FLAG;
    }
    
    
    return;
}

sub connect_to_db {
    my ($server, $db, $username, $password) = @_;
    #print Dumper([@_]);
    
    unless (defined $db) {
        confess "Error, require parameters: server, db, username, password\n"
            . " If you're using sqlite, you can specify dummy values for server, username, and password\n"
            . " but the 'db' parameter must be specified as a fully qualified path. \n"
            . " ie.   /path/to/my/sqlite.db\n\n";
    }
    
      
    unless($ENV{DBI_DRIVER} && $ENV{DBI_DRIVER} =~ /^(SQLite|mysql)$/) {
        &configure_db_driver($db);
    }
        
    
    if ($ENV{DBI_DRIVER} eq 'mysql') {
        if (! (defined ($server) && defined($db) && defined($username) && defined($password)) ) {
            confess "Error, need all method parameters (server, db, username, password) for mysql ";
        }
    }
    
    my $dbh = DBI->connect("dbi::database=$db;host=$server", $username, $password);
    unless (ref $dbh) {
        croak "Cannot connect to $db: $DBI::errstr";
    }
    $dbh->{RaiseError} = 1; #turn on raise error.  Must use exception handling now.
    

    
    ## add attributes so can reconnect later in case mysql server goes away.
    
    my $dbproc = new DB_connect(); ## temporary fix to deal with lost connections
    $dbproc->{dbh} = $dbh;
    
    if ($dbh->{Driver}->{Name} ne 'SQLite') {
        $dbproc->{__server} = $server;
        $dbproc->{__db} = $db;
        $dbproc->{__username} = $username;
        $dbproc->{__password} = $password;
    }

    unless ($CONNECTED_FLAG) {
        $CONNECTED_FLAG = 1;
    }
    
    return($dbproc);
}


####
sub new {
    my $packagename = shift;
    my $self = {};
    bless ($self, $packagename);
    return ($self);
}


####
sub disconnect {
    my $self = shift;
    $self->{dbh}->disconnect;
}

####
sub get_dbh {
    my $self = shift;
    return ($self->{dbh});
}


####
sub reconnect_to_server {
    my ($dbproc) = @_;

    return $dbproc if $dbproc->{dbh}->{Driver}->{Name} eq 'SQLite'; # shouldn't be needed for SQLite

    my $new_dbh;
    
    do {
        eval {
            $new_dbh = &connect_to_db($dbproc->{__server},  $dbproc->{__db}, $dbproc->{__username}, $dbproc->{__password});
        };
        if ($@) {
            $new_dbh = undef;
            sleep(30);
        }
    } until ($new_dbh);
    
    $dbproc->{dbh} = $new_dbh->{dbh};
    
    return ($dbproc);
}



## return results in 2-Dimensional array.
sub do_sql_2D {
    my ($dbproc,$query, @values) = @_;
    my ($statementHandle,@x,@results);
    my ($i,$result,@row);

    
    
    ## Use $QUERYFAIL Global variable to detect query failures.
    $QUERYFAIL = 0; #initialize
    print "QUERY: $query\tVALUES: @values\n" if($::DEBUG||$::DB_SEE);
    
    eval {
        $statementHandle = $dbproc->{dbh}->prepare($query);
    };
    if ($@) {
        confess "Error $@";
    }
    
    if ( !defined $statementHandle) {
        print "Cannot prepare statement: $DBI::errstr\nQUERY: $query\tVALUES: @values\n";
        $QUERYFAIL = 1;
    } else {
        
        # Keep trying to query thru deadlocks:
        do {
            $QUERYFAIL = 0; #initialize
            eval {
                {
                    lock $LOCKVAR;
                    
                    $statementHandle->execute(@values);
                    while ( @row = $statementHandle->fetchrow_array() ) {
                        push(@results,[@row]);
                    }
                }
            };
            ## exception handling code:
            if ($@) {
                
                ## check for mysql gone away:
                if ($DBI::errstr =~ /server has gone away|Lost connection/) {
                    ## reestablish connection and try again:
                    $dbproc = &reconnect_to_server($dbproc);
                    return (&do_sql_2D($dbproc, $query, @values));
                }
                else {
                    print STDERR "\n\n====\nFailed query: <$query>\tvalues: @values\nErrors: $DBI::errstr\n====\n";
                    $QUERYFAIL = 1;
                }
            }
            
        } while ($statementHandle->errstr() =~ /deadlock/);
        #release the statement handle resources
        $statementHandle->finish;
    }
    if ($QUERYFAIL) {
        confess "Failed query: <$query>\tvalues: @values\nErrors: $DBI::errstr\n";
    }
    return(@results);
}


sub RunMod {
    my ($dbproc,$query, @values) = @_;
    my ($result);

    if($::DEBUG||$::DB_SEE) {print "QUERY: $query\tVALUES: @values\n";}
    if($::DEBUG) {
        $result = "NOT READY";
    } else {
        eval {
            {
                lock $LOCKVAR;
                $dbproc->{dbh}->do($query, undef, @values);
            }
        };
        if ($@) { #error occurred
            
            ## check for mysql gone away:
            if ($DBI::errstr =~ /server has gone away|Lost connection/) {
                ## reestablish connection and try again:
                &reconnect_to_server($dbproc);
                return (&RunMod($dbproc, $query, @values));
            }
            else {
                confess "failed query: <$query>\tvalues: @values\nErrors: $DBI::errstr\n";
            }
        }
    }

    return;
}


sub first_result_sql {
    my ($dbproc, $query, @values) = @_;
    my @results = &do_sql_2D ($dbproc, $query, @values);
    return ($results[0]);
}

sub very_first_result_sql {
    my ($dbproc, $query, @values) = @_;
    my @results = &do_sql_2D ($dbproc, $query, @values);
    if ($results[0]) {
        return ($results[0]->[0]);
    } else {
        return (undef());
    }
}

sub get_last_insert_id {
    my ($dbproc) = @_;
    my $query = ($dbproc->{dbh}->{Driver}->{Name} eq 'SQLite') ? 'select last_insert_rowid()' : "select LAST_INSERT_ID()";
    return (&very_first_result_sql($dbproc, $query));
}

sub delete_table {
    my ($dbproc, $table) = @_;

    my $query = ($dbproc->{dbh}->{Driver}->{Name} eq 'SQLite') ? "DELETE FROM $table" : "TRUNCATE TABLE $table"; 
    &RunMod($dbproc, $query);

    return;
}

1; #EOM
