#!/usr/bin/env perl

package main;
our ($DB_SEE, $DEBUG);


package Mysql_connect;

require 5.6.0;
require Exporter;
use Carp;
use strict;
use Data::Dumper;

our @ISA = qw(Exporter);

## export the following for general use
our @EXPORT = qw (do_sql_2D connect_to_db RunMod first_result_sql very_first_result_sql get_last_insert_id );

my $QUERYFAIL;

############### DATABASE CONNECTIVITY ################################
####
sub connect_to_db {
    my ($server, $db, $username, $password) = @_;
    
    unless (defined ($server) && defined($db) && defined($username) && defined($password)) {
        confess "Error, need all method parameters (server, db, username, password) ";
    }
    
    my $dbh = DBI->connect("DBI:mysql:$db:$server", $username, $password);
    unless (ref $dbh) {
        croak "Cannot connect to $server: $DBI::errstr";
    }
    $dbh->{RaiseError} = 1; #turn on raise error.  Must use exception handling now.
    

    
    ## add attributes so can reconnect later in case mysql server goes away.
    
    my $dbproc = new Mysql_connect(); ## temporary fix to deal with lost connections
    $dbproc->{dbh} = $dbh;
    
    $dbproc->{__server} = $server;
    $dbproc->{__db} = $db;
    $dbproc->{__username} = $username;
    $dbproc->{__password} = $password;
    
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
                $statementHandle->execute(@values);
                while ( @row = $statementHandle->fetchrow_array() ) {
                    push(@results,[@row]);
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
                    print STDERR "failed query: <$query>\tvalues: @values\nErrors: $DBI::errstr\n";
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
            $dbproc->{dbh}->do($query, undef, @values);
        };
        if ($@) { #error occurred
            
            ## check for mysql gone away:
            if ($DBI::errstr =~ /server has gone away/) {
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
    my $query = "select LAST_INSERT_ID()";
    return (&very_first_result_sql($dbproc, $query));
}



1; #EOM
