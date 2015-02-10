package Pasa_conf;

use strict;
no strict qw (refs);
use warnings;
use ConfigFileReader;
use Carp;

my %conf;

BEGIN {
    unless ($ENV{PASAHOME}) {
        die "ERROR, env var PASAHOME not set to base PASA installation directory.\n";
    }
    
    my $confFile = "$ENV{PASAHOME}/pasa_conf/conf.txt";
    unless (-s $confFile) {
        die "ERROR, cannot find conf file $confFile\n";
    }
    
    %conf = &ConfigFileReader::readConfig($confFile);
}

sub getParam {
    my $param = shift;
    
    my $value = $conf{$param};
    
    if ($value) {
        $value =~ s/__PASAHOME__/$ENV{PASAHOME}/g;
    }
    return ($value);
}

sub dumpParams {
    foreach my $param (keys %conf) {
        print "$param = $conf{$param}\n";
    }
}



sub call_hook {
    my ($hook_name, @params) = @_;
    
    my $hook_method = &_get_hook($hook_name);

    my $extra_param = &_get_hook_extra_param($hook_method);

    if ($extra_param) {
        unshift (@params, $extra_param);
    }
    
    ## call the method:
    my %conf_copy = %conf;
    
    my @results;
    eval {
        print "Calling hook: $hook_method(" . join (",", \%conf_copy, @params) . "\n";
        @results = &$hook_method(\%conf_copy, @params);
    };

    if ($@) {
        confess "Error, calling hook method [$hook_method] failed: $@";
    }
    
    if (@results) {
        if (wantarray) {
            return @results;
        }
        else {
            return ($results[0]);
        }
    }
    else {
        return;
    }
    
}


#####################################################################
## private methods:
#####################################################################

sub _get_hook {
    my ($hook_name) = @_;
    
    ## ensure hook perl libs are loaded:
    &_load_hook_perl_libs();
    
    my $qualified_method_name = &getParam($hook_name);
    
    my @components = split (/::/, $qualified_method_name);
    
    my $method_name = pop @components;
    my $packagename = join ("::", @components);  ## could be a deeper package name with lots of :: ::'s.
    
    unless ($packagename && $method_name) {
        croak "Error, couldn't dissociate the packagename from the method name of supposed qualified hook method name [$qualified_method_name";
    }

    ## dynamically load the package.
    ## locate it in @INC
    
    my $found = 0;
    foreach my $lib (@INC) {
        my $path = "$lib/$packagename.pm";
        $path =~ s|::|/|g;
        if (-e $path) {
            # got it.  Lock and load
            require $path;
            import $packagename;
            print "Required and Imported $path\n";
            $found=1;
            last;
        }
    }
    
    unless ($found) {
        confess "Error, couldn't resolve path for $packagename\n";
    }

    return ($qualified_method_name);
}


sub _get_hook_extra_param {
    my ($hook_method_name) = @_;
    
    my $key_name = "$hook_method_name" . "~EXTRA_PARAM";
    
    return (&getParam($key_name));
}



####
sub _load_hook_perl_libs {
    
    my $lib_list = &getParam("HOOK_PERL_LIBS");
    unless ($lib_list) { 
        return;
    }

    my @libs = split (/,/, $lib_list);
    # remove whitespace
    foreach my $lib (@libs) {
        $lib =~ s/\s//g;
        print "examining lib: $lib\n";
        unless (grep { $_ eq $lib } @INC) {
            print "adding lib: $lib\n";
            push (@INC, $lib);
        }
    }
}




1;

