package Pasa_CGI;

use strict;
use warnings;
use Carp;
use File::Basename;


####
sub get_CSS {
    my ($script_name) = @_;
    
    return ( &get_common_CSS() . &get_page_specific_CSS($script_name));
}


####
sub get_common_CSS { 
    my $common_CSS_file = "CSS/common.css";
    
    return (&_get_text_from_file($common_CSS_file));
}

####
sub get_page_specific_CSS {
    my ($script_name) = @_;
    my $file = basename($script_name);
    
    my $css_file = "CSS/$file.css";

    return (&_get_text_from_file($css_file));
}

####
sub _get_text_from_file {
    my ($file) = @_;
    
    unless (-e $file) {
        confess "Error, cannot locate the CSS file at: $file";
    }
    
    my $text = `cat $file`;
    return ($text);
}

1; #EOM


