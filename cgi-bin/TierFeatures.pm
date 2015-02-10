package TierFeatures;
use strict;
use warnings;

sub new {
    my $packagename = shift;
    
    my $self = { tiers => [ [] ], # give it one empty tier to start
             };

    bless ($self, $packagename);

    return ($self);
}


sub tier_features {
    my $self = shift;
    my @features = @_;
            
    ## start at first tier:
    foreach my $feature (@features) {
        my ($feat_lend, $feat_rend) = ($feature->{lend}, $feature->{rend});
        my @tiers = @{$self->{tiers}};
        
        my $tiered_feature_flag = 0;
        
      tiers:
        foreach my $tier (@tiers) {
            my @tiered_feats = @$tier;
          feats:
            foreach my $feat (@tiered_feats) {
                my ($lend, $rend) = ($feat->{lend}, $feat->{rend});
                # check for overlap

                if ($lend <= $feat_rend && $rend >= $feat_lend) {
                    # got overlap
                    next tiers;
                }
            }
         
            # if got here, no overlap in current tier.  Just add it:
            push (@$tier, $feature);
            $tiered_feature_flag = 1;
            last tiers;
        }
        
        unless ($tiered_feature_flag) {
            # no current tier can accommodate it.  Add another tier with this element
            push (@{$self->{tiers}}, [$feature]);
        }
    }
 
    ## return tiers:

    return (@{$self->{tiers}});
    
}







package TierFeatures::Feature;
use strict;
use warnings;


sub new {
    my $packagename = shift;
    my ($lend, $rend, $feat_id) = @_;
    
    if ($lend > $rend) {
        ($lend, $rend) = ($rend, $lend);
    }

    my $self = { lend => $lend,
                 rend => $rend,
                 feat_id => $feat_id,
             };

    bless ($self, $packagename);
    return ($self);
}
                


1; #EOM
