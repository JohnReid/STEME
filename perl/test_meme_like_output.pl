#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::Matrix::PSM::IO;

my $psmIO = new Bio::Matrix::PSM::IO(	
	-file   => 'output/meme.out',
	-format => 'meme');

while (my $psm = $psmIO -> next_psm) {
   my $instances = $psm -> instances;
   foreach my $instance (@{ $instances }) {
	   my $start = $instance -> start;
	   my $score = $instance -> score;
	   printf("%5d : %.1e\n", $start, $score);
   }
}

