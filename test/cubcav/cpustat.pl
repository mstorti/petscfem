#!/usr/bin/perl -w
use strict;

my $file = shift();

open OUT,$file;

sub search {
    my ($tag) = @_;
    while (my $line = <OUT>) {
	return ($1,$line) if $line =~ /.*$tag.*\[\d.*\s+(\S*)\]/;
    }
    return undef;
}

my $pred_res_prev = 0;
my $line;
while (1) {
    my ($pred_res,$pred_sol,$poi_res,$poi_sol,$prj_res,$prj_sol);

    ($pred_sol,$line) = search("-PREDICTOR- Before solving");
    last unless $pred_sol;

    ($poi_res,$line) = search("-POISSON- Before residual ");
    last unless $poi_res;
    ($poi_sol,$line) = search("-POISSON- Before solving");
    last unless $poi_sol;

    ($prj_res,$line) = search("-PROJECTION- Before residual ");
    last unless $prj_res;
    ($prj_sol,$line) = search("-PROJECTION- Before solving");
    last unless $prj_sol;

    ($pred_res,$line) = search("-PREDICTOR- Before residual ");
    last unless $pred_res;

    if (0) {
	printf "%.3f %.3f %.3f %.3f %.3f %.3f\n",
	$pred_res_prev,
	$pred_sol,
	$poi_res,
	$poi_sol,
	$prj_res,
	$prj_sol,
	$pred_res;
    } else {
	printf "%.3f %.3f %.3f %.3f %.3f %.3f\n",
	$pred_sol-$pred_res_prev,
	$poi_res-$pred_sol,
	$poi_sol-$poi_res,
	$prj_res-$poi_sol,
	$prj_sol-$prj_res,
	$pred_res-$prj_sol;
    }

    $pred_res_prev = $pred_res;
}

close OUT;
