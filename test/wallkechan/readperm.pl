#!/usr/bin/perl -n

print "$1 $2\n" if /row (\d*): \((\d*) -> /;
