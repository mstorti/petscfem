#!/usr/bin/perl
use File::Copy;


my ($sub_dir,$node_list) = @ARGV;

if (-e $sub_dir) { system "rrm -rf $sub_dir"; }
mkdir $sub_dir;

open PROCTABLE,"./proctable";
my %proctable = ();
while (<PROCTABLE>) {
    next if /^\s*\#/;
    die "bad line: $_\n" unless /^\s*(node\d*)\s*(\S*)(|\s*server)/;
    my ($node,$line,$w) = ($1,$_,$2);
    chomp $line;
    chomp $w;
    $proctable{$node} = [$line,$w];
}

=cut
while (($k,$v) = each(%proctable)) {
    print "$k -> [$v->[0],$v->[1]]\n";
}
=cut

open SDPROCT,">$sub_dir/proctable";
print SDPROCT "node1 1.3 server\n";
my @node_list = split " ",$node_list;
for my $node (@node_list) {
    my $v =$proctable{$node};
    die "can't find $node\n" unless defined $v;
#    print "$node -> [$v->[0],$v->[1]]\n";
    print SDPROCT "$v->[0]\n";
}

copy "cubcav.epl",$sub_dir;
copy "Makefile.sub-dir","$sub_dir/Makefile";
