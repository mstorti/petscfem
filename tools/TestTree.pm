use tree;
package TestTree;

sub new {
    my $class = shift();
    my $node = new tree(shift());
    $test_tree = {
	'node' => $node,
	'root' => $node,
	};
    bless $test_tree, $class;
    return $test_tree;
}

sub begin_section {
    my ($test_tree,$title) = @_;
    ## OK - not OK - can't open
    $test_tree->{'node'}->add_brother([$title,0,0,0]);
}

sub expect {
    my ($test_tree,$title) = @_;
    print "succes? [o/n/c (ok/not_ok/can't open)] >";
    my $ans = <STDIN>;
    chomp $ans;
    if ($ans eq 'o') {
	$test_tree->{'node'}->[1] +=1;
    } elsif ($ans eq 'n') {
	$test_tree->{'node'}->[2] +=1;
    } else {
	$test_tree->{'node'}->[3] +=1;
    }
}

sub end_section {
    
	
	
