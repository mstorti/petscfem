# node = father - list of sons - user data
package tree;

$Lambda = undef;

sub left_son {
    $node = shift();
    return $node->[0];
}

sub right_brother {
    my $node = shift();
    return $node->[1];
}    

sub father {
    my $node = shift();
    return $node->[2];
}

sub data {
    my $node = shift();
    return $node->[3];
}

sub node_label {
    my $node = shift();
    return $node->[3];
}

sub add_son {
    my $node = shift();
    my $newson_data = shift();
    my $son = $node->[0];
    my $new_node = [$Lambda,$son,$node,$newson_data];
    $node->[0] = $new_node;
}

sub add_brother {
    my $node = shift();
    my $data = shift();
    my $brother = $node->[1];
    my $new_node = [$Lambda,$brother,$node->[2],$data];
    $node->[1] = $new_node;
}

sub post_print {
    my $node = shift();
    return unless defined $node;
    print data($node),"\n";
    post_print(left_son($node));
    post_print(right_brother($node));
}

sub pre_print {
    my $node = shift();
    return unless defined $node;
    pre_print(left_son($node));
    print data($node),"\n";
    pre_print(right_brother($node));
}

sub printn {
    my $node = shift();
    print "<",join("><",@{$node}),">";
}

$tree=[$Lambda,$Lambda,$Lambda,$Lambda];

$node = add_son($tree,'a');
$node = add_brother($node,'b');
$node = add_brother($node,'c');
$node = add_son($node,'c1');
$node = add_brother($node,'c2');
$node = add_brother($node,'c3');
$node = add_brother($node,'c4');
$node = add_brother($node,'c5');
$node = father($node);
$node = add_brother($node,'d');
$node = add_brother($node,'e');
$node = add_brother($node,'f');

pre_print($tree);
