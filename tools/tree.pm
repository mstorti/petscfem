# node = father - list of sons - user data
package tree;

$Lambda = [];
bless $Lambda,"tree";

sub new {
    my $class = shift();
    my $node = ($#_ == 3 ? [@_] : [$Lambda,$Lambda,$Lambda,undef]);
    bless $node, $class;
    return $node;
}

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
    my $new_node = new tree($Lambda,$son,$node,$newson_data);
    $node->[0] = $new_node;
}

sub add_brother {
    my $node = shift();
    my $data = shift();
    my $brother = $node->[1];
    my $new_node = new tree($Lambda,$brother,$node->[2],$data);
    $node->[1] = $new_node;
}

sub post_print {
    my $node = shift();
    return if $node==$Lambda;
    print $node->data(),"\n";
    left_son($node)->post_print();
    right_brother($node)->post_print();
}

sub pre_print {
    my $node = shift();
    return if $node==$Lambda;
    $node->left_son()->pre_print();
#    $node->print_data(),"\n";
    print $node->data(),"\n";
    $node->right_brother()->pre_print->();
}

sub printn {
    my $node = shift();
    print "<",join("><",@{$node}),">";
}

1;

