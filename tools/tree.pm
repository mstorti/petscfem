# node = [ left_son , right_brother , father, label ]
package tree;

$Lambda = [];
bless $Lambda,"tree";

sub new {
    my $class = shift();
#    my $node = ($#_ == 3 ? [@_] : [$Lambda,$Lambda,$Lambda,undef]);
    my $node = ($#_ == 3 ? [@_] : $#_==0 ? [$Lambda,$Lambda,$Lambda,$_[0]] : $Lambda);
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
    my $son = $node->left_son();
    while ($son!=$Lambda) { 
	$son->post_print();
	$son = $son->right_brother();
    }
    print $node->data(),"\n";
}

sub pre_print {
    my $node = shift();
    return if $node==$Lambda;
    print $node->data(),"\n";
    my $son = $node->left_son();
    while ($son!=$Lambda) { 
	$son->pre_print();
	$son = $son->right_brother();
    }
}

sub pre_order {
    my $node = shift();
    my $sub = shift();
    my $retval = shift();
    return if $node==$Lambda;
    $retval = &{$sub}($node->data(),$retval);
#    print $node->data(),"\n";
    my $son = $node->left_son();
    while ($son!=$Lambda) { 
       $retval = $son->pre_order($sub,$retval);
	$son = $son->right_brother();
    }
    return $retval;
}

sub printn {
    my $node = shift();
    print "<",join("><",@{$node}),">";
}

1;

