#!/usr/bin/perl

use tree;

$tree = new tree('Root');
$node = $tree->add_son('a');
$node = $node->add_brother('b');
$node = $node->add_brother('c');
$node = $node->add_son('c1');
$node = $node->add_brother('c2');
$node = $node->add_brother('c3');
$node = $node->add_brother('c4');
$node = $node->add_brother('c5');
$node = $node->father();
$node = $node->add_brother('d');
$node = $node->add_brother('e');
$node = $node->add_brother('f');

print "Pre order:\n";
$tree->pre_print();
print "Post order:\n";
$tree->post_print();

$p = sub { 
    my ($v,$r) = @_;
    $r .= " - $v";
    return $r;
};

$r = $tree->pre_order($p,"");
print "$r\n";

