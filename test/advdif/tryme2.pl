sub pp {
    my $a = shift();
    my $b = shift();
    print "a: $a, b: $b\n";
}

use subs 'pp';
