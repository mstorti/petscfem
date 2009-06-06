require "$ENV{'PETSCFEM_DIR'}/test/eperlini.pl"; # Initializes ePerl 

# Synchronize the time step with external frequency
$piecewise_fun = 1 if $spline_fun;
$N = 16 unless $N;   # number of time steps per period
$omega=2.*$PI*1000.; # pulsation
$T=2*$PI/$omega;     # period
$Dt = $T/$N;         # we put $N time steps per period
$Nperiods = 1 unless $Nperiods; # Numbeer of periods to be run
$alpha=1 unless $alpha; # trapezoidal rule parameter

sub heredoc {
    my $txt = shift();
    my $name = shift();
    $name = "tmp_file_$$.tmp" unless $name;
    open TMP,">$name";
    print TMP $txt;
    close TMP;
#    print $name;
}

if ($piecewise_fun) {
heredoc(<<EOF,"tmp_file_some.tmp");
1
11
EOF
} elsif ($spline_periodic_fun) {
heredoc(<<EOF,"tmp_file_some.tmp");
1
EOF
} else {
heredoc(<<EOF,"tmp_file_some.tmp");
11
EOF
}

1;
