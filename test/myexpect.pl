#__INSERT_LICENSE__

# usage: expect($file,$pattern);
#Advances $file until finds a little of "\n" delimited $pattern's. 
$DEBUG_EXPECT = 0;

sub expect {
    ($file,$descr,$pattern_list) = @_;
    open (SAL,$file) || do {print "can't open $file\n"; $cant_open++; return;};
    @sal=(<SAL>);
    close SAL;
    print "Testing: \"$descr\" on file \"$file\"...";
    print "\n" if $DEBUG_EXPECT;
    @pattern = split("\n",$pattern_list);
    $record=0;
    $inc=1;
    while ($pattern=shift @pattern) {
	if ($pattern =~ /^__REWIND__$/) {
	    $inc=1;
	    $record=0;
	    next;
	}
	if ($pattern =~ /^__BACKWARD__$/) {
	    $inc=-1;
	    next;
	}
	if ($pattern =~ /^__FORWARD__$/) {
	    $inc=+1;
	    next;
	}
	print "trying pattern: \"$pattern\"...  \n" if $DEBUG_EXPECT;
	while ($record<=$#sal) {
	    $_ = $sal[$record];
	    $record += $inc;
	    print "$record: $_" if $DEBUG_EXPECT;
	    do {chomp; 
		print "        -> found: \"$_\"\n" if $DEBUG_EXPECT; 
		goto NEXT;} if /$pattern/;
	}
	$notok++;
	print "not OK.\n        --->  Couldn't find \"$pattern\"\n";
	return 0;
      NEXT:;
    }
    $ok++;
    print "OK. \n";
}

sub final_check {
    $ok = 0 unless $ok;
    $notok = 0 unless $notok;
    $cant_open = 0 unless $cant_open;
    $total = $ok+$notok+$cant_open;
    print "Summary:  OK: $ok.  Not OK: $notok. ",
    "Couldn't open: $cant_open. Total: $total\n";
}

1;    
