#!/usr/bin/perl
# usage: verif.pl <dir> <file>
# Verifies that all headers in subdirectory rooted at
# `dir' are present in doc++ file `file'

$docxx = shift();

$dir = shift();
@files = split "\0",`find $dir -name '*.h' -print0`;
for $file (@files) {
#      if ($file =~ m|/([^/]*)$|) {
#  	$filedb{$file} = 1;
#    	print "file $file\n";
#      } else {
#    	print "nof found file in \"$file\"\n";
#      }
    $file =~ s|^../.||;
    print "file: \"$file\"\n";
    $filedb{$file} = 1;
}
print scalar @files," files in subdirectory tree\n";

open DOCXX,$docxx;
while (<DOCXX>) {
    if (m|../../(.*\.h)|) {
	my $f = "./$1";
	$docxx{$f} = 1;
	print "file: \"$f\"\n";
	if ($filedb{$f}) {
	    $filedb{$f} = 0;
	}
    }
}
print scalar (keys(%filedb))," entries in doc++ file\n";

print "Not in DOCXX:\n";
while (($k,$v) = each %filedb) {
    if ($v) { print "//\@Include: ../.$k\n"; }
}

# print "Files: \n",join ">\n<",@files;
