#!/usr/bin/env perl

sub usage {
    my $message = "@_";
    die "
Usage: ind2pheno.perl infile outfile

Required Arguments:
  infile  : input .ind file
  outfile : output .pheno file

$message

";
}

unless (@ARGV == 2) {usage("OOPS unexpected number of arguments")}

$in = $ARGV[0]; # .ind file
$out = $ARGV[1]; # .pheno file

open(IN,$in) || die("Cannot open file: $in");
open(OUT,">$out") || die("Cannot open file: $out");

while($line = <IN>)
{
  if($line =~ /Case/) { print OUT ("1"); $case=1; }
  elsif($line =~ /Control/) { print OUT ("0"); $control=1; }
  else { print OUT ("9"); }
}
print OUT ("\n");
unless($case) { print("WARNING: no cases\n"); }
unless($control) { print("WARNING: no controls\n"); }

