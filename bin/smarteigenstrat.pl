#!/usr/bin/perl

# perl wrapper for smarteigenstrat program.  Run smarteigenstrat.perl with no options for usage

use Getopt::Std;

my @flaglist = ("i","a","b","o","p","q","l","k");
$x = @ARGV;
for($n=0;$n<$x;$n++)  {
  foreach $flag (@flaglist)  {
    if ($ARGV[$n] eq "-$flag")  {
      $specified{$flag} = 1;
    }
  }
}

# check for mandatory options
foreach $flag ("i","a","b","p","o","l")  {
  unless ($specified{$flag})  {
    usage();
    die("Error:  -$flag not specified");
  }
}

# get opts from hash
getopts('i:a:b:p:o:q:k:l:', \%opts);
$genofilename = $opts{"i"};
$indfilename  = $opts{"b"};
$snpfilename  = $opts{"a"};
$pcafilename  = $opts{"p"};
$outfilename  = $opts{"o"};
$logfilename  = $opts{"l"};
$qtmode = "NO";
if ( $specified{"q"} )  {
  $qtmode = $opts{"q"};
}
$k = 1;
if ( $specified{"k"} )  {
  $k = $opts{"k"};
}

# write parameter file
open(PAR, ">$outfilename.par") || die("Error:  unable to open $outfilename.par\n");
print PAR "genotypename:  $genofilename\n";
print PAR "snpname:       $snpfilename\n";
print PAR "indivname:     $indfilename\n";
print PAR "pcaname:       $pcafilename\n";
print PAR "outputname:    $outfilename\n";
print PAR "numpc:         $k\n";
print PAR "qtmode:        $qtmode\n";
close(PAR);

# run smarteigenstrat
$cmd = "smarteigenstrat -p $outfilename.par >$logfilename";
print "$cmd\n";
system($cmd);

sub usage  {
  print "smarteigenstrat.perl -i <genotypefile>  -a <indivfile>  -b <snpfile>  -p <pcafile>   -o <outputfile>  ";
  print "    -l <logfile>  -k <numpca>  -q <qtlmode>";
  print "\n";
  print "-i  genotype file (PED, PACKEDPED, EIGENSTRAT, ANCESTRYMAP or PACKEDANCESTRYMAP format)";
  print "-o  output file (chisq)\n"; 
  print "-l  logfile (screen output,including error messages)\n";
  print "-q  YES for quantitative phenotype or NO otherwise\n";
  print "\n";
  print "For quantitative phenotype, sixth column of .ped file or third column of EIGENSTRAT .ind file\n";
  print "should be real numbers.  For non-quantitative phenotype, sixth column of .ped or third column\n";
  print "should be 'Case' or 'Control'\n";
}






