#!/usr/bin/env perl

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
    usage("OOPS -$flag not specified");
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
open(PAR, ">$outfilename.par") || die("OOPS unable to open $outfilename.par\n");
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
    my $message = "@_";
#       10       20       30       40       50       60       70       80       90
#---+----|---+----|---+----|---+----|---+----|---+----|---+----|---+----|---+----|
    die "
Usage: smarteigenstrat.perl [FLAGS]

This program is a PERL wrapper which calls the C program smarteigenstrat. 
  Note: the bin directory containing smarteigenstrat MUST be in your path. 
  See ./example.perl for a toy example.
We recommend smarteigenstrat.perl for users who prefer command-line flags.
However, users who prefer parameter files can run smarteigenstrat instead.

Required Flags:
  -i example.geno  : genotype file in any format (PED, PACKEDPED, EIGENSTRAT, 
                     ANCESTRYMAP or PACKEDANCESTRYMAP format)
  -a example.snp   : snp file in any format (see CONVERTF/README)
  -b example.ind   : individual file in any format (see CONVERTF/README).
                     We note that phenotype information will be contained in
                     example.ind, either as Case/Control labels or quantitative
                     phenotypes if -q set to YES.
  -q YES/NO        : If set to YES, use quantitative phenotypes in example.ind.
                     If -q is set to YES, the third column of the input individual
                     file in EIGENSTRAT format (or sixth column of input 
                     individual file in PED format) should be real numbers. The
                     value -100.0 signifies 'missing data'. If -q is set to NO,
                     these values should be 'Case' or 'Control'. The default value
                     for the -q parameter is NO.     
  -p example.pca   : input file of principal components (output of smartpca.perl)
  -k topk          : (Default is 10) number of principal components along which to
                     correct for stratification.  Note that l must be less than or
                     equal to the number of principal components reported in the 
                     file example.pca.
  -o example.chisq : chisq association statistics.  File contains log of flags to
                     eigenstrat program, followed by one line per SNP:
                     The first entry of each line is Armitage chisq statistic 
                       (Armitage, 1955) defined as NSAMPLES x (correlation between
                       genotype and phenotype)^2. If the set of individuals with 
                       genotype and phenotype both valid is monomorphic for either
                       genotype or phenotype, then NA is reported.
                     The second entry of each line is the EIGENSTRAT chisq statistic,
                       defined as (NSAMPLES-l-1) x (corr between adjusted_genotype 
                       and adjusted_phenotype)^2. If the set of individuals with 
                       genotype and phenotype both valid is monomorphic for either
                       genotype or phenotype, then NA is reported.
                     Note: even if l=0, there is a tiny difference between the two
                       statistics because Armitage uses NSAMPLES while we use 
                       NSAMPLES-1, which we consider to be appropriate.
  -l example.log   : standard output logfile (screen output,including error messages)

The running time of smarteigenstrat.perl is very fast compared to the running time
  of smartpca.perl.

For quantitative phenotype, sixth column of .ped file or third column of EIGENSTRAT
.ind file should be real numbers.  For non-quantitative phenotype, sixth column of
.ped or third column should be 'Case' or 'Control'

$message

"
}






