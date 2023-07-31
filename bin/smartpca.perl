#!/usr/bin/env perl

use Getopt::Std ;
use File::Basename ;

sub usage {
    my $message = "@_";
#       10       20       30       40       50       60       70       80       90
#---+----|---+----|---+----|---+----|---+----|---+----|---+----|---+----|---+----|
    die "
Usage: smartpca.perl [FLAGS]

This program calls the smartpca program (see POPGEN/README). 
For this to work, the bin directory containing smartpca MUST be in your path. 

Required Flags:
  -i example.geno  : genotype file in any format (see CONVERTF/README)
  -a example.snp   : snp file in any format (see CONVERTF/README)
  -b example.ind   : indiv file in any format (see CONVERTF/README)
  -k k             : (Default is 10) number of principal components to output
  -o example.pca   : output file of principal components.  Individuals removed
                     as outliers will have all values set to 0.0 in this file.
  -p example.plot  : prefix of output plot files of top 2 principal components.
                     (labeling individuals according to labels in indiv file)
  -e example.eval  : output file of all eigenvalues
  -l example.log   : output logfile
  -m maxiter       : (Default is 5) maximum number of outlier removal iterations.
                     To turn off outlier removal, set -m 0.
  -t topk          : (Default is 10) number of principal components along which 
                     to remove outliers during each outlier removal iteration.
  -s sigma         : (Default is 6.0) number of standard deviations which an
                     individual must exceed, along one of topk top principal
                     components, in order to be removed as an outlier.

Optional Flags:
  -w poplist       : compute eigenvectors using populations in poplist only,
                     where poplist is an ASCII file with one population per line
  -y plotlist      : output plot will include populations in plotlist only, 
                     where plotlist is an ASCII file with one population per line
  -z badsnpname    : list of SNPs which should be excluded from the analysis
  -q YES/NO        : If set to YES, assume that there is a single population and
                     the population field contains real-valued phenotypes.
                     (Corresponds to qtmode parameter in smartpca program.)
                     The default value for this parameter is NO.

Estimated running time of the smartpca program is 
  2.5e-12 * nSNP * NSAMPLES^2 hours            if not removing outliers.
  2.5e-12 * nSNP * NSAMPLES^2 hours * (1+m)    if m outlier removal iterations.
Thus, under the default of up to 5 outlier removal iterations, running time is 
  up to 1.5e-11 * nSNP * NSAMPLES^2 hours.

Recommendation: we advise after running pca to check for large correlations
between each principal component and all variables of interest.  For example,
large correlations with % missing data (per sample) could imply assay issues
large correlations with plate membership could imply assay issues
large correlations with phenotype indicate highly mismatched cases vs. controls
  which will lead to a loss of power upon applying eigenstrat correction.
  If input indiv file contains Case and Control labels only, then
  correlations between each principal component and Case/Control status will be
  listed at end of output logfile (-l flag).

$message

"
}
### process flags
# -w poplist is compute eigenvectors using populations in poplist, then project
# -y poplistplot is use populations in poplistplot for plot
# -z badsnpfile is use badsnpname: badsnpfile in call to smartpca
my @flaglist = ("i","a","b","k","o","p","e","l","m","q","t","s","w","y","z");
$x = @ARGV;
for($n=0; $n<$x; $n++)
{
  foreach $flag (@flaglist) 
  {
    if($ARGV[$n] eq "-$flag") { $specified{$flag} = 1; }
  }
}
foreach $flag ("i","a","b","o","p","e","l")
{
  unless($specified{$flag}) { usage("OOPS -$flag flag not specified"); }
}
getopts('i:a:b:k:o:p:e:l:m:t:s:w:y:z:q:',\%opts);
$i = $opts{"i"}; 
$a = $opts{"a"}; 
$b = $opts{"b"}; 
$k = 10; if($specified{"k"}) { $k = $opts{"k"}; }
$o = $opts{"o"}; 
$q = 0; if($specified{"q"}) { $q = $opts{"q"}; }
$p = $opts{"p"}; 
$e = $opts{"e"}; 
$l = $opts{"l"}; 
$m = 5; if($specified{"m"}) { $m = $opts{"m"}; }
$t = 10; if($specified{"t"}) { $t = $opts{"t"}; }
$s = 6.0; if($specified{"s"}) { $s = $opts{"s"}; }
if($specified{"w"}) { $w = $opts{"w"}; }
if($specified{"y"}) { $y = $opts{"y"}; }
if($specified{"z"}) { $z = $opts{"z"}; }

### run smartpca
$parfile = "$o.par"; 
$evec = "$o.evec";
open(PARFILE,">$parfile") || die("OOPS couldn't open file $parfile for writing");
print PARFILE ("genotypename: $i\n");
print PARFILE ("snpname: $a\n");
print PARFILE ("indivname: $b\n");
print PARFILE ("evecoutname: $evec\n");
print PARFILE ("evaloutname: $e\n");
print PARFILE ("altnormstyle: NO\n");
print PARFILE ("numoutevec: $k\n");
print PARFILE ("numoutlieriter: $m\n");
print PARFILE ("numoutlierevec: $t\n");
print PARFILE ("outliersigmathresh: $s\n");
print PARFILE ("qtmode: $q\n");
if($specified{"w"}) { print PARFILE ("poplistname: $w\n"); }
if($specified{"z"}) { print PARFILE ("badsnpname: $z\n"); }
close(PARFILE);
$command = "smartpca";          # MUST put bin directory in path
$command .= " -p $parfile >$l";
print("$command\n");
system("$command");

### make string of populations for ploteig
$popstring = "";
open(EVEC,$evec) || die("OOPS couldn't open file $evec for reading");
while($line = <EVEC>)
{
  chomp($line);
  my @array = split(/[\t ]+/,$line);
  $x = @array;
  if($array[1] =~ /eigvals/) { next; } # eigvals header line
  $pop = $array[$x-1];
  if($popfound{$pop}) { next; }
  $popstring = $popstring . "$pop:";
  $popfound{$pop} = 1;
}
close(EVEC);
chop($popstring); # remove last ":"

if($specified{"y"}) 
{
  ### make string of populations for ploteig based on -y flag input
  $popstring = "";
  open(Y,$y) || die("Cannot open file: $y");
  while($line = <Y>)
  {
    chomp($line);
    $popstring .= "$line:";
  }
  chop($popstring);
}

### cax ploteig
$command = "ploteig";           # MUST put bin directory in path
$command .= " -i $evec";
$command .= " -c 1:2 ";
$command .= " -p $popstring ";
$command .= " -x ";
$command .= " -y ";
$command .= " -o $p.xtxt "; # must end in .xtxt
print("$command\n");
system("$command");

### translate .evec file to .pca file expected by eigenstrat program
### Note: .evec file does not contain entries for outliers 
###       .pca  file does contain entries (set to all 0.0) for outliers

# ----- If this looks like a PLINK run, call the PLINK kludge
if ( $i =~ m/\.ped$/ || $i =~ m/\.PED/ )  {
  $command = "evec2pca-ped.perl $k $evec $b $o";
}
else  {
  $command = "evec2pca.perl $k $evec $b $o";
}
print("$command\n");
system("$command");

### If labels are Case and Control only, compute correlations between
### Case/Control status and each eigenvector.  Append to logfile.
if(($popstring eq "Case:Control") || ($popstring eq "Control:Case"))
{
  open(LOG,">>$l") || die("OOPS couldn't open file $l for appending");
  print LOG ("\n");
  for($x=0; $x<$k; $x++) # compute correlation for evec $x
  {
    open(EVEC,$evec) || die("OOPS couldn't open file $evec for reading");
    $sum1=0; $sumx=0; $sumxx=0; $sumy=0; $sumyy=0; $sumxy=0;
    $line = <EVEC>; chomp($line); # eigvals line
    while($line = <EVEC>)
    {
      chomp($line);
      my @array = split(/[\t ]+/,$line);
      $this = $array[2+$x];
      $sumy += $this;
      $sumyy += $this*$this;
      $sum1 += 1;
      if($line =~ /Case/) # Case is 1, Control is 0
      {
        $sumx += 1;
        $sumxx += 1;
        $sumxy += $this;
      }
    }
    close(EVEC);
    $meanx = $sumx/$sum1;
    $meany = $sumy/$sum1;
    if($sum1 == 0) { next; }
    $sdevx = sqrt($sumxx/$sum1 - $meanx*$meanx);
    $sdevy = sqrt($sumyy/$sum1 - $meany*$meany);
    if($sdevx * $sdevy == 0) { next; }
    $corr = ($sumxy/$sum1) / ($sdevx*$sdevy);
    $x1 = $x+1;
    printf LOG ("Correlation between eigenvector $x1 (of $k) and Case/Control status is %.03f\n",$corr);
  }
  close(LOG);
}
