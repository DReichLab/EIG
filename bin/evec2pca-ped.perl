#!/usr/bin/env perl

### translate .evec file to .pca file expected by eigenstrat program
### Note: .evec file does not contain entries for outliers
###       .pca  file does contain entries (set to all 0.0) for outliers

# ----- This is a new version for PLINK input files.  It differs from the
# ----- original in two ways.  (1) The indiv name is in the second column 
# ----- of the .fam file, but the first of a .ind file.  (2) If the 
# ----- indiv names are not found in the .evec file, try the 
# ----- familyname:indivname combination.  

$k = $ARGV[0];
$evec = $ARGV[1]; 
$ind = $ARGV[2];  
$pca = $ARGV[3];  
open(EVEC,$evec) || die("OOPS couldn't open file $evec for reading");
open(PCA,">$pca") || die("OOPS couldn't open file $pca for writing");

print PCA ("$k\n"); # number of output eigenvectors/eigenvalues
$line = <EVEC>; chomp($line); # eigvals line
my @array = split(/[\t ]+/,$line);
for($x=0; $x<$k; $x++) { printf PCA ("%.04f\n",$array[$x+2]); } # x-th eval
while($line = <EVEC>)
{
  chomp($line);
  $line = " " . $line;
  my @array = split(/[\t ]+/,$line);
  $l = @array; 
  unless($l == 3+$k) { die("OOPS #evec in $evec is different from $k"); }
  $sample = $array[1];
  for($x=0; $x<$k; $x++) { $evecarray{$sample}[$x] = $array[$x+2]; }
  $found{$sample} = 1;
}

# ----- Figure out which name convention to use
my $count1 = 0;
my $count2 = 0;
open(IND,$ind) || die("OOPS couldn't open indiv file $ind for reading");
while ( my $line = <IND> )  {
  chomp($line);
  $line =~ s/^[\s]+//;     # remove leading white-space
  my @E = split(/[\s]+/,$line);

  my $s = $E[1];
  my $t = $E[0] . ":" . $E[1];

  if ( exists $found{$s} )  {
    $count1++;
  }
  if ( exists $found{$t} )  {
    $count2++;
  }
}
close(IND);

open(IND,$ind) || die("OOPS couldn't open indiv file $ind for reading");
while($line = <IND>)
{
  chomp($line);
  $line =~ s/^[\s]+//;
  my @array = split(/[\s]+/,$line);
  $sample = ($count1 >= $count2 ? $array[1] : $array[0] . ":" . $array[1]);
  if($sample eq "") { $sample = $array[1]; }
  unless($found{$sample})
  {
    for($x=0; $x<$k; $x++) { $evecarray{$sample}[$x] = 0.0; }
  }
  for($x=0; $x<$k; $x++)
  {
    printf PCA (" ");
    if($evecarray{$sample}[$x] > 0) { printf PCA (" "); }
    printf PCA ("%.04f",$evecarray{$sample}[$x]);
  }
  printf PCA ("\n");
}



