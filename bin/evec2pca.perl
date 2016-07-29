#!/usr/bin/perl

### translate .evec file to .pca file expected by eigenstrat program
### Note: .evec file does not contain entries for outliers
###       .pca  file does contain entries (set to all 0.0) for outliers

$k = $ARGV[0];
$evec = $ARGV[1]; 
$ind = $ARGV[2];  
$pca = $ARGV[3];  

open(EVEC,$evec) || die("OOPS couldn't open file $evec for reading");
open(IND,$ind) || die("OOPS couldn't open indiv file $ind for reading");
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

while($line = <IND>)
{
  chomp($line);
  my @array = split(/[\t ]+/,$line);
  $sample = $array[0];
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



