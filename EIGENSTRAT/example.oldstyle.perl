#!/usr/bin/perl

$ENV{'PATH'} = "../bin:$ENV{'PATH'}"; 
# MUST put pca bin directory in path for smartpca.perl to work

$command = "pca";
$command .= " -i example.geno ";
$command .= " -o example.pca ";
$command .= " -e example.eval ";
$command .= " -l example.log ";
$command .= " -k 2 ";
$command .= " -m 5 ";
$command .= " -t 2 ";
$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "eigenstrat"; # or eigenstrat.big.perl for large data sets
$command .= " -i example.geno ";
$command .= " -j example.pheno ";
$command .= " -p example.pca ";
$command .= " -l 1 ";
$command .= " -o example.chisq ";
print("$command\n");
system("$command");

$command = "gc.perl example.chisq example.chisq.GC";
print("$command\n");
system("$command");
