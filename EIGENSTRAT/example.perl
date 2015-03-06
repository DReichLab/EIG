#!/usr/bin/perl

$ENV{'PATH'} = "../bin:$ENV{'PATH'}"; 
# MUST put smartpca bin directory in path for smartpca.perl to work

$command = "smartpca.perl";
$command .= " -i example.geno ";
$command .= " -a example.snp ";
$command .= " -b example.ind " ;
$command .= " -k 2 ";
$command .= " -o example.pca ";
$command .= " -p example.plot ";
$command .= " -e example.eval ";
$command .= " -l example.log ";
$command .= " -m 5 ";
$command .= " -t 2 ";
$command .= " -s 6.0 ";
print("$command\n");
system("$command");

$command = "smarteigenstrat.perl "; 
$command .= " -i example.geno ";
$command .= " -a example.snp ";
$command .= " -b example.ind ";
$command .= " -p example.pca ";
$command .= " -k 1 ";
$command .= " -o example.chisq ";
$command .= " -l example.log ";
print("$command\n");
system("$command");

$command = "gc.perl example.chisq example.chisq.GC";
print("$command\n");
system("$command");
