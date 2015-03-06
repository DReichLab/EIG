#!/usr/bin/perl

$ENV{'PATH'} = "../bin:$ENV{'PATH'}"; 
# MUST put smartpca bin directory in path for smartpca.perl to work

$command = "smartpca.perl";
$command .= " -i example.geno ";
$command .= " -a example.snp ";
$command .= " -b example.QTL.ind " ;
$command .= " -k 2 ";
$command .= " -o example.pca ";
$command .= " -p example.plot ";
$command .= " -e example.eval ";
$command .= " -l example.log ";
$command .= " -m 5 ";
$command .= " -t 2 ";
$command .= " -s 6.0 ";
$command .= " -q YES ";
print("$command\n");
system("$command");

$command = "smarteigenstrat.perl "; 
$command .= " -i example.geno ";
$command .= " -a example.snp ";
$command .= " -b example.QTL.ind ";
$command .= " -p example.pca ";
$command .= " -k 1 ";
$command .= " -o example.QTL.chisq ";
$command .= " -l example.log ";
$command .= " -q YES ";
print("$command\n");
system("$command");

$command = "gc.perl example.QTL.chisq example.QTL.chisq.GC";
print("$command\n");
system("$command");
