#!/usr/bin/perl

$ENV{'PATH'} = "../bin:$ENV{'PATH'}";

$command = "smartpca";
$command .= " -p par.example >example.log";
print("$command\n");
system("$command");

$command = "ploteig";
$command .= " -i example.evec ";
$command .= " -c 1:2 ";
$command .= " -p Case:Control ";
$command .= " -x ";
$command .= " -o example.plot.xtxt "; # must end in .xtxt
print("$command\n");
system("$command");

$command = "evec2pca.perl 2 example.evec example.ind example.pca";
print("$command\n");
system("$command");
