#!/usr/bin/perl

$command = "../bin/twstats";
$command .= " -t twtable ";
$command .= " -i twexample.eval ";
$command .= " -o twexample.out";
print("$command\n");
system("$command");

