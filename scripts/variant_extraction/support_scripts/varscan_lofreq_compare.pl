#!usr/bin/perl
use warnings;

my %hash;

open F1, $ARGV[0] or die $!; #open lofreq.snp
while(<F1>){
        chomp;
        @a=split "\t",$_;
        $hash{$a[0]}=1;
    }
                close F1;

open F2, $ARGV[1] or die $!;  #open varscan.snp
while(<F2>){
    chomp;
    @a=split "\t",$_;
    if(exists $hash{$a[0]}){
        print "$_\n";
        }
    }
close F2;