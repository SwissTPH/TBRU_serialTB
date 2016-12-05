#!usr/bin/perl
use warnings;

while(<>){
@a=split "\t",$_;
$a[4]=~s/%//;
if($a[4]>90){
@b=split "=",$a[7];
#print "$b[1]\n";
@c=split ":",$b[1];
if(($c[0]>=2)&&($c[1]>=2)){
if($a[8]=~m/Pass=1E0/){
print "$_";
}
}
}
}
