#!/usr/bin/perl
use strict;
use warnings;
use IO::File;

my @files=<*.fastq>;
#print join "\n", @files;
foreach my $f (@files){
    my($id,undef) = split /\./, $f;
    `rm ${id}_lane*`;
    print info(),"Reading $f ..\n";
    my %hash;
    open(IN,$f) or die $!;
    while(<IN>){
        my $lane = (split(/:/, (split(/\s+/,$_))[1]))[1];
        #print $lane,"\n";
        my $filename="${id}_lane${lane}.txt";
        if(!exists($hash{$filename})){
            my $fh=new IO::File;
            $fh->open(">$filename");
            $hash{$filename}=$fh;
        }
        my $fh=$hash{$filename};
        print $fh $_;
        my $l=<IN>;
        print $fh $l;
       $l= <IN>;
        print $fh $l;
        $l=<IN>;
        print $fh  $l;
    }
    close IN;
    foreach my $key (keys %hash){
        $hash{$key}->close;
        delete($hash{$key});
    }
}
sub info{
    return "[",scalar(localtime),"] @ ";
}
