#!/usr/bin/perl
use strict;
use warnings;

my %hash=("SRR002320"=>"kidney3pm_run1",
	  "SRR002325"=>"kidney3pm_run2",
	  "SRR002321"=>"liver3pm_run1",
	  "SRR002323"=>"liver3pm_run2",
	"SRR002324"=>"kidney1.5pm_run2",
	"SRR002322"=>"liver1.5pm_run2");

my @files=<*lane*>;
#print join "\n",@files;
foreach my $f (@files){
	my ($id) = split /_/,$f;
	print $id,"\t";
	if(exists($hash{$id})){
		print "OK\n";
		my $newfile=$f;
		$newfile=~s/$id/$hash{$id}/g;
	#	print $newfile ,"\n";
		`mv $f $newfile`;
	}else{
		print "NO\n";
	}
	
}
