#!/usr/bin/perl
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: Mon 02 Jul 2012 04:23:24 PM CDT
#Vanderbilt Center for Quantitative Sciences
#############################################
#RNA-Seq DE Pipeline
use strict;
use warnings;
use FindBin qw($Bin);
use Term::ANSIColor qw(:constants)
  ;    #print out message in console with different color
use Data::Dumper;    #To print out the variable structure in the debug mode
use Parallel::ForkManager;    #Run parallel
use Cwd;

#Begin block
BEGIN {
	open( MYLOG, ">run.log" );
	chdir(getcwd);
}

END {
	close MYLOG;
}

my $usage = <<USAGE;
################RNASeqDEPipeline#############################
#Author: Jiang Li
#Email:  riverlee2008\@gmail.com
#Usage:  DirofRNASeqDEPipeline/RNASeqDEPipeline.pl RNASeqDEPipeline.cfg
#URL:    http://github.com/riverlee/RNASeqDEPipeline
#Info:   Copy RNASeqDEPipeline.cfg from RNASeqDEPipeline directory to your 
#        working space and edit its variables
##############################################################
USAGE

my $debug = 0;    #For testing

if ( @ARGV == 0 ) {
	print $usage;
	exit(0);
}

#1) Load configure file

=config data structure
%config=(general=>{BowtieIndex=>string,
				   FileList=>string,
				   TranscriptGtf=>string,
				   Threads=>numeric}
		 tophat=>{-p=>numeric},
		 cuffdiff=>{-p=>numeric,
		 			--upper-quartile-norm=>
				 })
=cut

my $configureFile = $ARGV[0];
my %config;
my $global_start = time();

# Main
#Step 1) Loading config from configure file
_loadConfigure( $configureFile, \%config );
if ($debug) {
	print "After loading configure,data structure of \%config\n";
	print Dumper( \%config );
}

#Step 2) Checking config, not neccessary parameter will be removed
_checkingConfigure( \%config );
if ($debug) {
	print "After checking configure,data structure of \%config\n";
	print Dumper( \%config );
}

#Step 3) Run tophat
my ( $tophatCommandsRef, $samplesRef, $tophatRunningTime ) =
  runTophat( \%config );
if ($debug) {
	print "Sample structure\n";
	print Dumper($samplesRef);
}

#Step 4) Run cuffdiff
my ( $cuffdiffCommand, $cuffdiffRunningTime ) =
  runCuffidff( \%config, $samplesRef );

#Step 5) Replace content in the DEReport.Rmd.template
# Make necessary folder first
if ( !-d "DE" )   { mkdir "DE" }
if ( !-d "DE/R" ) { mkdir "DE/R" }
makeRMarkdown( \%config, $samplesRef, $tophatCommandsRef, $tophatRunningTime,
	$cuffdiffCommand, $cuffdiffRunningTime );

#Step 6) knit Rmd and generate html report
system(
	"R --no-save --slave --vanilla --args DEReport.Rmd <$Bin/RMarkdown2Html.R");

########################################################
#
sub makeRMarkdown {
	my ( $configRef, $samplesRef, $tophatCommandsRef, $tophatRunningTime,
		$cuffdiffCommand, $cuffdiffRunningTime )
	  = @_;
	open( IN, "$Bin/DEReport.Rmd.template" ) or die $!;
	my $codes = join "", <IN>;
	close IN;
	print $codes if ($debug);

	# Replace created time
	my $time = scalar(localtime);
	$codes =~ s/replacecreatedtime/$time/;

	# Replace in the template (tophat and cuffdiff)
	## Estimate count level expression
	my $tophatCommands = join "\n", @{$tophatCommandsRef};
	$codes =~ s/replacetophatcode/$tophatCommands/;
	$codes =~ s/replacetophatrunningtime/$tophatRunningTime/;
	$codes =~ s/replacecuffdiffcode/$cuffdiffCommand/;
	$codes =~ s/replacecuffdiffrunningtime/$cuffdiffRunningTime/;

	# Replace in the template (Differential gene expression in R)
	$codes =~ s/replacethreads/$configRef->{general}->{countThreads}/;
	$codes =~ s/replacebindir/$Bin/g;
	
	# Sort sample by its group
	my @samples =
	  sort { $samplesRef->{$a}->{group} <=> $samplesRef->{$b}->{group} }
	  keys %{$samplesRef};
	my $bamfiles     = "c(";
	my $bamfilenames = "c(";
	foreach my $s (@samples) {
		$bamfiles     .= "'tophat/$s/accepted_hits.bam',";
		$bamfilenames .= "'$s',";
	}
	chop($bamfilenames);
	chop($bamfiles);
	$bamfiles     .= ")";
	$bamfilenames .= ")";

	$codes =~ s/replacebamfiles/$bamfiles/;
	$codes =~ s/replacebamfilenames/$bamfilenames/;
	$codes =~ s/replacegtffile/"$configRef->{general}->{TranscriptGtf}"/;
   
    ## Differential gene expression by DESeq, edgeR, baySeq and TSPM
    my $group="c(";
    foreach my $s (@samples){
    	my $v = 1;
    	if($samplesRef->{$s}->{group} ne $samplesRef->{$samples[0]}->{group}){
    		$v = 2;
    	}
    	$group.="'$v',";
    }
    chop($group);
    $group.=")";
    $codes=~s/replacegroup/$group/g;
    
    
    
	#write output
	print $codes if ($debug);
	open( OUT, ">DEReport.Rmd" ) or die $!;
	print OUT $codes;
}

########################################################

########################################################
#Function to run Tophat and cuffdiff
sub runTophat {
	my ($configRef) = @_;
	my @commands;    #used to generate markdown file for reproducible research
	my $start_time = time();
	info("Running tophat ...");

	#Get samples
	my %samples = _readFileList( $configRef->{general}->{FileList} );

	print Dumper( \%samples ) if ($debug);
	if ( !-d "tophat" ) { mkdir "tophat" }

	my $pm = new Parallel::ForkManager( $configRef->{general}->{Threads} );
	foreach my $samplename ( sort keys %samples ) {

		#generat command
		my $command = "tophat ";

		#Add options
		if ( exists( $configRef->{tophat} ) ) {
			foreach my $param ( sort keys %{ $configRef->{tophat} } ) {
				if ( defined( $configRef->{tophat}->{$param} ) ) {
					$command .=
					  $param . " " . $configRef->{tophat}->{$param} . " ";
				}
				else {
					$command .= $param . " ";
				}
			}
		}

		#define output folder
		$command .=
		  " -o tophat/$samplename " . $configRef->{general}->{BowtieIndex};

		#if paired
		if ( $samples{$samplename}->{paired} == 1 ) {
			$command .= " "
			  . $samples{$samplename}->{file}->{first} . " "
			  . $samples{$samplename}->{file}->{second};
		}
		else {
			$command .= " " . $samples{$samplename}->{file};
		}

		#put output into a log file
		#$command.=" > tophat/$samplename/run.log";

		info($command);    # if ($debug);
		push @commands, $command;

		#Add bam file into %samples
		$samples{$samplename}->{bam} = "tophat/$samplename/accepted_hits.bam";

		#Parallel Run
		my $pid = $pm->start and next;

		system($command);
		$pm->finish;
	}
	$pm->wait_all_children;

	my $runningTime = sec_to_dhms( time() - $start_time );
	return ( \@commands, \%samples, $runningTime );
}

sub runCuffidff {
	my ( $configRef, $samplesRef ) = @_;
	my $start_time = time();

	#my @commands; #used to generate markdown file for reproducible research
	info("Running Cuffdiff ...");

	if ( !-d "DE" ) { mkdir "DE" }
	my %group2sample;
	map { $group2sample{ $samplesRef->{$_}->{group} }->{$_} = undef }
	  keys %{$samplesRef};

	print Dumper( \%group2sample ) if ($debug);

	my $command = "cuffdiff ";

	#Add options
	if ( exists( $configRef->{cuffdiff} ) ) {
		foreach my $param ( sort keys %{ $configRef->{cuffdiff} } ) {
			if ( defined( $configRef->{cuffdiff}->{$param} ) ) {
				$command .=
				  $param . " " . $configRef->{cuffdiff}->{$param} . " ";
			}
			else {
				$command .= $param . " ";
			}
		}
	}

	#define output folder
	$command .=
	  " -o DE/cuffdiff " . $configRef->{general}->{TranscriptGtf} . " ";

	#Add group
	foreach my $group ( sort keys %group2sample ) {
		foreach my $samplename ( sort keys %{ $group2sample{$group} } ) {
			$command .= $samplesRef->{$samplename}->{bam} . ",";
		}
		chop($command);
		$command .= "  ";
	}

	#put output into a log file
	#$command.=" > tophat/$samplename/run.log";

	info($command);    # if ($debug);
	                   #push @commands,$command;
	system($command);
	my $runningTime = sec_to_dhms( time() - $start_time );
	return ( $command, $runningTime );
}

#Invoke by runTophat function to read the file defined in the general section
sub _readFileList {
	my ($in) = @_;
	my %samples;
	open( IN, $in ) or error("File '$in' doesn't exist");
	while (<IN>) {
		next if (/^#|^$/);
		s/\r|\n//g;

		#whichend value in 0,1,2; 0 corresponding to $paired=0(single)
		my ( $samplename, $group, $paired, $whichend, $file ) = split /\s+/;
		$samples{$samplename}->{group}  = $group;
		$samples{$samplename}->{paired} = $paired;    #0 single; 1 paired
		if ($paired) {
			my $s = 'first';                          #5'
			if ( $whichend == 2 ) {
				$s = 'second';                        #3'
			}
			$samples{$samplename}->{file}->{$s} = $file;
		}
		else {
			$samples{$samplename}->{file} = $file;
		}
	}
	return %samples;
}

#End Function to run Tophat and cuffdiff
####################################################################

#####################################################################
#Functions related to load, check config
#Load the configure files into %config variable
sub _loadConfigure {
	my ( $configureFile, $ref ) = @_;
	info("Loading configure file ...");
	open( IN, $configureFile )
	  or error("Can't open configure file '$configureFile'");
	my $string = "";
	while (<IN>) {
		if (/\[\[(\w+)\]\]/) {
			if ( $string ne "" && $string =~ /\[\[(\w+)\]\]/ ) {

				#begin parse
				_parseConfigure( $string, $ref );
			}
			$string = $_;
		}
		else {
			$string .= $_;
		}
	}
	close IN;
	if ( $string =~ /\[\[(\w+)\]\]/ ) {
		_parseConfigure( $string, $ref );
	}
}

#parse each section information in the configure file, the  section begin with [[]]
sub _parseConfigure {
	my ( $string, $ref ) = @_;
	my $parentKey = "";
	if ( $string =~ /\[\[(\w+)\]\]/ ) {
		$parentKey = $1;
	}
	else {
		error(
			"parsing cofigure file failed, please refer to the template files 
		in the template directory for instruction"
		);
		exit 1;
	}

	if ( $parentKey eq "pbs" ) {
		my $t = $string;
		$t =~ s/\[\[pbs\]\]\s*\r?\n?//g;
		$ref->{$parentKey} = $t;
		return;
	}
	else {
		foreach ( split /\r|\n/, $string ) {
			next if (/#|\[\[\w+\]\]|^$/);
			chomp;
			my @GAcmd = split( /\s*=\s*/, $_ );
			@GAcmd = map { $_ =~ s/\s+//g; $_; }
			  @GAcmd;    #skip the blank in key and value
			$ref->{$parentKey}->{ $GAcmd[0] } = $GAcmd[1];
		}
	}
}

#After running '_loadConfigure' function to store information into %config variable,
#we will check
sub _checkingConfigure {
	my ($configRef) = @_;
	info("Checking configure file ...");

	my $flag = 0;        #flag changed to 1 if there is an error
	my $msg  = "";       #store the error message

	my $warningflag = 0;
	my $warningmsg  = "";

	#checking error
	if ( !exists( $configRef->{general} ) ) {
		$flag++;
		$msg .= "$flag) 'general section doesn't exists\n'";
	}
	else {

		#bowtie
		if ( !exists( $configRef->{general}->{BowtieIndex} ) ) {
			$flag++;
			$msg .= "$flag) 'BowtieIndex' is not defined \n";
		}

		#GTF
		if ( !exists( $configRef->{general}->{TranscriptGtf} ) ) {
			$flag++;
			$msg .= "$flag) 'TranscriptGtf' is not defined \n";
		}
		elsif ( !-e $configRef->{general}->{TranscriptGtf} ) {
			$flag++;
			$msg .=
"$flag) gtf file '$configRef->{general}->{TranscriptGtf}' doesn't exists \n";
		}

		#FileList
		if ( !exists( $configRef->{general}->{FileList} ) ) {
			$flag++;
			$msg .= "$flag) 'FileList' is not defined \n";
		}
		elsif ( !-e $configRef->{general}->{FileList} ) {
			$flag++;
			$msg .=
"$flag) FileList file '$configRef->{general}->{FileList}' doesn't exists \n";
		}

		#Threads
		if ( !exists( $configRef->{general}->{Threads} ) ) {
			$flag++;
			$msg .= "$flag) 'Threads' is not defined \n";
		}
		elsif ( $configRef->{general}->{Threads} !~ /^\d*$/ ) {
			$flag++;
			$msg .=
"$flag) 'Threads=$configRef->{general}->{Threads}' is not numeric \n";
		}
	}

	#check tophat warning
	if ( !exists( $configRef->{tophat} ) ) {
		$warningflag++;
		$warningmsg .=
		  "$warningflag) No parameters defined for tophat, will use default\n";
	}
	else {

		#get the tophat -- or -
		my @tophat_params = split "\n",
		  `tophat 2>&1 |awk '{print \$1}' | grep '^-' |tr '/' '\n'`;
		my %tophat_hash;
		@tophat_hash{@tophat_params} = undef;

		#Loop parameters
		foreach my $param ( keys %{ $configRef->{tophat} } ) {
			if (   $param eq '-o'
				|| $param eq '--output-dir'
				|| $param eq '-h'
				|| $param eq "--help"
				|| $param eq "-v"
				|| $param eq "--versbose" )
			{
				delete( $configRef->{tophat}->{$param} );
				next;
			}

			if ( !exists( $tophat_hash{$param} ) ) {
				$warningflag++;
				$warningmsg .=
"$warningflag) Option '$param' is not implemented in tophat\n";
				delete( $configRef->{tophat}->{$param} );
			}
		}
	}

	if ($flag) {
		error($msg);
	}

	#check cuffdiff warning
	if ( !exists( $configRef->{cuffdiff} ) ) {
		$warningflag++;
		$warningmsg .=
"$warningflag) No parameters defined for cuffdiff, will use default\n";
	}
	else {

		#get the cuffdiff -- or -
		my @cuffdiff_params = split "\n",
		  `cuffdiff 2>&1 |awk '{print \$1}' | grep '^-' |tr '/' '\n'`;
		my %cuffdiff_hash;
		@cuffdiff_hash{@cuffdiff_params} = undef;

		#Loop parameters
		foreach my $param ( keys %{ $configRef->{cuffdiff} } ) {
			if (   $param eq '-o'
				|| $param eq '--output-dir'
				|| $param eq '-h'
				|| $param eq "--help"
				|| $param eq "-v"
				|| $param eq "--versbose" )
			{
				delete( $configRef->{cuffdiff}->{$param} );
				next;
			}

			if ( !exists( $cuffdiff_hash{$param} ) ) {
				$warningflag++;
				$warningmsg .=
"$warningflag) Option '$param' is not implemented in cuffdiff\n";
				delete( $configRef->{tophat}->{$param} );
			}
		}
	}

	#mean it exists error when parsing inputed parameters
	if ($warningflag) {
		warning($warningmsg);
	}

}

#This is used for tophat parameters check
sub _tophat_param_check {
	my ($param) = @_;

	if ( $param eq '-o' || $param eq '--output-dir' ) {
		return 1;
	}

	#get the tophat -- or -
	my @tophat_params = split "\n",
	  `tophat 2>&1 |awk '{print \$1}' | grep '^-' |tr '/' '\n'`;
	my %hash;
	@hash{@tophat_params} = undef;

	#print Dumper(\%hash);
	if ( !exists( $hash{$param} ) ) {
		warning( $param . " option is not implemented in tophat" );
		return 1;
	}

	return 0;
}

#End Functions related to load, check config
#####################################################################

###################################################
#Functions to print out message and write to log
sub debug {
	my ($msg) = @_;
	my $str = "[Debug][" . scalar(localtime) . "]:" . $msg . "\n";
	print YELLOW, $str;
	logging($str);
}

#default will exist program
sub error {
	my ( $msg, $flag ) = @_;
	my $str = "[Error][" . scalar(localtime) . "]:" . $msg . "\n";
	print RED, $str;
	logging($str);
	unless ($flag) {
		exit 1;
	}
}

sub warning {
	my ($msg) = @_;
	my $str = "[Warning][" . scalar(localtime) . "]:" . $msg . "\n";
	print YELLOW, $str;
	logging($str);
}

sub info {
	my ($msg) = @_;
	my $str = "[Info][" . scalar(localtime) . "]:" . $msg . "\n";
	print GREEN, $str;
	logging($str);
}

sub logging {
	my ($str) = @_;
	print MYLOG $str;
}

sub sec_to_dhms {
	use integer;
	local $_ = shift;
	my ( $d, $h, $m, $s );
	$s = $_ % 60;
	$_ /= 60;
	$m = $_ % 60;
	$_ /= 60;
	$h = $_ % 24;
	$_ /= 24;
	$d = $_;
	return sprintf( "%02d:%02d:%02d:%02d", $d, $h, $m, $s );
}

#End Functions to print out message and write to log
###############################################

