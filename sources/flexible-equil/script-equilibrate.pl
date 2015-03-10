#!/usr/bin/perl

use strict;
use warnings;

#---logfile
our $logfile;
open $logfile, ">>", 'log-script-master' or die "$!";
open STDERR, ">>", 'log-script-master';

#---tools
use lib '../amx/';
use gmxwrap (
	'reporter',
	'timestring',
	'errorscan',
	'getbashstr',
	'syscatch',
	'fileparse',
	'getmaxforce',
	'finish',
	);

#---inputs
system("bash settings.sh");
my $stepstring = getbashstr("settings.sh","STEPSTRING");
my @steps = split(/\s+/, $stepstring);
my %steplist = map { $_ => 1 } @steps;
my %config = fileparse('../gmxpaths.conf');

#---reuse variables
my $cmd;
my $logname;

#-------------------------------------------------------------------------------------------------------------

reporter("EQUILIBRATION ROUTINE");
reporter timestring;

my $num;
my $conf;
my @files = glob("input-md-stage-*-in.mdp");
@files = sort {$a cmp $b} @files;
foreach my $file (@files) {
    print "$file\n";
}
foreach my $file ( @files ) {
	($num) = $file =~ /^input-md-stage-(\d+)/;
	if (-e "md-stage-" . $num . ".gro") { 
		reporter("found md-stage-" . $num . ".gro so skipping to the next step");
	} else {
		reporter "running: grompp, stage " . $num;
		$conf = ($num == "0") ? "system-input.gro" : "md-stage-" . ($num-1) . ".gro";
		$logname = "log-grompp-md-stage-" . $num;
		$cmd = $config{"grompp"};
		$cmd .= " -f input-md-stage-" . $num . "-in.mdp";
		$cmd .= " -c " . $conf;
		$cmd .= " -o md-stage-" . $num;
		$cmd .= " -p system.top";
		$cmd .= " -n system-groups.ndx";
		$cmd .= " -po md-stage-" . $num;
		$cmd .= " -maxwarn 10";
		$cmd .= " > " . $logname . " 2>&1";
		reporter "executing " . $cmd;
		syscatch($cmd);
		errorscan($logname)==0 or finish "$!";
		reporter "logfile = " . $logname;
		reporter timestring;
		reporter "running: mdrun, equilibration, stage " . $num;
		$logname = "log-mdrun-md-stage-" . $num;
		$cmd = $config{"mdrun"};
		$cmd .= " -v";
		$cmd .= " -deffnm md-stage-" . $num;
		$cmd .= " > " . $logname . " 2>&1";
		syscatch($cmd);
		errorscan($logname)==0 or finish "$!";
		reporter "logfile = " . $logname;
		reporter timestring;
	}
}

if (-e "md.part0001.tpr") { 
	reporter("found md.part0001.tpr so skipping to the next step");
} else {
	#---PRODUCTION MD RUN
	reporter "running: grompp, production MD run";
	$logname = "log-grompp-md-0001";
	$conf = ($num == "0") ? "system-input.gro" : "md-stage-" . $num . ".gro";
	$cmd = $config{"grompp"};
	$cmd .= " -f input-md-in.mdp";
	$cmd .= " -c " . $conf;
	$cmd .= " -o md.part0001";
	$cmd .= " -p system.top";
	$cmd .= " -n system-groups.ndx";
	$cmd .= " -po md.part0001";
	$cmd .= " -maxwarn 10";
	$cmd .= " > " . $logname . " 2>&1";
	syscatch($cmd);
	errorscan($logname)==0 or finish "$!";
	reporter "logfile = " . $logname;
	reporter timestring;
	reporter "running: mdrun, production MD run";
	$logname = "log-mdrun-md-0001";
	$cmd = $config{"mdrun"};
	$cmd .= " -v";
	$cmd .= " -s md.part0001";
	$cmd .= " -x md.part0001";
	$cmd .= " -g md.part0001";
	$cmd .= " -o md.part0001";
	$cmd .= " -c md.part0001";
	$cmd .= " -e md.part0001";
	$cmd .= " -cpo md.part0001";
	$cmd .= " > " . $logname . " 2>&1";
	syscatch($cmd);
	errorscan($logname)==0 or finish "$!";
	reporter "logfile = " . $logname;
	reporter timestring;
	#---job complete
	reporter timestring;
	if (-e "md.part0001.gro") { reporter "completed"; }
}
