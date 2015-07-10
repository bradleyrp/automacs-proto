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

#---NVT EQUILIBRATION, SHORT TIMESTEP
if ($steplist{"nvt-short"} == 1) {
	reporter "running: grompp, equilibration, NVT, short step";
	$logname = "log-grompp-md-nvt-shortstep";
	$cmd = $config{"grompp"};
	$cmd .= " -f input-md-nvt-eq-short-in.mdp";
	$cmd .= " -c system-input.gro";
	$cmd .= " -o md-nvt-eq-short";
	$cmd .= " -p system.top";
	$cmd .= " -n system-groups.ndx";
	$cmd .= " -po md-nvt-eq-short";
	$cmd .= " -maxwarn 10";
	$cmd .= " > " . $logname . " 2>&1";
	syscatch($cmd);
	errorscan($logname)==0 or finish "$!";
	reporter "logfile = " . $logname;
	reporter timestring;
	reporter "running: mdrun, equilibration, NVT, short step";
	$logname = "log-mdrun-md-nvt-shortstep";
	$cmd = $config{"mdrun"};
	$cmd .= " -v";
	$cmd .= " -deffnm md-nvt-eq-short";
	$cmd .= " > " . $logname . " 2>&1";
	syscatch($cmd);
	errorscan($logname)==0 or finish "$!";
	reporter "logfile = " . $logname;
	reporter timestring;
} else { 
	reporter "skipping nvt-short";
	if (-e "md-nvt-eq-short.gro") {
		reporter "using preexisting md-nvt-eq-short.gro"; 
	} else {
		$cmd = "cp system-input.gro md-nvt-eq-short.gro";
		reporter $cmd;
		system($cmd);
	}
}

#---NVT EQUILIBRATION
if ($steplist{"nvt"} == 1) {
	reporter "running: grompp, equilibration, NVT";
	$logname = "log-grompp-md-nvt";
	$cmd = $config{"grompp"};
	$cmd .= " -f input-md-nvt-eq-in.mdp";
	$cmd .= " -c md-nvt-eq-short.gro";
	$cmd .= " -o md-nvt-eq";
	$cmd .= " -p system.top";
	$cmd .= " -n system-groups.ndx";
	$cmd .= " -po md-nvt-eq";
	$cmd .= " -maxwarn 10";
	$cmd .= " > " . $logname . " 2>&1";
	syscatch($cmd);
	errorscan($logname)==0 or finish "$!";
	reporter "logfile = " . $logname;
	reporter timestring;
	reporter "running: mdrun, equilibration, NVT";
	$logname = "log-mdrun-md-nvt";
	$cmd = $config{"mdrun"};
	$cmd .= " -v";
	$cmd .= " -deffnm md-nvt-eq";
	$cmd .= " > " . $logname . " 2>&1";
	syscatch($cmd);
	errorscan($logname)==0 or finish "$!";
	reporter "logfile = " . $logname;
	reporter timestring;
} else { 
	reporter "skipping nvt";
	if (-e "md-nvt-eq.gro") {
		reporter "using preexisting md-nvt-eq.gro"; 
	} else {
		$cmd = "cp md-nvt-eq-short.gro md-nvt-eq.gro";
		reporter $cmd;
		system($cmd);
	}
}

#---NPT EQUILIBRATION
if ($steplist{"npt"} == 1) {
	reporter "running: grompp, equilibration, NPT";
	$logname = "log-grompp-md-npt";
	$cmd = $config{"grompp"};
	$cmd .= " -f input-md-npt-eq-in.mdp";
	$cmd .= " -c md-nvt-eq.gro";
	$cmd .= " -o md-npt-eq";
	$cmd .= " -p system.top";
	$cmd .= " -n system-groups.ndx";
	$cmd .= " -po md-npt-eq";
	$cmd .= " -maxwarn 10";
	$cmd .= " > " . $logname . " 2>&1";
	syscatch($cmd);
	errorscan($logname)==0 or finish "$!";
	reporter "logfile = " . $logname;
	reporter timestring;
	reporter "running: mdrun, equilibration, NPT";
	$logname = "log-mdrun-md-npt";
	$cmd = $config{"mdrun"};
	$cmd .= " -v";
	$cmd .= " -deffnm md-npt-eq";
	$cmd .= " > " . $logname . " 2>&1";
	syscatch($cmd);
	errorscan($logname)==0 or finish "$!";
	reporter "logfile = " . $logname;
	reporter timestring;
} else { 
	reporter "skipping npt";
	if (-e "md-npt-eq.gro") {
		reporter "using preexisting md-npt-eq.gro"; 
	} else {
		$cmd = "cp md-nvt-eq.gro md-npt-eq.gro";
		reporter $cmd;
		system($cmd);
	}
}
#---if this fails, particularly for the bilayer builder, we use the NVT structure
if (! -e "md-npt-eq.gro") {
	reporter "it looks like NPT crashed so continue with md-nvt-eq.gro"; 
	$cmd = "cp md-nvt-eq.gro md-npt-eq.gro";
	reporter $cmd;
	system($cmd);
}

#---PRODUCTION MD RUN
reporter "running: grompp, production MD run";
$logname = "log-grompp-md-0001";
$cmd = $config{"grompp"};
$cmd .= " -f input-md-in.mdp";
$cmd .= " -c md-npt-eq.gro";
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

