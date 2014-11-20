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

#---extension time
my $extend_time = getbashstr("settings.sh","EXTENDTIME");;

#-------------------------------------------------------------------------------------------------------------

reporter("SIMULATION RESTART");
reporter timestring;

#---get index of most recent cpt file
opendir(MYDIR,"./") or die "$!";
my @files = grep{ /md.part[0-9]{4}\.cpt/ } readdir(MYDIR);
close MYDIR;
my $last = 0;
foreach (@files) { 
	$_ =~ s/\D//g;
	$_ =~ s/^0+//g;
	if ($_ gt $last) { $last = $_ };
}
my $prun = sprintf("%04d", $last);
my $nrun = sprintf("%04d", $last+1);
reporter "running: grompp, production MD run restart";
$logname = "log-grompp-md-" . $nrun;
$cmd = $config{"grompp"};
$cmd .= " -f input-md-in.mdp";
if (-f "md.part" . $prun . ".gro") { $cmd .= " -c md.part" . $prun . ".gro"; }
elsif (-f "md.part" . $prun . ".tpr" . $prun . ".gro") { $cmd .= " -c md.part" . $prun . ".tpr"; }
elsif (-f "system-input.gro") { $cmd .= " -c system-input.gro"; }
$cmd .= " -o md.part" . $nrun . ".tpr";
$cmd .= " -po md.part" . $nrun . ".mdp";
$cmd .= " -p system.top";
$cmd .= " -maxwarn 10";
$cmd .= " -n system-groups.ndx" if -f "system-groups.ndx";
$cmd .= " > " . $logname . " 2>&1";
print $cmd . "\n";
syscatch($cmd);
errorscan($logname)==0 or finish "$!";
reporter "logfile = " . $logname;
reporter timestring;
reporter "running: mdrun, production MD run continuation";
$logname = "log-mdrun-md-" . $nrun;
$cmd = $config{"mdrun"};
$cmd .= " -s md.part" . $nrun . ".tpr";
$cmd .= " -o md.part" . $nrun . ".trr";
$cmd .= " -x md.part" . $nrun . ".xtc";
$cmd .= " -c md.part" . $nrun . ".gro";
$cmd .= " -g md.part" . $nrun . ".log";
$cmd .= " -e md.part" . $nrun . ".edr";
$cmd .= " -cpo md.part" . $nrun . ".cpt";
$cmd .= " > " . $logname . " 2>&1";
syscatch($cmd);
errorscan($logname)==0 or finish "$!";
reporter "logfile = " . $logname;
reporter timestring;

#---job complete
reporter timestring;
if (-e "md.part" . $nrun . ".gro") { reporter "completed"; }

