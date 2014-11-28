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

reporter("SIMULATION CONTINUATION");
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

#---CONTINUE PRODUCTION MD RUN
reporter "running: tpbconv, production MD run continuation";
$logname = "log-tpbconv-md-" . $nrun;
$cmd = $config{"tpbconv"};
$cmd .= " -s md.part" . $prun . ".tpr";
$cmd .= " -o md.part" . $nrun . ".tpr";
$cmd .= " -extend " . $extend_time;
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
$cmd .= " -deffnm md";
$cmd .= " -noappend";
$cmd .= " -cpi md.part" . $prun . ".cpt";
$cmd .= " -cpo md.part" . $nrun . ".cpt";
$cmd .= " > " . $logname . " 2>&1";
syscatch($cmd);
errorscan($logname)==0 or finish "$!";
reporter "logfile = " . $logname;
reporter timestring;

#---job complete
reporter timestring;
if (-e "md.part" . $nrun . ".gro") { reporter "completed"; }

