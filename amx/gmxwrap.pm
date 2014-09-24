#!/usr/bin/perl

package gmxwrap;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = (
	'reporter',
	'timestring',
	'errorscan',
	'getbashstr',
	'syscatch',
	'fileparse',
	'getmaxforce',
	'finish',
	);

#---error list
my @error_signifiers = (
	'File input/output error:',
	'command not found',
	'Fatal Error:',
	'Can not open file:',
);

#---return a string with a formatted timestamp
sub timestring { 
	my ($mday, $mon, $year, $hr, $min, $sec) = (localtime(time))[3, 4, 5, 2, 1, 0];
	$mon  += 1;
	$year += 1900;
	my $timestring = sprintf("time = %04d.%02d.%02d.%02d%02d.%02d", 
		$year, $mon, $mday, $hr, $min, $sec);
}

#---simultaneously write to stdout and logfile defined in main
sub reporter {
	print $main::logfile $_[0] . "\n";
	print $_[0] . "\n";
}

#---abort the program and write the time to std and logfile defined in main
sub finish {
	if ($_[0] ne "") {
		print $main::logfile $_[0] . "\n";
		print $_[0] . "\n";
	}
	my $endtime = timestring;
	print $main::logfile $endtime . "\n";
	print $endtime . "\n";
	exit;
}

#---scan an input file with path given as argument for typical GROMACS error messages
sub errorscan {
	my $error = 0;
	open(my $infile, '<', @_) or finish "error opening @_: $!";
	while (my $line = <$infile>) {
		foreach my $sig (@error_signifiers) {
			if ($line =~ /$sig/) {
				$line =~ s/\s+$//;
				reporter "caught gmx error: " . $line;
				reporter "fail"
				$error = 1;
			}
		}
	}
	return $error;
}

#---retrieve a bash variable from a bash script path if it is defined in double quotes
sub getbashstr {
	my $filename = shift;
	my $setting_name = shift;
	open(my $infile, '<', $filename) or finish "error opening @_: $!";
	while (my $line = <$infile>) {
		if ($line =~ /$setting_name/) {
			$line =~ /"(.*?)"/;
			return $1
		}
	}
}

#---wrapper for system which detects INT signals
sub syscatch {
	system(@_);
	if ($? & 127) {
	    my $sig = $? & 127;
	     finish "terminating because signal INT" if $sig == 2;
	}
}

#---parse the configuration file into a hash
sub fileparse {
	open(my $parsefile, '<', @_) or finish "error opening @_: $!";
	my %hash;
	while (my $line=<$parsefile>) {   
		$line =~ s/#.*|\s+#.*//;
		$line =~ s/^\s+|\s+$//g;
		if ($line ne "") {
			(my $word1,my $word2) = split /\s+/, $line, 2;
			$hash{$word1} = $word2;
		}
	}
	return %hash;
}

#---get maximum force from minimization log
sub getmaxforce {
	my $maxforce;
	open(my $parsefile, '<', @_) or die "error opening: $!";
	while (my $line=<$parsefile>) {   
		chomp $line;
		if ($line =~ /^Maximum/) {
			$line =~ s/Maximum\s+force\s+=\s+//;
			$line = (split /\s+/, $line)[0];
			#print $line . "\n";
			#$maxforce = $line;
			my $dec = sprintf("%.10g", $line);
			return $dec;
		}
	}
}
