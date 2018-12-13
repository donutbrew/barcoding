#!/usr/bin/env perl
#########
# -*- mode: cperl; tab-width: 8; indent-tabs-mode: nil; basic-offset: 2 -*-
# vim:ts=8:sw=2:et:sta:sts=2
#
# Copyright (c) 2014 Oxford Nanopore Technologies Ltd.
#
# Author:        dturner
# Last Modified: $12Feb15$
# Id:            $Id$
# $HeadURL$
#
# Usage: ./split_barcodes <filename> <stringency> (default = 13, lower = more stringent)
#
# Modified by: Clint Paden <cpaden@cdc.gov>, 30 Mar 2017
#    - Uses split and Parallel::ForkManager to speed things up.
#    - Added ability to read in barcode lists
#    - added --threads option for parallelization
#    - Modified default stringency
# Modified by Clint Paden Nov 2018
#    - Added require_both option
#    - Added check_hybrid and enforce_orientation options

# use strict;
use warnings;
use Text::LevenshteinXS qw(distance);
use Bio::SeqIO;
use Carp;
use English qw(-no_match_vars);
use Readonly;
use Getopt::Long;
use Parallel::ForkManager;
use File::Basename;
use File::Path;
use POSIX;

our $VERSION = q[1.0.0];

Readonly::Scalar my $SEQ_LENGTH   => 200;
Readonly::Scalar my $MATCH_WITHIN => 150;
Readonly::Scalar my $STRINGENCY   => 6; #cp change 13->6

my $opts = {};
GetOptions( $opts, qw(help barcodes=s threads=i stringency=i keep require_both enforce_orientation check_hybrid verbose) );

if ( $opts->{help} ) {	
    print <<"EOT" or croak "Error printing: $ERRNO";
Usage: \n$PROGRAM_NAME  --barcodes <barcodes-file> [--threads <threads>] [--stringency <stringency> (default = $STRINGENCY, lower = more stringent)] [--require_both] [--check_hybrid] [--enforce_orientation] <FASTQ/A file>
EOT
    exit;
}

### Load barcodes from file.
my @barcodes;
my %bclist;
if ( $opts->{barcodes} ) {
    
    if ($opts->{barcodes} =~ /\.(fa|fasta)$/) {
		my $bc_in        = Bio::SeqIO->new(
	        -file   => $opts->{barcodes},
	        -format => "Fasta",
	    ) or croak "Could not open barcode fasta file '$opts->{barcodes}' $ERRNO\n";
	    while ( my $b = $bc_in->next_seq() ) {
	        $bclist{$b->seq()} = $b->id();
	        push @barcodes, $b->seq();
	    }

    }
    else {	
	    open BCFILE, "$opts->{barcodes}";
	        # @barcodes = split /\s+/smx, <BCFILE>;
	        while (<BCFILE>) {
	        chomp;
	        push @barcodes, $_;
	        $bclist{$_} = $_;
	    }    
	    close BCFILE;
	}
}
else {
    # @barcodes   = split /\s+/smx, <DATA>; # this broke for some reason
    while (<DATA>) { chomp; push @barcodes, $_; $bclist{$_} = $_; }
}

my $infile     = $ARGV[0];
my $stringency = $opts->{stringency} // $STRINGENCY;
my $threads    = $opts->{threads} // 1;                #cp
local $RS = undef;
my ($basename) = $infile =~ m{([^\\/]+)$}smx;
$basename =~ s{[.]\S+$}{}smx;
my @tempfiles;
my $basedir = dirname($infile);

my %masterct;    #keeping track of number of things
my $filetype = "Fasta";
$filetype = "Fastq" if $infile =~ /(\.fastq$|\.fq$)/;
my $ft = tolower($filetype);
print STDERR "$filetype detected in $infile\n";

## Clint - split file

my $cmd      = "wc -l < $infile";
my $numlines = `$cmd`;
my $linesper = $numlines / $threads;


$linesper = ceil($linesper);                           #round up
while ( $linesper % 4 != 0 ) { $linesper++; }
print STDERR "Splitting into files with $linesper lines\n" if ($opts->{verbose});
$cmd = "split -l $linesper $infile temp.";   ### TO DO: WRITE TO LOCAL STORAGE
system($cmd);
opendir( my $fdir, $basedir );
@tempfiles = grep /^temp/, readdir($fdir);




$threads = 0 if ($threads == 1);  #turns off parallel for debugging
my $fork = new Parallel::ForkManager($threads);
###
for my $infile (@tempfiles) {    #cp generate new infile based on temp.*

    $fork->start and next;
    my $cmd = "mkdir $infile.d";
    system($cmd);
    my $basename = "$infile.d/$basename";
    for my $bc (@barcodes) {
        unlink "$basename-$bclist{$bc}.fasta";
    }

    #########
    # precompute a few things
    #
    my $bc_lengths  = {};
    my $revcomps    = {};
    my $out_handles = {};

    for my $bc (@barcodes) {
        $bc_lengths->{$bc} = length $bc;

        my $revcomp_bc = scalar reverse $bc;
        $revcomp_bc =~ tr/[A,T,G,C]/[T,A,C,G]/;
        $revcomps->{$bc} = $revcomp_bc;
        
        $out_handles->{$bc} = Bio::SeqIO->new(
            -file   => ">>$basename-$bclist{$bc}.$ft",
            -format => "$filetype",
            -width  => '999999999',
        );
    }
    $out_handles->{"unknown"} = Bio::SeqIO->new(
            -file   => ">>$basename-unknown.$ft",
            -format => "$filetype",
            -width  => '999999999',
        );

    my $sequences_in = 0;
    my $count        = {};
    my $io_in        = Bio::SeqIO->new(
        -file   => $infile,
        -format => "$filetype",
    ) or croak "Could not open '$infile' $ERRNO\n";
    while ( my $rec = $io_in->next_seq() ) {
        my $header  = $rec->id();
        my $seq     = $rec->seq();
        my $seq_len = length $seq;
        $sequences_in++;

        my %embedded_bc;

        if ( $seq_len <= $SEQ_LENGTH ) {
            #########
            # sequence too short to bother
            #
            next;
        }

        my $h_min_distance = $SEQ_LENGTH;    # arbitrary high number of changes
        my $t_min_distance = $SEQ_LENGTH;
        my $h_min_bc;
        my $t_min_bc;

        #########
        # iterate over the barcodes we have
        #
        for my $bc (@barcodes) {
            my $bc_length  = $bc_lengths->{$bc};
            my $revcomp_bc = $revcomps->{$bc};

            #########
# only try and match within the first $MATCH_WITHIN bases of the target sequence
#

            ###### Match both ends
            if ($opts->{require_both}) {
	            for my $scan ( 0 .. $MATCH_WITHIN ) {
                    # Note: this will actually match within $MATCH_WITHIN + $bc_lengt
	                my $head = substr $seq, $scan, $bc_length;
	                my $tail = substr $seq, $seq_len - $bc_length - $scan, $bc_length;
	            	my @searchbc = ($bc, $revcomp_bc) ;
                    @searchbc = ($revcomp_bc) if ($opts->{enforce_orientation});
	            	for my $barcode ( @searchbc ) {   # ( $bc, $revcomp_bc ) { # Mod to only look for proper barcode from 5'
	                    my $distance = distance( $barcode, $head );

	                    if ( $distance < $h_min_distance && $distance <= $stringency)
	                    {    ## no critic (ProhibitDeepNests)
	                        $h_min_distance = $distance;
	                        $h_min_bc       = $bc;
	                    }
	                    elsif ( $distance == $h_min_distance && $h_min_bc && $bc ne $h_min_bc)  #case where we can't tell one BC from another---> Issue is that as you scan, you pick up BC multiple times...
	                    {
	                        print STDERR "SDHBC\t$header\t$bclist{$h_min_bc}($h_min_distance),$bclist{$bc}($distance)q\n";
	                        $h_min_bc       = undef;

	                    }
	                }

                    
                    @searchbc = ($bc, $revcomp_bc);
                    @searchbc = ($bc) if ($opts->{enforce_orientation});
	                for my $barcode ( @searchbc ) { #( $bc, $revcomp_bc ) { # Mod to only loof for RC on tail.
	                    my $distance = distance( $barcode, $tail );

	                    if ( $distance < $t_min_distance && $distance <= $stringency)
	                    {    ## no critic (ProhibitDeepNests)
	                        $t_min_distance = $distance;
	                        $t_min_bc       = $bc;
	                    }
	                    elsif ( $distance == $h_min_distance && $t_min_bc && $bc ne $t_min_bc)  #case where we can't tell one BC from another
	                    {
	                        print STDERR "SDTBC\t$header\t$bclist{$t_min_bc}($t_min_distance),$bclist{$bc}($distance)\n";
                            $t_min_bc       = undef;
	                    }
	                }


	            }
        	} else {
	            for my $scan ( 0 .. $MATCH_WITHIN ) {
	                my $head = substr $seq, $scan, $bc_length;
	                my $tail = substr $seq, $seq_len - $bc_length - $scan, $bc_length;

	                for my $window ( $head, $tail ) {
	                    for my $barcode ( $bc, $revcomp_bc ) {
	                        my $distance = distance( $barcode, $window );

	                        if ( $distance < $h_min_distance && $distance <= $stringency)
	                        {    ## no critic (ProhibitDeepNests)
	                            $h_min_distance = $distance;
	                            $h_min_bc       = $bc;
	                        }
	                        elsif ( $distance == $h_min_distance && $h_min_bc && $bc ne $h_min_bc)  #case where we can't tell one BC from another
		                    {
		                        print STDERR "MBC\t$header\t$bclist{$h_min_bc}($h_min_distance),$bclist{$bc}($distance)\n";
                                $h_min_bc       = undef;
		                        
		                    }
	                    }
	                }
	            }

        	}
            if ($opts->{check_hybrid}) {
                @searchbc = ($bc, $revcomp_bc); 
                for my $scan (($MATCH_WITHIN + $bc_length) .. ($seq_len - $MATCH_WITHIN)) {
                    my $window = substr $seq, $scan, $bc_length;
                    for my $barcode ( @searchbc ) { 
                        my $distance = distance( $barcode, $window );
                        if ( $distance <= $stringency ) {
                            $embedded_bc{$bclist{$bc}} = 1;
                        }

                    }
                }
                # if (keys %embedded_bc) {
                #     print STDERR "$header\t$bclist{$h_min_bc}\t$bclist{$t_min_bc}\t";
                #     for my $i (sort keys %embedded_bc) { print STDERR "$i,";}
                #     print STDERR "\n";
                # }
            }    
        } #end for

        ### Write out 
        $bclist{"unknown"} = "unknown";
        $h_min_bc ||= "unknown"; 
        $t_min_bc ||= "unknown";
        my $exclude = scalar(%embedded_bc);
            print STDERR "$header   $bclist{$h_min_bc} ($h_min_distance) : $bclist{$t_min_bc} ($t_min_distance)  exclude=$exclude\n" if $opts->{verbose};
        if ((  $opts->{require_both} && $h_min_bc && $t_min_bc && $h_min_distance <= $stringency  && $t_min_distance <= $stringency && $h_min_bc eq $t_min_bc && ! $exclude) ||
        	(! $opts->{require_both} && $h_min_bc && $h_min_distance <= $stringency )) {
            my $io_out = $out_handles->{$h_min_bc};
            $io_out->write_seq($rec);
            $count->{$h_min_bc}++;
        } else {
        	my $io_out = $out_handles->{"unknown"};
        	$io_out->write_seq($rec);
            $count->{"unknown"}++;
        }
        if (keys %embedded_bc) {
                    print STDERR "EMB\t$header\t$bclist{$h_min_bc}\t$bclist{$t_min_bc}\t";
                    for my $i (sort keys %embedded_bc) { print STDERR "$i,";}
                    print STDERR "\n";
                }
        
    } #end while

    $io_in->close();

    my $sequences_out = 0;

# open COUNT, ">>", "count_temp";
# for my $key (sort keys %{$count}) {
#   printf COUNT "%s\t%d\n", $key, $count->{$key} or croak qq[Error printing: $ERRNO];
# #   $masterct{$key} += $count->{$key};
#   $sequences_out += $count->{$key};
# }
# close COUNT;
    printf "Sequences in:  %d\nSequences out: %d\n", $sequences_in,
        $sequences_out
        or croak qq[Error printing: $ERRNO];

    $fork->finish;
} #end for
$fork->wait_all_children;
print "All children finished\n";

# open my $tcount, "count_temp";
# while (<$tcount>) {
#   chomp;
#   print $_;
#   my @temp = split "\t", $_;
#   $masterct{$temp[0]} += $temp[1];
# }
# close $tcount;
# my $sequences_out = 0;
# print "======SUMMARY=======\n";
# for my $key (sort keys %masterct) {
#   printf "%s\t%d\n", $key, $masterct{$key};
#   $sequences_out += $masterct{$key};
#  }
#  print "Total Sequences out:\t$sequences_out\n";

##### stitch all files from temp.XX.d
$bclist{"unknown"} = "unknown"; #needed to set outside of parallel
for my $bc (@barcodes, "unknown") {


    open OUT, ">>", "$basename-$bclist{$bc}.$ft";
    for my $infile (@tempfiles) {
        open IN, "$infile.d/$basename-$bclist{$bc}.$ft";

        # select OUT;
        print OUT while <IN>;
        close IN;

    }
    close OUT;

}

# Clean up
if (! $opts->{keep}) {
	for my $infile (@tempfiles) {
    rmtree "$infile.d";
    unlink $infile;
}

unlink "count_temp";
}
1;

__DATA__
AAGAAAGTTGTCGGTGTCTTTGTG
TCGATTCCGTTTGTAGTCGTCTGT
GAGTCTTGTGTCCCAGTTACCAGG
TTCGGATTCTATCGTGTTTCCCTA
CTTGTCCAGGGTTTGTGTAACCTT
TTCTCGCAAAGGCAGAAAGTAGTC
GTGTTACCGTGGGAATGAATCCTT
TTCAGGGAACAAACCAAGTTACGT
AACTAGGCACAGCGAGTCTTGGTT
AAGCGTTGAAACCTTTGTCCTCTC
GTTTCATCTATCGGAGGGAATGGA
CAGGTAGAAAGAAGCAGAATCGGA
