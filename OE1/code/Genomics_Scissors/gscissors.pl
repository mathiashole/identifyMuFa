#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Spec;

# -----------------------------------------------------------------------------
# Extractor Script
# -----------------------------------------------------------------------------
# Original Author: Miguel Ponce de Leon
# Modified by: Mathias Mangino
# Additional modifications for format support and help/version features
# -----------------------------------------------------------------------------

# Get command line options
my ($fasta_file, $coord_file, $coord_format, $output_file, $flag_not_to_upper, $help, $version, $log);
GetOptions(
    'fasta=s'       => \$fasta_file,
    'coordinates=s' => \$coord_file,
    'format=s'      => \$coord_format,
    'output=s'      => \$output_file,
    'noupper!'      => \$flag_not_to_upper,
    'log'         => \$log,
    'help|h'        => \$help,
    'version|v'     => \$version,
) or die "Error in command line arguments\n";

# Show help or version if requested
show_help() if $help;
show_version() if $version;

# Check required parameters
die "Error: Missing multifasta file\n" unless $fasta_file;
die "Error: Missing coordinate file\n" unless $coord_file;
die "Error: Missing coordinate format\n" unless $coord_format;
die "Error: Missing output file\n" unless $output_file;

# Determine the log file name based on the output file name if --log is specified
my $log_file;
if ($log) {
    my ($vol, $dir, $file) = File::Spec->splitpath($output_file);
    $file =~ s/\.[^.]+$//;  # Remove the file extension
    $log_file = File::Spec->catpath($vol, $dir, "$file.log");
}

print "🔄 Processing started...\n";
# Main program
extractor($fasta_file, $coord_file, $coord_format, $output_file, $flag_not_to_upper, $log_file);
print "✅ Processing completed successfully. Output saved to $output_file\n";
print "📜 Log file saved to $log_file\n" if $log_file;

sub extractor {
    my ($fasta_file, $coord_file, $coord_format, $output_file, $flag_not_to_upper, $log_file) = @_;

    my %hash_sequence = read_fasta($fasta_file);
    open(my $output_fh, '>', $output_file) or die "\nError: Cannot open output file $output_file: $!\n";
    open(my $log_fh, '>', $log_file) or die "\nError: Cannot open log file $log_file: $!\n" if $log_file;

    my @coordinates = parse_coordinate_file($coord_file, $coord_format);
    foreach my $coord (@coordinates) {
        # my ($name, $start, $end, $strand, $seq_name, @rest) = @$coord;
        # my ($name, $start, $end, $strand, @rest) = @$coord; # Adjust for strand
        if ($coord_format eq 'id') {
            my $name = $coord->[0];

            # Exact match
            if (exists $hash_sequence{$name}) {
                my $sequence = $hash_sequence{$name};
                format_sequence(\$sequence);
                unless ($flag_not_to_upper) {
                    $sequence = uc($sequence);
                }
                print $output_fh ">$name\n$sequence\n";

            # Partial match
            } else {
                my $base_name = $name; 
                $base_name =~ s/[:;].*//;  # Normalize the base ID
                my ($match) = grep { /^$name\b/ } keys %hash_sequence;

                if ($match) {
                    my $sequence = $hash_sequence{$match};
                    format_sequence(\$sequence);
                    # $sequence = uc($sequence) unless $flag_not_to_upper;
                    unless ($flag_not_to_upper) {
                        $sequence = uc($sequence);
                    }
                    print $output_fh ">$match\n$sequence\n";

                # No match found
                } else {
                    warn "Warning: Sequence $name not found in FASTA file\n";
                    print $log_fh "Sequence $name not found in FASTA file\n" if $log_file;
                }
            }

        } else {
            my ($name, $start, $end, $strand, $seq_name, @rest) = @$coord;
            if (exists $hash_sequence{$name}) {
                my $contig = $hash_sequence{$name};
                # my $sequence = extract_sequence($contig, $start, $end);
                my $sequence = extract_sequence($contig, $start, $end, $strand);

                format_sequence(\$sequence);

                unless ($flag_not_to_upper) {
                    $sequence = uc($sequence);
                }

                my $fasta_string = ">$name" . "_" . ($seq_name // '') . " [$start $end] " . join(" ", @rest) . "\n$sequence\n";
                # my $fasta_string = ">$name" . "_" . ($strand // '') . " [$start $end] " . join(" ", @rest) . "\n$sequence\n";
                print $output_fh $fasta_string;
            } else {
                warn "Warning: Sequence $name not found in FASTA file\n";
                print $log_fh "Sequence $name\t$start\t$end not found in FASTA file\n" if $log_file;
            }
        }    
    }

    close($output_fh);
    close($log_fh) if $log_file;
}

sub format_sequence {
    my $sequence_ref = $_[0];
    $$sequence_ref =~ s/(.{80})/$1\n/g;
}

sub read_fasta {
    my ($fasta_file) = @_;
    my %sequences;

    open(my $fasta_fh, '<', $fasta_file) or die "Error opening FASTA file: $!\n";
    my $current_name = '';
    my $sequence = '';

    while (<$fasta_fh>) {
        chomp;
        if (/^>(\S+)/) {
            if ($current_name) {
                $sequences{$current_name} = $sequence;
            }
            $current_name = $1;
            $sequence = '';
        } else {
            $sequence .= $_;
        }
    }

    $sequences{$current_name} = $sequence if $current_name;
    close($fasta_fh);

    return %sequences;
}

sub extract_sequence {

    # return $result;
    my ($sequence, $start_coordinate, $end_coordinate, $strand) = @_;
    $start_coordinate--;
    $end_coordinate--;

    my $length = abs($end_coordinate - $start_coordinate) + 1;
    my $result;

    # Check if the start is less than or equal to the end
    if ($start_coordinate <= $end_coordinate) {
        $result = substr($sequence, $start_coordinate, $length);
    } else {
        # Case in which the start is greater than the end: extract and apply reverse complement
        my $seq_aux = substr($sequence, $end_coordinate, $length);
        $result = reverse_complement($seq_aux);
    }

    # Apply reverse complement if the strand is negative
    if (defined $strand && $strand eq '-') {
        $result = reverse_complement($result);
    }

    return $result;

}

sub parse_coordinate_file {
    my ($coord_file, $format) = @_;
    my @coordinates;

    open(my $coord_fh, '<', $coord_file) or die "\nError: This file does not exist $coord_file: $!\n";
    while (<$coord_fh>) {
        chomp;
        my @fields;
        if ($format eq 'txt') {
            @fields = split(/\s+/, $_);
        } elsif ($format eq 'gff') {
            @fields = (split(/\t/, $_))[0, 3, 4, 6, 8]; # Include strand as [6]
        } elsif ($format eq 'bed') {
            @fields = (split(/\t/, $_))[0, 1, 2, 5]; # Include strand as [5]
        } elsif ($format eq 'blast') {
            @fields = (split(/\t/, $_))[0, 8, 9];
        } else {
            die "Error: Unknown coordinate format $format\n";
        }
        push @coordinates, \@fields;
    }
    close($coord_fh);

    return @coordinates;
}


sub reverse_complement {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}

## Function to show help of the program
sub show_help {
    print << 'HELP';
  🧬 GScissors 🔪 - A Tool for Sequence Data Manipulation

  💻 Usage:
    GScissors [OPTIONS] --fasta <fasta_file> --coordinates <coord_file> --format <format> --output <output_file> [--noupper] [--log <log_file>] [--verbose <level>] [--output_format <format>]

  📄 Options:
    -h, --help        Show this help message and exit.
    -v, --version     Show the version of the program.
    -fasta, --fasta   Specify that your sequence data is in FASTA format.
    -coordinates, --coordinates Specify the coordinate file (TXT, GFF, BED, BLAST).
    -format, --format Specify the format of the coordinate file (txt, gff, bed, blast).
    -output, --output Specify the output file.
    --noupper         Specify if the output sequences should not be converted to uppercase.
    --log             Specify a log file to store IDs not found in the FASTA file.
    --output_format   Specify the format of the output file (fasta, csv, json).

  📂 Arguments:
    --fasta           The sequence in FASTA format.
    --coordinates     The input file in TXT, BED, GFF, or BLAST format.
    --format          The format of the coordinate file.
    --output          The output file in FASTA format.
    --log             The log file to store IDs not found in the FASTA file.

  📝 TXT Table Format:
    Your TXT file should have the following format:
      sequence_ID  start  end
      sequence_ID  start  end  output_ID

  📨 Contact:
    For more information or to report issues, visit:
    https://github.com/mathiashole

  MIT © Mathias Mangino
HELP
    exit;
}

# Function to show the version of the program
sub show_version {
    print "gscissors.pl v2.3.1\n";
    exit;
}