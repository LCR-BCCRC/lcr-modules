#!/usr/bin/env perl
use strict;
use warnings;

# Stream alleleCounts.tab -> minimal VCF for PLINK ROH.
# REF/ALT are faked (A/C). PLINK ROH only uses GT + positions.

# Usage:
#   perl battenberg_normal_to_vcf.pl \
#     --in 07-40648_CLC02055_alleleCounts.tab \
#     --sample 07-40648_CLC02055 \
#     --out 07-40648_CLC02055.normal.fake.vcf \
#     --minDP 12

my %arg = parse_args(@ARGV);
my $in     = required(\%arg, 'in');
my $sample = required(\%arg, 'sample');
my $out    = required(\%arg, 'out');

my $MIN_DP        = exists $arg{minDP} ? int($arg{minDP}) : 12;
my $HOM_MAX_MINOR = exists $arg{homMaxMinor} ? int($arg{homMaxMinor}) : 1;
my $MIN_MINCOUNT  = exists $arg{minMinorCount} ? int($arg{minMinorCount}) : 3;
my $HET_LO        = exists $arg{hetLo} ? $arg{hetLo} + 0.0 : 0.20;
my $HET_HI        = exists $arg{hetHi} ? $arg{hetHi} + 0.0 : 0.80;

open(my $IN,  '<', $in)  or die "ERROR: cannot open $in: $!\n";
my $OUT;
if ($out eq "-") {
    $OUT = *STDOUT;
} else {
    open($OUT, '>', $out) or die "ERROR: cannot write $out: $!\n";
}

# Minimal VCF header
print $OUT "##fileformat=VCFv4.2\n";
print $OUT "##source=alleleCounts_to_fakeREFALT_vcf.pl\n";
print $OUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print $OUT "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n";
print $OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample\n";

# Skip alleleCounts header
my $hdr = <$IN>;
die "ERROR: input seems empty\n" unless defined $hdr;

while (my $line = <$IN>) {
    chomp($line);
    next if $line eq '';

    my ($chr, $pos, $t1, $t2, $n1, $n2) = split(/\t/, $line);
    next unless defined $n2;

    $chr =~ s/^chr//i;
    $n1 = 0 + $n1;
    $n2 = 0 + $n2;
    my $dp = $n1 + $n2;

    # Conservative genotype calls from counts
    my $gt = "./.";
    if ($dp >= $MIN_DP) {
        my $minc = ($n1 < $n2) ? $n1 : $n2;
        my $maxc = ($n1 > $n2) ? $n1 : $n2;
        my $baf  = ($dp > 0) ? ($n1 / $dp) : 0.0;

        # Homozygous if one allele essentially absent
        if ($minc <= $HOM_MAX_MINOR && $maxc >= ($MIN_DP - $HOM_MAX_MINOR)) {
            $gt = ($n1 >= $n2) ? "0/0" : "1/1";
        }
        # Heterozygous if both alleles present + balanced
        elsif ($minc >= $MIN_MINCOUNT && $baf >= $HET_LO && $baf <= $HET_HI) {
            $gt = "0/1";
        }
    }

    # Fake REF/ALT and ID. Use CHR:POS as ID for uniqueness.
    my $id  = "$chr:$pos";
    my $ref = "A";
    my $alt = "C";

    print $OUT join("\t", $chr, $pos, $id, $ref, $alt, ".", "PASS", ".", "GT:DP", "$gt:$dp"), "\n";
}

close($IN);
if ($out ne "-") {
    close($OUT);
}
exit 0;

sub parse_args {
    my @argv = @_;
    my %a;
    while (@argv) {
        my $k = shift @argv;
        die "ERROR: expected --key\n" unless $k =~ /^--(.+)/;
        my $key = $1;
        my $v = shift @argv;
        die "ERROR: missing value for --$key\n" unless defined $v;
        $a{$key} = $v;
    }
    return %a;
}

sub required {
    my ($href, $k) = @_;
    die "ERROR: missing required --$k\n" unless exists $href->{$k};
    return $href->{$k};
}
