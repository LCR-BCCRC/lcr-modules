#!/usr/bin/perl -w
#########################################################
# Author:	Nancy F. Hansen
# Program:	"Shimmer"
# Function:	shimmer
#               Program to use hypothesis testing to 
#               predict single-nucleotide somatic
#               mutations using Fisher's Exact Test with
#               multiple testing correction.
##########################################################

use strict;
use Getopt::Long;
#use Pod::Usage;
use File::Temp qw /tempfile tempdir /;

use vars qw($VERSION);

our %Opt;
our $R_EXE = 'R';
our $ANNOVAR_EXE = 'annotate_variation.pl';

our @SOM_COLORS = ( {'power' => 0.0, 'color' => '255,0,0'}, # red
                    {'power' => 0.4, 'color' => '255,100,0'}, # orange
                    {'power' => 0.6, 'color' => '255,255,0'}, # yellow
                    {'power' => 0.8, 'color' => '0,255,0'}, # green
                    {'power' => 0.9, 'color' => '0,0,255'} ); # blue

$| = 1; # let's see output right away!

my $Id = q$Id:$;
$VERSION = sprintf "%.4f", substr(q$Rev:10000$, 4)/10000;

my $program_name = $0;
if ($program_name !~ /^\//) { # add path
    my $path = `pwd`;
    chomp $path;
    $program_name = "$path/$program_name";
}
my $short_program_name = $program_name;
$short_program_name =~ s:.*/::; # without path

my $Usage = "Usage: $short_program_name <--region chr1:1000-2000> <--bedfile bed file of regions> <--ref reference fasta> <bam file from normal sample> <bam file from mutated sample>\nFor more information, type \"perldoc $short_program_name\".";

my $printcounts_exe = "printCompCounts";

process_commandline();

($#ARGV == 1)
    or die "$Usage";

my $bam1 = $ARGV[0];
my $bam2 = $ARGV[1];
my $ref_fasta = $Opt{'ref'};
my $region = $Opt{'region'};
my $bedfile = $Opt{'bedfile'};
my $minqual = $Opt{'minqual'}; # disregard any base with quality lower than minqual
my $mapqual = $Opt{'mapqual'}; # disregard any read with mapping quality lower than mapqual
my $min_indel_reads = $Opt{'minindel'}; # disregard any indel without this total number of reads' coverage (in both samples)

if ($Opt{'counts'}) {

    # check files:
    foreach my $file ($ref_fasta, $bam1, $bam2) {
        if (!(-r $file)) {
            die "Can\'t read file $file!\n";
        }
    }
    my $som_file = $Opt{'som_file'};
    my $het_file = $Opt{'het_file'};
    my $indel_file = $Opt{'indel_file'};

    if (!$som_file || !$het_file || !(-e $som_file) || !(-e $het_file)) { # need to generate count files
        # print counts of different alleles at all interesting positions
        print_counts($ref_fasta, $bam1, $bam2, $region, $bedfile, $minqual, $mapqual, $som_file, $het_file, $indel_file, $printcounts_exe);
    }
    # test counts for significance with Fisher's exact test (without mult testing corr)
    test_counts($som_file, $het_file, $indel_file);
}
elsif ($Opt{'bh'}) {

    # apply Benjamini-Hochberg to determine cutoff for significance
    my $max_q = $Opt{'max_q'};
    my $test_fof = $Opt{'test_fof'};
    my $outfile = $Opt{'outfile'};
    my $vsfile = $Opt{'vs_file'};
    my $vcffile = $Opt{'vcf_file'};
    if (!(-r $test_fof)) {
        die "Cannot read file $test_fof to apply B-H correction!\n";
    }
    bh_correct_tests($test_fof, $max_q, $outfile);
    write_vs_file($outfile, $vsfile) if ($vsfile);
    write_vcf_file($outfile, $vcffile) if ($vcffile);
}
elsif ($Opt{'annotate'}) {
    my $annovar_db = $Opt{'annovardb'};
    my $vs_file = $Opt{'vs_file'};
    my $ann_vs_file = $Opt{'outfile'};
    my $buildver = $Opt{'buildver'};
    annotate_variants($annovar_db, $buildver, $vs_file, $ann_vs_file);
}
elsif ($Opt{'covg'}) {

    my $input_file = $Opt{'input'}
        or die "Must specify input file with --input option when running $0 -covg\n";

    my $outdir = $Opt{'outdir'}
        or die "Must specify output directory with --outdir option when running $0 --covg\n";

    calc_power_file($input_file, $outdir, $bam1, $bam2, $ref_fasta, $printcounts_exe);

}
else {
    run_shimmer($program_name, $ref_fasta, $bam1, $bam2, $region, $bedfile, $minqual);
}

## Subroutine to run all of the steps of Shimmer

sub run_shimmer {

    my $shimmer = shift;
    my $ref_fasta = shift;
    my $bam1 = shift;
    my $bam2 = shift;
    my $region = shift;
    my $bedfile = shift;
    my $minqual = shift;

    if ($Opt{'outdir'} && !(-e $Opt{'outdir'})) {
        mkdir $Opt{'outdir'}
            or die "Couldn\'t create $Opt{'outdir'}: $!\n";
    }

    my $tmpdir = $Opt{'outdir'} || tempdir( "run_shimmer_XXXXXX", DIR => '.' );
    if ($tmpdir !~ /^\//) {
        my $pwd = `pwd`;
        chomp $pwd;
        $tmpdir = "$pwd/$tmpdir";
    }
   
    if ((!$Opt{'som_file'}) || (!$Opt{'het_file'})) { 
        my $command = "$shimmer --counts --ref $ref_fasta $bam1 $bam2 --min_som_reads $Opt{min_som_reads} --som_file $tmpdir/som_counts.txt --het_file $tmpdir/het_counts.txt --indel_file $tmpdir/indel_counts.txt";
        $command .= " --region $region" if ($region);
        $command .= " --bedfile $bedfile" if ($bedfile);
        $command .= " --minqual $minqual" if ($minqual);
        $command .= " --mapqual $mapqual" if ($mapqual);
        $command .= " --minindel $min_indel_reads" if ($min_indel_reads);
        $command .= " --testall" if ($Opt{'testall'});
        print "$command\n";
        system($command) == 0
            or die "Failed to run $shimmer --counts!\n";
    }

    if (!$Opt{'skip_tests'}) {
        my $max_q = $Opt{'max_q'};
        open SOMFOF, ">$tmpdir/som_counts.fof"
            or die "Couldn\'t open $tmpdir/som_counts.fof for writing: $!\n";
        my $test_output = $Opt{'som_file'} || "$tmpdir/som_counts.txt";
        $test_output =~ s/\.txt$//;
        $test_output .= '.tests.txt'; 
        print SOMFOF "$test_output\n";
        close SOMFOF;

        my $command = "$shimmer --max_q $max_q --test_fof $tmpdir/som_counts.fof --bh --vs_file $tmpdir/somatic_diffs.vs --vcf_file $tmpdir/somatic_diffs.vcf --outfile $tmpdir/som_counts.bh.txt $bam1 $bam2";
        $command .= " --testall" if ($Opt{'testall'});
        print "$command\n";
        (system("$command") == 0)
            or die "Failed to run $command!\n";

        open INDELFOF, ">$tmpdir/indel_counts.fof"
            or die "Couldn\'t open $tmpdir/indel_counts.fof for writing: $!\n";
        my $indel_test_output = $Opt{'indel_file'} || "$tmpdir/indel_counts.txt";
        $indel_test_output =~ s/\.txt$//;
        $indel_test_output .= '.tests.txt'; 
        print INDELFOF "$indel_test_output\n";
        close INDELFOF;

        (system("$shimmer --max_q $max_q --test_fof $tmpdir/indel_counts.fof --bh --vs_file $tmpdir/somatic_indels.vs --outfile $tmpdir/indel_counts.bh.txt $bam1 $bam2") == 0)
            or die "Failed to run $shimmer --test_fof $tmpdir/indel_counts.fof!\n";

        if ($Opt{'annovardb'}) {
            my $annovardb = $Opt{'annovardb'};
            my $buildver = $Opt{'buildver'};
            (system("$shimmer --annotate --buildver $buildver --annovardb $annovardb --annovar $ANNOVAR_EXE --vs_file $tmpdir/somatic_diffs.vs --outfile $tmpdir/somatic_diffs.ANN.vs $bam1 $bam2")==0)
                or die "Failed to run $shimmer --annotate --buildver $buildver --annovardb $annovardb --annovar $ANNOVAR_EXE --vs_file $tmpdir/somatic_diffs.vs --outfile $tmpdir/somatic_diffs.ANN.vs $bam1 $bam2: $!\n";
            (system("$shimmer --annotate --buildver $buildver --annovardb $annovardb --annovar $ANNOVAR_EXE --vs_file $tmpdir/somatic_indels.vs --outfile $tmpdir/somatic_indels.ANN.vs $bam1 $bam2")==0)
                or die "Failed to run $shimmer --annotate --buildver $buildver --annovardb $annovardb --annovar $ANNOVAR_EXE --vs_file $tmpdir/somatic_indels.vs --outfile $tmpdir/somatic_indels.ANN.vs $bam1 $bam2: $!\n";
        }

    }

    if ($Opt{'power'}) {
        my $run_string = "$shimmer --input $tmpdir/som_counts.bh.txt --outdir $tmpdir --covg --ref $ref_fasta $bam1 $bam2";

        system("$run_string") == 0
            or die "Couldn\'t run $run_string\n";
    }
   
} ## end run_shimmer

sub print_counts {

    my $ref_fasta = shift;
    my $bam1 = shift;
    my $bam2 = shift;
    my $region = shift;
    my $bedfile = shift;
    my $minqual = shift;
    my $mapqual = shift;
    my $som_file = shift;
    my $het_file = shift;
    my $indel_file = shift;
    my $printcounts_exe = shift;

    my $insert = $Opt{'insert'};
    my $min_tumor_reads = $Opt{'min_som_reads'};

    my $ra_het_count_limits = read_min_max_data(); # store limits of different genotypes

    # call mpileup (via the c-script "printCompCounts"), and select sites that are independent, but have enough coverage/diversity to be informative
    my $printcounts_call = "$printcounts_exe -bam1 $bam1 -bam2 $bam2 -fasta $ref_fasta";
    $printcounts_call .= " -region $region" if ($region);
    $printcounts_call .= " -bedfile $bedfile" if ($bedfile);
    $printcounts_call .= " -minqual $minqual" if ($minqual);
    $printcounts_call .= " -mapqual $mapqual" if ($mapqual);
    
    print "Calling $printcounts_call\n";

    open COUNTS, "$printcounts_call | "
        or die "Couldn\'t execute $printcounts_call!\n";

    open SOM, ">$som_file"
        or die "Couldn\'t open $som_file for writing: $!\n";

    open INDELS, ">$indel_file"
        or die "Couldn\'t open $indel_file for writing: $!\n";

    open CNV, ">$het_file"
        or die "Couldn\'t open $het_file for writing: $!\n";

    my ($cur_chr, $cur_pos, $cur_win_start, $cur_win_end, $best_geno, $best_cov, $best_string);

    while (<COUNTS>) {
        chomp;
        if (/^#Indels/) {
            my ($indel, $chr, $pos, $ref, $total1, $total2, $indel_string, $indel1, $indel2) = split /\t/, $_;
            my $nonindel1 = $total1 - $indel1;
            my $nonindel2 = $total2 - $indel2;
            next if ($nonindel1 < 0 || $nonindel2 < 0); # skip ugly regions for now
            my $geno = call_genotype($indel1, $total1, $ra_het_count_limits);
            if (($total1 >= $min_indel_reads) && ($indel1 + $indel2 >= $min_indel_reads)) {
                print INDELS "$indel\t$chr\t$pos\t$ref\t$nonindel1\t$nonindel2\t$indel_string\t$indel1\t$indel2\t$geno\n";
            }
            next;
        }
        my ($chr, $pos, $ref_base, $base1, $normal1_count, $tumor1_count, $base2, $normal2_count, $tumor2_count) = split /\t/, $_;
        $ref_base = uc $ref_base;
        my $total_norm = $normal1_count + $normal2_count;
        my $total_tumor = $tumor1_count + $tumor2_count;
        my $total_alt = ($base2 eq $ref_base) ? $normal1_count + $tumor1_count : $normal2_count + $tumor2_count;
        my $normal_ref = ($base2 eq $ref_base) ? $normal2_count : $normal1_count;
        my $tumor_ref = ($base2 eq $ref_base) ? $tumor2_count : $tumor1_count;
        my $normal_alt = ($base2 eq $ref_base) ? $normal1_count : $normal2_count;
        my $tumor_alt = ($base2 eq $ref_base) ? $tumor1_count : $tumor2_count;
        my $alt_base = ($base2 eq $ref_base) ? $base1 : $base2;
        my $first_base = ($base2 eq $ref_base) ? $base2 : $base1;

        my $geno = call_genotype($normal_alt, $total_norm, $ra_het_count_limits);

        # first check for potential somatic alterations:
        if (($total_norm >= $min_tumor_reads) && ($total_tumor >= $min_tumor_reads) && ($total_alt >= $min_tumor_reads)) {
            print SOM "$chr\t$pos\t$ref_base\t$first_base\t$normal_ref\t$tumor_ref\t$alt_base\t$normal_alt\t$tumor_alt\t$geno\n";
        }

        next if ($geno eq 'und');
        # then check if this position is a potential "best position" in its window:
        if (!$cur_win_start) { 
            $cur_chr = $chr;
            $cur_win_start = $pos;
            $cur_win_end = $pos + $insert;
        }

        if (($chr ne $cur_chr) || ($pos > $cur_win_end)) { # process old window, create new
            print CNV "$best_string";
            $best_string = '';
            $best_geno = '';
            $best_cov = '';

            $cur_chr = $chr;
            $cur_win_start = $pos;
            $cur_win_end = $pos + $insert;
        }

        # replace best string, if appropriate:

        if ((!$best_string) || ($geno eq $best_geno && $total_norm > $best_cov) || ($best_geno ne 'het' && $geno eq 'het')) {
            $best_string = "$chr\t$pos\t$ref_base\t$first_base\t$normal_ref\t$tumor_ref\t$alt_base\t$normal_alt\t$tumor_alt\t$geno\n";
            $best_geno = $geno;
            $best_cov = $total_norm;
        }
    }

    # print final window:
    print CNV "$best_string";

    close COUNTS;
    close CNV;
    close SOM;
    close INDELS;

}

sub read_min_max_data {
    my $ra_min_max = [];

    while (<DATA>) {
        next if ((/^#/) || (/^COV/));
        chomp;
        if (/^(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)$/) {
            my ($cov, $max_hom, $min_het, $max_het, $min_hnr) = ($1, $2, $3, $4, $5);
            $ra_min_max->[$cov] = { 'max_hom' => $max_hom,
                                    'min_het' => $min_het,
                                    'max_het' => $max_het,
                                    'min_hnr' => $min_hnr }; 
        }
    } 
    return $ra_min_max;
}

sub call_genotype {
    my $alt_count = shift;
    my $total_count = shift;
    my $ra_count_limits = shift;

    my $rh_limits;
    if ($total_count > $#{$ra_count_limits}) {
        $rh_limits = {'max_hom' => 17,
                      'min_het' => 0.20*$total_count,
                      'max_het' => 0.80*$total_count };
    }
    else {
        $rh_limits = $ra_count_limits->[$total_count];
    }
    if ($rh_limits) {
        if ($alt_count <= $rh_limits->{'max_hom'}) {
            return 'hom';
        }
        elsif ($alt_count >= $rh_limits->{'min_het'} && 
                $alt_count <= $rh_limits->{'max_het'}) {
            return 'het';
        }
    }
    return 'und';
}

sub test_counts {

    # test counts for significance with Fisher's exact test (without mult testing corr)
    my $som_file = shift;
    my $het_file = shift;
    my $indel_file = shift;

    my $test_opt = ($Opt{'testall'}) ? 'all' : 'hom';

    foreach my $file ($som_file, $het_file, $indel_file) {
        my $type = ($file eq $som_file) ? 'som' :
                   ($file eq $het_file) ? 'het' : 'indel';

        my $r_command_file = "$file.r";
    
        open COM, ">$r_command_file"
            or die "Couldn\'t open $r_command_file for writing: $!\n";
    
        print COM <<"DOC";

library("statmod");
con <- file("$file", "r");
while (length(input <- readLines(con, n=1000)) > 0) {
    for (i in 1:length(input)) {
        line <- input[i];
        linevec <- strsplit(line, split="\\t");
        allele_counts <- c(linevec[[1]][[5]], linevec[[1]][[6]], linevec[[1]][[8]], linevec[[1]][[9]]);
        allele_counts <- as.numeric(allele_counts);
        geno <- linevec[[1]][[10]];
        dim(allele_counts) <- c(2,2);
        if ((("$type" == "som") && ((geno == "hom") || (geno == "und"))) || (("$type"=="som") && ("$test_opt"=="all")) || (("$type" == "het") && (geno == "het")) || (("$type" == "indel") && (geno == "hom"))) {
            exact_result <- fisher.test(allele_counts);
            pvalue <- exact_result["p.value"];
            output <- paste(linevec[[1]][[1]], linevec[[1]][[2]], linevec[[1]][[3]], linevec[[1]][[4]], linevec[[1]][[5]], linevec[[1]][[6]], linevec[[1]][[7]], linevec[[1]][[8]], linevec[[1]][[9]], linevec[[1]][[10]], pvalue, sep=":");
        }
        else {
            output <- paste(linevec[[1]][[1]], linevec[[1]][[2]], linevec[[1]][[3]], linevec[[1]][[4]], linevec[[1]][[5]], linevec[[1]][[6]], linevec[[1]][[7]], linevec[[1]][[8]], linevec[[1]][[9]], linevec[[1]][[10]], "NA", sep=":");
        }
        print(output, quote=FALSE, max.levels=0);
    }
}

DOC
 
        close COM
            or die "Couldn\'t close file $r_command_file: $!\n";
        
        my $r_pipe = "$R_EXE --file=$r_command_file | ";
        open ROUTPUT, "$r_pipe"
            or die "Couldn\'t open pipe to $r_pipe!\n";
       
        my $test_output = $file;
        $test_output =~ s/\.txt$//;
        $test_output .= '.tests.txt'; 

        open TEST, ">$test_output"
            or die "Couldn\'t open $test_output for writing: $!\n";

        my $last = 0;
        my @vals = ();
        while (<ROUTPUT>) {
            if (/^\[1\]\s(.+)$/) {
                 my $output = $1;
                 $output =~ s/:\s*/\t/g;
                 print TEST "$output\n";
            }
        }
        while (<ROUTPUT>) {
            next;
        }

        close TEST;
        close ROUTPUT;
    }
}

# subroutine to apply Benjamini-Hochberg procedure to p-values to obtain q values

sub bh_correct_tests {
    my $test_fof = shift;
    my $max_q = shift;
    my $outfile = shift;

    my @files = ();
    open FILES, "$test_fof"
        or die "Couldn\'t open $test_fof for reading: $!\n";
    while (<FILES>) {
        chomp;
        push @files, $_;
    }
    close FILES;

    my @lines = ();
    my @pvalues = ();
    my $no_tests = 0;
    my $line_index = 0;
    foreach my $test_file (@files) {
        open TESTS, $test_file
            or die "Couldn\'t open $test_file for reading: $!\n";
    
        while (<TESTS>) {
            my $line = $_;
            push @lines, $line;
            chomp $line;
            my @fields = split /\s/, $line;
            my $geno = $fields[$#fields - 1];
            my $pvalue = $fields[$#fields];
            $no_tests++ if (($pvalue ne 'NA') || ($geno eq 'und') || ($Opt{'acctests'}));
            if (($pvalue ne 'NA') && ((!$max_q) || ($pvalue <= $max_q))) {
                push @pvalues, {'pvalue' => $pvalue, 'line_index' => $line_index};
            }
            $line_index++;
        }
        close TESTS;
    }
   
    print "Applying Benjamini-Hochberg correction to somatic change predictions with $no_tests tests.  Results in $outfile.\n"; 

    # now order and assign q-values:
    
    my $index = 1;
   
    my @lines_to_print = ();
    my @sorted_pvalues = sort bypthenna @pvalues;
    foreach my $rh_line (@sorted_pvalues) {
        my $pvalue = $rh_line->{'pvalue'};
        my $qvalue = ($pvalue eq "NA") ? "NA" : $no_tests*$rh_line->{'pvalue'}/$index;
        $qvalue = 1.0 if ($qvalue > 1.0);
        last if (($max_q) && ($qvalue ne "NA") && ($qvalue > $max_q));
    
        $index++;
        $lines[$rh_line->{'line_index'}] =~ s/\n/\t$qvalue\n/;
        my $this_line = $lines[$rh_line->{'line_index'}];
        push @lines_to_print, {'line_index' => $rh_line->{'line_index'}, 'line' => $lines[$rh_line->{'line_index'}]};
    }

    sub bypthenna {
        my $avalue = $a->{'pvalue'};
        my $bvalue = $b->{'pvalue'};
        if ($avalue eq "NA") {
            return 1;
        }
        elsif ($bvalue eq "NA") {
            return -1;
        }
        else {
            return $avalue <=> $bvalue;
        }
    }
    
    # and write out the lines to stdout:
   
    open OUT, ">$outfile"
        or die "Couldn\'t open $outfile for writing: $!\n";
 
    foreach my $line (sort {$a->{'line_index'} <=> $b->{'line_index'}} @lines_to_print) {
        print OUT $line->{'line'};
    }
    close OUT;

} # end bh_correct_tests

sub write_vs_file {
    my $outfile = shift;
    my $vsfile = shift;

    my $indel_flag = ($vsfile =~ /somatic_indels/) ? 1 : 0;

    open SOM, "$outfile"
        or die "Couldn\'t open $outfile for reading: $!\n";

    open VS, ">$vsfile"
        or die "Couldn\'t open $vsfile for writing: $!\n";

    print VS "Index\tChr\tLeftFlank\tRightFlank\tref_allele\tvar_allele\tmuttype\tnormal_covg\ttumor_covg\tnormal_ratio\ttumor_ratio\tq_value\n";
    my $index = 1;
    while (<SOM>) {
        chomp;
        if (!$indel_flag) {
            my ($chr, $pos, $ref, $allele1, $norm1_count, $tumor1_count, $allele2, $norm2_count, $tumor2_count, $gen, $pvalue, $qvalue) = split /\t/, $_;
            $ref = uc $ref;
            my $normal_covg = $norm1_count + $norm2_count;
            my $tumor_covg = $tumor1_count + $tumor2_count;
            my $normal_ratio = $norm2_count/$normal_covg;
            my $tumor_ratio = $tumor2_count/$tumor_covg;
            my $lfe = $pos - 1;
            my $rfs = $pos + 1;
            print VS "$index\t$chr\t$lfe\t$rfs\t$ref\t$allele2\tSNP\t$normal_covg\t$tumor_covg\t$normal_ratio\t$tumor_ratio\t$qvalue\n";
            $index++;
        }
        else { # indel file has different format
            my ($indels, $chr, $pos, $ref, $normref_count, $tumorref_count, $indel_string, $normdiv_count, $tumordiv_count, $geno, $pvalue, $qvalue) = split /\t/, $_;
            my $indel_type = ($indel_string =~ /^\+/) ? 'ins' : 'del';
            $indel_string =~ s/^[+-]\d+//;
            my $indel_length = length($indel_string);
            #my $lfe = ($indel_type eq 'ins') ? $pos : $pos - 1;
            my $lfe = $pos;
            my $rfs = ($indel_type eq 'ins') ? $pos + 1 : $pos + $indel_length + 1;
            my $alt_string = ($indel_type =~ 'del') ? '*' : $indel_string;
            my $ref_string = ($indel_type =~ 'del') ? $indel_string : '*';
            my $normal_covg = $normref_count + $normdiv_count;
            my $tumor_covg = $tumorref_count + $tumordiv_count;
            my $normal_ratio = $normdiv_count/$normal_covg;
            my $tumor_ratio = $tumordiv_count/$tumor_covg;
            print VS "$index\t$chr\t$lfe\t$rfs\t$ref_string\t$alt_string\tINDEL\t$normal_covg\t$tumor_covg\t$normal_ratio\t$tumor_ratio\t$qvalue\n";
            $index++;
        }
    }

    close SOM;
    close VS;

} ## end write_vs_file
  
sub write_vcf_file {
    my $outfile = shift;
    my $vcffile = shift;

    open SOM, "$outfile"
        or die "Couldn\'t open $outfile for reading: $!\n";

    open VCF, ">$vcffile"
        or die "Couldn\'t open $vcffile for writing: $!\n";

    print VCF "##fileformat=VCFv4.1\n";
    my ($sec, $min, $hour, $mday, $mon, $year ) = localtime();
    $mon++;
    $year += 1900;
    printf VCF "##fileDate=%d%02d%02d\n", $year, $mon, $mday;
    print VCF "##source=Shimmer\n";

    # include info for RM flag for repeat-masked sequence:
    print VCF "##INFO=<ID=RM,Number=0,Type=Flag,Description=\"Lower-case reference (probably masked for repeat)\">\n";

    # included genotype id's:
    print VCF "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    #print VCF "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
    print VCF "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";

    # print required eight fields:
    print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n";

    my $index = 1;
    while (<SOM>) {
        chomp;
        my ($chr, $pos, $ref, $allele1, $norm1_count, $tumor1_count, $allele2, $norm2_count, $tumor2_count, $gen, $pvalue, $qvalue) = split /\t/, $_;
        my $normal_covg = $norm1_count + $norm2_count;
        my $tumor_covg = $tumor1_count + $tumor2_count;
        my $normal_ratio = $norm2_count/$normal_covg;
        my $tumor_ratio = $tumor2_count/$normal_covg;
        my $info_flag = ($ref eq uc $ref) ? '.' : 'RM';
        my $qual_score = ($qvalue == 0) ? 999 : int(-10.0*log($qvalue)/log(10.0));
        $qual_score = 0 if ($qual_score < 0);
        $ref = uc $ref;
        print VCF "$chr\t$pos\t\.\t$ref\t$allele2\t$qual_score\t.\t$info_flag\tGT:DP\t0/0:$normal_covg\t0/1:$tumor_covg\n";
        $index++;
    }

    close SOM;
    close VCF;

} ## end write_vcf_file

sub annotate_variants {
    my $annovar_db = shift;
    my $buildver = shift;
    my $vs_file = shift;
    my $ann_vs_file = shift;

    # write Annovar input file:
    open VS, "$vs_file"
        or die "Couldn\'t open $vs_file: $!\n";
  
    my $anv = "$vs_file.annovar.txt";
    open ANV, ">$anv"
        or die "Couldn\'t open $anv for writing: $!\n";

    while (<VS>) {
        chomp;
        my ($index, $chr, $lfe, $rfs, $ref, $var, $muttype, $rest) = split /\t/, $_;
        next if ($index eq 'Index');
        my $pos = $lfe + 1;
        if ($muttype eq 'SNP') {
            print ANV "$chr\t$pos\t$pos\t$ref\t$var\t$index\n";
        }
        elsif ($muttype eq 'INDEL') {
            $ref =~ s/\*/\-/;
            $var =~ s/\*/\-/;
            my ($left, $right);
            if ($ref eq '-' && ($lfe == $rfs - 1)) {
                $left = $lfe + 1;
                $right = $left;
            }
            else {
                $left = $lfe + 1;
                $right = $rfs - 1;
            }
            print ANV "$chr\t$left\t$right\t$ref\t$var\t$index\n";
        }
    }

    close ANV; 
    close VS;

    my $annovar_cmd = "$ANNOVAR_EXE --buildver $buildver --geneanno --dbtype knowngene --hgvs $anv $annovar_db";
    (system("$annovar_cmd") == 0)
        or die "Some problem running $ANNOVAR_EXE!\n";

    # parse output, add columns to VS file:
    my $rh_function = {};
    open VAR, "$anv.variant_function"
        or die "Couldn\'t open $anv.variant_function: $!\n";
    while (<VAR>) {
        chomp;
        my @fields = split /\t/;
        $rh_function->{$fields[7]} = {'loc_type' => $fields[0], 'gene' => $fields[1]};
    }
    close VAR;

    my $rh_cons = {};
    open EXON, "$anv.exonic_variant_function"
        or die "Couldn\'t open $anv.exonic_function: $!\n";
    while (<EXON>) {
        chomp;
        my @fields = split /\t/;
        my $type = $fields[1];
        $type =~ s/\s/_/g;
        $rh_cons->{$fields[8]} = {'type' => $type, 'cons' => $fields[2]};
    }
    close EXON;

    print "About to write annotated file $vs_file!\n";
    # now write annotated VS file with extra columns:
    open VS, "$vs_file"
        or die "Couldn\'t open $vs_file: $!\n";
  
    open ANN, ">$ann_vs_file"
        or die "Couldn\'t open $ann_vs_file for writing: $!\n";

    while (<VS>) {
        chomp;
        my @fields = split /\t/, $_;
        my @new_fields = @fields[0..6];
        if ($new_fields[0] eq 'Index') {
            push @new_fields, qw( type Gene_name loc_type consequence );
        }
        else {
            my $index = $new_fields[0];
            my $type = $rh_cons->{$index}->{'type'} || $rh_function->{$index}->{'loc_type'} || 'NA';
            my $gene = $rh_function->{$index}->{'gene'} || 'NA';
            my $loc_type = $rh_function->{$index}->{'loc_type'} || 'NA';
            my $consequence = $rh_cons->{$index}->{'cons'} || 'NA';
            push @new_fields, ($type, $gene, $loc_type, $consequence);
        }
        push @new_fields, @fields[7..$#fields];
        my $annotated_string = join "\t", @new_fields;
        print ANN "$annotated_string\n";
    }

    close ANN;
    close VS;

} ## end annotate_variants 

sub calc_power_file {
    my $input_file = shift;
    my $outdir = shift;
    my $bam1 = shift;
    my $bam2 = shift;
    my $ref_fasta = shift;
    my $printcounts_exe = shift;

    # call mpileup (via the c-script "printCompCounts"), and record bed regions for somatic and CNA power values.
    my $printcounts_call = "$printcounts_exe -bam1 $bam1 -bam2 $bam2 -fasta $ref_fasta";

    my $perc_sum = 0;
    my $perc_no = 0;
    open SOM, "$input_file"
        or die "Couldn\'t open $input_file: $!\n";
    while (<SOM>) {
        chomp;
        my ($chr, $pos, $ref, $base1, $norm1, $tumor1, $base2, $norm2, $tumor2, $geno) = split /\t/, $_;
        $perc_sum += ($norm1 + $norm2)/($norm1 + $tumor1 + $norm2 + $tumor2);
        $perc_no++;
    }
    close SOM;
    my $tumor_perc_estimate = ($perc_no) ? int(100*$perc_sum/$perc_no*2) : 100;
    $tumor_perc_estimate = 100 if ($tumor_perc_estimate > 100);

    my $purity = int($perc_sum/$perc_no*100);

    my $rh_somatic_power = read_somatic_power($purity);
    open COUNTS, "$printcounts_call | "
        or die "Couldn\'t execute $printcounts_call!\n";

    my ($current_chr, $current_start, $current_pos, $current_color);
    open SOMBED, ">$outdir/somatic_power.bed"
        or die "Couldn\'t open $outdir/somatic_power.bed for writing: $!\n";

    while (<COUNTS>) {
        my ($chr, $pos, $ref_base, $base1, $normal1_count, $tumor1_count, $base2, $normal2_count, $tumor2_count) = split /\t/, $_;
        my $total_norm = $normal1_count + $normal2_count;
        my $total_tumor = $tumor1_count + $tumor2_count;

        my $power_color = power_color($rh_somatic_power, $total_norm, $total_tumor);
        if (!$current_color || !$current_chr || $chr ne $current_chr || $pos != $current_pos + 1 || $power_color ne $current_color) {
            # write out old entry:
            if ($current_color) {
                print SOMBED "$current_chr\t$current_start\t$current_pos\t-\t-\t-\t-\t-\t$current_color\n";
            }
            $current_start = $pos - 1;
            $current_color = $power_color;
        }
        $current_chr = $chr;
        $current_pos = $pos;
    }
    close COUNTS;
    if ($current_color) {
        print SOMBED "$current_chr\t$current_start\t$current_pos\t-\t-\t-\t-\t-\t$current_color\n";
    }

    close SOMBED;
    close COUNTS;

} # end calc_power_file

sub read_somatic_power {
    my $purity = shift;
    my $som_power_file = "/home/nhansen/projects/shimmer/power_graphs/somatic_power_table.txt";
    open SOM, $som_power_file
        or die "Couldn\'t open $som_power_file: $!\n";
    my $rh_som_power = {};
    my $best_purity;
    while (<SOM>) {
        chomp;
        my ($perc_tumor, $normal_reads, $tumor_reads, $power) = split;
        next if ($perc_tumor > $purity);
        if (!$best_purity) {
            $best_purity = $perc_tumor;
        }
        elsif ($perc_tumor != $best_purity) {
            last;
        }
        $rh_som_power->{$normal_reads}->{$tumor_reads} = $power;
        $rh_som_power->{$normal_reads}->{'max'} = $tumor_reads;
        $rh_som_power->{'max'} = $normal_reads;
    }
    close SOM;

    return $rh_som_power;
}

sub power_color {
    my $rh_power = shift;
    my $normal_total = shift;
    my $tumor_total = shift;

    my $normal_tens = 10*int($normal_total/10); 
    my $tumor_tens = 10*int($tumor_total/10);

    $normal_tens = $rh_power->{'max'} if ($normal_tens > $rh_power->{'max'});

    my $power = 0;
    if ($rh_power->{$normal_tens} && $rh_power->{$normal_tens}->{$tumor_tens}) {
        $power = $rh_power->{$normal_tens}->{$tumor_tens};
    }
    elsif ($rh_power->{$normal_tens}) {
        if ($tumor_tens > $rh_power->{$normal_tens}->{'max'}) {
            $power = $rh_power->{$normal_tens}->{'max'};
        }
    }
    my $power_color;
    for (my $i=0; $i<= $#SOM_COLORS; $i++) {
        if ($power >= $SOM_COLORS[$i]->{'power'}) {
            $power_color = $SOM_COLORS[$i]->{'color'};
        }
    }
    return $power_color;
}

sub calculate_depth_norm {
    my $file = shift;

    open DIFFS, $file
        or die "Couldn\'t open $file: $!\n";

    my $sum_ratio = 0;
    my $sumsq_ratio = 0;
    my $total_points = 0;
    my @norm_ratios = ();
    while (<DIFFS>) {
        chomp;
        my @fields = split /\t/, $_;
        my $line = $_;
        my $no_fields = @fields;
        if ($no_fields == 12) { # data line 
            my $pvalue = $fields[10];
    
            if (($pvalue eq 'NA') || ($pvalue < 0.5)) {
                next;
            }

            my $normal_a = $fields[4];
            my $tumor_a = $fields[5];
            my $normal_b = $fields[7];
            my $tumor_b = $fields[8];
            #next if (!($normal_a + $normal_b)); # need some tumor reads to assess anything
            next if (($normal_a + $normal_b) < 100 || ($tumor_a + $tumor_b) < 100); # need some reads to assess anything
            my $this_ratio = ($normal_a + $normal_b) / ($tumor_a + $tumor_b);
            $sum_ratio += $this_ratio;
            $sumsq_ratio += $this_ratio**2;
            push @norm_ratios, $this_ratio;
            $total_points++;
        }
    }
    close DIFFS;
    if (!$total_points) {
        die "Unable to calculate a normalization ratio for depths--too few \"normal\" points!\n";
    }
    my @sorted_norm_ratios = sort {$a <=> $b} @norm_ratios;
    my $median_ratio = $sorted_norm_ratios[int($#sorted_norm_ratios/2)];
    
    my $avg_ratio = $sum_ratio/$total_points;
    print "Calculated average read depth ratio of $avg_ratio, median $median_ratio (normal divided by tumor).\n";

    return $avg_ratio;

} # end calculate_depth_norm

sub process_commandline {
    
    # Set defaults here
    %Opt = ( 
             max_q => 0.05, insert => 300, min_som_reads => 10, minindel => 10
           );
    GetOptions(\%Opt, qw(
                region=s bedfile=s ref=s counts som_file=s indel_file=s
                het_file=s bh vs_file=s vcf_file=s max_q=f test_fof=s
                outfile=s outdir=s input=s plots power covg minqual=i 
                min_som_reads=i minindel=i mapqual=i acctests testall
                insert=i annovar=s annovardb=s buildver=s skip_tests
                annotate help+ version verbose 
               )) || pod2usage(0);
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}); }
    if ($Opt{version}) { die "$0, ", q$Revision: $, "\n"; }

    # argument checking:
    if ($Opt{'annovardb'}) {
        $ANNOVAR_EXE = $Opt{'annovar'} || `which $ANNOVAR_EXE`;
        chomp $ANNOVAR_EXE;
        if (!$ANNOVAR_EXE || !(-e $ANNOVAR_EXE)) {
            die "Cannot find $ANNOVAR_EXE--ignoring opt --annovardb!\n";
        }
        if (!$Opt{'buildver'}) {
            die "Must specify a build version (e.g., hg18) for annovar with --buildver\n";
        }
    }
}

=pod

=head1 NAME

shimmer.pl - call somatic single base changes from matched tumor and normal 
next generation sequences.

=head1 SYNOPSIS

Tally counts of alleles in two BAM files, and write out data for sites with 
unexpectedly large deviations in allele frequencies, controlling the false
discovery rate (FDR) using the Benjamini-Hochberg procedure:

  shimmer.pl [options] <normal_bam_file> <tumor_bam_file> --ref <ref_fasta_file>

Run shimmer.pl -man for a detailed description of options and the output files.

=head1 DESCRIPTION

The script creates a randomly named directory (run_shimmer_XXXXXX) in the 
user's current working directory, and reads through the two BAM files
provided as options with samtools mpileup, recording the normal sample's
genotype, as well as the counts of the two most frequently seen alleles at
each site.  It then uses the R "statmod" library to calculate p-values with
the Fisher's exact test on each site where a minimum threshold of alternate
allele copies are seen (see --min_som_reads option).

Once all p-values have been calculated, shimmer.pl uses the Benjamini-Hochberg 
procedure to report only changes with a false discovery rate below
the specified maximum FDR (see --max_q option).  The single-nucleotide 
variants are reported in VarSifter and VCF formats in the files 
"somatic_diffs.vs" and "somatic_diffs.vcf", respectively.

=head1 INPUT

The first and second arguments to shimmer are the paths of two BAM-formatted 
files of aligned sequencing reads.  These files must be sorted
using samtools prior to running shimmer, and indexed if the --region option
wil be used.

The path of a valid, fasta-formatted file for the reference sequence must
be passed with the option --ref.  This fasta file must have a corresponding 
samtools index file with the same name except for an appended ".fai" if the
--region option will be used.

=head1 OPTIONS

=over 5

=item B<--region> I<chr> or

=item B<--region> I<chr:start-end>

This option specifies a region as a reference entry (chromosome), optionally 
followed by a position range, and causes the program to limit somatic call to
only that region.  By default, the program calls variants in all regions 
that covered by reads in both BAM files.

=item B<--bedfile> I<bedfilename>

This option specifies the path of a BED formatted file containing regions to
be tested for somatic variants.  Limiting regions with this option can increase
power to detect somatic variation by reducing the number of tests performed.

=item B<--ref> I<reference_fasta_file>

This option specifies the reference file to which the reads in the BAM files 
were aligned.  It is a required option.

=item B<--minqual> I<min_base_quality_score>

This option specifies a minimum phred quality score to be required for 
read bases to be included in the counts for the Fisher's exact tests.  By
default, all bases are included.

=item B<--mapqual> I<min_mapping_quality_score>

This option specifies a minimum read mapping quality score to be required for 
a read's bases to be included in the counts for the Fisher's exact tests.  By
default, all reads' bases are included.

=item B<--max_q> I<max_acceptable_FDR>

This option specifies the maximum FDR level to be set for the 
Benjamini-Hochberg procedure for multiple testing correction.  A value of 0
will cause shimmer to report all tests including those with q values equal to
1.  (Default=0.05)

=item B<--annovardb> I<path_to_annovar_db_directory>

This option allows the user to specify an ANNOVAR database
directory against which to annotate variants (see Wang et al, 
"ANNOVAR: functional annotations of genetic variants from 
high-throughput sequencing data". Nucl. Acids Res. 38, 2010).

=item B<--buildver> <hg18, hg19, etc.>

This option is passed directly to ANNOVAR.

=item B<--annovar> I<path_to_annovar>

This option allows the user to specify a particular path to the
"annotate_variants.pl" script from ANNOVAR.  As a default, if 
the annovardb option has been specified (see above), Shimmer 
will call the first copy of "annotate_variation.pl" in the user's
path.

=item B<--outdir> <path_to_output_directory>

This option specifies a directory in which to place result files for this run.  If the 
directory doesn't exist, it will be created.  By default, Shimmer will create a randomly named
directory called "run_shimmer_XXXXXX" within the current working directory (where "XXXXXX"
is a random string of length 6).

=item B<--testall>

This option causes the Fisher's Exact Test to be performed on all sites with the specified
minimum number of alternate alleles, not just sites at which the normal sample has homozygous
reference genotype.

=back

=head1 OUTPUT

The single nucleotide variants (sSNVs) are written in both VarSifter and VCF
formats.

The fields in the VarSifter file (see Teer et al., "VarSifter: Visualizing and
analyzing exome-scale sequence variation data on a desktop computer".
Bioinformatics 28, 2012) are as follows:

=over 5

=item B<Index>

A numerical identifier for each variant.

=item B<Chr>

The entry name of the reference sequence in the BAM files and reference fasta file.

=item B<LeftFlank>

Position one base to the left of the variant base.

=item B<RightFlank>

Position one base to the right of the variant base.

=item B<ref_allele>

Reference base.

=item B<var_allele>

Variant base.

=item B<muttype>

Type of mutation.  In this version of shimmer.pl, all variants are of type "SNP".
Future versions of shimmer.pl will also call somatic variants of type "INDEL".

=item B<normal_covg>

Number of reads covering this position in the normal BAM file.  If the
--minqual option has been used, only reads with the required base quality at this position
will be counted.

=item B<tumor_covg>

Number of reads covering this position in the tumor BAM file.  If the
--minqual option has been used, only reads with the required base quality at this position
will be counted.

=item B<normal_ratio>

Ratio of reads with the alternate base to total reads covering this 
position in the normal BAM file.  If the --minqual option has been used, only reads with 
the required base quality at this position will be counted.

=item B<tumor_ratio>

Ratio of reads with the alternate base to total reads covering this 
position in the tumor BAM file.  If the --minqual option has been used, only reads with 
the required base quality at this position will be counted.

=item B<q_value>

The expected value of the false discovery rate if all variants with q_value
values higher than this value were excluded.

=back

The VCF file conforms to the standards in VCFv4.0, but doesn't include as much information
as the VarSifter file does.  In particular, for each somatic SNV predicted, it reports the
reference and alternate allele, the q-value (as described above in the VarSifter file
description), and the genotype (always 0/0 for the normal and 0/1 for the tumor) and depth
of coverage for each sample.  Suggestions are welcome for how to best utilize VCF format
for these data, as it's a work in progress!

=head1 AUTHOR

 Nancy F. Hansen - nhansen@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation. 

=cut

__END__
__DATA__
COV	MAXHR	MINHET	MAXHET	MINHNR
6	NA	NA	NA	NA
7	NA	4	4	NA
8	NA	4	5	NA
9	NA	4	6	NA
10	NA	4	7	NA
11	0	4	7	NA
12	0	4	8	NA
13	0	4	9	NA
14	0	4	10	NA
15	0	5	11	15
16	0	5	12	16
17	0	5	13	17
18	0	5	14	18
19	1	5	14	19
20	1	5	15	20
21	1	6	16	21
22	1	6	17	22
23	1	6	18	23
24	1	6	18	24
25	1	6	19	25
26	1	6	20	26
27	1	7	21	27
28	1	7	22	28
29	1	7	22	29
30	1	7	23	30
31	1	8	24	31
32	1	8	25	31
33	1	8	26	32
34	1	8	26	33
35	1	8	27	34
36	1	9	28	35
37	1	9	29	36
38	1	9	30	37
39	1	9	30	38
40	1	9	31	39
41	1	10	32	40
42	1	10	33	41
43	1	10	33	42
44	1	10	34	43
45	1	11	35	44
46	2	11	36	45
47	2	11	37	46
48	2	11	37	47
49	2	11	38	48
50	2	11	39	49
51	2	12	40	50
52	2	12	41	51
53	2	12	41	52
54	2	12	42	53
55	2	13	43	54
56	2	13	44	55
57	2	13	45	56
58	2	13	45	57
59	2	13	46	58
60	2	14	47	59
61	2	14	48	60
62	2	14	49	60
63	2	14	49	61
64	2	14	50	62
65	2	15	51	63
66	2	15	52	64
67	2	15	52	65
68	2	15	53	66
69	2	16	54	67
70	2	16	55	68
71	2	16	56	69
72	2	16	56	70
73	2	16	57	71
74	2	16	58	72
75	2	17	59	73
76	3	17	60	74
77	3	17	60	75
78	3	17	61	76
79	3	18	62	77
80	3	18	63	78
81	3	18	63	79
82	3	18	64	80
83	3	18	65	81
84	3	19	66	82
85	3	19	67	83
86	3	19	68	84
87	3	19	68	85
88	3	19	69	86
89	3	20	70	87
90	3	20	71	88
91	3	20	71	89
92	3	20	72	89
93	3	21	73	90
94	3	21	74	91
95	3	21	75	92
96	3	21	75	93
97	3	21	76	94
98	3	22	77	95
99	3	22	78	96
100	3	22	79	97
101	3	22	79	98
102	3	22	80	99
103	3	23	81	100
104	3	23	82	101
105	3	23	82	102
106	4	23	83	103
107	4	23	84	104
108	4	24	85	105
109	4	24	86	106
110	4	24	87	107
111	4	24	87	108
112	4	24	88	109
113	4	25	89	110
114	4	25	90	111
115	4	25	90	112
116	4	25	91	113
117	4	26	92	114
118	4	26	93	115
119	4	26	94	116
120	4	26	94	117
121	4	26	95	118
122	4	27	96	118
123	4	27	97	119
124	4	27	98	120
125	4	27	98	121
126	4	27	99	122
127	4	28	100	123
128	4	28	101	124
129	4	28	101	125
130	4	28	102	126
131	4	28	103	127
132	4	29	104	128
133	4	29	105	129
134	4	29	105	130
135	4	29	106	131
136	5	29	107	132
137	5	30	108	133
138	5	30	109	134
139	5	30	109	135
140	5	30	110	136
141	5	31	111	137
142	5	31	112	138
143	5	31	113	139
144	5	31	113	140
145	5	31	114	141
146	5	32	115	142
147	5	32	116	143
148	5	32	117	144
149	5	32	117	145
150	5	32	118	146
151	5	33	119	147
152	5	33	120	147
153	5	33	120	148
154	5	33	121	149
155	5	34	122	150
156	5	34	123	151
157	5	34	124	152
158	5	34	124	153
159	5	34	125	154
160	5	34	126	155
161	5	35	127	156
162	5	35	128	157
163	5	35	128	158
164	5	35	129	159
165	5	36	130	160
166	6	36	131	161
167	6	36	132	162
168	6	36	132	163
169	6	36	133	164
170	6	37	134	165
171	6	37	135	166
172	6	37	136	167
173	6	37	136	168
174	6	37	137	169
175	6	38	138	170
176	6	38	139	171
177	6	38	139	172
178	6	38	140	173
179	6	39	141	174
180	6	39	142	175
181	6	39	143	176
182	6	39	143	176
183	6	39	144	177
184	6	39	145	178
185	6	40	146	179
186	6	40	147	180
187	6	40	147	181
188	6	40	148	182
189	6	41	149	183
190	6	41	150	184
191	6	41	151	185
192	6	41	151	186
193	6	41	152	187
194	6	42	153	188
195	6	42	154	189
196	7	42	155	190
197	7	42	155	191
198	7	42	156	192
199	7	43	157	193
200	7	43	158	194
201	7	43	158	195
202	7	43	159	196
203	7	44	160	197
204	7	44	161	198
205	7	44	162	199
206	7	44	162	200
207	7	44	163	201
208	7	44	164	202
209	7	45	165	203
210	7	45	166	204
211	7	45	166	205
212	7	45	167	205
213	7	46	168	206
214	7	46	169	207
215	7	46	169	208
216	7	46	170	209
217	7	46	171	210
218	7	47	172	211
219	7	47	173	212
220	7	47	174	213
221	7	47	174	214
222	7	47	175	215
223	7	48	176	216
224	7	48	177	217
225	7	48	177	218
226	8	48	178	219
227	8	49	179	220
228	8	49	180	221
229	8	49	181	222
230	8	49	181	223
231	8	49	182	224
232	8	50	183	225
233	8	50	184	226
234	8	50	185	227
235	8	50	185	228
236	8	50	186	229
237	8	51	187	230
238	8	51	188	231
239	8	51	188	232
240	8	51	189	233
241	8	51	190	234
242	8	52	191	234
243	8	52	192	235
244	8	52	193	236
245	8	52	193	237
246	8	52	194	238
247	8	53	195	239
248	8	53	196	240
249	8	53	196	241
250	8	53	197	242
251	8	54	198	243
252	8	54	199	244
253	8	54	200	245
254	8	54	200	246
255	8	54	201	247
256	9	55	202	248
257	9	55	203	249
258	9	55	204	250
259	9	55	204	251
260	9	55	205	252
261	9	56	206	253
262	9	56	207	254
263	9	56	207	255
264	9	56	208	256
265	9	56	209	257
266	9	57	210	258
267	9	57	211	259
268	9	57	211	260
269	9	57	212	261
270	9	57	213	262
271	9	58	214	263
272	9	58	215	263
273	9	58	215	264
274	9	58	216	265
275	9	59	217	266
276	9	59	218	267
277	9	59	219	268
278	9	59	219	269
279	9	59	220	270
280	9	60	221	271
281	9	60	222	272
282	9	60	223	273
283	9	60	223	274
284	9	60	224	275
285	9	61	225	276
286	10	61	226	277
287	10	61	226	278
288	10	61	227	279
289	10	62	228	280
290	10	62	229	281
291	10	62	230	282
292	10	62	230	283
293	10	62	231	284
294	10	62	232	285
295	10	63	233	286
296	10	63	234	287
297	10	63	234	288
298	10	63	235	289
299	10	64	236	290
300	10	64	237	291
301	10	64	238	292
302	10	64	238	292
303	10	64	239	293
304	10	65	240	294
305	10	65	241	295
306	10	65	242	296
307	10	65	242	297
308	10	65	243	298
309	10	66	244	299
310	10	66	245	300
311	10	66	245	301
312	10	66	246	302
313	10	67	247	303
314	10	67	248	304
315	10	67	249	305
316	11	67	249	306
317	11	67	250	307
318	11	67	251	308
319	11	68	252	309
320	11	68	253	310
321	11	68	253	311
322	11	68	254	312
323	11	69	255	313
324	11	69	256	314
325	11	69	256	315
326	11	69	257	316
327	11	69	258	317
328	11	70	259	318
329	11	70	260	319
330	11	70	261	320
331	11	70	261	321
332	11	70	262	321
333	11	71	263	322
334	11	71	264	323
335	11	71	264	324
336	11	71	265	325
337	11	72	266	326
338	11	72	267	327
339	11	72	268	328
340	11	72	268	329
341	11	72	269	330
342	11	72	270	331
343	11	73	271	332
344	11	73	272	333
345	11	73	272	334
346	12	73	273	335
347	12	74	274	336
348	12	74	275	337
349	12	74	275	338
350	12	74	276	339
351	12	74	277	340
352	12	75	278	341
353	12	75	279	342
354	12	75	280	343
355	12	75	280	344
356	12	75	281	345
357	12	76	282	346
358	12	76	283	347
359	12	76	283	348
360	12	76	284	349
361	12	77	285	350
362	12	77	286	350
363	12	77	287	351
364	12	77	287	352
365	12	77	288	353
366	12	78	289	354
367	12	78	290	355
368	12	78	291	356
369	12	78	291	357
370	12	78	292	358
371	12	79	293	359
372	12	79	294	360
373	12	79	294	361
374	12	79	295	362
375	12	79	296	363
376	13	80	297	364
377	13	80	298	365
378	13	80	299	366
379	13	80	299	367
380	13	80	300	368
381	13	81	301	369
382	13	81	302	370
383	13	81	302	371
384	13	81	303	372
385	13	82	304	373
386	13	82	305	374
387	13	82	306	375
388	13	82	306	376
389	13	82	307	377
390	13	83	308	378
391	13	83	309	379
392	13	83	310	379
393	13	83	310	380
394	13	83	311	381
395	13	84	312	382
396	13	84	313	383
397	13	84	313	384
398	13	84	314	385
399	13	84	315	386
400	13	85	316	387
401	13	85	317	388
402	13	85	317	389
403	13	85	318	390
404	13	85	319	391
405	13	86	320	392
406	14	86	321	393
407	14	86	321	394
408	14	86	322	395
409	14	87	323	396
410	14	87	324	397
411	14	87	325	398
412	14	87	325	399
413	14	87	326	400
414	14	88	327	401
415	14	88	328	402
416	14	88	329	403
417	14	88	329	404
418	14	88	330	405
419	14	89	331	406
420	14	89	332	407
421	14	89	332	408
422	14	89	333	408
423	14	90	334	409
424	14	90	335	410
425	14	90	336	411
426	14	90	336	412
427	14	90	337	413
428	14	90	338	414
429	14	91	339	415
430	14	91	340	416
431	14	91	340	417
432	14	91	341	418
433	14	92	342	419
434	14	92	343	420
435	14	92	344	421
436	15	92	344	422
437	15	92	345	423
438	15	93	346	424
439	15	93	347	425
440	15	93	348	426
441	15	93	348	427
442	15	93	349	428
443	15	94	350	429
444	15	94	351	430
445	15	94	351	431
446	15	94	352	432
447	15	95	353	433
448	15	95	354	434
449	15	95	355	435
450	15	95	355	436
451	15	95	356	437
452	15	95	357	437
453	15	96	358	438
454	15	96	359	439
455	15	96	359	440
456	15	96	360	441
457	15	97	361	442
458	15	97	362	443
459	15	97	362	444
460	15	97	363	445
461	15	97	364	446
462	15	98	365	447
463	15	98	366	448
464	15	98	367	449
465	15	98	367	450
466	16	98	368	451
467	16	99	369	452
468	16	99	370	453
469	16	99	370	454
470	16	99	371	455
471	16	100	372	456
472	16	100	373	457
473	16	100	374	458
474	16	100	374	459
475	16	100	375	460
476	16	100	376	461
477	16	101	377	462
478	16	101	378	463
479	16	101	378	464
480	16	101	379	465
481	16	102	380	466
482	16	102	381	466
483	16	102	381	467
484	16	102	382	468
485	16	102	383	469
486	16	103	384	470
487	16	103	385	471
488	16	103	386	472
489	16	103	386	473
490	16	103	387	474
491	16	104	388	475
492	16	104	389	476
493	16	104	389	477
494	16	104	390	478
495	16	105	391	479
496	17	105	392	480
497	17	105	393	481
498	17	105	393	482
499	17	105	394	483
500	17	106	395	484
