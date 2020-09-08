#!/usr/bin/perl -w
#
use strict;
use File::Basename;
my $dirname = dirname(__FILE__);
require "$dirname/enumerate_mm_edit.pl";

sub printUsage{
  print "$0 <whitelist_barcodes> <bam_file> <out_file_prefix> <nthreads (pass to samtools)> <dropEst_mode: R1_dT_N6_pairs.txt>\n";
  exit 1;
}

my $whitelist = $ARGV[0];
my $bam_file = $ARGV[1];
my $out_file = $ARGV[2];
my $ncpu = $ARGV[3];
my $dropEst_mode = $ARGV[4];
if(!$dropEst_mode){
  $dropEst_mode = 0;
}
my $max_distance = 2;

if(!$ARGV[3]){
   printUsage();
}

my %barcodes;
my %unplaced_barcodes;
my %dT_N6_pairs;
my $total_dT_reads = 0;
my $total_N6_reads = 0;

if($dropEst_mode){
  open(PAIRS, "$dropEst_mode") || die("Error reading $dropEst_mode");
  while(my $line = <PAIRS>){
    chomp($line);
    my @fields = split /\s+/, $line;
    my $dT = $fields[0];
    my $n6 = $fields[1];
    $dT_N6_pairs{$n6}->{"dT"} = $dT;
    $dT_N6_pairs{$n6}->{"n6"} = $n6;
    $dT_N6_pairs{$dT}->{"dT"} = $dT;
    $dT_N6_pairs{$dT}->{"n6"} = $n6;
  }
  close(PAIRS);
}

open(WHITELIST, "$whitelist") || die("Error reading $whitelist\n");
while(my $line = <WHITELIST>){
  chomp($line);
  my @fields = split " ", $line;
  foreach my $s1 (@fields){
    $barcodes{$s1} = 1;
    my @enum = enumerate($s1, $max_distance);
    foreach my $en (@enum){
      my ($string, $poss) = split "\t", $en;
      push(@{$unplaced_barcodes{$string}}, $poss);
    }
  }
}
close(WHITELIST);

my $m = 0;
open(BAMFILE, "samtools view -H $bam_file |") || die("Error reading $bam_file\n");
open(OUTFILE, ">$out_file.sam") || die("Error writing to $out_file\n");
while(my $header_line = <BAMFILE>){
  print OUTFILE $header_line;
}
close(BAMFILE);

open(BAMFILE, "samtools view $bam_file |") || die("Error reading $bam_file\n");
while(my $sam_line = <BAMFILE>){
  chomp($sam_line);
  $m++;
  if($m%1000000==0){
    print "$m reads in bam\n";
  }
  my @fields = split "\t", $sam_line;
  my ($cell, $read_id) = split ":", $fields[0];
  if($dropEst_mode){
    ($cell, $read_id) = split "#", $fields[0];
  }
  $cell =~ s/_2/\.2/;
  $cell =~ s/_N/\.N/;
  $cell =~ s/_S/\.S/;
  my ($pool_id, $cell_barcode) = split "_", $cell;
  if($dropEst_mode){
    ($pool_id, $cell_barcode) = split "!", $cell;
  }
  my $barcode1 = substr($cell_barcode, 0, 8);
  my $barcode2 = substr($cell_barcode, 8, 8);
  my $barcode3 = substr($cell_barcode, 16, 8);
  #print $barcode1, ",", $barcode2, ",", $barcode3, "\n";
  if(!exists($barcodes{$barcode1}) and exists($unplaced_barcodes{$barcode1})){
    my @Poss = @{$unplaced_barcodes{$barcode1}};
    next if(scalar(@Poss) > 1);
    my ($query, $edit) = split ":", $Poss[0];
    $barcode1 = $query;
  }
  if(!exists($barcodes{$barcode2}) and exists($unplaced_barcodes{$barcode2})){
    my @Poss = @{$unplaced_barcodes{$barcode2}};
    next if(scalar(@Poss) > 1);
    my ($query, $edit) = split ":", $Poss[0];
    $barcode2 = $query;
  }
  if(!exists($barcodes{$barcode3}) and exists($unplaced_barcodes{$barcode3})){
    my @Poss = @{$unplaced_barcodes{$barcode3}};
    next if(scalar(@Poss) > 1);
    my ($query, $edit) = split ":", $Poss[0];
    $barcode3 = $query;
  }
  if(exists($barcodes{$barcode1}) and exists($barcodes{$barcode2}) and exists($barcodes{$barcode3})){
    $cell = $pool_id . "_" . $barcode1 . $barcode2 . $barcode3;
    $fields[0] = $cell . ":" . $read_id;
    if($dropEst_mode){
      next if(!$dT_N6_pairs{$barcode3});
      my $dT = $dT_N6_pairs{$barcode3}->{"dT"};
      my $n6 = $dT_N6_pairs{$barcode3}->{"n6"};
      $total_dT_reads++ if($dT eq $barcode3);
      $total_N6_reads++ if($n6 eq $barcode3);
      $cell = $pool_id . "!" . $barcode1 . $barcode2 . $dT;
      $fields[0] = $cell . "#" . $read_id;
    }
    print OUTFILE join("\t", @fields), "\n";
  }
}
close(BAMFILE);
close(OUTFILE);
print "total reads in bam: ", $m, "\n";
print "total dT reads in bam: ", $total_dT_reads, "\n";
print "total N6 reads in bam: ", $total_N6_reads, "\n";
system("samtools view --threads $ncpu -b $out_file.sam | samtools sort -n --threads $ncpu - > $out_file.bam");
system("rm $out_file.sam");
#print "total unplaced: ", scalar(keys %unplaced_barcodes), "\n";
#print "total barcodes: ", scalar(keys %barcodes), "\n";

