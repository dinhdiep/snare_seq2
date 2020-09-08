#!/usr/bin/perl -w

use strict;

my %fragments;
my %uniqueFragments;
my $bam_file = $ARGV[0];

open(BAMFILE, "samtools view $bam_file |") or die ("Error reading $bam_file.\n");
while(my $line = <BAMFILE>){
  chomp($line);
  my @fields = split "\t", $line;
  my ($name, $chr, $start, $len) = ($fields[0], $fields[2], $fields[3], $fields[8]);
  if($fields[1] & 16){ # read is reverse complementary
    next;
  }else{
    next if(abs($len) > 1000 or $len == 0);
  }
  my $end = $start + $len;
  my $frag_info;
  $frag_info = $chr . ":" . $start . ":" . $end . ":" . $len if($start < $end);
  $frag_info = $chr . ":" . $end . ":" . $start . ":" . $len if($start > $end);
  if(!$fragments{$name}){
    $fragments{$name} = $frag_info;
  }else{
    my @cur_info = split ":", $fragments{$name};
    if($len > $cur_info[3]){
      $fragments{$name} = $frag_info;
    }
  }
}
close(BAMFILE);

foreach my $frag_name (keys %fragments){
  my ($cell_id, $umi) = split ":", $frag_name;
  my ($chr, $start, $end, $len) = split ":", $fragments{$frag_name};
  $uniqueFragments{$cell_id}->{$chr}->{$start}->{$end}++;
  delete $fragments{$frag_name};
}
undef %fragments;

foreach my $cell_id (sort keys %uniqueFragments){
  foreach my $chr (keys %{$uniqueFragments{$cell_id}}){
    foreach my $start (sort {$a <=> $b} keys %{$uniqueFragments{$cell_id}->{$chr}}){
      my @ends = keys %{$uniqueFragments{$cell_id}->{$chr}->{$start}};
      foreach my $end (@ends){
        print $chr, "\t", $start, "\t", $end, "\t", $cell_id, "\t", $uniqueFragments{$cell_id}->{$chr}->{$start}->{$end}, "\n";
      }
    }
  }
}
undef %uniqueFragments;
