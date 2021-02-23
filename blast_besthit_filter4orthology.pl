use strict ;
use warnings ;
use File::Basename ;

my$dir = $ARGV[0] ;

opendir(HITS, $dir) ;
my@all_files = readdir(HITS) ;

##only keep blast hit file names
my@files = () ;
foreach my$i (@all_files) {if ($i =~ m/_hits.txt/) {push @files, $i ;}}

foreach my$blast (@files) {
  my$sample = $blast ;
  my$output ;
  my%hits ;

  if ($sample =~ m/(.+)_aa_hits.txt/) {$output = $1 . "_aa_filtered_hits.txt" ;}
  elsif ($sample =~ m/(.+)_nt_hits.txt/) {$output = $1 . "_nt_filtered_hits.txt" ;}
  elsif ($sample =~ m/(.+)_tbx_hits.txt/) {$output = $1 . "_tbx_filtered_hits.txt" ;}

  open BLAST, "<$dir/$blast" or die "cannot open $dir/$blast" ;

  while (<BLAST>) {
    if ($_ =~ m/^#/) {next ;}
    chomp ;
    my@hit = split(/\t/, $_) ;

    my$qseqid = $hit[0] ;
    my$sallacc = $hit[1] ;
    my$salltitles = $hit[2] ;
    my$query_coverage = $hit[3] ;
    my$pident = $hit[4] ;
    my$length = $hit[5] ;
    my$mismatch = $hit[6] ;
    my$gapopen = $hit[7] ;
    my$qstart = $hit[8] ;
    my$qend = $hit[9] ;
    my$sstart = $hit[10] ;
    my$send = $hit[11] ;
    my$evalue = $hit[12] ;
    my$bitscore = $hit[13] ;

    # Note for some blast outputs: This matches up with the blast output, but not the requested output format - ssequid was specified and misspelled in some blasts, and results in no output (and no error)
    if (! exists $hits{$qseqid}) {
      if ($length > 100) {
        if ($pident > 30) {
          $hits{$qseqid}{"SUBJECT"} = $sallacc ;
          $hits{$qseqid}{"TITLES"} = $salltitles ;
          $hits{$qseqid}{"QCOV"} = $query_coverage ;
          $hits{$qseqid}{"PIDENTITY"} = $pident ;
          $hits{$qseqid}{"LENGTH"} = $length ;
          $hits{$qseqid}{"MISMATCH"} = $mismatch ;
          $hits{$qseqid}{"GAPS"} = $gapopen ;
          $hits{$qseqid}{"QSTART"} = $qstart ;
          $hits{$qseqid}{"QEND"} = $qend ;
          $hits{$qseqid}{"SSTART"} = $sstart ;
          $hits{$qseqid}{"SEND"} = $send ;
          $hits{$qseqid}{"Evalue"} = $evalue ;
          $hits{$qseqid}{"BITSCORE"} = $bitscore ;
        }
      }
    }

    else {
      # if more than 1 hit, only keep best hit
      if ($evalue < $hits{$qseqid}{"Evalue"}) {
        $hits{$qseqid}{"SUBJECT"} = $sallacc ;
        $hits{$qseqid}{"TITLES"} = $salltitles ;
        $hits{$qseqid}{"QCOV"} = $query_coverage ;
        $hits{$qseqid}{"PIDENTITY"} = $pident ;
        $hits{$qseqid}{"LENGTH"} = $length ;
        $hits{$qseqid}{"MISMATCH"} = $mismatch ;
        $hits{$qseqid}{"GAPS"} = $gapopen ;
        $hits{$qseqid}{"QSTART"} = $qstart ;
        $hits{$qseqid}{"QEND"} = $qend ;
        $hits{$qseqid}{"SSTART"} = $sstart ;
        $hits{$qseqid}{"SEND"} = $send ;
        $hits{$qseqid}{"Evalue"} = $evalue ;
        $hits{$qseqid}{"BITSCORE"} = $bitscore ;
      }
    }
  }

  close BLAST ;

  open OUT, ">$output" ;

  foreach my$query (sort {$a cmp $b} keys %hits) {
    print OUT $query, "\t", $hits{$query}{"SUBJECT"}, "\t", $hits{$query}{"TITLES"}, "\t", $hits{$query}{"QCOV"}, "\t", $hits{$query}{"PIDENTITY"}, "\t", $hits{$query}{"LENGTH"}, "\t", $hits{$query}{"MISMATCH"}, "\t", $hits{$query}{"GAPS"}, "\t", $hits{$query}{"QSTART"}, "\t", $hits{$query}{"QEND"}, "\t", $hits{$query}{"SSTART"}, "\t", $hits{$query}{"SEND"}, "\t", $hits{$query}{"Evalue"}, "\t", $hits{$query}{"BITSCORE"}, "\n" ;
  }

  close OUT ;
}
