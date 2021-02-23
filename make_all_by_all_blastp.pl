use strict ;
use warnings ;

##Usage: perl /home/shelbirussell/shelbilrussell@gmail.com/lab/projects/symbiont_coalescent/empirical_samples/cs_population_analyses/scripts/rbb/make_all_by_all_blasts.pl /home/shelbirussell/shelbilrussell@gmail.com/lab/projects/symbiont_coalescent/empirical_samples/cs_population_analyses/relative_coalescence_times/genomes_for_rbb/ /home/shelbirussell/shelbilrussell@gmail.com/lab/projects/symbiont_coalescent/empirical_samples/cs_population_analyses/relative_coalescence_times/scripts/rbb/

my$dir = $ARGV[0] ;
my$blast_dir = $ARGV[1] ;

opendir(GENOMES, $dir) ;
my@all_files = readdir(GENOMES) ;

##only keep fasta file names
my@files = () ;
foreach my$i (@all_files) {if ($i =~ m/.faa/) {push @files, $i ;}}

my$n = scalar @files ;
my$total_files = binom($n, 2) ;

print "\nThere are ", $n, " files in ${dir}\n" ;
print "So there will be ", $total_files, " pairwise comparisons for blast\n\n" ;

my%blast_scripts ;

foreach my$i (0..$#files) {
  foreach my$j (0..$#files) {
    ##skip same genome
    if ($files[$i] eq $files[$j]) {next ;}
    my$name1 = $files[$i] ;
    $name1 =~ s/_protein.faa// ;
    my$name2 = $files[$j] ;
    $name2 =~ s/_protein.faa// ;

    my$out_name = "${name1}_vs_${name2}" ;
    my$rev_out_name = "${name2}_vs_${name1}" ;

    ##skip reverse file name
    if ($blast_scripts{$out_name} || $blast_scripts{$rev_out_name}) {next ;}
    $blast_scripts{$out_name} = 1 ;

    ##print blast commands for this species/genome pair to file
    my$job = $out_name . "_rbb_blast" ;
    my$err = $out_name . "_rbb_blast.err" ;
    my$out = $out_name . "_rbb_blast.out" ;

    open JOB, ">${blast_dir}/${job}" or die "cannot open ${blast_dir}/${job}\n" ;

    print JOB "#!/bin/sh
#SBATCH -J ${job}
#SBATCH -o ${out}
#SBATCH -e ${err}
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 2880\n

module load ncbi-blast/2.2.29+-fasrc03

AADB=/n/scratchlfs/cavanaugh_lab/shelbirussell/Wolbachia_genomes/genomes_aa_seqs

blastp -subject \$AADB/$files[$i] -query \$AADB/$files[$j] -max_target_seqs 5 \\
-task blastp -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -evalue 1e-6 -num_threads 4 -show_gis \\
-outfmt \"6 qseqid ssequid sallacc salltitles qcovs pident length mismatch gapopen qstart qend sstart send evalue bitscore\" \\
-out ${out_name}_aa_hits.txt

blastp -subject \$AADB/$files[$j] -query \$AADB/$files[$i] -max_target_seqs 5 \\
-task blastp -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -evalue 1e-6 -num_threads 4 -show_gis \\
-outfmt \"6 qseqid ssequid sallacc salltitles qcovs pident length mismatch gapopen qstart qend sstart send evalue bitscore\" \\
-out ${rev_out_name}_aa_hits.txt" ;

    close JOB ;

    my@command = ("sbatch", ${job}) ;

    system(@command) ;

  }
}

my@blast_scripts = keys %blast_scripts ;
print scalar @blast_scripts, " blast scripts written to file\n" ;

sub binom {
    my( $n, $r ) = @_;
    return unless defined $n && $n =~ /^\d+$/ && defined $r && $r =~ /^\d+$/;
    my $product = 1;
    while( $r > 0 ) {
        $product *= $n--;
        $product /= $r--;
    }
    return $product;
}
