#! /usr/bin/perl 
use strict;
use vars qw(%OPT);
use FindBin qw( $RealBin );
use lib $RealBin;
use Exception;
use Vcf;
use VEP;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use Data::Dumper;
use Cwd 'abs_path';

# Command line arguments
GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "vcf_in=s",
	   "outdir=s",
	   "outfile=s",
	   "cell_annotation_file=s",
	   "cell_type=s",
	   "annotate_conf=s", 
	   "min_varcount_total=i", 
	   "min_varcount_celltype=i", 
	   "min_varportion_total=s", 
	   "min_varportion_celltype=s",
	   "min_datacount_total=i",
	   "min_datacount_celltype=i",
	   "min_dataportion_total=s",
	   "min_allelefreq_mean_total=s",
	   "min_totalportion_celltype=s",
       "zscore_cutoff=s",
	   "only_somatic",
	   "no_run", 
	   "sample_file=s",
	   "control_file=s",
       "chr=s",
       "celltype_sort_column=s",
       "all_sort_column=s",
       "plot_genes=s",
       "rare_cutoff=s",
       "skip_vep",
       "odds_ratio_cutoff=s", 
       "overwrite", 
       "debug"
   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{vcf_in});

=pod

=head1 SYNOPSIS

tapestri.pl 
	-vcf_in output_dir 
	-outdir outdir(default=cwd)
	-outfile output_summary_name(in_outdir)
	-cell_annotation_file cell_annotation_file(tsv_with_format 'sample celltype')
	-cell_type specific_celltype_to_use(default=all)
	-annotate_conf conf_file_if_running_vep_annotation_steps
	-min_varcount_total minimum_number_samples_with_variant
	-min_varcount_celltype minimum_number_samples_in_celltype_with_variant
	-min_varportion_total minimum_portion_total_samples_with_variant
	-min_varportion_celltype minimum_portion_celltype_samples_with_variant
	-min_datacount_total minimum_number_samples_with_variant
	-min_datacount_celltype minimum_number_samples_in_celltype_with_variant
	-min_dataportion_total minimum_portion_total_samples_with_variant
	-min_total_qual minimum_total_quality_score
	-min_total_persample_qual minimum_persample_quality_score
	-min_allelefreq_mean_total min_AF_average_total_cutoff
	-min_totalportion_celltype min_total_portion_of_variants_in_allsample
	-only_somatic only_celltype_specific_variants
	-no_run list_commands_and_just_generate_report 
	-sample_file run_for_sample_subset
	-control_file for_mixed_datasets_where_specific_controls_are_needed 
	-celltype_sort_column columnname_to_sortby_in_celltype_files(assumes_numerical_descending,default=zscore)
	-all_sort_column columnname_to_sortby_in_all_files(assumes_numerical_descending,default=average_qual_per_sample)
	-chr run_for_chr 
	-plot_genes file_of_genes_to_plot(requires_vep_for_priority_variants)
	-rare_cutoff cutoff_for_rare_allele_freq(default=0.02)
	-odds_ratio_cutoff cutoff_for_oddsratio(default=2)
	-zscore_cutoff cutoff_for_zscore(default=1.65)
	-overwrite overwrite_all_existing_files
	-skip_vep don't_rerun_vep_even_with_annotation_set
	-debug report_on_filtering

Required flags: -vcf_in 

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

tapestri.pl -> generate summary files from a tapestri vcf for all clones as well as for each celltype 

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

#All variants 
./quick_annotate_vcf.pl -outdir results/ -outfile patient2_WGS_all_variants  -vcf joint_calls.vcf 

#All variants sorted by specific column
./quick_annotate_vcf.pl -outdir results/ -outfile patient2_WGS_all_variants  -vcf joint_calls.vcf -sort_column average_qual_per_sample

#Variants in 50 samples (or cells)
./quick_annotate_vcf.pl -outdir results/ -outfile patient2_WGS_min50samples  -vcf joint_calls.vcf -min_sample_count 50

#Rerun with new filters on existing dir
./quick_annotate_vcf.pl -skip_vep -no_run -outdir results/ -outfile patient2_WGS_min10samples  -vcf joint_calls.vcf -min_sample_count 10

#Find somatic variants (variants exclusively found in sample_list only)
./quick_annotate_vcf.pl -outdir results/ -outfile celiac4_somatic_2samples -vcf joint_call_chr.vcf -sample_file celiac4_samples_rogue -somatic -min_sample_count 2

#In sample_list but not in controls
./quick_annotate_vcf.pl -outdir results_all/ -outfile patient2_WGS_nodonor  -vcf joint_calls.vcf -control_file sample_control -sample_file samples_of_interest

#Max no_data and min sample count cutoffs applied
./quick_annotate_vcf.pl -outdir results_all/ -outfile patient2_WGS_1nodata_2orMoreSamples  -vcf joint_calls.vcf -control_file sample_control -sample_file sample_noncontrols -max_nocall_count 1 -min_sample_count 2

#MB generate matrix
./quick_annotate_vcf.pl -outdir results_1912/  -outfile 1912_all -vcf 1912.cells.hg38.vcf

#MB for single cluster
./quick_annotate_vcf.pl -vcf_in 1912.cells.hg38.vcf -outdir analysis/ -outfile Celiac_1912_rogue -no_run -mb_matrix -sample_file Rogue_barcodes.txt -count_all -sort_column portion_in_cluster

#MB with all cluster groups and no sample zygosities and 25 samples and sorted by avg_variant_score
./quick_annotate_vcf.pl -vcf_in 1912.cells.hg38.vcf -outdir analysis/ -outfile Celiac_1912_cat1groups_min25samples -group_file group_file_cat1.tsv -no_run -no_zyg -min_sample_count 25 -sort_column average_qual_per_sample



=cut


my $rare_cutoff = defined $OPT{rare_cutoff}?$OPT{rare_cutoff}:0.02;
if ($rare_cutoff > 0.5 && $rare_cutoff <= 0) {
	Exception->throw("ERROR: Rare cutoff freq must be >1 and <0.5");
} 
my $odds_ratio_cutoff = defined $OPT{odds_ratio_cutoff}?$OPT{odds_ratio_cutoff}:-1000;

my $zscore_cutoff = defined $OPT{zscore_cutoff}?$OPT{zscore_cutoff}:0;

my $debug = defined $OPT{debug}?1:0;

my $vcf_file = $OPT{vcf_in};
$vcf_file = abs_path($vcf_file);

if ( !-e $vcf_file ) {
	Exception->throw("File $vcf_file doesn't exist");	
}

my $vcf = Vcf->new(-vcf=>$vcf_file);

#Keep track of total var/sample to allow potential filtering (rerun with reduced sample list)
my %sample_varcount = ();

#For output files
my ($vcf_short) = basename($vcf_file);
(my $vcf_out = $vcf_short) =~ s/.vcf/.txt/;

#Default run all chromosomes
my $chr_filter = defined $OPT{chr}?$OPT{chr}:"all";

#Default of current dir for output
my $outdir = defined $OPT{outdir}?$OPT{outdir}:`pwd`;
chomp $outdir;
$outdir = abs_path($outdir);

my $script_dir = dirname(__FILE__);
my $vep = defined $OPT{annotate_conf}?1:0;

my $only_somatic = defined $OPT{only_somatic}?1:0;


my $cell_anno = defined $OPT{cell_annotation_file}?1:0;
if ( $cell_anno && !-e $OPT{cell_annotation_file} ) {
	Exception->throw("File $OPT{cell_annotation_file} doesn't exist");
}

my $specific_cell_type = defined $OPT{cell_type}?$OPT{cell_type}:'All';

my $overwrite = defined $OPT{overwrite}?1:0;

if (!-d $outdir) {
        `mkdir $outdir`;
}

#For creating matrices for downstream analysis for vatrix and tapestri
my %mb_map = (
				"no_call" => 0,
				"ref" => 1,
				"hom" => 2,
				"het" => 3
			);


my $min_varcount_total = defined $OPT{min_varcount_total}?$OPT{min_varcount_total}:1;
my $min_varcount_celltype = defined $OPT{min_varcount_celltype}?$OPT{min_varcount_celltype}:1;
my $min_varportion_total = defined $OPT{min_varportion_total}?$OPT{min_varportion_total}:0;

if ($min_varportion_total > 1 || $min_varportion_total < 0) {
	Exception->throw("Error: Portion varaiable must be in range 0-1\n");
}

my $min_varportion_celltype = defined $OPT{min_varportion_celltype}?$OPT{min_varportion_celltype}:0;

if ($min_varportion_celltype > 1 || $min_varportion_celltype < 0) {
	Exception->throw("Error: Portion varaiable must be in range 0-1\n");
}

my $min_datacount_total = defined $OPT{min_datacount_total}?$OPT{min_datacount_total}:1;
my $min_datacount_celltype = defined $OPT{min_datacount_celltype}?$OPT{min_datacount_celltype}:1;
my $min_dataportion_total = defined $OPT{min_dataportion_total}?$OPT{min_datacount_total}:0;

if ($min_dataportion_total > 1 || $min_dataportion_total < 0) {
	Exception->throw("Error: Portion varaiable must be in range 0-1\n");
}


my $min_total_qual = defined $OPT{min_total_qual}?$OPT{min_total_qual}:0;
my $min_total_persample_qual = defined $OPT{min_total_persample_qual}?$OPT{min_total_persample_qual}:0;
my $min_odd_ratio = defined $OPT{min_odd_ratio}?$OPT{min_odd_ratio}:0;
my $min_allelefreq_mean_total = defined $OPT{min_allelefreq_mean_total}?$OPT{min_allelefreq_mean_total}:0;
my $min_allelefreq_mean_celltype = defined $OPT{min_allelefreq_mean_celltype}?$OPT{min_allelefreq_mean_celltype}:0;
my $min_totalportion_celltype = defined $OPT{min_totalportion_celltype}?$OPT{min_totalportion_celltype}:0;

#Handling of edge cases

if ($OPT{control_file} && !$OPT{sample_file}) {
	Exception->throw("ERROR: Can't use control option without sample_file\n");
}

if ($OPT{plot_genes} && !$vep) {
	Exception->throw("ERROR: Plot genes only works for priority variants annotated with vep\n");
}


my $vep_wrapper = "$script_dir/vep_wrapper.pl";
my $vep_conf = defined $OPT{vep_conf}?$OPT{vep_conf}:"$script_dir/tapestri.csv";

my @var_types = qw(snv indel);


#Look at specific samples only
my %samples = ();
if ($OPT{sample_file}) {
	open(SAMPLE,"$OPT{sample_file}") || Exception->throw("Can't open file $OPT{sample_file}\n");
	while (<SAMPLE>) {
		chomp;
		my ($sample) = $_ =~ /^(\S+)/;
		$samples{$sample} = 1;
	}
}

#Use specific samples as controls
my %controls = ();
if ($OPT{control_file}) {
	open(CONTROL,"$OPT{control_file}") || Exception->throw("Can't open file $OPT{control_file}\n");
	while (<CONTROL>) {
		chomp;
		my ($sample) = $_ =~ /(\S+)/;
		$controls{$sample} = 1;
	}
}

#Variant groups (i.e. cell annotations)
my %cellanno = ();
my %cellanno_counts = ();

#Generate per gene plots
my $plot = defined $OPT{plot_genes}?1:0;
my %genes_to_plot = ();

if ($plot) {
	open(PLOTGENES,$OPT{plot_genes}) || Exception->throw("Can't open file $OPT{plot_genes}\n");

	while (<PLOTGENES>) {
		chomp;
		my ($gene) = $_ =~ /(\S+)/;
		$genes_to_plot{$gene} = 1;
	}
}

if ($cell_anno) {
	open(ANNOFILE,$OPT{cell_annotation_file}) || Exception->throw("Can't open file $OPT{cell_anno_file}\n");
	while (<ANNOFILE>) {
		chomp;
		my ($sample,$celltype) = split();
		$celltype =~ s/ /_/g;
		$celltype =~ s/\//_/g;
		$cellanno{$sample} = $celltype;
		$cellanno_counts{$celltype}++;
	}
} 
	




#headers
my @common_headers = (
						'chr',
						'start',
						'end',
						'total_sample_quals',
						'average_qual_per_sample',
						'variant_count (het/hom)',
						'ref_count',
						'no_data_count',
						'mean_variant_af',
						'median_variant_af',
						'variant_samples',
						'var_type',
						'ref_base',
						'var_base'
						);


#headers for reporting single cell type
my @single_celltype_headers = (
						'odds_ratio',
						'z-score',
						'CELLTYPE_variant_count',
						'CELLTYPE_variant_freq',
						'nocluster_variant_count',
						'nocluster_variant_freq',
						'CELLTYPE_vs_rest_freq_diff',
						'total_portion_in_CELLTYPE'
						);



my @vep_headers = (			
						'gene_name',
						'ens_gene',
						'ens_trans',
						'dbsnp',
						'gnomad_AF',
						'gmaf',
						'variant_consequence',
						'aa_change',
						'poly_cat',
						'poly_score',
						'sift_cat',
						'sift_score',
						'domain',
						'pubmed',
						'clinical'
						);


#Now parse the files for the final report
my @samples = ();

open(VCF,"$vcf_file") || Exception->throw("Can't open file $vcf_file\n");
while (<VCF>) {
	chomp;
	next unless /^#CHROM/;
	my @fields = split("\t",$_);
	for my $field (@fields) { 
		next if $field eq '#CHROM';
		next if $field eq 'POS';
		next if $field eq 'ID';
		next if $field eq 'REF';
		next if $field eq 'ALT';
		next if $field eq 'QUAL';
		next if $field eq 'FILTER';
		next if $field eq 'INFO';
		next if $field eq 'FORMAT';
		push @samples, $field;
		if (!$cell_anno) {
			#Just create a generic single cell type
			$cellanno{$field} = 'Single_celltype';
			$cellanno_counts{'Single_celltype'}++;
		}
	}
  last;
}

#headers for reporting all celltypes
my @all_celltype_headers = ();

for my $celltype ( sort keys %cellanno_counts ) {
	if ($specific_cell_type ne 'All') {
		next unless $specific_cell_type eq $celltype;
	}
    push @all_celltype_headers, "$celltype var_count (het/hom)", "$celltype var_portion", "$celltype ref_count", "$celltype nodata_count";
}

my @all_headers = my @celltype_headers = ();
if ($vep) {
	@all_headers = (@common_headers, @all_celltype_headers, @vep_headers);
	@celltype_headers = (@common_headers, @single_celltype_headers, @vep_headers);
} else {
	@all_headers = (@common_headers, @all_celltype_headers);
	@celltype_headers = (@common_headers, @single_celltype_headers);
}


my $celltype_sort_column = defined $OPT{celltype_sort_column}?$OPT{celltype_sort_column}:'z-score';
my $all_sort_column = defined $OPT{all_sort_column}?$OPT{all_sort_column}:'average_qual_per_sample';

#Get the index for the sort command later
my $col_index_all = 0;
my $col_index_celltype = 0;


for ( my $colcount = 0 ; $colcount < @celltype_headers ; $colcount++ ) {
    if ($celltype_headers[$colcount] eq $celltype_sort_column) {
	    	$col_index_celltype = $colcount;
	}
}
if ( $col_index_celltype == 0 ) {
	Exception->throw("ERROR: Couldn't find sort_column for celltype \n");
}	

for ( my $colcount = 0 ; $colcount < @all_headers ; $colcount++ ) {
    if ($all_headers[$colcount] eq $all_sort_column) {
	    	$col_index_all = $colcount;
	}
	
}
if ( $col_index_all == 0 ) {
	Exception->throw("ERROR: Couldn't find sort_column for file\n");
}	




my $sample_count = @samples;


#First check and parse the vcf
if ($overwrite || !-e "$outdir/$vcf_out") {
	print STDERR "Normalising vcf...\n";
	$vcf->check_vcf(-vcf_file=>$vcf_file);
	$vcf->parse_vcf(-vcf_file=>$vcf_file);
	$vcf->write_normalised(-vcf_file=>$vcf_file,-vcf_out=>"$outdir/$vcf_out");
}

#Build up command list
my @commands = ();


if ($vep) {
	if (!$OPT{skip_vep}) {
		#Split by type (for vep input)
		push @commands, "grep SNV $outdir/$vcf_out > $outdir/$vcf_out.snv";
		push @commands, "grep -v SNV $outdir/$vcf_out > $outdir/$vcf_out.indel";
		#Generate VEP inputs for SNV/INS/DEL
		my $vep_in_command ="cat $outdir/$vcf_out.snv". ' | sed -e "s/:/ /g" -e "s/;/ /g" -e "s/->/ /" | awk \'{print $1,$2,$3,$7,$8,"+"}'."' > $outdir/$vcf_out.vep.in"; 
		my $vep_indel_command1 = "cat $outdir/$vcf_out.indel". ' | grep DEL |  sed -e "s/:/ /g" -e "s/;/ /g" -e "s/-/ /g" | awk \'{print $1,$2,$3,$8,"-","+"}'."' > $outdir/$vcf_out.vep.indel.in";
		my $vep_indel_command2 = "cat $outdir/$vcf_out.indel". ' | grep INS |  sed -e "s/:/ /g" -e "s/;/ /g" -e "s/+/ /g" -e "s/REF=//" | awk \'{print $1,$2,$3,$15,$15$7,"+"}'."' >> $outdir/$vcf_out.vep.indel.in";
		push @commands, $vep_in_command, $vep_indel_command1, $vep_indel_command2;
		
		
		if ($overwrite) {
			push @commands, "$vep_wrapper -vep_conf $vep_conf -vep_in $outdir/$vcf_out.vep.indel.in -all > $outdir/$vcf_out.vep.indel";
			push @commands, "$vep_wrapper -vep_conf $vep_conf -vep_in $outdir/$vcf_out.vep.in > $outdir/$vcf_out.vep.exon";
			push @commands, "$vep_wrapper -vep_conf $vep_conf -vep_in $outdir/$vcf_out.vep.in -all > $outdir/$vcf_out.vep.all";
		} else {
			if (!-e "$outdir/$vcf_out.vep.indel") {
				push @commands, "$vep_wrapper -vep_conf $vep_conf -vep_in $outdir/$vcf_out.vep.indel.in -all > $outdir/$vcf_out.vep.indel";
			}
			if (!-e "$outdir/$vcf_out.vep.exon") {
				push @commands, "$vep_wrapper -vep_conf $vep_conf -vep_in $outdir/$vcf_out.vep.in > $outdir/$vcf_out.vep.exon";
			}
			if (!-e "$outdir/$vcf_out.vep.all") {
				push @commands,"$vep_wrapper -vep_conf $vep_conf -vep_in $outdir/$vcf_out.vep.in -all > $outdir/$vcf_out.vep.all"
			}
			
		}
		push @commands, "rm -f $outdir/$vcf_out.snv";
		push @commands, "rm -f $outdir/$vcf_out.indel";
			
	}
}



#Run the commands
for my $command (@commands) {
	print STDERR "$command\n";
	`$command` unless $OPT{no_run};
}

my %data = ();

open(PARSED,"$outdir/$vcf_out") || Exception->throw("Can't open file $outdir/$vcf_out\n");

my $line_count = 0;

#multiple allele handling for counting
my %mult_allele = ();
my %total_alleles = (); #Use for generating average score per variant cell


print STDERR "Parsing vcf...\n";
while (<PARSED>) {
    $_ =~ s/^chr//;
	chomp;
    my @fields = split ("\t");
	my ($chr,$start,$end,$data,@genotypes) = split("\t");

  	if ($chr_filter =~ /[0-9X]/) {
    	next unless $chr_filter eq $chr;
  	}

	my ($var_type,$var_base_str,$qual,$allele_count,$zyg_count,$var_allele_total,$mean_af,$median_af,$var_read_count) = $data =~ /([A-Z]+);.*:(\S+);Q=(\S+);AC=(\d+);ZC=(\d+);ALLELE=(\d+).*MEANAF=(\S+);MEDAF=([0-9\.]+);VAR_READ_COUNTS=([0-9\/,]+)/;
	
	if ($var_type !~ /./) {
		print STDERR "ERROR: Data $data\n";
	}
	
	my $var_base = my $ref_base;
	if ($var_type eq 'SNV') {
		($ref_base,$var_base) = split('->',$var_base_str); 
	} else {
		$var_base = $var_base_str;
		$ref_base = 'N/A';
	} 
	my $key = "$chr:$start:$end:$var_base";
	
	

	$data{$key}{var_type} = $var_type;
	$data{$key}{ref} = $ref_base;
	$data{$key}{var} = $var_base; 
	$data{$key}{qual} = $qual;
	$data{$key}{mean_af} = $mean_af =~ /\d/?$mean_af:'N/A';
	$data{$key}{median_af} = $median_af =~ /\d/?$median_af:'N/A';
	$data{$key}{total_ac} = $var_allele_total;

	my $zyg = "N/A";
	my $sample;


	for (my $count = 0; $count < @genotypes; $count++) {
		my @geno_fields = split(':',$genotypes[$count]);
		$sample = defined $samples[$count]?$samples[$count]:0; #Mutect doesn't list samples
		
		#Don't count if sample not included (except with controls or count_all flag)
		if (keys %samples && !keys %controls) {
			next unless exists $samples{$sample} ;
		}
		
		my ($allele1,$allele2);
		if ($geno_fields[0] =~ /\//) {
			($allele1,$allele2) = split('/',$geno_fields[0]);
		} elsif ($geno_fields[0] =~ /\|/) {
			($allele1,$allele2) = split('\|',$geno_fields[0]);
		} else {
			#Exception->throw("ERROR: Can't handle genotype $geno_fields[0]\n");
			next;
		}
		if ($allele1 eq '0' && $allele2 eq '0') {
			$zyg = 'ref';
			$data{$key}{ref_count}++;
			$data{$key}{celltype}{$cellanno{$sample}}{ref_count}++;
		} elsif ($geno_fields[0] eq './.' || $geno_fields[0] eq '.|.') {
			$zyg = 'no_call';
			$data{$key}{no_data_count}++;
			$data{$key}{celltype}{$cellanno{$sample}}{nodata_count}++;
		} elsif ($allele1 == $allele2) {
			if ($allele1 == $zyg_count) {
				$zyg = 'hom';
				$data{$key}{hom_count}++;
				$data{$key}{var_count}++;
				$data{$key}{celltype}{$cellanno{$sample}}{var_count}++;
				$data{$key}{celltype}{$cellanno{$sample}}{hom_count}++;
				if ($sample) {
					push @{$data{$key}{var_samples}},$sample;
					$sample_varcount{$sample}++;
				}
			} else {
				#Handling for second varaint allele not reported
				$zyg = 'ref';
			}
		} elsif ($allele1 != $allele2) {
			if ($zyg_count == $allele1 || $zyg_count == $allele2) {
				$zyg = 'het';
				$data{$key}{het_count}++;
				$data{$key}{var_count}++;
				$data{$key}{celltype}{$cellanno{$sample}}{var_count}++;
				$data{$key}{celltype}{$cellanno{$sample}}{het_count}++;
				if ($sample) {
					push @{$data{$key}{var_samples}},$sample;
					$sample_varcount{$sample}++;
				}
			} else {
				#Handling for second varaint allele not reported
				$zyg = 'ref';
			} 
			
		}  else {
			Exception->throw("ERROR with $genotypes[$count]\n");
		}
		
		$data{$key}{zyg}{$sample} = $zyg;
		
	}
	
	
	
  	my $allele_add = 0;
	if (exists $data{$key}{var_count}){
		$allele_add = $data{$key}{var_count};
	}	
	$total_alleles{"$chr:$start:$end"} += $allele_count;
	$line_count++;

  if ($line_count % 10000 == 0) {
    print STDERR "Parsing vcf $chr $start\n";
  }
}

print STDERR "Parsed vcf...\n";


if ($vep) {
	open(VEPINDEL,"$outdir/$vcf_out.vep.indel") || Exception->throw("Can't open file $outdir/$vcf_out.vep.indel\n");
	
	while (<VEPINDEL>) {
	    $_ =~ s/^chr//;
	    next unless /^[0-9XY]+\s/;
	    
	    chomp;
	    my @fields = split("\t");
	    my $key = $fields[0].':'.$fields[1].':'.$fields[2] .':'.$fields[4];
	    
	    
	    if ($chr_filter =~ /[0-9X]/) {
	        next unless $fields[0] eq $chr_filter;
	    }
	    if (!exists $data{$key}) {
	    	Exception->throw("ERROR: Key $key doesn't exist\n");
	    	next;
	    }
	    $data{$key}{rs} = $fields[5];
	    $data{$key}{gmaf} = $fields[6];
	    $data{$key}{domain} = $fields[7];
	    $data{$key}{pubmed} = $fields[8];
	    $data{$key}{clin} = $fields[9];
	    $data{$key}{exon_str} = $fields[10]; 
	    $data{$key}{ens_gene} = $fields[11];
	    $data{$key}{ens_trans} = $fields[12];
	    $data{$key}{cadd_phred} = $fields[13];
	    $data{$key}{var_consequence} = $fields[14];
	    $data{$key}{gnomad_af} = $fields[15];
	    $data{$key}{genename} = $fields[16];
	}
	
	open(VEPALL,"$outdir/$vcf_out.vep.all") || Exception->throw("Can't open file $outdir/$vcf_out.vep.all\n");
	
	while (<VEPALL>) {
	    $_ =~ s/^chr//;
	    next unless /^[0-9XY]+\s/;
	    
	    chomp;
	    my @fields = split("\t");
	    my $key = $fields[0].':'.$fields[1].':'.$fields[1] .':'.$fields[4];
	    
	    
	    if ($chr_filter =~ /[0-9X]/) {
	        next unless $fields[0] eq $chr_filter;
	    }
	    if (!exists $data{$key}) {
	    	#next;
	    	Exception->throw("ERROR: Key $key doesn't exist\n");
	    }
	    $data{$key}{rs} = $fields[5];
	    $data{$key}{gmaf} = $fields[6];
	    $data{$key}{domain} = $fields[7];
	    $data{$key}{pubmed} = $fields[8];
	    $data{$key}{clin} = $fields[9];
	    $data{$key}{exon_str} = $fields[10]; 
	    $data{$key}{ens_gene} = $fields[11];
	    $data{$key}{ens_trans} = $fields[12];
	    $data{$key}{cadd_phred} = $fields[13];
	    $data{$key}{var_consequence} = $fields[14];
	    $data{$key}{gnomad_af} = $fields[15];
	    $data{$key}{genename} = $fields[16];
	}
	
	
	#print Dumper \%data;
	
	open(VEPEXON,"$outdir/$vcf_out.vep.exon") || Exception->throw("Can't open file $outdir/$vcf_out.vep.exon\n");
	
	while (<VEPEXON>) {
	    $_ =~ s/^chr//;
	    next unless /^[0-9XY]+\s/;
	    chomp;
	    my @fields = split("\t");
	    my $key = $fields[0].':'.$fields[1].':'.$fields[1] .':'.$fields[4];
	    if ($chr_filter =~ /[0-9X]/) {
	      	next unless $chr_filter eq $fields[0];
	   	}
		if (!exists $data{$key}) {
	    	next;
	    	#Exception->throw("ERROR: Key $key doesn't exist\n");
	    }
	    my ($poly_score) = $fields[9] =~ /([0-9\.]+)/;
	    my ($sift_score) = $fields[11] =~ /([0-9\.]+)/;
	    
	    $data{$key}{aa_change} = $fields[5];
	    $data{$key}{ens_gene} = $fields[6];
	    $data{$key}{ens_trans} = $fields[7];
	    $data{$key}{poly_cat} = $fields[8];
	    $data{$key}{poly_score} = $poly_score;
		$data{$key}{sift_cat} = $fields[10];
	    $data{$key}{sift_score} = $sift_score;
	}
	print STDERR "Parsed VEP...\n";
}


my $out = defined $OPT{outfile}?"$outdir/".$OPT{outfile}:"$outdir/${vcf_out}.all.annotated.tsv";

if ($out !~ /tsv$/) {
	$out .= '.tsv';
}



(my $out_short = $out) =~ s/.tsv//;

$out =~ s/.tsv/_all.tsv/;

my $sample_count_file = $out_short."_sample_variant_totals.tsv";
my $sample_group_count = $out_short."_variantcount_by_celltype.tsv";
my $plot_gene_pdf = $out_short."_priority_gene_plots.pdf";

#File containing number of variants per clone
open(VARCOUNT,">$sample_count_file") || Exception->throw("Can't open file to write $sample_count_file \n");
print VARCOUNT join("\t",
    					"Sample",
    					"Varcount",
    					"Celltype"
    					) . "\n";


for my $sample ( keys %sample_varcount ) {
    	print VARCOUNT join("\t",
    							$sample,
    							$sample_varcount{$sample},
    							$cellanno{$sample}
    						) . "\n";
}

#Contains variant count divided by group; write headers
open(GROUP,">$sample_group_count") || Exception->throw("Can't open file to write $sample_group_count \n");
print GROUP join("\t","Coord",sort keys %cellanno_counts) ."\tSamples\n";

#Same as above but only for priority variants (vep included)
if ($vep) {
	my $sample_group_count_priority = $out_short."_variantcount_by_celltype_priority.tsv";
	open(GROUPPRIORITY,">$sample_group_count_priority") || Exception->throw("Can't open file to write $sample_group_count_priority \n");
	print GROUPPRIORITY join("\t","Coord",sort keys %cellanno_counts) ."\n";
}


#If gene plot templates are needed
if ($plot) {
	open(GENES_RSCRIPT,">${out_short}_gene_plot.R") || Exception->throw("Can't open gene plot R script") if $plot; 
	print GENES_RSCRIPT join("\n",
							"library(dplyr)",
							"library(ggplot2)",
							"library(tidyr)",
							'pdf(file="'.$plot_gene_pdf.'")'
							) . "\n";
								
	for my $r_gene (sort keys %genes_to_plot) {
		print GENES_RSCRIPT $r_gene.' <- read.csv2("'.$r_gene.'_cell_count.tsv", header=TRUE,sep="\t")'."\n";
		my $font_size = 5;
		
		print GENES_RSCRIPT $r_gene.' %>% pivot_longer(-Coord,names_to="type") %>% ggplot(aes(Coord,value,fill=type)) + geom_col() + coord_cartesian(ylim=c(0,20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size='.$font_size.')) + ggtitle("'.$r_gene. ' Range_to_20")'."\n";
		print GENES_RSCRIPT $r_gene.' %>% pivot_longer(-Coord,names_to="type") %>% ggplot(aes(Coord,value,fill=type)) + geom_col() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size='.$font_size.')) + ggtitle("'.$r_gene. ' Full_Range")'."\n";
	}
		print GENES_RSCRIPT "dev.off()\n";
    
}


#File for all variants
open(OUT,">$out") || Exception->throw("Can't open file to write\n");

my @all_fhs = (*OUT);

#Create priority file if VEP annotations exist
if ($vep) {
	my $out_priority = $out_short . '_all_priority.tsv';
	open(PRIORITY,">$out_priority") || Exception->throw("Can't open file to write $out_priority\n");
	push @all_fhs,\*PRIORITY;
}

#Print the headers for the full variant files
for my $fh ( @all_fhs ) {
	print $fh join("\t",@all_headers) ."\n\n";
}

#Create all the output files and file handles
my %celltype_fhs = ();
my %celltype_fhs_priority = ();

for my $celltype ( keys %cellanno_counts ) {
	if ($specific_cell_type ne 'All') {
		next unless $specific_cell_type eq $celltype;
	}
	my $celltype_out =  $out_short.'_'.$celltype.'.tsv';
	local *FILE;
	open(FILE,">$celltype_out") || Exception->throw("Can't open file $celltype_out\n");
	push @{$celltype_fhs{$celltype}}, *FILE;
	
	if ($vep) {
		my $celltype_priority =  $out_short.'_'.$celltype.'_priority.tsv';
		local *PRIORITY;
		open(PRIORITY,">$celltype_priority") || Exception->throw("Can't open file $celltype_priority\n");
		push @{$celltype_fhs_priority{$celltype}}, *PRIORITY;
	}   
}



#Update headers for each file handle
for my $celltype ( keys %celltype_fhs ) {
	for my $fh (@{$celltype_fhs{$celltype}}) {
		my @specific_celltype_headers = ();
		for my $col (@celltype_headers) {
			(my $local_col = $col) =~ s/CELLTYPE/$celltype/;
			push @specific_celltype_headers, $local_col;
		}
		print $fh join("\t",@specific_celltype_headers) ."\n\n";
		
	}
	
}





#Matrix that can be input to R for further analysis
my $mb_out = $out_short . '_tapestri.tsv';
open(MB,">$mb_out") || Exception->throw("Can't open file $mb_out\n");

print MB join("\t",
				"ID",
				@samples
				) ."\n";


#Store all the data to write to the various files sorted by the relevant columns
my %all_files_to_write = ();
my %celltype_files_to_write = ();


my @keys = sort { my ($a_chr,$a_coord) = $a =~ /([0-9X]+):(\d+)/; my ($b_chr,$b_coord) = $b =~ /([0-9X]+):(\d+)/; $a_chr cmp $b_chr || $a_coord <=> $b_coord } keys(%data);

for my $key (@keys) {
	my ($chr,$start,$end,$var_base) = split(":",$key);
	if ($chr_filter =~ /[0-9X]/) {
    		next unless $chr_filter eq $chr;
  	}
  	#Flag to decide whether to print to the 'all' file
	my $print_all = 1;
  	next unless $chr =~ /^[0-9X]/;
	my $aa_change = exists $data{$key}{aa_change}?$data{$key}{aa_change}:'NO_AA_CHANGE';
	my $ens_trans = exists $data{$key}{ens_trans}?$data{$key}{ens_trans}:'NO_ENS_TRANS';
	my $ens_gene = exists $data{$key}{ens_gene}?$data{$key}{ens_gene}:'NO_ENS_GENE';
	my $poly_cat = exists $data{$key}{poly_cat}?$data{$key}{poly_cat}:'NO_POLY_CAT';
	my $poly_score = exists $data{$key}{poly_score}?$data{$key}{poly_score}:'NO_POLY_SCORE';
	my $sift_cat = exists $data{$key}{sift_cat}?$data{$key}{sift_cat}:'NO_SIFT_CAT';
	my $sift_score = exists $data{$key}{sift_score}?$data{$key}{sift_score}:'NO_SIFT_SCORE';
	my $rs = !exists $data{$key}{rs} || $data{$key}{rs} eq 'N/A'?'NO_DBSNP':$data{$key}{rs};
	my $gmaf = !exists $data{$key}{gmaf} || $data{$key}{gmaf} eq 'N/A'?'NO_GMAF':$data{$key}{gmaf};	
	my $var_mean_af = !exists $data{$key}{mean_af} || $data{$key}{mean_af} eq 'N/A'?'NO_MEAN_AF':$data{$key}{mean_af};	
	my $var_median_af = !exists $data{$key}{median_af} || $data{$key}{median_af} eq 'N/A'?'NO_MEDIAN_AF':$data{$key}{median_af};	
	my $var_read_count = !exists $data{$key}{var_read_count} || $data{$key}{var_read_count} eq 'N/A'?'NO_VAR_READ_COUNT': $data{$key}{var_read_count};
	my $domain = !exists $data{$key}{domain} || $data{$key}{domain} eq 'N/A'?'NO_DOMAIN':$data{$key}{domain};
	my $pubmed = !exists $data{$key}{pubmed} || $data{$key}{pubmed} eq 'N/A'?'NO_PUBMED':$data{$key}{pubmed};
	my $clin = !exists $data{$key}{clin} || $data{$key}{clin} eq 'N/A'?'NO_CLIN':$data{$key}{clin};
	my $gnomad = !exists $data{$key}{gnomad_af} || $data{$key}{gnomad_af} eq 'N/A'?'NO_GNOMAD':$data{$key}{gnomad_af};
	my $genename = !exists $data{$key}{genename} || $data{$key}{genename} eq 'N/A'?'NO_GENENAME':$data{$key}{genename};
	my $var_consequence = !exists $data{$key}{var_consequence} || $data{$key}{var_consequence} eq 'N/A'?'NO_VARIANT_CONSEQUENCE':$data{$key}{var_consequence};
	my $var_samples;
	
	#Find samples involved 
	if (!exists $data{$key}{var_samples}) {
		#Rare case with nested variant events
		$var_samples = "Complex overlapping event";
	} elsif (@{$data{$key}{var_samples}} > 100) {
		#Don't list more than 100 samples in the report
		$var_samples = ">100 samples";
	} else {
		$var_samples = join(",",@{$data{$key}{var_samples}});
	}
	
	
	#Total variant stats
	my $het_count = exists $data{$key}{het_count}?$data{$key}{het_count}:0;
	my $hom_count = exists $data{$key}{hom_count}?$data{$key}{hom_count}:0;
	my $ref_count = exists $data{$key}{ref_count}?$data{$key}{ref_count}:0;
	my $nd_count = exists $data{$key}{no_data_count}?$data{$key}{no_data_count}:0;
	my $data_count = $het_count + $hom_count + $ref_count;
	my $data_portion = sprintf("%.2f",$data_count/($nd_count+$data_count));
	
	my $var_count = 0;
	if (exists $data{$key}{var_count}) {
		$var_count =  $data{$key}{var_count};
	}
	
	my @celltype_columns = ();
	my %celltype_stats = ();
	my @celltype_vars = ();
	for my $celltype ( sort keys %cellanno_counts ) {
		if ($specific_cell_type ne 'All') {
			next unless $specific_cell_type eq $celltype;
		}
    	#push @group_headers, "$group var_count", "$group ref_count", "$group nodata_count";
    	if (!exists $data{$key}{celltype}{$celltype}) {
    		push @celltype_columns, 0, 0, 0, 0;
    		push @celltype_vars,0;
    	} else {
    		my $celltype_var_count = defined $data{$key}{celltype}{$celltype}{var_count}?$data{$key}{celltype}{$celltype}{var_count}:0;
    		my $celltype_het_count = defined $data{$key}{celltype}{$celltype}{het_count}?$data{$key}{celltype}{$celltype}{het_count}:0;
    		my $celltype_hom_count = defined $data{$key}{celltype}{$celltype}{hom_count}?$data{$key}{celltype}{$celltype}{hom_count}:0;
    		my $celltype_ref_count = defined $data{$key}{celltype}{$celltype}{ref_count}?$data{$key}{celltype}{$celltype}{ref_count}:0;
    		my $celltype_nodata_count = defined $data{$key}{celltype}{$celltype}{nodata_count}?$data{$key}{celltype}{$celltype}{nodata_count}:0;
    		my $celltype_data_count = $celltype_var_count+$celltype_ref_count;
			my $celltype_var_portion = $celltype_data_count != 0?sprintf("%.4f",$celltype_var_count/$celltype_data_count):0;

    		
    		my $other_varcount = $var_count - $celltype_var_count;
    		my $other_ref_count = $ref_count - $celltype_ref_count;
    		
    		
    		my $other_varportion = $data_count - $celltype_data_count != 0?sprintf("%.4f",($var_count - $celltype_var_count)/($data_count - $celltype_data_count)):0;

    		$celltype_stats{$celltype}{$key}{varcount} = $celltype_var_count;
    		$celltype_stats{$celltype}{$key}{varportion} = $celltype_var_portion;
    		$celltype_stats{$celltype}{$key}{datacount} = $celltype_data_count;
    		$celltype_stats{$celltype}{$key}{other_varcount} = $other_varcount;
    		$celltype_stats{$celltype}{$key}{other_varportion} =  $other_varportion;
    		$celltype_stats{$celltype}{$key}{portion_diff} = abs($other_varportion-$celltype_var_portion);
    		my $portion_in_cluster = $var_count != 0?sprintf("%.2f",$celltype_var_count/$var_count):0;
    		$celltype_stats{$celltype}{$key}{portion_in_cluster} = $portion_in_cluster;
    		
    		my $odds_ratio = my $zscore = 'N/A';
    		
    		if ($other_varcount == 0) {
    			$odds_ratio = 'N/A (SOMATIC TO '.$celltype.')';
    		} elsif ($celltype_var_count == 0) {
    			$odds_ratio = 'N/A (NOT IN '.$celltype.')';
    		} elsif ($celltype_ref_count == 0) {
    			$odds_ratio = 'N/A (NO REF IN '.$celltype.')';
    		} elsif ($other_ref_count == 0) {
    			$odds_ratio = 'N/A (NO REF IN OTHER)';
    		} else {
    			$odds_ratio = sprintf("%.2f",($celltype_var_count/$celltype_ref_count)/($other_varcount/$other_ref_count));
	    		
	    		#SE of OR
	    		my $sum_to_OR = 0;
	    		if ($celltype_var_count >0) {
	    			$sum_to_OR += 1/$celltype_var_count;
	    		}
	    		if ($celltype_ref_count >0) {
	    			$sum_to_OR += 1/$celltype_ref_count;
	    		}
	    		if ($other_varcount >0) {
	    			$sum_to_OR += 1/$other_varcount;
	    		}
	    		if ($other_ref_count >0) {
	    			$sum_to_OR += 1/$other_ref_count;
	    		}
	    		
	    		my $standard_OR = sqrt($sum_to_OR);
	    		#95% upper CI
	    		my $CI_lower = $odds_ratio*exp(-1.96*$standard_OR);
	    		#Zscore
	    		$zscore = sprintf("%.2f",log($odds_ratio)/(log($odds_ratio)-log($CI_lower))/1.96);
    		} 
    		$celltype_stats{$celltype}{$key}{zscore} = $zscore;
    		$celltype_stats{$celltype}{$key}{odds_ratio} = $odds_ratio;
	
    		push @celltype_columns, $celltype_var_count ." (${celltype_het_count}/${celltype_hom_count})", $celltype_var_portion, $celltype_ref_count, $celltype_nodata_count;
    		push @celltype_vars,$celltype_var_count;
    	}
	}
	
	#print Dumper \%celltype_stats;
	
	
	my $var_str = $var_count . '('.$het_count . '/'. $hom_count .')';
	my $average_score = 'N/A';

	my $alleles_key = "$chr:$start:$end";

	#First try allele field from vcf
	if (exists $data{$key}{total_ac}) {
		if ($data{$key}{total_ac} > 0) {
				$average_score = sprintf("%.2f",$data{$key}{qual} / $data{$key}{total_ac});
		} elsif (exists $total_alleles{$alleles_key}) {
			#otherwise use allele count from parsing (problems with complex snv/indel mixed events)
			if ($total_alleles{$alleles_key} > 0) {
				$average_score = sprintf("%.2f",$data{$key}{qual} / $total_alleles{$alleles_key});
			}
		}
		
	}
	
	#Start filtering here
	
	#If below min_mean_af
	if ($var_mean_af =~ /\d/ && $var_mean_af <  $min_allelefreq_mean_total) {
		print "$key Filtered ALL Mean_AF $var_mean_af <  $min_allelefreq_mean_total\n" if $debug;
		$print_all = 0;
	}
	
	if ($data_count <  $min_datacount_total) {
		print "$key Filtered ALL Min_datacount $data_count <  $min_datacount_total\n" if $debug;
		$print_all = 0;
	}
	
	if (exists $data{$key}{var_count} && $data{$key}{var_count} < $min_varcount_total) {
		print "$key Filtered ALL Min_varcount $data{$key}{var_count} < $min_varcount_total\n" if $debug;
		$print_all = 0;
	}
	
	
	if (exists $data{$key}{var_count}) {
		my $var_portion = sprintf("%.4f",$var_count/$data_count);
		if ($var_portion <  $min_varportion_total) {
			print "$key Filtered ALL Min_varportion $var_portion <  $min_varportion_total\n" if $debug;
			$print_all = 0;
		}
	}
	
	if ($data{$key}{qual} <  $min_total_qual) {
		print "$key Filtered ALL Min_total_qual $data{$key}{qual} <  $min_total_qual\n" if $debug;
		$print_all = 0;
	}
	
	if ($average_score <  $min_total_persample_qual) {
		print "$key Filtered ALL Min_persample_qual $average_score <  $min_total_persample_qual\n" if $debug;
		$print_all = 0;
	}
	
	if ($data_portion <  $min_dataportion_total) {
		print "$key Filtered ALL Min_dataportion $data_portion <  $min_dataportion_total\n" if $debug;
		$print_all = 0;
	}
	
	#Only keep rare nonsense/missense/frameshift
	my $priority_flag = 1;
	
	if ($vep) {
		if ($aa_change eq 'NO_AA_CHANGE') {
			$priority_flag = 0 unless $var_consequence =~ /frameshift/ || $var_consequence =~ /stop_gained/;
		} 
		
		if ($gmaf =~ /\d/) {
			if ($gmaf > $rare_cutoff) {
				$priority_flag = 0;
			}
		} 
		
		if ($gnomad =~ /\d/) {
			if ($gnomad > $rare_cutoff) {
				$priority_flag = 0;
			}
		} 
		
	}
	
	
	#Here we filter for specific samples; 
	if (keys %samples) {
		my $report = 0;
		for my $var_sample ( @{$data{$key}{var_samples}} ) {
			#Here we found the variant if at least once required sample
			if (exists $samples{$var_sample}) {
				$report = 1;
			}
		}
		#Means we didn't find any sample containing the variant
		next if $report == 0;
	}
	
	#Now we check if it's somatic
	if ($only_somatic) {
		my $somatic = 1;
		for my $var_sample ( @{$data{$key}{var_samples}} ) {
			#Here we found the variant in at least once control sample so don't keep report it
			if (!exists $samples{$var_sample}) {
				$somatic = 0;
			}
		}
		#Means a sample not in the list had the variant so it's not somatic
		next unless $somatic == 1;
	} elsif (%controls) {
		my $control = 0;
		for my $var_sample ( @{$data{$key}{var_samples}} ) {
			#Here we found the variant in at least once control sample so don't keep report it
			if (exists $controls{$var_sample}) {
				$control = 1;
			}
		}
		
		#Skip as found in control
		next unless $control == 0;
		
		#Special handling for sample/control and min_sample_count
		my $sample_include_count = 0;
		
		for my $var_sample ( @{$data{$key}{var_samples}} ) {
			#Here we found the variant if at least one required sample
			if (exists $samples{$var_sample}) {
				$sample_include_count++;
			}
		}
		
		if ($min_varcount_total > $sample_include_count) {
			next;
		}
	}
	
	if ($aa_change =~ /Stop/) {
		$poly_score = 'N/A';
	}
	
	if ($sift_cat eq 'N/A') {
		$sift_score = 'N/A';
	}
	
	#Here all filters have been applied so store results to print
	print GROUP join("\t", $key,@celltype_vars,$var_samples) . "\n";
	print GROUPPRIORITY join("\t", $key,@celltype_vars,$var_samples) . "\n" if $priority_flag;
	
	my @first_alllines = ($chr,$start,$end,$data{$key}{qual},$average_score,$var_str,$ref_count,$nd_count,$var_mean_af,$var_median_af,$var_samples,$data{$key}{var_type},$data{$key}{ref},$var_base);
	my @last_alllines = ($genename,$ens_gene,$ens_trans,$rs,$gnomad,$gmaf,$var_consequence,$aa_change,$poly_cat,$poly_score,$sift_cat,$sift_score,$domain,$pubmed,$clin);
	
	#Here we divide variants into cluster and non-cluster; only proceed with non-filtered variants
	if ($cell_anno && $print_all) {
		
		#New MB variables needed;  set defaults to account to few complex overlapping cases...	
		my $cluster_var_portion = 'N/A';
		my $noncluster_var_portion = 'N/A';
		my $cluster_var_count = 'N/A';
		my $noncluster_var_count = 'N/A';
		my $diff = 'N/A';
		my $OR = 'N/A';
		my $zscore = 'N/A';
		my $portion_cluster = 'N/A';
		my $celltype_datacount = 'N/A';
		
		for my $celltype (sort keys %celltype_stats) {
			my $celltype_printall = 1;
			my $celltype_priority = 1;		
			
			$cluster_var_count = $celltype_stats{$celltype}{$key}{varcount} if $celltype_stats{$celltype}{$key}{varcount} =~ /\d/;
    		$cluster_var_portion = $celltype_stats{$celltype}{$key}{varportion} if $celltype_stats{$celltype}{$key}{varportion} =~ /\d/;
    		$noncluster_var_count = $celltype_stats{$celltype}{$key}{other_varcount} if $celltype_stats{$celltype}{$key}{other_varcount} =~ /\d/;
    		$noncluster_var_portion = $celltype_stats{$celltype}{$key}{other_varportion} if $celltype_stats{$celltype}{$key}{other_varportion} =~ /\d/;
    		$diff = $celltype_stats{$celltype}{$key}{portion_diff} if $celltype_stats{$celltype}{$key}{portion_diff}  =~ /\d/;
    		$portion_cluster = $celltype_stats{$celltype}{$key}{portion_in_cluster}  if $celltype_stats{$celltype}{$key}{portion_in_cluster} =~ /\d/;
			$celltype_datacount = $celltype_stats{$celltype}{$key}{datacount} if $celltype_stats{$celltype}{$key}{datacount} =~ /\d/;
    		$zscore = $celltype_stats{$celltype}{$key}{zscore} if $celltype_stats{$celltype}{$key}{zscore} =~ /\d/;
			$OR = $celltype_stats{$celltype}{$key}{odds_ratio} if $celltype_stats{$celltype}{$key}{odds_ratio} =~ /\S/;
			
			#Apply filters
			if ($cluster_var_count !~ /N\/A/ && $cluster_var_count < $min_varcount_celltype ) {
				print "$key FilteredCell $celltype clus_var_count $cluster_var_count < $min_varcount_celltype\n" if $debug;
				$celltype_printall = 0;
			} 	
			
			if ($cluster_var_portion !~ /N\/A/ && $cluster_var_portion < $min_varportion_celltype) {
				print "$key FilteredCell $celltype clus_var_portion $cluster_var_portion < $min_varportion_celltype\n" if $debug;
				$celltype_printall = 0;
			}
			
			if ($celltype_datacount !~ /N\/A/ && $celltype_datacount < $min_datacount_celltype ) {
				print "$key FilteredCell $celltype celltype_data_count $celltype_datacount < $min_datacount_celltype\n" if $debug;
				$celltype_printall = 0;
			}
			
			if ($zscore !~ /N\/A/ && $zscore < $zscore_cutoff) {
				print "$key FilteredCell $celltype zscore $zscore < $zscore_cutoff\n" if $debug;
				$celltype_printall = 0;
			}
			
			if ($OR !~ /N\/A/ && $OR < $odds_ratio_cutoff) {
				print "$key FilteredCell $celltype oddsratio $OR < $odds_ratio_cutoff\n" if $debug;
				$celltype_printall = 0;
			}
			
			if ($portion_cluster < $min_totalportion_celltype) {
				print "$key FilteredCell $celltype min_celltype_portion $portion_cluster < $min_totalportion_celltype\n" if $debug;
				$celltype_printall = 0;
			}
			
			
			#Apply filters
			if ($OR =~ /N\A/ && $OR !~ /SOMATIC/) {
				$celltype_priority = 0;
			}
			
			#Here every read is variant 
			if ($ref_count == 0) {
				$celltype_priority = 0;
			}


			if ($celltype_printall) {
				my @full_line = ();
				my @celltype_specific = ();
				push @celltype_specific, $OR, $zscore, $cluster_var_count, $cluster_var_portion, $noncluster_var_count, $noncluster_var_portion, $diff, $portion_cluster;
				@full_line = (@first_alllines,@celltype_specific,@last_alllines);
				push @{$celltype_files_to_write{$celltype}{'all'}},join ("\t",@full_line);
				
				#Check this passes both sets of filters
				if ($celltype_priority && $priority_flag) {
					push @{$celltype_files_to_write{$celltype}{'priority'}},join("\t",@full_line);
				}
				
			}
			
			
			
		}
	}  


	my @full_alllines = (@first_alllines,@celltype_columns,@last_alllines);
	push @{$all_files_to_write{'all'}}, join("\t",@full_alllines) if $print_all;
	
	
	if ($vep && $priority_flag) {
		push @{$all_files_to_write{'priority'}},join("\t",@full_alllines) if $print_all;
	}

	 
	
	
	
	
	my @mb_values = ();	
	push @mb_values, $key;
	for my $sample (@samples) {
		my $zyg = '?';
		if (defined $data{$key}{zyg}{$sample}) {
			$zyg = $data{$key}{zyg}{$sample};
		}
		
		if (exists $mb_map{$zyg}) {
			push @mb_values, $mb_map{$zyg};
		} else {
			push @mb_values, "?";
		}
	}
	print MB join("\t",
				@mb_values
				) ."\n"
}



#Write the files for all variants; sort the lines by numerical descending
my @all_lines = sort {my @afields = split("\t",$a); my @bfields = split("\t",$b); $bfields[$col_index_all] <=> $afields[$col_index_all]} @{$all_files_to_write{'all'}};
print OUT join ("\n",@all_lines);

if ($vep) {
	my @priority_lines = sort {my @afields = split("\t",$a); my @bfields = split("\t",$b); $bfields[$col_index_all] <=> $afields[$col_index_all]} @{$all_files_to_write{'priority'}};
  	print PRIORITY join ("\n",@priority_lines);
}

#Write the files for each celltype
if ($cell_anno) {
	
	for my $celltype ( keys %celltype_files_to_write ) {
	    my @celltype_all_lines = sort {my @afields = split("\t",$a); my @bfields = split("\t",$b); $bfields[$col_index_celltype] <=> $afields[$col_index_celltype]} @{$celltype_files_to_write{$celltype}{'all'}};
	    
	    for my $fh (@{$celltype_fhs{$celltype}}) {
		    print $fh join ("\n",@celltype_all_lines);
	    }
	    
	    if ($vep) {
	    	#Only create the file if there are entries
	    	if (exists $celltype_files_to_write{$celltype}{'priority'}) {
	    		#First the headers
	    		for my $fh (@{$celltype_fhs_priority{$celltype}}) {
					my @specific_celltype_headers = ();
					for my $col (@celltype_headers) {
						(my $local_col = $col) =~ s/CELLTYPE/$celltype/;
						push @specific_celltype_headers, $local_col;
					}
					print $fh join("\t",@specific_celltype_headers) ."\n\n";
				}
	    		#then add the sorted data
		    	my @celltype_priority_lines = sort {my @afields = split("\t",$a); my @bfields = split("\t",$b); $bfields[$col_index_celltype] <=> $afields[$col_index_celltype]} @{$celltype_files_to_write{$celltype}{'priority'}};
		    	for my $fh (@{$celltype_fhs_priority{$celltype}}) {
			    	print $fh join ("\n",@celltype_priority_lines);
	    		}
	    		
	    	}
	    }
	}
	
}




