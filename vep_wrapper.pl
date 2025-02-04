#! /usr/bin/perl -w

use strict;
use FindBin qw( $RealBin );
use lib $RealBin;
use Exception;
use VEP;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use Data::Dumper;

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
		"all",
		"vep_in=s",
		"vep_conf=s",
		"parse_file=s"
    	);

	   
pod2usage(1) if ($OPT{help} || !$OPT{vep_in});


=pod

=head1 SYNOPSIS

./vep_wrapper.pl -all run_on_all -indel run_for_indels(default=SNVs) -vep_in pass_in_batch_field -parse_file parse_file_name(don't run)

Required flags: -vep_in

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

vep_wrapper.pl -> Get vep info from variant outside pipeline

=head1 DESCRIPTION

Mar 30, 2011

It's ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./sample_list.pl -sample_type normal_cell_line

=cut

my $vep_conf = defined $OPT{vep_conf}?$OPT{vep_conf}:'./tapestri.conf';

my %vep_confs  = ();
my $workdir = `pwd`;
chomp $workdir;

open(CONF,"$vep_conf") || Exception->throw("Can't open file $vep_conf\n");


while (<CONF>) {
	chomp;
	my ($arg,$value) = split(",",$_);
	$vep_confs{$arg} =$value;
}

my $vep_db_dir = 0;

if (exists $vep_confs{'vep_db'}) {
	$vep_db_dir = $vep_confs{'vep_db'};
}

my $vep_executable = $vep_confs{'vep_bin'};
my $exon = defined $OPT{all}?0:1;

my $vep_args;
if ($exon) {
	$vep_args = $vep_confs{'vep_args_exon'};
} else {
	$vep_args = $vep_confs{'vep_args_all'};
}

my $vep_input = 'vep.input'. $$;
open(INPUT,">$vep_input") || Exception->throw("Can't open file for input");

my %ref_base_map;


open(FILE,"$OPT{vep_in}") || Exception->throw("Can't open file to write $OPT{vep_in}\n");
while (<FILE>) {
    chomp;
    next unless $_ =~ /^[0-9XY]/;
    my ($chr,$start_coord,$end_coord,$ref_base,$var_base) = split;
	$ref_base_map{"$chr:$start_coord-$end_coord"} = $ref_base;
	print INPUT join("\t",
				$chr,
				$start_coord,
				$end_coord,
				$ref_base .'/'.$var_base,
				'+'
				) . "\n";
	    
}

	
my $outfile;

if ($OPT{parse_file}) {
	$outfile = $OPT{parse_file};
	if ( !-e $outfile ) {
		Exception->throw("File $outfile doesn't exist");
	}
} else {
	$outfile = 'vep.out'.$$;
}


my $vep;

if (!$vep_db_dir) {
	$vep = VEP->new(-input_file      => $vep_input,
							-executable_path => $vep_executable,
							-working_dir     => $workdir,
							-exon => $exon,
							-vep_args => $vep_args,
							-output_file => $outfile
							);
} else {
	$vep = VEP->new(-input_file      => $vep_input,
							-executable_path => $vep_executable,
	  						-db_dir          => $vep_db_dir,
							-working_dir     => $workdir,
							-exon => $exon,
							-vep_args => $vep_args,
							-output_file => $outfile
							);
}
							
$vep->run unless $OPT{parse_file};
	
my $results = $vep->parse_result;

foreach my $result (@$results){
	if ($exon) {
		my ($snv_chr,
		    $snv_start,
		    $snv_end,
		    $var_type,
			$se_type,
		    $var_base,
		    $snv_gene,
		    $snv_exon,
		    $snv_ref_aa,
		    $snv_var_aa,
		    $snv_aa_pos,
		    $poly_predict,
		    $poly_score,
		    $sift_predict,
		    $sift_score,
	    	    $cadd_phred) = @$result;


		my $ref_base = $ref_base_map{"$snv_chr:$snv_start-$snv_end"};
		my $aa_string = $snv_ref_aa.$snv_aa_pos.$snv_var_aa;
		print join("\t",
					$snv_chr,
					$snv_start,
					$snv_end,
					$ref_base,
					$var_base,
					$aa_string,
					$snv_gene,
					$snv_exon,
					$poly_predict,
	    			$poly_score,
	    			$sift_predict,
	    			$sift_score,
				$cadd_phred
					) ."\n";
	    
	}  else {
		my ($snv_chr,
		    $snv_start,
		    $snv_end,
			$var_type,
			$vep_category,
		    $var_base,
		 	$rs, 
		 	$gmaf, 
		 	$domain, 
		 	$pubmed, 
		 	$clinical, 
		 	$exon_str,
		 	$snv_gene,
		 	$snv_trans,
		 	$cadd_phred,
	    	    $gnomad_af,
	    	    $gene_name) = @$result;
		 my $ref_base = $ref_base_map{"$snv_chr:$snv_start-$snv_end"};
		 print join("\t",
					$snv_chr,
					$snv_start,
					$snv_end,
					$ref_base,
					$var_base,
					$rs, 
		 			$gmaf, 
		 			$domain, 
		 			$pubmed, 
		 			$clinical, 
		 			$exon_str,
		 			$snv_gene,
		 			$snv_trans,
					$cadd_phred,
					$vep_category,
	    	    $gnomad_af,
	    	    $gene_name
					) ."\n";
	}
	
		
	
} 










