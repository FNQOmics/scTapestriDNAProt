package Vcf;

use strict;
use Data::Dumper;
use Exception;
use Statistics::Descriptive::Full;

sub new {
    my ($class, @args) = @_;

	my @required_args = (
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    Exception->throw("Required argument [$required_arg] not set");
		}
    }


    my $self = bless {}, $class;
    
    my $organism = defined $args{-organism}?$args{-organism}:'human';
    
    $self->organism($organism);
    
    return $self;
}




#Get the organism
sub organism {
    my ($self, $organism) = @_;

    if (defined $organism) {
		$self->{'organism'} = $organism;
    } elsif (! defined $self->{'organism'}) {
		Exception->throw("organism not set");
    }

    return $self->{'organism'};
}

#Parse a vcf file
sub parse_vcf {
    my ($self, @args) = @_;

    my @required_args = (
			             -vcf_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	Exception->throw("File $args{-vcf_file} doesn't exist");	
    }
    
    
    my %vcf_data = ();
    

    open(VCF,$args{-vcf_file}) || Exception->throw("Can't open file $args{-vcf_file}\n");
    
   	#Genotype quality index
   	my $gq_index = 0;
   	my $sample_index = 0;


    while (<VCF>) {
	    next if /^#/;
    	next unless /\w/;
    	chomp;
    	my ($chr,$first_coord,undef,$ref,$var_str,$qual,undef,$rest,$gt_fields,@alleles) = split;
    	my $line = $_;
   		my @fields = split;
        my @gt_fields = split(':',$gt_fields) if $gt_fields;
		my @details = split(';',$rest);	

    	
    	if ($chr eq 'MT') {
			$chr = 'M';
		}
    	

    	if ($qual !~ /\d/ && $qual ne '.') {
    		Exception->throw("ERROR: Error with qual $qual format at line $_");
    	}

		if ($ref =~ /N/ || $var_str =~ /N/) {
			next;
		}

		if ($var_str eq '.') {
			next;
		}


		my @vars = split(",",$var_str);

				

		my @ac_fields =  ();
		my $total_ac = 0;

		for my $detail (@details) {
			if ($detail =~ /^AC=/) {
				$detail =~ s/AC=//;
				@ac_fields = split(",",$detail);
				for my $count (@ac_fields) {
					$total_ac += $count;
				} 
			}
			
			
			
		}

		my $zyg_num = 1;	
		for my $var ( @vars ) {
			next if $var eq '*'; #Due to upstream deletion
			next if $var eq $ref; #Bug with Tapestri / very rare
			my ($var_key,$var_type,$ref_base) = _get_variant_key(-type=>'vcf',-chrom=>$chr,-first=>$first_coord,-ref_seq=>$ref,-var_seq=>$var);
			next unless $var_key =~ /(\d+)-(\d+):.+/;
			my ($start,$end) = $var_key =~ /(\d+)-(\d+)/;

			$vcf_data{$var_key}{var_count} = shift @ac_fields;
			$vcf_data{$var_key}{zyg_count} = $zyg_num;
			
			if ($qual eq '.') {
				$vcf_data{$var_key}{qual} = "N/A";				
			} else {
				$vcf_data{$var_key}{qual} = $qual;
			}
			$vcf_data{$var_key}{type} = $var_type;
			$vcf_data{$var_key}{ref_base} = $ref_base;
			$vcf_data{$var_key}{zyg} = \@alleles;

			#Get the number of variant samples
			my ($allele_total) = $rest =~ /AN=(\d+)/;
			$allele_total = 0 if !defined $allele_total;
			$vcf_data{$var_key}{allele} = $total_ac . '/' . $allele_total;
			
			#Get the frequency across all the samples containing the variant
				
			my $stats = Statistics::Descriptive::Full->new();
			
			my @var_counts;
			
			
			for my $allele_str (@alleles) {
				my @gt_fields_data = split(':',$allele_str);
				#print Dumper \@gt_fields_data;
				my $var_count_str = my $sample_af = 0;
					
				#Some strange GTs don't have multiple values listed....
				if (@gt_fields_data > 1) {
    				my @ads = split(',',$gt_fields_data[1]);
    				my $sum = 0;
    				for my $element (@ads) {
    					$sum += $element if $element =~ /\d/;
					}	
						
					my $var_count;
					if (defined $ads[$zyg_num] && $ads[$zyg_num] =~ /\d/) {
						$var_count = $ads[$zyg_num];
					} else {
						$var_count = 0;	
					}
    				if ($sum == 0) {
    					$sample_af = 0;
    				} elsif ($var_count == 0) {
    					$sample_af = 0;
    				} else {
    					$sample_af = sprintf("%.4f",$ads[$zyg_num]/$sum)
    				}
    				$var_count_str =  "$var_count/$sum";
				}
    				
    			push @var_counts,$var_count_str;
    				
    			$stats->add_data($sample_af);
			}
				
			my $mean = defined $stats->mean?$stats->mean:0;
			my $median = defined $stats->median?$stats->median:0;
			
			$vcf_data{$var_key}{mean_af} = sprintf("%.3f",$mean);
			$vcf_data{$var_key}{median_af} = sprintf("%.3f",$median);
			$vcf_data{$var_key}{var_read_counts} = join(",",@var_counts);
			
			$zyg_num++;
			push @{$self->{vcf_order}}, $var_key;
		}
    }
    $self->{data}{$args{-vcf_file}} = \%vcf_data;
}

#Write normalised file
sub write_normalised {
	my ($self, @args) = @_;

    my @required_args = (
			             -vcf_out,
			             -vcf_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    
    
    
    open(FILE,">$args{-vcf_out}") || Exception->throw("Can't open file to write $args{-vcf_out}\n");
	
	
	my $vcf_data = $self->get_vcf(-vcf_file=>$args{-vcf_file});
	
	for my $var_key (@{$self->{vcf_order}}) {	
		
		my ($chr,$coord,$event) = split(':',$var_key);
		my ($start,$end) = split('-',$coord);
		
		$chr =~ s/chr//;
	
		my $var_count = 0;
		if (exists $vcf_data->{$var_key} && defined $vcf_data->{$var_key}{var_count}) {
			$var_count = $vcf_data->{$var_key}{var_count};
		}	
	
		my $zyg_count = 0;
	    if (exists $vcf_data->{$var_key} && defined $vcf_data->{$var_key}{zyg_count}) {
	    	$zyg_count = $vcf_data->{$var_key}{zyg_count};
	    }
	
	
		my $vcf_str  = $vcf_data->{$var_key}{type}.';'.$var_key.';Q='.$vcf_data->{$var_key}{qual}. ';AC='.$var_count . ';ZC='.$zyg_count.";ALLELE=".$vcf_data->{$var_key}{allele}.".MEANAF=".$vcf_data->{$var_key}{mean_af}.";MEDAF=".$vcf_data->{$var_key}{median_af}.";VAR_READ_COUNTS=".$vcf_data->{$var_key}{var_read_counts};
		

		if ($vcf_data->{$var_key}{type} eq 'INS') {
			$vcf_str .= ";REF=".$vcf_data->{$var_key}{ref_base};
		}
		
		
		my $zyg_str = join("\t",@{$vcf_data->{$var_key}{zyg}});
		print FILE join("\t", 
							$chr,
							$start,
							$end,
							$vcf_str,
							$zyg_str
							)."\n";
		
	}
	close FILE;
    
}


#Get vcf data for a file
sub get_vcf {
	my ($self, @args) = @_;

    my @required_args = (
			             -vcf_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	Exception->throw("File $self->{data}{$args{-vcf_file}} doesn't exist");	
    }
    return $self->{data}{$args{-vcf_file}};
}



#Check that a vcf file is 'complete'
sub check_vcf {
	my ($self,@args) = @_;
	my @required_args = (
			             -vcf_file,
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    Exception->throw("Required argument [$required_arg] not set");
		}
    }
	
	my $vcf = $args{-vcf_file};
	if ( !-e $vcf) {
		Exception->throw("File $vcf doesn't exist");
	}
	my $organism = $self->{organism};
	
	if ($organism ne 'human' && $organism ne 'mouse') {
		Exception->throw("ERROR: organism must be human or mouse");
	}
	

	
	

	return 1;
}


#Gets variant key from either vcf or vep; standardises naming for loading into data structure
sub _get_variant_key {
	 my @args = @_;
	
	 my %args = @args;


    my @required_args = (
    					-chrom,
    					-first,
    					-ref_seq,
    					-var_seq,
    					-type
    					);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    Exception->throw("Required argument [$required_arg] not set $args{$required_arg}");
		} 
    }
    
    my $ref = $args{-ref_seq};
    my $var = $args{-var_seq};
    my $first_coord = $args{-first};
    my $chr = $args{-chrom};
    my $type = $args{-type};
    
    my $start_coord = my $end_coord = my $bases;
    my $length_ref = length($ref);
    my $length_var = length($var);
    my $var_type;
    my $ref_base = "N/A";
    
    if ($type eq 'vcf') {
		if ($length_ref > $length_var) {
			$var_type = 'DEL';
			$start_coord = $first_coord + 1;
			$end_coord = $start_coord + $length_ref - $length_var - 1;				
			my $del_length = $length_ref - $length_var;
			#print "VCF R $ref L $del_length\n";
			
			$bases = '-'. substr($ref,1,$del_length);
		} elsif ($length_ref < $length_var) {
			#Add the ref length and var length difference to the coord 
			#$start_coord = $end_coord = $first_coord + 1;
			$var_type = 'INS';
			$start_coord = $end_coord = $first_coord;
			my $ins_length = $length_var - $length_ref;
			$bases = '+'.substr($var,1,$ins_length);
			$ref_base = substr($var,0,1);
		} elsif ($length_ref == $length_var && $length_ref != 1) {
			#Handling for cases like AT->AC ot ATATA->CTATA; turn first different base into an SNV
			$var_type = 'SNV';
			$start_coord = $end_coord = $first_coord;
			
			my @ref_bases = split("",$ref);
			my @var_bases = split("",$var);

			if ($ref_bases[0] ne $var_bases[0]) {
				#Simplest case where first base is different
				$bases = $ref_bases[0] . '->' .$var_bases[0];
			} else {
				my $until_snp_count = 0; #How far away from start we are
				my $length = @ref_bases;
				while ($length > 0) {
					if ($ref_bases[$until_snp_count] ne $var_bases[$until_snp_count]) {
						#print "Before $chr : $start_coord - $end_coord : $ref -> $var\n";
						$start_coord += $until_snp_count;
						$end_coord += $until_snp_count;
						$bases = $ref_bases[$until_snp_count] . '->' .$var_bases[$until_snp_count];
						#print "After $chr : $start_coord - $end_coord : $bases\n\n";
						last;
					}
					$until_snp_count++;
					$length--;
				}
				
			}
			
		} else {
			$var_type = 'SNV';
			$start_coord = $end_coord = $first_coord;
			$bases = $ref . '->' .$var;
		}
    	
    } else {
    	if ($ref eq '-' || $length_ref < $length_var) {
	    	$start_coord = $end_coord = $first_coord - 1;
			$var_type = 'INS';	
			my $ins_length = $length_var - $length_ref;
			$ins_length++ if $ref eq '-';
			$bases = '+'.substr($var,0,$ins_length);
		}  elsif ($var eq '-' || $length_ref > $length_var) {
			$var_type = 'DEL';
			$start_coord = $first_coord;
			my $del_length = $length_ref - $length_var;
			$end_coord = $start_coord + $length_ref - $length_var - 1;
			$del_length++ if $var eq '-';
			$end_coord++ if $var eq '-';
			$bases = '-'. substr($ref,0,$del_length);	
		} elsif ($length_ref == $length_var && $length_ref == 1) {
			#single snvs
			$var_type = 'SNV'; 
			$bases = $ref .'->'.$var;
			$start_coord = $end_coord = $first_coord;
		} else {
			Exception->warning("ERROR: Can't identify var_type doesn't match any var type\n");
			#next;
		}
    }
    
	my $var_key = $chr . ':'.$start_coord .'-'.$end_coord .':' .$bases;
	return($var_key,$var_type,$ref_base);
	
}


return 1;
