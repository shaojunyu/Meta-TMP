#!/usr/bin/perl

## author yushaojun
## last update:
## Huazhong Universitiy of Science and Technology

use utf8;
use strict;
use Pod::Usage;
use Getopt::Long;
use Data::Dumper;
use v5.10.1;
use File::Path qw(make_path remove_tree);
use File::Basename;
no if ($] >= 5.018), 'warnings' => 'experimental';

use constant VERSION => "1.0";

## Parse options
my %opt;
my $config_file;
my %runtime_config;
my %meta_data;
my $Meta_TMP_Home;

## Check for config file
die("$0: Too many files given.\n")  if (@ARGV > 1);
pod2usage("$0: Please input the config_file.\n")  if (@ARGV == 0);
$config_file = $ARGV[0];
pod2usage("$0: Invalid config_file: $config_file\n") unless(-e $config_file);

## check for ENV
my $Meta_TMP_Home = $ENV{"MetaTMP"};
pod2usage("Invalid MetaTMP ENV!") unless(-e $Meta_TMP_Home);

##init config
%runtime_config =load_runtime_config($config_file);
check_config(%runtime_config);

## Parallel-META

open (IN,"<",$runtime_config{"sequence_list_file"});
open (OUT_TAXA, ">", $runtime_config{"single_sample_list_dir"}."/taxa.list") || die($!);
open (OUT_FUNC, ">", $runtime_config{"single_sample_list_dir"}."/func.list") || die($!);
while (<IN>) {
	chomp;
	my $sample_path = $_;
	my $sample_name = basename($sample_path);
	my $output_path = $runtime_config{"single_sample_dir"}."/".$sample_name;
	print OUT_TAXA $sample_name."\t".$output_path."/classification.txt\n";
	print OUT_FUNC $sample_name."\t".$output_path."/functions.txt\n";
	# make_path($output_path);
	if ($runtime_config{"sequence_type"} eq "F") {#16s rRNA
		#system("$Meta_TMP_Home/bin/parallel-meta -r $sample_path -o $output_path -e $runtime_config{\"alignment_mode\"} -f F -L $runtime_config{\"rRNA_length_threshold\"} -k $runtime_config{\"sequence_format_check\"}" );
	}else{
		#system("$Meta_TMP_Home/bin/parallel-meta -m $sample_path -o $output_path -e $runtime_config{\"alignment_mode\"} -f F -L $runtime_config{\"rRNA_length_threshold\"} -k $runtime_config{\"sequence_format_check\"}" );
	}
}
close IN;
close OUT_FUNC;
close OUT_TAXA;

if ($runtime_config{"functional_analysis"} eq "T") {
	# system("$Meta_TMP_Home/bin/class-func -l $runtime_config{\"single_sample_list_dir\"}/taxa.list -o $runtime_config{\"single_sample_dir\"}");
}

##Sample visualization
# system("$Meta_TMP_Home/bin/class-tax -l $runtime_config{\"single_sample_list_dir\"}/taxa.list -o $runtime_config{\"single_sample_views_dir\"} ");

###Feature Selection
my @taxonomical_levels = split(/ /,$runtime_config{"taxonomical_levels"});
my @T_level = ("phylum", "class", "order", "family", "genus", "species", "OTU");
# print Dumper($runtime_config{"taxonomical_levels"});
# print Dumper(@levels);
foreach my $level (@taxonomical_levels){
	system("$Meta_TMP_Home/bin/taxa-sel -L $level -l $runtime_config{\"single_sample_list_dir\"}/taxa.list -o $runtime_config{\"abundance_dir\"}/taxa.$T_level[$level-1]  -P T");
}
## OTU
system("$Meta_TMP_Home/bin/taxa-sel -L 7 -l $runtime_config{\"single_sample_list_dir\"}/taxa.list -o $runtime_config{\"abundance_dir\"}/taxa.OTU  -P T");

if ($runtime_config{"functional_analysis"} eq "T") {
	my @functional_levels = split(/ /, $runtime_config{"functional_levels"});
	my @F_level = ("l1", "l2", "l3", "KO");
	foreach my $level (@functional_levels){
		system("$Meta_TMP_Home/bin/func-sel -L $level -l $runtime_config{\"single_sample_list_dir\"}/func.list -o $runtime_config{\"abundance_dir\"}/func.$F_level[$level-1]");
	}
	system("$Meta_TMP_Home/bin/class-func-nsti -l $runtime_config{\"single_sample_list_dir\"}/func.list -o $runtime_config{\"abundance_dir\"}/func.nsti");
}

#Dist Calculation
system("$Meta_TMP_Home/bin/comp-sam -T $runtime_config{\"abundance_dir\"}/taxa.OTU.Abd -o $runtime_config{\"distance_matrix_dir\"}/taxa.dist -d T -c 2");
if ($runtime_config{"functional_analysis"} eq "T") {
	system("$Meta_TMP_Home/bin/comp-sam -T $runtime_config{\"abundance_dir\"}/func.KO.Abd -o $runtime_config{\"distance_matrix_dir\"}/func.dist -d T -c 2");
}

#Correlation Calculation
foreach my $level (@taxonomical_levels){
	system("$Meta_TMP_Home/bin/comp-corr -i $runtime_config{\"abundance_dir\"}/taxa.$T_level[$level-1].Abd -o $runtime_config{\"network_dir\"}/taxa.$T_level[$level-1] -N F -f 1");
}






############################################################################
## sub functions
############################################################################
sub load_runtime_config {
	pod2usage("parameters error...") unless (@_ == 1);
	my %config_hash = (
		#basic para
		"sequence_list_file"=>"",
		"meta_data_file"=>"",
		"output_path"=>"",

		# Profiling parameters
		"sequence_type"=>"",
		"alignment_mode"=>"",
		"pair_end_sequence_orientation"=>"",
		"rRNA_length_threshold"=>"",
		"sequence_format_check"=>"",
		"functional_analysis"=>"",

		# Statistic parameters
		"taxonomical_levels"=>"",
		"functional_levels"=>"",
		"sequence_number_normalization_depth"=>"",
		"bootstrap_for_sequence_number_normalization"=>"",
		"rarefaction_curve"=>"",
		"paird_samples"=>"",
		"cluter_number"=>"",
		"network_analysis_edge_threshold"=>""
	);
	my $config_file = $_[0];
	open(IN,"<$config_file") || die("cant not open config file!");
	while (<IN>) {
		chomp;
		unless ($_ =~ /^\w/) {
			next;
		}
		unless ($_ =~ /([^=]+)=([^=]+)/) {
			next;
		}
		my $key = $1;
		my $value = $2;
		$key = trim($key);
		$value = trim($value);
		# print $key."\n";
		# print $value."\n";
		if (exists($config_hash{$key})) {
			$config_hash{$key} = $value;
		}
	}
	close(IN);

	$config_hash{"single_sample_dir"} = $config_hash{"output_path"}."/"."Single_Sample";
	$config_hash{"single_sample_list_dir"} = $config_hash{"output_path"}."/"."Single_Sample_List";
	$config_hash{"single_sample_views_dir"} = $config_hash{"output_path"}."/"."Sample_Views";
	$config_hash{"abundance_dir"} = $config_hash{"output_path"}."/"."Abundance_Tables";
	$config_hash{"distance_matrix_dir"} = $config_hash{"output_path"}."/"."Distance_Matrix";
	$config_hash{"network_dir"} = $config_hash{"output_path"}."/"."Network";

	return %config_hash;
}

sub trim { 
	my $s = shift; 
	$s =~ s/^\s+|\s+$//g;
	return $s 
}

sub check_config {
	my (%config_hash) = @_;
	# print Dumper(%config_hash);
	foreach my $key(keys %config_hash){
		my $value = $config_hash{$key};
		given ($key){
			########################### basic para
			when ("sequence_list_file"){
				unless (-e $value){ pod2usage("Meta-TMP Error: Invalid sequence list file"); }
				my $seq_file_count = `wc -l < $value`;
				chomp($seq_file_count);
				my $sample_count = `wc -l < $config_hash{"meta_data_file"}`;
				chomp($sample_count);
				unless ($sample_count - 1 == $seq_file_count){ pod2usage("Meta-TMP Error: Sequence files (pairs) and meta data should have the same sample number and order"); }
			}
			when ("meta_data_file"){
				unless (-e $value){ pod2usage("Meta-TMP Error: Invalid meta data file"); }
			}
			when ("output_path"){
				if (-e $value) {
					# check for writeable
					unless (-d $value && -w $value) {
						pod2usage("Meta-TMP Error: Can not write to  $value!");
					}
					# check for empty
					opendir my $dir, $value or die $!;
					if( grep ! /^\.\.?$/, readdir $dir ){
					  	# pod2usage("Meta-TMP Error: output directory is not empty!");
					}
				}else{
					make_path($value, {error => \my $err});
					if (@$err) {
						pod2usage("Meta-TMP Error: Fail to create output path!");
					}
				}
			}
			###################
			when ("sequence_type"){
				if ($value =~ /[^TF]/) {#0: very fast, 1: fast, 2: sensitive, 3: very-sensitive
					pod2usage("Invalid sequence_type parameter");
				}
			}
			when ("alignment_mode"){
				if ($value =~ /[^0123]/) {
					pod2usage("Invalid alignment mode parameter");
				}
			}
			when ("pair_end_sequence_orientation"){ #0: Fwd & Rev, 1: Fwd & Fwd, 2: Rev & Fwd, default is 0
				if ($value =~ /[^012]/) {
					pod2usage("Invalid pair-end_sequence_orientation parameter");
				}
			}
			when ("rRNA_length_threshold"){
				if ($value =~ /[^\d]/) {
					pod2usage("Invalid rRNA_length_threshold parameter");
				}
			}
			when ("sequence_format_check") {
				if ($value =~ /[^TF]/) {
					pod2usage("Invalid sequence_format_check parameter");
				}
			}
			when ("functional_analysis") {
				if ($value =~ /[^TF]/) {
					pod2usage("Invalid functional_analysis parameter");
				}
			}
			#################################Statistic parameters
			when ("taxonomical_levels"){
				if ($value =~ /[^123456\s]/) {
					pod2usage("Invalid taxonomical_levels parameter");
				}
			}
			when ("functional_levels"){
				if ($value =~ /[^1234\s]/) {
					pod2usage("Invalid taxonomical_levels parameter");
				}
			}
			when ("sequence_number_normalization_depth"){
				if ($value =~ /[^\d]/) {
					pod2usage("Invalid sequence_number_normalization_depth parameter");
				}
			}
			when ("bootstrap_for_sequence_number_normalization"){
				if ($value =~ /\d+/ and $value <= 800 and $value >= 200) {
				}else{
					pod2usage("Invalid bootstrap_for_sequence_number_normalization parameter");
				}
			}
			when ("rarefaction_curve"){
				if ($value =~ /[^TF]/) {
					pod2usage("Invalid rarefaction_curve parameter");
				}
			}
			when ("paird_samples"){
				if ($value =~ /[^TF]/) {
					pod2usage("Invalid paird_samples parameter");
				}
			}
			when ("cluter_number"){
				if ($value =~ /\d+/ and $value > 0) {
				}else{
					pod2usage("Invalid cluter_number parameter");
				}
			}
			when ("network_analysis_edge_threshold"){
				if ($value) {
				}else{
					pod2usage("Invalid network_analysis_edge_threshold parameter");
				}
			}
			when (/_dir$/){
				unless (-e $value) {
					make_path($value);
				}
			}

		}#end given $key
	}# end foreach
}#end sub check_config