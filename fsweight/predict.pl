#!/usr/bin/perl
use strict;
use PP::PP_Engine;
use PP::PP_Validation;
use PP::PP_Weighting;
use Getopt::Std;

my %argopts;
if (! getopts('s:a:i:w:f:t:o:m:n:y:p:x:b:r:v:lkgue', \%argopts)) { printusage(); }

#get params

# program to run
if (!defined($argopts{p})) { $argopts{p} = 'p'; }
my $action = lc($argopts{p});
if ($action ne 'p' && $action ne 'c' && $action ne 'w' && $action ne 'v'  && $action ne 'd')
{
  die("Invalid action. Must be (c, p or w)\n");
}

# check if annotation scheme is supplied
if ($action ne 'w' && !defined($argopts{s})) { printusage(); } 
my $annot_scheme = $argopts{s};

# check if annotations are supplied
if ($action ne 'w' && !defined($argopts{a})) { printusage(); }
my $annot_file = $argopts{a};

# check if interactions are supplied
if (!defined($argopts{i})) { printusage(); }
my $int_file = $argopts{i};

# check if outputfile is supplied
my $output_file = 'output.txt';
if (defined($argopts{o})) { 
 $output_file = $argopts{o};
}

# number of folds, default is 1 (LOOCV)
if (!defined($argopts{f})) { $argopts{f} = 1; }
my $folds = $argopts{f};
if ($folds <1 || $folds > 10)
{
 die("Number of folds muste be between 1 to 10\n");
}

# number of bins, default is 20. Only applicable to prob_combine
if (!defined($argopts{b})) { $argopts{b} = 20; }
my $bins = $argopts{b};
if ($bins <1 || $bins > 100)
{
 die("Number of bins muste be between 1 to 100\n");
}

# vague list
my $vaguelist = '0000004,0005554,0008372,0003674,0005575,0008150';
if (defined($argopts{v})) 
{ $vaguelist = $argopts{v}; }

# annotation level for reliability estimation
my $simlevel = 4;
if (defined($argopts{r})) 
{ $simlevel = $argopts{r}; } 

# global method
my $global = defined($argopts{g})?1:0;

# weighting scheme
my $weighting = defined($argopts{w})?$argopts{w}:'FSWEIGHT';

# prediction method
my $method = defined($argopts{m})?$argopts{m}:'PROB_COMBINE';

print "# Annot Scheme: $annot_scheme\n";
print "# Annotations: $annot_file\n";

# Init engine
my $p_engine = PP_Engine->new();
$p_engine->{VAGUE_LIST} = $vaguelist;
$p_engine->{SIM_LEVEL} = $simlevel;

if (defined($argopts{k}))
{
 $p_engine->{IGNOREISOLATED} = 1;
}
$p_engine->{NUMBINS} = $bins;

# Read Annotations
if (defined($argopts{s}) && defined($argopts{a}))
{
 $p_engine->ReadAnnotations($annot_scheme, $annot_file);
}
# Read interactions
my @int_files = split(/,/, $int_file);
foreach (@int_files)
{
 $_ =~ s/^\s+//;
 $_ =~ s/\s+$//;
 print "# Interactions: $_\n";
 $p_engine->ReadDataSource($_, 'Protein Interaction');
}

# use indirect interactions fof prob_combine
if ($method == 'PROB_COMBINE')
{
 my $level2 = defined($argopts{l})?1:0;
 if ($level2)
 {
  $p_engine->{USE_INDIRECT} = 1;
#  $p_engine->ComputeIndirect();
 }
}

if (defined($argopts{t}))
{
 $p_engine->ReadTestList($argopts{t}); # read testlist
}
else
{
 foreach (keys %{$p_engine->{ENTITIES}})
 {
  $p_engine->{ENTITIES}{$_}{TEST} = 1;
 }
}

# Compute Informative
my $informative = 30;
if (defined($argopts{x})) { $informative = $argopts{x}; }
$p_engine->ComputeInformative($informative);

if ($action eq 'c')
{
 print "# ACTION: Cross Validate\n";
 print "# METHOD: $method\n";
 print "# FOLDS: $folds\n";
 print "# WEIGHT: $weighting\n\n";
 $p_engine->PrintSummary();
 $p_engine->MultiLevelCrossValidate($method, $weighting, $folds, $global, $output_file);
}
elsif ($action eq 'w')
{
 print "# ACTION: Compute Interaction Weights\n";
 print "# WEIGHT: $weighting\n\n";
 $p_engine->PrintSummary();
 $p_engine->AnalyzeWeights($weighting, $output_file);
}
elsif ($action eq 'p')
{
 print "# ACTION: Predict Novel Annotations\n";
 print "# METHOD: $method\n";
 print "# WEIGHT: $weighting\n\n";
 $p_engine->PrintSummary();
 $p_engine->PredictForUnknown($method, $weighting, $output_file);
}
elsif ($action eq 'v')
{
 print "# ACTION: Validate Weight\n";
 print "# WEIGHT: $weighting\n\n";
 $p_engine->PrintSummary();
 $p_engine->ValidateFunctionalWeight($weighting);
}
elsif ($action eq 'd')
{
 print "# ACTION: Find L2 Coverage\n";
 $p_engine->PrintSummary();
 $p_engine->DetectL2OnlySimilarity();
}


sub printusage()
{
 print "\n";
 print "usage: predict -i interactions [-s annotscheme] [-a annotations] [-o outputfile] [-p action] [-f folds] [-x informative] [-k] [-m prediction method] [-w weighting scheme] [-t protein list]\n";
 print "options:\n";
 print "\t-s\tAnnotation Scheme\n";
 print "\t-a\tAnnotations\n";
 print "\t-i\tInteractions\n";
 print "\t-o\tOutput File\n";
 print "\t-p\tAction (Default = p)\n";
 print "\t-r\tMinimum annotation level for reliability estimation (Default = 4, Use 1 for annotation schemes other than GO)\n";
 print "\t-v\tComma-delimited list of vague annotation terms (Default = 0000004,0005554,0008372,0003674,0005575,0008150)\n";
 print "\t-f\tFolds (Default = 1)\n";
 print "\t\tp = prediction of novel annotations\n";
 print "\t\tc = cross validation\n";
 print "\t\tw = weight interactions\n";
 print "\t-x\tThreshold for Informative GO Terms (Default = 30)\n";
 print "\t-k\tIgnore isolated proteins (no neighbours) (Default = Off)\n";
# print "\t-l\tUse indirect interactions\n";
# print "\t-g\tIterative\n";
 print "\t-m\tPrediction Method (MAJORITY_VOTE, CHI_SQUARE, WEIGHTED_AVG, PROB_COMBINE) (Default = PROB_COMBINE)\n";
 print "\t-w\tWeighting Scheme (CD_DIST, GEOMETRIC, FSWEIGHT) (Default = FSWEIGHT)\n";
 print "\t-o\tProtein List File (Perform prediction/validation only on selected proteins)\n";
 print "\n";
 exit;
}
