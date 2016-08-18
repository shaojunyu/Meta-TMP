package PP_Engine; 
use strict;

use PP::PP_Weighting;
use PP::PP_Scoring; 
use PP::PP_DataSource; 
use PP::PP_Validation; 
#use PP::FloatMatrix; 
#use PP::OLS; 
#use Algorithm::SVM; 
#use Algorithm::SVM::DataSet; 
#use GD::Graph::points;

########################################################################################################### # Configurations ###########################################################################################################

use constant DEBUG_MODE => 0; 
use constant YEAST_GENES_ONLY => 0; 
#srand(101);
srand ();

########################################################################################################### # Constants ###########################################################################################################

use constant ENTITIES => 1; 
use constant ANNOTATIONS => 2; 
use constant INTERACTIONS => 3; 
use constant DATASOURCE => 4; 
use constant DATATYPE => 5; 
use constant WEIGHTS => 6; 
use constant RELIABILITY => 7; 
use constant INDEX => 8; 
use constant PREDICTIONS => 9; 
use constant NEIGHBOURS => 10; 
use constant LAST_NEIGHB => 11; 
use constant ENTITIES_TABLE => 12; 
use constant NAME => 13; 
use constant BASE_SIMILARITY => 14; 
use constant FOLDLABEL => 15; 
use constant INF_ANNOTS => 17; 
use constant TEST => 18; 
use constant ANNOTATED => 19; 
use constant ANNOTATIONS_SAFE => 20; 
use constant PREDICTIONS_SAFE => 21;
use constant LEVEL => 22;
use constant RANGE => 23;
use constant DELTA => 24;
use constant MINVAL => 25;
use constant CHILDREN => 26;
use constant INF => 27;
use constant PARENTS => 28;
use constant IGNOREISOLATED => 29;
use constant VAGUE_LIST => 30;
use constant SIM_LEVEL => 31;
use constant NUMBINS => 32;
use constant USE_INDIRECT => 33;
 
# PP_Engine->new() 
sub new() 
{
 my($type) = $_[0];
 my($self) = {};

 $self->{ENTITIES} = ();
 $self->{ANNOTATIONS} = ();
 $self->{DATASOURCE} = ();
 $self->{DATATYPE} = ();
 $self->{INTERACTIONS} = ();
 $self->{WEIGHTS} = ();
 $self->{NUMBINS} = 20;
 $self->{IGNOREISOLATED} = 0;
 $self->{VAGUE_LIST} = "0000004,0005554,0008372,0003674,0005575,0008150"; 
 $self->{SIM_LEVEL} = 4;
 $self->{USE_INDIRECT} = 0;

 bless($self, $type);
 return $self; 
}

########################################################################################################### # Utilities ###########################################################################################################

sub GetInteractionKey() 
{
 my ($self, $a, $b) = @_;
 return ($a lt $b)?"$a|$b":"$b|$a"; 
}

########################################################################################################### # File parsing ###########################################################################################################

sub ReadAnnotations() 
{
 ###########################################################
 # ReadAnnotations(
 # annotation_scheme_filepath,
 # annotation_filepath
 # )
 # Reads functions from scheme and annotation file
 ###########################################################

 my($self, $annotation_scheme_filepath, $annotation_filepath) = @_;
 my $line;

 ####################################################
 # clear existing annotations
 ####################################################

 $self->{ANNOTATIONS} = ();
 foreach (keys(%{$self->{ENTITIES}}))
 {
  $self->{ENTITIES}{$_}{ANNOTATIONS} = ();
 }

 ####################################################
 # read annotation scheme
 ####################################################
 
 open(DAT, "<$annotation_scheme_filepath")
 || die ("Cannot read file $annotation_scheme_filepath");

 while($line = <DAT>)
 {
  chomp($line);
  $line =~ s/^\s+//;
  $line =~ s/\s+$//;
  if (substr($line, 0, 1) eq '#' || $line eq '')
  {
   next;
  }
  my @params = split(/\s/, $line);
  my $id = shift(@params);
  my $name = join(' ', @params);

  ####################################################
  # filter vague annotations
  ####################################################

  if ($self->{VAGUE_LIST} =~ m/($id),/ || $self->{VAGUE_LIST} =~ m/,($id)$/ || $self->{VAGUE_LIST }=~ m/^($id)$/)
  { 
# print "$id $name\n";
   next;
  } 
# print "$id $name\n";
  $self->{ANNOTATIONS}{$id}{NAME} = $name;
  my @cols = split(/\|/, $name);
  if (scalar(@cols)>1)
  { $self->{ANNOTATIONS}{$id}{LEVEL} = $cols[1]; }
  else
  { $self->{ANNOTATIONS}{$id}{LEVEL} = 1; }
  if (scalar(@cols)>2)
  {
   my @parents = split(/\,/, $cols[3]);
   foreach (@parents)
   {
    if ($self->{VAGUE_LIST} =~ m/($_),/ || $self->{VAGUE_LIST} =~ m/,($_)$/ || $self->{VAGUE_LIST} =~ m/^($_)$/) { next; }
#    print "$_\n";
    $self->{ANNOTATIONS}{$_}{CHILDREN}{$id} = 1;
    $self->{ANNOTATIONS}{$id}{PARENTS}{$_} = 1;
   }
  }
#  print "$name $cols[1]\n";
 }
 close(DAT);

 if (DEBUG_MODE)
 {
  foreach(sort {$a cmp $b} keys %{$self->{ANNOTATIONS}})
  {
   print "$_ ".$self->{ANNOTATIONS}{$_}{NAME}."\n";
  }
 }

 ####################################################
 # read annotations
 ####################################################
 
 open(DAT, "<$annotation_filepath") || die ("Cannot read file $annotation_filepath");

 my $lines = 0;
 while($line = <DAT>)
 {
  chomp($line);
  $line =~ s/^\s+//;
  $line =~ s/\s+$//;
  if (substr($line, 0, 1) eq '#' || $line eq '')
  {
   next;
  }

  $lines++;
  my $id;
  my $annotation;
  ($id, $annotation) = split(/\|/, $line);

  $id = uc($id);
  if ($id eq '')
  {
   die ($line);
  }

  if (!defined($self->{ANNOTATIONS}{$annotation}))
  { 
#   print "#$id -> $annotation not in scheme!\n";
   next;
  }
  $self->{ENTITIES}{$id}{ANNOTATIONS}{$annotation} = 1;
  if (!defined($self->{ENTITIES}{$id}{LEVEL})) { $self->{ENTITIES}{$id}{LEVEL} = 0; }
  if ($self->{ANNOTATIONS}{$annotation}{LEVEL} > $self->{ENTITIES}{$id}{LEVEL}) { $self->{ENTITIES}{$id}{LEVEL} = $self->{ANNOTATIONS}{$annotation}{LEVEL}; }
  $self->{ENTITIES}{$id}{ANNOTATED} = 1;
 }
 close(DAT);

 if (DEBUG_MODE)
 {
  foreach(sort {$a cmp $b} keys %{$self->{ENTITIES}})
  {
   my $annotations = join(',', keys %{$self->{ENTITIES}{$_}{ANNOTATIONS}});
   print "$_ $annotations\n";
  }
  foreach(sort {$a <=> $b} keys %{$self->{ANNOTATIONS}})
  {
   print "$_ $self->{ANNOTATIONS}{$_}{NAME}\n";
  }
 }

 my $annotation_count = 0;
 my $annotated_entity_count = 0;
 foreach(keys(%{$self->{ENTITIES}}))
 {
  if (scalar(keys %{$self->{ENTITIES}{$_}{ANNOTATIONS}}))
  {
   $annotated_entity_count ++;
   $annotation_count += scalar(keys %{$self->{ENTITIES}{$_}{ANNOTATIONS}});
  }
 }

 print "#$lines lines\n";
 print "#$annotated_entity_count annotated entities\n";
 print "#$annotation_count annotations\n"; 
}

sub ReadDataSource() 
{
 ###########################################################
 # ReadDataSource(
 # datasource_filepath,
 # datasource_type
 # )
 # Reads binary from datasource
 ###########################################################
 my($self)= $_[0];

 my($datasource_filepath) = $_[1];
 my($datasource_type) = $_[2];
 my $line;

 ####################################################
 # read annotation scheme
 ####################################################
 open(DAT, "<$datasource_filepath")
 || die ("Cannot read file $datasource_filepath");

 my $paircounts = 0;
 my $repeats = 0;
 while($line = <DAT>)
 {
  chomp($line);
  $line =~ s/^\s+//;
  $line =~ s/\s+$//;
  my $id_a;
  my $id_b;
  my $subtype;
  my $score;
  my $key;

  ($id_a, $id_b, $subtype, $score) = split(/\t/, $line);
  if ($subtype eq '')
  {
   $subtype = 'Undefined';
  }
  if ($score eq '')
  {
   $score = 1;
  }
  $id_a = uc($id_a);
  $id_b = uc($id_b);
  if ($id_a eq $id_b) # || (!defined($self->{ENTITIES}{$id_a}{TEST}) && !defined($self->{ENTITIES}{$id_b}{TEST})))
  {
   next;
  }
  if ($id_a eq '' || $id_b eq '')
  {
   die ($line);
  }

#  print "$id_a, $id_b, $subtype, $score\n";
  $self->{DATATYPE}{$datasource_type}{$subtype}{'occurrences'} = 1;
  $self->{ENTITIES}{$id_a}{NEIGHBOURS}{$id_b}=1;
  $self->{ENTITIES}{$id_b}{NEIGHBOURS}{$id_a}=1;
  $key = $self->GetInteractionKey($id_a, $id_b);
  $score = defined($score)?$score:1;
  $paircounts ++;
  if (defined($self->{INTERACTIONS}{$key}{DATASOURCE}{$datasource_type}))
  {
    $repeats ++;
  }
  push @{$self->{INTERACTIONS}{$key}{DATASOURCE}{$datasource_type}{$subtype}{'scores'}}, $score;
 }
 close(DAT);

 print "$paircounts entity pairs\n".($paircounts - $repeats)." unique entity pairs\n";

 if (DEBUG_MODE)
 {
  my $key;
  my $datasource;
  my $subtype;
  my $id_a;
  my $id_b;

  foreach $key(keys %{$self->{INTERACTIONS}})
  {
   ($id_a, $id_b) = split(/\|/, $key);

   foreach $datasource(keys %{$self->{INTERACTIONS}{$key}{DATASOURCE}})  # data source types
   {
     foreach $subtype(keys %{$self->{INTERACTIONS}{$key}{DATASOURCE}{$datasource}})  # data source types
     {
      print "$id_a, $id_b, $datasource, $subtype, $self->{INTERACTIONS}{$key}{DATASOURCE}{$datasource}{$subtype}{'score'}\n";
     }
   }
  }

  my $id;
  my $neighbour;
  foreach $id(keys %{$self->{ENTITIES}})
  {
   print "$id\t";
   foreach $neighbour(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})
   {
    print ",$neighbour";
   }
   print "\n";
  }
 }

}

sub ReadTestList() 
{
 ###########################################################
 # ReadTestList(
 # testlist_filepath,
 # )
 # Reads test list
 ###########################################################
 my($self)= $_[0];
 my($testlist_filepath) = $_[1];
 my $line;

 foreach (keys %{$self->{ENTITIES}})
 {
  delete $self->{ENTITIES}{$_}{TEST};
 }

 if (open (DAT, "<$testlist_filepath"))
 {
  while($line = <DAT>)
  {
   $line =~ s/\s+$//; 
 # print "$line\n";
   $self->{ENTITIES}{uc($line)}{TEST} = 1;
  }
  close DAT; 
 }
 else
 {
  print "Cannot open $testlist_filepath\n";
  foreach (keys %{$self->{ENTITIES}})
  {
   $self->{ENTITIES}{$_}{TEST} = 1;
  }
 }
}

sub HasSimilarity() 
{
 my($self)= $_[0];
 my($id_a)= $_[1];
 my($id_b)= $_[2];
 my ($annot_level) = defined($_[3])?$_[3]:1;

 my $sim = 0;
 if (!scalar(keys %{$self->{ENTITIES}{$id_a}{ANNOTATIONS}}) || !scalar(keys %{$self->{ENTITIES}{$id_b}{ANNOTATIONS}}))
 {
  return 0;
 }

 foreach (keys %{$self->{ENTITIES}{$id_a}{ANNOTATIONS}})
 {
  if ($self->{ANNOTATIONS}{$_}{LEVEL} < $annot_level) { next; }
  if (defined($self->{ENTITIES}{$id_b}{ANNOTATIONS}{$_}))
  {
   return 1;
  }
 }

 return 0; 
}

sub GetSimilarity() 
{
 my($self)= $_[0];
 my($id_a)= $_[1];
 my($id_b)= $_[2];
 my ($annot_level) = defined($_[3])?$_[3]:1;

 if (!scalar(keys %{$self->{ENTITIES}{$id_a}{ANNOTATIONS}}) || !scalar(keys %{$self->{ENTITIES}{$id_b}{ANNOTATIONS}}))
 {
  return 0;
 }

 my $count = 0;
 my $count_a = 0;
 foreach (keys %{$self->{ENTITIES}{$id_a}{ANNOTATIONS}})
 {
  if ($self->{ANNOTATIONS}{$_}{LEVEL} != $annot_level) { next; }
  $count_a++;
  if (defined($self->{ENTITIES}{$id_b}{ANNOTATIONS}{$_}))
  {
    $count ++;
  }
 }
 my $count_b = 0;
 foreach (keys %{$self->{ENTITIES}{$id_b}{ANNOTATIONS}})
 {
  if ($self->{ANNOTATIONS}{$_}{LEVEL} != $annot_level) { next; }
  $count_b++;
 }

 return $count/($count_a + $count_b - $count); 
}

sub HasSharedNeighbours() 
{
 my($self)= $_[0];
 my($id_a)= $_[1];
 my($id_b)= $_[2];

 my $neighbour;
 foreach $neighbour(keys %{$self->{ENTITIES}{$id_a}{NEIGHBOURS}})
 {
  if (defined($self->{ENTITIES}{$id_b}{NEIGHBOURS}{$neighbour}))
  {
   return 1;
  }
 }
 return 0; }

sub HasL2Similarity() 
{
 my($self)= $_[0];
 my($id_a)= $_[1];
 my($id_b)= $_[2];

 my $similarity = 0;
 if (!scalar(keys %{$self->{ENTITIES}{$id_a}{ANNOTATIONS}}) || !scalar(keys %{$self->{ENTITIES}{$id_b}{ANNOTATIONS}}))
 {
  return 0;
 }

 foreach (keys %{$self->{ENTITIES}{$id_a}{ANNOTATIONS}})
 {
  if (defined($self->{ENTITIES}{$id_b}{ANNOTATIONS}{$_}))
  {
   return 1;
  }
 }

 my $neighbour;
 foreach $neighbour(keys %{$self->{ENTITIES}{$id_a}{NEIGHBOURS}})
 {
  if ($neighbour eq $id_b)
  {
   next;
  }
  foreach (keys %{$self->{ENTITIES}{$neighbour}{ANNOTATIONS}})
  {
   if (defined($self->{ENTITIES}{$id_b}{ANNOTATIONS}{$_}))
   {
    return 1;
   }
  }
 }
 foreach $neighbour(keys %{$self->{ENTITIES}{$id_b}{NEIGHBOURS}})
 {
  if ($neighbour eq $id_a)
  {
   next;
  }
  foreach (keys %{$self->{ENTITIES}{$neighbour}{ANNOTATIONS}})
  {
   if (defined($self->{ENTITIES}{$id_a}{ANNOTATIONS}{$_}))
   {
    return 1;
   }
  }
 }

 return 0 }

sub computeFunctionalSimilarity() 
{
 my($self)= $_[0];
 my($id_a)= $_[1];
 my($id_b)= $_[2];

 my $similarity = 0;
 my $Au = scalar(keys %{$self->{ENTITIES}{$id_a}{ANNOTATIONS}});
 if ($Au == 0) {return 0;}

 my $Av = scalar(keys %{$self->{ENTITIES}{$id_b}{ANNOTATIONS}});
 if ($Av == 0) {return 0;}

 foreach (keys %{$self->{ENTITIES}{$id_a}{ANNOTATIONS}})
 {
  if (defined($self->{ENTITIES}{$id_b}{ANNOTATIONS}{$_}))
  {
   $similarity ++;
  }
 }
 if ($similarity > 0)
 {
  return $similarity / ($Au + $Av);
 }
 else
 {
  return 0;
 } 
}

sub computeWeights() 
{
 ###########################################################
 # computeWeights(
 # weighting_technique
 # )
 # compute functional weight between protein pairs
 ###########################################################

 my($self) = @_;
 my $level = $self->{SIM_LEVEL};
 my $process_data = 0;
 if (defined($_[1])) { $process_data = $_[1]; }
 my $annotation;
 my $key;
 my $id_a;
 my $id_b;
 my %completed = ();

 # clear weights;
 $self->{WEIGHTS} = ();

 foreach $key(keys %{$self->{INTERACTIONS}})
 {
  delete($self->{INTERACTIONS}{$key}{RELIABILITY});
 }

 # compute background frequency
 %completed = ();
 my $sim = 0;
 my $annotated_count = 1;
# my $paircount = 0;
# my %sharedandannot = ();
 foreach $annotation(keys %{$self->{ANNOTATIONS}})
 {
  $self->{ANNOTATIONS}{$annotation}{'frequency'} = 0;
 }
 foreach $id_a(keys %{$self->{ENTITIES}})
 {
  if (!scalar(keys %{$self->{ENTITIES}{$id_a}{ANNOTATIONS}})) # || !defined($self->{ENTITIES}{$id_a}{TEST}))
  {
   next;
  }

  $annotated_count ++;
  foreach $annotation(keys %{$self->{ENTITIES}{$id_a}{ANNOTATIONS}})
  {
   $self->{ANNOTATIONS}{$annotation}{'frequency'} ++;
  }
 }

 foreach $annotation(keys %{$self->{ANNOTATIONS}})
 {
  $self->{ANNOTATIONS}{$annotation}{'frequency'} /= $annotated_count;
 }

 #compute Average neighbours
 $self->{'g_lambda'} = 0;
 foreach $id_a(keys %{$self->{ENTITIES}})
 {
  $self->{'g_lambda'} += scalar(keys %{$self->{ENTITIES}{$id_a}{NEIGHBOURS}});
 }
 $self->{'g_lambda'} /= scalar(keys %{$self->{ENTITIES}});

 # compute datasources reliability
 my $datasource_manager = PP_DataSource->new($self);

 my $datatype;
 my $subtype;
 my $range;
 my $annotation;

 if ($process_data) {
  foreach $datatype(keys %{$self->{DATATYPE}})
  {
   $datasource_manager->ProcessDataSource($datatype, 0);
  }
 }

 foreach $datatype(keys %{$self->{DATATYPE}})
 {
  foreach $subtype(keys %{$self->{DATATYPE}{$datatype}})
  {
   $self->{DATATYPE}{$datatype}{$subtype}{'occurrence'} = 1;
   if (scalar(keys %{$self->{ANNOTATIONS}})) {
   $self->{DATATYPE}{$datatype}{$subtype}{RELIABILITY} = 0.5; 
   }else{
   $self->{DATATYPE}{$datatype}{$subtype}{RELIABILITY} = 1; 
   }
#  foreach $annotation(keys %{$self->{DATATYPE}{$datatype}{$subtype}{'functions'}}) 
#  { 
#   foreach $range(sort {$a <=> $b} keys %{$self->{DATATYPE}{$datatype}{$subtype}{'functions'}{$annotation}}) 
#   { 
#    print "$datatype -> $subtype -> $annotation ($self->{DATATYPE}{$datatype}{$subtype}{'functions'}{$annotation}{$range}{'bound'}): $self->{DATATYPE}{$datatype}{$subtype}{'functions'}{$annotation}{$range}{RELIABILITY}\n"; 
#   } 
#  }
  }
 }

 my $interactions = 0;
 my $known = 0;
 my $similar = 0;
 my $key;
#  print "l$level\n";
 foreach $key(keys %{$self->{INTERACTIONS}})
 {
  my $id_a;
  my $id_b;
  ($id_a, $id_b) = split(/\|/, $key);
   
  if (!scalar(keys %{$self->{ENTITIES}{$id_a}{ANNOTATIONS}}) || !scalar(keys %{$self->{ENTITIES}{$id_b}{ANNOTATIONS}})) 
  {
   next; 
  }
  if ($self->{ENTITIES}{$id_a}{LEVEL} < $level || $self->{ENTITIES}{$id_b}{LEVEL} < $level) { next; }
  $interactions++;

  my $sim = $self->HasSimilarity($id_a, $id_b, $level);
  $known ++;
  $similar += $sim;
   
  foreach $datatype(keys %{$self->{INTERACTIONS}{$key}{DATASOURCE}})
  {
   foreach $subtype(keys %{$self->{INTERACTIONS}{$key}{DATASOURCE}{$datatype}})
   { 
# print "$key -> $datatype -> $subtype\n";
    $self->{DATATYPE}{$datatype}{$subtype}{'occurrence'} ++;
    $self->{DATATYPE}{$datatype}{$subtype}{RELIABILITY} += $sim;
   }
  }
 }

# print "$interactions interactions\n";
# print "$interactions1 interactions with known proteins\n";
# exit;
 if ($similar>0)
 {
  $self->{'g_reliability'} = $similar/$known;
 }
 else
 {
  $self->{'g_reliability'} = 0;
 }
 $self->{'g_lambda'} *= $self->{'g_reliability'};
 
 foreach $datatype(keys %{$self->{DATATYPE}})
 {
  foreach $subtype(keys %{$self->{DATATYPE}{$datatype}})
  {
   $self->{DATATYPE}{$datatype}{$subtype}{RELIABILITY} /=
   $self->{DATATYPE}{$datatype}{$subtype}{'occurrence'};
   print "$datatype -> $subtype: ".$self->{DATATYPE}{$datatype}{$subtype}{RELIABILITY}."\n";
  }
 }


 # compute links reliability
 foreach $key(keys %{$self->{INTERACTIONS}})
 {
   my $weight = 1;
   foreach $datatype(keys %{$self->{INTERACTIONS}{$key}{DATASOURCE}})
   {
    foreach $subtype(keys %{$self->{INTERACTIONS}{$key}{DATASOURCE}{$datatype}})
    {
     my $score;
     foreach $score(@{$self->{INTERACTIONS}{$key}{DATASOURCE}{$datatype}{$subtype}{'scores'}})   
     {
      $weight *= (1 - $self->{DATATYPE}{$datatype}{$subtype}{RELIABILITY});
     }
    }
   }
   $self->{INTERACTIONS}{$key}{RELIABILITY} = 1-$weight;
 }

 print "#OK\n";

 if (DEBUG_MODE)
 {
  foreach $annotation(keys %{$self->{ANNOTATIONS}})
  {
   print "$annotation\t$self->{ANNOTATIONS}{$annotation}{'frequency'}\n";
  }

  foreach $datatype(keys %{$self->{DATATYPE}})
  {
   foreach $subtype(keys %{$self->{DATATYPE}{$datatype}})
   {
    foreach $range(keys %{$self->{DATATYPE}{$datatype}{$subtype}{'ranges'}})
    {
     print "$datatype -> $subtype : $self->{DATATYPE}{$datatype}{$subtype}{'ranges'}{$range}{RELIABILITY}\n";
    }
   }
  }

  my $key;
  foreach $key(keys %{$self->{INTERACTIONS}})
  {
   ($id_a, $id_b) = split(/\|/, $key);
   print "$id_a, $id_b, $self->{WEIGHTS}{$id_a}{$id_b}\n";
  }

 }

}

sub RemoveIsolated() 
{
 my($self)= $_[0];

 my $id;
 foreach $id(keys %{$self->{ENTITIES}})
 {
  my $neighbour;
  my $num_annotated = 0;
  foreach $neighbour(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})
  {
   if (scalar(keys %{$self->{ENTITIES}{$neighbour}{ANNOTATIONS}})>0)
   {
    $num_annotated ++;
   }
  }
  if (!$num_annotated || ($num_annotated == 1 && !scalar(keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})))
  {
   #print "a\n";
   # remove from neighbours' neighbourlist and interactionlist
   foreach $neighbour(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})
   {
    delete $self->{ENTITIES}{$neighbour}{NEIGHBOURS}{$id};
    my $key = $self->GetInteractionKey($id, $neighbour);
    delete $self->{INTERACTIONS}{$key};
   }
   # remove self
   delete $self->{ENTITIES}{$id};
  }
 } 
}

sub RemoveHubs() 
{
 my($self)= $_[0];
 my($limit)= $_[1];

 my $id;
 foreach $id(keys %{$self->{ENTITIES}})
 {
  my $neighbour;
  my $num_annotated = scalar(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}});
  if ($num_annotated > $limit)
# || ($num_annotated == 1 && !scalar(keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})))
  {
   #print "a\n";
   # remove from neighbours' neighbourlist and interactionlist
   foreach $neighbour(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})
   {
    delete $self->{ENTITIES}{$neighbour}{NEIGHBOURS}{$id};
    my $key = $self->GetInteractionKey($id, $neighbour);
    delete $self->{INTERACTIONS}{$key};
   }
   # remove self
   delete $self->{ENTITIES}{$id};
  }
 } 
}

sub SetGOAnnotationLevel()
{
 my($self)= $_[0];
 my($level)= $_[1];
 my $annotation;
 my $id;
 
 foreach $id(keys %{$self->{ENTITIES}})
 {
  foreach $annotation(keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})
  {
   if ($self->{ANNOTATIONS}{$annotation}{LEVEL} != $level)
   {
    delete $self->{ENTITIES}{$id}{ANNOTATIONS}{$annotation};
   }
  }
 }
}

sub SetAnnotationLevel() 
{
 my($self)= $_[0];

 my $selectedlevel = $_[1];
 my $id;
 my $annotation;

 my @annotations;
 foreach $annotation(keys %{$self->{ANNOTATIONS}})
 {
  my $numlevel = (length($annotation)-2)/3;
  push(@{$annotations[$numlevel]}, $annotation);
 }

 my $level = scalar(@annotations)-1;
 while($level>=0)
 {
  if ($level != $selectedlevel)
  {
   foreach $annotation(@{$annotations[$level]})
   {
    if (defined($self->{ANNOTATIONS}{$annotation}))
    {
     delete($self->{ANNOTATIONS}{$annotation});
    }
   }
  }else{
   $self->{ANNOTATIONS}{$annotation}{'frequency'} = 0;
  }
  $level --;
 }

 foreach $id(keys %{$self->{ENTITIES}})
 {
  foreach $annotation(keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})
  {
   delete $self->{ENTITIES}{$id}{ANNOTATIONS}{$annotation};
   my $numlevel = (length($annotation)-2)/3;
   while($numlevel>=0)
   {
    my $annot = substr($annotation, 0, 2+($numlevel*3));
    if (defined($self->{ANNOTATIONS}{$annot}))
    {
     $self->{ENTITIES}{$id}{ANNOTATIONS}{$annot} = 1;
     $self->{ANNOTATIONS}{$annot}{'frequency'} ++;
     last;
    }
    $numlevel--;
   }
  }
 }

 foreach $annotation(keys %{$self->{ANNOTATIONS}})
 {
  if (!$self->{ANNOTATIONS}{$annotation}{'frequency'})
  {
   delete $self->{ANNOTATIONS}{$annotation};
  }
 }

# foreach (sort {$a cmp $b} keys %{$self->{ANNOTATIONS}}) 
# { 
# print "$_\t$self->{ANNOTATIONS}{$_}{NAME}\n"; 
# } 
}

sub FindInformative()
{
 my ($self, $annotation, $min) = @_;
 my $ok = 1;
# print "$annotation\n";
 foreach (keys %{$self->{ANNOTATIONS}{$annotation}{CHILDREN}})
 {
  if ($self->{ANNOTATIONS}{$_}{'occ'} >= $min)
  {
   $self->FindInformative($_, $min); 
   $ok = 0;
  }
 }

 if ($ok && $self->{ANNOTATIONS}{$annotation}{'occ'} >= $min)
 {
  $self->{ANNOTATIONS}{$annotation}{INF} = 1;
 }

}

sub ComputeInformative() 
{
 my($self, $min)= @_;
 
 my $id;
 my $annotation;
 my @level1s;
 foreach $annotation(keys %{$self->{ANNOTATIONS}})
 { 
  $self->{ANNOTATIONS}{$annotation}{'occ'} = 0;
  delete $self->{ANNOTATIONS}{$annotation}{INF};
  if ($self->{ANNOTATIONS}{$annotation}{LEVEL}==1)
  {
   push @level1s, $annotation;
  }
 }
 foreach $id(keys %{$self->{ENTITIES}})
 {
  if (
#($self->{IGNOREISOLATED} && !scalar(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})) || 
!scalar(keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}}) || !defined($self->{ENTITIES}{$id}{TEST}))
  {
   next;
  }
  foreach $annotation(keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})
  {
   $self->{ANNOTATIONS}{$annotation}{'occ'} ++;
  }
 }
 foreach (@level1s)
 {
  $self->FindInformative($_, $min);
 }
}

sub PrintSummary() 
{
 my($self)= $_[0];
 my $annotated = 0;
 my $connected = 0;
 my $id;
 my $key;

 my $informative = 0;
 foreach (keys %{$self->{ANNOTATIONS}})
 { 
# print "$self->{ANNOTATIONS}{$_}{NAME}\n";
  if (defined($self->{ANNOTATIONS}{$_}{INF}))
  {
   $informative ++;
  }
 }
 print "Annotations: ".scalar(keys %{$self->{ANNOTATIONS}})."\n";
 print "Informative Annotations: $informative\n";
 print "Entities: ".scalar(keys %{$self->{ENTITIES}})."\n";

 $informative = 0;
 foreach $id(keys %{$self->{ENTITIES}})
 {
  if (scalar(keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}}))
  {
   $annotated ++;
   foreach (keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})
   {
    if (defined($self->{ANNOTATIONS}{$_}{INF}))
    {
     $informative ++;
     last;
    }
   }
  }
  if (scalar(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}}))
  {
   $connected ++;
  }
 }

 print "Annotated Entities: $annotated\n";
 print "Inf Annotated Entities: $informative\n";
 print "Data Sources: ".scalar(keys %{$self->{DATATYPE}})."\n";
 print "Interactions Pairs: ".scalar(keys %{$self->{INTERACTIONS}})."\n";
 print "Interacting Entities: $connected\n"; 
}

sub AnalyzeWeights() 
{
 ###########################################################
 # AnalyzeWeights(
 # weighting_technique
 # level
 # )
 # performs cross validation for annotation prediction
 ###########################################################

 my($self)= $_[0];
 my($weighting_technique)= $_[1];
 my($outfile)= $_[2];

 open (DAT, ">$outfile") || die ("Cannot write to $outfile\n");
 my $weight_manager = PP_Weighting->new($self);

 $self->computeWeights();

 my $id;
 my %completed = ();
 foreach $id (keys %{$self->{ENTITIES}})
 {
  my %neighbours;
  my $neighbour;
  foreach $neighbour(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})
  {
   # l1 neighbours
   if (!defined($completed{$neighbour}))
   {
    if (!defined($neighbours{$neighbour}))
    {
     $neighbours{$neighbour}{'weight'} = $weight_manager->ComputeLocalWeight($id, $neighbour, $weighting_technique);
    }
    $neighbours{$neighbour}{'l1'} = 1;
   }
   next;

   my $l2neighbour;
   foreach $l2neighbour(keys %{$self->{ENTITIES}{$neighbour}{NEIGHBOURS}})
   {
    # l2 neighbours
    if ($l2neighbour ne $id && !defined($completed{$l2neighbour}))
    {
     if (!defined($neighbours{$l2neighbour}))
     {
      $neighbours{$l2neighbour}{'weight'} = $weight_manager->ComputeLocalWeight($id, $l2neighbour, $weighting_technique);
     }
     $neighbours{$l2neighbour}{'l2'} = 1;
    }
   }
  }

  foreach $neighbour(keys %neighbours)
  {
   if ($neighbours{$neighbour}{'l1'} )
   {
    if ($neighbours{$neighbour}{'l2'})
    {
     print DAT "$id\t$neighbour\tL12\t$neighbours{$neighbour}{'weight'}\n";
    }
    else
    {
     print DAT "$id\t$neighbour\tL1\t$neighbours{$neighbour}{'weight'}\n";
    }
   }
   elsif ($neighbours{$neighbour}{'l2'} && $neighbours{$neighbour}{'weight'}) # >= 0.05)
   {
    print DAT "$id\t$neighbour\tL2\t$neighbours{$neighbour}{'weight'}\n";
   }
  }
  $completed{$id} = 1;
 } 
 close DAT;
}

sub ComputeIndirect() 
{
 ###########################################################
 # ComputeIndirect(
 # weighting_technique
 # )
 # compute Indirect Interactions
 ###########################################################

 my($self)= $_[0];
 my($weighting_technique)= $_[1];
 my $datasource_type = 'Protein Interaction';
 my $subtype = 'INDIRECT';
 my $weight_manager = PP_Weighting->new($self);

 my $paircounts = 0;
 my $id;
 my %completed = ();
 foreach $id (keys %{$self->{ENTITIES}})
 {
  my %neighbours;
  my $neighbour;
  foreach $neighbour(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})
  {
   my $l2neighbour;
   foreach $l2neighbour(keys %{$self->{ENTITIES}{$neighbour}{NEIGHBOURS}})
   {
    # l2 neighbours
    if ($l2neighbour ne $id && !defined($completed{$l2neighbour}))
    {
     if (!defined($neighbours{$l2neighbour}) && !defined($self->{ENTITIES}{$id}{NEIGHBOURS}{$l2neighbour}))
     {
      #my $score = $weight_manager->ComputeLocalWeight($id, $l2neighbour, $weighting_technique);
      my $score = $weight_manager->ComputeFSWeight($id, $l2neighbour, 1);
      if ($score < 0.05) { next; }
      $self->{DATATYPE}{$datasource_type}{$subtype}{'occurrences'} = 1;
      $self->{ENTITIES}{$id}{NEIGHBOURS}{$l2neighbour}=1;
      $self->{ENTITIES}{$l2neighbour}{NEIGHBOURS}{$id}=1;
      my $key = $self->GetInteractionKey($id, $l2neighbour);
      push @{$self->{INTERACTIONS}{$key}{DATASOURCE}{$datasource_type}{$subtype}{'scores'}}, $score;
      $paircounts ++;
     }
    }
   }
  }
  $completed{$id} = 1;
 } 
 print "$paircounts entity pairs\n";
}

sub SVMCrossValidate() 
{
 my($self)= $_[0];
 my ($folds) = $_[1];
 my $id;

 my $i = 1;
 if ($folds > 1)
 {
  foreach $id(keys %{$self->{ENTITIES}})
  {
   $self->{ENTITIES}{$id}{'index'} = $i++;
   if (scalar(keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}}) && defined($self->{ENTITIES}{$id}{TEST}))
   {
    $self->{ENTITIES}{$id}{FOLDLABEL} = int(rand($folds));
   }
   else
   {
    $self->{ENTITIES}{$id}{FOLDLABEL} = -1;
   }
  }
 }

 my $p_validation = PP_Validation->new($self);
 my $annotation;
 for ($i=0; $i<$folds; $i++)
 {
  my $ds;
  foreach $annotation(keys %{$self->{ANNOTATIONS}})
  {
   print "$annotation\n";

   open DAT, ">train.dat";
   # train SVM
   foreach $id(keys %{$self->{ENTITIES}})
   {
    if ($self->{ENTITIES}{$id}{FOLDLABEL} != $i && $self->{ENTITIES}{$id}{FOLDLABEL} != -1)
    {
     my $label = defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$annotation})?1:0;
     print DAT "$label";
     foreach (keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})
     {
      print DAT " $self->{ENTITIES}{$_}{'index'}:1";
     }
     print DAT "\n";
    }
   }
   close DAT;
   system ("/home/kenny/libsvm-2.82/svm-train -b 1 -s 0 -t 2 train.dat svm/$annotation.model");

   # test SVM
   open DAT, ">test.dat";
   foreach $id(keys %{$self->{ENTITIES}})
   {
    if ($self->{ENTITIES}{$id}{FOLDLABEL} == $i && $self->{ENTITIES}{$id}{FOLDLABEL} != -1)
    {
     my $label = defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$annotation})?1:0;
     print DAT "$label";
     foreach (keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})
     {
      print DAT " $self->{ENTITIES}{$_}{'index'}:1";
     }
     print DAT "\n";
    }
   }
   close DAT;
   system ("/home/kenny/libsvm-2.82/svm-predict -b 1 test.dat svm/$annotation.model out.dat");

   open DAT, "<out.dat";
   my $line = <DAT>;
   foreach $id(keys %{$self->{ENTITIES}})
   {
    if ($self->{ENTITIES}{$id}{FOLDLABEL} == $i && $self->{ENTITIES}{$id}{FOLDLABEL} != -1)
    {
     $line = <DAT>;
     $line =~ s/\s+$//;
     my @cols = split(/\s+/, $line);
     $self->{ENTITIES}{$id}{PREDICTIONS}{$annotation} = $cols[2];
    }
   }
   close DAT;
  }

  $p_validation->Compute_ROC($i);
 }
 $p_validation->Compute_PrecisionVSRecall();
 $p_validation->Compute_ROC(-1); 
}

sub ReportStats()
{
 my($self)= $_[0];
 my($message)= $_[1];

 my $annotated = 0;
 my $predicted = 0;
 my $predictions = 0;
 foreach (keys %{$self->{ENTITIES}})
 {
  if (scalar(keys %{$self->{ENTITIES}{$_}{ANNOTATIONS}})) { $annotated++; }
  if (scalar(keys %{$self->{ENTITIES}{$_}{PREDICTIONS}})) { $predicted++; }
  $predictions += scalar(keys %{$self->{ENTITIES}{$_}{PREDICTIONS}});
 }
 print "###########################################################################\n";
 print "## >>> [$message] <<<\n";
 print "## >>> $annotated annotated, $predicted predicted, $predictions predictions\n";
 print "###########################################################################\n";
}

sub SplitData
{
 my($self)= $_[0];
 my $folds = $_[1];

 # Label all instances
 my $id;
 foreach $id(keys %{$self->{ENTITIES}})
 {
  $self->{ENTITIES}{$id}{FOLDLABEL} = -1;
 }

 if ($folds > 1)
 {
  foreach $id(keys %{$self->{ENTITIES}})
  {
   if (defined($self->{ENTITIES}{$id}{ANNOTATED}) && defined($self->{ENTITIES}{$id}{TEST}))
   {
    $self->{ENTITIES}{$id}{FOLDLABEL} = int(rand($folds));
   }
  }
 }

}

sub MultiLevelCrossValidate() 
{
 ###########################################################
 # MultiLevelCrossValidate(
 # prediction_technique
 # weighting_technique
 # )
 # performs cross validation for annotation prediction
 ###########################################################

 my ($self)= $_[0];
 my ($prediction_technique)= $_[1];
 my ($weighting_technique)= $_[2];
 my ($folds) = $_[3];
 my ($multi) = $_[4];
 my ($outputfile) = $_[5];

 $self->SplitData($folds);

 my $levels = 1;
 if ($multi) { $levels = 10; }

 my $score_manager = PP_Scoring->new($self);
 my $p_validation = PP_Validation->new($self);
 my $level = 0;
 my $id;

# open DAT, ">$outputfile";
 $self->ReportStats("Validation Start");
 my $i;
 for ($i=0; $i<$folds; $i++)
 {
#  print DAT "### FOLD $i ###\n";
  
  # clear weights
  foreach (keys %{$self->{WEIGHTS}} )
  {
   delete $self->{WEIGHTS}{$_};
  }
  # hide annotations for all instances in current fold
  if ($folds > 1)
  {
   foreach $id(keys %{$self->{ENTITIES}})
   {
    if ($self->{ENTITIES}{$id}{FOLDLABEL} == $i && defined($self->{ENTITIES}{$id}{ANNOTATIONS}))
    {
     %{$self->{ENTITIES}{$id}{ANNOTATIONS_SAFE}} = %{$self->{ENTITIES}{$id}{ANNOTATIONS}};
     %{$self->{ENTITIES}{$id}{ANNOTATIONS}} = ();
    }
   }

   # clear all predictions
   if (!$multi)
   {
    foreach $id(keys %{$self->{ENTITIES}})
    {
     %{$self->{ENTITIES}{$id}{PREDICTIONS}} = ();
    }
   }
  }

  if ($folds > 1)
  {
   $self->ReportStats("Fold $i start");
  }

  # compute prior probabilities
  $self->computeWeights($prediction_technique eq 'PROB_COMBINE'?1:0);
 
  for ($level = 0; $level <= $levels; $level++)
  {

   # clear all predictions
   if ($multi)
   {
    foreach $id(keys %{$self->{ENTITIES}})
    {
     %{$self->{ENTITIES}{$id}{PREDICTIONS}} = ();
    }
    $self->ReportStats("Iteration $level start");
   }

   foreach $id(keys %{$self->{ENTITIES}})
   {

    # reject unnecessary instances if not iterative
    if (!$multi && (
   !defined($self->{ENTITIES}{$id}{ANNOTATED}) ||            # not annotated
   !defined($self->{ENTITIES}{$id}{TEST}) ||                 # not test
   ($folds > 1 && $self->{ENTITIES}{$id}{FOLDLABEL} != $i)   # not in current fold
    ))
    {
     next;
    }
    my %tempannots;
    if ($multi || $folds == 1)
    {
     %tempannots = %{$self->{ENTITIES}{$id}{ANNOTATIONS}};
     %{$self->{ENTITIES}{$id}{ANNOTATIONS}} = ();
    }

    if ($prediction_technique eq 'MAJORITY_VOTE')
    {
     $score_manager->DoMajorityVote($id);
    }
    elsif ($prediction_technique eq 'PROB_COMBINE')
    {
     $score_manager->DoProbCombine($id);
    }
    elsif ($prediction_technique eq 'PROB_COMBINE_SIMPLE')
    {
     $score_manager->DoProbCombineSimple($id);
    }
    elsif ($prediction_technique eq 'PROB_COMBINE_NOWEIGHT')
    {
     $score_manager->DoProbCombineNoWeight($id);
    }
    elsif ($prediction_technique eq 'PROB_COMBINE_BAYES')
    {
     $score_manager->DoProbCombineBayes($id);
    }
    elsif ($prediction_technique eq 'WEIGHTED_MAJORITY_VOTE')
    {
     $score_manager->DoWeightedMajorityVote($id, $weighting_technique);
    }
    elsif ($prediction_technique eq 'WEIGHTED_MAJORITY_VOTE_L2')
    {
     $score_manager->DoWeightedMajorityVoteL2($id, $weighting_technique);
    }
    elsif ($prediction_technique eq 'CHI_SQUARE')
    {
     $score_manager->DoChiSquare($id);
    }
    elsif ($prediction_technique eq 'WEIGHTED_AVG')
    {
     $score_manager->DoWeightedAvg($id, $weighting_technique);
    }
    elsif ($prediction_technique eq 'WEIGHTED_AVG_L1')
    {
     $score_manager->DoWeightedAvgL1($id, $weighting_technique);
    }
    else
    {
     die ('Invalid Prediction Technique');
    }
    # print prediction
#    foreach (keys %{$self->{ENTITIES}{$id}{PREDICTIONS}})
#    {
#     print DAT "$id\t$_\t$self->{ENTITIES}{$id}{PREDICTIONS}{$_}\n";
#    }
    
    if ($multi || $folds == 1)
    {
     %{$self->{ENTITIES}{$id}{ANNOTATIONS}} = %tempannots;
    }
   }

   # assign new functions
   if (!$multi || !$self->AssignNewFunctions($level, $folds>1?$i:-1))
   {
    last;
   }

   if ($multi)
   {
    $self->ReportStats("Iteration $level ended");
   }

   if ($multi && $folds == 1)
   {
#    $p_validation->Compute_ROC(-1); 
#    $p_validation->Compute_PrecisionVSRecall();
   }

  } # end level

  # restore annotations for instances in current folds
  if ($folds > 1)
  {

   $self->ReportStats("Fold $i ended");

   foreach $id(keys %{$self->{ENTITIES}})
   {
    if ($self->{ENTITIES}{$id}{FOLDLABEL} == $i && defined($self->{ENTITIES}{$id}{ANNOTATIONS_SAFE}))
    {
     %{$self->{ENTITIES}{$id}{ANNOTATIONS}} = %{$self->{ENTITIES}{$id}{ANNOTATIONS_SAFE}};
     %{$self->{ENTITIES}{$id}{ANNOTATIONS_SAFE}} = ();
    }
   }

   # remove any temporary annotations
   if ($multi)
   {
    foreach $id(keys %{$self->{ENTITIES}})
    {
     my $annotation;
     foreach $annotation (keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})
     {
      if ($self->{ENTITIES}{$id}{ANNOTATIONS}{$annotation} > 1)
      {
       delete $self->{ENTITIES}{$id}{ANNOTATIONS}{$annotation};
      }
     }
    }
   }
   
#   $p_validation->Compute_ROC($i);

   # keep predictions for all instances in current fold
   foreach $id(keys %{$self->{ENTITIES}})
   {
    if ($self->{ENTITIES}{$id}{FOLDLABEL} == $i)
    {
     %{$self->{ENTITIES}{$id}{PREDICTIONS_SAFE}} = %{$self->{ENTITIES}{$id}{PREDICTIONS}};
    }
   }
  }

 } # end fold

 # restore predictions for instances in all fold and remove predictions in non test instances
 if ($folds > 1)
 {
  foreach $id(keys %{$self->{ENTITIES}})
  {
   if ($self->{ENTITIES}{$id}{FOLDLABEL} != -1)
   {
    %{$self->{ENTITIES}{$id}{PREDICTIONS}} = %{$self->{ENTITIES}{$id}{PREDICTIONS_SAFE}};
    %{$self->{ENTITIES}{$id}{PREDICTIONS_SAFE}} = ();
   }
   else
   {
    %{$self->{ENTITIES}{$id}{PREDICTIONS}} = ();
   }
  }
 }

# if (!$multi && $folds == 1)
 {
#  $p_validation->Compute_MultiLevelPrecisionVSRecall($outputfile);
   $p_validation->Compute_PrecisionVSRecall();
   $p_validation->Compute_ROC(-1);
 }
 
# close DAT;
}

sub AssignNewFunctions() 
{
 my ($self)= $_[0];
 my ($level)= $_[1]+2;
 my ($fold) = $_[2];

 print "$level/$fold\n";

 my $span = 10;
 my $edge = $span*2 + 1;
 my %corrects;
 my %preds;
 my %threshold;
 my %index;

 my $id;
 my %completed;
 my %predictions;
 foreach $id(keys %{$self->{ENTITIES}})
 {
  foreach (keys %{$self->{ENTITIES}{$id}{PREDICTIONS}})
  {
   if ($self->{ENTITIES}{$id}{PREDICTIONS}{$_}>0)
   {
    $predictions{"$id|$_"} = $self->{ENTITIES}{$id}{PREDICTIONS}{$_};
    $completed{$_} = 1;
   }
  }
 }
 my $predfunctions = scalar(keys %completed);

 foreach (sort {$predictions{$b} <=> $predictions{$a}} keys %predictions)
 {
  my @cols = split(/\|/, $_);
  my $id = $cols[0];
  my $function = $cols[1];

  if (defined($self->{ENTITIES}{$id}{ANNOTATED}) && ($fold == -1 || $self->{ENTITIES}{$id}{FOLDLABEL} != $fold) && !defined($threshold{$function}))
  {
   $preds{$function}[$index{$function}][0] = $predictions{$_};
   $preds{$function}[$index{$function}][1] = (defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$function}) && $self->{ENTITIES}{$id}{ANNOTATIONS}{$function} == 1)?1:0;
   $corrects{$function} += $preds{$function}[$index{$function}++][1];
     
   if ($index{$function} >= $span)
   {
    if ($index{$function} > $edge)
    {
     $corrects{$function} -= $preds{$function}[$index{$function}-$edge-1][1];
    }
    my $base = $index{$function} >= $edge?$edge:$index{$function};
    my $precision = $corrects{$function}/$base;
    if ($precision < 0.8)
    {
     if ($index{$function} == $span)
     {
      $threshold{$function} = 1;
     }
     else
     {
      $threshold{$function} = $preds{$function}[$index{$function}-$span-1][0];
     }
    }
   }
  }
 }

 foreach (keys %{$self->{ANNOTATIONS}})
 {
  if (!defined($threshold{$_}))
  {
   $threshold{$_} = 1;
  }
 }

 my $newannots = 0;
 %completed = ();
 foreach (sort {$predictions{$b} <=> $predictions{$a}} keys %predictions)
 {
  my @cols = split(/\|/, $_);
  my $id = $cols[0];
  my $function = $cols[1];
 
  if ($predictions{$_} < $threshold{$function})
  {
   $completed{$function} = 1;
   if (scalar(keys %completed) == $predfunctions) { last; }
   next;
  }

 if (
(!defined($self->{ENTITIES}{$id}{ANNOTATED}) || 
($fold != -1 && $self->{ENTITIES}{$id}{FOLDLABEL} == $fold)) && 
!defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$function})
)
  {
   $self->{ENTITIES}{$id}{ANNOTATIONS}{$function} = $level;
   $newannots ++;
  }
 }
 print "$newannots new annotations assigned.\n";
 return $newannots; 
}

sub ShowPrediction
{
 my($self)= $_[0];
 my($id)= $_[1];
 my($annot)= $_[2];
 my($weighting_technique)= $_[3];
 my($output_file)= $_[4];

 $self->computeWeights();
 my $score_manager = PP_Scoring->new($self);
 $score_manager->DoWeightedAvgS($id, $annot, $weighting_technique, $output_file);
 print ($self->{ENTITIES}{$id}{PREDICTIONS}{$annot})." OK";
}

sub PredictForUnknown() 
{
 ###########################################################
 # PredictForUnknown(
 # prediction_technique
 # weighting_technique
 # )
 # performs cross validation for annotation prediction
 ###########################################################

 my($self)= $_[0];
 my($prediction_technique)= $_[1];
 my($weighting_technique)= $_[2];
 my($output_file)= $_[3];
 my %predictions;

 open DAT, ">$output_file";
 $self->computeWeights($prediction_technique eq 'PROB_COMBINE'?1:0);
# exit;

 my $score_manager = PP_Scoring->new($self);
 my $validation = PP_Validation->new($self);

 my %predictions;
 my $id;
 foreach $id(sort {$a cmp $b} keys %{$self->{ENTITIES}})
 {
  $self->{ENTITIES}{$id}{PREDICTIONS} = ();
  if (!defined($self->{ENTITIES}{$id}{TEST}) || defined($self->{ENTITIES}{$id}{ANNOTATED}))
  {
   next;
  }

  if ($prediction_technique eq 'MAJORITY_VOTE')
  {
   $score_manager->DoMajorityVote($id);
  }
  elsif ($prediction_technique eq 'WEIGHTED_MAJORITY_VOTE')
  {
   $score_manager->DoWeightedMajorityVote($id, $weighting_technique);
  }
  elsif ($prediction_technique eq 'WEIGHTED_MAJORITY_VOTE_L2')
  {
   $score_manager->DoWeightedMajorityVoteL2($id, $weighting_technique);
  }
  elsif ($prediction_technique eq 'PROB_COMBINE')
  {
   $score_manager->DoProbCombine($id, -1);
  }
  elsif ($prediction_technique eq 'CHI_SQUARE')
  {
   $score_manager->DoChiSquare($id);
  }
  elsif ($prediction_technique eq 'WEIGHTED_AVG')
  {
   $score_manager->DoWeightedAvg($id, $weighting_technique);
  }
  elsif ($prediction_technique eq 'WEIGHTED_AVG_L1')
  {
   $score_manager->DoWeightedAvgL1($id, $weighting_technique);
  }
  else
  {
   die ('Invalid Prediction Technique');
  }

  foreach (keys %{$self->{ENTITIES}{$id}{PREDICTIONS}})
  {
#   if ($self->{ANNOTATIONS}{$_}{LEVEL} != 4) { next; }
#   @{$predictions{$_}[scalar(@{$predictions{$_}})]} = ($id, $self->{ENTITIES}{$id}{PREDICTIONS}{$_});
   print DAT "$id\t$_\t$self->{ENTITIES}{$id}{PREDICTIONS}{$_}\n";
   delete $self->{ENTITIES}{$id}{PREDICTIONS}{$_};
  }
 }
 close DAT;
 return;

 # compute rocs
 my %rocs;
 my $annotation;
 foreach $annotation(keys %predictions)
 {
  $rocs{$annotation} = $validation->ComputeTermROC($annotation, -1, 30);
#  print "$annotation\t$rocs{$annotation}\n";
 }

 my $annot;
 foreach $annot(sort {$rocs{$b} <=> $rocs{$a}} keys %rocs)
 {
#  if ($rocs{$annot} < 0.5) { last; }

  
  my %scores;
  my $lastval = -1;
  my $truestotal = 0;
  my $trues = 0;
  foreach (sort {$b->[1] <=> $a->[1]} @{$predictions{$annot}})
  {
   my ($id, $score) = @{$_};
   if ($score == 0) { last; }
   if (defined($self->{ENTITIES}{$id}{ANNOTATED}))
   {
    if (defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$annot})) { $trues++; }
#    else { $falses++; }
   }
   if ($score != $lastval && $lastval != -1)
   {
    $scores{$lastval} = $trues;
   }
   $lastval = $score;
  }
  if ($lastval > 0)
  {
   $scores{$lastval} = $trues;
  }
  $truestotal = $trues;

  foreach (sort {$b->[1] <=> $a->[1]} @{$predictions{$annot}})
  {
   my ($id, $score) = @{$_};
   if (defined($self->{ENTITIES}{$id}{ANNOTATED})) { next; }
   $trues = $scores{$score};
   print DAT "$id\t$annot\t".$self->{ANNOTATIONS}{$annot}{LEVEL}."\t$trues\t$truestotal\t$score\t$self->{ANNOTATIONS}{$annot}{NAME}\t$rocs{$annot}\n";   
  }
 }

 close DAT;
 return;
 my $span = 10;
 my $edge = $span*2 + 1;
 my %corrects;
 my %preds;
 my %threshold;
 my %index;

 foreach (sort {$predictions{$b} <=> $predictions{$a}} keys %predictions)
 {
  my @cols = split(/\|/, $_);
  my $id = $cols[0];
  my $function = $cols[1];

  if (defined($self->{ENTITIES}{$id}{ANNOTATED}) && !defined($threshold{$function}))
  {
   $corrects{$function} += defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$function})?1:0;
   $index{$function}++;
   if ($index{$function} >= 10)
   {
    my $recall = $corrects{$function}/$index{$function};
    if ($index{$function} == 10 && $recall < 0.8)
    {
     $threshold{$function} = 1;
    }
    elsif ($recall < 0.8)
    {
     $threshold{$function} = $preds{$function};
    }
   }
   $preds{$function} = $predictions{$_};
   next;

   $preds{$function}[$index{$function}][0] = $predictions{$_};
   $preds{$function}[$index{$function}][1] = defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$function})?1:0;
   $corrects{$function} += $preds{$function}[$index{$function}++][1];

   if ($index{$function} >= $span)
   {
    if ($index{$function} > $edge)
    {
     $corrects{$function} -= $preds{$function}[$index{$function}-$edge-1][1];
    }
    my $base = $index{$function} >= $edge?$edge:$index{$function};
    my $precision = $corrects{$function}/$base;
    if ($precision < 0.6)
    {
     if ($index{$function} == $span)
     {
      $threshold{$function} = 1;
     }
     else
     {
      $threshold{$function} = $preds{$function}[$index{$function}-$span-1][0];
     }
    }
   }
  }
 }

 my %fpredictions;
 my %fcorrects;

 foreach (keys %{$self->{ANNOTATIONS}})
 {
  $fpredictions{$_} = $fcorrects{$_} = 0;
  if (!defined($threshold{$_}))
  {
   $threshold{$_} = 1;
  } 
# print "$threshold{$_}\n";
 }
 
 foreach (sort {$predictions{$b} <=> $predictions{$a}} keys %predictions)
 {
  my @cols = split(/\|/, $_);
  my $id = $cols[0];
  my $function = $cols[1];

  if (defined($self->{ENTITIES}{$id}{ANNOTATED}))
  {
   my $correct = defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$function})?1:0;
   $fpredictions{$function} ++;
   $fcorrects{$function} += $correct;
  }

  if ($predictions{$_} < $threshold{$function}) { next; }

  if (!defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$function}) || $self->{ENTITIES}{$id}{ANNOTATIONS}{$function} > 1)
  {
   my %neighbours;
   my $neighbour;
   foreach $neighbour(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})
   {
    $neighbours{$neighbour} = $score_manager->GetPairWeight($id, $neighbour, $function);
   }
   my $sources;
   foreach $neighbour(sort {$neighbours{$b} <=> $neighbours{$a}} keys %neighbours)
   {
    if (!defined($self->{ENTITIES}{$neighbour}{ANNOTATIONS}{$function}))
    {
     next;
    }
    $sources .= "$neighbour|$neighbours{$neighbour}|";

    my $key = $self->GetInteractionKey($id, $neighbour);
    my $datatype;
    foreach $datatype(keys %{$self->{INTERACTIONS}{$key}{DATASOURCE}})
    {
     my $subtype;
     foreach $subtype(keys %{$self->{INTERACTIONS}{$key}{DATASOURCE}{$datatype}})
     {
      my $weight = 1;
      my $score;
      foreach $score(sort {$b <=> $a} @{$self->{INTERACTIONS}{$key}{DATASOURCE}{$datatype}{$subtype}{'scores'}})
      {
       my $reliability = 0;
       my $range;
       my $ranges = $self->{DATATYPE}{$datatype}{$subtype}{'functions'}{$function};
       foreach $range(sort {$ranges->{$b}{'bound'} <=> $ranges->{$a}{'bound'}} keys %{$ranges})
       {
        if ($score >= $ranges->{$range}{'bound'})
        {
         $reliability = $ranges->{$range}{RELIABILITY};
         last;
        }
       }
       $weight *= (1-$reliability);
      }
      $weight = 1 - $weight;
      $sources .= "$subtype=$weight;";
     }
    }
   }
   print "$id\t$function\t$predictions{$_}\t".($fcorrects{$function}>0?($fcorrects{$function}/$fpredictions{$function}):0)."\t$fpredictions{$function}\t$sources\t$self->{ANNOTATIONS}{$function}{NAME}\n";
  }
 }

}

sub ValidateFunctionalWeight() 
{
 my($self)= $_[0];
 my($weighting_technique) = $_[1];

 my $directcount = 0;
 my $indirectcount = 0;
 my @predictions;
 my @l2predictions;
 my @allpredictions;
 my $index = 0;
 my $l2index = 0;
 my $allindex = 0;

 $self->computeWeights();

 my $weight_manager = PP_Weighting->new($self);

 my $id;
 my %completed = ();
 foreach $id (keys %{$self->{ENTITIES}})
 {  
  if (!defined($self->{ENTITIES}{$id}{ANNOTATED}) || $self->{ENTITIES}{$id}{LEVEL} < $self->{SIM_LEVEL}) 
  {
   $completed{$id} = 1;  
   next; 
  }
  my %neighbours;
  my $neighbour;
  foreach $neighbour(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})
  {
   # l1 neighbours
   if (!defined($completed{$neighbour}) && defined($self->{ENTITIES}{$neighbour}{ANNOTATED}) && $self->{ENTITIES}{$neighbour}{LEVEL} >= $self->{SIM_LEVEL})
   {
    if (!defined($neighbours{$neighbour}))
    {
     $neighbours{$neighbour}{'weight'} = $weight_manager->ComputeLocalWeight($id, $neighbour, $weighting_technique);
    }   
    $neighbours{$neighbour}{'l1'} = 1;
   }
   next;    
   my $l2neighbour;
   foreach $l2neighbour(keys %{$self->{ENTITIES}{$neighbour}{NEIGHBOURS}})
   {
    if (!defined($self->{ENTITIES}{$l2neighbour}{ANNOTATED})) { next; }
    # l2 neighbours
    if ($l2neighbour ne $id && !defined($completed{$l2neighbour}) && $self->{ENTITIES}{$l2neighbour}{LEVEL} >= $self->{SIM_LEVEL})
    {
     if (!defined($neighbours{$l2neighbour}))
     {
      $neighbours{$l2neighbour}{'weight'} = $weight_manager->ComputeLocalWeight($id, $l2neighbour, $weighting_technique);
     }
     $neighbours{$l2neighbour}{'l2'} = 1;
    }
   }
  }
 
  foreach $neighbour(keys %neighbours)
  {
   my $sim = $self->HasSimilarity($id, $neighbour, $self->{SIM_LEVEL});
   my $weight = $neighbours{$neighbour}{'weight'};
#   if ($weight < 0.05) { next; }

   $allpredictions[$allindex][0] = $weight;
   $allpredictions[$allindex++][1] = $sim>0?1:0;

   if ($neighbours{$neighbour}{'l1'} )
   {
#    if ($neighbours{$neighbour}{'l2'})
#     $x += $sim;
#     $y += $weight;
#     $xx += $sim*$sim;
#     $yy += $weight*$weight;
#     $xy += $sim*$weight;
     $predictions[$index][0] = $weight;
     $predictions[$index++][1] = $sim>0?1:0;
#     $directcount++;
   }
   elsif ($neighbours{$neighbour}{'l2'})
   {
     $l2predictions[$l2index][0] = $weight;
     $l2predictions[$l2index++][1] = $sim>0?1:0;
#     $l2x += $sim;
#     $l2y += $weight;
#     $l2xx += $sim*$sim;
#     $l2yy += $weight*$weight;
#     $l2xy += $sim*$weight;
#     $indirectcount++;
   }
  }  
  $completed{$id} = 1;
 }

# print "$directcount\t".scalar(@predictions)."\t".scalar(@predictions)/$directcount."\n";
 $self->GetPearson(\@predictions);
 $self->GetPearson(\@l2predictions);
 $self->GetPearson(\@allpredictions);
}

sub GetPearson
{
 my($self, $predictions) = @_;

 my $x = 0;
 my $xx = 0;
 my $y = 0;
 my $yy = 0;
 my $xy = 0;
 my $count = 0;
 my $predicts = 0;
 my $corrects= 0;
 my $threshold = 0.02;
 my $recall = 0;
 my $lastvalue = -1;
 my $total = scalar(@{$predictions});
 foreach (sort {$$b[0] <=> $$a[0]} @{$predictions})
 {
  if ($lastvalue != $$_[0])
  {
   if ($lastvalue != -1)
   {
     my $sim = $corrects/$predicts;
     $x += $sim;
     $y += $lastvalue;
     $xx += $sim*$sim;
     $yy += $lastvalue*$lastvalue;
     $xy += $sim*$lastvalue;
     $count++;
   }
   if ($lastvalue != -1 && ($recall = $predicts/$total) > $threshold)
   {
    print "$lastvalue\t$predicts\t$recall\t".$corrects/$predicts."\n";
    while($threshold < $recall)
    {
     $threshold += 0.02;
    }
   }
   $lastvalue = $$_[0];
  }

  if ($$_[1] == 1)
  {
   $corrects++;
  }
  $predicts++;
 }
   if ($lastvalue != -1)
   {
     my $sim = $corrects/$predicts;
     $x += $sim;
     $y += $lastvalue;
     $xx += $sim*$sim;
     $yy += $lastvalue*$lastvalue;
     $xy += $sim*$lastvalue;
     $count ++;
   }
 
# if (($recall = $predicts/$total) >= $threshold)
 {
  print "$lastvalue\t$predicts\t$recall\t".$corrects/$predicts."\n";
 }
 my $pearson = ($xy - $x*$y/$count)/sqrt(($xx - $x*$x/$count)*($yy - $y*$y/$count));
 print "$pearson\n";


}

sub DetectL2OnlySimilarity() 
{
 my($self)= $_[0];
 my $annot_level = $self->{SIM_LEVEL};

 my $id;
 my $total = 0;
 my $L1count = 0;
 my $L2count = 0;
 my $L3count = 0;
 my %functions;
 my @functioncount;
 my @count;
 my @sim;
 my @fcount;
 my @fsim;
 my %completed;
 my $weight_manager = PP_Weighting->new($self);

 $self->computeWeights();

 foreach $id(keys %{$self->{ENTITIES}})
 {
  if (!scalar(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})
|| $self->{ENTITIES}{$id}{LEVEL} < $annot_level 
|| !defined($self->{ENTITIES}{$id}{TEST}) 
  )
  {
   next;
  }
  $total ++;
  my @similarity;
  my $neighbour;
  my %neighbours;
  my %nextlevelneighbours;
  my @queue = ();
  my $level = 0;

  unshift(@queue, $id);
  while(scalar(@queue))
  {
   my $orf = pop @queue;

   foreach $neighbour(keys %{$self->{ENTITIES}{$orf}{NEIGHBOURS}})
   {
    if ($neighbour eq $id)
    {
     next;
    }
    $nextlevelneighbours{$neighbour} = 1;
    $neighbours{$neighbour}[$level] = 1;
   }

   if (!scalar(@queue))
   {
 # print "=>$level ";
    if (++$level == 2)
    {
     last;
    }
    @queue = keys %nextlevelneighbours;
    %nextlevelneighbours = ();
   }
  }

  foreach $neighbour (keys %neighbours)
  { 
   # ignore unannotated neighbours
   if ($self->{ENTITIES}{$neighbour}{LEVEL} < $annot_level) { next; }
   if (!defined($self->{ENTITIES}{$neighbour}{ANNOTATED})) { next; }
   my $filtered = ($weight_manager->ComputeLocalWeight($id, $neighbour, 'FSWEIGHT') < 0.2)?0:1;

   if ($neighbours{$neighbour}[0] && !$neighbours{$neighbour}[1])
   {
    if (!defined($completed{$neighbour}))
    {
     $count[0]++;
     $fcount[0]+=$filtered;     
    }
   
    if ($self->HasSimilarity($id, $neighbour, $annot_level))
    {
     if (!defined($completed{$neighbour}))
     { 
      $sim[0]++; 
      if ($filtered) { $fsim[0]++; }
     }
     $similarity[0]++;

     foreach (keys %{$self->{ENTITIES}{$neighbour}{ANNOTATIONS}})
     {
      if (defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$_}))
      {
       $functions{$_}[0]++;
       $functioncount[0]++;
      }
     }
    }
   }
   elsif ($neighbours{$neighbour}[1] && !$neighbours{$neighbour}[0])
   {
    if (!defined($completed{$neighbour}))
    {
     $count[1]++;
     $fcount[1]+=$filtered;     
    }

    if ($self->HasSimilarity($id, $neighbour, $annot_level))
    {
     if (!defined($completed{$neighbour}))
     { 
      $sim[1]++; 
      if ($filtered) { $fsim[1]++; }
     }
     $similarity[1]++;

     foreach (keys %{$self->{ENTITIES}{$neighbour}{ANNOTATIONS}})
     {
      if (defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$_}))
      {
       $functions{$_}[1]++;
       $functioncount[1]++;
      }
     }
    }
   }
   elsif ($neighbours{$neighbour}[0] && $neighbours{$neighbour}[1])
   {
    if (!defined($completed{$neighbour}))
    {
     $count[2]++;
     $fcount[2]+=$filtered;     
    }

    if ($self->HasSimilarity($id, $neighbour, $annot_level))
    {
     if (!defined($completed{$neighbour}))
     { 
       $sim[2]++; 
       if ($filtered) { $fsim[2]++; }
     }
     $similarity[2]++;

     foreach (keys %{$self->{ENTITIES}{$neighbour}{ANNOTATIONS}})
     {
      if (defined($self->{ENTITIES}{$id}{ANNOTATIONS}{$_}))
      {
       $functions{$_}[2]++;
       $functioncount[2]++;
      }
     }
    }
   }
  }

  foreach (keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})
  {
    $functions{$_}[3]++;
  }

  if ($similarity[0] && !$similarity[1] && !$similarity[2])
  {
   $L1count++;
#   print "$id\tL1-L2\t".scalar(keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})."\t".scalar(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})."\n";
   foreach (keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})
   {
    $functions{$_}[4]++;
   }
  }
  elsif (!$similarity[0] && $similarity[1] && !$similarity[2])
  {
   $L2count++;
#   print "$id\tL2-L1\t".scalar(keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})."\t".scalar(keys %{$self->{ENTITIES}{$id}{NEIGHBOURS}})."\n";
   foreach (keys %{$self->{ENTITIES}{$id}{ANNOTATIONS}})
   {
    $functions{$_}[5]++;
   }
  }elsif (($similarity[0] && $similarity[1]) || $similarity[2])
  {
   $L3count++;
  # print "$id\n";
  }

  $completed{$id} = 1;
 }

 print "All Pairs:\n";
 print "L1 only pairs with similarity: $sim[0] / $count[0] = ".($count[0]>0?($sim[0]/$count[0]):0)."\n";
 print "L2 only pairs with similarity: $sim[1] / $count[1] = ".($count[1]>0?($sim[1]/$count[1]):0)."\n";
 print "L1 and L2 pairs with similarity: $sim[2] / $count[2] = ".($count[2]>0?($sim[2]/$count[2]):0)."\n";
 print "Pairs with FS-Weight >= 0.2:\n";
 print "L1 only pairs with similarity: $fsim[0] / $fcount[0] = ".($fcount[0]>0?($fsim[0]/$fcount[0]):0)."\n";
 print "L2 only pairs with similarity: $fsim[1] / $fcount[1] = ".($fcount[1]>0?($fsim[1]/$fcount[1]):0)."\n";
 print "L1 and L2 pairs with similarity: $fsim[2] / $fcount[2] = ".($fcount[2]>0?($fsim[2]/$fcount[2]):0)."\n";
 print "Entities with neighbours and annotation:\n";
 print "Fraction of annotated entities with L1 only similarity: $L1count / $total = ".($total>0?($L1count/$total):0)."\n";
 print "Fraction of annotated entities with L2 only similarity: $L2count / $total = ".($total>0?($L2count/$total):0)."\n";
 print "Fraction of annotated entities with L1 and L2 similarity: $L3count / $total = ".($total>0?($L3count/$total):0)."\n";
 print "Fraction of annotated entities with L1 or L2 similarity: ($L1count + $L2count + $L3count) / $total = ".($total>0?(($L1count + $L2count + $L3count)/$total):0)."\n";
 return;
 foreach (keys %functions)
 {
  print "$_";
  print "\t$self->{ANNOTATIONS}{$_}{NAME}";
  print "\t".(defined($functions{$_}[0])?$functions{$_}[0]:0);
  print "\t".(defined($functions{$_}[1])?$functions{$_}[1]:0);
  print "\t".(defined($functions{$_}[2])?$functions{$_}[2]:0);
  print "\t".(defined($functions{$_}[0])?$functions{$_}[0]/$functioncount[0]:0);
  print "\t".(defined($functions{$_}[1])?$functions{$_}[1]/$functioncount[1]:0);
  print "\t".(defined($functions{$_}[2])?$functions{$_}[2]/$functioncount[2]:0);

  print "\t".(defined($functions{$_}[3])?$functions{$_}[3]:0);
  print "\t".(defined($functions{$_}[4])?$functions{$_}[4]:0);
  print "\t".(defined($functions{$_}[5])?$functions{$_}[5]:0);
  print "\t".(defined($functions{$_}[4])?$functions{$_}[4]/$functions{$_}[3]:0);
  print "\t".(defined($functions{$_}[5])?$functions{$_}[5]/$functions{$_}[3]:0);

  print "\n";
 } 
}

return 1;
