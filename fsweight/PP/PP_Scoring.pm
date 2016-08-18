package PP_Scoring;
use strict;
use PP::PP_Weighting;
use constant INFINITY          => 9999999;
use POSIX;

sub new()
{
 ###########################################################
 # new(
 #  pp_engine
 # )
 # weighting functions
 ###########################################################

 my($type) = $_[0];
 my($self) = {};

 $self->{'engine'} = $_[1];

 bless($self, $type);
 return($self);
}

sub DoMajorityVote()
{
 ###########################################################
 # DoMajorityVote(
 # )
 # predict annotation of entities
 ###########################################################

 my($self) = $_[0];
 my($id) = $_[1];
 my($engine) = $self->{'engine'};
 my %neighbours;
 my $neighbour;
 my $annotation;
 my $rank = 0;
 my $lastscore = 0;
 my $repeats = 1;

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})
 {
  foreach $annotation(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}})
  {
   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} += 1;
  }
 }

 foreach $annotation(sort
        {$engine->{ENTITIES}{$id}{PREDICTIONS}{$b} <=>
        $engine->{ENTITIES}{$id}{PREDICTIONS}{$a}}
        keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
 {
  if ($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} < $lastscore || $lastscore == 0)
  {
   $rank += $repeats;
   $repeats = 1;
   $lastscore = $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation};
  }else
  {
   $repeats++;
  }
  $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} = 1/$rank;
 }
}

sub GetPairWeight()
{
 my ($self, $id_a, $id_b, $annotation) = @_;
 my($engine) = $self->{'engine'};

# return 1;

 my $key = $engine->GetInteractionKey($id_a, $id_b); 
 my $weight = 1;
 my $datatype;
 my $subtype;  
 my $score;

 foreach $datatype(keys %{$engine->{INTERACTIONS}{$key}{DATASOURCE}})
 {
  foreach $subtype(keys %{$engine->{INTERACTIONS}{$key}{DATASOURCE}{$datatype}})
  {
   my $score_delta = $engine->{DATATYPE}{$datatype}{$subtype}{DELTA};
   my $minvalue = $engine->{DATATYPE}{$datatype}{$subtype}{MINVAL};
   foreach $score(sort {$b <=> $a} @{$engine->{INTERACTIONS}{$key}{DATASOURCE}{$datatype}{$subtype}{'scores'}})
   {
     if ($score<$minvalue) { next; }
     my $reliability = 0;
     my $bin = 1;
     if ($score_delta > 0) { $bin = floor(($score-$minvalue)/$score_delta)+1; }
     if ($bin > $engine->{NUMBINS}) { $bin = $engine->{NUMBINS}; }
     $reliability = $engine->{DATATYPE}{$datatype}{$subtype}{RANGE}{$annotation}{$bin};
     $weight *= (1-$reliability);
   }
  }
 }
 return 1-$weight;
}

sub DoProbCombine()
{
 my($self, $id) = @_;
 my($engine) = $self->{'engine'};
 
 my $annotation;
 my %base = ();
 my $neighbour;
 my $l2neighbour;
 my @neighbours;
 my @l2neighbours;
 my $total = 1;

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})
 {  
  if (!scalar(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}}))
  {   
   next;
  }

  push @neighbours, $neighbour;
  foreach $annotation(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}})
  {
   my $weight = $self->GetPairWeight($id, $neighbour, $annotation);
   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} += $weight;
  }
 }

 foreach $neighbour(@neighbours)
 {  
  foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
  {
#   if (!defined($engine->{ANNOTATIONS}{$annotation}{INF})) { next; }
   my $weight = $self->GetPairWeight($id, $neighbour, $annotation);
   $base{$annotation} += $weight;
  }
 }
  
 foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
 {
  if ($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} > 0)
  {
   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} /= ($base{$annotation} + 1);
  } 
  else
  { 
   delete $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation};
  }
 }

}

sub DoProbCombineNoWeight()
{
 my($self, $id) = @_;
 my($engine) = $self->{'engine'};
 
 my $annotation;
 my %base = ();
 my $neighbour;
 my @neighbours;

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})
 {  
  if (!scalar(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}}))
  {   
   next;
  }

  push @neighbours, $neighbour;

  foreach $annotation(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}})
  {
#   if (!defined($engine->{ANNOTATIONS}{$annotation}{INF})) { next; }
   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} ++;
  }  
 }

 foreach $neighbour(@neighbours)
 {  
  foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
  {
#   if (!defined($engine->{ANNOTATIONS}{$annotation}{INF})) { next; }
   $base{$annotation} ++;
  }
 }
  
 foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
 {
  if ($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} > 0)
  {
   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} /= ($base{$annotation} + 1);
  } 
  else
  { 
   delete $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation};
  }
 }

}

sub DoProbCombineSimple()
{
 my($self, $id) = @_;
 my($engine) = $self->{'engine'};
 
 my $annotation;
 my %base = ();
 my $neighbour;
 my @neighbours;

# foreach $annotation(keys %{$engine->{ANNOTATIONS}})
# {
#  if (!defined($engine->{ANNOTATIONS}{$annotation}{INF})) { next; }
#  if ($engine->{ANNOTATIONS}{$annotation}{'frequency'})
#  {
#   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} = $engine->{ANNOTATIONS}{$annotation}{'frequency'};#*$engine->{'g_reliability'};
#  }
# }

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})
 {  
  if (!scalar(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}}))
  {   
   next;
  }

  push @neighbours, $neighbour;
  my $key = $engine->GetInteractionKey($id, $neighbour);
  my $weight = $engine->{INTERACTIONS}{$key}{RELIABILITY};

  foreach $annotation(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}})
  {
   if (!defined($engine->{ANNOTATIONS}{$annotation}{INF})) { next; }
   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} += $weight;
  }  
 }

 foreach $neighbour(@neighbours)
 {  
  my $key = $engine->GetInteractionKey($id, $neighbour);
  my $weight = $engine->{INTERACTIONS}{$key}{RELIABILITY};
  foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
  {
   $base{$annotation} += $weight;
  }
 }
  
 foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
 {
  if ($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} > 0)
  {
   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} /= ($base{$annotation} + 1);
  } 
  else
  { 
   delete $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation};
  }
 }

}

sub DoProbCombineBayes()
{
 my($self, $id) = @_;
 my($engine) = $self->{'engine'};
 
 my $annotation;
 my %base = ();
 my $neighbour;
 my @neighbours;

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})
 {  
  if (!scalar(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}}))
  {   
   next;
  }

  push @neighbours, $neighbour;
#  my $key = $engine->GetInteractionKey($id, $neighbour);
#  my $weight = $engine->{INTERACTIONS}{$key}{RELIABILITY};

  foreach $annotation(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}})
  {
   if (!defined($engine->{ANNOTATIONS}{$annotation}{INF})) { next; }
   my $weight = $self->GetPairWeight($id, $neighbour, $annotation);

   if (!defined($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation}))
   {
    $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} = 1;
   }
   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} *= (1-$weight);
  }  
 }
  
 foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
 {
  $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} = 1 - $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation};
  if ($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} > 0)
  {
  } 
  else
  { 
   delete $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation};
  }
 }
}

sub DoWeightedMajorityVote()
{
 ###########################################################
 # DoWeightedMajorityVote(
 # )
 # predict annotation of entities
 ###########################################################

 my($self) = $_[0];
 my($id) = $_[1];
 my($weighting_technique) = $_[2];
 my($engine) = $self->{'engine'};
 my %neighbours;
 my $neighbour;
 my $annotation;
 my $rank = 0;
 my $lastscore = 0;
 my $repeats = 1;
 my $weight_manager = PP_Weighting->new($engine);

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})
 {
  my $weight = $weight_manager->ComputeLocalWeight($id, $neighbour, $weighting_technique);
  foreach $annotation(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}})
  {
   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} += $weight;  
  }
 }

 foreach $annotation(sort
        {$engine->{ENTITIES}{$id}{PREDICTIONS}{$b} <=>
        $engine->{ENTITIES}{$id}{PREDICTIONS}{$a}}
        keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
 {
  if ($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} < $lastscore || $lastscore == 0)
  {
   $rank += $repeats;
   $repeats = 1;
   $lastscore = $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation};
  }else
  {
   $repeats++;
  }
  $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} = 1/$rank;
 }

}

sub DoWeightedMajorityVoteL2()
{
 ###########################################################
 # DoWeightedMajorityVote(
 # )
 # predict annotation of entities
 ###########################################################

 my($self) = $_[0];
 my($id) = $_[1];
 my($weighting_technique) = $_[2];
 my($engine) = $self->{'engine'};
 my %neighbours;
 my $neighbour;
 my $annotation;
 my $rank = 0;
 my $lastscore = 0;
 my $repeats = 1;
 my $weight_manager = PP_Weighting->new($engine);

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})
 {
  $neighbours{$neighbour}++;
  my $neighbour1;
  foreach $neighbour1(keys %{$engine->{ENTITIES}{$neighbour}{NEIGHBOURS}})
  {
   if ($neighbour1 eq $id) { next; }
   $neighbours{$neighbour1}++;
  }  
 }

 foreach $neighbour(keys %neighbours)
 {
  my $weight = $weight_manager->ComputeLocalWeight($id, $neighbour, $weighting_technique);
  foreach $annotation(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}})
  {
   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} += $weight*$neighbours{$neighbour};  
  }
 }

 foreach $annotation(sort
        {$engine->{ENTITIES}{$id}{PREDICTIONS}{$b} <=>
        $engine->{ENTITIES}{$id}{PREDICTIONS}{$a}}
        keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
 {
  if ($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} < $lastscore || $lastscore == 0)
  {
   $rank += $repeats;
   $repeats = 1;
   $lastscore = $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation};
  }else
  {
   $repeats++;
  }
  $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} = 1/$rank;
 }

}

sub DoChiSquare()
{
 ###########################################################
 # DoChiSquare(
 # )
 # predict annotation of entities
 ###########################################################

 my($self) = $_[0];
 my($id) = $_[1];
 my($engine) = $self->{'engine'};
 my $neigbour;
 my $annotation;
 my $rank = 0;
 my $lastscore = 0;
 my $repeats = 1;
 my $Nu = 0; #= scalar(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}});

 foreach $neigbour(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})
 {
  if (!scalar(keys %{$engine->{ENTITIES}{$neigbour}{ANNOTATIONS}})) { next; }
  $Nu++;
  foreach $annotation(keys %{$engine->{ENTITIES}{$neigbour}{ANNOTATIONS}})
  {
   if ($engine->{ANNOTATIONS}{$annotation}{'frequency'}==0) { next; }
   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation}++;
  }
 }
 foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
 {
  my $expected_occ = $Nu * $engine->{ANNOTATIONS}{$annotation}{'frequency'};
  $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} -= $expected_occ;
  $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} *= $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation};
  $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} /= $expected_occ;
 }
 foreach $annotation(sort
        {$engine->{ENTITIES}{$id}{PREDICTIONS}{$b} <=>
        $engine->{ENTITIES}{$id}{PREDICTIONS}{$a}}
        keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
 {
  if ($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} < $lastscore || $lastscore == 0)
  {
   $rank += $repeats;
   $repeats = 1;
   $lastscore = $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation};
  }else
  {
   $repeats++;
  }
  $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} = 1/$rank;
 }
}

sub DoWeightedAvg()
{
 ###########################################################
 # DoWeightedAvg(
 # )
 # predict annotation of entities
 ###########################################################

 my($self) = $_[0];
 my($id) = $_[1];
 my($weighting_technique) = $_[2];
 my($engine) = $self->{'engine'};

 my $annotation;
 my %neighbours =();
 my %scores =();
 my $neighbour;

 my %weights;
 my $neighbour1;
 my $weight_manager = PP_Weighting->new($engine);

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})
 {
  if (!defined($neighbours{$neighbour}))
  {
   $weights{$id}{$neighbour} = $scores{$neighbour} = $weight_manager->ComputeLocalWeight($id, $neighbour,$weighting_technique );
  }
  $neighbours{$neighbour}++;

  my $weight1 = $weights{$id}{$neighbour};

  #next;

  foreach $neighbour1(keys %{$engine->{ENTITIES}{$neighbour}{NEIGHBOURS}})
  {
    if ($neighbour1 eq $id) {next;}
    if (!defined($neighbours{$neighbour1}))
    {
     $weights{$id}{$neighbour1} = $scores{$neighbour1} = $weight_manager->ComputeLocalWeight($id, $neighbour1, $weighting_technique);
    }
    $neighbours{$neighbour1}++;

    my $weight2 = $weight1*$weight_manager->ComputeLocalWeight($neighbour, $neighbour1, $weighting_technique);
    if ($scores{$neighbour1} < $weight2)
    {
     $scores{$neighbour1} = $weight2;
    }
  }
 }

 my $base = 1;
# foreach $annotation(keys %{$engine->{ANNOTATIONS}})
# foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
# {
#  if ($engine->{ANNOTATIONS}{$annotation}{'frequency'})
#  {
#   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} = $engine->{ANNOTATIONS}{$annotation}{'frequency'}#*$engine->{'g_reliability'};
#  }
# }

 foreach $neighbour(keys %neighbours)
 {
   foreach $annotation(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}})
   {
    $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} += $scores{$neighbour}*$neighbours{$neighbour};
   }
   if (scalar(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}}))
   {   
    $base += $scores{$neighbour}*$neighbours{$neighbour};
   }
 }

 foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
 {
  $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} /= $base;
 }

}

sub DoWeightedAvgL1()
{
 ###########################################################
 # DoWeightedAvgL1(
 # )
 # predict annotation of entities
 ###########################################################

 my($self) = $_[0];
 my($id) = $_[1];
 my($weighting_technique) = $_[2];
 my($engine) = $self->{'engine'};

 my $annotation;
 my %neighbours =();
 my %scores =();
 my $neighbour;

 my $neighbour1;
 my $weight_manager = PP_Weighting->new($engine);

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})
 {
  if (!defined($neighbours{$neighbour}))
  {
   $scores{$neighbour} = $weight_manager->ComputeLocalWeight($id, $neighbour,$weighting_technique );
  }
  $neighbours{$neighbour}++;

  my $weight1 = $engine->{WEIGHTS}{$id}{$neighbour};
  if ($scores{$neighbour} < $weight1)
  {
   $scores{$neighbour} = $weight1;
  }
 }

 my $base = 1;
# foreach $annotation(keys %{$engine->{ANNOTATIONS}})
# foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
# {
#  if ($engine->{ANNOTATIONS}{$annotation}{'frequency'})
#  {
#   $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} = $engine->{ANNOTATIONS}{$annotation}{'frequency'};#*$engine->{'g_reliability'};
#  }
# }

 foreach $neighbour(keys %neighbours)
 {
   foreach $annotation(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}})
   {
    $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} += $scores{$neighbour}*$neighbours{$neighbour};;
   }
   if (scalar(keys %{$engine->{ENTITIES}{$neighbour}{ANNOTATIONS}}))
   {   
    $base += $scores{$neighbour}*$neighbours{$neighbour};
   }
 }

 foreach $annotation(keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
 {
  $engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation} /= $base;
 }

}

sub DoFunctionalFlow()
{
 my($self) = $_[0];
 my($id) = $_[1];
 my($engine) = $self->{'engine'};

 my $i;
 my $entity;
 my $annotation;

 foreach $entity(keys %{$engine->{ENTITIES}})
 {
  delete $engine->{ENTITIES}{$entity}{'reservoirs'};
  foreach $annotation(keys %{$engine->{ENTITIES}{$entity}{ANNOTATIONS}})
  {
   $engine->{ENTITIES}{$entity}{'reservoirs'}{$annotation} = INFINITY;
  }
 }

 for ($i=0; $i=6; $i++)
 {


 }

}

return 1;
