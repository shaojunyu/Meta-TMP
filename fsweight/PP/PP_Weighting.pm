package PP_Weighting;
use strict;
use PP::PP_Math;
#use Time::HiRes qw( usleep ualarm gettimeofday tv_interval );


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

sub ComputeLocalWeight()
{
 my($self) = $_[0];
 my ($id_a) = $_[1];
 my ($id_b) = $_[2];
 my($weighting_technique) = $_[3];
 my $engine = $self->{'engine'};

 my $key = $id_a lt $id_b?"$id_a|$id_b":"$id_b|$id_a";
 if (defined($engine->{WEIGHTS}{$key}))
 {
  return $engine->{WEIGHTS}{$key};
 }
 elsif ($weighting_technique eq 'GEOMETRIC')
 {
  $engine->{WEIGHTS}{$key} = $self->ComputeGeometric($id_a, $id_b);
 }
 elsif ($weighting_technique eq 'CD_DIST')
 {
  $engine->{WEIGHTS}{$key} = $self->ComputeCDDist($id_a, $id_b);
 }
 elsif ($weighting_technique eq 'FSWEIGHT')
 {
  $engine->{WEIGHTS}{$key} = $self->ComputeFSWeight($id_a, $id_b);
 }
 elsif ($weighting_technique eq 'FSLOCAL')
 {
  $engine->{WEIGHTS}{$key} = $self->ComputeFSLocal($id_a, $id_b);
 }
 else
 {
  die ('Invalid Weighting Technique - '.$weighting_technique);
 }
 return $engine->{WEIGHTS}{$key};
}

sub ComputeNeighbourWeights()
{
 my($self) = $_[0];
 my ($id) = $_[1];
 my($weighting_technique) = $_[2];

 if ($weighting_technique eq 'PATHRATIO')
 {
  return $self->ComputePathRatio($id);
 }
 elsif ($weighting_technique eq 'PATHRATIONEW')
 {
  return $self->ComputePathRatioNew($id);
 }
 else
 {
  die ('Invalid Weighting Technique - '.$weighting_technique);
 }
}

sub ComputeWeight()
{
 my($self) = $_[0];
 my($weighting_technique) = $_[1];

 my $engine = $self->{'engine'};
 my $weight = 0;

# my $starttime = [gettimeofday];

 if ($weighting_technique eq 'FSGLOBAL')
 {
  $self->ComputeFSGlobal();
 }
 elsif ($weighting_technique eq 'PATHRATIONEW'
    || $weighting_technique eq 'PATHRATIO')
 {
  my $id;
  foreach $id(sort {$a cmp $b} keys %{$engine->{ENTITIES}})
  {
   if (scalar(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}}))
   {
    $self->ComputeNeighbourWeights($id, $weighting_technique);
#     print "$id\n";
   }
  }
 }
 else
 {
  my $key;
  my $id_a;
  my $id_b;

  print "a";
  foreach $id_a(keys %{$engine->{ENTITIES}})
  {
   foreach $id_b(keys %{$engine->{ENTITIES}{$id_a}{NEIGHBOURS}})
   {
    $engine->{WEIGHTS}{$id_a}{$id_b} = $engine->{WEIGHTS}{$id_b}{$id_a} = $self->ComputeLocalWeight($id_a, $id_b, $weighting_technique);
   }
  }

 }
 return;
# my $elapsed = tv_interval($starttime, [gettimeofday]);
 #print "$elapsed\n";

 my $key;
 my $id_a;
 my $id_b;
 open(DAT, ">weights.txt");
 foreach $key(keys %{$engine->{INTERACTIONS}})
 {
   ($id_a, $id_b) = split(/\|/, $key);
   print DAT "$id_a\t$id_b\t$engine->{WEIGHTS}{$id_a}{$id_b}\n";
 }
 close(DAT);
}

sub ComputeCDDist()
{
 ###########################################################
 # ComputeCDDist(
 # )
 # compute 1 - CD-Distance between pairs of entities
 ###########################################################

 my($self) = $_[0];
 my($id_a) = $_[1];
 my($id_b) = $_[2];

 my $engine = $self->{'engine'};

 my $neighbour;
 my $Na = 0;
 my $Nb = 0;
 my $shared = 0;

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id_a}{NEIGHBOURS}})
 {
  if ($neighbour eq $id_b) { next; }
  my $key = $engine->GetInteractionKey($id_a, $neighbour);

  if (defined($engine->{ENTITIES}{$id_b}{NEIGHBOURS}{$neighbour}))
  {
   my $key2 = $engine->GetInteractionKey($id_b, $neighbour);
   $shared ++; #= $engine->{INTERACTIONS}{$key}{RELIABILITY}*$engine->{INTERACTIONS}{$key2}{RELIABILITY};
  }else{
   $Na ++; #+= $engine->{INTERACTIONS}{$key}{RELIABILITY};
  }
 }

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id_b}{NEIGHBOURS}})
 {
  if ($neighbour eq $id_a) { next; }
  my $key = $engine->GetInteractionKey($id_b, $neighbour);
  if (!defined($engine->{ENTITIES}{$id_a}{NEIGHBOURS}{$neighbour}))
  { 
   $Nb ++; #+= $engine->{INTERACTIONS}{$key}{RELIABILITY};
  }
 }

 if (defined($engine->{ENTITIES}{$id_a}{NEIGHBOURS}{$id_b}))
 {
  my $key = $engine->GetInteractionKey($id_a, $id_b);
  $shared ++; #+= $engine->{INTERACTIONS}{$key}{RELIABILITY};
 }else{
  $Na ++;
  $Nb ++;
 }

 if ($shared == 0)
 {
  return 0;
 }
 return 1 - ($Na+$Nb)/($Na+$Nb+$shared*2);
}

sub ComputeGeometric()
{
 ###########################################################
 # ComputeGeometric(
 # )
 # compute geometric weight between pairs of entities
 ###########################################################

 my($self) = $_[0];
 my($id_a) = $_[1];
 my($id_b) = $_[2];

 my $engine = $self->{'engine'};

 my $neighbour;
 my $Na = 0;
 my $Nb = 0;
 my $shared = 0;

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id_a}{NEIGHBOURS}})
 {
  if ($neighbour eq $id_b)
  {
   next;
  }
  my $key = $engine->GetInteractionKey($id_a, $neighbour);
  $Na ++; #+= $engine->{INTERACTIONS}{$key}{RELIABILITY};
  if (defined($engine->{ENTITIES}{$id_b}{NEIGHBOURS}{$neighbour}))
  {
   my $key2 = $engine->GetInteractionKey($id_b, $neighbour);
   $shared ++; #+= $engine->{INTERACTIONS}{$key}{RELIABILITY}*$engine->{INTERACTIONS}{$key2}{RELIABILITY};
  }
 }
 foreach $neighbour(keys %{$engine->{ENTITIES}{$id_b}{NEIGHBOURS}})
 {
  if ($neighbour eq $id_a)
  {
   next;
  }
  my $key = $engine->GetInteractionKey($id_b, $neighbour);
  $Nb ++; #+= $engine->{INTERACTIONS}{$key}{RELIABILITY};
 }
 if ($shared == 0)
 {
  return 0;
 }
 return ($shared*$shared)/($Na*$Nb);
}

sub ComputeFSWeight()
{
 ###########################################################
 # ComputeFSWeight(
 # )
 # compute FS weight between pairs of entities
 ###########################################################

 my($self) = $_[0];
 my($id_a) = $_[1];
 my($id_b) = $_[2];
 my($ignoreweights) = defined($_[3])?$_[3]:0;

 my $engine = $self->{'engine'};

 my $neighbour;
 my $base1 = 1;
 my $base2 = 1;
 my $sim = 0;
 my $weight;
 my $weight2;
 my $key;

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id_a}{NEIGHBOURS}})
 {
  $key = $engine->GetInteractionKey($id_a, $neighbour);
  $weight = $ignoreweights?1:$engine->{INTERACTIONS}{$key}{RELIABILITY};

  if ($neighbour eq $id_b)
  {
   $base1 += (1 - $weight);
   $sim += $weight;
  }
  elsif (defined($engine->{ENTITIES}{$id_b}{NEIGHBOURS}{$neighbour}))
  {
   $key = $engine->GetInteractionKey($id_b, $neighbour);
   $weight2 = $ignoreweights?1:$engine->{INTERACTIONS}{$key}{RELIABILITY};

   $sim += $weight*$weight2;
   $base1 += $weight*(1-$weight2);
   $base2 += $weight2*(1-$weight);
  }
  else
  {
   $base1 += $weight;
  }
 }

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id_b}{NEIGHBOURS}})
 {
  $key = $engine->GetInteractionKey($id_b, $neighbour);
  $weight = $engine->{INTERACTIONS}{$key}{RELIABILITY};

  if ($neighbour eq $id_a)
  {
   $base2 += (1 - $weight);
   $sim += $weight;
  }
  elsif (!defined($engine->{ENTITIES}{$id_a}{NEIGHBOURS}{$neighbour}))
  {
   $base2 += $weight;
  }
 }

 if (!defined($engine->{ENTITIES}{$id_a}{NEIGHBOURS}{$id_b}))
 {
  $base1 ++;
  $base2 ++;
 }

# $base1 += $sim;
# $base2 += $sim;
# $base1 = $engine->{'g_lambda'}>$base1?$engine->{'g_lambda'}:$base1;
# $base2 = $engine->{'g_lambda'}>$base2?$engine->{'g_lambda'}:$base2;
# $base1 -= $sim;
# $base2 -= $sim;

 $base1 += $sim*2;
 $base2 += $sim*2;

 if (!$sim)
 {
  return 0;
 }
 elsif (!$base1 || !$base2)
 {
  die($sim." ".$id_a." ".$id_b);
 }

 return ($sim*2/$base1) * ($sim*2/$base2);
}

# pathratio subroutines

sub _computeMaxPathStrength()
{
 my($self) = $_[0];
 my($id_a) = $_[1];
 my($id_b) = $_[2];
 my($k) = $_[3];

 my($engine) = $self->{'engine'};

 #compute Max Path Strength
 if ($k<2)
 {
  return 0;
 }
 elsif ($k==2)
 {
  return sqrt(scalar(keys %{$engine->{ENTITIES}{$id_a}{NEIGHBOURS}})
              *scalar(keys %{$engine->{ENTITIES}{$id_b}{NEIGHBOURS}}));
 }
 elsif ($k==3)
 {
  return scalar(keys %{$engine->{ENTITIES}{$id_a}{NEIGHBOURS}})
              *scalar(keys %{$engine->{ENTITIES}{$id_b}{NEIGHBOURS}});
 }
 else
 {
  my $neighbourA;
  my $neighbourB;
  my $maxpathstrength = 0;

  foreach $neighbourA(keys %{$engine->{ENTITIES}{$id_a}{NEIGHBOURS}})
  {
   foreach $neighbourB(keys %{$engine->{ENTITIES}{$id_b}{NEIGHBOURS}})
   {
    $maxpathstrength += $self->_computeMaxPathStrength($neighbourA, $neighbourB, $k-2);
   }
  }
  return $maxpathstrength;
 }
}

sub _recurse_computePathRatio()
{
 my($self) = $_[0];
 my($id_a) = $_[1];
 my($id_b) = $_[2];
 my($pathstrengths) = $_[3];
 my($visited) = $_[4];
 my($targets) = $_[5];
 my($limit) = $_[6];
 my($length) = $_[7];
 my($strength) = $_[8];
 my($engine) = $self->{'engine'};
 my $neighbour;

 $visited->{$id_a} = 1;

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id_a}{NEIGHBOURS}})
 {

  my $newtargets = $targets;
  my $key = $engine->GetInteractionKey($id_a, $neighbour);
  if (!defined($engine->{INTERACTIONS}{$key}{RELIABILITY}))
  {
   die($key);
  }
  my $weight = $strength*$engine->{INTERACTIONS}{$key}{RELIABILITY};

  if (!defined($visited->{$neighbour}))
  {
    if (defined($engine->{ENTITIES}{$id_b}{NEIGHBOURS}{$neighbour}))
    {
      $pathstrengths->{$neighbour}[$length] += $weight;
#      $newtargets --;
    }

    if ($length < $limit && $newtargets)
    {
     $self->_recurse_computePathRatio($neighbour, $id_b, $pathstrengths, $visited, $newtargets, $limit, $length+1, $weight);
    }
  }
 }

 delete $visited->{$id_a};
}

sub ComputePathRatio()
{
 my($self) = $_[0];
 my($id) = $_[1];
 my($engine) = $self->{'engine'};

 my $i;
 my %pathstrengths;
 my $neighbour;
 my %visited = ();
 my $limit = 5;

 $self->_recurse_computePathRatio($id, $id, \%pathstrengths, \%visited, scalar(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}}), $limit, 1, 1);
# print "done\n";

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})
 {
  my $score = 0;
  for ($i=2; $i<6; $i++)
  {
   $score += $pathstrengths{$neighbour}[$i] / $self->_computeMaxPathStrength($id, $neighbour, $i);
  }
  $engine->{WEIGHTS}{$id}{$neighbour} = $engine->{WEIGHTS}{$neighbour}{$id} = $score;
 }
}

sub _recurse_computePathRatioNew()
{
 my($self) = $_[0];
 my($id_a) = $_[1];
 my($pathstrengths) = $_[2];
 my($visited) = $_[3];
 my($limit) = $_[4];
 my($length) = $_[5];
 my($strength) = $_[6];
 my($engine) = $self->{'engine'};
 my $neighbour;

 $visited->{$id_a} = 1;

 foreach $neighbour(keys %{$engine->{ENTITIES}{$id_a}{NEIGHBOURS}})
 {
  my $key = $engine->GetInteractionKey($id_a, $neighbour);
  if (!defined($engine->{INTERACTIONS}{$key}{RELIABILITY}))
  {
   die($key);
  }
  my $weight = $strength*$engine->{INTERACTIONS}{$key}{RELIABILITY};

  if (!defined($visited->{$neighbour}))
  {
    $pathstrengths->{$neighbour}[$length] += $weight;

    if ($length < $limit)
    {
     $self->_recurse_computePathRatioNew($neighbour, $pathstrengths, $visited, $limit, $length+1, $weight);
    }
  }
 }

 delete $visited->{$id_a};
}

sub ComputePathRatioNew()
{
 my($self) = $_[0];
 my($id_a) = $_[1];
 my($engine) = $self->{'engine'};

 my $i;
 my $neighbour;
 my %pathstrengths;
 my %visited = ();
 my $limit = 3;

 %visited = ();
 $self->_recurse_computePathRatioNew($id_a, \%pathstrengths, \%visited, $limit, 1, 1);

 foreach $neighbour(keys %pathstrengths)
 {
  if ($neighbour eq $id_a)
  {
   next;
  }
  my $key = $engine->GetInteractionKey($id_a, $neighbour);
  $engine->{WEIGHTS}{$id_a}{$neighbour} = $engine->{WEIGHTS}{$neighbour}{$id_a} = 0;
  for ($i=2; $i<6; $i++)
  {
   my $maxpath = $self->_computeMaxPathStrength($id_a, $neighbour, $i) + 2;
   if (!$maxpath)
   {
    last;
   }

   if (defined($engine->{ENTITIES}{$id_a}{NEIGHBOURS}{neighbour}))
   {
    $pathstrengths{$neighbour}[$i] += 2;
   }

   if ($pathstrengths{$neighbour}[$i] > 0)
   {
     $engine->{WEIGHTS}{$id_a}{$neighbour} = $engine->{WEIGHTS}{$neighbour}{$id_a} = $pathstrengths{$neighbour}[$i] / $maxpath;
     last;
   }
  }
 }

 print "$id_a\n";
}

return 1;
