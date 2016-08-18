package PP_Validation; 
use strict;

sub new() 
{
 ###########################################################
 # new(
 # pp_engine
 # )
 # validation methods
 ###########################################################

 my($type) = $_[0];
 my($self) = {};

 $self->{'engine'} = $_[1];

 bless($self, $type);
 return($self); 
}

sub ComputeTermROC() 
{
 ###########################################################
 # Compute_ROC(
 # )
 # compute roc of preditions for each annotation
 ###########################################################
	
 my($self, $annotation, $fold, $min) = @_;
 my $engine = $self->{'engine'};
 my $id;

 my $annotated = 0;
 my @predictions;
 my $index = 0;
 foreach $id(keys %{$engine->{ENTITIES}})
 {
  if (($fold != -1 && $engine->{ENTITIES}{$id}{FOLDLABEL} != $fold) || 
       !defined($engine->{ENTITIES}{$id}{ANNOTATED}) || 
       ($engine->{IGNOREISOLATED} && !scalar(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})) || 
       !defined($engine->{ENTITIES}{$id}{TEST}))
  { next; }
  
  if (defined($engine->{ENTITIES}{$id}{ANNOTATIONS}{$annotation})) 
  { $annotated++; }

  my $score = defined($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation})?$engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation}:0;
  my $label = (defined($engine->{ENTITIES}{$id}{ANNOTATIONS}{$annotation}) && $engine->{ENTITIES}{$id}{ANNOTATIONS}{$annotation} == 1)?1:0; 
  $predictions[$index][0] = $score;
  $predictions[$index++][1] = $label;
 }
 if ($annotated < $min) { return -1; }

 my $roc = 0;
 my $tp = 0;
 my $fp = 0;
 my $lasttp = 0;
 my $numfp = 0;
 my $lastscore = -1;
 foreach (sort {$b->[0] <=> $a->[0]} @predictions)
 {
   my $score = $_->[0];
   my $label = $_->[1];

   if ($score != $lastscore && $lastscore != 1)
   {
      $roc += ($lasttp + ($tp-$lasttp)/2)*$numfp;
      $fp += $numfp;
      $numfp = 0;
      $lasttp = $tp;
   }
   $lastscore = $score;
  
   if ($label == 1)
   {
     $tp ++;
   }
   else
   {
     $numfp ++;
   }
 }
 $roc += ($lasttp + ($tp-$lasttp)/2)*$numfp;   
 $fp += $numfp;
 if ($roc > 0)
 {
  $roc /= ($fp*$tp);
 }
 return $roc;
}

sub ComputePredictedROC() 
{
 ###########################################################
 # Compute_ROC(
 # )
 # compute roc of preditions for each annotation
 ###########################################################
	
 my($self, $annotation, $fold, $min) = @_;
 my $engine = $self->{'engine'};
 my $id;

 my $annotated = 0;
 my @predictions;
 my $index = 0;
 foreach $id(keys %{$engine->{ENTITIES}})
 {
  if (($fold != -1 && $engine->{ENTITIES}{$id}{FOLDLABEL} != $fold) || 
       !defined($engine->{ENTITIES}{$id}{ANNOTATED}) || 
       ($engine->{IGNOREISOLATED} && !scalar(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}})) || 
       !defined($engine->{ENTITIES}{$id}{TEST}))
  { next; }

  if (!defined($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation})) { next; }
  
  if (defined($engine->{ENTITIES}{$id}{ANNOTATIONS}{$annotation})) 
  { $annotated++; }

  my $score = defined($engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation})?$engine->{ENTITIES}{$id}{PREDICTIONS}{$annotation}:0;
  my $label = (defined($engine->{ENTITIES}{$id}{ANNOTATIONS}{$annotation}) && $engine->{ENTITIES}{$id}{ANNOTATIONS}{$annotation} == 1)?1:0; 
  $predictions[$index][0] = $score;
  $predictions[$index++][1] = $label;
 }
 if ($annotated < $min) { return -1; }

 my $roc = 0;
 my $tp = 0;
 my $fp = 0;
 my $lasttp = 0;
 my $numfp = 0;
 my $lastscore = -1;
 foreach (sort {$b->[0] <=> $a->[0]} @predictions)
 {
   my $score = $_->[0];
   my $label = $_->[1];

   if ($score != $lastscore && $lastscore != 1)
   {
      $roc += ($lasttp + ($tp-$lasttp)/2)*$numfp;
      $fp += $numfp;
      $numfp = 0;
      $lasttp = $tp;
   }
   $lastscore = $score;
  
   if ($label == 1)
   {
     $tp ++;
   }
   else
   {
     $numfp ++;
   }
 }
 $roc += ($lasttp + ($tp-$lasttp)/2)*$numfp;   
 $fp += $numfp;
 if ($roc > 0)
 {
  $roc /= ($fp*$tp);
 }
 return $roc;
}


sub Compute_ROC() 
{
 ###########################################################
 # Compute_ROC(
 # )
 # compute roc of preditions for each annotation
 ###########################################################
	
 my($self, $fold) = @_;
 my $engine = $self->{'engine'};
 my $annotation;
 my %rocs;

 foreach $annotation(sort {$a <=> $b} keys %{$engine->{ANNOTATIONS}})
 {
  if (!defined($engine->{ANNOTATIONS}{$annotation}{INF})) { next; }
#  if ($engine->{ANNOTATIONS}{$annotation}{LEVEL} != 4) { next; }

  my $roc = $self->ComputeTermROC($annotation, $fold, 1);
  if ($roc < 0) { next; }
  $rocs{$annotation} = $roc;
  print $engine->{ANNOTATIONS}{$annotation}{NAME}."\t$annotation\t $roc\n"; 
 }

 my $i = 0;
 foreach ($i=0; $i<=1; $i+=0.1)
 {
  my $count = 0;
  foreach (values %rocs)
  {
   if ($_ >= $i)
   {
    $count++;
   }
  }
  print "$i\t$count\n";
 }
}

sub Compute_PrecisionVSRecall() 
{
 ###########################################################
 # Compute_PrecisionVSRecall(
 # )
 # compute Precision VS Recall of preditions for entities
 ###########################################################

 my($self) = $_[0];
 my $engine = $self->{'engine'};

 # Precision vs Recall
 my $id;
 my @predictions = ();
 my $index = 0;
 my $annotation_count = 0;
 my $total = 0;
 my $annotated = 0;

 foreach $id(keys %{$engine->{ENTITIES}})
 {
  if (!defined($engine->{ENTITIES}{$id}{TEST}) || 
      !defined($engine->{ENTITIES}{$id}{ANNOTATED}) 
||     ($engine->{IGNOREISOLATED} && !scalar(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}}))
  )
  {
   next;
  }
  $annotated++;

  foreach (keys %{$engine->{ENTITIES}{$id}{ANNOTATIONS}})
  {
   if (!defined($engine->{ENTITIES}{$id}{ANNOTATIONS}{$_}) 
    || $engine->{ENTITIES}{$id}{ANNOTATIONS}{$_} != 1 
    || !defined($engine->{ANNOTATIONS}{$_}{INF})
   )
   {
    next;
   }
   $annotation_count++;
  }
  foreach (keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
  {
   if (!defined($engine->{ANNOTATIONS}{$_}{INF}))
   {
    next;
   }
#   if ($engine->{ENTITIES}{$id}{LEVEL} >= $engine->{ANNOTATIONS}{$_}{LEVEL})
   {
    push @{$predictions[$index]}, $engine->{ENTITIES}{$id}{PREDICTIONS}{$_};
    push @{$predictions[$index++]}, (defined($engine->{ENTITIES}{$id}{ANNOTATIONS}{$_}) && $engine->{ENTITIES}{$id}{ANNOTATIONS}{$_} == 1)?1:0;
   }
  }
 }
 print "$annotation_count\t".scalar(@predictions)."\t$annotated\n";

 my $predicts = 0;
 my $corrects= 0;
 my $lastvalue = -1;
 my $threshold = 0.02;
 my $recall;

 foreach (sort {$$b[0] <=> $$a[0]} @predictions)
 {
# print "$predicts\t$$_[0]\t$$_[1]\t".$predicts/$annotation_count."\n";
  if ($lastvalue != $$_[0])
  {
   $recall = $corrects/$annotation_count;
   if ($$_[0] != $lastvalue && $lastvalue != -1 && $recall > $threshold)
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
 print "$lastvalue\t$predicts\t$recall\t".$corrects/$predicts."\n";
}

sub Compute_MultiLevelPrecisionVSRecall() 
{
 ###########################################################
 # Compute_PrecisionVSRecall(
 # )
 # compute Precision VS Recall of preditions for entities
 ###########################################################

 my($self) = $_[0];
 my($outfile) = $_[1];
 my $engine = $self->{'engine'};

 my $maxlevel = 0;
 foreach (keys %{$engine->{ANNOTATIONS}})
 {
  if ($engine->{ANNOTATIONS}{$_}{LEVEL} > $maxlevel)
  {
   $maxlevel = $engine->{ANNOTATIONS}{$_}{LEVEL};
  }
 }

 open DAT, ">$outfile";
 my $i=0;
 for ($i=1; $i<=$maxlevel; $i++)
 {
  print DAT "###### LEVEL $i ######\n";
  # Precision vs Recall
  my $id;
  my @predictions = ();
  my $index = 0;
  my $annotation_count = 0;
  my $total = 0;
  my $annotated = 0;

  foreach $id(keys %{$engine->{ENTITIES}})
  {
   if (!defined($engine->{ENTITIES}{$id}{TEST}) || 
       !defined($engine->{ENTITIES}{$id}{ANNOTATED}) || 
       $engine->{ENTITIES}{$id}{LEVEL} < $i || 
       ($engine->{IGNOREISOLATED} && !scalar(keys %{$engine->{ENTITIES}{$id}{NEIGHBOURS}}))
   )
   {
    next;
   }
   $annotated++;

   foreach (keys %{$engine->{ENTITIES}{$id}{ANNOTATIONS}})
   {
    if ($engine->{ENTITIES}{$id}{ANNOTATIONS}{$_} != 1 || $engine->{ANNOTATIONS}{$_}{LEVEL} != $i)
    {
     next;
    }
    $annotation_count++;
   }
   foreach (keys %{$engine->{ENTITIES}{$id}{PREDICTIONS}})
   {
    if ($engine->{ANNOTATIONS}{$_}{LEVEL} != $i)
    {
     next;
    }
#    if ($engine->{ENTITIES}{$id}{LEVEL} >= $engine->{ANNOTATIONS}{$_}{LEVEL})
    {
     push @{$predictions[$index]}, $engine->{ENTITIES}{$id}{PREDICTIONS}{$_};
     push @{$predictions[$index++]}, (defined($engine->{ENTITIES}{$id}{ANNOTATIONS}{$_}) && $engine->{ENTITIES}{$id}{ANNOTATIONS}{$_} == 1)?1:0;
    }
   }
  }
  print DAT "$annotation_count\t".scalar(@predictions)."\t$annotated\n";

  my $predicts = 0;
  my $corrects= 0;
  my $lastvalue = -1;
  my $threshold = 0.02;
  my $recall;

  foreach (sort {$$b[0] <=> $$a[0]} @predictions)
  {
#   print "$predicts\t$$_[0]\t$$_[1]\t".$predicts/$annotation_count."\n";
   if ($lastvalue != $$_[0])
   {
    $recall = $annotation_count>0?($corrects/$annotation_count):0;
    if ($$_[0] != $lastvalue && $lastvalue != -1 && $recall > $threshold)
    {
     print DAT "$lastvalue\t$predicts\t$recall\t".$corrects/$predicts."\n";
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
  print DAT "$lastvalue\t$predicts\t$recall\t".($predicts>0?($corrects/$predicts):0)."\n";
 } 
 close DAT;
}

return 1;
