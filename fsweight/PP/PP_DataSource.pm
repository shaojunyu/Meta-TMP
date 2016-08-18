use strict; 
package PP_DataSource; 
use POSIX; 

sub new() {
 ###########################################################
 # new(
 # pp_engine
 # )
 ###########################################################

 my($type) = $_[0];
 my($self) = {};

 $self->{'engine'} = $_[1];

 bless($self, $type);
 return($self); }

sub ProcessDataSource() {
 my($self) = $_[0];
 my $datatype = $_[1];
 my $minprecision = $_[2];

 my $engine = $self->{'engine'};

 my $subtype;
 my $range;
 my $key;
 my $annotation;

 foreach $subtype(keys %{$engine->{DATATYPE}{$datatype}})
 {
   my $link;
   my %indexes;
   my %links;
   my $maxvalue = -1;
   my $minvalue = -1;
   
   # determine range based on precision
   foreach $key(keys %{$engine->{INTERACTIONS}})
   {
    my $id_a;
    my $id_b;
    ($id_a, $id_b) = split(/\|/, $key);
    
    if (!scalar(keys %{$engine->{ENTITIES}{$id_a}{ANNOTATIONS}}) || !scalar(keys %{$engine->{ENTITIES}{$id_b}{ANNOTATIONS}})
    )
    {
     next;
    }

    if (defined($engine->{INTERACTIONS}{$key}{DATASOURCE}{$datatype}{$subtype}))
    {
     my %functions;
     foreach $annotation(keys %{$engine->{ENTITIES}{$id_a}{ANNOTATIONS}})
     {
      $functions{$annotation} = 0;
     }
     foreach $annotation(keys %{$engine->{ENTITIES}{$id_b}{ANNOTATIONS}})
     {
      if (defined($functions{$annotation}))
      {
       $functions{$annotation} = 1;
      }
      else
      {
       $functions{$annotation} = 0;
      }
     }

     my $score;
     delete $engine->{INTERACTIONS}{$key}{DATASOURCE}{$datatype}{$subtype}{RANGE};
     foreach $score(sort {$b <=> $a} @{$engine->{INTERACTIONS}{$key}{DATASOURCE}{$datatype}{$subtype}{'scores'}})
     {
       if ($maxvalue == -1 || $maxvalue < $score) { $maxvalue = $score; }
       if ($minvalue == -1 || $minvalue > $score) { $minvalue = $score; }
       foreach $annotation(keys %functions)
       {
        if ($engine->{ANNOTATIONS}{$annotation}{NAME} =~ m/\|NIF_/) { next; }
        push @{$links{$annotation}[$indexes{$annotation}]}, $score;
        push @{$links{$annotation}[$indexes{$annotation}]}, $functions{$annotation};
        $indexes{$annotation} ++;
       } 
       # last;
     }
    }
   }

   my $score_delta = ($maxvalue-$minvalue)/$engine->{NUMBINS};
   $engine->{DATATYPE}{$datatype}{$subtype}{DELTA} = $score_delta;
   $engine->{DATATYPE}{$datatype}{$subtype}{MINVAL} = $minvalue;
   foreach $annotation(keys %links)
   {
    my $total = 0;
    my $corrects = 0;
    my %bins;
    foreach $link(sort {$$a[0] <=> $$b[0]} @{$links{$annotation}})
    {
     my $bin = 1;
     if ($score_delta > 0) { $bin = floor(($$link[0]-$minvalue)/$score_delta)+1; }
     if ($bin > $engine->{NUMBINS}) { $bin = $engine->{NUMBINS}; }
#     print "$$link[0] / $score_delta - $bin\n";
     $bins{$bin}[0] ++;
     $bins{$bin}[1] += $$link[1];
     $total ++;
     $corrects += $$link[1];
    }
#    my $globalprecision = $corrects/($total+1);
    foreach (keys %bins)
    {
     my $precision = ($bins{$_}[1]+0.01)/($bins{$_}[0]+1);
#     if ($bins{$_}[0] < 10) { $precision = $globalprecision; }
     if ($precision < $minprecision) { next; }
     my $bound = 0;
     if ($score_delta > 0) { $bound = $score_delta*($_-1)+$minvalue; }
#     print "$datatype -> $subtype -> $annotation -> $_ -> $bound ($bins{$_}[1]/$bins{$_}[0]) = $precision\n";
     $engine->{DATATYPE}{$datatype}{$subtype}{RANGE}{$annotation}{$_} = $precision;
    }
  } # end each annot
 } # end each subtype 
} # end processData()

return 1;
