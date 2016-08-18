use strict;
package PP_Math;

my $factsize = 2;
my @factcache = ();
$factcache[0] = $factcache[1] = 0;

sub new()
{
 my($type) = $_[0];
 my($self) = {};

 bless($self, $type);
 return $self;
}

sub logfactorial()
{
 my ($self, $k) = @_;

 if ($k < $factsize)
 {
  return $factcache[$k];
 }

 while ($factsize <= $k)
 {
  $factcache[$factsize] =
  $factcache[$factsize-1] + log($factsize);
  $factsize ++;
 }
 return $factcache[$k];
}

sub Choose()
{
 my ($self, $N,$k) = @_;
 return $self->logfactorial($N) - ($self->logfactorial($k) + $self->logfactorial($N-$k));
}

return 1
