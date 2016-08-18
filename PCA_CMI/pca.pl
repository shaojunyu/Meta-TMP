open IN,"taxa.phylum.Abd";

%OTU;
$line_num = 1;
while (<IN>) {
	chomp;
	@a = split(/\t/, $_);
	print scalar @a."\n";
	# print $_;
}