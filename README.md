# fastq_stats
Count the number of each base and quality at each base position across multiple fastq files

```
fastq_stats *.fastq.gz > fastq_stats.output
```

To summarise the output in to 5kb bins and 1 large bin for positions higher than 65kb:

```
perl -Mstrict -we 'open(my $fh, shift); my $s = 0;  my $p = 0; my @counts; print "bases (A,C,G,T %):\n"; while (<$fh>) { if ($_ =~ /^quals/) { last } if ($_ !~ /^\d/) { next } my @vals = split(",", $_); $p = shift @vals; for my $i (0 .. $#vals) { $counts[$i] += $vals[$i] } if ($s == 0) { $s = $p } if ($p < 65001 && ($p-4999 == $s)) { printbases($s, $p, \@counts); @counts = (); $s = 0; } } printbases($s, $p, \@counts); print "quals (avg Q):\n"; $s = 0; while (<$fh>) { if ($_ !~ /^\d/) { next } my @vals = split(",", $_); $p = shift @vals; for my $i (0 .. $#vals) { $counts[$i] += $vals[$i] } if ($s == 0) { $s = $p } if ($p < 65001 && ($p-4999 == $s)) { printquals($s, $p, \@counts); @counts = (); $s = 0; } } printquals($s, $p, \@counts); sub printbases { my ($s, $p, $counts) = @_; my $t = 0; foreach my $c (@$counts) { $t += $c } printf("%05d..%06d", $s, $p); foreach my $c (@$counts) { printf("\t%0.2f", ($c/$t)*100) } print "\n" } sub printquals { my ($s, $p, $counts) = @_; my $t = 0; my $n = 0; for my $i (0 .. $#counts) { my $c = $counts[$i]; $n += $c; $t += ($c*$i); } printf("%05d..%06d\t%0.2f\n", $s, $p, $t/$n) }' fastq_stats.output > fastq_stats.output.binned
```

To compare 2 different summary files:

```
perl -Mstrict -we 'open(my $fh1, shift); open(my $fh2, shift); while (my $l1 = <$fh1>) { my $l2 = <$fh2>; if ($l1 !~ /^\d/) { print "change in $l1"; next } chomp $l1; chomp $l2; compare($l1, $l2); } sub compare { my ($l1, $l2) = @_; my ($p, @v1) = split("\t", $l1); my @v2 = split("\t", $l2); shift @v2; my @diffs; for my $i (0 .. $#v1) { push(@diffs, sprintf("%+0.2f", $v2[$i]-$v1[$i])) } printf("%s\t%s\n", $p, join("\t", @diffs)) }' sample1.fastq_stats.output.binned sample2.fastq_stats.output.binned > diff
```
