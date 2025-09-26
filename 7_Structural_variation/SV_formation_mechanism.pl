#!/bin/perl
# meerkat.pl
# Identify SVs by read pairs and split reads, give precise break points by local alignment
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
# Email: lixing_yang@hms.harvard.edu
# System requirements:
# samtools 0.1.5 or above
# BWA 0.5.7 or above
# NCBI blast 2.2.10 or above
# i.e. perl scripts/meerkat.pl -F /db/hg18/hg18_fasta/ -W /opt/bwa/ -B /opt/blast/bin/ -b

use strict;
use Getopt::Std;
use FindBin '$Bin';
my $version = 'v.0.174';

my %opts = (o=>1, t=>100000, z=>1000000000);
getopts("o:t:z:R:I:D:h", \%opts);

my $intra_refine_type = $opts{I};
my $outfile = $opts{D};
my $include_other = $opts{o};
my $rmskfile = $opts{R};
my $te_size_max = $opts{t};
my $sv_size_cutoff = $opts{z};
my $del_ins_size_cutoff_d = 0.8; # size ratio of del and ins in del_ins events
my $del_ins_size_cutoff_u = 1.2;
my $ovl = 0.8; # overlap of a predicted events and a annotated TE


foreach my $key (keys %opts) {
    print "Option $key = $opts{$key}\n";
}


my $intra_refine_type = $opts{I} // die "Input file not specified\n";
my $outfile = $opts{D} // die "Output file not specified\n";


print "Input file: $intra_refine_type\n";
print "Output file: $outfile\n";

my (%te, %tei, %sr, %sri);
my $newline;
open RMSK, "<$rmskfile";
while ($newline = <RMSK>)
{
	chomp $newline;
	my @data = split (/\t/, $newline);
	if ($data[3] eq 'LINE' or $data[3] eq 'SINE' or $data[3] eq 'LTR' or $data[3] eq 'DNA')
	{
		$tei{$data[0]} = 0 unless ($tei{$data[0]});
		$te{$data[0]}[$tei{$data[0]}][0] = $data[1];# start
		$te{$data[0]}[$tei{$data[0]}][1] = $data[2];# end
		$te{$data[0]}[$tei{$data[0]}][2] = $data[3];# class
		$te{$data[0]}[$tei{$data[0]}][3] = $data[3];# name
		$tei{$data[0]}++;
	}
	if ($include_other and $data[3] eq 'Other')
	{
		$tei{$data[0]} = 0 unless ($tei{$data[0]});
		$te{$data[0]}[$tei{$data[0]}][0] = $data[1];# start
		$te{$data[0]}[$tei{$data[0]}][1] = $data[2];# end
		$te{$data[0]}[$tei{$data[0]}][2] = $data[3];# class
		$te{$data[0]}[$tei{$data[0]}][3] = $data[3];# name
		$tei{$data[0]}++;
	}
	if ($data[3] eq 'Satellite' or $data[3] eq 'Simple_repeat' or $data[3] eq 'Low_complexity')
	{
		$sri{$data[0]} = 0 unless ($sri{$data[0]});
		$sr{$data[0]}[$sri{$data[0]}][0] = $data[1];# start
		$sr{$data[0]}[$sri{$data[0]}][1] = $data[2];# end
		$sr{$data[0]}[$sri{$data[0]}][2] = $data[3];# class
		$sri{$data[0]}++;
	}
}
close RMSK;

my @variant;
my $i = 0;
open FILE, "<$intra_refine_type";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split (/\t/, $newline);
	$variant[$i] = \@data;
	#print "@{$variant[$i]}\n";
	$i++;
}
close FILE;

open OUT, ">$outfile";
LOOP: foreach (@variant)
{
	my @data = @$_;
	my $type = $data[0];
	my ($mechanism, @bp_annotation);
	# print "@data\n";
	# call mechanism
	if ($type eq 'DEL' or $type eq 'DUP')
	{
		# large events
		if ($data[5] > $sv_size_cutoff)
		{
			;
		}
		else
		{
			my ($te_class, $te_name) = &tei($data[2], $data[3], $data[4]);
			$mechanism = 'TEI_'.$te_class.'_'.$te_name if ($te_class);
			$mechanism = 'TEI_complex' if ($te_class eq 'complex');
			goto LEND if ($mechanism);
			my $vntr = &vntr($data[2], $data[3], $data[4]);
			$mechanism = 'VNTR' if ($vntr);
			goto LEND if ($mechanism);
			$mechanism = 'NAHR' if ($data[6] > 100);
			goto LEND if ($mechanism);
			$mechanism = 'alt-EJ' if ($data[6] >= 2 and $data[6] <= 100);
			goto LEND if ($mechanism);
			$mechanism = 'NHEJ';
			goto LEND if ($mechanism);
		}
	}
LEND:	$mechanism = 'NA' unless ($mechanism);
	shift(@data);
	unshift(@data, $mechanism);
	unshift(@data, $type);
	
	# annotate break points
	if ($type eq 'DEL' or $type eq 'DUP')
	{
		if (&sr_overlap($data[3], $data[4]))
		{
			$bp_annotation[0] = 'SR';
		}
		elsif (&te_overlap($data[3], $data[4]))
		{
			$bp_annotation[0] = 'TE';
		}
		else
		{
			$bp_annotation[0] = '';
		}
		if (&sr_overlap($data[3], $data[5]))
		{
			$bp_annotation[1] = 'SR';
		}
		elsif (&te_overlap($data[3], $data[5]))
		{
			$bp_annotation[1] = 'TE';
		}
		else
		{
			$bp_annotation[1] = '';
		}
	}
	my $bp_annotation = join ("_", @bp_annotation);
	$bp_annotation = 'BP:'.$bp_annotation;
	push @data, $bp_annotation;
	my $toprint = join ("\t", @data);
	print OUT "$toprint\n";
}
close OUT;

# TE insertion only
sub tei
{
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	return 0 if (abs($end-$start) > $te_size_max);
	my @overlap;
	foreach (@{$te{$chr}})
	{
		my $te = $_;
		next if ($$te[0] < $end and $$te[1] < $end and $$te[0] < $start and $$te[1] < $start);
		last if ($$te[0] > $end and $$te[1] > $end and $$te[0] > $start and $$te[1] > $start);
		if (&covered($start, $end, $$te[0], $$te[1]))
		{
			my ($overlap1, $overlap2) = &overlap($start, $end, $$te[0], $$te[1]);
			if (abs($overlap2-$overlap1)/(abs($$te[1]-$$te[0])+0.5)>=$ovl)
			{
				if (abs($overlap2-$overlap1)/(abs($end-$start)+0.5)>=$ovl)
				{
					return ($$te[2], $$te[3]);
				}
				else
				{
					push @overlap, $overlap2-$overlap1;
				}
			}
		}
	}
	if ($overlap[0])
	{
		my $overlap;
		foreach (@overlap)
		{
			$overlap += $_;
		}
		if ($overlap/(abs($end-$start)+0.5)>=$ovl)
		{
			return ('complex', 1);
		}
	}
	return 0;	
}


sub vntr
{
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	return 0 if (abs($end-$start) > $te_size_max);
	my @overlap;
	foreach (@{$sr{$chr}})
	{
		my $sr = $_;
		next if ($$sr[0] < $end and $$sr[1] < $end and $$sr[0] < $start and $$sr[1] < $start);
		last if ($$sr[0] > $end and $$sr[1] > $end and $$sr[0] > $start and $$sr[1] > $start);
		if (&covered($start, $end, $$sr[0], $$sr[1]))
		{
			my ($overlap1, $overlap2) = &overlap($start, $end, $$sr[0], $$sr[1]);
			if (abs($overlap2-$overlap1)/(abs($end-$start)+0.5)>=$ovl)
			{
				return 1;
			}
		}
	}
	return 0;	
}

# overlap satellite repeat, simple repeat, low complexity repeat
sub sr_overlap
{
	my $chr = shift;
	my $bp = shift;
	foreach (@{$sr{$chr}})
	{
		next if ($$_[0] < $bp and $$_[1] < $bp);
		last if ($$_[0] > $bp and $$_[1] > $bp);
		if ($bp >= $$_[0] and $bp <= $$_[1])
		{
			return 1;
		}
	}
	return 0;
}

# overlap TE
sub te_overlap
{
	my $chr = shift;
	my $bp = shift;
	foreach (@{$te{$chr}})
	{
		next if ($$_[0] < $bp and $$_[1] < $bp);
		last if ($$_[0] > $bp and $$_[1] > $bp);
		if ($bp >= $$_[0] and $bp <= $$_[1])
		{
			return 1;
		}
	}
	return 0;	
}


sub covered
{
	my $a1 = shift;
	my $a2 = shift;
	my $b1 = shift;
	my $b2 = shift;
	return 0 if ($a1 =~ /\D/ or $a2 =~ /\D/ or $b1 =~ /\D/ or $b2 =~ /\D/);
	($a1, $a2) = sort{ $a <=> $b } ($a1, $a2);
	($b1, $b2) = sort { $a <=> $b } ($b1, $b2);
	if (($a1 >= $b1 and $a1 <= $b2) or ($a2 >= $b1 and $a2 <= $b2) or ($b1 >= $a1 and $b1 <= $a2) or ($b2 >= $a1 and $b2 <= $a2))
	{
		return 1;
	}
}


sub overlap
{
	my $a1 = shift;
	my $a2 = shift;
	my $b1 = shift;
	my $b2 = shift;
	return ($a1, $a2) unless ($b1 and $b2);
	($a1, $a2) = sort { $a <=> $b } ($a1, $a2);
	($b1, $b2) = sort { $a <=> $b } ($b1, $b2);
	if (($a1 >= $b1 and $a1 <= $b2) or ($a2 >= $b1 and $a2 <= $b2) or ($b1 >= $a1 and $b1 <= $a2) or ($b2 >= $a1 and $b2 <= $a2))
	{
		my $a = ($a1>$b1)?$a1:$b1;
		my $b = ($a2>$b2)?$b2:$a2;
		return ($a, $b);
	}
	else
	{
		return (0, 0);
	}
}
