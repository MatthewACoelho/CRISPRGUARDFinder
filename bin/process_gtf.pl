#!/usr/bin/env perl

use strict;

my $type_codes =
{
	exon				=> 0x1,
	CDS					=> 0x2,
	five_prime_utr		=> 0x4,
	three_prime_utr		=> 0x8,
	intron				=> 0x10,
	five_prime_flank	=> 0x20,
	three_prime_flank	=> 0x40,
	enhancer			=> 0x80,
	promoter			=> 0x100,
	tss					=> 0x200
};

my $flank_size = 1000;

my $out_dir = $ARGV[0] || ".";
my $last_chr = '';
my $genes = {};
my $gene_names = {};
my $exons = {};
my $cds = {};
my $tags = {};
my $tx = {};
my $features = {};
while (<STDIN>)
{
	chomp;
	my ($chr, $source, $type, $start, $end, $x, $strand, $y, $data) = split /\t/;

	next if not $type;

	if ($chr ne $last_chr)
	{
		if ($last_chr)
		{
			output($last_chr, $genes, $features);
		}
		$last_chr = $chr;
		$genes = {};
		$tx = {};
		$tags = {};
		$exons = {};
		$cds = {};
		$features = {};
	}

	if ($type eq 'gene')
	{
		my $gene_id = get('gene_id', $data);
		my $gene_name = get('gene_name', $data);

		$gene_names->{$gene_name} = { id => $gene_id, pos => "$chr:$start-$end" };;

		my $left_type;
		my $right_type;
		if ($strand eq '+')
		{
			$left_type = 'five_prime_flank';
			$right_type = 'three_prime_flank';
		}
		else
		{
			$left_type = 'three_prime_flank';
			$right_type = 'five_prime_flank';
		}
		my $left_start = ($start <= $flank_size ? 1 : ($start-$flank_size));
		my $left_end = $start-1;
		$features->{$left_start}->{$left_end}->{$gene_id}->{type}->{$left_type}++;
		$features->{$left_start}->{$left_end}->{$gene_id}->{type_code} |= $type_codes->{$left_type};

		my $right_start = $end+1;
		my $right_end = $end+$flank_size;
		$features->{$right_start}->{$right_end}->{$gene_id}->{type}->{$right_type}++;
		$features->{$right_start}->{$right_end}->{$gene_id}->{type_code} |= $type_codes->{$right_type};
	}
	elsif ($type eq 'transcript')
	{
	}
	elsif ($type eq 'exon' or $type eq 'CDS' or $type =~ /utr/)
	{
		my $gene_id = get('gene_id', $data);
		my $tx_id = get('transcript_id', $data);
		my $exon_number = get('exon_number', $data);
		my $gene_name = get('gene_name', $data);
		my $gene_biotype = get('gene_biotype', $data);
		if (exists $genes->{$gene_id})
		{
			$genes->{$gene_id}->{start} = $start if $start < $genes->{$gene_id}->{start};
			$genes->{$gene_id}->{end} = $end if $genes->{$gene_id}->{end} < $end;
		}
		else
		{
			$genes->{$gene_id} = { name => $gene_name || $gene_id, biotype => $gene_biotype, strand => $strand, start => $start, end => $end };
		}
		$tx->{$gene_id}->{$tx_id}++;
		$features->{$start}->{$end}->{$gene_id}->{type}->{$type}++;
		$features->{$start}->{$end}->{$gene_id}->{type_code} |= $type_codes->{$type};
		$features->{$start}->{$end}->{$gene_id}->{tx}->{$tx_id}++;
		$tags->{$gene_id}->{$tx_id} = get_tags($tags->{$gene_id}->{$tx_id}, $data);

		if ($type eq 'exon')
		{
			push @{$exons->{$gene_id}->{$tx_id}}, { start => $start, end => $end };
		}
		elsif ($type eq 'CDS')
		{
			push @{$cds->{$gene_id}->{$tx_id}}, { start => $start, end => $end };
		}
	}
}
output($last_chr, $genes, $features);

open(OUT, ">$out_dir/gene_names.txt");
for my $gene_name (sort keys %$gene_names)
{
	print OUT "$gene_name\t$gene_names->{$gene_name}->{id}\t$gene_names->{$gene_name}->{pos}\n";
}
close(OUT);

sub get
{
	my ($name, $data) = @_;

	if ($data =~ /$name "([^"]*)"/)
	{
		return $1;
	}
	return "";
}

sub get_tags
{
	my ($tags, $data) = @_;

	while ($data =~/tag "([^"]*)"/g)
	{
		$tags->{$1}++;
	}
	return $tags;
}

sub output
{
	my ($chr, $genes, $features) = @_;

	for my $gene_id (keys %$exons)
	{
		for my $tx_id (keys %{$exons->{$gene_id}})
		{
			my @e = sort { $a->{start} <=> $b->{start} } @{$exons->{$gene_id}->{$tx_id}};
			for(my $i = 1; $i < scalar @e; $i++)
			{
				my $intron_start = $e[$i-1]->{end} + 1;
				my $intron_end = $e[$i]->{start} - 1;
				$features->{$intron_start}->{$intron_end}->{$gene_id}->{type}->{intron}++;
				$features->{$intron_start}->{$intron_end}->{$gene_id}->{type_code} |= $type_codes->{intron};
				$features->{$intron_start}->{$intron_end}->{$gene_id}->{tx}->{$tx_id}++;
			}
		}
	}

	open(CHR, ">$out_dir/chr$chr.info") || die;
	for my $start (sort {$a <=> $b} keys %{$features})
	{
		for my $end (sort {$a <=> $b} keys %{$features->{$start}})
		{
			my $gs = $features->{$start}->{$end};
			for my $gene_id (sort keys %$gs)
			{
				my $g = $genes->{$gene_id};
				my $ntx = scalar keys %{$tx->{$gene_id}};
				my $e = $gs->{$gene_id};
				my $etype = join(",", sort keys %{$e->{type}});
				my $etype_code = $e->{type_code};
				my $entx = scalar keys %{$e->{tx}};

				#print "$gene_id\t$etype\t$etype_code\t$chr\t$start\t$end\t$g->{strand}\t";
				#print "$entx/$ntx\t$g->{name}\t$g->{biotype}\n";

				print CHR "$gene_id\t$etype_code\t$start\t$end\t$g->{strand}\t$entx\n";
			}
		}
	}
	close(CHR);

	open(CHR, ">$out_dir/chr$chr.genes") || die;
	for my $gene_id (sort keys %$genes)
	{
		my $g = $genes->{$gene_id};
		my $ntx = scalar keys %{$tx->{$gene_id}};
		print CHR "$gene_id\t$chr\t$g->{start}\t$g->{end}\t$g->{strand}\t$ntx\t$g->{name}\t$g->{biotype}\n";
	}
	close(CHR);

	open(CDS, ">$out_dir/chr$chr.cds") || die;
	for my $gene_id (sort keys %$cds)
	{
		for my $tx_id (sort keys %{$cds->{$gene_id}})
		{
			my $tag_s = join(",", sort keys %{$tags->{$gene_id}->{$tx_id}});
			my $g = $genes->{$gene_id};
			print CDS "$gene_id\t$g->{name}\t$g->{strand}\t$tx_id\t";
			my @c = sort { $a->{start} <=> $b->{start} } @{$cds->{$gene_id}->{$tx_id}};
			my $sep = "";
			for my $feature (@c)
			{
				print CDS "$sep$feature->{start}:$feature->{end}";
				$sep = ",";
			}
			print CDS "\t$tag_s";
			print CDS "\n";
		}
	}
	close(CDS);
}

1;
