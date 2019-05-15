#!/bin/env Rscript

suppressMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))

init_options =
function(options)
{
	if (options$genome == 'hg38')
	{
		suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
		options$genome_seq = BSgenome.Hsapiens.UCSC.hg38
	}
	else if (options$genome == 'mm10')
	{
		suppressPackageStartupMessages(library("BSgenome.Mmusculus.UCSC.mm10"))
		options$genome_seq = BSgenome.Mmusculus.UCSC.mm10
	}
	else if (file.exists(sprintf("%s/%s.fa", options$genome, options$genome)))
	{
		options$genome_seq = readDNAStringSet(sprintf("%s/%s.fa", options$genome, options$genome))
	}
	else
	{
		stop('unknown genome')
	}

	options$pam_len = nchar(options$pam)

	return (options)
}

pam_side =
function(pam)
{
	if (pam == "TTN" || pam == "TTTN")
	{
		return ('left')
	}
	return ('right')
}

pam_extend_grna =
function(pam, len, pam_seq, strand)
{
	#extend the PAM to find the gRNA seq
	side = pam_side(pam)
	if (side == 'right')
	{
		# extend to the left
		if (strand == '+')
		{
			start(pam_seq) = start(pam_seq)-len
		}
		else
		{
			end(pam_seq) = end(pam_seq)+len
		}
	}
	else
	{
		# extend to the right
		if (strand == '+')
		{
			end(pam_seq) = end(pam_seq)+len
		}
		else
		{
			start(pam_seq) = start(pam_seq)-len
		}
	}
	return (pam_seq)
}

calc_overlap =
function(region_a, region_b)
{
	if (region_a[1] > region_b[2])
	{
		return (0)
	}
	if (region_a[2] < region_b[1])
	{
		return (0)
	}
	return (min(region_a[2], region_b[2])-max(region_a[1], region_b[1])+1)
}

calc_gc_content =
function(seq)
{
	return (nchar(gsub("[^GC]", "", seq))/nchar(seq))
}

find_guards_in_seq =
function(options, sequence, guard_strand="+")
{
	pam = options$pam
	pam_seq = DNAString(pam)
	if (guard_strand == "-")
	{
		pam_seq = reverseComplement(pam_seq)
	}

	pam_data = matchPattern(pam_seq, sequence, fixed=FALSE)
	# Check that no PAMs includes Ns
	toremove  = grep("N", as.character(pam_data))
	if (length(toremove) > 0) { pam_data = pam_data[-toremove] }

	# extend the PAM to find the gRNA seq
	guard_view = pam_extend_grna(pam, options$guard_len, pam_data, guard_strand)

	# remove the hit that cannot be extended 
	toremove = which(start(guard_view) <= 0 | end(guard_view) > length(sequence))
	if (length(toremove) > 0) { guard_view = guard_view[-toremove] }

	if (length(guard_view) == 0)
	{
		return (NULL)
	}

	hits = NULL
	for(i in 1:length(guard_view))
	{
		guard_with_pam = as.character(guard_view[i])
		forward_guard_with_pam = guard_with_pam
		if (guard_strand == '-') forward_guard_with_pam = as.character(reverseComplement(guard_view[i]))

		forward_guard = toupper(forward_guard_with_pam)
		forward_pam = substring(forward_guard, options$guard_len+1)
		forward_guard = substring(forward_guard, 1, options$guard_len)

		guard_region = c(start(guard_view[i])+options$target_region[1]-1,
						 end(guard_view[i])+options$target_region[1]-1)

		row = c(ID=options$id,
				GuardChr=options$chr,
				GuardStart=guard_region[1],
				GuardStop=guard_region[2],
				GuardStrand=guard_strand,
				GuardWithPAM=guard_with_pam,
				GuardGC=calc_gc_content(forward_guard),
				ForwardGuardPAM=forward_pam,
				ForwardGuardWithPAM=forward_guard_with_pam,
				OffGuideOverlap=calc_overlap(options$region, guard_region),
				OffGuideSeedOverlap=calc_overlap(options$seed_region, guard_region),
				OffGuidePAMOverlap=calc_overlap(options$pam_region, guard_region),
				OffGuide=options$off_guide,
				OffGuideGC=options$off_guide_gc,
				Guard=forward_guard)
		hits = rbind(hits, row)
	}
	hits = as.data.frame(hits, stringsAsFactors=FALSE)
	return (hits)
}

find_guards_in_seqs =
function(options, seqs)
{
	out = NULL
	for(i in 1:length(seqs))
	{
		out = rbind(out, find_guards_in_seq(options, seqs[[i]], guard_strand="+"))
		out = rbind(out, find_guards_in_seq(options, seqs[[i]], guard_strand="-"))
	}
	return (out)
}

find_guard_seqs =
function(options, coord)
{
	options$id = coord$id
	options$chr = coord$chr
	options$strand = coord$strand
	options$region = c(coord$start, coord$end)
	options$target_region = c(coord$start-options$guard_len+1, coord$end+options$guard_len-1)

	if (coord$strand == '+')
	{
		pam_start = coord$end-options$pam_len+1
		options$pam_region = c(pam_start, pam_start+options$pam_len-1)
		options$seed_region = c(pam_start-options$seed_len, pam_start-1)
	}
	else
	{
		pam_end = coord$start+options$pam_len-1
		options$pam_region = c(coord$start, pam_end)
		options$seed_region = c(pam_end+1, pam_end+options$seed_len)
	}

	seqs = Views(options$genome_seq[[coord$chr]], options$target_region[1], options$target_region[2])
	return (find_guards_in_seqs(options, seqs))
}

find_guards =
function(options, on_coord, off_coord)
{
	id = paste0(on_coord$id, "_", off_coord$id);
	out_prefix = paste0(options$out_dir, "/guards_", id)
	out_file = paste0(out_prefix, "_final.txt")

	if (file.exists(out_file))
	{
		return (read.table(out_file, sep="\t", header=TRUE, stringsAsFactors=FALSE))
	}

	tmp_file = paste0(out_prefix, "_tmp.txt")

	if (file.exists(tmp_file))
	{
		off_guards = read.table(tmp_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
	}
	else
	{
		#on_guards = find_guard_seqs(options, on_coord)
		off_guards = find_guard_seqs(options, off_coord)
		write.table(off_guards, file=tmp_file, sep="\t", row.names=FALSE, quote=FALSE)
	}

	guides_file = paste0(out_prefix, "_guides.txt")
	writeLines(off_guards$Guard, guides_file)

	if (options$dry_run)
	{
		return (off_guards)
	}

	off_guards_ot_stats = read.table(paste0(guides_file, ".stderr"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
	off_guards_ot_hits = read.table(paste0(guides_file, ".stdout"), sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)

	off_guards = merge(off_guards, off_guards_ot_stats, by.x="Guard", by.y="Query", all=TRUE)
	off_guards = off_guards[, !(colnames(off_guards) %in% c('Query', 'QueryNo', 'QuerySize'))]

	on_hits = unique(off_guards_ot_hits[ off_guards_ot_hits$Chr == on_coord$chr &
					 # !(off_guards_ot_hits$End < on_coord$start | on_coord$end < off_guards_ot_hits$Start) # partial overlap
					 (on_coord$start <= off_guards_ot_hits$Start & off_guards_ot_hits$End <= on_coord$end) # total overlap
					, c("QuerySeq", "HitSeq", "HitMismatches", "SeedMismatches") ])
	colnames(on_hits) = c("Guard", "OnTargetHit", "OnTargetMismatches", "OnTargetSeedMismatches")
	off_guards = merge(off_guards, on_hits, by="Guard", all=TRUE)

	write.table(off_guards, file=paste0(out_prefix, "_final.txt"), sep="\t", row.names=FALSE, quote=FALSE)

	return (off_guards)
}

get_guide =
function(options, coord)
{
	guide = Views(options$genome_seq[[coord$chr]], coord$start, coord$end)
	if (coord$strand == "-")
	{
		guide = reverseComplement(guide)
	}
	guide = as.character(guide)
	guide = substring(guide, 1, nchar(guide)-options$pam_len)
	message(guide, " ", coord$id)
	return(guide)
}

mask_mismatches =
function(a, b)
{
	a = strsplit(a, split="")[[1]]
	b = strsplit(b, split="")[[1]]
	a[ a != b ] = "N"
	return (paste(a, collapse=""))
}

find_all =
function(options, on_coord, off_coords)
{
	options$on_guide = get_guide(options, on_coord)

	guards = NULL
	for(off_coord in off_coords)
	{
		options$off_guide = get_guide(options, off_coord)
		options$off_guide_gc = calc_gc_content(mask_mismatches(options$off_guide, options$on_guide))
		guards = rbind(guards, find_guards(options, on_coord, off_coord))
	}

	if (options$dry_run)
	{
		return (guards)
	}

	#guards$PAMScore = ifelse(guards$OffGuidePAMOverlap > 0, 1, 0)
	#guards$SeedScore = ifelse(guards$OffGuideSeedOverlap > 0, 1, 0)
	guards$PAMAndSeedScore = ifelse(guards$OffGuidePAMOverlap > 0 | guards$OffGuideSeedOverlap > 0, 1, 0)
	guards$OnTargetScore = ifelse(is.na(guards$OnTargetHit), 1, 0)
	guards$OffTargetScore = (100*(guards$X0-1) + 10*guards$X1 + guards$X2)/10000
	guards$GCScore = (guards$GuardGC - guards$OffGuideGC) / guards$GuardGC

	guards$Score = guards$PAMAndSeedScore +
				   guards$OnTargetScore -
				   guards$OffTargetScore +
				   guards$GCScore

	write.table(guards, file=paste0(options$out_dir, "/", on_coord$id, "_final.txt"), sep="\t", row.names=FALSE, quote=FALSE)
	return (guards)
}

find_from_on_target =
function(options, id="G", on_coord=NULL, on_guide=NULL, guides_file="guide.txt")
{
	if (is.null(on_guide))
	{
		on_guide = get_guide(options, on_coord)
	}

	out_prefix = paste0(options$out_dir, "/guide_", id)

	ot_stats = read.table(paste0(guides_file, ".stderr"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
	ot_hits = read.table(paste0(guides_file, ".stdout"), sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)

	if (is.null(on_coord))
	{
		exact_hits = unique(ot_hits[ ot_hits$HitMismatches == 0, c("Chr", "Start", "End", "Strand"), drop=FALSE])
		if (nrow(exact_hits) == 1)
		{
			hit = exact_hits[1,]
			on_coord = list(id=id, chr=hit$Chr, start=hit$Start, end=hit$End, strand=hit$Strand)
		}
		else if (nrow(exact_hits) > 1)
		{
			message("Too many (", nrow(exact_hits), ") exact hits for ", on_guide, ":")
			for(i in 1:nrow(exact_hits))
			{
				hit = exact_hits[i, ]
				message("    ", i, ". chr=\"", hit$Chr, "\" start=", hit$Start, " end=", hit$End, " strand=\"", hit$Strand, "\"")
			}
			return ()
		}
		else if (nrow(exact_hits) == 0)
		{
			message("No exact hits for ", on_guide)
			return ()
		}
	}

	message("All off-targets: ", nrow(ot_hits))
	# Remove on-target
	ot_hits = ot_hits[ ot_hits$HitMismatches > 0,, drop=FALSE ]
	# Remove OTs based on min p-value
	ot_hits = ot_hits[ ot_hits$pOffTarget >= options$crrna_min_pvalue &
					   ot_hits$HitMismatches <= options$crrna_mismatches,, drop=FALSE ]
	message("Filtered off-targets: ", nrow(ot_hits))

	if (nrow(ot_hits) == 0)
	{
		message("No off-targets")
		return ()
	}

	off_coords = list()
	for(i in 1:nrow(ot_hits))
	{
		row = ot_hits[i, ]
		id = sprintf("OT%04d", i)
		if (row$GeneSymbol != "" && row$GeneSymbol != "-")
		{
			id = paste0(id, "-", row$GeneSymbol)
		}
		off_coords[[ id ]] =
			list(id=id,
				 chr=row$Chr,
				 start=row$Start,
				 end=row$End,
				 strand=row$Strand)
	}
	find_all(options, on_coord, off_coords)
}

main =
function()
{
	option_list = list(
		make_option(c("-n", "--dry_run"), action="store_true", default=FALSE,
					help="Dry-run (default %default)"),
		make_option(c("-g", "--genome"), action="store", default="hg38", type='character',
					help="Genome: hg38, mm10, cho, cho1 (default %default)"),
		make_option(c("-p", "--pam"), action="store", default="NRG", type='character',
					help="PAM (default %default)"),
		make_option(c("-l", "--crrna_len"), action="store", default=20, type='integer',
					help="Guide length (default %default)"),
		make_option(c("--crrna_min_pvalue"), action="store", default=0.0, type='numeric',
					help="Minimum pOffTarget (default %$default)"),
		make_option(c("--crrna_mismatches"), action="store", default=5, type='integer',
					help="Guide mismatches (default %$default)"),
		make_option(c("--guard_len"), action="store", default=14, type='integer',
					help="Guard length (default %$default)"),
		make_option(c("--guard_mismatches"), action="store", default=3, type='integer',
					help="Guard mismatches (default %$default)"),
		make_option(c("--seed_len"), action="store", default=12, type='integer',
					help="Seed length (default %$default)"),
		make_option(c("--out_dir"), action="store", default=".", type='character',
					help="Output directory, (default %$default)"),
		make_option(c("--id"), action="store", default="guide", type='character',
					help="Identifier, prefix for output files (default %$default)"),
		make_option(c("--crrna"), action="store", default=NA, type='character',
					help="Guide"),
		make_option(c("--chr"), action="store", default="", type='character',
					help="Chromosome"),
		make_option(c("--start"), action="store", default="", type='integer',
					help="Start position (including PAM)"),
		make_option(c("--end"), action="store", default="", type='integer',
					help="End position (including PAM)"),
		make_option(c("--strand"), action="store", default="", type='character',
					help="Strand (+ or -)")
	)

	opt = parse_args2(OptionParser(option_list=option_list))
	options = opt$options

	options = init_options(options)

	if (options$chr != "" &&
		options$start != "" &&
		options$end != "" &&
		options$strand != "")
	{
		on_coord = list(id=options$id, chr=options$chr, start=options$start, end=options$end, strand=options$strand)
		x = find_from_on_target(options, options$id, on_coord=on_coord)
	}
	else if (!is.na(options$crrna))
	{
		x = find_from_on_target(options, options$id, on_guide=options$crrna)
	}
}

main()
