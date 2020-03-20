#include <stdlib.h>
#include <stdio.h>
#define _GNU_SOURCE
#include <pthread.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

#define MAX_BUFF 1024
#define MAX_INDEXES 10240
#define MAX_SIZE 20
#define MAX_THREADS 64
#define MAX_SEQS 3000000
#define MAX_QUERIES 100000
#define PAM_RIGHT 0x1ull << 63
#define TRUE 1
#define FALSE 0
#define FLIP(x) x = (x == '+' ? '-' : '+');

typedef enum side_e { LEFT = 0, RIGHT = 1 } side_e;

#define NPAMS		20
typedef enum pam_e
{

	NGG			= 0,
	NGCG		= 1,
	NGAG		= 2,
	NGAN		= 3,
	NNGRRT		= 4,
	NNGRRN		= 5,
	NNAGAAW		= 6,
	NNNNGHTT	= 7,

	NGA			= 8,
	NAG			= 9,
	NGNG		= 10,
	NNGRRW		= 11,
	NGGRRT		= 12,
	NNNNGATT	= 13,
	NNNNGMTT	= 14,

	/* These left-end PAMs are indexed as NAA and NAAA, and the results flipped strand on output */
	TTN			= 15,
	TTTN		= 16,

	NG			= 17,
	YG			= 18,
	NNG			= 19

} pam_e;
int pam_sizes[NPAMS]			= { 3, 4, 4, 4, 6, 6, 7, 8, 3, 3, 4, 6, 6, 8, 8, 3, 4, 2, 2, 3 };
int pam_distance_to_cut[NPAMS]	= { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 22, 22, 3, 3, 3 };
side_e pam_sides[NPAMS]			= { RIGHT, RIGHT, RIGHT, RIGHT, RIGHT, RIGHT, RIGHT, RIGHT, RIGHT, RIGHT, RIGHT, RIGHT, RIGHT, RIGHT, RIGHT, LEFT, LEFT, RIGHT, RIGHT, RIGHT };

#define EXON				0x1
#define CDS					0x2
#define UTR5				0x4
#define UTR3				0x8
#define INTRON				0x10
#define FLANK5				0x20
#define FLANK3				0x40
#define ENHANCER			0x80
#define PROMOTER			0x100
#define TSS					0x200
#define TTS					0x400
#define SPLICE				0x800
#define FIRST_EXON			0x1000
#define LAST_EXON			0x2000
#define CONSTITUTIVE_EXON	0x4000

typedef enum format_e { NO_FORMAT, GFF_FORMAT, TABLE_FORMAT, ANNOTATED_FORMAT, STATS_FORMAT, SHORT_ANNOTATED_FORMAT, EXACT_ANNOTATED_FORMAT } format_e;
typedef enum cut_region_e { NONCODING = 0, CODING = 1 } cut_region_e;

/* Priors from "Lei S. Qi Cell 2013" 
 * Results for a crRNA of 20 !!! The values are the activities of repression still remaining despite
 * one mismatch in a specific position.
 * The last position is the one which is the closest from the PAM motif (reverse in the figure of the article).
 */
double priors[MAX_SIZE] = {0.79, 0.415, 0.45, 0.64, 0.66, 0.695, 0.67, 0.64, 0.3, 0.26,
						   0.25, 0.3, 0.245, 0.01, 0.12, 0.11, 0.07, 0.06, 0.09, 0.065};

typedef unsigned long int uint64_t;
typedef unsigned int uint32_t;
typedef unsigned short int uint16_t;
typedef unsigned char uint8_t;
typedef struct entry_t
{
	uint32_t pos;
	uint64_t seq;
} entry_t;

typedef struct feature_entry_t
{
	char		gene[30];
	uint16_t	type;
	uint8_t		strand;
	uint32_t	start, end, ntx;
} feature_entry_t;

typedef struct gene_entry_t
{
	char		gene[30];
	char		name[40];
	char		type[40];
	uint8_t		strand;
	uint32_t	start, end, ntx;
} gene_entry_t;

typedef struct info_t
{
	feature_entry_t	features[1000000];
	int				nfeatures;
	gene_entry_t	genes[1000000];
	int				ngenes;
} info_t;

char *prog;
bool count_only = FALSE;
format_e format = NO_FORMAT;
uint8_t cmap[256];
uint64_t ERROR_STR = -1;
entry_t seqs[MAX_SEQS];

typedef struct search_data_t
{
	bool		is_free;
	int			thread_no;
	int			index_no;
/*
	int			query_no;
	uint64_t	query;
	int			query_size;
	int			max_mismatches;
*/
	char		*name;
	int			nseqs;
	entry_t		seqs[MAX_SEQS];
	int			max_size;
	pam_e		pam;
	info_t		*info;
	int			mismatches[MAX_SIZE+1];
	int			off_targets;
} search_data_t;

pthread_t threads[MAX_THREADS];
search_data_t search_data[MAX_THREADS];

int nqueries = 0;
int query_sizes[MAX_QUERIES];
uint64_t queries[MAX_QUERIES];
int max_mismatches[MAX_SIZE+1];
uint32_t query_mismatches[MAX_QUERIES][MAX_SIZE+1];
double query_prob[MAX_QUERIES][2];
int query_count[MAX_QUERIES][2];
int query_exact[MAX_QUERIES];
int seed_length = 12;
int max_seed_mismatches = 2;
int bad_seed[MAX_QUERIES];

info_t *read_gene_info(FILE *fp, info_t *info)
{
	char *line, buffer[MAX_BUFF];
	info->ngenes = 0;
	while ((line = fgets(buffer, MAX_BUFF, fp)))
	{
		size_t len = strlen(line);
		line[--len] = '\0';
		char *gene = strtok(line, "\t");
		char *chr = strtok(NULL, "\t");
		int start = atoi(strtok(NULL, "\t"));
		int end = atoi(strtok(NULL, "\t"));
		char *strand = strtok(NULL, "\t");
		int ntx = atoi(strtok(NULL, "\t"));
		char *name = strtok(NULL, "\t");
		char *type = strtok(NULL, "\t");

		gene_entry_t *inf = &info->genes[info->ngenes++];
		strcpy(inf->gene, gene);
		strcpy(inf->name, name);
		strcpy(inf->type, type);
		inf->start = start;
		inf->end = end;
		inf->strand = *strand;
		inf->ntx = ntx;
	}

	return (info);
}

info_t *read_feature_info(FILE *fp)
{
	info_t *info = (info_t *) malloc(sizeof(info_t));
	char *line, buffer[MAX_BUFF];
	info->nfeatures = 0;
	while ((line = fgets(buffer, MAX_BUFF, fp)))
	{
		size_t len = strlen(line);
		line[--len] = '\0';
		char *gene = strtok(line, "\t");
		int type = atoi(strtok(NULL, "\t"));
		int start = atoi(strtok(NULL, "\t"));
		int end = atoi(strtok(NULL, "\t"));
		char *strand = strtok(NULL, "\t");
		int ntx = atoi(strtok(NULL, "\t"));

		feature_entry_t *inf = &info->features[info->nfeatures++];
		strcpy(inf->gene, gene);
		inf->type = type;
		inf->start = start;
		inf->end = end;
		inf->strand = *strand;
		inf->ntx = ntx;
	}

	return (info);
}

gene_entry_t *find_gene(info_t *info, char *gene_id)
{
	int j;
	gene_entry_t *inf;

	for(j = 0; j < info->ngenes; j++)
	{
		inf = &info->genes[j];
		if (strcmp(inf->gene, gene_id) == 0)
		{
			return (inf);
		}
	}
	return (NULL);
}

cut_region_e find_features(int query_no, char *query_seq, char *name,
						   uint32_t hit_start, uint32_t hit_end, char hit_strand, char *hit_seq,
						   int hit_mm, int seed_mm, double p_ot,
						   info_t *info, uint32_t pos)
{
	int i, j;
	feature_entry_t hits[100];
	int nhits = 0;
	feature_entry_t *inf;

	if (info)
	{
	for(j = 0; j < info->nfeatures; j++)
	{
		inf = &info->features[j];

		if (pos < inf->start)
		{
			break;
		}
		if (inf->start <= pos && pos <= inf->end)
		{
			int hit = -1;
			for(i = 0; i < nhits; i++)
			{
				if (strcmp(hits[i].gene, inf->gene) == 0)
				{
					hit = i;
					break;
				}
			}
			if (hit == -1)
			{
				hit = nhits++;
				strcpy(hits[i].gene, inf->gene);
				hits[i].type = 0;
				hits[i].ntx = 0;
				hits[i].strand = inf->strand;
			}

			hits[hit].type |= inf->type;
			if ((inf->type & INTRON) && (pos <= inf->start+2 || inf->end-2 <= pos))
			{
				hits[hit].type |= SPLICE;
			}
			hits[hit].ntx += inf->ntx;
		}
	}
	}

	if (format == ANNOTATED_FORMAT || format == SHORT_ANNOTATED_FORMAT || format == EXACT_ANNOTATED_FORMAT)
	{
		if (nhits == 0)
		{
			if (format == ANNOTATED_FORMAT || (format == EXACT_ANNOTATED_FORMAT && hit_mm == 0))
			{
				printf("%d\t%s\t%s\t%u\t%u\t%c\t%s\t%d\t%d\t%f\t-\t-\t-\n", query_no, query_seq,name, hit_start, hit_end, hit_strand, hit_seq, hit_mm, seed_mm, p_ot);
			}
			else if (format == SHORT_ANNOTATED_FORMAT)
			{
				printf("%s\t%u\t%u\t%c\t%s\t-\t-\t-\n", name, hit_start, hit_end, hit_strand, hit_seq);
			}
			return (NONCODING);
		}
		else
		{
			uint16_t type = 0;
			for(i = 0; i < nhits; i++)
			{
				inf = &hits[i];
				if (inf->type == FLANK5) continue;
				if (inf->type == FLANK3) continue;
				gene_entry_t *gene = find_gene(info, inf->gene);
				char regions[MAX_BUFF];
				regions[0] = '\0';
				regions[1] = '\0';

				if (inf->type & EXON) strcat(regions, ",exon");
				if (inf->type & CDS) strcat(regions, ",cds");
				if (inf->type & UTR5) strcat(regions, ",5'utr");
				if (inf->type & UTR3) strcat(regions, ",3'utr");
				if (inf->type & INTRON) strcat(regions, ",intron");
				if (inf->type & SPLICE) strcat(regions, ",splice");
				if (inf->type & FLANK5) strcat(regions, ",5'flank");
				if (inf->type & FLANK3) strcat(regions, ",3'flank");
				if (inf->type & ENHANCER) strcat(regions, ",enhancer");
				if (inf->type & PROMOTER) strcat(regions, ",promoter");
				if (inf->type & TSS) strcat(regions, ",tss");
				if (inf->type & TTS) strcat(regions, ",tts");
				if (inf->type & FIRST_EXON) strcat(regions, ",first_exon");
				if (inf->type & LAST_EXON) strcat(regions, ",last_exon");
				if (inf->type & CONSTITUTIVE_EXON) strcat(regions, ",constitutive_exon");

				if (format == ANNOTATED_FORMAT || (format == EXACT_ANNOTATED_FORMAT && hit_mm == 0))
				{
					printf("%d\t%s\t%s\t%u\t%u\t%c\t%s\t%d\t%d\t%f\t%s\t%s\t%s\n",
						   query_no, query_seq,
						   name, hit_start, hit_end, hit_strand, hit_seq, hit_mm, seed_mm, p_ot,
						   inf->gene, gene->name, regions+1);
				}
				else if (format == SHORT_ANNOTATED_FORMAT)
				{
					printf("%s\t%u\t%u\t%c\t%s\t%s\t%s\t%s\n",
						   name, hit_start, hit_end, hit_strand, hit_seq,
						   inf->gene, gene->name, regions+1);
				}
				type |= inf->type;
			}
			if (type == 0) /* must be all flanking */
			{
				if (format == ANNOTATED_FORMAT || (format == EXACT_ANNOTATED_FORMAT && hit_mm == 0))
				{
					printf("%d\t%s\t%s\t%u\t%u\t%c\t%s\t%d\t%d\t%f\t-\t-\t-\n", query_no, query_seq,name, hit_start, hit_end, hit_strand, hit_seq, hit_mm, seed_mm, p_ot);
				}
				else if (format == SHORT_ANNOTATED_FORMAT)
				{
					printf("%s\t%u\t%u\t%c\t%s\t-\t-\t-\n", name, hit_start, hit_end, hit_strand, hit_seq);
				}
			}
			return ((type == 0 || type == INTRON) ? NONCODING : CODING);
		}
	}
	else
	{
		if (nhits == 0)
		{
			return (NONCODING);
		}
		else
		{
			uint16_t type = 0;
			for(i = 0; i < nhits; i++)
			{
				inf = &hits[i];
				if (inf->type & FLANK5) continue;
				if (inf->type & FLANK3) continue;
				type |= inf->type;
			}
			return ((type == 0 || type == INTRON) ? NONCODING : CODING);
		}
	}
}

void populate_cmap()
{
    /* create an array for all possible char values, with 4 as the default value */
	int i;

	for(i = 0; i < 256; i++)
	{
		cmap[i] = 4;
	}

    /* fill in the 8 entries we're actually interested in with their two bit values */
    cmap['a'] = cmap['A'] = 0; /* 00 */
    cmap['c'] = cmap['C'] = 1; /* 01 */
    cmap['g'] = cmap['G'] = 2; /* 10 */
    cmap['t'] = cmap['T'] = 3; /* 11 */
}


char *bits_to_string(uint64_t text, uint64_t match, int size, char *s)
{
	/* have to & 0x3 to turn off all bits but what we actually want. */
	int shift = 2*(size-1); /* there are twice as many bits as there are characters */
	int i;

	/* fill with N if its an error string (all bits set to 1) */
	if (text == ERROR_STR)
	{
		for(i = 0; i < size; i++)
		{
			s[i] = 'N';
		}
		s[i] = '\0';
	}

	/* extract each character from the text */
	for (i = 0; i < size; i++, shift -= 2)
	{
		/* put the character we're interested in at the very end */
		/* of the integer, and switch all remaining bits to 0 with & 0x3 */
		uint8_t character = (text >> shift) & 0x3;
		uint8_t m = (match >> shift) & 0x3;
		switch ( character )
		{
		case 0: s[i] = m ? 'a' : 'A'; break;
		case 1: s[i] = m ? 'c' : 'C'; break;
		case 2: s[i] = m ? 'g' : 'G'; break;
		case 3: s[i] = m ? 't' : 'T'; break;
		default: break;
		}
	}
	s[i] = '\0';

	return s;
}

uint64_t string_to_bits(const char *seq, int len, int base)
{
	/* loop through each character */
	/* optionally allow the user to prepend the sequence with a short */
	uint64_t bits = base;
	int j;

	for(j = 0; j < len; j++)
	{
		uint8_t const c = seq[j]; /* get char */

		if ( cmap[c] == 4 )
		{
			bits = ERROR_STR; /* set to all 1s */
			break;
		}
		else
		{
			bits <<= 2; /* shift left to make room for new char */
			bits |= cmap[c]; /* add our new char to the end */
		}
	}

	return bits;
}

uint64_t rev(uint64_t text, int size)
{
	uint64_t reversed = 0;
	int shift = 0;

	/* now reverse the sequence, 2 bits at a time */
	/* could just do shift = 0; shift < size*2; shift += 2 */
	int i;
	for (i = 0; i < size; i++, shift += 2)
	{
		reversed <<= 2;
		reversed |= ( text >> shift ) & 0x3;
	}

	return reversed;
}

uint64_t revcom(uint64_t text, int size)
{
	/* for a size of 23 we & with 23 1s to undo the complement on the part of */
	/* the integer we aren't interested in, for consistency */
	/* i.e. set the unused bits back to 0. */
	/* assumes size is smaller than sizeof(text)... */
	unsigned int num_bits = sizeof(text) * 8;
	/* we have to -1 from here to account for the pam_right bit which we DO want flipped */
	uint64_t mask = 0xFFFFFFFFFFFFFFFFull >> ( (num_bits - (size * 2)) - 1 );
	/* bit complement here maps A -> T, G -> C etc. */
	text = ~text & mask;

	uint64_t reversed = 0;
	int shift = 0;

	/* now reverse the sequence, 2 bits at a time */
	/* could just do shift = 0; shift < size*2; shift += 2 */
	int i;
	for (i = 0; i < size; i++, shift += 2)
	{
		reversed <<= 2;
		reversed |= ( text >> shift ) & 0x3;
	}

	return reversed;
}

int pop_count(uint64_t x)
{
	/* as everything is two bit we must convert them all to one bit, */
	/* to do this we must turn off all MSBs, but before we can do that */
	/* we need to ensure that when an MSB is set to 1, the LSB is also set. */
	/* the 4 will change pam_right to 0 so it doesnt get counted. */
	/* 5 is 0101, 4 is 0100 */
	/* x = (x | (x >> 1)) & (0x5555545555555555ull); */
	x = (x | (x >> 1)) & (0x5555555555555555ull);

	x = x-( (x>>1) & 0x5555555555555555ull);
	x = (x & 0x3333333333333333ull) + ((x >> 2) & 0x3333333333333333ull);
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0full;
	return  (0x0101010101010101ull*x >> 56);
}

double pam_to_prob(uint64_t u, pam_e pam)
{
	int pam_size = pam_sizes[pam];

	switch (pam)
	{
	case NGG:
		/* xNN = 00 1111 = 0xF */
		u &= 0xFull;
		/* xRG = 00 xx10 */
		/* xGG = 00 1010 = 0xA */
		/* xAG = 00 0010 = 0x2 */
		return (u == 0xAull ? 1.0 : 0.2); /* GG = 1.0, AG = 0.2 */
	case NGGRRT:
		/* xNNNNN = 0011 1111 1111 = 0x3FF */
		u &= 0x3FFull;
		/* xGG R R T = 0010 10xx xx11 */
		/* xGG G G T = 0010 1010 1011 = 0x2AB */
		/* xGG A G T = 0010 1000 1011 = 0x28B */
		/* xGG G A T = 0010 1010 0011 = 0x2A3 */
		/* xGG A A T = 0010 1000 0011 = 0x283 */
		/* return (u == 0x2ABull || u == 0x28Bull || u == 0x2A3ull || u == 0x283ull); */
		return (1.0);
	case NG:
	case YG:
	case NNG:
	case NGCG:
	case NGAG:
	case NGA:
	case NGNG:
	case NNGRRT:
	case NNGRRW:
	case NNNNGATT:
	case NNNNGMTT:
	case NNAGAAW:
	case NGAN:
	case NNGRRN:
	case NNNNGHTT:
	case TTN:
	case TTTN:
		return (1.0);
	case NAG:
		return (0.2);
	default:
		fprintf(stderr, "# Unimplemented PAM: %d\n", pam);
		exit(1);
	}
}

double match_to_prob(uint64_t match, int size, uint64_t seq, pam_e pam)
{
	if (size > 20) size = 20;
	int i, shift = 2*(size-1);
	double p = 1.0, *ps = &priors[20-size];

	for (i = 0; i < size; i++, shift -= 2)
	{
		uint8_t m = (match >> shift) & 0x3;
		if (m)
		{
			p *= ps[i];
		}
	}
	p *= pam_to_prob(seq, pam);

	return (p);
}

void record_match_prob(int query_no, cut_region_e cut_region, double p)
{
	/* If this is the first exact match (p=1.0) then do not include p=1.0 - assume this is the intended target */
	if (p == 1.0 && !query_exact[query_no])
	{
		query_exact[query_no] = 1;
		return;
	}
	query_prob[query_no][cut_region] *= 1.0-p;
	query_count[query_no][cut_region]++;
	/* printf("%d %d %d %f %f\n", query_no, cut_region, query_mismatches[query_no][0], p, query_prob[query_no][cut_region]); */
}

int search(int query_no, uint64_t query, int query_size,
		   char *name, int num_seqs, entry_t *seqs, int max_size, pam_e pam,
		   info_t *info,
		   int max_mismatches, int *mismatches)
{
	int off_targets = 0;
	int j;
	int delta = max_size-query_size;
	uint64_t mask = 0xFFFFFFFFFFull >> (2*delta);
	uint64_t seed_mask = 0xFFFFFFFFFFull >> (2*(max_size-seed_length));
	char query_buff[MAX_BUFF];
	char buff[MAX_BUFF];
	int pam_size = pam_sizes[pam];
	int distPAMcut = pam_distance_to_cut[pam];
	side_e pam_side = pam_sides[pam];
	char *query_seq = bits_to_string(query, 0ull, query_size, query_buff);

	for(j = 0; j < num_seqs; j++)
	{
		uint64_t s = seqs[j].seq;

		/* xor the two bit strings */
		if (s & PAM_RIGHT)
		{
			char strand = (pam_side == LEFT ? '-' : '+');
			s >>= 2*pam_size;
			s &= mask;
			uint64_t match = query ^ s;
			int mm = pop_count(match);
			if (mm <= max_mismatches)
			{
				uint64_t seed = match & seed_mask;
				int seed_mm = pop_count(seed);
				if (!count_only)
				{
					s = seqs[j].seq;
					char *hit_seq;
					if (pam_side == RIGHT)
					{
						hit_seq = bits_to_string(s, match << (2*pam_size), query_size+pam_size, buff);
					}
					else
					{
						uint64_t match_rev = rev(match, query_size);
						uint64_t seq_rev = revcom(s, query_size+pam_size);
						hit_seq = bits_to_string(seq_rev, match_rev, query_size+pam_size, buff);
					}
					uint32_t hit_start = seqs[j].pos+delta;
					uint32_t hit_end = seqs[j].pos+max_size+pam_size-1;
					uint32_t cut_pos = seqs[j].pos+max_size-distPAMcut;

					if (format == TABLE_FORMAT)
					{
						printf("%s\t%u\t%u\t+\t%s\t%d\t%d\n", name, hit_start, hit_end, hit_seq, mm, seed_mm);
					}
					else if (format == ANNOTATED_FORMAT || format == EXACT_ANNOTATED_FORMAT || format == STATS_FORMAT)
					{
						double p = match_to_prob(match, query_size, seqs[j].seq, pam);
						cut_region_e cut_region = find_features(query_no+1, query_seq, name, hit_start, hit_end, strand, hit_seq, mm, seed_mm, p, info, cut_pos);
						record_match_prob(query_no, cut_region, p);
					}
					else
					{
						printf("%d\t%s\t%u\t+\t%d\t%s\n", query_no, name, hit_start, mm, hit_seq);
					}
					fflush(stdout);
				}
				mismatches[mm]++;
				query_mismatches[query_no][mm]++;
				if (seed_mm <= max_seed_mismatches) bad_seed[query_no]++;
				off_targets++;
			}
		}
		else
		{
			char strand = (pam_side == LEFT ? '+' : '-');
			uint64_t query_rev = revcom(query, query_size);

			s >>= 2*delta;
			s &= mask;
			uint64_t match = query_rev ^ s;
			int mm = pop_count(match);

			if (mm <= max_mismatches)
			{
				uint64_t match_rev = rev(match, query_size);
				uint64_t seed = match_rev & seed_mask;
				int seed_mm = pop_count(seed);
				if (!count_only)
				{
					s = seqs[j].seq >> (2*delta);
					uint64_t seq_rev = revcom(s, query_size+pam_size);
					char *hit_seq;
					if (pam_side == RIGHT)
					{
						hit_seq = bits_to_string(seq_rev, match_rev << (2*pam_size), query_size+pam_size, buff);
					}
					else
					{
						hit_seq = bits_to_string(s, match, query_size+pam_size, buff);
					}
					uint32_t hit_start = seqs[j].pos;
					uint32_t hit_end = seqs[j].pos+max_size+pam_size-delta-1;
					uint32_t cut_pos = seqs[j].pos+pam_size+distPAMcut;

					if (format == TABLE_FORMAT)
					{
						printf("%s\t%u\t%u\t-\t%s\t%d\t%d\n", name, hit_start, hit_end, hit_seq, mm, seed_mm);
					}
					else if (format == ANNOTATED_FORMAT || format == EXACT_ANNOTATED_FORMAT || format == STATS_FORMAT)
					{
						double p = match_to_prob(match_rev, query_size, seq_rev, pam);
						cut_region_e cut_region = find_features(query_no+1, query_seq, name, hit_start, hit_end, strand, hit_seq, mm, seed_mm, p, info, cut_pos);
						record_match_prob(query_no, cut_region, p);
					}
					else
					{
						printf("%d\t%s\t%u\t-\t%d\t%s\n", query_no, name, hit_start, mm, hit_seq);
					}
					fflush(stdout);
				}
				mismatches[mm]++;
				query_mismatches[query_no][mm]++;
				if (seed_mm <= max_seed_mismatches) bad_seed[query_no]++;
				off_targets++;
			}
		}
	}

	return (off_targets);
}

void *search_thread(void *p)
{
	search_data_t *d = (search_data_t *) p;
	int query_no;

	d->off_targets = 0;
	bzero(d->mismatches, sizeof(d->mismatches));
	for(query_no = 0; query_no < nqueries; query_no++)
	{
		d->off_targets += search(query_no, queries[query_no], query_sizes[query_no],
			   					 d->name, d->nseqs, d->seqs, d->max_size, d->pam, d->info,
			   					 max_mismatches[query_sizes[query_no]], d->mismatches);
	}
/*
	d->off_targets = search(d->query_no, d->query, d->query_size,
		   					d->name, d->nseqs, d->seqs, d->max_size, d->pam,
		   					d->max_mismatches, d->mismatches);
*/
}

void init_threads()
{
	int i;

	for(i = 0; i < MAX_THREADS; i++)
	{
		search_data[i].is_free = TRUE;
	}
}

void join_data(int thread_no, int *off_targets, int *mismatches)
{
	int index_no = search_data[thread_no].index_no;
	off_targets[index_no] += search_data[thread_no].off_targets;
	int k;
	for(k = 0; k < MAX_SIZE; k++)
	{
		mismatches[k] += search_data[thread_no].mismatches[k];
	}
	search_data[thread_no].is_free = TRUE;
}

void join_threads(int nthreads, uint32_t *off_targets, int *mismatches)
{
	int thread_no;
	for(thread_no = 0; thread_no < nthreads; thread_no++)
	{
		if (search_data[thread_no].is_free)
		{
			continue;
		}
		if (pthread_join(threads[thread_no], NULL) == 0)
		{
			join_data(thread_no, off_targets, mismatches);
		}
	}
}

int find_thread(int nthreads, uint32_t *off_targets, int *mismatches)
{
	int thread_no = nthreads+1;
	while (thread_no >= nthreads)
	{
		for(thread_no = 0; thread_no < nthreads; thread_no++)
		{
			if (search_data[thread_no].is_free)
			{
				return (thread_no);;
			}
		}
		for(thread_no = 0; thread_no < nthreads; thread_no++)
		{
			if (pthread_tryjoin_np(threads[thread_no], NULL) == 0)
			{
				join_data(thread_no, off_targets, mismatches);
			}
		}
	}

	return (thread_no);
}

void create_thread(int thread_no, search_data_t *d)
{
	if (pthread_create(&threads[thread_no], NULL, search_thread, d) != 0)
	{
		fprintf(stderr, "# Failed to create thread %d\n", thread_no);
		exit(1);
	}
}

void search_index(int index_no, FILE *fp, char *chr,
				  int *mismatches, uint32_t *index_off_targets,
				  int nthreads, pam_e pam,
				  info_t *info)
{
	int nread;

	if (nthreads > 1)
	{
		/* Find free slot */
		int thread_no = find_thread(nthreads, index_off_targets, mismatches);
		search_data_t *d = &search_data[thread_no];
		while ((d->nseqs = fread(d->seqs, sizeof(entry_t), MAX_SEQS, fp)) > 0)
		{
			d->is_free = FALSE;
			d->thread_no = thread_no;
			d->index_no = index_no;
			d->name = chr;

			d->max_size = MAX_SIZE;
			d->pam = pam;
			d->info = info;

			create_thread(thread_no, d);

			thread_no = find_thread(nthreads, index_off_targets, mismatches);
			d = &search_data[thread_no];
		}
	}
	else
	{
		while ((nread = fread(seqs, sizeof(entry_t), MAX_SEQS, fp)) > 0)
		{
			int query_no;
			for(query_no = 0; query_no < nqueries; query_no++)
			{
				index_off_targets[index_no] +=
					search(query_no, queries[query_no], query_sizes[query_no],
						   chr, nread, seqs, MAX_SIZE, pam, info,
						   max_mismatches[query_sizes[query_no]], mismatches);
			}
		}
	}
}

bool pam_left(uint64_t u, pam_e pam)
{
	switch (pam)
	{
	case NG:
		/* Nx = 1100 = 0xC */
		u &= 0xCull;
		/* Cx = 0100 = 0x4 */
		return (u == 0x4ull);
	case YG:
		/* NN = 1111 = 0xF */
		u &= 0xFull;
		/* CG = 0110 = 0x6 */
		/* CA = 0100 = 0x4 */
		return (u == 0x6ull || u == 0x4ull);
	case NNG:
		/* Nxx = 11 0000 = 0x30 */
		u &= 0x30ull;
		/* Cxx = 01 0000 = 0x10 */
		return (u == 0x10ull);
	case NGG:
		/* NNx = 11 1100 = 0x3C */
		u &= 0x3Cull;
		/* CCx = 01 0100 = 0x14 */
		/* CTx = 01 1100 = 0x1C */
		return (u == 0x14ull || u == 0x1Cull);
	case NGGRRT:
		/* NNNNNNx = 1111 1111 1100 = 0xFFC */
		u &= 0xFFCull;
		/* AYYCCx = 00xx xx01 0100 */
		/* ACCCCx = 0001 0101 0100 = 0x154 */
		/* ACTCCx = 0001 1101 0100 = 0x1D4 */
		/* ATCCCx = 0011 0101 0100 = 0x354 */
		/* ATTCCx = 0011 1101 0100 = 0x3D4 */
		return (u == 0x154ull || u == 0x1D4ull || u == 0x354ull || u == 0x3D4ull);
	case NGCG:
		/* NNNx = 1111 1100 = 0xFC */
		u &= 0xFCull;
		/* CGCx = 0110 0100 = 0x64 */
		return (u == 0x64ull);
	case NGAG:
		/* NNNx = 1111 1100 = 0xFC */
		u &= 0xFCull;
		/* CTCx = 0111 0100 = 0x74 */
		return (u == 0x74ull);
	case NGA:
		/* NNx = 11 1100 = 0x3C */
		u &= 0x3Cull;
		/* TCx = 11 0100 = 0x34 */
		return (u == 0x34ull);
	case NGNG:
		/* NxNx = 1100 1100 = 0xCC */
		u &= 0xCCull;
		/* CxCx = 0100 0100 = 0x44 */
		return (u == 0x44ull);
	case NNGRRT:
		/* NNNNxx = 1111 1111 0000 = 0xFF0 */
		u &= 0xFF0ull;
		/* AYYCxx = 00xx xx01 0000 */
		/* ACCCxx = 0001 0101 0000 = 0x150 */
		/* ACTCxx = 0001 1101 0000 = 0x1D0 */
		/* ATCCxx = 0011 0101 0000 = 0x350 */
		/* ATTCxx = 0011 1101 0000 = 0x3D0 */
		return (u == 0x150ull || u == 0x1D0ull || u == 0x350ull || u == 0x3D0ull);
	case NNGRRW:
		/* NNNNxx = 1111 1111 0000 = 0xFF0 */
		u &= 0xFF0ull;
		/* WYYCxx = 00xx xx01 0000 */
		/* ACCCxx = 0001 0101 0000 = 0x150 */
		/* ACTCxx = 0001 1101 0000 = 0x1D0 */
		/* ATCCxx = 0011 0101 0000 = 0x350 */
		/* ATTCxx = 0011 1101 0000 = 0x3D0 */
		/* TCCCxx = 1101 0101 0000 = 0xD50 */
		/* TCTCxx = 1101 1101 0000 = 0xDD0 */
		/* TTCCxx = 1111 0101 0000 = 0xF50 */
		/* TTTCxx = 1111 1101 0000 = 0xFD0 */
		return (u == 0x150ull || u == 0x1D0ull || u == 0x350ull || u == 0x3D0ull ||
				u == 0xD50ull || u == 0xDD0ull || u == 0xF50ull || u == 0xFD0ull);
	case NNNNGATT:
		/* NNNNxxxx = 1111 1111 0000 0000 = 0xFF00 */
		u &= 0xFF00ull;
		/* AATCxxxx = 0000 1101 0000 0000 = 0xD00 */
		return (u == 0xD00ull);
	case NNNNGMTT:
		/* NNNNxxxx = FFFF FFFF 0000 0000 = 0xFF00 */
		u &= 0xFF00ull;
		/* AAKCxxxx = 0000 1101 0000 0000 = 0xD00 */
		/* AATCxxxx = 0000 1101 0000 0000 = 0xD00 */
		/* AAGCxxxx = 0000 1001 0000 0000 = 0x900 */
		return (u == 0xD00ull || u == 0x900ull);
	case NNAGAAW:
		/* NNNNNxx = 11 1111 1111 0000 = 0x3FF0 */
		u &= 0x3FF0ull;
		/* WTTCTxx = xx 1111 0111 0000 */
		/* ATTCTxx = 00 1111 0111 0000 = 0xF70 */
		/* TTTCTxx = 11 1111 0111 0000 = 0x3F70 */
		return (u == 0xF70ull || u == 0x3F70ull);
	case NAG:
		/* NNx = 11 1100 = 0x3C */
		u &= 0x3Cull;
		/* CTx = 01 1100 = 0x1C */
		return (u == 0x1Cull);
	case NNNNGHTT:
		/* NNNNxxxx = FFFF FFFF 0000 0000 = 0xFF00 */
		u &= 0xFF00ull;
		/* AADCxxxx = 0000 xx01 0000 0000 = 0xD00 */
		/* AATCxxxx = 0000 1101 0000 0000 = 0xD00 */
		/* AAGCxxxx = 0000 1001 0000 0000 = 0x900 */
		/* AAACxxxx = 0000 0001 0000 0000 = 0x100 */
		return (u == 0xD00ull || u == 0x900ull || u == 0x100ull);
	case NGAN:
		/* xNNx = 0011 1100 = 0x3C */
		u &= 0x3Cull;
		/* xTCx = 0011 0100 = 0x34 */
		return (u == 0x34ull);
	case NNGRRN:
		/* xNNNxx = 0011 1111 0000 = 0x3F0 */
		u &= 0x3F0ull;
		/* xYYCxx = 00xx xx01 0000 */
		/* xCCCxx = 0001 0101 0000 = 0x150 */
		/* xCTCxx = 0001 1101 0000 = 0x1D0 */
		/* xTCCxx = 0011 0101 0000 = 0x350 */
		/* xTTCxx = 0011 1101 0000 = 0x3D0 */
		return (u == 0x150ull || u == 0x1D0ull || u == 0x350ull || u == 0x3D0ull);
	case TTN:
		/* NNx = 11 1100 = 3C */
		u &= 0x3Cull;
		/* TTx = 11 1100 = 0x3C */
		return (u == 0x3Cull);
	case TTTN:
		/* NNNx = 1111 1100 = FC */
		u &= 0xFCull;
		/* TTTx = 1111 1100 = 0xFC */
		return (u == 0xFCull);
	default:
		fprintf(stderr, "# Unimplemented PAM: %d\n", pam);
		exit(1);
	}
}

bool pam_right(uint64_t u, pam_e pam)
{
	switch (pam)
	{
	case NG:
		/* xN = 0011 = 0x3 */
		u &= 0x3ull;
		/* xG = 00 10 = 0x2 */
		return (u == 0x2ull);
	case YG:
		/* NN = 1111 = 0xF */
		u &= 0xFull;
		/* CG = 01 10 = 0x6 */
		/* TG = 11 10 = 0xE */
		return (u == 0x6ull || u == 0xEull);
	case NNG:
		/* xxN = 00 0011 = 0x3 */
		u &= 0x3ull;
		/* xxG = 00 0010 = 0x2 */
		return (u == 0x2ull);
	case NGG:
		/* xNN = 00 1111 = 0xF */
		u &= 0xFull;
		/* xRG = 00 xx10 */
		/* xGG = 00 1010 = 0xA */
		/* xAG = 00 0010 = 0x2 */
		return (u == 0xAull || u == 0x2ull);
	case NGGRRT:
		/* xNNNNN = 0011 1111 1111 = 0x3FF */
		u &= 0x3FFull;
		/* xGG R R T = 0010 10xx xx11 */
		/* xGG G G T = 0010 1010 1011 = 0x2AB */
		/* xGG A G T = 0010 1000 1011 = 0x28B */
		/* xGG G A T = 0010 1010 0011 = 0x2A3 */
		/* xGG A A T = 0010 1000 0011 = 0x283 */
		return (u == 0x2ABull || u == 0x28Bull || u == 0x2A3ull || u == 0x283ull);
	case NGCG:
		/* xNNN = 0011 1111 = 0x3F */
		u &= 0x3Full;
		/* xGCG = 0010 0110 = 0x26 */
		return (u == 0x26ull);
	case NGAG:
		/* xNNN = 0011 1111 = 0x3F */
		u &= 0x3Full;
		/* xGAG = 0010 0010 = 0x22 */
		return (u == 0x22ull);
	case NGA:
		/* xNN = 00 1111 = 0xF */
		u &= 0xFull;
		/* xGA = 00 1000 = 0x8 */
		return (u == 0x8ull);
	case NGNG:
		/* xNxN = 0011 0011 = 0x33 */
		u &= 0x33ull;
		/* xGxG = 0010 0010 = 0x22 */
		return (u == 0x22ull);
	case NNGRRT:
		/* xxNNNN = 0000 1111 1111 = 0xFF */
		u &= 0xFFull;
		/* xxG R R T = 0000 10xx xx11 */
		/* xxG G G T = 0000 1010 1011 = 0xAB */
		/* xxG A G T = 0000 1000 1011 = 0x8B */
		/* xxG G A T = 0000 1010 0011 = 0xA3 */
		/* xxG A A T = 0000 1000 0011 = 0x83 */
		return (u == 0xABull || u == 0x8Bull || u == 0xA3ull || u == 0x83ull);
	case NNGRRW:
		/* xxNNNN = 0000 1111 1111 = 0xFF */
		u &= 0xFFull;
		/* xxG R R W = 0000 10xx xxxx */
		/* xxG G G A = 0000 1010 1000 = 0xA8 */
		/* xxG A G A = 0000 1000 1000 = 0x88 */
		/* xxG G A A = 0000 1010 0000 = 0xA0 */
		/* xxG A A A = 0000 1000 0000 = 0x80 */
		/* xxG G G T = 0000 1010 1011 = 0xAB */
		/* xxG A G T = 0000 1000 1011 = 0x8B */
		/* xxG G A T = 0000 1010 0011 = 0xA3 */
		/* xxG A A T = 0000 1000 0011 = 0x83 */
		return (u == 0xA8ull || u == 0x88ull || u == 0xA0ull || u == 0x80ull ||
				u == 0xABull || u == 0x8Bull || u == 0xA3ull || u == 0x83ull);
	case NNNNGATT:
		/* xxxxNNNN = 0000 0000 1111 1111 = 0xFF */
		u &= 0xFFull;
		/* xxxxGATT = 0000 0000 1000 1111 = 0x8F */
		return (u == 0x8Full);
	case NNNNGMTT:
		/* xxxxNNNN = 0000 0000 1111 1111 = 0xFF */
		u &= 0xFFull;
		/* xxxx G M T T = 0000 0000 10xx 1111 = 0x8F */
		/* xxxx G A T T = 0000 0000 1000 1111 = 0x8F */
		/* xxxx G C T T = 0000 0000 1001 1111 = 0x9F */
		return (u == 0x8Full || u == 0x9Full);
	case NNAGAAW:
		/* xxNNNNN = 00 0011 1111 1111 = 0x3FF */
		u &= 0x3FFull;
		/* xxA G A A W = 00 0000 1000 00xx */
		/* xxA G A A A = 00 0000 1000 0000 = 0x80 */
		/* xxA G A A T = 00 0000 1000 0011 = 0x83 */
		return (u == 0x80ull || u == 0x83ull);
	case NAG:
		/* xNN = 00 1111 = 0xF */
		u &= 0xFull;
		/* xAG = 00 0010 = 0x2 */
		return (u == 0x2ull);
	case NGAN:
		/* xNNx = 0011 1100 = 0x3C */
		u &= 0x3Cull;
		/* xGAx = 0010 0000 = 0x20 */
		return (u == 0x20ull);
	case NNGRRN:
		/* xxNNNx = 0000 1111 1100 = 0xFC */
		u &= 0xFCull;
		/* xxG R R x = 0000 10xx xx00 */
		/* xxG G G x = 0000 1010 1000 = 0xA8 */
		/* xxG A G x = 0000 1000 1000 = 0x88 */
		/* xxG G A x = 0000 1010 0000 = 0xA0 */
		/* xxG A A x = 0000 1000 0000 = 0x80 */
		return (u == 0xA8ull || u == 0x88ull || u == 0xA0ull || u == 0x80ull);
	case NNNNGHTT:
		/* xxxxNNNN = 0000 0000 1111 1111 = 0xFF */
		u &= 0xFFull;
		/* xxxx G H T T = 0000 0000 10xx 1111 = 0x8F */
		/* xxxx G A T T = 0000 0000 1000 1111 = 0x8F */
		/* xxxx G C T T = 0000 0000 1001 1111 = 0x9F */
		/* xxxx G T T T = 0000 0000 1011 1111 = 0xBF */
		return (u == 0x8Full || u == 0x9Full || u == 0xBFull );
	case TTN:
		/* xNN = 00 1111 = 0xF */
		u &= 0xFull;
		/* xAA = 00 0000 = 0x0 */
		return (u == 0x0ull);
	case TTTN:
		/* xNNN = 0011 1111 = 0x3F */
		u &= 0x3Full;
		/* xAAA = 0000 0000 = 0x0 */
		return (u == 0x0ull);
	default:
		fprintf(stderr, "# Unimplemented PAM: %d\n", pam);
		exit(1);
	}
}

pam_e string_to_pam(char *s)
{
	if      (strcasecmp(s, "NGG") == 0)			{ return (NGG); }
	else if (strcasecmp(s, "NGCG") == 0)		{ return (NGCG); }
	else if (strcasecmp(s, "NGAG") == 0)		{ return (NGAG); }
	else if (strcasecmp(s, "NGAN") == 0)		{ return (NGAN); }
	else if (strcasecmp(s, "NNGRRT") == 0)		{ return (NNGRRT); }
	else if (strcasecmp(s, "NNGRRN") == 0)		{ return (NNGRRN); }
	else if (strcasecmp(s, "NNAGAAW") == 0)		{ return (NNAGAAW); }
	else if (strcasecmp(s, "NNNNGHTT") == 0)	{ return (NNNNGHTT); }

	else if (strcasecmp(s, "NAG") == 0)			{ return (NAG); }
	else if (strcasecmp(s, "NGA") == 0)			{ return (NGA); }
	else if (strcasecmp(s, "NGNG") == 0)		{ return (NGNG); }
	else if (strcasecmp(s, "NGGRRT") == 0)		{ return (NGGRRT); }
	else if (strcasecmp(s, "NNGRRW") == 0)		{ return (NNGRRW); }
	else if (strcasecmp(s, "NNNNGATT") == 0)	{ return (NNNNGATT); }
	else if (strcasecmp(s, "NNNNGMTT") == 0)	{ return (NNNNGMTT); }

	else if (strcasecmp(s, "TTN") == 0)			{ return (TTN); }
	else if (strcasecmp(s, "TTTN") == 0)		{ return (TTTN); }

	else if (strcasecmp(s, "NG") == 0)			{ return (NG); }
	else if (strcasecmp(s, "YG") == 0)			{ return (YG); }
	else if (strcasecmp(s, "NNG") == 0)			{ return (NNG); }

	else
	{
		fprintf(stderr, "# %s: unimplemented PAM\n", prog, s);
		exit(1);
	}
}

char *pam_to_string(pam_e pam)
{
	switch (pam)
	{
	case NGG:		return ("NGG");
	case NGGRRT:	return ("NGGRRT");
	case NGCG:		return ("NGCG");
	case NGAG:		return ("NGAG");
	case NGA:		return ("NGA");
	case NGNG:		return ("NGNG");
	case NNGRRT:	return ("NNGRRT");
	case NNGRRW:	return ("NNGRRW");
	case NNNNGATT:	return ("NNNNGATT");
	case NNNNGMTT:	return ("NNNNGMTT");
	case NNAGAAW:	return ("NNAGAAW");
	case NAG:		return ("NAG");
	case NGAN:		return ("NGAN");
	case NNGRRN:	return ("NNGRRN");
	case NNNNGHTT:	return ("NNNNGHTT");
	case TTN:		return ("TTN");
	case TTTN:		return ("TTTN");
	case NG:		return ("NG");
	case YG:		return ("YG");
	case NNG:		return ("NNG");
	default:
		fprintf(stderr, "# Unimplemented PAM: %d\n", pam);
		exit(1);
	}
}

void index_file(FILE *fp, char *out_dir, int max_size, pam_e pam)
{
	FILE *out_fp = NULL;
	char *line;
	char buffer[MAX_BUFF];
	char name[MAX_BUFF];
	char filename[MAX_BUFF];
	int n = 0;
	uint32_t pos = 0;
	uint64_t u = 0;
	int nu = 0;
	int nseqs = 0;
	int tot_seqs = 0;
	int pam_size = pam_sizes[pam];
	uint64_t all_mask = 0xFFFFFFFFFFFFFFFFull >> (64 - 2*(max_size+pam_size-1));

	name[0] = '\0';
	while (line = fgets(buffer, MAX_BUFF, fp))
	{
		size_t len = strlen(line);
		line[--len] = '\0';
		if (line[0] == '>')
		{
			if (out_fp)
			{
				if (nseqs > 0)
				{
					fwrite(&seqs, sizeof(entry_t), nseqs, out_fp);
				}
				fclose(out_fp);
				fprintf(stderr, "%d bases, %d sequences\n", n, tot_seqs);
			}

			strcpy(name, line+1);
			char *p = strchr(name, ' ');
			if (p) { *p = '\0'; }
			p = strchr(name, '\t');
			if (p) { *p = '\0'; }

			if (strncmp(name, "chr", 3) == 0)
			{
				sprintf(filename, "%s/%s.inx", out_dir, name);
			}
			else
			{
				sprintf(filename, "%s/chr%s.inx", out_dir, name);
			}
			out_fp = fopen(filename, "w");
			if (!out_fp)
			{
				fprintf(stderr, "can't open index file %s for output\n", filename);
				exit(1);
			}
			tot_seqs = 0;
			nseqs = 0;
			pos = 0;
			n = 0;
			fprintf(stderr, "%s ... ", name);

			/* Header */
			fwrite(&pam, sizeof(pam_e), 1, out_fp);
		}
		else
		{
			int i;

			for(i = 0; i < len; i++)
			{
				char c = line[i];
				if (cmap[c] == 4)
				{
					u = 0;
					nu = 0;
				}
				else
				{
					u <<= 2;
					u |= cmap[c];
					nu++;
				}
				if (nu == max_size + pam_size)
				{
					/*
					fprintf(stderr, "%d %20.20s %ld\n", pos+i+1, bits_to_string(u, max_size), u);
					*/
					if (pam_left(u >> (2*max_size), pam))
					{
						/* fprintf(stderr, "L %d %26.26s %ld\n", pos+i-max_size-pam_size+2, bits_to_string(u, max_size+2*pam_size), u); */
						seqs[nseqs].pos = pos+i-max_size-pam_size+2;
						seqs[nseqs].seq = u;
						nseqs++;
						tot_seqs++;
						if (nseqs == MAX_SEQS)
						{
							fwrite(&seqs, sizeof(entry_t), MAX_SEQS, out_fp);
							nseqs = 0;
						}
					}
					if (pam_right(u, pam))
					{
						/* fprintf(stderr, "R %d %26.26s %ld\n", pos+i-max_size-pam_size+2, bits_to_string(u, max_size+2*pam_size), u); */
						seqs[nseqs].pos = pos+i-max_size-pam_size+2;
						seqs[nseqs].seq = u | PAM_RIGHT;
						nseqs++;
						tot_seqs++;
						if (nseqs == MAX_SEQS)
						{
							fwrite(&seqs, sizeof(entry_t), MAX_SEQS, out_fp);
							nseqs = 0;
						}
					}
					u &= all_mask;
					nu--;
				}
			}
			pos += len;
			n += len;
		}
	}
	if (nseqs > 0)
	{
		fwrite(&seqs, sizeof(entry_t), nseqs, out_fp);
	}
	fclose(out_fp);
	fprintf(stderr, "%d bases, %d sequences\n", n, tot_seqs);
}

void usage()
{
	fprintf(stderr,
"Usage: %s command [options]\n\
\n\
%s index [-out D] F.fa\n\
    -out D\n\
        Output directory for .inx files. Default current directory.\n\
    -pam (NGG|NGGRRT|...)\n\
        PAM sequence to index.\n\
\n\
%s query [options] F.inx ...\n\
    -seq Q\n\
        Query sequence - max 20bp - without PAM.\n\
    -list F\n\
        File of query sequences - one per line.\n\
    -count_only\n\
        Only count numbers of hits per query sequence.\n\
    -mismatch L N\n\
        Allow N mismatches for queries of length L.\n\
        May be repeated for different lengths. Default 0 mismatches.\n\
    -seed L N\n\
        Count hits with <= N mismatches in the seed region of length L as \"bad\".\n\
        May be repeated for different lengths. Default 8 and 2.\n\
    -n N\n\
        Use N cores to do the search. Default 1.\n\
    -format (table|gff|annotated|stats)\n\
        Output format. Default gff.\n\
    -info D\n\
        Direcory containing .info and .genes files.\n\
    -quiet\n\
        Don't output progress.\n\
\n\
%s region -chr C -start N -end N -index D [-info D]\n\
    -index D\n\
        Direcory containing .inx files.\n\
    -info D\n\
        Direcory containing .info and .genes files.\n\
    -chr C\n\
        Chromosome.\n\
    -start N\n\
        Start position.\n\
    -end N\n\
        End position.\n\
    -format (table|gff|annotated)\n\
        Output format. Default gff.\n\
", prog, prog, prog, prog);
	exit(1);
}

int index_main(int argc, char *argv[])
{
	char *out_dir = ".";
	int i, j;
	int nindexes = 0;

	pam_e pam = NGG;

	for(i = 2; i < argc; i++)
	{
		char *arg = argv[i];

		if (strcmp(arg, "-out") == 0 && i+1 < argc)
		{
			struct stat sb;
			out_dir = argv[++i];
			if (stat(out_dir, &sb) == 0 && S_ISDIR(sb.st_mode))
			{
				/* exists already */
			}
			else if (mkdir(out_dir, 0755))
			{
				fprintf(stderr, "%s: %s can't be created\n", prog, out_dir);
				perror(out_dir);
				exit(1);
			}
			if (access(out_dir, W_OK) != 0)
			{
				fprintf(stderr, "%s: %s is not writable\n", prog, out_dir);
				exit(1);
			}
		}
		else if (strcmp(arg, "-pam") == 0 && i+1 < argc)
		{
			pam = string_to_pam(argv[++i]);
		}
		else if (arg[0] == '-')
		{
			fprintf(stderr, "%s: bad index argument %s\n", prog, arg);
			usage();
		}
		else
		{
			FILE *fp = fopen(arg, "r");
			if (fp)
			{
				index_file(fp, out_dir, MAX_SIZE, pam);
				fclose(fp);
				nindexes++;
			}
			else
			{
				fprintf(stderr, "%s: can't open fasta file: %s\n", prog, arg);
				exit(1);
			}
		}
	}

	if (nindexes == 0)
	{
		usage();
	}
}

info_t *read_info(char *info_dir, char *index_file, char *chr)
{
	char info_file[MAX_BUFF];

	if (strncmp(chr, "chr", 3) == 0)
	{
		chr += 3;
	}

	if (info_dir)
	{
		sprintf(info_file, "%s/chr%s.info", info_dir, chr);
	}
	else
	{
		strcpy(info_file, index_file);
		strcpy(&info_file[strlen(info_file)-1], "fo");
	}

	info_t *info = NULL;
	FILE *info_fp = fopen(info_file, "r");
	if (info_fp)
	{
		info = read_feature_info(info_fp);
		fclose(info_fp);

		if (info_dir)
		{
			sprintf(info_file, "%s/chr%s.genes", info_dir, chr);
		}
		else
		{
			strcpy(info_file, index_file);
			strcpy(&info_file[strlen(info_file)-3], "genes");
		}

		info_fp = fopen(info_file, "r");
		if (info_fp)
		{
			read_gene_info(info_fp, info);
			fclose(info_fp);
		}
	}

	return (info);
}

int query_main(int argc, char *argv[])
{
	int i, j;
	char *info_dir = NULL;
	int quiet = FALSE, header_printed = FALSE;

	uint32_t mismatches[MAX_SIZE+1];
	int nthreads = 1;

	int nindexes = 0;
	char *indexes[MAX_INDEXES];
	uint32_t index_off_targets[MAX_INDEXES];

	for(i = 0; i <= MAX_SIZE; i++)
	{
		mismatches[i] = 0;
		max_mismatches[i] = 0;
	}
	bzero(query_mismatches, sizeof(query_mismatches));
	bzero(query_exact, sizeof(query_exact));
	bzero(bad_seed, sizeof(bad_seed));

	init_threads();

	for(i = 2; i < argc; i++)
	{
		char *arg = argv[i];

		if (strcmp(arg, "-count_only") == 0)
		{
			count_only = TRUE;
		}
		else if (strcmp(arg, "-quiet") == 0)
		{
			quiet = TRUE;
		}
		else if (strcmp(arg, "-n") == 0 && i+1 < argc)
		{
			int n = atoi(argv[++i]);
			if (n > 0 && n <= MAX_THREADS)
			{
				nthreads = n;
			}
		}
		else if (strcmp(arg, "-seed") == 0 && i+2 < argc)
		{
			int n = atoi(argv[++i]);
			int m = atoi(argv[++i]);
			if (n >= 0 && n <= MAX_SIZE)
			{
				if (m >= 0 && m <= n)
				{
					seed_length = n;
					max_seed_mismatches = m;
				}
				else
				{
					fprintf(stderr, "# %s: bad seed mismatches (%s), option ignored\n", prog, argv[i]);
				}
			}
			else
			{
				fprintf(stderr, "# %s: bad seed length (%s), option ignored\n", prog, argv[i-1]);
			}
		}
		else if (strcmp(arg, "-mismatch") == 0 && i+2 < argc)
		{
			int n = atoi(argv[++i]);
			int m = atoi(argv[++i]);
			if (n >= 0 && n <= MAX_SIZE && m >= 0 && m < MAX_SIZE)
			{
				if (n == 0)
				{
					int k;
					for(k = 1; k <= MAX_SIZE; k++)
					{
						max_mismatches[k] = m;
					}
				}
				else
				{
					max_mismatches[n] = m;
				}
			}
			else if (n > MAX_SIZE && m >=0 && m < n)
			{
				max_mismatches[MAX_SIZE] = m;
			}
		}
		else if (strcmp(arg, "-format") == 0 && i+1 < argc)
		{
			arg = argv[++i];
			if (strcmp(arg, "gff") == 0)
			{
				format = GFF_FORMAT;
			}
			else if (strcmp(arg, "table") == 0)
			{
				format = TABLE_FORMAT;
			}
			else if (strcmp(arg, "annotated") == 0)
			{
				format = ANNOTATED_FORMAT;
			}
			else if (strcmp(arg, "exact_annotated") == 0)
			{
				format = EXACT_ANNOTATED_FORMAT;
			}
			else if (strcmp(arg, "stats") == 0)
			{
				format = STATS_FORMAT;
			}
			else
			{
				fprintf(stderr, "# %s: %s is not a valid format\n", prog, arg);
				exit(1);
			}
		}
		else if (strcmp(arg, "-list") == 0 && i+1 < argc)
		{
			char *query_file = argv[++i];
			FILE *fp = fopen(query_file, "r");

			if (fp)
			{
				char *line;
				char buffer[MAX_BUFF];
				while (nqueries < MAX_QUERIES && (line = fgets(buffer, MAX_BUFF, fp)))
				{
					size_t len = strlen(line);
					line[--len] = '\0';
					if (len > MAX_SIZE)
					{
						line += len-MAX_SIZE;
						len = MAX_SIZE;
					}
					query_sizes[nqueries] = len;
					query_prob[nqueries][NONCODING] = 1.0;
					query_prob[nqueries][CODING] = 1.0;
					query_count[nqueries][NONCODING] = 0;
					query_count[nqueries][CODING] = 0;
					queries[nqueries++] = string_to_bits(line, len, 0);
				}
				fclose(fp);
				if (nqueries == MAX_QUERIES)
				{
					fprintf(stderr, "# %s: warning: maximum number of queries (%d) reached\n", prog, MAX_QUERIES);
				}
			}
			else
			{
				fprintf(stderr, "# %s: can't open query file: %s\n", prog, query_file);
				exit(1);
			}
		}
		else if (strcmp(arg, "-seq") == 0 && i+1 < argc)
		{
			char *query_seq = argv[++i];
			int len = strlen(query_seq);
			if (len > MAX_SIZE)
			{
				query_seq += len-MAX_SIZE;
				len = MAX_SIZE;
			}
			if (nqueries < MAX_QUERIES)
			{
				query_sizes[nqueries] = len;
				query_prob[nqueries][NONCODING] = 1.0;
				query_prob[nqueries][CODING] = 1.0;
				query_count[nqueries][NONCODING] = 0;
				query_count[nqueries][CODING] = 0;
				queries[nqueries++] = string_to_bits(query_seq, len, 0);
			}
		}
		else if (strcmp(arg, "-info") == 0 && i+1 < argc)
		{
			struct stat sb;
			info_dir = argv[++i];
			if (stat(info_dir, &sb) == 0 && S_ISDIR(sb.st_mode))
			{
				/* exists already */
			}
			else
			{
				fprintf(stderr, "# %s: %s doesn't exist or is not a directory\n", prog, info_dir);
				exit(1);
			}
			if (access(info_dir, R_OK) != 0)
			{
				fprintf(stderr, "# %s: %s is not readable\n", prog, info_dir);
				exit(1);
			}
		}
		else if (strcmp(&arg[ strlen(arg)-4], ".inx") == 0)
		{
			if (nqueries == 0)
			{
				fprintf(stderr, "# %s: no queries: %s skipped\n", prog, arg);
				continue;
			}
			if (nindexes >= MAX_INDEXES)
			{
				fprintf(stderr, "# %s: too many indexes: %s skipped\n", prog, arg);
				continue;
			}
			if ((format == ANNOTATED_FORMAT || format == EXACT_ANNOTATED_FORMAT) && !info_dir)
			{
				fprintf(stderr, "# %s: warning annotated format specified but no -info\n", prog);
			}
			if (!header_printed && !count_only)
			{
				if (format == TABLE_FORMAT)
				{
					printf("Chr\tStart\tEnd\tStrand\tHitSeq\tHitMismatches\tSeedMismatches\n");
					fflush(stdout);
				}
				else if (format == ANNOTATED_FORMAT || format == EXACT_ANNOTATED_FORMAT)
				{
					printf("QueryNo\tQuerySeq\tChr\tStart\tEnd\tStrand\tHitSeq\tHitMismatches\tSeedMismatches\tpOffTarget\tGeneID\tGeneSymbol\tGeneRegion\n");
					fflush(stdout);
				}
				header_printed = TRUE;
			}

			char *s, chr[MAX_BUFF];
			if (s = strrchr(arg, '/'))
			{
				strcpy(chr, s+1);
			}
			else
			{
				strcpy(chr, arg);
			}
			chr[ strlen(chr)-4 ] = '\0';

			int index_no = nindexes++;
			/* s = (strncmp(chr, "chr", 3) == 0 ? &chr[3] : chr);
			indexes[index_no] = strdup(s); */
			indexes[index_no] = strdup(chr);
			index_off_targets[index_no] = 0;

			FILE *fp = fopen(arg, "r");

			if (fp)
			{
				/* Header */
				pam_e pam;
				fread(&pam, sizeof(pam_e), 1, fp);

				if (!quiet)
				{
					fprintf(stderr, "# Searching %s (%s, PAM %s) ...\n", arg, indexes[index_no], pam_to_string(pam));
				}

				info_t *info = read_info(info_dir, arg, chr);

				search_index(index_no, fp, indexes[index_no],
							 mismatches, index_off_targets,
							 nthreads, pam, info);
			}
			else
			{
				fprintf(stderr, "# %s: can't open index file: %s\n", prog, arg);
				exit(1);
			}
		}
		else
		{
			fprintf(stderr, "# %s: bad query argument %s\n", prog, arg);
			exit(1);
		}
	}

	/*if (nqueries > 0)*/
	{
		if (nthreads > 1)
		{
			join_threads(nthreads, index_off_targets, mismatches);
		}

		FILE *out_fp = (count_only ? stdout : stderr);

		uint32_t tot_off_targets = 0;
		for(i = 0; i < nindexes; i++)
		{
			if (index_off_targets[i] > 0)
			{
				if (!quiet)
				{
					fprintf(out_fp, "# %s: %u off targets\n", indexes[i], index_off_targets[i]);
				}
				tot_off_targets += index_off_targets[i];
			}
		}

		int max_mm = 0;
		uint32_t tot_mismatches[MAX_SIZE+1];
		char buff[MAX_BUFF];
		for(j = 0; j <= MAX_SIZE; j++)
		{
			if (max_mismatches[j] > max_mm) max_mm = max_mismatches[j];
			tot_mismatches[j] = 0;
		}
		if (!quiet)
		{
			fprintf(out_fp, "# Total: %u off targets\n", tot_off_targets);
		}
		fprintf(out_fp, "QueryNo\tQuery\tQuerySize\tpCoding\tnCoding\tpNonCoding\tnNonCoding\tnBadSeed");
		for(j = 0; j <= max_mm; j++)
		{
			fprintf(out_fp, "\t%u", j);
		}
		fprintf(out_fp, "\n");
		for(i = 0; i < nqueries; i++)
		{
			fprintf(out_fp, "%d\t%s\t%d\t%g\t%d\t%g\t%d\t%d", i+1,
					bits_to_string(queries[i], 0ull, query_sizes[i], buff),
					query_sizes[i],
					1-query_prob[i][CODING],
					query_count[i][CODING],
					1-query_prob[i][NONCODING],
					query_count[i][NONCODING],
					bad_seed[i]);
			for(j = 0; j <= max_mm; j++)
			{
				fprintf(out_fp, "\t%u", query_mismatches[i][j]);
				tot_mismatches[j] += query_mismatches[i][j];
			}
			fprintf(out_fp, "\n");
		}
		if (!quiet)
		{
			fprintf(out_fp, "\tTotal\t\t\t\t\t\t");
			for(i = 0; i <= max_mm; i++)
			{
				fprintf(out_fp, "\t%u", mismatches[i]);
				if (mismatches[i] != tot_mismatches[i])
				{
					fprintf(out_fp, "/%u", tot_mismatches[i]);
				}
			}
			fprintf(out_fp, "\n");
		}
	}
	/*
	else
	{
		usage(prog);
	}
	*/
}

int region_main(int argc, char *argv[])
{
	int i, j;

	int nindexes = 0;
	char *indexes[MAX_INDEXES];
	int index_off_targets[MAX_INDEXES];

	char *chr = NULL;
	uint32_t start = 1, end = 1000000000;
	char *index_dir = NULL, *info_dir = NULL;

	format = GFF_FORMAT;

	for(i = 2; i < argc; i++)
	{
		char *arg = argv[i];

		if (strcmp(arg, "-index") == 0 && i+1 < argc)
		{
			struct stat sb;
			index_dir = argv[++i];
			if (stat(index_dir, &sb) == 0 && S_ISDIR(sb.st_mode))
			{
				/* exists already */
			}
			else
			{
				fprintf(stderr, "%s: %s doesn't exist or is not a directory\n", prog, index_dir);
				exit(1);
			}
			if (access(index_dir, R_OK) != 0)
			{
				fprintf(stderr, "%s: %s is not readable\n", prog, index_dir);
				exit(1);
			}
		}
		else if (strcmp(arg, "-info") == 0 && i+1 < argc)
		{
			struct stat sb;
			info_dir = argv[++i];
			if (stat(info_dir, &sb) == 0 && S_ISDIR(sb.st_mode))
			{
				/* exists already */
			}
			else
			{
				fprintf(stderr, "%s: %s doesn't exist or is not a directory\n", prog, info_dir);
				exit(1);
			}
			if (access(info_dir, R_OK) != 0)
			{
				fprintf(stderr, "%s: %s is not readable\n", prog, info_dir);
				exit(1);
			}
		}
		else if (strcmp(arg, "-format") == 0 && i+1 < argc)
		{
			arg = argv[++i];
			if (strcmp(arg, "gff") == 0)
			{
				format = GFF_FORMAT;
			}
			else if (strcmp(arg, "table") == 0)
			{
				format = TABLE_FORMAT;
				printf("Chr\tStart\tEnd\tStrand\tHitSeq\n");
				fflush(stdout);
			}
			else if (strcmp(arg, "annotated") == 0)
			{
				format = SHORT_ANNOTATED_FORMAT;
				printf("Chr\tStart\tEnd\tStrand\tHitSeq\tGeneID\tGeneSymbol\tGeneRegion\n");
				fflush(stdout);
			}
			else
			{
				fprintf(stderr, "%s: %s is not a valid format\n", prog, arg);
				exit(1);
			}
		}
		else if (strcmp(arg, "-chr") == 0 && i+1 < argc)
		{
			chr = argv[++i];
			if (strncmp(chr, "chr", 3) == 0)
			{
				chr += 3;
			}
		}
		else if (strcmp(arg, "-start") == 0 && i+1 < argc)
		{
			start = atoi(argv[++i]);
		}
		else if (strcmp(arg, "-end") == 0 && i+1 < argc)
		{
			end = atoi(argv[++i]);
		}
	}

	if (chr && start > 0 && end > start && index_dir)
	{
		if (format == SHORT_ANNOTATED_FORMAT && !info_dir)
		{
			fprintf(stderr, "# %s: warning annotated format specified but no -info\n", prog);
		}

		char index_file[MAX_BUFF];
		sprintf(index_file, "%s/chr%s.inx", index_dir, chr);

		FILE *fp = fopen(index_file, "r");

		if (fp)
		{
			info_t *info = read_info(info_dir, index_file, chr);

			pam_e pam;
			fread(&pam, sizeof(pam_e), 1, fp);
			int pam_size = pam_sizes[pam];
			side_e pam_side = pam_sides[pam];

			int nseqs = 0;
			while ((nseqs = fread(seqs, sizeof(entry_t), MAX_SEQS, fp)) > 0)
			{
				int first_pos = seqs[0].pos;
				int last_pos = seqs[nseqs-1].pos;

				if (start < last_pos && end > first_pos)
				{
					int j;
					for(j = 0; j < nseqs; j++)
					{
						uint32_t pos = seqs[j].pos;
						uint64_t seq = seqs[j].seq;
						char buff[MAX_BUFF];
						if (start <= pos+MAX_SIZE+pam_size && pos <= end)
						{
							char strand = (seq & PAM_RIGHT) ? '+' : '-';
							int cr_start = pos, cr_end = cr_start+MAX_SIZE+pam_size-1;

							if (format == GFF_FORMAT)
							{
								int pam_start, pam_end, c_start, c_end;
								if (strand == '+')
								{
									pam_start = cr_end-pam_size+1;
									pam_end = cr_end;
									c_start = cr_start;
									c_end = cr_end-pam_size;
								}
								else if (strand == '-')
								{
									pam_start = cr_start;
									pam_end = cr_start+pam_size-1;
									c_start = cr_start+pam_size;
									c_end = cr_end;
								}
								if (pam_side == LEFT)
								{
									FLIP(strand);
								}

								if (strand == '-')
								{
									seq = revcom(seq, MAX_SIZE+pam_size);
								}
							   	char *s = bits_to_string(seq, 0, MAX_SIZE+pam_size, buff);
								printf("chr%s\tWGE\tCrispr\t%d\t%d\t.\t%c\t.\tID=C%d;Name=%d;PAM=%s;Seq=%s;OT_Summary={}\n",
									   chr, cr_start, cr_end, strand, j, j, pam_to_string(pam), s);
								printf("chr%s\tWGE\tCDS\t%d\t%d\t.\t%c\t.\tID=Cr%d;Parent=C%d;Name=%d;color=#45A825\n",
									   chr, c_start, c_end, strand, j, j, j);
								printf("chr%s\tWGE\tCDS\t%d\t%d\t.\t%c\t.\tID=PAM%d;Parent=C%d;Name=%d;color=#1A8599\n",
									   chr, pam_start, pam_end, strand, j, j, j);
							}
							else if (format == TABLE_FORMAT)
							{
								if ((strand == '-' && pam_side == RIGHT) ||
									(strand == '+' && pam_side == LEFT))
								{
									/* Flip so the PAM is always forward */
									seq = revcom(seq, MAX_SIZE+pam_size);
								}
								if (pam_side == LEFT)
								{
									FLIP(strand);
								}
								printf("chr%s\t%d\t%d\t%c\t%s\n",
									   chr, cr_start, cr_end, strand,
							   		   bits_to_string(seq, 0, MAX_SIZE+pam_size, buff));
								fflush(stdout);
							}
							else if (format == SHORT_ANNOTATED_FORMAT)
							{
								uint32_t distPAMcut = pam_distance_to_cut[pam];
								uint32_t cut_pos = cr_start+MAX_SIZE-distPAMcut;
								if ((strand == '-' && pam_side == RIGHT) ||
									(strand == '+' && pam_side == LEFT))
								{
									/* Flip so the PAM is always forward */
									seq = revcom(seq, MAX_SIZE+pam_size);
									cut_pos = cr_start+pam_size+distPAMcut;
								}
								if (pam_side == LEFT)
								{
									FLIP(strand);
								}
							   	char *hit_seq = bits_to_string(seq, 0, MAX_SIZE+pam_size, buff);
								find_features(-1, "", chr, cr_start, cr_end, strand, hit_seq, -1, -1, -1.0, info, cut_pos);
								fflush(stdout);
							}
						}
					}
				}
				else if (first_pos > end)
				{
					/* printf("break %d %d - %d\n", nseqs, first_pos, last_pos); */
					break;
				}
			}
			fclose(fp);
		}
	}
	else
	{
		usage();
	}
}

int main(int argc, char *argv[])
{
	int i, j;

	prog = argv[0];

	populate_cmap();

	if (argc < 2)
	{
		usage();
	}
	else if (strcmp(argv[1], "index") == 0)
	{
		index_main(argc, argv);
	}
	else if (strcmp(argv[1], "query") == 0)
	{
		query_main(argc, argv);
	}
	else if (strcmp(argv[1], "region") == 0)
	{
		region_main(argc, argv);
	}
	else
	{
		fprintf(stderr, "# %s: unknown command: %s\n", prog, argv[1]);
		usage();
	}
	exit(0);
}
