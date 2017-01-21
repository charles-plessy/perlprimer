#!/usr/bin/perl -w

# PerlPrimer
# Designs primers for PCR, Bisulphite PCR, QPCR (Realtime), and Sequencing

# version 1.1.1 (15/3/2004)
# Copyright © 2003-2004, Owen Marshall

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 
# USA

use strict;
no strict 'refs';

#---------------#
# Usage message #
#---------------#

my ($version, $commandline);
BEGIN {
	$version = "1.1.1";
	($commandline) = @ARGV;

	if ($commandline && $commandline =~ /^-[\w-]+/) {
		print <<EOT;
PerlPrimer v$version
Designs primers for PCR, Bisulphite PCR, QPCR (Realtime), and Sequencing

Copyright © 2003-2004 Owen Marshall\n
Usage: perlprimer.pl [file.ppr]\n
EOT
		exit 0;
	}
}


#---------------#
# Load packages #
#---------------#

my $failed_packages;
BEGIN {
	# A modular package checking routine ...
	# ($failed_packages is used by the check_packages() subroutine)
	
	$failed_packages = " "; # check_packages gives errors if undefined
	my $warning;
	my $eval_sub = sub {
		eval $_;
		if ($@) {
			unless ($warning) {
				print "PerlPrimer v$version\nCopyright © 2003-2004 Owen Marshall\n\n";
				$warning = 1;
			}
			/(use|require)\s*([\w:-]+)\;/;
			my $package = $2;
			
			# Special circumstances!
			return if $package eq "Benchmark";
			if ($package eq "Tk") {
				print "Error: Perl/Tk not found!\nPerlPrimer requires Perl/Tk to run - please download from CPAN (http://cpan.org) and install before running.\n\n";
				exit 0;
			} elsif ($1 eq "require") {
				print "Error: $package not found!\nPerlPrimer requires $package to run - please download from CPAN (http://cpan.org) and install before running.\n\n";
				exit 0;
			}		
					
			print "Warning: PerlPrimer requires $package for full functionality\n- some features may be disabled.\n$package should be obtainable from CPAN (http://cpan.org)\n\n";
			$failed_packages .= " $package ";
		}
	};	
	
	foreach (split("\n","	
use Tk;
use Benchmark;
use HTTP::Request;
use LWP::UserAgent;
require Tk::NoteBook;
require Tk::HList;
require Tk::ItemStyle;
require Tk::Dialog;
require Tk::LabFrame;
require Tk::ROText;
require Tk::Balloon;
require Tk::BrowseEntry;
	")) {&$eval_sub($_)};
}


#------------------------------------------------------------------#
# Primer pair storage code for two dimensional array @primer_pair: #
#------------------------------------------------------------------#

# @primer_pair[index][see below]:
#
# 0. primer_f
# 1. pos_f
# 2. window_length_f
# 3. Tm_f
#
# 4. primer_r
# 5. pos_r
# 6. window_length_r
# 7. Tm_r
# 8. real_pos_r
#
# 9. amplicon_size
# 10. deltaG of most stable extensible primer-dimer
# 11. [primer_f sequence pre-conversion for Bisulphite PCR]
# 12. [primer_r sequence pre-conversion for Bisulphite PCR]


#-----
#
# NN thermodynamics hashes (AA = 5' AA 3'/3' TT 5') derived from ...
# 
# Allawi HT, SantaLucia J Jr.  Thermodynamics and NMR of internal G.T mismatches in DNA.
# 	Biochemistry. 1997 Aug 26;36(34):10581-94
#
# SantaLucia J Jr.  A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics.
# 	Proc Natl Acad Sci U S A. 1998 Feb 17;95(4):1460-5. 
# 
# ... with mismatch dG data (AGTG = 5' AG 3'/3' TG 5') derived from ...
# 
# Peyret N, Seneviratne PA, Allawi HT, SantaLucia J Jr.  Nearest-neighbor thermodynamics and NMR of DNA sequences with internal A.A, C.C, G.G, and T.T mismatches.
# 	Biochemistry. 1999 Mar 23;38(12):3468-77. 
# 
# Allawi HT, SantaLucia J Jr.  Nearest-neighbor thermodynamics of internal A.C mismatches in DNA: sequence dependence and pH effects.
# 	Biochemistry. 1998 Jun 30;37(26):9435-44.
# 
# Allawi HT, SantaLucia J Jr.  Thermodynamics of internal C.T mismatches in DNA.
# 	Nucleic Acids Res. 1998 Jun 1;26(11):2694-701. 
# 
# Allawi HT, SantaLucia J Jr.  Nearest neighbor thermodynamic parameters for internal G.A mismatches in DNA.
# 	Biochemistry. 1998 Feb 24;37(8):2170-9.
# 
# Allawi HT, SantaLucia J Jr.  Thermodynamics and NMR of internal G.T mismatches in DNA.
# 	Biochemistry. 1997 Aug 26;36(34):10581-94
# 
#-----

#-------------------#
# deltaH (kcal/mol) #
#-------------------#

my %oligo_dH=qw(
	AA -7.9 TT -7.9 
	AT -7.2 TA -7.2 
	CA -8.5 TG -8.5 
	GT -8.4 AC -8.4 
	CT -7.8 AG -7.8 
	GA -8.2 TC -8.2 
	CG -10.6 GC -9.8 
	GG -8.0 CC -8.0 
	initC 0.1 initG 0.1 
	initA 2.3 initT 2.3
);

my %oligo_dH_full=(
	qw(AATT -7.9 	TTAA -7.9 
	ATTA -7.2 	TAAT -7.2 
	CAGT -8.5 	TGAC -8.5 
	GTCA -8.4 	ACTG -8.4 
	CTGA -7.8 	AGTC -7.8 
	GACT -8.2 	TCAG -8.2 
	CGGC -10.6 	GCCG -9.8 
	GGCC -8.0 	CCGG -8.0
		
	initC 0.1 	initG 0.1 
	initA 2.3 	initT 2.3),
	
	# Like pair mismatches 
		
	qw(AATA 1.2 	ATAA 1.2
	CAGA -0.9 	AGAC -0.9
	GACA -2.9 	ACAG -2.9
	TAAA 4.7 	AAAT 4.7 
	
	ACTC 0.0 	CTCA 0.0 
	CCGC -1.5 	CGCC -1.5
	GCCC 3.6 	CCCG 3.6 
	TCAC 6.1 	CACT 6.1 
	
	AGTG -3.1 	GTGA -3.1
	CGGG -4.9 	GGGC -4.9
	GGCG -6.0 	GCGG -6.0
	TGAG 1.6 	GAGT 1.6 
	
	ATTT -2.7 	TTTA -2.7
	CTGT -5.0 	TGTC -5.0
	GTCT -2.2 	TCTG -2.2
	TTAT 0.2 	TATT 0.2  ),
	
	# G.T mismatches 
	
	qw(AGTT 1.0  	TTGA 1.0
	ATTG  -2.5 	GTTA  -2.5
	CGGT  -4.1 	TGGC  -4.1
	CTGG  -2.8 	GGTC  -2.8
	GGCT  3.3 	TCGG  3.3
	GGTT  5.8 	TTGG  5.8
	GTCG  -4.4 	GCTG  -4.4
	GTTG  4.1 	GTTG  4.1
	TGAT  -0.1 	TAGT  -0.1
	TGGT  -1.4 	TGGT  -1.4
	TTAG  -1.3 	GATT  -1.3), 
	
	# G.A mismatches 
	
	qw(AATG  -0.6 	GTAA  -0.6
	AGTA  -0.7 	ATGA  -0.7
	CAGG  -0.7 	GGAC  -0.7
	CGGA  -4.0 	AGGC  -4.0
	GACG  -0.6 	GCAG  -0.6
	GGCA  0.5 	ACGG  0.5
	TAAG  0.7 	GAAT  0.7
	TGAA  3.0 	AAGT  3.0), 
	
	# C.T mismatches 
	
	qw(ACTT  0.7 	TTCA  0.7
	ATTC  -1.2 	CTTA  -1.2
	CCGT  -0.8 	TGCC  -0.8
	CTGC  -1.5 	CGTC  -1.5
	GCCT  2.3 	TCCG  2.3 
	GTCC  5.2 	CCTG  5.2 
	TCAT  1.2 	TACT  1.2 
	TTAC  1.0 	CATT  1.0), 
	
	# A.C mismatches 
	
	qw(AATC  2.3	CTAA  2.3
	ACTA  5.3 	ATCA  5.3 
	CAGC  1.9 	CGAC  1.9 
	CCGA  0.6 	AGCC  0.6 
	GACC  5.2 	CCAG  5.2 
	GCCA  -0.7 	ACCG  -0.7
	TAAC  3.4  	CAAT  3.4 
	TCAA  7.6 	AACT  7.6),

);

#--------------------#
# deltaS (cal/K.mol) #
#--------------------#

my %oligo_dS=qw(
	AA -22.2 TT -22.2 
	AT -20.4 TA -21.3 
	CA -22.7 TG -22.7 
	GT -22.4 AC -22.4 
	CT -21.0 AG -21.0 
	GA -22.2 TC -22.2 
	CG -27.2 GC -24.4 
	GG -19.9 CC -19.9 
	initC -2.8 initG -2.8 
	initA 4.1 initT 4.1 
	sym -1.4
);

my %oligo_dS_full=(
	qw(AATT -22.2 	TTAA -22.2 
	ATTA -20.4 	TAAT -21.3 
	CAGT -22.7 	TGAC -22.7 
	GTCA -22.4 	ACTG -22.4 
	CTGA -21.0 	AGTC -21.0 
	GACT -22.2 	TCAG -22.2 
	CGGC -27.2 	GCCG -24.4 
	GGCC -19.9 	CCGG -19.9
		
	initC -2.8 	initG -2.8 
	initA 4.1 	initT 4.1
	sym -1.4),
	
	# Like pair mismatches
		
	qw(AATA 1.7 	ATAA 1.7
	CAGA -4.2 	AGAC -4.2 
	GACA -9.8 	ACAG -9.8 
	TAAA 12.9 	AAAT 12.9 
	
	ACTC -4.4 	CTCA -4.4 
	CCGC -7.2 	CGCC -7.2 
	GCCC 8.9 	CCCG 8.9 
	TCAC 16.4 	CACT 16.4 
	
	AGTG -9.5 	GTGA -9.5 
	CGGG -15.3 	GGGC -15.3
	GGCG -15.8 	GCGG -15.8
	TGAG 3.6 	GAGT 3.6 
	
	ATTT -10.8 	TTTA -10.8
	CTGT -15.8 	TGTC -15.8
	GTCT -8.4 	TCTG -8.4 
	TTAT -1.5 	TATT -1.5),
	
	# G.T mismatches
	
	qw(AGTT 0.9 	TTGA 0.9
	ATTG  -8.3 	GTTA  -8.3
	CGGT  -11.7 	TGGC  -11.7
	CTGG  -8.0 	GGTC  -8.0
	GGCT  10.4 	TCGG  10.4
	GGTT  16.3 	TTGG  16.3
	GTCG  -12.3 	GCTG  -12.3
	GTTG  9.5 	GTTG  9.5
	TGAT  -1.7 	TAGT  -1.7
	TGGT  -6.2 	TGGT  -6.2
	TTAG  -5.3 	GATT  -5.3), 
	
	# G.A mismatches
	
	qw(AATG  -2.3 	GTAA  -2.3
	AGTA  -2.3 	ATGA  -2.3
	CAGG  -2.3 	GGAC  -2.3
	CGGA  -13.2 	AGGC  -13.2
	GACG  -1.0 	GCAG  -1.0
	GGCA  3.2 	ACGG  3.2
	TAAG  0.7 	GAAT  0.7
	TGAA  7.4 	AAGT  7.4), 
	
	# C.T mismatches
	
	qw(ACTT  0.2 	TTCA  0.2
	ATTC  -6.2 	CTTA  -6.2
	CCGT  -4.5 	TGCC  -4.5
	CTGC  -6.1 	CGTC  -6.1
	GCCT  5.4 	TCCG  5.4 
	GTCC  13.5 	CCTG  13.5
	TCAT  0.7 	TACT  0.7 
	TTAC  0.7 	CATT  0.7), 
	
	# A.C mismatches
	
	qw(AATC  4.6 	CTAA  4.6
	ACTA  14.6 	ATCA  14.6
	CAGC  3.7 	CGAC  3.7 
	CCGA  -0.6 	AGCC  -0.6
	GACC  14.2 	CCAG  14.2
	GCCA  -3.8 	ACCG  -3.8
	TAAC  8.0  	CAAT  8.0 
	TCAA  20.2 	AACT  20.2),

);


# Genetic code hash
my %genetic_code=qw(
		TTT F TTC F TTA L TTG L
		CTT L CTC L CTA L CTG L
		ATT I ATC I ATA I ATG M
		GTT V GTC V GTA V GTG V
		TCT S TCC S TCA S TCG S
		CCT P CCC P CCA P CCG P
		ACT T ACC T ACA T ACG T
		GCT A GCC A GCA A GCG A
		TAT Y TAC Y TAA * TAG *
		CAT H CAC H CAA Q CAG Q
		AAT N AAC N AAA K AAG K
		GAT D GAC D GAA E GAG E
		TGT C TGC C TGA * TGG W
		CGT R CGC R CGA R CGG R
		AGT S AGC S AGA R AGG R
		GGT G GGC G GGA G GGG G
);

		
# Starting ionic concentration variables
my $oligo_conc = 200; #in nM
my $mg_conc=1.5; #in mM
my $monovalent_cation_conc=50; #in mM
my $dntp_conc=0.2; #in mM


# Default primer variables
my $max_tm_pr=63; # standard PCR
my $min_tm_pr=57;
my $max_diff_pr=3;
my $pri_win_min_pr=20;
my $pri_win_max_pr=24;
my $exclude_gc=1;
my $exclude_clamp=1;

my $max_tm_bs=65; # bisulphite PCR
my $min_tm_bs=55;
my $max_diff_bs=5;
my $pri_win_min_bs=25;
my $pri_win_max_bs=30;
my $bisul_min_c=30;
my $pre_bs=1;
my $exclude_3c=1;
my $exclude_cpg=0;

my $max_tm_q=62; # QPCR
my $min_tm_q=58;
my $max_diff_q=2;
my $pri_win_min_q=20;
my $pri_win_max_q=24;
my $min_ampsize_q=100;
my $max_ampsize_q=300;
my $exclude_ie=7;

my $min_tm_seq = 56; # Sequencing
my $max_tm_seq = 64;
my $pri_win_min_seq = 20;
my $pri_win_max_seq = 24;
my $seq_spacing_min = 500;
my $seq_spacing_max = 700;
my $exclude_gc_seq = 1;
my $exclude_clamp_seq = 1;
my $exclude_pd_seq = 1;
my $seq_pd_min = 5;

# Establish default varibles for opening a new file - perhaps a bit verbose
# but it is rather neat ...
my %default_variables = (
	pd => {
		\$min_tm_pr => $min_tm_pr,
		\$max_tm_pr => $max_tm_pr,
		\$max_diff_pr => $max_diff_pr,
		\$pri_win_min_pr => $pri_win_min_pr,
		\$pri_win_max_pr => $pri_win_max_pr,
		\$exclude_gc => $exclude_gc,
		\$exclude_clamp => $exclude_clamp,
	},
	seq => {
		\$min_tm_seq => $min_tm_seq,
		\$max_tm_seq => $max_tm_seq,             
		\$pri_win_min_seq => $pri_win_min_seq,        
		\$pri_win_max_seq => $pri_win_max_seq,        
		\$seq_spacing_min => $seq_spacing_min,
		\$seq_spacing_max => $seq_spacing_max,
		\$exclude_gc_seq => $exclude_gc_seq,
		\$exclude_clamp_seq => $exclude_clamp_seq,
		\$exclude_pd_seq => $exclude_pd_seq,
		\$seq_pd_min => $seq_pd_min,  
	},
	bis => {
		\$min_tm_bs => $min_tm_bs,
		\$max_tm_bs => $max_tm_bs,
		\$max_diff_bs => $max_diff_bs,
		\$pri_win_min_bs => $pri_win_min_bs,
		\$pri_win_max_bs => $pri_win_max_bs,
		\$exclude_cpg => $exclude_cpg,
		\$pre_bs => $pre_bs,
		\$exclude_3c => $exclude_3c,
		\$bisul_min_c => $bisul_min_c,
	},
	qpcr => {
		\$min_tm_q => $min_tm_q,
		\$max_tm_q, => $max_tm_q,
		\$max_diff_q => $max_diff_q,
		\$pri_win_min_q => $pri_win_min_q,
		\$pri_win_max_q => $pri_win_max_q,
		\$min_ampsize_q => $min_ampsize_q,
		\$max_ampsize_q => $max_ampsize_q,
		\$exclude_gc => $exclude_gc,
		\$exclude_clamp => $exclude_clamp,
		\$exclude_ie => $exclude_ie,
	},
);

# flags
my $bs=0;
my $qpcr_flag=0;
my $cancel=0;
my $defer_to_caps=0;

# benchmarking - for code optimisation
my $benchmark=0;

# DNA canvas sizing:
my $dna_canvas_offset=15;
my $dna_canvas_height=21;
my $dna_canvas_middle=int($dna_canvas_height/2)+2;

my $half_dna_size=3;
my $dc_dna_y1 = $dna_canvas_middle - $half_dna_size;
my $dc_dna_y2 = $dna_canvas_middle + $half_dna_size;

my $dc_selection_offset = 3;
my $dc_sel_y1 = $dc_dna_y1-$dc_selection_offset;
my $dc_sel_y2 = $dc_dna_y2+$dc_selection_offset;
my $dc_sel_offset2 = 2;

# Repeats/runs
my $exclude_rr = 1;
my $repeat = 3;
my $run = 4;

my $exclude_rr_bs = 1;
my $repeat_bs = 3;
my $run_bs = 4;

# %GC exclusions
my $max_gc = 60;
my $min_gc = 40;

# CpG prediction variables
my $cpg_window = 200;
my $min_cpg_island = 200;
my $cpg_gc = 50;
my $cpg_oe = 0.6;
my $cpgplot_method = 0;

# Http proxy
my $use_proxy = 0;
my ($http_proxy, $http_port);

# Ensembl data
my $ensembl; # GUI page
my $ensembl_gene;
my $ensembl_organism = 'Homo_sapiens';
my $ensembl_type = 'cdna';
my @ensembl_species = split("\n",
"Homo_sapiens
Mus_musculus
Rattus_norvegicus
Danio_rerio
Fugu_rubripes
Caenorhabditis_elegans
Drosophila_melanogaster
Anopheles_gambiae");

my @ensembl_types = split("\n",
"genomic
cdna
coding
utr5
utr3");

# BLAST search parameters
my $blast_expect = 10;
my $blast_database = 'nr';
my @blast_database_array = split("\n",
"nr
est
est_human
est_mouse
est_others
gss
htgs
pat
yeast
mito
vector
pdb
month
alu
dbsts
chromosome
E. coli
Drosophila genome");
 
my $blast_entrez_query = 'none';
my @blast_entrez_array = split("\n",
"none
-------------------
Viruses [ORGN]
Archaea [ORGN]
Bacteria [ORGN]
Eukaryota [ORGN]
Viridiplantae [ORGN]
Fungi [ORGN]
Metazoa [ORGN]
Arthropoda [ORGN]
Vertebrata [ORGN]
Mammalia [ORGN]
Rodentia [ORGN]
Primates [ORGN]
-------------------
Aeropyrum pernix [ORGN]
Aquifex aeolicus [ORGN]
Arabidopsis thaliana [ORGN]
Bacillus subtilis [ORGN]
Bos taurus [ORGN]
Caenorhabditis elegans [ORGN]
Danio rerio [ORGN]
Dictyostelium discoideum [ORGN]
Drosophila melanogaster [ORGN]
Escherichia coli [ORGN]
Gallus gallus [ORGN]
Homo sapiens [ORGN]
Human immunodeficiency virus type 1 [ORGN]
Methanococcus jannaschii [ORGN]
Mus musculus [ORGN]
Oryctolagus cuniculus [ORGN]
Oryza sativa [ORGN]
Ovis aries [ORGN]
Plasmodium falciparum [ORGN]
Rattus norvegicus [ORGN]
Saccharomyces cerevisiae [ORGN]
Schizosaccharomyces pombe [ORGN]
Simian immunodeficiency virus [ORGN]
Synechocystis PCC6803 [ORGN]
Takifugu rubripes [ORGN]
Xenopus laevis [ORGN]
Zea mays [ORGN]");

my $blast_count;
my $http_abort;

# Primer-dimer parameters
my $pd_full=0;
my $pd_extensible=1;
my $pd_temperature=37;

# Restriction enzyme cloning sequence parameters
my $cloning_anchor = "GCGCGC";
my $simple_sites = 1;
my $exclude_found_sites = 1;

# Mouse wheel
my $scroll_factor = 2;

# Home environment variables and OS-specific tweaking:
# (HOME defaults to home folder on *nix and C drive on Win32)
### Need a tmp folder search here for spidey files to be placed into!!

# These variables are the total changes required for nice cross-platform compatibility 
# (and most of these are cosmetic - strange that the Tk look&feel is not consistent with
# things such as widget spacing, most notably the checkbuttons)
my ($HOME, $os, $gui_font_face, $gui_font_size, $font_override, $menu_relief, $frame_pady, $check_pady, $button_pady, $button_pack_padx, $button_pack_pady);
unless ($^O =~ /win/i) {
	$HOME = (getpwuid($<))[7].'/';
	$gui_font_size = 10;
	$os = 'nix';
	$gui_font_face = "Helvetica";
	$menu_relief = 'raised';
	$button_pady = 2;
	$button_pack_padx = 0;
	$button_pack_pady = 0;
	$frame_pady = 1;
	$check_pady = 3;
} else {
	$HOME = 'c:\\';
	$os = 'win';
	$gui_font_size = 8;
	$gui_font_face = "MS Sans Serif";
	$menu_relief = 'flat';
	$button_pady = -2;
	$button_pack_padx = 2;
	$button_pack_pady = 2;
	$frame_pady = 2;
	$check_pady = 0;
}

# Win32 only - for button pre-lights
my $activebackground_color="#ececec";

# Spidey details
my (%file_pre, %spidey_exec);
$file_pre{win} = "$HOME";
$file_pre{nix} = "$HOME";
$spidey_exec{win} = 'Spidey.exe';
$spidey_exec{nix} = 'spidey.linux';

my $file_pre = $file_pre{$os};
my $spidey_exec = $file_pre.$spidey_exec{$os};		
my $spidey_command = $spidey_exec." -i ".$file_pre.".dna_tmp -m ".$file_pre.".mrna_tmp ";
my $spidey_out;

# More global variables
my (
	$primer_f, $primer_r, $pos,
	$rprimer_r,
	$pd, @score_sort, $dna_canvas_size, $reverse,
	$pfkeys, $pkeys, @PF, @PR, %primer_hash,
	$gc_exclude, $gc_percent_ex, @primer_pairs,
	$prl, $pfl, $pl, @bind_string, %rating_hash, @score,
	
	$save_seq, $save_seq2,
	@primer_pairs_pr_s, @primer_pairs_seq_s, @primer_pairs_bs_s, @primer_pairs_q_s, @save_selection,
	@intron_exon_bounds,
	
	$primer_seq_5f, $primer_seq_5r, $primer_seq_5f_frame, $primer_seq_5r_frame,
	$primer_seq_5f_atg, $primer_seq_5r_atg,
);

my (
	$min_ampsize_pr, $max_ampsize_pr, $max_range_5p_pr, $min_range_pr, $max_range_pr, $max_range_3p_pr,
	$min_range_seq, $max_range_seq,	
	$min_ampsize_bs, $max_ampsize_bs, $max_range_5p_bs, $min_range_bs, $max_range_bs, $max_range_3p_bs,
	$max_range_5p_q, $min_range_q, $max_range_q, $max_range_3p_q,
);

# GUI globals
my (
	%prr, %prc, %pre, %prb, %prl, %prf,
	$fprimer, $fprimer_tm, $fprimer_len, $fprimer_ds, $fprimer_dh, $fprimer_dg,
	$rprimer, $rprimer_tm, $rprimer_len, $rprimer_ds, $rprimer_dh, $rprimer_dg,
	$old_reference, $text_widget_ref,
	$style_primer, $style_tm,
	%min_amp_canvas, %max_amp_canvas, %min_range_canvas, %max_range_canvas,
	$flag, $rid_get, @blast_results_1, @blast_results_2, @blast_results_3, $rptid, $blast_status,
	$popup_sort, $popup_sort_seq,
	$stored_page,
	$forward_re_site, $reverse_re_site,
);

# GUI dialogues
my (
	$blast_d, $view_base, $prefs, $ack_d, $canvas_info_d, $info_d,
	$cloning_d,
);

my $balloon_help = 1;		

# Recently used files
my @mru;
my $mru_number = 8;

		
# Define the perlprimer recognised filetypes
my $file_types = [
	['PerlPrimer Files', '.ppr'],
	['Fasta Files', '.fasta'],
	['All Files', '*']
	];

# Open file hash (for writing reports)
my %open_file;
for my $i (qw(pd seq bis qpcr)) {
	$open_file{$i} = "File not saved"
}

# Variable/array hashes for opening and saving
my %variables = (
	pd => {
		mintm => \$min_tm_pr,
		maxtm => \$max_tm_pr,
		maxdiff => \$max_diff_pr,
		minwin => \$pri_win_min_pr,
		maxwin => \$pri_win_max_pr,
		minamp => \$min_ampsize_pr,
		maxamp => \$max_ampsize_pr,
		maxrange5p => \$max_range_5p_pr,
		minrange => \$min_range_pr,
		maxrange => \$max_range_pr,
		maxrange3p => \$max_range_3p_pr,
		exclude_gc => \$exclude_gc,
		gc_clamp => \$exclude_clamp,
		seq => \$save_seq,
		primer_seq_5f => \$primer_seq_5f,
		primer_seq_5f_frame => \$primer_seq_5f_frame,
		primer_seq_5r => \$primer_seq_5r,
		primer_seq_5r_frame => \$primer_seq_5r_frame,
	},
	seq => {
		mintm => \$min_tm_seq,
		maxtm => \$max_tm_seq,
		minwin => \$pri_win_min_seq,
		maxwin => \$pri_win_max_seq,
		minrange => \$min_range_seq,
		maxrange => \$max_range_seq,
		seq_spacing_min => \$seq_spacing_min,
		seq_spacing_max => \$seq_spacing_max,
		exclude_gc => \$exclude_gc_seq,
		gc_clamp => \$exclude_clamp_seq,
		exclude_pd => \$exclude_pd_seq,
		seq_pd_min => \$seq_pd_min,
		seq => \$save_seq,
	},
	bis => {
		mintm => \$min_tm_bs,
		maxtm => \$max_tm_bs,
		maxdiff => \$max_diff_bs,
		minwin => \$pri_win_min_bs,
		maxwin => \$pri_win_max_bs,
		minamp => \$min_ampsize_bs,
		maxamp => \$max_ampsize_bs,
		maxrange5p => \$max_range_5p_bs,
		minrange => \$min_range_bs,
		maxrange => \$max_range_bs,
		maxrange3p => \$max_range_3p_bs,
		exclude_cpg => \$exclude_cpg,
		rrpostconv => \$pre_bs,
		exclude_cs => \$exclude_3c,
		min_c => \$bisul_min_c,
		seq => \$save_seq,
	},
	qpcr => {
		mintm => \$min_tm_q,
		maxtm => \$max_tm_q,
		maxdiff => \$max_diff_q,
		minwin => \$pri_win_min_q,
		maxwin => \$pri_win_max_q,
		minamp => \$min_ampsize_q,
		maxamp => \$max_ampsize_q,
		exclude_gc => \$exclude_gc,
		gc_clamp => \$exclude_clamp,
		mrna_seq => \$save_seq,
		dna_seq => \$save_seq2,
		exclude_ie => \$exclude_ie,
	},
);

my %arrays = (
	pd => {
		res => \@primer_pairs_pr_s,
		selection => \@save_selection,
	},
	seq => {
		res => \@primer_pairs_seq_s,
		selection => \@save_selection,
	},
	bis => {
		res => \@primer_pairs_bs_s,
		selection => \@save_selection,
	},
	qpcr => {
		res => \@primer_pairs_q_s,
		selection => \@save_selection,
		intron_exon => \@intron_exon_bounds,
	},
);
		
# Hash for opening/saving prefs
my %pref_variables = (
	home => \$HOME,
	spidey_path => \$file_pre{$os},
	repeats => \$repeat,
	runs => \$run,
	exclude_rr => \$exclude_rr,
	repeats_bs => \$repeat_bs,
	runs_bs => \$run_bs,
	exclude_rr_bs => \$exclude_rr_bs,
	max_gc => \$max_gc,
	min_gc => \$min_gc,
	mg_conc => \$mg_conc,
	monovalent_cation_conc => \$monovalent_cation_conc,
	dNTP_conc => \$dntp_conc,
	oligo_conc => \$oligo_conc,
	cpg_window => \$cpg_window,
	min_cpg_island => \$min_cpg_island,
	cpg_gc_percent => \$cpg_gc,
	cpg_observed_expected => \$cpg_oe,
	balloon_help => \$balloon_help,
	blast_database => \$blast_database,
	blast_entrez_query => \$blast_entrez_query,
	blast_expect => \$blast_expect,
	gui_font_face => \$gui_font_face,
	gui_font_size => \$gui_font_size,
	font_override => \$font_override,
	scroll_factor => \$scroll_factor,
	ensembl_organism => \$ensembl_organism,
	ensembl_type => \$ensembl_type,
	use_proxy => \$use_proxy,
	http_proxy => \$http_proxy,
	http_proxy_port => \$http_port,
	cpgplot_method => \$cpgplot_method,
	defer_to_caps => \$defer_to_caps,
	mru_number => \$mru_number,
	re_simple_sites => \$simple_sites,
	cloning_anchor_sequence => \$cloning_anchor,
	pd_temperature => \$pd_temperature,
	exclude_found_sites => \$exclude_found_sites,
);

my %pref_arrays = (
	mru => \@mru,
);


# Hash for variable references
my $null = undef;
my %page_specific_vars = (
	pd => {
		min_tm => \$min_tm_pr,
		max_tm => \$max_tm_pr,
		max_diff => \$max_diff_pr,
		pri_win_win => \$pri_win_min_pr,
		pri_win_max => \$pri_win_max_pr,
		min_ampsize => \$min_ampsize_pr,
		max_ampsize => \$max_ampsize_pr,
		max_range_5p => \$max_range_5p_pr,
		min_range => \$min_range_pr,
		max_range => \$max_range_pr,
		max_range_3p => \$max_range_3p_pr,
		seq => \$prf{seq},
		primers => \@primer_pairs_pr_s,
		hlist => \$prf{res},
		canvas => \$prf{primer_canvas},
		subroutine => \&get_gene,
		find_sub => \&find_orf,
		primer_sub => \&get_primers,
		popup => \$popup_sort,
	},
	seq => {
		min_tm => \$min_tm_seq,
		max_tm => \$max_tm_seq,
		pri_win_win => \$pri_win_min_seq,
		pri_win_max => \$pri_win_max_seq,
		min_ampsize => \$null,
		max_ampsize => \$null,
		max_range_5p => \$null,
		min_range => \$min_range_seq,
		max_range => \$max_range_seq,
		max_range_3p => \$null,
		seq => \$prf{seq_seq},
		primers => \@primer_pairs_seq_s,
		hlist => \$prf{seq_res},
		canvas => \$prf{seq_canvas},
		subroutine => \&get_gene,
		find_sub => \&find_orf,
		primer_sub => \&get_seq_primers,
		popup => \$popup_sort_seq,
	},
	bis => {
		min_tm => \$min_tm_bs,
		max_tm => \$max_tm_bs,
		max_diff => \$max_diff_bs,
		pri_win_win => \$pri_win_min_bs,
		pri_win_max => \$pri_win_max_bs,
		min_ampsize => \$min_ampsize_bs,
		max_ampsize => \$max_ampsize_bs,
		max_range_5p => \$max_range_5p_bs,
		min_range => \$min_range_bs,
		max_range => \$max_range_bs,
		max_range_3p => \$max_range_3p_bs,
		seq => \$prf{bisul_seq},
		primers => \@primer_pairs_bs_s,
		hlist => \$prf{bisul_res},
		canvas => \$prf{bisul_canvas},
		subroutine => \&get_cpg,
		find_sub => \&find_cpg,
		primer_sub => \&get_bisulphite,
		popup => \$popup_sort,
	},
	qpcr => {
		min_tm => \$min_tm_q,
		max_tm => \$max_tm_q,
		max_diff => \$max_diff_q,
		pri_win_win => \$pri_win_min_q,
		pri_win_max => \$pri_win_max_q,
		min_ampsize => \$min_ampsize_q,
		max_ampsize => \$max_ampsize_q,
		max_range_5p => \$null,
		min_range => \$null,
		max_range => \$null,
		max_range_3p => \$null,
		seq => \$prf{qmrna_seq},
		primers => \@primer_pairs_q_s,
		hlist => \$prf{qpcr_res},
		canvas => \$prf{qprimer_canvas},
		find_sub => \&find_orf,
		popup => \$popup_sort,
	},
);

# prefs_file
my $pref_file = $HOME.'.perlprimer';
read_prefs();

# Balloon help messages
my %balloonmsg = (
	'primer_getgene', "Find the longest ORF within the sequence and set the selected range",
	'primer_reset', "Reset the range to the ORF boundaries",
	'primer_stepin', "Reduce the inner range by 10bp on each side",
	'primer_stepout', "Increase the outer range by 10bp on each side",
	'exclude_gc', "Exclude primers with greater than 60% GC content\nor less than 40% GC content",
	'gc_clamp', "Require primers to have a 3' GC clamp\n(two of the last three bases G or C)",
	'primer_seq_5f', "Sequence to add to the 5' end of the forward primer\nPlace an underscore (_) before restriction enzyme site\n\nUse the 'Add cloning sequences' menu option to configure this automatically",
	'primer_seq_5r', "Sequence to add to the 5' end of the reverse primer\nPlace an underscore (_) before restriction enzyme site\n\nUse the 'Add cloning sequences' menu option to configure this automatically",
	'primer_seq_5f_frame', "Frame of the restriction enzyme site in the cloning plasmid",
	'primer_seq_5r_frame', "Frame of the restriction enzyme site in the cloning plasmid",
	'seq', "DNA sequence: right-click for menu with options\nincluding opening and saving files",
	'res', "Matching primer pairs are displayed here:\nDouble click to see primer-dimers\nRight-click for option menu",
	'primerbutton', "Find primer pairs",
	'autobuttonin', "Successively reduce the inner range by 10bp increments\nuntil primer pairs are found",
	'autobuttonout', "Successively increase the outer range by 10bp increments\nuntil primer pairs are found",
	'primer_cancel', "Cancel the current task",
	'primer_view', "Copy the selected primer pairs to the clipboard\n(format is tab-delimited text)",
	'primer_canvas', "Primer canvas",
	
	'bisul_getcpg', "Find CpG islands within the sequence and set the selected range",
	'bisul_reset', "Reset the range to the ORF boundaries",
	'bisul_stepin', "Reduce the inner range by 10bp on each side",
	'bisul_stepout', "Increase the outer range by 10bp on each side",
	'bisul_seq', "DNA sequence: right-click for menu with options\nincluding opening and saving files",
	'bisul_res', "Matching primer pairs are displayed here:\nDouble click to see primer-dimers\nRight-click for option menu",
	'bisul_exclude_cpg', "Exclude primers with CpG residues",
	'bisul_pre_bs', "Exclude primers with repeats / runs after bisulphite conversion",
	'bisul_exclude_cs', "Exclude primers without 3' C content\n(either last base as C, or two of last three bases C)",
	'bisul_min_c', "Minimum \%C content for primers",
	'bisul_button', "Find primer pairs",
	'bisul_autobuttonin', "Successively reduce the inner range by 10bp increments\nuntil primer pairs are found",
	'bisul_autobuttonout', "Successively increase the outer range by 10bp increments\nuntil primer pairs are found",
	'bisul_cancel', "Cancel the current task",
	'bisul_view', "Copy the selected primer pairs to the clipboard\n(format is tab-delimited text)",
	
	'qexclude_gc', "Exclude primers with greater than 60% GC content\nor less than 40% GC content",
	'qgc_clamp', "Require primers to have a 3' GC clamp\n(two of the last three bases G or C)",
	'qdna_seq', "Genomic DNA sequence: right-click for menu with options\nincluding opening and saving files",
	'qmrna_seq', "mRNA sequence: right-click for menu with options\nincluding opening and saving files",
	'qpcr_res', "Matching primer pairs are displayed here:\nDouble click to see primer-dimers\nRight-click for option menu",
	'qprimerbutton', "Find primer pairs\n(at least one primer must span an intron-exon boundary)",
	'qprimer_cancel', "Cancel the current task",
	'qprimer_view', "Copy the selected primer pairs to the clipboard\n(format is tab-delimited text)",
	'qprimer_spidey', "Detailed intron/exon boundary information",
	'qexclude_ie', "Minimum number of bases that an intron/exon boundary\nmust lie within at least one primer of a primer pair",
	
	'sspacingmin', "Minimum distance between sequencing primers",
	'sspacingmax', "Maximum distance between sequencing primers",
	'seqminrange', "The region of the DNA to be sequenced",
	'seqmaxrange', "The region of the DNA to be sequenced",
	'seq_getgene', "Set the sequencing region to the boundaries of the largest ORF present",
	'exclude_gc_seq', "Exclude primers with greater than 60% GC content\nor less than 40% GC content",
	'gc_clamp_seq', "Require primers to have a 3' GC clamp\n(two of the last three bases G or C)",
	'exclude_pd_seq', "Do not consider primers which can form dimers\nwith a stability greater than this value",
	'spdmin', "Do not consider primers which can form dimers\nwith a stability greater than this value",
	'seq_button', "Find primer pairs\n(at least one primer must span an intron-exon boundary)",
	'seq_cancel', "Cancel the current task",
	'seq_view', "Copy the selected primer pairs to the clipboard\n(format is tab-delimited text)",
	
	'tmbutton', "Calculate Tm, thermodynamic properties and primer-dimers from the primer sequences",
	'blastbutton', "Perform a BLAST search on the primers (requires an internet connection)",
	
	'prefs_defer', "Do not attempt to find an ORF or CpG Islands if\na captalised region is already present in the\nDNA sequence - use that region to set the range instead",
	'prefs_cpgplot_method', "Emulate the behaviour of the program cpgplot when\nfinding CpG islands (please note that this will cause the \%GC and O/E\nto be an overestimate - see the documentation for details)",
	'prefs_gui_override', "Use system default fonts for most widgets\n(requires program restart)",
	'prefs_gui_family', "Changes will not take effect until program restart",
	'prefs_simple_sites', "Limit restriction enzyme database to 6-base cutting enzymes\n(Also excludes degenerate cutting enzymes)",
	'prefs_exclude_found_sites', "When displaying the list of restriction enzymes,\nonly list those that will not cut the input sequence\n(Highly recommended!)",
	
	'cloning_anchor', "Sequence to add 5' of the restiction enzyme sites\n(recommended for successful digestion after amplification)",
	
	'blast_f', "Display results for the forward primer",
	'blast_r', "Display results for the reverse primer",
	'blast_search_string', "Enter a string to limit results",
	'blast_search', "Limit results",
	
	'view_base_copy', "Copies the layout to the clipboard\nas 80-column wrapped text",
	'view_base_pf', "Jump to the forward primer",
	'view_base_pr', "Jump to the reverse primer", 
);


#--------------------#
# deltaG (kcal/mol)  #
#--------------------#

# Recalculate oligo_dG (kcal/mol) for PCR salt conditions
# (the initiation values are salt independent)
my %oligo_dG=(
	qw(initC 0.98 	initG 0.98 
	initA 1.03 	initT 1.03), 
);
recalculate_dG();


# Load pixmap icon data
my (
	$perlprimer_icon, $icon_open, $icon_save, $icon_clear,
	$icon_info, $icon_magnify, $icon_new, $dna_canvas_pixmap,
	$icon_separator, $icon_prefs, $icon_report, $icon_save_as,
	$icon_ensembl, $icon_dna_open, $icon_dna_save,
);
load_icon_data();


		#-----------------#
		#  		  #
		#  GUI Interface  #
		#  		  #
		#-----------------#

#-------------#
# Main Window #
#-------------#

my $top = MainWindow->new(-title=>"PerlPrimer v$version");
$top->withdraw();

#--------------#
# Balloon help #
#--------------#

# Set the balloon background colour (nobody actually *likes* the default sick yellow
# colour, do they??) and the font if we're using Win32
$top->optionAdd("*Balloon.Background", "#ececec");
		
# configure fonts
my $gui_font = "{$gui_font_face} $gui_font_size";
my $font = $gui_font;
my $gui_font_bold = $gui_font . ' bold';
my $text_font = "Courier $gui_font_size";
unless ($font_override) {
	$top->optionAdd("*font", $gui_font);
	$top->optionAdd("*ROText.font", $text_font);
	$top->optionAdd("*Text.font", $text_font);
	$top->optionAdd("*HList.font", "{Verdana} ".($gui_font_size-1)) if $os eq 'win';
}
		
my $Balloon = $top->Balloon();
$Balloon->configure(-state => 'none') if $balloon_help == 0;
$Balloon->configure(-state => 'balloon') if $balloon_help == 1;

# Notebook reference hash
my %nb_page_ref = (
		1 => 'pd',
		2 => 'bis',
		3 => 'qpcr',
		4 => 'seq',
		);

# Load user-defined default variables
load_defaults();


#---------#
# Menubar #
#---------#

# It seems impossible to get the same menu behaviour across platforms with the 
# same code!  This is a compromise which seems to work on my systems.  Hopefully
# this will hold for other systems too ...

my $menu = $top->Menu(-bd => 1, -relief => $menu_relief, -type=>'menubar');

# setting $menu as the -menu for $top causes geometry problems under *nix ...
# packing the menubar as a widget causes alignment problems under Win32 ...
# (sometimes I think it was easier with menubuttons ....)
$menu->pack(-fill=>'x') if $os eq 'nix';
$top->configure(-menu => $menu) if $os eq 'win';

# recently used items
my $menu_mru = $menu->Menu(-title => "Recently opened", -menuitems => [
		[command => "- none available -", -state => "disabled" ]
		]);
recently_used_files();
 
my $menu_file = $menu->cascade(-label => "File", -menuitems => [
		[command => "New file", -command =>\&new_file, -accelerator=>'Ctrl-N' ],
		[command => "Open ...", -command =>\&pp_file_open, -accelerator=>'Ctrl-O' ],
		[cascade => "Open Previous", -menu=>$menu_mru ],
		"-",
		[command => "Save ...", -command =>\&pp_file_save, -accelerator=>'Ctrl-S' ],
		[command => "Save as ...", -command =>[sub {pp_file_save(1)}], ],
		"-",
		[command => "Retrieve gene from Ensembl", -command =>\&get_ensembl, -accelerator=>'Ctrl-E' ],
		"-",
		[command => "Exit", -command => \&end_prog, -accelerator=>'Ctrl-Q' ]
		]);		
		
my $menu_tools = $menu->cascade(-label => 'Tools', -menuitems => [
		[command => "Find primers for cloning", -command => \&get_primers_cloning, ],
		[command => "Add cloning sequences", -command => \&find_re_sites, ],
		"-",
		[command => "Generate report", -command => \&generate_report, -accelerator=>'Ctrl-R' ],
		"-",
		[command => "Save default values for this page", -command => \&save_defaults, ],
		[command => "Restore in-built default values", -command => \&restore_defaults, ],
		"-",
		[command => "Preferences", -command => \&prefs, -accelerator=>'Ctrl-P' ],
		]);
my $menu_help = $menu->cascade(-label => "Help", -menuitems => [
		['Checkbutton' => "Balloon help", -variable => \$balloon_help, -command => \&balloon_toggle],
		"-",
		[command => "Acknowledgements", -command => \&acknowledgements ],
		[command => "About ...", -command => \&info, -accelerator=>'Ctrl-Shift-A' ],
		]);


#--------------#
# Key bindings #
#--------------#

# If you're wondering about the anonymous sub for file opening - it's because
# Perl/Tk sends the calling widget as the first argument ...

$top->bind('<Control-n>' => \&new_file);
$top->bind('<Control-o>' => [sub {pp_file_open()}]);
$top->bind('<Control-s>' => [sub {pp_file_save()}]);
$top->bind('<Control-q>' => \&end_prog);
$top->bind('<Control-r>' => \&generate_report);
$top->bind('<Control-p>' => \&prefs);
$top->bind('<Control-A>' => \&info);
$top->bind('<Control-e>' => \&get_ensembl);


#---------#
# Toolbar #
#---------#

my $separator_menu = $top->Frame(-relief => "groove", -height=>2, -bd => 2);
$separator_menu->pack(-fill=> 'x') if $os eq 'win';

	my $toolbar = $top->Frame(-relief => $menu_relief, -bd => 1) ->pack(-fill=> 'x');
	my $tool_new = pack_button($toolbar, $top->Pixmap(-data => $icon_new), \&new_file)->pack(-side=>'left');
		$Balloon->attach($tool_new, -balloonposition => 'mouse', -balloonmsg => "New file");
	my $tool_open = pack_button($toolbar, $top->Pixmap(-data => $icon_open), \&pp_file_open)->pack(-side=>'left');
		$Balloon->attach($tool_open, -balloonposition => 'mouse', -balloonmsg => "Open PerlPrimer file");
	my $tool_save = pack_button($toolbar, $top->Pixmap(-data => $icon_save), \&pp_file_save)->pack(-side=>'left');
		$Balloon->attach($tool_save, -balloonposition => 'mouse', -balloonmsg => "Save file");
	my $tool_save_as = pack_button($toolbar, $top->Pixmap(-data => $icon_save_as), [sub {pp_file_save(1)}])->pack(-side=>'left');
		$Balloon->attach($tool_save_as, -balloonposition => 'mouse', -balloonmsg => "Save file as ...");

my $tool_sep1 = $toolbar->Label(-image => $top->Pixmap(-data => $icon_separator))->pack(-side=>'left', -padx=>2);

	my $tool_ensembl = pack_button($toolbar, $top->Pixmap(-data => $icon_ensembl), \&get_ensembl)->pack(-side=>'left'); 
		$Balloon->attach($tool_ensembl, -balloonposition => 'mouse', -balloonmsg => "Retrieve gene from Ensembl");

my $tool_sep2 = $toolbar->Label(-image => $top->Pixmap(-data => $icon_separator))->pack(-side=>'left', -padx=>2);

	my $tool_prefs = pack_button($toolbar, $top->Pixmap(-data => $icon_prefs), \&prefs)->pack(-side=>'left'); 
		$Balloon->attach($tool_prefs, -balloonposition => 'mouse', -balloonmsg => "Preferences");
	
	my $tool_report = pack_button($toolbar, $top->Pixmap(-data => $icon_report), \&generate_report)->pack(-side=>'left');
		$Balloon->attach($tool_report, -balloonposition => 'mouse', -balloonmsg => "Generate report");

my $separator_tools = $top->Frame(-relief => "groove", -height=>2, -bd => 2);
$separator_tools->pack(-fill=> 'x') if $os eq 'win';
	
#---------------#
# Draw Notebook	#
#---------------#

# notebook inactivebackground calculation
# (uses a 90% shade of the default widget colour - as per Perl/Tk *nix defaults.
# Win32 does not use this by default, but I feel it improves usability)
my $nb_colour = $top->cget(-bg);
my ($red, $green, $blue) = $top->rgb($nb_colour);
$red *= 0.9;
$green *= 0.9;
$blue *= 0.9;
$nb_colour = sprintf "%lx%lx%lx", $red, $green, $blue;

my $nb = $top->NoteBook(
		-relief => 'raised',
		-inactivebackground => "#$nb_colour",
		-bd => 1,
	)->pack(
		-expand=>1,
		-padx=>2,
		-pady=>4,
		-fill=>'both',
	);

# Notebook pages
my $page_primer_design = $nb->add("pd", -label=>"Standard PCR", -anchor=>"nw");
my $page_bisul_seq = $nb->add("bis", -label=>"Bisulphite PCR", -anchor=>"nw");
my $page_qpcr = $nb->add("qpcr", -label=>"Real-time PCR", -anchor=>'nw');
my $page_seq = $nb->add("seq", -label=>"Sequencing", -anchor=>"nw");
my $page_primer = $nb->add("primer", -label=>"Primers", -anchor=>"nw");

# Event binding to allow name display
$page_primer_design->bind('<Visibility>' => \&update_title);
$page_bisul_seq->bind('<Visibility>' => \&update_title);
$page_qpcr->bind('<Visibility>' => \&update_title);
$page_seq->bind('<Visibility>' => \&update_title);

# Have to use frames to get the nw anchor - grid plonks everything in the centre ...
my $page_primerf = $page_primer->Frame()->pack(-anchor=>'nw', -expand=>'1', -fill=>'both');
my $page_primerfb = $page_primer->Frame()->pack(-side=>'bottom', -anchor=>'sw');

my $page_primer_designf = $page_primer_design->Frame()->pack(-anchor=>'nw', -expand=>'1', -fill=>'both');
my $page_primer_designfb = $page_primer_design->Frame()->pack(-side=>'bottom', -anchor=>'sw');

my $page_sequencingf = $page_seq->Frame()->pack(-anchor=>'nw', -expand=>'1', -fill=>'both');
my $page_sequencingfb = $page_seq->Frame()->pack(-side=>'bottom', -anchor=>'sw');

my $page_bisul_seqf = $page_bisul_seq->Frame()->pack(-anchor=>'nw', -expand=>'1', -fill=>'both');
my $page_bisul_seqfb = $page_bisul_seq->Frame()->pack(-side=>'bottom', -anchor=>'sw');

my $page_qpcrf = $page_qpcr->Frame()->pack(-anchor=>'nw', -expand=>'1', -fill=>'both');
my $page_qpcrfb = $page_qpcr->Frame()->pack(-side=>'bottom', -anchor=>'sw');

#------------#
# Status bar #
#------------#

my $status_bar_f = $top->Frame()->pack(-side=>'bottom', -fill=>'x', -padx=>0,);
my $status_bar = $status_bar_f->ROText(-height=>1,
		-bg=>'#eeeeee',
		-relief=>'sunken',
		)->pack(-fill=>'x');

$status_bar->tagConfigure('red',
	-foreground => 'red');
$status_bar->tagConfigure('blue',
	-foreground => 'blue');
$status_bar->tagConfigure('grey',
	-foreground => 'grey20');

	
tie (*SBAR, 'Tk::Text', $status_bar);

# Packing system:
# This started off as an extremely simple system and has just kept on growing!
# It's now pretty robust and provides a good looking layout of widgets, all of
# which can be nicely arranged within labframes.  Another advantage is providing
# a consistent look-and-feel.

# NB: This has really become a complicated GUI - and the performance hit when drawing/
# resizing is quite noticable on old machines.  I'm not too sure of a way around this, 
# though, without sacrificing good looks or usability ...
# (Suggestions welcome!!)

# Format is:  
# (well, it *was* ... this has changed somewhat now - needs re-documenting!  Sometime ...)
#
# IMPORTANT UPDATE: reliance on grid() has now been greatly reduced (only used for labframe
# widgets) - each widget row is now separately defined by the nr() subroutine and widgets
# are packed from the left until nr() is called again.  Much simpler and easier!
# (NB - GUI code still has lots of relic options in the pack_gui() calls - needs to be
# tidied up a lot!)
#
# <label> <variable/subroutine> <widget hash reference> <frame reference>
# with optional <length of entry widget> <column>
#		<r = radiobutton> <option1> <option2> ... 
#               <t = text> <width> <height> <colspan> <wrap> <optional: weight> <optional: labframe> <optional: column> <optional: filebutton> <fb orientation>
#               <b = button> <state>
#		<c = checkbutton> <offvalue> <onvalue>
#		<h = hlist> <width> <height> <colspan> <columns> <optional: weight> <browse command sub> <optional: labframe>
#		<v = canvas> <colspan> <optional: in labframe?>
#		<l = labframe> <optional: column> <optional: rowspan>
#		<lt = label text> <string1> <$ (for textvariable marker)> <string2> ...

my $gui_count=0;
my $guilab_count=0;
my @row_counter;

#-------------#
# Primer page #
#-------------#

pack_gui('Forward primer', '', 'forward_l', \$page_primerf, 'l');
	nr(\$prf{forward_l});
		pack_gui('Sequence', $fprimer, 'fprimerseq', \$prf{forward_f1}, 40, 0 , 2);
	nr();
		pack_gui('Tm: ', $fprimer_tm, 'fprimertm', \$prl{forward_l}, 'ltg', 0, $fprimer_tm, '°C');
		pack_gui('Length: ', $fprimer_len, 'fprimerlen', \$prl{forward_l}, 'ltg', 1, $fprimer_len, 'bases');
		pack_gui('dS°: ', $fprimer_ds, 'fprimerds', \$prl{forward_l}, 'ltg', 0, $fprimer_ds, 'eu');
		pack_gui('dH°: ', $fprimer_dh, 'fprimerdh', \$prl{forward_l}, 'ltg', 1, $fprimer_dh, 'kcal/mol');
		pack_gui("dG°$pd_temperature: ", $fprimer_dg, 'fprimerds', \$prl{forward_l}, 'ltg', 0, $fprimer_dg, 'kcal/mol');

pack_gui('Reverse primer', '', 'reverse_l', \$page_primerf, 'l');
	nr(\$prf{reverse_l});
		pack_gui('Sequence', $rprimer, 'rprimerseq', \$prl{reverse_l}, 40, 0, 2);
	nr();
		pack_gui('Tm: ', $rprimer_tm, 'rprimertm', \$prl{reverse_l}, 'ltg', 0, $rprimer_tm, '°C');
		pack_gui('Length: ', $rprimer_len, 'rprimerlen', \$prl{reverse_l}, 'ltg', 1, $rprimer_len, 'bases');
		pack_gui('dS°: ', $rprimer_ds, 'rprimerds', \$prl{reverse_l}, 'ltg', 0, $rprimer_ds, 'eu');
		pack_gui('dH°: ', $rprimer_dh, 'rprimerdh', \$prl{reverse_l}, 'ltg', 1, $rprimer_dh, 'kcal/mol');
		pack_gui("dG°$pd_temperature: ", $rprimer_dg, 'fprimerds', \$prl{reverse_l}, 'ltg', 0, $rprimer_dg, 'kcal/mol');

pack_gui('Calculate Tm', \&get_tm, 'tmbutton', \$page_primerfb, 'b', 'active', 0);
pack_gui('BLAST primers', \&blast_primers, 'blastbutton', \$page_primerfb, 'b', '', 1);

pack_gui('Dimers', '', 'dim', \$page_primerf, 't', 60, 15, 1, 0, 2);
	$prf{dim}->bind('<Key-Right>' => \&jump_to_tm);
	$prf{dim}->bind('<Key-Left>' => \&jump_back);


$prf{dim}->configure(-fg=>'grey30');
$prf{dim}->tagConfigure('blue',
	-foreground => 'midnightblue');
$prf{dim}->tagConfigure('red',
	-foreground => 'red');
$prf{dim}->tagConfigure('black',
	-foreground => 'black');

# Layout tweak to make the grid column take up all possible space:
$page_primerf->gridColumnconfigure(0,-weight=>1);

tie (*DIMER, 'Tk::Text',$prf{'dim'});

#--------------------#
# Primer Design page #
#--------------------#

# $gui_count=0;
# $guilab_count=0;
# 
pack_gui('Primer Tm', '', 'primer_tm_l', \$page_primer_designf, 'l');
	nr(\$prl{primer_tm_l});
		pack_gui('', $min_tm_pr, 'mintm', \$prl{primer_tm_l}, 4, '', '', '');
		pack_gui('-', $max_tm_pr, 'maxtm', \$prl{primer_tm_l}, 4, 1, '', '°C');
		pack_gui(' Difference ', $max_diff_pr, 'maxdiff', \$prl{primer_tm_l}, 3, 2, '', '°C');
	
pack_gui('Primer Length', '', 'primer_len_l', \$page_primer_designf, 'l', 1);
	nr(\$prl{primer_len_l});
		pack_gui('', $pri_win_min_pr, 'minwin', \$prl{primer_len_l}, 4);
		pack_gui('-', $pri_win_max_pr, 'maxwin', \$prl{primer_len_l}, 4, 1, '', 'bases');
			
pack_gui('Amplified range', '', 'primer_range_l', \$page_primer_designf, 'l');
	nr(\$prl{primer_range_l}, 2);
		pack_gui('5\'', $max_range_5p_pr, 'maxrange5p', \$prf{primer_range_tf}, 5);
  		pack_gui('-', $min_range_pr, 'minrange', \$prf{primer_range_tf}, 5, 1);
  		pack_gui(' 3\'', $max_range_pr, 'maxrange', \$prf{primer_range_tf}, 5, 2);
		pack_gui('-', $max_range_3p_pr, 'maxrange3p', \$prf{primer_range_tf}, 5, 3);
	nr('');
		pack_gui('Amplicon size: ', '', 'minamp', \$prf{primer_ampsize_tf}, 'lt', 0, $min_ampsize_pr, '-', $max_ampsize_pr, 'bases');
	nr('', 0);
  		pack_gui('Set from ORF', \&reset_bounds, 'primer_getgene', '', 'b', 'normal', 0);
		pack_gui('-10', \&step_in, 'primer_stepin', '', 'b', 'normal', 1);
		pack_gui('+10', \&step_out, 'primer_stepout','' , 'b', 'normal', 2);
	
pack_gui('Options', '', 'primer_options_l', \$page_primer_designf, 'l', 1);
 	nr(\$prl{primer_options_l});
		pack_gui('Exclude %GC', $exclude_gc, 'exclude_gc', \$prl{primer_options_l}, 'c', 0, 1);
		pack_gui('GC clamp', $exclude_clamp, 'gc_clamp', \$prl{primer_options_l}, 'c', 0, 1, 1);
 	nr();
		pack_gui('Add 5\' F seq ', $primer_seq_5f, 'primer_seq_5f', \$prf{primer_options_seqadd}, 12, '','', ' ');
		pack_gui('Frame', $primer_seq_5f_frame, 'primer_seq_5f_frame', \$prf{primer_options_seqadd}, 2, 1,'', ' ');
		# not currently implemented and may never be ...
		# pack_gui('ATG', $primer_seq_5f_atg, 'primer_seq_5f_atg', '', 'c', 0, 1);
	nr();
		pack_gui('Add 5\' R seq ', $primer_seq_5r, 'primer_seq_5r', \$prf{primer_options_seqadd}, 12, '','', ' ');
		pack_gui('Frame', $primer_seq_5r_frame, 'primer_seq_5r_frame', \$prf{primer_options_seqadd}, 2, 1, '', ' ');
		# pack_gui('ATG', $primer_seq_5r_atg, 'primer_seq_5r_atg', '', 'c', 0, 1);

pack_gui('Sequence', '', 'seq', \$page_primer_designf, 't', 60, 6, 2, 0, 1, '', '', 1);

pack_gui('Results', '', 'res', \$page_primer_designf, 'h', 60, 10, 3, 10, 2, \&browse_primer, 1);
pack_gui('Canvas', '', 'primer_canvas', \$prl{res}, 'v', 4, 1);

pack_gui('Find primers', \&get_primers, 'primerbutton', \$page_primer_designfb, 'b', 'active', 0);
pack_gui('Find inwards', \&get_primers_auto_in, 'autobuttonin', \$page_primer_designfb, 'b', 'normal', 1);
pack_gui('Find outwards', \&get_primers_auto_out, 'autobuttonout', \$page_primer_designfb, 'b', 'normal', 2);
pack_gui('Cancel', \&cancel, 'primer_cancel', \$page_primer_designfb, 'b', 'normal', 3);
pack_gui('Copy selected', \&copy_selected_primers, 'primer_view', \$page_primer_designfb, 'b', 'normal', 4);

# We need to fiddle to get the two columns looking neat:
$page_primer_designf->gridColumnconfigure(0,-weight=>1);
$page_primer_designf->gridColumnconfigure(1,-weight=>1);


#-----------------#
# Sequencing page #
#-----------------#

pack_gui('Primer Tm', '', 'seq_tm_l', \$page_sequencingf, 'l');
	nr(\$prl{seq_tm_l});
		pack_gui('', $min_tm_seq, 'smintm', \$prl{seq_tm_l}, 4, '', '', '');
		pack_gui('-', $max_tm_seq, 'smaxtm', \$prl{seq_tm_l}, 4, 1, '', '°C');

pack_gui('Primer Length', '', 'seq_len_l', \$page_sequencingf, 'l', 1);
	nr(\$prl{seq_len_l});
		pack_gui('', $pri_win_min_seq, 'sminwin', \$prl{seq_len_l}, 4);
		pack_gui('-', $pri_win_max_seq, 'smaxwin', \$prl{seq_len_l}, 4, 1, '', 'bases');

pack_gui('Spacing / Coverage', '', 'seq_amp_l', \$page_sequencingf, 'l');
	nr(\$prl{seq_amp_l}, 2);
		pack_gui('Primers every ', $seq_spacing_min, 'sspacingmin', \$prf{primer_seq_amp_f1}, 5, '', '', '');
		pack_gui('-', $seq_spacing_max, 'sspacingmax', \$prf{primer_seq_amp_f1}, 5, 2, '', 'bases');
	nr('');
		pack_gui('Range:  5\' ', $min_range_seq, 'seqminrange', \$prf{primer_seq_amp_f2}, 5);
		pack_gui('-', $max_range_seq, 'seqmaxrange', \$prf{primer_seq_amp_f2}, 5, 2, '', ' 3\'');
	nr('', 0);
		pack_gui('Set range from ORF', \&reset_bounds, 'seq_getgene', '', 'b', 'normal', 0);
		
pack_gui('Options', '', 'seq_options_l', \$page_sequencingf, 'l', 1);  
	nr(\$prl{seq_options_l}, 0);
		pack_gui('Exclude %GC', $exclude_gc_seq, 'exclude_gc_seq', \$prf{primer_seq_options_f1}, 'c', 0, 1);
		pack_gui('GC clamp', $exclude_clamp_seq, 'gc_clamp_seq', \$prf{primer_seq_options_f1}, 'c', 0, 1, 1);
	nr('', 0);
		pack_gui('Exclude self-complimentarity >', $exclude_pd_seq, 'exclude_pd_seq', \$prf{primer_seq_options_f2}, 'c', 0, 1);
		pack_gui('-', $seq_pd_min, 'spdmin', \$prf{primer_seq_options_f2}, 3, 1, '', 'dG°37');

pack_gui('Sequence', '', 'seq_seq', \$page_sequencingf, 't', 60, 6, 2, 0, 1, '', '', 1);

pack_gui('Results', '', 'seq_res', \$page_sequencingf, 'h', 60, 10, 3, 5, 2, \&browse_primer, 1);
pack_gui('Canvas', '', 'seq_canvas', \$prl{seq_res}, 'v', 4, 1);

pack_gui('Find primers', \&get_seq_primers, 'seq_button', \$page_sequencingfb, 'b', 'active', 0);
pack_gui('Cancel', \&cancel, 'seq_cancel', \$page_sequencingfb, 'b', 'normal', 1);
pack_gui('Copy selected', \&copy_selected_primers, 'seq_view', \$page_sequencingfb, 'b', 'normal', 2);
	
$page_sequencingf->gridColumnconfigure(0,-weight=>1);
$page_sequencingf->gridColumnconfigure(1,-weight=>1);

	
#----------------------------#
# Bisulphite sequencing page #
#----------------------------#

# $gui_count=$guilab_count=0;
pack_gui('Primer Tm', '', 'bisul_tm_l', \$page_bisul_seqf, 'l');
	nr(\$prl{bisul_tm_l});
		pack_gui('', $min_tm_bs, 'mintm', \$prl{bisul_tm_l}, 4);
		pack_gui('-', $max_tm_bs, 'maxtm', \$prl{bisul_tm_l}, 4, 2, '', '°C');
		pack_gui(' Difference ', $max_diff_bs, 'maxdiff', \$prl{bisul_tm_l}, 3, 3, '', '°C');

pack_gui('Primer Length', '', 'bisul_len_l', \$page_bisul_seqf, 'l', 1);
	nr(\$prl{bisul_len_l});
		pack_gui('', $pri_win_min_bs, 'minwin', \$prl{bisul_len_l}, 4);
		pack_gui('-', $pri_win_max_bs, 'maxwin', \$prl{bisul_len_l}, 4, 2, '', 'bases');

pack_gui('Amplified range', '', 'bisul_range_l', \$page_bisul_seqf, 'l');
	nr(\$prl{bisul_range_l}, 2);
		pack_gui('5\'', $max_range_5p_bs, 'bisul_maxrange5p', \$prf{bisul_range_tf}, 5);
 		pack_gui('-', $min_range_bs, 'bisul_minrange', \$prf{bisul_range_tf}, 5, 1);
  		pack_gui(' 3\'', $max_range_bs, 'bisul_maxrange', \$prf{bisul_range_tf}, 5, 2);
		pack_gui('-', $max_range_3p_bs, 'bisul_maxrange3p', \$prf{bisul_range_tf}, 5, 3);
	nr('');
		pack_gui('Amplicon size: ', '', 'bisul_minamp', \$prf{bisul_ampsize_tf}, 'lt', 0, $min_ampsize_bs, '-', $max_ampsize_bs, 'bases');
	nr('', 0);
  		pack_gui('Set from CpG island', \&reset_bounds, 'bisul_getcpg', '', 'b', 'normal');
		pack_gui('-10', \&step_in, 'primer_stepin', '', 'b', 'normal');
		pack_gui('+10', \&step_out, 'primer_stepout', '', 'b', 'normal');

pack_gui('Options', '', 'bisul_options_l', \$page_bisul_seqf, 'l', 1);
	nr(\$prl{bisul_options_l}, 0);
		pack_gui('Exclude CpG', $exclude_cpg, 'bisul_exclude_cpg', \$prl{bisul_options_l}, 'c', 0, 1);
		pack_gui('Require 3\' C', $exclude_3c, 'bisul_exclude_cs', \$prl{bisul_options_l}, 'c', 0, 1);
	nr('', 0);
		pack_gui('Repeats / runs post conversion', $pre_bs, 'bisul_pre_bs', \$prl{bisul_options_l}, 'c', 0, 1, 2);
	nr();	
		pack_gui('Min C%', $bisul_min_c, 'bisul_min_c', \$prl{bisul_options_l}, 4, 2);

pack_gui('Sequence', '', 'bisul_seq', \$page_bisul_seqf, 't', 60, 6, 2, 0, 1, '', '', 1);

pack_gui('Results', '', 'bisul_res', \$page_bisul_seqf, 'h', 60, 10, 3, 10, 2, \&browse_bisulphite, 1);
pack_gui('Canvas', '', 'bisul_canvas', \$prl{bisul_res}, 'v', 4, 1);

pack_gui('Find primers', \&get_bisulphite, 'bisul_button', \$page_bisul_seqfb, 'b', 'active', 0);
pack_gui('Find inwards', \&get_primers_auto_in, 'bisul_autobuttonin', \$page_bisul_seqfb, 'b', 'normal', 1);
pack_gui('Find outwards', \&get_primers_auto_out, 'bisul_autobuttonout', \$page_bisul_seqfb, 'b', 'normal', 2);
pack_gui('Cancel', \&cancel, 'bisul_cancel', \$page_bisul_seqfb, 'b', 'normal', 3);
pack_gui('Copy selected', \&copy_selected_primers, 'bisul_view', \$page_bisul_seqfb, 'b', 'normal', 4);

# More column tidying ...
$page_bisul_seqf->gridColumnconfigure(0,-weight=>1);
$page_bisul_seqf->gridColumnconfigure(1,-weight=>1);

#-----------#
# QPCR page #
#-----------#

# $gui_count=0;
# $guilab_count=0;
pack_gui('Primer Tm', '', 'qprimer_tm_l', \$page_qpcrf, 'l');
	nr(\$prl{qprimer_tm_l});
		pack_gui('', $min_tm_q, 'qmintm', \$prl{qprimer_tm_l}, 4);
		pack_gui('-', $max_tm_q, 'qmaxtm', \$prl{qprimer_tm_l}, 4, 2, '', '°C');
		pack_gui(' Difference ', $max_diff_q, 'qmaxdiff', \$prl{qprimer_tm_l}, 3, 3, '', '°C');

pack_gui('Primer Length', '', 'qprimer_len_l', \$page_qpcrf, 'l', 1);
	nr(\$prl{qprimer_len_l});
		pack_gui('', $pri_win_min_q, 'qminwin', \$prl{qprimer_len_l}, 4);
		pack_gui('-', $pri_win_max_q, 'qmaxwin', \$prl{qprimer_len_l}, 4, 2, '', 'bases');

pack_gui('Amplicon size', '', 'qprimer_amp_l', \$page_qpcrf, 'l');
	nr(\$prl{qprimer_amp_l});
		pack_gui('', $min_ampsize_q, 'qminamp', \$prl{qprimer_amp_l}, 5);
		pack_gui('-', $max_ampsize_q, 'qmaxamp', \$prl{qprimer_amp_l}, 5, 2, '', 'bases');
	
pack_gui('Options', '', 'qprimer_options_l', \$page_qpcrf, 'l', 1);
	nr(\$prl{qprimer_options_l}, 0);
		pack_gui('Exclude %GC', $exclude_gc, 'qexclude_gc', \$prf{qprimer_options_tf}, 'c', 0, 1);
		pack_gui('GC clamp', $exclude_clamp, 'qgc_clamp', \$prf{qprimer_options_tf}, 'c', 0, 1, 2);
	nr();
		pack_gui('Minimum Intron/Exon boundary: ', $exclude_ie, 'qexclude_ie', \$prl{qprimer_options_l}, 3, '', '', 'bases');

pack_gui('Genomic sequence', '', 'qdna_seq', \$page_qpcrf, 't', 25, 5, 1, 0, 1, '', '', 1, 1);
pack_gui('mRNA sequence', '', 'qmrna_seq', \$page_qpcrf, 't', 25, 5, 1, 0, 1, '', 1, 1, 1);

pack_gui('Results', '', 'qpcr_res', \$page_qpcrf, 'h', 60, 10, 3, 10, 2, \&browse_primer, 1);
pack_gui('Canvas', '', 'qprimer_canvas', \$prl{qpcr_res}, 'v', 4, 1);

pack_gui('Find primers', \&get_qprimers, 'qprimerbutton', \$page_qpcrfb, 'b', 'active', 0);
pack_gui('Cancel', \&cancel, 'qprimer_cancel', \$page_qpcrfb, 'b', 'normal', 1);
pack_gui('Copy selected', \&copy_selected_primers, 'qprimer_view', \$page_qpcrfb, 'b', 'normal', 2);
pack_gui('Spidey output', \&view_spidey_out, 'qprimer_spidey', \$page_qpcrfb, 'b', 'normal', 3);

# We need to fiddle to get the two columns looking neat:
$page_qpcrf->gridColumnconfigure(0,-weight=>1);
$page_qpcrf->gridColumnconfigure(1,-weight=>1);

#-------------#
# Window icon #
#-------------#

my $pixmap = $top->Pixmap(-data => $perlprimer_icon);

# set icon for main window ...
$top->Icon(-image => $pixmap);


#-------------#
# Popup Menus #
#-------------#

# A big problem with perl/Tk:
# HList does not have any way to bind to column headers or even
# to identify in which column a click was - this is the only way to
# provide sorting ...  It's very clumsy!

# (update: after using this for a while, I actually think it works rather well
# for this system - since not all the Hlist column headers are generally visible
# at any one time, it's nice to have this system in place)

my $sort_reverse=0;
my $sort_prev = 10;

$popup_sort = $top->Menu(-menuitems => [
		[cascade => 'Sort primers by ...', -menuitems => [
			# [command => 'Sort primers:', -state=>'disabled'],
			# '-',
			[command => 'Forward position', -command => [\&sort_primers, 1]],
			[command => '  ... length', -command => [\&sort_primers, 2]],
			[command => '  ... Tm', -command => [\&sort_primers, 3]],
			'-',
			[command => 'Reverse position', -command => [\&sort_primers, 8]],
			[command => '  ... length', -command => [\&sort_primers, 6]],
			[command => '  ... Tm', -command => [\&sort_primers, 7]],
			'-',
			[command => 'Amplicon size', -command => [\&sort_primers, 9]],
			[command => 'Dimer dG', -command => [\&sort_primers, 10]],
			'-',
			[Checkbutton => 'Reversed', -variable => \$sort_reverse, -command =>[\&sort_primers, "r"]],
		]],
		'-',
		[command => 'Select all', -command => [\&select_all_primers]],
		[command => 'Copy', -command => [\&copy_selected_primers]],
		'-',
		[command => 'Set range from selected', -command => [\&primer_take_range]],
		]);

$popup_sort_seq = $top->Menu(-menuitems => [
		[Button => 'Select all', -command => [\&select_all_primers]],
		'-',
		[Button => 'Copy', -command => [\&copy_selected_primers]],
		]);

# Text menu
my $popup_text = $top->Menu(-menuitems => [
		[Button => 'Open sequence ...', -command => [\&open_seq, \$text_widget_ref]],
		[Button => 'Save sequence ...', -command => [\&save_seq, \$text_widget_ref]],
		'-',
		[Button => 'Select All', -command => [\&select_all_text, \$text_widget_ref]],
		[Button => 'Clear', -command => [\&clear_text, \$text_widget_ref]],
		[Button => 'Lower case', -command => [\&lc_text, \$text_widget_ref]],
		'-',
		[Button => 'Reverse complement', -command => [\&reverse_complement, \$text_widget_ref]],
		]);

# Show the main window
$top->deiconify();

# prevent the user from resizing the window any smaller (which makes things hidden/messy)
my $min_width = $top->geometry;
my ($min_x, $min_y) = ($min_width =~ /(\d.*)x(\d.*?)\+/);
$top->minsize($min_x, $min_y);

# Open file if specified on the command line
if ($commandline) {
	$top->update;
	pp_file_open($commandline);
}

# Override HList bindings for left/right and other commands...
my $class = "Tk::HList";
$top->bind($class,'<Right>', \&jump_to_tm);
$top->bind($class,'<Left>', \&jump_back);
$top->bind($class,'<2>', \&copy_selected_primers);
$top->bind($class,'<Control-a>', \&select_all_primers);
$top->bind($class,'<Control-c>', \&copy_selected_primers);
$top->bind($class,'<3>', \&menu_popup);

MainLoop();

# Program end
exit_program();

		#------------------------------#
		#  			       #
		#  Subroutines start here ...  #
		#  			       #
		#------------------------------#


#------------------#
# Post-GUI routine #
#------------------#

sub exit_program {
	# Save updated prefs, (including Balloon help)
	my $file_data = "";
	foreach my $i (keys %pref_variables) {
		my $pointer = $pref_variables{$i};
		$file_data .= "$i = $$pointer\n";
	}
	
	foreach my $i (keys %pref_arrays) {
		my $pointer = $pref_arrays{$i};
		$file_data .= "$i = [".join(",",@$pointer)."]\n";
	}
	
	open (PREFS, ">$pref_file") || die "Error: could not open $pref_file for writing: $!\n";
	print PREFS $file_data;
	close (PREFS);
	exit 0;
}


#-------------------------#
# GUI building subroutine #
#-------------------------#

# The messy guts that allow all the other GUI parts of the program to be so neat
# and tidy.  It's nothing fancy - in fact, it's downright inelegant - but it gets
# the job done, keeps things consistent, and writing new gui code (eg, the
# preferences dialogue) is a matter of minutes to do ...

sub nr {
	# New row signal
	my ($reference, $pady) = @_;
	$reference ||= $old_reference;
	
	$pady = $frame_pady unless defined($pady);
	
	push @row_counter, $$reference->Frame(
			)->pack(-side=>'top', -anchor=>'w', -pady=>$pady, -padx=>0);
	$old_reference = $reference;
	return;
}

sub pack_gui {		
	my ($label, $reference, $widget_name, $widget_reference) = @_; 
	my ($col, $colspan, $len, $wrap, $scrollbars, $widget_ref, $pady, $pack);
	
	# Check for signals ...
	if ($_[4] eq "t") {
		
		#------------#
		# Text entry #
		#------------#
		
		if ($_[8] == 1) {
			$wrap = "none";
			$scrollbars = "osoe";
		} elsif ($_[8] == 2) {
			$wrap = "none";
			$scrollbars = "s";
		} elsif ($_[8] == 3) {
			$wrap = "word";
			$scrollbars = "oe";
		} else {
			$wrap = "char";
			$scrollbars = "oe";
		}
		
		$col = ($_[11] ? $_[11] : 0);
		if ($_[11]) {
			$gui_count--
		}

		$colspan = ($_[7] ? $_[7] : 1);
		
		my $side =  'left'; 
		my $bside = 'top';
		my $text_ref;
		
		if ($_[9]) {
			$prl{"$widget_name"} = $$widget_reference->LabFrame(
					-label=>"$label",
					-labelside=>'acrosstop',
				)->grid(-column=>$col, -row=>$gui_count, -columnspan=>$colspan, -sticky=>'nsew', -padx=>4);
			$text_ref = \$prl{"$widget_name"};
		} else {
			$text_ref = $widget_reference;
		}
		
		my $widget_type = ($_[12]) ? 'Text' : 'ROText';
		
		$prf{"$widget_name"}=$$text_ref->Scrolled($widget_type,
				-scrollbars => $scrollbars,
				-relief => 'groove',
				-bg => '#eeeeee',
				-wrap => $wrap,
				-width => $_[5],
				-height => $_[6],
			)->pack(-expand=>1, -fill=>'both', -side=>$side, -padx=>2, -pady=>2);
		# button flag?
		if ($_[12]) {
			#open/save/clear button stuff ...
			$prb{"$widget_name open"}=pack_button($prl{$widget_name}, $top->Pixmap(-data => $icon_dna_open), [\&open_seq, \$prf{$widget_name}])
					->pack(-side=>$bside, -anchor=>'w', -fill=>'x');
			$Balloon->attach($prb{"$widget_name open"}, -balloonposition => 'mouse', -balloonmsg => "Open DNA sequence");
			
			$prb{"$widget_name save"}=pack_button($prl{"$widget_name"}, $top->Pixmap(-data => $icon_dna_save), [\&save_seq, \$prf{"$widget_name"}])
					->pack(-side=>$bside, -anchor=>'w', -fill=>'x');
			$Balloon->attach($prb{"$widget_name save"}, -balloonposition => 'mouse', -balloonmsg => "Save DNA sequence");
			
			$prb{"$widget_name clear"}=pack_button($prl{"$widget_name"}, $top->Pixmap(-data => $icon_clear), [\&clear_text, \$prf{"$widget_name"}])
					->pack(-side=>$bside, -anchor=>'w', -fill=>'x');
			$Balloon->attach($prb{"$widget_name clear"}, -balloonposition => 'mouse', -balloonmsg => "Clear sequence");
			
			my $text_widget = \$prf{"$widget_name"};
			$$text_widget->menu(undef);
			$$text_widget->bind("<3>", [\&menu_text, \$prf{"$widget_name"}, $widget_name]);
		}
		
		if ($_[9]) {
			$$widget_reference->gridRowconfigure($gui_count, -weight=>$_[9]);
		}
		$widget_ref=\$prf{"$widget_name"};
		bind_mousewheel($$widget_ref);
	} elsif ($_[4] eq "h") {
		
		#-------------#
		# hlist entry #
		#-------------#
		
		my $sub=$_[10];
		
		if ($_[8]) {
			$wrap = "none";
			$scrollbars = "osoe";
		} else {
			$wrap = "char";
			$scrollbars = "oe";
		}

		$colspan = ($_[7] ? $_[7] : 1);
		
		$prl{$widget_name} = $$widget_reference->LabFrame(
				-label=>"$label",
				-labelside=>'acrosstop',
			)->grid(-column=>0, -row=>$gui_count, -columnspan=>$colspan, -sticky=>'nsew', -padx=>4);
		
		$prf{$widget_name}=$prl{$widget_name}->Scrolled('HList',
				-scrollbars => $scrollbars,
				-relief => 'groove',
				-bg => '#eeeeee',
				-width => $_[5],
				-height => $_[6],
				# -font => $hlist_font,
				-header => 1,
				-selectmode => 'extended',
				-columns => $_[8],
				-browsecmd => \&$sub,
				-command => \&hlist_command,
			)->pack(-expand=>1, -fill=>'both', -padx=>2, -pady=>2);
		
		if ($_[9]) {
			$$widget_reference->gridRowconfigure($gui_count, -weight=>$_[9]);
		}
		
		# various primer bindings ...
		# (assume the only use of Hlist is for primer display)
		my $hlist = \$prf{"$widget_name"};
		bind_mousewheel($$hlist);
		header_create($hlist, $widget_name);	
		$widget_ref=\$prf{$widget_name};		
	} elsif ($_[4] eq "b") {
		
		#--------#
		# Button #
		#--------#
		
		$prb{"$widget_name"}=($widget_reference ? $$widget_reference : $row_counter[-1])->Button(
				-text=>"$label",
				# -font => $font,
				-padx=>5,
				-pady=>$button_pady,
				-command=>\$_[1],
				-state=>($_[9] ? $_[9] : 'normal'),
				-default=>($_[5] ? $_[5] : 'normal'),
			)->pack(-side=>'left', -anchor=>'w', -pady=>$button_pack_pady, -padx=>$button_pack_padx);
		$widget_ref=\$prb{"$widget_name"};
	} elsif ($_[4] eq "v") {
		
		#--------#
		# Canvas #
		#--------#
		
		$prb{"$widget_name info"}=pack_button($$widget_reference, $top->Pixmap(-data => $icon_info), \&canvas_info)
				->pack(-side=>'right', -expand=>0, -anchor=>'e', -fill=>'none');
			$Balloon->attach($prb{"$widget_name info"}, -balloonposition => 'mouse', -balloonmsg => "Graphical display help");
		$prb{"$widget_name magnify"}=pack_button($$widget_reference, $top->Pixmap(-data => $icon_magnify), [\&dna_magnify, 0])
				->pack(-side=>'right', -expand=>0, -anchor=>'e', -fill=>'none');
			$Balloon->attach($prb{"$widget_name magnify"}, -balloonposition => 'mouse', -balloonmsg => "Magnified view");
		
		$prf{"$widget_name"}=$$widget_reference->Canvas(
				-height=>$dna_canvas_height,
			)->pack(-expand=>1, -fill=>'x', -pady=>2, -padx=>6, -side=>'left');
		
		# Bindings for DNA canvas
		my $canvas = \$prf{"$widget_name"};
		$$canvas->CanvasBind('<B1-Motion>' => sub {items_drag($Tk::event->x, $canvas)});
		$$canvas->CanvasBind('<B2-Motion>' => sub {amplicon_drag($Tk::event->x, $canvas)});
		$$canvas->CanvasBind('<Control-B1-Motion>' => sub {amplicon_drag($Tk::event->x, $canvas)});
		$$canvas->CanvasBind('<3>' => sub {dna_magnify($Tk::event->x)});
		$$canvas->CanvasBind('<Control-3>', \&draw_dna);

		$widget_ref=\$prf{"$widget_name"};
	} elsif ($_[4] eq "c") {
		
		#--------------#
		# Check button #
		#--------------#

		$prc{"$widget_name"}=$row_counter[-1]->Checkbutton(
				-offvalue=>"$_[5]",
				-onvalue=>"$_[6]",
				# -font => $font,
				-text=>$label,
				-pady=>$check_pady,
				-variable=>\$_[1]
			)->pack(-side=>'left', -anchor=>'w', -pady=>0, -ipady=>0);
		$widget_ref=\$prc{"$widget_name"};
	} elsif ($_[4] eq "l") {
		
		#----------#
		# LabFrame #
		#----------#
		
		$col = ($_[5] ? $_[5] : 0);
		if ($_[5]) {
			$guilab_count--
		}
		
		my $rowspan = ($_[6] ? $_[6] : 1);

		$prf{"$widget_name"}=$$widget_reference->LabFrame(
					-label=>"$label",
					-labelside=>'acrosstop',
				)->grid(-column=>$col, -row=>$guilab_count, -rowspan=>$rowspan, -sticky=>'nswe', -padx=>4);
		$prl{"$widget_name"}=$prf{"$widget_name"}->Frame(
				)->pack(-anchor=>'nw', -expand=>'1', -fill=>'none');
		$guilab_count++;
		$widget_ref=\$prf{"$widget_name"};
	} elsif ($_[4] eq "lt") {
		
		#------------#
		# label text #
		#------------#
		
		$len = @_ -1;

		$prl{$widget_name}=$row_counter[-1]->Label(
				-text=>$label,
				-font => $font,
				-justify=>'left'
			)->pack(-side=>'left', -padx=>0, -ipadx=>0);
		
		for my $i (6 .. $len) {
			# Very messy, but it works ...
			if (defined($_[$i]) && ($_[$i] =~ /^[^\d]/)) {
				$prl{"$widget_name$i"}=$row_counter[-1]->Label(
						-text=>$_[$i],
						# -font => $font,
						-justify=>'left',
					)->pack(-side=>'left', -padx=>0, -ipadx=>0);
			} else {
				$prl{"$widget_name$i"}=$row_counter[-1]->Label(
						# -font => $font,
						-textvariable=>\$_[$i],
						-justify=>'right',
					)->pack(-side=>'left', -padx=>0, -ipadx=>0);
			}
		}
		$widget_ref=\$prl{"$widget_name"};
	} elsif ($_[4] eq "ltg") {
		
		#-----------------#
		# label text grid #
		#-----------------#
		
		$col = ($_[5] ? $_[5] : 0);
		$gui_count-- unless $col == 0;
		$len = @_ -1;
		
		$prf{$widget_name}=$row_counter[-1]->Frame(
				)->grid(-column=>$col, -row=>$gui_count, -sticky=>'w', -pady=>0, -padx=>4);

		$prl{$widget_name}=$prf{$widget_name}->Label(
				-text=>$label,
				# -font => $font,
				-justify=>'left'
			)->pack(-side=>'left', -padx=>0, -ipadx=>0);
		
		for my $i (6 .. $len) {
			# Very messy, but it works ...
			if (defined($_[$i]) && ($_[$i] =~ /^[^\d]/)) {
				$prl{"$widget_name$i"}=$prf{$widget_name}->Label(
						-text=>$_[$i],
						# -font => $font,
						-justify=>'left',
					)->pack(-side=>'left', -padx=>0, -ipadx=>0);
			} else {
				$prl{"$widget_name$i"}=$prf{$widget_name}->Label(
						# -font => $font,
						-textvariable=>\$_[$i],
						-justify=>'right',
					)->pack(-side=>'left', -padx=>0, -ipadx=>0);
			}
		}
		$widget_ref=\$prf{"$widget_name"};

	} elsif ($_[4] eq "f") {
		
		#-------#
		# Frame #
		#-------#
		
		$pady = (defined($_[5]) ? $_[5] : 2);
		$prf{$widget_name}=$$widget_reference->Frame(
				)->pack(-side=>'top', -anchor=>'w', -pady=>$pady, -padx=>0);
		$widget_ref=\$prf{$widget_name};
	} elsif ($_[4] eq "be") {
		
		#--------------#
		# Browse Entry #
		#--------------#
		
		my $array_ref = $_[6];
		my $width = ($_[7] ? $_[7] : 20);
		my $listwidth = ($_[8] ? $_[8] : undef);
		$prl{$widget_name}=$row_counter[-1]->Label(
				-text=>$label,
				-font => $font,
				-justify=>'left'
			)->pack(-side=>'left', -padx=>0, -ipadx=>0);
		$pre{$widget_name}=$row_counter[-1]->BrowseEntry(
				-variable=>$reference,
				-choices=> [@$array_ref],
				-listwidth => $listwidth,
				-width => $width,
				# -bg => '#eeeeee',
			)->pack(-side=>'left', -padx=>0, -ipadx=>0);
		$widget_ref=\$pre{$widget_name};
	} else {
		
		#---------------------#
		# label:entry(:label) #
		#---------------------#
		
		$prl{"$widget_name"}=$row_counter[-1]->Label(
				-text=>"$label",
				# -font => $font,
				-justify=>'left'
			)->pack(-side=>'left');
		
		$pre{"$widget_name"}=$row_counter[-1]->Entry(
				-relief => 'groove',
				-bg => '#eeeeee',
				-width => ($_[4] ? $_[4] : ($_[1] ? length($_[1])+1 : 4)),
				-textvariable=>\$_[1]
			)->pack(-side=>'left');
		if ($_[7]) {
			$prl{"$widget_name-r"}=$row_counter[-1]->Label(
					-text=>"$_[7]",
					# -font => $font,
					-justify=>'left'
				)->pack(-side=>'left');
		}
		$widget_ref=\$pre{"$widget_name"};
	}
	$gui_count++;
	
	#--------------#
	# Balloon help #
	#--------------#
		
	if ($balloonmsg{$widget_name}) {
		$Balloon->attach($$widget_ref, -balloonposition => 'mouse', -balloonmsg => $balloonmsg{$widget_name});
	}
}

sub pack_button {
	my ($widget, $image, $sub_ref) = @_;
	return $widget->Button(
			-relief => 'flat',
			-image => $image,
			-command=> $sub_ref,
			-activebackground=>$activebackground_color,
		) if $os eq 'win';
	return $widget->Button(
			-relief => 'flat',
			-image => $image,
			-command=> $sub_ref,
		)->pack(-side=>'left');
}

#---------------------#
# MouseWheel bindings #
#---------------------#

sub bind_mousewheel {
	# Thanks to Slaven Rezic, Steve Liddie and "Mastering Perl/Tk" for this routine (slightly modified) ...
    	my ($w) = @_;
		
    	if ($^O eq 'MSWin32') {
		# Windows bindings
        	$w->bind('<MouseWheel>' => [ sub { 
				$_[0]->yview('scroll', -($_[1] / 120) * $scroll_factor, 'units') 
			}, Ev('D') ]
        	);
    	} else {
		# *nix bindings
        	$w->bind('<4>' => sub {
            			$_[0]->yview('scroll', -$scroll_factor, 'units') unless $Tk::strictMotif;
        		});
        	$w->bind('<5>' => sub {
                  		$_[0]->yview('scroll', $scroll_factor, 'units') unless $Tk::strictMotif;
        		});
    	}
}

sub update_title {
	my $nb_page = which_nb_page();
	if ($open_file{$nb_page}) {
		if ($open_file{$nb_page} eq 'File not saved') {
			$top->configure(-title=>"PerlPrimer v$version");
			return;
		}
		my $file = $open_file{$nb_page};
		$top->configure(-title=>"PerlPrimer v$version - $file");
	}
	return;
}

#-------------------#
# Read DNA sequence #
#-------------------#

sub read_windows {
	my ($sub_ref, $dnaseq_f, $dnaseq_f_len, $dnaseq_r) = @_;
	
	$reverse=0;
	sbarprint("\nWindowing forward seq ...");
	@PF = &$sub_ref($dnaseq_f, $dnaseq_f_len);
	$pfkeys=@PF;
	
	if ($dnaseq_r) {
		$reverse=1;
		sbarprint("\nWindowing reverse seq ...");
		@PR = &$sub_ref($dnaseq_r, $dnaseq_f_len);
	}
}


#-------------------------------------------------#
# basic complementation and bisulphite conversion #
#-------------------------------------------------#

sub complement {
	$_ = shift;
	tr/AGCTagct/TCGAtcga/;
	return $_;
}

sub bisulphite {
	$_ = shift;
	tr/C/T/;
	return $_;
}


#----------------------------#
# primer window calculations #
#----------------------------#

sub primer_window {
	my ($primer_seq, $dnaseq_len) = @_;
	my $key;
	my @primer_list=();
	my $repeat_real = $repeat-1;
	
	my ($min_tm, $max_tm, $max_diff, $pri_win_min, $pri_win_max, $min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p) 
		= get_variables(qw(min_tm max_tm max_diff pri_win_win pri_win_max min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p)); 

	# old method - using amplicon range only ...
	my ($seqbound5,$seqbound3)=($$min_range,$$max_range);

	$seqbound5 ||= 0;	
	$seqbound5 = $$max_range-$$max_ampsize if ($$max_range) && ($$max_ampsize) && ($reverse==0) && ($$max_range-$$max_ampsize>0);
	$seqbound5 = $dnaseq_len-$$min_range-$$max_ampsize if ($$min_range) && ($reverse==1) && ($dnaseq_len-$$min_range-$$max_ampsize > 0);
	
	if (($$max_ampsize) && ($$min_ampsize)) {	
		if ($reverse==0) {
			unless ($$min_range) {
				$seqbound3 = $dnaseq_len-$$min_ampsize;
			} else {
				$seqbound3 = ($dnaseq_len-$$min_ampsize>$$min_range ? $$min_range : $dnaseq_len-$$min_ampsize);
			}
		} else {
			$seqbound3 = ($$max_range ? $dnaseq_len-$$max_range+$$pri_win_max : $dnaseq_len)
		}
	}
	
	# the above covers the unlikely possibility that $max_ampsize is set
	# and $max_range_5p and $max_range_3p are not.  It might be desirable at some
	# point, and it was the original way of doing things.  It doesn't hurt!
	# ... and is now used for qpcr!
	
	# To use without specific range boundaries but only with a minimum and 
	# maximum amplicon distance, keep $max_range_5p and $min_range equal,
	# and likewise at the 3' end

	if (defined($$max_range_5p)&&defined($$max_range_3p)) {
		unless (($$max_range_5p==$$min_range)&&($$max_range_3p==$$max_range)) {
			if ($reverse==0) {
				$seqbound5 = $$max_range_5p;
				$seqbound3 = $$min_range + $$pri_win_max;
			}
			
			if ($reverse==1) {
				$seqbound5 = $dnaseq_len-$$max_range_3p-1;
				$seqbound3 = $dnaseq_len-$$max_range + $$pri_win_max;
			}
		}
	}
		
	for my $i ($seqbound5 .. $seqbound3) {
		for my $primer_win ($$pri_win_min .. $$pri_win_max) {
			next if $i+$primer_win>$dnaseq_len;
			
			my $currwindow=substr($primer_seq, $i, $primer_win);
		
			if ($exclude_gc) {
				# exclude based on %GC content
				my $gc=gc($currwindow);
				next if $gc < $min_gc or $gc > $max_gc;
			}
			
			if ($exclude_clamp) {
				# calculate GC clamp at 5' end: two out of last three residues required...
				my $gc_clamp=gc(substr($currwindow, -3, 3));
				next unless $gc_clamp > 50;
			}
			
			if ($exclude_rr) {
				# simple run exclusion:
				# ($run is number of consecutive bases to exclude)
				next if ($currwindow =~ /(C{$run,}|A{$run,}|G{$run,}|T{$run,})/);
	
				# simple repeat exclusion:
				# ($repeat is number of consecutive repeats of two or more bases)
				next if ($currwindow =~ /(.{2,})\1{$repeat_real,}/);
			}
			
			# All possible primers have now been excluded
			# Get the tm and check this against the accepted range:
			my ($tm)=tm($currwindow);
			if (($$min_tm < $tm) && ($tm < $$max_tm)) {
				# Looks good, so let's add the primer to the array
				push @primer_list, [$i, $primer_win, $currwindow, $tm];
			}
		}		 
	}
	return @primer_list;
}


#---------------------------------#
#  bisulphite window calculations #
#---------------------------------#

sub bisul_window {
	# See the primer_window subroutine above for general comments
	my ($primer_seq, $dnaseq_len)=@_;
	my $key;
	my $repeat_bs_real = $repeat_bs-1;

	my ($seqbound5,$seqbound3) = 0;
	my @primer_list=();
	
	my ($min_tm, $max_tm, $max_diff, $pri_win_min, $pri_win_max, $min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p) 
		= get_variables(qw(min_tm max_tm max_diff pri_win_win pri_win_max min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p)); 


	$seqbound5 = 0;
	$seqbound5 = $$max_ampsize-$$max_range if ($$max_range) && ($reverse==0) && ($$max_ampsize-$$max_range>0);
	$seqbound5 = $dnaseq_len-$$min_range-$$max_ampsize if ($$min_range) && ($reverse==1) && ($dnaseq_len-$$min_range-$$max_ampsize > 0);
	
	if (defined($$max_range_5p)&&defined($$max_range_3p)) {
		unless (($$max_range_5p==$$min_range)&&($$max_range_3p==$$max_range)) {
			if ($reverse==0) {
				$seqbound5 = $$max_range_5p;
				$seqbound3 = $$min_range + $$pri_win_max;
			}
			
			if ($reverse==1) {
				$seqbound5 = $dnaseq_len-$$max_range_3p-1;
				$seqbound3 = $dnaseq_len-$$max_range + $$pri_win_max;
			}
		}
	} else {
		if ($reverse==0) {
			unless ($$min_range) {
				$seqbound3 = $dnaseq_len-$$min_ampsize;
			} else {
				$seqbound3 = ($dnaseq_len-$$min_ampsize>$$min_range ? $$min_range : $dnaseq_len-$$min_ampsize);
			}
		} else {
			$seqbound3 = ($$max_range ? $dnaseq_len-$$max_range+$$pri_win_max : $dnaseq_len)
		}
	}
		
	return unless defined($seqbound5) && defined($seqbound3);
	
	my $count;
	for my $i ($seqbound5 .. $seqbound3) {
		for my $primer_win ($$pri_win_min .. $$pri_win_max) {
			$count++;
			next if $i+$primer_win>$dnaseq_len;

			my $currwindow=substr($primer_seq, $i, $primer_win);
			
			# Check for %C content
			next if cc($currwindow) < $bisul_min_c;

			#next if ga($currwindow) < 50;
			
			# Check for exclude_cpg
			# If exclude_cpg is set, no primers with CpG residues are include
			# Otherwise, primers are allowed to have up to a maximum of 2 CpG residues
			if ($exclude_cpg) {
				next if $currwindow =~ /CG/
			} else {
				my $count_cpg = 0;
				while ($currwindow =~ /CG/g) {
					$count_cpg++
				}
				next if $count_cpg >2;
			}
			
			# Check for 3' C content
			if ($exclude_3c) {
				# It is important that primer specificity to the converted template is given
				# by having C residues at the 3' end
				
				# The following gives the most flexibility: a primer is accepted if
				# it either ends in a C or has two of the last three bases as C's
				next unless ((substr($currwindow, -1, 1) eq "C") or (cc(substr($currwindow, -3, 3)) > 50));
			}
			
			# now convert the primer ...
			my $currwindowbs = bisulphite($currwindow);
			
			# do we exclude repeats and runs before or after bisulphite conversion?
			my $checkref = ($pre_bs ? \$currwindowbs : \$currwindow);
			
			if ($exclude_rr_bs) {
				# Run exclusion:
				next if ($$checkref =~ /(C{$run_bs,}|A{$run_bs,}|G{$run_bs,}|T{$run_bs,})/);
	
				# Repeat exclusion:
				next if ($$checkref =~ /(.{2,})\1{$repeat_bs_real,}/);
			}
			
			# Need to complement if we're creating reverse primers,
			# as the two converted strands of DNA are non-complementary
			if ($reverse==1) {
				$currwindowbs=complement($currwindowbs);
			}

			my ($tm)=tm($currwindowbs);
			
			if (($$min_tm < $tm) && ($tm < $$max_tm)) {
				# Add the original, unconverted sequence to the array as well:
				# Note - this will not be complemented if reversed seq
				push @primer_list, [$i, $primer_win, $currwindowbs, $tm, $currwindow];
			}
		}		 
	}
	return @primer_list;
}


#-----------------------------------------#
# Calculation of amplicon lengths and     #
# compatible temperature matching primers #
#-----------------------------------------#

sub calc_amplicon {
	my $dnaseq_r_len = shift;
	my $count_amplicon=0;
	@primer_pairs=();
	my $count=0;
	
	my ($min_tm, $max_tm, $max_diff, $pri_win_min, $pri_win_max, $min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p) 
		= get_variables(qw(min_tm max_tm max_diff pri_win_win pri_win_max min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p)); 
		
	# (the reverse sequence was reversed to make it 5'->3')
	# So I either have to recalculate the revese numbers or
	# not reverse the sequence and count backwards!
	#
	# Solution:
	#
	# realpos_start = total_length - calc_pos
	# and closest_real_base = realpos_start - window
	#               thus    = total_length - calc_pos - window
	
	foreach my $key_f (@PF) {
		return if ($cancel==1);
		
		my ($pos_f, $len_f, $seq_f, $tm_f, $unconverted_f ) = @$key_f;
				
		foreach my $key_r (@PR) {
			return if ($cancel==1);
			my ($pos_r, $len_r, $seq_r, $tm_r, $unconverted_r ) = @$key_r;
			
			my $realpos_r = $dnaseq_r_len - $pos_r - 1;
			my $amp_size=$realpos_r-$pos_f;
			
			# Skip all pairs where forward > reverse ...			
			next unless ($realpos_r-$len_r>$pos_f+$len_f);
			
			# Range check (messy):
			# Now should be covered in windowing! but need for qpcr
			if (($$min_range) && ($$max_range)) {
				next unless ($pos_f < $$min_range && $$max_range < $realpos_r);
			}
						
			# Temperature check:
			### could be more efficient
			next unless ($tm_r > ($tm_f-$$max_diff) && $tm_f > ($tm_r-$$max_diff));
			
			# QPCR specific:
			if ($qpcr_flag == 1) {
				# Amplicon size check
				next if ($amp_size < $$min_ampsize || $amp_size > $$max_ampsize);

				#check amp range spans i/e boundary
				my $qpcr_check=0;
				foreach my $i (@intron_exon_bounds) {
					$qpcr_check=1 if $pos_f<$i && $i<$realpos_r
				}
				next unless $qpcr_check==1;
				
				#check at least one primer falls across i/e boundary
				$qpcr_check=0;
				foreach my $i (@intron_exon_bounds) {
					$qpcr_check=1 if $pos_f<($i-$exclude_ie) && ($i+$exclude_ie)<($pos_f+$len_f);
					$qpcr_check=1 if ($realpos_r-$len_r)<($i-$exclude_ie) && ($i+$exclude_ie)<$realpos_r;
				}
				next unless $qpcr_check==1;
			}
                        
			# all OK up to here, so let's check primer-dimers:
			# NB: this is a big speed hit ...
			my (@pd_score, @pd_sorted)=();

             		push(@pd_score,primer_dimer($seq_f,$seq_f));
             		push(@pd_score,primer_dimer($seq_f,$seq_r));
             		push(@pd_score,primer_dimer($seq_r,$seq_r));
             		@pd_sorted=sort{$a <=> $b} @pd_score;

			$tm_f = sprintf("%.2f", $tm_f);
			$tm_r = sprintf("%.2f", $tm_r);
			
			push @primer_pairs, [ $seq_f, $pos_f, $len_f, $tm_f,
					$seq_r, $pos_r, $len_r, $tm_r, $realpos_r, $amp_size, $pd_sorted[0], ( $unconverted_f ?  $unconverted_f : 0), ( $unconverted_r ?  $unconverted_r : 0) ];
			
		}
		$count++;
		my $percent = sprintf("%.f", ($count/$pfkeys)*100);
		sbarprint("\nCalculating amplicons ... $percent% completed ...");
	}
	sbarprint("\nFinished ... found ".($#primer_pairs+1)." primer pairs");
}


#-----------------------------------------#
# Calculation of amplicon lengths and     #
# compatible temperature matching primers #
#-----------------------------------------#

sub calc_seq_primers {
	# scaled down calc_amplicon routine
	my $count_amplicon=0;
	@primer_pairs=();
	my $count=0;

	my $pdcheck=1;
	
	my ($min_range, $max_range) 
		= get_variables(qw(min_range max_range)); 
	
	my $last_seq;
	
	$pd_full = 1;
	# $pd_extensible = 0;
		
	foreach (@PF) {
		return if ($cancel==1);
		
		$count++;
		my $percent = sprintf("%.f", ($count/$pfkeys)*100);
		sbarprint("\nFinding primers ... $percent% completed ...");
		
		my ($pos, $len, $seq, $tm_f ) = @$_;
		
		# take the first one that fits if last_seq hasn't been set ...
		unless (defined($last_seq)) {
			next if $pos < $$min_range;
			if ($pos > $seq_spacing_max+$$min_range) {
				dialogue("Contiguous primers could not be found within the set parameters");
				last;
			}
		} else {
			next if $pos > $$max_range;
			next if $pos < $last_seq + $seq_spacing_min;
			if ($pos > $last_seq + $seq_spacing_max) {
				dialogue("Contiguous primers could not be found within the set parameters");
				last;
			}
		}		
		
		# exclude all primer-dimers if asked ...
		my $pd;
		if ($exclude_pd_seq == 1) {
             		next if ($pd = primer_dimer($seq,$seq)) < -$seq_pd_min;
		}
					
		$tm_f = sprintf("%.2f", $tm_f);
		
		push @primer_pairs, [ $seq, $pos, $len, $tm_f, $pd ];
		$last_seq = $pos;		
	}
	
	# $pd_extensible = 1;
	sbarprint("\nFinished ... found ".($#primer_pairs+1)." primer pairs");
}

# sub calc_seq_primers {
	# # scaled down calc_amplicon routine
	# my $count_amplicon=0;
	# @primer_pairs=();
	# my $count=0;
# 
	# my $pdcheck=1;
	# 
	# my ($min_range, $max_range) 
		# = get_variables(qw(min_range max_range)); 
	# 
	# my $last_seq;
	# my $seq_spacing_mid = int(($seq_spacing_max+$seq_spacing_min)/2);
	# 
	# $pd_full = 1;
	# # $pd_extensible = 0;
	# 
	# # build hash table
	# my $i;
	# my $flag=0;
	# my $skip;
	# my $interval=1;
	# my %seq_primers;
	# foreach (@PF) {
		# $i++;
		# my ($pos) = @$_;
		# $seq_primers{$pos} = $i;
	# }
	# 
	# my $j;	
	# for ($j=$$min_range; $j < $$max_range;) {
		# my ($pos, $len, $seq, $tm_f);
		# return if ($cancel==1);
		# 
		# if ($skip) {
			# $flag = 1-$flag;
			# if ($flag) {
				# $j+=$interval;
			# } else {
				# $j-=$interval;
				# $interval++;
			# }
			# $skip=0;
		# }
		# 
		# $count++;
		# # my $percent = sprintf("%.f", ($count/$pfkeys)*100);
		# sbarprint("\nFinding primers ... $count primers checked ...");
		# 
		# 
		# # take the first one that fits if last_seq hasn't been set ...
		# unless (defined($last_seq)) {
			# if ($seq_primers{$j}) {
				# ($pos, $len, $seq, $tm_f ) = @{ $PF[$seq_primers{$j}] };
			# } else {
				# # we want the first primer to be as near to the start as possible
				# print "$j\n";
				# $j++;
				# next;
			# }
			# print "$j: OK! pos was $pos\n";
			# 
			# next if $pos < $$min_range;
			# if ($pos > $seq_spacing_max+$$min_range) {
				# dialogue("Contiguous primers could not be found within the set parameters");
				# last;
			# }
		# } else {
			# if ($seq_primers{$j}) {
				# ($pos, $len, $seq, $tm_f ) = @{ $PF[$seq_primers{$j}] };
			# } else {
				# $skip=1;
				# next;
			# }
			# if ($pos > $$max_range) {
				# $skip=1;
				# next;
			# }
			# if ($pos < $last_seq + $seq_spacing_min) {
				# $skip=1;
				# next;
			# }
			# if ($pos > $last_seq + $seq_spacing_max) {
				# dialogue("Contiguous primers could not be found within the set parameters");
				# last;
			# }
		# }		
		# 
		# # exclude all primer-dimers if asked ...
		# my $pd;
		# if ($exclude_pd_seq == 1) {
             		# if ($pd = primer_dimer($seq,$seq) < -$seq_pd_min) {
				# $skip=1;
				# next;
			# }
		# }
					# 
		# $tm_f = sprintf("%.2f", $tm_f);
		# 
		# push @primer_pairs, [ $seq, $pos, $len, $tm_f, $pd ];
		# $interval = 1;
		# $last_seq = $pos;
		# $j = $pos + $seq_spacing_mid;		
	# }
	# 
	# # $pd_extensible = 1;
	# sbarprint("\nFinished ... found ".($#primer_pairs+1)." primer pairs");
# }

#-------------------------#
# calculate %base content #
#-------------------------#

# %GC	
sub gc {
	$_ = $_[0];
	my($gc,$countgc,$counttotal);
	$gc=0;
	$countgc=0;
		
	$countgc = tr/GCgc/GCgc/;
	$counttotal = length();
	
	$gc = $countgc/$counttotal*100;
	return $gc;
}

# %C
sub cc {
	$_ = $_[0];
	my($cc,$countc,$counttotal);
	$cc=0;
	$countc=0;

	$countc = tr/C/C/;
	$counttotal = length();

	$cc = $countc/$counttotal*100;
	return $cc;
}
		

#-------------------------------#
# forward and reverse sequences #
#-------------------------------#

sub getseq {
	my ($min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p) 
		= get_variables(qw(min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p)); 

	if (($$min_range) && ($$max_range)) {
		$$min_ampsize=$$max_range - $$min_range;
		$$max_ampsize=$$max_range - $$min_range if ($$max_range - $$min_range)>$$max_ampsize;
	}

	# forward				
	$_ = uc($_[0]);
	s/[\s|\n]//g;
	s/[^ACGT]//g;
	my $dnaseq_f = $_;
	
	#reverse						
	my $dnaseq_r_top = reverse($dnaseq_f);
	my $dnaseq_r = complement($dnaseq_r_top);
				
	my $dnaseq_f_len = length($dnaseq_f);
	my $dnaseq_r_len = length($dnaseq_r);
	
	return ($dnaseq_f, $dnaseq_f_len, $dnaseq_r_top, $dnaseq_r, $dnaseq_r_len);
}


#-----
# subroutine for calculating Tm as per SantaLucia(1998), cited above, with Mg++
# (and K+, Tris++) concentration calculated via equations 7 and 8 from
# von Ahsen, et al, 2001, Clinical Chemistry 47(11):1956-1961
#-----

sub tm {
	my ($primer) = @_;
	$primer = uc($primer); # if user enters primer directly as lower-case
	my ($i, $nn, $initterm, $endterm);
	my $primer_len = length($primer);
	my ($deltaH, $deltaS);
		
	#-----------------------------#
	# calculate deltaH and deltaS #
	#-----------------------------#

	for ($i=0; $i<$primer_len-1; $i++) {
		$nn = substr($primer, $i, 2);
		
		$deltaH+= $oligo_dH{$nn};
		$deltaS+= $oligo_dS{$nn};
	}
		
	#-------------------------#
	# initial term correction #
	#-------------------------#

	$initterm="init" . substr($primer, 0, 1);
	$deltaH+= $oligo_dH{$initterm};
	$deltaS+= $oligo_dS{$initterm};
	
	$endterm="init" . substr($primer, -1, 1);
	$deltaH+= $oligo_dH{$endterm};
	$deltaS+= $oligo_dS{$endterm};
				
	# Tm at 1M NaCl
	# $tm= ($deltaH * 1000) / ($deltaS + (1.987 * log($oligo_conc / 4))) - 273.15;
	
	#------------------------------------------#
	# correct for salt concentration on deltaS #
	#------------------------------------------#
	
	my $na_eq=$monovalent_cation_conc + 120 * sqrt($mg_conc - $dntp_conc);
	$na_eq /= 1000;
	
	# deltaS correction:
	$deltaS += (0.368 * ($primer_len - 1) * log($na_eq));
	# print "na_eq was $na_eq\nentropy correct was ".(0.368 * ($primer_len - 1) * log($na_eq))."\n";
	
	my $oligo_conc_mols = $oligo_conc / 1000000000;

	# Salt corrected Tm
	# NB - for PCR I'm assuming for the moment that the [strand target] << [oligo]
	# and that therefore the C(t) correction term approx equals [oligo]
	my $corrected_tm=(($deltaH * 1000) / ($deltaS + (1.987 * log($oligo_conc_mols)))) - 273.15;
	return ($corrected_tm, $deltaH, $deltaS);
}


#----------------------------------------------------------#
# Recalculate the oligo_dG hash on current salt conditions #
#----------------------------------------------------------#

sub recalculate_dG {
	# because dG = dH - TdS, and dS is dependent on the salt concentration ...
	my $temperature = shift || $pd_temperature;
	my $na_eq=$monovalent_cation_conc + 120 * sqrt($mg_conc - $dntp_conc);
	$na_eq /= 1000;
	
	# the length of each NN dimer is 2, therefore the modifier is 1
	my $entropy_adjust = (0.368 * log($na_eq));
		
	foreach my $key (keys(%oligo_dH_full)) {
		next if $key =~ /init/; # the length of each monomer is 1, thus the modifier of dS is 0 and the values are precalulated
		
		my $dS = $oligo_dS_full{$key} + $entropy_adjust;
		my $dG = $oligo_dH_full{$key}-((273.15+$temperature)*($dS/1000));
		$oligo_dG{$key} = $dG;
	}
}

#--------------------------#
# Primer-dimer calculation #
#--------------------------#

sub primer_dimer {
	# This is my second attempt at implementing a primer-dimer system:
	# The first attempt aligned each combination together explicitly; this attempt
	# creates a binding matrix and then reads each primer combination from the
	# matrix.  It's not significantly faster, but it does have the advantage of
	# extending this algorithm to allow for unequal loops (although this has not
	# been implemented as yet) and providing the possiblity of reading all
	# primer-dimers (not just 3') by adjusting a single variable (which is used when
	# displaying primer information.
	
	### TODO:
	# It is easily possible to calculate *all* primers when searching for primer
	# pairs ... I really should add a preferences option to do this ...
	
	my ($primer_f, $primer_r) = @_;
	return unless ($primer_f) && ($primer_r);
			
	my ($k, $l);
	@score=();
	%primer_hash=();
	@score_sort=();
	@bind_string=();
	%rating_hash=();
		
	# $pl = greatest length
	$pfl=length($primer_f);
	$prl=length($primer_r);
	$pl = ($pfl>$prl ? $pfl : $prl);
	
	my $rcompr = reverse(complement($primer_r));
	my $rcomprlc = lc($rcompr);
	my $fprimer_r=lc(reverse($primer_f));
	$rprimer_r=reverse($primer_r);
	
	# Scan the primers against each other:
	# The default is to only consider 3' primer-dimers, for speed concerns (5'
	# pd's will reduce the primer population, but won't cause extendible dimers) -
	# however, setting $pd_full to 1 will calculate all primer-dimers.  This is
	# currently used only when viewing individual primers, for both speed concerns
	# and because it's 3' primer-dimers that are the real problem in PCR.
	
	# create a binding array for each of the four bases
	for $l (0 .. $pfl-1) {
		my $mbase = substr($fprimer_r, $l, 1);
		$primer_hash{$mbase}[$l]=1;
		for $k qw(a g c t) {
			$primer_hash{$k}[$l] ||=0;
		}
	}
		
	# create the primer matrix
	my @primer_comp;
	for $k (0 .. $prl-1) {
		$primer_comp[$k]=$primer_hash{substr($rcomprlc, $k, 1)};
	}
	
	# # print the matrix - for debugging
	# for $k (0 .. $prl-1) {
		# print "$k:\t@{$primer_comp[$k]}\n";
	# }
	
	# read each combination from the matrix, calculate dG for each dimer
	my $pd_len = ($pd_full ? $pfl+$prl-1 : $pl-2);
	for $k (0 .. $pd_len) {
		$score[$k]=0;
		my $bind;
		my $score_p=0;
		
		# extensible primer short-circuit - ignore all primers that will
		# not create extensible (i.e. amplifiable) dimers
		my $start = $k>$pfl-1 ? $pfl-1 : $k;
		my $end = $k>$prl-1 ? $prl-1 : $k;
		if ($pd_extensible && !$pd_full) {
			next unless $primer_comp[0][$start] == 1;
			next unless $primer_comp[$end][$start-$k] == 1;
		}
		
		# elsif ($pd_extensible) {
			# # no point reconsidering them!
			# next if $primer_comp[0][$start] == 1 && $primer_comp[$end][$start-$k] == 1;
		# }
		
		# read the binding data
		for $l (0 .. $prl-1) {
			if (($k-$l)<$pfl) {
				$bind .= $primer_comp[$l][$k-$l] if ($k-$l)>=0;
			} else {
				# spacer
				$bind .= "2";
			}
		}
		
		# Single matched bases surrounded by mismatches are unstable,
		# so we remove them with the regexp (look ahead is needed otherwise
		# strings of consecutive match/mismatches are not caught)
		$bind =~ s/01(?=[^1])/00/gx;
		
		# Short circuit if there's nothing to bind
		next unless $bind =~ /[1]/;
		
		# Find start and end of similarity
		my ($pb_init,$pb_end);
		for $l (0 .. length($bind)-1) {
			# at first I tried finding the initiating terminal bases with
			# regexps, but that was much slower ...
			if (substr($bind, $l, 1) eq "1") {
				defined($pb_init) || ($pb_init = $l);
				$pb_end=$l;
			}
		}
				
		if (defined($pb_init)) {
			# deltaG calculation
			for $l ($pb_init .. $pb_end-1) {
				next if substr($bind, $l, 2) eq "00";
				next if substr($bind, $l, 1) eq "2";
				$score_p+=$oligo_dG{substr($primer_f, $pfl-$k+$l-1, 2).substr($rprimer_r, $l, 2)};
			}
			
			# init term corrections
			my $initterm="init" . substr($rprimer_r, $pb_init, 1);
			$score_p+= $oligo_dG{$initterm};
			
			my $endterm="init" . substr($rprimer_r, $pb_end, 1);
			$score_p+= $oligo_dG{$endterm};
			
			# add to the hash ...
			$score[$k]=sprintf("%.2f",$score_p);
			$bind_string[$k]=$bind;
			$rating_hash{$score[$k]}=$k;
		}
	}
	
	# sort the dimers to give the most stable:	
	@score_sort = sort { $a <=> $b } @score;
		
	# Returns the most stable dimer
	return $score_sort[0];
}


#--------------------------------#
# Rountine to draw primer-dimers #
#--------------------------------#

sub draw_dimer {
	# This all seems a bit cumbersome!!
	my ($primer_f, $primer_r, $pos, $FH) = @_;
	
	my $rprimer_r=reverse($primer_r);
	my $dimer_binding="";
	my $pr_space="";
	my $fspace="";
	my $rspace="";
			
	my $fspace_def = $pl-$pfl>0 ? $pl-$pfl : 0;
	$fspace=" "x($fspace_def+($pos>$pl-1?$pos-$pl+1:0));
	
	if ($pos+1>=$pfl) {
		$rspace=" "x($pl-$pos-1);
	} else {
		$rspace=$fspace;
	}
	
	$pr_space=" "x($pfl-$pos-1);
	
	for my $j (0 .. $pos) {
		next unless $j < $prl;
		if (substr($bind_string[$pos],$j,1)==1) {
			$dimer_binding=$dimer_binding."|"
		} elsif (substr($bind_string[$pos],$j,1)==0) {
			$dimer_binding=$dimer_binding."."
		} else {
			$dimer_binding=$dimer_binding." "
		}
	}
				
	print $FH "$fspace"."5' "."$primer_f"." 3'\n".
		"$rspace"."   "."$pr_space"."$dimer_binding\n".
		"$rspace"."$pr_space"."3' "."$rprimer_r"." 5'\n\n";
}	


sub end_prog {
	# This sends the program to the pref writing routine after Mainloop() ...
	$top->destroy;
}

		
#---------------------#
# GUI button routines #
#---------------------#

sub get_tm {
	my ($report) = @_;
	# Lower case bug fix:
	$fprimer = uc($fprimer);
	$rprimer = uc($rprimer);
	
	my ($deltaG, $deltaH, $deltaS);
	if ($fprimer) {
		($fprimer_tm, $deltaH, $deltaS) = tm($fprimer);
		$fprimer_tm = sprintf("%.2f", $fprimer_tm);
		
		# since dG = dH - TdS, we don't need to calculate based on NN's ...
		$deltaG = $deltaH-((273.15+$pd_temperature)*($deltaS/1000));
		
		$fprimer_ds = sprintf("%.2f", $deltaS);
		$fprimer_dh = sprintf("%.2f", $deltaH);
		$fprimer_dg = sprintf("%.2f", $deltaG);
		$fprimer_len = length($fprimer);
	}
	
	if ($rprimer) {
		($rprimer_tm, $deltaH, $deltaS) = tm($rprimer);
		$rprimer_tm = sprintf("%.2f", $rprimer_tm);
		
		$deltaG = $deltaH-((273.15+$pd_temperature)*($deltaS/1000));
		
		$rprimer_ds = sprintf("%.2f", $deltaS);
		$rprimer_dh = sprintf("%.2f", $deltaH);
		$rprimer_dg = sprintf("%.2f", $deltaG);
		$rprimer_len = length($rprimer);
	}
	
	my $repeat_real = $repeat-1;
	$prf{dim}->delete(0.1,"end");
	$prf{dim}->insert('end', "Warning: forward primer run found\n", 'red') if ($fprimer =~ /(C{$run,}|A{$run,}|G{$run,}|T{$run,})/);
	$prf{dim}->insert('end', "Warning: forward primer repeat found : $1\n", 'red') if ($fprimer =~ /(.{2,})\1{$repeat_real,}/);

	$prf{dim}->insert('end', "Warning: reverse primer run found\n", 'red') if ($rprimer =~ /(C{$run,}|A{$run,}|G{$run,}|T{$run,})/);
	$prf{dim}->insert('end', "Warning: reverse primer repeat found\n", 'red') if ($rprimer =~ /(.{2,})\1{$repeat_real,}/);

	$prf{dim}->insert('end', "Most stable 3' extensible primer-dimers (at $pd_temperature°C)\n\n", 'blue');
	
	my ($pd1, $pd2, $pd3);
	if ($fprimer) {
		$pd1 = primer_dimer($fprimer,$fprimer);
		$pos=$rating_hash{$score_sort[0]};
		$prf{dim}->insert('end', "Forward vs. Forward: $score_sort[0] kcal/mol\n\n", 'black');
		draw_dimer($fprimer, $fprimer, $pos, \ *DIMER) unless ($score_sort[0]==0);
		
		if ($rprimer) {
			$pd2 = primer_dimer($fprimer,$rprimer);
			$pos=$rating_hash{$score_sort[0]};
			$prf{dim}->insert('end', "Forward vs. Reverse: $score_sort[0] kcal/mol\n\n", 'black');
			draw_dimer($fprimer, $rprimer, $pos, \ *DIMER) unless ($score_sort[0]==0);
		}
	}
	
	if ($rprimer) {		
		$pd3 = primer_dimer($rprimer,$rprimer,1);
		$pos=$rating_hash{$score_sort[0]};
		$prf{dim}->insert('end', "Reverse vs. Reverse: $score_sort[0] kcal/mol\n\n", 'black');
		draw_dimer($rprimer, $rprimer, $pos, \ *DIMER) unless ($score_sort[0]==0);
	}
	
	$pd_full = 1;
	$prf{dim}->insert('end', "\nMore stable non-extensible primer-dimers (at $pd_temperature°C), if any\n\n", 'blue');
	
	if ($fprimer) {
		primer_dimer($fprimer,$fprimer);
		$pos=$rating_hash{$score_sort[0]};
		if ($score_sort[0]<$pd1) {
			$prf{dim}->insert('end', "Forward vs. Forward: $score_sort[0] kcal/mol\n\n", 'black');
			draw_dimer($fprimer, $fprimer, $pos, \ *DIMER);
		}
		
		if ($rprimer) {
			primer_dimer($fprimer,$rprimer);
			$pos=$rating_hash{$score_sort[0]};
			if ($score_sort[0]<$pd2) {
				$prf{dim}->insert('end', "Forward vs. Reverse: $score_sort[0] kcal/mol\n\n", 'black');
				draw_dimer($fprimer, $rprimer, $pos, \ *DIMER);
			}
		}
	}
	if ($rprimer) {
		primer_dimer($rprimer,$rprimer);
		$pos=$rating_hash{$score_sort[0]};
		if ($score_sort[0]<$pd3) {
			$prf{dim}->insert('end', "Reverse vs. Reverse: $score_sort[0] kcal/mol\n\n", 'black');
			draw_dimer($rprimer, $rprimer, $pos, \ *DIMER);
		}
	}
	
	$pd_full = 0;
	
	# again ... messy tie in for report generating routine ...
	# Yes, I know this is a really, really, ugly way to do things!!
	if ($report) {
		my $dimer_text = $prf{dim}->get('0.1','end');
		my $primer_text = <<EOT;
Forward primer: $fprimer

Tm: $fprimer_tm°C\t\tLength: $fprimer_len bases
dS°: $fprimer_ds eu\t\tdH°: $fprimer_dh kcal/mol
dG°$pd_temperature: $fprimer_dg kcal/mol

		
Reverse primer: $rprimer

Tm: $rprimer_tm°C\t\tLength: $rprimer_len bases
dS°: $rprimer_ds eu\t\tdH°: $rprimer_dh kcal/mol
dG°$pd_temperature: $rprimer_dg kcal/mol

		
$dimer_text\n
EOT
		return $primer_text;
	}
	
	$prf{dim}->see(0.1);
}

sub blast_primers {
	return if check_packages("HTTP::Request", "LWP::UserAgent");
	# I've tried to make this the way blast *should* be, with colour output
	# and the ability to limit the results ... 
	if (($fprimer eq "") && ($rprimer eq "")) {
		dialogue("Please enter primers to blast");
		return;
	}
	my ($blast_out_f,$blast_out_r)=("Waiting on Blast Server for forward primer results ...","Waiting on Blast Server for reverse primer results ...");
	my ($blast_header_f, @blast_summary_f, @blast_results_f1, @blast_results_f2, @blast_results_f3);
	my ($blast_header_r, @blast_summary_r, @blast_results_r1, @blast_results_r2, @blast_results_r3);
	my $display='f';
	my $blast_search_string;
	$blast_status = "";
	
	my $blast_display = sub {
		# I don't know whether this is messy or beautiful ...
		# But I really love the fact that I can do it!!
		my ($blast_header_ref, $blast_summary_ref, $blast_results_1_ref, $blast_results_2_ref, $blast_results_3_ref, $search_string);
		if ($display eq 'f') {
			$blast_header_ref = \$blast_header_f;
			$blast_summary_ref = \@blast_summary_f;
			$blast_results_1_ref = \@blast_results_f1;
			$blast_results_2_ref = \@blast_results_f2;
			$blast_results_3_ref = \@blast_results_f3;
		} else {
			$blast_header_ref = \$blast_header_r;
			$blast_summary_ref = \@blast_summary_r;
			$blast_results_1_ref = \@blast_results_r1;
			$blast_results_2_ref = \@blast_results_r2;
			$blast_results_3_ref = \@blast_results_r3;
		}
		
		$prf{'blast_text'}->delete('0.1', 'end');
		unless ($$blast_header_ref) {
			return
		}
		$prf{'blast_text'}->insert('0.1', $$blast_header_ref, 'black');
		
		### TODO: hyperlinks??
		
		# Users might get confused about "."'s in the string box ...
		($search_string = $blast_search_string) ||=".";
		for my $i (0 .. $#$blast_summary_ref) {
			$_ = $$blast_summary_ref[$i];
			next unless /$search_string/;
			/^(\w*\|[\w\|\.\d]+)(.*?)(\d+\s+[\d\.\w\-]+\s*$)/;
			$prf{'blast_text'}->insert('end', $1);
			$prf{'blast_text'}->insert('end', $2, 'blue');
			$prf{'blast_text'}->insert('end', $3."\n");
			
		}
		for my $i (0 .. $#$blast_results_1_ref) {
			next unless $$blast_results_1_ref[$i] =~ /$search_string/;
			$prf{'blast_text'}->insert('end', "\n\n".$$blast_results_1_ref[$i]."\n\n ", 'blue');
			$prf{'blast_text'}->insert('end', $$blast_results_2_ref[$i]."\n\n", 'black');
			$prf{'blast_text'}->insert('end', $$blast_results_3_ref[$i]."\n");
		}
	};

	# GUI code
	# destroy old window if one exists (and cancel last blast search ...)
	if (Exists($blast_d)) {
		$top->afterCancel($rptid);
		$blast_d->destroy;
	}	
		
	$blast_d = $top->Toplevel(-title=>'BLAST Search');
	$blast_d->withdraw();
	my $blast_d_f = $blast_d->Frame()->pack(-expand=>1, -side=>'top', -fill=>'both');
	my $blast_d_fb = $blast_d->Frame()->pack(-side=>'left', -anchor=>'sw');
	my $blast_d_fs = $blast_d->Frame()->pack(-side=>'right', -anchor=>'se');
	$blast_d->Icon(-image => $pixmap);

	nr(\$blast_d_f);
	pack_gui('', '', 'blast_text', \$blast_d_f, 't', 95, 30, 1, 1);
	
	nr(\$blast_d_fb);

	pack_gui("OK", sub{
			$top->afterCancel($rptid);
			$blast_d->destroy;
		}, "blast_full", '', "b", "active", 0);
	pack_gui("Forward", sub {
			$display='f';
			&$blast_display;
			
		}, "blast_f", '', "b", '', 1, '', '', 'disabled');
		
	pack_gui("Reverse", sub {
			$display='r';
			&$blast_display;
		}, "blast_r", '', "b", '', 2, '', '', 'disabled');

	
	pack_gui('Status: ', $blast_status, 'blast_status', '', 'lt', 3, $blast_status);

	# text search button
	nr(\$blast_d_fs);
	pack_gui('String:', $blast_search_string, 'blast_search_string', \$blast_d_fs, 20);
	pack_gui('Search', sub {
			&$blast_display;
		}, "blast_search", '', "b", '', 1, '', '', 'disabled');
	# $top->update;
	# $blast_d->deiconify();
	$prf{blast_text}->configure(-fg=>'grey30');
	$prf{blast_text}->tagConfigure('blue',
		-foreground => 'midnightblue');
	$prf{blast_text}->tagConfigure('black',
		-foreground => 'black');
	
	$blast_status = "Sending BLAST request ...";
	$blast_d->update;
	$blast_d->deiconify();
	$top->update;
	
	my $blast_summary;
	
	if ($fprimer) {
		# Blast forward primer
		$blast_status="Blasting forward primer ...";
		
		($blast_header_f, $blast_summary) = blast($fprimer);
		if ($blast_summary) {
			@blast_summary_f = split("\n",$blast_summary);
			@blast_results_f1 = @blast_results_1;
			@blast_results_f2 = @blast_results_2;
			@blast_results_f3 = @blast_results_3;
			
			$display='f';
			&$blast_display;
			
			$prb{blast_f}->configure(-state=>'normal');
			$prb{blast_search}->configure(-state=>'normal');
		} else {
			$blast_status="Unable to connect to server";
		}
	}
	
	if ($rprimer) {
		#Blast reverse primer
		$blast_status="Blasting reverse primer ...";
		
		($blast_header_r, $blast_summary) = blast($rprimer);
		if ($blast_summary) {
			@blast_summary_r = split("\n",$blast_summary);
			@blast_results_r1 = @blast_results_1;
			@blast_results_r2 = @blast_results_2;
			@blast_results_r3 = @blast_results_3;
			
			$prb{blast_r}->configure(-state=>'normal');
			$prb{blast_search}->configure(-state=>'normal');
			
			$blast_status="Blast search complete";
		} else {
			$blast_status="Unable to connect to server";
		}
	}
}

sub get_primers {
	my ($max_ampsize, $max_range_5p, $max_range_3p) 
		= get_variables(qw(max_ampsize max_range_5p max_range_3p)); 

	$cancel=0;
	$bs=0;
	my $seq = get_seq();
	check_range();
	
	# set the max_ampsize in case the user has modified
	# either max_range_5p or max_range_3p by hand ...
	if (($$max_range_5p)&&($$max_range_3p)) {
		$$max_ampsize = $$max_range_3p-$$max_range_5p
	}
	
	@PF=@PR=();
	my ($gene_5p, $gene_3p)=find_gene($seq);
	
	sbarprint("\nMoving forward and reverse into arrays ...");
	my ($dnaseq_f, $dnaseq_f_len, $dnaseq_r_top, $dnaseq_r, $dnaseq_r_len) = getseq($seq);
	return if ($cancel==1);
	
	# Check that there was actually a DNA sequence found!
	if (length($dnaseq_f)==0) {
		sbarprint("\nNo DNA sequence");
		$cancel = 1;
		return;
	}

	read_windows(\&primer_window, $dnaseq_f, $dnaseq_f_len, $dnaseq_r);
	return if ($cancel==1);

	sbarprint("\nCalculating amplicons ...");
	calc_amplicon($dnaseq_r_len);

	@primer_pairs_pr_s = @primer_pairs;	
	sort_primers('10');
	
	sbarprint("\nFinished ... found ".($#primer_pairs+1)." primer pairs");	
	# draw_dna(\$prf{primer_canvas},\$prf{seq});
	draw_dna();
}


sub get_seq_primers {
	my ($min_range, $max_range, $slist) = get_variables(qw(min_range max_range primers)); 

	$cancel=0;
	$bs=0;
	my $seq = get_seq();
	# check_range();
	
	# unless boundaries are set, set them to sequence limits	
	unless (($$max_range)&&($$max_range)) {
		$$min_range = 0;
		$$max_range = length($seq);
	}
	
	$$max_range = length($seq) if $$max_range > length($seq);
	
	@PF=@PR=();
	my ($gene_5p, $gene_3p)=find_gene($seq);
	
	sbarprint("\nMoving forward and reverse into arrays ...");
	my ($dnaseq_f, $dnaseq_f_len) = getseq($seq);
	return if ($cancel==1);
	
	# Check that there was actually a DNA sequence found!
	if (length($dnaseq_f)==0) {
		sbarprint("\nNo DNA sequence");
		return;
	}

	read_windows(\&primer_window, $dnaseq_f, $dnaseq_f_len);
	return if ($cancel==1);

	sbarprint("\nFinding primers ...");
	calc_seq_primers();

	@$slist = @primer_pairs;	
	sort_primers('1');
	
	sbarprint("\nFinished ... found ".($#primer_pairs+1)." primer pairs");	
	# draw_dna(\$prf{primer_canvas},\$prf{seq});
	draw_dna();
}


sub get_qprimers {
	$cancel=0;
	$bs=0;
		
	my $mrna_seq;
	($spidey_out, $mrna_seq) = run_spidey(1);
	return if $cancel == 1;
			
	@intron_exon_bounds = ();
	while ($spidey_out =~ /(\d{1,})-(\d{1,}) \(mRNA\)/g) {
		push (@intron_exon_bounds, $2);
	}
	
	# last boundary is end of mRNA sequence: remove
	pop @intron_exon_bounds;
			
	@PF=@PR=();	
	my ($gene_5p, $gene_3p)=find_gene($mrna_seq);
	
	sbarprint("\nMoving forward and reverse into arrays ...");
	my ($dnaseq_f, $dnaseq_f_len, $dnaseq_r_top, $dnaseq_r, $dnaseq_r_len) = getseq($mrna_seq);
	return if ($cancel==1);
	
	# Check that there was actually a DNA sequence found!
	if (length($dnaseq_f)==0) {
		sbarprint("\nNo DNA sequence");
		return;
	}
	
	read_windows(\&primer_window, $dnaseq_f, $dnaseq_f_len, $dnaseq_r);			
	return if ($cancel==1);
	
	sbarprint("\nCalculating amplicons ...");
	# set qpcr flag for amplicon and draw_dna subroutines:
	$qpcr_flag=1;
	calc_amplicon($dnaseq_r_len);
	
	@primer_pairs_q_s = @primer_pairs;
	sort_primers('10');
		
	sbarprint("\nFinished ... found ".($#primer_pairs+1)." primer pairs");	
	# draw_dna(\$prf{qprimer_canvas},\$prf{mrna_seq});
	draw_dna();
	
	# unset flag
	$qpcr_flag=0;
}


sub find_re_sites {
	my $page = which_nb_page();
	unless ($page eq "pd") {
		dialogue("This feature is only relevant for the Standard PCR page");
		return;
	}

	my $seq = get_seq();
	my ($max_range_5p, $max_range_3p) 
		= get_variables(qw(max_range_5p max_range_3p));
	$seq = substr($seq, $$max_range_5p, $$max_range_3p-$$max_range_5p);
	
	# find the rebase file (different versions are released constantly)
	my @gcg_paths = glob("$HOME"."gcg.*");
	if (@gcg_paths == 0) {
		# search for the file in the program directory
		my $program_directory = $0;
		$program_directory =~ s/([^\\\/]*$)//;
		@gcg_paths = glob("$program_directory"."gcg.*");
	}
	
	my $gcg_path;
	foreach (@gcg_paths) {
		if (/.*gcg\.\d*$/) {
			$gcg_path = $_;
		}
	}
	
	unless ($gcg_path) {
		dialogue("Error: Cannot find the Restriction Enzyme database.\n\nYou need the GCG database from REBASE (http://rebase.neb.com/rebase/rebase.html) for this function to work.  The file 'gcg.###' (where ### is the version) should be included in the PerlPrimer distribution - please place that file either in the directory where the program is located, or in the directory $HOME, and try again.\n\nAlternatively, you can download the latest version of the database from the REBASE ftp site.");
		return;
	}
	
	# open rebase
	unless (open (REBASE, "<$gcg_path")) {
		dialogue("Error: Cannot open the Restriction Enzyme database: $!");
		return;
	}		
	my @re_data = <REBASE>;
	close REBASE;
	
	# parse data
	my $i;
	my %enzymes;
	my %enzyme_sites;
	my %enzymes_full;
	foreach (@re_data) {
		s/[\n\r]//g;
		my ($enzyme, $re_site) = /^\;*([A-Z|a-z|0-9|\;]+)\s+\d+\s+([A-z\'_]+)\s+/;
		
		next unless $enzyme && $re_site;
		my $original_site = $re_site;
		$re_site =~ s/\'//g;
		$re_site =~ s/_//g;
		
		if ($enzymes_full{$re_site}) {
			$enzymes_full{$re_site}.=" | $enzyme";
		} else {
			$enzymes_full{$re_site}=$enzyme;
		}
		
		$enzyme_sites{$enzyme}=$original_site;
				
		if ($simple_sites) {
			next if $re_site =~ /[rymkswbdhvn]/i;
			next if length($re_site) > 6;
			next if length($re_site) < 6;
		} else {
			$re_site =~ s/r/[ga]/ig;
			$re_site =~ s/y/[ct]/ig;
			$re_site =~ s/m/[ac]/ig;
			$re_site =~ s/k/[gt]/ig;
			$re_site =~ s/s/[gc]/ig;
			$re_site =~ s/w/[at]/ig;
			$re_site =~ s/b/[cgt]/ig;
			$re_site =~ s/d/[agt]/ig;
			$re_site =~ s/h/[act]/ig;
			$re_site =~ s/v/[acg]/ig;
			$re_site =~ s/n/[acgt]/ig;
		}
		
		$enzymes{$enzyme}=$re_site;
	}
		
	my @no_matches;
	if ($exclude_found_sites) {
		foreach my $key (keys(%enzymes)) {
			if ($seq !~ /$enzymes{$key}/i) {
				push @no_matches, "$key - $enzyme_sites{$key}";
			}
		}
	} else {
		foreach my $key (keys(%enzymes)) {
			push @no_matches, "$key - $enzyme_sites{$key}";
		}
	}
		
	unless (Exists($cloning_d)) {
		$cloning_d = $top->Toplevel(-title=>'Cloning site configuration');
		my $cloning_f = $cloning_d->Frame()->pack(-fill=>'both', -pady=>7);
		my $cloning_fb = $cloning_d->Frame()->pack(-side=>'bottom', -fill=>'none');
		$gui_count=0;
		
		my $re_enzymes = [sort @no_matches];
		
		# deconstruct current info if it exists
		my ($tmp_site, $enzyme, $key);
		if ($primer_seq_5f) {
			$_ = uc($primer_seq_5f);
			($cloning_anchor, $tmp_site) = /([A-Z]*)_([A-Z]*)/;
			if ($enzymes_full{$tmp_site}) {
				$enzyme = $enzymes_full{$tmp_site};
				($key) = ($enzyme =~ /([A-z0-9]*)/); 
				$forward_re_site = "$enzyme - $enzyme_sites{$key}";
			} else {
				$forward_re_site = $tmp_site;
			}
		}
		
		if ($primer_seq_5r) {
			$_ = uc($primer_seq_5r);
			my $tmp_site;
			($cloning_anchor, $tmp_site) = /([A-Z]*)_([A-Z]*)/;
			if ($enzymes_full{$tmp_site}) {
				$enzyme = $enzymes_full{$tmp_site};
				($key) = ($enzyme =~ /([A-z0-9]*)/);
				$reverse_re_site = "$enzyme - $enzyme_sites{$key}";
			} else {
				$reverse_re_site = $tmp_site;
			}
		}
				
		nr(\$cloning_f);		
			pack_gui('Forward restriction enzyme site: ', \$forward_re_site, 'cloning_re_f', \$cloning_f, 'be', 0, $re_enzymes);
		nr();
			pack_gui('Reverse restriction enzyme site: ', \$reverse_re_site, 'cloning_re_r', \$cloning_f, 'be', 0, $re_enzymes);
		nr();
			pack_gui("5' anchor sequence", $cloning_anchor, "cloning_anchor", \$cloning_f, 15);
		
		pack_gui('OK', sub {
				my ($enzyme) = ($forward_re_site =~ /([A-z0-9]*)/); 
				$primer_seq_5f = lc("$cloning_anchor\_$enzymes{$enzyme}");
				($enzyme) = ($reverse_re_site =~ /([A-z0-9]*)/);
				$primer_seq_5r = lc("$cloning_anchor\_$enzymes{$enzyme}");
				$cloning_d->destroy;
				return;
			}, 'cloning_ok', \$cloning_fb, "b", "active", 0);
		pack_gui('Cancel', sub {$cloning_d->destroy;}, 'cloning_cancel', \$cloning_fb, "b", "normal", 1);
		
		$cloning_d->Icon(-image => $pixmap);
	} else {
		$cloning_d->deiconify;
		$cloning_d->raise;
	}
	
}


sub run_spidey {
	my $mrna_seq = $prf{"qmrna_seq"}->get(0.1,"end");
	my $dna_seq = $prf{"qdna_seq"}->get(0.1,"end");
	
	# Spidey will only accept input as files, so we need to move the two
	# sequences to temporary files
	### Only they're not temporary at present!  Really should delete them at the end ...
	### or at the very least, save them to a user-definable temp directory (/tmp or c:\tmp)
	### perhaps run a check on a default temp directory's existence, and otherwise use $HOME
	unless (open (MRNA, ">".$file_pre.".mrna_tmp")) {
		dialogue("Error: could not write mRNA temp file: $!\n");
		$cancel=1;
		return;
	}
	print MRNA $mrna_seq;
	close (MRNA);
	
	unless (open (DNA, ">".$file_pre.".dna_tmp")) {
		dialogue("Error: could not write DNA temp file: $!\n");
		$cancel=1;
		return;
	}
	print DNA $dna_seq;
	close (DNA);

	# run spidey ...
	sbarprint("\nRunning spidey ...");
	$_ = `$spidey_command -p $_[0]`;
	
	# abort and complain if we don't see the --SPIDEY signature
	unless (/--SPIDEY/) {
		dialogue("Cannot run Spidey executable\n(check that $spidey_exec exists?)");
		sbarprint("\nCancelled - Spidey executable not found");
		$cancel=1;
		return;
	}
	
	# Warn if problems occur ...
	print "\n\nWarning: incomplete sequence data or large intron present\n" unless /overall percent identity: 100\.0%/;
	print "\n\nWarning: mRNA coverage incomplete\n" unless /mRNA coverage: 100%/;
	
	# ... and try using the large intron option if they have
	unless (/overall percent identity: 100\.0%/ && /mRNA coverage: 100%/) {
		print "Trying spidey again using large intron option ...\n";
		$_ = `$spidey_command -p $_[0] -X`;
		
		# if things are still crap, display a dialogue with spidey's output,
		# so the user can work out what's gone wrong
		dialogue("Alignment error:\nSequence data is incomplete\n\n$_") unless /overall percent identity: 100\.0%/;
		dialogue("Alignment error:\nmRNA coverage is incomplete\n\n$_") unless /mRNA coverage: 100%/;
	}
	
	# delete those tmp files
	unlink("$file_pre.mrna_tmp") || dialogue("Warning: could not delete temporary file $file_pre.mrna_tmp");
	unlink("$file_pre.dna_tmp") || dialogue("Warning: could not delete temporary file $file_pre.dna_tmp");
	return ($_, $mrna_seq);
}


sub get_bisulphite {
	# Criteria:
	# Nested or heminested PCR
	# Primer amplicon not greater than 450bp
	# Selection for bisulphite-converted templates
	# Even distribution of bases prior to conversion
	# As high GA content as posible [currently not used]
	# Lack of CpG residues or degeneracy
	# Long (25-30mer) primer
	$cancel=0;
	
	# old benchmarking code - used for optimising routines
	my $t_old = new Benchmark if $benchmark;
	$bs = 1;
	check_range();
	
	my ($max_ampsize, $max_range_5p, $max_range_3p) 
		= get_variables(qw(max_ampsize max_range_5p max_range_3p));
	
	# set the max_ampsize in case the user has modified
	# either max_range_5p or max_range_3p by hand ...
	if (($$max_range_5p)&&($$max_range_3p)) {
		$$max_ampsize = $$max_range_3p-$$max_range_5p
	}

	my $seq = $prf{"bisul_seq"}->get(0.1,"end");
	my ($dnaseq_f, $dnaseq_f_len, $dnaseq_r_top, $dnaseq_r, $dnaseq_r_len) = getseq($seq);
	
	# Check that there was actually a DNA sequence found!
	if (length($dnaseq_f)==0) {
		sbarprint("\nNo DNA sequence");
		return;
	}

	@PF=@PR=();
	return if ($cancel==1);
	
	# Only use top strand for primer design (because conversion makes
	# the two templates non-complementary), so use $dnaseq_r_top
	read_windows(\&bisul_window, $dnaseq_f, $dnaseq_f_len, $dnaseq_r_top);
	return if ($cancel==1);

	sbarprint("\nCalculating amplicons ...");
	calc_amplicon($dnaseq_r_len);
	
	@primer_pairs_bs_s = @primer_pairs;
	sort_primers('10');
		
	$bs = 0;
	
	# benchmarking code
	if ($benchmark) {
		my $t_new = new Benchmark;
		my $diff = timediff ($t_new, $t_old);
		my $str = timestr ($diff, "all", "5.3f");
		print "Time for bisulphite was $str\n\n";
	}
}

sub check_range {
	my ($max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p) 
		= get_variables(qw(max_ampsize max_range_5p min_range max_range max_range_3p)); 

	$$max_range_5p = undef if defined($$max_range_5p) && $$max_range_5p eq "";
	$$max_range_3p = undef if defined($$max_range_3p) && $$max_range_3p eq "";
	unless (defined($$max_range_5p)&&defined($$max_range_3p)) {
		return unless defined($$min_range)&&defined($$max_range);
		$$max_range_5p = $$min_range unless $$max_range_5p;
		$$max_range_3p = $$max_range unless $$max_range_3p;
		$$max_ampsize = $$max_range_3p - $$max_range_5p;
	}
}

sub copy_selected_primers {
	# Get variable refs
	my ($ref, $slist) = get_variables(qw(hlist primers));
	my $page = which_nb_page();
	
	# get selected primers
	my @sel = $$ref->selectionGet;
	
	# Clipboard copying routine	
	my $clip = "Forward Primer\tPos\tLen\tTm\tReverse Primer\tPos\tLen\tTm\tAmp\tdG\n";
	
	my (@gene_array, $gene_frame, $seq);
	if ($page eq "pd" && (($primer_seq_5f) || ($primer_seq_5r))) {
		$seq = get_seq();
		@gene_array=find_gene($seq);
		$gene_frame=$gene_array[0][0]%3;
	}
	
	foreach (@sel) {
		my $sel = $_;
		for my $j ( 0 .. 4, 8, 6, 7, 9, 10 ) {
	    		if ($page eq "bis") {
	    			# Primer redundancy for CpG residues:
				# Replaces T with Y (pyrimidine) for forward
				# Replaces A with R (purine) for reverse
		    		
				my $temp_seq;
				if ($j ==  0) {
					$temp_seq = $$slist[$sel][0];
					my $f_original = $$slist[$sel][11];
                       			for my $i (0 .. length($f_original)) {
                             				if (substr($f_original, $i, 2) eq "CG") {
                      	       				substr($temp_seq, $i, 1) = 'Y';
                             				} 
                     				}
	       				$clip .= "$temp_seq\t";
				} elsif ($j == 4) {
					$temp_seq = $$slist[$sel][4];
					my $r_original = $$slist[$sel][12];
                       			for my $i (0 .. length($r_original)) {
      	                       			if (substr($r_original, $i, 2) eq "CG") {
      	                	       			substr($temp_seq, $i, 1) = 'R';
      	                       			} 
      	               			}
		        		$clip .= "$temp_seq\t";
				} else {
					$clip .= "$$slist[$sel][$j]\t";
				}
			} elsif ($page eq "pd") {
				if ($j == 0 && $primer_seq_5f) {
					my $fprimer = $$slist[$sel][$j];
					my ($primer_seq_5f_real, $insert_f) = add_cloning($seq, $fprimer, $$slist[$sel][1]);
					$clip .= uc("$primer_seq_5f_real$insert_f\_$fprimer\t");
				} elsif ($j == 4 && $primer_seq_5r) {
					my $rprimer = $$slist[$sel][$j];
					my ($primer_seq_5r_real, $insert_r) = (add_cloning($seq, '', '', $rprimer, $$slist[$sel][8]))[5,6];
					$clip .= uc("$primer_seq_5r_real$insert_r\_$rprimer\t");
				} else {
					$clip .= "$$slist[$sel][$j]\t";
				}
			} else {
				$clip .= "$$slist[$sel][$j]\t";
			}
		}
		$clip .= "\n";
	}
	
	$top->clipboardClear;
	$top->clipboardAppend($clip);
}

sub select_all_primers {
	my ($ref, $slist) = get_variables(qw(hlist primers));
	$$ref->selectionSet(0,$#$slist);
}

sub primer_take_range {
	# Get variable refs
	my ($max_range_5p, $min_range, $max_range, $max_range_3p, $hlist, $slist) 
		= get_variables(qw(max_range_5p min_range max_range max_range_3p hlist primers)); 

	# get selected primer
	my @sel = $$hlist->selectionGet;
	my $primer = shift @sel;
	
	$$max_range_5p = $$min_range = $$slist[$primer][1]+$$slist[$primer][2];
	$$max_range_3p = $$max_range = $$slist[$primer][8]-$$slist[$primer][6];
	
	draw_dna();
	browse_primer($primer);
}

sub menu_popup {
	my ($popup_ref) = get_variables(qw(popup));
	$$popup_ref->Popup(
			-popanchor  => 'nw',
			-popover => 'cursor');
}

sub menu_text {
	$text_widget_ref = $_[0];
	$popup_text->Popup(
			-popanchor  => 'nw',
			-popover => 'cursor');
}

sub sort_primers {
	my ($ref, $slist) = get_variables(qw(hlist primers));
	my @old_slist = @$slist;
	
	$$ref->delete('all');
	
	# Re-sort forwards or backwards
	# allow the ability to simply redisplay the list (when opening files)
	if ($_[0]) {
		$_[0] = $sort_prev if $_[0] eq "r";
		
		# We want dG to be sorted from highest to lowest by default:
		# this is really inelegant!
		my $sort_order = $sort_reverse;
		$sort_order = 1-$sort_reverse if ($_[0] eq '10');
		
		@$slist = sort {@$a[$_[0]] <=> @$b[$_[0]]} @old_slist if $sort_order==0;
		@$slist = sort {@$b[$_[0]] <=> @$a[$_[0]]} @old_slist if $sort_order==1;
	}
	
	my $num_elements = @{ @$slist[0] } if @$slist;
	$num_elements ||=0;
	
	# Re-draw
	if ($num_elements == 5) { # sequencing quick fix!!
		for my $i ( 0 .. $#$slist ) {
			$$ref->add($i);
			$$ref->itemCreate($i, 0, -text=>"$$slist[$i][0]", -style=> $style_primer);
			$$ref->itemCreate($i, 1, -text=>"$$slist[$i][1]");
			$$ref->itemCreate($i, 2, -text=>"$$slist[$i][2]");
			$$ref->itemCreate($i, 3, -text=>"$$slist[$i][3]", -style=> $style_tm);
			
			$$ref->itemCreate($i, 4, -text=>"$$slist[$i][4]");
		}
	} else {
		for my $i ( 0 .. $#$slist ) {
			$$ref->add($i);
			$$ref->itemCreate($i, 0, -text=>"$$slist[$i][0]", -style=> $style_primer);
			$$ref->itemCreate($i, 1, -text=>"$$slist[$i][1]");
			$$ref->itemCreate($i, 2, -text=>"$$slist[$i][2]");
			$$ref->itemCreate($i, 3, -text=>"$$slist[$i][3]", -style=> $style_tm);
			
			$$ref->itemCreate($i, 4, -text=>"$$slist[$i][4]", -style=> $style_primer);
			$$ref->itemCreate($i, 5, -text=>"$$slist[$i][8]");
			$$ref->itemCreate($i, 6, -text=>"$$slist[$i][6]");
			$$ref->itemCreate($i, 7, -text=>"$$slist[$i][7]", -style=> $style_tm);
			
			$$ref->itemCreate($i, 8, -text=>"$$slist[$i][9]");
			$$ref->itemCreate($i, 9, -text=>"$$slist[$i][10]");
		}
	}
	
	# save last value so that "reverse" can be applied if neccessary
	$sort_prev = $_[0];
}

sub cancel {
	# Set cancel flag
	# (this stops the primer-pairs calculation,
	# which can go on for some time ...)
	$cancel = 1;
	sbarprint("\nOperation cancelled");
}

sub get_primers_auto_in {
	my ($sub_ref, $primers_ref) = get_variables(qw(primer_sub primers));
	&$sub_ref();
	while ($#$primers_ref < 0) {
		step_in();
		return if $cancel==1;
		&$sub_ref();
	}
}

sub get_primers_auto_out {
	my ($sub_ref, $primers_ref) = get_variables(qw(primer_sub primers));
	&$sub_ref();
	while ($#$primers_ref < 0) {
		step_out();
		return if $cancel==1;
		&$sub_ref();
	}
}

sub get_primers_cloning {
	my $page = which_nb_page();
	unless ($page eq 'pd') {
		dialogue("This feature is only relevant for the Standard PCR page");
		return;
	}
	my $seq = get_seq();
	unless ($seq) {
		dialogue("No DNA sequence found!");
		return;
	}
	
	my ($sub_ref, $primers_ref) = get_variables(qw(primer_sub primers));
	&$sub_ref();
	while ($#$primers_ref < 0) {
		step_cloning();
		return if $cancel==1;
		&$sub_ref();
	}
}

sub new_file {
	my $page = which_nb_page();
	if ($page eq "primer") {
		dialogue("You can't open a new file from the primer information page - please switch to the relevant project page first");
		return;
	}

	my ($seq_ref, $hlist) = get_variables(qw(seq hlist));
	my $seq = $$seq_ref->search(-regexp => "[a|g|t|c|A|T|G|C]", "0.1");
	
	# If there's a sequence present, prompt to make sure the user wants to
	# lose all changes ...
	if (defined($seq)) {
		my $answer = dialogue("Data will be lost!", "OKCancel");
		if ($answer =~ /cancel/i) {
			return;
		}
	}
	
	if ($page eq "qpcr") {
		# only for qpcr
		$prf{qdna_seq}->delete(0.1,"end");
	}
	
	# Clear variables
	foreach my $key (keys %{ $variables{$page}}) {
		my $pointer = $variables{$page}{$key};
		if (eval {$$pointer}) {
			undef($$pointer);
			if ($default_variables{$page}{$pointer}) {
				$$pointer = $default_variables{$page}{$pointer};
			}
		}
	}
	
	# Clear arrays
	foreach my $key (keys %{ $arrays{$page}}) {
		my $pointer = $arrays{$page}{$key};
		if (eval {@$pointer}) {
			undef(@$pointer);
		}
	}
	
	# Clear sequence and hlist displays		
	$$seq_ref->delete(0.1,"end");
	$$hlist->delete('all');
	draw_dna();
	
	# Clear title
	$top->configure(-title=>"PerlPrimer v$version");
	$open_file{$page} = 'File not saved';
}

sub save_defaults {
	# Save user-defined default values
	my $page = which_nb_page();
	my $file_defaults = "$HOME.perlprimer.$page";
	my $defaults;
	# get list of variables and values
	foreach my $key (keys %{ $variables{$page}}) {
		my $pointer = $variables{$page}{$key};
		if (eval {defined($$pointer)}) {
			if (defined($default_variables{$page}{$pointer})) {
				$default_variables{$page}{$pointer} = $$pointer;
				$defaults .= "$key = $$pointer\n";
			}
		}
	}
	unless (open (DEFAULTS, ">$file_defaults")) {
		dialogue("Error: could not open $file_defaults for writing: $!");
		return;
	}
	print DEFAULTS $defaults;
	close DEFAULTS;
}

sub restore_defaults {
	# Remove user-specified preferences file if exists
	my $page = which_nb_page();
	my $file_defaults = "$HOME.perlprimer.$page";
	
	if (-e $file_defaults) {
		my $answer = dialogue("Warning: This will permanently remove the user-defined default values for this page and restore the built-in values (requires program restart)",'OkCancel');
		return if $answer =~ /cancel/i;
		unlink($file_defaults) || dialogue("Error: cannot delete file $file_defaults: $1");
		dialogue("Built-in values restored.  Please restart the program for changes to take effect");
	} else {
		dialogue("Error - User-defined defaults have not been saved for this page.\n\nIf you are trying to reset the parameters for this tab and remove all sequence data, use the \"New File\" command from the File menu");
	}
}

sub load_defaults {
	# load user-defined default values (same idea as save defaults above, but in reverse)
	foreach my $key (keys (%nb_page_ref)) {
		my $page = $nb_page_ref{$key};
		my $file_defaults = "$HOME.perlprimer.$nb_page_ref{$key}";
		if (-e $file_defaults) {
			print "Loading user-defined defaults file: $file_defaults\n";
			unless (open (DEFAULTS, "<$file_defaults")) {
				dialogue("Error: could not open $file_defaults for reading: $!");
				return;
			}
			while (<DEFAULTS>) {
				chomp;
				my ($variable, $value) = split / = /;
				my $pointer = $variables{$page}{$variable};
				if (defined($default_variables{$page}{$pointer})) {
					$$pointer = $value;
					$default_variables{$page}{$pointer} = $value;
				}
			}				
			close DEFAULTS;
		}
	}
}

#----------------#
# Blast routines #
#----------------#

sub blast {	
	# for more info on using BLAST via http requests, see:
	# http://www.ncbi.nlm.nih.gov/BLAST/Doc/
	# which details the request syntax
	
	unless (Exists($blast_d)) {
		# this routine should not continue running if the user has destroyed
		# the blast search dialogue.  Hopefully this will stop it!
		$top->afterCancel($rptid);
		return;
	}
	
	return unless (my $query=shift);
	my $blast_put = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?QUERY=$query&DATABASE=$blast_database&ENTREZ_QUERY=($blast_entrez_query)&EXPECT=$blast_expect&FORMAT_TYPE=Text&PROGRAM=blastn&SERVICE=plain&CMD=Put";		
			
	# send off the query
	$_ = http_get($blast_put);
	
	# abort if unsuccessfull
	return if $_ eq " ";
	
	# find the RID string
	/QBlastInfoBegin(.*)QBlastInfoEnd/sm;
	$_ = $1;
	/RID = ([\d\w\-.]*)/;

	$rid_get = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?RID=$1&FORMAT_TYPE=Text&CMD=Get";
	
	# Hah! How about this for a simple way to unblock your GUI ... the only way I
	# know of to use non-blocking, sleeping subroutines without forks - which we
	# can't use with ActivePerl for Win32 :(	
	$flag=undef;
	$blast_count = 0;
	$rptid = $blast_d->repeat(15000, \&blast_wait);
	$blast_d->waitVariable(\$flag);
	$blast_d->afterCancel($rptid);
	
	# Having asked for text output, blast still provides some HTML formatting!
	s/\<!\-\-.*\-\-\>//sg;
	s/\<[\w\\]*\>//g;
			
	# OK, let's try something really nifty ... let's pull each entry out:
	my $blast_out = $_;
	undef(@blast_results_1);
	undef(@blast_results_2);
	undef(@blast_results_3);
	
	/(BLASTN.*?Value.*?)(\w+.*?)ALIGNMENTS/s;
	my ($blast_header, $blast_summary) = ($1,$2);
	
	# search for each sequence identifier and find corresponding blast entry
	while (/^(\w*\|[\w\|\.\d]+)/mg) {
			my $found_result = $1;
			if ($blast_out =~ /\n\>(\Q$found_result\E.*?Length = \d*)[\n\s]*(Score = .*?Strand = \w* \/ \w*)[\n\s]*(Query.*?Sbjct.*?)\n/s) {
				push @blast_results_1, $1;
				push @blast_results_2, $2;
				push @blast_results_3, $3;
			}
			
	}
		
	return ($blast_header, $blast_summary);
}

sub blast_wait {
	#print "Waiting for server response ...\n";
	$blast_count++;
	my $seconds = $blast_count*15;
	my $seconds_remainder = $seconds % 60;
	my $minutes = ($seconds-$seconds_remainder)/60;
	$seconds_remainder = "00" if $seconds_remainder == 0;
	
	$blast_status = "Waiting for server response ... $minutes:$seconds_remainder elapsed";
	$top->update;

	# print "status: $blast_status\n";
	$_ = http_get($rid_get);
	unless (/Status=WAITING/) {
		$flag=0;
	} else {
		#### debugging ... really need to use a timeout with repeat above ...
		# print "output was $_\n";
	}
}


#------------------#
# Ensembl routines #
#------------------#

sub get_ensembl {
	return if check_packages("HTTP::Request", "LWP::UserAgent");
	# retrieve a gene sequence (or a gene and genomic sequence) from ensembl.org ...
	
	# Why do we use Ensembl?  Simply because it happens to provide easy sequence 
	# retrieval, including the genomic sequence of a gene (i.e. not just the cDNA
	# sequence).  However, unlike NCBI, Ensembl provides no documentation on their
	# cgi scripts; this was all worked out from going through the page sources.
	my $nb_page = which_nb_page();
	if ($nb_page eq "primer") {
		dialogue("Please switch to the project page that you wish to enter the Ensembl data into");
		return;
	}

	unless (Exists($ensembl)) {
		$ensembl = $top->Toplevel(-title=>'Retrieve gene from Ensembl');
		my $ensembl_f = $ensembl->Frame()->pack(-fill=>'both', -pady=>7);
		my $ensembl_fb = $ensembl->Frame()->pack(-side=>'bottom', -fill=>'none');
		$gui_count=0;
		
		nr(\$ensembl_f);		
			pack_gui('Gene name ', $ensembl_gene, "ensembl_gene", \$ensembl_f, 15);
		nr();
			pack_gui('Organism', \$ensembl_organism, 'ensembl_organsim', \$ensembl_f, 'be', 0, \@ensembl_species, 20);
		nr();
			pack_gui('Retrieve', \$ensembl_type, 'ensembl_type', \$ensembl_f, 'be', 0, \@ensembl_types, 10);
		
		pack_gui('OK', \&fetch_ensembl, 'ensembl_ok', \$ensembl_fb, "b", "active", 0);
		pack_gui('Cancel', sub {$ensembl->destroy;}, 'ensembl_cancel', \$ensembl_fb, "b", "normal", 1);
		
		$ensembl->Icon(-image => $pixmap);
	} else {
		$ensembl->deiconify;
		$ensembl->raise;
	}
}

sub fetch_ensembl {
	my $page = which_nb_page();
	unless (($ensembl_gene) && ($ensembl_organism) && ($ensembl_type)) {
		dialogue("Please make sure gene name, organism and retrieval type are filled in");
		return;
	}
	
	# Search for the gene:
	
	# Ensembl's textview cgi will work regardless of the species origin given in the
	# address - the Homo_sapiens in the address will not limit the search.
	$_ = http_get("http://www.ensembl.org/Homo_sapiens/textview?species=$ensembl_organism&idx=Gene&q=$ensembl_gene");
	
	# find the Ensembl gene ID, and count the number - if there's more than one
	# we'll have to ask the user to be more specific
	my $gene_id;
	my @gene_names;
	my $name;
	my %ids;
	my $count=0;
	while (/Ensembl Gene: \<\/[bB]\>\<[aA].*?\>(.*?)\<.*?\<br \/\>(.*?)\.\s*\[/mg) {
		$gene_id = $1;
		$name = $2;
		$name =~ s/\<.*?\>//g;
		push @gene_names, $name;
		$ids{$name}=$gene_id;
		$count++;
	}
	
	if ($count>1) {
		# multiple/ambiguous matches
		@gene_names = sort(@gene_names);
		my $ensembl_mm = $top->Toplevel(-title=>"More than one match found ...");
		my $ensembl_mm_f = $ensembl_mm->Frame()->pack(-fill=>'both', -pady=>7);
		my $ensembl_mm_fb = $ensembl_mm->Frame()->pack(-side=>'bottom', -fill=>'none');
		nr(\$ensembl_mm_f);		
			pack_gui('Please select gene of interest','',"ensemble_mm_d_note",'', 'lt');
		nr();
			pack_gui('', \$name, 'ensembl_mm_d_genes', '', 'be', 0, \@gene_names, 40);
		
		my $cancel=1;
		pack_gui('OK', sub {
				$cancel=undef;
				$ensembl_mm->destroy;
				$gene_id = $ids{$name};
			}, 'ensembl_ok', \$ensembl_mm_fb, "b", "active", 0);
		pack_gui('Cancel', sub {
				$ensembl_mm->destroy;
			}, 'ensembl_cancel', \$ensembl_mm_fb, "b", "normal", 1);
		
		$ensembl_mm->Icon(-image => $pixmap);
		
		# this seems a bit messy, but at least it works
		# (we need it to freeze execution at this point, since the user may
		# wish to cancel and refine their choice)
		$ensembl_mm->waitWindow;
		return if $cancel;

		# dialogue("More than one match found:\n$multiple\nPlease enter a more specific query");
		# return;
	}
		
	unless ($gene_id) {
		# no matches
		if (/Your query found no matches/) {
			dialogue("Your query found no matches");
		} elsif ($_ eq " ") {
			# returned if response->is_error below
			return;
		} else {
			dialogue("Unable to find gene_id in response from server");
		}
		
		return;
	}
	# print "gene_id was $gene_id\n";
	
	if ($page eq 'qpcr') {
		# retrieve both gene and transcript sequences - retrieval type is ignored
		$_ = http_get("http://www.ensembl.org/$ensembl_organism/exportview?tab=fasta&embl_format=fasta&type=feature&ftype=gene&id=$gene_id&fasta_option=genomic&out=text&btnsubmit=Export&embl_format=embl");
		$prf{qdna_seq}->delete(0.1,"end");
		$prf{qdna_seq}->insert(0.1,$_);
				
		$_ = http_get("http://www.ensembl.org/$ensembl_organism/exportview?tab=fasta&embl_format=fasta&type=feature&ftype=gene&id=$gene_id&fasta_option=cdna&out=text&btnsubmit=Export&embl_format=embl");
		$prf{qmrna_seq}->delete(0.1,"end");
		$prf{qmrna_seq}->insert(0.1,$_);
	} else {
		# retrieve requested sequence
		my ($seq_ref) = get_variables('seq');
			
		$_ = http_get("http://www.ensembl.org/$ensembl_organism/exportview?tab=fasta&embl_format=fasta&type=feature&ftype=gene&id=$gene_id&fasta_option=$ensembl_type&out=text&btnsubmit=Export&embl_format=embl");
		$$seq_ref->delete(0.1,"end");
		$$seq_ref->insert(0.1,$_);
	}
	
	# update title
	$top->configure(-title=>"PerlPrimer v$version - $ensembl_gene");
	$open_file{$page} = $ensembl_gene;
	
	# destroy dialogue, draw the sequence
	$ensembl->destroy;
	draw_dna();
	reset_bounds();
	return;	
}


#------------------------#
# http retrieval routine #
#------------------------#

sub http_get {
	return if check_packages("HTTP::Request", "LWP::UserAgent");
	# simple html retrieval ... used by BLAST and Ensembl routines
	# I'd love to fork this and stop the GUI blocking here - but I cannot find
	# an effective way to do this ... :(
	my ($address, $method) = @_;
	my $ua = LWP::UserAgent->new();
	
	# timeout - default value of 180 seconds is too long
	$ua->timeout(60);
	
	# proxy server
	if ($use_proxy) {
		# allow for users using "http://" in the proxy address
		$http_proxy =~ /(http:\/\/)*(.*)/;
		my $proxy = $2;
		$ua->proxy('http', "http://$proxy:$http_port");
	} else {
		# still allow env_proxy override
		$ua->env_proxy();
	}
	
	my $req = HTTP::Request->new(GET => "$address");

	my $response = $ua->request($req);
	if ($response->is_error) {
		dialogue("HTTP request unsuccessful");
		return " ";
	}
	
	$top->update;
	my $http = $response->content;
	# print "...\n";
	
	return $http;
}


#---------------------#
# Status bar updating # 
#---------------------#

sub sbarprint {
	print SBAR $_[0];
	$top->update;
}


#------------#
# HList subs #
#------------#

sub header_create {
	my ($ref, $widget_name) = @_;	
	unless ($widget_name eq 'seq_res') {
		$$ref->header('create', 0, -text => 'Forward Primer');
		$$ref->header('create', 1, -text => 'Pos');
		$$ref->header('create', 2, -text => 'Len');
		$$ref->header('create', 3, -text => 'Tm');
		$$ref->header('create', 4, -text => 'Reverse Primer');
		$$ref->header('create', 5, -text => 'Pos');
		$$ref->header('create', 6, -text => 'Len');
		$$ref->header('create', 7, -text => 'Tm');
		$$ref->header('create', 8, -text => 'Amp');
		$$ref->header('create', 9, -text => 'Primer-Dimer dG');
	} else {
		$$ref->header('create', 0, -text => 'Forward Primer');
		$$ref->header('create', 1, -text => 'Pos');
		$$ref->header('create', 2, -text => 'Len');
		$$ref->header('create', 3, -text => 'Tm');
		$$ref->header('create', 4, -text => 'Primer-dimer dG');
	}
	
	$style_primer = $$ref->ItemStyle("text",
		# -font => $hlist_font,
		-background => '#eeeeee',
		-foreground => 'red');
	$style_tm = $$ref->ItemStyle("text",
		# -font => $hlist_font,
		-background => '#eeeeee',
		-foreground => 'blue');
}


sub browse_bisulphite {
	my $seq = get_seq();
	# Believe it or not, when a primer-pair is double-clicked to send the user
	# to the hlist-command sub, the browse-command sub is still executed! In
	# our case, that's pretty bad since we're now requesting a sequence on a
	# page that doesn't have one.  So there's a failsafe here:
	return if $seq eq "1";

	my $hlist_sel = $_[0];
	
	# This routine colours original C residues red and original CpG C's blue
	if ($primer_pairs_bs_s[$hlist_sel][10]) {
		# brute force!!		
		my $f_original = $primer_pairs_bs_s[$hlist_sel][11];
		my $f_converted = $primer_pairs_bs_s[$hlist_sel][0];
		
		my $r_original = $primer_pairs_bs_s[$hlist_sel][12];
		my $r_comp_original = complement($r_original);
		my $r_converted = $primer_pairs_bs_s[$hlist_sel][4];
		
		$status_bar->insert('end',"\nF: 5' ", 'grey');
		
		for my $i (0 .. length($f_original)) {
			if (substr($f_original, $i, 2) eq "CG") {
				$status_bar->insert('end',substr($f_converted, $i, 1), 'blue');
			} elsif (compare_bs($f_original, $f_converted, $i)) {
				$status_bar->insert('end',substr($f_converted, $i, 1), 'red')
			} else {
				$status_bar->insert('end',substr($f_converted, $i, 1))
			}
		}
		
		$status_bar->insert('end'," 3'  R: 5' ", 'grey');
		
		for my $i (0 .. length($r_comp_original)) {
			if (substr($r_comp_original, $i, 2) eq "GC") {
				$status_bar->insert('end',substr($r_converted, $i, 1), 'blue');
			} elsif (compare_bs($r_comp_original, $r_converted, $i)) {
				$status_bar->insert('end',substr($r_converted, $i, 1), 'red')
			} else {
				$status_bar->insert('end',substr($r_converted, $i, 1))
			}
		}
		
		$status_bar->insert('end'," 3'", 'grey');
		$status_bar->see('end');

		browse_primer($hlist_sel);
	}
}


sub browse_primer {
	my ($hlist, $slist, $canv) = get_variables(qw(hlist primers canvas));
	my $seq = get_seq();
	return if $seq eq "1";
	
	my $hlist_sel;
		
	if (defined($_[0])) {
		$hlist_sel = $_[0];
	} else {		
		my @sel = $$hlist->selectionGet;
		$hlist_sel = shift @sel;
	}
	
	my $width = $$canv->width;
	my $height= $$canv->height;
	$dna_canvas_size=($width-($dna_canvas_offset*2))/length($seq);
	
	$$canv->delete('primerf','primerr');
	
	my $fprimerpos = $$slist[$hlist_sel][1];
	my $rprimerpos = $$slist[$hlist_sel][8];
	
	my $fprimerposl = $$slist[$hlist_sel][2];
	my $rprimerposl = $$slist[$hlist_sel][6];
	
	my ($fprimer_x1,$fprimer_x2,$rprimer_x1,$rprimer_x2);
	if (defined($fprimerpos)) {
		$fprimer_x1 = $fprimerpos*$dna_canvas_size + $dna_canvas_offset;
		$fprimer_x2 = ($fprimerpos+$fprimerposl)*$dna_canvas_size + $dna_canvas_offset;
	}
	
	if (defined($rprimerpos)) {
		$rprimer_x1 = $rprimerpos*$dna_canvas_size + $dna_canvas_offset;
		$rprimer_x2 = ($rprimerpos-$rprimerposl)*$dna_canvas_size + $dna_canvas_offset;
	}
	
	my $line_arrow_h = 4;
	my $line_arrow_w = 4;
	
	$$canv->createLine($fprimer_x1,$dc_sel_y1,$fprimer_x2,$dc_sel_y1,$fprimer_x2-$line_arrow_w,$dc_sel_y1-$line_arrow_h, -fill=>'red', -width=>2, -tag=>'primerf') if $fprimerpos;
	$$canv->createLine($rprimer_x1,$dc_sel_y2,$rprimer_x2,$dc_sel_y2,$rprimer_x2+$line_arrow_w,$dc_sel_y2+$line_arrow_h, -fill=>'red', -width=>2, -tag=>'primerr') if $rprimerpos;
}


sub compare_bs {
	return 1 if substr($_[0], $_[2], 1) ne substr($_[1], $_[2], 1);
}


sub jump_to_tm {
	# an arrow key shortcut ... jumps from the primer listings to detailed info
	# including primer-dimers
	my $nb_page = which_nb_page();
	return if $nb_page eq "primer";
	$stored_page = $nb_page;
	
	my ($hlist, $slist) = get_variables(qw(hlist primers));
	my @sel = $$hlist->selectionGet;
	my $hlist_sel = $sel[0];
	
	return unless @$slist;
	
	$fprimer = $$slist[$hlist_sel][0];
	$rprimer = $$slist[$hlist_sel][4] if $$slist[$hlist_sel][5];

	get_tm();
	
	$nb->raise('primer');
	$prf{dim}->focus;
}


sub jump_back {
	# jumps back from the primer information page to the project page
	my $nb_page = which_nb_page();
	return unless $nb_page eq "primer";
	return unless $stored_page;
	
	$nb->raise($stored_page);
	my ($hlist) = get_variables(qw(hlist));
	$$hlist->focus;
}


sub hlist_command {
	my $hlist_sel = $_[0];
	my ($hlist_ref, $slist) = get_variables(qw(hlist primers));
	
	my $nb_page = which_nb_page();
	$stored_page = $nb_page;
		
	$fprimer = $$slist[$hlist_sel][0];
	if ($$slist[$hlist_sel][5]) {
		$rprimer = $$slist[$hlist_sel][4];
	} else {
		$rprimer = "";
	}

	get_tm();
	
	$nb->raise('primer');
}


#----------------------------#
# Generate Report subroutine #
#----------------------------#

sub generate_report {
	my ($hlist_ref, $slist) = get_variables(qw(hlist primers));
	my $nb_page = which_nb_page();
	if ($nb_page eq "primer") {
		dialogue("You can't generate a report from the primer information page - please switch to the project page first");
		return;
	}
	
	my ($hlist_sel) = $$hlist_ref->selectionGet;
	unless (defined(@$slist) && defined($hlist_sel)) {
		dialogue("The Generate Report function saves the statistics and alignment of a particular primer pair - please select a primer pair first");
		return;
	}
	
	my $fprimer_mod = $fprimer = $$slist[$hlist_sel][0];
	my $rprimer_mod = $rprimer = $$slist[$hlist_sel][4] unless $nb_page eq 'seq';
	
	my $fprimerpos = $$slist[$hlist_sel][1];
	my $rprimerpos = $$slist[$hlist_sel][8] unless $nb_page eq 'seq';
	my $amplicon_size = $$slist[$hlist_sel][9] unless $nb_page eq 'seq';
	
	# Generate modified primers:
	my (@gene_array, $gene_frame, $seq);
	if ($nb_page eq "pd" && (($primer_seq_5f) || ($primer_seq_5r))) {
		$seq = get_seq();
		@gene_array=find_gene($seq);
		$gene_frame=$gene_array[0][0]%3;
	}
	
	for my $primer ($fprimer, $rprimer) {
    		if ($nb_page eq "bis") {
    			# Primer redundancy for CpG residues:
			# Replaces T with Y (pyrimidine) for forward
			# Replaces A with R (purine) for reverse
    			
			my $f_original = $$slist[$hlist_sel][11];
                       	for my $i (0 .. length($f_original)) {
                             	if (substr($f_original, $i, 2) eq "CG") {
                             		substr($fprimer_mod, $i, 1) = 'Y';
                             	} 
                     	}
			my $r_original = $$slist[$hlist_sel][12];
                       	for my $i (0 .. length($r_original)) {
                             	if (substr($r_original, $i, 2) eq "CG") {
                             		substr($rprimer_mod, $i, 1) = 'R';
                             	} 
                     	}
		} elsif ($nb_page eq "pd") {
			if ($primer_seq_5f) {
				my ($primer_seq_5f_real, $insert_f) = add_cloning($seq, $fprimer, $fprimerpos);
				$fprimer_mod = uc("$primer_seq_5f_real$insert_f\_$fprimer");
			}
			
			if ($primer_seq_5r) {
				my ($primer_seq_5r_real, $insert_r) = (add_cloning($seq, '', '', '', $rprimer, $rprimerpos))[5,6];
				$rprimer_mod = uc("$primer_seq_5r_real$insert_r\_$rprimer");
			}
		}
	}

	my $tm_text = get_tm(1);
	my $text = dna_magnify('',1);
	
	my $time = localtime;
	my $amplicon = "$amplicon_size bases ($fprimerpos - $rprimerpos)" unless $nb_page eq 'seq';
	$amplicon ||= "";
	
	my $output = <<EOT;
[$open_file{$nb_page} - Report generated at $time]

PRIMERS
-------

Forward: 5' $fprimer_mod 3' 
Reverse: 5' $rprimer_mod 3'
		
$amplicon

PRIMER DETAILS
--------------

$tm_text

SEQUENCE MAP
------------

$text
EOT
	my $file;
	
	unless ($open_file{$nb_page} eq 'File not saved') {
		my $filename = $open_file{$nb_page};
		$filename =~ s/\.ppr//g;
		$filename .= "_report";
		$file = $top->getSaveFile(-defaultextension=>'.txt', -initialfile=>$filename);
	} else {
		my $file = $top->getSaveFile(-defaultextension=>'.txt');
	}
	
	if (defined($file)) {
		open (REPORT, ">$file") || dialogue("Could not open file: $!");
		print REPORT $output;
		close (REPORT);
		sbarprint("\n$file report generated");
	}
}


#---------------#
# Gui dialogues #
#---------------#

sub dialogue {
	my ($message, $type) = @_;
	$type ||= 'OK';
	my $info_d = $top->messageBox(-title=>'Warning',
				-message=> $message,
				-type => $type,
				-icon=>'info'
				);
	return $info_d;
}


sub info {
	if (Exists($info_d)) {
		$info_d->deiconify();
		$info_d->raise();
		return;
	}	
	
	my $text = <<EOT;
PerlPrimer v$version
Copyright © 2003-2004 Owen Marshall\n
EOT
	my $text2 = <<EOT;
An application to design primers for PCR, Bisulphite PCR, Real-time PCR and Sequencing.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.\n
EOT
	my $address = <<EOT;	
http://perlprimer.sourceforge.net
EOT
	$info_d = $top->Toplevel(-title=>'About PerlPrimer ...');
	my $info_d_f = $info_d->Frame()->pack(-padx=>4, -pady=>4, -expand=>1, -fill=>'both');
	my $info_icon = $info_d_f->Label(-image => $pixmap)->pack(-side=>'left', -anchor=>'n');
	
	# Because of a few tweaks and the custom layout, this is the one window that
	# doesn't use pack_gui()
	
	# Want the text widget to have the same background colour as the window
	# (so the user can't tell it's a text widget!)  We're using a text widget
	# rather than a frame to get some formatting options (i.e. bold text)
	my $colour = $info_d->cget(-bg);	
	my $info_text = $info_d_f->ROText(
				-relief => 'flat',
				-wrap => 'word',
				-bg => $colour,
				-width => 50,
				-height => 12,
				-font => $gui_font,
			)->pack( -expand=>0, -fill=>'both', -anchor=>'nw');
	$info_text->tagConfigure('bold',
		-font => "$gui_font bold");
	$info_text->tagConfigure('center',
		-justify => 'center');
	$info_text->insert('end', $text, 'bold');
	$info_text->insert('end', $text2);
	$info_text->insert('end', $address);
	my $info_OK = $info_d->Button(
			-text => 'OK',
			-font => $gui_font,
			-padx=>4,
			-pady=>$button_pady,
			-command => sub {$info_d->destroy},
			)->pack(-pady=>4, -side=>'bottom');
	$info_d->Icon(-image => $pixmap);
}


sub canvas_info {
	if (Exists($canvas_info_d)) {
		$canvas_info_d->deiconify();
		$canvas_info_d->raise();
		return;
	}	
	
	$canvas_info_d = $top->Toplevel(-title=>'PerlPrimer Help');
	my $canvas_info_d_f = $canvas_info_d->Frame()->pack(-expand=>1, -fill=>'both');
	my $canvas_info_d_fb = $canvas_info_d->Frame()->pack(-side=>'bottom', -fill=>'x');
	pack_gui('', '', 'canvas_info_text', \$canvas_info_d_f, 't', 80, 35, 1, 3);
	pack_gui("OK", sub {$canvas_info_d->destroy}, "canvas_info_OK", \$canvas_info_d_fb, "b", "active", 0);
	$prf{canvas_info_text}->configure(-font=>"$gui_font");
	$prf{canvas_info_text}->tagConfigure('bold',
		-font => "$gui_font bold");
	$prf{canvas_info_text}->tagConfigure('center',
		-justify => 'center');
	
	my $info_text_dna=<<EOT;
\nThe grey line represents the DNA sequence, with genes [Standard PCR and Real-time PCR tabs] or CpG islands [Bisulphite PCR tab] marked in dark blue.  Intron/exon boundaries are displayed in the Real-time PCR tab.
EOT

	my $info_text_ranges=<<EOT;
\nThe light blue rectangle [Standard PCR and Bisulphite PCR tabs only] marks the minimum amplified area.  It can be automatically set to the gene boundaries by clicking the "Set from ORF" or "Set from CpG Islands" buttons, or can be resized using the left mouse button.  Alternatively, you can set it precisely in the "Amplified range" section.
		
The orange rectangle marks the maximum amplified area - no primers will be taken from outside of this range.  You can resize it by using the middle mouse button or Ctrl-left mouse button.
EOT

	my $info_text_primers=<<EOT;
\nClicking on a primer pair in the Results list displays the primers in red, drawn to scale.
EOT
	
	my $info_text_other=<<EOT;	
\nRight-clicking on the DNA will open a detailed view aligning the primers, DNA, translated ORF [Standard PCR and Real-time PCR tabs] or CpG island boundaries [Bisulphite PCR tab].  Intron/exon boundaries [Real-time PCR tab] or CpG residues [Bisulphite PCR tab] are also highlighted.  The "Copy printable" button in this dialogue window will copy the entire layout to the clipboard wrapped at 80 characters, which can then be pasted into a word processor/text editor and printed.
EOT
	
	$prf{canvas_info_text}->insert('1.0', " ", 'center');
	$prf{canvas_info_text}->imageCreate('end', -image=>$top->Pixmap(-data => $dna_canvas_pixmap));
	
	$prf{canvas_info_text}->insert('end', "\n\n\nDNA sequence\n", 'bold');
	$prf{canvas_info_text}->insert('end', $info_text_dna);
	$prf{canvas_info_text}->insert('end', "\n\nAmplified ranges\n", 'bold');
	$prf{canvas_info_text}->insert('end', $info_text_ranges);
	$prf{canvas_info_text}->insert('end', "\n\nPrimers\n", 'bold');
	$prf{canvas_info_text}->insert('end', $info_text_primers);
	$prf{canvas_info_text}->insert('end', "\n\nOther features\n", 'bold');
	$prf{canvas_info_text}->insert('end', $info_text_other);
}


sub acknowledgements {
	# just raise the dialogue if it already exists ....
	if (Exists($ack_d)) {
		$ack_d->deiconify();
		$ack_d->raise();
		return;
	}	
	
	$ack_d = $top->Toplevel(-title=>'Acknowledgements');
	my $ack_d_f = $ack_d->Frame()->pack(-expand=>1, -fill=>'both');
	my $ack_d_fb = $ack_d->Frame()->pack(-side=>'bottom', -fill=>'none');
	pack_gui('', '', 'ack_text', \$ack_d_f, 't', 70, 35, 1, 3);
	pack_gui("OK", sub {$ack_d->destroy}, "ack_OK", \$ack_d_fb, "b", "active", 0);
	$prf{ack_text}->configure(-font=>"$gui_font");
	$prf{ack_text}->tagConfigure('bold',
		-font =>"$gui_font bold");
	my $text_thermo=<<EOT;
Thermodynamic parameters are based on the following papers:

Allawi HT, SantaLucia J Jr.  Thermodynamics and NMR of internal G.T mismatches in DNA.  Biochemistry. 1997 Aug 26;36(34):10581-94

SantaLucia J Jr.  A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics.  Proc Natl Acad Sci U S A. 1998 Feb 17;95(4):1460-5. 

Allawi HT, SantaLucia J Jr.  Nearest neighbor thermodynamic parameters for internal G.A mismatches in DNA.  Biochemistry. 1998 Feb 24;37(8):2170-9.

Allawi HT, SantaLucia J Jr.  Thermodynamics of internal C.T mismatches in DNA.  Nucleic Acids Res. 1998 Jun 1;26(11):2694-701. 

Allawi HT, SantaLucia J Jr.  Nearest-neighbor thermodynamics of internal A.C mismatches in DNA: sequence dependence and pH effects.  Biochemistry. 1998 Jun 30;37(26):9435-44.

Peyret N, Seneviratne PA, Allawi HT, SantaLucia J Jr.  Nearest-neighbor thermodynamics and NMR of DNA sequences with internal A.A, C.C, G.G, and T.T mismatches.  Biochemistry. 1999 Mar 23;38(12):3468-77. 
EOT
	
	my $text_entropy=<<EOT;
Entropy corrections for PCR salt conditions are based on:

von Ahsen N, Wittwer CT, Schütz E.  Oligonucleotide Melting Temperatures under PCR Conditions: Nearest-Neighbor Corrections for Mg2+, Deoxynucleotide Triphosphate, and Dimethyl Sulfoxide Concentrations with Comparison to Alternative Empirical Formulas.  Clinical Chemistry.  2001; 47(11):1956-1961
EOT
		
	my $text_cpg=<<EOT;
Parameters for Bisulphate PCR primer design are based upon:
		
Warnecke PM, Stirzaker C, Song J, Grunau C, Melki JR, Clark SJ.  Identification and resolution of artifacts in bisulfite sequencing.  Methods. 2002 Jun;27(2):101-7.
EOT
			
	my $text_thanks=<<EOT;
Chris Vega
Katrina Bell
Richard Saffery
Nick Wong
Karl Billeter
EOT
	
	$prf{ack_text}->insert('0.1', "Thermodynamic parameters\n\n", 'bold');
	$prf{ack_text}->insert('end', $text_thermo);
	
	$prf{ack_text}->insert('end', "\n\nEntropy corrections\n\n", 'bold');
	$prf{ack_text}->insert('end', $text_entropy);
	
	$prf{ack_text}->insert('end', "\n\nBisulphate PCR primer design\n\n", 'bold');
	$prf{ack_text}->insert('end', $text_cpg);
	
	$prf{ack_text}->insert('end', "\n\nSpecial thanks for comments, testing and suggestions\n(in no particular order)\n\n", 'bold');
	$prf{ack_text}->insert('end', $text_thanks);
	
	# icon	
	$ack_d->Icon(-image => $pixmap);
}
	
		
sub view_spidey_out {
	my $spidey_d = $top->Toplevel(-title=>'Spidey output');
	my $spidey_d_f = $spidey_d->Frame()->pack(-expand=>1, -fill=>'both');
	my $spidey_d_fb = $spidey_d->Frame()->pack(-side=>'bottom', -fill=>'x');
	$spidey_d->Icon(-image => $pixmap);
	
	pack_gui('', '', 'spidey_text', \$spidey_d_f, 't', 60, 15, 1, '');
	pack_gui("OK", sub {$spidey_d->destroy}, "spidey_OK", \$spidey_d_fb, "b", "active", 0);
	pack_gui("Full output", sub {
			($spidey_out,) = run_spidey(0);
			$prf{'spidey_text'}->delete('0.1', 'end');
			$prf{'spidey_text'}->insert('0.1', $spidey_out)
		}, "spidey_full", \$spidey_d_fb, "b", "normal", 1);
	
	$prf{spidey_text}->insert('0.1', $spidey_out);
}	

#-------------#
# Preferences #
#-------------#

sub prefs {
	# Preferences dialogue
	
	# Preference file is $pref_file; all variables saved/restored are given in
	# the hash %pref_variables.  It should be fairly clear that this makes the
	# preference system easily extensible and yet also extremely simple ...
	
	# just raise the dialogue if it already exists ....
	if (Exists($prefs)) {
		$prefs->deiconify();
		$prefs->raise();
		return;
	}
	
	$prefs = $top->Toplevel(-title=>'Preferences');
	my $prefs_nb = $prefs->NoteBook(
			-inactivebackground=>"#$nb_colour",
			-relief => 'raised',
			-bd => 1,
			# -font => $gui_font,
		)->pack(
			-expand=>1,
			-padx=>4,
			-pady=>4,
			-fill=>'both',
		);
	
	my $prefs_page_general = $prefs_nb->add('files', -label=>'General', -anchor=>'nw');
	my $prefs_page_repeats = $prefs_nb->add('repeats', -label=>'Exclusions', -anchor=>'nw');
	my $prefs_page_bioinformatics = $prefs_nb->add('blast', -label=>'Bioinformatics', -anchor=>'nw');
	my $prefs_page_cloning = $prefs_nb->add('cloning', -label=>'Cloning', -anchor=>'nw');
	my $prefs_page_dimers = $prefs_nb->add('dimers', -label=>'Dimers', -anchor=>'nw');
	my $prefs_page_connection = $prefs_nb->add('connection', -label=>'Connection', -anchor=>'nw');
	my $prefs_page_gui = $prefs_nb->add('gui', -label=>'GUI', -anchor=>'nw');
	
	my $prefs_page_general_f = $prefs_page_general->Frame()->pack(-anchor=>'nw', -expand=>0, -fill=>'none');
	my $prefs_page_repeats_f = $prefs_page_repeats->Frame()->pack(-anchor=>'nw', -expand=>0, -fill=>'none');
	my $prefs_page_bioinformatics_f = $prefs_page_bioinformatics->Frame()->pack(-anchor=>'nw', -expand=>0, -fill=>'none');
	my $prefs_page_cloning_f = $prefs_page_cloning->Frame()->pack(-anchor=>'nw', -expand=>0, -fill=>'none');
	my $prefs_page_dimers_f = $prefs_page_dimers->Frame()->pack(-anchor=>'nw', -expand=>0, -fill=>'none');
	my $prefs_page_connection_f = $prefs_page_connection->Frame()->pack(-anchor=>'nw', -expand=>0, -fill=>'none');
	my $prefs_page_gui_f = $prefs_page_gui->Frame()->pack(-anchor=>'nw', -expand=>0, -fill=>'none');
	
	my $prefs_fb = $prefs->Frame()->pack(-side=>'bottom', -padx=>2, -fill=>'none');

	nr(\$prefs_page_general_f, 2);
		nr();
			$font = $gui_font_bold;
			pack_gui("Files/executables", "", "prefs_files_l", \$prf{prefs_gen}, 'lt');
			$font = $gui_font;
		nr();
			pack_gui("Home directory", $HOME, "prefs_files_home", \$prf{prefs_gen}, 20);
		nr();
			pack_gui("Path to Spidey executable", $file_pre{$os}, "prefs_files_spidey", \$prf{prefs_gen}, 20);
		
		nr('', 7);
		
		nr();	
			$font = $gui_font_bold;
			pack_gui("PCR component concentrations", "", "prefs_oligo", \$prf{prefs_gen}, 'lt');
			$font = $gui_font;
		nr();
			pack_gui("Mg++: ", $mg_conc, "prefs_oligo_runs", \$prf{prefs_gen}, 5, "", "", "mM");
		nr();
			pack_gui("Oligos: ", $oligo_conc, "prefs_oligo_repeats", \$prf{prefs_gen}, 5, '', "", "nM");
		nr();
			pack_gui("dNTPs: ", $dntp_conc, "prefs_oligo_runs", \$prf{prefs_gen}, 5, "", "", "mM");
		nr();
			pack_gui("Monovalent cations: ", $monovalent_cation_conc, "prefs_oligo_runs", \$prf{prefs_gen}, 5, '', "", "mM");
			
		nr('', 7);
		
		nr();	
			$font = $gui_font_bold;
			pack_gui("ORF / CpG island finding behaviour", "", "prefs_oligo", \$prf{prefs_gen}, 'lt');
			$font = $gui_font;
		nr('',0);
			pack_gui('Defer to capitalised regions', $defer_to_caps, "prefs_defer", \$prf{prefs_gen}, 'c', 0, 1);		
		nr('',2);

	nr(\$prefs_page_repeats_f, 2);
		nr();
			$font = $gui_font_bold;
			pack_gui("PCR excluded repeats / runs", "", "prefs_rrgen", \$prf{prefs_rr}, 'lt');
			$font = $gui_font;
		nr('',0);
			pack_gui('Exclude primers containing more than', $exclude_rr, "prefs_exclude_rr", '', 'c', 0, 1);
		nr();
			pack_gui('', $repeat, "prefs_rr_repeats", \$prf{prefs_rr}, 3, '','', 'Repeats or');
		# nr();
			pack_gui('', $run, "prefs_rr_runs", \$prf{prefs_rr}, 3, '', '', 'Runs');
	
		nr('', 7);
	
		nr();	
			$font = $gui_font_bold;
			pack_gui("Bisulphite PCR excluded repeats / runs", "", "prefs_rrbs_l", \$prf{prefs_rr}, 'lt');
			$font = $gui_font;
		nr('',0);
			pack_gui('Exclude primers containing more than', $exclude_rr_bs, "prefs_exclude_rr_bs", '', 'c', 0, 1);
		nr();
			pack_gui('', $repeat_bs, "prefs_rrbs_repeats", \$prf{prefs_rr}, 3, '', '', 'Repeats or');
		# nr();		
			pack_gui('', $run_bs, "prefs_rrbs_runs", \$prf{prefs_rr}, 3, '', '', 'Runs');
		nr('', 7);
		
		nr();
			$font = $gui_font_bold;
			pack_gui('Exclude %GC content', "", "prefs_gc", \$prf{prefs_rr}, 'lt');
			$font = $gui_font;
		nr();
			pack_gui("Only consider primers with ", $min_gc, "prefs_gc_max", \$prf{prefs_rr}, 3);
			pack_gui('-', $max_gc, "prefs_gc_min", \$prf{prefs_rr}, 3, '', '', '% GC content');
		nr();
			pack_gui('(Used when "Exclude %GC" is checked on a project page)','',"prefs_gc_note",\$prf{prefs_rr}, 'lt');
		
		nr('',2);
			
	nr(\$prefs_page_bioinformatics_f, 2);
		nr();
			$font = $gui_font_bold;
			pack_gui("BLAST search parameters", "", "prefs_bioinf_l", \$prf{prefs_bioinf}, 'lt');
			$font = $gui_font;
		nr();
			pack_gui('Expect value: ', $blast_expect, "prefs_blast_expect", \$prf{prefs_bioinf}, 4);
		nr();	
			pack_gui('Database: ', \$blast_database, 'prefs_blast_database', \$prf{prefs_bioinf}, 'be', 0, \@blast_database_array, 10);
		nr();	
			pack_gui('Limit to organism (Entrez query): ', \$blast_entrez_query, 'prefs_blast_entrez', \$prf{prefs_bioinf}, 'be', 0, \@blast_entrez_array, 20);
		
		nr('',7);
				
		nr();		
			$font = $gui_font_bold;
			pack_gui("CpG island prediction", "", "prefs_bioinf_l2", \$prf{prefs_bioinf}, 'lt');
			$font = $gui_font;
		nr();	
			pack_gui("Window size: ", $cpg_window, "prefs_bioinf_window", \$prf{prefs_bioinf}, 5,'','','bases');
		nr();	
			pack_gui("Minimum island size: ", $min_cpg_island, "prefs_bioinf_island", \$prf{prefs_bioinf}, 5, '','','bases');
		nr();	
			pack_gui("Minimum obs/exp: ", $cpg_oe, "prefs_bioinf_oe", \$prf{prefs_bioinf}, 5);
		nr();	
			pack_gui("Minimum GC content: ", $cpg_gc, "prefs_bioinf_gc", \$prf{prefs_bioinf}, 5, '','', '%');
			
		nr('',5);
		
		nr('',0);
			pack_gui('Emulate cpgplot', $cpgplot_method, "prefs_cpgplot_method", '', 'c', 0, 1);
		nr('',2);
		
	nr(\$prefs_page_cloning_f, 2);
		nr();
			$font = $gui_font_bold;
			pack_gui("Restriction enzyme cloning sequences", "", "prefs_bioinf_re_l", \$prf{prefs_bioinf}, 'lt');
			$font = $gui_font;
		nr('',0);
			pack_gui('Use 6-base cutting enzymes only', $simple_sites, "prefs_simple_sites", \$prf{prefs_bioinf}, 'c', 0, 1);
		nr('',0);
			pack_gui('Only list enzymes that do not cut sequence', $exclude_found_sites, "prefs_exclude_found_sites", \$prf{prefs_bioinf}, 'c', 0, 1);
		
		nr('',7);
	
	nr(\$prefs_page_dimers_f, 2);
		nr();
			$font = $gui_font_bold;
			pack_gui("Primer-dimer parameters", "", "prefs_dimer_l",'', 'lt');
			$font = $gui_font;
		nr();
			pack_gui('Calculate primer-dimer dG at ', $pd_temperature, 'prefs_dimer_temp', '', 3, '', '', '°C');
		
	nr(\$prefs_page_connection_f, 2);
		nr();
			$font = $gui_font_bold;
			pack_gui("Proxy server", "", "prefs_connection_font", \$prf{prefs_connection}, 'lt');
			$font = $gui_font;
		nr('',0);	
			pack_gui('Use http proxy server', $use_proxy, "prefs_connection_proxy", \$prf{prefs_connection}, 'c', 0, 1);
		nr();	
			pack_gui('Address: ', $http_proxy, 'prefs_connection_address', \$prf{prefs_connection}, 30);
			pack_gui('Port', $http_port, 'prefs_connection_port', \$prf{prefs_connection}, 5, 2);
		nr('',2);
	
			
	nr(\$prefs_page_gui_f, 3);
		nr();
			$font = $gui_font_bold;
			pack_gui("Fonts", "", "prefs_gui_font", \$prf{prefs_gui}, 'lt');
			$font = $gui_font;
		nr();	
			my $font_families = [sort $top->fontFamilies];
			pack_gui('Font:', \$gui_font_face, 'prefs_gui_family', \$prf{prefs_gui}, 'be', 0, $font_families);
			pack_gui(' Size:', \$gui_font_size, 'prefs_gui_size', \$prf{prefs_gui}, 'be', 1, [(3 .. 32)], 3);
		nr();
			pack_gui('Use OS font defaults', $font_override, 'prefs_gui_override', \$prf{prefs_gui}, 'c', 0, 1);
		nr('',7);
		
		nr();	
			$font = $gui_font_bold;
			pack_gui("Previously opened files", "", "prefs_gui_mru", \$prf{prefs_gui}, 'lt');
			$font = $gui_font;
		nr();
			pack_gui('List the last ', $mru_number, "prefs_gui_mru_number", \$prf{prefs_gui}, 3, '','', 'files');
		nr('',7);
		
		nr();	
			$font = $gui_font_bold;
			pack_gui("Mouse Wheel", "", "prefs_gui_mouse", \$prf{prefs_gui}, 'lt');
			$font = $gui_font;
		nr();
			pack_gui('Mouse wheel scrolls ', $scroll_factor, "prefs_gui_wheel", \$prf{prefs_gui}, 3, '','', 'lines');
		nr('',2);
		
					
	pack_gui("OK", sub {
			# if the salt concentration has been changed, we need to recalculate %oligo_dG
			recalculate_dG();
			
			# write the data to the pref_file ...
			my $file_data = "";
			foreach my $i (keys %pref_variables) {
				my $pointer = $pref_variables{$i};
				$file_data .= "$i = $$pointer\n";
			}
			foreach my $i (keys %pref_arrays) {
				my $pointer = $pref_arrays{$i};
				$file_data .= "$i = [".join(",",@$pointer)."]\n";
			}
			open (PREFS, ">$pref_file") || dialogue("Could not open file: $!");
			print PREFS $file_data;
			close (PREFS);
			
			# You really don't want to know why we have to read the prefs file
			# just after writing it ...
			
			# (hint: if we don't, and we open and OK the prefs window twice,
			# perl starts saying that 1.5-0.2 = 1 ...  I don't know why ... )
			
			read_prefs();
			$prefs->destroy;
		}, "prefs_OK", \$prefs_fb, "b", "active", 0);
	
	pack_gui("Revert", sub {
			# This is perhaps a bit inefficient, re-reading the pref_file each time
			# the dialogue is cancelled.  However, it does not seem to cause a
			# noticable time-lag and it's certainly the simplest way to do things
			# (otherwise we'd have to save each value to a temporary hash and update
			# each variable on prefs_OK ... to make this extensible we'd have to make
			# a hash with the same keys as %pref_variables ...
			# it's definitely possible, but a bit of a fuss!!
			read_prefs();
			$prefs->destroy;
		}, "prefs_cancel", \$prefs_fb, "b", "normal", 1);

	$prefs->Icon(-image => $pixmap);
}


sub read_prefs {
	return unless -e $pref_file;
	open (PREFS, "<$pref_file") || dialogue("Could not open prefs file for reading: $!");
	my $pointer;
	while (<PREFS>) {
		s/[\n\r]//g;
		next if /^$/;
		my ($key, $value) = split / = /;
		if ($value =~ /\[(.+)\]/) {
			# an array
			$pointer = $pref_arrays{$key};
			@$pointer = split(",",$1);
		} else {
			$pointer = $pref_variables{$key};
			$$pointer = $value;
		}
	}
	close (PREFS);
}


#------------------------#
# DNA graphical routines #
#------------------------#

sub find_gene {
	# Find genes (or CpG islands, or whatever else you might fancy) based on
	# case (upper-case marks genes) (can have more than one ORF in the same
	# sequence, which the ORF-finding sub doesn't allow; also speeds things
	# up a bit by not having to ORF/CpG-find each time)
	$_ = $_[0];
	return if length($_) == 0;
	my ($subroutine) = get_variables(qw(find_sub));

	s/\>.*\n//g; #remove FASTA formatting if it exists
	s/[\s|\n]//g; #remove spaces/tabs/new lines
	s/[^ACGTXN]//ig; #remove any other characters that shouldn't be there
	my @array = ();
	my $prev = 0;
	my $nb_page = which_nb_page();
	
	# if the user wishes to defer to capitalised regions, do so ...
	# i.e. a sequence with lowercase and capitalised regions will be assumed to
	# use caps to mark the gene/cpg island
	# NB - this is no longer the default behaviour
	if ($defer_to_caps) {
		my $i;	
		while (/([a-z]*)([A-Z]*)(?=[a-z]*)/xg) {
			$i++;
			push @array, [$prev+length($1), $prev+length($1)+length($2)] if (length($2)>0 && length($1)>0);
			$prev += length($1)+length($2);
		}
		
		return @array if @array;
	}
	
	# if we didn't find any marked regions, let's hand it over to our built in,
	# trusty orf and cpg finding routines ...
	return  &$subroutine($_);
	
	# my ($sel5, $sel3);
	# ($sel5, $sel3) = &$subroutine($_);
	# return ([$sel5, $sel3]);
}

sub find_orf {
	# no gene marked - try to find ORF since user has requested it
	$_ = lc($_[0]);
	s/[\s|\n]//g;
	s/[^ACGTXN]//ig;
	my ($seq_ref) = get_variables('seq');

	my @orf=();
	my $seq = $_;
	for my $i (0 .. 2) {
		$orf[$i][0] = 0;
		my $met_init = 0;
		my $orf_count = 0;
		my $orf_start = -1;
		my $j;
		for ($j = $i; $j<(length($seq)-3); $j+=3) {
			my $codon = uc(substr($seq, $j, 3));
			my $aa = $genetic_code{$codon};
			unless ($aa eq "*") {
				# only take ORF from initiating Met codon
				next if ($aa ne "M") && ($met_init==0);
				$met_init=1;
				$orf_start = $j if $orf_start==-1;
				$orf_count++;
				next;
			}
			# stop codon
			if ($orf_count > $orf[$i][0]) {
				$orf[$i][0]=$orf_count+1;
				$orf[$i][1]=$orf_start;
				# end orf
				$orf[$i][2]=$j+3;
				# frame
				$orf[$i][3]=$i;
			}
			$orf_count=$met_init=0;
			$orf_start=-1;
		}
		if ($orf_count > $orf[$i][0]) {
			$orf[$i][0]=$orf_count+1;
			$orf[$i][1]=$orf_start;
			# end orf
			$orf[$i][2]=$j+3;
			# frame
			$orf[$i][3]=$i;
		}
		$orf_start=$orf_count=$met_init=0;
	}
		
	# sort by orf length:
	@orf = sort {@$b[0] <=> @$a[0]} @orf;
	
	return unless $orf[0] > 1;
		
	# Recreate DNA sequence based on ORF
	my $dna3 ="";
	my $dna1 = substr($seq, 0, $orf[0][1]);
	my $dna2 = uc(substr($seq, $orf[0][1], $orf[0][0]*3));
	$dna3 = substr($seq, $orf[0][0]*3+$orf[0][1], ) unless $orf[0][0]*3+$orf[0][1]>=length($seq);
	
	push my @orf_return, [$orf[0][1], $orf[0][2]];
	
	# Re-enter DNA sequence using capitalised ORF, lowercase upstream and downstream regions
	$$seq_ref->delete(0.1,"end");
	$$seq_ref->insert(0.1,$dna1.$dna2.$dna3);
	
	return @orf_return;
	# return $orf[0][1], $orf[0][2];
}

sub find_cpg {
	# no CpG marked - try to find CpG islands from sequence
	# two methods: correct and cpgplot emulation
	$_ = lc($_[0]);
	s/[\s|\n]//g;
	s/[^ACGT]//ig;

	my $seq = $_;
	my $seq_len = length($seq);
	my @cpg_island =();
		
	my $cpg_flag = 0;
	my $cpg_start = 0;
	my $cpg_count = 0;
	my $island_max = 0;
	my $last_cpg_pos = 0;
	my $win_count =0;
	
	# when cpgplot method is 1, we're emulating cpgplot behaviour (see comment
	# below)
	my $cpg_av = ($cpgplot_method ? 10 : 9);
	
	# CpG finding routine
	#
	# This is based on the method first described in Gardiner-Garden and Frommer
	# (1987) and used in many other programs since.  The most common utility is
	# cpgplot, but cpgplot does not give the same results as this routine! 
	# Analysing the source, it's apparent that cpgplot is buggy and in fact does
	# not use the parameters it  claims to, but slightly an underestimate. 
	#
	# The worst bug of all in cpgplot is that the avg readings for each window are
	# inflated: they are not an average at all; rather, 11 frames are read and the
	# totals divided by 10!! The second bug comes in calculating the number of CpG
	# residues in a window - cpgplot actually adds in an extra base when
	# performing this calculation, but does not add in an extra base when counting
	# the c's and g's in a sequence, or the %GC content!!
	# 
	# PerlPrimer has the option to deliberately allow for this and return the
	# same results as cpgplot, simply because cpgplot has been around for so long 
	# seems to be the de facto standard for CpG island prediction.  However, by
	# default PerlPrimer uses my algorithmically correct approach.
	
	my (@cpg_oe, @cpg_pgc)=();
	
	# We use an average of 10 (or 11 in cpgplot's case) window values at each
	# point, so to save recalculating the same value 10 times we save it in an
	# array ...
	for my $i (0 .. $seq_len-$cpg_window) {
		$_ = substr($seq, $i, $cpg_window);
		
		# calculate O/E
		my $g = tr/Gg/Gg/;
		my $c = tr/Cc/Cc/;
		my $exp_gc = ($c*$g)/(length($_)-1);
		if ($cpgplot_method) {
			$_ = substr($seq, $i, $cpg_window+1);
		}
		my $cg = s/cg/cg/ig;
		
		# save data
		$cpg_oe[$i] = $cg/$exp_gc;
		$cpg_pgc[$i] = gc($_);
	}
	
	# Calculate the averages and find the islands
	for my $i (0 .. $seq_len-$cpg_window-10) {
		# the real position is in the middle of the window...
		my $cpg_real_pos = int($i+($cpg_window/2));
		
		# of course, this is (10 or 11)/2 bases upstream from the real mid-point,
		# since we're averaging 10 frames from here (another bug in cpgplot!)
		$cpg_real_pos += int($cpg_av/2) unless $cpgplot_method;
		
		my $sum_cg = 0;
		my $sum_obs_exp = 0;
		
		# Cpgplot uses an 11 window "average" (it thinks it's only using 10 ...)
		for my $j ($i .. $i+$cpg_av) {
			$sum_cg += $cpg_pgc[$j];
			$sum_obs_exp += $cpg_oe[$j];
		}
		
		my $av_obs_exp = $sum_obs_exp/10;
		my $av_cg = $sum_cg/10;
		
		if (($av_cg > $cpg_gc) && ($av_obs_exp>$cpg_oe)) {
			# A possible island: set the start point if it's the first
			$cpg_start=$cpg_real_pos if $cpg_flag==0;
			$cpg_flag=1;
			next;
		}
		unless ($cpg_flag==0) {
			if ($i-$cpg_start < $min_cpg_island) {
				# don't save the island if it's less than the minimum size
				$cpg_flag=0;
				$cpg_count=0;
				$cpg_start=0;
				next;
			}
			
			# Save the island
			# This should have an end point of $cpg_real_pos-1, as we're
			# now a base ahead of the last positive, but cpgplot has it
			# as $cpg_real_pos+1 (?!?)
			my $cpg_island_end = ($cpgplot_method ? $cpg_real_pos-1 : $cpg_real_pos+1);
			push @cpg_island, [$cpg_start, $cpg_island_end];
			
			# reset flags
			$cpg_flag=0;
			$cpg_count=0;
			$cpg_start=0;
		}
	}
	
	$prf{"bisul_seq"}->delete(0.1,"end");
	my $cpg_pos = 0;

	for my $i (0 .. $#cpg_island) {
		print "Island No. $i: start $cpg_island[$i][0], length $cpg_island[$i][1]\n";
		
		# Recreate DNA sequence based on ORF
		my $dna1 = substr($seq, $cpg_pos, $cpg_island[$i][0]-$cpg_pos);
		# print "pos, len = $cpg_pos, $cpg_island[$i][0]\n";
		my $dna2 = uc(substr($seq, $cpg_island[$i][0], $cpg_island[$i][1]-$cpg_island[$i][0]));
		# print "dna2 pos, len = $cpg_island[$i][0], $cpg_island[$i][1]\n";
		$cpg_pos = $cpg_island[$i][1];
		# print "pos = $cpg_pos\n";
	 	
		# Re-enter DNA sequence using capitalised ORF, lowercase 5' and 3' regions
		$prf{"bisul_seq"}->insert("1.end",$dna1.$dna2);
	}
	
	# Last part of the reconstructed sequence:
	my $dna3 = substr($seq, $cpg_pos);
	$prf{"bisul_seq"}->insert("1.end",$dna3);
	
	return (-1) unless @cpg_island;
	return @cpg_island;
	# return the start of the first island and the end of the last island
	return $cpg_island[0][0], $cpg_island[$#cpg_island][0]+$cpg_island[$#cpg_island][1];
}

sub draw_dna {
	my $old_defer = $defer_to_caps; # save state
	$defer_to_caps = 1 if @_ && $_[0] == 1;
	
	my ($min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p, $canvas) 
		= get_variables(qw(min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p canvas)); 

	# my $canvas = which_canvas_ref();
	my $seq = get_seq();
	
	# delete old canvas elements
	$$canvas->delete('dna', 'amp_range', 'gene_range', 'dna_range','primerf','primerr', 'direction', 'ie_boundary');
	
	return unless $seq =~ /[a|t|c|g]/i;
	
	# get dimensions of canvas
	my $width = $$canvas->width;
	my $height= $$canvas->height;
	$dna_canvas_size=($width-($dna_canvas_offset*2))/length($seq);
	
	my $canv = which_nb_page();
	
	# Amplicon size
	if (($$max_ampsize)&&($$min_ampsize)&&($$max_range)) {
		# print "max_range_5p = $$max_range_5p, max_range_3p = $$max_range_3p\n";
		$min_amp_canvas{$canv} = $$max_range_5p * $dna_canvas_size + $dna_canvas_offset if defined($$max_range_5p);
		$max_amp_canvas{$canv} = $$max_range_3p * $dna_canvas_size + $dna_canvas_offset if defined($$max_range_3p);
		$$canvas->createRectangle($min_amp_canvas{$canv},$dc_sel_y1-$dc_sel_offset2,$max_amp_canvas{$canv},$dc_sel_y2+$dc_sel_offset2, -fill=>'orange', -outline=>'orange4', -width=>1, -tag=>'amp_range');
	}
	
	# range
	if (defined($$max_range)&&defined($$min_range)) {
		$min_range_canvas{$canv} = $$min_range * $dna_canvas_size + $dna_canvas_offset;
		$max_range_canvas{$canv} = $$max_range * $dna_canvas_size + $dna_canvas_offset;
		$$canvas->createRectangle($min_range_canvas{$canv},$dc_sel_y1,$max_range_canvas{$canv},$dc_sel_y2, -fill=>'dodgerblue', -outline=>'dodgerblue4', -tag=>'dna_range');
	}
	
	# max_range_5p and 3p variables set
	check_range();
	
	# draw DNA
	$$canvas->createRectangle($dna_canvas_offset,$dc_dna_y1,$width-($dna_canvas_offset+1),$dc_dna_y2, -fill=>'grey50', -tag=>'dna');
	$$canvas->createText(5,$dna_canvas_middle, -text=>"5'", -justify=>'right', -tag=>'direction');
	$$canvas->createText($width-6,$dna_canvas_middle, -text=>"3'", -justify=>'right', -tag=>'direction');
	
	# gene (from capitals, if present)
	my @gene_array = find_gene($seq);

	for my $i (0 ... $#gene_array) {
		my($gene_5p,$gene_3p)=($gene_array[$i][0], $gene_array[$i][1]);
		if (defined($gene_5p)&&defined($gene_3p)) {
			$min_range_canvas{$canv} = $gene_5p * $dna_canvas_size + $dna_canvas_offset;
			$max_range_canvas{$canv} = $gene_3p * $dna_canvas_size + $dna_canvas_offset;
			$$canvas->createRectangle($min_range_canvas{$canv},$dc_dna_y1,$max_range_canvas{$canv},$dc_dna_y2, -fill=>'midnightblue', -tag=>'gene_range');
		}
	}
	
	# QPCR specific:
	if ($qpcr_flag==1) {
		foreach my $i (@intron_exon_bounds) {
			my $iex = $i * $dna_canvas_size + $dna_canvas_offset;
			$$canvas->createLine($iex, $dc_dna_y1, $iex, $dc_dna_y2, -fill=>'white', -tag=>'ie_boundary');
		}
	}
	
	$defer_to_caps = $old_defer; # restore state
}

sub items_drag {
	my ($min_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p) 
		= get_variables(qw(min_ampsize max_range_5p min_range max_range max_range_3p)); 

	my $canv = which_nb_page();
	return if $canv eq 'qpcr';
	my $ref = $_[1];
	my $x = $_[0];
	
	my $width = $$ref->width;
	my $dna_max = $width-$dna_canvas_offset;
	my $dna_min = $dna_canvas_offset;
	
	if (($$max_range)&&($$min_range)) {
		$min_range_canvas{$canv} = $$min_range * $dna_canvas_size + $dna_canvas_offset;
		$max_range_canvas{$canv} = $$max_range * $dna_canvas_size + $dna_canvas_offset;
	}

	$x = sprintf("%.0f", $$ref->canvasx($x));
	$x = $dna_min if $x < $dna_min;
	$x = $dna_max if $x > $dna_max;

	if (($min_amp_canvas{$canv}) && ($max_amp_canvas{$canv})) {
		$x = $min_amp_canvas{$canv} if $x < $min_amp_canvas{$canv};
		$x = $max_amp_canvas{$canv} if $x > $max_amp_canvas{$canv};
	}

	if ($x-$min_range_canvas{$canv}<$max_range_canvas{$canv}-$x) {
		$$ref->coords('dna_range', $x, $dc_sel_y1, $max_range_canvas{$canv},$dc_sel_y2);
		$$min_range = sprintf("%.0f", ($x-$dna_canvas_offset) / $dna_canvas_size);
		$$min_range = $$max_range_5p if $$max_range_5p && $min_range<$$max_range_5p;

	} else {
		$$ref->coords('dna_range', $min_range_canvas{$canv}, $dc_sel_y1, $x, $dc_sel_y2);
		$$max_range = sprintf("%.0f", ($x-$dna_canvas_offset) / $dna_canvas_size);
		$$max_range = $$max_range_3p if $$max_range_3p && $$max_range>$$max_range_3p;
	}
	
	$$min_ampsize = $$max_range - $$min_range;
}

sub amplicon_drag {
	my ($min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p) 
		= get_variables(qw(min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p)); 

	my $canv = which_nb_page();
	return if $canv eq 'qpcr';
	my $ref = $_[1];
	my $x = $_[0];

	my $seq = get_seq();
	$seq =~ s/[\n|\s]//g;
	my $seq_len = length($seq)-1;
	
	my $width = $$ref->width;
	my $dna_max = $width-$dna_canvas_offset;
	my $dna_min = $dna_canvas_offset;
	my $dna_max2 = $$max_range*$dna_canvas_size + $dna_canvas_offset;
	my $dna_min2 = $$min_range*$dna_canvas_size + $dna_canvas_offset;

	# grab default values if not assigned
	if (($min_ampsize)&&($max_ampsize)) {
		$min_amp_canvas{$canv} = $$min_range * $dna_canvas_size + $dna_canvas_offset unless $min_amp_canvas{$canv};
		$max_amp_canvas{$canv} = ($$max_ampsize+$$min_range) * $dna_canvas_size + $dna_canvas_offset unless $max_amp_canvas{$canv} ;
	}

	$x = sprintf("%.0f", $$ref->canvasx($x));
	$x = $dna_min if $x < $dna_min;
	$x = $dna_max if $x > $dna_max;

	
	if ($x-$min_amp_canvas{$canv}<$max_amp_canvas{$canv}-$x)  {
		# Changing amp_size at 5' end
		$x = $dna_min2 if $x > $dna_min2;
		$$ref->coords('amp_range', $x, $dc_sel_y1-$dc_sel_offset2, $max_amp_canvas{$canv},$dc_sel_y2+$dc_sel_offset2);
		$min_amp_canvas{$canv}=$x;
		$$max_range_5p = sprintf("%.0f", ($x-$dna_canvas_offset) / $dna_canvas_size);
	} else {
		# Changing amp_size at 3' end
		$x = $dna_max2 if $x < $dna_max2;
		$$ref->coords('amp_range', $min_amp_canvas{$canv}, $dc_sel_y1-$dc_sel_offset2, $x, $dc_sel_y2+$dc_sel_offset2);
		$max_amp_canvas{$canv}=$x;
		$$max_range_3p = sprintf("%.0f", ($x-$dna_canvas_offset) / $dna_canvas_size);
		$$max_range_3p = $seq_len if $$max_range_3p>$seq_len;

	}
	$$max_ampsize = sprintf("%.0f", ($max_amp_canvas{$canv}-$min_amp_canvas{$canv}) / $dna_canvas_size);
}

sub get_gene {
	my ($min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p) 
		= get_variables(qw(min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p)); 
	my $seq = get_seq();
	my @gene_array = find_gene($seq);
	
	return unless @gene_array;		
	($$min_range, $$max_range)=($gene_array[0][0], $gene_array[$#gene_array][1]-1);
	
					
	$$max_range_5p ||= $$min_range;
	$$max_range_3p ||= $$max_range;
	$$max_range_5p = $$min_range if $$min_range < $$max_range_5p;
	$$max_range_3p = $$max_range if $$max_range > $$max_range_3p;
	
	$$min_ampsize=$$max_range - $$min_range;
	$$max_ampsize=$$max_range - $$min_range; # if ($max_range - $min_range)>$max_ampsize;
		
	draw_dna(1);
}

sub get_cpg {
	my ($min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p) 
		= get_variables(qw(min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p)); 
	my $seq = get_seq();
	my @gene_array = find_gene($seq);
	if ($gene_array[0] == -1) {
		sbarprint("\nNo CpG Island found!");
		return;
	}
	
	($$min_range, $$max_range)=($gene_array[0][0], $gene_array[$#gene_array][1]-1);
	
	$$max_range_5p ||= $$min_range;
	$$max_range_3p ||= $$max_range;
	$$max_range_5p = $$min_range if $$min_range < $$max_range_5p;
	$$max_range_3p = $$max_range if $$max_range > $$max_range_3p;
	
	$$min_ampsize=$$max_range - $$min_range;
	$$max_ampsize=$$max_range - $$min_range; # if ($max_range - $min_range)>$max_ampsize;
	
	draw_dna(1);
}

sub reset_bounds {
	my $page = which_nb_page();
	return if $page eq "qpcr";
	my ($min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p, $sub_ref) 
		= get_variables(qw(min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p subroutine)); 
	$$max_range_5p = $$max_range_3p = $$min_range = $$max_range = $$min_ampsize = $$max_ampsize = undef;
	&$sub_ref();
}

sub step_in {
	my ($min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p, $sub_ref, $canv_ref, $seq_ref) 
		= get_variables(qw(min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p subroutine canvas seq)); 
	
	&$sub_ref() unless defined($$max_range_5p) && defined($$max_range_3p);
	$$min_range += 10;
	$$max_range -= 10;
	if ($$min_range >= $$max_range) {
		sbarprint("\nNo primers found");
		$cancel=1;
		return;
	}

	$$min_ampsize=$$max_range - $$min_range;
	$$max_ampsize=$$max_range_3p - $$max_range_5p;

	draw_dna(1);
}

sub step_out {
	my ($min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p, $sub_ref, $canv_ref, $seq_ref) 
		= get_variables(qw(min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p subroutine canvas seq));
	my $seq_len = length(get_seq());
	&$sub_ref() unless defined($$max_range_5p) && defined($$max_range_3p);
	
	$cancel=1 if ($$max_range_5p==0) && ($$max_range_3p==$seq_len);
	
	$$max_range_5p -= 10;
	$$max_range_3p += 10;
	$$max_range_5p = 0 if $$max_range_5p < 0;
	$$max_range_3p = $seq_len if $$max_range_3p>$seq_len;
		
	$$min_ampsize=$$max_range - $$min_range;
	$$max_ampsize=$$max_range_3p - $$max_range_5p;

	draw_dna(1);
}

sub step_cloning {
	my ($min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p, $sub_ref, $canv_ref, $seq_ref) 
		= get_variables(qw(min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p subroutine canvas seq));
	my $seq_len = length(get_seq());
	&$sub_ref() unless defined($$max_range_5p) && defined($$max_range_3p);
	
	$cancel=1 if ($$max_range_5p==0) && ($$max_range_3p==$seq_len);
	
	$$max_range_5p -= 10;
	$$max_range -= 10;
	$$max_range_5p = 0 if $$max_range_5p < 0;
	$$max_range = $$min_range if $$max_range < $$min_range;
	
	return if $$max_range_5p == 0 && $$max_range == $$min_range;
		
	$$min_ampsize=$$max_range - $$min_range;
	$$max_ampsize=$$max_range_3p - $$max_range_5p;

	draw_dna(1);
}

sub add_cloning {
	my ($seq, $fprimer, $fprimerpos, $fprimerposl, $rprimer, $rprimerpos, $rprimerposl) = @_;
	
	# Adding an extra sequence is all very well.
	# But we want to automagically keep it in frame when creating
	# fusion constructs ...
	
	# Warning!! This is clumsy!
	
	my @gene_array=find_gene($seq);
	my $gene_frame=$gene_array[0][0]%3;
	
	my $insert_f = "";
	my $insert_r = "";
	my ($primer_seq_5f_real, $rprimer_seq_5r_real);
	
	if ($primer_seq_5f && $fprimer) {
		# calculate frame position for forward primer + insert
		$_ = $primer_seq_5f;
		
		my ($fspl1, $fspl2)=split('_');
		my $renz_offset_f = 3-length($fspl2);
		
		s/[\|_]//g;
		$primer_seq_5f_real = $_;
		
		if (defined($primer_seq_5f_frame) && ($primer_seq_5f_frame ne "")) {
			my $fjoinframe = 3-$primer_seq_5f_frame;
			my $primer_frame_f= (($fprimerpos - $gene_frame)%3 + $fjoinframe + $renz_offset_f)%3;
	
			# Add spacers if needed between primer sequences and the inserts,
			# and adjust the spacing of the primerpos and primerposl variables
			# so that the primers are still displayed in the correct positions			
			$insert_f = "A" x $primer_frame_f;
		}
		
		$fprimer = "$primer_seq_5f_real$insert_f$fprimer";
		$fprimerpos -= length($insert_f)+length($primer_seq_5f_real);
		$fprimerposl += length($insert_f)+length($primer_seq_5f);
	}
	
	if ($primer_seq_5r && $rprimer) {
		# calculate frame position for reverse primer + insert
		$_ = $primer_seq_5r;
		
		my ($rspl1, $rspl2)=split('_');
		# my $renz_offset_r = 3-length($rspl2);
		
		s/[\|_]//g;
		my $primer_seq_5r_real = $_;
		$rprimer_seq_5r_real = reverse($primer_seq_5r_real);
		
		if (defined($primer_seq_5r_frame) && ($primer_seq_5r_frame ne "")) {
			my $rjoinframe = $primer_seq_5r_frame;
			my $primer_frame_r= (2-($rprimerpos - $gene_frame)%3 + $rjoinframe)%3;			
			$insert_r = "A" x $primer_frame_r;
		}
		
		$rprimer = "$rprimer$insert_r$rprimer_seq_5r_real";
		$rprimerpos += length($insert_r)+length($primer_seq_5r_real);
		$rprimerposl += length($insert_r)+length($primer_seq_5r_real);
	}
	
	return ($primer_seq_5f_real, $insert_f, $fprimer, $fprimerpos, $fprimerposl, $rprimer_seq_5r_real, $insert_r, $rprimer, $rprimerpos, $rprimerposl);
}

sub dna_magnify {
	# This is one long, convoluted subroutine that attempts
	# to please all of the people all of the time :)
	# A bit ungainly, but it just kept on growing ...
	# Design could still do with a bit of tweaking ...
	my ($x, $report) = @_;
	my $seq = get_seq();
	my $page = which_nb_page();
	my ($max_range_5p, $min_range, $max_range, $max_range_3p, $primer_array, $list_ref, $ref) 
		= get_variables(qw(max_range_5p min_range max_range max_range_3p primers hlist canvas));
	my ($fprimer, $rprimer, $fprimerpos, $rprimerpos, $fprimerposl, $rprimerposl, $primer_seq_5f_real, $insert_f, $primer_seq_5r_real, $insert_r);
	
	# get dimensions of canvas
	my $width = $$ref->width;
	$dna_canvas_size=($width-($dna_canvas_offset*2))/length($seq);
	
	# This sub used to jump to the exact position right clicked; however, all I ever
	# want to do is look at the primers, so I feel the routine would be better jumping
	# straight to the forward primer ...
	# my $base_pos = sprintf("%.0f",($x-$dna_canvas_offset)/$dna_canvas_size);
	
	my $text_ref = \$prf{view_base_text};
	
	# Local subroutines:	
	my $copy_report = sub {
		# Copies the content into 80-column wrapped text
		my ($clip, $last_i);
		my $length = 80;
		for (my $i=$length; $i<=length($seq)+$length; $i+=$length) {
			$last_i ||= 0;
			for my $j (0 .. 4) {
				$clip .= $$text_ref->get("$j.$last_i","$j.$i");
				$clip .= "\n";
			}
			$clip .= "\n\n";
			$last_i=$i;
		}
		return $clip;
	};
	
	my $copy_printable = sub {
		# copy report for the clipboard ...
		$top->clipboardClear;
		my $clip = &$copy_report();
		$top->clipboardAppend($clip);
	};
	
	my $view_fprimer = sub {
		$$text_ref->see("1.".($fprimerpos+10));
	};
		
	my $view_rprimer = sub {
		$$text_ref->see("1.".($rprimerpos-10));
	};
	
	###  In progress!!
	# my $multiline_pos = sub {
		# my ($x,$y) = @_;
		# # if ($multiline) {
			# # # single lines to multiple lines
			# # 
		# # }
	# };	
			
	
	# if previous window exists we destroy it first, to avoid confusion
	# (trying to simply erase the contents in the text box and 
	# re-configuring the buttons has about the same effect with Perl/Tk,
	# so there's not much point in doing so ...)
	if (Exists($view_base)) {
		$view_base->destroy;
	}	
	
	$view_base = $top->Toplevel(-title=>'Details ...');
	my $view_basef = $view_base->Frame()->pack(-expand=>1, -fill=>'both');
	my $view_basefb = $view_base->Frame()->pack(-side=>'left');
	my $view_basefbr = $view_base->Frame()->pack(-side=>'right');
	
	$view_base->withdraw if $report;
		
	pack_gui('', '', 'view_base_text', \$view_basef, 't', 80, 4, 1, 2);
		
	pack_gui('OK', sub {$view_base->destroy}, 'view_base_ok', \$view_basefb, 'b', 'active', 0);
	pack_gui('Copy printable', \&$copy_printable, 'view_base_copy', \$view_basefb, 'b', '', 1);
	pack_gui('Forward', \&$view_fprimer, 'view_base_pf', \$view_basefbr, 'b', '', 0, '', '', 'disabled' );
	pack_gui('Reverse', \&$view_rprimer, 'view_base_pr', \$view_basefbr, 'b', '', 1, '', '', 'disabled' );
	
	$$text_ref->tagConfigure('red',
		-foreground => 'red');
	$$text_ref->tagConfigure('blue',
		-foreground => 'midnightblue',
		-background => '#cee3ee');
	$$text_ref->tagConfigure('blue_cpg',
		-foreground => 'royalblue',
		-background => 'royalblue');
	$$text_ref->tagConfigure('numbers',
		-foreground => 'grey75',
		-background => 'grey35');
	$$text_ref->tagConfigure('grey60',
		-foreground => 'grey50');
	$$text_ref->tagConfigure('grey40',
		-foreground => 'grey40');
	$$text_ref->tagConfigure('black',
		-foreground => 'black',
		-background => 'grey89');
	$$text_ref->tagConfigure('codon',
		-foreground => 'black',
		-background => 'grey81');
	$$text_ref->tagConfigure('orange',
		-foreground => 'orange');
	$$text_ref->tagConfigure('dodgerblue',
		-foreground => 'dodgerblue');
	$$text_ref->tagConfigure('ie',
		-foreground => 'white',
		-background => 'red');
			
	for my $i (1 .. 3) {
		$$text_ref->insert("$i.0","\n")
	}
	
	# my $list_ref = which_hlist();
	my @sel= $$list_ref->selectionGet;
	my $sel=shift(@sel);
	
	# Selected primers, if any ...
	if (defined($sel)) {
		$fprimer = $$primer_array[$sel][0];
		$rprimer = reverse($$primer_array[$sel][4]) if $$primer_array[$sel][5];
		
		$fprimerpos = $$primer_array[$sel][1];
		$rprimerpos = $$primer_array[$sel][8];
		
		$fprimerposl = $$primer_array[$sel][2];
		$rprimerposl = $$primer_array[$sel][6];
		
		# need to do this before adding sequences ...
		my $spacerrp = " "x($rprimerpos-$rprimerposl-$fprimerpos-$fprimerposl-5) if $rprimer;
		
		# Added sequences?
		if ($page eq 'pd' && (($primer_seq_5f) || ($primer_seq_5r))) {
			($primer_seq_5f_real, $insert_f, $fprimer, $fprimerpos, $fprimerposl, $primer_seq_5r_real, $insert_r, $rprimer, $rprimerpos, $rprimerposl) = add_cloning($seq, $fprimer, $fprimerpos, $fprimerposl, $rprimer, $rprimerpos, $rprimerposl);
		}
		
		my $spacerfp = " "x($fprimerpos-3);
		# print "($rprimerpos-$rprimerposl-$fprimerpos-$fprimerposl-5) = ".($rprimerpos-$rprimerposl-$fprimerpos-$fprimerposl-5);
		$$text_ref->insert('2.0',$spacerfp."5' ", 'grey60');
		$$text_ref->insert('2.end',$fprimer, 'red');
		
		if ($rprimer) {
			$$text_ref->insert('2.end'," 3'".$spacerrp."3' ", 'grey60');		
			$$text_ref->insert('2.end',$rprimer, 'red');
			$$text_ref->insert('2.end'," 5'", 'grey60');
		} else {
			$$text_ref->insert('2.end'," 3'", 'grey60');
		}
		
		# Make buttons active or inactive:
		$prb{view_base_pf}->configure(-state=>($fprimer ? 'active' : 'disabled'));
		$prb{view_base_pr}->configure(-state=>($rprimer ? 'active' : 'disabled'));
	}
	
	# Numbering along the top (starts at 0)
	my $numbers="0.........";
	for my $i(1 .. length($seq)/10) {
		$numbers.=($i*10)."."x(9-length($i));
	}
	
	# insert the numbers and the sequence
	$$text_ref->insert('1.0',$numbers,'numbers');
	$$text_ref->insert('3.0',$seq, 'grey40');
	
	# translation of ORFs
	# (or if Bisulphite sequencing highlights CpG island)
	my @gene_array=find_gene($seq);
	my $prev_end = 0;
	for my $i (0 .. $#gene_array) {
		my $peptide="";
		my ($start, $end)=($gene_array[$i][0],$gene_array[$i][1]);
		$$text_ref->tagAdd('black', "3.$start", "3.$end");
		if ($page eq 'bis') {
			$peptide = "-"x($end-$start);
		} else {
			my ($j,$k)=(0,0);
			for ($j=$start; $j<$end; $j+=3) {
				my $codon = uc(substr($seq, $j, 3));
				$peptide .= $genetic_code{$codon}."  ";
				if ($k == 1) {
					my $tend = $j+3;
					$$text_ref->tagAdd('codon', "3.$j", "3.$tend");
				}
				$k = 1 - $k;
			}
		}
		my $spacer=" "x($start-$prev_end);
		$prev_end=$end;
		$$text_ref->insert("4.0",$spacer);
		
		unless ($page eq 'bis') {
			$$text_ref->insert("4.end",$peptide,'blue');
		} else {
			$$text_ref->insert("4.end",$peptide,'blue_cpg');
		}
	}
	
	# highlight intron-exon boundaries for QPCR
	# or ranges otherwise
	if ($page eq 'qpcr') {
		foreach my $i (@intron_exon_bounds) {
			my $tag_start=$i-1;
			my $tag_end=$i+1;
			$$text_ref->tagAdd('ie', "3.$tag_start", "3.$tag_end");
		}
	} else {
		my $sel_range=$$max_range+1;
		my $sel_range_3p=$$max_range_3p+1 if $$max_range_3p;
		$$text_ref->tagAdd('orange', "3.$$max_range_5p", "3.$sel_range_3p") if $$max_range_3p;
		$$text_ref->tagAdd('dodgerblue', "3.$$min_range", "3.$sel_range");
	}
	
	if ($page eq 'bis') {
		# CpG labeling
		for my $i (0 .. length($seq)) {
			if (substr($seq, $i, 2) =~ /CG|cg/ ) {
				my $tagcpg = $i+2;
				$$text_ref->tagAdd('ie', "3.$i", "3.$tagcpg");
			}
		}
		
	}
	
	# this is very inefficient and messy, but gets the job done ...
	if ($report) {
		my $sequence_report = &$copy_report;
		$view_base->destroy;
		return $sequence_report;
	}
	
	# icon	
	$view_base->Icon(-image => $pixmap);
	
	# Show the position where the user clicked
	# my $center_base_pos = $base_pos-30;
	# $$text_ref->see("1.$center_base_pos");
	
	# Show the forward primer
	my $center_base_pos = $fprimerpos-30;
	$$text_ref->see("1.$center_base_pos");
}


#---------------#
# File commands #
#---------------#

sub pp_file_open {
	# NB - there's several reasons why we're not using XML as a file format here. 
	# The first is simplicity, both for the programmer and for the end-user (who
	# would have to download and install a CPAN module otherwise - which might be
	# daunting to many).  The second is that I can't abide XML.
	
	# The system that follows is both simple and highly extensible (adding
	# support for another varible simple entails a single line in either the
	# %variables or %arrays hashes - everything else is automatic ... )
	
	my ($file) = @_;
	$file ||= $top->getOpenFile(-filetypes=>$file_types);
	return unless defined($file);
		
	# clear previous selections:
	#$max_range=$min_range=$max_range_5p=$max_range_3p=undef;
	
	unless (open (SEQ, "<$file")) {
		dialogue("Could not open file $file: $!");
		return;
	}
	sbarprint("\nOpening $file ...");
	my @file_data = <SEQ>;
	close SEQ;
	
	my ($nb_page, $name) = open_file_type($file, @file_data);
	
	# redraw dna
	draw_dna();
				
	# reset qpcr_flag
	$qpcr_flag = 0;
	sbarprint("\n$file opened successfully");
	$file =~ s/\//\\/g if $os eq 'win'; # path bug in Win32 Perl/Tk
	
	# since we're taking the name from FASTA files directly ...
	my $title_filename = ($name ? $name : $file);
	
	# my $full_path = $file;
	# $file =~ s/.*[\/\\]//g;
	$top->configure(-title=>"PerlPrimer v$version - $title_filename");
	$open_file{$nb_page} = $title_filename;
	
	recently_used_files($file);
}

sub open_file_type {
	my $file = shift;
	my @file_data = @_;
	my ($type, $page, $name);
	foreach (@file_data) {
		# This loop is in case the first line is blank.
		# It's probably unnecessary ..
		if (/nb = .+/) {
			# file appears to be a perlprimer file
			$type = 'ppr';
			return ($page, $name) = open_ppr($file, @file_data);
		} elsif (/^\>.+/) {
			# file appears to be FASTA format
			$type = 'fasta';
			return ($page, $name) = open_fasta($file, @file_data);
		} elsif (/[\w|\d]+/) {
			# file is unknown format
			last;
		}
	}
	dialogue("File does not appear to be a PerlPrimer file or in FASTA format.\nIf you are trying to open a DNA sequence file,\nplease use the open icon next to the sequence entry field") unless ($type);
}

sub open_ppr {
	my $file = shift;
	my @file_data = @_;
	
	my $nb_page;
	my $array_flag=0;
	@save_selection=();
	$save_seq=$save_seq2="";
	my %flag=();

	my $pointer;
	my ($key, $value, $lines, $open_percent);
	my $total_lines = @file_data + 1;
	foreach (@file_data) {
		# we're supporting both *nix and win32 ...
		s/[\n\r]//g;
	
		$open_percent = sprintf("%.f", $lines++/$total_lines*100);
		sbarprint("\nOpening $file ... $open_percent\%");

		# get nb_page entry (NB: must come before all data entry)
		# or look for [array] section
		if (/nb = (.*)/) {
			$nb->raise($1);
			$top->update;
			$nb_page = $1;
			next;
		} elsif (/\[arrays\]/) {
			$array_flag=1;
			next;
		}
		
		# ignore blank lines
		next if /^$/;
		
		# grab values
		# my ($key, $value) = split / = /;
		if (/ = /) {
			($key, $value) = split / = /;
		} elsif ($array_flag==0) {
			# allow multi-line variables - eg for genomic DNA
			$$pointer .= "\n$_";
			next;
		}
		
		# set variables
		unless ($array_flag==1) {
			$pointer = $variables{$nb_page}{$key};
			$$pointer = $value;
			next;
		}
		
		# array parsing
		$pointer = $arrays{$nb_page}{$key};
		unless ($flag{$key}) {
			# clear the array
			@$pointer=();
			$flag{$key}=1;
		}
		# this is so easy!
		my @temp = split(/ /,$value);
		if ($#temp == 0) {
			push (@$pointer, $temp[0]);
		} else {
			push (@$pointer, [@temp]);
		}
	}
	
	# restore sequence(s)
	my ($seq_ref, $hlist) = get_variables(qw(seq hlist));
	$$seq_ref->delete(0.1,"end");
	$$seq_ref->insert(0.1,$save_seq);
	if ($save_seq2) {
		# only for qpcr
		$prf{qdna_seq}->delete(0.1,"end");
		$prf{qdna_seq}->insert(0.1,$save_seq2);
	}
	
	# restore results
	sort_primers();
		
	# restore intron-exon boundaries if qpcr
	$qpcr_flag = 1 if $nb_page eq 'qpcr';
		
	# restore selection
	for my $i (0 .. $#save_selection) {
		$$hlist->selectionSet("$save_selection[$i]");
	}
	if (@save_selection) {
		browse_primer($save_selection[$#save_selection]);
		$$hlist->see($save_selection[$#save_selection]);
	}
	
	return $nb_page;
}

sub open_fasta {
	my $file = shift;
	my @file_data = @_;
	my ($flag, $dna, $name, $range_5a, $range_5b, $range_3a, $range_3b, $page, $nb_page);
	my ($key, $value, $lines, $open_percent);
	my $total_lines = @file_data + 1;
	foreach (@file_data) {
		$open_percent = sprintf("%.f", $lines++/$total_lines*100);
		sbarprint("\nOpening $file ... $open_percent\%");

		next if /^$/;
		if ($flag) {
			$dna .= "$_";
			next;
		} elsif (/^\>/) {
			($name, $range_5a, $range_5b, $range_3a, $range_3b, $page) = /^\>\s*(.+?)\s*(?:5prime_region\[(\d+)-(\d+)\])?\s*(?:3prime_region\[(\d+)-(\d+)\])?\s*(?:page\[(\d+)\])?\s*$/;
			$page = 1 unless ($page && $nb_page_ref{$page});
			$nb_page = $nb_page_ref{$page};
			$nb->raise($nb_page);
			$top->update;
			
			my ($min_ampsize, $max_ampsize, $max_range_5p, $min_range, $max_range, $max_range_3p) 
				= get_variables(qw(min_ampsize max_ampsize max_range_5p min_range max_range max_range_3p));
			
			if ($range_5a && $range_5b && $range_3a && $range_3b) {
				$$max_range_5p = $range_5a;
				$$min_range = $range_5b;
				$$max_range = $range_3a;
				$$max_range_3p = $range_3b;
				
				$$min_ampsize=$$max_range - $$min_range;
				$$max_ampsize=$$max_range_3p - $$max_range_5p;
			}
			
			$flag = 1;
		}
	}
	
	# insert FASTA sequence
	my ($seq_ref) = get_variables(qw(seq));
	$$seq_ref->delete(0.1,"end");
	$$seq_ref->insert(0.1,$dna);
	
	return ($nb_page, $name);
}

sub pp_file_save {
	my $nb_page = which_nb_page();
	if ($nb_page eq "primer") {
		dialogue("You can't save a file from the primer information page - please save from the project page instead");
		return;
	}
	
	# Either prompt for file name if file not saved or "Save as ..." was
	# called, or save without prompting if file has been saved before
	# (It's amazing how complicated "simple" behaviour can be ...)
	my ($flag) = @_ || 0;
	my $file;
	if ($open_file{$nb_page} eq 'File not saved') {
		$file = $top->getSaveFile(-filetypes=>$file_types, -defaultextension=>'.ppr');
	} elsif (($flag == 1) || ($open_file{$nb_page} !~ /\.ppr/)) {
		$file = $top->getSaveFile(-filetypes=>$file_types, -defaultextension=>'.ppr', -initialfile=>$open_file{$nb_page});
	} else {
		$file = $open_file{$nb_page};
	}
	
	return unless defined($file);
	$file .= '.ppr' unless $file =~ /\.ppr/;
	
	sbarprint("\nSaving $file ...");
	
	my $total_keys = keys(%{ $variables{$nb_page} }) + keys(%{ $arrays{$nb_page} }) + 1;
	my ($keys_saved, $saved_percent);
	
	# Now the guts of the saving routine - see the pp_file_open routine for general
	# comments ... this is all straightforward stuff
		
	# save sequence data
	$save_seq = get_seq();
	$save_seq2 = $prf{qdna_seq}->get(0.1,"end") if $nb_page eq "qpcr";
	
	# save selection data
	my ($hlist) = get_variables('hlist');;
	@save_selection = $$hlist->selectionGet;
	
	# save variables
	my $file_data = "nb = $nb_page\n";
	my $pointer;
	foreach my $i (keys %{ $variables{$nb_page} }) {
		$pointer = $variables{$nb_page}{$i};
		$file_data .= "$i = $$pointer\n" if defined($$pointer);
		$saved_percent = sprintf("%.f", $keys_saved++/$total_keys*100);
		sbarprint("\nSaving $file ... $saved_percent\%");
	}
	
	# save arrays
	$file_data .= "[arrays]\n";
	foreach my $i (keys %{ $arrays{$nb_page} }) {
		$pointer = $arrays{$nb_page}{$i};
		for my $j (0 .. $#$pointer) {
			unless (@{$$pointer[$j]}) {
				$file_data .= "$i = $$pointer[$j]\n";
			} else {
				$file_data .= "$i = @{$$pointer[$j]}\n";
			}
		}
		$saved_percent = sprintf("%.f", $keys_saved++/$total_keys*100);
		sbarprint("\nSaving $file ... $saved_percent\%");
	}
	
	# write file
	open (SEQ, ">$file") || dialogue("Could not open file for saving: $!\n");
	print SEQ $file_data;
	close (SEQ);
	
	sbarprint("\n$file saved successfully");
	
	# Set the program title to reflect the new file name
	# $file =~ s/.*\/.*\///g;
	$file =~ s/\//\\/g if $os eq 'win'; # path bug in Win32 Perl/Tk
	$top->configure(-title=>"PerlPrimer v$version - $file");
	$open_file{$nb_page} = $file;
	
	recently_used_files($file);
}

sub recently_used_files {
	# recently used files
	my ($file) = @_;
	if ($file) {
		for my $i (0 .. $#mru) {
			$mru[$i] = "" if $mru[$i] eq $file;
		}
		push(@mru, $file);
	}
	
	# the subroutine gets called when the program loads, so we need an escape		
	return if @mru == 0;
	
	# clean up the space from before if it exists
	my @new_mru;
	for my $i (0 .. $#mru) {
		push(@new_mru, $mru[$i]) if ($mru[$i]);
	}
	
	# take only the last requested files so the list doesn't grow too big
	my $start = $#new_mru-($mru_number-1)<0 ? 0 : $#new_mru-($mru_number-1);
	my $end = $#new_mru;
	@mru = @new_mru[$start .. $end];
	
	# clear the old menu
	$menu_mru->delete(1,'end');
	
	# ... and insert the new list
	for my $i (0 .. $#mru) {
		$menu_mru->insert(0,'command', -label=>"$mru[$i]", -command =>[sub {pp_file_open($mru[$i])}] );
	}
}

sub open_seq {
	# opens file and writes contents to text widget
	my $ref = $_[0];
	my $file = $top->getOpenFile();
	if (defined($file)) {
		open (SEQ, "<$file") || dialogue("Could not open file: $!");
		$$ref->delete(0.1,'end');
		my $seq;
		while (<SEQ>) {
			s/^>.*//g; # remove fasta
			# s/[\n\r]//g; # Remove lines (is this necessary??)
			s/[\s\t]//g; # remove white space
			# s/^.*[^AGCTNX].*$//ig; # not sure about this!
			$seq .= "$_\n";
		}
		
		close (SEQ);
		$$ref->insert('end',$seq);
		sbarprint("\n$file opened successfully");
	}
}

sub save_seq {
	# saves sequence in text widget to file
	my $ref = $_[0];
	my $file = $top->getSaveFile(-defaultextension=>'.txt');
	my $seq = $$ref->get(0.1, 'end');
	if (defined($file)) {
		open (SEQ, ">$file") || dialogue("Could not open file: $!");
		print SEQ $seq;
		close (SEQ);
		sbarprint("\n$file saved successfully");
	}
}


#------------------#
# Text subroutines #
#------------------#

sub clear_text {
	my $ref = $_[0];
	$$ref->delete(0.1,'end');
}

sub select_all_text {
	my $ref = $_[0];
	$$ref->selectAll;
}

sub lc_text {
	my $ref = $_[0];
	my $text = $$ref->get(0.1,"end");
	my $text_lc = lc($text);
	$$ref->delete(0.1,'end');
	$$ref->insert(0.1,$text_lc);
}

sub reverse_complement {
	my $ref = $_[0];
	my $text = $$ref->get(0.1,"end");
	$text =~ s/[\n\s]//g;
	my $text_comp_rev = reverse(complement($text));
	$$ref->delete(0.1,'end');
	$$ref->insert(0.1,$text_comp_rev);
}


#----------------#
# Reference subs #
#----------------#

sub get_seq {
	my ($seq_ref) = get_variables('seq');
	return 1 unless $$seq_ref;
	$_ = $$seq_ref->get(0.1,"end");
	
	s/\>.*\n//g; #remove FASTA formatting if it exists
	s/[\n\r]//g; # remove line breaks
	return $_;
}

sub which_nb_page {
	return $nb->raised;
}

sub get_variables {
	my $note_page = $nb->raised;
	$null = undef;
	my @return_args;
	foreach my $var (@_) {
		if ($page_specific_vars{$note_page}{$var}) {
			push @return_args, $page_specific_vars{$note_page}{$var};
		} else {
			push @return_args, \undef;
		}
	}
	return @return_args;
}

sub check_packages {
	my $failed;
	foreach (@_) {
		if ($failed_packages =~ / $_ /) {
			dialogue("PerlPrimer requires the package $_ for this feature to work\n\nPlease download and install $_ from CPAN (http://cpan.org)");
			$failed = 1;
		}
	}
	return 1 if $failed;
}

#----------------------#
# Balloon state toggle #
#----------------------#

sub balloon_toggle {
	$Balloon->configure(-state => 'none') if $balloon_help == 0;
	$Balloon->configure(-state => 'balloon') if $balloon_help == 1;
}


sub load_icon_data {

#----------------------#
# Icon and pixmap data #
#----------------------#

# ... just a little bit of bloat :)

$perlprimer_icon = <<'end_of_pixmap';
/* XPM */
static char * dna_icon2_xpm[] = {
"32 32 25 1",
" 	c #FFFFFF",
".	c #402E33",
"+	c #A0302D",
"@	c #918FC4",
"#	c #B89590",
"$	c #CAC7EE",
"%	c #E4B82A",
"&	c #969096",
"*	c #E3382B",
"=	c #924E32",
"-	c #E76029",
";	c #EFD4C7",
">	c #69647F",
",	c #F9FBF8",
"'	c #DD6C69",
")	c #80696A",
"!	c #C0817E",
"~	c #9F6160",
"{	c #E08B27",
"]	c #F7B6B3",
"^	c #D2D1CF",
"/	c #B3AEB2",
"(	c #F92C2B",
"_	c #AEACF1",
":	c #E7A19D",
"     *(;,,;:>@_>##!((;=((+      ",
"      -;,,/>@_@/^,;(-;)=.       ",
"      +*~&#/__&,,,-(({          ",
"      .!,::#.>*',;*(**          ",
"       ),;~&,'--,,~=++          ",
"        &/&,!(*%:^)             ",
"         +#,'*({(#&             ",
"        .(:,^+=**/^.            ",
"       +(%%],,#&&&)++           ",
"       +({*(!&>>....~           ",
"      +#+{(()$@@@>)!^^          ",
"     +%;^/))/$____::,,!(+       ",
"     ({(*;/&&&@)>)!#,,/*+       ",
"     +-((~^^,__@:],@$^,;=       ",
"      .=~&,,,___:]___,#+*+      ",
"      =*#/&&&&>.!]$_$^!(***     ",
"      +**+#__#:@@/^//,,]'%*     ",
"       ++./__::__,,,^+'#/{      ",
"          &/@#!@_$,^'*(-(+      ",
"          !,##>&&^/,]({{*+      ",
"          &,:@_@~'/^,,{{+       ",
"           /~>@@**-(~,:+        ",
"             ).&(%-(*,)         ",
"           )#^^,;%{*~/&)        ",
"           (-*+,,'(~))&/.       ",
"           {((*;,;+))!))+       ",
"           %-(',,,>)!&^^'=      ",
"       =*=!](*!,/>>)~#,,/*+     ",
"       **(;,''~.>>>&:];':(*     ",
"       *(*,,,^>>@@_$/^'**++     ",
"      +=*{,,/)~@__>^,,,-*       ",
"     *(''*!)#~.>@$]];,;((+      "};
end_of_pixmap
		
$icon_open = <<'end_of_pixmap';
/* XPM */
static char * open_20_xpm[] = {
"20 20 27 1",
" 	c None",
".	c #020501",
"+	c #37372E",
"@	c #5B5D5A",
"#	c #8E8F7E",
"$	c #BABCAC",
"%	c #B8BB98",
"&	c #ACB0A3",
"*	c #1A1B18",
"=	c #43433B",
"-	c #CDD0AD",
";	c #737561",
">	c #A1A68A",
",	c #898F73",
"'	c #D7D9C6",
")	c #525548",
"!	c #AEB293",
"~	c #C4C7A4",
"{	c #242521",
"]	c #61624E",
"^	c #808673",
"/	c #C2C5B4",
"(	c #959B7E",
"_	c #D5D8B6",
":	c #B1B39D",
"<	c #999B8E",
"[	c #6C6D59",
"                    ",
"                    ",
"  .....             ",
" +#<<<^=            ",
" @____~[{           ",
"=$_~~~-(+           ",
"=/%!!%%>;======+    ",
"+->(>>>>>>:::::#.   ",
"{/,)+++++++++++{..  ",
"{:]<$$$///$////$$$@ ",
"{$=&__---_-----~~!) ",
"{<<'_---------~%!,+ ",
"{<<'-~--~~~~%~%:>]* ",
"*&/~%%%%%%%!!>>>^+  ",
"*$/!!>>>>>>>((,,].  ",
"{&,^;;;;[[[]]]))+   ",
" ...............    ",
"                    ",
"                    ",
"                    "};
end_of_pixmap

$icon_save = <<'end_of_pixmap';
/* XPM */
static char * save_20_xpm[] = {
"20 20 27 1",
" 	c None",
".	c #1B2226",
"+	c #728BA0",
"@	c #ADC9E0",
"#	c #424B54",
"$	c #CA766B",
"%	c #5F7688",
"&	c #E3E5E1",
"*	c #8DAAC4",
"=	c #2A3138",
"-	c #AEAEAC",
";	c #878A89",
">	c #656764",
",	c #DF9A8E",
"'	c #4E3F3B",
")	c #55534E",
"!	c #DCDDDA",
"~	c #7D8690",
"{	c #3B424A",
"]	c #F8FAF7",
"^	c #C7C6C1",
"/	c #6E757D",
"(	c #4A5E6D",
"_	c #7F98AD",
":	c #757879",
"<	c #9BBAD0",
"[	c #F5E8E8",
"                    ",
"                    ",
" .{='''''''''''={.  ",
" #@+,,,,,,,,,,,/*#  ",
" {<%,$$$$$$$$$$/_=  ",
" #<~![[[[[[[[[[/+=  ",
" {<~!]]]]][&[[[:_=  ",
" {<~^!&!!&!&&!!/+=  ",
" {</[]]]]]]]]]]~+=  ",
" {<~^!!!!&!&!!!/+=  ",
" {</!]]]]]]]]]]:+=  ",
" {<_~//////////~+.  ",
" {*__++++++_+++_+.  ",
" {*_%:;----;;((+_.  ",
" {*+#^^;:^^^-:(+_.  ",
" #*+#[^))^^^^;(%_.  ",
" {*+#&-)'-^^[/%%_.  ",
" .%%#^->:-!&&:(%~.  ",
"  .==))>>>::>=={(.  ",
"                    "};
end_of_pixmap
		
$icon_clear = <<'end_of_pixmap';
/* XPM */
static char * clear_20_xpm[] = {
"20 20 27 1",
" 	c None",
".	c #000100",
"+	c #504B3D",
"@	c #989981",
"#	c #332B22",
"$	c #6D7163",
"%	c #ADB8A6",
"&	c #1E1E1B",
"*	c #CEBF94",
"=	c #8D745E",
"-	c #554537",
";	c #A6AA95",
">	c #847B60",
",	c #3D3E32",
"'	c #5A5C4F",
")	c #BDA793",
"!	c #B5AB88",
"~	c #E3D3A4",
"{	c #CFBAAA",
"]	c #48382B",
"^	c #282017",
"/	c #6A634C",
"(	c #797E6C",
"_	c #6D5644",
":	c #A38E79",
"<	c #C9BE9F",
"[	c #11100E",
"                    ",
"                    ",
"  ^-&               ",
"  -)_.              ",
"  ]))^              ",
"  ^={_.             ",
"   #=:^             ",
"    #=_^  ..        ",
"     ]_- .['&       ",
"      ]_[(@%,   .   ",
"       -;%%$+#     .",
"      &$%%$>*' . .  ",
"     &(%%$>~<!&..   ",
"     [%%'>~~>@>..  .",
"     [$'>~<~!'(/..  ",
"      .'*~)@;@^',.. ",
"       #>!!/]'+.    ",
"        [,+,&..     ",
"                    ",
"                    "};
end_of_pixmap

$icon_info = <<'end_of_pixmap';	
/* XPM */
static char * info_20_xpm[] = {
"20 20 147 2",
"  	c None",
". 	c #ABABB1",
"+ 	c #D6D6DE",
"@ 	c #E5E5EF",
"# 	c #E5E5F2",
"$ 	c #D8D8E6",
"% 	c #BBBBC9",
"& 	c #6F7078",
"* 	c #CECED4",
"= 	c #F4F4F9",
"- 	c #F8F8FB",
"; 	c #F0F0F7",
"> 	c #E9E9F4",
", 	c #E4E4F2",
"' 	c #DEDFF0",
") 	c #D7D8EA",
"! 	c #8B8B97",
"~ 	c #C4C4CA",
"{ 	c #F6F6FB",
"] 	c #FCFCFD",
"^ 	c #F9F9FC",
"/ 	c #F1F1F8",
"( 	c #ECECF6",
"_ 	c #E6E6F3",
": 	c #E0E1F1",
"< 	c #DADDEF",
"[ 	c #D6D9EB",
"} 	c #6C6D77",
"| 	c #EEEEF3",
"1 	c #FBFBFD",
"2 	c #F3F3F9",
"3 	c #EEEEF7",
"4 	c #E8E8F5",
"5 	c #E2E2F2",
"6 	c #D7DBED",
"7 	c #C5C8DA",
"8 	c #BFBFC4",
"9 	c #FAFAFC",
"0 	c #F6F6FA",
"a 	c #F2F3F9",
"b 	c #E7E8F5",
"c 	c #DCE0EF",
"d 	c #D6DCEC",
"e 	c #D5DBED",
"f 	c #55575E",
"g 	c #CBCBD2",
"h 	c #F2F2F9",
"i 	c #F7F7FB",
"j 	c #F9F9FB",
"k 	c #F7F8FB",
"l 	c #F4F4FA",
"m 	c #E0E4F1",
"n 	c #D6DFEC",
"o 	c #D4DCED",
"p 	c #6A6C77",
"q 	c #C6C6CD",
"r 	c #EFEFF7",
"s 	c #F6F6F8",
"t 	c #F9F9FA",
"u 	c #F6F7FA",
"v 	c #ECECF2",
"w 	c #E1E6F1",
"x 	c #D6E3EB",
"y 	c #D4DEEC",
"z 	c #676973",
"A 	c #A8A8AF",
"B 	c #EAEAF4",
"C 	c #EDEEF5",
"D 	c #F1F1F7",
"E 	c #E7E7E9",
"F 	c #F2F2F8",
"G 	c #EBECF3",
"H 	c #E3E6ED",
"I 	c #DCE5EC",
"J 	c #D5E8EA",
"K 	c #D3E0EC",
"L 	c #404248",
"M 	c #DADBE6",
"N 	c #E3EAEF",
"O 	c #E5EBF1",
"P 	c #E0E3E6",
"Q 	c #E5EAF2",
"R 	c #E3E8F2",
"S 	c #D8DDE6",
"T 	c #CED6DC",
"U 	c #D6E9E9",
"V 	c #D3E9E9",
"W 	c #B2BBC8",
"X 	c #8D8F96",
"Y 	c #DEE3ED",
"Z 	c #DCE7EB",
"` 	c #C6D0D0",
" .	c #D8E2E9",
"..	c #D7E0E9",
"+.	c #CED5D4",
"@.	c #D4E4E6",
"#.	c #D3EBE9",
"$.	c #C8D6DF",
"%.	c #43454B",
"&.	c #8D8E98",
"*.	c #D5D8E6",
"=.	c #CBCFD2",
"-.	c #C5CBC5",
";.	c #CFD5CE",
">.	c #C1C9C2",
",.	c #D1DDDE",
"'.	c #C4D0DA",
").	c #55585F",
"!.	c #847F75",
"~.	c #C9C3B4",
"{.	c #C7C5C0",
"].	c #CDC9BE",
"^.	c #D0CAB5",
"/.	c #C1B99B",
"(.	c #2E2E2C",
"_.	c #F0E4BA",
":.	c #EBE0B7",
"<.	c #E1D29A",
"[.	c #D5C58C",
"}.	c #AE9E66",
"|.	c #13110A",
"1.	c #F1E5BC",
"2.	c #F1E7C5",
"3.	c #E6D69E",
"4.	c #D8C890",
"5.	c #AFA06C",
"6.	c #15130C",
"7.	c #EFE4BF",
"8.	c #E7DEBE",
"9.	c #DACD9C",
"0.	c #CDBE8A",
"a.	c #B0A16D",
"b.	c #14120B",
"c.	c #E7DDB6",
"d.	c #E2D8B4",
"e.	c #D1C28C",
"f.	c #C3B47D",
"g.	c #9A8D5F",
"h.	c #0F0D08",
"i.	c #898060",
"j.	c #BCB18A",
"k.	c #AD9F6C",
"l.	c #82764C",
"m.	c #312C1A",
"n.	c #2D2D2D",
"o.	c #404040",
"p.	c #060606",
"                                        ",
"            . + @ # $ % &               ",
"          * = - ; > , ' ) !             ",
"        ~ { ] ^ / ( _ : < [ }           ",
"        | ^ 1 - 2 3 4 5 < 6 7           ",
"      8 = ^ ] 9 0 a 3 b c d e f         ",
"      g h i j j ^ k l 3 m n o p         ",
"      q r 2 s t 9 9 u v w x y z         ",
"      A B C D E F ; G H I J K L         ",
"        M N O P Q R S T U V W           ",
"        X Y Z `  ...+.@.#.$.%.          ",
"          &.*.=.-.;.>.,.'.).            ",
"            !.~.{.].^./.(.              ",
"              _.:.<.[.}.|.              ",
"              1.2.3.4.5.6.              ",
"              7.8.9.0.a.b.              ",
"              c.d.e.f.g.h.              ",
"              i.j.k.l.m.                ",
"                n.o.p.                  ",
"                                        "};
end_of_pixmap

$icon_magnify = <<'end_of_pixmap';	
/* XPM */
static char * stock_zoom_in_20_xpm[] = {
"20 20 108 2",
"  	c None",
". 	c #3E3E3E",
"+ 	c #535353",
"@ 	c #585858",
"# 	c #565656",
"$ 	c #525252",
"% 	c #3F3F3F",
"& 	c #484848",
"* 	c #797979",
"= 	c #A7A7A7",
"- 	c #BDBDBD",
"; 	c #C1C1C1",
"> 	c #B5B5B5",
", 	c #8C8C8C",
"' 	c #2B2B2B",
") 	c #4D4D4D",
"! 	c #9A9A9A",
"~ 	c #CBCBCB",
"{ 	c #D9D9D9",
"] 	c #E2E2E2",
"^ 	c #E6E6E6",
"/ 	c #E4E4E4",
"( 	c #D5D5D5",
"_ 	c #6E6E6E",
": 	c #272727",
"< 	c #D1D1D1",
"[ 	c #F3F3F3",
"} 	c #F8F8F8",
"| 	c #F9F9F9",
"1 	c #F7F7F7",
"2 	c #F1F1F1",
"3 	c #C6C6C6",
"4 	c #626262",
"5 	c #1A1A1A",
"6 	c #757575",
"7 	c #CFCFCF",
"8 	c #FAFAFA",
"9 	c #FCFCFC",
"0 	c #A1A1A1",
"a 	c #454545",
"b 	c #DDDDDD",
"c 	c #F0F0F0",
"d 	c #E8E8E8",
"e 	c #C3C3C3",
"f 	c #343434",
"g 	c #434343",
"h 	c #AAAAAA",
"i 	c #E0E0E0",
"j 	c #FDFDFD",
"k 	c #A2A2A2",
"l 	c #F2F2F2",
"m 	c #EAEAEA",
"n 	c #DFDFDF",
"o 	c #808080",
"p 	c #080808",
"q 	c #4C4C4C",
"r 	c #C9C9C9",
"s 	c #F4F4F4",
"t 	c #9F9F9F",
"u 	c #737373",
"v 	c #909090",
"w 	c #9D9D9D",
"x 	c #B7B7B7",
"y 	c #060606",
"z 	c #D4D4D4",
"A 	c #E9E9E9",
"B 	c #444444",
"C 	c #424242",
"D 	c #7A7A7A",
"E 	c #DBDBDB",
"F 	c #C0C0C0",
"G 	c #050505",
"H 	c #EDEDED",
"I 	c #D7D7D7",
"J 	c #8E8E8E",
"K 	c #D2D2D2",
"L 	c #A3A3A3",
"M 	c #030303",
"N 	c #929292",
"O 	c #979797",
"P 	c #414141",
"Q 	c #D6D6D6",
"R 	c #C8C8C8",
"S 	c #BCBCBC",
"T 	c #D8D8D8",
"U 	c #DEDEDE",
"V 	c #CDCDCD",
"W 	c #949494",
"X 	c #171717",
"Y 	c #181818",
"Z 	c #6A6A6A",
"` 	c #CECECE",
" .	c #D3D3D3",
"..	c #3B3B3B",
"+.	c #ACACAC",
"@.	c #C2C2C2",
"#.	c #898989",
"$.	c #393939",
"%.	c #000000",
"&.	c #0A0A0A",
"*.	c #262626",
"=.	c #656565",
"-.	c #A0A0A0",
";.	c #878787",
">.	c #474747",
",.	c #141414",
"'.	c #101010",
").	c #313131",
"                                        ",
"          . + @ # $ %                   ",
"        & * = - ; > , # '               ",
"      ) ! ~ { ] ^ / ( > _ :             ",
"    % ! < / [ } | 1 2 ^ 3 4 5           ",
"  ' 6 7 ^ 8 9 0 a b 1 c d e f           ",
"  g h i [ 8 j k a b } l m n o p         ",
"  q r d s t 0 u a v w x m n > y         ",
"  ) z A l g B B B B C D ^ E F G         ",
"  a 3 ^ H z I J g - K I i ( L M         ",
"  ' N E / m H O P 7 ^ i Q R #           ",
"    ) S T U ] h u V E z ~ W X           ",
"    Y Z F `  .z  .K ` 3 0 ..y           ",
"      Y # +.; @.; ; x #.$.p %.          ",
"        &.*.=.N -.;.>.,.y %.%.%.%.      ",
"            %.%.%.%.%.      5 *.y %.    ",
"                            %.>.$.&.%.  ",
"                              '.B ).G   ",
"                                ,.X G   ",
"                                        "};
end_of_pixmap

$icon_new = <<'end_of_pixmap';
/* XPM */
static char * new_20_xpm[] = {
"20 20 160 2",
"  	c None",
". 	c #2D2D2D",
"+ 	c #4E4E4E",
"@ 	c #545454",
"# 	c #555555",
"$ 	c #505050",
"% 	c #3D3D3D",
"& 	c #727272",
"* 	c #F3F3F3",
"= 	c #FEFEFE",
"- 	c #FFFFFF",
"; 	c #EEEEEE",
"> 	c #B2B2B2",
", 	c #838383",
"' 	c #7E7E7E",
") 	c #EAEAEA",
"! 	c #CFCFCF",
"~ 	c #E3E3E3",
"{ 	c #7D7D7D",
"] 	c #808080",
"^ 	c #FDFDFC",
"/ 	c #FDFDFD",
"( 	c #FCFCFB",
"_ 	c #CDCDCC",
": 	c #F6F6F6",
"< 	c #DADADA",
"[ 	c #666666",
"} 	c #FEFEFD",
"| 	c #FCFCFC",
"1 	c #FBFBFB",
"2 	c #FAFAF9",
"3 	c #FBFBFA",
"4 	c #FAFAF8",
"5 	c #DCDCDB",
"6 	c #B6B6B5",
"7 	c #D0D0D0",
"8 	c #CDCDCD",
"9 	c #A6A6A6",
"0 	c #626262",
"a 	c #F9F9F8",
"b 	c #F8F8F7",
"c 	c #F9F8F7",
"d 	c #F8F7F6",
"e 	c #E3E3E2",
"f 	c #828282",
"g 	c #484747",
"h 	c #3F3E3E",
"i 	c #323232",
"j 	c #242423",
"k 	c #141412",
"l 	c #FAFAFA",
"m 	c #F7F7F6",
"n 	c #F6F6F4",
"o 	c #F7F6F5",
"p 	c #F6F5F4",
"q 	c #F4F4F2",
"r 	c #EAEAE8",
"s 	c #DAD9D7",
"t 	c #CACAC8",
"u 	c #B0AFAD",
"v 	c #878582",
"w 	c #3B3A35",
"x 	c #F5F4F3",
"y 	c #F4F3F2",
"z 	c #F2F2F0",
"A 	c #EDEDEB",
"B 	c #E8E7E5",
"C 	c #E4E3E1",
"D 	c #CFCECC",
"E 	c #ADABA7",
"F 	c #464440",
"G 	c #F3F3F1",
"H 	c #F1F1EF",
"I 	c #F0EFED",
"J 	c #EDECEA",
"K 	c #EBEAE8",
"L 	c #E6E5E3",
"M 	c #DFDEDB",
"N 	c #BAB7B2",
"O 	c #47443F",
"P 	c #F3F2F1",
"Q 	c #F0F0EE",
"R 	c #F1F0EE",
"S 	c #EEEDEB",
"T 	c #EAE9E7",
"U 	c #E6E5E2",
"V 	c #E4E3E0",
"W 	c #E2E1DE",
"X 	c #C9C6C1",
"Y 	c #4D4A45",
"Z 	c #F2F1F0",
"` 	c #EFEFED",
" .	c #EFEEEC",
"..	c #E9E8E6",
"+.	c #E5E4E2",
"@.	c #E4E2DF",
"#.	c #E2E0DD",
"$.	c #CCC9C3",
"%.	c #524F49",
"&.	c #F4F4F3",
"*.	c #EEEEEC",
"=.	c #ECEBE9",
"-.	c #EBEAE7",
";.	c #E7E5E2",
">.	c #E1E0DC",
",.	c #DFDEDA",
"'.	c #CAC8C1",
").	c #514E48",
"!.	c #E8E7E4",
"~.	c #E9E8E5",
"{.	c #E4E3DF",
"].	c #E2E1DD",
"^.	c #DDDCD8",
"/.	c #C9C6C0",
"(.	c #E7E6E3",
"_.	c #E6E4E1",
":.	c #E5E4E0",
"<.	c #DBD9D5",
"[.	c #C6C3BC",
"}.	c #504D47",
"|.	c #7C7C7C",
"1.	c #E3E2DF",
"2.	c #E3E2DE",
"3.	c #DDDBD7",
"4.	c #DCDAD6",
"5.	c #D9D7D3",
"6.	c #C5C1BB",
"7.	c #4F4C46",
"8.	c #E2E1DF",
"9.	c #E1DFDC",
"0.	c #E0DFDB",
"a.	c #E0DEDB",
"b.	c #DFDDD9",
"c.	c #DEDCD8",
"d.	c #D7D6D1",
"e.	c #C3C0B9",
"f.	c #4E4B45",
"g.	c #5B5A59",
"h.	c #BEBBB5",
"i.	c #C2BEB7",
"j.	c #C1BDB6",
"k.	c #C1BDB5",
"l.	c #C0BDB5",
"m.	c #C0BCB4",
"n.	c #BFBBB3",
"o.	c #BEBAB2",
"p.	c #BBB7AF",
"q.	c #BAB6AE",
"r.	c #A6A299",
"s.	c #413E39",
"t.	c #20201F",
"u.	c #353431",
"v.	c #373531",
"w.	c #363431",
"x.	c #363430",
"y.	c #353330",
"z.	c #35332F",
"A.	c #34322F",
"B.	c #2E2C29",
"C.	c #151412",
"                                        ",
"    . + @ # # # # # # $ %               ",
"    & * = - - - - - = ; > ,             ",
"    ' = - - - - - - = ) ! ~ {           ",
"    ] - - - - = ^ / ( ~ _ : < [         ",
"    ] - - } | 1 2 3 4 5 6 7 8 9 0       ",
"    ] = | 3 2 a b c d e f g h i j k     ",
"    ] | l a b m n o p q r s t u v w     ",
"    ] 3 b m n x q x y z A B C D E F     ",
"    ] a p q y G z z H I J K L M N O     ",
"    ] d P z Q Q I R I S T U V W X Y     ",
"    ] o Z H ` `  .I  .J ..+.@.#.$.%.    ",
"    ] &.*.S J J K =.-...;.V >.,.'.).    ",
"    ] z K T ....!.~.!.;.{.].,.^./.).    ",
"    ] ` (.L U U _._.:.{.>.,.^.<.[.}.    ",
"    |.T 1.2.2.2.].].].>.,.3.4.5.6.7.    ",
"    & 8.9.0.0.0.0.a.0.b.c.4.5.d.e.f.    ",
"    g.h.i.j.j.j.k.l.k.m.n.o.p.q.r.s.    ",
"    t.u.v.w.x.x.x.x.x.x.y.z.A.A.B.C.    ",
"                                        "};
end_of_pixmap
		
$icon_report = <<'end_of_pixmap';
/* XPM */
static char * report_20_xpm[] = {
"20 20 95 2",
"  	c None",
". 	c #323232",
"+ 	c #282828",
"@ 	c #292929",
"# 	c #2C2C2C",
"$ 	c #E8E8E8",
"% 	c #ECECEC",
"& 	c #EDEDED",
"* 	c #EEEEEE",
"= 	c #EFEFEF",
"- 	c #F0F0F0",
"; 	c #F1F1F1",
"> 	c #F2F2F2",
", 	c #F3F3F3",
"' 	c #F4F4F4",
") 	c #000000",
"! 	c #E1E1E1",
"~ 	c #E2E2E2",
"{ 	c #E3E3E3",
"] 	c #E5E4E5",
"^ 	c #E6E6E6",
"/ 	c #E7E6E7",
"( 	c #E8E7E8",
"_ 	c #E9E9E9",
": 	c #EAEAEA",
"< 	c #EBECEC",
"[ 	c #ECEDED",
"} 	c #EFEEEF",
"| 	c #F0F0F1",
"1 	c #EDEDEC",
"2 	c #E4E5E5",
"3 	c #777777",
"4 	c #EDEEED",
"5 	c #E4E4E4",
"6 	c #E6E5E5",
"7 	c #E7E7E7",
"8 	c #EEEDEE",
"9 	c #E6E7E6",
"0 	c #EEEEEF",
"a 	c #E5E5E6",
"b 	c #E6E6E7",
"c 	c #E7E7E8",
"d 	c #E9E9E8",
"e 	c #EBEBEB",
"f 	c #EEEDED",
"g 	c #F5F5F5",
"h 	c #E8E8E9",
"i 	c #F7F7F7",
"j 	c #EBECEB",
"k 	c #F6F6F6",
"l 	c #F6F7F7",
"m 	c #F8F8F7",
"n 	c #F8F9F9",
"o 	c #F4F4F5",
"p 	c #F9F9F8",
"q 	c #FAFAFA",
"r 	c #FBFBFB",
"s 	c #EFEFF0",
"t 	c #F4F3F4",
"u 	c #F8F8F8",
"v 	c #F9F8F9",
"w 	c #FAF9F9",
"x 	c #FBFBFC",
"y 	c #FCFCFC",
"z 	c #FDFDFE",
"A 	c #F5F4F4",
"B 	c #F1F0F1",
"C 	c #F1F2F1",
"D 	c #F2F3F3",
"E 	c #F3F3F4",
"F 	c #F7F7F8",
"G 	c #F9F9F9",
"H 	c #FAF9FA",
"I 	c #FEFDFD",
"J 	c #FEFEFE",
"K 	c #F0EFF0",
"L 	c #F1F1F2",
"M 	c #F2F2F3",
"N 	c #F3F4F3",
"O 	c #F4F5F4",
"P 	c #F6F6F5",
"Q 	c #F7F7F6",
"R 	c #FAFAFB",
"S 	c #FCFBFC",
"T 	c #FDFDFD",
"U 	c #9A9A9A",
"V 	c #A2A2A2",
"W 	c #A4A4A4",
"X 	c #A5A5A5",
"Y 	c #A6A6A6",
"Z 	c #A7A7A7",
"` 	c #A8A8A8",
" .	c #A9A9A9",
"..	c #AAAAAA",
"+.	c #979797",
"  . + + + + + @ @ @ @ @ @ @ @ @ #       ",
". $ % & & * * = - - ; > > , ' ' & )     ",
"+ % ! ~ { ] ^ / ( _ : < [ * } - | )     ",
"+ 1 ) ) 2 3 3 3 3 3 3 3 3 3 3 3 > )     ",
"+ 4 { 5 6 7 $ _ : % [ 8 = - ; > , )     ",
"+ * ) ) 9 3 3 3 3 3 3 3 3 3 3 3 ' )     ",
"@ 0 a b c d : e % f 0 - ; > , ' g )     ",
"@ = ) ) h 3 3 3 3 3 3 3 3 3 3 3 i )     ",
"@ - c h : j % & * - ; > , ' k l m )     ",
"@ ; ) ) e 3 3 3 3 3 3 3 3 3 3 3 n )     ",
"@ ; : e % & * - ; > , o g k m p q )     ",
"@ > ) ) & 3 3 3 3 3 3 3 3 3 3 3 r )     ",
"@ , % & * s | > , t g k u v w x y )     ",
"@ ' ) ) = 3 3 3 3 3 3 3 3 3 3 3 z )     ",
"@ A * = B C D E g k F G H r y I J )     ",
"@ % K L M N O P Q u G R S y T J ' )     ",
"# U V W W X X Y Y Z `  .......` +.)     ",
"  ) ) ) ) ) ) ) ) ) ) ) ) ) ) ) )       ",
"                                        ",
"                                        "};
end_of_pixmap

$icon_prefs = <<'end_of_pixmap';
/* XPM */
static char * preferences_20_xpm[] = {
"20 20 138 2",
"  	c None",
". 	c #201F1E",
"+ 	c #242422",
"@ 	c #1B1B19",
"# 	c #97948D",
"$ 	c #B6B3AC",
"% 	c #5E5C57",
"& 	c #2E2D2C",
"* 	c #9F9C98",
"= 	c #C9C6BE",
"- 	c #20201E",
"; 	c #404040",
"> 	c #676767",
", 	c #888786",
"' 	c #D6D4CE",
") 	c #222120",
"! 	c #272727",
"~ 	c #BBBBBB",
"{ 	c #BEBEBE",
"] 	c #373737",
"^ 	c #7C7B77",
"/ 	c #434341",
"( 	c #000000",
"_ 	c #3E3E3E",
": 	c #CDCCCB",
"< 	c #CECBC4",
"[ 	c #72706B",
"} 	c #262626",
"| 	c #D2D2D2",
"1 	c #696969",
"2 	c #111111",
"3 	c #6C6966",
"4 	c #BCB8B1",
"5 	c #CFCBC4",
"6 	c #E7E6E2",
"7 	c #D6D2CD",
"8 	c #BFBBB1",
"9 	c #BBB6AC",
"0 	c #161514",
"a 	c #A8A8A8",
"b 	c #6B6863",
"c 	c #827F79",
"d 	c #878580",
"e 	c #8E8C87",
"f 	c #C7C5C0",
"g 	c #D0CDC6",
"h 	c #8B8781",
"i 	c #272624",
"j 	c #3D3D3D",
"k 	c #959595",
"l 	c #5B5B5A",
"m 	c #DEDDDC",
"n 	c #D4D2CD",
"o 	c #716F6A",
"p 	c #2C2C2C",
"q 	c #616161",
"r 	c #1A1A1A",
"s 	c #494949",
"t 	c #D2D2D1",
"u 	c #545352",
"v 	c #585858",
"w 	c #131313",
"x 	c #393939",
"y 	c #6B6B6B",
"z 	c #3A3937",
"A 	c #4B4A45",
"B 	c #0D0D0C",
"C 	c #171B20",
"D 	c #191D21",
"E 	c #040505",
"F 	c #7C7C7C",
"G 	c #232323",
"H 	c #717170",
"I 	c #898885",
"J 	c #85827C",
"K 	c #151413",
"L 	c #4B5967",
"M 	c #9CB0C6",
"N 	c #99AFC5",
"O 	c #343F4A",
"P 	c #3C3C3C",
"Q 	c #AAAAA9",
"R 	c #CAC8C3",
"S 	c #242321",
"T 	c #53606E",
"U 	c #A9BACC",
"V 	c #99A8B8",
"W 	c #94ACC3",
"X 	c #6B839D",
"Y 	c #575757",
"Z 	c #D9D9D8",
"` 	c #D6D3CE",
" .	c #87847E",
"..	c #292725",
"+.	c #5A6673",
"@.	c #BBC9D8",
"#.	c #8C9DAF",
"$.	c #7A90A6",
"%.	c #859FBB",
"&.	c #637A93",
"*.	c #484848",
"=.	c #EFEEED",
"-.	c #CBC9C5",
";.	c #74726F",
">.	c #1D1C1B",
",.	c #5F6A76",
"'.	c #D1DAE2",
").	c #899BAF",
"!.	c #748BA4",
"~.	c #87A1BC",
"{.	c #69819B",
"].	c #1D232B",
"^.	c #474747",
"/.	c #C0BFBE",
"(.	c #6A6968",
"_.	c #7D7D7C",
":.	c #4B5661",
"<.	c #D2DCE5",
"[.	c #A3AEBC",
"}.	c #7289A2",
"|.	c #8AA5C0",
"1.	c #5E758D",
"2.	c #1B2128",
"3.	c #3A3A3A",
"4.	c #969695",
"5.	c #8F8F8E",
"6.	c #4D5864",
"7.	c #BFCFDF",
"8.	c #A0ACB9",
"9.	c #92ABC4",
"0.	c #5C728A",
"a.	c #202831",
"b.	c #212932",
"c.	c #6A8097",
"d.	c #93A9C0",
"e.	c #6D8299",
"f.	c #242C36",
"g.	c #1D242C",
"            . +                         ",
"          @ # $ %                       ",
"            & * = -             ; >     ",
"              , ' )           ! ~ { ]   ",
"      ^ / ( _ : < [ (         } | 1 2   ",
"      3 4 5 6 7 8 9 0       } a ( (     ",
"      ( b c d e f g h i   j k (         ",
"          ( ( ( l m n o p q r           ",
"                  s t u v w             ",
"                    x y z A B           ",
"              C D E F G H I J K         ",
"            L M N O (   P Q R h S       ",
"          T U V W X (     Y Z `  ...    ",
"        +.@.#.$.%.&.(       *.=.-.;.>.  ",
"      ,.'.).!.~.{.].          ^./.(._.r ",
"    :.<.[.}.|.1.2.              3.4.5.r ",
"    6.7.8.9.0.a.(                 ( (   ",
"    b.c.d.e.a.(                         ",
"      ].f.g.(                           ",
"                                        "};
end_of_pixmap
		
$icon_save_as = <<'end_of_pixmap';
/* XPM */
static char * stock_save_as_20_xpm[] = {
"20 20 274 2",
"  	c None",
". 	c #272309",
"+ 	c #29200C",
"@ 	c #3F3B0F",
"# 	c #C2A930",
"$ 	c #CA9F41",
"% 	c #000000",
"& 	c #2A2E34",
"* 	c #5B6976",
"= 	c #43494F",
"- 	c #6D5C58",
"; 	c #78645F",
"> 	c #78615C",
", 	c #775E5A",
"' 	c #77605B",
") 	c #614C48",
"! 	c #2A270A",
"~ 	c #BCA22F",
"{ 	c #CAA13E",
"] 	c #44381C",
"^ 	c #282A2C",
"/ 	c #4A555F",
"( 	c #2D353B",
"_ 	c #3F474F",
": 	c #A9C6DE",
"< 	c #76818F",
"[ 	c #D0968C",
"} 	c #E69989",
"| 	c #E49586",
"1 	c #E29384",
"2 	c #E19284",
"3 	c #DC9488",
"4 	c #5C4C36",
"5 	c #C1A72F",
"6 	c #C59B3C",
"7 	c #493B1E",
"8 	c #492820",
"9 	c #53565D",
"0 	c #7891A4",
"a 	c #34424D",
"b 	c #3E474E",
"c 	c #A1BED6",
"d 	c #727B8B",
"e 	c #BB7E79",
"f 	c #CE7D73",
"g 	c #CE7B72",
"h 	c #CE7C72",
"i 	c #CD7C71",
"j 	c #3F2623",
"k 	c #D6B932",
"l 	c #CB9A3C",
"m 	c #45361D",
"n 	c #552D27",
"o 	c #8B4E4B",
"p 	c #5E626C",
"q 	c #6D899F",
"r 	c #293742",
"s 	c #3D464D",
"t 	c #9FBED4",
"u 	c #7A858D",
"v 	c #E8E8E8",
"w 	c #FFFFFF",
"x 	c #E9E9E9",
"y 	c #625F4D",
"z 	c #D5B631",
"A 	c #C29643",
"B 	c #433C26",
"C 	c #6E6E6E",
"D 	c #CBCBCB",
"E 	c #E3E3E3",
"F 	c #72777B",
"G 	c #718DA4",
"H 	c #293640",
"I 	c #3E464D",
"J 	c #9FBBD3",
"K 	c #79838C",
"L 	c #DADADA",
"M 	c #EEEEEE",
"N 	c #D7D7D7",
"O 	c #7E785E",
"P 	c #BE9C2C",
"Q 	c #C8963A",
"R 	c #403A23",
"S 	c #5C5C5C",
"T 	c #BCBCBC",
"U 	c #D3D3D3",
"V 	c #E2E2E2",
"W 	c #757A7E",
"X 	c #728EA4",
"Y 	c #26333D",
"Z 	c #9EBAD2",
"` 	c #737D86",
" .	c #D2D2D2",
"..	c #E6E6E6",
"+.	c #989898",
"@.	c #73602A",
"#.	c #BC9341",
"$.	c #42371D",
"%.	c #6B6B6B",
"&.	c #B6B6B6",
"*.	c #CDCDCD",
"=.	c #DEDEDE",
"-.	c #DFDFDF",
";.	c #767B7F",
">.	c #728DA3",
",.	c #24303A",
"'.	c #3D444D",
").	c #9BB7D0",
"!.	c #747F87",
"~.	c #E7E7E7",
"{.	c #4A473D",
"].	c #3E331A",
"^.	c #2D2617",
"/.	c #6D6C6A",
"(.	c #CACACA",
"_.	c #F6F6F6",
":.	c #FDFDFD",
"<.	c #718CA1",
"[.	c #232E37",
"}.	c #3A444C",
"|.	c #96B6CE",
"1.	c #76828A",
"2.	c #C6C6C6",
"3.	c #B0B0B0",
"4.	c #040401",
"5.	c #474642",
"6.	c #828282",
"7.	c #A3A3A3",
"8.	c #BEBEBE",
"9.	c #D4D4D4",
"0.	c #D5D5D5",
"a.	c #D0D0D0",
"b.	c #758C9F",
"c.	c #293035",
"d.	c #3A454C",
"e.	c #98B8CD",
"f.	c #797880",
"g.	c #D5D3D3",
"h.	c #E4E4E4",
"i.	c #C4C4C4",
"j.	c #D6D6D7",
"k.	c #757B80",
"l.	c #71899C",
"m.	c #272D31",
"n.	c #39444C",
"o.	c #91B3CA",
"p.	c #7B95A8",
"q.	c #7D8992",
"r.	c #788289",
"s.	c #747E86",
"t.	c #727B84",
"u.	c #707C87",
"v.	c #768087",
"w.	c #79828A",
"x.	c #79838A",
"y.	c #757E86",
"z.	c #798A98",
"A.	c #6D879A",
"B.	c #1C272E",
"C.	c #39424B",
"D.	c #8EADC8",
"E.	c #7696B0",
"F.	c #7A99B3",
"G.	c #637684",
"H.	c #647686",
"I.	c #6C8091",
"J.	c #708596",
"K.	c #6B8192",
"L.	c #708494",
"M.	c #6E8291",
"N.	c #6D8190",
"O.	c #637B8D",
"P.	c #627585",
"Q.	c #829FB8",
"R.	c #778D9F",
"S.	c #1B252B",
"T.	c #8EADC7",
"U.	c #7090A9",
"V.	c #5B6F7E",
"W.	c #828380",
"X.	c #ADAFAF",
"Y.	c #B8BBBC",
"Z.	c #BDC0C2",
"`.	c #BBBFC0",
" +	c #B1B3B4",
".+	c #A6A9A8",
"++	c #939491",
"@+	c #5F696E",
"#+	c #465A69",
"$+	c #6C8AA3",
"%+	c #8299AA",
"&+	c #1A2328",
"*+	c #39434B",
"=+	c #8BACC4",
"-+	c #6D8DA5",
";+	c #49545D",
">+	c #D3D0CB",
",+	c #C2C1BC",
"'+	c #74716B",
")+	c #6C6962",
"!+	c #C9C7C4",
"~+	c #CAC7C2",
"{+	c #C0BCB5",
"]+	c #C1BDB7",
"^+	c #71777A",
"/+	c #4C6478",
"(+	c #607F94",
"_+	c #879BAB",
":+	c #192025",
"<+	c #38424B",
"[+	c #89AAC3",
"}+	c #6F8EA4",
"|+	c #47525A",
"1+	c #E4E3DF",
"2+	c #C1BFBB",
"3+	c #5A5850",
"4+	c #524F47",
"5+	c #B8B5B0",
"6+	c #C2BEB8",
"7+	c #C6C2BC",
"8+	c #D6D3CF",
"9+	c #7C8286",
"0+	c #4D6578",
"a+	c #607E92",
"b+	c #8A9CAB",
"c+	c #191F24",
"d+	c #323B42",
"e+	c #809DB4",
"f+	c #6C89A0",
"g+	c #3F4950",
"h+	c #DEDCD8",
"i+	c #B3B1AC",
"j+	c #4F4C44",
"k+	c #4E4B43",
"l+	c #B1AEA7",
"m+	c #C8C5BE",
"n+	c #DAD8D3",
"o+	c #E7E6E3",
"p+	c #7C8386",
"q+	c #4B6272",
"r+	c #5B798D",
"s+	c #7F96A8",
"t+	c #1C1F21",
"u+	c #1E2226",
"v+	c #505F6C",
"w+	c #5F7889",
"x+	c #465159",
"y+	c #C0BCB7",
"z+	c #B3AFA9",
"A+	c #88847D",
"B+	c #8C8982",
"C+	c #C1BEB8",
"D+	c #DDDAD7",
"E+	c #DAD8D5",
"F+	c #717679",
"G+	c #485D6F",
"H+	c #597385",
"I+	c #7C95A6",
"J+	c #1B1E20",
"K+	c #14181B",
"L+	c #1D2327",
"M+	c #1B1F22",
"N+	c #383735",
"O+	c #3F3E3B",
"P+	c #413F3D",
"Q+	c #444341",
"R+	c #474645",
"S+	c #4C4C4B",
"T+	c #4D4C4B",
"U+	c #222425",
"V+	c #181F25",
"W+	c #222B32",
"X+	c #3A4853",
"Y+	c #0B0D0E",
"                          . +           ",
"                        @ # $ %         ",
"  & * = - ; > , ' ' ) ! ~ { ] ^ / (     ",
"  _ : < [ } | 1 2 3 4 5 6 7 8 9 0 a     ",
"  b c d e f g h i j k l m n o p q r     ",
"  s t u v w w x y z A B C D E F G H     ",
"  I J K L M N O P Q R S T U V W X Y     ",
"  s Z `  ...+.@.#.$.%.&.*.=.-.;.>.,.    ",
"  '.).!.~.V {.].^./.(.E _.:.w ;.<.[.    ",
"  }.|.1.2.3.4.5.6.7.8.*.9.0.a.;.b.c.    ",
"  d.e.f.g.h.9.i.(.j.V ~.v v ..k.l.m.    ",
"  n.o.p.q.r.s.t.s.u.v.w.x.w.y.z.A.B.    ",
"  C.D.E.F.G.H.I.J.K.L.M.N.O.P.Q.R.S.    ",
"  C.T.U.V.W.X.Y.Z.`. +.+++@+#+$+%+&+    ",
"  *+=+-+;+>+,+'+)+!+~+{+]+^+/+(+_+:+    ",
"  <+[+}+|+1+2+3+4+5+6+7+8+9+0+a+b+c+    ",
"  d+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t+    ",
"  u+v+w+x+y+z+A+B+C+D+o+E+F+G+H+I+J+    ",
"    K+L+M+N+O+P+Q+R+S+T+Q+U+V+W+X+Y+    ",
"                                        "};	
end_of_pixmap
		
$icon_separator = <<'end_of_pixmap';
/* XPM */
static char * sep_20_xpm[] = {
"4 20 3 1",
" 	c None",
".	c #000000",
"+	c #FFFFFF",
"    ",
"    ",
"    ",
"    ",
" .+ ",
" .+ ",
" .+ ",
" .+ ",
" .+ ",
" .+ ",
" .+ ",
" .+ ",
" .+ ",
" .+ ",
" .+ ",
" .+ ",
"    ",
"    ",
"    ",
"    "};
end_of_pixmap

$icon_ensembl = <<'end_of_pixmap';
/* XPM */
static char * ensembl_20_xpm[] = {
"20 20 236 2",
"  	c None",
". 	c #2D2D2D",
"+ 	c #4E4E4E",
"@ 	c #545454",
"# 	c #555555",
"$ 	c #505050",
"% 	c #3D3D3D",
"& 	c #CC0000",
"* 	c #950000",
"= 	c #727272",
"- 	c #F3F3F3",
"; 	c #FEFEFE",
"> 	c #FFFFFF",
", 	c #EEEEEE",
"' 	c #B2B2B2",
") 	c #838383",
"! 	c #DC0000",
"~ 	c #AB0000",
"{ 	c #6A0000",
"] 	c #7E7E7E",
"^ 	c #EAEAEA",
"/ 	c #CFCFCF",
"( 	c #E3E3E3",
"_ 	c #867373",
": 	c #9B0000",
"< 	c #530000",
"[ 	c #808080",
"} 	c #FDFDFC",
"| 	c #FAFCFD",
"1 	c #F4F8FB",
"2 	c #DDE0E3",
"3 	c #CDCDCC",
"4 	c #F6F6F6",
"5 	c #DDB8B8",
"6 	c #B60C0C",
"7 	c #850000",
"8 	c #FEFEFD",
"9 	c #F9FBFC",
"0 	c #DBEAF9",
"a 	c #B2CEF1",
"b 	c #90B4EB",
"c 	c #7EA5E5",
"d 	c #7598D7",
"e 	c #7D93BB",
"f 	c #C7C9CD",
"g 	c #D49F9F",
"h 	c #B70B0B",
"i 	c #620C0C",
"j 	c #FCFCFC",
"k 	c #D6E8F9",
"l 	c #97BBEB",
"m 	c #5F8DD7",
"n 	c #4168AF",
"o 	c #4D68A0",
"p 	c #5A73A8",
"q 	c #3E66B8",
"r 	c #234A9B",
"s 	c #313A4E",
"t 	c #6F2C2C",
"u 	c #9F0202",
"v 	c #360E0E",
"w 	c #131311",
"x 	c #B2CFF2",
"y 	c #6594DB",
"z 	c #3D63AF",
"A 	c #61759F",
"B 	c #9EA3AF",
"C 	c #D6D6D7",
"D 	c #C7D7E6",
"E 	c #5F89D2",
"F 	c #37589E",
"G 	c #9399A4",
"H 	c #D38281",
"I 	c #840D0D",
"J 	c #665C5A",
"K 	c #3B3A35",
"L 	c #B3D0F2",
"M 	c #5385D4",
"N 	c #3A5EA7",
"O 	c #7D8BA7",
"P 	c #D6D5D4",
"Q 	c #EFF2F2",
"R 	c #CADEF3",
"S 	c #80A6DF",
"T 	c #395EA6",
"U 	c #7080A0",
"V 	c #D5D4D3",
"W 	c #DC7E7D",
"X 	c #843A3A",
"Y 	c #A19F9B",
"Z 	c #464440",
"` 	c #718292",
" .	c #6695DB",
"..	c #416BBF",
"+.	c #5573AE",
"@.	c #ABBCD4",
"#.	c #ADC7E7",
"$.	c #83A4D3",
"%.	c #557BB4",
"&.	c #566C9A",
"*.	c #8D95A8",
"=.	c #D6D5D3",
"-.	c #EBE9E7",
";.	c #C76969",
">.	c #9A7C7A",
",.	c #BAB7B2",
"'.	c #47443F",
").	c #5978A1",
"!.	c #4572C7",
"~.	c #3A5DA6",
"{.	c #5E729A",
"].	c #7589A8",
"^.	c #808FA9",
"/.	c #939BAC",
"(.	c #B4B6BD",
"_.	c #D8D7D6",
":.	c #EEEDEB",
"<.	c #EAE9E7",
"[.	c #E6E0DD",
"}.	c #B76968",
"|.	c #C0B9B6",
"1.	c #C9C6C1",
"2.	c #4D4A45",
"3.	c #446DA9",
"4.	c #416CC0",
"5.	c #516B9E",
"6.	c #D0D0D1",
"7.	c #EAEAE8",
"8.	c #EEEEEC",
"9.	c #EFEEEC",
"0.	c #F0EFED",
"a.	c #EDECEA",
"b.	c #E9E8E6",
"c.	c #E3D5D3",
"d.	c #BC8684",
"e.	c #DCDAD7",
"f.	c #CCC9C3",
"g.	c #524F49",
"h.	c #3862AB",
"i.	c #416ABF",
"j.	c #5570A6",
"k.	c #DBDBDB",
"l.	c #EBEAE8",
"m.	c #ECEBE9",
"n.	c #E5E8E7",
"o.	c #D5D8DD",
"p.	c #E7E5E2",
"q.	c #E2D9D6",
"r.	c #CFC1BD",
"s.	c #DFDEDA",
"t.	c #CAC8C1",
"u.	c #514E48",
"v.	c #3960A7",
"w.	c #406BBE",
"x.	c #3E65B4",
"y.	c #93A7D0",
"z.	c #CED7E6",
"A.	c #D5DDE7",
"B.	c #C9D7E5",
"C.	c #ABC4E3",
"D.	c #7B95B6",
"E.	c #A3AEBF",
"F.	c #E8908D",
"G.	c #D55E5C",
"H.	c #CA918E",
"I.	c #DDDCD8",
"J.	c #C9C6C0",
"K.	c #4F6694",
"L.	c #3E67B8",
"M.	c #3D67BC",
"N.	c #406BC0",
"O.	c #547FCD",
"P.	c #5883CB",
"Q.	c #4B6FAC",
"R.	c #637A9F",
"S.	c #9A9EA8",
"T.	c #E0B6B3",
"U.	c #CC2221",
"V.	c #A40101",
"W.	c #912827",
"X.	c #C4BAB6",
"Y.	c #C6C3BC",
"Z.	c #504D47",
"`.	c #71747A",
" +	c #8595B7",
".+	c #4C669C",
"++	c #405E99",
"@+	c #475F92",
"#+	c #687695",
"$+	c #979BA4",
"%+	c #D0CFCB",
"&+	c #E2E1DD",
"*+	c #DFB3B0",
"=+	c #B11211",
"-+	c #9C0101",
";+	c #6F1716",
">+	c #B7B1AE",
",+	c #C5C1BB",
"'+	c #4F4C46",
")+	c #DFDEDC",
"!+	c #C9C8C8",
"~+	c #BDBDBE",
"{+	c #C7C6C4",
"]+	c #D9D8D4",
"^+	c #E0DFDB",
"/+	c #E0DEDB",
"(+	c #DCD5D1",
"_+	c #AA6463",
":+	c #712323",
"<+	c #836765",
"[+	c #D1D0CC",
"}+	c #C3C0B9",
"|+	c #4E4B45",
"1+	c #5B5A59",
"2+	c #BEBBB5",
"3+	c #C2BEB7",
"4+	c #C1BDB6",
"5+	c #C1BDB5",
"6+	c #C0BDB5",
"7+	c #C0BCB4",
"8+	c #BFBBB3",
"9+	c #BEBAB2",
"0+	c #BBB7AF",
"a+	c #BAB6AE",
"b+	c #A6A299",
"c+	c #413E39",
"d+	c #20201F",
"e+	c #353431",
"f+	c #373531",
"g+	c #363431",
"h+	c #363430",
"i+	c #353330",
"j+	c #35332F",
"k+	c #34322F",
"l+	c #2E2C29",
"m+	c #151412",
"                                        ",
"    . + @ # # # # # # $ %       & *     ",
"    = - ; > > > > > ; , ' )   ! ~ {     ",
"    ] ; > > > > > > ; ^ / ( _ & : <     ",
"    [ > > > > ; } | 1 2 3 4 5 6 7       ",
"    [ > > 8 9 0 a b c d e f g h i       ",
"    [ ; j k l m n o p q r s t u v w     ",
"    [ 9 x y z A B C D E F G H I J K     ",
"    [ L M N O P Q R S T U V W X Y Z     ",
"    `  ...+.@.#.$.%.&.*.=.-.;.>.,.'.    ",
"    ).!.~.{.].^./.(._.:.<.[.}.|.1.2.    ",
"    3.4.5.6.7.8.9.0.9.a.b.c.d.e.f.g.    ",
"    h.i.j.k.a.a.l.m.n.o.p.q.r.s.t.u.    ",
"    v.w.x.y.z.A.B.C.D.E.F.G.H.I.J.u.    ",
"    K.L.M.N.O.P.Q.R.S.T.U.V.W.X.Y.Z.    ",
"    `. +.+++@+#+$+%+&+*+=+-+;+>+,+'+    ",
"    = )+!+~+{+]+^+/+^+(+_+:+<+[+}+|+    ",
"    1+2+3+4+4+4+5+6+5+7+8+9+0+a+b+c+    ",
"    d+e+f+g+h+h+h+h+h+h+i+j+k+k+l+m+    ",
"                                        "};
end_of_pixmap
		
$icon_dna_save = <<'end_of_pixmap';
/* XPM */
static char * save_dna_20_xpm[] = {
"20 20 247 2",
"  	c None",
". 	c #2A3034",
"+ 	c #525860",
"@ 	c #4E545A",
"# 	c #756966",
"$ 	c #7E7370",
"% 	c #877D7A",
"& 	c #908784",
"* 	c #9A918F",
"= 	c #A29A98",
"- 	c #ACA5A3",
"; 	c #B6AFAE",
"> 	c #BEB8B7",
", 	c #C8C3C2",
"' 	c #D1CDCC",
") 	c #D2D4D5",
"! 	c #E1E2E3",
"~ 	c #E6E7E8",
"{ 	c #4E575F",
"] 	c #B6CFE3",
"^ 	c #899EB0",
"/ 	c #E6B0A7",
"( 	c #E7B5AC",
"_ 	c #E9BAB2",
": 	c #EBC0B8",
"< 	c #ECC5BE",
"[ 	c #EECAC4",
"} 	c #F0CFCA",
"| 	c #F1D5D0",
"1 	c #F3DAD5",
"2 	c #F5DFDB",
"3 	c #F6E4E1",
"4 	c #E0E2E3",
"5 	c #EDF2F5",
"6 	c #EAEBEC",
"7 	c #484E56",
"8 	c #A6C2D5",
"9 	c #798D9C",
"0 	c #D89B93",
"a 	c #DBA29B",
"b 	c #DDA9A2",
"c 	c #E0B0AA",
"d 	c #E3B8B2",
"e 	c #E6BFB9",
"f 	c #E8C6C1",
"g 	c #EBCDC9",
"h 	c #EED4D1",
"i 	c #F1DBD8",
"j 	c #EAEEF2",
"k 	c #E9E9EA",
"l 	c #A7C2D5",
"m 	c #929AA2",
"n 	c #E3E4E2",
"o 	c #F7EEEE",
"p 	c #F8EFEF",
"q 	c #F8F0F0",
"r 	c #F9F1F1",
"s 	c #F9F3F3",
"t 	c #FAF4F4",
"u 	c #FAF5F5",
"v 	c #FBF6F6",
"w 	c #FBF7F7",
"x 	c #FCF9F9",
"y 	c #E8ECF0",
"z 	c #E8E9E9",
"A 	c #F9FBF9",
"B 	c #FAFBF9",
"C 	c #FAFBFA",
"D 	c #7C8BDE",
"E 	c #6E6F72",
"F 	c #726F6F",
"G 	c #6F6F6F",
"H 	c #737070",
"I 	c #737171",
"J 	c #766765",
"K 	c #862A21",
"L 	c #B6382B",
"M 	c #E4D3D2",
"N 	c #939AA3",
"O 	c #D3D2CE",
"P 	c #E5E6E4",
"Q 	c #ECEDEA",
"R 	c #E9E9E7",
"S 	c #7483E0",
"T 	c #5769DE",
"U 	c #DCDFEA",
"V 	c #F0F1EF",
"W 	c #F2F3F2",
"X 	c #F2F2F1",
"Y 	c #F3F4F3",
"Z 	c #DEE0E1",
"` 	c #DB675A",
" .	c #DB7E77",
"..	c #878C93",
"+.	c #F7EDED",
"@.	c #F1F2F8",
"#.	c #5E6FE1",
"$.	c #485CD3",
"%.	c #484F80",
"&.	c #757675",
"*.	c #747576",
"=.	c #757676",
"-.	c #696A6B",
";.	c #6D6E71",
">.	c #BF8B87",
",.	c #E7E8E6",
"'.	c #EAEBE9",
").	c #8F9BE4",
"!.	c #455ADD",
"~.	c #5568DE",
"{.	c #8A8FB1",
"].	c #707070",
"^.	c #676869",
"/.	c #CD7C6C",
"(.	c #DD958C",
"_.	c #868C92",
":.	c #FAFCFA",
"<.	c #FBFCFA",
"[.	c #FBFCFB",
"}.	c #E2E5F6",
"|.	c #5768DF",
"1.	c #4156DC",
"2.	c #6F75DA",
"3.	c #DA775D",
"4.	c #DA6544",
"5.	c #E6DEDD",
"6.	c #95A9BB",
"7.	c #99A0A8",
"8.	c #959AA0",
"9.	c #9DA1A7",
"0.	c #A4A8AD",
"a.	c #ABAFB4",
"b.	c #B3B7BB",
"c.	c #BBBEC2",
"d.	c #C2C5C8",
"e.	c #CACACD",
"f.	c #9A6F98",
"g.	c #C85F55",
"h.	c #A05B80",
"i.	c #ACACDF",
"j.	c #E7E8E8",
"k.	c #9AB4CA",
"l.	c #94A9BA",
"m.	c #9BAFBF",
"n.	c #98AABA",
"o.	c #9FB0BE",
"p.	c #A7B6C3",
"q.	c #AEBCC8",
"r.	c #B6C2CD",
"s.	c #BDC8D2",
"t.	c #D5AEA1",
"u.	c #DF7F5D",
"v.	c #D16846",
"w.	c #8E6163",
"x.	c #3645AC",
"y.	c #3647B4",
"z.	c #C6CAE5",
"A.	c #8294A2",
"B.	c #9A9D9D",
"C.	c #AEB0AF",
"D.	c #CCCCCB",
"E.	c #D0D0CF",
"F.	c #D8B7AD",
"G.	c #E08768",
"H.	c #E28664",
"I.	c #DE8B6C",
"J.	c #BCB3B4",
"K.	c #9DA1A4",
"L.	c #ACB0B3",
"M.	c #6471BF",
"N.	c #909BE2",
"O.	c #8A9FB0",
"P.	c #6C737A",
"Q.	c #D6D5D2",
"R.	c #D9D8D5",
"S.	c #B4B6B5",
"T.	c #C78B7D",
"U.	c #DF6F51",
"V.	c #DF8468",
"W.	c #796058",
"X.	c #666666",
"Y.	c #616161",
"Z.	c #5F6164",
"`.	c #67696D",
" +	c #5960B6",
".+	c #A0A6E3",
"++	c #4D565E",
"@+	c #959390",
"#+	c #D9654B",
"$+	c #C98373",
"%+	c #C6C7C8",
"&+	c #CBCCCA",
"*+	c #CDCDCC",
"=+	c #C9CACA",
"-+	c #969BDB",
";+	c #5960DF",
">+	c #5962DF",
",+	c #E0E1E7",
"'+	c #6B7279",
")+	c #EAECE9",
"!+	c #C8C8C7",
"~+	c #C1665B",
"{+	c #747473",
"]+	c #7A7B7C",
"^+	c #7C7D7E",
"/+	c #5E5E90",
"(+	c #4F50BD",
"_+	c #535ADF",
":+	c #7378E0",
"<+	c #E1E6F0",
"[+	c #292F33",
"}+	c #718696",
"|+	c #9FA09E",
"1+	c #C05452",
"2+	c #7A3434",
"3+	c #786D89",
"4+	c #7777E4",
"5+	c #6565E1",
"6+	c #6D6FDE",
"7+	c #B0B7DA",
"8+	c #DDE2E6",
"9+	c #EAECED",
"0+	c #353C3F",
"a+	c #4D5359",
"b+	c #585E63",
"c+	c #83817D",
"d+	c #8C8A87",
"e+	c #9EA09E",
"f+	c #A99EA0",
"g+	c #9E5796",
"h+	c #7770DF",
"i+	c #9363A9",
"j+	c #B38FA9",
"k+	c #BCBEC0",
"l+	c #C7C8CA",
"m+	c #D6D7D9",
"n+	c #E2E5E8",
"o+	c #6360E1",
"p+	c #6B67E2",
"q+	c #6A4E8B",
"r+	c #942D22",
"s+	c #D7402D",
"t+	c #D34434",
"u+	c #5F5EE0",
"v+	c #D7402F",
"w+	c #D63B2F",
"x+	c #D34441",
"                                        ",
"  . + @ # $ % & * = - ; > , ' ) ! ~     ",
"  { ] ^ / ( _ : < [ } | 1 2 3 4 5 6     ",
"  7 8 9 / 0 a b c d e f g h i 4 j k     ",
"  { l m n o p q r s t u v w x 4 y z     ",
"  7 8 m n A B C D E F G H I J K L M     ",
"  7 8 N O P Q R S T U V W X Y Z `  .    ",
"  7 8 ..+.A B C @.#.$.%.&.*.=.-.;.>.    ",
"  7 8 N O P ,.R '.V ).!.~.{.].^./.(.    ",
"  7 8 _.n A B C :.<.[.}.|.1.2.3.4.5.    ",
"  7 8 6.7.8.9.0.a.b.c.d.e.f.g.h.i.j.    ",
"  7 k.l.m.n.o.p.q.r.s.t.u.v.w.x.y.z.    ",
"  7 k.l.A.B.C.D.E.F.G.H.I.J.K.L.M.N.    ",
"  7 k.O.P.Q.R.S.T.U.V.W.X.Y.Z.`. +.+    ",
"  ++k.^ P.o R.@+#+$+%+&+*+=+-+;+>+,+    ",
"  7 k.^ '+)+!+@+~+{+]+^+/+(+_+:+<+~     ",
"  [+}+9 '+Q.!+|+1+2+3+4+5+6+7+8+9+j.    ",
"    0+a+b+c+d+e+f+g+h+i+j+k+l+m+n+j.    ",
"                o+p+q+r+s+t+            ",
"                u+        v+w+x+        "};
end_of_pixmap

$icon_dna_open = <<'end_of_pixmap';
/* XPM */
static char * open_dna_20_iii_xpm[] = {
"20 20 214 2",
"  	c None",
". 	c #1C1F1B",
"+ 	c #292C28",
"@ 	c #373936",
"# 	c #444643",
"$ 	c #515350",
"% 	c #414138",
"& 	c #999A8B",
"* 	c #A9AA9F",
"= 	c #AEAFA5",
"- 	c #B3B5AB",
"; 	c #A7AB9E",
"> 	c #888883",
", 	c #5569D6",
"' 	c #0E0F15",
") 	c #000000",
"! 	c #210B08",
"~ 	c #84251B",
"{ 	c #B53528",
"] 	c #646563",
"^ 	c #D9DCBD",
"/ 	c #DBDEC1",
"( 	c #DDE0C5",
"_ 	c #E0E2C9",
": 	c #D6D8C0",
"< 	c #A2A296",
"[ 	c #80817F",
"} 	c #4B60DE",
"| 	c #485CDD",
"1 	c #D74130",
"2 	c #D5463A",
"3 	c #43433B",
"4 	c #BDBFB0",
"5 	c #CDD0B2",
"6 	c #D0D2B6",
"7 	c #D3D5BB",
"8 	c #DCDEC7",
"9 	c #BCBFAD",
"0 	c #8A8A85",
"a 	c #4B5EDE",
"b 	c #485CD3",
"c 	c #212A65",
"d 	c #040407",
"e 	c #030405",
"f 	c #020102",
"g 	c #A44F48",
"h 	c #C5C8B8",
"i 	c #BFC2A2",
"j 	c #BBBEA4",
"k 	c #BEC2A9",
"l 	c #CACCB3",
"m 	c #CED0B8",
"n 	c #C3C6B5",
"o 	c #AEAFA3",
"p 	c #9B9B97",
"q 	c #A6A6A2",
"r 	c #7580C9",
"s 	c #4459DC",
"t 	c #5063DA",
"u 	c #787D9E",
"v 	c #616160",
"w 	c #010102",
"x 	c #C65E48",
"y 	c #D65848",
"z 	c #37372E",
"A 	c #CFD2B1",
"B 	c #AAAF95",
"C 	c #A5AA92",
"D 	c #B4B8A2",
"E 	c #B9BDA8",
"F 	c #BEC1AE",
"G 	c #C8CBBB",
"H 	c #CDD0C1",
"I 	c #DADBD0",
"J 	c #DEDFD5",
"K 	c #CBCEDA",
"L 	c #5466DC",
"M 	c #4156DC",
"N 	c #6B71D6",
"O 	c #D8755B",
"P 	c #DA5F3D",
"Q 	c #242521",
"R 	c #C5C7B7",
"S 	c #959A81",
"T 	c #6D6F64",
"U 	c #61615A",
"V 	c #6B6B64",
"W 	c #767670",
"X 	c #80807B",
"Y 	c #8B8B86",
"Z 	c #959591",
"` 	c #A0A09B",
" .	c #ABABA7",
"..	c #B5B5B1",
"+.	c #C0BEBC",
"@.	c #9A6E97",
"#.	c #C85F55",
"$.	c #A05B80",
"%.	c #AAA8D8",
"&.	c #B5B7A2",
"*.	c #717260",
"=.	c #C8CABD",
"-.	c #CCCDC1",
";.	c #CFD1C6",
">.	c #D8DACF",
",.	c #DBDDD3",
"'.	c #DEE0D7",
").	c #DEDFD7",
"!.	c #E5E6DF",
"~.	c #E5B8A5",
"{.	c #E07F5E",
"].	c #D16846",
"^.	c #926465",
"/.	c #3645AC",
"(.	c #3647B4",
"_.	c #D1D5F0",
":.	c #56564F",
"<.	c #B9BCB1",
"[.	c #DFE1CB",
"}.	c #E1E3CF",
"|.	c #E8EAD8",
"1.	c #E4C4B0",
"2.	c #E18968",
"3.	c #E28664",
"4.	c #E28F6F",
"5.	c #DED1C7",
"6.	c #B7B8B3",
"7.	c #BBBBB7",
"8.	c #6773BF",
"9.	c #97A1E8",
"0.	c #9EA093",
"a.	c #A3A599",
"b.	c #DDDECE",
"c.	c #DADCC2",
"d.	c #E0A58D",
"e.	c #DF6F51",
"f.	c #DF8568",
"g.	c #7A6257",
"h.	c #6D6D69",
"i.	c #6D6E69",
"j.	c #6E6E6B",
"k.	c #6E6E6F",
"l.	c #5B61B6",
"m.	c #A7ACE9",
"n.	c #9EA094",
"o.	c #D7D9BE",
"p.	c #DCDEC6",
"q.	c #DCDECA",
"r.	c #DC684D",
"s.	c #C98470",
"t.	c #C8CABF",
"u.	c #C9CAC0",
"v.	c #CED0C7",
"w.	c #DADBD2",
"x.	c #A6A8E4",
"y.	c #5A61DF",
"z.	c #5962DE",
"A.	c #EBECF2",
"B.	c #1A1B18",
"C.	c #B0B4A7",
"D.	c #C7C9AD",
"E.	c #D2D4BD",
"F.	c #D5D7C3",
"G.	c #CE7566",
"H.	c #787970",
"I.	c #767873",
"J.	c #797A77",
"K.	c #595A89",
"L.	c #4F50BD",
"M.	c #535ADF",
"N.	c #767AE0",
"O.	c #E0E1E7",
"P.	c #BFC2A9",
"Q.	c #C65A55",
"R.	c #7A3433",
"S.	c #716680",
"T.	c #7575E1",
"U.	c #6565E1",
"V.	c #6F71DF",
"W.	c #C1C4E0",
"X.	c #E6E6E3",
"Y.	c #E4E4E4",
"Z.	c #B0B4A8",
"`.	c #939888",
" +	c #909282",
".+	c #97998A",
"++	c #9EA092",
"@+	c #A6A79B",
"#+	c #A9AA9E",
"$+	c #B2A8A0",
"%+	c #9E5796",
"&+	c #7770DF",
"*+	c #9363A9",
"=+	c #B592A9",
"-+	c #D1D2CE",
";+	c #DADBD8",
">+	c #DFDFDE",
",+	c #0E110D",
"'+	c #363835",
")+	c #454744",
"!+	c #525451",
"~+	c #5F615E",
"{+	c #6B6C70",
"]+	c #6562D7",
"^+	c #6B67DF",
"/+	c #77648D",
"(+	c #942F24",
"_+	c #D7402D",
":+	c #D05142",
"<+	c #C9A7A6",
"[+	c #5F5EE0",
"}+	c #D7402F",
"|+	c #D63B2F",
"1+	c #D34441",
"                                        ",
"                                        ",
"    . + @ # $                           ",
"  % & * = - ; >   , ' ) ) ) ) ! ~ {     ",
"  ] ^ / ( _ : < [ } |             1 2   ",
"3 4 ^ 5 6 7 8 9 0   a b c ) d e ) f g   ",
"3 h i j k l m n o p q r s t u v w x y   ",
"z A B C D E F n G H I J K L M N O P     ",
"Q R S T U V W X Y Z `  ...+.@.#.$.%.    ",
"Q &.*.* =.-.;.>.,.'.).!.~.{.].^./.(._.  ",
"Q 4 :.<.( _ 8 [.}.|.1.2.3.4.5.6.7.8.9.  ",
"Q 0.a.b.( c.8 [.}.d.e.f.g.h.i.j.k.l.m.  ",
"Q n.a.b.o.7 p.[.q.r.s.t.u.v.w.x.y.z.A.  ",
"B.C.G 5 D.l m E.F.G.H.I.J.K.L.M.N.O.    ",
"B.4 G j P.E F n G Q.R.S.T.U.V.W.X.Y.    ",
"Q Z.S `. +.+++@+#+$+%+&+*+=+-+;+>+      ",
"  ,+. + '+)+!+~+{+]+^+/+(+_+:+<+        ",
"                  [+        }+|+1+      ",
"                                        ",
"                                        "};
end_of_pixmap

									
$dna_canvas_pixmap = <<'end_of_pixmap';		
/* XPM */
static char * perlprimer_dna_canvas_2_aa_xpm[] = {
"528 52 306 2",
"  	c None",
". 	c #EEEEEE",
"+ 	c #000000",
"@ 	c #1B1B1B",
"# 	c #1C1C1C",
"$ 	c #545454",
"% 	c #DBDBDB",
"& 	c #DEDEDE",
"* 	c #C8C8C8",
"= 	c #484848",
"- 	c #1A1A1A",
"; 	c #181818",
"> 	c #424242",
", 	c #C4C4C4",
"' 	c #E8E8E8",
") 	c #D0D0D0",
"! 	c #191919",
"~ 	c #9A9A9A",
"{ 	c #ECECEC",
"] 	c #A2A2A2",
"^ 	c #171717",
"/ 	c #E3E3E3",
"( 	c #212121",
"_ 	c #595959",
": 	c #7C7C7C",
"< 	c #858585",
"[ 	c #474747",
"} 	c #EEE7E7",
"| 	c #F2AFAF",
"1 	c #525252",
"2 	c #555555",
"3 	c #DFDFDF",
"4 	c #4F4F4F",
"5 	c #AEAEAE",
"6 	c #9B9B9B",
"7 	c #060606",
"8 	c #E5E5E5",
"9 	c #5D5D5D",
"0 	c #8B8B8B",
"a 	c #C7C7C7",
"b 	c #373737",
"c 	c #151515",
"d 	c #777777",
"e 	c #333333",
"f 	c #626262",
"g 	c #2A2A2A",
"h 	c #C3C3C3",
"i 	c #A5A5A5",
"j 	c #6E6E6E",
"k 	c #232323",
"l 	c #131313",
"m 	c #878787",
"n 	c #676767",
"o 	c #DDDDDD",
"p 	c #4E4E4E",
"q 	c #6B6B6B",
"r 	c #EBEBEB",
"s 	c #D2D2D2",
"t 	c #070707",
"u 	c #7E7E7E",
"v 	c #727272",
"w 	c #141414",
"x 	c #848484",
"y 	c #353535",
"z 	c #5E5E5E",
"A 	c #EFD6D6",
"B 	c #F67B7B",
"C 	c #FE0000",
"D 	c #D7D7D7",
"E 	c #616161",
"F 	c #EAEAEA",
"G 	c #606060",
"H 	c #505050",
"I 	c #909090",
"J 	c #767676",
"K 	c #585858",
"L 	c #161616",
"M 	c #ACACAC",
"N 	c #EDEDED",
"O 	c #5F5F5F",
"P 	c #969696",
"Q 	c #C9C9C9",
"R 	c #E4E4E4",
"S 	c #444444",
"T 	c #7F7F7F",
"U 	c #888888",
"V 	c #707070",
"W 	c #8E8E8E",
"X 	c #B1B1B1",
"Y 	c #ADADAD",
"Z 	c #949494",
"` 	c #464646",
" .	c #8A8A8A",
"..	c #030303",
"+.	c #D8D8D8",
"@.	c #E9E9E9",
"#.	c #818181",
"$.	c #BABABA",
"%.	c #E6E6E6",
"&.	c #494949",
"*.	c #7A7A7A",
"=.	c #F1C4C4",
"-.	c #F86262",
";.	c #FB3030",
">.	c #F39B9B",
",.	c #0A0A0A",
"'.	c #CFCFCF",
").	c #090909",
"!.	c #D3D3D3",
"~.	c #272727",
"{.	c #1F1F1F",
"].	c #808080",
"^.	c #6C6C6C",
"/.	c #B9B9B9",
"(.	c #383838",
"_.	c #A9A9A9",
":.	c #565656",
"<.	c #DCDCDC",
"[.	c #0B0B0B",
"}.	c #E2E2E2",
"|.	c #0F0F0F",
"1.	c #C6C6C6",
"2.	c #D5D5D5",
"3.	c #D6D6D6",
"4.	c #F2B5B5",
"5.	c #CDCDCD",
"6.	c #0D0D0D",
"7.	c #A6A6A6",
"8.	c #414141",
"9.	c #CECECE",
"0.	c #D1D1D1",
"a.	c #E7E7E7",
"b.	c #E0E0E0",
"c.	c #5C5C5C",
"d.	c #7B7B7B",
"e.	c #DADADA",
"f.	c #F48F8F",
"g.	c #FD0E0E",
"h.	c #6D6D6D",
"i.	c #252525",
"j.	c #828282",
"k.	c #C0C0C0",
"l.	c #8C8C8C",
"m.	c #4C4C4C",
"n.	c #787878",
"o.	c #919191",
"p.	c #3B3B3B",
"q.	c #4B4B4B",
"r.	c #E1E1E1",
"s.	c #303030",
"t.	c #F0CCCC",
"u.	c #9C9C9C",
"v.	c #A7A7A7",
"w.	c #111111",
"x.	c #AAAAAA",
"y.	c #8F8F8F",
"z.	c #8D8D8D",
"A.	c #B8B8B8",
"B.	c #323232",
"C.	c #242424",
"D.	c #B0B0B0",
"E.	c #393939",
"F.	c #1D1D1D",
"G.	c #CCCCCC",
"H.	c #EEE0E0",
"I.	c #FF0000",
"J.	c #959595",
"K.	c #F49494",
"L.	c #636363",
"M.	c #BCBCBC",
"N.	c #F0CDCD",
"O.	c #F1BABA",
"P.	c #1E1E1E",
"Q.	c #2D2D2D",
"R.	c #B6B6B6",
"S.	c #F76666",
"T.	c #8B5A00",
"U.	c #A44600",
"V.	c #CF2400",
"W.	c #9A4D00",
"X.	c #FFA500",
"Y.	c #FFA100",
"Z.	c #FF7200",
"`.	c #FF6600",
" +	c #FF9400",
".+	c #104E8B",
"++	c #FF6A00",
"@+	c #FF5A00",
"#+	c #FF9B00",
"$+	c #1E90FF",
"%+	c #FFA400",
"&+	c #FF9900",
"*+	c #FF9C00",
"=+	c #191970",
"-+	c #916D6D",
";+	c #391560",
">+	c #9D6060",
",+	c #F60707",
"'+	c #551252",
")+	c #F00107",
"!+	c #956868",
"~+	c #E41919",
"{+	c #B14C4C",
"]+	c #827B7B",
"^+	c #8B7272",
"/+	c #CA3434",
"(+	c #070000",
"_+	c #650000",
":+	c #CB0000",
"<+	c #2D0000",
"[+	c #430000",
"}+	c #EF0000",
"|+	c #F1BEBE",
"1+	c #FE0E0E",
"2+	c #596ABC",
"3+	c #F00910",
"4+	c #248CF8",
"5+	c #77569A",
"6+	c #3481E6",
"7+	c #F85E5E",
"8+	c #EFDDDD",
"9+	c #536DC2",
"0+	c #FC0103",
"a+	c #AF325A",
"b+	c #4B73CC",
"c+	c #2D85ED",
"d+	c #4E3966",
"e+	c #F00408",
"f+	c #553762",
"g+	c #A01E37",
"h+	c #F10407",
"i+	c #3F3E6F",
"j+	c #204881",
"k+	c #FF7900",
"l+	c #FF0A00",
"m+	c #FF8400",
"n+	c #FF4100",
"o+	c #FF0900",
"p+	c #A94200",
"q+	c #F70500",
"r+	c #935300",
"s+	c #A24800",
"t+	c #D12300",
"u+	c #F80400",
"v+	c #EEEDED",
"w+	c #F0D2D2",
"x+	c #FE0D0D",
"y+	c #F76B6B",
"z+	c #FB2F2F",
"A+	c #EEE8E8",
"B+	c #F1B9B9",
"C+	c #F49999",
"D+	c #FE0505",
"E+	c #F3A5A5",
"F+	c #FE0303",
"G+	c #FE0404",
"H+	c #F85555",
"I+	c #F48E8E",
"J+	c #F85C5C",
"K+	c #979797",
"L+	c #C1C1C1",
"M+	c #515151",
"N+	c #202020",
"O+	c #050505",
"P+	c #C2C2C2",
"Q+	c #3C3C3C",
"R+	c #343434",
"S+	c #282828",
"T+	c #939393",
"U+	c #363636",
"V+	c #5B5B5B",
"W+	c #D4D4D4",
"X+	c #0E0E0E",
"Y+	c #B7B7B7",
"Z+	c #747474",
"`+	c #989898",
" @	c #BEBEBE",
".@	c #9F9F9F",
"+@	c #222222",
"@@	c #2F2F2F",
"#@	c #535353",
"$@	c #BDBDBD",
"%@	c #CBCBCB",
"&@	c #D9D9D9",
"*@	c #010101",
"=@	c #717171",
"-@	c #BFBFBF",
";@	c #292929",
">@	c #3E3E3E",
",@	c #313131",
"'@	c #9E9E9E",
")@	c #9D9D9D",
"!@	c #5A5A5A",
"~@	c #080808",
"{@	c #B2B2B2",
"]@	c #040404",
"^@	c #BBBBBB",
"/@	c #A4A4A4",
"(@	c #6F6F6F",
"_@	c #CACACA",
":@	c #404040",
"<@	c #646464",
"[@	c #0C0C0C",
"}@	c #2B2B2B",
"|@	c #868686",
"1@	c #A8A8A8",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + @ # # # $ . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . # & . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . * = - ; > , . . . . . . . . . . ' . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ) ! ~ { . ] ^ * . . . . . . . . / ( . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . _ : . . . . < [ . . . . . . . . % + . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . } | . . . + % . . . . & 1 ! - 2 3 . . + 4 ; + 5 . 6 7 8 . 9 0 a b # c d . . + 4 ; < e # f + % . . . . + f # g h . . + 4 ; i + % . + j k l m $ c n . o p @ - q r . + 4 ; i . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ; , . . . . s t r . + % . . + % @ + b o p @ - q r . + 4 ; i . . . + 4 ; u b # c d . . + v k w x . s y # z + % o p @ - q r . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . A B C | . . + @ # # # D 4 E F F G H . . + I . J q . K L M N ; D O P . Q + R . + I . S T . U + % . . . . + V . W e . . + I . { + % . + : . X + Y s + R H Z . 8 `  .. + I . { . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..+.. . . . @.+ 3 . + % . . + % % + . H Z . 8 `  .. + I . { . . . + I . 9 P . Q + R . + #.. $...%.&.*.. U + % H Z . 8 `  .. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . =.-.;.>.. . . + % . . . . ,.'.. . '.).. . + !.. , ~.. {.].^./.(.. r _.u :.+ <.. + !.. [.s . }.+ % . . . . + ) . R + 8 . + !.. . + % . + s . % + F . + % [.@ # # |.H . + !.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ; 1.. . . . 2.7 r . + % . . + % % + . [.@ # # |.H . + !.. . . . . + !.. r _.u :.+ <.. + 3.. N + <.[.s . / + % [.@ # # |.H . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . } 4.;.-.=.. . . . + % . . . . ,.'.. . 5.6.. . + % . . ( 7.8.* g V < . ` :.] 9.+ % . + % . [.0.. 8 + % . . . . + s . / + a.. + % . . + % . + % . % + . . + % ).9.. . b.8 . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . c.d.. . . .  .` . . + e.. %.+ % % + . ).9.. . b.8 . + % . . . . . + % . ` :.] 9.+ % . + % . . + % ).s . R + % ).9.. . b.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . A f.g.f.A . . . . . + % . . . . p E F F G 4 . . + % . . h.i.j.. # @ s . t k.. l.+ 3.. + % . m.n.. o.+ % . . . . + u . 0 p.. . + % . . + % . + % . % + . . + % q.f { r.s.U . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . } t.. !.- u.N . v.^ a . . w.x.. y.+ % & + @.q.f { r.s.U . + % . . . . . + % . t k.. l.+ 3.. + % . . + % = u . z.+ <.q.f { r.s.U . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . =.-.;.4.} . . . . . . + % . . . . <.H ! # 2 o . . + % . . A.+ h . 9 B.. . ].c C.d ~.D.. + % . 2.E.# O + % . . . . + _ F.e G.. . + % . . + % . + % . % + . . + % % m.! L n { . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . H.f.I.t.. 5.= ! - > , . . . J.; k u + % N s.C.% m.! L n { . + % . . . . . + % . ].c C.d ~.D.. + % . . + % s e @ :.+ r.% m.! L n { . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4.;.-.=.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . A K.C 4.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . L.M.. 5 @ . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . t.g.f.A . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . N.-.C O.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . l.P.F.Q.R.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . I.I.. . . . t.} . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . N.S.C =.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . I.I.I.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . =.C C N.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . T.T.T.T.I.I.I.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.U.C V.W.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . T.X.X.I.I.I.I.I.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.Y.Z.I.`. +T.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . T.X.X.I.I.I.I.I.X.X.X.X.X.X.X..+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+X.X.X.X.X.X.X.X.Y.++I.@+#+X.T.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . + + + + + . . + . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . T.X.X.X.X.X.X.X.X.X.X.X.X.X.X..+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+.+X.X.X.X.X.X.X.X.%+&+++*+X.X.T.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + + + . . . + . . . . . ",
". . . . . . + . . . . . . + . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . T.X.X.X.X.X.X.X.X.X.X.X.X.X.X..+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+.+X.X.X.X.X.X.X.X.X.%+Y.X.X.X.T.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + . . . + . . + . . . . . ",
". . . . . . + . . . . . . + . . . . . . + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + . . . . . . . . . + . . + . . . . . ",
". . . . . . + + + + . . . . . . . . . . + T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T + =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=++ T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T + . . . . . . . + + . . . . . . . . . ",
". . . . . . . . . . + . . . . . . . . . + T T T T T T T T T T T T T T T T T T T T T T T T T T T T T -+T T T T T T T T T T T T T T T T T T T T T + =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+;+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=++ T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T + . . . . . . . . . + . . . . . . . . ",
". . . . . . . . . . + . . . . . . . . . + T T T T T T T T T T T T T T T T T T T T T T T T T T T T >+,+-+T T T T T T T T T T T T T T T T T T T T + =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+'+)+;+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=++ T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T + . . . . . . . . . + . . . . . . . . ",
". . . . . . + . . . + . . . . . . . . . + T T T T T T T T T T T T T T T T T T T T T T T T T T T !+~+{+]+T T T T T T T T T T T T T T T T T T T T + =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+'+)+'+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=++ T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T + . . . . . + . . . + . . . . . . . . ",
". . . . . . + . . . + . . . . . . . . . + T T T T T T T T T T T T T T T T T T T T T T T T T T ^+/+/+^+T T T T T T T T T T T T T T T T T T T T T + =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+'+)+'+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=++ T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T T + . . . . . + . . . + . . . . . . . . ",
". . . . . . . + + + . . . . . . . . . . + + + + + + + + + + + + + + + + + + + + + + + + + + (+_+:+<++ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + [+}+[++ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + (++ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + . . . . . . + + + . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . |+1+4.. . . . . . . . . . . . . . . . . . . . T.X.X.X.X.X.X.X.X.X.X.X.X.X.X..+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+2+3+2+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+4+5+6+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+.+X.X.X.X.X.X.X.X.X.X.X.X.X.X.T.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . N.S.7+8+. . . . . . . . . . . . . . . . . . . . T.X.X.X.X.X.X.X.X.X.X.X.X.X.X..+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+2+3+2+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+9+0+a+b+c+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+$+I.I.I.I.I.X.X.X.X.X.X.X.X.X.X.X.X.X.X.T.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 8+7+S.N.. . . . . . . . . . . . . . . . . . . . . T.X.X.X.X.X.X.X.X.X.X.X.X.X.X..+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+d+e+d+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+f+g+h+g+i+j+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+I.I.I.I.I.I.X.X.X.X.X.X.X.X.X.X.X.X.X.X.T.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4.1+|+. . . . . . . . . . . . . . . . . . . . . . T.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.k+l+k+X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.&+m+n+o+n+m+&+X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.I.I.I.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.X.T.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . =.;.f.} . . . . . . . . . . . . . . . . . . . . . . T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.p+q+p+T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.r+s+t+u+t+s+r+T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.I.I.I.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.T.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . v+w+-.-.A . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . } | g.| . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 8+|+7+x+7+|+8+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . I.I.I.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . v+A y+z+=.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . } B g.| . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 8+|+7+x+7+|+8+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . A+B+8+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . A+C+D+E+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4.F+>.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 8+|+7+x+7+|+8+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . A+C+G+H+|+8+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . A+B+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 8+|+7+x+7+|+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . v+8+I+J+x+7+|+8+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 8+|+7+1+t.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . v+A+8+|+7+x+7+|+8+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + @ # # ~.K+. . . + n.. . . % + . . . { - [ . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . * = - ; > , . . . + @ # F.L y a . . + @ # # # $ . . . . . 1 L+. . . . 5.M+@ # k u { . . . . . . . . . . 5.M+@ # k u { . . . . . # & . . . . . . . . + % . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 8+t.. . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 8+|+7+x+7+|+. . . + @ # F.L y a . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . # & . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . 5.N+v.. . + O+P+. . % + . . . Y Q+w 3.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ) ! ~ { . ] ^ * . . + % . . N *.R+. . + % . . . . . . . . { S+N . . . s ; T+F . s U+*.. . . . . . . . . s ; T+F . s U+*.. . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 8+|+7+1+t.. . + % . . N *.R+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . ] U+. . + V+g ' . % + . . . 2 U G : . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . _ : . . . . < [ . . + % . . . s O+. . + % . . . . . . . . M.2 . . . . O V . . . . 1.0 . . . . . . . . . O V . . . . 1.0 . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 8+t.. . . + % . . . s O+. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . & ..r . + W+~.v . % + . . }.X+e.Y+N+N . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ; , . . . . s t r . + % . . { Z+B.. . + % . . . . . . . . d.`+. . . . !  @. . . . . . . . + f # g h . . !  @. . . . . . . . . . . + % . .@+@@ @@k.. . + % a b # c d . . + v k w x . 9.e # f + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . + v k w x . . + v k w x . o p @ - q r . + 4 ; i . . . + 4 ; u b # c d . . + v k w x . s y # z + % o p @ - q r . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . { Z+B.. o p @ - q r #@#.. . l.H o p @ - q r . + 4 ; i .@+@@ @@k.. o p @ - q r . . . . + f # g h . . + 4 ; i + % . + j k l m $ c n . o p @ - q r . + 4 ; i . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . r + b.. + % $@t $@% + . . P [ . . +@D.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..+.. . . . @.+ 3 . + @ # # w ` %@. . + @ # # # D . . . . E.&@. . . . *@3.. . =@# # # . . + V . W e . . *@3.. . =@# # # . . . . . + % . 6.-@. `+: . . + % O P . Q + R . + #.. $...%.` T . U + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . + #.. $...%.. + #.. $...%.H Z . 8 `  .. + I . { . . . + I . 9 P . Q + R . + #.. $...%.&.*.. U + % H Z . 8 `  .. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + @ # # w ` %@. H Z . 8 `  .X ;@. . B.5 H Z . 8 `  .. + I . { 6.-@. `+: . H Z . 8 `  .. . . . + V . W e . . + I . { + % . + : . X + Y s + R H Z . 8 `  .. + I . { . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . &@O+r . + % . v ~.W++ . . >@X+# # ,.#@. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ; 1.. . . . 2.7 r . + % N $.! m . . . + % . . . . . . . r.,@. . . . . c L+. . . . % + . . + ) . R + 8 . c L+. . . . % + . . . . . + % . :.# G '@%.. . + % r _.u :.+ <.. + 3.. N + <.[.s . }.+ % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . + 3.. N + <.. + 3.. N + <.[.@ # # |.H . + !.. . . . . + !.. r _.u :.+ <.. + 3.. N + <.[.s . / + % [.@ # # |.H . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % N $.! m . . [.@ # # |.H N ( $@1.P.N [.@ # # |.H . + !.. . :.# G '@%.. [.@ # # |.H . . . . + ) . R + 8 . + !.. . + % . + s . % + F . + % [.@ # # |.H . + !.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . )@E.. . + % . ' g !@+ . !.X+a.. . s ~@& . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . c.d.. . . .  .` . . + % . . {@]@^@. . + % . . . . . . . /@(@. . . . . _ =@. . . . % + . . + s . / + a.. _ =@. . . . % + . . . . . + % . . _@l.:@y . . + % ` :.] 9.+ % . + % . . + % [.0.. 8 + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . + % . . + % . + % . . + % ).9.. . b.8 . + % . . . . . + % . ` :.] 9.+ % . + % . . + % ).s . R + % ).9.. . b.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . {@]@^@. ).9.. . b.8 . u <@^.d.. ).9.. . b.8 . + % . . . _@l.:@y . ).9.. . b.8 . . . . + s . / + a.. + % . . + % . + % . % + . . + % ).9.. . b.8 . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . _@N+_.. . + % . . L+]@+ . u H . . . . y m . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . !.- u.N . v.^ a . . + % . . . f @@{ . + % . . . . . . . f X . . . . . s ! l.a.. ) p [@. . + u . 0 p.. . s ! l.a.. ) p [@. . . . . + % / |. @. A.w.. . + % t k.. l.+ 3.. + % . . + % m.n.. o.+ % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . + % . . + % . + % . . + % q.f { r.s.U . + % . . . . . + % . t k.. l.+ 3.. + % . . + % = u . z.+ <.q.f { r.s.U . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . f @@{ q.f { r.s.U . D # # 2.. q.f { r.s.U . + % . / |. @. A.w.. q.f { r.s.U . . . . + u . 0 p.. . + % . . + % . + % . % + . . + % q.f { r.s.U . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + @ # @ }@u.. . . + % . . . n.+ . ~.6 . . . . |@g . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 5.= ! - > , . . . + % . . . 3 w z.. + % . . . . . . . g ' . . . . . . s !@# - P.<@s . . + _ F.e G.. . . s !@# - P.<@s . . . . . + % . `+F.@ ~.1@. . + % ].c C.d ~.D.. + % . . + % 2.E.# O + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . + % . . + % . + % . . + % % m.! L n { . + % . . . . . + % . ].c C.d ~.D.. + % . . + % s e @ :.+ r.% m.! L n { . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . 3 w z.% m.! L n { . . q.[ . . % m.! L n { . + % . . `+F.@ ~.1@. % m.! L n { . . . . + _ F.e G.. . + % . . + % . + % . % + . . + % % m.! L n { . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . L.M.. 5 @ . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . l.P.F.Q.R.. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . + % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . "};

end_of_pixmap
}
