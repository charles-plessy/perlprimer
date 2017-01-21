PerlPrimer
----------

Copyright © 2003-2005, Owen Marshall (owenjm@users.sourceforge.net)


Contents
--------

1. Introduction
2. Installation / Using the program
3. Links to external programs
4. Suggestions/patches/bugs/support


1. Introduction
---------------

Perlprimer is a GUI application that designs primers for standard PCR, bisulphite PCR, Real-time PCR (QPCR) and sequencing.

Perlprimer's current features include the following:

* Calculation of possible primer-dimers
* Retrieval of genomic or cdna sequences from Ensembl (including both sequences automatically for QPCR)
* Ability to BLAST search primers using either the NCBI server or a local server
* Results can be saved or optionally exported in a tab-delimited format that is compatible with most spreadsheet applications.
* ORF and CpG island detection algorithms
* Ability to add cloning sequences to primers, automatically adjusted to be in-frame
* QPCR primer design without manual intron-exon boundary entry

Perlprimer calculates primer melting temperature using J. SantaLucia's extensive nearest-neighbour thermodynamic parameters.  To adjust for the salt conditions of the PCR, PerlPrimer uses the empirical formula derived by von Ahsen, et al. (2001) and allows the user to specify the concentration of Mg2+, dNTPs and primers, or use default, standard, PCR conditions. The result is a highly accurate prediction of primer melting temperature, giving rise to a maximum yeild of product when amplified.

Perlprimer is written in Perl and requires Perl/Tk.  In addition, for QPCR functionality perlprimer requires the open-source Spidey executable from NCBI.  The program is designed to be cross-platform and has been tested on both Microsoft Windows and GNU/Linux-based operating systems.  Users have also had success using the program under Mac OS X.

Restriction enzyme data is provided by the REBASE project (http://rebase.neb.com/)


2. Installation / Using the program
-----------------------------------

PerlPrimer requires Perl and Perl/Tk.  These should be provided by your distribution if using a UNIX-based operating system; if using Windows I recommend using ActiveState's distribution of Perl 5.8, which includes Perl/Tk, freely available from

	http://www.activestate.com/Products/ActivePerl/

After installing ActivePerl under Windows, the file "perlprimer.pl" should run if double-clicked.

PerlPrimer can also be started from the commandline, and this would be the normal mode of operation for UNIX-based OSes (including Mac OS X):

	$ perl perlprimer.pl [file to open]

(the file to open parameter is optional)


Users under UNIX-based OSes may need to install Perl/Tk manually if it is not included in your distribution.  Perl/Tk can be downloaded from 

	http://search.cpan.org/~ni-s/Tk/

and can be installed by following the instructions included in the archive.  Users may also need the Perl modules HTTP::Request and LWP::UserAgent for BLAST searching and gene retieval from Ensembl; these are part of libwww-perl-5.76 which can be found at

	http://search.cpan.org/~gaas/libwww-perl-5.76/

In addition to installing Perl/Tk (and libwww-perl if required), Mac OS X users will generally require OS 10.3 (Panther) or later and an X-server.  (Please note that as I do not have access to a Mac OS X system I cannot guarantee compatibility, although I will try to fix any issues that are reported).


Using PerlPrimer should be fairly self-explanatory, with extensive "balloon help" (turned on by default) and a separate help window detailing the operation of the graphical display of the DNA sequence, selection ranges and primers. A tutorial is provided with this distribution, and is also available at

	http://perlprimer.sourceforge.net/tutorial.html

which covers the most commonly used features of the program.


3. Links to external programs
-----------------------------

For Real-time PCR (QPCR) functionality, PerlPrimer uses the program Spidey, released freely by the NCBI, to calculate intron/exon boundaries.  This may be obtained from

	http://www.ncbi.nlm.nih.gov/spidey/

The program by default expects the spidey executable to be in your home directory on UNIX-based systems, or at C:\Spidey.exe for Windows.  The location may be changed in the preferences window.


4.  Suggestions/patches/bugs/support
------------------------------------

The latest release of PerlPrimer can be obtained from
	
	http://perlprimer.sourceforge.net/

Please check the lastest version before reporting bugs.

Bugs, feature suggestions and support requests can be placed on the web forums provided at the above address, or alternatively by emailing me directly.