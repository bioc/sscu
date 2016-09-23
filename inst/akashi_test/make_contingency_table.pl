#!/usr/bin/perl -w

use strict;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Perl;
use Getopt::Std;
use Data::Dumper;

# usage example: perl make_contingency_table.pl -t 25 -r 35 -f ~/Data/codon_usage/bee_endosymbionts/akashi_test/lactos_all -o CAC,GAC,AAC,ATC,CGT,GGT,GTT
# -t the number of the target sequence in the alignment
# -r the number of the reference sequence in the alignment
# -f the folder/directory for the alignment files
# -o the optimal codons. The optimal codons could be either identified previously, or identified by functions in the sscu package.

getopts('t:r:f:o:');
our($opt_t, $opt_r, $opt_f, $opt_o);

# Gvag 19, Bbif 7, num1 is the target sequence, num2 is the reference sequence
my $aln_num1 = $opt_t;  
my $aln_num2 = $opt_r;

my $dir = $opt_f;
my @files = <$dir/*>;

my $optimal_codon = $opt_o;
my @optimal_codons = split(/,/, $optimal_codon);

#print "@optimal_codons\n";

#my @files = <~/Data/codon_usage/bee_endosymbionts/akashi_test/bifidos_ribo/*.mafft>;

my $contingency_data;

foreach my $file (@files){
    #if($file eq '/home/yu/Data/codon_usage/bee_endosymbionts/akashi_test/bifidos_ribo/0082.fas.ffn.mafft'){
    my $str = Bio::AlignIO->new(-file => $file,
				-format => 'fasta');
    my $aln = $str->next_aln();
    my $obj1 = $aln->get_seq_by_pos($aln_num1);
    my $seq1 = $obj1->seq;
    
    my $obj2 = $aln->get_seq_by_pos($aln_num2);
    my $seq2 = $obj2->seq;
    
    my $length = $aln->length();
    
    my %conserved_codons  = (   
    'TCA' => 0,'TCC' => 0,'TCG' => 0,'TCT' => 0,'TTC' => 0,'TTT' => 0,'TTA' => 0,'TTG' => 0,'TAC' => 0,'TAT' => 0,'TAA' => 0,'TAG' => 0,'TGC' => 0,'TGT' => 0,'TGA' => 0,'TGG' => 0,'CTA' => 0,'CTC' => 0,'CTG' => 0,'CTT' => 0,'CCA' => 0,'CCC' => 0,'CCG' => 0,'CCT' => 0,'CAC' => 0,'CAT' => 0,'CAA' => 0,'CAG' => 0,'CGA' => 0,'CGC' => 0,'CGG' => 0,'CGT' => 0,'ATA' => 0,'ATC' => 0,'ATT' => 0,'ATG' => 0,'ACA' => 0,'ACC' => 0,'ACG' => 0,'ACT' => 0,'AAC' => 0,'AAT' => 0,'AAA' => 0,'AAG' => 0,'AGC' => 0,'AGT' => 0,'AGA' => 0,'AGG' => 0,'GTA' => 0,'GTC' => 0,'GTG' => 0,'GTT' => 0,'GCA' => 0,'GCC' => 0,'GCG' => 0,'GCT' => 0,'GAC' => 0,'GAT' => 0,'GAA' => 0,'GAG' => 0,'GGA' => 0,'GGC' => 0,'GGG' => 0,'GGT' => 0    
	);
    my %variable_codons = (   
    'TCA' => 0,'TCC' => 0,'TCG' => 0,'TCT' => 0,'TTC' => 0,'TTT' => 0,'TTA' => 0,'TTG' => 0,'TAC' => 0,'TAT' => 0,'TAA' => 0,'TAG' => 0,'TGC' => 0,'TGT' => 0,'TGA' => 0,'TGG' => 0,'CTA' => 0,'CTC' => 0,'CTG' => 0,'CTT' => 0,'CCA' => 0,'CCC' => 0,'CCG' => 0,'CCT' => 0,'CAC' => 0,'CAT' => 0,'CAA' => 0,'CAG' => 0,'CGA' => 0,'CGC' => 0,'CGG' => 0,'CGT' => 0,'ATA' => 0,'ATC' => 0,'ATT' => 0,'ATG' => 0,'ACA' => 0,'ACC' => 0,'ACG' => 0,'ACT' => 0,'AAC' => 0,'AAT' => 0,'AAA' => 0,'AAG' => 0,'AGC' => 0,'AGT' => 0,'AGA' => 0,'AGG' => 0,'GTA' => 0,'GTC' => 0,'GTG' => 0,'GTT' => 0,'GCA' => 0,'GCC' => 0,'GCG' => 0,'GCT' => 0,'GAC' => 0,'GAT' => 0,'GAA' => 0,'GAG' => 0,'GGA' => 0,'GGC' => 0,'GGG' => 0,'GGT' => 0 
	);

    for(my $i=0; $i<$length; $i=$i+3){
	my $codon1 = substr $seq1, $i, 3;
	my $codon2 = substr $seq2, $i, 3;
	if($codon1 ne '---' and $codon2 ne '---'){
	    my $aa1 = translate_as_string($codon1);
	    my $aa2 = translate_as_string($codon2);
	    if($aa1 eq $aa2){
		#print "$aa1";

		$conserved_codons{$codon1}++;
	    }
	    else{
		#print "$aa1";

		$variable_codons{$codon1}++;
	    }
	    
	}
    }
   
    #print Dumper(\%conserved_codons);

    my @codons_set_ref = (
	['TCA','TCC','TCG','TCT','AGC','AGT'], #ser
	['CGA','CGT','AGA','CGG','CGC','AGG'], #arg
	['CTA','CTT','TTA','CTG','CTC','TTG'], #LEU
	['ATA','ATT','ATC'], #ILE
	['GTA','GTT','GTG','GTC'], #VAL
	['CCA','CCT','CCG','CCC'], #PRO
	['ACA','ACT','ACG','ACC'], #THR
	['GCA','GCT','GCG','GCC'], #ALA
	['GGA','GGT','GGG','GGC'], #GLY
	['AAT','AAC'], #ASN
	['GAT','GAC'], #ASP
	['TGT','TGC'], #CYS
	['CAA','CAG'], #GLN
	['GAA','GAG'], #GLU
	['CAT','CAC'], #HIS
	['AAA','AAG'], #LYS
	['TTT','TTC'], #PHE
	['TAT','TAC'] #TYR
	);


    foreach my $codon_set_ref (@codons_set_ref){
	$contingency_data .= contingency_table($codon_set_ref,\@optimal_codons,\%conserved_codons,\%variable_codons);
    }
    chop $contingency_data;
    $contingency_data .= "\n";
    
    
    #}
}
chop $contingency_data;
print "$contingency_data\n";


sub contingency_table
{
    my $codons_ref = shift;
    my $optimal_codons_ref = shift;
    my $conserved_codons_ref = shift;
    my $variable_codons_ref = shift;

    my @optimal_codons = @$optimal_codons_ref;
    my $con_op = 0; 
    my $con_nop = 0;
    my $var_op = 0;
    my $var_nop = 0;
    foreach my $codon (@$codons_ref){
	if($codon ~~ @optimal_codons){
	    $con_op += $conserved_codons_ref->{$codon};
	    $var_op += $variable_codons_ref->{$codon};
	}
	else{
	    $con_nop += $conserved_codons_ref->{$codon};
	    $var_nop += $variable_codons_ref->{$codon};
	}
    }
    if($con_op+$con_nop+$var_op+$var_nop>1){
	return("$con_op,$var_op,$con_nop,$var_nop\n");
    }
}
