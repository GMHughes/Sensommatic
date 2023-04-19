#!/usr/bin/perl
use strict;
use warnings;

###########################################################################
# OPTIONS AND PARAMETERS: 
###########################################################################

############
# ESSENTIAL:
############

my $receptorfile="MammalReceptors.fa"; #Change this to suit the class of vertebrates you are interest in
my $augustus_species = "human"; #change this to suit your query species. See README file for suitable augustus species.

############
# OPTIONAL:
############

#Options
my $block_region = 750; #Mask 750 nt at both ends of prediction, reduces occurance of single loci being reported as two predictions
my $threads=20; #number of threads for blast
my $hmm_filter = "yes"; #if "yes", use hmmsearch to ensure prediction contains 7TM domain. Minimises off target hits. May reduce recovery slightly.
my $keep_same_seq= "yes"; #if yes, you may get multiple copies of the same gene if they occur at different positions (assuming 100% identical sequence)
my $pseudogenome_esl = "yes"; #write out contigs with hmmer esl-sfetch, requires hmmer esl suite. Improves speed, does not impact or alter predictions
my $repeat_preds = 0; #Repeat predictions at difficult regions #0: Fast, #1: Slow (may pick up some extra hits)
my $remove_intermediate= "yes"; #if yes, remove all intermediate output files
my $print_datestring = "no"; #if datesting is "yes", prints start and end date/time
my $minlength=300; #Minimum length for sequence to be considered a gene in the output file. Predictions lower than this length will not be included in the output.

#Augustus pararmeters
my $multi_exon_range = 50000; #blast hit +/- this range is fed into augustus - multiexon receptors
my $single_exon_range = 2500; #blast hit +/- this range is fed into augustus - single exon receptors

#BLAST parameters: Thresholds for a hit to be considered
my $cover_threshold = 300; #minimum coverage (number of nucleotides) required for blast hits to be considered
my $identity_threshold = 70; #minimum identity required for blast hits to be considered
my $evalue_inclusion = "1e-30"; #Blast hits with evalues lower than this value will be considered
my $qcov =30; #percentage cover required for single exon blast hits to be considered

#Pseudogene length parameters: Predictions with lengths lower than these thresholds will be labelled as 'pseudogene_short'
my $single_exon = 860; #single exon pseudogene length
my $opsin_multiexon = 1000; #opsin (multiexon) pseudogene length
my $taste_multiexon = 1000; #taste (TAS1R - multiexon) pseudogene length
my $vomeronasl_multiexon = 2000; #vomeronasal (V2R -multiexon) pseudogene length

#Sensommatic Parameters (only applies when augustus fails to predict a single exon gene):
my $max_length_external = 1300; #external extension limit; if missing start or stop, this is the length that the sequence can expand to when searching for start and stop codons.
my $min_length_internal  = 500; #internal extension limit; this is the length that the sequence can contract to when searching internally for start and stop codons.
my $type_length_switch = 0; #1: ON, overwrites $max_length_external and uses family specific extension lengths, #0: OFF, uses $max_length_external for all receptors regardless of query type.
my $OLFR_max_length = 1300; #external extension limit: Olfactory Receptors
my $TAAR_max_length = 1300; #external extension limit: Trace amine associated receptors
my $TASR_max_length = 1300; #external extension limit: Taste receptors 
my $VOMR_max_length = 1300; #external extension limit: Vomeronasal receptors
my $percent = 0.25; #percent sequence is contracted prior to external search

#Blastback parameters: Prediction is blasted back to reference receptor in a validation step to avoid over-extension or over-trimming.
my $blastback_verification = 1; #1:ON, blastback verification step on; #0:OFF, blastback verification step off.
my $percentage_cover_threshold = 0.6; #0.6; #number of nucleotides aligned with reference receptor divided by prediction length (cover/prediction length) 
my $percentage_identity_threshold = 0.5; #0.5; #number of identical nucleotides between prediction and reference receptor divided by prediction length (identity/prediction length)
my $evalue_threshold ="1e-30";#"1e-30"; #evalue threshold for blastback step
my $include_score = 0; #Include blastback scores (percentage cover, percentage identity and evalue) in fasta headers. #1: ON, include scores; #0: OFF, don't include scores

#Annotations: Functional vs Pseudogene
my $annotation_short = "pseudogene_short"; #If prediction is shorter than $pseudogene_length, gene will be annotated as pseudogene regardless of conditions (1-6).
my $annotation_1 = "functional"; #START codon && no in frame stop codons...........................: ATG -----------
my $annotation_2 = "functional"; #no START codon && no stop codons in any frame....................: ---------------
my $annotation_3 = "functional"; #no START codon && no in frame stop codons........................: ---------------
my $annotation_4 = "pseudogene_nonsense"; #START codon && in frame stop codon.............: ATG------TGA---
my $annotation_5 = "pseudogene_nonsense"; #no START codon && stop codons in all frames....: ---------TGA---
my $annotation_6 = "pseudogene_nonsense"; #no START codon && in frame stop codon..........: ---------TGA---

#Truncated option: for predictions with no end stop codon
my $truncated = 1; #1: adds truncated annotation; 0: truncated annotations off
my $annotation_truncated = "truncated"; # gets added to above annotations if no end STOP codon and if $truncated option equals 1.


###########################################################################
# SENSOMMATIC:                                         
###########################################################################

#Declare variables
my $datestring = localtime();##Get local time
if ($print_datestring eq "yes"){
    print "Starting at: ".$datestring."\n";
}
my @ranset = ('0' ..'9', 'A' .. 'F');#create a random set of variables
my $randomstr = join '' => map $ranset[rand @ranset], 1 .. 8;
my $genomefile=$ARGV[0];#take in genome file
chomp $genomefile;#remove whitespace
my $backup_genomefile=$genomefile;#create a copy of name of genomefile for end of script
my $refDB=$randomstr."_refDB";#create a series of output files using random string:
my $prot_frame_db =$randomstr."protein_frame_DB"; 
my $target_nuc=$randomstr."_target_nuc.fa"; 
my $tmp_X_file = $randomstr."_tmp_X_file.fa";
my $fullrun_out=$randomstr."_fullrun_out"; 
my $frame_out= $randomstr."_protein_frame_out";
my $tmp_seq=$randomstr."_tmp_seq.fa";
my $blat_output=$randomstr."_blat_output";
my $renamed_seqfile=$randomstr."_Reference_receptors_unique_names.fa";
my $prot_frames = $randomstr."_protein_translated_frames.fa"; 
my $query_prot = $randomstr."_query_receptor_prot.fa";
my $prot_ref = $randomstr."_tblastn_reference_receptor.fa";
my $tblastn_out = $randomstr."_tblastn_out";
my $Tblastn_set = $randomstr."_tblastn_database.fa";
my $Tblastn_db = $randomstr."_tblastn.db";
my $ref_file_db = $randomstr."_reference_file_database.db";
my $multiout = $randomstr."_multiout.fa";
my $multiblatout = $randomstr."_multiblatout.fa";
my $max_length; 
my @coordinates=();#array stores coordinates.
my @final_unique_coords=();#array stores final coordinates to prevent overlapping predictions.
my $dna=""; 
my $contig_name="";
my $protseq="";
my $first=-1;
my @exonfiles;
my $plusminus = "";
my $multiexon = 0;
my $count = 0;
my $pseudogene_length = "";
my $output_format = "6 qseqid qlen sseqid slen qstart qend sstart send evalue length nident qcovhsp pident";
my $profilehmms = "ProfileHmms.hmm";

# Split ProfileHmms file to individual profile hmms
my $profile_hmms = "";
{
    local $/;
    open HMMS, $profilehmms;
    $profile_hmms = <HMMS>;
}
close HMMS;
my @hmms = split(/SPLIT/, $profile_hmms);
my $hmm1 = $hmms[0];
my $hmm2 = $hmms[1];
my $hmm3 = $hmms[2];

my $classA_hmm = $randomstr."_GPCRClassA.hmm";
my $classC_hmm = $randomstr."_GPCRClassC.hmm";
my $v2r_hmm = $randomstr."_vn2r.hmm";

open(HMM, ">>$classA_hmm");
print HMM $hmm1;
close HMM;

open(HMM, ">>$classC_hmm");
print HMM $hmm2;
close HMM;

open(HMM, ">>$v2r_hmm");
print HMM $hmm3;
close HMM;

#Check option
my $pseudogenome_perl = "no"; #write out contigs with hits to pseudogenome using perl, does not require hmmer esl suite, slower
if($pseudogenome_esl ne "yes" || $pseudogenome_esl ne "Yes"){
    $pseudogenome_perl = "yes";
    $pseudogenome_esl = "no";
}
if($pseudogenome_perl eq "yes" && $pseudogenome_esl eq"yes"){ #can't be both, reset to perl if both are yes
    $pseudogenome_esl = "no";
}
my $reference_file =  $receptorfile;
if($blastback_verification == 0){ #if blastback is off, $include_score must be off
    $include_score = 0;
}

#Genome information
my @genome_info = ();
push(@genome_info, $genomefile);
push(@genome_info, $receptorfile);
push(@genome_info, $randomstr);
push(@genome_info, $threads);
push(@genome_info, $pseudogenome_perl);
push(@genome_info, $pseudogenome_esl); 
my $pseudogenome = "";

#Output pseudogenome - only contigs with hits in output file (speed)
if($pseudogenome_perl eq "yes" || $pseudogenome_esl eq "yes"){
    $pseudogenome =&pseudogenome(@genome_info); 
    push (@genome_info, $pseudogenome);
}
my @names =&get_contigs(@genome_info); #get contig names
my @seqs =&parse_fasta($receptorfile); #parse reference file
my @receptors =&rename_refs(@seqs); #assign unique names to reference file
my %receptor_key =&gene_key(@receptors); #convert to key

#Split single exon and multi exon reference files
my $multiexon_file = $randomstr."_multiexon_refs.fa";
my $singleexon_file = $randomstr."_singleexon_refs.fa";
open(ME, ">$multiexon_file");
open(SE, ">$singleexon_file");
foreach my $rec(@receptors){
    if($rec =~ m/OPN/i ||
       $rec =~ m/VN2/i ||
       $rec =~ m/Vom2/i ||
       $rec =~ m/Vmn2/i ||
       $rec =~ m/V2/i ||
       $rec =~ m/V2R/i ||
       $rec =~ m/TAS1R/i ||
       $rec =~ m/T1R/i ||
       $rec =~ m/RHO/i || 
       $rec =~ m/Rho/i ||
       $rec =~ m/Rh[0-9]/i){          
	print ME $rec;
    }
    else{
	print SE $rec;
    }
}
close ME;
close SE;

#Open genome or pseudogenome
if ($pseudogenome_perl eq "yes" || $pseudogenome_esl eq "yes"){
    open(GENOME, "$pseudogenome"); #quicker version
}
else{
    open(GENOME, "$genomefile"); #slowest version
}

#Generate output file for hits filtered out with hmmer
my $suffix_hmm = "";
if($genomefile =~ m/(GC[A-Z]?\_[^\_]+)/){
    $suffix_hmm = $1;
}
my $hmm_filter_outfile = "";
if($suffix_hmm){
    $hmm_filter_outfile = "Hmm_filtered_out_".$suffix_hmm.".fa";
}
else{
    $hmm_filter_outfile = "Hmm_filtered_out.fa";
}
if($hmm_filter eq "yes" || $hmm_filter eq "Yes"){
    open(HMMOUT, ">>$hmm_filter_outfile");
}


#Sensommatic:
`makeblastdb -in $receptorfile -dbtype="nucl" -out $ref_file_db`; #make blastdb from the target contig
while(<GENOME>){
    if($_=~m/\>([\S\s]+)/){
	$contig_name=$1; #store contig name
	my @final_unique_coords = (); #Resets at each contig, coordinates considered overlapping only if on same contig
	if(length($dna)==0){#if no contig sequence has been held in memory yet
	    $first++;#take note of it 
	}
       	else{
	    $contig_name=$names[$first];#otherwise take name of contig #CONTIG NAME
	    $first++;
	    chomp $contig_name;
	    $contig_name=~s/\>//g;
	    $dna=~s/[\n\r]//g; #Remove line breaks from contig data. If lines are <80 format, fails without this statement.
	    open(OUT, ">$target_nuc");
	    print OUT ">".$contig_name."\n".$dna; #Print current target contig to tmp file
	    close OUT;
	    my $opsin_type = 0;
	    my $fullrun_out1 = $fullrun_out."_single_exon.out";
	    my $fullrun_out2 = $fullrun_out."_multi_exon.out";
	    `makeblastdb -in $target_nuc -dbtype="nucl" -out $refDB`; #make blastdb from the target contig
	    `blastn -task 'megablast' -db $refDB -query $singleexon_file -out $fullrun_out1 -num_threads $threads -outfmt \"$output_format\" -qcov_hsp_perc $qcov`;
	    `blastn -task 'megablast' -db $refDB -query $multiexon_file -out $fullrun_out2 -num_threads $threads -outfmt \"$output_format\" -qcov_hsp_perc 15` ;
	    `cat $fullrun_out2 >> $fullrun_out1`;
	    my @input_firstsweep = ();
	    push (@input_firstsweep, $fullrun_out1);
	    push(@input_firstsweep, $cover_threshold);
	    push(@input_firstsweep, $qcov);
	    my @touse=();#array to hold sensory genes that map to contig of interest
	    @touse=&firstsweep(@input_firstsweep);
	    foreach my $ordered_gene(@touse){ #foreach receptor that had a hit (ordered by increasing e-value, and decreasing identity to ensure most signifigant reference receptors are used)
		my $gene_seq = $receptor_key{$ordered_gene};
		my $gene_name = $ordered_gene;
		my $sense_gene = ">".$gene_name."\n".$gene_seq;
		my $prediction_included = 0;
		if($type_length_switch ==1){ #if receptor specific length option is on, check what type of receptor your reference is
		    if($gene_name =~ m/OR/ || $gene_name =~ m/Olr/i || $gene_name =~ m/Olfr/i){ #olfactory receptors
			$max_length = $OLFR_max_length;
		    }
		    elsif($gene_name =~m/TAAR/i){ #trace amine accosiated receptors
			$max_length = $TAAR_max_length;
		    }
		    elsif($gene_name =~ m/TAS/ || $gene_name =~ m/T2R/){ #taste receptors
			$max_length = $TASR_max_length;
		    }
		    elsif($gene_name =~ m/VN/ || $gene_name =~ m/Vom/ || $gene_name =~ m/Vmn/ || $gene_name =~ m/V1R/){ #vomeronasal receptors
			$max_length = $VOMR_max_length;
		    }
		    else{
			$max_length = $max_length_external; #if doesnt match any of the above annotations, will just use the $max_length_external
		    }
		}
		else{ #if $type_length_switch is off, use single max_length_external
		    $max_length = $max_length_external;
		}
		if($gene_name =~ m/OPN/i ||
		   $gene_name =~ m/VN2/i || 
		   $gene_name =~ m/Vom2/i || 
		   $gene_name =~ m/Vmn2/i || 
		   $gene_name =~ m/V2/i || 
		   $gene_name =~ m/TAS1R/i || 
		   $gene_name =~ m/T1R/i ||
		   $gene_name =~ m/OPN1S/i ||
		   $gene_name =~ m/OPN1M/i ||
		   $gene_name =~ m/OPN1L/i ||
		   $gene_name =~ m/OPN2/i ||
		   $gene_name =~ m/RHO/i ||
		   $gene_name =~ m/Rhol/i ||
		   $gene_name =~ m/rh[0-9]/i){
		    $multiexon = 1;
		    $plusminus = $multi_exon_range;
		    if($gene_name =~ m/OPN/i ||
		       $gene_name =~ m/OPN1S/i ||
		       $gene_name =~ m/OPN1M/i ||
		       $gene_name =~ m/OPN1L/i ||
		       $gene_name =~ m/OPN2/i ||
		       $gene_name =~ m/RHO/i ||
		       $gene_name =~ m/Rho/i || 
		       $gene_name =~ m/Rh[0-9]/i){
			$pseudogene_length = $opsin_multiexon;
			$opsin_type = 1;
		    }
		    elsif($gene_name =~ m/TAS1R/i ||
			  $gene_name =~ m/T1R/i){
			$pseudogene_length =$taste_multiexon;
			$opsin_type = 1;
		    }
		    elsif($gene_name =~ m/VN2/i ||
			  $gene_name =~ m/Vom2/i ||
			  $gene_name =~ m/Vmn2/i ||
			  $gene_name =~ m/V2/i){
			$pseudogene_length =$vomeronasl_multiexon;
		    }			      
		}
		else{
		    $plusminus = $single_exon_range;
		    $multiexon = 0;
		    $pseudogene_length =  $single_exon;
		}
		open(OUTGENE, ">$tmp_seq"); #write out reference receptor to $tmp_seq
		print OUTGENE $sense_gene;
		close OUTGENE;
		`blastn -task 'megablast' -db $refDB -query $tmp_seq -out $blat_output -num_threads $threads -max_hsps 10`;#blast the reference receptor to the target contig
		my @details=();
		my @input_blat = ();
		push(@input_blat, $blat_output);
		push(@input_blat, $max_length);
		push(@input_blat, $min_length_internal);
		push(@input_blat, $plusminus);
		push(@input_blat, $pseudogene_length);
		push(@input_blat, $percent);
		push(@input_blat, $cover_threshold);
		push(@input_blat, $identity_threshold);
		push(@input_blat, $evalue_inclusion);
		push(@input_blat, $max_length_external);
		push(@input_blat, $blastback_verification);
		push(@input_blat, $percentage_cover_threshold);
		push(@input_blat, $percentage_identity_threshold);
		push(@input_blat, $evalue_threshold);
		push(@input_blat, $multiexon);
		push(@input_blat, $minlength);
		push(@input_blat, $ref_file_db);
		push(@input_blat, $multiout);
		push(@input_blat,  $multiblatout);
		push(@input_blat, $opsin_type);
		push(@input_blat, $tmp_seq);
		push(@input_blat, $classA_hmm);
		push(@input_blat, $hmm_filter);
		push(@input_blat, $classC_hmm);
		push(@input_blat, $v2r_hmm);
		@details =&blatout(@input_blat);
		my @already=();
		my $infoheader;
		my $infoheader_original;
		if (@details){ #if details is populated
		    my $info = $details[0]; #prediction header and sequence
		    if($info ne ""){
			my $coord1 = $details[1]; #blank coord 1
			my $coord2 = $details[2]; #blank coord 2
			my $maximum_contig_bound = $details[3];
			my $homolog = $details[4];
			if($info ne "NOPREDICTION" && $info ne ""){
			    if($info=~m/(\>Seq[\S]+)\n([\S]+)/){ 
				$infoheader=$1; my $infoseq=$2; #store prediction header and sequence
				$infoheader_original = $infoheader;
				$infoheader_original =~ s/Seq//g;
				$info=~s/Seq//g;
				my $hmm_remove_status = 0;
				if($infoheader =~ m/hmm\_remove/){
				    $hmm_remove_status = 1;
				}
				if($infoheader =~ m/([^\|]+?\|)/){
				    $infoheader = $1; #gets rid of scores from infoheader for coordinate check - important step.
				    $infoheader =~ s/\|//g;
				}
				my $overlap="N";
				#Take coordinates, if not already stored in @final_unique_coords, keep overlap as N, add coords. If coordinates overlap, set overlap as Y.
				if($coord2<$coord1){ #if reverse, flip so coordinates are ordered small to big.
				    my $tmpc=$coord1;
				    $coord1=$coord2;
				    $coord2=$tmpc;
				}
				my $coordspan=$coord1."\.\.\.".$coord2;
				my $coordlen=$coord2-$coord1; #prediction length 
				my $maxcoord="";my $mincoord="";
				foreach my $uniquecoord(@final_unique_coords){ #loop through coordinates already stored
				    if($uniquecoord=~m/([0-9]+)\.\.\.([0-9]+)/){
					my $refcoord1=$1;my $refcoord2=$2; #store coordinates for each prediction already in output 
					if($refcoord2<$refcoord1){
					    my $tmpc=$refcoord1; #if reverse, flip so coordinates are ordered small to big
					    $refcoord1=$refcoord2;
					    $refcoord2=$tmpc;
					}
					my $refcoordlen=$refcoord2-$refcoord1; #length of prediction already in output
					if($refcoord2>$coord2){ #find maximum coordinate - comparing new prediction to prediction already in output
					    $maxcoord=$refcoord2;
					}
					else{
					    $maxcoord=$coord2;
					}
					if($refcoord1<$coord1){ #find minimum coordinate - comparing new prediction to prediction already in output
					    $mincoord=$refcoord1;
					}
					else{
					    $mincoord=$coord1;
					}
					my $maxspan=$maxcoord-$mincoord; #coordinate span: maximum coordinate - minimum coordinate
					my $totallength=$coordlen+$refcoordlen; #total length of new prediction and prediction aldready in output
					if($maxspan<$totallength){ #if span is greater than total length, the coordinates overlap
					    $overlap="Y"; #set overlap to YES
					}
				    }
				}
				if($overlap ne "Y"){ #if coordinates do not overlap, store coordinates
				    push @final_unique_coords,$coordspan; #only push nonoverlapping coordinates
				}
				if($overlap eq "Y" || $info~~@already || $infoheader~~@already || $keep_same_seq eq "no" && $infoseq~~@already){ 
				    #if coordinates overlap, or exact sequence is already stored: exclude prediction from output
				    #if $keep_same_seq is "no", and sequences are identical, despite unique locations (e.g on different contigs), exclude duplicated sequence
				}
				else{ #otherwise - include prediction in output
				    my $infoseq_wrap = $infoseq;
				    $infoseq_wrap =~ s/.{80}\K/\n/g;
				    $infoseq_wrap = uc($infoseq_wrap);
				    if($hmm_remove_status == 0){
					print $infoheader_original."\n".$infoseq_wrap."\n\n";
				    }
				    else{
					print HMMOUT $infoheader_original."\n".$infoseq_wrap."\n\n";
				    }
				    push @already,$info;
				    push @already,$infoheader;
				    $prediction_included = 1;
				    if($keep_same_seq eq "no"){
					push @already,$infoseq;
				    }
				}
			    }
			}
			#block out prediction region with "N"s for next blast round to avoid repeating blast hits
			if($multiexon ==0 && $prediction_included == 1){
			    $coord1 -=$block_region;
			    if($coord1 < 1){
				$coord1 = 1;
			    }
			    $coord2 +=$block_region;
			    if($coord2 > $maximum_contig_bound){
				$coord2 = $maximum_contig_bound;
			    }
			}
			my $block_out_region = 0;
			if($repeat_preds == 1){
			    if($prediction_included ==1){
				$block_out_region = 1; # if prediction found, block out
			    }
			    else{
				$block_out_region = 0; #if no prediction, dont block out
			    }
			}
			elsif($repeat_preds ==0){
			    $block_out_region = 1;
			    if($prediction_included == 0){
				if($homolog =~ m/rho/i || $homolog =~ m/opn/i || $homolog =~ m/Tas1/i || $homolog =~ m/rh[0-9]/i){
				    $block_out_region = 0; #if opsin or a taste 1 receptor, don't block out
				}
				else{
				    $block_out_region = 1; #otherwise, block region regardless of prediction found or not
				}
			    }
			}
			if($block_out_region == 1){
			    my $cordlen=($coord2-$coord1)+1;
			    open(INFILE, "$target_nuc");
			    my $tmpname="";my $tmpseq="";
			    while(<INFILE>){
				if($_=~m/\>/){
				    $tmpname=$_;
				}
				else{
				    $tmpseq.=$_;
				}
			    }
			    close INFILE;
			    $coord1 -=1;
			    for(my $j=0;$j<$cordlen;$j++){
				my $locus=$coord1+$j;
				substr($tmpseq,$locus,1)="N"; #block out with Ns
			    }
			    open(OUTFILE, ">$target_nuc");
			    print OUTFILE $tmpname.$tmpseq;
			    $tmpname="";$tmpseq="";
			    close OUTFILE;
			    `makeblastdb -in $target_nuc -dbtype="nucl" -out $refDB`;
			} 
		    }
		}
	    }
	}
	$dna="";
    }
    else{
	chomp $_;
	$dna.=$_;
    }
}
if($remove_intermediate eq "yes"){ #remove intermediate files
    `rm $randomstr*`;
}
open(INX, "$backup_genomefile");#to remove the >X appended to genome file
open(OUTX, ">$tmp_X_file");
while(<INX>){
    if($_=~m/\>PSEUDOCONTIG/){
    }
    else{
	print OUTX $_;
    }
}
system("mv $tmp_X_file $backup_genomefile");
$datestring = localtime();
if ($print_datestring eq "yes"){ #print end date and time
    print "Ending at: ".$datestring."\n";
}
if ($pseudogenome_esl eq "yes"){
    `rm *ssi`; #remove index file 
}



###########################################################################
# SUBROUTINES:
###########################################################################

#blatout subroutine: Takes blast hit and searches for start and stop codons
#Extension loops look for start and stop codons by extending outwards in 5' and 3' directions
#If no START and STOP found with extension - sequence is sent to internal loops
#Internal loops look for start and stop codons by trimming inwards in 3' and 5' directions
sub blatout{
    my @inputs = @_;
    my $blatoutputfile=$inputs[0];
    my $max_length =$inputs[1];
    my $min_length_internal =$inputs[2];
    my $plusminus = $inputs[3];
    my $pseudogene_length = $inputs[4];
    my $percent = $inputs[5];
    my $cover_thresh = $inputs[6];
    my $identity_threshold = $inputs[7];
    my $evalue_inclusion = $inputs[8];
    my $max_length_external = $inputs[9];
    my $blastback_verification = $inputs[10];
    my $percentage_cover_threshold = $inputs[11];
    my $percentage_identity_threshold = $inputs[12]; 
    my $evalue_threshold = $inputs[13];
    my $multiexon = $inputs[14];
    my $minlength = $inputs[15];
    my $reff_db = $inputs[16];
    my $multiout = $inputs[17];
    my $multiblatout = $inputs[18];
    my $opsin_switch = $inputs[19];
    my $opsin_ref = $inputs[20];
    my $classA_HMM = $inputs[21];
    my $HMM_filter = $inputs[22];
    my $classC_HMM = $inputs[23];
    my $vom2r_hmm = $inputs[24];
    my $block_start = "";
    my $block_end = "";
    my @scoords;my @ecoords;#series of variables to be used
    my @finalarray=();
    my $extend_checkpoint = 0; #controlls extension loops
    my $Frame;
    my $tool = "";
    my $cover_threshold;
    my $bu_sstart="";my $bu_send="";
    open(BLATRES, "$blatoutputfile"); #open blast file
    my $tar_gene="";my $receptor_name="";
    my @blatres=<BLATRES>;
    @blatres=split(/Query\=/,join('',@blatres));
    my $EV; my $COVER_rounded; my $ID_rounded; my $blast_program;
    if($multiexon == 1){
	if ($cover_thresh > 200){
	    $cover_threshold = 200;
	}
	else{
	    $cover_threshold = $cover_thresh;
	}
    }
    else{
	$cover_threshold = $cover_thresh;
    }
    my $top_hit=0; #sets to 1 when tophit has been used for prediction, keeps references ordered by signifigance
    foreach my $x(@blatres){
	my $qlen=0;my $qstart=0;my $qend=0;my $sstart=0;my $send=0;my $strand="";my $startdiff=0;my $enddiff=0;
	my @hits=split(/\>/,join('',$x));
	if($hits[0]=~m/Length\=([0-9]+)/){
	    $qlen=$1; 
	}
	if($hits[0]=~m/\_(O[A-Za-z0-9]+)/ ||
	   $hits[0]=~m/\_(V[A-Za-z0-9]+)/ ||
	   $hits[0]=~m/\_(T[A-Za-z0-9]+)/ ||
           $hits[0]=~m/\_(RHO)/ ||
           $hits[0]=~m/\_(rho)/i ||
           $hits[0]=~ m/\_(rhol)/i ||
           $hits[0] =~ m/\_(Taar[A-Za-z0-9]+)/i ||
           $hits[0] =~ m/\_(TAS[0-9A-Za-z]+)/i ||
           $hits[0] =~ m/\_(Vom[0-9A-Za-z]+)/i ||
           $hits[0] =~ m/\_(VN[\S]+)/i ||
           $hits[0] =~ m/\_(V2[\S]+)/i ||
           $hits[0] =~ m/\_(Olf[0-9A-Za-z]+)/i ||
           $hits[0] =~ m/\_(OPN[0-9A-Za-z]+)/i ||
           $hits[0] =~ m/\_(OR[0-9A-Za-z]+)/i ||
           $hits[0] =~ m/\_(Rh[0-9])/i){
	    $receptor_name=$1;
	}
	for(my $i=1;$i<scalar(@hits);$i++){#for all blast HSP
	    my $CONTIG_MAX_LENGTH = ""; #contig length stops extension beyond contig limits
	    if($hits[$i] =~m/Length\=([0-9]+)/){ 
		$CONTIG_MAX_LENGTH = $1; #contig length
	    }
	    my @loci=split(/Score/,$hits[$i]);#split each query chunk by hits
	    foreach my $y(@loci){ 
		if ($top_hit == 1){ #if tophit has already been used, move on to next reference receptor
		}
		else{ #if top hit, proceed to prediction
		    my $ident=0;
		    my $copy=$y;
		    my $evalue_score;
		    if($y=~m/Identities[\s]+\=[\s]+([0-9]+)\/([0-9]+)[\s]+\(([0-9]+)\%/){
			my $cover=$1;my $ident=$3; #get cover and identity scores
			if($y =~ m/Expect[\s]+\=[\s]+([0-9]+.*[0-9]+)/){
			    $evalue_score = $1; #get evalue
			}
			if($cover>=$cover_threshold && $ident>=$identity_threshold && $evalue_score <= $evalue_inclusion){#if at least 200 nucs covered
			    $top_hit = 1; #if meets user determined thresholds, proceed to predict the top hit
			    if($y=~m/Query[\s]+([0-9]+)/){#Get start and end coordinates for query and hit
				$qstart=$1;
			    }
			    if($y=~m/Sbjct[\s]+([0-9]+)/){
				$sstart=$1;
			    }   
			    while($copy=~s/Query[\s]+[0-9]+[\s]+[\S]+[\s]+([0-9]+)//){
				$qend=$1;
			    }
			    while($copy=~s/Sbjct[\s]+[0-9]+[\s]+[\S]+[\s]+([0-9]+)//){
				$send=$1;
			    }
			    if($sstart>$send){#determine if hit is fwd or rev relative to query
				$strand="neg"; #reverse 
			    }
			    else{
				$strand="pos"; #forward 
			    }
			    open(FILE, "$target_nuc");#read contig into memory
			    my $refseq="";my $reffasta="";
			    while(<FILE>){
				if($_=~m/\>/){
				    $reffasta=$_;
				}
				else{
				    $refseq.=$_;
				}
			    }
			    close FILE;
			    my $ref_name=""; #Contig name
			    if($reffasta=~m/\>([\S]+)/){
				$ref_name=$1; #Contig name
			    }
			    my $y1 = "";
			    if($strand eq "pos"){#if the strand was positive relative to the reference
				my $targene_len=0;
				my $genomeseq=uc(substr($refseq,($sstart-1),(($send-$sstart)+1))); #store sequence in $genomeseq using blast coordinates
				my $length_gs = length($genomeseq);
				$bu_sstart=$sstart; $bu_send=$send;#make backups of starting and ending coordinates
				#####
				my $check_opsin_ref = $randomstr."_opsin_reference.fa";
				my $opsin_blast_out = $randomstr."_opsin_blast_out.fa";
				if($opsin_switch == 1){
				    $opsin_ref = $randomstr."_opsin_hint_refs.fa";
				    if(-e $opsin_ref){
					`rm $opsin_ref`;
				    }
				    $sstart -=1000;
				    if($sstart < 1){
					$sstart = 1;
				    }
				    $send +=1000;
				    if($send > $CONTIG_MAX_LENGTH){
					$send = $CONTIG_MAX_LENGTH;
				    }
				    $genomeseq=uc(substr($refseq,($sstart-1),(($send-$sstart)+1))); #store sequence in $genomeseq using blast coordinates
				    $sstart = $bu_sstart;$send = $bu_send;#make backups of starting and ending coordinates
				    open(SS, ">$check_opsin_ref");
				    my $ss = ">opsin\n".$genomeseq;
				    print SS $ss;
				    close SS;
				    `blastn -task 'megablast' -db $reff_db -query $check_opsin_ref -out $opsin_blast_out -num_threads $threads`;
				    my @opsin_hits =&parse_blast_hits($opsin_blast_out);
				    my $count_hits = 0;
				    foreach my $opsinhit(@opsin_hits){
					my $ohit = $opsin_hits[0];
					shift(@opsin_hits);
					if($count_hits <=5){
					    if ($ohit){
						if($ohit =~m/([\S]+)/){
						    my $opsin_name = $1;
						    my $opsin_sequence =  $receptor_key{$opsin_name};
						    open(OPSIN, ">>$opsin_ref");
						    print OPSIN ">".$opsin_name."\n".$opsin_sequence."\n";
						    $count_hits +=1;
						}
					    }
					}
				    }				    
				}
				my $framecheck=0;#frame check variable
				my $frame;
				my $cds;
				$sstart = $sstart - $plusminus;
				$send = $send + $plusminus;
				my $hit_start = "";
				my $hit_end = "";
				my $x1 = "";
				if($sstart < 1){
				    $x1 = 1 - $sstart;
				    $y1 = $plusminus - $x1;
				    $hit_start =  $y1;
				    $sstart = 1;
				}
				elsif($sstart >=1){
				    $hit_start =  $plusminus;
				}
				if($send > $CONTIG_MAX_LENGTH){
				    $hit_end = $hit_start + $length_gs;
				    $send = $CONTIG_MAX_LENGTH;
				}
				elsif($send <= $CONTIG_MAX_LENGTH){
				    $hit_end = $hit_start + $length_gs;
				}
				$genomeseq=uc(substr($refseq,($sstart-1),(($send-$sstart)+1)));
				my $max_bound = length($genomeseq);
				my $Nblock= "N" x 30; #10 codons of Ns. If extension reaches this, it is likely extending into a blocked out gene.
				$count +=1;
				my $header = ">Hit".$count."\n";
				my $genome_seq = $header.$genomeseq."|".$hit_start."|".$hit_end."|".$receptorfile."|".$augustus_species."|".$multiexon."|".$randomstr."|".$opsin_switch."|".$opsin_ref;
			        my @augustus_preds =&augustus_predict($genome_seq);
				my $aug_avail = 0;
				my $aug_prediction = $augustus_preds[0];
				if($aug_prediction){
				    my @split_aug_info = split(/\|/, $aug_prediction);
				    my $aug_start = $split_aug_info[1];
				    my $aug_end = $split_aug_info[2];
				    my $augustus_prediction = $split_aug_info[0];
				    $genomeseq = $augustus_prediction;
				    if(length($genomeseq) < $minlength){
					$aug_avail = 0;
				    }
				    elsif($genomeseq =~ m/$Nblock/i){
					$aug_avail = 0;
				    }
				    else{
					$aug_avail = 1;
					$tool = "augustus";
				    }
				    if($aug_avail ==1){
					my $startdiff = "";
					my $enddiff = "";
					if($aug_start > $hit_start){
					    $startdiff = $aug_start - $hit_start;
					    $startdiff -=1;
					    $sstart = $bu_sstart;
					    $sstart += $startdiff;
					}
					if($aug_start < $hit_start){
					    $startdiff = $hit_start - $aug_start;
					    $startdiff +=1;
					    $sstart = $bu_sstart;
					    $sstart -= $startdiff;
					}
					if($aug_start == $hit_start){
					    $sstart = $bu_sstart;
					}
					if($aug_end > $hit_end){
					    $enddiff = $max_bound - $aug_end;
					    $send -= $enddiff;
					}
					if($aug_end < $hit_end){
					    $enddiff = $hit_end - $aug_end;
					    $send = $bu_send;
					    $send -= $enddiff;
					}
					if($aug_end == $hit_end){
					    $send = $bu_send;
					}
				    }
				}
				if($aug_avail == 0 && $multiexon == 1){
				    $genomeseq = "NOPREDICTION";
				    $sstart = $bu_sstart;
                                    $send = $bu_send;
				}
				if($aug_avail == 0 && $multiexon == 0){
				    $tool = "sensommatic";
				    $sstart = $bu_sstart;
				    $send = $bu_send;
				    $genomeseq=uc(substr($refseq,($sstart-1),(($send-$sstart)+1)));
				    if(length($genomeseq) % 3 ==0){ #if divisible by 3, proceed
					$cds = $genomeseq;
				    }else{ #if not divisible by three
					$cds = substr($genomeseq, 1, length($genomeseq));
					$sstart +=1; #cut in 1, now divisible by 3?
					if (length($cds) % 3 != 0) {
					    $cds = substr($genomeseq, 2, length($genomeseq));
					    $sstart +=1; #cut in 1, now divisible by 3?
					}
				    }
				    my $bu_start_if = $sstart; my $bu_end_if = $send; #back up coordinates in frame
				    my $Frame=&getframe_tblastn($cds);
				    if ($Frame eq "Frame1"){ 
					$genomeseq = $cds; #Frame 1, no adjustments needed
				    }
				    elsif($Frame eq "Frame2"){
					$genomeseq = substr($cds, 1, -2); #Frame 2, start +1, end -2
					$sstart += 1;
					$send -=2;
				    }elsif($Frame eq "Frame3"){
					$genomeseq = substr($cds, 2, -1); #Frame3 3, start +2, end -1
					$sstart +=2;
					$send -=1;
				    }
				    my $percent = 0.25;
				    my $initial_contract = length($genomeseq) * $percent;
				    my $either_end = $initial_contract / 2;
				    $either_end = sprintf("%.0f", $either_end);
				    if ($either_end % 3 !=0){
					$either_end-=1;
					if ($either_end % 3!=0){
					    $either_end -=1;
					}
				    }
				    ####################################
				    # Positive Extension loops         #
				    # search for START and STOP codons #
				    ####################################
				    while($framecheck==0){ #while start and stop not found
					if($genomeseq=~m/^ATG/ && $genomeseq=~m/TGA$/ || $genomeseq=~m/^ATG/ && $genomeseq=~m/TAG$/ || $genomeseq=~m/^ATG/ && $genomeseq=~m/TAA$/){
					    $framecheck=1;#if start and stop codon, exit extension loops
					}
					elsif($genomeseq =~ m/^$Nblock/ || $genomeseq =~ m/$Nblock$/){ #if sequence extends into a blocked out gene (NNNN), stop extending
					    $framecheck =1; #if extension reaches 10 codons of Ns in either direction, stop extension
					}
					##########################
					# Search for STOP codon  #
					##########################
					elsif($genomeseq=~m/^ATG/){#if seq has START codon, but no STOP codon
					    if(length($genomeseq) % 3 ==0){#if the sequence is divisible by 3
						$send+=3; #extend in the 3' direction to find STOP codon
						if ($send > $CONTIG_MAX_LENGTH) { #if you exceed contig bound in extension, revert back and exit extension
						    $send -=3;#Go back to max length in frame. No stop will be found.
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1)));
						    $framecheck =1; # In this case, no stop will be found, exit extension and send to internal search.
						}
						else{ #if does not exceed contig bound, proceed with search for STOP codon
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1)));
						}
					    }
					    else{ #if the sequence is not divisible by 3 (this condition should never be met - here as buffer)
						$send+=1;#if the sequences is not divisible by 3, extend by one nucleotide until divisible by 3.
						if ($send <= $CONTIG_MAX_LENGTH){ #if doesn't cross contig end, all is good, continue and check frame
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1)));
						}
						elsif($send > $CONTIG_MAX_LENGTH) { #otherwise, if crosses contig end.
						    $send -=2; #Go back two nucleotides until divisible by three. No stop will be found.
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1))); #get sequence
						    if(length($genomeseq) % 3 !=0){ #if still not divisible by three
							$send -=1; #go inward one more step. Now will be divisible by 3.
							$genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1))); # Take new substring. Divisible by 3.
							$framecheck =1; #gets sent to internal search
						    }
						    elsif(length($genomeseq) % 3 ==0){ #if divisible by three 
							$framecheck =1; #send to internal search
						    }
						}
					    }
					}
					##########################
					# Search for START codon #
					###########################
					elsif($genomeseq=~m/TAA$/ || $genomeseq=~m/TGA$/ || $genomeseq=~m/TAG$/){#alternatively, if no start codon, but has a stop codon
					    if(length($genomeseq) % 3 ==0){# if divisisble by 3, extend in 5' direction to search for START codon
						$sstart-=3; #extend in 5' direction
						if ($sstart < 1) { #if extension crosses contig bound
						    $sstart +=3; #revert back (+3)   
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1)));
						    $framecheck =1; #stop extension and send to to internal search                                                                                                                                     
						}
						elsif($sstart >=1){ #if does not exceed contig bound
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1))); #keep extending and search for START
						}
					    }
					    else{
						$sstart-=1;#if not divisible by 3 (condition should not be met - included as buffer), go back by 1 until it is
						if ($sstart >=1){
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1)));
						}elsif($sstart <1){
						    $sstart +=2; #If crosses the start of contig, go inward one nucleotide (-1 + 2 = +1)
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1)));
						    if(length($genomeseq) % 3 !=0){ #If not divisible by 3
							$sstart +=1; #go inwards one more nucleotide. Now divisible by three.
							$genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1))); #Take new substring.
							$framecheck =1; #send to internal search
						    }
						    elsif(length($genomeseq) % 3 ==0){ #if divisible by three
							$framecheck =1; #send to internal search
						    }
						}
					    }
					}
					#############################
					# Search for START and STOP #
					##############################
					else{#if the sequence has neither start or stop codon
					    if(length($genomeseq) % 3 ==0){#if divisible by 3
						$extend_checkpoint = 0;
						$send+=3;$sstart-=3; #extend in both directions (3' and 5' directions)
						if($send > $CONTIG_MAX_LENGTH){ #if crosses end of contig, go back to longest transcript in frame
						    $send -=3; #go back to longest transcript in frame
						    $extend_checkpoint += 0.5; #for framecheck statement
						}
						if($sstart < 1){ #if cross the start position, revert back to longest transcript in frame
						    $sstart +=3; #revert back to longest transcript in frame
						    $extend_checkpoint += 0.5; #for framecheck statement
						}
						$genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1))); #take the coordinates (adjusted or reverted)
						if ($extend_checkpoint == 1){ # if both ends cross the contig, stop extending in both directions and set framecheck to 1.
						    $framecheck =1; #send to internal loops
						}
						elsif($extend_checkpoint < 1){ #if at least one end does not exceed contig bound, continue extending in that direction
						}
					    }
					    else{ #if not divisible by 3 (should not meet this condition - included as buffer)
						$sstart-=1;$send+=1;#extend by 1 site until divisible by three
						if($send > $CONTIG_MAX_LENGTH && $sstart < 1){ #if crosses both ends
						    $send -=1; #revert end coord
						    $sstart +=1; #revert start coord
						    $send -=1; #take 1 base from end since we don't have a stop
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1))); #take sequence
						    if (length($genomeseq) % 3 ==0){ #if removing base results in divisible by 3, move on.
							$framecheck =1;
						    }
						    elsif (length($genomeseq) % 3 !=0) { #else if not divisible by 3
							if ($genomeseq =~m/TAA$/ || $genomeseq=~m/TGA$/ || $genomeseq=~m/TAG$/){ #if seq ends in stop, start removing from beginning
							    $sstart +=1;
							    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1)));
							    if (length($genomeseq) % 3 ==0) { #now divisible by three
								$framecheck =1; #send to internal search
							    } 
							}
							else{ #if no STOP, remove one more nucelotide from the end (5')
							    $send -=1;
							    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1))); #take sequence
							    if (length($genomeseq) % 3 ==0){ #now divisible by three
								$framecheck=1; #send to internal search
							    }
							}
						    }
						}
						elsif($send > $CONTIG_MAX_LENGTH && $sstart >=1){ #if crosses contig end (not start), adjust only in 3' direction
						    $send -=1; #restore end
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1))); #take sequence
						    if(length($genomeseq) % 3 ==0){ #if divisible by three
							$framecheck=1; #send to internal search
						    }
						}
						elsif($sstart <1 && $send <= $CONTIG_MAX_LENGTH){ #if crosses contig start (not end), adjust only in 5' direction
						    $sstart +=1; #restore start
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1))); #take sequence
						    if(length($genomeseq) % 3 ==0){ #if divisible by three
							$framecheck=1; #send to internal search
						    }
						}
						elsif($sstart >=1 && $send <= $CONTIG_MAX_LENGTH){ #if does not exceed either contig bound
						    $sstart +=1; #adjust by 1 until divisible by 3
						    $genomeseq=uc(substr($refseq,$sstart-1,(($send-$sstart)+1))); #take sequence
						}
					    }
					}
					if(length($genomeseq) > $max_length && $framecheck ==0){ #keep extending until START and STOP found OR exceeds user determined maximum length
					    $framecheck=1; #exit extension loops
					}
				    }
				    #################################
				    # BLASTBACK verification step:  #
				    # Prevents over extension       #
				    #################################
				    my $revert=0;
				    my @blastbackin = ();
				    push(@blastbackin, $genomeseq);
				    push(@blastbackin, $blastback_verification);
				    push(@blastbackin, $percentage_cover_threshold);
				    push(@blastbackin, $percentage_identity_threshold);
				    push(@blastbackin, $evalue_threshold);
				    if($blastback_verification ==1){
					my @blast_values=&blastback(@blastbackin);
				        $revert = $blast_values[0];
				        if ($revert ==1){ #If after extension loops, blast hit is too diluted
					    $genomeseq = uc(substr($refseq,($bu_start_if -1),(($bu_end_if-$bu_start_if)+1))); #revert to backup coordinates and send to internal search
					    $sstart = $bu_start_if; #adjust start coordinate accordingly
					    $send = $bu_end_if; #adjust end coordinate accordingly
					}
					else{ #if prediction is still within cover and identity thresholds, keep extension
					}
				    }
				    ############################
				    # Positive internal search #
				    ############################
				    if($genomeseq=~m/^ATG/ && $genomeseq=~m/TGA$/ || $genomeseq=~m/^ATG/ && $genomeseq=~m/TAG$/ || $genomeseq=~m/^ATG/ && $genomeseq=~m/TAA$/){ #if start and stop, don't enter internal search
				    } 
				    elsif($genomeseq!~m/^ATG/ || $genomeseq!~m/TGA$/ || $genomeseq!~m/TAG$/ || $genomeseq!~m/TAA$/){ #if either no START codon or no STOP codon
					if($genomeseq!~m/^ATG/ && $genomeseq!~m/TGA$/ && $genomeseq!~m/TAG$/ && $genomeseq!~m/TAA$/){#if no start AND no stop codon
					    $genomeseq=uc(substr($refseq,($bu_start_if -1),(($bu_end_if-$bu_start_if)+1))); #revert to frame adjusted back up coordinates
					    $sstart = $bu_start_if;
					    $send = $bu_end_if;
					}
					my $subseq="";#look internally
					if(length($genomeseq) % 3 !=0){#if not divisible by three, adjust to be divisible by 3 (condition should never be met - included as buffer)
					    my $newseq=uc(substr($genomeseq,0,((length($genomeseq))-1))); #trim end by 1
					    if (length($newseq) % 3 ==0){
						$genomeseq = $newseq; #if now divisible by 3, good to go
						$send -=1;
					    }
					    else{ #if still not divisible by 3
						my $newseq=uc(substr($genomeseq,0,((length($genomeseq))-2))); #trim end by 1
						$genomeseq = $newseq; #now divisible by 3
						$send -=2;
					    }
					}
					if(length($genomeseq) % 3 ==0){ #if divisble by 3
					    my $Backup = $genomeseq;
					    my $startfound=0;
					    my $start_reduce =0; #0 + 3 will give 3rd index --> 4th position, we want 4th nucelotide to access second codon
					    if ($genomeseq =~ m/^ATG/){ #if start is found
						$subseq = $genomeseq; 
						$startfound =1; #exit start trim loops
					    }
					    else{ #if no start codon
						while($startfound ==0){ #trim internally in 3 -> 5' direction until START codon is found
						    $start_reduce +=3; 
						    $sstart +=3;
						    $subseq = uc(substr($genomeseq,$start_reduce,length($genomeseq))); #get sequence
						    if ($subseq =~ m/^ATG/){ #if START is found 
							$startfound =1; #exit start trim loops
						    }
						    if (length($subseq) < $min_length_internal){ #if internal trim is too short
							$startfound =1; #exit start trim loops
						    }
						}
					    }  
					}
					if (length($subseq) < $min_length_internal || $subseq !~ m/^ATG/){ #if too short OR start codon not found
					    $subseq = uc(substr($refseq,($bu_start_if -1),(($send-$bu_start_if)+1))); #revert to in frame start and keep extended end
					    $sstart = $bu_start_if; #update coordinates accordingly
					}
					my $backup = $subseq;
					my $codonstring="";
					my $end_reduce = length($subseq);
					my $end_trim = 0; #trim internally from end coordinate in 3' -> 5' to search for STOP codon
					if ($subseq =~ m/TAG$/ || $subseq =~ m/TAA$/ || $subseq =~ m/TGA$/){ #if STOP codon found
					    $codonstring = $subseq; 
					    $end_trim = 1; #exit end trim loops
					}
					else{ #if no STOP codon
					    while($end_trim == 0){ #while no STOP codon
						$end_reduce -= 3; #trim end by 3
						$send -=3; # adjust coordinate accordingly
						$codonstring= uc(substr($subseq,0,$end_reduce)); #get sequence 
						if($codonstring =~ m/TAG$/ || $codonstring =~ m/TAA$/ || $codonstring =~ m/TGA$/){ #if STOP codon found
						    $end_trim = 1; #exit end trim loops
						}
						if (length($codonstring) < $min_length_internal){ #if internal trim is too short
						    $end_trim = 1; #exit end trim loops
						}
					    }
					}
					if($codonstring =~ m/TAG$/ || $codonstring =~ m/TAA$/ || $codonstring =~ m/TGA$/){ #if STOP is found
					    if(length($codonstring)<$min_length_internal){ #but sequence is too short
						$genomeseq = uc(substr($refseq,($sstart -1),(($bu_end_if-$sstart)+1))); #revert end back to backup, but retain work done from START search
						$send = $bu_end_if; #adjust end coordinate accordingly
					    }
					    else{
						$genomeseq = $codonstring; #if STOP found AND is long enough, keep internal end search
					    }
					}
					else{ #if no STOP codon found
					    $genomeseq = uc(substr($refseq,($sstart -1),(($bu_end_if-$sstart)+1))); #revert to back up end coordinare regardless of length
					    $send = $bu_end_if; #adjust end coordinate accordingly
					}
				    } 
				    ################################
				    # BLASTBACK verification step: #
				    # Prevents over trimming       #
				    ################################
				    $revert=0; my $COV =0; my $IDY =0; $EV=1; my $exclude=0; #my $COVER_rounded; my $ID_rounded; my $blast_program;
				    @blastbackin = ();
				    push(@blastbackin, $genomeseq);
				    push(@blastbackin, $blastback_verification);
				    push(@blastbackin, $percentage_cover_threshold);
				    push(@blastbackin, $percentage_identity_threshold);
				    push(@blastbackin, $evalue_threshold);				    
				    if($blastback_verification ==1){
					my @blast_values=&blastback(@blastbackin);
					$revert = $blast_values[0]; #revert status 
					$COV = $blast_values[1]; #percentage cover score
					$IDY = $blast_values[2]; #percentage identity score
					$EV = $blast_values[3]; #evalue
					$blast_program = $blast_values[4]; #megablast or nblast (if no hit with megablast, nblast is called)
					$exclude =0;
					if ($revert ==1){ #if prediction does not meet cover and identity thresholds after internal/extension loops
					    $genomeseq = uc(substr($refseq,($bu_start_if -1),(($bu_end_if-$bu_start_if)+1))); #revert to in-frame adjusted back up coordinates
					    $sstart = $bu_start_if; #adjust start accordingly
					    $send = $bu_end_if; #adjust end accordingly
					    @blastbackin = ();
					    push(@blastbackin, $genomeseq);
					    push(@blastbackin, $blastback_verification);
					    push(@blastbackin, $percentage_cover_threshold);
					    push(@blastbackin, $percentage_identity_threshold);
					    push(@blastbackin, $evalue_threshold);
					    @blast_values=&blastback(@blastbackin); #do one more blastback step, with backup coorinates 
					    $revert = $blast_values[0]; #if revert is 1 at this point, sequence will be excluded
					    $COV = $blast_values[1]; #new percentage cover score
					    $IDY = $blast_values[2]; #new percentage identity score
					    $EV = $blast_values[3]; #new evalue 
					    $blast_program = $blast_values[4]; #megablast or nblast
					    if ($revert ==1){ #if revert is still 1 
						$exclude = 1; #exclude sequence
					    }
					}
					$COVER_rounded = sprintf("%.4f", $COV); #round percentage cover score to 4 decimals for annotation
					$ID_rounded = sprintf("%.4f", $IDY); #round percentage identity score to 4 decimals for annotation
				    }
				}
				##########################################################
				#Hmm filter - remove prediction if no 7transmembrane hit #
				##########################################################
				my $hmmremove = "no";
				if($HMM_filter eq "Yes" || $HMM_filter eq "yes"){
				    if($genomeseq ne "NOPREDICTION"){
					my $Protein_prediction = ">Prediction_in\n";
					for (my $d=0;$d<(length($genomeseq) -2);$d+=3)
					{
					    my $Qcodon = substr($genomeseq,$d,3);
					    $Protein_prediction.=&dna2prot($Qcodon); #translate reference receptor using dna2protein
					}
					my $prediction_out =  $randomstr."_hmm_input_protein.fa";
					open(PPOUT, ">$prediction_out");
					print PPOUT $Protein_prediction;
					close PPOUT;
					my @hmm_inputs = ();
					push(@hmm_inputs, $classA_HMM);
					push(@hmm_inputs, $prediction_out);
					push(@hmm_inputs, $randomstr);
					my $hmm_result =&hmmsearch(@hmm_inputs);
					if($hmm_result == 0){
                                            my @hmm_inputs = ();
                                            push (@hmm_inputs, $classC_hmm);
                                            push (@hmm_inputs, $prediction_out);
					    push(@hmm_inputs, $randomstr);
                                            my $hmm_result =&hmmsearch(@hmm_inputs);
                                            if($hmm_result == 0){
						my @hmm_inputs = ();
						push (@hmm_inputs, $vom2r_hmm);
						push (@hmm_inputs, $prediction_out);
						push(@hmm_inputs, $randomstr);
						my $hmm_result =&hmmsearch(@hmm_inputs);
						if($hmm_result == 0){
						    $hmmremove = "yes";
						}
                                            }
                                        }	
				    }   
				}
				#########################################################
				#Blast prediction back to reference file to re-annotate #
				#########################################################
				if($genomeseq ne "NOPREDICTION"){
				    open(mOUT, ">$multiout");
				    print mOUT ">multi_out\n".$genomeseq;
				    close mOUT;
				    `blastn -task 'megablast' -db $reff_db -query $multiout -out $multiblatout -num_threads $threads`;
				    my @qhits =&parse_blast_hits($multiblatout);
				    my $thit = $qhits[0];
				    if ($thit){
					if($thit =~ m/\_(O[A-Za-z0-9]+)/ ||
					   $thit =~m/\_(V[A-Za-z0-9]+)/ ||
					   $thit =~m/\_(T[A-Za-z0-9]+)/ ||
					   $thit =~m/\_(RHO)/ ||
					   $thit =~m/\_(rho)/i ||
					   $thit =~ m/\_(rhol)/i ||
					   $thit =~ m/\_(Taar[A-Za-z0-9]+)/i ||
					   $thit =~ m/\_(TAS[0-9A-Za-z]+)/i ||
					   $thit =~ m/\_(Vom[0-9A-Za-z]+)/i ||
					   $thit =~ m/\_(VN[\S]+)/i ||
					   $thit =~ m/\_(V2[\S]+)/i ||
					   $thit =~ m/\_(Olf[0-9A-Za-z]+)/i ||
					   $thit =~ m/\_(OPN[0-9A-Za-z]+)/i ||
					   $thit =~ m/\_(OR[0-9A-Za-z]+)/i ||
					   $thit =~ m/\_(Rh[0-9])/i){
					    my $name = $1;
					    $name =~ s/\_//g;
					    $receptor_name = $name;
					}
				    }			    
				}
				################################################
				#checkframe subroutine to annotate prediction  #
				################################################		
				my @status=&checkframe($genomeseq); #use checkframe subroutineto get appropriate annotation status
				my $stat=$status[0]; #get annotation status
				if($truncated == 1){ #if truncated option is ON
				    if($genomeseq =~ m/TAG$/ || $genomeseq =~ m/TAA$/ || $genomeseq =~ m/TGA$/){ #if STOP codon
					if($include_score ==1 && $tool eq "sensommatic"){ #if blastback annotation score option is ON
					    $tar_gene=">Seq".$ref_name."_".$receptor_name."_status_fwd_".$sstart."_".$send."|".$COVER_rounded."|".$ID_rounded."|".$EV."|".$blast_program."|".$tool."\n".uc($genomeseq)."\n"; #add blastback scores to annotation
					}
					elsif($include_score ==0 || $tool eq "augustus"){ #if blastback annotation score option is OFF
					    $tar_gene=">Seq".$ref_name."_".$receptor_name."_status_fwd_".$sstart."_".$send."|".$tool."\n".uc($genomeseq)."\n"; #dont include blastback scores to annotation
					}
				    }
				    else{ #if no STOP codon, add truncated annotation
					if($include_score ==1 && $tool eq "sensommatic"){ #if blastback annotation score option is ON
					    $tar_gene=">Seq".$ref_name."_".$receptor_name."_status_".$annotation_truncated."_fwd_".$sstart."_".$send."|".$COVER_rounded."|".$ID_rounded."|".$EV."|".$blast_program."|".$tool."\n".uc($genomeseq)."\n"; #add blastback scores to annotation
					}
					elsif($include_score ==0 || $tool eq "augustus"){ #if blastback annotation score option is OFF
					    $tar_gene=">Seq".$ref_name."_".$receptor_name."_status_".$annotation_truncated."_fwd_".$sstart."_".$send."|".$tool."\n".uc($genomeseq)."\n"; #dont include blastback scores to annotation
					}
				    }  
				}
				else{ #if truncated option is OFF
				    if($include_score ==1 && $tool eq "sensommatic"){ #if blastback annotation score option is ON
					$tar_gene=">Seq".$ref_name."_".$receptor_name."_status_fwd_".$sstart."_".$send."|".$COVER_rounded."|".$ID_rounded."|".$EV."|".$blast_program."|".$tool."\n".uc($genomeseq)."\n"; #add blastback scores to annotation
				    }
				    elsif($include_score ==0 || $tool eq "augustus"){ #if blastback annotation score option is OFF
					$tar_gene=">Seq".$ref_name."_".$receptor_name."_status_fwd_".$sstart."_".$send."|".$tool."\n".uc($genomeseq)."\n"; #dont include blastback scores to annotation
				    }				
				}
				my $length_check = 0;
				$targene_len=length($genomeseq); #get length of prediction
				if ($targene_len < $pseudogene_length){ #if length is below user determined threshold, annotate as pseudogene
				    $tar_gene=~s/status/$annotation_short/; #If prediction is shorter than $pseudogene_length, gene will be annotated as pseudogene
				    $length_check = 1;
				}
				if($stat eq "1"){
				    unless($length_check ==1){
					$tar_gene=~s/status/$annotation_1/; #START codon && no in frame stop codons
				    }
				}
				if($stat eq "2"){
				    unless($length_check ==1){
					$tar_gene=~s/status/$annotation_2/; #no START codon && no stop codons in any frame   
				    }
				}
				if($stat eq "3"){
				    unless($length_check ==1){
					$tar_gene=~s/status/$annotation_3/; #no START codon && no in frame stop codons   
				    }
				}
				if($stat eq "4"){
				    if($length_check ==1){
					my $new_annotation = $annotation_4."_short";
					$tar_gene =~ s/$annotation_short/$new_annotation/;
				    }
				    else{
					$tar_gene=~s/status/$annotation_4/; #START codon && in frame stop codon
				    }  
				}
				if($stat eq "5"){
				    if($length_check ==1){
					my $new_annotation = $annotation_5."_short";
					$tar_gene =~ s/$annotation_short/$new_annotation/;
				    }
				    else{
					$tar_gene=~s/status/$annotation_5/; #START codon && in frame stop codon                                                                                                                                                             
				    }
				}
				if($stat eq "6"){
				    if($length_check ==1){
					my $new_annotation = $annotation_6."_short";
					$tar_gene =~ s/$annotation_short/$new_annotation/;
				    }
				    else{
					$tar_gene=~s/status/$annotation_6/; #START codon && in frame stop codon                                                                                                                                                             
				    }
				}
				if($stat eq "7"){
				    $tar_gene=~s/\>Seq/>Ignore/;
				}
				if($targene_len>=$minlength){ #finally, if prediction is greater than user determined threshold
				    if($hmmremove ne "yes"){
					push @finalarray,$tar_gene; #include prediction in final output
				    }
				    else{
					if($tar_gene =~ m/([^\n]+)\n([\S]+)/){
					    my $hmmheader = $1;
					    my $hmmseq = $2;
					    $hmmheader .="_hmm_remove";
					    my $hmm_targene = $hmmheader."\n".$hmmseq;
					    push @finalarray, $hmm_targene;
					}
				    }
				}
				else{#otherwise, prediction will be exlcuded
				    push @finalarray, "NOPREDICTION";
				}
				$targene_len=0;
			    }
			    #################
			    #Reverse strand #
			    #################
			    else{#If the reference was reversed relative to query,do the same as above but acocunting for this
				my $targene_len=0;
				my $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #store sequence in $rev_seq using blast coordinates
				$rev_seq=reverse($rev_seq); #reverse sequence
				$rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
				my $length_rs = length($rev_seq);
				my $framecheck=0; #frame check variable
				$bu_sstart=$sstart;$bu_send=$send; #make backups of starting and ending coordinates
				my $check_opsin_ref = $randomstr."_opsin_reference.fa";
				my $opsin_blast_out = $randomstr."_opsin_blast_out.fa";
				if($opsin_switch == 1){
				    $opsin_ref = $randomstr."_opsin_hint_refs.fa";
				    if(-e $opsin_ref){
					`rm $opsin_ref`;
				    }
				    $send -=1000;
				    if($send < 1){
					$send = 1;
				    }
				    $sstart +=1000;
				    if($sstart > $CONTIG_MAX_LENGTH){
					$sstart = $CONTIG_MAX_LENGTH;
				    }
				    my $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #store sequence in $rev_seq using blast coordinates
				    $rev_seq=reverse($rev_seq); #reverse sequence
				    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
				    $sstart = $bu_sstart;$send = $bu_send;#make backups of starting and ending coordinates
				    open(SS, ">$check_opsin_ref");
				    my $ss = ">opsin\n".$rev_seq;
				    print SS $ss;
				    close SS;
				    `blastn -task 'megablast' -db $reff_db -query $check_opsin_ref -out $opsin_blast_out -num_threads $threads`;
				    my @opsin_hits =&parse_blast_hits($opsin_blast_out);
				    my $count_hits = 0;
				    foreach my $opsinhit(@opsin_hits){
					my $ohit = $opsin_hits[0];
					shift(@opsin_hits);
					if($count_hits <=5){
					    if ($ohit){
						if($ohit =~m/([\S]+)/){
						    my $opsin_name = $1;
						    my $opsin_sequence =  $receptor_key{$opsin_name};
						    open(OPSIN, ">>$opsin_ref");
						    print OPSIN ">".$opsin_name."\n".$opsin_sequence."\n";
						    $count_hits +=1;
						}
					    }
					}
				    }				    
				}
				my $frame;
				my $cds;
				$sstart = $sstart + $plusminus;
				$send = $send - $plusminus;
				my $hit_start = "";
				my $hit_end = "";
				my $y1 = "";
				my $x1 = "";
				if($send < 1){
				    $x1 = 1 - $sstart;
				    $y1 = $plusminus - $x1;
				    $hit_start = $y1;
				    $send = 1;
				}
				elsif($send >=1){
				    $hit_start = $plusminus;
				}
				if($sstart > $CONTIG_MAX_LENGTH){
				    $hit_end = $hit_start + $length_rs;
				    $sstart = $CONTIG_MAX_LENGTH;
				}
				elsif($sstart <= $CONTIG_MAX_LENGTH){
				    $hit_end = $hit_start + $length_rs;
				}	
				$rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1)));
				$rev_seq=reverse($rev_seq); #reverse sequence
				$rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
				my $max_bound = length($rev_seq);
				my $Nblock= "N" x 30; #10 codons of Ns. If extension reaches this, it is likely extending into a blocked out gene.
				$count +=1;
				my $header = ">Hit".$count."\n";
				my $reverse_seq = $header.$rev_seq."|".$hit_start."|".$hit_end."|".$receptorfile."|".$augustus_species."|".$multiexon."|".$randomstr."|".$opsin_switch."|".$opsin_ref;;
			        my @augustus_preds =&augustus_predict($reverse_seq);
				my $aug_avail = 0;
				my $aug_rev_prediction = $augustus_preds[0];			
				if($aug_rev_prediction){
				    my @split_aug_info = split(/\|/, $aug_rev_prediction);
				    my $aug_start = $split_aug_info[1];
				    my $aug_end = $split_aug_info[2];
				    my $augustus_rev_prediction = $split_aug_info[0];
				    $rev_seq = $augustus_rev_prediction;
				    if(length($rev_seq) < $minlength){
					$aug_avail = 0;
				    }
				    elsif($rev_seq =~ m/$Nblock/){
					$aug_avail =0;
				    }
				    else{
					$aug_avail = 1;
					$tool = "augustus";
				    }
				    if($aug_avail == 1){
					my $startdiff = "";
					my $enddiff = "";
					if($aug_start > $hit_start){ 
					    $startdiff = $aug_start - $hit_start;
					    $startdiff -=1;
					    $sstart = $bu_sstart;
					    $sstart -= $startdiff;
					}
					if($aug_start < $hit_start){
					    $startdiff = $hit_start - $aug_start;
					    $startdiff +=1;
					    $sstart = $bu_sstart;
					    $sstart += $startdiff;
					}
					if($aug_start == $hit_start){
					    $sstart = $bu_sstart;
					}
					if($aug_end > $hit_end){
					    $enddiff = $aug_end - $hit_end;
					    $send = $bu_send;
					    $send -= $enddiff;
					}
					if($aug_end < $hit_end){ 
					    $enddiff = $hit_end - $aug_end;
					    $send = $bu_send;
					    $send += $enddiff;
					}
					if($aug_end == $hit_end){
					    $send = $bu_send;
					}
				    }
				}
				if($aug_avail == 0 && $multiexon == 1){
                                    $rev_seq = "NOPREDICTION";
				    $sstart = $bu_sstart;
                                    $send = $bu_send;
                                }
				if($aug_avail == 0 && $multiexon == 0){
				    $tool = "sensommatic";
				    $sstart = $bu_sstart;
				    $send = $bu_send;
				    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #store sequence in $rev_seq using blast coordinates
				    $rev_seq=reverse($rev_seq); #reverse sequence
				    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
				    if(length($rev_seq) % 3 ==0){ #if divisible by 3, proceed
					$cds = $rev_seq;
				    }else{ #if not divisible by three
					$cds = substr($rev_seq, 1, length($rev_seq)); 
					$sstart -=1; #cut in 1, now divisible by 3?
					if (length($cds) % 3 != 0) {
					    $cds = substr($rev_seq, 2, length($rev_seq));
					    $sstart -=1; #cut in 1, now divisible by 3? 
					}
				    }
				    my $Frame=&getframe_tblastn($cds);
				    if ($Frame eq "Frame1"){
					$rev_seq = $cds; #Frame 1, no adjustments needed
				    }
				    elsif($Frame eq "Frame2"){
					$rev_seq = substr($cds, 1, -2); #Frame 2, start +1, end -2
					$sstart -=1;
					$send += 2; 
				    }elsif($Frame eq "Frame3"){
					$rev_seq = substr($cds, 2, -1); #Frame3 3, start +2, end -1
					$sstart -=2; 
					$send +=1; 
				    }
				    my $percent = 0.25;
				    my $initial_contract = length($rev_seq) * $percent;
				    my $either_end = $initial_contract / 2;
				    $either_end = sprintf("%.0f", $either_end);
				    if ($either_end % 3 !=0){
					$either_end-=1;
					if ($either_end % 3!=0){
					    $either_end -=1;
					}
				    }
				    $sstart -= $either_end;
				    $send += $either_end;
				    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1)));
				    $rev_seq=reverse($rev_seq); #reverse sequence
				    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
				    my $bu_start_if = $sstart; my $bu_end_if = $send; #back up coordinates in frame
				    my $Nblock= "N" x 30; #10 codons of Ns. If extension reaches this, it is likely extending into a blocked out gene.
				    ###################################
				    # Reverse Extension loops          #
				    # search for START and STOP codons #
				    ####################################
				    while($framecheck==0){ #while start and stop not found
					if($rev_seq=~m/^ATG/ && $rev_seq=~m/TGA$/ || $rev_seq=~m/^ATG/ && $rev_seq=~m/TAG$/ || $rev_seq=~m/^ATG/ && $rev_seq=~m/TAA$/){ 
					    $framecheck=1;#if start and stop codon, exit extension loops
					}
					elsif($rev_seq =~ m/^$Nblock/ || $rev_seq =~ m/$Nblock$/){ #if sequence extends into a blocked out gene (NNNN), stop extending
					    $framecheck =1; #if extension reaches 10 codons of Ns in either direction, stop extension
					}
					############################
					# Search for STOP codon    #
					############################
					elsif($rev_seq =~m/^ATG/){ #if seq has START codon, but no STOP codon
					    if(length($rev_seq) % 3 ==0){ #if the sequence is divisible by 3
						$send -=3; #extend in the 3' direction to find STOP codon  
						if($send <1 ){ #if extension crosses contig bound
						    $send+=3; #revert back (+3)
						    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1)));
						    $rev_seq=reverse($rev_seq); #reverse 
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						    $framecheck =1; #stop extension and send to to internal search 
						}
						elsif($send >=1){ #if does not exceed contig bound
						    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #keep extending and search for STOP
						    $rev_seq=reverse($rev_seq); #reverse
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						}
					    }
					    else{ 
						$send-=1; #if not divisible by 3, go back by 1 until it is
						if($send<1){ 
						    $send+=2;#If crosses the start of contig, go inward one nucleotide (-1 + 2 = +1) 
						    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #take new sequence
						    $rev_seq=reverse($rev_seq); #reverse it
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						    if(length($rev_seq) % 3 !=0){ #if not divisible by 3
							$send +=1; #go inwards one more nucleotide. Now divisible by three.
							$rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #take new sequence
							$rev_seq=reverse($rev_seq); #reverse it
							$rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
							$framecheck =1; #send to internal search
						    }
						    elsif(length($rev_seq) % 3 ==0){ #if divisible by 3
							$framecheck =1; #send to internal search
						    }
						}
						elsif($send >=1){ #if does not exceed contig bound
						    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #keep extending and search for STOP
						    $rev_seq=reverse($rev_seq); #reverse it
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						}
					    }
					}
					##########################
					# Search for START codon #
					##########################
					elsif($rev_seq=~m/TAA$/ || $rev_seq=~m/TGA$/ || $rev_seq=~m/TAG$/){ #alternatively, if no start codon, but has a stop codon
					    if(length($rev_seq) % 3 ==0){ #if divisible by 3, extend in 5' direction to search for START codon
						$sstart+=3; #extend in 5' direction
						if($sstart > $CONTIG_MAX_LENGTH){ #if extension crosses contig bound
						    $sstart -=3; #revert back (-3)
						    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #Take $sstart -3 coordinate. Retrieve seq again as below frame-adjust loop also feeds into this loop.
						    $rev_seq=reverse($rev_seq); #reverse it
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						    $framecheck =1; # In this case, no stop will be found, exit extension and send to internal search.
						}
						elsif($sstart <= $CONTIG_MAX_LENGTH){ #if does not exceed contig bound, proceed with search for STOP codon
						    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1)));
						    $rev_seq=reverse($rev_seq);#reverse it
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						}
					    }
					    else{ #if the sequence is not divisible by 3 (this condition should never be met - here as buffer)
						$sstart+=1; #if the sequences is not divisible by 3, extend by one nucleotide until divisible by 3.
						if($sstart > $CONTIG_MAX_LENGTH){ #if exceeds contig bound
						    $sstart -=2; ##Go back two nucleotides (-1 site inwards).
						    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #take new sequence 
						    $rev_seq=reverse($rev_seq); #reverse it
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						    if(length($rev_seq) % 3 ==0){ #if divisible by three
							$framecheck =1; #send to internal search 
						    }
						    elsif(length($rev_seq) % 3 !=0){ #if still not divisible by 3
							$sstart -=1; #go back inwards one more step. Now will be divisible by 3.
							$rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #take new sequence 
							$rev_seq=reverse($rev_seq); #reverse it
							$rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
							$framecheck =1; #gets sent to internal search
						    }
						}
						elsif($sstart <= $CONTIG_MAX_LENGTH){ #if does not exceed contig bound, proceed with search for START codon
						    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1)));
						    $rev_seq=reverse($rev_seq);#reverse it
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						}
					    }
					}
					#############################
					# Search for START and STOP #
					#############################
					else{ #if the sequence has neither start or stop codon 
					    if(length($rev_seq) % 3 ==0){ #if divisible by 3
						$extend_checkpoint = 0;
						$send-=3;$sstart+=3; #extend in both directions (3' and 5' directions)
						if($sstart > $CONTIG_MAX_LENGTH){ #if crosses end of contig, go back to longest transcript in frame 
						    $sstart -=3; #revert back to max length in frame
						    $extend_checkpoint += 0.5; #for framecheck statement
						}
						if($send < 1){ #if cross contig bound, revert 
						    $send +=3; #revert back to max length in frame
						    $extend_checkpoint +=0.5; #for framecheck statement
						}
						$rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #take coordinates, (adjusted or reverted) 
						$rev_seq=reverse($rev_seq); #reverse it
						$rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						if ($extend_checkpoint ==1) { #if both ends cross the contig, stop extending in both directions and set framecheck to 1
						    $framecheck =1;
						} 
						elsif($extend_checkpoint <1){ #if at least one end does not exceed contig bound, continue extending in that direction
						}
					    }
					    elsif(length($rev_seq) % 3 !=0){ #if not divisible by 3 (should not meet this condition - included as buffer)
						$send-=1;$sstart+=1; #extend by 1 site until divisible by three
						if($sstart > $CONTIG_MAX_LENGTH && $send < 1){ #If crosses both contig bounds
						    $sstart-=1; #revert start coord
						    $send+=1; #revert end coord
						    $send +=1; #take one base from end since we don't have a stop
						    $rev_seq =uc(substr($refseq,($send -1), (($sstart-$send)+1))); #take sequence
						    $rev_seq=reverse($rev_seq); #reverse it
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						    if(length($rev_seq) % 3 ==0) { #if removing base results in divisible by 3, move on.
							$framecheck =1;
						    }
						    elsif(length($rev_seq) % 3 !=0){ #else not divisible by 3
							if ($rev_seq =~m/TAA$/ || $rev_seq=~m/TGA$/ || $rev_seq=~m/TAG$/){ #if seq ends in stop, start removing from beginning
							    $sstart -=1; 
							    $rev_seq =uc(substr($refseq,($send -1), (($sstart-$send)+1))); #take sequence
							    $rev_seq=reverse($rev_seq); #reverse it
							    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
							    if (length($rev_seq) % 3 ==0){ #now divisible by three
								$framecheck =1; #send to internal loops
							    }
							}
							else{ #if no STOP, remove one more nucelotide from the end (5') 
							    $send +=1; 
							    $rev_seq =uc(substr($refseq,($send -1), (($sstart-$send)+1))); #take sequence
							    $rev_seq=reverse($rev_seq); #reverse it
							    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
							    if (length($rev_seq) % 3 ==0){ #now divisible by three 
								$framecheck =1; #send to internal loops
							    }
							}   
						    }
						}
						elsif($sstart > $CONTIG_MAX_LENGTH && $send >=1){ #if crosses contig end (not start), adjust only in 3' direction 
						    $sstart -=1; #restore start 
						    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #take sequence 
						    $rev_seq=reverse($rev_seq); #reverse it
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						    if(length($rev_seq) % 3 ==0){ #if divisible by 3
							$framecheck =1; #send to internal loops
						    }
						}
						elsif($send <1 && $sstart <= $CONTIG_MAX_LENGTH){ #if crosses contig start (not end), adjust only in 5' direction
						    $send +=1; #restore end
						    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #take sequence 
						    $rev_seq=reverse($rev_seq); #reverse it
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						    if(length($rev_seq) % 3 ==0){ #if divisible by 3
							$framecheck =1; #send to internal search
						    }
						}
						elsif($send >=1 && $sstart <= $CONTIG_MAX_LENGTH){ #if does not exceed either contig bound
						    $sstart -=1; #adjust by 1 until divisible by 3
						    $rev_seq=uc(substr($refseq,($send -1), (($sstart-$send)+1))); #take sequence
						    $rev_seq=reverse($rev_seq); #reverse it
						    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						}
					    }
					}
					if(length($rev_seq)>$max_length && $framecheck==0){ #keep extending until START and STOP found OR exceeds user determined maximum length
					    $framecheck=1; #exit extension loops
					}
				    }
				    #################################
				    # BLASTBACK verification step:  #
				    # Prevents over extension       #
				    #################################
				    my $revert=0;
				    my @blastbackin = ();
				    push(@blastbackin, $rev_seq);
				    push(@blastbackin, $blastback_verification);
				    push(@blastbackin, $percentage_cover_threshold);
				    push(@blastbackin, $percentage_identity_threshold);
				    push(@blastbackin, $evalue_threshold);
				    if($blastback_verification ==1){
					my @blast_values=&blastback(@blastbackin);
					$revert = $blast_values[0];
					if ($revert ==1){ #If after extension loops, blast hit is too diluted 
					    $rev_seq=uc(substr($refseq,($bu_end_if-1), (($bu_start_if-$bu_end_if)+1))); #revert to backup coordinates and send to internal search
					    $rev_seq=reverse($rev_seq); 
					    $rev_seq=~tr/ATGCatgc/TACGtacg/;
					    $sstart = $bu_start_if; #adjust start coordinate accordingly
					    $send = $bu_end_if; #adjust end coordinate accordingly
					}
					else{ #if prediction is still within cover and identity thresholds, keep extension 
					}
				    }
				    ###########################
				    # Reverse internal search #
				    ###########################
				    if($rev_seq=~m/^ATG/ && $rev_seq=~m/TGA$/ || $rev_seq=~m/^ATG/ && $rev_seq=~m/TAG$/ || $rev_seq=~m/^ATG/ && $rev_seq=~m/TAA$/){ #if start and stop, don't enter internal search
				    }
				    elsif($rev_seq!~m/^ATG/ || $rev_seq!~m/TGA$/ || $rev_seq!~m/TAG$/ || $rev_seq!~m/TAA$/){ #if either no START codon or no STOP codon
					if ($rev_seq!~m/^ATG/ && $rev_seq!~m/TGA$/ &&  $rev_seq!~m/TAG$/ && $rev_seq!~m/TAA$/){ #if no start AND no stop codon
					    $rev_seq=uc(substr($refseq,($bu_end_if-1), (($bu_start_if-$bu_end_if)+1))); #revert to frame adjusted back up coordinates
					    $rev_seq=reverse($rev_seq); 
					    $rev_seq=~tr/ATGCatgc/TACGtacg/;
					    $sstart = $bu_start_if;
					    $send = $bu_end_if;
					}
					my $subseq=""; #look internally
					if(length($rev_seq) % 3 !=0){ #if not divisible by three, adjust to be divisible by 3 (condition should never be met - included as buffer)
					    my $newseq=uc(substr($rev_seq,0,((length($rev_seq))-1))); #trim end by 1
					    if (length($newseq) % 3 ==0){
						$rev_seq = $newseq; #if now divisible by 3, good to go
						$send +=1;
					    }
					    else{ #if still not divisible by 3 
						my $newseq=uc(substr($rev_seq,0,((length($rev_seq))-2))); #trim end by 1
						$rev_seq = $newseq; #now divisible by 3
						$send +=2;
					    }
					}
					if(length($rev_seq) % 3 ==0){ #is divisible by 3
					    my $Backup = $rev_seq;
					    my $startfound=0;
					    my $start_reduce =0; #0 + 3 will give 3rd index --> 4th position, we want 4th nucelotide to access second codon
					    if ($rev_seq =~ m/^ATG/){ #if START is found
						$subseq = $rev_seq;
						$startfound =1; #exit start trim loops
					    }
					    else{ #if no start codon 
						while($startfound ==0){ #trim internally in 3 -> 5' direction until START codon is found 
						    $start_reduce +=3; 
						    $sstart -=3;
						    $subseq = uc(substr($rev_seq,$start_reduce,length($rev_seq))); #get sequence
						    if ($subseq =~ m/^ATG/){ #if START is found
							$startfound =1; #exit start trim loops
						    }
						    if (length($subseq) < $min_length_internal){ #if internal trim is too short 
							$startfound =1; #exit start trim loops
						    }
						}
					    }
					}
					if (length($subseq) < $min_length_internal || $subseq !~ m/^ATG/){ #if too short OR start codon not found
					    $subseq =uc(substr($refseq,($send -1), (($bu_start_if - $send)+1))); #revert to in frame start and keep extended end
					    $subseq=reverse($subseq); #reverse
					    $subseq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
					    $sstart = $bu_start_if; #update coordinates accordingly
					}
					my $codonstring="";
					my $backup = $subseq;
					my $end_reduce = length($subseq);
					my $end_trim = 0; #trim internally from end coordinate in 3' -> 5' to search for STOP codon
					if ($subseq =~ m/TAG$/ || $subseq =~ m/TAA$/ || $subseq =~ m/TGA$/){ #if STOP codon is found
					    $codonstring = $subseq;
					    $end_trim = 1; #exit end trim loops
					}
					else{ #if no STOP codon
					    while($end_trim == 0){ #while no STOP codon 
						$end_reduce -= 3; #trim end by 3
						$send +=3; # adjust coordinate accordingly
						$codonstring= uc(substr($subseq,0,$end_reduce)); #get sequence
						if($codonstring =~ m/TAG$/ || $codonstring =~ m/TAA$/ || $codonstring =~ m/TGA$/){ #if STOP codon found
						    $end_trim = 1; #exit end trim loops
						}
						if (length($codonstring) < $min_length_internal){ #if internal trim is too short 
						    $end_trim = 1; #exit end trim loops 
						}
					    }
					}
					if($codonstring =~ m/TAG$/ || $codonstring =~ m/TAA$/ || $codonstring =~ m/TGA$/){ #if STOP codon is found
					    if(length($codonstring)<$min_length_internal){ #but sequence is too short
						$rev_seq = uc(substr($refseq,($bu_end_if -1),(($sstart - $bu_end_if)+1)));#revert end back to backup, but retain work done from START search 
						$rev_seq=reverse($rev_seq); #reverse
						$rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
						$send = $bu_end_if; #adjust end coordinate accordingly
					    }
					    else{
						$rev_seq = $codonstring; #if matches end AND is long enough, keep internal end search  
					    }
					}
					else{ #if no STOP codon found
					    $rev_seq =uc(substr($refseq,($bu_end_if-1), (($sstart-$bu_end_if)+1))); #revert to back up end coordinare regardless of length
					    $rev_seq=reverse($rev_seq); #reverse
					    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
					    $send = $bu_end_if; #adjust end coordinate accordingly
					}
				    }
				    ################################
				    # BLASTBACK verification step: #
				    # Prevents over trimming       #
				    ################################
				    $revert =0; my $COV=0; my $IDY=0; $EV=1; my $exclude =0;
				    @blastbackin = ();
				    push(@blastbackin, $rev_seq);
				    push(@blastbackin, $blastback_verification);
				    push(@blastbackin, $percentage_cover_threshold);
				    push(@blastbackin, $percentage_identity_threshold);
				    push(@blastbackin, $evalue_threshold);
				    if($blastback_verification ==1){
					my @blast_values=&blastback(@blastbackin);
					#my @blast_values=&blastback("$rev_seq");
					$revert = $blast_values[0]; #revert status
					$COV = $blast_values[1]; #percentage cover score
					$IDY = $blast_values[2]; #percentage identity score
					$EV = $blast_values[3]; #evalue
					$blast_program = $blast_values[4]; #megablast or nblast (if no hit with megablast, nblast is called)
					$exclude =0;
					if ($revert ==1){ #if prediction does not meet cover and identity thresholds after internal/extension loops
					    $rev_seq=uc(substr($refseq,($bu_end_if-1), (($bu_start_if-$bu_end_if)+1))); #revert to in-frame adjusted back up coordinates
					    $rev_seq=reverse($rev_seq); #reverse 
					    $rev_seq=~tr/ATGCatgc/TACGtacg/; #get reverse compliment
					    $sstart = $bu_start_if; #adjust start accordingly 
					    $send = $bu_end_if; #adjust end accordingly
					    my @blastbackin = ();
					    push(@blastbackin, $rev_seq);
					    push(@blastbackin, $blastback_verification);
					    push(@blastbackin, $percentage_cover_threshold);
					    push(@blastbackin, $percentage_identity_threshold);
					    push(@blastbackin, $evalue_threshold);
					    @blast_values=&blastback(@blastbackin);
					   # @blast_values=&blastback("$rev_seq"); #do one more blastback step, with backup coorinates
					    $revert = $blast_values[0]; #if revert is 1 at this point, sequence will be excluded
					    $COV = $blast_values[1]; #new percentage cover score
					    $IDY = $blast_values[2]; #new percentage identity score
					    $EV = $blast_values[3]; #new evalue
					    $blast_program = $blast_values[4]; #megablast or nblast
					    if ($revert ==1){ #if revert is still 1
						$exclude =1; #exclude sequence
					    }
					}
					$COVER_rounded = sprintf("%.4f", $COV); #round percentage cover score to 4 decimals for annotation 
					$ID_rounded = sprintf("%.4f", $IDY); #round percentage identity score to 4 decimals for annotation 
				    }
				}
				##########################################################
				#Hmm filter - remove prediction if no 7transmembrane hit #
				##########################################################
				my $hmmremove = "no";
				if($HMM_filter eq "Yes" || $HMM_filter eq "yes"){
				    if($rev_seq ne "NOPREDICTION"){
					my $Protein_prediction = ">Prediction_in\n";
					for (my $d=0;$d<(length($rev_seq) -2);$d+=3)
					{
					    my $Qcodon = substr($rev_seq,$d,3);
					    $Protein_prediction.=&dna2prot($Qcodon); #translate reference receptor using dna2protein
					}
					my $prediction_out =  $randomstr."_hmm_input_protein.fa";
					open(PPOUT, ">$prediction_out");
					print PPOUT $Protein_prediction;
					close PPOUT;
					my @hmm_inputs = ();
					push(@hmm_inputs, $classA_HMM);
					push(@hmm_inputs, $prediction_out);
					push(@hmm_inputs, $randomstr);
					my $hmm_result =&hmmsearch(@hmm_inputs);
					if($hmm_result == 0){
					    my @hmm_inputs = ();
					    push (@hmm_inputs, $classC_hmm);
					    push (@hmm_inputs, $prediction_out);
					    push(@hmm_inputs, $randomstr);
					    my $hmm_result =&hmmsearch(@hmm_inputs);
					    if($hmm_result == 0){
                                                my @hmm_inputs = ();
                                                push (@hmm_inputs, $vom2r_hmm);
                                                push (@hmm_inputs, $prediction_out);
                                                push(@hmm_inputs, $randomstr);
						my $hmm_result =&hmmsearch(@hmm_inputs);
                                                if($hmm_result == 0){
                                                    $hmmremove = "yes";
						}
					    }
					}
				    }   
				}
				#########################################################
				#Blast prediction back to reference file to re-annotate #
				#########################################################
				if($rev_seq ne "NOPREDICTION"){
				    open(mOUT, ">$multiout");
				    print mOUT ">multi_out\n".$rev_seq;
				    close mOUT;
				    `blastn -task 'megablast' -db $reff_db -query $multiout -out $multiblatout -num_threads $threads`;
				    my @qhits =&parse_blast_hits($multiblatout);
				    my $thit = $qhits[0];
				    if ($thit){
					if($thit =~ m/\_(O[A-Za-z0-9]+)/ ||
                                           $thit=~m/\_(V[A-Za-z0-9]+)/ ||
                                           $thit=~m/\_(T[A-Za-z0-9]+)/ ||
                                           $thit=~m/\_(RHO)/ ||
                                           $thit=~m/\_(rho)/i ||
                                           $thit=~ m/\_(rhol)/i ||
                                           $thit =~ m/\_(Taar[A-Za-z0-9]+)/i ||
                                           $thit =~ m/\_(TAS[0-9A-Za-z]+)/i ||
                                           $thit =~ m/\_(Vom[0-9A-Za-z]+)/i ||
                                           $thit =~ m/\_(VN[\S]+)/i ||
                                           $thit =~ m/\_(V2[\S]+)/i ||
                                           $thit =~ m/\_(Olf[0-9A-Za-z]+)/i ||
                                           $thit =~ m/\_(OPN[0-9A-Za-z]+)/i ||
                                           $thit =~ m/\_(OR[0-9A-Za-z]+)/i ||
                                           $thit =~ m/\_(Rh[0-9])/i){
					    my $name = $1;
					    $name =~ s/\_//g;
					    $receptor_name = $name;
					}
				    }			    
				}
				#################################################
				# checkframe subroutine to annotate prediction  #
				#################################################
				my @status=&checkframe($rev_seq); #use checkframe subroutineto get appropriate annotation status 
				my $stat=$status[0]; #get annotation status
				if($truncated == 1){ #if truncated option is ON
				    if($rev_seq =~ m/TAG$/ || $rev_seq =~ m/TAA$/ || $rev_seq =~ m/TGA$/){ #if STOP codon
					if($include_score ==1 && $tool eq "sensommatic"){ #if blastback annotation score option is ON
					    $tar_gene=">Seq".$ref_name."_".$receptor_name."_status_rev_".$send."_".$sstart."|".$COVER_rounded."|".$ID_rounded."|".$EV."|".$blast_program."|".$tool."\n".uc($rev_seq)."\n"; #add blastback scores to annotation
					}
					elsif($include_score==0 || $tool eq "augustus"){ #if blastback annotation score option is OFF
					    $tar_gene=">Seq".$ref_name."_".$receptor_name."_status_rev_".$send."_".$sstart."|".$tool."\n".uc($rev_seq)."\n"; #dont include blastback scores to annotation
					}
				    }
				    else{
					if($include_score ==1 && $tool eq "sensommatic"){ #if blastback annotation score option is ON
					   $tar_gene=">Seq".$ref_name."_".$receptor_name."_status_".$annotation_truncated."_rev_".$send."_".$sstart."|".$COVER_rounded."|".$ID_rounded."|".$EV."|".$blast_program."|".$tool."\n".uc($rev_seq)."\n"; #include blastback scores to annotation
					}
					elsif($include_score==0 || $tool eq "augustus"){ #if blastback annotation score option is OFF
					    $tar_gene=">Seq".$ref_name."_".$receptor_name."_status_".$annotation_truncated."_rev_".$send."_".$sstart."|".$tool."\n".uc($rev_seq)."\n"; #dont include blastback scores to annotation
					}
				    }
				}
				else{
				    if($include_score ==1 && $tool eq "sensommatic"){ #if blastback annotation score option is ON
					$tar_gene=">Seq".$ref_name."_".$receptor_name."_status_rev_".$send."_".$sstart."|".$COVER_rounded."|".$ID_rounded."|".$EV."|".$blast_program."|".$tool."\n".uc($rev_seq)."\n"; #include blastback scores to annotation
				    }
				    elsif($include_score ==0 || $tool eq "augustus"){ #if blastback annotation score option is OFF
					$tar_gene=">Seq".$ref_name."_".$receptor_name."_status_rev_".$send."_".$sstart."|".$tool."\n".uc($rev_seq)."\n"; #don't include blastback scores to annotation
				    }
				}
				my $length_check = 0;
				$targene_len=length($rev_seq); #get length of prediction
				if ($targene_len < $pseudogene_length){ #if length is below user determined threshold, annotate as pseudogene
				    $tar_gene=~s/status/$annotation_short/; #If prediction is shorter than $pseudogene_length, gene will be annotated as pseudogene   
                                    $length_check = 1;
                                }
				if($stat eq "1"){
				    unless($length_check ==1){
					$tar_gene=~s/status/$annotation_1/; #START codon && no in frame stop codons
				    }
				}
				if($stat eq "2"){
				    unless($length_check ==1){
					$tar_gene=~s/status/$annotation_2/; #no START codon && no stop codons in any frame
				    }
				}
				if($stat eq "3"){
				    unless($length_check ==1){
					$tar_gene=~s/status/$annotation_3/; #no START codon && no in frame stop codons
				    }
				}
				if($stat eq "4"){
				    if($length_check ==1){
					my $new_annotation = $annotation_4."_short";
					$tar_gene =~ s/$annotation_short/$new_annotation/;
				    }
				    else{
					$tar_gene=~s/status/$annotation_4/; #START codon && in frame stop codon
				    }
				}
				if($stat eq "5"){
				    if($length_check ==1){
					my $new_annotation = $annotation_5."_short";
					$tar_gene =~ s/$annotation_short/$new_annotation/;
				    }
				    else{
					$tar_gene=~s/status/$annotation_5/; #START codon && in frame stop codon                                                                                                                                                             
				    }
				}
				if($stat eq "6"){
				    if($length_check ==1){
					my $new_annotation = $annotation_6."_short";
					$tar_gene =~ s/$annotation_short/$new_annotation/;
				    }
				    else{
					$tar_gene=~s/status/$annotation_6/; #START codon && in frame stop codon                                                                                                                                                             
				    }
				}
				if($stat eq "7"){
				    $tar_gene=~s/\>Seq/>Ignore/;
				}
				if($targene_len>=$minlength){ #finally, if prediction is greater than user determined threshold
				    if($hmmremove ne "yes"){
					push @finalarray,$tar_gene; #include prediction in final output                                                                                          
				    }
				    else{
					if($tar_gene =~ m/([^\n]+)\n([\S]+)/){
                                            my $hmmheader = $1;
                                            my $hmmseq = $2;
                                            $hmmheader .="_hmm_remove";
                                            my $hmm_targene = $hmmheader."\n".$hmmseq;
                                            push @finalarray, $hmm_targene;
					}
				    }
				}
				else{#otherwise, prediction will be excluded
				    push @finalarray, "NOPREDICTION";
				}
				$targene_len=0;
			    }
			    #If contracted, take the blat upper and lower limits to blank out. Avoid genes being split into two predictions.
			    if ($sstart > $send){
				#rev
				if ($sstart < $bu_sstart){ #start contracted inwards 
				    $block_start = $bu_sstart; #block the upper limit
				}
				else{
				    $block_start = $sstart;
				}
				if($send > $bu_send){ #end contracted inwards     
				    $block_end = $bu_send; #block the lower limit
				}
				else{
				    $block_end = $send;
				}
				push(@finalarray, $block_end); #end first since rev
				push(@finalarray, $block_start); #start second since rev
				
			    }
			    elsif($sstart < $send){
				#pos
				if ($sstart > $bu_sstart){
				    $block_start = $bu_sstart;
				}
				else{
				    $block_start = $sstart;
				}
				if($send < $bu_send){
				    $block_end = $bu_send;
				}
				else{
				    $block_end = $send;
				}
				push(@finalarray, $block_start); #start first since pos
				push(@finalarray, $block_end); #end second since pos
			    }
			    push(@finalarray, $CONTIG_MAX_LENGTH);
			    push(@finalarray, $receptor_name);
			}
		    }
		}
	    }
	} 
    }
    return @finalarray; #retrun prediction to main code
}


#Firstsweep subroutine: Order reference receptors by signifigance:
#Sort each hit for each query by evalue and identity, so that each reference receptor can be stored twice, for two hits with different evals ..eg:
#QueryA, 0.0
#QueryB, 0.0001
#QueryA, 0.001,
#QueryC, 0.1
#hence Query A maps to two regions of the contig, but each hit has different e-vals. 
#If QueryB maps to the same region as second QueryA hit, but is more signifigant, we want Query B to predict the hit, not QueryA. 
sub firstsweep{#subroutine to take all sequences with a blast hit
    my @input_sweep = @_;
    my $blatoutputfile = $input_sweep[0];
    my $cover_thresh = $input_sweep[1];
    my $qcover = $input_sweep[2];
    my @finalseqs=();
    my %hash_query_evals;
    my $query_r;
    my $percent_id; my $inverse_percent;
    my $p = 0;
    my @Hits =&parse_blastn($blatoutputfile);
    my $multi = 0;
    foreach my $hit(@Hits){
	$multi = 0;
	$p +=1;
	my @val =&get_blastn_info($hit);
	my $name = $val[0];
	my $cover = $val[10];
	$cover = sprintf("%.0f", $cover);
	my $percent_id = $val[11];
	$percent_id = sprintf("%.0f", $percent_id);
	my $Evalue = $val[8];
	my $query_receptor = $val[0];
	my $inverse_percent = 100 - $percent_id; #find 100 - identity, so that identity is ordered correctly
	if($inverse_percent < 10){ #to ensure numeric sorting
	    $inverse_percent = "0".$inverse_percent;
	}
	if($name =~ m/OPN/i ||
	   $name =~ m/VN2/i || 
	   $name =~ m/Vom2/i || 
	   $name =~ m/Vmn2/i || 
	   $name =~ m/V2R/i || 
	   $name =~ m/TAS1R/i || 
	   $name =~ m/T1R/i ||
	   $name =~ m/OPN1S/i ||
	   $name =~ m/OPN1M/i ||
	   $name =~ m/OPN1L/i ||
	   $name =~ m/OPN2/i ||
	   $name =~ m/RHO/i ||
	   $name =~ m/Rho/i){
	    $multi = 1;
	}
	else{
	    $multi = 0;
	}
	if($multi == 1){
	    if ($cover_thresh > 200){
		$cover_threshold = 200;
	    }
	    else{
		$cover_threshold = $cover_thresh;
	    }
	}
	else{
	    $cover_threshold = $cover_thresh;
	}
	if($cover>$cover_threshold){#if meets cover threshold
	    $query_r = $inverse_percent."|".$query_receptor."|".$p; #store header with identity appended
	    $hash_query_evals{$query_r} = $Evalue; #store header as key and evalue as value in hash as pair
	}
    }
    #sort reference receptors by 1) increasing evalue and 2) decreasing percent query cover and 3) alphabetically for reproducability
    if(%hash_query_evals){
	foreach my $element(sort { $hash_query_evals{$a} <=> $hash_query_evals{$b} or $a cmp $b } keys %hash_query_evals){
	    $element =~ s/\|[0-9]+//g; #remove appended number
	    $element =~ s/[0-9]+\|//g; #remove number identity
	    push(@finalseqs, $element);
	}
    }
    return @finalseqs; #return array with reference receptors ordered by signifigance
}


#getframe subroutine:
#Returns the frame for extension loops
#translates prediction in three frames, and uses blastp with translated reference receptor to find best frame
sub getframe_tblastn{
    my $F1=$_[0];
    my $reference_receptor = "";
    my $header = "";
    my $Reference = "";
    my $Query_Protein = "";
    open(OUT, ">$Tblastn_set");
    print OUT ">QUERY\n$F1";
    close OUT;
    open(RR, "$tmp_seq"); #open reference receptor file
    {
        local $/; #read in seq as one chunk
        $reference_receptor = <RR>; #store reference receptor sequence
    }
    close RR;
    if($reference_receptor=~m/(\>[\S]+)\n([\S]+)/){
        $header = $1;
    }
    $reference_receptor =~ s/$header//; #remove header from sequence
    $reference_receptor =~s/[\n\r]//g; #remove line breaks for translation
    $Reference = $reference_receptor; #reference now just a sequence to be translated
    for (my $d=0;$d<(length($Reference) -2);$d+=3)
    {
        my $Qcodon = substr($Reference,$d,3);
        $Query_Protein.=&dna2prot($Qcodon); #translate reference receptor using dna2protein
    }
    open(OUT, ">$prot_ref"); #print out translated reference receptor to file                                                                                                               
    print OUT ">Reference"."\n".$Query_Protein;
    close OUT;
    `makeblastdb -in $Tblastn_set -dbtype nucl -out $Tblastn_db`; #make database with translated predictions in three frames   
    `tblastn -task 'tblastn-fast' -db $Tblastn_db -query $prot_ref -out $tblastn_out -num_threads $threads`;#blast the reference receptor to the target contig
    my @tblastn_hits =&parse_blast_hits($tblastn_out);
    my $top_hit = $tblastn_hits[0];
    my $frame = "";
    if($top_hit){
	my @query_hits =&parse_blast_subhits($top_hit);
	my $top_hit_alignment = $query_hits[1];
	my @blast_info =&get_blast_info($top_hit_alignment);
	my $tblastn_start = $blast_info[5];
	my $tblastn_end = $blast_info[6];
	my $frame_tblastn = $blast_info[8];
	if($frame_tblastn =~m/\+1/){
	    $frame = "Frame1";
	}
	elsif($frame_tblastn =~m/\+2/){
	    $frame = "Frame2";
	}
	elsif($frame_tblastn =~m/\+3/){
	    $frame = "Frame3";
	}
	return $frame;
    }
    else{
	$frame = "Frame1";
	return $frame;
    }
}

#parse blast output format 6 
sub parse_blastn{  
    my $tsv_file = $_[0];
    open(BLSTN, $tsv_file);
    my @HSPs = (<BLSTN>);
    return @HSPs;
    close BLSTN;
}

#Get blast information (blast output format 6)
sub get_blastn_info{
    my $HSP = $_[0];
    my @vals = split(/\t/, $HSP);
    return @vals;
}

#read in blast out file and output hits per subject
sub parse_blast_hits{ 
    my $blastoutputfile=$_[0];
    my @blast_results = ();
    my @contig_hits = ();
    open(BLAST, $blastoutputfile);
    {
	local $/ = "Query\=";
	while(<BLAST>){
	    my $qc = $_;
	    if($qc =~ m/Query\=/){
		$qc =~ s/Query\=//g;
		$qc = "Query\=".$qc;
	    }
	    else{
		$qc = "Query\=".$qc;
	    }
	    push (@blast_results, $qc);
	}
    }
    close BLAST;
    shift(@blast_results);
    my $query_results = $blast_results[0];
    @contig_hits = split (/\>/, $query_results);
    shift @contig_hits;
    return @contig_hits;
}


#Split each blast hit by "score" for local alignments
sub parse_blast_subhits{ #read in array of > hits and split each by score
    my $query_contig = $_[0];
    my @score_results = ();
    if($query_contig){
	@score_results = split (/Score/, $query_contig);
    }
    return @score_results;
}


#Get blast info: read in blast alignment and output information
sub get_blast_info{ 
    my $alignment = $_[0];
    my $copy = $alignment;
    my $cover = "";
    my $identity = "";
    my $evalue = "";
    my $subj_length = "";
    my $subj_start = "";
    my $subj_end = "";
    my $query_start = "";
    my $query_end = "";
    my $tblastn_frame = "";
    $alignment = "Score".$alignment;
    my @alignment_info = ();
    if($alignment=~m/Identities[\s]+\=[\s]+([0-9]+)\/[0-9]+[\s]+\(([0-9]+)\%/){
	$cover=$1;
	$identity=$2;
    }
    if($alignment =~ m/Expect[\s]+\=[\s]+([0-9]+.*[0-9]+)/){
	$evalue= $1;
    }
    if($alignment =~m/Length\=([0-9]+)/){ 
	$subj_length = $1;
    }
    if($alignment=~m/Query[\s]+([0-9]+)/){
	$query_start=$1;
    }
    if($alignment=~m/Sbjct[\s]+([0-9]+)/){
	$subj_start=$1;
    }   
    while($copy=~s/Query[\s]+[0-9]+[\s]+[\S]+[\s]+([0-9]+)//){
	$query_end=$1;
    }
    while($copy=~s/Sbjct[\s]+[0-9]+[\s]+[\S]+[\s]+([0-9]+)//){
	$subj_end=$1;
    }
    if($alignment =~ m/Frame[\s]\=[\s]([^\s][0-3])/){
	$tblastn_frame = $1;
    }
    push (@alignment_info, $cover); #................cover: 0
    push (@alignment_info, $identity); #..........identity: 1
    push (@alignment_info, $evalue); #..............evalue: 2
    push (@alignment_info, $query_start); #....Query start: 3
    push (@alignment_info, $query_end); #........Query end: 4
    push (@alignment_info, $subj_start); #...Subject start: 5
    push (@alignment_info, $subj_end); #.......Subject end: 6
    push (@alignment_info, $subj_length); #..contig length: 7
    push (@alignment_info, $tblastn_frame); #........frame: 8
    return @alignment_info;
}


#checkframe subroutine: Checks for in frame stop codons to annotate functional vs pseudogene status
sub checkframe{
    my $status=0; #holds annotation status 
    my $dnaseqf1=uc($_[0]); #takes prediction as input (frame 1)
    my @framedata=();
    my $ignore_seq=0;
    if($dnaseqf1=~m/NOHITS/){ 
	$ignore_seq++;
	$dnaseqf1=~s/NOHITS//g;
    }
    if($dnaseqf1=~m/TGA$/){
	$dnaseqf1=~s/TGA$//; #remove stop codon
    }
    elsif($dnaseqf1=~m/TAA$/){
	$dnaseqf1=~s/TAA$//; #remove stop codon
    }
    elsif($dnaseqf1=~m/TAG$/){
	$dnaseqf1=~s/TAG$//; #remove stop codon
    }
    else{ #if no stop codon do nothing
    }
    my $dnaseqf2=substr($dnaseqf1,1,(length($dnaseqf1)-1)); #frame2
    my $dnaseqf3=substr($dnaseqf1,2,(length($dnaseqf1)-2)); #frame3
    my $f1stop=0; my $f2stop=0; my $f3stop=0; #declare variables to track stop codons in each frame
    my @framecheck=();
    $framecheck[0]=$dnaseqf1;$framecheck[1]=$dnaseqf2;$framecheck[2]=$dnaseqf3; #store counts in @framecheck array
    for(my $k=0;$k<scalar(@framecheck);$k++){
	for(my $n=0;$n<length($framecheck[$k])-2;$n+=3){ #read in codons for each frame
	    my $codon=substr($framecheck[$k],$n,3);
	    if($codon eq "TGA" || $codon eq "TAA" || $codon eq "TAG"){
		if($k==0){
		    $f1stop++; #stop codon in frame 1
		}
		elsif($k==1){
		    $f2stop++; #stop codon in frame 2
		}
		elsif($k==2){
		    $f3stop++; #stop codon in frame 3
		}
	    }
	}
    }
    if($dnaseqf1=~m/^ATG/){ #if start codon
	if($f1stop==0){
	    $status=1; ##1 means start codon, with no in frame stop codons
	}
	else{
	    $status=4; #4 means start codon, but contains at least 1 in frame stop codon
	}
    }
    else{ #no start codon
	if($f1stop==0 && $f2stop==0 && $f3stop==0){
	    $status=2; #2 means no start but no stop codons in any frame
	}
	elsif($f1stop>0 && $f2stop>0 && $f3stop>0){
	    $status=5; #5 means no start and stop codons in in all frames 
	}
	elsif($f1stop>0){ #at least frame1 has stop codon
	    $status =6; #6 means no START codon and in frame stop codon
	}
	else{ #only frames 2 and 3 have stop codon
	    $status=3; #3 means no START codon and no in frame stop codon (stop codon in frame 2 and 3) 
	}
    }
    if($ignore_seq>0){
	$status=7; #no hit
    }
    push @framedata,$status; #push status to @framedata array
    return @framedata; #subroutine returns this array with annotation status
}


#dna2prot subroutine: translation: nucelotide sequence -> amino acid sequence 
sub dna2prot{
    my $codon=$_[0]; #takes codon as input
    $codon=uc $codon;
    my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G','NNN'=>'X','GCN'=>'A','RAY'=>'B','TGY'=>'C','GAY'=>'D','GAR'=>'E','TTY'=>'F','GGN'=>'G','CAY'=>'H','ATH'=>'I','AAR'=>'K','TTR'=>'L','CTN'=>'L','YTR'=>'L','AAY'=>'N','CCN'=>'P','CAR'=>'Q','CGN'=>'R','AGR'=>'R','MGR'=>'R','TCN'=>'S','AGY'=>'S','ACN'=>'T','GUN'=>'V','TAY'=>'Y','SAR'=>'Z','TAR'=>'*','TRA'=>'*');
    if(exists $g{$codon})
    {
	return $g{$codon}; #translate codon and return amino acid
    }
    else
    {
	return "X";
    }
}

#blastback subroutine:
#Blasts prediction back to reference to ensure over-extension or over-trimming does not occur. 
#Returns array with scores from blastback step
sub blastback{ #verification subroutine #Is prediction actually a receptor? If <X% of prediction macthes query with BLAST, revert to original blast hit (inframe coordinates)
    my @invals = @_;
    my $prediction = $invals[0];
    my $blastback_verification = $invals[1];
    my $percentage_cover_threshold = $invals[2];
    my $percentage_identity_threshold =$invals[3];
    my $evalue_threshold = $invals[4];
    my @blast_back_values = (); #array to hold revert status and blast scores 
    my $BB; 
    my $bb_identity_nt=0;
    my $identity_ratio=0;
    my $revert=0;
    my $bb_cover=0;
    my $PRED_length=0;
    my $cover_ratio=0;
    my $bb_identity=0;
    my $bb_evalue=1;
    my $prediction_file=$randomstr."_prediction_seq.fa"; #prediction sequence file
    open(Prediction_OUT, ">$prediction_file"); #open new file
    print Prediction_OUT ">Prediction\n".$prediction."\n"; #print out prediction sequence
    close Prediction_OUT; #close prediction file
    my $reference_receptor_DB=$randomstr."_refRDB"; #blast database (preediction is database)
    my $blastback_out=$randomstr."_blastback_out";#blast output file
    `makeblastdb -in $prediction_file -dbtype="nucl" -out $reference_receptor_DB`; #use prediction as database
    `blastn -task 'dc-megablast' -db $reference_receptor_DB -query $tmp_seq -out $blastback_out -num_threads $threads`; #blast reference receptor against prediction
    if(-e $prediction_file){ #if prediction file exists
	`rm $prediction_file`; #remove for next round
    }
    open(BBIN, "$blastback_out"); #open the blast output file
    {
	local $/; #set delimiter to nothing - allows file to be read as one chunk
	$BB = <BBIN>; #store blast output file in $BB
    }
    close BBIN; #close file
    my @BBOUT = split(/\>/,$BB); #Split at each query
    my $discard = $BBOUT[0]; #discard rubbish stored in element 0 of array
    my $bbouthit = $BBOUT[1]; #store element 1 as tophit
    my $blast_type = "mb"; #megablast used
    if ($discard =~ m/[\S]+\*\*\*[\S]+/){ #if first element matches no hit
	`blastn -db $reference_receptor_DB -query $tmp_seq -out $blastback_out -num_threads $threads`; #try using nblast
	open(BBIN, "$blastback_out"); #open nblast out file
	{ 
	    local $/; #set delimiter to nothing - allows file to be read as one chunk 
	    $BB = <BBIN>; #store blast output file in $BB
	} 
	close BBIN; #close file
	@BBOUT = split(/\>/,$BB); #split at each query
	$discard = $BBOUT[0]; #discard rubbish stored in element 0 of array
	$bbouthit = $BBOUT[1]; #store element 1 as tophit
	$blast_type = "nb"; #nblast used
	
    }
    if($discard =~ m/[\S]+\*\*\*[\S]+/){ #if discard still matches no hits - exclude
    }
    else{
	if ($bbouthit =~ m/Identities[\s]+\=[\s]+([0-9]+)\/([0-9]+)[\s]+\(([0-9]+).*\)/){ 
	    $bb_identity = $3; #store identity percentage
	    $bb_cover = $2; #store cover value (number of nucleotides aligned to the reference receptor)
	    $bb_identity_nt = $1; #store identity value (number of nucleotides identical to the reference receptor)
	}
	if ($bbouthit =~ m/Expect[\s]+\=[\s]+([0-9]+.*[0-9]+)/){
	    $bb_evalue = $1; #store evalue
	}
    }
    $PRED_length = length($prediction); #get length of prediction
    $cover_ratio = $bb_cover / $PRED_length; #number of nucleotides aligned with reference divided by total prediction length: (value usually between 0 and 1 -> can exceed 1 if the blast alignments introduces gaps)
    $identity_ratio = $bb_identity_nt / $PRED_length; #number of identical nucelotides aligned with reference divided by total prediction length: (value between 0 and 1)
    if ($cover_ratio < $percentage_cover_threshold || $identity_ratio < $percentage_identity_threshold|| $bb_evalue > $evalue_threshold ){ #if hit does not exceed user defined thresholds
	$revert = 1; #revert =1 (sequence will be reverted back to backup coords in blatout routine)
	push (@blast_back_values, $revert); #element 0 of retruned array
    }
    else{ #if thresholds are met
	$revert = 0; #don't revert
	push (@blast_back_values, $revert); #element 0 of returned array
    }
    push (@blast_back_values, $cover_ratio); #element 1 of returned array
    push (@blast_back_values, $identity_ratio); #element 2 of returned array
    push (@blast_back_values, $bb_evalue); #element 3 of returned array 
    push (@blast_back_values, $blast_type); #element 4 of returned array
    return @blast_back_values; #return array with revert status and blast scores
}


#Run augustus and return prediction and prediction info
sub augustus_predict{
    my $seq_info_in = $_[0];
    my @splits = split(/\|/, $seq_info_in);
    my $sequence = $splits[0];
    my $seq_start = $splits[1];
    my $seq_end = $splits[2];
    my $reference_file = $splits[3];
    my $augustus_species = $splits[4];
    my $multiexon_option = $splits[5];
    my $random_str = $splits[6];
    my $opsin_pred = $splits[7];
    my $opsin_reference = $splits[8];
    my $header = "";
    my $seq = "";
    my $count = 0;
    my @preds = ();
    my $tmp_out = $random_str."_tmp.fa";
    open(OUT, ">$tmp_out");
    print OUT $sequence;
    close OUT;
    my $minidentity = 60;
    my $psl =  $random_str."_ref.psl";
    my $hints =  $random_str."_hints.gff";
    my $reference_outseq =  $random_str."_ref.fa";
    my $prediction_gff =  $random_str."_prediction_out.gff";
    if($opsin_pred == 1){
	`blat -minIdentity=$minidentity $tmp_out $opsin_reference $psl`;
    }
    else{
	`blat -minIdentity=$minidentity $tmp_out $reference_file $psl`;
    }
    `perl blat2hints.pl --in=$psl --out=$hints`;
    if($multiexon_option == 0){ #single exon, genemodel==intronless
	system("augustus --species=$augustus_species --strand=forward --genemodel=intronless --codingseq=on --softmasking=0 --hintsfile=$hints --extrinsicCfgFile=extrinsic.ME.cfg $tmp_out > $prediction_gff");
    }
    elsif($multiexon_option ==1){ #multiexon
	system("augustus --species=$augustus_species --strand=forward --codingseq=on --softmasking=0 --hintsfile=$hints --extrinsicCfgFile=extrinsic.ME.cfg $tmp_out > $prediction_gff")
    }
    my $prediction_in = "";
    open(AUG, $prediction_gff);
    {
	local $/;
	$prediction_in = <AUG>;
    }
    close AUG;
    my @intron_coords = ();
    my @exon_coords = ();
    my @predictions = split(/#[\s]start[\s]gene/, $prediction_in);
    shift @predictions;
    foreach my $pred(@predictions){
	my @pred_split = split (/#[\s]protein[\s]sequence/, $pred);
	my $gene_cds_details = $pred_split[0];
	my $protein_details = $pred_split[1];
	my @split_again = split(/#[\s]coding[\s]sequence/, $gene_cds_details);
	my $gene_details = $split_again[0];
	my $cds_details = $split_again[1];
	my $pred_start = "";
	my $pred_end = "";
	$gene_details =~ s/[\s]+/\|/g;
	my $copy = $gene_details;
	if($gene_details =~ m/\|CDS\|([0-9]+)\|([0-9]+)/i){
	    $pred_start = $1;
	}
	my @end = ();
	while($copy =~s/\|CDS\|([0-9]+)\|([0-9]+)/hello/){
	    $pred_end = $2;
	    push(@end, $pred_end);
	}
	if(@end){
	    $pred_end = pop(@end);
	}
	my $new_copy = $gene_details;
	my $stop_codon_start = "";
	my $stop_codon_end = "";
	if($gene_details =~ m/\|stop\_codon\|([0-9]+)\|([0-9]+)/i){
	    $stop_codon_start = $1;
	}
	while($new_copy =~ s/\|stop\_codon\|([0-9]+)\|([0-9]+)/hello/){
	    $stop_codon_end = $2;
	}
	if($stop_codon_end){
	    if($stop_codon_end > $pred_end){
		$pred_end = $stop_codon_end;
	    }
	}
	my $maximum_end = "";
	my $minimum_start = "";
	my $prediction_found = "";
	if($pred_end && $pred_start){
	    if($pred_end > $seq_end){
		$maximum_end = $pred_end;
	    }else{
		$maximum_end = $seq_end;
	    }
	    if($pred_start < $seq_start){
		$minimum_start = $pred_start;
	    }else{
		$minimum_start = $seq_start;
	    }
	    my $prediction_length = ($pred_end - $pred_start) +1;
	    my $seq_length = ($seq_end - $seq_start) +1;
	    my $total_length = $prediction_length + $seq_length;
	    my $max_span = ($maximum_end - $minimum_start) +1;
	    if ($max_span < $total_length){
		$prediction_found = 1;
		if ($cds_details =~ m/\[([^\]]+)\]/){
		    my $cds_seq = $1;
		    $cds_seq =~ s/\#//g;
		    $cds_seq =~ s/\s//g;
		    $cds_seq = uc($cds_seq);
		    $cds_seq = $cds_seq."|".$pred_start."|".$pred_end;
		    push @preds, $cds_seq;
		}
	    }
	    else{
		$prediction_found = 0;
	    }
	}
	else{
	    $prediction_found = 0;
	}
    }
    return @preds;
}


#Parse fasta file
sub parse_fasta{ #returns sequences stored in an array     
    my $sequence_file = $_[0];
    my @sequence_array = ();
    open(SEQS, "$sequence_file");
    {
        local $/ = ">";
        while(<SEQS>){
            my $gene_sequence = $_;
            if($gene_sequence =~ m/\>/){
                $gene_sequence =~ s/\>//g;
                $gene_sequence = ">".$gene_sequence;
            }
            else{
                $gene_sequence = ">".$gene_sequence;
            }
            push (@sequence_array, $gene_sequence);
        }
    }
    close SEQS;
    shift @sequence_array;
    return @sequence_array;
}

#Convert fasta array to key
sub gene_key{ #takes seq array and outputs hash where headers are key and seqs are values     
    my @genes = @_;
    my %seq_key;
    foreach my $gene(@genes){
        if($gene=~m/(\>[\S]+)\n([\S\n]+)/){
            my $header = $1;
            my $seq = $2;
            $seq =~ s/\n//g;
            $header =~ s/\>//g;
            $seq_key{$header} = $seq;
        }
    }
    return %seq_key;
}


#Pseudogenome subroutine:
#Blast reference receptors against genome and write out contigs to pseudogenome file
sub pseudogenome{
    my @input_vals = @_;
    my $genome_assembly = $input_vals[0];
    my $queries_file = $input_vals[1];
    my $rand = $input_vals[2];
    my $threads = $input_vals[3];
    my $pseudogenome_perl = $input_vals[4];
    my $pseudogenome_esl = $input_vals[5];
    my @contigs = ();
    my $nblast_contig_out=$randomstr."_contig_blast_out";              
    my $pseudogenome=$randomstr."_pseudogenome";
    my $genomedb = $randomstr."_genomeDB";
    my $PseudoGenome = "";
    `makeblastdb -in $genome_assembly -dbtype="nucl" -out $genomedb`; #make blastdb from the target contig     
    `blastn -task 'megablast' -db $genomedb -query $queries_file -out $nblast_contig_out -num_threads $threads -qcov_hsp_perc 15`;#blast all sensory genes against target contig
    open(BLAT, "$nblast_contig_out");
    my $Contig_Hit="";
    my @contigblast = ();
    {
	local $/ = "Query\=";
	while(<BLAT>){
	    my $qc = $_;
	    $qc =~ s/Query\=//g;
	    push (@contigblast, $qc);
	}
    }
    foreach my $query_chunk(@contigblast){#loop through each result
	my $query_chunk = "Query\=".$query_chunk;
	my @hits=split(">",$query_chunk);
	foreach my $hit(@hits){
	    unless($hit =~ m/BLASTN\s2\.9\.0\+/ ||  $hit =~ m/Query\=.*/i) { #Only reappend the > to contig IDs
		$hit = ">".$hit; #Reappend the fasta header to contig chunk
	    }
	    if ($hit =~ m/(\>.*?\s)/i){ #if contig has a hit
		$Contig_Hit = $1;
		unless($Contig_Hit ~~ @contigs){
		    push(@contigs, $Contig_Hit); #store contig in @contigs array
		}
	    }
	}
    }
    open(PSEUDO, ">>$pseudogenome");
    if($pseudogenome_perl eq "yes"){ #write out contigs using perl
	foreach my $C(@contigs) {
	    $C =~ s/\>//;
	    $C =~ s/\s//;
	    {
		local $/ = ">"; #  change line delimter to read in file by contig
		open(GENOME, "$genome_assembly");
		while(<GENOME>) {
		    my $ContigSequence = $_; #store each contig of the genome file in $contig_seq
		    if ($ContigSequence =~ m/$C/i){
			$ContigSequence =~ s/\>//g;
			my $PseudoGenome = ">".$ContigSequence."\n";
			print(PSEUDO $PseudoGenome); #write out contig to pseudogenome
		    }
		}
	    }
	}
    }
    elsif($pseudogenome_perl ne "yes" && $pseudogenome_esl eq "yes"){ #write out contigs with esl-sfetch
	`esl-sfetch --index $genome_assembly`; #index genome (makes parsing genome much faster)
	foreach my $C(@contigs) {
	    $C =~ s/\>//;
	    $C =~ s/\s//;
	    `esl-sfetch $genome_assembly $C >> $pseudogenome`; #write out contigs to pseudogenome with esl-sfetch
	}
    }
    close GENOME;
    close PSEUDO;
    return $pseudogenome;
}


#Get contig names, add 'PSEUDOCONTIG' to mark end of genome file
sub get_contigs{
    my @contig_info = @_;
    my $genome_name = $contig_info[0];
    my $pseudogenome_perl = $contig_info[4];
    my $pseudogenome_esl = $contig_info[5];
    my $pseudo_file = $contig_info[6];
    my @names = ();
    if ($pseudogenome_perl eq "yes" || $pseudogenome_esl eq "yes"){
	open(PSEUDO, ">>$pseudo_file");
	print PSEUDO "\n>PSEUDOCONTIG\n";#This will be undone at the end of the script
	close PSEUDO;
	open (PSEUDO, $pseudo_file);
	while(<PSEUDO>){
	    if($_=~m/\>/){
		push @names,$_; 
	    }
	}
	push @names, ">PSEUDOCONTIG";
	close PSEUDO;
    }
    else{
	open(GENOME, ">>$genome_name");
	print GENOME "\n>PSEUDOCONTIG\n";
	close GENOME;
	open (GENOME, $genome_name);
	while(<GENOME>){
	    if($_=~m/\>/){
		push @names,$_;
	    }
	}
	push @names, ">PSEUDOCONTIG";
	close GENOME;
    }
    return @names;
}


#If reference receptors do not have unique names, assign unique names
sub rename_refs{
    my @refs = @_;
    my @receptor_array = ();
    my @receptor_list=();
    my $alpha = "a";
    foreach my $R(@refs){ #foreach sequence in reference file
	my $g_name="";
	my $g_seq="";
	my $renamed_seq;
	if($R=~m/(\>[\S]+)\n([\S\n]+)/){
	    $g_name = $1; #store header
	    $g_seq = $2; #store sequence
	}
	if($g_name ~~ @receptor_list){ #if header is duplicated
	    my $new_name = $g_name."_".$alpha; #rename with unique header by appending letters 
	    $alpha++; #move up alpahbet for next sequence
	    $renamed_seq = $new_name."\n".$g_seq; #rename sequence
	    push(@receptor_array, $renamed_seq); #store sequence with unique name in array
	    push(@receptor_list, $new_name); #store new name
	}
	else{
	    push(@receptor_list, $g_name); #if not duplicated, do not rename
	    $renamed_seq = $g_name."\n".$g_seq; #original sequence
	    push(@receptor_array, $renamed_seq); #store original sequence
	}
    }
    return @receptor_array;
}


#Hmmer filter: Takes prediction and scans for 7 transmembrane domain.
sub hmmsearch{
    my @hmminputs = @_;
    my $hmm_file = $hmminputs[0];
    my $prediction_file = $hmminputs[1]; #protein prediction file
    my $hmmrandom = $hmminputs[2];
    my $hmm_out = $hmmrandom."_hmmfilter.out";
    `hmmsearch --incE 0.05 $hmm_file $prediction_file > $hmm_out`;
    my $results_file;
    open(IN, $hmm_out);
    {
        local $/;
        $results_file = <IN>;
    }
    close IN;
    my $regex = "\-";
    my $expression = $regex x 37;
    my @hmmhits = split(/$expression/, $results_file);
    my $summary = $hmmhits[1];
    my @stats = split(/\n/, $summary);
    my $target_seqs = $stats[2];
    my $passed_hits = $stats[8];
    my $targets = "";
    my $hits = "";
    if($target_seqs){
	if ($target_seqs =~ m/\:[\s]+([0-9]+)\s/){
	    $targets = $1;
	}
    }
    if($passed_hits){
	if ($passed_hits =~ m/\:[\s]+([0-9]+)\s/){
	    $hits = $1;
	}
    }
    my $hmm_status = 0;
    if($targets){
	if($targets == 1){
	    if ($hits ==1){
		$hmm_status = 1;
            return $hmm_status;
	    }
	    else{
		$hmm_status = 0;
		return $hmm_status;
	    }
	}
	else{
	    $hmm_status = 0;
	    return $hmm_status;
	}
    }
    else{
	$hmm_status = 0;
	return $hmm_status;   
    }
}

