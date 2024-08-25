#!/usr/bin/perl

###################################################################################################################################
#
#
#
# 	bagpipe_phylo.pl
#		  
#    	Copyright (C) 2013-2024  Douglas Chesters and Anna Papadopoulou
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	contact address dc0357548934@live.co.uk
#
#
#
#
##########################################################################################################################################
#
#
#
#	rewrite of the script parse_clades.pl, 
#	which was part of the BAGpipe taxonomic annotation pipeline (Papadopoulou et al. 2014).
#	bagpipe_phylo.pl is written to be a general, standalone script for 
#	outputting likely taxonomies for queries, 
#	where the input is a phylogeny containing both queries and references (latter must have binomial species IDs).
#	
#	you need to download the NCBI taxonomy database, go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/
#	download the file taxdump.tar.gz, and unzip in the working directory. 
#	the script will read 2 of the resulting files: names.dmp and nodes.dmp.
# 	it will get from these the higher taxonomies to be used for assignment to species given in your tree.
#	only other input is a fasta file of the queries (-seqfile), the script simply identifies query names from this file 
#
# 	to run:
#
# 	perl bagpipe_phylo.pl -treefile RAxML_bestTree.nwk -seqfile seqs.fas -support 0.7
#
# 	use the command above to run this script.
# 	-treefile = ROOTED newick tree file name. Tested using raxml trees, so these will work best.
# 	-seqfile = the name for the fasta file of queries.
#	-support = cutoff for bootstrap support, can be proportion (0.0-1.0) or percentage (0-100)
#		if you have no support values, just put 0
#	-references = name of a file containing a list of the reference species to be used for taxonomic assignment.
#		this argument is not mandatory. if you dont give it, 
#		the program will assume all terminals of the tree which arnt queries, must be references.
#		however you may wish to restrict your analysis to a subset of all references e.g. regional taxa only.
#		if so, make a file list the names of reference species to be used. 
#		one species per line (Genus_species), same format as used in the tree.
#
#
#	Change Log:
# 	v1.02: 		included option for keyfile, so it can print out full taxon names.
#			uses different method for input of options in command line, a bit more flexable.
# 	v1.03: 		trying to read differently-rooted trees, reads root trifucation and root bifucation with no support.
# 	v1.04: 		reads newick trees with proportional support values in addition to the usual percentages
# 	v1.05: 		bugfix, spaces in query list
# 	20130923: 	bugfix in translation of common substrings ending in underscore
#	20130924:	and problem translating new format species strings, where common string ends with \d
#	20131009: 	added gpl. published with Papadopoulou et al. (2014).
#	20131112: 	another bugfix in rare problematic string translation
# 	20140617: 	option to run the script where full species names are used as ID's in the trees.
#	20140806 (v2.01): big re-write. inserted all code for directly reading the NCBI taxonomy DB, instead of the custom DB format.
#			so users can just download the files from NCBI and proceed.
#			discontinued use of 'tobycodes', use standard binomials instead.
#			new simple code for newick parsing.
#			changed program name. now works as a standalone program for taxonomic assignment. 
#			not neccessary to run from within the bagpipe framework. 
#			test_for_monospecific_clade discontinued, it was unreliable
#			published with Chesters, Zheng, Zhu (2015).
#	20140901: 	reading trees with no bootstrap, all nodes read as 100. 
#	20140928 (v2.02): user can give list of reference species subset, 
#			for using taxonomic information. other references are ignored if this is specified.
#	20141001: 	option to remove accessions from reference species
#	20150528: 	prints full assigned lineage on last column of results file.
#	20150718: 	bug fix, one column of results file had extra rubbish printed by mistake
#			script assumes strict use of binomials (required for getting taxonomies). 
#			however, if lineage is not returned for a reference ID, 
#			script now tries to parse just genus name and get taxonomy for that.
#	20151215: 	wont crash is a taxonomy is not found, 
#			in such cases just tests for whatevers on the reference label (which should be genus name at least)	
#	20151219: 	test to shortest reference were just calculated from a parent node of query (inner or outer)
#			not also adding branches from query to parent. now corrected.
#			note the reference terminal with shortest distance to query is calcualted
#			looking at both 'inner' and 'outer' nodes
#	20160130: 	quick implementation to read parsed jplace files (produced by evolutionary placement software)
#			needs further work, idealy to read jplace files directly, and perhaps write that foramt also
#	20160520: 	for CX and DY, queries are output to seperate fasta files, named by the assigned taxon
#			given by command -sorted_fasta [results column, 4 or 6] [filename_prefix]
#	20240825: 	Established Github repository. Tool being used in current Lepidoptera study.
#	
#	
#	
#	
#	
#	
#
#
##########################################################################################################################################





$arg_string  = join ' ', @ARGV;



# additional user options:

$verbose 		= 1; # 0 simple screen output, 1=verbose, 2=ulta-verbose

$minimum_required_number_of_non_diet_terminals_inner_clade = 1;
$minimum_required_number_of_non_diet_terminals_outer_clade = 2;

$remove_accessions_from_reference_species = 1;

##########################################################################################################################################


# globals:

$treefile;
$keyfile;
$fas_file;
$reference_file;
$output_filename;# 	= "$treefile.diet_clades";
%species_names;
$starting_node;
$support_cutoff = "NA";
%filter_these_assignment;

#####################################
read_command_arguments($arg_string);#
#####################################


###############################################################################################################################################

# from parse_ncbi_tax_db ...

$ignore_subspecies			= 1;	# 0=read taxonomies beyond species level (e.g. include subspecies, varietas, and unnamed bacteria ranks)
						# if 1, these extra ranks are ignored
						# this applies to full species name, assigned to given ncbi tax number
# globals

%ncbi_nodes;	# ncbi taxonomy structure is recorded in this object




###############
store_nodes();#
###############



###################
parse_namesfile();#
###################


$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;
print "name for the starting node ($starting_node) which has been specified:$starting_name\n";
print "traversing NCBI taxonomy tree from this node. recording lineage for each taxon.\n\n";
$starting_name =~ s/^(\w).+$/$1/;

#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################



###############################################################################################################################################



my %terminals = ();
my $root_identity = "";
my $date = localtime time;

%query_IDs;	# the fasta file contains names of queries. 
		# each are stored in the hash %query_IDs.


############
read_fas();#
############


# if the user wants to limit taxonomic work to a subset of the references:

%reference_species;
$how_many_references_will_be_used = 0;

if($reference_file =~ /[\w\d]/)
	{
	open(IN7, $reference_file) || die "\n\nerror. 
you have specified to use a file with list of reference species ($reference_file). but this file cannot be found.\n";
	while(my $line = <IN7>)
		{
		#print "line:$line";
		$line =~ s/\n//;$line =~ s/\r//;
		if($line =~ /(.+)/)
			{
			my $sp = $1;#print "\tsp:$sp\n";
			$reference_species{$sp} = 1;$how_many_references_will_be_used++;
			};
		};
	close IN7;

	print "\nreference_file:$reference_file has been read.\n";
	print "how_many_references_will_be_used:$how_many_references_will_be_used\n";

	my @testkeys = keys %reference_species;

	if($#testkeys <= 1)
		{die "\nwarning. you specified a file with references to be used. but non were recorded. quitting.\n\n"}

	};


# if you wish to filter those assignments with bad ML scores

#$raxml_ML_classification = "/home/douglas/scripted_analyses/COIevol/classification/RAxML_classification.Yu2012.raxmlEPA.G-0.12";






%nodes;		# user phylogeny stored in this hash
$root_node;	# this global will keep the value of the last in the loop, which will be the root



# read tree into hash. node identities for hash key, then parent/child identities as hash entry

#########################
record_tree2($treefile);#
#########################



my @test_query_IDs = keys %query_IDs;
foreach my $id(@test_query_IDs)
	{
	#print "query:$id\n";
	unless (exists($terminals{$id}))
		{print "\n\nWARNING:member found in fasta file ($id) is absent from tree. are you sure you are using the correct pair of files?\n\n"}
	};


$terminals_belonging_to_current_node;


# analysis start:

open(OUT_RESULTS, ">$output_filename") || die "\nerror 63.\n";

print OUT_RESULTS 	"input\tquery\touter_support\touter_taxon_id\t" ,
			"inner_support\tinner_taxon_id\t" , 
			"distance_to_reference_tree\tshortest_distance_to_reference_leaf\treference_terminal_with_shortest_dist\tshared_lineage\n";




###################################
get_clades_for_each_query();#
###################################

close(OUT_RESULTS);













print "\nnumber queries assigned according to rank ....\n";
my @ranks = keys %number_queries_assigned_to_this_rank;
my @taxas = keys %store_assigned_tax;
				$number_queries_assigned_to_this_rank{$rank}++;
				$number_queries_assigned_to_this_tax{$rank}{$assigned_tax}++;
				${$assigned_tax}=1;

foreach my $rank(@ranks)
	{
	my $count = $number_queries_assigned_to_this_rank{$rank};
	print "$rank\t$count\n";
	foreach my $tax(@taxas)
		{
		my $count2= $number_queries_assigned_to_this_tax{$rank}{$tax};
		if($count2 >= 1){print "\t$tax\t$count2\n"}
		};
	}

if($sorted_fasta == 1)
	{
	######################################################
	write_taxonomically_sorted_fasta_files_for_queries();#
	######################################################
	};
	
#	$sorted_fasta_column = $1; $sorted_fasta_prefix = $2;


print "
user taxon:$starting_name_ff
";

print "\n\n\nend of script\n";
exit;





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub record_tree2
{


my $tree1= "";
open(IN, $treefile) || die "\n\nerror 1408 cant open file $treefile\n\n";
while (my $line = <IN>)
	{
	if($line =~ /(.+\(.+)/){$tree1 = $1}
	}
close(IN);

print "\nlooking at tree:$treefile\n";

$tree1 =~ s/ //g;



$tree_parse = 0; 	# assuming sensible tree.
if($tree_parse ==1)
	{
	# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
	$tree1 =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/g; 
	# remove regular branchlengths: 0.02048
	$tree1 =~ s/\:\-*\d+\.\d+//g; 
	# remove 0 length branchlengths
	$tree1 =~ s/\:\d+//g;
	}


my $newick_string 	= $tree1;
my $interal_node	= 0;



# new newick parser ....

while ($newick_string =~ s/\(([^\(\)]+)\)(\d*)/INTERNAL_NODE_$interal_node/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
	{
	my $node = $1;my $boot = $2; my $nodeID = "INTERNAL_NODE_$interal_node"; #print "nodeID:$nodeID node:$node\n";

	# found a pair of parentheses (with no parentheses inside).
	# this string is taken, and replaced with a single node (INTERNAL_NODE_\d)
	# node is numbered (\d in variable $interal_node), 
	# the number being incremented at the end of the loop
	# find what is at the current node (decendents), by splitting at the commas

	if($boot =~ /\d/){$nodes{$nodeID}{support} = $boot}else{$nodes{$nodeID}{support} = 100};#print "\tnodeID:$nodeID boot:$boot\n";

	my @child_nodes = split /\,/ , $node;
	$child_counts{$nodeID} = $#child_nodes;

	for $i(0 .. $#child_nodes)
		{
		#print "$child_nodes[$i]\n";
		my $bl = "NA"; 	
		if($child_nodes[$i] =~ /\:(.+)/)	# branchlength found
			{$bl = $1};
		$child_nodes[$i] =~ s/\:(.+)//;
		#print "node:$child_nodes[$i]\tbl:$bl\n ";

		if($remove_accessions_from_reference_species == 1)
			{$child_nodes[$i] =~ s/^([A-Z][a-z]+_[a-z]+)_.+$/$1/}

		# record each decendent of the current node, in a hash nodes{ the current node }{ child number }

		$nodes{$nodeID}{$i} 			= $child_nodes[$i];
		$nodes{$child_nodes[$i]}{parent} 	= $nodeID;
		$nodes{$child_nodes[$i]}{branchlength} 	= $bl;

		# and record whether the current node is a terminal

		unless($child_nodes[$i] =~ /INTERNAL_NODE_/){$terminals{$child_nodes[$i]} =1}

		}
	#print "node:$interal_node\n\tchild nodes:@child_nodes\n";

	$root_node = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.

	$interal_node++;
	}


print "your tree has been read, it has $interal_node nodes.\n";
#print "newick_string:$newick_string\n";die;


unless($interal_node >= 2){die "\nerror reading your phylogeny.\n"}



}#sub record_tree2




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub read_fas
{

open(IN, $fas_file) || die "\nerror ... fasta file you have given ($fas_file) cannot be opened. 
make sure it is named correctly and in the working directory. quitting\n";

my $count_diet_IDs = 0;
while (my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ />(.+)/)
		{
		my $id = $1;#print "id:$id\n";
		$count_diet_IDs++;
		$query_IDs{$id}= 1;
		}
	}
close(IN);


print "\n$count_diet_IDs fasta IDs in query file $fas_file. will be looking for these in the tree ...\n\n";

unless($count_diet_IDs >= 1)
	{die "\nerror ... no fasta ID's found in the file ($fas_file). quitting.\n"}

}







#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub get_clades_for_each_query
{

my @all_queries = keys %query_IDs;@all_queries = sort @all_queries;

my $count_diet_member =0;

foreach my $diet_member(@all_queries)
	{
	my $outer_support = "NA";my $outer_taxon_id = "NA";my $outer_included_tax = "NA";my $outer_nonqueries = "NA";my $outer_queries = "NA";
	my $inner_support = "NA";my $inner_taxon_id = "NA";my $inner_included_tax = "NA";my $inner_nonqueries = "NA";my $inner_queries = "NA";
	my $outer_species = "NA";my $inner_species = "NA";
	my $outer_taxon_string_is_fullname =0;my $inner_taxon_string_is_fullname =0;
	my $shared_lineage = "";

	my $distance_to_reference_tree = $nodes{$diet_member}{branchlength};# length not enough, need to proceed until reference members are included in clade
		# distance to closest leaf is common measure.

	my $parent_node = $diet_member;

	if($verbose >= 1){print "\n\n\n****** TESTING NEW QUERY, number $count_diet_member of $#all_queries ******\n";$count_diet_member++;
		print "\nquery:($diet_member). subtending_branchlength:$distance_to_reference_tree. heading towards root ......\n";}

	my $at_root = 0;my $check_count=0;my $found_node =0;my $sum_branch_dist_to_current_parent =0;

	while ($at_root == 0)
		{
		# progress from diet terminal, one node at a time, towards root

		$parent_node = $nodes{$parent_node}{parent};
		my $support = $nodes{$parent_node}{support};
		my $parent_branchlength = $nodes{$parent_node}{branchlength};
		if($parent_branchlength =~ /[\d\.]+/){$sum_branch_dist_to_current_parent += $parent_branchlength};

		$count_non_diet_terminals=0;$terminals_belonging_to_current_node = "";$query_terminals_belonging_to_current_node = "";
		#@distances_to_reference_leaves_from_this_node = ();
		$shortest_distance_to_a_reference_terminal = 100;$reference_terminal_with_shortest_dist = "NA";

		############################################
		get_terminals_from_this_node($parent_node , $sum_branch_dist_to_current_parent);#
		############################################

		print "test1:reference_terminal_with_shortest_dist:$reference_terminal_with_shortest_dist\n";

		if($verbose >= 1){print "\tnode ($parent_node) count non-diet terminals from this node ($count_non_diet_terminals)\n";
			print "\tbootstrap support for this node is ($support) , user defined support cutoff ($support_cutoff)\n";
		print "node($parent_node) boot($support) dist_to_ref($distance_to_reference_tree)\n";
		print "sum dist to this parent node from query:$sum_branch_dist_to_current_parent\n";
}
		#$support =1;
		if($support >= $support_cutoff && $count_non_diet_terminals >= $minimum_required_number_of_non_diet_terminals_inner_clade)
			{

			#print "distance to reference leaves = @distances_to_reference_leaves_from_this_node\n";
			#if($count_diet_member >= 6){die};

			# stop progressing towards root, relevent inner node has been found.

			$found_node  = 1;$inner_support = $support;$inner_nonqueries=$terminals_belonging_to_current_node;
			$inner_queries = $query_terminals_belonging_to_current_node;
			if($verbose >= 1){print "found appropriate node! \n";
				#print "members derived from this node:$terminals_belonging_to_current_node\n\n";
				}

			my $shared_tobycode = "";
			my $sub_results;


			##################################################
			$sub_results = get_shared_taxonomic_name($terminals_belonging_to_current_node);#
			##################################################
		
			if($sub_results =~ /^([^\t]+)\t(.+)/){$shared_tobycode = $1; $shared_lineage = $2};


			$inner_taxon_id = $shared_tobycode;print "\ninner_taxon_id:$inner_taxon_id\n";
			if($verbose >= 1){print "\ntaxonstring shared by these terminals is ($shared_tobycode)\n\n";		
				print "now finding OUTER CLADE\n"}

			my $at_root2 = 0;my $check_count2=0;my $found_node2 =0;
			while ($at_root2 == 0)
				{

				# starting at the inner node, again move towards the root

				$parent_node = $nodes{$parent_node}{parent};
				my $support = $nodes{$parent_node}{support};
				my $parent_branchlength = $nodes{$parent_node}{branchlength};
				if($parent_branchlength =~ /[\d\.]+/){$sum_branch_dist_to_current_parent += $parent_branchlength};

				$count_non_diet_terminals=0;$terminals_belonging_to_current_node = "";$query_terminals_belonging_to_current_node = "";

				############################################
				get_terminals_from_this_node($parent_node , $sum_branch_dist_to_current_parent);#
				############################################

		print "test2:reference_terminal_with_shortest_dist:$reference_terminal_with_shortest_dist\n";


				if($verbose >= 1){
					print "\tpotential outer node ($parent_node) count non-diet terminals from this node ($count_non_diet_terminals)\n";
					print "\tsupport (bootstrap) for this node is ($support) , user defined support cutoff ($support_cutoff)\n";
					print "sum dist to this parent node from query:$sum_branch_dist_to_current_parent\n";
					}

				if($support >= $support_cutoff && $count_non_diet_terminals >= $minimum_required_number_of_non_diet_terminals_outer_clade)
					{

					# outer node has been found, stop the second round of root-bound moves

					$found_node2  = 1;$outer_support = $support;$outer_nonqueries=$terminals_belonging_to_current_node;
					$outer_queries = $query_terminals_belonging_to_current_node;
				#	$outer_species = get_species_list($terminals_belonging_to_current_node);

					if($verbose >= 1){
						print "found OUTER node! \n";
					#	print "members derived from this node:$terminals_belonging_to_current_node\n\n";
						}

					my $shared_tobycode = "";my $sub_results;

					#################################################################################
					$sub_results = get_shared_taxonomic_name($terminals_belonging_to_current_node);#
					#################################################################################
			
					if($sub_results =~ /^([^\t]+)\t(.+)/){$shared_tobycode = $1; $shared_lineage = $2};


					if($verbose >= 1){print "\nOUTER NODE. taxonstring shared by these terminals is ($shared_tobycode)\n\n"}
		
					$outer_taxon_id = $shared_tobycode;
					last;

					}#if($support >= $support_cutoff ...... 

				unless($parent_node =~ /INTERNAL_NODE_/)	{$at_root2 =1;print ""}
				if($parent_node == $root_identity)	{$at_root2 =1;print ""}	
				$check_count2++;if($check_count2>=100000){last}

				}

			last;
			}else{#if($support >= $support_cutoff .......
	
			if($distance_to_reference_tree =~ /\d/ && $parent_branchlength =~ /\d/)
				{
				$distance_to_reference_tree += $parent_branchlength;
				}else{
				$distance_to_reference_tree = "NA";
				};

			};

	
		unless($parent_node =~ /INTERNAL_NODE_/)		
			{$at_root =1;#print "test1\n"
			};#print "no parent. skipping .. \n"}
		if($parent_node eq $root_node)
			{$at_root =1;#print "test2 root_node:$root_node\n"
			};#print "parent node ($parent_node) == root ($root_identity) ... skipping\n"}	
		$check_count++;
		if($check_count>=100000)
			{#print "test3\n";
			last}

		}#while ($at_root == 0)


$outer_nonqueries =~ s/\t/ /g;$outer_queries =~ s/\t/ /g;$inner_nonqueries =~ s/\t/ /g;$inner_queries =~ s/\t/ /g;
my $print_fas_file = $fas_file;$print_fas_file =~ s/.+(....................)$/"...." . $1/e;

	
	if($filter_these_assignment{$diet_member} == 1){$found_node = 0;print "ignoring current assignment, it has low raml score\n"};

	if($arg_string =~ /raxml_jplace_parsed/)
		{
		#print "arguments =~ raxml_jplace_parsed\n";
		unless($keep_these_assignment{$diet_member} == 1)
			{$found_node = 0;#print "ignoring current assignment, it has low raml score\n";
			};
		};

	# if($diet_member =~ /411159039.+JX828806.+Alysia/){print "Query:$diet_member, found_node:$found_node arguments:$arguments";die};

	if($found_node == 1)# || $is_monospecific == 1)
		{

		my $shortest_distance_to_reference_leaf = $shortest_distance_to_a_reference_terminal+$distance_to_reference_tree;

	#	print "shortest_distance_to_reference_leaf:$shortest_distance_to_reference_leaf\n";

	#	print	"($print_fas_file)\t($diet_member)\t($outer_support)\t($outer_taxon_id)\t" ,
	#		"($inner_support)\t($inner_taxon_id)\t" , 
	#		"distance_to_reference_tree:($distance_to_reference_tree)\tshortest_distance_to_reference_leaf:($shortest_distance_to_reference_leaf)\t($reference_terminal_with_shortest_dist)\t($shared_lineage)\n";


		if($outer_taxon_id eq ""){$outer_taxon_id = "NA"};if($inner_taxon_id eq ""){$inner_taxon_id = "NA"};

		print OUT_RESULTS 	
			"$print_fas_file\t$diet_member\t$outer_support\t$outer_taxon_id\t" ,
			"$inner_support\t$inner_taxon_id\t" , 
			"$distance_to_reference_tree\t$shortest_distance_to_reference_leaf\t$reference_terminal_with_shortest_dist\t$shared_lineage\n";

		#print "shared_lineage:$shared_lineage\n";
		my @ranks = ("order", "family", "genus", "species");
		foreach my $rank(@ranks)	
			{
			if($shared_lineage =~ / $rank\:([a-z]+)/i)
				{
				my $assigned_tax = $1;
				$number_queries_assigned_to_this_rank{$rank}++;
				$number_queries_assigned_to_this_tax{$rank}{$assigned_tax}++;
				$store_assigned_tax{$assigned_tax}=1;
				}
			};
		};

#	if($diet_member =~ /prospera/){die ""};
	}#foreach my $diet_member(@all_queries)


}#sub get_clades_for_each_diet_member




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub get_terminals_from_this_node
{
my $next = $_[0];my $sum_branchlength = $_[1];

my $child_nodes = $child_counts{$next};
my @next1 = ();
for $i(0 .. $child_nodes)
	{
	push @next1 , $nodes{$next}{$i};	#print "node:$next sum:$sum_branchlength i:$i child:$nodes{$next}{$i}\n";
	}

for my $index(0 .. $#next1)
	{
	my $test = @next1[$index];#print "test:$test\n";
	my $branchlength_to_child = $nodes{$test}{branchlength};
	if($nodes{$test}{branchlength} =~ /\d/)
		{
		#print "\nerror 567. ($test) \n"
		#$branchlength_to_child = $nodes{$test}{branchlength};
		}else{
		$branchlength_to_child = 0;
		}

	if($test =~ /^INTERNAL_NODE_/)
		{

		#####################################
		get_terminals_from_this_node($test , $sum_branchlength+$branchlength_to_child );#	# recurse
		#####################################

		}else{

		if(exists($query_IDs{$test}))
			{
			$query_terminals_belonging_to_current_node .= $test . " ";#print "exists:$test\n";
			}else{

			# here user may require a subset of the references testing,
			# so determine whether current terminal of of this subset, 
			# only if so, increment $count_non_diet_terminals

			#print "reference terminal($test)\n";

			my $include_this_reference = 0;
			if($how_many_references_will_be_used == 0 )
				{
				$include_this_reference = 1;#print "\tno list, so including\n"
				}else{
				if($reference_species{$test} == 1)
					{
					$include_this_reference = 1;
					#print "\tfound in list\n";
					};
				}

			if($include_this_reference == 1)
				{

				# reference tip reached
			#	print "reference tip reached sum:$sum_branchlength bl child:$branchlength_to_child\n";
				$new_val = $sum_branchlength+$branchlength_to_child;

				if($new_val =~ /\d/ && $new_val<$shortest_distance_to_a_reference_terminal)
					{
					$shortest_distance_to_a_reference_terminal = $new_val;
					$reference_terminal_with_shortest_dist = "$test";
					};


				$terminals_belonging_to_current_node .= $test . "\t";$count_non_diet_terminals++;#print "does not exist:$test\n";

				};

			}

		}
	}

return();

}






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################






sub get_shared_taxonomic_name
{
my $tax_list = shift;
my @tax_array = split /\t/ , $tax_list;

my $shared_tax_substring = "";

if($verbose >= 1){
print "\nfinding largest common substring  ....\n\n";
}


my $index=0;

my $test_species = $tax_array[0];
#$test_species =~ s/Ceratina_nGU321539/Ceratina/;
my $test_taxonomy = $complete_lineage_for_this_species{$test_species};

unless($test_taxonomy =~ /\w+/)
	{
	$test_species =~ s/^([A-Z][a-z]+)_.+/$1/;# try genus only ...
	$test_taxonomy = $complete_lineage_for_this_species{$test_species};
	unless($test_taxonomy =~ /\w+/)
		{
		print "\nerror 770, no taxonomy found for this:($test_species)\n";
		$test_taxonomy = " no_rank:cellular_organisms genus:$test_species";
		}
	};

#print "test_taxonomy:$test_taxonomy\n";
my $most_inclusive_name = "";
my $most_inclusive_lineage = "";

while($test_taxonomy =~ s/^\s*([^:]+):(\w+)//)
	{
	my $current_rank = $1;my $current_taxname = $2;
	#print "current_rank:$current_rank current_taxname:$current_taxname\n";

	my $all_members_have_this_name =1;
	foreach my $tax(@tax_array)
		{	
		my $test_taxonomy2 = $complete_lineage_for_this_species{$tax};
	#	print "\tspecies:$tax complete taxonomy:$test_taxonomy2\n";
		if($test_taxonomy2 =~ /\:$current_taxname\s/)
			{
			#print "1";
			}else{
			$all_members_have_this_name =0;#print "0";
			}

		}
	#print "\nall_members_have_this_name:$all_members_have_this_name\n";

	if($all_members_have_this_name == 1)
		{
		$most_inclusive_name = $current_taxname;
		$most_inclusive_lineage .= "$current_rank:$current_taxname ";
		};

	}

my $returnstring = $most_inclusive_name . "	" .  $most_inclusive_lineage;
return($returnstring);


}





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub test_for_monospecific_clade
{

my $diet_terminal = shift;

my $parent_node 	= $nodes{$diet_terminal}{parent};
my $nodesupport		= $nodes{$diet_terminal}{bootstrap};
my $grandparent_node 	= $nodes{$parent_node}{parent};

$count_non_diet_terminals=0;$terminals_belonging_to_current_node = "";


#################################################
get_terminals_from_this_node($grandparent_node);#
#################################################



my $is_diet_seq_monospecific = 0;my $monospecies = "NA";

if($terminals_belonging_to_current_node =~ /^([^_]+)_.+\t([^_]+)_/)
	{
	my $tobycodeA=$1;my $tobycodeB=$2;
	if($tobycodeA eq $tobycodeB)
		{
		$monospecies = $tobycodeB;$is_diet_seq_monospecific = 1;
		}
	}

if($verbose == 1)
	{
	print "\n\ntest_for_monospecific_clade.\n";		
	print "diet_terminal:$diet_terminal parent_node:$parent_node grandparent_node:$grandparent_node\n\t non diet terminals from grandparent:$terminals_belonging_to_current_node\n";

	if($is_diet_seq_monospecific == 1 && $nodesupport >= $support_cutoff)
		{
		print "\n$diet_terminal IS monospecific with species $monospecies\n\n"
		}else{
		print "\n$diet_terminal is NOT monospecific\n\n";
		}

	}

return($is_diet_seq_monospecific);


}



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub get_species_list
{
my $list_of_terminals = shift;

my %hash_of_terminals = ();

if($verbose == 1){print "list_of_terminals:$list_of_terminals\n"};

my @array_of_terminals = split /\t/, $list_of_terminals;

foreach my $termnal(@array_of_terminals)
	{
	$termnal =~ s/^([^_]+)_.+/$1/;
	$hash_of_terminals{$termnal}= 1;
	}

my @array_of_terminals2 = keys %hash_of_terminals;

if($verbose == 1){
#print "list_of_species:@array_of_terminals2\n";
}

my $list_of_species = join ' ', @array_of_terminals2;

return($list_of_species);

}



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;

print "
\n\n
     ***** bagpipe_phylo.pl *****
		  
\n";


# perl parse_clades.pl -treefile RAxML_bipartitions.diet_group.0.reformatted -seqfile diet_group.0.fas -keyfile key_Mar2013_Magnoliophyta
#$treefile;
#$fas_file;
#$output_filename 	= "$treefile.diet_clades";



if($arguments =~ /-treefile\s+(\S+)/)
	{
	$treefile = $1;
	open(INTEST , $treefile) || die "\ntree file not found in current directlry:$treefile\n";
	close INTEST;
	}else{
	print "\nerror reading command arguments  (-treefile)\n";die""
	}

if($arguments =~ /-support\s+(\S+)/)
	{
	$support_cutoff = $1;
	}else{
	$support_cutoff = 50;
	print "user did not give cutoff for acceptable boot support (not relevent if you are not using bootstrapped trees)\nusing default ($support_cutoff).\n"
	}

if($arguments =~ /-seqfile\s+(\S+)/)
	{
	$fas_file = $1;
	}else{
	}

if($arguments =~ /-references\s+(\S+)/)
	{
	$reference_file = $1;
	}else{
	}

if($arguments =~ /-raxml_classification\s+(\S+)/)
	{
	$raxml_classification = $1;
open(INTEST, $raxml_classification) || die "\nerror 1064, cant open file ($raxml_classification)\n";
my $raxml_dists_read = 0;
while (my $line = <INTEST>)
	{
	if($line =~ /^(\S+)\s(I\d+)\s([\d\.]+)\s([\d\.]+)/)
		{my $query=$1; my $score = $4;#print "query:$query score:$score\n";
		if($score > 0.07){$filter_these_assignment{$query}=1;$poorassign++;};
		}
	};
print "assignments with low ML score to be ignored:$poorassign\n";

	}else{
	}


if($arguments =~ /-raxml_jplace_parsed\s+(\S+)/)
	{
	$raxml_classification = $1;
	open(INTEST, $raxml_classification) || die "\nerror 1082, cant open ($raxml_classification)\n";
	my $raxml_dists_read = 0;
	while (my $line = <INTEST>)
		{
#$query\t$edge_num\t$likelihood\t$like_weight_ratio\t$distal_length\t$pendant_length\n";
		if($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/)
			{my $query=$1; my $like_weight_ratio = $4;#print "query:$query score:$score\n";
			$keep_these_assignment{$query}=1;$assignmens_read_from_parsed_jpalce++;
			}
		};

	print "
	assignmens_read_from_parsed_jpalce:$assignmens_read_from_parsed_jpalce
	";
	};




if($arguments =~ /-node\s+(\d+)/)
	{
	$starting_node = $1;
	}else{
	print "user did not given NCBI taxonomic number. using default 33208 (Metazoa)
if you have references outside this default, you need to get the appropriate NCBI taxonomy number from http://www.ncbi.nlm.nih.gov/taxonomy/
and input this here with the option -node taxon_number
";
	$starting_node = 33208;
	}

if($arguments =~ /-sorted_fasta\s+([^\-\s]+)\s([^\-\s]+)/)
	{
	$sorted_fasta = 1;$sorted_fasta_column = $1; $sorted_fasta_prefix = $2;
	unless($sorted_fasta_prefix =~ /[a-z0-9]/i){die "\nunexpected specification of output fasta prefix. quitting\n"};
	}elsif($arguments =~ /-sorted_fasta/){
	die "\ncommand error, syntax: -sorted_fasta [column] [filename_prefix]\n"
	};




$output_filename	= "$treefile.query_clades";


print "\n
user options have been read.....
support_cutoff:			$support_cutoff
taxon node encompassing refs:	$starting_node
treefile:			$treefile
query fasta file:		$fas_file
reference_file:			$reference_file

output will be written to file:	$output_filename


";



}#sub read_command_arguments



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub read_keyfile
{
open(IN, $keyfile) || die "\n\nerror 91. cant open $keyfile\n";
my $species_strings_read=0;

while (my $line = <IN>)
	{
	#print "\n$line";

	#MeeRBBegrac Berberis_gracilis 258166 species no_rank:eudicotyledons no_rank:eudicotyledons order:Ranunculales family:Berberidaceae genus:Berberis species:gracilis

	if($line =~ /^(\S+)\s(\S+)\s(\S+)\s(\S+)\s(.+)/)
		{
		my $tobycode = $1;my $species_name = $2;my $ncbi_number = $3; my $rank = $4; my $taxonomic_path = $5;
		#print "\ttobycode:$tobycode species_name:$species_name ncbi_number:$ncbi_number rank:$rank taxonomic_path:$taxonomic_path\n";
		$species_names{$tobycode} = $species_name;$taxonomic_paths{$tobycode} = $taxonomic_path;
		$species_strings_read++;
		$tobycodes{$species_name} = $tobycode;

	#	if($tobycode eq "MLAA4"){print "$tobycode $species_name ... quit\n\n";die ""}

		}
	}

close(IN);

if($species_strings_read == 0){die "\n\nerror 112. seems to be problem reading keyfile, lines in that file dont match regex.\n"}

print "$species_strings_read species strings read from file, plus associated taxonomic information\n";

}




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub store_nodes
{

my %rank_hash;

# nodes.dmp contains each taxon, described in a tree-like structure. 
# so each node has a name (eg polyphaga) and the higher group node to which it belongs (beetles). 

open(NODES , "nodes.dmp") || die "cant open nodes.dmp
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";

print "\nreading files from NCBI taxonomy database .... nodes.dmp .... ";


my $line_counter=0;
while (my $line = <NODES>)
	{
	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]+)\t[\|]\t/)
		{
		my $tax_id = $1;my $parent_tax_id = $2;my $rank = $3;
		$rank_hash{$rank}++;
#		print "tax_id:$tax_id parent_tax_id:$parent_tax_id rank:$rank\n";

			$ncbi_nodes{$tax_id}{rank} = $rank;
		#	$ncbi_nodes{$tax_id}{rank_code} = $current_rankcode;
			$ncbi_nodes{$tax_id}{parent} = $parent_tax_id;
			$ncbi_nodes{$parent_tax_id}{child_nodes} .= "\t" . $tax_id;

		}else{
		print "line_counter:$line_counter line:$line";
		die "UNEXPECTED LINE:$line\nquitting\n";
		}
	$line_counter++;
	}

close(NODES);

#my @ranks = keys %rank_hash;@ranks = sort @ranks;
#print "ranks found in nodes.dmp:\n";
#print LOG "ranks found in nodes.dmp:\n";

#foreach my $rank(@ranks){print "$rank\t" , $rank_hash{$rank} , "\n";print LOG "$rank\t" , $rank_hash{$rank} , "\n"};

my @all_nodes = keys %ncbi_nodes;@all_nodes = sort @all_nodes;

print scalar @all_nodes , " nodes have been read.\n";


}




#####################################################################################################
#
#
#
#####################################################################################################


sub parse_namesfile
{

# here just parse the scientific name of each node. ignore synonyms etc

open(NAMES , "names.dmp") || die "cant open names.dmp
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";


print "\nnames.dmp, parsing 'scientific name', ignoring others ... ";

my $names_line_counter=0;
while (my $line = <NAMES>)
	{
# 24	|	Shewanella putrefaciens	|		|	scientific name	|

	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]*)\t[\|]\tscientific name/)
		{
		my $tax_id = $1;my $name = $2;#my $rank = $3;
		# print "tax_id:$tax_id name:$name\n";

		# if you want to remove non-alphanumerical characters from assigned species names:
		$name =~ s/[\(\)\,\[\]\'\#\&\/\:\.\-]/ /g;
		$name =~ s/\s\s+/ /g;$name =~ s/\s+$//;$name =~ s/^\s+//;
		$ncbi_nodes{$tax_id}{name} = $name;

		$names_line_counter++;#print "$names_line_counter\n";


		}else{
		if($line =~ /^(\d+).+scientific name/){die "UNEXPECTED LINE:\n$line\nquitting\n"}
		}

	}

close(NAMES);

print "$names_line_counter names parsed.\n";


}



#####################################################################################################
#
#
#
#####################################################################################################




sub traverse_nodes
{
my $current_node = $_[0];my $current_node_taxstring = $_[1];# $current_node_taxstring to be deprecated

my $child_nodes = $ncbi_nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;
my @child_nodes_array = split(/\t/, $child_nodes);



if($current_node == $starting_node)
	{
	my $rank = $ncbi_nodes{$starting_node}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $current_node_taxstring = substr $current_node_taxstring, 0,1;
	my $originalname = $ncbi_nodes{$starting_node}{name};$originalname =~ s/[\s\t]/_/g;
	my $child_complete_lineage = $ncbi_nodes{$starting_node}{complete_lineage} . "$rank:$originalname ";
	}



foreach my $child(@child_nodes_array)
	{
	unless($current_node == $child)# one the child nodes of the root node (1), is also 1 
	{
	my $rank = $ncbi_nodes{$child}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	
	
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/; ####################### sep2013
	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name
	
	$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;
	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;

	$ncbi_tax_number_for_this_species{$name_string}=$child;

#	print "complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";

	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;

	my $name_assignment_to_taxnumber = "";

	if($ignore_subspecies == 1)
		{
		if ( $rank eq "subspecies" || $rank eq "varietas"   || $nodes{$current_node}{complete_lineage} =~ / species:/)# 
			{
			my $parentoriginalname = $ncbi_nodes{$current_node}{name};
			if($parentoriginalname =~ /^[A-Z][a-z]+\s[a-z]+$/)
				{
				$parentoriginalname =~ s/[\s\t]/_/g;
				$name_assignment_to_taxnumber = $parentoriginalname;
			#	print "node:$child named:$originalname rank:$rank appears to be subspecfic\n";
			#	print "\tassigning parent name instead:$parentoriginalname\n\n";
				}else{$name_assignment_to_taxnumber = $originalname}
			}else{
			$name_assignment_to_taxnumber = $originalname
			}

		}else{
		$name_assignment_to_taxnumber = $originalname
		}

	#print "$name_assignment_to_taxnumber $child $rank $child_complete_lineage\n";


		###########################################
		traverse_nodes($child , $child_taxstring);#
		###########################################
	}}


	
}#sub traverse_nodes





#####################################################################################################
#
#
#
#####################################################################################################



sub write_taxonomically_sorted_fasta_files_for_queries
{


print "
you have specified to output queries in seperate fasta files for each assigned taxon.
you have specified that taxa are read from this column of the results file of this program:$sorted_fasta_column
	( should probably be either columns 4 or 6, 4 has conservative taxonomic assignments, 6 liberal)
and output fasta files have this prefix:$sorted_fasta_prefix
";
#	$sorted_fasta_column = $1; $sorted_fasta_prefix = $2;


open(IN3, "$output_filename") || die "\nerror 1392, cant open file which should have just been written.\n";
print "opened results file $output_filename\n";
my $lines_read = 0;
while (my $line = <IN3>)
	{
	$lines_read++;print "$line";
	$line =~ s/\n//;$line =~ s/\r//;

# input	query	outer_support	outer_taxon_id	inner_support	inner_taxon_id	distance_to_reference_tree	shortest_distance_to_reference_leaf	reference_terminal_with_shortest_dist	shared_lineage
# ....6Svert.fasta.pynast2	uniques_1000864	100	Spalacidae	100	Rhizomys_pruinosus	0.111417	0.477840004	Rhizomys_pruinosus	no_rank:Theria no_rank:Eutheria no_rank:Boreoeutheria superorder:Euarchontoglires no_rank:Glires order:Rodentia suborder:Sciurognathi no_rank:Muroidea family:Spalacidae 

	my @split = split /\t/ , $line;
	my $tax = $split[$sorted_fasta_column - 1];print "\tTAX:$tax\n";
	if($tax =~ /\w/)
		{
		my $query_ID = $split[1];
	#	unless($queries_assigned_to_each_taxon{$tax}=~ /\t$query_ID\t/)
	#		{
	#		$queries_assigned_to_each_taxon{$tax} .= "\t$query_ID\t";
	#		};
		# this way instead:
		$which_taxon_is_this_query_assigned_to{$query_ID} = $tax;
		};

	};
close IN3;
print "lines_read from results file:$lines_read
for each taxon in column $sorted_fasta_column, the queries assigned to it have been stored.
next will read query sequences (you have specified file name:$fas_file), and organize seqs taxonomically ....
\n";


##################################


# open(IN5 , $fas_file) || die "\nerror 1429, cant open query fasta file named $fas_file\n";

open(FASTA_IN, $fas_file) || die "Cant open $database_file.\n";
my $fasta_entry = "";
while(my $fasta_line= <FASTA_IN>)
	{
	if($fasta_line =~ /^>.+/)
		{
		unless(length($fasta_entry)<=2)
			{
			$fasta_entry =~ s/^>//;
			########################################
			process_entry($fasta_entry);#
			########################################
			}
		$fasta_entry = $fasta_line;
		}else{
		$fasta_entry .= $fasta_line;
		}
	};
close(FASTA_IN);

# do the last entry!:
unless(length($fasta_entry)<=2)
	{
	$fasta_entry =~ s/^>//;
	########################################
	process_entry($fasta_entry);#
	########################################
	}



sub process_entry
	{
	my $line = shift;
	my $current_id = "";
	if ($line =~ /^(.+)\n/ )
		{
		$current_id = $1;
		$line =~ s/^.+\n//;#$line =~ s/\n//;$line =~ s/\r//;my $current_seq_length = length($line);
		}else{
		die "\nerror 5631.\n"
		}

	# $which_taxon_is_this_query_assigned_to{$query_ID} = $tax;

	$current_id =~ s/\s+$//;
	my $query_assigned_to_this = $which_taxon_is_this_query_assigned_to{$current_id};
	if($query_assigned_to_this =~ /\w/)
		{
		$query_assigned_to_this =~ s/\W/_/g;
	#	print "sequence ID $current_id from file $fas_file is assigned taxon $query_assigned_to_this\n";

		open(OUT8, ">>$sorted_fasta_prefix.$query_assigned_to_this") || die "";
		print OUT8 ">$current_id\n$line";
		close OUT8;
		};

	}



};



#####################################################################################################
#
#
#
#####################################################################################################










