#!/usr/bin/perl

# TGMM2SVF_XMLfinalResult_folder_fix_cell_NaNs.pl 
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# This script is used to fix the TGMM XML data to alleviate the following error during SVF:
#	Traceback (most recent call last):
#	  File "/home/martin/Downloads/TGMM2SVF.app/SVF-prop.py", line 627, in <module>
#	    VF = build_VF_propagation_backward(LT_main, t_b = te, t_e = tb, neighb_size = 10, nb_proc=24)
#	  File "/home/martin/Downloads/TGMM2SVF.app/SVF-prop.py", line 99, in build_VF_propagation_backward
#	    C_LT = to_check_LT[idx3d.query(LT.VF.pos[C])[1]]
#	IndexError: list index out of range

# Please "cd" to make working directory XML_finalResult_lht; make a backup of all XML files!
# Then run "perl /path/to/the/script/TGMM2SVF_XMLfinalResult_folder_fix_cell_NaNs.pl" as many times as needed until the error goes away

use Cwd qw( cwd );
my $path = Cwd::cwd(); #. "/interestpoints";

#==================
#Main subroutine
#==================
sub main {
	#my $z_column = 3;#default
	#my $id_column = 0; #default
	
	if ( opendir( my $dh, $path ) ) {
		my @files_ = readdir( $dh );
		@files_ = sort { $a cmp $b } @files_;
		
		#my @z_values = ();
		#my @z_values_transformed = ();
		my @timepoints_;
		my @files_to_fix;
		my @timepoints_to_fix;
		my $html, $html1;
		my $cell_text, $blank_text;#, $new_cell_text;
		my @split_cell_text;
		#my @altered_text = ();
		my @altered_index;
		
		for ( my $i=scalar(@files_)-1; $i>=0; $i-- ) {
			if ( $files_[$i] =~ /^GMEMfinalResult_frame(\d+)\.xml$/ && open(FILE, "<$path/$files_[$i]" ) ) {
				my $tp = $1;
				$timepoints_[$i] = $tp;
				flock( FILE, LOCK_EX );
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /nan/i ) {
						push( @files_to_fix, $files_[$i] );
						push( @timepoints_to_fix, $tp );
						last;
					}
				}
				flock( FILE, LOCK_UN );
				close( FILE );
			}
		}
		#print scalar(@files_to_fix) . "\n";
		for ( my $i=0; $i<scalar(@files_to_fix); $i++ ) {
			if ( open(FILE, "<$path/$files_to_fix[$i]" ) ) {
				flock(FILE, LOCK_EX);
				$_ = $/;
				$/ = undef;
				$html = <FILE>; 
				flock(FILE, LOCK_UN);
				close(FILE);
				$/ = $_;
			} else {
				print "Error opening $path/$files_to_fix[$i]!\n";
				return;
			}
			
			my @orig_cells;
			my @replace_cells;
			
			while( $html =~ /<GaussianMixtureModel (.*?)>(.*?)<\/GaussianMixtureModel>/igs ) {
				$cell_text = $1;
				$blank_text = $2;
				
				next unless ( $cell_text =~ /nan/i );
				
				my @this_cell_data = split( /\ (?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)/, $cell_text ); #parse actual key/value pairs, in plaintext
				my @parent_cell_data;
				my @child_cell_data;
				my $this_id = undef;
				my $this_parent_id = undef;
				#@altered_text = ();
				@altered_index = ();
				
				for ( my $q=0; $q<scalar(@this_cell_data); $q++ ) {
					if ( $this_cell_data[$q] =~ /id=\"(\w+)\"/i ) {
						$this_id = $1;
					} elsif ( $this_cell_data[$q] =~ /parent=\"(.*?)\"/i ) {
						$this_parent_id = $1
					} elsif ( $this_cell_data[$q] =~ /nan/i ) {
						push( @altered_index, $q );
						#push( @altered_text, $this_cell_data[$q] )
					}
				}
				#print "NaN encountered in $files_to_fix[$i] $this_parent_id $this_id\n";
				next unless ( scalar(@altered_index) > 0 && $this_id =~ /\w/ && $this_parent_id =~ /\w/ );
				
				#gather child and parent data here
				for ( my $qq=-1; $qq<=1; $qq++ ) { #-1 is minus one timepoint, or parent; +1 is child
					next if ( $qq==0);
					next if ( $qq< 0 && $timepoints_to_fix[$i] <= 0 ); #don't gather parent if we are at first timepoint
					next if ( $qq > 0 && $timepoints_to_fix[$i] >= $#timepoints_ ); #don't gather child if we are at last timepoint

					my $find_file = $timepoints_to_fix[$i] + $qq;
					#print "NaN encountered in $files_to_fix[$i]:$qq\n";
					my $found = -1;
					for ( my $q=0; $q<scalar(@timepoints_); $q++ ) {
						if ( $find_file == $timepoints_[$q] ) {
							$found = $q;
							last;
						}
					}
					
					if ( $found >= 0 ) {
						#print "NaN encountered in $files_to_fix[$i]:$files_[$found]\n";
						if ( open(FILE, "<$path/$files_[$found]" ) ) {
							flock(FILE, LOCK_EX);
							$_ = $/;
							$/ = undef;
							$html1 = <FILE>; 
							flock(FILE, LOCK_UN);
							close(FILE);
							$/ = $_;
						} else {
							print "Error opening $path/$files_[$found]!\n";
							return;
						}
						
						$found = -1;
						while( $html1 =~ /<GaussianMixtureModel (.*?)>.*?<\/GaussianMixtureModel>/igs ) {
							#$new_cell_text = $1;
							#print "here";
							my @new_cell_data = split( /\ (?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)/, $1 ); #parse actual key/value pairs, in plaintext
							for ( my $q=0; $q<scalar(@new_cell_data); $q++ ) {
								if ( $new_cell_data[$q] =~ /id=\"(\w+)\"/i ) {
									if ( $qq<0 && $this_parent_id eq $1 ) {
										#found parent
										@parent_cell_data = @new_cell_data;
										$found = 1;
									}
								} elsif ( $new_cell_data[$q] =~ /parent=\"(\w+)\"/i ) {
									if ( $qq>0 && $this_id eq $1 ) {
										#found parent
										@child_cell_data = @new_cell_data;
										$found = 1;
									}
								} elsif ( $new_cell_data[$q] =~ /nan/i ) {
									#error here
									$found = -1;
									last;
								}
							}
							
							last if ( $found > 0 );
						}
					}
				}
				
				#print "NaN encountered in $files_to_fix[$i]: " . scalar(@parent_cell_data) . "," . scalar(@child_cell_data) . " $parent_cell_data[0],$child_cell_data[0]\n";
				#now, iterate over each field needing to be altered
				#print "NaN encountered in $files_to_fix[$i]: " . scalar(@altered_index) . "\n";
				for ( my $k=0; $k<scalar(@altered_index); $k++ ) {
					#print "here";
					if ( $this_cell_data[$altered_index[$k]] =~ /(\w+)=\"(.*)\"/i ) { #print "here";
						my $find_field = $1;
						#my $find_values = $2;
						my $child_data = "";
						my $parent_data = "";
						
						for ( my $l=0; $l<scalar(@child_cell_data); $l++ ) {
							if ( $child_cell_data[$l] =~ /$find_field=\"(.*)\"/i ) {
								$parent_data = $1;
								last;
							}
						}
						for ( my $l=0; $l<scalar(@parent_cell_data); $l++ ) {
							if ( $parent_cell_data[$l] =~ /$find_field=\"(.*)\"/i ) {
								$parent_data = $1;
								last;
							}
						}
						
						if ($parent_data eq "" && $child_data eq "") {
							#no data at all!
							print "Problem collecting parent and/or child data for cell id $this_id in file $files_to_fix[$i]!\n";
						} elsif ( $parent_data eq "" ) {
							#no parent data, just copy child's
							$this_cell_data[$altered_index[$k]] = $find_field . "=\"" . $child_data . "\"";
						} elsif ( $child_data eq "" ) {
							#no child data, just copy parent's
							$this_cell_data[$altered_index[$k]] = $find_field . "=\"" . $parent_data . "\"";
						} elsif ( $parent_data =~ /\s/ || $child_data =~ /\s/ ) { #compound field
							my @parent_data_values = split( /\s/, $parent_data );
							my @child_data_values = split( /\s/, $child_data );
						
							if ( scalar(@parent_data_values) != scalar(@child_data_values)  ) {
								print "Problem with mismatched parent and child data for cell id $this_id in file $files_to_fix[$i]!\n";
							}
							my $final_data_field = "";
							for ( my $m=0; $m<scalar(@parent_data_values); $m++ ) {
								$final_data_field .= (($parent_data_values[$m]+$child_data_values[$m])/2) . " ";
							}
							
							$this_cell_data[$altered_index[$k]] = $find_field . "=\"" . $final_data_field . "\"";
						} else {
							if ( $child_data =~ /\d/ && $parent_data =~ /\d/ ) {
								#numerical value, so average
								$this_cell_data[$altered_index[$k]] = $find_field . "=\"" . (($parent_data+$child_data)/2) . "\"";		
							} else {
								#just copy parent data
								$this_cell_data[$altered_index[$k]] = $find_field . "=\"" . $parent_data . "\"";
							}
						}
					}
				}
				
				push( @orig_cells, "<GaussianMixtureModel $cell_text>$blank_text<\/GaussianMixtureModel>" );
				push( @replace_cells, "<GaussianMixtureModel " . join(' ', @this_cell_data) . ">$blank_text<\/GaussianMixtureModel>" );
			}
			
			
			
			#next; #TODO:comment this line
			for ( my $i=0; $i<scalar(@orig_cells); $i++ ) {
				#print "replacing $orig_cells[$i]\n\nwith\n\n$replace_cells[$i]\n\n\n\n";
				$html =~ s/$orig_cells[$i]/$replace_cells[$i]/is;			
			}
			
			rename( "$path/$files_to_fix[$i]", "$path/$files_to_fix[$i]." . time() . ".bak" );
			if ( open(FILE, ">$path/$files_to_fix[$i]" ) ) {
				flock(FILE, LOCK_EX);
				print FILE $html;
				flock(FILE, LOCK_UN);
				close(FILE);
			} else {
				print "Error writing to $path/dataset\.xml!\n";
				return;
			}
		}
	

	} else {
		print "Can't perform opendir $path: $!\n";
	}	
}

#==================
#Array remove duplicates
#==================
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

sub set_stats {

	my $avg = 0;
	my $std = 0;
	my $med = 0;
	#my $mdn;
	#use Data::Dumper;	print Dumper( \@_ );
	for ( my $i = 0; $i < scalar(@_); $i++ ) {
		#find the outliers of theta, and the average theta -> the average/outlier differential will determine how excessively to adjust the exponent
		$avg += @_[$i];#abs(@_[$i]);
	}
	$avg /= scalar(@_);
	for ( my $i = 0; $i < scalar(@_); $i++ ) {
		#find the outliers of theta, and the average theta -> the average/outlier differential will determine how excessively to adjust the exponent
		$std += (@_[$i] - $avg) ** 2;
	}
	$std = sqrt( $std / (scalar(@_)-1) );
	
	#TODO: fix this median-determining function!
	@_ = sort { $a <=> $b } @_;#{ $a cmp $b };
	#@_ = sort(@_);#{ $a cmp $b };
	if ( scalar(@_) & 1 ) { #odd number of numbers
		$med = int( scalar(@_)/2);
		$med = @_[$med];
		#$med = "me!";
	} else {
		$med = scalar(@_) / 2;
		$med = ( @_[$med] + @_[$med -1] ) / 2;
	}
	
	return ( $avg, $med, $std );
}

sub return_counts_statistics {

	#this sub does the main heavy lifting of this script: it finds the counts for every item in an array, then returns the avg,std,med,mode,count_of_mode of those counts
	my @unique = uniq(@_);
	my @unique_counts = ( 0 ) x scalar(@unique);

	my $max_index = 0;
	
	for ( my $i=0; $i<scalar(@unique); $i++ ) {
		#$unique_counts[$i] = 0;
		map { $unique_counts[$i]++ if ( $unique[$i] eq $_ ); } @_;
			$max_index = $i if ( $unique_counts[$i] > $unique_counts[$max_index] );
		}
			
	my ( $avg, $med, $std ) = set_stats(@unique_counts); #find the average number of counts and the standard deviation
	
	return ($avg,$std,$med,$unique[$max_index],$unique_counts[$max_index]);
}




main();


