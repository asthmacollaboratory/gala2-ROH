#!/usr/bin/perl -w 

use strict;

if($#ARGV != 2){
	die "Usage: ./roh_overlap.pl <roh bed file> <lanc bed file> <chr:pos query>";
}

my $rohbedfile = $ARGV[0];
my $lancbedfile = $ARGV[1];
my $querysite = $ARGV[2];

my $qchr = "";
my $qpos = "";

if($querysite =~ m/((chr)?\d+):(\d+)/){
	if(! defined $2){
		$qchr = "chr";
	}
	$qchr .= $1;
	$qpos = $3;
}
else{
	print STDERR "ERROR: Did not recognize query. Please provide a genomic location formatted chr:pos, e.g. chr1:12345";
}

my @indlist;
my %output;

open(FIN,"<",$rohbedfile) or die $!;
my $ind;
my $pop;
my $n = -1;

my $mstart = -1;
my $mend = 3000000000;
while(defined(my $line = <FIN>)){
	chomp $line;
	if ( $line =~ m/^track .+Ind: (.+) Pop:(.+) ROH.+/ ) {
        if($n == 0){
        	$output{$ind} = "$ind\t$pop\tNA\tNA\tNA\tNA";
        }
        $ind = $1;
        push(@indlist, $ind);
        $pop = $2;
        $n = 0;
    }
    else {
        my ( $chr, $start, $end, $class, $size, $anc, $junk ) = split( /\s+/, $line, 7 );
        if($chr ne $qchr){
        	next;
        }
        if($qpos >= $start and $qpos <= $end){
        	$output{$ind} = "$ind\t$pop\t$chr\t$start\t$end\t$class";

        	if($start >= $mstart){
        		$mstart = $start;
        	}
        	if($mend >= $end){
        		$mend = $end;
        	}

        	$n++;
        }
    }
}

if($n == 0){
	$output{$ind} = "$ind\t$pop\tNA\tNA\tNA\tNA";
}

open( LANC, "<", $lancbedfile ) or die $!;

$n = -1;

my $nind = 0;
my %lancoutput;
my $bpoverlap = "NA";
my $ancoverlap = "NA";
my $overlap;

#for $ind (@indlist){
#	$lancoutput{$ind} = "$mstart\t$mend";
#}

while ( my $line = <LANC> ) {
    chomp $line;

    if ( $line =~ m/^track .+Ind: (.+) Pop:(.+) Admixture.+/ ) {
        if($n == 0){
        	$output{$ind} .= "\tNA";
        }

        if($nind > 0){
        	$lancoutput{$ind} = "\t$bpoverlap\t$ancoverlap"; 
        }

        $ind = $1;
        $pop = $2;

        $nind++;
        $n = 0;
        $bpoverlap = "NA";
        $ancoverlap = "NA";
        $overlap = 0;
    }
    else {
        my ( $chr, $start, $end, $anc, $size, $junk ) = split( /\s+/, $line, 6 );
        if($chr ne $qchr){
        	next;
        }
        $anc =~ tr/123/EAN/;
        if($qpos >= $start and $qpos <= $end){
        	if(exists $output{$ind}){

        		my $eur = 0;
        		my $afr = 0;
        		my $nam = 0;

        		if($anc eq "EE"){
        			$eur = 2;
        		}
        		elsif($anc eq "AA"){
        			$afr = 2;
        		}
        		elsif($anc eq "NN"){
        			$nam = 2;
        		}
        		elsif($anc eq "EN"){
        			$eur = 1;
        			$nam = 1;
        		}
        		elsif($anc eq "EA"){
        			$eur = 1;
        			$afr = 1;
        		}
        		elsif($anc eq "AN"){
        			$afr = 1;
        			$nam = 1;
        		}

        		$output{$ind} .= "\t$anc\t$eur\t$afr\t$nam";
        		$n++;
        	}
        }

        $overlap = getOverlap($start,$end,$mstart,$mend);

        if($overlap > 0 and $bpoverlap eq "NA"){
        	$bpoverlap = $overlap;
        	$ancoverlap = $anc;
        }
        elsif($overlap > 0){
        	$bpoverlap .= ",$overlap";
        	$ancoverlap .= ",$anc";
        }
    }
}

if($n == 0){
	$output{$ind} .= "\tNA";
}
if($nind > 0){
	$lancoutput{$ind} = "\t$bpoverlap\t$ancoverlap"; 
}

close(LANC);

print "pop\tchr\trohStart\trohEnd\tclass\tqueryAnc\tqueryEurAlleles\tqueryAfrAlleles\tqueryNamAlleles\tmaxSpanStart\tmaxSpanEnd\tancOverlap\tancTypes\n";

for $ind (@indlist){
	if(exists $output{$ind} and exists $lancoutput{$ind}){
		print $output{$ind}, "\t$mstart\t$mend\t", $lancoutput{$ind}, "\n";
	}
}


sub getOverlap{
	my $start = $_[0];
	my $end = $_[1];
	my $qstart = $_[2];
	my $qend = $_[3];

	if ( $qstart >= $start and $qstart <= $end ) {
    	if($qend >= $end){
    		return ($end-$qstart)+1;
    	}
    	else{
    		return ($qend-$qstart)+1;
    	}
    }
    elsif ( $qend >= $start and $qend <= $end ) {
        return ($qend-$start)+1;
    }
    elsif ( $qstart <= $start and $qend >= $end ) {
        return ($start-$end)+1;
    }

    return 0;

}
