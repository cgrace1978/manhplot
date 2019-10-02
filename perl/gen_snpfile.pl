use strict;
use Getopt::Long qw(GetOptions);
my %locus=();
my $chrcol="";
my $snpcol="";
my $poscol="";
my $pcol="";
my $betacol="";
my $eafcol="";
my $chridx=0;
my $snpidx=0;
my $posidx=0;
my $pidx=0;
my $betaidx=0;
my $eafidx=0;

my $debug=0;
my $file="";
my $genefile="glist-hg19";#insert the file for gene coordinates depending on the build
my $snpgenedist=500000;#if gene is beyond this distance from SNP either side, the gene annotation is "gene desert"
my $gwasp=5e-8;#pvalue cut-off to define a peak
my $boundary=500000;#distance between 2 peaks to be classified as a locus. 
                    #Can be changed to 0.5 or 1 if interested in cM boundaries 
		    #provided the posidx is pointing to cM positions
my $totalchrs=0;
my $outfile="";
GetOptions(
    'chrom=s' => \$chrcol,
    'var=s'   => \$snpcol,
    'pos=s'   => \$poscol,
    'pval=s'  => \$pcol,
    'beta=s'  => \$betacol,
    'eaf=s'   => \$eafcol,
    'gwasfile=s' => \$file,
    'genefile=s' => \$genefile,
    'gwaspcut=s' => \$gwasp,
    'locusbounds=s' => \$boundary,
    'snpgenebounds=s'=> \$snpgenedist,
    'out=s'   => \$outfile,
    'help' => \$debug,
) or die "Usage: $0 --help\n";

if ($debug || $file eq ""){
print "\n\n";
  print <<HELP;
######################################################################
Help Documentation to generate locus file from GWAS summary statistics
Contact: Anuj Goel
Version: 1.0
Date: 1 Oct 2019
Please cite Manhattan++ software if you use this script.

Usage:
  --gwasfile       : [Required] The input GWAS summary file name (tab delimited). 
                   : Can be a gzipped file. Must have .gz extension.
		   : The input file must have the below 6 columns.
  --chrom          : [Required] Column name containing chromosome numbers (numeric)
  --var            : [Required] Column name containing variant names
  --pos            : [Required] Column name having variant positions (basepairs)
  --pval           : [Required] Column name having variant p-values (not log10 transformed)
  --beta           : [Required] Column name having variant beta values (not Odds Ratios)
  --eaf            : [Required] Column name having variant effect allele frequency
  --genefile       : [Required] The input gene database file name
                   : Please download from PLINK/GitHub. (eg. glist-hg19, glist-hg38)
  --gwaspcut       : The p-value cut-off for locus (default: 5e-8)
  --locusbounds    : The distance between 2 loci (default: 500000)
  --snpgenebounds  : The max distance between SNP & gene to consider (default: 500000)
  --out            : Output file name prefix. (default: gwasfile.script.txt)
  --help           : 
#######################################################################
HELP
print "\n\n";
exit;
}

print "\n\n";
if ($outfile eq ""){
  $outfile="$file.script";
}
my $datestring = localtime();
print "Analysis started: $datestring\n";
print "Parameters:\n";
print "  GWAS Filename: $file\n";
print "  Variant column name: $snpcol\n";
print "  Variant chromosome column name: $chrcol\n";
print "  Variant position column name: $poscol\n";
print "  Variant beta column name: $betacol\n";
print "  Variant allele frequency column name: $eafcol\n";
print "  Variant pvalue column name: $pcol\n";
print "  GWAS threshold: $gwasp\n";
print "  Gene database to use: $genefile\n";
print "  Distance between 2 loci: $boundary\n";
print "  Max. distance between variant and gene to consider: $snpgenedist\n";
print "  Output file: $outfile.txt\n\n\n";


my %orsnps=();
my %eafsnps=();
my %psnps=();
#Take care of gzip input summary statistics file
if ($file=~/\.gz$/){
  print "Opening gzipped file ...\nMight not work if running this script in an OS with no zcat installed\n\n";
  open (IF, "zcat $file |") or die "Cannot open $file\n;exit;";
}else{
  open (IF, $file) or die "Cannot open $file\n";
}


open (OF ,">$outfile.txt") or die "Cannot write output file\n";
print OF "chr markername pos Pvalue OR eaf gene novel\n";
#chr snp pos p beta eaf
my $header=<IF>;
chomp $header;
#figure out column names and column numbers
my @colnames=split(/\t/,$header);
my $colx=0;#should add up to 6 as we need 6 column names
for (my $z=0;$z<=$#colnames;$z++){
  if    ($colnames[$z] eq $chrcol){
    $chridx=$z;
    $colx++;
  }elsif($colnames[$z] eq $snpcol){
    $snpidx=$z;
    $colx++;
  }elsif($colnames[$z] eq $poscol){
    $posidx=$z;
    $colx++;
  }elsif($colnames[$z] eq $pcol){
    $pidx=$z;
    $colx++;
  }elsif($colnames[$z] eq $betacol){
    $betaidx=$z;
    $colx++;
  }elsif($colnames[$z] eq $eafcol){
    $eafidx=$z;
    $colx++;
  }
}
if($colx<6){
  print "$colx ERROR: Not all column names specified. Please use --help for usage.\n";
  exit;
}

print "Reading GWAS summary file ... \n\n";
my %sig=();
my $co=1;
#store information of the lead SNPs only in hash.
my $maxchr=0;
while (my $line=<IF>){
  chomp $line;
  my @cells=split('\t',$line);
  
  if ($cells[$chridx]=~m/\D/){
    print "ERROR: Chromosomes need to be numeric\n";
    exit;
  }elsif($cells[$chridx]>$maxchr){
    $maxchr=$cells[$chridx];
  }
  if ($cells[$pidx]<$gwasp){ 
    $sig{$line}=$cells[$pidx];#store whole line as key. pval as value
    $psnps{$cells[$snpidx]}=$cells[$pidx];
    $eafsnps{$cells[$snpidx]}=sprintf("%.2f", $cells[$eafidx]);
    my $or=exp($cells[$betaidx]);
    $orsnps{$cells[$snpidx]}=sprintf("%.2f", $or);
  }
}

$totalchrs=$maxchr;
close IF;
print "Reading Gene database file ... \n\n";
open (IF, $genefile) or die "Cannot open $genefile\n";
#19 58858171 58864865 A1BG
#19 58863335 58866549 A1BG-AS1
my %genesbychr=();#chr as key and array of gene details as value
my %chr=("X"=>23,"Y"=>24);
while (my $line=<IF>){
  my ($c,$s,$e,$g)=split(/\s+/,$line);
  if ($c eq "X"){
    $c=23;
  }elsif($c eq "Y"){
    $c=24;
  }
  if (exists $genesbychr{$c}){#if chr exists
    my %genes=%{$genesbychr{$c}};#gene hash of genes for this chr
    if (!exists $genes{$g}){#store this gene as its not in the hash
       my @coors=($s,$e);
       $genes{$g}=\@coors;
       $genesbychr{$c}=\%genes;
    }
  }else{#create chr entry in the hash
    my @coors=($s,$e);
    my %genes=();
    $genes{$g}=\@coors;
    $genesbychr{$c}=\%genes;
  }
}

close IF;
print "Identifying loci ... \n";
for (my $chr=1;$chr<=$totalchrs;$chr++){#process chr by chr
  my @sortedchr=();
  my $c=0;
  foreach my $line (sort { $sig{$a} <=> $sig{$b} } keys  %sig){#sort rows by p-values
    my @cells=split('\t',$line);
    #print "$line\n";
    if ($cells[$chridx] == $chr){
      $c++;
      push(@sortedchr,$line);
    }
  }
  #print "Total in this chr=$c\n";
  
  #Do the filtering in the subroutine
  my %list=%{&getloci(@sortedchr)};
  
  #print the output for the chromosome
  foreach my $name (sort { $list{$a} <=> $list{$b} } keys %list){
    my $snppos=$list{$name};
    #get gene info
    my %genelist=%{$genesbychr{$chr}};
    my $dist5=$snpgenedist;
    my $dist3=$snpgenedist;
    my $gene5="ABC";
    my $gene3="DEF";
    my $genein="GEH";
    foreach my $genes (keys %genelist){
      my ($s,$e)=@{$genelist{$genes}};
      if ($snppos>=$s && $snppos<=$e){#simple situation
         $genein=$genes;
         last;
      }elsif($snppos<$s){#get 3' gene
        my $dist=$s-$snppos;
        if ($dist<=$snpgenedist && $dist<=$dist3){#gene within defined distance
          $dist3=$dist;#overwrite this to get the closest gene 3'
          $gene3=$genes;
        }
      }elsif($snppos>$e){#get 5' gene
       my $dist=$snppos-$e;
       if ($dist<=$snpgenedist && $dist<=$dist5){#gene within defined distance
         $dist5=$dist;
         $gene5=$genes;
       }
      }
    }
    #format final gene output
    my $geneout="gene desert";#genes are far apart as defined from SNP
    if ($genein ne "GEH"){#best situation
       $geneout=$genein;
    }elsif($gene5 ne "ABC" && $gene3 ne "DEF"){#found 2 genes either side of SNP & within defined dist
       $geneout="$gene5/$gene3";
    }elsif($gene5 ne "ABC" && $gene3 eq "DEF"){#found only 1 gene 5' of SNP. 3' gene is too far
       $geneout=$gene5;
    }elsif($gene5 eq "ABC" && $gene3 ne "DEF"){#found only 1 gene 3' of SNP. 5' gene is too far
       $geneout=$gene3;
    }

    print OF $chr," ",$name," ",$snppos," ",$psnps{$name}," ",$orsnps{$name}," ",$eafsnps{$name}," $geneout FALSE\n";
    $co++;
  }
  
}#process next chromosome

print "\n************************************************************************************\n\n\n";
print "WARNING: Please remember to update the \"novel\" column in the output file ($outfile.txt)\n\n";
print "\n************************************************************************************\n";
close OF;
$datestring = localtime();
print "Analysis finished: $datestring\n";

exit;

sub getloci(){
  my (@sorted)=@_;
  my $total=scalar @sorted;
  my @sortedtmp=();
  my %list=();
  @sortedtmp=@sorted;
  my $spliced=0;
  my $loci=0;
  my @cells1=();
  
  #keep on going through the array till the time the 
  #number of rows dropped and numb of loci found add 
  #up to the total number of rows in this chromosome
  
  while($spliced+$loci<=$total){
    my $line1=shift @sorted;#the top row is the locus
    @sortedtmp=@sorted;#make a copy as you dont want to iterate the array which you are processing too
    $loci++;#increment locus by 1
    @cells1=split('\t',$line1);
    if (defined $cells1[$snpidx]){
      #print "Y:$cells1[1] $cells1[2]\n";
      $list{$cells1[$snpidx]}=$cells1[$posidx]; 
    }
    #print "XX $loci $cells1[1] $cells1[2]\n";
    my @coord=();
    #iterate through the list of loci to drop rows (cells in array) which are nearby (+/-1 Mb)
    foreach my $key(keys %list){
      my $specific=0;
      #iterate through the rows (@sortedtmp)
      for (my $x=0;$x<=$#sortedtmp;$x++){
        #foreach $line(@sorted){
        my $line=$sortedtmp[$x];
        # print "$line\n";
        my @cells=split('\t',$line);
      
        if (($cells[$posidx]>(($list{$key})-$boundary)) && ($cells[$posidx]<(($list{$key})+$boundary)) && $key ne $cells[$snpidx]){
          #print "Splicing $cells[2] $list{$key} $key $cells[1] ",$cells[2]-$list{$key},"\n";
	  #splice(@sorted,$x,1);
	  push(@coord,$x);
	  $specific++;#keep track of number of rows filtered for this locus
	  $spliced++;#keep track of all rows filtered for this chromosome
        }
      }
      if($specific>=1){#if there are any to be filtered for this locus
        #print "Array before=",scalar @sorted, " ";
        my @tmp=@{&splicearray(\@sorted,\@coord)};
        #print "Splice=$specific, Array after=",scalar @tmp,"\n";
        @sortedtmp=@tmp;#update the arrays
        @sorted=@tmp;#update the arrays
        #print "$loci $spliced $total\n";
      }
    }
  }
  #sort out the last element. 
  my $line1=shift @sorted;
  my @cells1=split('\t',$line1);
  #print "Y:$cells1[1] $cells1[2]\n";
  if (defined $cells1[$snpidx]){  
    $list{$cells1[$snpidx]}=$cells1[$posidx];
  }
  return \%list;
}

#subroutine to delete elements of the array and then tidy up the array
sub splicearray(){
  my ($array,$coord)=@_;
  my @tmp=@{$array};
  foreach my $co(@{$coord}){
    delete $tmp[$co];
  }
  my @new=();
  foreach(@tmp){
    if( ( defined $_) and !($_ =~ /^$/ )){
        push(@new, $_);
    }
  }
  return \@new;
}
