use strict;
use warnings;
use YAML qw/Dump Load/;

use CAG::SampleSheet;

my $usage = "perl create_vcfs.pl /path/to/config.yml";
 
if ( @ARGV != 1 ) {
    die $usage;
}

my $path=shift @ARGV;

open my $fh, '<', "$path/config.yml" 
  or die "can't open config file: $path/config.yml $!";
my $yml = do { local $/; <$fh> };

my $config = Load($yml);

#print join("\n",@{$config->{'projects'}}),"\n";

my $HISEQ=$ENV{'HISEQ'};

die "$HISEQ doesn't exist. Exit!" if ( ! -d $HISEQ );

#find all samples for projects
my $sheet=CAG::SampleSheet->new();
my $SIDs=$sheet->find_samples_by_projects($config->{'projects'});


#retrieve BAMs for each sample
#system("bash $HISEQ/NGS/shell/hiseq_retrieve_SID.sh @{$SIDs}");


open GATK,">$path/gatk.sh" or die "Error open $path/gatk.sh";
print GATK "bamlist=(".join(" ",@{$SIDs}).")\n";
print GATK "VQSR=$config->{'VQSR'}\n";
print GATK "run=$config->{'run_name'}\n";
print GATK "target=$config->{'target'}\n";
print GATK "platform=$config->{'platform'}\n";

system("bash $HISEQ/NGS/shell/hiseq_gatk_resolute.sh $p");
#call variants
#system("bash $HISEQ/NGS/shell/hiseq_genome_resolute.sh @{$SIDs}");

#annotate
#system("bash $HISEQ/NGS/shell/hiseq_annovar.sh");
