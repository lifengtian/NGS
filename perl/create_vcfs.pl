use strict;
use warnings;
use YAML qw/Dump Load/;

use CAG::SampleSheet;

open my $fh, '<', 'config.yml' 
  or die "can't open config file: $!";
my $yml = do { local $/; <$fh> };

my $config = Load($yml);

print join("\n",@{$config->{'projects'}}),"\n";

my $HISEQ=$ENV{'HISEQ'};


#find all samples for projects
my $sheet=CAG::SampleSheet->new();
my $SIDs=$sheet->find_samples_by_projects($config->{'projects'});



#retrieve BAMs for each sample
#system("bash $HISEQ/NGS/shell/hiseq_retrieve_SID.sh @{$SIDs}");


open GATK,">gatk.sh" or die "Error open gatk.sh";
print GATK "bamlist=(".join(" ",@{$SIDs}).")\n";
print GATK "VQSR=$config->{'VQSR'}\n";
print GATK "run=$config->{'run_name'}\n";
print GATK "target=$config->{'target'}\n";
print GATK "platform=$config->{'platform'}\n";

system("bash $HISEQ/NGS/shell/hiseq_gatk_resolute.sh @{$SIDs}");
#call variants
#system("bash $HISEQ/NGS/shell/hiseq_genome_resolute.sh @{$SIDs}");

#annotate
#system("bash $HISEQ/NGS/shell/hiseq_annovar.sh");
