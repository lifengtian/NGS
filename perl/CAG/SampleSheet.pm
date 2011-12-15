package CAG::SampleSheet;

use strict;
use warnings;
use Text::CSV;

sub new {
my $thing = shift;
my $self = {};
bless $self, ref($thing) || $thing;
$self->init(@_);
return $self;
}



sub init {
my $self = shift;
# Extract various interesting things from
# @_ and use them to create a data structure
# that represents a customer.
}



sub validate {
my $self = shift;
# Call a number of methods, each of which validates
# one data item in the customer record.
return $self->is_valid_sales_ref
&& $self->is_valid_other_attr
&& $self->is_valid_another_attr;
}


=head1 find_samples_by_projects 
right now, we retrieve from files
in near future, we will do it from a rdbms
=cut

sub find_samples_by_projects {
	my ($self, $projects) = @_ ;
	
	my $HISEQ=$ENV{'HISEQ'};
	my $CSVparser=Text::CSV->new();
	my %h; #accumulate SIDs
	if ( ! -d $HISEQ ) {
		die "$HISEQ not there!";
	}

	opendir my $dir, "$HISEQ/SampleSheet/" or die "Error open $HISEQ/SampleSheet/";
	my @csv=grep(/\.csv$/, readdir($dir));
	closedir($dir);

	for my $csv (@csv){
		open my $fh , "$HISEQ/SampleSheet/$csv" or die "Error open $HISEQ/SampleSheet/$csv";
		<$fh>;
		while(<$fh>){
			 if ( $CSVparser->parse($_) ) {
                		my @cols=$CSVparser->fields();	
				push @{$h{$cols[9]}}, $cols[2];
			 	#print $cols[9],"\n";#print join("\t",@cols),"\n";	
			}
		}
		close $fh;
	}
	my %sids = ();
	for my $proj (@{$projects}){
		#print $proj,"\t",join(";",@{$h{$proj}}),"\n";
		foreach ( @{$h{$proj}} ) {		
			$sids{$_}++;
		}
	}
	return [keys %sids];
}




1;
