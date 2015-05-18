use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'mysql-ensembl-mirror.ebi.ac.uk', # alternatively 'useastdb.ensembl.org'
    -port => 4240,
    -user => 'anonymous'
);

open(INFILE, @ARGV[0]);

# my ( $species, $object_type, $db_type ) = $registry->get_species_and_object_type($stable_id);

# print $species, $object_type, $db_type, "\n";

# my $adaptor = $registry->get_adaptor( $species, $object_type, $db_type );
my $adaptor = $registry->get_adaptor( 'Human', 'Core', 'Translation' );

for my $stable_id (<INFILE>) {
    chomp $stable_id;
    if ($stable_id eq "") {
	next;
    }
    my $translation = $adaptor->fetch_by_stable_id($stable_id);

    if (! $translation) {
	print STDERR "Problem fetching annotations for", $stable_id, "\n";
	next;
    }
    my $pfeatures = $translation->get_all_ProteinFeatures();
    while ( my $pfeature = shift @{$pfeatures} ) {
	my $logic_name = $pfeature->analysis()->logic_name();

	printf(
	    "%s\t%d\t%d\t%s\t%s\t%s\t%s\n",
	    $stable_id,
	    $pfeature->start(), $pfeature->end(), $logic_name,
	    $pfeature->interpro_ac(),
	    $pfeature->hseqname(),
	    $pfeature->idesc()
	    );
    }
}
