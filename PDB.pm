package Chemistry::File::PDB;

$VERSION = '0.10';

use base "Chemistry::File";
use Chemistry::MacroMol;
use Chemistry::Domain;
use Carp;
use strict;
use warnings;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

=head1 NAME

Chemistry::File::PDB

=head1 SYNOPSIS

    use Chemistry::File::PDB 'pdb_read';

    my $macro_mol = pdb_read("myfile.pdb");

=cut


require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw(  );

@EXPORT_OK = qw( pdb_read );

%EXPORT_TAGS = (
   all  => [@EXPORT, @EXPORT_OK]
);

=head1 DESCRIPTION

This module reads PDB files. The PDB file format is commonly used to describe
proteins, particularly those stored in the Protein Data Bank
(L<http://www.rcsb.org/pdb/>). The current version of this module only reads
the ATOM records, ignoring everything else.

This module automatically registers the 'pdb' format with Chemistry::Mol,
so that PDB files may be identified and read by Chemistry::Mol::read_mol().

The PDB reader returns a Chemistry::MacroMol object, by default, but the
user is free to give any object to the pdb_read subroutine, as long as it
implements the same interface as Chemistry::Mol.

=cut

Chemistry::Mol->register_format(pdb => __PACKAGE__);

=head1 FUNCTIONS

=over 4

=item $my_mol = pdb_read($fname, option => value...)

Returns a Chemistry::MacroMol object (by default) from the specified PDB file.
The only option so far is mol => $my_mol, which would use a previously created
object instead of creating a new MacroMol object. The object should be a 
Chemistry::Mol object or a derived class.

=cut

sub parse_file {
    my $class = shift;
    my $fname = shift;
    my %options = @_; # a molecule
    my @mols; 
    my ($n_mol, $n_atom);
    my $n_res = 0;
    my $domain;

    open F, $fname or croak "Could not open file $fname";

    my $mol = $options{mol} || Chemistry::MacroMol->new(id => "mol". ++$n_mol);
    my $is_macro = $mol->isa('Chemistry::MacroMol');
    while (<F>) {
	if (/^TER/) {
#	    $mol->{name} = $name;
	    #push @mols, $mol;
	    #$mol = new Chemistry::Mol(id => "mol". ++$n_mol);
	    #$n_atom = 0;
	} elsif (/^ATOM/) {
	    my ($symbol, $suff, $res_name, $seq_n, $x, $y, $z) = 
		unpack "x12A2A2x1A3x2A4x4A8A8A8", $_;
	    #print "S:$symbol; N:$name; x:$x; y:$y; z:$z\n";
            $seq_n =~ s/ //g;
            if ($seq_n != $n_res) {
                $domain = Chemistry::Domain->new(parent=>$mol, name=>$res_name,
                    type => $res_name, id => "d".$seq_n);
                $n_res = $seq_n;
                $mol->add_domain($domain);
            }
            my $atom_name = $symbol.$suff;
            $atom_name =~ s/ //g;
	    my $a = $domain->new_atom(
		symbol => $symbol, 
		coords => [$x, $y, $z], 
		#id    => "$mol->{id}-$res_name-a".++$n_atom,
		id    => "a".++$n_atom,
                name => $atom_name,
	    );
            $a->attr('pdb/residue', $domain->name.$seq_n);
	}
    }
    close F;

    return $mol;
}

=item is_pdb($fname)

Returns true if the specified file is a PDB file. The test is not incredibly
smart; it assumes that a file is a PDB if the name ends with .pdb or if it has
any line beginning with "ATOM  ".

=cut

sub file_is {
    my $class = shift;
    my $fname = shift;
    
    return 1 if $fname =~ /\.pdb$/i;

    open F, $fname or croak "Could not open file $fname";
    
    while (<F>){
	if (/^ATOM  /) {
	    close F;
	    return 1;
	}
    }

    return 0;
}

1;


=back

=head1 SEE ALSO

L<Chemistry::MacroMol>, L<Chemistry::Mol>

The PDB format description at 
L<http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html>

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>

=cut

