package Chemistry::File::PDB;

$VERSION = '0.15';
# $Id: PDB.pm,v 1.8 2004/06/18 19:45:04 itubert Exp $

use base qw(Chemistry::File);
use Chemistry::MacroMol;
use Chemistry::Domain;
use Carp;
use strict;
use warnings;

=head1 NAME

Chemistry::File::PDB

=head1 SYNOPSIS

    use Chemistry::File::PDB;

    my $macro_mol = Chemistry::MacroMo->read("myfile.pdb");

=cut

=head1 DESCRIPTION

This module reads PDB files. The PDB file format is commonly used to describe
proteins, particularly those stored in the Protein Data Bank
(L<http://www.rcsb.org/pdb/>). The current version of this module only reads
the ATOM and HETATM records, ignoring everything else.

This module automatically registers the 'pdb' format with Chemistry::Mol,
so that PDB files may be identified and read by Chemistry::Mol->read(). For 
autodetection purpuses, it assumes that files ending in .pdb or having 
a line matching /^(ATOM  |HETATM)/ are PDB files.

The PDB reader is designed for generating Chemistry::MacroMol objects, but
it can also create Chemistry::Mol objects by throwing some information away.

=head2 Properties

When reading the file, this PDB reades stores some of the information in the
following places:

=over

=item $domain->type

The residue type, such as "ARG".

=item $domain->name

The type and sequence number, such as "ARG114". The system doesn't deal with
chains yet.

=item $atom->name

The PDB atom name, such as "CA".

=item $atom->attr("pdb/residue_name")

The name of the residue, as discussed above.

=back

=cut

Chemistry::Mol->register_format(pdb => __PACKAGE__);

sub parse_file {
    my $class = shift;
    my $fname = shift;
    my %options = @_;
    my @mols; 
    my ($n_mol, $n_atom);
    my $n_res = 0;
    my $domain;

    open F, $fname or croak "Could not open file $fname";

    my $mol_class = $options{mol_class} || "Chemistry::MacroMol";
    my $mol = $mol_class->new(id => "mol". ++$n_mol);
    my $is_macro = $mol->isa('Chemistry::MacroMol');
    $domain = $mol unless $is_macro;
    while (<F>) {
	if (/^TER/) {
	    #$mol->{name} = $name;  # create multiple molecules
	    #push @mols, $mol;
	    #$mol = new Chemistry::Mol(id => "mol". ++$n_mol);
	    #$n_atom = 0;
	} elsif (/^(HETATM|ATOM)/) {
	    my ($symbol, $suff, $res_name, $seq_n, $x, $y, $z) = 
		unpack "x12A2A2x1A3x2A4x4A8A8A8", $_;
	    #print "S:$symbol; N:$name; x:$x; y:$y; z:$z\n";
            $seq_n =~ s/ //g;
            if (!$domain || $seq_n != $n_res) {
                if ($is_macro) {
                    $domain = Chemistry::Domain->new(
                        parent => $mol, name => "$res_name$seq_n",
                        type => $res_name, id => "d".$seq_n);
                    $mol->add_domain($domain);
                }
                $n_res = $seq_n;
            }
            my $atom_name = $symbol.$suff;
            $atom_name =~ s/ //g;
	    my $a = $domain->new_atom(
		symbol => $symbol, 
		coords => [$x, $y, $z], 
		id    => "a".++$n_atom,
                name => $atom_name,
	    );
            $a->attr('pdb/residue_name',    "$res_name$seq_n");
            $a->attr('pdb/sequence_number', $n_atom);
	}
    }
    close F;

    return $mol;
}

sub file_is {
    my $class = shift;
    my $fname = shift;
    
    return 1 if $fname =~ /\.pdb$/i;

    open F, $fname or croak "Could not open file $fname";
    
    while (<F>){
	if (/^ATOM  / or /^HETATM/) {
	    close F;
	    return 1;
	}
    }

    return 0;
}

1;

=head1 VERSION

0.15

=head1 SEE ALSO

L<Chemistry::MacroMol>, L<Chemistry::Mol>, L<http://www.perlmol.org/>.

The PDB format description at 
L<http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html>

=head1 AUTHOR

Ivan Tubert-Brohman <itub@cpan.org>

=cut

