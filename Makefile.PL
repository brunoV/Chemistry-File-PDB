use 5.006001;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'		=> 'Chemistry::File::PDB',
    'VERSION_FROM'	=> 'PDB.pm', # finds $VERSION
    'PREREQ_PM'		=> {
        'Chemistry::MacroMol' => 0.05,
        'Chemistry::Mol' => '0.10',
    }, # e.g., Module::Name => 1.1
);
