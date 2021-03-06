# File : Makefile.PL
use ExtUtils::MakeMaker;
use Config;

# unfortunately MakeMaker thinks it is clever
# with forcing several compiler and/or linker flags
# However, this behavior is a real pain in the a**
# when crosscompiling and therefore we need to adjust
# some things prior to the call of WriteMakefile()

my $CCFLAGS = $Config{'ccflags'};
my $OPTIMIZE = $Config{'optimize'};

# Perl is built with -Wdeclaration-after-statement on RHEL5 - this isn't
# meaningful for C++ - it only emits a warning but it's easy to fix.
$CCFLAGS =~ s/(?:^|\s+)-Wdeclaration-after-statement(?:\s+|$)/ /;

# The generated code causes "variable may be used uninitialized" warnings
# if Perl was built with -Wall.
$CCFLAGS =~ s/(^|\s+)-Wall(\s+|$)/$1-Wall -Wno-uninitialized$2/;

# add CPP/CXX flags if they exist
$CCFLAGS .= ' ' . $var{CPPFLAGS} if exists $var{CPPFLAGS};
$CCFLAGS .= ' ' . $var{CXXFLAGS} if exists $var{CXXFLAGS};

# add CXXFLAGS and AM_CXXFLAGS as shell variables
$CCFLAGS .= ' $(AM_CXXFLAGS) $(CXXFLAGS) -I../../H';

# we also need to strip some compiler flags from 'optimize'
# since architecture features of the host system may interfere
# with the target architecture when crooscompiling
$OPTIMIZE =~ s/-march=[a-zA-Z0-9_\-]+//g;
$OPTIMIZE =~ s/-mtune=[a-zA-Z0-9_\-]+//g;

# print "CCFLAGS:\t", $CCFLAGS, "\n";
# print "OPTIMIZE:\t", $OPTIMIZE, "\n";

WriteMakefile(
              MAKEFILE    =>  "Makefile.perl",
              NAME        =>  "RNA",
              LIBS        =>  ["-lm"],
              CCFLAGS     =>  $CCFLAGS,
              OPTIMIZE    =>  $OPTIMIZE,
              MYEXTLIB    =>  "../../lib/libRNA.a",
              LD          =>  '${CXX}',
              PM          =>  {"RNA.pm", '$(INST_LIBDIR)/RNA.pm'},
              OBJECT      =>  "RNA_wrap.o",
              dynamic_lib =>  {OTHERLDFLAGS => '$(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -fopenmp'},
              AUTHOR      =>  'Ivo Hofacker <ivo@tbi.univie.ac.at>',
              VERSION     =>  '2.1.8'
);
