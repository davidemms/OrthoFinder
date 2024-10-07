#!/usr/bin/perl

#####################################################################
# Author:  KM Amada (kmamada@ifrec.osaka-u.ac.jp)
#
# Ver. Date      Changelog
#####################################################################
# 1.0  07.26.13  Initial release
# 2.0  09.03.13  Added extensive warnings and error messages
# 3.0  10.28.13  Fix for retrieving large files. Added STDERR logs
# 3.1  11.08.13  Added LWP failsafe. Made hat3 not a required output
# 3.2  12.08.14  Removed 5-char restriction for own structure files
#
#####################################################################

use strict;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use LWP::Simple;
use LWP::UserAgent;

# to prevent error 'Header line too long (limit is 8192)' [v3.1]
use LWP::Protocol::http;
push(@LWP::Protocol::http::EXTRA_SOCK_OPTS, MaxLineLength => 0);



my $BASEURL = "http://sysimm.ifrec.osaka-u.ac.jp/MAFFTash/REST/service.cgi/premafft";

my ( $WORKDIR, $PDBLIST, $OWNLIST, $HAT3FILE, $INSTRFILE );

GetOptions
(
    'd=s' => \$WORKDIR,
    'p=s' => \$PDBLIST,
    'o=s' => \$OWNLIST,
    'h=s' => \$HAT3FILE,
    'i=s' => \$INSTRFILE,
);

print STDERR "[MAFFTash-premafft]\n";

# set temp directory
my $TMP = "/tmp/mapremafft$$";
make_path($TMP) unless -d $TMP;



######
# validation
&help("Required parameter : atleast one of either '-p' or '-o'") unless ( defined $PDBLIST || defined $OWNLIST);
&help("Required parameter : '-d'") if defined $OWNLIST && ! defined $WORKDIR;

$HAT3FILE = "hat3" unless defined $HAT3FILE;
$INSTRFILE = "instr" unless defined $INSTRFILE;
chop $WORKDIR if defined $WORKDIR && $WORKDIR =~ m/\/$/g;


######
# prepare inputs
print STDERR "Preparing inputs for service request...\n";

my @files = ();
push(@files, "strweight" => "0.5");
push(@files, "premafft" => "1");


# pdb entries
if ( defined $PDBLIST )
{
    print STDERR "PDB List defined!\n";
    &bail("Error: Input file $PDBLIST does not exists!") unless -e $PDBLIST;
    my $listfile = "$TMP/pdblist.inp";


    open(INPF,"<$PDBLIST") or &bail("Error: Cannot open file $PDBLIST for reading!");
    open(OUTF,">$listfile") or &bail("Error: Cannot open temporary file $listfile for writing!");

    while(<INPF>)
    {
        chomp;
        if ( /^(\w{5})$/ )
        {
            print OUTF ">PDBID\n$1\n";
        }
    }

    close OUTF;
    close INPF;

    push(@files, "inputfile" => ["$listfile"]);
}



# upload own structures
my %ownids = ();

if ( defined $OWNLIST  )
{
    print STDERR "OWN List defined!\n";
    &bail("Error: Input file $OWNLIST does not exists!") unless -e $OWNLIST;


    open(OWNINPF,"<$OWNLIST") or &bail("Error: Cannot open file $OWNLIST for reading!");

    while(<OWNINPF>)
    {
        chomp;

        if ( /^(\S+)$/ )
        {
            my $fileref = "$WORKDIR/$1.pdb";

            unless (-e $fileref)
            {
                close OWNINPF;
                &bail("Error: File $fileref does not exists!");
            }

            push(@files, "inputownfile[]" => ["$fileref"]);
            $ownids{$1} = 1;
        }
    }

    close OWNINPF;
}



######
# start rest service
print STDERR "Sending service request...\n";

my $browser = LWP::UserAgent->new;
$browser->timeout(0);


# post: running a mafftash job
my $postResponse = $browser->post( $BASEURL, \@files, 'Content_Type' => 'form-data' );
&bail(sprintf("[%d] %s\n", $postResponse->code, &parseError($postResponse->content))) unless($postResponse->is_success);


# get response from post request
my ($status, $mafftashid) = &parseResponse($postResponse->content);



my $MAXTRIES = 3;
my $STIMER = 4;
my $longtimer = 0;

print STDERR "Request sent! Waiting for response...[$mafftashid]\n";


# wait for results until it becomes available
while(1)
{
    $longtimer = $longtimer <= ($STIMER*3) ? $longtimer+$STIMER : $STIMER;
    sleep $longtimer;


    # get: get results for mafftash job
    my $getResponse = $browser->get("$BASEURL/$mafftashid");

    if ( $getResponse->is_success )
    {

        # get response from get request
        ($status, $mafftashid) = &parseResponse($getResponse->content);
        next unless ( $status eq "done" );


        # if job is finished and ready
        print STDERR "Results found!\n";
        my $csfile = "$TMP/checksum.tar.gz";
        my $try1 = 1;


        while(1)
        {
            print STDERR "Fetching Results... [Trial $try1]\n";

            if ( is_success(getstore("$BASEURL/getmdlist/$mafftashid", $csfile)) && -e $csfile && -s $csfile )
            {
                # get response from get request
                my $checklist = &extractchecksum($csfile);
                &bail("Error retrieving list of compressed files!") unless ( scalar %$checklist > 0 );


                foreach my $id ( keys %$checklist )
                {
                    my $checkfile = "$TMP/$id";
                    my $checkid = $checklist->{$id};
                    my $try2 = 1;

                    while(1)
                    {
                        unlink $checkfile if -e $checkfile;

                        if ( is_success(getstore("$BASEURL/get/$mafftashid/$id", $checkfile)) && -e $checkfile && -s $checkfile )
                        {
                            my $hashid = &getchecksum($checkfile);
                            #print STDERR "[hashid]$hashid [checkid]$checkid\n";

                            if ($hashid ne "" && $hashid ne $checkid )
                            {
                                unlink $checkfile if -e $checkfile;
                                &bail("Error retrieving compressed file from server! [Checksum Failed]") if $try2 >= $MAXTRIES;
                                $try2++;
                                sleep $STIMER;
                            }
                            else
                            {
                                last;
                            }
                        }
                        else
                        {
                            &bail("Error retrieving compressed file from server!") if $try2 >= $MAXTRIES;
                            $try2++;
                            sleep $STIMER;
                        }
                    }
                }

                last;
            }
            else
            {
                &bail("Error retrieving list of compressed files from server!") if $try1 >= $MAXTRIES;
                $try1++;
                sleep $STIMER;
            }
        }

        last;

    }
    else
    {
        &bail(sprintf("[%d] %s\n", $getResponse->code, &parseError($getResponse->content)));
    }

}


# make sure outputs were generated
# decompress
print STDERR "Assembling final results...\n";

&backticks("cat $TMP/archive.tar.gz* | tar -zxf - -C $TMP/");
&backticks("mv -f $TMP/instr $INSTRFILE") if -e "$TMP/instr";
&backticks("mv -f $TMP/hat3 $HAT3FILE") if -e "$TMP/hat3";

# sometimes no hat3 file is generated [v3.1]
#&bail("Error: Output file $HAT3FILE not found!") unless -e $HAT3FILE;
&bail("Error: Output file $INSTRFILE not found!") unless -e $INSTRFILE;


# warn if some ownids were ommitted
if ( scalar keys(%ownids) > 0 )
{
    my %instrids = ();

    open(INSTRF,"<$INSTRFILE") or &bail("Error: Cannot open file $INSTRFILE for reading!");

    while(<INSTRF>)
    {
        chomp;

        if ( /^>\d+_(\S+)$/ )
        {
            $instrids{$1} = 1;
        }
    }

    close INSTRF;

    foreach my $id ( keys %ownids )
    {
        warn "Warning: Own structure $id was excluded from instr/hat3.\n" unless $instrids{$id};
    }

}



&cleanup();



####################
####################



sub parseResponse
{
    my $response = shift;

    #"status":"wait","mafftashid":"Ma8211432R"

    my $status = "";
    my $mafftashid = "";

    if ( $response =~ /^([^\s:]+):([^\s:]+)$/ )
    {
        $mafftashid = $1;
        $status = $2;
    }

    return ($status, $mafftashid);

}


sub extractchecksum
{
    my $infile = shift;
    my %dataset = ();

    open CSUM, "tar -zxf $infile -O|" or return \%dataset;

    while(<CSUM>)
    {
        chomp;
        if ( /^(\S+)\s+(\S+)$/ )
        {
            $dataset{$2} = $1;
        }

    }

    close CSUM;

    return \%dataset;

}


sub parseError
{
    my $response = shift;

    #"error":"Invalid number of inputs found."
    my $errorstr = ( $response =~ /\"error\"\s*:\s*\"([^\"]+)\"/ ) ? $1 : "";
    return $errorstr;
}


sub getchecksum
{
    my $infile = shift;

    # md5 binary check
    my $MD5BIN = "";

    if ( -x "/usr/bin/md5sum" )
    {
        $MD5BIN = "/usr/bin/md5sum";
    }
    elsif ( -x "/sbin/md5" )
    {
        $MD5BIN = "/sbin/md5 -q";
    }

    return "" if $MD5BIN eq "";


    my $checksum = "";
    open MD5EXE, "$MD5BIN $infile|" or return "";

    while(<MD5EXE>)
    {
        if (/^(\S+)\s+(\S+)$/)
        {
            $checksum = $1;
            last;
        }
        elsif (/^(\S+)$/)
        {
            $checksum = $1;
            last;
        }
    }

    close MD5EXE;

    return $checksum;

}


sub backticks
{
    my $command = shift;

    `$command`;
    return ($? == -1) ? 0 : 1;
}


sub bail
{
    my $str = shift;
    print STDERR "$str\n" if defined $str;

    &cleanup();
    exit(1);
}


sub cleanup
{
    return if ($TMP eq "" || !-d $TMP);

    opendir(MAINDIR, $TMP);
    my @files = readdir(MAINDIR);
    closedir(MAINDIR);

    foreach my $file (@files)
    {
        unlink "$TMP/$file" if -e "$TMP/$file";
    }

    remove_tree($TMP);

}


sub help
{
    my $str = shift;

    print <<'HELPME';

USAGE
  ./mafftash_premafft.pl -p [FILE]
  ./mafftash_premafft.pl -o [FILE] -d [DIRECTORY]
  ./mafftash_premafft.pl -p [FILE] -o [FILE] -d [DIRECTORY]


PARAMETERS
  -p [FILE]
     FILE contains a list of PDBIDs (one entry per line); make sure that the PDBIDs are in the standard 5-character pdbid+chain naming format

  -o [FILE] -d [DIRECTORY]
     FILE contains a list of IDs from your own structure/pdb files (one entry per line)
     for each ID in the list make sure that a corresponding structure file (same ID with .pdb extension) is stored in DIRECTORY

  -h [HATFILE]
     save the output hat3 file in HATFILE; if not set, the output is written to a file named 'hat3' in your current directory

  -i [INSTRFILE]
     save the output instr file in INSTRFILE; if not set, the output is written to a file named 'instr' in your current directory

HELPME

    &bail($str);
}



