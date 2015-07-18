use strict;
use warnings;
use Storable qw(dclone);

#open file and read data
sub get_file_data 
{

    my($filename) = @_;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}

#parse all the data absed on the record types
sub parsePDBrecordtypes {

    my @file = @_;
    
    my %recordtypes = (  );
    
    foreach my $line (@file) {
    
        # Get the record type name which begins at the
        # start of the line and ends at the first space

        # The pattern (\S+) is returned and saved in $recordtype
        my($recordtype) = ($line =~ /^(\S+)/);
    
        # .= fails if a key is undefined, so we have to
        # test for definition and use either .= or = depending
        if(defined $recordtypes{$recordtype} ) {
            $recordtypes{$recordtype} .= $line;
        }else{
            $recordtypes{$recordtype} = $line;
        }
    }    
    return %recordtypes;
}

#get helix information
sub parseHELIX
{
    my($helixrecord) = @_;
    my %results = ( );

    my(@helixrecord) = split(/\n/, $helixrecord);
    my $i = 1;
    foreach my $record (@helixrecord)
    {
	my $number		= $i;
	my $startPos 	= substr($record, 22, 3);
	my $endPos 		= substr($record, 34, 3);

	$startPos =~ s/^\s*//;
	$endPos =~ s/^\s*//;

	$results{$number} = "$startPos $endPos";
	$i++;
    }
	return %results;
}

#get sheet information
sub parseSHEET
{
    my($sheetrecord) = @_;
    my %results = ( );

    my(@sheetrecord) = split(/\n/, $sheetrecord);
    my $i = 1;
    foreach my $record (@sheetrecord)
    {
        my $number 		= $i;
        my $startPos            = substr($record, 23, 3);
        my $endPos 		= substr($record, 34, 3);

        $startPos =~ s/^\s*//;
        $endPos =~ s/^\s*//;
	
		$results{$number} = "$startPos $endPos";
		$i++;
    }
	return %results;
}

#get Atom information
sub parseATOM 
{
    my($atomrecord, $startPos, $endPos, $outFileName) = @_;

    open(my $fh, '>', $outFileName) or die "Could not open file $outFileName $!";

    my(@atomrecord) = split(/\n/, $atomrecord);
    my $i = 1;
    foreach my $record (@atomrecord) 
    {
	my $location = substr($record, 23, 3);
	$location =~ s/^\s*//;
	if($location >= $startPos && $location <= $endPos)
	{
	    print $fh $record;
	    print $fh "\n";
	}
    }
    close $fh;
    print "written to the file $outFileName \n";
}


#-------->start the program
#get file name from user
print "Enter the file name: \n";
chomp(my $name = <>);

#get desired length of the loop from user
print "Enter the length: \n";
chomp(my $length = <>);

my @file = get_file_data($name);

#Parse the record types of the pdb files

my %recordtypes = parsePDBrecordtypes(@file);

#extract the atoms of all chains
my %helix = parseHELIX($recordtypes{'HELIX'});
my %sheet = parseSHEET($recordtypes{'SHEET'});

#get all sheets and sort them by their locations
my @sheetIntervals = ();


#parse all the sheet information and sort it
for my $i(1..keys %sheet)
{
    my $str = $sheet{$i};
    my @temp = split /s+/, $str;
    push @sheetIntervals, [ @temp ];
}

my $notcomplete = 1;
my $notpresent = 1;

while($notcomplete)
{
    $notcomplete = 0;
    for my $index(0..scalar(@sheetIntervals)-2)
    {
	my $temp1 = $sheetIntervals[$index][0];
	my $temp2 = $sheetIntervals[$index+1][0];
	($temp1) = ($temp1 =~ /\A(.*?) /);
	($temp2) = ($temp2 =~ /\A(.*?) /);
	if($temp1 > $temp2)
	{
	    my $temp = $sheetIntervals[$index+1][0];
	    $sheetIntervals[$index+1][0] = $sheetIntervals[$index][0];
	    $sheetIntervals[$index][0] = $temp;
	    $notcomplete = 1;
	}   
    }
}

#parse all the helix information and parse it
my @helixIntervals = ();

for my $i(1..keys %helix)
{
    my $str = $helix{$i};
    my @temp = split /s+/, $str;
    push @helixIntervals, [ @temp ];
}

$notcomplete = 1;

while($notcomplete)
{
    $notcomplete = 0;
    for my $index(0..scalar(@helixIntervals)-2)
    {
        my $temp1 = $helixIntervals[$index][0];
        my $temp2 = $helixIntervals[$index+1][0];
        ($temp1) = ($temp1 =~ /\A(.*?) /);
        ($temp2) = ($temp2 =~ /\A(.*?) /);
        if($temp1 > $temp2)
        {
            my $temp = $helixIntervals[$index+1][0];
            $helixIntervals[$index+1][0] = $helixIntervals[$index][0];
            $helixIntervals[$index][0] = $temp;
            $notcomplete = 1;
	}   
    }
}

#print all the sorted Sheet and Helix structure locations
print "=========Sheet=========\n";
for my $i(0..scalar(@sheetIntervals)-1)
{
    print $sheetIntervals[$i][0], "\n";
}

print "=========Helix=========\n";
for my $i(0..scalar(@helixIntervals)-1)
{
    print $helixIntervals[$i][0], "\n";
}

#Check for the desired length of SLS
for my $i(0..scalar(@sheetIntervals)-2)
{
    my $sheet1 = substr($sheetIntervals[$i][0], index($sheetIntervals[$i][0], ' ')+1);
    my $sheet2 = substr($sheetIntervals[$i+1][0], 0, index($sheetIntervals[$i+1][0], ' '));

    my $distance = $sheet2-$sheet1;
    
#compare the length
    if($distance == $length)
    {
		my $flag = 0;
		#check for the helix between two sheet structures
		for my $j(0..scalar(@helixIntervals)-1)
		{
			my $helixStartPos = substr($helixIntervals[$j][0], 0, index($sheetIntervals[$j][0], ' ')+1);
			$helixStartPos =~ s/^\s+|\s+$//g;
			if(($helixStartPos >= $sheet1) && ($helixStartPos <= $sheet2))
			{
				$flag = 1;
				last;
			}
			else
			{
				$flag = 0;
			}
		}
		if($flag == 0)
		{
			my $startPos = substr($sheetIntervals[$i][0], 0, index($sheetIntervals[$i][0], ' '));
			my $endPos = substr($sheetIntervals[$i+1][0], index($sheetIntervals[$i+1][0], ' ')+1);
			
			my $outFileName = substr($name, 0, index($name, '.'))."_".$sheet1."_".$distance."."."pdb";
			
			parseATOM($recordtypes{'ATOM'}, $startPos, $endPos, $outFileName);
			$notpresent = 0;
		}
    }
}
if($notpresent  == 1)
{
	print "No such SLS present..\n";
}
exit;
