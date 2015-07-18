use strict;
use warnings;
use Storable qw(dclone);
use Math::Trig;
no warnings ('uninitialized', 'substr', 'numeric');

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

sub parseATOM1 
{
    my($atomrecord, $startPos, $endPos, $outFileName) = @_;

#    my $filename = 'report.pdb';
    open(my $fh, '>', $outFileName) or die "Could not open file $outFileName $!";

    # Turn the scalar into an array of ATOM lines
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

#get Atom information
sub parseATOM
{
    my($atomrecord, $startPos, $endPos) = @_;
	my @results = (  );
    my(@atomrecord) = split(/\n/, $atomrecord);
    my $i = 0;
    foreach my $record (@atomrecord) 
    {
		if(substr($record, 13, 2) eq 'CA')
		{
			my $number = $i;
			my $location = substr($record, 23, 3);
			$location =~ s/^\s*//;
			if($location >= $startPos && $location <= $endPos)
			{
				$results[$number] = $record;
				$i++;
			}
		}
    }
	# Return the hash
	
    return @results;
}

sub addSpaces
{
	my($tmx1, $sign) = @_;
	my $endOfLoop;
	if($sign == 1)
	{
		$endOfLoop = 11;
	}
	if($sign == 2)
	{
		$endOfLoop = 7;
	}
	if($sign == 3)
	{
		$endOfLoop = 7;
	}
	my @line = split("", $tmx1);
	my $tmx11 = $tmx1;
	for my $i(reverse 0..$endOfLoop)
	{
		if($i > length($tmx1))
		{
			$tmx11 = " ".$tmx11;
		}
	}
	return $tmx11;
}

#-------->start the program
#get file names from user
print "Enter the file names (seperated by comma): \n";
chomp(my $names = <>);

my @names = split(',',$names);

print "Do you want to align the structures? (Enter 'y' or 'n'): ";
chomp(my $sign = <>);

#get desired length of the loop from user
print "Enter the length: \n";
chomp(my $length = <>);

my @allAtoms =();
my @allAtoms2 = ();
my @atoms = ();
my @atoms2 = ();
my $base_CA1;
my $x1;
my $y1;
my $z1;
my $base_CA2;
my $x2;
my $y2;
my $z2;
my $diffx1;my $diffy1;my $diffz1;
#my $aaPosition = " ";
my $str = " ";
my $rho;
my $theta;
my $phi;

my $CA2_x;
my $CA2_y;
my $CA2_z;

for my $protein(0..scalar(@names)-1)
{
    print $names[$protein], "\n";
    my @file = get_file_data($names[$protein]);

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
	my $current_CA1;
	my $current_CA2;
	
	for my $i(0..scalar(@sheetIntervals)-2)
	{
		my $sheet1 = substr($sheetIntervals[$i][0], index($sheetIntervals[$i][0], ' ')+1);
		my $sheet2 = substr($sheetIntervals[$i+1][0], 0, index($sheetIntervals[$i+1][0], ' '));
		if($protein == 0 && $i == 0)
		{
			$base_CA1 = $sheet1;
			$base_CA2 = $sheet2;
		}
		else
		{
			$current_CA1 = $sheet1;
			$current_CA2 = $sheet2;
		}	
		#print "diif1 = ", $base_CA1 - $current_CA1, "\ndiff2 = ", $base_CA2 - $current_CA2, "<<<<<<< interesting\n";
 		my $distance = $sheet2-$sheet1-1;
		
		if($sign eq 'n')
		{
			if($distance == $length)
			{
				my $flag = 0;
				for my $j(0..scalar(@helixIntervals)-1)
				{
					my $helixStartPos = substr($helixIntervals[$j][0], 0, index($helixIntervals[$j][0], ' ')+1);
								print $helixStartPos;
					$helixStartPos =~ s/^\s+|\s+$//g;
					print $helixStartPos, ",",$sheet1, ",",$sheet2, "\n";
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
					
					my $outFileName = substr($protein, 0, index($protein, '.'))."_".$sheet1."_".$distance."."."pdb";
					
					parseATOM1($recordtypes{'ATOM'}, $startPos, $endPos, $outFileName);
					$notpresent = 0;
				}
			}
		}
		if($sign eq 'y')
		{
		#compare the length
		if($distance == $length)
		{
			my $flag = 0;
			#check for the helix between two sheet structures
			for my $j(0..scalar(@helixIntervals)-1)
			{
				my $helixStartPos = substr($helixIntervals[$j][0], 0, index($helixIntervals[$j][0], ' ')+1);
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
				@atoms = ();
				@atoms2 = ();
				@atoms = parseATOM($recordtypes{'ATOM'}, $startPos, $endPos);
				$str = " ";
				foreach $str (@atoms)
				{
					my $aaPosition = substr($str, 22,4);
					$aaPosition =~ s/^\s*//;
					if($protein == 0 && $i == 0)
					{
						###########CA1
						if($aaPosition eq $base_CA1)
						{
							$x1       = substr($str, 31, 8);  # columns 31-38
							$y1       = substr($str, 38, 8);  # columns 39-46
							$z1       = substr($str, 46, 8);  # columns 47-54
							$x1 =~ s/^\s*//;$y1 =~ s/^\s*//;$z1 =~ s/^\s*//;
						}
						if($aaPosition eq $base_CA2)
						{
							$x2       = substr($str, 31, 8);  # columns 31-38
							$y2       = substr($str, 38, 8);  # columns 39-46
							$z2       = substr($str, 46, 8);  # columns 47-54
							$x2 =~ s/^\s*//;$y2 =~ s/^\s*//;$z2 =~ s/^\s*//;
						}
					}
					if($i != 0)
					{
						if($aaPosition eq $current_CA1)
						{
							my $mx1       = substr($str, 31, 8);  # columns 31-38
							my $my1       = substr($str, 38, 8);  # columns 39-46
							my $mz1       = substr($str, 46, 8);  # columns 47-54
							$mx1 =~ s/^\s*//;$my1 =~ s/^\s*//;$mz1 =~ s/^\s*//;
							
							$diffx1 = $mx1 - $x1;
							$diffy1 = $my1 - $y1;
							$diffz1 = $mz1 - $z1;
						}
					}
				}
				foreach my $line (@atoms)
				{
					
					if($i != 0)
					{
						my $mx1       = substr($line, 31, 8);  # columns 31-38
						my $my1       = substr($line, 38, 8);  # columns 39-46
						my $mz1       = substr($line, 46, 8);  # columns 47-54
						$mx1 =~ s/^\s*//;$my1 =~ s/^\s*//;$mz1 =~ s/^\s*//;
						
						my $tmx1 = sprintf("%.3f", $mx1 - $diffx1);
						my $tmy1 = sprintf("%.3f", $my1 - $diffy1);
						my $tmz1 = sprintf("%.3f", $mz1 - $diffz1);
						
						$tmx1 = addSpaces($tmx1, 1);
						$tmy1 = addSpaces($tmy1, 2);
						$tmz1 = addSpaces($tmz1, 3);
						
						$line       =~ substr($line, 27, 11, $tmx1);  # columns 31-38
						$line       =~ substr($line, 39, 7, $tmy1);  # columns 39-46
						$line       =~ substr($line, 47, 7, $tmz1);  # columns 47-54
						#end of transformation of CA1
						
						my $aaPosition = substr($line, 22,4);
						$aaPosition =~ s/^\s*//;
						############CA2
						
						if($aaPosition == $base_CA2)
						{
							$CA2_x       = substr($line, 31, 8);  # columns 31-38
							$CA2_y       = substr($line, 38, 8);  # columns 39-46
							$CA2_z       = substr($line, 46, 8);  # columns 47-54
							$CA2_x =~ s/^\s*//;$CA2_y =~ s/^\s*//;$CA2_z =~ s/^\s*//;
						}	
					}
				}
				#############spherical coords
				#($rho, $theta, $phi)   = cartesian_to_spherical($x2, $y2, $z2);
				
				#print $rho, "  ", $theta,"  ", $phi, "<<<<<<<<<\n";
				#############
				
				#calculate vector A same for all the iterations
				my $vectorX = $x2 - $x1;
				my $vectorY = $y2 - $y1;
				my $vectorZ = $z2 - $z1;
				
				#calculate other vector B 
				my $newVectorX = $CA2_x - $x1;
				my $newVectorY = $CA2_y - $y1;
				my $newVectorZ = $CA2_z - $z1;
				########################
				
				#my $magnitude = sqrt(($x2 - $x1)**2 + ($y2 - $y1)**2 + ($z2 - $z1)**2);
				#my $unitVectorX = $vectorX / $magnitude;
				#my $unitVectorY = $vectorY / $magnitude;
				#my $unitVectorZ = $vectorZ / $magnitude;
				
				#print $x1, $x2, $y1, $y2, $z1, $z2, "\n";
				#print $newVectorX, " ", $newVectorY, " ", $newVectorZ, " \n";
				#print $magnitude, "<<<<<<<<<<<<<<<<<<<<<<\n";
				#######################
				
				my $cosPsi = (($vectorY*$newVectorY) + ($vectorZ*$newVectorZ))/(sqrt($vectorY**2 + $vectorZ**2) * sqrt($newVectorY**2 + $newVectorZ**2));
				my $cosTheta = (($vectorX*$newVectorX) + ($vectorZ*$newVectorZ))/(sqrt($vectorX**2 + $vectorZ**2) * sqrt($newVectorX**2 + $newVectorZ**2));
				my $cosPhi = (($vectorY*$newVectorY) + ($vectorX*$newVectorX))/(sqrt($vectorY**2 + $vectorX**2) * sqrt($newVectorY**2 + $newVectorX**2));
				
				#angles
				my $psi = acos($cosPsi);
				my $theta = acos($cosTheta);
				my $phi = acos($cosPhi);
				
				#print $theta, $beta, $gaama, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
						
						#rotation matrix
						my $R11 = cos($theta)*cos($phi);
						my $R12 = (sin($psi)*sin($theta)*cos($phi))-(cos($psi)*sin($phi));
						my $R13 = (cos($psi)*sin($theta)*cos($phi)) + (sin($psi)*sin($phi));
						my $R21 = cos($theta)*sin($phi);
						my $R22 = (sin($psi)*sin($theta)*sin($phi))+(cos($psi)*cos($phi));
						my $R23 = (cos($psi)*sin($theta)*sin($phi)) - (sin($psi)*cos($phi));
						my $R31 = -sin($theta);
						my $R32 = sin($psi)*cos($theta);
						my $R33 = cos($psi)*cos($theta);
						########################
						if($i != 0)
						{
						foreach my $frag(@atoms)
						{
							my $temp_x       = substr($frag, 31, 8);  # columns 31-38
							my $temp_y       = substr($frag, 38, 8);  # columns 39-46
							my $temp_z       = substr($frag, 46, 8);  # columns 47-54
							$temp_x =~ s/^\s*//;$temp_y =~ s/^\s*//;$temp_z =~ s/^\s*//;
							
							my $tempVectorX = $temp_x - $x1;
							my $tempVectorY = $temp_y - $y1;
							my $tempVectorZ = $temp_z - $z1;
							
							#transformation
							my $transformedCA2_x = $R11*$tempVectorX + $R12*$tempVectorY + $R13*$tempVectorZ;
							my $transformedCA2_y = $R21*$tempVectorX + $R22*$tempVectorY + $R23*$tempVectorZ;
							my $transformedCA2_z = $R31*$tempVectorX + $R32*$tempVectorY + $R33*$tempVectorZ;
							########################
							
							#new coordinates
							my $transCoordinatesX = sprintf("%.3f", $transformedCA2_x + $x1);
							my $transCoordinatesY = sprintf("%.3f", $transformedCA2_y + $y1);
							my $transCoordinatesZ = sprintf("%.3f", $transformedCA2_z + $z1);	
							
							$transCoordinatesX = addSpaces($transCoordinatesX, 1);
							$transCoordinatesY = addSpaces($transCoordinatesY, 2);
							$transCoordinatesZ = addSpaces($transCoordinatesZ, 3);
							
							$frag       =~ substr($frag, 27, 11, $transCoordinatesX);  # columns 31-38
							$frag       =~ substr($frag, 39, 7, $transCoordinatesY);  # columns 39-46
							$frag       =~ substr($frag, 47, 7, $transCoordinatesZ);  # columns 47-54
						}
						}
						#print $transCoordinatesX, " ",$transCoordinatesY, " ", $transCoordinatesZ,"\n";
				
				#print $CA2_x, " ", $CA2_y, " ", $CA2_z, "\n";
				
				
				push @allAtoms2, [@atoms];
				push @allAtoms, [@atoms];
				$notpresent = 0;
			}
		}
		}
	}
	if($notpresent  == 1)
	{
		print "No such SLS present..\n";
	}
}

open(my $fh, '>>', 'coordinates.pdb') or die "Could not open file coordinates.pdb $!";
foreach my $a (@allAtoms2)
	{
		foreach my $b (@$a)
		{
			#print $b, "\n";
			print $fh $b;
			print $fh "\n";
			#print $b, "\n";
		}
	}
exit;
