#!/usr/bin/perl

@program_names = ("caseq");
%program_nodes = ("caseq", 4);

$program_to_run = $ARGV[0];
$number_lines = $ARGV[1];
$number_its = $ARGV[2];

if (!$program_to_run || !$program_nodes{$program_to_run}) {
  die "Must enter program name to run. Possible programs are: " .
      "\n@program_names\n";
} else {
  if ($ENV{"MPIRUN"}) {
    $mpirun = $ENV{"MPIRUN"};
  } else {
    $mpirun = "mpirun";
  }
  if ($ENV{"MPI_HOSTS"}) {
    $hosts = "-f " . $ENV{"MPI_HOSTS"};
  } else {
    $hosts = "";
  }

  print "$mpirun -n $program_nodes{$program_to_run} $hosts ./$program_to_run $number_lines $number_its\n";
  system("$mpirun -n $program_nodes{$program_to_run} $hosts ./$program_to_run $number_lines $number_its");
}
