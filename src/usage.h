
char const usage[] =
"crossArea: program to find average tumbling area.\n"
"\n\n"
"  Usage:\n"
"  -gasradius    <g> (double)  radius of probe molecule/atom.\n"
"  -radlib       <r> (string) filename of library containing pdb atomic radii definitions.\n"
"  -infile       <i> (string) filename containing atomic coordinates in pdb or pqr format.\n"
"  -anglesteps   <a> (int)    number of latitude orientations to average tumbling area over,\n"
"                             number of longitude orientations follows from this to give even sampling.\n"
"\n"
"  Optionally:\n"
"  -seed    <s> (int)    to ensure repeatability by seeding the random number generator.\n"
"  -verbose <v>          to engage talkative mode.\n"
"  -log     <l> (string) to redirect output to a logfile\n";

