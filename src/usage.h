
char const usage[] =
"crossArea: program to find average tumbling area."
""
"  Usage:"
"  -gasradius    <g> (float)  radius of probe molecule/atom."
"  -radlib       <r> (string) filename of library containing pdb atomic radii definitions."
"  -infile       <i> (string) filename containing atomic coordinates in pdb or pqr format."
"  -anglesteps   <a> (int)    number of latitude orientations to average tumbling area over,"
"                             number of longitude orientations follows from this to give even sampling."
""
"  Optionally:"
"  -seed    <s> (int)    to ensure repeatability by seeding the random number generator."
"  -verbose <v>         to engage talkative mode."
"  -log     <l> (string) to redirect output to a logfile";

