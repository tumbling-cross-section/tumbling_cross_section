
char const usage[] =
"crossArea: program to find average tumbling area.\n"
"\n\n"
"  Usage:\n"
"  -gasradius    <g> (double) radius of probe molecule/atom.\n"
"  -radlib       <r> (string) filename of library containing pdb atomic radii definitions.\n"
"  -infile       <i> (string) filename containing atomic coordinates in pdb or pqr format.\n"
"  -nsteps       <n> (int)    number of tumbling orientations to sample (from a Fibbonacci grid).\n"
"\n"
"  Optionally:\n"
"  -seed     <s> (int)    to ensure repeatability by seeding the random number generator.\n"
"  -norients <o> (int)    number of random grid orientations to sample over (default 10)\n"
"  -verbose  <v>          to engage talkative mode.\n"
"  -log      <l> (string) to redirect output to a logfile\n";

