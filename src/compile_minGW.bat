REM ARK: run this command from the mingW command prompt
REM ARK: eg: start->programs->MinGW->MinGW Command Prompt
REM ARK: eg: C:\crossArea> compile_minGW.bat


gcc -o crossArea_WIN32 -mno-cygwin -O2 -Wall -D _NO_LUT crossArea.c getPdbStructure.c -lm


