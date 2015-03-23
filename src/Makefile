

crossArea: crossArea.c getPdbStructure.c usage.h crossArea.h getPdbStructure.h
	$(CC) -o crossArea -O2 -Wall -D _NO_LUT crossArea.c getPdbStructure.c -lm	
	cp crossArea ../.

install: crossArea
	cp crossArea /usr/bin/.

.PHONY: clean
clean:
	rm -Rf *~ *\#




