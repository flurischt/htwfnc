default:
	gcc-4.9 fw.c -o fw -O1 -mno-abm -fno-tree-vectorize

clean:
	rm -f ./fw

test: default
	./fw 1000
